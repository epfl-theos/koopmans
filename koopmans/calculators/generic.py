"""

Calculator module for python_KI

Written by Edward Linscott Jan 2020

Major modifications
May 2020: replaced Extended_Espresso_cp with QE_calc, a calculator class agnostic to the underlying ASE machinery
Sep 2020: moved individual calculators into calculators/

"""

import os
import sys
import copy
import numpy as np
import ase.io as ase_io
from ase.io import espresso_kcp as kcp_io
from ase.build import make_supercell
from koopmans import defaults
import koopmans.utils as utils


class QE_calc:

    '''
    A quantum espresso calculator object that...
     - stores kcp.x/pw.x keywords in self._settings, but can be accessed like direct attributes
       e.g. self.<keyword> will return self._settings[<keyword>]
     - runs a kcp.x/pw.x calculation upon self.calculate()
     - the calculation input/output files are 'self.directory/self.name.(cp/pw)(i/o)'
     - stores the results of this calculation in self.results

    Under the hood. it uses ASE to manage I/O and executing calculations.
      self.calc -> self._ase_calc
      self.directory -> self._ase_calc.directory
      self.name -> self._ase_calc.prefix
      self.results -> self._ase_calc.results
    This could be changed in the future without affecting the rest of the code

    From this generic class we will later define
      CP_calc for calculations using kcp.x
      PW_calc for calculations using pw.x
      ... and others as listed in calculators/
    These are differentiated by self._io = kcp_io/pw_io/...

    Arguments:
        calc       an ASE calculator object to initialize the QE calculation settings
        qe_files   an alternative to calc, initialize the settings using qe input file(s)
        kwargs     any valid quantum espresso keywords to alter
    '''

    def __init__(self, calc=None, qe_files=[], skip_qc=False, **kwargs):

        # If qe_input/output_files are provided, use them instead of calc
        if len(qe_files) > 0 and calc is not None:
            raise ValueError(
                f'Please only provide either the calc or qe_files argument to {self.__class__}')

        # Interpreting the qe_files argument
        if isinstance(qe_files, str):
            # If no extension is provided, automatically read both input and output files
            if '.' not in qe_files.split('/')[-1]:
                qe_files = [qe_files + self.ext_in, qe_files + self.ext_out]
            else:
                qe_files = [qe_files]

        # Read qe input file
        for qe_file in [f for f in qe_files if self.ext_in in f]:
            calc = ase_io.read(qe_file).calc

        # Read qe output file
        for qe_file in [f for f in qe_files if self.ext_out in f]:
            if calc is None:
                calc = ase_io.read(qe_file).calc
            else:
                # Copy over the results
                calc.results = ase_io.read(qe_file).calc.results

        # Initialise the calculator object
        if isinstance(calc, QE_calc):
            self.calc = copy.deepcopy(calc.calc)
        elif calc is None:
            self.calc = self._ase_calc_class()
        else:
            self.calc = copy.deepcopy(calc)

        # Initialise a dictionary to store QE settings in
        self._settings = {}

        # Copy over settings from calc object into self._settings
        self._update_settings_dict()

        # Handle any recognised QE keywords passed as arguments
        for key, val in kwargs.items():
            if key not in self._recognised_keywords:
                raise ValueError(
                    f'{key} is not a recognised Quantum Espresso keyword')
            setattr(self, key, val)

        # Set up preprocessing flags
        self.preprocessing_flags = ['']

        # Check the corresponding program is installed
        self.check_code_is_installed()

        # Load defaults
        self.load_defaults()

        # Parse any algebraic expressions used for keywords
        self.parse_algebraic_settings()

        # Initialise quality control variables
        self.skip_qc = skip_qc
        self.results_for_qc = []

    # By default, use kcp.x
    _io = kcp_io

    # extensions for i/o files
    ext_out = ''
    ext_in = ''

    @property
    def calc(self):
        return self._ase_calc

    @calc.setter
    def calc(self, value):
        self._ase_calc = value

    @property
    def directory(self):
        return self._ase_calc.directory

    @directory.setter
    def directory(self, value):
        self._ase_calc.directory = value

    @property
    def name(self):
        return self._ase_calc.prefix

    @name.setter
    def name(self, value):
        self._ase_calc.prefix = value

    @property
    def results(self):
        return self._ase_calc.results

    @results.setter
    def results(self, value):
        self._ase_calc.results = value

    @property
    def preprocessing_flags(self):
        return self._preprocessing_flags

    @preprocessing_flags.setter
    def preprocessing_flags(self, value):
        if not hasattr(self.calc, 'original_command'):
            self.calc.original_command = self.calc.command
        if not isinstance(value, list):
            value = [value]
        self._preprocessing_flags = value

        # Updating self._ase_calc.command
        command = self.calc.original_command
        if command[:6] == 'mpirun':
            mpirun, npflag, npval, exe, after = command.split(' ', 4)
            before = [mpirun, npflag, npval, exe]
        elif command[:4] == 'srun':
            srun, exe, after = command.split(' ', 2)
            before = [srun, exe]
        else:
            before, after = command.split(' ', 1)
            before = [before]
        after = [after]
        self.calc.command = ' '.join(before + value + after)

    def calculate(self):
        # Generic function for running a calculation

        # If pseudo_dir is a relative path then make sure it accounts for self.directory
        if getattr(self, 'pseudo_dir', None) is not None and self.pseudo_dir[0] != '/':
            directory_depth = self.directory.strip('./').count('/') + 1
            self.pseudo_dir = '../' * directory_depth + self.pseudo_dir

        self._ase_calculate()

    def _ase_calculate(self):
        # Perform the calculation
        self._ase_calc.calculate()

    def _update_settings_dict(self):
        # Updates self._settings based on self._ase_calc
        self._settings = self._ase_calc.parameters

    def parse_algebraic_setting(self, expr):
        # Checks if self._settings[expr] is defined algebraically, and evaluates them
        if not isinstance(expr, str):
            return expr
        if all([c.isalpha() or c in ['_'] for c in expr]):
            return expr

        expr = expr.replace('/', ' / ').replace('*', ' * ').split()
        for i, term in enumerate(expr):
            if term in ['*', '/']:
                continue
            elif all([c.isalpha() for c in term]):
                if getattr(self, term, None) is None:
                    raise ValueError('Failed to parse ' + ''.join(map(str, expr)))
                else:
                    expr[i] = getattr(self, term)
            else:
                expr[i] = float(term)

        value = float(expr[0])
        for op, term in zip(expr[1::2], expr[2::2]):
            if op == '*':
                value *= float(term)
            elif op == '/':
                value /= float(term)
            else:
                raise ValueError('Failed to parse '
                                 ''.join([str(e) for e in expr]))

        return value

    _universal_settings_to_not_parse = set(['outdir'])

    @property
    def settings_to_not_parse(self):
        # A list of settings that typically contain mathematical operators (e.g...
        # )
        return getattr(self, '_settings_to_not_parse', self._universal_settings_to_not_parse)

    @settings_to_not_parse.setter
    def settings_to_not_parse(self, val):
        self._settings_to_not_parse = self._universal_settings_to_not_parse.union(val)

    def parse_algebraic_settings(self):
        # Checks self._settings for keywords defined algebraically, and evaluates them
        for key, value in self._settings.items():
            if key in self.settings_to_not_parse:
                continue
            self._settings[key] = self.parse_algebraic_setting(value)

    def is_converged(self):
        raise ValueError(
            f'is_converged() function has not been implemented for {self.__class__}')

    def is_complete(self):
        raise ValueError(
            f'is_complete() function has not been implemented for {self.__class__}')

    def transform_to_supercell(self, matrix, **kwargs):
        # Converts to a supercell as given by a 3x3 transformation matrix
        assert np.shape(matrix) == (3, 3)
        self._ase_calc.atoms = make_supercell(
            self._ase_calc.atoms, matrix, **kwargs)
        self._ase_calc.atoms.calc = self._ase_calc

    def check_code_is_installed(self):
        # Checks the corresponding code is installed
        command = self.calc.command
        if command[:4] == 'srun':
            i_exe = 1
        elif command[:6] == 'mpirun':
            i_exe = 3
        else:
            i_exe = 0
        executable = command.split()[i_exe]

        executable_with_path = utils.find_executable(executable)

        if executable_with_path is None:
            raise OSError(f'{executable} is not installed')

        return executable_with_path

    def write_alphas(self):
        raise NotImplementedError(f'{self.__class__}.write_alphas() has not been implemented/should not be called')

    def read_alphas(self):
        raise NotImplementedError(f'{self.__class__}.read_alphas() has not been implemented/should not be called')

    def load_defaults(self):
        calc_class = self.__class__

        if calc_class not in defaults.defaults:
            utils.warn(f'No defaults found for {calc_class} calculator')
            return

        for key, value in defaults.defaults[calc_class].items():
            if getattr(self, key, value) not in [None, value]:
                # If a setting has already been set, keep that value but print a warning
                utils.warn(
                    f'Suggested value for {key} is being overwritten; do this with caution')
            else:
                setattr(self, key, value)
        return
