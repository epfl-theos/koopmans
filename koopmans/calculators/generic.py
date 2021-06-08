"""

Calculator module for python_KI

Written by Edward Linscott Jan 2020

Major modifications
May 2020: replaced Extended_Espresso_cp with GenericCalc, a calculator class agnostic to the underlying ASE machinery
Sep 2020: moved individual calculators into calculators/
Feb 2021: Split calculators further into GenericCalc and EspressoCalc
"""

import os
import sys
import copy
import numpy as np
import ase.io as ase_io
from ase.io.espresso import koopmans_cp as kcp_io
from ase.build import make_supercell
from ase.calculators.calculator import FileIOCalculator
import koopmans.utils as utils


class GenericCalc:

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

    # By default, use kcp.x
    _io = kcp_io

    # extensions for i/o files
    ext_out = ''
    ext_in = ''

    _valid_settings = None
    _settings_that_are_paths = []

    def __init__(self, calc=None, qe_files=[], skip_qc=False, dct={}, **kwargs):

        # Initialise a dictionary to store QE settings in
        self._settings = {}

        # Create settings shortcuts
        self.create_settings_shortcuts()

        # Construct from dct if this is provided
        if dct:
            self.fromdict(dct)
            return

        # If qe_input/output_files are provided, use them instead of calc
        if len(qe_files) > 0 and calc is not None:
            raise ValueError(f'Please only provide either the calc or qe_files argument to {self.__class__}')

        # Interpreting the qe_files argument
        if isinstance(qe_files, str):
            # If no extension is provided, automatically read both input and output files
            if '.' not in qe_files.split('/')[-1]:
                qe_files = [qe_files + self.ext_in, qe_files + self.ext_out]
            else:
                qe_files = [qe_files]

        # Read qe input file
        for qe_file in [f for f in qe_files if self.ext_in in f]:
            calc = self.read_input_file(qe_file)

        # Read qe output file
        for qe_file in [f for f in qe_files if self.ext_out in f]:
            if calc is None:
                calc = self.read_output_file(qe_file)
            else:
                # Copy over the results
                calc.results = self.read_output_file(qe_file).results

        # Initialise the calculator object
        if isinstance(calc, GenericCalc):
            self.calc = copy.deepcopy(calc.calc)
        elif calc is None:
            self.calc = self._ase_calc_class()
        elif isinstance(calc, FileIOCalculator):
            self.calc = copy.deepcopy(calc)
        else:
            raise ValueError(
                f'Unrecognizable object "{calc.__class__.__name__}" provided to QE_calc() as the "calc" argument')

        # Copy over settings from calc object into self._settings
        self._update_settings_dict()

        # If we initialised from qe_files, update self.directory and self.name
        if len(qe_files) > 0:
            self.directory, name = qe_files[0].rsplit('/', 1)
            self.name = name.rsplit('.', 1)[0]

        # Handle any recognised QE keywords passed as arguments
        for key, val in kwargs.items():
            if key not in self._valid_settings:
                raise ValueError(f'{key} is not a recognised keyword for {self.__class__.__name__}')
            setattr(self, key, val)

        # Load defaults
        self.load_defaults()

        # Parse any algebraic expressions used for keywords
        self.parse_algebraic_settings()

        # Initialise quality control variables
        self.skip_qc = skip_qc
        self.results_for_qc = []
        self.qc_results = {}

    def create_settings_shortcuts(self):
        # Dynamically add keywords as decorated properties of the child class being initialised. This means one can
        # set and get keywords as self.<keyword> but internally they are stored as self._settings['keyword'] rather
        # than self.<keyword>
        assert self._valid_settings is not None, \
            f'self._valid_settings has not been defined for {self.__class__.__name__}'

        for k in self._valid_settings:
            if hasattr(self.__class__, k):
                continue

            # We need to use these make_get/set functions so that get/set_k are
            # evaluated immediately (otherwise we run into late binding and 'k'
            # is not defined when get/set_k are called)

            def make_get_and_set(key):
                def get_k(self):
                    # Return 'None' rather than an error if the keyword has not
                    # been defined
                    return self._settings.get(key, None)

                if k in self._settings_that_are_paths:
                    # Insist on paths being stored as absolute paths
                    def set_k(self, value):
                        if value is None:
                            self._settings[key] = None
                        elif value.startswith('/'):
                            self._settings[key] = value
                        else:
                            self._settings[key] = os.path.abspath(self.directory + '/' + value)
                else:
                    def set_k(self, value):
                        self._settings[key] = value
                return get_k, set_k

            get_k, set_k = make_get_and_set(k)

            # Add the keyword to self.__class__ (not self)
            setattr(self.__class__, k, property(get_k, set_k))

    @property
    def calc(self):
        return self._ase_calc

    @calc.setter
    def calc(self, value):
        self._ase_calc = value

    @property
    def directory(self):
        # Ensure directory is an absolute path
        if not self.calc.directory.startswith('/'):
            self.calc.directory = os.path.abspath(self.calc.directory)
        return self.calc.directory

    @directory.setter
    def directory(self, value):
        # Insist on directory being an absolute path
        self.calc.directory = os.path.abspath(value)

    @property
    def name(self):
        return self._ase_calc.prefix

    @name.setter
    def name(self, value):
        self._ase_calc.prefix = value

    @property
    def results(self):
        return self.calc.results

    @results.setter
    def results(self, value):
        self.calc.results = value

    @property
    def settings(self):
        return self._settings

    def calculate(self):
        # Generic function for running a calculation

        # First, check the corresponding program is installed
        self.check_code_is_installed()

        # If pseudo_dir is a relative path then make sure it accounts for self.directory
        if getattr(self, 'pseudo_dir', None) is not None and self.pseudo_dir[0] != '/':
            directory_depth = self.directory.strip('./').count('/') + 1
            self.pseudo_dir = '../' * directory_depth + self.pseudo_dir

        self._ase_calculate()

    def _ase_calculate(self):
        # ASE expects self.command to be a string
        command = copy.deepcopy(self._ase_calc.command)
        self._ase_calc.command = str(command)

        # Perform the calculation
        self._ase_calc.calculate()

        # Restore self._ase_calc.command
        self._ase_calc.command = command

    def _update_settings_dict(self):
        # Points self._settings to self.calc.parameters
        self._settings = self.calc.parameters

    def write_input_file(self, input_file=None):
        # By default, use ASE
        if input_file is None:
            directory = self.directory
            fname = self.name + self.ext_in
        else:
            directory, fname = input_file.rsplit('/', 1)
        with utils.chdir(directory):
            ase_io.write(fname, self.calc.atoms)

    def read_input_file(self, input_file=None):
        # By default, use ASE
        if input_file is None:
            input_file = self.directory + '/' + self.name + self.ext_in
        return ase_io.read(input_file).calc

    def load_input_file(self, input_file=None):
        if input_file is None:
            input_file = self.directory + '/' + self.name + self.ext_in
        elif '.' not in input_file.split('/')[-1]:
            # Add extension if necessary
            input_file += self.ext_in

        # Save directory, name, and command
        directory = self.directory
        name = self.name
        command = self.calc.command

        # Load calculator from input file and update self._settings
        self.calc = self.read_input_file()
        self._update_settings_dict()

        # Restore directory, name, and command
        self.directory = directory
        self.name = name
        self.calc.command = command

    def read_output_file(self, output_file=None):
        # By default, use ASE
        if output_file is None:
            output_file = self.directory + '/' + self.name + self.ext_out
        return ase_io.read(output_file).calc

    def parse_algebraic_setting(self, expr):
        # Checks if self._settings[expr] is defined algebraically, and evaluates them
        if not isinstance(expr, str):
            return expr
        if all([c.isalpha() or c in ['_', '"', "'"] for c in expr]):
            return expr.strip('"').strip("'")

        expr = expr.replace('/', ' / ').replace('*', ' * ').split()
        for i, term in enumerate(expr):
            if term in ['*', '/']:
                continue
            elif all([c.isalpha() for c in term]):
                if getattr(self, term, None) is None:
                    raise ValueError('Failed to parse ' +
                                     ''.join(map(str, expr)))
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
        default = self._universal_settings_to_not_parse.union(self._settings_that_are_paths)
        return getattr(self, '_settings_to_not_parse', default)

    @settings_to_not_parse.setter
    def settings_to_not_parse(self, val):
        default = self._universal_settings_to_not_parse.union(self._settings_that_are_paths)
        self._settings_to_not_parse = default.union(val)

    def parse_algebraic_settings(self):
        # Checks self._settings for keywords defined algebraically, and evaluates them
        for key, value in self._settings.items():
            if key in self.settings_to_not_parse:
                continue
            setattr(self, key, self.parse_algebraic_setting(value))

    def is_converged(self):
        raise ValueError(
            f'is_converged() function has not been implemented for {self.__class__.__name__}')

    def is_complete(self):
        raise ValueError(
            f'is_complete() function has not been implemented for {self.__class__.__name__}')

    def transform_to_supercell(self, matrix, **kwargs):
        # Converts to a supercell as given by a 3x3 transformation matrix
        assert np.shape(matrix) == (3, 3)
        self._ase_calc.atoms = make_supercell(
            self._ase_calc.atoms, matrix, **kwargs)
        self._ase_calc.atoms.calc = self._ase_calc

    def check_code_is_installed(self):
        # Checks the corresponding code is installed
        if self.calc.command.path is '':
            executable_with_path = utils.find_executable(self.calc.command.executable)
            if executable_with_path is None:
                raise OSError(f'{self.calc.command.executable} is not installed')
            self.calc.command.path = executable_with_path.rsplit('/', 1)[0] + '/'
        else:
            assert os.path.isfile(self.calc.command.path + self.calc.command.executable)
        return

    def write_alphas(self):
        raise NotImplementedError(
            f'{self.__class__.__name__}.write_alphas() has not been implemented/should not be called')

    def read_alphas(self):
        raise NotImplementedError(
            f'{self.__class__.__name__}.read_alphas() has not been implemented/should not be called')

    @property
    def defaults(self):
        raise NotImplementedError(f'{self.__class__.__name__} has not got defined defaults')

    def load_defaults(self):
        for key, value in self.defaults.items():
            if getattr(self, key, None) is None:
                setattr(self, key, value)
        return

    def todict(self):
        dct = self.__dict__
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    def fromdict(self, dct):
        for k, v in dct.items():
            setattr(self, k, v)


class EspressoCalc(GenericCalc):

    def __init__(self, *args, **kwargs):
        self._valid_settings = [k for sublist in self._io.KEYS.values() for k in sublist]
        super().__init__(*args, **kwargs)

    @property
    def calc(self):
        # First, update the param block
        self._ase_calc.parameters['input_data'] = self.construct_namelist()
        return self._ase_calc

    @calc.setter
    def calc(self, value):
        self._ase_calc = value

    def _ase_calculate(self):
        # Before running the calculation, update the keywords for the ASE calculator object
        self._ase_calc.parameters['input_data'] = self.construct_namelist()
        super()._ase_calculate()

    def construct_namelist(self):
        # Returns a namelist of settings, grouped by their Quantum Espresso headings
        return self._io.construct_namelist(**self._settings, warn=True)

    def _update_settings_dict(self):
        # Updates self._settings based on self._ase_calc
        self._settings = {}
        for namelist in self._ase_calc.parameters.get('input_data', {}).values():
            for key, val in namelist.items():
                setattr(self, key, val)


class KCWannCalc(EspressoCalc):
    # Parent class for kc_ham.x, kc_screen.x and wann2kc.x calculators
    defaults = {'outdir': 'TMP',
                'kc_iverbosity': 1,
                'kc_at_ks': False,
                'homo_only': False,
                'read_unitary_matrix': True,
                'check_ks': True,
                'have_empty': True,
                'has_disentangle': True}

    _settings_that_are_paths = ['outdir']

    def __init__(self, *args, **kwargs):
        self.settings_to_not_parse = ['assume_isolated']

        super().__init__(*args, **kwargs)

        self.results_for_qc = []

    @property
    def name(self):
        return self._ase_calc.prefix

    @name.setter
    def name(self, value):
        self._ase_calc.prefix = value
        self.prefix = value

    def is_complete(self):
        return self.calc.results.get('job_done', False)

    @property
    def filling(self):
        return [[True for _ in range(self.num_wann_occ)] + [False for _ in range(self.num_wann_emp)]]
