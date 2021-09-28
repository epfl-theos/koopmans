"""

Calculator utilities for koopmans

Written by Edward Linscott Jan 2020

Major modifications
May 2020: replaced Extended_Espresso_cp with GenericCalc, a calculator class agnostic to the underlying ASE machinery
Sep 2020: moved individual calculators into calculators/
Feb 2021: Split calculators further into GenericCalc and EspressoCalc
Sep 2021: Reshuffled files to make imports cleaner
"""

import copy
import numpy as np
from typing import Union, Optional
from pathlib import Path
import ase.io as ase_io
from ase.build import make_supercell
from ase.calculators.calculator import FileIOCalculator
from koopmans import utils, settings, pseudopotentials

# Directories of the various QE calculators
qe_parent_directory = Path(__file__).parents[2] / 'quantum_espresso'
qe_bin_directory = qe_parent_directory / 'qe_koopmans/bin/'
kcp_bin_directory = qe_parent_directory / 'cp_koopmans/bin/'


class ExtendedCalculator:

    '''
    This generic class is designed to be a parent class of a calculator that also inherits from an ASE calculator

    Arguments:
        calc       an ASE calculator object to initialize the QE calculation settings
        qe_files   an alternative to calc, initialize the settings using qe input file(s)
        kwargs     any valid quantum espresso keywords to alter
    '''

    results: dict
    ext_in: str = ''
    ext_out: str = ''

    def __init__(self, calc=None, qe_files=[], skip_qc=False, dct={}, **kwargs):
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
                try:
                    outcalc = self.read_output_file(qe_file)
                except:
                    # Calculation could not be read; must have been incomplete
                    continue
                # Copy over the results
                calc.results = outcalc.results
                if hasattr(outcalc, 'kpts'):
                    calc.kpts = outcalc.kpts

        # Initialise the calculator object
        if isinstance(calc, ExtendedCalculator):
            self.__dict__ = copy.deepcopy(calc.__dict__)
        elif isinstance(calc, FileIOCalculator):
            # We must convert from an ASE Calculator to an ExtendedCalculator
            for k, v in calc.__dict__.items():
                if k == 'parameters':
                    # Use our custom parameters class instead of ASE's Parameters
                    self.parameters.update(calc.parameters)
                else:
                    self.__dict__[k] = v
            self.atoms.calc = self
        else:
            raise ValueError(
                f'Unrecognizable object "{calc.__class__.__name__}" provided to QE_calc() as the "calc" argument')

        # If we initialised from qe_files, update self.directory and self.parameters.prefix
        if len(qe_files) > 0:
            if '/' in qe_files[0]:
                self.directory, prefix = qe_files[0].rsplit('/', 1)
            else:
                self.directory, prefix = '.', qe_files[0]
            self.parameters.prefix = prefix.rsplit('.', 1)[0]

        # Handle any recognised QE keywords passed as arguments
        self.parameters.update(**kwargs)

        # Extract nelec from the pseudos if it has not been specified explicitly
        if 'pseudopotentials' in self.parameters and 'nelec' not in self.parameters:
            self.parameters.nelec = pseudopotentials.nelec_from_pseudos(self)

        # Parse any algebraic expressions used for keywords
        self.parse_algebraic_settings()

        # Initialise quality control variables
        self.skip_qc = skip_qc
        self.results_for_qc = []
        self.qc_results = {}

    @property
    def parameters(self):
        if not hasattr(self, '_parameters'):
            raise ValueError(f'{self}.parameters has not yet been set')
        return self._parameters

    @parameters.setter
    def parameters(self, value: Union[settings.SettingsDict, dict]):
        if isinstance(value, settings.SettingsDict):
            self._parameters = value
        else:
            # If setting with a standard dictionary, retain all of the information about valid keywords etc
            self._parameters.data = {}
            self._parameters.update(**value)

    @property
    def directory(self) -> Path:
        return self._directory

    @directory.setter
    def directory(self, value: Union[Path, str]):
        if not isinstance(value, Path):
            value = Path(value)
        # Insist on directory being an absolute path
        self._directory = value.resolve()

        # Update parameters' record of self.directory
        self.parameters.directory = self._directory

    def calculate(self):
        # Generic function for running a calculation

        # First, check the corresponding program is installed
        self.check_code_is_installed()

        # If pseudo_dir is a relative path then make sure it accounts for self.directory
        if getattr(self.parameters, 'pseudo_dir', None) is not None and self.parameters['pseudo_dir'][0] != '/':
            directory_depth = self.directory.strip('./').count('/') + 1
            self.parameters.pseudo_dir = '../' * directory_depth + self.parameters.pseudo_dir

        self._ase_calculate()

    def _ase_calculate(self):
        # ASE expects self.command to be a string
        command = copy.deepcopy(self.command)
        self.command = str(command)

        # Perform the calculation
        super().calculate()

        # Restore self.command
        self.command = command

    def write_input_file(self, input_file=None):
        # By default, use ASE
        if input_file is None:
            directory = self.directory
            fname = self.parameters.prefix + self.ext_in
        else:
            directory, fname = input_file.rsplit('/', 1)
        with utils.chdir(directory):
            ase_io.write(fname, self.atoms)

    def read_input_file(self, input_file: Optional[Path] = None):
        # By default, use ASE
        if input_file is None:
            input_file = self.directory / (self.parameters.prefix + self.ext_in)
        return ase_io.read(input_file).calc

    def load_input_file(self, input_file: Optional[Path] = None):
        if input_file is None:
            input_file = self.directory / (self.parameters.prefix + self.ext_in)
        elif not input_file.suffix:
            # Add extension if necessary
            input_file = input_file.with_suffix(self.ext_in)

        # Save command
        command = self.command

        # Load calculator from input file and update self.parameters
        calc = self.read_input_file()
        self.parameters = calc.parameters
        self.atoms = calc.atoms
        self.atoms.calc = self

        # Restore command
        self.command = command

    def read_output_file(self, output_file=None):
        # By default, use ASE
        if output_file is None:
            output_file = self.directory + '/' + self.prefix + self.ext_out
        return ase_io.read(output_file).calc

    def parse_algebraic_setting(self, expr):
        # Checks if expr is defined algebraically, and evaluates them
        if not isinstance(expr, str):
            return expr
        if all([c.isalpha() or c in ['_', '"', "'"] for c in expr]):
            return expr.strip('"').strip("'")

        expr = expr.replace('/', ' / ').replace('*', ' * ').split()
        for i, term in enumerate(expr):
            if term in ['*', '/']:
                continue
            elif all([c.isalpha() for c in term]):
                if getattr(self.parameters, term, None) is None:
                    raise ValueError('Failed to parse ' + ''.join(map(str, expr)))
                else:
                    expr[i] = getattr(self.parameters, term)
            else:
                expr[i] = float(term)

        value = float(expr[0])
        for op, term in zip(expr[1::2], expr[2::2]):
            if op == '*':
                value *= float(term)
            elif op == '/':
                value /= float(term)
            else:
                raise ValueError('Failed to parse ' + ''.join([str(e) for e in expr]))

        return value

    def parse_algebraic_settings(self):
        # Checks self.parameters for keywords defined algebraically, and evaluates them
        for key in list(self.parameters.keys()):
            if key in self.parameters.to_not_parse:
                continue
            self.parameters[key] = self.parse_algebraic_setting(self.parameters[key])

    def is_converged(self):
        raise ValueError(
            f'is_converged() function has not been implemented for {self.__class__.__name__}')

    def is_complete(self):
        raise ValueError(
            f'is_complete() function has not been implemented for {self.__class__.__name__}')

    def transform_to_supercell(self, matrix, **kwargs):
        # Converts to a supercell as given by a 3x3 transformation matrix
        assert np.shape(matrix) == (3, 3)
        self.atoms = make_supercell(self.atoms, matrix, **kwargs)
        self.atoms.calc = self

    def check_code_is_installed(self):
        # Checks the corresponding code is installed
        if self.command.path == '':
            executable_with_path = utils.find_executable(self.command.executable)
            if executable_with_path is None:
                raise OSError(f'{self.command.executable} is not installed')
            self.command.path = executable_with_path.rsplit('/', 1)[0] + '/'
        else:
            assert (self.command.path / self.command.executable).is_file
        return

    def write_alphas(self):
        raise NotImplementedError(
            f'{self.__class__.__name__}.write_alphas() has not been implemented/should not be called')

    def read_alphas(self):
        raise NotImplementedError(
            f'{self.__class__.__name__}.read_alphas() has not been implemented/should not be called')

    def todict(self):
        # Shallow copy of self.__dict__
        dct = dict(self.__dict__)

        # Remove keys that we don't need to reconstruct the calculator
        for k in ['_ase_calc_class']:
            dct.pop(k, None)

        # Add additional information required by the json decoder
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    def fromdict(self, dct):
        for k, v in dct.items():
            setattr(self, k, v)


class KCWannCalculator(ExtendedCalculator):
    # Parent class for kc_ham.x, kc_screen.x and wann2kc.x calculators

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.results_for_qc = []

    def is_complete(self):
        return self.results.get('job_done', False)

    @property
    def filling(self):
        return [[True for _ in range(self.parameters.num_wann_occ)] + [False for _ in range(self.parameters.num_wann_emp)]]


kc_wann_defaults = {'outdir': 'TMP',
                    'kc_iverbosity': 1,
                    'kc_at_ks': False,
                    'homo_only': False,
                    'read_unitary_matrix': True,
                    'check_ks': True,
                    'have_empty': True,
                    'has_disentangle': True}
