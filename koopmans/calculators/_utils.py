"""

Calculator utilities for koopmans

The central objects defined in this submodule are CalculatorABC and CalculatorExt, designed to extend
ASE calculators to have several additional useful features.

We can create a new 'extended' version of a preexisting ASE calculator via
    class ExtendedCalc(CalculatorExt, ASECalc, CalculatorABC):
        pass

Written by Edward Linscott Jan 2020

Major modifications
May 2020: replaced Extended_Espresso_cp with GenericCalc, a calculator class agnostic to the underlying ASE machinery
Sep 2020: moved individual calculators into calculators/
Feb 2021: Split calculators further into GenericCalc and EspressoCalc
Sep 2021: Reshuffled files to make imports cleaner
"""

import copy
from typing import Union, Optional, List, TypeVar, Generic, Type
from pathlib import Path
from abc import ABC, abstractclassmethod, abstractmethod, abstractproperty
from ase import Atoms
import ase.io as ase_io
from ase.calculators.calculator import FileIOCalculator
from koopmans import utils, settings, pseudopotentials

# Directories of the various QE calculators
qe_parent_directory = Path(__file__).parents[2] / 'quantum_espresso'
qe_bin_directory = qe_parent_directory / 'qe_koopmans/bin/'
kcp_bin_directory = qe_parent_directory / 'cp_koopmans/bin/'


def sanitise_filenames(filenames: Union[str, Path, List[str], List[Path]], ext_in: str, ext_out: str) -> List[Path]:
    # Generic function for sanitising the input of CalculatorExt.fromfile()
    if isinstance(filenames, List):
        sanitised_filenames = [Path(f) for f in filenames]
    else:
        if isinstance(filenames, str):
            filenames = Path(filenames)
        # If the input is a single string...
        if filenames.suffix in [ext_in, ext_out]:
            # ... and it has a valid suffix, convert it to a list and proceed
            sanitised_filenames = [filenames]
        elif filenames.suffix == '':
            # ... and it has no suffix, automatically add the expected suffixes for both the input and output files
            sanitised_filenames = [filenames.with_suffix(ext_in), filenames.with_suffix(ext_out)]
        else:
            raise ValueError(f'Unrecognised file format {filenames.suffix}')
    return sanitised_filenames


T = TypeVar('T', bound='CalculatorExt')


class CalculatorABC(ABC, Generic[T]):

    '''
    This abstract base class defines various functions we expect any Calculator to possess

    '''

    ext_in: str
    ext_out: str

    def __init__(self, atoms: Atoms) -> None:
        self.prefix: str = ''
        pass

    @abstractmethod
    def read_input(self, input_file: Optional[Path] = None):
        ...

    @abstractmethod
    def read_results(self) -> None:
        ...

    @abstractproperty
    def directory(self) -> Path:
        ...

    @directory.setter
    def directory(self, value: Union[Path, str]) -> None:
        ...

    @abstractmethod
    def is_converged(self) -> bool:
        ...

    @abstractmethod
    def is_complete(self) -> bool:
        ...

    @abstractmethod
    def todict(self) -> dict:
        ...

    @abstractclassmethod
    def fromdict(cls, dct: dict) -> T:
        ...

    @classmethod
    def fromfile(cls, filenames: Union[str, Path, List[str], List[Path]]):
        sanitised_filenames = sanitise_filenames(filenames, cls.ext_in, cls.ext_out)

        # Initialise a new calc object
        calc = cls(atoms=Atoms())

        # Read qe input file
        for filename in [f for f in sanitised_filenames if f.suffix == cls.ext_in]:
            calc.read_input(input_file=filename)

        # Read qe output file
        for filename in [f for f in sanitised_filenames if f.suffix == cls.ext_out]:
            calc.directory = filename.parent
            calc.prefix = filename.stem
            try:
                calc.read_results()
            except:
                # Calculation could not be read; must have been incomplete
                pass

        # Update calc.directory and calc.parameters.prefix
        calc.directory = sanitised_filenames[0].parent
        calc.prefix = sanitised_filenames[0].stem

        # Return the new calc object
        return calc


class CalculatorExt():

    '''
    This generic class is designed to be a parent class of a calculator that also inherits from an ASE calculator and
    CalculatorABC

    '''

    prefix: str = ''
    results: dict
    ext_in: str = ''
    ext_out: str = ''

    def __init__(self, skip_qc=False, **kwargs):
        # Handle any recognised QE keywords passed as arguments
        self.parameters.update(**kwargs)

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
    def parameters(self, value: Union[settings.SettingsDict, dict, None]):
        if isinstance(value, settings.SettingsDict):
            self._parameters = value
        else:
            # If setting with a standard dictionary or None, retain all of the information about valid keywords etc
            self._parameters.data = {}
            if value is not None:
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

        # Then call the relevant ASE calculate() function
        self._ase_calculate()

    def _ase_calculate(self):
        # ASE expects self.command to be a string
        command = copy.deepcopy(self.command)
        self.command = str(command)

        # Perform the calculation
        super().calculate()

        # Restore self.command
        self.command = command

    def read_input(self, input_file: Optional[Path] = None):
        # Auto-generate the appropriate input file name if required
        if input_file is None:
            input_file = self.directory / (self.prefix + self.ext_in)
        elif not input_file.suffix:
            # Add extension if necessary
            input_file = input_file.with_suffix(self.ext_in)

        # Load calculator from input file
        calc = ase_io.read(input_file).calc

        # Update self based off the input file
        self.parameters = calc.parameters
        if calc.atoms is not None:
            # Some calculators (e.g. wann2kc) can't reconstruct atoms from an input file
            self.atoms = calc.atoms
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

    @classmethod
    def fromdict(cls, dct: dict) -> 'CalculatorExt':
        calc = cls(dct.pop('atoms'))
        for k, v in dct.items():
            setattr(calc, k.lstrip('_'), v)
        return calc


class KCWannCalculator(CalculatorExt):
    # Parent class for kc_ham.x, kc_screen.x and wann2kc.x calculators

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.results_for_qc = []

    def is_complete(self):
        return self.results.get('job_done', False)

    @property
    def filling(self):
        return [[True for _ in range(self.parameters.num_wann_occ)]
                + [False for _ in range(self.parameters.num_wann_emp)]]
