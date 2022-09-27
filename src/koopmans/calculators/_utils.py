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

from __future__ import annotations

import copy
import os
from abc import ABC, abstractmethod, abstractproperty
from pathlib import Path
from typing import Any, Dict, Generic, List, Optional, Type, TypeVar, Union

import ase.io as ase_io
import numpy as np
from ase import Atoms
from ase.calculators.calculator import CalculationFailed, Calculator
from ase.dft.kpoints import BandPath
from ase.spectrum.band_structure import BandStructure
from numpy import typing as npt

from koopmans import settings, utils


def sanitize_filenames(filenames: Union[str, Path, List[str], List[Path]], ext_in: str, ext_out: str) -> List[Path]:
    # Generic function for sanitizing the input of CalculatorExt.fromfile()
    if isinstance(filenames, List):
        sanitized_filenames = [Path(f) for f in filenames]
    else:
        if isinstance(filenames, str):
            filenames = Path(filenames)
        # If the input is a single string...
        if filenames.suffix in [ext_in, ext_out]:
            # ... and it has a valid suffix, convert it to a list and proceed
            sanitized_filenames = [filenames]
        elif filenames.suffix == '':
            # ... and it has no suffix, automatically add the expected suffixes for both the input and output files
            sanitized_filenames = [filenames.with_suffix(ext_in), filenames.with_suffix(ext_out)]
        else:
            raise ValueError(f'Unrecognized file format {filenames.suffix}')
    return sanitized_filenames


TCalc = TypeVar('TCalc', bound='CalculatorExt')
TCalcABC = TypeVar('TCalcABC', bound='CalculatorABC')


class CalculatorExt():

    '''
    This generic class is designed to be a parent class of a calculator that also inherits from an ASE calculator and
    CalculatorABC

    '''

    prefix: str = ''
    results: Dict[str, Any]
    ext_in: str = ''
    ext_out: str = ''

    def __init__(self, skip_qc: bool = False, **kwargs: Any):
        # Remove arguments that should not be treated as QE keywords
        kwargs.pop('directory', None)

        # Handle any recognized QE keywords passed as arguments
        self.parameters.update(**kwargs)

        # Some calculations we don't want to check their results for when performing tests; for such calculations, set
        # skip_qc = True
        self.skip_qc = skip_qc

    def __repr__(self):
        entries = []

        # prefix
        entries.append(f'prefix={self.prefix}')

        # directory
        entries.append(f'directory={os.path.relpath(self.directory, ".")}')

        entries.append(f'parameters={self.parameters.briefrepr()}')

        return self.__class__.__name__ + '(' + ',\n   '.join(entries) + ')'

    @property
    def parameters(self) -> settings.SettingsDict:
        if not hasattr(self, '_parameters'):
            raise ValueError(f'{self}.parameters has not yet been set')
        return self._parameters

    @parameters.setter
    def parameters(self, value: Union[settings.SettingsDict, Dict[str, Any], None]):
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
        self._calculate()

        # Then check if the calculation completed
        if not self.is_complete():

            raise CalculationFailed(
                f'{self.directory}/{self.prefix} failed; check the Quantum ESPRESSO output file for more details')

        # Then check convergence
        self.check_convergence()

    def _calculate(self):
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
        calc: Calculator = ase_io.read(input_file).calc

        # Update self based off the input file, first updating self.directory in order to ensure any settings that are
        # relative paths are appropriately stored
        self.directory = input_file.parent
        self.parameters = calc.parameters
        if isinstance(calc.atoms, Atoms):
            # Some calculators (e.g. wann2kc) can't reconstruct atoms from an input file
            self.atoms = calc.atoms
            self.atoms.calc = self

    def check_code_is_installed(self):
        # Checks the corresponding code is installed
        if self.command.path == Path():
            executable_with_path = utils.find_executable(self.command.executable)
            if executable_with_path is None:
                raise OSError(f'{self.command.executable} is not installed')
            self.command.path = executable_with_path.parent
        else:
            if not (self.command.path / self.command.executable).is_file():
                raise OSError(f'{self.command.executable} is not installed')
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
    def fromdict(cls: Type[TCalc], dct: Any) -> TCalc:
        calc = cls(dct.pop('atoms'))
        for k, v in dct.items():
            setattr(calc, k.lstrip('_'), v)
        return calc


class CalculatorABC(ABC, Generic[TCalc]):

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

    def check_convergence(self) -> None:
        # Default behaviour is to check self.is_converged(), and raise an error if this returns False. Override
        # this function if this behaviour is undesired
        if not self.is_converged():
            raise CalculationFailed(f'{self.directory}/{self.prefix} did not converge; check the Quantum ESPRESSO '
                                    'output file for more details')

    @abstractmethod
    def todict(self) -> Dict[str, Any]:
        ...

    @classmethod
    @abstractmethod
    def fromdict(cls: Type[TCalcABC], dct: Dict[str, Any]) -> TCalc:
        ...

    @classmethod
    def fromfile(cls, filenames: Union[str, Path, List[str], List[Path]]):
        sanitized_filenames = sanitize_filenames(filenames, cls.ext_in, cls.ext_out)

        # Initialize a new calc object
        calc = cls(atoms=Atoms())

        # Read qe input file
        for filename in [f for f in sanitized_filenames if f.suffix == cls.ext_in]:
            calc.read_input(input_file=filename)

        # Read qe output file
        for filename in [f for f in sanitized_filenames if f.suffix == cls.ext_out]:
            calc.directory = filename.parent
            calc.prefix = filename.stem
            try:
                calc.read_results()
            except Exception:
                # Calculation could not be read; must have been incomplete
                pass

        # Update calc.directory and calc.parameters.prefix
        calc.directory = sanitized_filenames[0].parent
        calc.prefix = sanitized_filenames[0].stem

        # Return the new calc object
        return calc


class ReturnsBandStructure(ABC):
    """
    Abstract base class to be used for calculators that return bandstructures. These classes implement a
    self.generate_band_structure() which is called after self.calculate(). This is done after self.calculate()
    (and not during) because we require access to the band path

    Putting the band structure in self.results is very un-ASE-y, so we might want to ultimately align all of this
    with the more general self.band_structure() method of ASE
    """

    @abstractmethod
    def eigenvalues_from_results(self) -> npt.NDArray[np.float_]:
        ...

    @abstractmethod
    def vbm_energy(self) -> float:
        ...

    def generate_band_structure(self):
        if isinstance(self.parameters.kpts, BandPath):
            path = self.parameters.kpts
            if len(path.kpts) > 1:
                # Fetch bandstructure from results
                eigenvalues_np = self.eigenvalues_from_results()
                self.results['band structure'] = BandStructure(path, eigenvalues_np, reference=self.vbm_energy())
        return


class KCWannCalculator(CalculatorExt):
    # Parent class for KCWHam, KCWScreen and Wann2KCW calculators
    def is_complete(self) -> bool:
        return self.results.get('job_done', False)

    @property
    def filling(self):
        return [[True for _ in range(self.parameters.num_wann_occ)]
                + [False for _ in range(self.parameters.num_wann_emp)]]


class CalculatorCanEnforceSpinSym(ABC):
    # Abstract base class for calculators that can run a sequence of calculations in order to enforce spin symmetry
    # (with the goal of avoiding spin contamination)
    @abstractproperty
    def from_scratch(self) -> bool:
        ...

    @abstractmethod
    def convert_wavefunction_2to1(self):
        ...

    @abstractmethod
    def nspin1_dummy_calculator(self) -> CalculatorCanEnforceSpinSym:
        ...

    @abstractmethod
    def nspin1_calculator(self) -> CalculatorCanEnforceSpinSym:
        ...

    @abstractmethod
    def convert_wavefunction_1to2(self):
        ...

    @abstractmethod
    def nspin2_dummy_calculator(self) -> CalculatorCanEnforceSpinSym:
        ...

    @abstractmethod
    def prepare_to_read_nspin1(self):
        ...
