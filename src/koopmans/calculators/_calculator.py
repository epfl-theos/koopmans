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
from typing import (TYPE_CHECKING, Any, Dict, Generic, List, Optional, Tuple,
                    Type, TypeVar, Union)
from uuid import UUID, uuid4

import ase_koopmans.io as ase_io
import numpy as np
from ase_koopmans import Atoms
from ase_koopmans.calculators.calculator import CalculationFailed, Calculator
from ase_koopmans.dft.kpoints import BandPath
from ase_koopmans.spectrum.band_structure import BandStructure
from numpy import typing as npt

from koopmans import settings, utils
from koopmans.files import File

if TYPE_CHECKING:
    from koopmans.workflows import Workflow


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
            raise ValueError(f'Unrecognized file format `{filenames.suffix}`')
    return sanitized_filenames


TCalc = TypeVar('TCalc', bound='CalculatorExt')
TCalcABC = TypeVar('TCalcABC', bound='CalculatorABC')


class CalculatorExt(utils.HasDirectory):

    '''
    This generic class is designed to be a parent class of a calculator that also inherits from an ASE calculator and
    CalculatorABC

    '''

    prefix: str = ''
    results: Dict[str, Any]
    ext_in: str = ''
    ext_out: str = ''
    parent_process: Workflow | None

    def __init__(self, parent_process=None, engine=None, skip_qc: bool = False, **kwargs: Any):
        super().__init__(parent_process=parent_process, engine=engine)

        # Remove arguments that should not be treated as QE keywords
        kwargs.pop('directory', None)

        # Handle any recognized QE keywords passed as arguments
        self.parameters.update(**kwargs)

        # Some calculations we don't want to check their results for when performing tests; for such calculations, set
        # skip_qc = True
        self.skip_qc = skip_qc

        # Prepare a dictionary to store a record of linked files
        self.linked_files: Dict[str, Tuple[File, bool, bool, bool]] = {}

    def __repr__(self):
        entries = []

        # prefix
        entries.append(f'prefix={self.prefix}')

        # directory
        if self.directory_has_been_set():
            entries.append(f'directory={os.path.relpath(self.directory, ".")}')

        entries.append(f'parameters={self.parameters.briefrepr()}')

        return self.__class__.__name__ + '(' + ',\n   '.join(entries) + ')'

    @property
    def parameters(self) -> settings.SettingsDict:
        if not hasattr(self, '_parameters'):
            raise ValueError(f'`{self}.parameters` has not yet been set')
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

    def run(self):
        # Alias for self.calculate so that calculators follow the Process protocol
        self.calculate()

    def _pre_run(self):
        # Alias for self._pre_calculate so that calculators follow the Process protocol
        self._pre_calculate()

    def _run(self):
        # Alias for self._calculate so that calculators follow the Process protocol
        self._calculate()

    def _post_run(self):
        # Alias for self._post_calculate so that calculators follow the Process protocol
        self._post_calculate()

    @property
    def name(self) -> str:
        # Alias for self.prefix so that calculators follow the Process protocol
        return self.prefix

    def calculate(self):
        """Generic function for running a calculator"""

        # First run any pre-calculation steps
        self._pre_calculate()

        # Then call the relevant ASE calculate() function
        self._calculate()

        # Then run any post-calculation steps
        self._post_calculate()

    def _post_calculate(self):
        """Perform any necessary post-calculation steps after running the calculation"""

        # Check if the calculation completed
        if not self.is_complete():
            raise CalculationFailed(
                f'`{self.directory}/{self.prefix}` failed; check the `Quantum ESPRESSO` output file for more details')

        # Check convergence
        self.check_convergence()

        return

    def _fetch_linked_files(self):
        """Link all files provided in self.linked_files

        This function is called in _pre_calculate() i.e. immediately before a calculation is run.
        """
        for dest_filename, (src_file, symlink, recursive_symlink, overwrite) in self.linked_files.items():
            # Convert to a File object
            dest_file = File(self, dest_filename)

            # Create the containing folder if it doesn't exist
            if dest_file.parents:
                containing_folder = dest_file.parent
                containing_folder.mkdir(parents=True, exist_ok=True)

            # Copy/link the file
            if recursive_symlink:
                assert src_file.is_dir(), 'recursive_symlink=True requires src to be a directory'
                dest_file.symlink_to(src_file, recursive=True)
            elif symlink:
                dest_file.symlink_to(src_file, overwrite=overwrite)
            else:
                src_file.copy_to(dest_file, exist_ok=overwrite)

    def _pre_calculate(self):
        """Perform any necessary pre-calculation steps before running the calculation"""

        # First, remove the directory (all files that the calculation will use must be linked, not manually
        # copied to the calculation directory)
        if self.directory.exists():
            utils.remove(self.directory)

        # By default, check the corresponding program is installed
        self.check_code_is_installed()

        # Copy over all files linked to this calculation
        self._fetch_linked_files()

        return

    def _calculate(self):
        """Run the calculation using the ASE calculator's calculate() method

        This method should NOT be overwritten by child classes. Child classes should only modify _pre_calculate() and
        _post_calculate() to perform any necessary pre- and post-calculation steps."""

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
            assert self.directory is not None
            input_file = self.directory / (self.prefix + self.ext_in)
        elif not input_file.suffix:
            # Add extension if necessary
            input_file = input_file.with_suffix(self.ext_in)

        # Load calculator from input file
        calc: Calculator = ase_io.read(input_file).calc

        # Update self based off the input file
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
                raise OSError(f'`{self.command.executable}` is not installed')
            self.command.path = executable_with_path.parent
        else:
            if not (self.command.path / self.command.executable).is_file():
                raise OSError(f'`{self.command.executable}` is not installed')
        return

    def write_alphas(self):
        raise NotImplementedError(
            f'`{self.__class__.__name__}.write_alphas()` has not been implemented/should not be called')

    def read_alphas(self):
        raise NotImplementedError(
            f'`{self.__class__.__name__}.read_alphas()` has not been implemented/should not be called')

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

    def link(self, src: File, dest_filename: File | Path | str | None = None,
             symlink: bool = False, recursive_symlink: bool = False, overwrite: bool = False):
        if dest_filename is None:
            dest_filename = src.name
        if isinstance(dest_filename, File):
            if dest_filename.parent_process != self:
                raise ValueError(
                    'When linking to a calculator, destination `File` objects must have that calculator as their `parent_process`')
            dest_filename = dest_filename.name
        self.linked_files[str(dest_filename)] = (src, symlink, recursive_symlink, overwrite)


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

    @property
    @abstractmethod
    def directory(self) -> Path | None:
        ...

    @directory.setter
    def directory(self, value: Path) -> None:
        ...

    @abstractmethod
    def is_converged(self) -> bool:
        ...

    @abstractmethod
    def is_complete(self) -> bool:
        ...

    def check_convergence(self) -> None:
        # Default behavior is to check self.is_converged(), and raise an error if this returns False. Override
        # this function if this behavior is undesired
        if not self.is_converged():
            raise CalculationFailed(f'`{self.directory}/{self.prefix}` did not converge; check the `Quantum ESPRESSO` '
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
            calc.prefix = filename.stem
            calc.directory = filename.parent
            try:
                calc.read_results()
            except Exception:
                # Calculation could not be read; must have been incomplete
                pass

        # Update calc.directory and calc.parameters.prefix
        assert hasattr(calc, 'parent_process')
        if calc.parent_process is None:
            base_directory = sanitized_filenames[0].parents[1]
        else:
            base_directory = calc.parent_process.base_directory
        calc.directory = Path(os.path.relpath(sanitized_filenames[0].parent, base_directory))
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
    def eigenvalues_from_results(self) -> npt.NDArray[np.float64]:
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
    @property
    @abstractmethod
    def from_scratch(self) -> bool:
        ...

    @property
    @abstractmethod
    def files_to_convert_with_spin2_to_spin1(self) -> Dict[str, List[File] | List[Path]]:
        ...

    @property
    @abstractmethod
    def files_to_convert_with_spin1_to_spin2(self) -> Dict[str, List[File] | List[Path]]:
        ...

    @abstractmethod
    def nspin1_dummy_calculator(self) -> CalculatorCanEnforceSpinSym:
        ...

    @abstractmethod
    def nspin1_calculator(self) -> CalculatorCanEnforceSpinSym:
        ...

    @abstractmethod
    def nspin2_dummy_calculator(self) -> CalculatorCanEnforceSpinSym:
        ...
