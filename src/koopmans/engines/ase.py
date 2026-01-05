"""An engine for running a workflow locally using the Atomic Simulation Environment."""

from __future__ import annotations

import contextlib
import shutil
from itertools import chain
from pathlib import Path
from time import sleep
from typing import (TYPE_CHECKING, Any, Generator, Generic, Literal, Optional,
                    overload)

import yaml
from ase.calculators.calculator import CalculationFailed
from ase.calculators.espresso import Espresso as ASEPWCalculator
from ase.calculators.espresso import EspressoProfile
from ase.io.formats import register_external_io_formats
from ase.spectrum.band_structure import BandStructure
from ase_pw2wannier90_plugin.calculator import \
    PW2Wannier90 as ASEPW2Wannier90Calculator
from ase_pw2wannier90_plugin.calculator import PW2Wannier90Profile
from ase_wannier90_plugin.calculator import Wannier90 as ASEWannier90Calculator
from ase_wannier90_plugin.calculator import Wannier90Profile
from upf_tools.pydantic import UPFModel

from koopmans import utils
from koopmans.calculators.calculator import CalcType
from koopmans.calculators.pw2wannier90 import (PW2Wannier90Inputs,
                                               PW2Wannier90Outputs)
from koopmans.calculators.wannier90 import Wannier90Inputs, Wannier90Outputs
from koopmans.files import File, LocalFile
from koopmans.pseudopotentials import (element_from_pseudo_filename,
                                       local_libraries)
from koopmans.status import Status

from .engine import Engine, Response

if TYPE_CHECKING:
    from koopmans.calculators.pw import PWInputs, PWOutputs
    from koopmans.processes.process import Process

config_file = Path.home() / '.koopmans/config.yaml'


class LocalASEEngine(Engine):
    """Engine for running a workflow locally using the Atomic Simulation Environment."""

    def run(self, process: Process[Any, Any], **kwargs):
        """Run a process."""
        from koopmans.calculators.calculator import CalculatorProcess

        self._process_running_message(process)

        self.processes[process.uid] = process

        assert process.directory is not None
        process_directory = LocalFile(path=process.directory)
        if self.file_exists(process_directory):
            self.rmdir(process_directory)

        if isinstance(process, CalculatorProcess):
            process._pre_run()
            self.run_calculator(process, **kwargs)
            process._post_run()
        else:
            process.run()

        self.set_status(process, Status.COMPLETED)
        self._process_completed_message(process)

        # If we reached here, all future processes should be performed from scratch
        self.from_scratch = True

        return

    def model_post_init(self, __context):
        """Perform post-initialization logic."""
        super().model_post_init(__context)
        # Register the required ASE plugins
        register_external_io_formats('ase_wannier90_plugin.formats')
        register_external_io_formats('ase_pw2wannier90_plugin.formats')

    def load_old_calculator(self, calc: CalcType):
        """Load a calculation that has already run."""
        return load_old_calculator(calc)

    def load_results(self, process: Process):
        """Load the results of a process.

        Note that for this class, process.run() also loads the results of the calculator so this function is not needed.
        """
        pass

    def _get_status(self, process: Process) -> Status:
        """Get the status of a process."""
        if not self.from_scratch:
            raise NotImplementedError(
                'Need to replace this with new functionality now that Process covers more classes')
            # to_run = True
            # if isinstance(process, ImplementedCalc):
            #     assert process.directory is not None
            #     calc_file = process.directory / process.prefix

            #     if calc_file.with_suffix(process.ext_out).is_file():
            #         loaded_process = self.load_old_calculator(process)

            #         if loaded_process.is_complete():
            #             to_run = False
            # elif isinstance(process, Process):
            #     if process.is_complete():
            #         to_run = False
            #         process.load_outputs()
            # else:
            #     raise ValueError(f'Unknown process type: {type(process)}')

            # if to_run:
            #     # This and subsequent calculations should be performed from scratch
            #     self.from_scratch = True
            # else:
            #     if process.uid not in self.statuses:
            #         self._process_skipped_message(process)
            #     self.set_status(process, Status.COMPLETED)

        if process.uid not in self.statuses:
            return Status.NOT_STARTED

        return self.statuses[process.uid]

    def _set_status(self, process: Process, status: Status):
        """Set the status of a process."""
        self.statuses[process.uid] = status

    def _update_statuses(self) -> None:
        """Update the statuses of all processes."""
        pass

    def _processes_are_running(self) -> bool:
        """Check if any processes are currently running.

        Because localhost blocks until the process is complete, there are never any processes running when python
        code is being executed. This method therefore always returns False.
        """
        return False

    def get_pseudopotential(self,
                            library: str,
                            element: Optional[str] = None,
                            filename: Optional[str] = None) -> UPFModel:
        """Return the content of a pseudopotential."""
        from koopmans.files import LocalFile
        if element is None and filename is None:
            raise ValueError('Either `element` or `filename` must be provided')

        matching_pseudo = None
        if library in local_libraries:
            pseudo_dir = LocalFile(path=Path(__file__).parents[1] / 'pseudopotentials' / library)
            if filename is not None:
                pseudo_name: str | None = filename
            elif element is not None:
                pseudo_name = None
                for pseudo_file in chain(self.glob(pseudo_dir, '*.upf'), self.glob(pseudo_dir, '*.UPF')):
                    if element_from_pseudo_filename(pseudo_file.stem) == element:
                        pseudo_name = str(pseudo_file.name)
                        break
                if pseudo_name is None:
                    raise ValueError(f'Pseudo for {element} not found in library {library}')
            else:
                raise ValueError('Either `element` or `filename` must be provided')
            matching_pseudo = pseudo_dir / pseudo_name
        else:
            try:
                with open(config_file, 'r') as f:
                    config = yaml.safe_load(f)
            except FileNotFoundError:
                raise FileNotFoundError('Could not find the pseudopotential library')

            if library not in config['installed_pseudopotentials']:
                raise FileNotFoundError(f'Could not find the custom pseudopotential library `{library}`')

            for pseudo_filepath in config['installed_pseudopotentials'][library]:
                assert pseudo_filepath is not None
                if filename is not None:
                    if Path(pseudo_filepath).name == filename:
                        matching_pseudo = LocalFile(pseudo_filepath)
                        break
                elif element is not None:
                    if element_from_pseudo_filename(Path(pseudo_filepath).name) == element:
                        matching_pseudo = LocalFile(pseudo_filepath)
                        break
            if matching_pseudo is None:
                raise FileNotFoundError(f'Could not find the requested pseudopotential in library {library}')

        if not self.file_exists(matching_pseudo):
            raise FileNotFoundError(f'Could not find the pseudopotential file `{matching_pseudo}`')

        return UPFModel.from_upf(self._file_to_path(matching_pseudo))

    def install_pseudopotential(self, file: Path, library: str) -> None:
        """Install a pseudopotential file."""
        try:
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
        except FileNotFoundError:
            config = {}

        if 'installed_pseudopotentials' not in config:
            config['installed_pseudopotentials'] = {}

        if library not in config['installed_pseudopotentials']:
            config['installed_pseudopotentials'][library] = set([str(file.resolve())])
        else:
            config['installed_pseudopotentials'][library].add(str(file.resolve()))

        with utils.chdir(config_file.parent):
            with open(config_file.name, 'w') as f:
                yaml.dump(config, f)

    def uninstall_pseudopotential_library(self, library):
        """Uninstall a pseudopotential library."""
        try:
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
        except FileNotFoundError:
            raise ValueError(f'No custom pseudopotential libraries are installed; skipping uninstallation of {library}')

        if library in config['installed_pseudopotentials']:
            del config['installed_pseudopotentials'][library]
        else:
            raise ValueError(f'Could not find the custom pseudopotential library `{library}`')

        with utils.chdir(config_file.parent):
            with open(config_file.name, 'w') as f:
                yaml.dump(config, f)

    @overload
    def read_file(self, file: File, binary: Literal[True]) -> bytes:
        ...

    @overload
    def read_file(self, file: File, binary: Literal[False]) -> str:
        ...

    @overload
    def read_file(self, file: File, binary: bool = False) -> bytes | str:
        ...

    def read_file(self, file: File, binary: bool = False) -> bytes | str:
        """Read the contents of a file."""
        assert file.parent_process.absolute_directory is not None
        full_path = file.parent_process.absolute_directory / file.relative_path
        fstring = 'rb' if binary else 'r'
        with open(full_path, fstring) as f:
            flines = f.read()
        return flines

    def _write_to_file(self, file: File, content: str | bytes) -> None:
        """Write content to a file."""
        path = self._file_to_path(file)
        fstring = 'wb' if isinstance(content, bytes) else 'w'
        if path.parents:
            parent = File(parent_process_uid=file.parent_process_uid, path=path.parent)
            self.mkdir(parent, parents=True, exist_ok=True)
        with open(path, fstring) as f:
            f.write(content)

    def copy_file(self, source: File, destination: File, exist_ok: bool = False) -> None:
        """Copy one file to another."""
        utils.copy_file(self._file_to_path(source), self._file_to_path(destination), exist_ok=exist_ok)

    def link_file(self, source: File, destination: File, recursive: bool = False, overwrite: bool = False) -> None:
        """Link one file to another."""
        source_path = self._file_to_path(source)
        destination_path = self._file_to_path(destination)
        if recursive:
            utils.symlink_tree(source_path, destination_path)
        else:
            utils.symlink(source_path, destination_path, exist_ok=overwrite, force=overwrite)

    def _file_to_path(self, file: File | LocalFile) -> Path:
        if isinstance(file, LocalFile):
            return file.path
        else:
            parent_process = self.processes[file.parent_process_uid]
            directory = parent_process.absolute_directory
            return directory / file.path

    def _file_exists(self, file: File | LocalFile) -> bool:
        return self._file_to_path(file).exists()

    def _file_is_dir(self, file: File) -> bool:
        """Return True if the file is a directory."""
        return self._file_to_path(file).is_dir()

    @contextlib.contextmanager
    def chdir(self, directory: Path):
        """Return a context manager that changes the current working directory."""
        return utils.chdir_logic(directory)

    def _rmdir(self, directory: File) -> None:
        """Remove a directory."""
        path = self._file_to_path(directory)
        if path.is_dir():
            shutil.rmtree(path)
        else:
            raise ValueError(f'{path} is not a directory')

    def _unlink_file(self, file: File) -> None:
        """Remove a file."""
        self._file_to_path(file).unlink()

    def _mkdir(self, directory: File, parents: bool = False, exist_ok: bool = False) -> None:
        """Make a directory."""
        self._file_to_path(directory).mkdir(parents=parents, exist_ok=exist_ok)

    def _glob_files(self, directory: File, pattern: str, recursive: bool = False) -> Generator[File, None, None]:
        """Generate files that match the pattern."""
        base_path = self._file_to_path(directory)
        for path in self._glob_path(base_path, pattern, recursive=recursive):
            yield File(path=path.relative_to(base_path), parent_process_uid=directory.parent_process_uid)

    def available_pseudo_libraries(self) -> set[str]:
        """List the available pseudopotential libraries."""
        try:
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
            installed_libraries = set(config['installed_pseudopotentials'].keys())
        except FileNotFoundError:
            installed_libraries = set()
        return local_libraries.union(installed_libraries)

    def wait(self, seconds: float) -> None:
        """Wait for a specified number of seconds."""
        sleep(seconds)

    def _teardown(self) -> None:
        pass

    def run_pw_calculator(self, uid: str, inputs: PWInputs, directory: Path, command: str) -> Response[PWOutputs]:
        """Run a PW calculation."""
        from koopmans.calculators.pw import PWOutputs
        from koopmans.files import File, LocalFile

        # Create an EspressoProfile
        command, _ = command.split('-in')
        profile = EspressoProfile(command=command, pseudo_dir='pseudopotentials')

        # Create an ASE PW calculator
        ase_calc = ASEPWCalculator(directory=directory, profile=profile)

        # Set up the calculator parameters
        dct = inputs.parameters.model_dump(exclude_none=True, exclude_defaults=True)
        dct['system'].pop('celldm', None)  # ASE does not support celldms
        dct['system']['ibrav'] = 0
        if dct['system'].get('nspin', 1) == 1:
            del dct['system']['tot_magnetization']
        ase_calc.parameters['input_data'] = dct

        # Set up the k-points
        if inputs.parameters.control.calculation == 'bands':
            ase_calc.parameters['kpts'] = inputs.kpoints.path
        else:
            ase_calc.parameters['kpts'] = inputs.kpoints.grid

        # Set up the pseudopotentials
        ase_calc.parameters['pseudopotentials'] = inputs.pseudopotentials

        # Copy over the pseudopotentials
        pseudo_folder = LocalFile(path=directory / 'pseudopotentials')
        self.mkdir(pseudo_folder, parents=True, exist_ok=True)
        for pseudo in set(inputs.pseudopotentials.values()):
            src = LocalFile(path=Path(__file__).parents[1]
                            / 'pseudopotentials' / inputs.pseudo_family / pseudo)
            self.link_file(src, pseudo_folder / pseudo, overwrite=True)

        ase_calc.prefix = 'espresso'

        outputs = None
        error_message = None
        error_type = None
        try:
            ase_calc.calculate(atoms=inputs.atoms, properties='energy', system_changes=None)
            status = Status.COMPLETED
        except CalculationFailed as error:
            status = Status.FAILED
            self._process_failed_message_by_id(uid)
            error_message = str(error)
            error_type = type(error)

        fermi_level = ase_calc.results.get('fermi_level', None)
        if fermi_level is None:
            fermi_level = []
        elif not isinstance(fermi_level, list):
            fermi_level = [fermi_level]

        if inputs.parameters.control.calculation == 'bands':
            band_structure = BandStructure(path=inputs.kpoints.path,
                                           energies=ase_calc.results['eigenvalues'],
                                           reference=max(fermi_level) if fermi_level else 0.0)
        else:
            band_structure = None

        if status == Status.COMPLETED:
            outdir = File(parent_process_uid=uid, path=inputs.parameters.control.outdir)
            outputs = PWOutputs(total_energy=ase_calc.results.get('energy', None),
                                band_structure=band_structure,
                                eigenvalues=ase_calc.results.get('eigenvalues', None),
                                fermi_level=fermi_level,
                                outdir=outdir)

        return Response(status=status, outputs=outputs, error_message=error_message, error_type=error_type)

    def run_wannier90_calculator(self,
                                 uid: str,
                                 inputs: Wannier90Inputs,
                                 directory: Path,
                                 command: str) -> Response[Wannier90Outputs]:
        """Run a Wannier90 calculation."""
        from koopmans.calculators.wannier90 import Wannier90Outputs

        command = command.split('PREFIX')[0]
        profile = Wannier90Profile(command=command)

        ase_calc = ASEWannier90Calculator(directory=directory, profile=profile)

        dct = inputs.parameters.model_dump(exclude_none=True, exclude_defaults=True)
        dct.update(inputs.projections.wannier90_kwargs)
        dct['mp_grid'] = inputs.kpoints.grid
        if dct.get('bands_plot', False):
            dct['kpoint_path'] = inputs.kpoints.path
        ase_calc.parameters = dct

        error_message = None
        error_type = None
        outputs = None
        try:
            ase_calc.calculate(atoms=inputs.atoms, properties=['spread'], system_changes=None)
            status = Status.COMPLETED
        except CalculationFailed as error:
            status = Status.FAILED
            self._process_failed_message_by_id(uid)
            error_message = str(error)
            error_type = type(error)

        if status == Status.COMPLETED:
            nnkp_file = ase_calc.results.pop('nnkp_file', None)
            nnkp_file_obj = File(parent_process_uid=uid, path=nnkp_file) if nnkp_file is not None else None

            outputs = Wannier90Outputs(nnkp_file=nnkp_file_obj, **ase_calc.results)

        return Response(status=status, outputs=outputs, error_message=error_message, error_type=error_type)

    def run_pw2wannier90_calculator(self,
                                    uid: str,
                                    inputs: PW2Wannier90Inputs,
                                    directory: Path,
                                    command: str) -> Response[PW2Wannier90Outputs]:
        """Run a PW2Wannier90 calculation using ASE."""
        from koopmans.calculators.pw2wannier90 import PW2Wannier90Outputs

        command = command.split('PREFIX')[0]
        profile = PW2Wannier90Profile(command=command)
        ase_calc = ASEPW2Wannier90Calculator(directory=directory, profile=profile)
        dct = inputs.parameters.model_dump(exclude_none=True, exclude_defaults=True)
        ase_calc.parameters = dct

        error_message = None
        error_type = None
        outputs = None
        try:
            ase_calc.calculate(atoms=inputs.atoms, properties=[], system_changes=None)
            status = Status.COMPLETED
        except CalculationFailed as error:
            status = Status.FAILED
            self._process_failed_message_by_id(uid)
            error_message = str(error)
            error_type = type(error)

        if status == Status.COMPLETED:
            outputs = PW2Wannier90Outputs(**ase_calc.results)

        return Response(status=status, outputs=outputs, error_message=error_message, error_type=error_type)


def load_old_calculator(calc: Generic[CalcType]) -> CalcType:
    """Load a calculation that has already run."""
    # This is a separate function so that it can be imported by other engines
    loaded_calc = calc.__class__.fromfile(calc.directory / calc.prefix)

    if loaded_calc.is_complete():
        # If it is complete, load the results
        calc.results = loaded_calc.results

        # Check the convergence of the calculation
        calc.check_convergence()

        raise NotImplementedError('Need to revisit')

        # # Load k-points if relevant
        # if hasattr(loaded_calc, 'kpts'):
        #     calc.kpts = loaded_calc.kpts

        # if isinstance(calc, ReturnsBandStructure):
        #     calc.generate_band_structure()

        # if isinstance(calc, ProjwfcCalculator):
        #     calc.generate_dos()

        # if isinstance(calc, PhCalculator):
        #     calc.read_dynG()

    return loaded_calc
