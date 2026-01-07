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
from ase.spectrum.doscollection import GridDOSCollection
from ase_projwfc_plugin.calculator import Projwfc as ASEProjwfcCalculator
from ase_projwfc_plugin.calculator import ProjwfcProfile
from ase_pw2wannier90_plugin.calculator import \
    PW2Wannier90 as ASEPW2Wannier90Calculator
from ase_pw2wannier90_plugin.calculator import PW2Wannier90Profile
from ase_wannier90_plugin.calculator import Wannier90 as ASEWannier90Calculator
from ase_wannier90_plugin.calculator import Wannier90Profile
from upf_tools.pydantic import UPFModel
from wannier90_input.models.parameters import Projection

from koopmans import utils
from koopmans.calculators.calculator import CalcType
from koopmans.calculators.projwfc import ProjwfcInputs, ProjwfcOutputs
from koopmans.calculators.pw2wannier90 import (PW2Wannier90Inputs,
                                               PW2Wannier90Outputs)
from koopmans.calculators.wannier90 import (Wannier90Inputs, Wannier90Outputs,
                                            Wannier90PreprocInputs,
                                            Wannier90PreprocOutputs)
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
        register_external_io_formats('ase_projwfc_plugin.formats')

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

        if process.uid not in self.processes:
            self.processes[process.uid] = process

        if not self.from_scratch:
            if process.is_complete():
                if process.uid not in self.statuses:
                    self._process_skipped_message(process)
                process.load_outputs()
                self.set_status(process, Status.COMPLETED)
            else:
                # This and subsequent calculations should be performed from scratch
                self.from_scratch = True

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
        path = self._file_to_path(file)
        fstring = 'rb' if binary else 'r'
        with open(path, fstring) as f:
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
        elif inputs.parameters.system.nosym:
            # Make sure we follow specific ordering
            ase_calc.parameters['kpts'] = inputs.kpoints.explicit_grid
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

    def _run_ase_wannier90_calculator(self, command: str, directory: Path, inputs: Wannier90Inputs | Wannier90PreprocInputs,
                                      ase_to_output: Callable[[dict], Wannier90Outputs | Wannier90PreprocOutputs]) -> Response:
        command = command.split('PREFIX')[0]
        profile = Wannier90Profile(command=command)

        ase_calc = ASEWannier90Calculator(directory=directory, profile=profile)

        dct = inputs.parameters.model_dump(exclude_none=True, exclude_defaults=True)
        # Avoid serialization of projections
        dct['projections'] = [Projection(**p) for p in dct.get('projections', [])]

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
            outputs = ase_to_output(ase_calc.results)

        return Response(status=status, outputs=outputs, error_message=error_message, error_type=error_type)

    def run_wannier90_pp_calculator(self,
                                    uid: str,
                                    inputs: Wannier90PreprocInputs,
                                    directory: Path,
                                    command: str) -> Response[Wannier90PreprocOutputs]:
        """Run a Wannier90 preprocessing calculation."""
        from koopmans.calculators.wannier90 import Wannier90PreprocOutputs

        def ase_to_output(results: dict) -> Wannier90PreprocOutputs:
            nnkp_file = results.pop('nnkp_file')
            nnkp_file_obj = File(parent_process_uid=uid, path=nnkp_file)

            # Remove empty values
            results = {k: v for k, v in results.items() if v and k not in ['job_done', 'converged']}

            return Wannier90PreprocOutputs(nnkp_file=nnkp_file_obj, **results)

        return self._run_ase_wannier90_calculator(command=command, directory=directory, inputs=inputs,
                                                  ase_to_output=ase_to_output)

    def run_wannier90_calculator(self,
                                 uid: str,
                                 inputs: Wannier90Inputs,
                                 directory: Path,
                                 command: str) -> Response[Wannier90Outputs]:
        """Run a Wannier90 calculation."""
        from koopmans.calculators.wannier90 import Wannier90Outputs

        def ase_to_output(results: dict) -> Wannier90Outputs:
            results = {k: v for k, v in results.items() if v and k not in ['job_done', 'converged']}

            # Record the files that were generated
            exts = {
                'hr_file': '_hr.dat',
                'u_file': '_u.mat',
                'centers_file': '_centres.xyz',
                'u_dis_file': '_u_dis.dat',
                'chk_file': '.chk',
            }
            filenames = {key: 'wannier90' + ext for key, ext in exts.items()}
            files = {key: File(parent_process_uid=uid, path=path) for key, path in filenames.items()}
            files_that_exist = {key: val for key, val in files.items() if self._file_exists(val)}
            return Wannier90Outputs(**files_that_exist, **results)

        return self._run_ase_wannier90_calculator(command=command, directory=directory, inputs=inputs,
                                                  ase_to_output=ase_to_output)

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
            results = ase_calc.results
            results.pop('job_done')
            if results.pop('mmn_calculated'):
                results['mmn_file'] = File(parent_process_uid=uid, path=inputs.parameters.seedname + '.mmn')
            if results.pop('amn_calculated'):
                results['amn_file'] = File(parent_process_uid=uid, path=inputs.parameters.seedname + '.amn')
            results['eig_file'] = File(parent_process_uid=uid, path=inputs.parameters.seedname + '.eig')

            outputs = PW2Wannier90Outputs(**results)

        return Response(status=status, outputs=outputs, error_message=error_message, error_type=error_type)

    def run_projwfc_calculator(self, uid: str, inputs: ProjwfcInputs, directory: Path,
                               command: str) -> Response[ProjwfcOutputs]:
        """Run a projwfc.x calculation."""
        from koopmans.calculators.projwfc import ProjwfcOutputs

        command = command.split('PREFIX')[0]
        profile = ProjwfcProfile(command=command)
        ase_calc = ASEProjwfcCalculator(directory=directory, profile=profile)
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
            results = ase_calc.results
            is_complete = results.pop('job_done')

            pdos_dict = results.pop('pdos_dict')
            pdos_list = []
            for atom in inputs.atoms:
                # pdos filenames do not encode the principal quantum number n. In order to recover this number, we compare
                # the reported angular momentum quantum number l against the list of expected orbitals, and infer n
                # assuming only that the file corresponding to nl will come before (n+1)l
                label = atom.symbol + str(atom.tag) if atom.tag > 0 else atom.symbol
                filenames = sorted([f for f in pdos_dict if f"pdos_atm#{atom.index + 1}(" in f])
                orbitals = inputs.expected_subshells[label].copy()
                for filename in filenames:
                    ang_mtm = filename[-2]
                    orbital = [o for o in orbitals if o[1] == ang_mtm][0]
                    orbitals.remove(orbital)
                    for pdos in pdos_dict[filename]:
                        pdos.info['n'] = int(orbital[0])
                        pdos_list.append(pdos)

            outputs = ProjwfcOutputs(dos=GridDOSCollection(pdos_list), **results)

        return Response(status=status, outputs=outputs, error_message=error_message, error_type=error_type)
