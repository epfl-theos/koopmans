"""An engine for running a workflow locally."""

import contextlib
import shutil
from itertools import chain
from pathlib import Path
from typing import Generator, Literal, Optional, overload

import yaml
from ase_koopmans.calculators.calculator import CalculationFailed
from upf_tools import UPFDict

from koopmans import utils
from koopmans.calculators import (Calc, ImplementedCalc, PhCalculator,
                                  ProjwfcCalculator, ReturnsBandStructure)
from koopmans.commands import ParallelCommandWithPostfix
from koopmans.files import File, LocalFile
from koopmans.processes import Process, ProcessProtocol
from koopmans.pseudopotentials import (element_from_pseudo_filename,
                                       local_libraries)
from koopmans.status import Status

from .engine import Engine

config_file = Path.home() / '.koopmans/config.yaml'


class LocalhostEngine(Engine):
    """Engine for running a workflow locally."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.statuses = {}

    def run(self, step: ProcessProtocol):
        """Run a step."""
        self._step_running_message(step)

        assert step.directory is not None
        if step.directory.exists():
            assert isinstance(step, utils.HasDirectory)
            self.rmdir(File(step, ''))

        # Update postfix if relevant
        if self.npool and isinstance(getattr(step, 'command', None), ParallelCommandWithPostfix):
            assert hasattr(step, 'command')
            step.command.postfix = f'-npool {self.npool}'
        try:
            step.run()
        except CalculationFailed:
            self.set_status(step, Status.FAILED)
            self._step_failed_message(step)
            raise

        self.set_status(step, Status.COMPLETED)
        self._step_completed_message(step)

        # If we reached here, all future steps should be performed from scratch
        self.from_scratch = True

        return

    def load_old_calculator(self, calc: Calc):
        """Load a calculation that has already run."""
        return load_old_calculator(calc)

    def load_results(self, step: ProcessProtocol):
        """Load the results of a step.

        Note that for this class, step.run() also loads the results of the calculator so this function is not needed.
        """
        pass

    def get_status(self, step: ProcessProtocol) -> Status:
        """Get the status of a step."""
        if not self.from_scratch:
            to_run = True
            if isinstance(step, ImplementedCalc):
                assert step.directory is not None
                calc_file = step.directory / step.prefix

                if calc_file.with_suffix(step.ext_out).is_file():
                    loaded_step = self.load_old_calculator(step)

                    if loaded_step.is_complete():
                        to_run = False
            elif isinstance(step, Process):
                if step.is_complete():
                    to_run = False
                    step.load_outputs()
            else:
                raise ValueError(f'Unknown step type: {type(step)}')

            if to_run:
                # This and subsequent calculations should be performed from scratch
                self.from_scratch = True
            else:
                if step.uid not in self.statuses:
                    self._step_skipped_message(step)
                self.set_status(step, Status.COMPLETED)

        if step.uid not in self.statuses:
            return Status.NOT_STARTED

        return self.statuses[step.uid]

    def set_status(self, step: ProcessProtocol, status: Status):
        """Set the status of a step."""
        self.statuses[step.uid] = status

    def update_statuses(self) -> None:
        """Update the statuses of all steps."""
        pass

    def get_pseudopotential(self,
                            library: str,
                            element: Optional[str] = None,
                            filename: Optional[str] = None) -> UPFDict:
        """Return the content of a pseudopotential."""
        if element is None and filename is None:
            raise ValueError('Either `element` or `filename` must be provided')

        matching_pseudo = None
        if library in local_libraries:
            pseudo_dir = LocalFile(Path(__file__).parents[1] / 'pseudopotentials' / library)
            if filename is not None:
                pseudo_name: str | None = filename
            elif element is not None:
                pseudo_name = None
                for pseudo_file in chain(self.glob(pseudo_dir, '*.upf'), self.glob(pseudo_dir, '*.UPF')):
                    if element_from_pseudo_filename(pseudo_file.name.name) == element:
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

        return UPFDict.from_upf(matching_pseudo.aspath())

    def install_pseudopotential(self, file: Path, library: str) -> None:
        """Install a pseudopotential file."""
        try:
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
        except FileNotFoundError:
            config = {'installed_pseudopotentials': {}}

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
        full_path = file.parent_process.absolute_directory / file.name
        fstring = 'rb' if binary else 'r'
        with open(full_path, fstring) as f:
            flines = f.read()
        return flines

    def write_file(self, content: str | bytes, file: File) -> None:
        """Write content to a file."""
        fstring = 'wb' if isinstance(content, bytes) else 'w'
        if file.parents:
            self.mkdir(file.parent, parents=True, exist_ok=True)
        with open(file.aspath(), fstring) as f:
            f.write(content)

    def copy_file(self, source: File, destination: File, exist_ok: bool = False) -> None:
        """Copy one file to another."""
        utils.copy_file(source.aspath(), destination.aspath(), exist_ok=exist_ok)

    def link_file(self, source: File, destination: File, recursive: bool = False, overwrite: bool = False) -> None:
        """Link one file to another."""
        if recursive:
            utils.symlink_tree(source.aspath(), destination.aspath())
        else:
            utils.symlink(source.aspath(), destination.aspath(), exist_ok=overwrite, force=overwrite)

    def file_exists(self, file: File) -> bool:
        """Return True if the file exists."""
        return file.aspath().exists()

    def file_is_dir(self, file: File) -> bool:
        """Return True if the file is a directory."""
        return file.aspath().is_dir()

    @contextlib.contextmanager
    def chdir(self, directory: Path):
        """Return a context manager that changes the current working directory."""
        return utils.chdir_logic(directory)

    def rmdir(self, directory: File) -> None:
        """Remove a directory."""
        path = directory.aspath()
        if isinstance(path, str):
            path = Path(path)
        if path.is_dir():
            shutil.rmtree(path)
        else:
            raise ValueError(f'{path} is not a directory')

    def unlink_file(self, file: File) -> None:
        """Remove a file."""
        file.aspath().unlink()

    def mkdir(self, directory: File, parents: bool = False, exist_ok: bool = False) -> None:
        """Make a directory."""
        directory.aspath().mkdir(parents=parents, exist_ok=exist_ok)

    def glob(self, directory: File, pattern: str, recursive: bool = False) -> Generator[File, None, None]:
        """Generate files that match the pattern."""
        assert directory.parent_process is not None
        assert directory.parent_process.directory is not None
        if recursive:
            generator = directory.aspath().rglob(pattern)
        else:
            generator = directory.aspath().glob(pattern)
        for path in generator:
            yield File(parent_process=directory.parent_process,
                       name=path.relative_to(directory.parent_process.directory))

    def available_pseudo_libraries(self) -> set[str]:
        """List the available pseudopotential libraries."""
        try:
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
            installed_libraries = set(config['installed_pseudopotentials'].keys())
        except FileNotFoundError:
            installed_libraries = set()
        return local_libraries.union(installed_libraries)


def load_old_calculator(calc):
    """Load a calculation that has already run."""
    # This is a separate function so that it can be imported by other engines
    loaded_calc = calc.__class__.fromfile(calc.directory / calc.prefix)

    if loaded_calc.is_complete():
        # If it is complete, load the results
        calc.results = loaded_calc.results

        # Check the convergence of the calculation
        calc.check_convergence()

        # Load k-points if relevant
        if hasattr(loaded_calc, 'kpts'):
            calc.kpts = loaded_calc.kpts

        if isinstance(calc, ReturnsBandStructure):
            calc.generate_band_structure()

        if isinstance(calc, ProjwfcCalculator):
            calc.generate_dos()

        if isinstance(calc, PhCalculator):
            calc.read_dynG()

    return loaded_calc
