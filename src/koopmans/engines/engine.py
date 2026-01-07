"""Abstract base class for any engine that wants to run a workflow."""

from __future__ import annotations

import logging
import shutil
import sys
from abc import ABC, abstractmethod
from pathlib import Path
from typing import (TYPE_CHECKING, Any, Generator, Generic, Literal, Optional,
                    TypeVar, overload)

from pydantic import Field
from upf_tools.pydantic import UPFModel

from koopmans import utils
from koopmans.commands import CommandConfig, CommandConfigs
from koopmans.files import File, LocalFile
from koopmans.protocol import ProcessProtocol as Process
from koopmans.pydantic import BaseModel
from koopmans.status import Status

if TYPE_CHECKING:
    from koopmans.calculators.calculator import CalculatorProcess
    from koopmans.calculators.projwfc import ProjwfcInputs, ProjwfcOutputs
    from koopmans.calculators.pw import PWInputs, PWOutputs
    from koopmans.calculators.pw2wannier90 import (PW2Wannier90Inputs,
                                                   PW2Wannier90Outputs)
    from koopmans.calculators.wannier90 import (Wannier90Inputs,
                                                Wannier90Outputs,
                                                Wannier90PreprocInputs,
                                                Wannier90PreprocOutputs)


class Engine(BaseModel, ABC):
    """An abstract base class for representing engines that run workflows."""

    from_scratch: bool = True
    commands: CommandConfigs = Field(default_factory=lambda: CommandConfigs())
    keep_tmpdirs: bool = True
    statuses: dict[str, Status] = Field(default_factory=dict)
    processes: dict[str, Process] = Field(default_factory=dict)
    workflows: dict[str, Process] = Field(default_factory=dict)

    def print(self, text: str, **kwargs):  # noqa: A003
        """Print text to the console."""
        utils.indented_print(text, **kwargs, wrap=False)

    def _process_message(self, uid, symbol, suffix, **kwargs):
        self.print(f'- {symbol} `{uid}` {suffix}', **kwargs)

    def _process_completed_message(self, process):
        self._process_completed_message_by_uid(process.uid)

    def _process_failed_message(self, process):
        self._process_failed_message_by_uid(process.uid)

    def _process_running_message(self, process, end: str = '\r'):
        self._process_running_message_by_uid(process.uid, end=end)

    def _process_skipped_message(self, process):
        self._process_skipped_message_by_uid(process.uid)

    def _process_completed_message_by_uid(self, uid):
        self._process_message(uid, '✅', 'completed  ')

    def _process_failed_message_by_uid(self, uid):
        self._process_message(uid, '❌', 'failed     ')

    def _process_running_message_by_uid(self, uid, end: str = '\r'):
        if sys.stdout.isatty():
            self._process_message(uid, '🖥️ ', 'running...', end=end, flush=True)

    def _process_skipped_message_by_uid(self, uid):
        self._process_message(uid, '⏭️ ', 'already complete  ')

    @abstractmethod
    def run(self, process: Process, **kwargs) -> None:
        """Run a process of the workflow.

        The command used to run the process can be modified by providing additional keyword arguments
        matching the fields of the CommandConfig for the code being run.
        """
        ...

    @abstractmethod
    def _get_status(self, process: Process) -> Status:
        """Get the status of a process."""
        ...

    def get_status(self, process: Process) -> Status:
        """Get the status of a process with logging."""
        status = self._get_status(process)
        logger = logging.getLogger(__name__)
        logger.debug(""f"Fetching status of process {process.uid}; it is {status.name}")
        return status

    @abstractmethod
    def _set_status(self, process: Process, status: Status):
        """Set the status of a process."""
        ...

    def set_status(self, process: Process, status: Status) -> None:
        """Set the status of a process with logging."""
        logger = logging.getLogger(__name__)
        logger.debug(f"Setting status of process {process.uid} to {status.name}")
        self._set_status(process, status)
        return

    @abstractmethod
    def _update_statuses(self) -> None:
        """Check the statuses of all processes and update them if necessary."""

    def update_statuses(self) -> None:
        """Check the statuses of all processes and update them if necessary with logging."""
        logger = logging.getLogger(__name__)
        logger.debug("Updating statuses of all processes")
        self._update_statuses()

    @abstractmethod
    def load_results(self, process: Process) -> None:
        """Load the results of a completed process."""
        ...

    def processes_are_running(self) -> bool:
        """Check if any processes are currently running.

        Implementations of this function must update the statuses of all processes before checking.
        """
        logger = logging.getLogger(__name__)
        are_running = self._processes_are_running()
        logger.debug(f"Checking if any processes are running: {are_running}")
        return are_running

    @abstractmethod
    def _processes_are_running(self) -> bool:
        """Check if any processes are currently running.

        Implementations of this function must update the statuses of all processes before checking.
        """
        ...

    @abstractmethod
    def get_pseudopotential(self, library: str,
                            element: Optional[str] = None,
                            filename: Optional[str] = None) -> UPFModel:
        """Load a pseudopotential.

        The pseudopotential is specified by its library and either by its species (element) or the name of the
        pseudopotential (filename).
        """
        ...

    @abstractmethod
    def install_pseudopotential(self, file: Path, library: str) -> None:
        """Install a local file so that it can be accessed by the engine via self.get_pseudopotential()."""
        ...

    @abstractmethod
    def uninstall_pseudopotential_library(self, library: str) -> None:
        """Uninstall a pseudopotential library."""
        ...

    @abstractmethod
    def available_pseudo_libraries(self) -> set[str]:
        """Return a set of available pseudopotential libraries."""
        ...

    @overload
    @abstractmethod
    def read_file(self, file: File | LocalFile, binary: Literal[True]) -> bytes:
        ...

    @overload
    @abstractmethod
    def read_file(self, file: File | LocalFile, binary: Literal[False]) -> str:
        ...

    @overload
    @abstractmethod
    def read_file(self, file: File | LocalFile, binary: bool = False) -> bytes | str:
        ...

    @abstractmethod
    def read_file(self, file: File | LocalFile, binary: bool = False) -> bytes | str:
        """Read content from file; should mimic Path.write_text() or Path.write_bytes."""
        ...

    def write_to_file(self, file: File | LocalFile, content: str | bytes) -> None:
        """Write content to a file; should mimic Path.write_text() or Path.write_bytes."""
        if isinstance(file, LocalFile):
            self._write_to_local_file(file, content)
        else:
            self._write_to_file(file, content)

        ...

    def _write_to_local_file(self, file: LocalFile, content: bytes | str, binary: bool = False) -> None:
        """Write content to a local file; should mimic Path.write_text() or Path.write_bytes."""
        fstring = 'wb' if isinstance(content, bytes) else 'w'
        self.mkdir(file.parent, parents=True, exist_ok=True)
        with open(file.path, fstring) as f:
            f.write(content)

    @abstractmethod
    def _write_to_file(self, file: File, content: str | bytes) -> None:
        """Write content to a file; should mimic Path.write_text() or Path.write_bytes."""

    @abstractmethod
    def link_file(self, source: File, destination: File, recursive: bool = False, overwrite: bool = False) -> None:
        """Create a symbolic link at destination that points to source; should mimic Path.symlink_to.

        Additionally must support recursive = True, which, if source is a directory, will link all the files
        individually rather than creating a single link to the entire directory.
        """
        ...

    def glob(self,
             directory: File | LocalFile,
             pattern: str,
             recursive: bool = False) -> Generator[File | LocalFile, None, None]:
        """Recurse through a directory and yield files that match the pattern; should mimic Path.glob."""
        if isinstance(directory, LocalFile):
            yield from self._glob_local_files(directory, pattern, recursive=recursive)
        else:
            yield from self._glob_files(directory, pattern, recursive=recursive)

    def _glob_local_files(self,
                          directory: LocalFile,
                          pattern: str,
                          recursive: bool = False) -> Generator[LocalFile, None, None]:
        for path in self._glob_path(directory.path, pattern, recursive=recursive):
            yield LocalFile(path=path)

    @abstractmethod
    def _glob_files(self, directory: File, pattern: str, recursive: bool = False) -> Generator[File, None, None]:
        ...

    def _glob_path(self, directory: Path, pattern: str, recursive: bool = False) -> Generator[Path, None, None]:
        if recursive:
            yield from directory.rglob(pattern)
        else:
            yield from directory.glob(pattern)

    @abstractmethod
    def copy_file(self, source: File, destination: File, exist_ok: bool = False) -> None:
        """Copy a file; should mimic shutil.copy."""
        ...

    @abstractmethod
    def _unlink_file(self, file: File) -> None:
        """Remove a file; should mimic Path.unlink."""
        ...

    def unlink_file(self, file: File | LocalFile) -> None:
        """Remove a file; should mimic Path.unlink."""
        if isinstance(file, LocalFile):
            file.path.unlink()
        else:
            self._unlink_file(file)
        return

    def file_exists(self, file: File | LocalFile) -> bool:
        """Check if a file exists; should mimic Path.exists."""
        if isinstance(file, LocalFile):
            return file.path.exists()
        else:
            return self._file_exists(file)

    @abstractmethod
    def _file_exists(self, file: File) -> bool:
        """Check if a file exists; should mimic Path.exists."""
        ...

    def file_is_dir(self, file: File | LocalFile) -> bool:
        """Check if a file is a directory; should mimic Path.is_dir."""
        if isinstance(file, LocalFile):
            return file.path.is_dir()
        else:
            return self._file_is_dir(file)

    @abstractmethod
    def _file_is_dir(self, file: File) -> bool:
        """Check if a file is a directory; should mimic Path.is_dir."""
        ...

    @abstractmethod
    def chdir(self, directory: Path):
        """Return a context manager that changes directory (returning to the original directory when finished)."""
        ...

    @abstractmethod
    def _mkdir(self, directory: File, parents: bool = False, exist_ok: bool = False) -> None:
        """Create a directory; should mimic Path.mkdir."""
        ...

    def mkdir(self, directory: File | LocalFile, parents: bool = False, exist_ok: bool = False) -> None:
        """Create a directory; should mimic Path.mkdir."""
        if isinstance(directory, LocalFile):
            directory.path.mkdir(parents=parents, exist_ok=exist_ok)
        else:
            self._mkdir(directory, parents=parents, exist_ok=exist_ok)
        return

    @abstractmethod
    def _rmdir(self, directory: File) -> None:
        """Remove a directory."""

    def rmdir(self, directory: File | LocalFile) -> None:
        """Remove a directory."""
        if isinstance(directory, LocalFile):
            shutil.rmtree(directory.path)
        else:
            self._rmdir(directory)
        return

    def rename_file(self, src: File | LocalFile, dst: File | LocalFile) -> None:
        """Rename a file; should mimic Path.rename."""
        src.rename(dst)
        raise NotImplementedError('Need to check behaviour of this')

    @abstractmethod
    def wait(self, seconds: float) -> None:
        """Wait for a specified amount of time, possibly performing background tasks while doing so."""
        ...

    @abstractmethod
    def _teardown(self) -> None:
        """Perform any necessary engine-specific cleanup actions when the workflow is completed."""
        ...

    def remove_tmpdirs(self):
        """Remove temporary directories created during the workflow."""
        for process in self.processes:
            outdir_path = get(process, 'parameters.control.outdir', None)
            if outdir_path is None:
                continue
            outdir = File(parent_process=process, relative_path=outdir_path)
            if outdir.exists():
                outdir.unlink()

    def run_calculator(self, calc: CalculatorProcess[Any, Any], **kwargs):
        """Run a calculator.

        The command used to run the calculator can be modified by providing additional keyword arguments
        """
        from koopmans.calculators.projwfc import ProjwfcCalculator
        from koopmans.calculators.pw import PWCalculator
        from koopmans.calculators.pw2wannier90 import PW2Wannier90Calculator
        from koopmans.calculators.wannier90 import (Wannier90Calculator,
                                                    Wannier90PreprocCalculator)

        if isinstance(calc, PWCalculator):
            caller = self.run_pw_calculator
        elif isinstance(calc, Wannier90PreprocCalculator):
            kwargs['pp'] = True
            caller = self.run_wannier90_pp_calculator
        elif isinstance(calc, Wannier90Calculator):
            caller = self.run_wannier90_calculator
        elif isinstance(calc, PW2Wannier90Calculator):
            caller = self.run_pw2wannier90_calculator
        elif isinstance(calc, ProjwfcCalculator):
            caller = self.run_projwfc_calculator
        else:
            raise NotImplementedError(
                f'Engine {self.__class__.__name__} does not implement running calculator of type '
                f'{calc.__class__.__name__}')

        # Determine the command to run the process based on the configuration of the code in question
        command_config: CommandConfig | None = getattr(self.commands, calc.code, None)
        if command_config is None:
            raise ValueError(f'The configuration for {calc.code} is missing from {self.__class__.__name__}.codes')

        # Update using the kwargs
        this_command_config = command_config.model_copy()
        for k, v in kwargs.items():
            setattr(this_command_config, k, v)
        command = this_command_config.command()

        response = caller(calc.uid, calc.inputs, calc.directory, command)

        if response.status == Status.COMPLETED:
            calc.outputs = response.outputs

        return

    @abstractmethod
    def run_pw_calculator(self, uid: str, inputs: PWInputs, directory: Path, command: str) -> Response[PWOutputs]:
        """Run a pw.x calculation."""
        ...

    @abstractmethod
    def run_wannier90_calculator(self, uid: str, inputs: Wannier90Inputs, directory: Path,
                                 command: str) -> Response[Wannier90Outputs]:
        """Run a wannier90.x calculation."""
        ...

    @abstractmethod
    def run_wannier90_pp_calculator(self, uid: str, inputs: Wannier90PreprocInputs, directory: Path,
                                    command: str) -> Response[Wannier90PreprocOutputs]:
        """Run a wannier90.x -pp calculation."""
        ...

    @abstractmethod
    def run_pw2wannier90_calculator(self, uid: str, inputs: PW2Wannier90Inputs, directory: Path,
                                    command: str) -> Response[PW2Wannier90Outputs]:
        """Run a pw2wannier90.x calculation."""
        ...

    @abstractmethod
    def run_projwfc_calculator(self, uid: str, inputs: ProjwfcInputs, directory: Path,
                               command: str) -> Response[ProjwfcOutputs]:
        """Run a projwfc.x calculation."""
        ...


OutputModel = TypeVar('OutputModel', bound=BaseModel)


class Response(BaseModel, Generic[OutputModel]):
    """Response of a run_calculator call."""

    status: Status
    outputs: OutputModel | None
    error_message: str | None = None
    error_type: type | None = None
