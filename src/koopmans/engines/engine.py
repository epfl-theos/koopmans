import sys
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generator, Literal, Optional, overload

from upf_tools import UPFDict

from koopmans import utils
from koopmans.files import File
from koopmans.status import Status
from koopmans.processes import ProcessProtocol


class Engine(ABC):

    def __init__(self, from_scratch: bool = True, npool: Optional[int] = None):
        self.from_scratch = from_scratch
        self.npool = npool

    def print(self, text: str, **kwargs):
        utils.indented_print(text, **kwargs, wrap=False)

    def _step_message(self, uid, symbol, suffix, **kwargs):
        self.print(f'- {symbol} `{uid}` {suffix}', **kwargs)

    def _step_completed_message(self, step):
        self._step_completed_message_by_uid(step.uid)

    def _step_failed_message(self, step):
        self._step_failed_message_by_uid(step.uid)

    def _step_running_message(self, step):
        self._step_running_message_by_uid(step.uid)

    def _step_skipped_message(self, step):
        self._step_skipped_message_by_uid(step.uid)

    def _step_completed_message_by_uid(self, uid):
        self._step_message(uid, '✅', 'completed  ')

    def _step_failed_message_by_uid(self, uid):
        self._step_message(uid, '❌', 'failed     ')

    def _step_running_message_by_uid(self, uid):
        if sys.stdout.isatty():
            self._step_message(uid, '🖥️ ', 'running...', end='\r', flush=True)

    def _step_skipped_message_by_uid(self, uid):
        self._step_message(uid, '⏭️ ', 'already complete  ')

    @abstractmethod
    def run(self, step: ProcessProtocol) -> None:
        # Run a step of the workflow
        ...

    @abstractmethod
    def get_status(self, step: ProcessProtocol) -> Status:
        # Get the status of a step
        ...

    @abstractmethod
    def set_status(self, step: ProcessProtocol, status: Status):
        # Set the status of a step
        ...

    @abstractmethod
    def update_statuses(self) -> None:
        # Check the statuses of all steps and update them if necessary
        ...

    @abstractmethod
    def load_results(self, step: ProcessProtocol) -> None:
        # Load the results of a completed step
        ...

    @abstractmethod
    def get_pseudopotential(self, library: str, element: Optional[str] = None, filename: Optional[str] = None) -> UPFDict:
        # Load a pseudopotential from a library either by the species of the atom (element) or the name of the pseudopotential (filename)
        ...

    @abstractmethod
    def install_pseudopotential(self, file: Path, library: str) -> None:
        # Install a local pseudopotential file so that it can be accessed by the engine via self.get_pseudopotential()
        ...
    
    @abstractmethod
    def uninstall_pseudopotential_library(self, library: str) -> None:
        # Uninstall a pseudopotential library
        ...

    @abstractmethod
    def available_pseudo_libraries(self) -> set[str]:
        # Return a set of available pseudopotential libraries
        ...

    @overload
    @abstractmethod
    def read_file(self, file: File, binary: Literal[True]) -> bytes: ...

    @overload
    @abstractmethod
    def read_file(self, file: File, binary: Literal[False]) -> str: ...

    @overload
    @abstractmethod
    def read_file(self, file: File, binary: bool = False) -> bytes | str: ...

    @abstractmethod
    def read_file(self, file: File, binary: bool = False) -> bytes | str: ...
        # Read content from file; should mimic Path.write_text() or Path.write_bytes()

    @abstractmethod
    def write_file(self, content: str | bytes, file: File) -> None:
        # Write content to a file; should mimic Path.write_text() or Path.write_bytes()
        ...

    @abstractmethod
    def link_file(self, source: File, destination: File, recursive: bool = False, overwrite: bool = False) -> None:
        # Create a symbolic link at destination that points to source; should mimic Path.symlink_to()
        # Additionally must support recursive = True, which, if source is a directory, will link all the files
        # individually rather than creating a single link to the entire directory
        ...

    @abstractmethod
    def glob(self, directory: File, pattern: str, recursive: bool = False) -> Generator[File, None, None]:
        ...

    @abstractmethod
    def copy_file(self, source: File, destination: File, exist_ok: bool = False) -> None:
        # Remove a file; should mimic shutil.copy()
        ...
    
    @abstractmethod
    def unlink_file(self, file: File) -> None:
        # Remove a file; should mimic Path.unlink()
        ...

    @abstractmethod
    def file_exists(self, file: File) -> bool:
        # Check if a file exists; should mimic Path.exists()
        ...

    @abstractmethod
    def file_is_dir(self, file: File) -> bool:
        # Check if a file is a directory; should mimic Path.is_dir()
        ...

    @abstractmethod
    def chdir(self, directory: Path):
        # A context manager that changes the working directory (and returns to the original directory when done)
        ...

    @abstractmethod
    def mkdir(self, directory: File, parents: bool = False, exist_ok: bool = False) -> None:
        # Create a directory; should mimic Path.mkdir()
        ...

    @abstractmethod
    def rmdir(self, directory: File) -> None:
        # Remove a directory; should mimic Path.rmdir()
        ...
