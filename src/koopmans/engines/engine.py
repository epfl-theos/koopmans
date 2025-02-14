import sys
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generator, Literal, Optional, overload

from upf_tools import UPFDict

from koopmans import utils
from koopmans.files import File
from koopmans.status import Status
from koopmans.step import Step


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
    def run(self, step: Step) -> None:
        ...

    @abstractmethod
    def get_status(self, step: Step) -> Status:
        ...

    @abstractmethod
    def set_status(self, step: Step, status: Status):
        ...

    @abstractmethod
    def update_statuses(self) -> None:
        ...

    @abstractmethod
    def load_results(self, step: Step) -> None:
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

    @abstractmethod
    def write_file(self, content: str | bytes, file: File) -> None:
        ...

    @abstractmethod
    def link_file(self, source: File, destination: File, recursive: bool = False, overwrite: bool = False) -> None:
        ...

    @abstractmethod
    def glob(self, directory: File, pattern: str, recursive: bool = False) -> Generator[File, None, None]:
        ...

    @abstractmethod
    def copy_file(self, source: File, destination: File, exist_ok: bool = False) -> None:
        ...

    @abstractmethod
    def file_exists(self, file: File) -> bool:
        ...

    @abstractmethod
    def file_is_dir(self, file: File) -> bool:
        ...

    @abstractmethod
    def chdir(self, directory: Path):
        # A context manager that changes the working directory (and returns to the original directory when done)
        ...

    @abstractmethod
    def mkdir(self, directory: File, parents: bool = False, exist_ok: bool = False) -> None:
        ...

    @abstractmethod
    def rmdir(self, directory: File) -> None:
        ...

    @abstractmethod
    def available_pseudo_families(self) -> set[str]:
        ...
