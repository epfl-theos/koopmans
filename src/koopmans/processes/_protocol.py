from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Protocol, runtime_checkable

from koopmans.utils import HasDirectory

if TYPE_CHECKING:
    from koopmans.engines import Engine


@runtime_checkable
class ProcessProtocol(Protocol):
    """Protocol for a Process.

    Ultimately to be merged with Process once that has a more general definition
    """

    parent_process: HasDirectory | None
    name: str
    engine: Engine | None

    @property
    def uid(self) -> str:
        """Return the unique identifier for the process."""
        ...

    @property
    def directory(self) -> Path | None:
        """Return the directory in which the process runs."""
        ...

    @directory.setter
    def directory(self, value: Path | str | None) -> None:
        ...

    def run(self) -> None:
        """Run the Process."""
        ...

    def _pre_run(self) -> None:
        ...

    def _run(self) -> None:
        ...

    def _post_run(self) -> None:
        ...

    def directory_has_been_set(self) -> bool:
        """Return true if the directory has been set."""
        ...

    @property
    def absolute_directory(self) -> Path | None:
        """Return the absolute path to the directory in which the process is contained."""
        ...

    @property
    def base_directory(self) -> Path | None:
        """Return the base directory in which the process is contained."""
        ...
