from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Protocol, runtime_checkable

from koopmans.utils import HasDirectory

if TYPE_CHECKING:
    from koopmans.engines import Engine


@runtime_checkable
class Step(Protocol):
    # Ultimately to be merged with Process once that has a more general definition
    parent: HasDirectory | None
    name: str

    @property
    def uid(self) -> str:
        ...

    @property
    def engine(self) -> Engine:
        ...

    @engine.setter
    def engine(self, value: Engine | None) -> None:
        ...

    @property
    def directory(self) -> Path | None:
        ...

    @directory.setter
    def directory(self, value: Path | str | None) -> None:
        ...

    def run(self) -> None:
        ...

    def _pre_run(self) -> None:
        ...

    def _run(self) -> None:
        ...

    def _post_run(self) -> None:
        ...

    def directory_has_been_set(self) -> bool:
        ...

    @property
    def absolute_directory(self) -> Path | None:
        ...

    @property
    def base_directory(self) -> Path | None:
        ...
