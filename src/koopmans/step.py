from __future__ import annotations

from pathlib import Path
from typing import Protocol, runtime_checkable

from koopmans.utils import HasDirectory


@runtime_checkable
class Step(Protocol):
    # Ultimately to be merged with Process once that has a more general definition
    parent: HasDirectory | None
    name: str
    uid: str

    @property
    def directory(self) -> Path | None:
        ...

    @directory.setter
    def directory(self, value: Path | str | None) -> None:
        ...

    def run(self) -> None:
        ...

    def directory_has_been_set(self) -> bool:
        ...

    @property
    def absolute_directory(self) -> Path | None:
        ...

    @property
    def base_directory(self) -> Path | None:
        ...
