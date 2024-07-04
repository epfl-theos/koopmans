from pathlib import Path
from typing import NamedTuple

from koopmans.utils import HasDirectoryAttr


class FilePointer(NamedTuple):
    parent: HasDirectoryAttr
    name: Path
