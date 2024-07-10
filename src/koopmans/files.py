import os
from pathlib import Path
from typing import NamedTuple

from koopmans.utils import HasDirectoryAttr


class FilePointer(NamedTuple):
    parent: HasDirectoryAttr
    name: Path | str

    def __repr__(self):
        return f'FilePointer({os.path.relpath(self.parent.directory, Path.cwd())}/{self.name})'
