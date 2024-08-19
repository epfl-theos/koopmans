import os
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np

from koopmans.utils import (HasDirectoryAttr, get_binary_content, get_content,
                            write_binary_content, write_content)


class FilePointer(NamedTuple):
    """ An abstract way of representing a file

    Because a file may not exist locally (specifically, when koopmans is run with AiiDA), we need a way of
    referring to a file that is more general than an absolute path. This class achieves this by storing a
    file as a parent (which is a Process, Calculator, or some other object that exists in a directory known
    to koopmans/AiiDA) and a name (which is the path of the file relative to the parent's directory).

    """
    parent: HasDirectoryAttr
    name: Path

    def __repr__(self):
        assert self.parent.directory is not None
        relpath: str = os.path.relpath(self.parent.directory, Path.cwd())
        return f'FilePointer({relpath}/{self.name})'

    def aspath(self):
        assert self.parent.directory is not None
        return self.parent.directory / self.name

    def copy(self, dst: Path, binary=False):
        if binary:
            binary_content = get_binary_content(self.parent, self.name)
            write_binary_content(dst, binary_content)
        else:
            content = get_content(self.parent, self.name)
            write_content(dst, content)

    def read(self, binary: bool = False, numpy: bool = False) -> Any:
        if binary:
            binary_content = get_binary_content(self.parent, self.name)
            if numpy:
                return np.frombuffer(binary_content)
            else:
                return binary_content
        else:
            content = get_content(self.parent, self.name)
            if numpy:
                return np.array(content)
            else:
                return content

    def __eq__(self, other):
        if not isinstance(other, FilePointer):
            return False
        return self.parent.directory == other.parent.directory and self.name == other.name


class AbsolutePath:
    """ A class that can stand in as a parent in a FilePointer when the file is unattached to a Calculator or Process """

    def __init__(self, directory: Path | str | None):
        directory = Path(directory) if isinstance(directory, str) else directory
        self.directory: Path | None = directory

    def __eq__(self, other):
        if not isinstance(other, AbsolutePath):
            return False
        return self.directory == other.directory


def AbsoluteFilePointer(path: Path | str) -> FilePointer:
    path = path if isinstance(path, Path) else Path(path)
    path = path.resolve()
    parent = AbsolutePath(directory=path.parent)
    return FilePointer(parent=parent, name=Path(path.name))
