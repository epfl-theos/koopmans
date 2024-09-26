import os
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np

from koopmans.utils import (HasDirectoryInfo, get_binary_content, get_content,
                            warn, write_binary_content, write_content)


class FilePointer(NamedTuple):
    """ An abstract way of representing a file

    Because a file may not exist locally (specifically, when koopmans is run with AiiDA), we need a way of
    referring to a file that is more general than an absolute path. This class achieves this by storing a
    file as a parent (which is a Process, Calculator, or some other object that exists in a directory known
    to koopmans/AiiDA) and a name (which is the path of the file relative to the parent's directory).

    """
    parent: HasDirectoryInfo
    name: Path

    def __repr__(self):
        assert self.parent.absolute_directory is not None
        return f'FilePointer({self.parent.absolute_directory}/{self.name})'

    def aspath(self, absolute: bool = True) -> Path:
        if absolute:
            assert self.parent.absolute_directory is not None
            return self.parent.absolute_directory / self.name
        else:
            return self.name

    def copy(self, dst: Path, binary=False):
        if binary:
            binary_content = get_binary_content(self.parent, self.name)
            write_binary_content(dst, binary_content)
        else:
            content = get_content(self.parent, self.name)
            write_content(dst, content)

    def exists(self):
        return self.aspath().exists()

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

    def rglob(self, pattern: str):
        for f in self.aspath().rglob(pattern):
            assert self.parent.absolute_directory is not None
            yield FilePointer(parent=self.parent, name=f.relative_to(self.parent.absolute_directory))

    def is_dir(self):
        return self.aspath().is_dir()

    def __eq__(self, other):
        if not isinstance(other, FilePointer):
            return False
        # Note that we don't check self.parent.absolute_directory
        return self.parent.directory == other.parent.directory and self.name == other.name

    def __reduce__(self):
        # We don't want to store the entire parent object in the database; we only need the directory information
        abs_dir = self.parent.absolute_directory
        dummy_parent = ParentPlaceholder(directory=self.parent.directory,
                                         absolute_directory=self.parent.absolute_directory)
        return (FilePointer, (dummy_parent, self.name))


class ParentPlaceholder:
    # Move into the test suite?
    def __init__(self, directory, absolute_directory):
        self.directory = directory
        self._absolute_directory = absolute_directory

    @property
    def absolute_directory(self):
        # This is a property in order to follow the HasDirectoryInfo protocol
        return Path(__file__).parents[2] / self._absolute_directory


class AbsolutePath:
    """ A class that can stand in as a parent in a FilePointer when the file is unattached to a Calculator or Process """

    def __init__(self, directory: Path | str | None):
        directory = Path(directory) if isinstance(directory, str) else directory
        self.directory: Path | None = directory

    def __eq__(self, other):
        if not isinstance(other, AbsolutePath):
            return False
        return self.directory == other.directory

    @property
    def absolute_directory(self) -> Path | None:
        return self.directory


def AbsoluteFilePointer(path: Path | str) -> FilePointer:
    path = path if isinstance(path, Path) else Path(path)
    path = path.resolve()
    parent = AbsolutePath(directory=path.parent)
    return FilePointer(parent=parent, name=Path(path.name))
