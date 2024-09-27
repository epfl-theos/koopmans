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
        dummy_parent = ParentPlaceholder.fromobj(self.parent)
        return (FilePointer, (dummy_parent, self.name))


class ParentPlaceholder:
    # Placeholder parent for FilePointers that don't have a Workflow/Process/Calculator as a parent
    def __init__(self, parent, directory, base_directory=None):
        self.parent = parent
        self.directory = directory
        if self.parent is None:
            self.base_directory = base_directory
        else:
            self._base_directory = None

    def __repr__(self):
        return f'ParentPlaceholder(directory={self.absolute_directory})'

    @classmethod
    def fromobj(cls, obj, replace_parents_with_placeholders=False):
        if replace_parents_with_placeholders:
            if obj.parent is None:
                parent = None
            else:
                parent = cls.fromobj(obj, replace_parents_with_placeholders=True)
        else:
            parent = obj.parent

        if parent is None:
            base_directory = obj.base_directory
        else:
            base_directory = None
        directory = obj.directory
        new_obj = cls(parent, directory, base_directory)

        # Sanity checks
        assert new_obj.directory == obj.directory
        assert new_obj.absolute_directory == obj.absolute_directory
        assert new_obj.base_directory == obj.base_directory

        return new_obj

    @classmethod
    def frompath(cls, path: Path):
        return cls(None, Path(), path)

    @property
    def base_directory(self) -> Path:
        if self.parent is not None:
            return self.parent.base_directory
        else:
            if self._base_directory is None:
                raise ValueError(f'{self.__class__.__name__}.base_directory has not been set')
            return self._base_directory

    @base_directory.setter
    def base_directory(self, value: Path):
        if self.parent is not None:
            raise ValueError(f'{self.__class__.__name__}.base_directory should not be set for processes with parents')
        self._base_directory = value.resolve()

    @property
    def absolute_directory(self) -> Path:
        assert self.directory is not None
        if self.parent is None:
            abs_dir = Path(self.base_directory / self.directory)
            assert abs_dir.is_absolute()
            return abs_dir
        path = self.parent.directory / self.directory

        # Recursive through the parents, adding their directories to path (these are all relative paths)...
        obj = self.parent
        while getattr(obj, 'parent', None):
            assert obj.parent is not None
            path = obj.parent.directory / path
            obj = obj.parent

        # Finally, 'path' is relative to self.base_directory
        return self.base_directory / path


def AbsoluteFilePointer(path: Path | str) -> FilePointer:
    path = path if isinstance(path, Path) else Path(path)
    parent = ParentPlaceholder.frompath(path.parent)
    return FilePointer(parent=parent, name=Path(path.name))
