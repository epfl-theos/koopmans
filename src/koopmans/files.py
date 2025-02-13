from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Generator, Union

import numpy as np

from koopmans.utils import HasDirectory

if TYPE_CHECKING:
    from koopmans.engines import Engine


class File:
    """ An abstract way of representing a file

    Because a file may not exist locally (specifically, when koopmans is run with AiiDA), we need a way of
    referring to a file that is more general than an absolute path. This class achieves this by storing a
    file as a parent (which is a Process, Calculator, or some other object that exists in a directory known
    to koopmans/AiiDA) and a name (which is the path of the file relative to the parent's directory).

    We also need to delegate file creation/modification/deletion to the engine

    """

    def __init__(self, parent: HasDirectory, name: Union[str, Path]):
        self.parent = parent
        self.name = Path(name)

    @property
    def _engine(self) -> Engine | None:
        return self.parent.engine

    def __repr__(self):
        return f'File({self.aspath()})'

    def aspath(self) -> Path:
        if self.parent.directory is None:
            return self.name
        else:
            return self.parent.directory / self.name

    def copy(self, dst: File, binary=False):
        assert self._engine is not None
        self._engine.copy_file(self, dst, binary)

    def exists(self):
        assert self._engine is not None
        return self._engine.file_exists(self)

    def read(self, binary: bool = False, numpy: bool = False) -> Any:
        assert self._engine is not None
        content = self._engine.read_file(self, binary=binary)
        if numpy:
            if binary:
                assert isinstance(content, bytes)
                return np.frombuffer(content)
            else:
                assert isinstance(content, str)
                return np.array(content)
        else:
            return content

    def rglob(self, pattern: str) -> Generator[File, None, None]:
        assert self._engine is not None
        yield from self._engine.glob(self, pattern, recursive=True)

    def glob(self, pattern: str) -> Generator[File, None, None]:
        assert self._engine is not None
        yield from self._engine.glob(self, pattern, recursive=False)

    def mkdir(self, *args, **kwargs):
        assert self._engine is not None
        self._engine.mkdir(self, *args, **kwargs)

    def is_dir(self) -> bool:
        assert self._engine is not None
        return self._engine.file_is_dir(self)

    def __eq__(self, other):
        if not isinstance(other, File):
            return False
        # Note that we only check the parent's directory, not the parent details
        return self.parent.absolute_directory == other.parent.absolute_directory and self.name == other.name

    def __reduce__(self):
        # We don't want to store the entire parent object in the database; we only need the directory information
        if self.parent is None:
            dummy_parent = None
        else:
            dummy_parent = ParentPlaceholder.fromobj(self.parent)
        return (File, (dummy_parent, self.name))

    def __gt__(self, other):
        if not isinstance(other, File):
            raise TypeError(f'Cannot compare File with {type(other)}')
        return self.aspath() > other.aspath()

    def __lt__(self, other):
        return not self > other

    def __truediv__(self, other):
        assert isinstance(other, Path) or isinstance(other, str)
        return File(self.parent, self.name / other)


class ParentPlaceholder(HasDirectory):
    # Placeholder parent for Files that don't have a Workflow/Process/Calculator as a parent OR for when we
    # don't want to store the parent in the database
    def __init__(self, parent, directory, engine):
        super().__init__(parent, directory, engine=engine, _directory_must_be_relative=False)

    def __repr__(self):
        return f'ParentPlaceholder(directory={self.absolute_directory})'

    @classmethod
    def fromobj(cls, obj, replace_parents_with_placeholders=True):
        if replace_parents_with_placeholders:
            if obj.parent is None:
                parent = None
            else:
                parent = cls.fromobj(obj.parent, replace_parents_with_placeholders=True)
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
    def frompath(cls, path: Path, engine: Engine | None = None):
        return cls(None, path, engine)


def LocalFile(path: Path | str) -> File:
    from koopmans.engines.localhost import LocalhostEngine
    path = path if isinstance(path, Path) else Path(path)
    engine = LocalhostEngine()
    if path.is_dir():
        parent = ParentPlaceholder.frompath(path.resolve(), engine=engine)
        return File(parent=parent, name='')
    else:
        parent = ParentPlaceholder.frompath(path.parent.resolve(), engine=engine)
        return File(parent=parent, name=Path(path.name))
