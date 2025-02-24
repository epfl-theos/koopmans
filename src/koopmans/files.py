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
    file as a parent_process (which is a Process, Calculator, or some other object that exists in a directory known
    to koopmans/AiiDA) and a name (which is the path of the file relative to the parent_process's directory).

    We also need to delegate file creation/modification/deletion to the engine

    """

    def __init__(self, parent_process: HasDirectory, name: Union[str, Path]):
        self.parent_process = parent_process
        self.name = Path(name)

    @property
    def _engine(self) -> Engine | None:
        return self.parent_process.engine

    def __repr__(self):
        return f'File({self.aspath()})'

    def aspath(self) -> Path:
        if self.parent_process.directory is None:
            return self.name
        else:
            return self.parent_process.directory / self.name

    @property  
    def parent(self) -> File:
        return File(self.parent_process, self.name.parent)

    @property
    def parents(self) -> Generator[File, None, None]:
        parent = self.parent
        while parent != self:
            yield parent
            parent = parent.parent

    def copy_to(self, dst: File, exist_ok=False):
        assert self._engine is not None
        self._engine.copy_file(self, dst, exist_ok=exist_ok)

    def exists(self):
        assert self._engine is not None
        return self._engine.file_exists(self)

    def read_text(self) -> str:
        assert self._engine is not None
        return self._engine.read_file(self, binary=False)

    def read_bytes(self) -> bytes:
        assert self._engine is not None
        return self._engine.read_file(self, binary=True)
    
    def write_text(self, content: str):
        assert self._engine is not None
        self._engine.write_file(content, self)
    
    def write_bytes(self, content: bytes):
        assert self._engine is not None
        self._engine.write_file(content, self)

    def rglob(self, pattern: str) -> Generator[File, None, None]:
        assert self._engine is not None
        yield from self._engine.glob(self, pattern, recursive=True)

    def glob(self, pattern: str) -> Generator[File, None, None]:
        assert self._engine is not None
        yield from self._engine.glob(self, pattern, recursive=False)
    
    def symlink_to(self, target: File, overwrite=False, recursive=False):
        # Create a symbolic link at self that points to target
        assert self._engine is not None
        self._engine.link_file(target, self, overwrite=overwrite, recursive=recursive)
    
    def unlink(self):
        assert self._engine is not None
        if self.is_dir():
            self._engine.rmdir(self)
        else:
            self._engine.unlink_file(self)

    def mkdir(self, *args, **kwargs):
        assert self._engine is not None
        self._engine.mkdir(self, *args, **kwargs)

    def is_dir(self) -> bool:
        assert self._engine is not None
        return self._engine.file_is_dir(self)

    def __eq__(self, other):
        if not isinstance(other, File):
            return False
        # Note that we only check the parent_process's directory, not the parent_process details
        return self.parent_process.absolute_directory == other.parent_process.absolute_directory and self.name == other.name

    def __reduce__(self):
        # We don't want to store the entire parent_process object in the database; we only need the directory information
        if self.parent is None:
            dummy_parent = None
        else:
            dummy_parent = ParentProcessPlaceholder.fromobj(self.parent_process)
        return (File, (dummy_parent, self.name))

    def __gt__(self, other):
        if not isinstance(other, File):
            raise TypeError(f'Cannot compare File with {type(other)}')
        return self.aspath() > other.aspath()

    def __lt__(self, other):
        return not self > other

    def __truediv__(self, other):
        assert isinstance(other, Path) or isinstance(other, str)
        return File(self.parent_process, self.name / other)


class ParentProcessPlaceholder(HasDirectory):
    # Placeholder parent_process for Files that don't have a Workflow/Process/Calculator as a parent_process OR for when we
    # don't want to store the parent in the database
    def __init__(self, parent_process, directory, engine, **kwargs):
        super().__init__(parent_process, directory, engine=engine, _directory_must_be_relative=False, **kwargs)

    def __repr__(self):
        return f'ParentProcessPlaceholder(directory={self.absolute_directory})'

    @classmethod
    def fromobj(cls, obj, replace_parents_with_placeholders=True):
        if replace_parents_with_placeholders:
            if obj.parent_process is None:
                parent_process = None
            else:
                parent_process = cls.fromobj(obj.parent_process, replace_parents_with_placeholders=True)
        else:
            parent_process = obj.parent_process

        if parent_process is None:
            base_directory = obj.base_directory
        else:
            base_directory = None
        new_obj = cls(parent_process, directory=obj.directory, engine=obj.engine, base_directory=base_directory)

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
        parent_process = ParentProcessPlaceholder.frompath(path.resolve(), engine=engine)
        return File(parent_process=parent_process, name='')
    else:
        parent_process = ParentProcessPlaceholder.frompath(path.parent.resolve(), engine=engine)
        return File(parent_process=parent_process, name=Path(path.name))
