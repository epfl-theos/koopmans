"""A class for representing files in a general way by tethering them to a parent process."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Generator, Union

from koopmans.utils import HasDirectory

if TYPE_CHECKING:
    from koopmans.engines import Engine


class File:
    """An abstract way of representing a file.

    Because a file may not exist locally (specifically, when koopmans is run with AiiDA), we need a way of
    referring to a file that is more general than an absolute path. This class achieves this by storing a
    file as a parent_process (which is a Process, Calculator, or some other object that exists in a directory known
    to koopmans/AiiDA) and a name (which is the path of the file relative to the parent_process's directory).

    We also need to delegate file creation/modification/deletion to the engine.
    """

    def __init__(self, parent_process: HasDirectory, name: Union[str, Path]):
        self.parent_process = parent_process
        self.name = Path(name)

    @property
    def _engine(self) -> Engine | None:
        return self.parent_process.engine

    def __repr__(self):
        return f'File({self.aspath()})'

    def __str__(self):
        return str(self.aspath())

    def aspath(self) -> Path:
        """Return the File as a Path object."""
        if self.parent_process.directory is None:
            return self.name
        else:
            return self.parent_process.directory / self.name

    @property
    def parent(self) -> File:
        """Return the parent directory."""
        return File(self.parent_process, self.name.parent)

    @property
    def parents(self) -> Generator[File, None, None]:
        """Return a generator of all parent directories."""
        parent = self.parent
        while parent != self:
            yield parent
            parent = parent.parent

    @property
    def suffix(self) -> str:
        """Return the file suffix, replicating the behavior of Path.suffix."""
        return self.name.suffix

    def copy_to(self, dst: File, exist_ok=False):
        """Copy this file to another."""
        assert self._engine is not None
        logger = logging.getLogger(__name__)
        logger.info(f'Copying {self.aspath()} to {dst.aspath()}')
        self._engine.copy_file(self, dst, exist_ok=exist_ok)

    def exists(self) -> bool:
        """Return true if the file exists."""
        assert self._engine is not None
        return self._engine.file_exists(self)

    def read_text(self) -> str:
        """Read text from this file."""
        assert self._engine is not None
        logger = logging.getLogger(__name__)
        logger.info(f'Reading text from {self.aspath()}')
        return self._engine.read_file(self, binary=False)

    def read_bytes(self) -> bytes:
        """Read bytes from this file."""
        assert self._engine is not None
        logger = logging.getLogger(__name__)
        logger.info(f'Reading bytes from {self.aspath()}')
        return self._engine.read_file(self, binary=True)

    def write_text(self, content: str):
        """Write text to this file."""
        assert self._engine is not None
        logger = logging.getLogger(__name__)
        logger.info(f'Writing text to {self.aspath()}')
        self._engine.write_file(content, self)

    def write_bytes(self, content: bytes):
        """Write bytes to this file."""
        assert self._engine is not None
        logger = logging.getLogger(__name__)
        logger.info(f'Writing bytes to {self.aspath()}')
        self._engine.write_file(content, self)

    def rglob(self, pattern: str) -> Generator[File, None, None]:
        """Iterate over all files within this directory and its subdirectories that match the pattern."""
        assert self._engine is not None
        yield from self._engine.glob(self, pattern, recursive=True)

    def glob(self, pattern: str) -> Generator[File, None, None]:
        """Iterate over all files within this directory that match the pattern."""
        assert self._engine is not None
        yield from self._engine.glob(self, pattern, recursive=False)

    def symlink_to(self, target: File, overwrite=False, recursive=False):
        """Symbolically link this file to another file."""
        # Create a symbolic link at self that points to target
        assert self._engine is not None
        logger = logging.getLogger(__name__)
        logger.info(f'Creating symlink from {self.aspath()} to {target.aspath()}')
        self._engine.link_file(target, self, overwrite=overwrite, recursive=recursive)

    def unlink(self):
        """Remove this file/directory."""
        assert self._engine is not None
        if self.is_dir():
            self._engine.rmdir(self)
        else:
            self._engine.unlink_file(self)
        logger = logging.getLogger(__name__)
        logger.info(f'Removing {self.aspath()}')

    def mkdir(self, *args, **kwargs):
        """Create this directory."""
        assert self._engine is not None
        logger = logging.getLogger(__name__)
        logger.info(f'Creating directory {self.aspath()}')
        self._engine.mkdir(self, *args, **kwargs)

    def is_dir(self) -> bool:
        """Return true if the file is a directory."""
        assert self._engine is not None
        return self._engine.file_is_dir(self)

    def __eq__(self, other):
        if not isinstance(other, File):
            return False
        # Note that we only check the parent_process's directory, not the parent_process details
        return self.parent_process.absolute_directory == other.parent_process.absolute_directory \
            and self.name == other.name

    def __reduce__(self):
        # We don't want to store the entire parent_process object in the database; we only need the directory
        # information
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

    def with_suffix(self, suffix: str) -> File:
        """Return a new File with the given suffix."""
        new_name = self.name.with_suffix(suffix)
        return File(self.parent_process, new_name)


class ParentProcessPlaceholder(HasDirectory):
    """Placeholder parent_process for Files that don't have a Workflow/Process/Calculator as a parent_process.

    Also used when we don't want to store the parent in the database.
    """

    def __init__(self, parent_process, directory, engine, **kwargs):
        super().__init__(parent_process, directory, engine=engine, _directory_must_be_relative=False, **kwargs)

    def __repr__(self):
        return f'ParentProcessPlaceholder(directory={self.absolute_directory})'

    @classmethod
    def fromobj(cls, obj, replace_parents_with_placeholders=True):
        """Create a ParentProcessPlaceholder from a Process/Calculator object."""
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
        """Create a ParentProcessPlaceholder from a path."""
        return cls(None, path, engine)


def LocalFile(path: Path | str) -> File:
    """Return a file object that does not have a parent_process, and is assumed to be on the local filesystem."""
    from koopmans.engines.localhost import LocalhostEngine
    path = path if isinstance(path, Path) else Path(path)
    engine = LocalhostEngine()
    if path.is_dir():
        parent_process = ParentProcessPlaceholder.frompath(path.resolve(), engine=engine)
        return File(parent_process=parent_process, name='')
    else:
        parent_process = ParentProcessPlaceholder.frompath(path.parent.resolve(), engine=engine)
        return File(parent_process=parent_process, name=Path(path.name))
