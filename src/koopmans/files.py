"""A class for representing files in a general way by tethering them to a parent process."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Annotated, Any, Generator

from pydantic import BeforeValidator
from typing_extensions import Self

from koopmans.pydantic import BaseModel


class FileABC(BaseModel, ABC):
    """An abstract base class for representing a file."""

    path: Path

    # @abstractmethod
    # def aspath(self) -> Path:
    #     """Return a `Path` object represetning this object."""

    # def __repr__(self):
    #     return f'File({self.aspath()})'

    # def __str__(self):
    #     return str(self.aspath())

    @property
    def stem(self) -> str:
        """Return the stem of the file, replicating the behavior of Path.stem."""
        return self.path.stem

    @property
    def suffix(self) -> str:
        """Return the file suffix, replicating the behavior of Path.suffix."""
        return self.path.suffix

    # @property
    # @abstractmethod
    # def parent(self) -> Self:
    #     """Return the parent directory of this file."""

    # @property
    # def parents(self) -> Generator[Self, None, None]:
    #     """Return a generator of all parent directories."""
    #     parent = self.parent
    #     while parent != self:
    #         yield parent
    #         parent = parent.parent

    # @abstractmethod
    # def copy_to(self, dst: FileABC, exist_ok: bool = False):
    #     """Copy this file to dst."""

    # @abstractmethod
    # def exists(self) -> bool:
    #     """Return true if this file exists."""

    # @abstractmethod
    # def read_text(self) -> str:
    #     """Read text from this file."""

    # @abstractmethod
    # def read_bytes(self) -> bytes:
    #     """Read bytes from this file."""

    # @abstractmethod
    # def write_text(self, content: str):
    #     """Write text to this file."""

    # @abstractmethod
    # def write_bytes(self, content: bytes):
    #     """Write bytes to this file."""

    # @abstractmethod
    # def rglob(self, pattern: str) -> Generator[Self, None, None]:
    #     """Iterate over all files within this directory and its subdirectories that match the pattern."""

    # @abstractmethod
    # def glob(self, pattern: str) -> Generator[Self, None, None]:
    #     """Iterate over all files within this directory that match the pattern."""

    # @abstractmethod
    # def symlink_to(self, target: FileABC, overwrite: bool = False, recursive: bool = False):
    #     """Symbolically link this file to another file."""

    # @abstractmethod
    # def unlink(self):
    #     """Remove this file/directory."""

    # @abstractmethod
    # def mkdir(self, *args, **kwargs):
    #     """Create this directory."""

    # @abstractmethod
    # def is_dir(self) -> bool:
    #     """Return true if the file is a directory."""

    # @abstractmethod

    def rename(self, dst: Self) -> None:
        """Rename this file, replicating the behavior of Path.rename."""
        self.path.rename(dst.path)

    @property
    def name(self) -> str:
        """Return the name of the file, replicating the behavior of Path.name."""
        return self.path.name

    @abstractmethod
    def with_stem(self, stem: str) -> Self:
        """Return a new object with the given stem."""

    @abstractmethod
    def with_suffix(self, suffix: str) -> Self:
        """Return a new File with the given suffix."""


# def accept_subclasses(v: Any) -> Process[Any, Any]:
#     if isinstance(v, Process):
#         return v
#     return Process.model_validate(v)


class File(FileABC):
    """An abstract way of representing a file.

    Because a file may not exist locally (specifically, when koopmans is run with AiiDA), we need a way of
    referring to a file that is more general than an absolute path. This class achieves this by storing a
    file as a parent_process (which is a Process, Calculator, or some other object that exists in a directory known
    to koopmans/AiiDA) and a name (which is the path of the file relative to the parent_process's directory).

    We also need to delegate file creation/modification/deletion to the engine.
    """

    parent_process_uid: str

    # def aspath(self) -> Path:
    #     """Return the File as a Path object."""
    #     if self.parent_process.directory is None:
    #         return self.relative_path
    #     else:
    #         return self.parent_process.directory / self.relative_path

    # @property
    # def parent(self) -> File:
    #     """Return the parent directory."""
    #     return File(parent_process=self.parent_process, relative_path=self.relative_path.parent)

    # def copy_to(self, dst: FileABC, exist_ok: bool = False):
    #     """Copy this file to another."""
    #     if not isinstance(dst, File):
    #         raise NotImplementedError('Can only copy File type to File type at this time')
    #     logger = logging.getLogger(__name__)
    #     logger.info(f'Copying {self.aspath()} to {dst.aspath()}')
    #     self._engine.copy_file(self, dst, exist_ok=exist_ok)

    # def exists(self) -> bool:
    #     """Return true if the file exists."""
    #     return self._engine.file_exists(self)

    # def read_text(self) -> str:
    #     """Read text from this file."""
    #     logger = logging.getLogger(__name__)
    #     logger.info(f'Reading text from {self.aspath()}', stacklevel=2)
    #     return self._engine.read_file(self, binary=False)

    # def read_bytes(self) -> bytes:
    #     """Read bytes from this file."""
    #     logger = logging.getLogger(__name__)
    #     logger.info(f'Reading bytes from {self.aspath()}', stacklevel=2)
    #     return self._engine.read_file(self, binary=True)

    # def write_text(self, content: str):
    #     """Write text to this file."""
    #     logger = logging.getLogger(__name__)
    #     logger.info(f'Writing text to {self.aspath()}', stacklevel=2)
    #     self._engine.write_file(content, self)

    # def write_bytes(self, content: bytes):
    #     """Write bytes to this file."""
    #     logger = logging.getLogger(__name__)
    #     logger.info(f'Writing bytes to {self.aspath()}', stacklevel=2)
    #     self._engine.write_file(content, self)

    # def rglob(self, pattern: str) -> Generator[File, None, None]:
    #     """Iterate over all files within this directory and its subdirectories that match the pattern."""
    #     yield from self._engine.glob(self, pattern, recursive=True)

    # def glob(self, pattern: str) -> Generator[File, None, None]:
    #     """Iterate over all files within this directory that match the pattern."""
    #     yield from self._engine.glob(self, pattern, recursive=False)

    # def symlink_to(self, target: FileABC, overwrite: bool = False, recursive: bool = False):
    #     """Symbolically link this file to another file."""
    #     # Create a symbolic link at self that points to target
    #     if not isinstance(target, File):
    #         raise NotImplementedError('Can only symlink File type to File type at this time')
    #     logger = logging.getLogger(__name__)
    #     logger.info(f'Creating symlink from {self.aspath()} to {target.aspath()}', stacklevel=2)
    #     self._engine.link_file(target, self, overwrite=overwrite, recursive=recursive)

    # def unlink(self):
    #     """Remove this file/directory."""
    #     if self.is_dir():
    #         self._engine.rmdir(self)
    #     else:
    #         self._engine.unlink_file(self)
    #     logger = logging.getLogger(__name__)
    #     logger.info(f'Removing {self.aspath()}')

    # def mkdir(self, *args, **kwargs):
    #     """Create this directory."""
    #     logger = logging.getLogger(__name__)
    #     logger.info(f'Creating directory {self.aspath()}')
    #     self._engine.mkdir(self, *args, **kwargs)

    # def is_dir(self) -> bool:
    #     """Return true if the file is a directory."""
    #     return self._engine.file_is_dir(self)

    # def rename(self, dst: File) -> None:
    #     """Rename this file, replicating the behavior of Path.rename."""
    #     return self._engine.rename_file(self, dst)

    def with_stem(self, stem: str) -> File:
        """Return a new File with the given stem."""
        new = self.path.with_stem(stem)
        return File(parent_process_uid=self.parent_process_uid, path=new)

    def with_suffix(self, suffix: str) -> File:
        """Return a new File with the given suffix."""
        new = self.path.with_suffix(suffix)
        return File(parent_process_uid=self.parent_process_uid, path=new)

    # def __gt__(self, other: Any) -> bool:
    #     if not isinstance(other, File):
    #         raise TypeError(f'Cannot compare File with {type(other)}')
    #     return self.path > other.path

    # def __lt__(self, other: Any) -> bool:
    #     return not self > other

    def __truediv__(self, other: Any) -> File:
        assert isinstance(other, Path) or isinstance(other, str)
        return File(parent_process_uid=self.parent_process_uid, path=self.path / other)


class LocalFile(FileABC):
    """A local file.

    This class simply wraps the Path class.
    """

    def with_stem(self, stem: str) -> LocalFile:
        """Return a new object with the given stem."""
        new = self.path.with_stem(stem)
        return LocalFile(path=new)

    def with_suffix(self, suffix: str) -> LocalFile:
        """Return a new File with the given suffix."""
        new = self.path.with_suffix(suffix)
        return LocalFile(path=new)

    def __truediv__(self, other: Any) -> LocalFile:
        assert isinstance(other, Path) or isinstance(other, str)
        return LocalFile(path=self.path / other)

    @property
    def parent(self) -> LocalFile:
        """Return the parent directory of this file."""
        return LocalFile(path=self.path.parent)

    @property
    def parents(self) -> Generator[LocalFile, None, None]:
        """Return all parent directories of this file."""
        for parent in self.path.parents:
            yield LocalFile(path=parent)

    # def aspath(self) -> Path:
    #     return self.relative_path

    # @property
    # def parent(self) -> LocalFile:
    #     return LocalFile(relative_path=self.relative_path.parent)

    # def copy_to(self, dst: FileABC, exist_ok: bool = False):
    #     """Copy this file to dst."""
    #     if isinstance(dst, LocalFile):
    #         raise NotImplementedError()
    #     elif isinstance(dst, File):
    #         raise NotImplementedError()

    # def exists(self) -> bool:
    #     return self.relative_path.exists()

    # def read_text(self) -> str:
    #     """Read text from this file."""
    #     with open(self.relative_path, 'r') as f:
    #         content = f.read()
    #     return content

    # def read_bytes(self) -> bytes:
    #     """Read bytes from this file."""
    #     with open(self.relative_path, 'rb') as f:
    #         content = f.read()
    #     return content

    # def write_text(self, content: str):
    #     """Write text to this file."""
    #     with open(self.relative_path, 'w') as f:
    #         f.write(content)

    # def write_bytes(self, content: bytes):
    #     """Write bytes to this file."""
    #     with open(self.relative_path, 'wb') as f:
    #         f.write(content)

    # def rglob(self, pattern: str) -> Generator[LocalFile, None, None]:
    #     yield from [LocalFile(relative_path=p) for p in self.relative_path.rglob(pattern)]

    # def glob(self, pattern: str) -> Generator[LocalFile, None, None]:
    #     yield from [LocalFile(relative_path=p) for p in self.relative_path.glob(pattern)]

    # def symlink_to(self, target: FileABC, overwrite: bool = False, recursive: bool = False):
    #     """Symbolically link this file to another file."""
    #     if isinstance(target, LocalFile):
    #         raise NotImplementedError()
    #     elif isinstance(target, File):
    #         raise NotImplementedError()
    #     else:
    #         raise ValueError('target must be a LocalFile or File')

    # def unlink(self):
    #     self.relative_path.unlink()

    # def mkdir(self, *args, **kwargs):
    #     self.relative_path.mkdir(*args, **kwargs)

    # def is_dir(self) -> bool:
    #     return self.relative_path.is_dir()

    # def rename(self, dst: LocalFile) -> None:
    #     """Rename this file, replicating the behavior of Path.rename."""
    #     self.relative_path.rename(dst.relative_path)

    # def with_stem(self, stem: str) -> LocalFile:
    #     """Return a new object with the given stem."""
    #     new = self.relative_path.with_stem(stem)
    #     return LocalFile(relative_path=new)

    # def with_suffix(self, suffix: str) -> LocalFile:
    #     new = self.relative_path.with_suffix(suffix)
    #     return LocalFile(relative_path=new)

    # def __truediv__(self, other: Any) -> LocalFile:
    #     assert isinstance(other, Path) or isinstance(other, str)
    #     return LocalFile(relative_path=self.relative_path / other)


class LinkedFile(BaseModel):
    """A link from a source to a destination file.

    This can be either a symlink or a copy (see `FileViaSymlink`, `FileViaRecursiveSymlink`, and `FileViaCopy` below).
    """

    dest: Path
    src: File
    symlink: bool = False
    recursive_symlink: bool = False
    overwrite: bool = False


def as_linkedfile(*, symlink: bool, recursive_symlink: bool, overwrite: bool):
    """Return a validator that allows conversion of File to LinkedFile.

    Allows specification of symlink, recursive_symlink, and overwrite parameters.
    """

    def convert(v: Any):
        if v is None or isinstance(v, LinkedFile):
            return v
        if isinstance(v, File):
            return LinkedFile(dest=v.path, src=v, symlink=symlink, recursive_symlink=recursive_symlink,
                              overwrite=overwrite)
        raise TypeError("must be File or LinkedFile")
    return BeforeValidator(convert)


# Shortcuts for common LinkedFile types
FileViaSymlink = Annotated[LinkedFile, as_linkedfile(symlink=True, recursive_symlink=False, overwrite=False)]
FileViaRecursiveSymlink = Annotated[LinkedFile, as_linkedfile(symlink=True, recursive_symlink=True, overwrite=False)]
FileViaCopy = Annotated[LinkedFile, as_linkedfile(symlink=False, recursive_symlink=False, overwrite=False)]
