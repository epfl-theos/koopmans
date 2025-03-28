"""System I/O functions for `koopmans.utils`."""

from __future__ import annotations

import contextlib
import os
import shutil
from glob import glob
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

if TYPE_CHECKING:
    from koopmans.engines import Engine


def system_call(command: str, check_ierr: bool = True):
    """Make a system call and check the exit code."""
    raise ValueError('System calls are no longer allowed')


def symlink(src: Union[str, Path], dest: Union[str, Path], relative: bool = True, exist_ok: bool = False,
            force: bool = False):
    """Create a symlink from "src" to "dest"."""
    if isinstance(src, str) and '*' in src:
        # Follow the syntax of ln, whereby ln -s src/* dest/ will create multiple links
        for src_file in glob(src):
            symlink(src_file, dest, relative, exist_ok, force)
    else:
        # Sanitize input
        if isinstance(src, str):
            src = Path(src)
        if isinstance(dest, str):
            dest = Path(dest)

        if dest.is_dir():
            dest /= src.name
        dest = dest.absolute()
        src = src.absolute()

        if dest == src:
            raise OSError('Cannot symlink a file to itself')

        # Check if the src exists
        if not src.exists():
            raise FileNotFoundError(src)

        if relative:
            # The equivalent of ln -sr
            src = Path(os.path.relpath(src, dest.parent))
        else:
            # The equivalent of ln -s
            pass

        if force and dest.exists():
            # The equivalent of ln -sf
            dest.unlink()

        if exist_ok:
            try:
                dest.symlink_to(src)
            except FileExistsError:
                pass
        else:
            if relative:
                dest.symlink_to(src)


def symlink_tree(src: Union[str, Path], dest: Union[str, Path], exist_ok: bool = False, force: bool = False):
    """Create a symlink tree from "src" to "dest"."""
    if isinstance(src, str):
        src = Path(src)
    if isinstance(dest, str):
        dest = Path(dest)

    if not src.exists():
        raise FileNotFoundError(src)
    if not src.is_dir():
        raise NotADirectoryError(src)

    if not dest.exists():
        dest.mkdir(parents=True, exist_ok=exist_ok)

    for f in src.rglob('*'):
        if f.is_dir():
            (dest / f.relative_to(src)).mkdir(exist_ok=exist_ok)
        else:
            symlink(f, dest / f.relative_to(src), exist_ok=exist_ok, force=force)


def copy_file(src: Union[str, Path], dest: Union[str, Path], exist_ok: bool = False):
    """Copy a file from "src" to "dest"."""
    if '*' in str(src) or '?' in str(src):
        raise ValueError('Do not use wildcards in `utils.copy()`')

    # Sanitize input
    if isinstance(src, str):
        src = Path(src)
    if isinstance(dest, str):
        dest = Path(dest)

    if dest.is_dir() and not src.is_dir():
        dest /= src.name
    dest = dest.absolute()
    src = src.absolute()

    # Check if the src exists
    if not src.exists():
        raise FileNotFoundError(src)

    if dest.exists() and not exist_ok:
        raise FileExistsError(dest)
    else:
        if src.is_file():
            shutil.copy(src, dest)
        else:
            shutil.copytree(src, dest)


def chdir_logic(path: Union[Path, str]):
    """Change the working directory.

    Allows for the context "with chdir(path)". All code within this
    context will be executed in the directory "path"
    """
    # Ensure path is a Path object
    if not isinstance(path, Path):
        path = Path(path)

    this_dir = Path.cwd()

    # Create path if it does not exist
    path.mkdir(parents=True, exist_ok=True)

    # Move to the directory
    os.chdir(path)
    try:
        yield
    finally:
        # Return to the original directory
        os.chdir(this_dir)


@contextlib.contextmanager
def chdir(path: Union[Path, str]):
    """Return a context that changes the working directory (returns to the original directory when done)."""
    return chdir_logic(path)


@contextlib.contextmanager
def set_env(**environ):
    """Temporarily set the process environment variables.

    :type environ: dict[str, unicode]
    :param environ: Environment variables to set
    """
    old_environ = dict(os.environ)
    os.environ.update(environ)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(old_environ)


def find_executable(program: Union[Path, str]) -> Optional[Path]:
    """Find where the executable 'program' is located."""
    if isinstance(program, str):
        program = Path(program)

    # Equivalent to the unix command "which"
    def is_exe(fpath: Path):
        return fpath.is_file() and os.access(fpath, os.X_OK)

    fpath = program.parent

    if fpath.samefile('.'):
        if is_exe(Path(fpath) / program):
            return Path(fpath) / program

    for path in os.environ["PATH"].split(os.pathsep):
        exe_file = Path(path) / program
        if is_exe(exe_file):
            return exe_file

    return None


class HasDirectory:
    """A class that has a directory and a base directory.

    This class will eventually be merged with the Process class
    For the moment it only contains information related to parent processes and directories. Once calculators and
    workflows have been transformed to have pydantic inputs and outputs then those classes will be able to inherit
    directly from Process.
    """

    __slots__ = ['parent_process', '_directory', '_base_directory', 'engine', '_directory_must_be_relative']

    def __init__(self, parent_process: Optional[HasDirectory] = None, directory=None, base_directory=Path(),
                 engine: Optional[Engine] = None,
                 _directory_must_be_relative=False):
        self._base_directory: Optional[Path] = None
        self._directory: Optional[Path] = None
        self.engine: Optional[Engine] = engine
        self.parent_process = parent_process
        self._directory_must_be_relative = _directory_must_be_relative

        if not self.parent_process:
            self.base_directory = base_directory
        self.directory = directory

    @property
    def directory(self) -> Path | None:
        """Return the directory of the process."""
        return self._directory

    @directory.setter
    def directory(self, value: Path | str | None):
        if value is None:
            return

        # Sanitize input
        if isinstance(value, str):
            value = Path(value)

        # Sanity checks
        if value.is_absolute() and self._directory_must_be_relative:
            raise ValueError(
                f'{self.__class__.__name__} directory must be a relative path '
                f'(relative to {self.__class__.__name__},base directory)'
            )

        self._directory = value

    @property
    def uid(self) -> str:
        """Return the UID of the process."""
        return str(self.directory)

    @property
    def base_directory(self) -> Path | None:
        """Return the base (i.e. root) directory of the entire workflow."""
        if self.parent_process:
            return self.parent_process.base_directory
        else:
            return self._base_directory

    @base_directory.setter
    def base_directory(self, value: Path | str | None):
        if self.parent_process is not None:
            raise ValueError('Do not directly set `base_directory` for objects with a parent_process')
        if isinstance(value, str):
            value = Path(value)
        self._base_directory = None if value is None else value.resolve()

    @property
    def absolute_directory(self) -> Path | None:
        """Return the absolute path of the directory."""
        if self.directory is None:
            return None
        if self.base_directory is None:
            return None
        return self.base_directory / self.directory

    def directory_has_been_set(self) -> bool:
        """Check if the directory has been set."""
        return self._directory is not None
