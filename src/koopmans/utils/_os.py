'''

System I/O functions for koopmans.utils

Written by Edward Linscott May 2020

'''

import contextlib
import os
import shutil
import subprocess
from glob import glob
from pathlib import Path
from typing import List, Optional, Protocol, Union, runtime_checkable


def system_call(command: str, check_ierr: bool = True):
    '''
    Make a system call and check the exit code
    '''

    raise ValueError('System calls are no longer allowed')
    ierr = subprocess.call(command, shell=True)
    if ierr > 0 and check_ierr:
        raise OSError(f'{command} exited with exit code {ierr}')


def symlink(src: Union[str, Path], dest: Union[str, Path], relative: bool = True, exist_ok: bool = False,
            force: bool = False):
    # Create a symlink of "src" at "dest"
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


def copy(src: Union[str, Path], dest: Union[str, Path], exist_ok: bool = False):
    # Copy a file from "src" to "dest"
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


@contextlib.contextmanager
def chdir(path: Union[Path, str]):
    # Allows for the context "with chdir(path)". All code within this
    # context will be executed in the directory "path"

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
def set_env(**environ):
    """
    Temporarily set the process environment variables.

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
    # This class will eventually be merged with the Process class
    # For the moment it only contains information related to parents and directories. Once calculators and workflows
    # have been transformed to have pydantic inputs and outputs then those classes will be able to inherit directly
    # from Process

    __slots__ = ['parent', '_directory', '_base_directory']

    def __init__(self, parent=None, directory=None, base_directory=None):
        self._base_directory: Optional[Path] = None
        self._directory: Optional[Path] = None
        self.parent: Optional[HasDirectory] = parent

        if not self.parent:
            self.base_directory = base_directory
        self.directory = directory

    @property
    def directory(self) -> Path | None:
        return self._directory

    @property
    def uid(self) -> str:
        return str(self.directory)

    @directory.setter
    def directory(self, value: Path | str | None):
        if value is None:
            return

        # Sanitize input
        if isinstance(value, str):
            value = Path(value)

        # Sanity checks
        if value.is_absolute():
            raise ValueError(
                f'{self.__class__.__name__} directory must be a relative path (relative to {self.__class__.__name__},base directory)')

        self._directory = value

    @property
    def base_directory(self) -> Path | None:
        if self.parent:
            return self.parent.base_directory
        else:
            return self._base_directory

    @base_directory.setter
    def base_directory(self, value: Path | str | None):
        if self.parent is not None:
            raise ValueError('Do not directly set `base_directory` for objects with a parent')
        if isinstance(value, str):
            value = Path(value)
        self._base_directory = None if value is None else value.resolve()

    @property
    def absolute_directory(self) -> Path | None:
        if self.directory is None:
            return None
        if self.base_directory is None:
            return None
        return self.base_directory / self.directory

    def directory_has_been_set(self) -> bool:
        return self._directory is not None


def get_binary_content(source: HasDirectory, relpath: Path | str) -> bytes:
    if isinstance(relpath, str):
        relpath = Path(relpath)
    assert source.absolute_directory is not None
    with open(source.absolute_directory / relpath, "rb") as f:
        flines = f.read()
    return flines


def write_binary_content(dst_file: Path | str, merged_filecontents: bytes):
    if isinstance(dst_file, str):
        dst_file = Path(dst_file)
    dst_file.parent.mkdir(parents=True, exist_ok=True)
    with open(dst_file, "wb") as f:
        f.write(merged_filecontents)


def get_content(source: HasDirectory | None, relpath: Path | str) -> List[str]:
    if isinstance(relpath, str):
        relpath = Path(relpath)

    if source is None or source.absolute_directory is None:
        full_path = relpath
    else:
        full_path = source.absolute_directory / relpath

    with open(full_path, "r") as f:
        flines = [l.strip('\n') for l in f.readlines()]
    return flines


def write_content(dst_file: Path | str, merged_filecontents: List[str]):
    if isinstance(dst_file, str):
        dst_file = Path(dst_file)
    dst_file.parent.mkdir(parents=True, exist_ok=True)
    with open(dst_file, "w") as f:
        f.write('\n'.join(merged_filecontents))


def remove(path: Union[Path, str]):
    if isinstance(path, str):
        path = Path(path)
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()
