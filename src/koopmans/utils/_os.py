'''

System I/O functions for koopmans.utils

Written by Edward Linscott May 2020

'''

import contextlib
import os
import subprocess
from glob import glob
from pathlib import Path
from typing import Optional, Union


def system_call(command: str, check_ierr: bool = True):
    '''
    Make a system call and check the exit code
    '''
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
            dest.symlink_to(src)


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
