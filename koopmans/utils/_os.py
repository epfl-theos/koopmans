'''

System I/O functions for koopmans.utils

Written by Edward Linscott May 2020

'''

import os
from pathlib import Path
from typing import Union
import subprocess
import contextlib


def system_call(command, check_ierr=True):
    '''
    Make a system call and check the exit code
    '''
    ierr = subprocess.call(command, shell=True)
    if ierr > 0 and check_ierr:
        raise OSError(f'{command} exited with exit code {ierr}')


def mkdir(path: Path):
    # Creates a (possibly nested) directory
    relpath = Path(os.path.relpath(path, Path.cwd()))
    for p in reversed(relpath.parents):
        if not p.is_dir():
            p.mkdir()


@contextlib.contextmanager
def chdir(path: Union[Path, str]):
    # Allows for the context "with chdir(path)". All code within this
    # context will be executed in the directory "path"

    # Ensure path is a Path object
    if isinstance(path, str):
        path = Path(path)

    this_dir = Path.cwd()

    # Create path if it does not exist
    mkdir(path)

    # Move to the directory
    os.chdir(path)
    try:
        yield
    finally:
        # Return to the original directory
        os.chdir(this_dir)


def find_executable(program: Path):
    # Equivalent to the unix command "which"
    def is_exe(fpath: Path):
        return fpath.is_file() and os.access(fpath, os.X_OK)

    fpath = program.parent
    if fpath.samefile('.'):
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = Path(path) / program
            if is_exe(exe_file):
                return exe_file

    return None
