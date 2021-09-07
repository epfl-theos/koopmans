'''

System I/O functions for koopmans.utils

Written by Edward Linscott May 2020

'''

import os
import subprocess
import contextlib


def system_call(command, check_ierr=True):
    '''
    Make a system call and check the exit code
    '''
    ierr = subprocess.call(command, shell=True)
    if ierr > 0 and check_ierr:
        raise OSError(f'{command} exited with exit code {ierr}')


def mkdir(path):
    # Creates a (possibly nested) directory
    relpath = os.path.relpath(path, os.getcwd())
    split_relpath = relpath.split('/')
    for i in range(len(split_relpath)):
        subdir = '/'.join(split_relpath[:i + 1])
        if not os.path.isdir(subdir):
            system_call(f'mkdir {subdir}')


@contextlib.contextmanager
def chdir(path):
    # Allows for the context "with chdir(path)". All code within this
    # context will be executed in the directory "path"
    this_dir = os.getcwd()

    # Create path if it does not exist
    mkdir(path)

    # Move to the directory
    os.chdir(path)
    try:
        yield
    finally:
        # Return to the original directory
        os.chdir(this_dir)


def find_executable(program):
    # Equivalent to the unix command "which"
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    if program[0] == '~':
        program = program.replace('~', os.environ["HOME"], 1)

    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
