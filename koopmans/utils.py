'''

utils module for python_KI

Written by Edward Linscott May 2020

'''

import os
import warnings
import subprocess


def _warning(message, category=UserWarning, filename='', lineno=-1, file=None, line=None):
    '''
    Monkey-patching warnings.warn
    '''
    print(f'{category.__name__}: {message}')


warnings.showwarning = _warning


def warn(message):
    '''
    Allowing the monkey-patched warnings.warn to be imported as utils.warn
    '''
    warnings.warn(message)


def system_call(command, check_ierr=True):
    '''
    Make a system call and check the exit code
    '''
    ierr = subprocess.call(command, shell=True)
    if ierr > 0 and check_ierr:
        raise OSError(f'{command} exited with exit code {ierr}')


def find_executable(program):
    # Equivalent to the unix command "which"
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
