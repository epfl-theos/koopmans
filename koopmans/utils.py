'''

utils module for koopmans

Written by Edward Linscott May 2020

'''

import os
import warnings
from collections import namedtuple
import subprocess
import contextlib
from ase.units import create_units


# Quantum Espresso -- and koopmans -- uses CODATA 2006 internally
units = create_units('2006')


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


def calc_diff(calcs, silent=False):
    # Returns the differences in the settings of a list of calculators

    # If calcs is a dict, convert it to a list (we only need the values)
    if isinstance(calcs, dict):
        calcs = calcs.values()

    diffs = []

    settings = [c._settings for c in calcs]

    keys = set([k for s in settings for k in s.keys()])
    for key in sorted(keys):
        vals = [s.get(key, None) for s in settings]
        if len(set(vals)) > 1:
            if not silent:
                print(f'{key}: ' + ', '.join(map(str, vals)))
            diffs.append(key)

    return diffs


Setting = namedtuple('Setting', ['name', 'description', 'type', 'default', 'options'])


def check_settings(settings, valid_settings, mandatory_settings=[], physicals=[], do_not_lower=[]):
    '''
    Checks settings (a dict) against a list of Settings (namedtuples), and adding defaults where they are not present
    Will additionally check that all mandatory settings are present, and arse any units for those settings named in
    the 'physicals' list
    '''

    valid_settings_dict = {s.name: s for s in valid_settings}

    # Check all mandatory keys are provided:
    for mandatory_setting in mandatory_settings:
        if mandatory_setting not in settings:
            raise ValueError(f'The mandatory key {mandatory_setting} has not been provided')

    # Populate checked_settings with the default values
    checked_settings = {s.name: s.default for s in valid_settings}

    for key, value in settings.items():

        # Check key is a valid keyword
        if key in valid_settings_dict:
            valid_setting = valid_settings_dict[key]

            # Lowers any uppercase strings
            if isinstance(value, str) and key not in do_not_lower:
                value = value.lower()

            # Check value is the correct type
            if not isinstance(value, valid_setting.type) and value is not None:
                if isinstance(valid_setting.type, tuple):
                    raise ValueError(
                        f'{type(value).__name__} is an invalid type for "{key}" (must be '
                        'one of ' + '/'.join([t.__name__ for t in valid_setting.type]) + ')')
                else:
                    raise ValueError(
                        f'{type(value).__name__} is an invalid type for "{key}" (must be '
                        f'{valid_setting.type.__name__})')

            # Check value is among the valid options
            if valid_setting.options is not None and valid_setting.default is not None and value not in \
                    valid_setting.options:
                raise ValueError(
                    f'"{value}" is an invalid value for "{key}" (options are {"/".join(valid_setting.options)})')

            checked_settings[key] = value
        else:
            raise ValueError(f'"{key}" is not a recognised setting')

    # Parse physicals
    for physical in physicals:
        checked_settings[physical] = parse_physical(checked_settings[physical])

    return checked_settings


def parse_physical(value):
    '''
    Takes in a value that potentially has a unit following a float,
    converts the value to ASE's default units (Ang, eV), and returns
    that value
    '''

    if isinstance(value, float):
        return value
    elif isinstance(value, int):
        return float(value)
    elif isinstance(value, str):
        splitline = value.strip().split()
        if len(splitline) == 1:
            return float(splitline[0])
        elif len(splitline) == 2:
            [value, val_units] = splitline
            value = float(value)
            matching_units = [
                u for u in units if u.lower() == val_units.lower()]
            if len(matching_units) == 1:
                return value * units[matching_units[0]]
            elif len(matching_units) > 1:
                raise ValueError(
                    f'Multiple matches for {val_units} found; this should not happen')
            else:
                raise NotImplementedError(
                    f'{val_units} not implemented in koopmans.utils.parse_physical')
