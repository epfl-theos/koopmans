'''

Module koopmans for dealing with settings

Written by Edward Linscott May 2020

'''

import os
from collections import UserDict
from typing import Union, Type, Tuple, NamedTuple, Dict, Any, Optional, List
from koopmans.utils import units


class SettingsDict(UserDict):
    '''
    A dictionary-like class that has a few extra checks that are performed when setting values (e.g. when setting
    variables identified as paths it will convert them to absolute paths) as well as a few extra useful attributes

    Modelled off ase.calculators.Parameters which allows us to refer to "self.key", which returns "self['key']"
    '''

    def __init__(self, valid: List[str], defaults: Dict[str, Union[int, str, float, bool]] = {}, are_paths: List[str] = [], to_not_parse: List[str] = [], directory='', **kwargs):
        super().__init__(**kwargs)
        self.valid = valid
        self.defaults = defaults
        self.update(**defaults)
        self.are_paths = are_paths
        self.to_not_parse = to_not_parse
        self.directory = directory
        self.update(**kwargs)

    def __getattr__(self, key):
        if key in ['data', 'valid', 'defaults', 'update', 'are_paths', 'to_not_parse', '_to_not_parse', 'directory']:
            return self.__dict__[key]
        elif key in self.valid:
            return self.data.get(key, None)
        else:
            if key not in self.data:
                return dict.__getattribute__(self.data, key)
            return self.data[key]

    def __setattr__(self, key, value):
        if key in ['data', 'valid', 'defaults', 'update', 'are_paths', 'to_not_parse', '_to_not_parse', 'directory']:
            self.__dict__[key] = value
        else:
            self.data[key] = value

    def __getitem__(self, key: str):
        if key not in self.data:
            if key in self.defaults:
                self.data[key] = self.defaults[key]
            else:
                raise KeyError(key)
        return self.data[key]

    def __setitem__(self, key: str, value: Union[int, str, float, bool]):
        # Insisting that all values corresponding to paths are absolute
        if key in self.are_paths and value.startswith('/'):
            value = os.path.abspath(self.directory + '/' + self.value)

        super().__setitem__(key, value)

    def update(self, *args, **kwargs):
        if args:
            if len(args) > 1:
                raise TypeError(f"update expected at most 1 arguments, got {len(args)}")
            other = dict(args[0])
            for key in other:
                self.data[key] = other[key]
        for key in kwargs:
            self.data[key] = kwargs[key]

    def setdefault(self, key: str, value: Optional[Any] = None):
        if key not in self:
            self.data[key] = value
        return self.data[key]

    @property
    def to_not_parse(self):
        if not '_to_not_parse' in self.__dict__:
            self._to_not_parse = set(self.are_paths)
        return self._to_not_parse

    @to_not_parse.setter
    def to_not_parse(self, value: Union[list, set]):
        self._to_not_parse = self.to_not_parse.union(value)


class Setting(NamedTuple):
    name: str
    description: str
    kind: Union[Type, Tuple[Type, ...]]
    default: Union[str, bool, float, list, None]
    options: Union[tuple, None]


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
            if not isinstance(value, valid_setting.kind) and value is not None:
                if isinstance(valid_setting.kind, tuple):
                    raise ValueError(
                        f'{type(value).__name__} is an invalid type for "{key}" (must be '
                        'one of ' + '/'.join([t.__name__ for t in valid_setting.kind]) + ')')
                else:
                    raise ValueError(
                        f'{type(value).__name__} is an invalid type for "{key}" (must be '
                        f'{valid_setting.kind.__name__})')

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
