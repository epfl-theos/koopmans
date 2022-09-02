'''

Module koopmans for dealing with settings

Written by Edward Linscott May 2020

'''

import os
from collections import UserDict
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Optional, Tuple, Type, Union

import numpy as np

from koopmans.utils import units


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


class Setting(NamedTuple):
    name: str
    description: str
    kind: Union[Type, Tuple[Type, ...]]
    default: Union[str, bool, float, list, None]
    options: Union[tuple, None]


class SettingsDict(UserDict):
    '''
    A dictionary-like class that has a few extra checks that are performed when setting values (e.g. when setting
    variables identified as paths it will convert them to absolute paths) as well as a few extra useful attributes

    Modelled off ase.calculators.Parameters which allows us to refer to "self.key", which returns "self['key']"

    Arguments:
    valid -- list of valid settings
    defaults -- dict of defaults for each setting
    are_paths -- list of settings that correspond to paths
    to_not_parse -- list of settings that should not be parsed algebraically
    directory -- the directory in which the calculation is being run (used to enforce all path settings to be absolute
                 paths)
    physicals -- list of keywords that have accompanying units from input

    By altering SettingsDict.use_relative_paths = True/False you can change if settings that are paths are returned as
    absolute or relative paths (absolute paths are returned by default)
    '''

    # Need to provide these here to allow copy.deepcopy to perform the checks in __getattr__
    valid: List[str] = []
    data: Dict[str, Any] = {}
    are_paths: List[str] = []
    to_not_parse: List[str] = []
    physicals: List[str] = []

    def __init__(self, valid: List[str], defaults: Dict[str, Any] = {}, are_paths: List[str] = [],
                 to_not_parse: List[str] = [], directory='', physicals: List[str] = [], **kwargs):
        super().__init__()
        self.valid = valid + self._other_valid_keywords
        self.are_paths = are_paths
        self.defaults = {k: v for k, v in defaults.items() if k not in self.are_paths}
        self.defaults.update(**{k: Path(v) for k, v in defaults.items() if k in self.are_paths})
        self.to_not_parse = to_not_parse
        self.directory = directory if directory is not None else Path.cwd()
        self.physicals = physicals
        self.update(**defaults)
        self.update(**kwargs)
        self.use_relative_paths = False

    def __getattr__(self, name):
        if name != 'valid' and self.is_valid(name):
            if name in self.data:
                return self.__getitem__(name)
            else:
                return None
        else:
            try:
                super().__getattr__(name)
            except AttributeError:
                raise AttributeError(name)

    def __setattr__(self, name, value):
        if self.is_valid(name):
            self.__setitem__(name, value)
        else:
            super().__setattr__(name, value)

    def __getitem__(self, key: str):
        if key not in self.data:
            if key in self.defaults:
                self.data[key] = self.defaults[key]
            else:
                raise KeyError(key)
        if key in self.are_paths and self.use_relative_paths:
            return Path(os.path.relpath(self.data[key], self.directory))
        else:
            return self.data[key]

    def __setitem__(self, key: str, value: Any):
        # If we set something to "None", simply remove it from the dictionary
        if value is None:
            self.pop(key, None)
            return

        # Insisting that all values corresponding to paths are absolute and are Path objects
        if key in self.are_paths:
            if isinstance(value, str):
                value = Path(value)
            elif not isinstance(value, Path):
                raise ValueError(f'{key} must be either a string or a Path')
            if not value.is_absolute():
                value = (self.directory / value).resolve()

        # Parse any units provided
        if key in self.physicals:
            value = parse_physical(value)

        # Perform additional checks that the key and corresponding value is valid
        self._check_before_setitem(key, value)

        # Set the item
        super().__setitem__(key, value)

    def is_valid(self, name: str) -> bool:
        # Check if a keyword is valid. This is a separate subroutine to allow child classes to overwrite it
        # e.g. QE calculators want to be able to set keywords such as Hubbard(i) where i is an arbitrary integer
        return name in self.valid

    def update(self, *args: Any, **kwargs: Any) -> None:
        if args:
            if len(args) > 1:
                raise TypeError(f"update expected at most 1 arguments, got {len(args)}")
            other = dict(args[0])
            for key in other:
                self.__setitem__(key, other[key])
        for key in kwargs:
            self.__setitem__(key, kwargs[key])

    def setdefault(self, key: str, value: Optional[Any] = None):
        if key not in self:
            self.data[key] = value
        return self.data[key]

    def _check_before_setitem(self, key, value):
        if not self.is_valid(key):
            raise KeyError(f'{key} is not a valid setting')
        return

    def parse_algebraic_setting(self, expr, **kwargs):
        # Checks if expr is defined algebraically, and evaluates them
        if not isinstance(expr, str):
            return expr
        if all([c.isalpha() or c in ['_', '"', "'"] for c in expr]):
            return expr.strip('"').strip("'")

        for replace in [('/', ' / '), ('*', ' * '), ('-', ' - '), ('+', ' + '), ('e - ', 'e-'), ('e + ', 'e+')]:
            expr = expr.replace(*replace)
        expr = expr.split()

        for i, term in enumerate(expr):
            if term in ['*', '/', '+', '-']:
                continue
            elif all([c.isalpha() for c in term]):
                if term in kwargs:
                    expr[i] = kwargs[term]
                elif getattr(self, term, None) is None:
                    raise ValueError('Failed to parse ' + ''.join(map(str, expr)))
                else:
                    expr[i] = getattr(self, term)
            elif len(expr) == 1 and any([c.isalpha() for c in term]):
                return term.strip('"').strip("'")
            else:
                expr[i] = float(term)

        value = float(expr[0])
        for op, term in zip(expr[1::2], expr[2::2]):
            if op == '*':
                value *= float(term)
            elif op == '/':
                value /= float(term)
            elif op == '+':
                value += float(term)
            elif op == '-':
                value -= float(term)
            else:
                raise ValueError('Failed to parse ' + ''.join([str(e) for e in expr]))

        return value

    def parse_algebraic_settings(self, **kwargs):
        # Checks self.parameters for keywords defined algebraically, and evaluates them.
        # kwargs can be used to provide values for keywords such as 'nelec', which the
        # SettingsDict may not have access to but is a useful keyword for scaling extensive
        # parameters
        for key in self.data.keys():
            if key in self.to_not_parse:
                continue
            self.data[key] = self.parse_algebraic_setting(self.data[key], **kwargs)

    @property
    def _other_valid_keywords(self):
        return ['pseudopotentials', 'gamma_only', 'kpts', 'koffset']

    def todict(self):
        # Construct a minimal representation of this dictionary. Most of the requisite information
        # (defaults, valid, are_paths, etc) is contained in the class itself so we needn't store this
        dct = {}
        for k in self.data:
            v = self[k]
            if k in self.defaults:
                if isinstance(v, np.ndarray) and np.all(v == self.defaults[k]):
                    continue
                elif v == self.defaults[k]:
                    continue
            dct[k] = v

        # Make sure if a default is missing entirely (which means it must have been manually wiped) we store this as
        # key: None
        for k in self.defaults:
            if k not in self.data:
                dct[k] = None

        # Adding information required by the json decoder
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__

        return dct

    @classmethod
    def fromdict(cls, dct):
        return cls(**dct)

    def briefrepr(self) -> str:
        entries = str(self)[1:-1]
        try:
            comma_index = entries[:80].rindex(',')
            body = entries[:comma_index + 2]
        except ValueError:
            body = ''
        return '{' + body + '...}'


class SettingsDictWithChecks(SettingsDict):
    def __init__(self, settings: List[Setting], **kwargs):
        self.settings = settings
        super().__init__(valid=[s.name for s in settings],
                         defaults={s.name: s.default for s in settings if s.default is not None},
                         are_paths=[s.name for s in settings if s.kind == Path], **kwargs)

    def _check_before_setitem(self, key, value):
        super()._check_before_setitem(key, value)

        # Always accept pseudopotentials and k-point data
        if key in self._other_valid_keywords:
            return

        # Fetch the record of the setting in question
        [setting] = [s for s in self.settings if s.name == key]

        # Check the value is the valid type
        if isinstance(setting.kind, tuple):
            if not any([isinstance(value, k) for k in setting.kind]):
                raise ValueError(f'{setting.name} must be a ' + '/'.join([str(k) for k in setting.kind]))
        else:
            if not isinstance(value, setting.kind):
                raise ValueError(f'{setting.name} must be a {setting.kind}')

        # Check the value is among the valid options
        if setting.options is not None and value not in setting.options:
            raise ValueError(f'{setting.name} may only be set to ' + '/'.join([str(o) for o in setting.options]))


class IbravDict():
    def __setitem__(self, key: str, value: Any) -> None:
        if key == 'celldms':
            if not isinstance(value, dict):
                raise ValueError('celldms should be a dictionary')
            for k, v in value.items():
                self[f'celldm({k})'] = v
            return
        else:
            return super().__setitem__(key, value)  # type: ignore


kc_wann_defaults = {'outdir': 'TMP',
                    'kcw_iverbosity': 1,
                    'kcw_at_ks': False,
                    'homo_only': False,
                    'read_unitary_matrix': True,
                    'lrpa': False,
                    'check_ks': True,
                    'have_empty': True,
                    'has_disentangle': True}
