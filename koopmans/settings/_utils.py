'''

Module koopmans for dealing with settings

Written by Edward Linscott May 2020

'''

from koopmans import pseudopotentials
from pathlib import Path
from collections import UserDict
from typing import Union, Type, Tuple, NamedTuple, Dict, Any, Optional, List
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
    directory -- the directory in which the calculation is being run (used to enforce all path settings to be absolute paths)
    physicals -- list of keywords that have accompanying units from input
    '''

    # Need to provide these here to allow copy.deepcopy to perform the checks in __getattr__
    valid: List[str] = []
    data: Dict[str, Any] = {}
    are_paths: List[str] = []
    to_not_parse: List[str] = []
    physicals: List[str] = []

    def __init__(self, valid: List[str], defaults: Dict[str, Any] = {}, are_paths: List[str] = [], to_not_parse: List[str] = [], directory='', physicals: List[str] = [], **kwargs):
        super().__init__()
        self.valid = valid + self._other_valid_keywords
        self.defaults = defaults
        self.are_paths = are_paths
        self.to_not_parse = to_not_parse
        self.directory = directory
        self.physicals = physicals
        self.update(**defaults)
        try:
            self.update(**kwargs)
        except:
            raise ValueError()

    def __getattr__(self, name):
        if name != 'valid' and self.is_valid(name):
            return self.data.get(name, None)
        else:
            try:
                super().__getattr__(name)
            except AttributeError:
                raise AttributeError(name)
        # if key == 'directory':
        #     return self._directory
        # elif key == 'to_not_parse':
        #     if not '_to_not_parse' in self.__dict__:
        #         self.__dict__['_to_not_parse'] = set(self.are_paths)
        #     return self.__dict__['_to_not_parse']
        # elif key in ['attributes'] + self.attributes:
        #     return self.__dict__[key]
        # elif key in self.valid:
        #     return self.data.get(key, None)
        # else:
        #     if key not in self.data:
        #         return dict.__getattribute__(self.data, key)
        #     return self.data[key]

    def __setattr__(self, name, value):
        if self.is_valid(name):
            self.data[name] = value
        else:
            super().__setattr__(name, value)
    # def __setattr__(self, key, value):
    #     if key == 'directory':
    #         self.__dict__['_directory'] = value
    #     elif key == 'to_not_parse':
    #         self.__dict__['_to_not_parse'] = self.to_not_parse.union(value)
    #     elif key in ['attributes'] + self.attributes:
    #         self.__dict__[key] = value
    #     else:
    #         self.data[key] = value

    def __getitem__(self, key: str):
        if key not in self.data:
            if key in self.defaults:
                self.data[key] = self.defaults[key]
            else:
                raise KeyError(key)
        return self.data[key]

    def __setitem__(self, key: str, value: Any):
        # Insisting that all values corresponding to paths are absolute and are Path objects
        if key in self.are_paths:
            if isinstance(value, str):
                value = Path(value)
            elif not isinstance(value, Path):
                raise ValueError(f'{key} must be either a string or a Path')
            value = value.resolve()

        # Parse any units provided
        if key in self.physicals:
            value = parse_physical(value)

        # Perform additional checks that the key and corresponding value is valid
        self._check_before_setitem(key, value)

        # Set the item
        super().__setitem__(key, value)

    def is_valid(self, name):
        # Check if a keyword is valid. This is a separate subroutine to allow child classes to overwrite it
        # e.g. QE calculators want to be able to set keywords such as Hubbard(i) where i is an arbitrary integer
        return name in self.valid

    def update(self, *args, **kwargs):
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

        expr = expr.replace('/', ' / ').replace('*', ' * ').split()
        for i, term in enumerate(expr):
            if term in ['*', '/']:
                continue
            elif all([c.isalpha() for c in term]):
                if term in kwargs:
                    expr[i] = kwargs[term]
                elif getattr(self, term, None) is None:
                    raise ValueError('Failed to parse ' + ''.join(map(str, expr)))
                else:
                    expr[i] = getattr(self, term)
            else:
                expr[i] = float(term)

        value = float(expr[0])
        for op, term in zip(expr[1::2], expr[2::2]):
            if op == '*':
                value *= float(term)
            elif op == '/':
                value /= float(term)
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
        return ['pseudopotentials', 'kpts', 'koffset', 'kpath']

    def todict(self):
        # Shallow copy
        dct = dict(self.__dict__)

        # Adding information required by the json decoder
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    def fromdict(self, dct):
        for k, v in dct.items():
            setattr(self, k, v)


class SettingsDictWithChecks(SettingsDict):
    def __init__(self, settings: List[Setting], **kwargs):
        self.settings = settings
        super().__init__(valid=[s.name for s in settings], defaults={
            s.name: s.default for s in settings if not s.default is None}, **kwargs)

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


kc_wann_defaults = {'outdir': 'TMP',
                    'kc_iverbosity': 1,
                    'kc_at_ks': False,
                    'homo_only': False,
                    'read_unitary_matrix': True,
                    'check_ks': True,
                    'have_empty': True,
                    'has_disentangle': True}