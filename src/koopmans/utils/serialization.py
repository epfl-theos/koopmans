"""

Encoder and decoder for koopmans objects.

"""

import functools
import inspect
import json
import os
from importlib import import_module
from pathlib import Path

from ase.io import jsonio as ase_json

from koopmans import utils


def convert_keys_to_strings(dct):
    if isinstance(dct, dict):
        return {str(key): convert_keys_to_strings(value) for key, value in dct.items()}
    elif isinstance(dct, list):
        return [convert_keys_to_strings(value) for value in dct]
    else:
        return dct


class KoopmansEncoder(ase_json.MyEncoder):
    def default(self, obj) -> dict:
        if isinstance(obj, set):
            return {'__set__': list(obj)}
        elif isinstance(obj, Path):
            return {'__path__': os.path.relpath(obj, '.')}
        elif inspect.isclass(obj):
            return {'__class__': {'__name__': obj.__name__, '__module__': obj.__module__}}
        elif isinstance(obj, functools.partial):
            return {'__partial__': {'func': obj.func, 'args': obj.args, 'keywords': obj.keywords}}
        elif inspect.isfunction(obj):
            module = import_module(obj.__module__)
            if hasattr(module, obj.__name__):
                return {'__function__': {'__name__': obj.__name__, '__module__': obj.__module__}}
            else:
                return {'__function__': {}}
        elif hasattr(obj, 'todict'):
            d = obj.todict()
            if '__koopmans_name__' in d:
                return convert_keys_to_strings(d)

        # If none of the above, use ASE's encoder
        return super().default(obj)


encode = KoopmansEncoder(indent=1).encode


def object_hook(dct):
    if '__koopmans_name__' in dct:
        return create_koopmans_object(dct)
    elif '__set__' in dct:
        return set(dct['__set__'])
    elif '__path__' in dct:
        return Path(dct['__path__'])
    elif '__class__' in dct:
        subdct = dct['__class__']
        module = import_module(subdct['__module__'])
        return getattr(module, subdct['__name__'])
    elif '__partial__' in dct:
        subdct = dct['__partial__']
        return functools.partial(subdct['func'], **subdct['keywords'])
    elif '__function__' in dct:
        subdct = dct['__function__']
        if '__module__' in subdct:
            module = import_module(subdct['__module__'])
            if hasattr(module, subdct['__name__']):
                return getattr(module, subdct['__name__'])
        utils.warn('Function was not able to be serialized')
        return f'<placeholder for function {subdct["__name__"]} lost during serialization>'
    else:
        # Patching bug in ASE where allocating an np.empty(dtype=str) will assume a particular length for each
        # string. dtype=object allows for individual strings to be different lengths
        if '__ndarray__' in dct:
            dtype = dct['__ndarray__'][1]
            if 'str' in dtype:
                dct['__ndarray__'][1] = object

        # Making it possible for "None" to be a key in a dictionary
        if 'null' in dct:
            dct[None] = dct.pop('null')

        return ase_json.object_hook(dct)


def create_koopmans_object(dct: dict):
    # Load object class corresponding to this dictionary
    name = dct.pop('__koopmans_name__')
    module = import_module(dct.pop('__koopmans_module__'))
    objclass = getattr(module, name)

    # Reconstruct the class from the dictionary
    return objclass.fromdict(dct)


decode = json.JSONDecoder(object_hook=object_hook).decode
