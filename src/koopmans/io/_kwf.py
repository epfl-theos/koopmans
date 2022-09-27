"""

kwf (Koopmans WorkFlow) I/O for koopmans

Written by Edward Linscott Mar 2021, largely modelled off ase.io.jsonio

"""

import inspect
import json
import os
from importlib import import_module
from pathlib import Path
from typing import TextIO, Union

from ase.io import jsonio as ase_json

import koopmans.workflows as workflows


class KoopmansEncoder(ase_json.MyEncoder):
    def default(self, obj) -> dict:
        if isinstance(obj, set):
            return {'__set__': list(obj)}
        elif isinstance(obj, Path):
            return {'__path__': os.path.relpath(obj, '.')}
        elif inspect.isclass(obj):
            return {'__class__': {'__name__': obj.__name__, '__module__': obj.__module__}}
        elif hasattr(obj, 'todict'):
            d = obj.todict()
            if '__koopmans_name__' in d:
                return d

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


def read_kwf(fd: TextIO):
    return decode(fd.read())


def write_kwf(obj: Union[workflows.Workflow, dict], fd: TextIO):
    if isinstance(obj, workflows.Workflow):
        use_relpath = obj.parameters.use_relative_paths
        obj.parameters.use_relative_paths = True
    fd.write(encode(obj))
    if isinstance(obj, workflows.Workflow):
        obj.parameters.use_relative_paths = use_relpath
