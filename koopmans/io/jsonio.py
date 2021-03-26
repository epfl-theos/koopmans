"""

JSON I/O for python_KI

Written by Edward Linscott Mar 2021, largely modelled off ase.io.jsonio

"""


from importlib import import_module
import inspect
import json
from ase.io import jsonio as ase_json
from ase.calculators.calculator import Calculator
from koopmans.calculators.generic import GenericCalc
from koopmans.workflows.generic import Workflow
from koopmans.workflows.singlepoint import SinglepointWorkflow


class KoopmansEncoder(ase_json.MyEncoder):
    def default(self, obj):

        if isinstance(obj, set):
            return {'__set__': list(obj)}
        elif isinstance(obj, Calculator):
            # ASE only stores the calculator parameters, with Atoms being the more fundamental object
            # Because we store calculators as the primary object, we need to make sure the atoms are also stored
            d = {'__calculator__': obj.todict(),
                 '__name__': obj.__class__.__name__,
                 '__module__': obj.__class__.__module__,
                 '__results__': obj.results,
                 '__atoms__': super().default(obj.atoms)}
            return d
        elif isinstance(obj, (GenericCalc, Workflow)):
            d = obj.todict()
            if d.get('__koopmans_name__', None):
                return d
        elif inspect.isclass(obj):
            return {'__class__': {'__name__': obj.__name__, '__module__': obj.__module__}}
        return super().default(obj)


encode = KoopmansEncoder().encode


def object_hook(dct):
    if '__koopmans_name__' in dct:
        dct = ase_json.numpyfy(dct)
        return create_koopmans_object(dct)
    elif '__calculator__' in dct:
        return create_ase_calculator(dct)
    elif '__set__' in dct:
        return set(dct['__set__'])
    elif '__class__' in dct:
        subdct = dct['__class__']
        module = import_module(subdct['__module__'])
        return getattr(module, subdct['__name__'])
    else:
        # Patching bug in ASE
        if '__ndarray__' in dct:
            dtype = dct['__ndarray__'][1]
            if 'str' in dtype:
                dct['__ndarray__'][1] = 'str'

        return ase_json.object_hook(dct)


def create_ase_calculator(dct):
    module = import_module(dct['__module__'])
    calc_class = getattr(module, dct['__name__'])
    calc = calc_class()
    calc.atoms = ase_json.object_hook(dct.pop('__atoms__'))
    calc.atoms.calc = calc
    calc.parameters = dct['__calculator__']
    calc.results = dct['__results__']
    return calc


def create_koopmans_object(dct):
    # Load object class corresponding to this dictionary
    name = dct.pop('__koopmans_name__')
    module = import_module(dct.pop('__koopmans_module__'))
    objclass = getattr(module, name)

    # Reconstruct the class from the dictionary
    return objclass(dct=dct)


decode = json.JSONDecoder(object_hook=object_hook).decode


def read_json(fd):
    return decode(fd.read())


def write_json(fd, obj):
    fd.write(encode(obj))
