"""

Generic I/O functions for koopmans

Written by Edward Linscott Jan 2020

"""

from typing import Union, List
from ._json import read_json, write_json
from ._kwf import read_kwf, write_kwf
from ._calculators import read_calculator
from koopmans.calculators import ExtendedCalculator
from koopmans.workflows import Workflow


def read(filename: Union[str, List[str]], **kwargs) -> Union[Workflow, ExtendedCalculator]:

    # Generic "read" function

    if isinstance(filename, str) and filename.endswith('kwf'):
        with open(filename, 'r') as fd:
            out = read_kwf(fd)
        return out
    elif isinstance(filename, str) and filename.endswith('.json'):
        return read_json(filename, **kwargs)
    else:
        try:
            return read_calculator(filename)
        except ValueError:
            raise ValueError(f'Unrecognised file type for {filename}')


def write(obj: Workflow, filename: str):

    # Generic "write" function

    if filename.endswith('kwf'):
        with open(filename, 'w') as fd:
            write_kwf(obj, fd)
    elif filename.endswith('.json'):
        write_json(obj, filename)
    else:
        raise ValueError(f'Unrecognised file type for {filename}')
