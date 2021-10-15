"""

Generic I/O functions for koopmans

Written by Edward Linscott Jan 2020

"""

from typing import Union, List
from pathlib import Path
from ._json import read_json, write_json
from ._kwf import read_kwf, write_kwf
from ._calculators import read_calculator
from koopmans.utils import chdir
from koopmans.calculators import CalculatorExt
from koopmans.workflows import Workflow


def read(filename: Union[str, Path, List[str], List[Path]], **kwargs) -> Union[Workflow, CalculatorExt]:
    if isinstance(filename, str):
        filename = Path(filename)
    elif isinstance(filename, list):
        filename = [Path(f) for f in filename]

    # Generic "read" function

    if isinstance(filename, Path) and filename.suffix == '.kwf':
        with open(filename, 'r') as fd:
            out = read_kwf(fd)
        return out
    elif isinstance(filename, Path) and filename.suffix == '.json':
        with chdir(filename.parent):
            with open(filename.name, 'r') as fd:
                out = read_json(fd, **kwargs)
        return out
    else:
        try:
            return read_calculator(filename)
        except ValueError:
            raise ValueError(f'Unrecognised file type for {filename}')


def write(obj: Workflow, filename: Union[str, Path]):
    if isinstance(filename, str):
        filename = Path(filename)

    # Generic "write" function

    if filename.suffix == '.kwf':
        with open(filename, 'w') as fd:
            write_kwf(obj, fd)
    elif filename.suffix == '.json':
        write_json(obj, filename)
    else:
        raise ValueError(f'Unrecognised file type for {filename}')
