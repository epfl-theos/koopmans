"""

Generic I/O functions for koopmans

Written by Edward Linscott Jan 2020

"""

from ._json import read_json, write_json
from ._kwf import read_kwf, write_kwf
from koopmans.workflows.generic import Workflow


def read(filename: str, **kwargs) -> Workflow:

    # Generic "read" function

    if filename.endswith('kwf'):
        with open(filename, 'r') as fd:
            out = read_kwf(fd)
        return out
    elif filename.endswith('.json'):
        return read_json(filename, **kwargs)
    else:
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
