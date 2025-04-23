"""Generic I/O functions for koopmans."""

from pathlib import Path
from typing import Any, List, Union

from koopmans.utils import chdir

from ._calculators import read_calculator
from ._dill import read_pkl, write_pkl
from ._json import read_json, write_json


def read(filename: Union[str, Path, List[str], List[Path]], **kwargs) -> Any:
    """Read an object from file, with the behavior depending on the file extension."""
    if isinstance(filename, str):
        filename = Path(filename)
    elif isinstance(filename, list):
        filename = [Path(f) for f in filename]

    # Generic "read" function

    if isinstance(filename, Path) and filename.suffix == '.pkl':
        with chdir(filename.parent):
            out = read_pkl(filename.name)
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
            raise ValueError(f'Unrecognized file type for `{filename}`')


def write(obj: Any, filename: Union[str, Path]):
    """Write an object to file. Behaves differently depending on the file extension."""
    if isinstance(filename, str):
        filename = Path(filename)

    # Generic "write" function

    if filename.suffix == '.pkl':
        write_pkl(obj, filename)
    elif filename.suffix == '.json':
        write_json(obj, filename)
    else:
        raise ValueError(f'Unrecognized file type for `{filename}`')
