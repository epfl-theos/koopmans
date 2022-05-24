'''

Miscellaneous functions for koopmans.utils

Written by Edward Linscott May 2020

'''


from typing import List, Union, Generator, Iterable, Any, Dict
import numpy as np
from ase.cell import Cell
from ase.dft.kpoints import BandPath


def flatten(l: Union[List[Any], Iterable[Any]]) -> Generator[Any, None, None]:
    # Converts a list of any kind of object (numbers, arrays, lists, strings, ecc.)
    # to a generator
    for item in l:
        if isinstance(item, Iterable) and not isinstance(item, str):
            for x in flatten(item):
                yield x
        else:
            yield item


def convert_kpath_str_to_bandpath(path: str, cell: Cell, density: int = 10) -> BandPath:
    npoints = density * len(path) - density + 1 - (3 * density - 1) * path.count(',')
    special_points: Dict[str, np.ndarray] = cell.bandpath().special_points
    return BandPath(cell=cell, path=path, special_points=special_points).interpolate(npoints=npoints)
