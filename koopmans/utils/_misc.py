'''

Miscellaneous functions for koopmans.utils

Written by Edward Linscott May 2020

'''


from typing import List, Union, Generator, Iterable, Any, Dict, Optional
import numpy as np
from ase.cell import Cell
from ase.dft.kpoints import BandPath, resolve_kpt_path_string, kpoint_convert


def flatten(l: Union[List[Any], Iterable[Any]]) -> Generator[Any, None, None]:
    # Converts a list of any kind of object (numbers, arrays, lists, strings, ecc.)
    # to a generator
    for item in l:
        if isinstance(item, Iterable) and not isinstance(item, str):
            for x in flatten(item):
                yield x
        else:
            yield item


def convert_kpath_str_to_bandpath(path: str, cell: Cell, density: Optional[int] = None) -> BandPath:
    special_points: Dict[str, np.ndarray] = cell.bandpath().special_points
    return BandPath(cell=cell, path=path, special_points=special_points).interpolate(density=density)


def kpath_length(path: BandPath, cell: Cell) -> float:
    _, paths = resolve_kpt_path_string(path.path, path.special_points)
    points = np.concatenate(paths)
    dists = points[1:] - points[:-1]
    lengths = [np.linalg.norm(d) for d in kpoint_convert(cell, skpts_kc=dists)]

    i = 0
    for path in paths[:-1]:
        i += len(path)
        lengths[i - 1] = 0

    return sum(lengths)


def kpath_to_dict(path: BandPath, cell: Cell) -> Dict[str, Any]:
    dct = {}
    dct['kpath'] = path.path
    dct['kpath_density'] = len(path.kpts) / kpath_length(path, cell)
    return dct
