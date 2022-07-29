"""
Defines a "kpoints" class that stores information related to k-points

Written by Edward Linscott, July 2022
"""


from __future__ import annotations

from typing import Any, Dict, List, Optional, Union

import numpy as np

from ase.cell import Cell
from ase.dft.kpoints import BandPath, kpoint_convert, resolve_kpt_path_string


class Kpoints:
    """
    A class for storing k-point information

    Attributes
    ----------
    grid : List[int]
        a list of three integers specifying the shape of the regular grid of k-points
    offset : List[int]
        a list of three integers, either zero or one. If one, the regular k-point grid is offset by half a grid step
        in that dimension
    path : ase.dft.kpoints.BandPath
        an ASE ``BandPath`` object specifying the k-path as defined by the special points of the Bravais
        lattice
    gamma_only : bool
        True if the calculation is only sampling the gamma point
    """

    _grid: Optional[List[int]]
    _offset: Optional[List[int]]
    _path: Optional[BandPath]
    gamma_only: bool

    def __init__(self, grid: Optional[List[int]] = [1, 1, 1], offset: Optional[List[int]] = [0, 0, 0], path: Optional[Union[str, BandPath]] = None, gamma_only: bool = False, cell: Optional[Cell] = None, density: float = 10.0):
        """
        Initialize a Kpoint object. The path can be initialized using a string, but if so the additional requirements
        density and cell are required, where...

        cell : ase.cell.Cell
            the simulation cell
        density : float
            k-points per inverse Angstrom along the k-path
        """
        self.gamma_only = gamma_only

        if gamma_only:
            if grid:
                raise ValueError(f'gamma_only = {gamma_only} and grid != None are incompatible')
            self.grid = None
            self.offset = [0, 0, 0]
        else:
            self.grid = grid
            self.offset = offset

        self.set_path(path, cell, density)

    def __repr__(self):
        return 'Kpoints(' + ', '.join(f'{k.lstrip("_")}={v}' for k, v in self.__dict__.items() if v is not None) + ')'

    @property
    def grid(self) -> Optional[List[int]]:
        return self._grid

    @grid.setter
    def grid(self, value: Optional[List[int]]):
        if value is not None:
            if len(value) != 3:
                raise ValueError('"grid" must be a list of three integers')
        self._grid = value

    @property
    def offset(self) -> Optional[List[int]]:
        return self._offset

    @offset.setter
    def offset(self, value: Optional[List[int]]):
        if isinstance(value, list):
            if len(value) != 3:
                raise ValueError('"offset" must be a list of three integers')
            if any([x not in [0, 1] for x in value]):
                raise ValueError('"offset" must only contain either 0 or 1s')
        self._offset = value

    @property
    def path(self) -> Optional[BandPath]:
        return self._path

    @path.setter
    def path(self, value: Optional[BandPath]):
        if isinstance(value, str):
            raise ValueError('To set Kpoints.path with a string, please use Kpoints.set_path()')
        self._path = value

    def set_path(self, path: Optional[Union[str, BandPath]], cell: Optional[Cell] = None, density: float = 10.0):
        # A function for setting self.path with a string
        if isinstance(path, str):
            # Convert the string to a BandPath object
            if cell is None:
                raise ValueError(
                    'To set the path of a Kpoints object with a string, please provide a value for "cell"')
            path = convert_kpath_str_to_bandpath(path, cell, density)
        self.path = path

    def tojson(self) -> Dict[str, Union[str, bool, List[int], Dict[str, Any]]]:
        dct: Dict[str, Union[str, bool, List[int], Dict[str, Any]]] = {}
        for k, v in self.__dict__.items():
            k = k.lstrip('_')
            if isinstance(v, BandPath):
                # Store the path and the density, but not the cell because that is stored elsewhere
                dct.update(**kpath_to_dict(v))
                dct.pop('cell')
            else:
                dct[k] = v
        return dct

    def todict(self) -> Dict[str, Union[str, bool, List[int], BandPath]]:
        dct = self.tojson()

        # We also need the cell
        assert isinstance(self.path, BandPath)
        dct['cell'] = self.path.cell

        # Adding information required by the json decoder
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct

    @classmethod
    def fromdict(cls, dct) -> Kpoints:
        # Remove name and module if they're there
        dct.pop('__koopmans_name__', None)
        dct.pop('__koopmans_module__', None)

        return cls(**dct)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Kpoints):
            return False
        for k, v in self.__dict__.items():
            vother = getattr(other, k, None)
            if isinstance(v, np.ndarray):
                assert isinstance(vother, np.ndarray)
                if v.shape != vother.shape:
                    return False
                if not np.allclose(v, vother):
                    return False
            else:
                if v != vother:
                    return False
        return True


def convert_kpath_str_to_bandpath(path: str, cell: Cell, density: Optional[float] = None) -> BandPath:
    special_points: Dict[str, np.ndarray] = cell.bandpath().special_points
    bp = BandPath(cell=cell, path=path, special_points=special_points)
    if len(path) > 1:
        bp = bp.interpolate(density=density)
    return bp


def kpath_length(path: BandPath) -> float:
    _, paths = resolve_kpt_path_string(path.path, path.special_points)
    points = np.concatenate(paths)
    dists = points[1:] - points[:-1]
    lengths: List[float] = [float(np.linalg.norm(d)) for d in kpoint_convert(path.cell, skpts_kc=dists)]

    i = 0
    for path in paths[:-1]:
        i += len(path)
        lengths[i - 1] = 0.0

    return np.sum(lengths)


def kpath_to_dict(path: BandPath) -> Dict[str, Any]:
    dct = {}
    dct['path'] = path.path
    dct['cell'] = path.cell.todict()
    if len(path.path) > 1:
        dct['density'] = len(path.kpts) / kpath_length(path)
    return dct


def dict_to_kpath(dct: Dict[str, Any]) -> BandPath:
    density = dct.pop('density', None)
    cell = Cell(dct.pop('cell')['array'])
    return BandPath(cell=cell, special_points=cell.bandpath().special_points, **dct).interpolate(density=density)
