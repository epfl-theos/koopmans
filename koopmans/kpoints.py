"""
Defines a "kpoints" class that stores information related to k-points

Written by Edward Linscott, July 2022
"""


from dataclasses import dataclass
from typing import Dict, List, Optional, Union

from ase.cell import Cell
from ase.dft.kpoints import BandPath
from koopmans import utils


@dataclass
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

    def __init__(self, grid: Optional[List[int]] = None, offset: Optional[List[int]] = None, path: Optional[Union[str, BandPath]] = None, gamma_only: bool = False, cell: Optional[Cell] = None, density: float = 10.0):
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
        if value is not None:
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
        self._path = value

    def set_path(self, path: Union[str, BandPath], cell: Optional[Cell] = None, density: float = 10.0):
        # A function for setting self.path with a string
        if isinstance(path, str):
            # Convert the string to a BandPath object
            if cell is None:
                raise ValueError(
                    'To set the path of a Kpoints object with a string, please provide a value for "cell"')
            path = utils.convert_kpath_str_to_bandpath(path, cell, density)
        self.path = path

    def tojson(self) -> Dict[str, Union[str, bool, List[int], BandPath]]:
        return {k.lstrip('_'): v for k, v in self.__dict__.items()}

    def todict(self) -> Dict[str, Union[str, bool, List[int], BandPath]]:
        dct = self.tojson()

        # Adding information required by the json decoder
        dct['__koopmans_name__'] = self.__class__.__name__
        dct['__koopmans_module__'] = self.__class__.__module__
        return dct
