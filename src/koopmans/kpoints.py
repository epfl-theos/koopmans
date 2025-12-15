"""Defines a "kpoints" class that stores information related to k-points.

Written by Edward Linscott, July 2022
"""


from __future__ import annotations

import copy
from typing import (Annotated, Any, ClassVar, Dict, List, Literal, Optional,
                    Union)

import numpy as np
from ase.cell import Cell
from ase_koopmans.dft.kpoints import (BandPath, kpoint_convert,
                                      parse_path_string,
                                      resolve_kpt_path_string)
from pydantic import Field, field_serializer, field_validator

from koopmans.pydantic import BaseModel
from koopmans.utils import remove_numpy_from_dict

ZeroOrOne = Literal[0, 1]
Fractional = Annotated[float, Field(gt=-1, le=1)]


class GammaOnlyKpoints(BaseModel):
    gamma_only: Literal[True] = True
    grid: ClassVar[tuple[int, int, int]] = (1, 1, 1)
    offset: ClassVar[tuple[ZeroOrOne, ZeroOrOne, ZeroOrOne]] = (0, 0, 0)
    offset_nscf: ClassVar[tuple[Fractional, Fractional, Fractional]] = (0.0, 0.0, 0.0)
    path: ClassVar[None] = None


class Kpoints(BaseModel):
    """A class for storing k-point information.

    Attributes
    ----------
    grid : List[int]
        a list of three integers specifying the shape of the regular grid of k-points
    offset : List[int]
        a list of three integers, either zero or one. If one, the regular k-point grid is offset by half a grid step
        in that dimension
    offset_nscf : Union[List[Union[int, float]]], int, float]
        a list of three numbers, within the interval (-1, 1].
        It has the same meaning of offset, but it is applied only to the grid of PWnscf calculations. Additionally it
        allows to provide explicitly any fractional offset.
    path : ase.dft.kpoints.BandPath
        an ASE ``BandPath`` object specifying the k-path as defined by the special points of the Bravais
        lattice
    gamma_only : bool
        True if the calculation is only sampling the gamma point
    """

    grid: tuple[int, int, int]
    offset: tuple[ZeroOrOne, ZeroOrOne, ZeroOrOne] = (0, 0, 0)
    offset_nscf: tuple[float, float, float] = (0.0, 0.0, 0.0)
    path: BandPath
    gamma_only: Literal[False] = False

    @field_serializer("path")
    def serialize_path(self, v: BandPath) -> dict[str, Any]:
        dct = v.todict()
        dct.pop('labelseq')
        dct['cell'] = dct['cell'].todict()
        return remove_numpy_from_dict(dct)

    @field_validator("path", mode="before")
    @classmethod
    def coerce_path(cls, value: Any) -> BandPath:
        if isinstance(value, dict):
            value['cell'] = Cell(**value['cell'])
            return BandPath(**value)
        return value

    @classmethod
    def from_path_as_string(cls, path: str | None, cell: Cell, density: float = 10.0, **kwargs) -> Kpoints:
        """Create a Kpoints object using a string for the path instead of a `BandPath` object.

        cell : ase.cell.Cell
            the simulation cell
        density : float
            k-points per inverse Angstrom along the k-path
        """
        bandpath = convert_kpath_str_to_bandpath(path, cell, density)
        return cls(path=bandpath, **kwargs)

    def __eq__(self, other: object) -> bool:
        raise NotImplementedError()
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
    """Convert a string of high-symmetry k-points to a BandPath object."""
    special_points: Dict[str, np.ndarray] = cell.bandpath(eps=1e-10).special_points
    special_points_on_path = set([x for y in parse_path_string(path) for x in y])
    for sp in special_points_on_path:
        if sp not in special_points.keys():
            raise KeyError(
                f'The path provided to convert_kpath_str_to_bandpath contains a special point ({sp}) '
                f'that is not in the set of special points for this cell ({", ".join(special_points.keys())})')
    bp = BandPath(cell=cell, path=path, special_points=special_points)
    if len(path) > 1:
        bp = bp.interpolate(density=density)
    return bp


def kpath_length(path: BandPath) -> float:
    """Calculate the length of a k-path."""
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
    """Convert a BandPath object to a dictionary."""
    dct = {}
    dct['path'] = path.path
    dct['cell'] = path.cell.todict()
    if len(path.path) > 1:
        attempts = 0
        density_guess = len(path.kpts) / kpath_length(path)
        density_max = None
        density_min = None

        # Super-dumb bisection approach to obtain a density that gives the correct nunber of kpoints.
        # We have to resort to this because the length of paths produced by BandPath().interpolate(density=...)
        # is unreliable
        while attempts < 100:
            dct['density'] = density_guess
            new_path = dict_to_kpath(dct)
            if len(new_path.kpts) == len(path.kpts):
                break
            elif len(new_path.kpts) < len(path.kpts):
                density_min = density_guess
                density_guess = density_min + 1.0
            else:
                density_max = density_guess
                density_guess = density_max - 1.0
            if density_min and density_max:
                density_guess = (density_max + density_min) / 2
            attempts += 1
        if attempts >= 100:
            raise ValueError('Search failed')
    return dct


def dict_to_kpath(dct: Dict[str, Any]) -> BandPath:
    """Convert a dictionary to a BandPath object."""
    dct = copy.deepcopy(dct)
    density = dct.pop('density', None)
    cell = Cell(dct.pop('cell')['array'])
    return BandPath(cell=cell, special_points=cell.bandpath(eps=1e-10).special_points,
                    **dct).interpolate(density=density)
