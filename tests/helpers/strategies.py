"""Hypothesis strategies for generating objects (cells, kpoints, etc) for testing purposes."""

from typing import Callable, List, Optional

import numpy as np
from ase_koopmans.cell import Cell
from ase_koopmans.dft.kpoints import BandPath
from ase_koopmans.lattice import (BCC, BCT, CUB, FCC, HEX, MCL, MCLC, ORC,
                                  ORCC, ORCF, ORCI, RHL, TET, TRI,
                                  BravaisLattice, UnconventionalLattice)
from hypothesis import assume, note
from hypothesis.strategies import (booleans, composite, decimals, floats,
                                   integers, lists)


@composite
def ase_cells(draw: Callable) -> Cell:
    """Generate ASE cells."""
    # Define the strategies for drawing lengths and angles
    def lengths(min_value: float = 1.0):
        return decimals(min_value=min_value, max_value=min_value + 5.0, places=3).map(float)

    # For the lengths, we will impose a < b < c
    delta = 0.1
    a: float = draw(lengths())
    b: float = draw(lengths(min_value=a + delta))
    c: float = draw(lengths(min_value=b + delta))

    # For the angles, we will impose...
    # - alpha < beta < gamma (by convention)
    # - alpha + beta + gamma <= 355 (strictly, we must have alpha + beta + gamma < 360)
    # - alpha + beta >= gamma + 4 (strictly, we must have alpha + beta > gamma)
    # - all angles >= 5 (just so the cell is not too problematic)
    alphas = integers(min_value=5, max_value=352 // 3)
    alpha: int = draw(alphas)
    note(f'Drew alpha = {alpha}')
    betas = integers(min_value=alpha + 1, max_value=(354 - alpha) // 2)
    beta: int = draw(betas)
    gammas = integers(min_value=beta + 1,
                      max_value=min(355 - alpha - beta, alpha + beta - 4))
    gamma: int = draw(gammas)
    parameters = {"a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma}

    # Randomly select a Bravais lattice class
    lat_classes: List[BravaisLattice] = [CUB, BCC, FCC, TET, BCT, ORC, ORCF, ORCI, ORCC, HEX, RHL, MCL, MCLC, TRI]
    index: int = draw(integers(min_value=0, max_value=len(lat_classes) - 1))
    lat_class = lat_classes[index]

    # Create a Bravais lattice instance with the relevant parameters
    lat: Optional[BravaisLattice] = None
    try:
        lat = lat_class(**{k: v for k, v in parameters.items() if k in lat_class.parameters})
    except UnconventionalLattice:
        pass
    # Throw away any cases which gave an UnconventionalLattice
    assume(lat is not None)

    # Convert the lattice to a Cell
    cell = lat.tocell()

    # Filter out any examples which are not Niggli-reduced
    new_cell, op = cell.niggli_reduce()
    assume(np.all(op == np.identity(3, dtype=int)))

    # Log these parameters so hypothesis can report them upon a failure
    parameter_strs = [f'{k}={v}' for k, v in parameters.items() if k in lat_class.parameters]
    note('Generated a {lat_class} cell with ' + ', '.join(parameter_strs))

    return new_cell


@composite
def bandpaths(draw: Callable) -> BandPath:
    """Generate band paths."""
    cell = draw(ase_cells())
    return cell.bandpath(eps=1e-12)


_grids = lists(integers(), min_size=3, max_size=3)
_offsets = lists(integers(min_value=0, max_value=1), min_size=3, max_size=3)
_offsets_nscf = lists(floats(min_value=0, max_value=1), min_size=3, max_size=3)


@composite
def _kpoints_via_bandpath(draw):
    """Generate k-points via a `BandPath`."""
    from koopmans.kpoints import Kpoints
    gamma_only = draw(booleans())
    if gamma_only:
        grid = None
        offset = None
        offset_nscf = None
    else:
        grid = draw(_grids)
        offset = draw(_offsets)
        offset_nscf = draw(_offsets_nscf)
    path = draw(bandpaths())
    return Kpoints(grid, offset, offset_nscf, path, gamma_only)


@composite
def _kpoints_via_pathstr(draw):
    """Generate k-points via a path string."""
    from koopmans.kpoints import Kpoints
    gamma_only = draw(booleans())
    if gamma_only:
        grid = None
        offset = None
        offset_nscf = None
    else:
        grid = draw(_grids)
        offset = draw(_offsets)
        offset_nscf = draw(_offsets_nscf)
    cell = draw(ase_cells())
    density = draw(floats(min_value=1.0, max_value=50.0))
    path = cell.bandpath().path
    return Kpoints(grid, offset, offset_nscf, path, gamma_only, cell, density)


kpoints = _kpoints_via_bandpath() | _kpoints_via_pathstr()
