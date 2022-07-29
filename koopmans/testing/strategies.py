from multiprocessing.sharedctypes import Value
from typing import List

from hypothesis import given, settings, strategies

from ase.cell import Cell
from ase.dft.kpoints import BandPath
from ase.lattice import (BCC, BCT, CUB, FCC, HEX, MCL, MCLC, ORC, ORCC, ORCF,
                         ORCI, RHL, TET, TRI, BravaisLattice,
                         UnconventionalLattice, tri_angles_explanation)
from koopmans.kpoints import Kpoints


@strategies.composite
def ase_cells(draw) -> Cell:
    # Randomly select a bravais lattice class
    lat_classes: List[BravaisLattice] = [CUB, BCC, FCC, TET, BCT, ORC, ORCF, ORCI, ORCC, HEX, RHL, MCL, MCLC, TRI]
    index: int = draw(strategies.integers(min_value=0, max_value=len(lat_classes) - 1))
    lat_class = lat_classes[index]

    # Define the strategies for drawing lengths and angles
    def lengths(min_value: float = 1.0):
        return strategies.decimals(min_value=min_value, max_value=min_value + 5.0, places=3).map(float)

    # For the lengths, we will impose a < b < c
    delta = 0.1
    a = draw(lengths())
    b = draw(lengths(min_value=a + delta))
    c = draw(lengths(min_value=b + delta))
    '''
    For the angles, we will impose...
    - alpha < beta < gamma (by convention)
    - alpha + beta + gamma <= 355 (strictly, we must have alpha + beta + gamma < 360)
    - alpha + beta >= gamma + 4 (strictly, we must have alpha + beta > gamma)
    - all angles >= 5 (just so the cell is not too problematic)
    '''
    alphas = strategies.integers(min_value=5, max_value=352 // 3).map(float)
    attempts = 0
    while attempts < 1000:
        alpha = draw(alphas)
        betas = strategies.integers(min_value=alpha + 1, max_value=(354 - alpha) // 2).map(float)
        beta = draw(betas)
        gammas = strategies.integers(min_value=beta + 1,
                                     max_value=min(355 - alpha - beta, alpha + beta - 4)).map(float)
        gamma = draw(gammas)

        parameters = {"a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma}

        # Return a lattice with the relevant parameters
        try:
            lat = lat_class(**{k: v for k, v in parameters.items() if k in lat_class.parameters})
            break
        except UnconventionalLattice:
            attempts += 1
    assert attempts < 1000, 'Failed to generate a set of valid cell parameters'

    return lat.tocell()


@strategies.composite
def bandpaths(draw) -> BandPath:
    cell = draw(ase_cells())
    return cell.bandpath(eps=1e-12)


_grids = strategies.lists(strategies.integers(), min_size=3, max_size=3)
_offsets = strategies.lists(strategies.integers(min_value=0, max_value=1), min_size=3, max_size=3)


@strategies.composite
def _kpoints_via_bandpath(draw) -> Kpoints:
    gamma_only = draw(strategies.booleans())
    if gamma_only:
        grid = None
        offset = None
    else:
        grid = draw(_grids)
        offset = draw(_offsets)
    path = draw(bandpaths())
    return Kpoints(grid, offset, path, gamma_only)


@strategies.composite
def _kpoints_via_pathstr(draw) -> Kpoints:
    gamma_only = draw(strategies.booleans())
    if gamma_only:
        grid = None
        offset = None
    else:
        grid = draw(_grids)
        offset = draw(_offsets)
    cell = draw(ase_cells())
    density = draw(strategies.floats(min_value=1.0, max_value=50.0))
    path = cell.bandpath().path
    return Kpoints(grid, offset, path, gamma_only, cell, density)


kpoints = _kpoints_via_bandpath() | _kpoints_via_pathstr()
