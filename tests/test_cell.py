import numpy as np
import pytest
from hypothesis import given, settings, strategies

from ase.cell import Cell
from ase.lattice import (BCC, BCT, CUB, FCC, HEX, MCL, MCLC, ORC, ORCC, ORCF,
                         ORCI, RHL, TET, TRI, UnconventionalLattice)
from koopmans.cell import cell_to_parameters, parameters_to_cell


@strategies.composite
def ase_cells(draw) -> Cell:
    # Randomly select a bravais lattice class
    lat_classes = [CUB, BCC, FCC, TET, BCT, ORC, ORCF, ORCI, ORCC, HEX, RHL, MCL, MCLC, TRI]
    index = draw(strategies.integers(min_value=0, max_value=len(lat_classes) - 1))
    lat_class = lat_classes[index]

    # Define the strategies for drawing lengths and angles
    lengths = strategies.floats(min_value=1.0, max_value=5.0)
    # For the angles, we will impose...
    # - alpha <= beta <= gamma (by convention)
    # - alpha + beta + gamma < 355.0 (alpha + beta + gamma > 360.0 is impossible)
    # - alpha + beta >= gamma + 4.0 (alpha + beta < gamma is impossible)
    alphas = strategies.floats(min_value=5.0, max_value=355.0 / 3)
    attempts = 0
    while attempts < 1000:
        alpha = draw(alphas)
        betas = strategies.floats(min_value=alpha, max_value=(355.0 - alpha) / 2)
        beta = draw(betas)
        gammas = strategies.floats(min_value=max(alpha, beta),
                                   max_value=min(355.0 - alpha - beta, alpha + beta - 4.0))
        gamma = draw(gammas)

        parameters = {"a": draw(lengths), "b": draw(lengths), "c": draw(lengths),
                      "alpha": alpha, "beta": beta, "gamma": gamma}

        # Return a lattice with the relevant parameters
        try:
            lat = lat_class(**{k: v for k, v in parameters.items() if k in lat_class.parameters})
            break
        except UnconventionalLattice:
            attempts += 1
    assert attempts < 1000, 'Failed to generate a set of valid cell parameters'

    return lat.tocell()


# Reference list of ASE Bravais lattice classes
bravais_lattices = {1: CUB, 2: FCC, 3: BCC, 4: HEX, 5: RHL, 6: TET, 7: BCT,
                    8: ORC, 9: ORCC, 10: ORCF, 11: ORCI, 12: MCL, 13: MCLC, 14: TRI}


@given(cell=ase_cells())
@settings(max_examples=100, report_multiple_bugs=False)
def test_cell_to_parameters(cell):
    try:
        params = cell_to_parameters(cell)
        bl_class = bravais_lattices[params['ibrav']]
        assert isinstance(cell.get_bravais_lattice(), bl_class)
    except ValueError as e:
        if 'You have provided a cell that is not Niggli-reduced' not in str(e):
            raise ValueError(e)


@given(cell=ase_cells())
@settings(max_examples=100, report_multiple_bugs=False)
def test_roundtrip(cell):
    try:
        params = cell_to_parameters(cell)
        new_cell = parameters_to_cell(**params)

        # Check is volume-conserving
        assert abs(cell.volume - new_cell.volume) < 1e-8

        # Check is idempotent
        new_params = cell_to_parameters(cell)
        assert params == new_params
    except ValueError as e:
        if 'You have provided a cell that is not Niggli-reduced' not in str(e):
            raise ValueError(e)
