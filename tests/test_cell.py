import numpy as np
import pytest
from hypothesis import given, settings, strategies

from ase.cell import Cell
from ase.lattice import (BCC, BCT, CUB, FCC, HEX, MCL, MCLC, ORC, ORCC, ORCF,
                         ORCI, RHL, TET, TRI, UnconventionalLattice)
from koopmans.cell import cell_to_parameters, parameters_to_cell
from koopmans.testing import strategies as kst

# Reference list of ASE Bravais lattice classes
bravais_lattices = {1: CUB, 2: FCC, 3: BCC, 4: HEX, 5: RHL, 6: TET, 7: BCT,
                    8: ORC, 9: ORCC, 10: ORCF, 11: ORCI, 12: MCL, 13: MCLC, 14: TRI}


@given(cell=kst.ase_cells())
@settings(max_examples=100, report_multiple_bugs=False)
def test_cell_to_parameters(cell: Cell):
    try:
        params = cell_to_parameters(cell)
        bl_class = bravais_lattices[params['ibrav']]
        assert isinstance(cell.get_bravais_lattice(), bl_class)
    except ValueError as e:
        if 'You have provided a cell that is not Niggli-reduced' not in str(e):
            raise


@given(cell=kst.ase_cells())
@settings(max_examples=100, report_multiple_bugs=False)
def test_roundtrip_cell_parameters(cell: Cell):
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
            raise
