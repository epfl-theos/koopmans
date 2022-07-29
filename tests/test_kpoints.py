import numpy as np
from hypothesis import given, settings

from ase.dft.kpoints import BandPath
from koopmans.kpoints import Kpoints, dict_to_kpath, kpath_to_dict
from koopmans.testing import strategies
from koopmans.testing import strategies as kst


@given(bp=kst.bandpaths())
def test_roundtrip_cell_bandpath(bp: BandPath):
    bp_roundtrip = bp.cell.bandpath()
    assert bp_roundtrip.special_points.keys() == bp.special_points.keys()


@given(bp=kst.bandpaths())
@settings(report_multiple_bugs=False)
def test_roundtrip_kpath_dict(bp: BandPath):
    bp_roundtrip = dict_to_kpath(kpath_to_dict(bp))
    assert np.allclose(bp_roundtrip.cell, bp.cell)
    assert bp_roundtrip.path == bp.path
    assert len(bp.kpts) == len(bp_roundtrip.kpts)
    assert bp == bp_roundtrip


@given(kpts=kst.kpoints)
def test_roundtrip_Kpoints_dict(kpts):
    kpts_roundtrip = Kpoints.fromdict(kpts.todict())
    assert kpts == kpts_roundtrip
