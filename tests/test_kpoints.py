"""Testing koopmans.kpooints."""

import numpy as np
import pytest  # noqa: F401
from ase_koopmans.dft.kpoints import BandPath
from hypothesis import given, settings

from koopmans.kpoints import Kpoints, dict_to_kpath, kpath_to_dict
from tests.helpers import strategies as kst


@given(bp=kst.bandpaths())
@settings(report_multiple_bugs=False, deadline=None)
def test_roundtrip_kpath_dict(bp: BandPath):
    """Test that BandPath can be converted to a dictionary and back."""
    # TODO work out why bp.special_points changes
    bp_roundtrip = dict_to_kpath(kpath_to_dict(bp))
    assert np.allclose(bp_roundtrip.cell, bp.cell)
    assert bp_roundtrip.path == bp.path
    assert len(bp.kpts) == len(bp_roundtrip.kpts)


@given(kpts=kst.kpoints)
@settings(report_multiple_bugs=False, deadline=None)
def test_roundtrip_Kpoints_dict(kpts):
    """Test that Kpoints can be converted to a dictionary and back."""
    # TODO work out why kpts.path.special_points changes
    kpts_roundtrip = Kpoints.fromdict(kpts.todict())
    assert kpts.gamma_only == kpts_roundtrip.gamma_only
    assert kpts.grid == kpts_roundtrip.grid
    assert np.allclose(kpts.path.cell, kpts_roundtrip.path.cell)
