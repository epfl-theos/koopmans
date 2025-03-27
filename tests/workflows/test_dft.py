"""Test the DFT workflows."""

import pytest  # noqa

from koopmans import utils, workflows


def test_dftbands_si(silicon, workflow_patch, tmp_path, sys2file):
    """Calculate the band structure of Si."""
    with utils.chdir(tmp_path):
        silicon['pseudo_library'] = 'PseudoDojo/0.4/PBEsol/SR/standard/upf'
        wf = workflows.DFTBandsWorkflow(
            name='si', **silicon)
        wf.run()


def test_dftph_tio2(tio2, workflow_patch, tmp_path, sys2file):
    """Calculate the dielectric permittivity of TiO2."""
    with utils.chdir(tmp_path):
        tio2['calculator_parameters']['ph'] = {'tr2_ph': 1.0e-14}
        tio2['pseudo_library'] = 'PseudoDojo/0.4/PBEsol/SR/standard/upf'
        wf = workflows.DFTPhWorkflow(
            name='tio2', **tio2)
        wf.run()
