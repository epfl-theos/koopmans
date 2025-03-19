import pytest

from koopmans import utils, workflows
from koopmans.kpoints import Kpoints


def test_dftbands_si(silicon, workflow_patch, tmp_path, sys2file):
    with utils.chdir(tmp_path):
        silicon['pseudo_library'] = 'PseudoDojo/0.4/PBEsol/SR/standard/upf'
        wf = workflows.DFTBandsWorkflow(
            name='si', **silicon)
        wf.run()


def test_dftph_tio2(tio2, workflow_patch, tmp_path, sys2file):
    with utils.chdir(tmp_path):
        tio2['calculator_parameters']['ph'] = {'tr2_ph': 1.0e-14}
        tio2['pseudo_library'] = 'PseudoDojo/0.4/PBEsol/SR/standard/upf'
        wf = workflows.DFTPhWorkflow(
            name='tio2', **tio2)
        wf.run()
