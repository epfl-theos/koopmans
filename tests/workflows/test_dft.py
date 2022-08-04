import pytest

from koopmans import utils, workflows
from koopmans.kpoints import Kpoints


def test_pwbandstructure_si(silicon, workflow_patch, tmp_path, sys2file):
    with utils.chdir(tmp_path):
        wf = workflows.PWBandStructureWorkflow(
            parameters={'pseudo_library': 'pseudo_dojo_standard', 'base_functional': 'pbesol', 'keep_tmpdirs': False},
            name='si', **silicon)
        wf.run()


def test_dftph_tio2(tio2, workflow_patch, tmp_path, sys2file):
    with utils.chdir(tmp_path):
        tio2['calculator_parameters']['ph'] = {'tr2_ph': 1.0e-14}
        wf = workflows.DFTPhWorkflow(
            parameters={'pseudo_library': 'pseudo_dojo_standard', 'base_functional': 'pbesol', 'keep_tmpdirs': False},
            name='tio2', **tio2)
        wf.run()
