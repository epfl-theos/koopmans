import pytest
from koopmans import workflows, utils


def test_pwbandstructure_si(silicon, workflow_patch, tmp_path, sys2file):
    with utils.chdir(tmp_path):
        wf = workflows.PWBandStructureWorkflow(
            parameters={'pseudo_library': 'pseudo_dojo_standard', 'base_functional': 'pbesol', 'keep_tmpdirs': False},
            name='si', **silicon)
        wf.run()
