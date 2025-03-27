"""Test the Convergence workflow."""

import pytest  # noqa

from koopmans import utils, workflows


def test_convergence_h2o(water, workflow_patch, tmp_path, sys2file):
    """Test convergence of the HOMO energy of water as a function of ecutwfc and celldm1."""
    with utils.chdir(tmp_path):
        parameters = {
            "functional": "dft",
            "task": "convergence",
            "converge": True}

        subwf = workflows.DFTCPWorkflow(parameters=parameters, **water)
        wf = workflows.ConvergenceWorkflowFactory(subworkflow=subwf,
                                                  observable='homo energy',
                                                  threshold=0.1,
                                                  variables=['ecutwfc', 'celldm1'])
        wf.run()
