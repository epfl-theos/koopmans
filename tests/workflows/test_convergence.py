import pytest

from koopmans import utils, workflows


def test_convergence_h2o(water, workflow_patch, tmp_path, sys2file):
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
