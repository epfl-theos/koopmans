import pytest

from koopmans import utils, workflows


def test_convergence_h2o(water, workflow_patch, tmp_path, sys2file):
    with utils.chdir(tmp_path):
        parameters = {
            "functional": "dft",
            "task": "convergence",
            "keep_tmpdirs": False,
            "convergence_observable": "homo energy",
            "convergence_threshold": "0.1 eV",
            "convergence_parameters": ["ecutwfc", "cell_size"]}

        wf = workflows.ConvergenceWorkflow(parameters=parameters, **water)
        wf.run()
