import pytest
from koopmans import workflows, utils


def test_convergence_h2o(water, mockable, tmp_path, sys2file):
    with utils.chdir(tmp_path):
        parameters = {
            "functional": "dft",
            "task": "convergence",
            "from_scratch": True,
            "convergence_observable": "homo energy",
            "convergence_threshold": "0.1 eV",
            "convergence_parameters": ["ecutwfc", "cell_size"]}

        wf = workflows.ConvergenceWorkflow(parameters=parameters, **water)
        wf.run()
