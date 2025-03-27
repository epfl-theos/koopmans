import pytest  # noqa

from ._benchmark import monkeypatch_bench
from ._check import monkeypatch_check
from ._mock import monkeypatch_mock
from ._utils import benchmark_filename  # noqa: F401


def monkeypatch_stumble(monkeypatch):
    """Patch the workflows to stumble."""
    from ._stumble import (StumblingConvergenceWorkflow,
                           StumblingDeltaSCFWorkflow, StumblingDFTCPWorkflow,
                           StumblingDFTPhWorkflow, StumblingDFTPWWorkflow,
                           StumblingFoldToSupercellWorkflow,
                           StumblingKoopmansDFPTWorkflow,
                           StumblingKoopmansDSCFWorkflow,
                           StumblingSinglepointWorkflow,
                           StumblingTrajectoryWorkflow,
                           StumblingUnfoldAndInterpolateWorkflow,
                           StumblingWannierizeWorkflow)
    monkeypatch.setattr('koopmans.workflows.WannierizeWorkflow', StumblingWannierizeWorkflow)
    monkeypatch.setattr('koopmans.workflows.KoopmansDSCFWorkflow', StumblingKoopmansDSCFWorkflow)
    monkeypatch.setattr('koopmans.workflows.SinglepointWorkflow', StumblingSinglepointWorkflow)
    monkeypatch.setattr('koopmans.workflows.ConvergenceWorkflow', StumblingConvergenceWorkflow)
    monkeypatch.setattr('koopmans.workflows.FoldToSupercellWorkflow', StumblingFoldToSupercellWorkflow)
    monkeypatch.setattr('koopmans.workflows.DFTCPWorkflow', StumblingDFTCPWorkflow)
    monkeypatch.setattr('koopmans.workflows.DFTPhWorkflow', StumblingDFTPhWorkflow)
    monkeypatch.setattr('koopmans.workflows.DFTPWWorkflow', StumblingDFTPWWorkflow)
    monkeypatch.setattr('koopmans.workflows.DeltaSCFWorkflow', StumblingDeltaSCFWorkflow)
    monkeypatch.setattr('koopmans.workflows.KoopmansDFPTWorkflow', StumblingKoopmansDFPTWorkflow)
    monkeypatch.setattr('koopmans.workflows.UnfoldAndInterpolateWorkflow',
                        StumblingUnfoldAndInterpolateWorkflow)
    monkeypatch.setattr('koopmans.workflows.TrajectoryWorkflow',
                        StumblingTrajectoryWorkflow)
    # When running with stumble mode, we want to check our results against the benchmarks by using CheckCalcs
    monkeypatch_check(monkeypatch)


@pytest.fixture
def tutorial_patch(monkeypatch, pytestconfig):
    """Patch the tutorials according to the user's request."""
    # For the tutorials...
    if pytestconfig.getoption('generate_benchmark'):
        # when generating benchmarks, use BenchCalcs
        monkeypatch_bench(monkeypatch)
    elif pytestconfig.getoption('stumble'):
        # when testing recovery from a crash, use StumblingWorkflows
        monkeypatch_stumble(monkeypatch)
    else:
        # we use MockCalcs when running our tests on github, OR if the user is running locally
        monkeypatch_mock(monkeypatch)


@pytest.fixture
def workflow_patch(monkeypatch, pytestconfig):
    """Patch the workflows according to the user's request."""
    # For tests involving the workflow...
    if pytestconfig.getoption('generate_benchmark'):
        # when generating benchmarks, use BenchCalcs
        monkeypatch_bench(monkeypatch)
    elif pytestconfig.getoption('stumble'):
        # when testing recovery from a crash, use StumblingWorkflows
        monkeypatch_stumble(monkeypatch)
    else:
        # we use MockCalcs when running our tests on github, OR if the user is running locally
        monkeypatch_mock(monkeypatch)


@pytest.fixture
def check_patch(monkeypatch, pytestconfig):
    """Patch any calculations that involve python routines only according to the users's request."""
    # For calculations that involve python routines only (such as Processes)...
    if pytestconfig.getoption('generate_benchmark'):
        # when generating benchmarks, use BenchCalcs
        monkeypatch_bench(monkeypatch)
    elif pytestconfig.getoption('stumble'):
        # when testing recovery from a crash, use StumblingWorkflows
        monkeypatch_stumble(monkeypatch)
    else:
        # we rely on CheckProcesses to compare the results
        monkeypatch_check(monkeypatch)


@pytest.fixture
def espresso_patch(monkeypatch, pytestconfig):
    """Patch the espresso calculations according to the user's request."""
    # For tests involving Quantum ESPRESSO...
    if pytestconfig.getoption('generate_benchmark'):
        # when generating benchmarks, use BenchCalcs
        monkeypatch_bench(monkeypatch)
    elif pytestconfig.getoption('stumble'):
        # when testing recovery from a crash, use StumblingWorkflows
        monkeypatch_stumble(monkeypatch)
    elif pytestconfig.getoption('ci'):
        # when running our tests on github, these tests should not be called!
        raise ValueError('These tests cannot be run with --ci')
    else:
        # when the user is running locally, use CheckCalcs
        monkeypatch_check(monkeypatch)
