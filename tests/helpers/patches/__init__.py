import pytest

from ._benchmark import monkeypatch_bench
from ._check import monkeypatch_check
from ._mock import monkeypatch_mock
from ._utils import benchmark_filename


@pytest.fixture(autouse=True)
def patch_path_comparison(monkeypatch):
    def path_eq(self, other):
        if not isinstance(other, type(self)):
            return False
        if self.is_absolute() == other.is_absolute():
            return str(self) == str(other)
        else:
            # In some cases we store paths relative to a base directory. Without knowing
            # the base directory, the best we can do is check that the relative path
            # ends the same as the absolute path
            return str(self.resolve()).endswith(str(other)) or str(other).endswith(str(self))
    monkeypatch.setattr('pathlib.Path.__eq__', path_eq)


def monkeypatch_stumble(monkeypatch):
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
