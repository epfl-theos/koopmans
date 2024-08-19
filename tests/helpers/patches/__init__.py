import pytest

from ._benchmark import (BenchGenBin2XMLProcess,
                         BenchGenComputePowerSpectrumProcess,
                         BenchGenConvertFilesFromSpin1To2,
                         BenchGenConvertFilesFromSpin2To1,
                         BenchGenEnvironCalculator, BenchGenExtendProcess,
                         BenchGenExtractCoefficientsFromXMLProcess,
                         BenchGenKoopmansCPCalculator,
                         BenchGenKoopmansHamCalculator,
                         BenchGenKoopmansScreenCalculator,
                         BenchGenMergeProcess, BenchGenPhCalculator,
                         BenchGenProjwfcCalculator,
                         BenchGenPW2WannierCalculator, BenchGenPWCalculator,
                         BenchGenUnfoldAndInterpolateProcess,
                         BenchGenWann2KCCalculator, BenchGenWann2KCPCalculator,
                         BenchGenWannier90Calculator)
from ._check import (CheckBin2XMLProcess, CheckComputePowerSpectrumProcess,
                     CheckConvertFilesFromSpin1To2,
                     CheckConvertFilesFromSpin2To1, CheckEnvironCalculator,
                     CheckExtendProcess,
                     CheckExtractCoefficientsFromXMLProcess,
                     CheckKoopmansCPCalculator, CheckKoopmansHamCalculator,
                     CheckKoopmansScreenCalculator, CheckMergeProcess,
                     CheckPhCalculator, CheckProjwfcCalculator,
                     CheckPW2WannierCalculator, CheckPWCalculator,
                     CheckUnfoldAndInterpolateProcess, CheckWann2KCCalculator,
                     CheckWann2KCPCalculator, CheckWannier90Calculator)
from ._mock import (MockBin2XMLProcess, MockComputePowerSpectrumProcess,
                    MockConvertFilesFromSpin1To2, MockConvertFilesFromSpin2To1,
                    MockEnvironCalculator, MockExtendProcess,
                    MockExtractCoefficientsFromXMLProcess,
                    MockKoopmansCPCalculator, MockKoopmansDSCFWorkflow,
                    MockKoopmansHamCalculator, MockKoopmansScreenCalculator,
                    MockMergeProcess, MockPhCalculator, MockProjwfcCalculator,
                    MockPW2WannierCalculator, MockPWCalculator,
                    MockUnfoldAndInterpolateProcess, MockWann2KCCalculator,
                    MockWann2KCPCalculator, MockWannier90Calculator,
                    MockWannierizeWorkflow)
from ._stumble import (StumblingConvergenceWorkflow, StumblingDeltaSCFWorkflow,
                       StumblingDFTCPWorkflow, StumblingDFTPhWorkflow,
                       StumblingDFTPWWorkflow,
                       StumblingFoldToSupercellWorkflow,
                       StumblingKoopmansDFPTWorkflow,
                       StumblingKoopmansDSCFWorkflow,
                       StumblingSinglepointWorkflow,
                       StumblingTrajectoryWorkflow,
                       StumblingUnfoldAndInterpolateWorkflow,
                       StumblingWannierizeWorkflow)
from ._utils import benchmark_filename


def monkeypatch_bench(monkeypatch):
    # After each calculation is run, store the results
    monkeypatch.setattr('koopmans.calculators.Wannier90Calculator', BenchGenWannier90Calculator)
    monkeypatch.setattr('koopmans.calculators.PW2WannierCalculator', BenchGenPW2WannierCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCPCalculator', BenchGenWann2KCPCalculator)
    monkeypatch.setattr('koopmans.calculators.PhCalculator', BenchGenPhCalculator)
    monkeypatch.setattr('koopmans.calculators.PWCalculator', BenchGenPWCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', BenchGenKoopmansCPCalculator)
    monkeypatch.setattr('koopmans.calculators.EnvironCalculator', BenchGenEnvironCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCCalculator', BenchGenWann2KCCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansScreenCalculator', BenchGenKoopmansScreenCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansHamCalculator', BenchGenKoopmansHamCalculator)
    monkeypatch.setattr('koopmans.calculators.ProjwfcCalculator', BenchGenProjwfcCalculator)

    # Processes
    monkeypatch.setattr('koopmans.processes.power_spectrum.ExtractCoefficientsFromXMLProcess',
                        BenchGenExtractCoefficientsFromXMLProcess)
    monkeypatch.setattr('koopmans.processes.power_spectrum.ComputePowerSpectrumProcess',
                        BenchGenComputePowerSpectrumProcess)
    monkeypatch.setattr('koopmans.processes.bin2xml.Bin2XMLProcess', BenchGenBin2XMLProcess)
    monkeypatch.setattr('koopmans.processes.koopmans_cp.ConvertFilesFromSpin1To2', BenchGenConvertFilesFromSpin1To2)
    monkeypatch.setattr('koopmans.processes.koopmans_cp.ConvertFilesFromSpin2To1', BenchGenConvertFilesFromSpin2To1)
    monkeypatch.setattr('koopmans.processes.wannier.ExtendProcess', BenchGenExtendProcess)
    monkeypatch.setattr('koopmans.processes.wannier.MergeProcess', BenchGenMergeProcess)
    monkeypatch.setattr('koopmans.processes.ui.UnfoldAndInterpolateProcess', BenchGenUnfoldAndInterpolateProcess)


def monkeypatch_mock(monkeypatch):
    # Replace calculators with mock versions that obtain results from the database
    monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', MockKoopmansCPCalculator)
    monkeypatch.setattr('koopmans.calculators.Wannier90Calculator', MockWannier90Calculator)
    monkeypatch.setattr('koopmans.calculators.PW2WannierCalculator', MockPW2WannierCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCPCalculator', MockWann2KCPCalculator)
    monkeypatch.setattr('koopmans.calculators.PhCalculator', MockPhCalculator)
    monkeypatch.setattr('koopmans.calculators.PWCalculator', MockPWCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', MockKoopmansCPCalculator)
    monkeypatch.setattr('koopmans.calculators.EnvironCalculator', MockEnvironCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCCalculator', MockWann2KCCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansScreenCalculator', MockKoopmansScreenCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansHamCalculator', MockKoopmansHamCalculator)
    monkeypatch.setattr('koopmans.calculators.ProjwfcCalculator', MockProjwfcCalculator)

    # Workflows
    monkeypatch.setattr('koopmans.workflows.KoopmansDSCFWorkflow', MockKoopmansDSCFWorkflow)
    monkeypatch.setattr('koopmans.workflows.WannierizeWorkflow', MockWannierizeWorkflow)

    # Processes
    monkeypatch.setattr('koopmans.processes.power_spectrum.ExtractCoefficientsFromXMLProcess',
                        MockExtractCoefficientsFromXMLProcess)
    monkeypatch.setattr('koopmans.processes.power_spectrum.ComputePowerSpectrumProcess',
                        MockComputePowerSpectrumProcess)
    monkeypatch.setattr('koopmans.processes.bin2xml.Bin2XMLProcess', MockBin2XMLProcess)
    monkeypatch.setattr('koopmans.processes.koopmans_cp.ConvertFilesFromSpin1To2', MockConvertFilesFromSpin1To2)
    monkeypatch.setattr('koopmans.processes.koopmans_cp.ConvertFilesFromSpin2To1', MockConvertFilesFromSpin2To1)
    monkeypatch.setattr('koopmans.processes.wannier.ExtendProcess', MockExtendProcess)
    monkeypatch.setattr('koopmans.processes.wannier.MergeProcess', MockMergeProcess)
    monkeypatch.setattr('koopmans.processes.ui.UnfoldAndInterpolateProcess', MockUnfoldAndInterpolateProcess)


def monkeypatch_check(monkeypatch):
    # Replace calculators with versions that double-check their results against
    monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', CheckKoopmansCPCalculator)
    monkeypatch.setattr('koopmans.calculators.Wannier90Calculator', CheckWannier90Calculator)
    monkeypatch.setattr('koopmans.calculators.PW2WannierCalculator', CheckPW2WannierCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCPCalculator', CheckWann2KCPCalculator)
    monkeypatch.setattr('koopmans.calculators.PhCalculator', CheckPhCalculator)
    monkeypatch.setattr('koopmans.calculators.PWCalculator', CheckPWCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', CheckKoopmansCPCalculator)
    monkeypatch.setattr('koopmans.calculators.EnvironCalculator', CheckEnvironCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCCalculator', CheckWann2KCCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansScreenCalculator', CheckKoopmansScreenCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansHamCalculator', CheckKoopmansHamCalculator)
    monkeypatch.setattr('koopmans.calculators.ProjwfcCalculator', CheckProjwfcCalculator)

    # Processes
    monkeypatch.setattr('koopmans.processes.power_spectrum.ExtractCoefficientsFromXMLProcess',
                        CheckExtractCoefficientsFromXMLProcess)
    monkeypatch.setattr('koopmans.processes.power_spectrum.ComputePowerSpectrumProcess',
                        CheckComputePowerSpectrumProcess)
    monkeypatch.setattr('koopmans.processes.bin2xml.Bin2XMLProcess', CheckBin2XMLProcess)
    monkeypatch.setattr('koopmans.processes.koopmans_cp.ConvertFilesFromSpin1To2', CheckConvertFilesFromSpin1To2)
    monkeypatch.setattr('koopmans.processes.koopmans_cp.ConvertFilesFromSpin2To1', CheckConvertFilesFromSpin2To1)
    monkeypatch.setattr('koopmans.processes.wannier.ExtendProcess', CheckExtendProcess)
    monkeypatch.setattr('koopmans.processes.wannier.MergeProcess', CheckMergeProcess)
    monkeypatch.setattr('koopmans.processes.ui.UnfoldAndInterpolateProcess', CheckUnfoldAndInterpolateProcess)


def monkeypatch_stumble(monkeypatch):
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
