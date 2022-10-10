from ._check import (CheckEnvironCalculator, CheckKoopmansCPCalculator,
                     CheckKoopmansHamCalculator, CheckKoopmansScreenCalculator,
                     CheckMLFittingWorkflow, CheckPhCalculator,
                     CheckProjwfcCalculator, CheckPW2WannierCalculator,
                     CheckPWCalculator, CheckUnfoldAndInterpolateCalculator,
                     CheckWann2KCCalculator, CheckWann2KCPCalculator,
                     CheckWannier90Calculator, compare)
from ._generate_benchmarks import (BenchGenEnvironCalculator,
                                   BenchGenKoopmansCPCalculator,
                                   BenchGenKoopmansHamCalculator,
                                   BenchGenKoopmansScreenCalculator,
                                   BenchGenMLFittingWorkflow,
                                   BenchGenPhCalculator,
                                   BenchGenProjwfcCalculator,
                                   BenchGenPW2WannierCalculator,
                                   BenchGenPWCalculator,
                                   BenchGenUnfoldAndInterpolateCalculator,
                                   BenchGenWann2KCCalculator,
                                   BenchGenWann2KCPCalculator,
                                   BenchGenWannier90Calculator)
from ._mock import (MockEnvironCalculator, MockKoopmansCPCalculator,
                    MockKoopmansDSCFWorkflow, MockKoopmansHamCalculator,
                    MockKoopmansScreenCalculator, MockMLFittingWorkflow,
                    MockPhCalculator, MockProjwfcCalculator,
                    MockPW2WannierCalculator, MockPWCalculator,
                    MockUnfoldAndInterpolateCalculator, MockWann2KCCalculator,
                    MockWann2KCPCalculator, MockWannier90Calculator,
                    MockWannierizeWorkflow)
from ._stumble import (StumblingConvergenceMLWorkflow,
                       StumblingConvergenceWorkflow, StumblingDeltaSCFWorkflow,
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
    # After each calculation is run, store the results in a json (one json per calculation)
    monkeypatch.setattr('koopmans.calculators.Wannier90Calculator', BenchGenWannier90Calculator)
    monkeypatch.setattr('koopmans.calculators.PW2WannierCalculator', BenchGenPW2WannierCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCPCalculator', BenchGenWann2KCPCalculator)
    monkeypatch.setattr('koopmans.calculators.PhCalculator', BenchGenPhCalculator)
    monkeypatch.setattr('koopmans.calculators.PWCalculator', BenchGenPWCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', BenchGenKoopmansCPCalculator)
    monkeypatch.setattr('koopmans.calculators.EnvironCalculator', BenchGenEnvironCalculator)
    monkeypatch.setattr('koopmans.calculators.UnfoldAndInterpolateCalculator',
                        BenchGenUnfoldAndInterpolateCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCCalculator', BenchGenWann2KCCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansScreenCalculator', BenchGenKoopmansScreenCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansHamCalculator', BenchGenKoopmansHamCalculator)
    monkeypatch.setattr('koopmans.calculators.ProjwfcCalculator', BenchGenProjwfcCalculator)
    monkeypatch.setattr('koopmans.workflows.MLFittingWorkflow', BenchGenMLFittingWorkflow)


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
    monkeypatch.setattr('koopmans.calculators.UnfoldAndInterpolateCalculator',
                        MockUnfoldAndInterpolateCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCCalculator', MockWann2KCCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansScreenCalculator', MockKoopmansScreenCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansHamCalculator', MockKoopmansHamCalculator)
    monkeypatch.setattr('koopmans.calculators.ProjwfcCalculator', MockProjwfcCalculator)

    # Workflows
    monkeypatch.setattr('koopmans.workflows.KoopmansDSCFWorkflow', MockKoopmansDSCFWorkflow)
    monkeypatch.setattr('koopmans.workflows.WannierizeWorkflow', MockWannierizeWorkflow)
    monkeypatch.setattr('koopmans.workflows.MLFittingWorkflow', MockMLFittingWorkflow)


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
    monkeypatch.setattr('koopmans.calculators.UnfoldAndInterpolateCalculator',
                        CheckUnfoldAndInterpolateCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCCalculator', CheckWann2KCCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansScreenCalculator', CheckKoopmansScreenCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansHamCalculator', CheckKoopmansHamCalculator)
    monkeypatch.setattr('koopmans.calculators.ProjwfcCalculator', CheckProjwfcCalculator)
    monkeypatch.setattr('koopmans.workflows.MLFittingWorkflow', CheckMLFittingWorkflow)


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
    monkeypatch.setattr('koopmans.workflows.ConvergenceMLWorkflow',
                        StumblingConvergenceMLWorkflow)
    # When running with stumble mode, we want to check our results against the benchmarks by using CheckCalcs
    monkeypatch_check(monkeypatch)
