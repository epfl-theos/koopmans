import pytest
from koopmans import testing


def pytest_addoption(parser):
    parser.addoption("--mock", action="store_true", default=False, help="TODO")
    parser.addoption("--generate_benchmark", action="store_true", default=False, help="TODO")


@pytest.fixture(autouse=True)
def benchmark_gen(monkeypatch, pytestconfig):
    '''
    After each calculation is run, store the results in a json (one json per calculation)

    Only perform this patching if pytest is run with the mark "benchmark_gen"
    '''

    if pytestconfig.getoption('generate_benchmark'):
        monkeypatch.setattr('koopmans.calculators.Wannier90Calculator', testing.BenchGenWannier90Calculator)
        monkeypatch.setattr('koopmans.calculators.PW2WannierCalculator', testing.BenchGenPW2WannierCalculator)
        monkeypatch.setattr('koopmans.calculators.Wann2KCPCalculator', testing.BenchGenWann2KCPCalculator)
        monkeypatch.setattr('koopmans.calculators.PWCalculator', testing.BenchGenPWCalculator)
        monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', testing.BenchGenKoopmansCPCalculator)
        monkeypatch.setattr('koopmans.calculators.EnvironCalculator', testing.BenchGenEnvironCalculator)
        monkeypatch.setattr('koopmans.calculators.UnfoldAndInterpolateCalculator',
                            testing.BenchGenUnfoldAndInterpolateCalculator)
        monkeypatch.setattr('koopmans.calculators.Wann2KCCalculator', testing.BenchGenWann2KCCalculator)
        monkeypatch.setattr('koopmans.calculators.KoopmansScreenCalculator', testing.BenchGenKoopmansScreenCalculator)
        monkeypatch.setattr('koopmans.calculators.KoopmansHamCalculator', testing.BenchGenKoopmansHamCalculator)
        monkeypatch.setattr('koopmans.calculators.ProjwfcCalculator', testing.BenchGenProjwfcCalculator)


@pytest.fixture(autouse=True)
def mock(monkeypatch, pytestconfig):
    '''
    Instead of running a calculation, access the results from the benchmarks directory
    '''

    if pytestconfig.getoption('mock'):
        # Calculators
        monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', testing.MockKoopmansCPCalculator)
        monkeypatch.setattr('koopmans.calculators.Wannier90Calculator', testing.MockWannier90Calculator)
        monkeypatch.setattr('koopmans.calculators.PW2WannierCalculator', testing.MockPW2WannierCalculator)
        monkeypatch.setattr('koopmans.calculators.Wann2KCPCalculator', testing.MockWann2KCPCalculator)
        monkeypatch.setattr('koopmans.calculators.PWCalculator', testing.MockPWCalculator)
        monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', testing.MockKoopmansCPCalculator)
        monkeypatch.setattr('koopmans.calculators.EnvironCalculator', testing.MockEnvironCalculator)
        monkeypatch.setattr('koopmans.calculators.UnfoldAndInterpolateCalculator',
                            testing.MockUnfoldAndInterpolateCalculator)
        monkeypatch.setattr('koopmans.calculators.Wann2KCCalculator', testing.MockWann2KCCalculator)
        monkeypatch.setattr('koopmans.calculators.KoopmansScreenCalculator', testing.MockKoopmansScreenCalculator)
        monkeypatch.setattr('koopmans.calculators.KoopmansHamCalculator', testing.MockKoopmansHamCalculator)
        monkeypatch.setattr('koopmans.calculators.ProjwfcCalculator', testing.MockProjwfcCalculator)

        # Workflows
        monkeypatch.setattr('koopmans.workflows.WannierizeWorkflow', testing.MockWannierizeWorkflow)


@pytest.fixture
def sys2file(capsys, tmp_path):
    # Run the test
    yield

    # Write stdout and stderr to file
    out, err = capsys.readouterr()
    if out or err:
        with open(tmp_path.with_suffix('.stdout'), 'w') as f:
            f.write(out)
            f.write(err)
