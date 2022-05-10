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
        monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', testing.BenchGenKoopmansCPCalculator)


@pytest.fixture(autouse=True)
def mock(monkeypatch, pytestconfig):
    '''
    Instead of running a calculation, access the results from the benchmarks directory
    '''

    if pytestconfig.getoption('mock'):
        monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', testing.MockKoopmansCPCalculator)
