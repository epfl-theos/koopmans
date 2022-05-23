import pytest
from ase.build import molecule, bulk
from ase.spacegroup import crystal
from koopmans import testing, base_directory
from koopmans.projections import ProjectionBlocks


def pytest_addoption(parser):
    parser.addoption("--mock", action="store_true", default=False, help="TODO")
    parser.addoption("--generate_benchmark", action="store_true", default=False, help="TODO")


@pytest.fixture
def datadir():
    # Returns the directory where various reference QE files are stored
    return base_directory / 'tests' / 'data'


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


@pytest.fixture
def mockable(monkeypatch, pytestconfig):
    '''
    Marker identifying tests where we either...
     a) instead of running a calculation, access the results from the benchmarks directory
    or
     b) after running a calculation. check the results against the benchmark
    '''

    if pytestconfig.getoption('mock'):
        # Replace calculators with mock versions that obtain results from the database
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

    elif not pytestconfig.getoption('generate_benchmark'):
        # Replace calculators with versions that double-check their results against
        monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', testing.CheckKoopmansCPCalculator)
        monkeypatch.setattr('koopmans.calculators.Wannier90Calculator', testing.CheckWannier90Calculator)
        monkeypatch.setattr('koopmans.calculators.PW2WannierCalculator', testing.CheckPW2WannierCalculator)
        monkeypatch.setattr('koopmans.calculators.Wann2KCPCalculator', testing.CheckWann2KCPCalculator)
        monkeypatch.setattr('koopmans.calculators.PWCalculator', testing.CheckPWCalculator)
        monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', testing.CheckKoopmansCPCalculator)
        monkeypatch.setattr('koopmans.calculators.EnvironCalculator', testing.CheckEnvironCalculator)
        monkeypatch.setattr('koopmans.calculators.UnfoldAndInterpolateCalculator',
                            testing.CheckUnfoldAndInterpolateCalculator)
        monkeypatch.setattr('koopmans.calculators.Wann2KCCalculator', testing.CheckWann2KCCalculator)
        monkeypatch.setattr('koopmans.calculators.KoopmansScreenCalculator', testing.CheckKoopmansScreenCalculator)
        monkeypatch.setattr('koopmans.calculators.KoopmansHamCalculator', testing.CheckKoopmansHamCalculator)
        monkeypatch.setattr('koopmans.calculators.ProjwfcCalculator', testing.CheckProjwfcCalculator)


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


@pytest.fixture
def water():
    # water
    return {'atoms': molecule('H2O', vacuum=5.0, pbc=False),
            'master_calc_params': {'kcp': {'ecutwfc': 20.0, 'nbnd': 5}}}


@pytest.fixture
def silicon():
    # bulk silicon
    si = bulk('Si')
    pdict = [{'fsite': [0.25, 0.25, 0.25], 'ang_mtm': 'sp3'}]
    si_projs = ProjectionBlocks.fromprojections([pdict, pdict], fillings=[True, False], spins=[None, None], atoms=si)
    return {'atoms': si,
            'master_calc_params': {'kcp': {'ecutwfc': 40.0},
                                   'pw': {'ecutwfc': 40.0, 'nbnd': 10},
                                   'w90_occ': {'conv_window': 5, },
                                   'w90_emp': {'conv_window': 5, 'dis_froz_max': 10.6, 'dis_win_max': 16.9},
                                   'ui': {'smooth_int_factor': 2},
                                   'plot': {'Emin': -10, 'Emax': 4, 'degauss': 0.5}
                                   },
            'projections': si_projs}


@pytest.fixture
def ozone():
    # ozone
    return {'atoms': molecule('O3', vacuum=5.0, pbc=False),
            'master_calc_params': {'pw': {'ecutwfc': 20.0, 'nbnd': 10}}}


@pytest.fixture
def tio2():
    # rutile TiO2
    a = 4.6
    c = 2.95
    atoms = crystal(['Ti', 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                    spacegroup=136, cellpar=[a, a, c, 90, 90, 90])

    projs = ProjectionBlocks.fromprojections([["Ti:l=0"], ["Ti:l=1"], ["O:l=0"], ["O:l=1"], ["Ti:l=0", "Ti:l=2"]],
                                             fillings=[True, True, True, True, False],
                                             spins=[None, None, None, None, None],
                                             atoms=atoms)

    return {'atoms': atoms,
            'projections': projs,
            'master_calc_params': {'pw': {'ecutwfc': 40, 'nbnd': 36}}}
