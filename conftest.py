from typing import Any, Dict

import pytest

from ase import Atoms
from ase.build import bulk, molecule
from ase.spacegroup import crystal
from koopmans import base_directory, testing
from koopmans.kpoints import Kpoints
from koopmans.projections import ProjectionBlocks


def pytest_addoption(parser):
    parser.addoption("--ci", action="store_true", default=False,
                     help="Run only those tests that do not require an installation of Quantum ESPRESSO")
    parser.addoption("--generate_benchmark", action="store_true", default=False, help="Generate new benchmark files")
    parser.addoption("--stumble", action="store_true", default=False,
                     help="Make the workflows deliberately crash and restart after each calculation, in order to test "
                     "that the code can restart at any stage")


@pytest.fixture
def datadir():
    # Returns the directory where various reference QE files are stored
    return base_directory / 'tests' / 'data'


def monkeypatch_bench(monkeypatch):
    # After each calculation is run, store the results in a json (one json per calculation)
    monkeypatch.setattr('koopmans.calculators.Wannier90Calculator', testing.BenchGenWannier90Calculator)
    monkeypatch.setattr('koopmans.calculators.PW2WannierCalculator', testing.BenchGenPW2WannierCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCPCalculator', testing.BenchGenWann2KCPCalculator)
    monkeypatch.setattr('koopmans.calculators.PhCalculator', testing.BenchGenPhCalculator)
    monkeypatch.setattr('koopmans.calculators.PWCalculator', testing.BenchGenPWCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', testing.BenchGenKoopmansCPCalculator)
    monkeypatch.setattr('koopmans.calculators.EnvironCalculator', testing.BenchGenEnvironCalculator)
    monkeypatch.setattr('koopmans.calculators.UnfoldAndInterpolateCalculator',
                        testing.BenchGenUnfoldAndInterpolateCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCCalculator', testing.BenchGenWann2KCCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansScreenCalculator', testing.BenchGenKoopmansScreenCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansHamCalculator', testing.BenchGenKoopmansHamCalculator)
    monkeypatch.setattr('koopmans.calculators.ProjwfcCalculator', testing.BenchGenProjwfcCalculator)


def monkeypatch_mock(monkeypatch):
    # Replace calculators with mock versions that obtain results from the database
    monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', testing.MockKoopmansCPCalculator)
    monkeypatch.setattr('koopmans.calculators.Wannier90Calculator', testing.MockWannier90Calculator)
    monkeypatch.setattr('koopmans.calculators.PW2WannierCalculator', testing.MockPW2WannierCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCPCalculator', testing.MockWann2KCPCalculator)
    monkeypatch.setattr('koopmans.calculators.PhCalculator', testing.MockPhCalculator)
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
    monkeypatch.setattr('koopmans.workflows.KoopmansDSCFWorkflow', testing.MockKoopmansDSCFWorkflow)
    monkeypatch.setattr('koopmans.workflows.WannierizeWorkflow', testing.MockWannierizeWorkflow)


def monkeypatch_check(monkeypatch):
    # Replace calculators with versions that double-check their results against
    monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', testing.CheckKoopmansCPCalculator)
    monkeypatch.setattr('koopmans.calculators.Wannier90Calculator', testing.CheckWannier90Calculator)
    monkeypatch.setattr('koopmans.calculators.PW2WannierCalculator', testing.CheckPW2WannierCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCPCalculator', testing.CheckWann2KCPCalculator)
    monkeypatch.setattr('koopmans.calculators.PhCalculator', testing.CheckPhCalculator)
    monkeypatch.setattr('koopmans.calculators.PWCalculator', testing.CheckPWCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', testing.CheckKoopmansCPCalculator)
    monkeypatch.setattr('koopmans.calculators.EnvironCalculator', testing.CheckEnvironCalculator)
    monkeypatch.setattr('koopmans.calculators.UnfoldAndInterpolateCalculator',
                        testing.CheckUnfoldAndInterpolateCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCCalculator', testing.CheckWann2KCCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansScreenCalculator', testing.CheckKoopmansScreenCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansHamCalculator', testing.CheckKoopmansHamCalculator)
    monkeypatch.setattr('koopmans.calculators.ProjwfcCalculator', testing.CheckProjwfcCalculator)


def monkeypatch_stumble(monkeypatch):
    monkeypatch.setattr('koopmans.workflows.WannierizeWorkflow', testing.StumblingWannierizeWorkflow)
    monkeypatch.setattr('koopmans.workflows.KoopmansDSCFWorkflow', testing.StumblingKoopmansDSCFWorkflow)
    monkeypatch.setattr('koopmans.workflows.SinglepointWorkflow', testing.StumblingSinglepointWorkflow)
    monkeypatch.setattr('koopmans.workflows.ConvergenceWorkflow', testing.StumblingConvergenceWorkflow)
    monkeypatch.setattr('koopmans.workflows.FoldToSupercellWorkflow', testing.StumblingFoldToSupercellWorkflow)
    monkeypatch.setattr('koopmans.workflows.DFTCPWorkflow', testing.StumblingDFTCPWorkflow)
    monkeypatch.setattr('koopmans.workflows.DFTPhWorkflow', testing.StumblingDFTPhWorkflow)
    monkeypatch.setattr('koopmans.workflows.DFTPWWorkflow', testing.StumblingDFTPWWorkflow)
    monkeypatch.setattr('koopmans.workflows.DeltaSCFWorkflow', testing.StumblingDeltaSCFWorkflow)
    monkeypatch.setattr('koopmans.workflows.KoopmansDFPTWorkflow', testing.StumblingKoopmansDFPTWorkflow)
    monkeypatch.setattr('koopmans.workflows.UnfoldAndInterpolateWorkflow',
                        testing.StumblingUnfoldAndInterpolateWorkflow)
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
def ui_patch(monkeypatch, pytestconfig):
    # For tests involving the UI python routines only...
    if pytestconfig.getoption('generate_benchmark'):
        # when generating benchmarks, use BenchCalcs
        monkeypatch_bench(monkeypatch)
    elif pytestconfig.getoption('stumble'):
        # when testing recovery from a crash, use StumblingWorkflows
        monkeypatch_stumble(monkeypatch)
    else:
        # we can run the calculations directly when running our tests on github, OR if the user is running locally
        pass


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
        # when running our tests on github, these tests shold not be called!
        raise ValueError('These tests cannot be run with --ci')
    else:
        # when the user is running locally, use CheckCalcs
        monkeypatch_check(monkeypatch)


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
def water() -> Dict[str, Any]:
    # water
    return {'atoms': molecule('H2O', vacuum=5.0, pbc=False),
            'ecutwfc': 20.0,
            'nbnd': 5}  # calculator_parameters': {'kcp': {'ecutwfc': 20.0, 'nbnd': 5}}}


@pytest.fixture
def silicon() -> Dict[str, Any]:
    # bulk silicon
    si: Atoms = bulk('Si')
    pdict = [{'fsite': [0.25, 0.25, 0.25], 'ang_mtm': 'sp3'}]
    si_projs = ProjectionBlocks.fromprojections([pdict, pdict], fillings=[True, False], spins=[None, None], atoms=si)
    kpoints = Kpoints(grid=[2, 2, 2], path='GXG', cell=si.cell)
    return {'atoms': si,
            'calculator_parameters': {'pw': {'nbnd': 10},
                                      'w90_emp': {'dis_froz_max': 10.6, 'dis_win_max': 16.9}
                                      },
            'plotting': {'Emin': -10, 'Emax': 4, 'degauss': 0.5},
            'projections': si_projs,
            'kpoints': kpoints,
            'ecutwfc': 40.0,
            'smooth_int_factor': 2}


@pytest.fixture
def ozone() -> Dict[str, Any]:
    # ozone
    return {'atoms': molecule('O3', vacuum=5.0, pbc=False),
            'calculator_parameters': {'pw': {'ecutwfc': 20.0, 'nbnd': 10}}}


@pytest.fixture
def tio2() -> Dict[str, Any]:
    # rutile TiO2
    a = 4.6
    c = 2.95
    atoms: Atoms = crystal(['Ti', 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                           spacegroup=136, cellpar=[a, a, c, 90, 90, 90])

    projs = ProjectionBlocks.fromprojections([["Ti:l=0"], ["Ti:l=1"], ["O:l=0"], ["O:l=1"], ["Ti:l=2"]],
                                             fillings=[True, True, True, True, False],
                                             spins=[None, None, None, None, None],
                                             atoms=atoms)

    kpoints = Kpoints(grid=[2, 2, 2], path='GXG')
    return {'atoms': atoms,
            'calculator_parameters': {'pw': {'nbnd': 34}},
            'projections': projs,
            'kpoints': kpoints,
            'ecutwfc': 40.0}


@pytest.fixture
def gaas() -> Dict[str, Any]:
    # bulk gallium arsenide
    atoms: Atoms = bulk('GaAs', crystalstructure='zincblende', a=5.6536)
    gaas_projs = ProjectionBlocks.fromprojections([["Ga: d"], ["As: sp3"], ["Ga: sp3"]],
                                                  fillings=[True, True, False],
                                                  spins=[None, None, None],
                                                  atoms=atoms)
    kpoints = Kpoints(grid=[2, 2, 2])
    return {'atoms': atoms,
            'calculator_parameters': {'pw': {'nbnd': 45},
                                      'w90_emp': {'dis_froz_max': 14.6, 'dis_win_max': 18.6}
                                      },
            'ecutwfc': 40.0,
            'smooth_int_factor': 4,
            'plotting': {'degauss': 0.5},
            'projections': gaas_projs,
            'kpoints': kpoints}
