from pathlib import Path
from typing import Any, Dict

import pytest
from ase import Atoms
from ase.build import bulk, molecule
from ase.spacegroup import crystal

from koopmans.kpoints import Kpoints
from koopmans.projections import ProjectionBlocks
from tests import patches


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
    return Path(__file__).parent / 'tests' / 'data'


@pytest.fixture
def tutorial_patch(monkeypatch, pytestconfig):
    # For the tutorials...
    if pytestconfig.getoption('generate_benchmark'):
        # when generating benchmarks, use BenchCalcs
        patches.monkeypatch_bench(monkeypatch)
    elif pytestconfig.getoption('stumble'):
        # when testing recovery from a crash, use StumblingWorkflows
        patches.monkeypatch_stumble(monkeypatch)
    else:
        # we use MockCalcs when running our tests on github, OR if the user is running locally
        patches.monkeypatch_mock(monkeypatch)


@pytest.fixture
def workflow_patch(monkeypatch, pytestconfig):
    # For tests involving the workflow...
    if pytestconfig.getoption('generate_benchmark'):
        # when generating benchmarks, use BenchCalcs
        patches.monkeypatch_bench(monkeypatch)
    elif pytestconfig.getoption('stumble'):
        # when testing recovery from a crash, use StumblingWorkflows
        patches.monkeypatch_stumble(monkeypatch)
    else:
        # we use MockCalcs when running our tests on github, OR if the user is running locally
        patches.monkeypatch_mock(monkeypatch)


@pytest.fixture
def ui_patch(monkeypatch, pytestconfig):
    # For tests involving the UI python routines only...
    if pytestconfig.getoption('generate_benchmark'):
        # when generating benchmarks, use BenchCalcs
        patches.monkeypatch_bench(monkeypatch)
    elif pytestconfig.getoption('stumble'):
        # when testing recovery from a crash, use StumblingWorkflows
        patches.monkeypatch_stumble(monkeypatch)
    else:
        # we can run the calculations directly when running our tests on github, OR if the user is running locally
        pass


@pytest.fixture
def espresso_patch(monkeypatch, pytestconfig):
    # For tests involving Quantum ESPRESSO...
    if pytestconfig.getoption('generate_benchmark'):
        # when generating benchmarks, use BenchCalcs
        patches.monkeypatch_bench(monkeypatch)
    elif pytestconfig.getoption('stumble'):
        # when testing recovery from a crash, use StumblingWorkflows
        patches.monkeypatch_stumble(monkeypatch)
    elif pytestconfig.getoption('ci'):
        # when running our tests on github, these tests shold not be called!
        raise ValueError('These tests cannot be run with --ci')
    else:
        # when the user is running locally, use CheckCalcs
        patches.monkeypatch_check(monkeypatch)


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
            'nbnd': 5}


@pytest.fixture
def silicon() -> Dict[str, Any]:
    # bulk silicon
    si: Atoms = bulk('Si')
    pdict = [{'fsite': [0.25, 0.25, 0.25], 'ang_mtm': 'sp3'}]
    si_projs = ProjectionBlocks.fromlist([pdict, pdict], spins=[None, None], atoms=si)
    kpoints = Kpoints(grid=[2, 2, 2], path='GXG', cell=si.cell)
    return {'atoms': si,
            'calculator_parameters': {'pw': {'nbnd': 10},
                                      'w90': {'dis_froz_max': 10.6, 'dis_win_max': 16.9}
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

    projs = ProjectionBlocks.fromlist([["Ti:l=0"], ["Ti:l=1"], ["O:l=0"], ["O:l=1"], ["Ti:l=2"]],
                                      spins=[None, None, None, None, None],
                                      atoms=atoms)

    kpoints = Kpoints(grid=[2, 2, 2], path='GXG', cell=atoms.cell)
    return {'atoms': atoms,
            'calculator_parameters': {'pw': {'nbnd': 34}},
            'projections': projs,
            'kpoints': kpoints,
            'ecutwfc': 40.0}


@pytest.fixture
def gaas() -> Dict[str, Any]:
    # bulk gallium arsenide
    atoms: Atoms = bulk('GaAs', crystalstructure='zincblende', a=5.6536)
    gaas_projs = ProjectionBlocks.fromlist([["Ga: d"], ["As: sp3"], ["Ga: sp3"]],
                                           spins=[None, None, None],
                                           atoms=atoms)
    kpoints = Kpoints(grid=[2, 2, 2])
    return {'atoms': atoms,
            'calculator_parameters': {'pw': {'nbnd': 45},
                                      'w90': {'dis_froz_max': 14.6, 'dis_win_max': 18.6}
                                      },
            'ecutwfc': 40.0,
            'smooth_int_factor': 4,
            'plotting': {'degauss': 0.5},
            'projections': gaas_projs,
            'kpoints': kpoints}
