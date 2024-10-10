import copy
from typing import Any, Dict

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk, molecule
from ase.spacegroup import crystal

np.random.seed(0)


@pytest.fixture
def water() -> Dict[str, Any]:
    # water
    return {'atoms': molecule('H2O', vacuum=5.0, pbc=False),
            'ecutwfc': 20.0,
            'nbnd': 5}


@pytest.fixture
def water_snapshots(water) -> Dict[str, Any]:
    atoms = water['atoms']
    atoms.pbc = True
    water['snapshots'] = []
    for _ in range(5):
        new_atoms = copy.deepcopy(atoms)
        new_atoms.positions += np.random.normal(0, 0.05, (3, 3))
        water['snapshots'].append(new_atoms)

    return water


@pytest.fixture
def silicon() -> Dict[str, Any]:
    # bulk silicon
    from koopmans.kpoints import Kpoints
    from koopmans.projections import ProjectionBlocks

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
    from koopmans.kpoints import Kpoints
    from koopmans.projections import ProjectionBlocks

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
    from koopmans.kpoints import Kpoints
    from koopmans.projections import ProjectionBlocks
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
