import pytest
from ase.build import molecule, bulk
from koopmans.workflows import SinglepointWorkflow
from koopmans.projections import ProjectionBlocks
from koopmans.io import write
from koopmans.utils import chdir


def test_singlepoint_h2o_all_dscf(tmp_path, sys2file):
    with chdir(tmp_path):
        atoms = molecule('H2O', vacuum=5.0, pbc=False)
        wf = SinglepointWorkflow(atoms=atoms, parameters={'functional': 'all', 'n_max_sc_steps': 2, 'from_scratch': True, 'orbital_groups_self_hartree_tol': 100.0},
                                 master_calc_params={'kcp': {'ecutwfc': 20.0, 'nbnd': 5}})
        wf.run()


def test_singlepoint_si_ki_dscf(tmp_path, sys2file):
    with chdir(tmp_path):
        atoms = bulk('Si')
        master_calc_params = {'kcp': {'ecutwfc': 40.0},
                              'pw': {'ecutwfc': 40.0, 'nbnd': 10},
                              'w90_occ': {'conv_window': 5, },
                              'w90_emp': {'conv_window': 5, 'dis_froz_max': 10.6, 'dis_win_max': 16.9},
                              'ui': {'smooth_int_factor': 2},
                              'plot': {'Emin': -10, 'Emax': 4, 'degauss': 0.5}
                              }
        p = [{'fsite': [0.25, 0.25, 0.25], 'ang_mtm': 'sp3'}]
        projs = ProjectionBlocks.fromprojections([p, p], fillings=[True, False], spins=[None, None], atoms=atoms)
        wf = SinglepointWorkflow(
            atoms, parameters={'functional': 'ki', 'method': 'dscf', 'from_scratch': True,
                               'mp_correction': True, 'eps_inf': 13.02, 'init_orbitals': 'mlwfs', 'alpha_guess': 0.077, 'orbital_groups_self_hartree_tol': 100.0},
            kgrid=[2, 2, 2], kpath='GL', projections=projs, master_calc_params=master_calc_params)
        wf.run()
