import pytest
from ase.build import molecule, bulk
from koopmans.workflows import SinglepointWorkflow
from koopmans.projections import ProjectionBlocks
from koopmans.io import write
from koopmans.utils import chdir

# water
water = molecule('H2O', vacuum=5.0, pbc=False)
water_master_calc_params = {'kcp': {'ecutwfc': 20.0, 'nbnd': 5}}

# bulk silicon
si = bulk('Si')
pdict = [{'fsite': [0.25, 0.25, 0.25], 'ang_mtm': 'sp3'}]
si_projs = ProjectionBlocks.fromprojections([pdict, pdict], fillings=[True, False], spins=[None, None], atoms=si)
si_master_calc_params = {'kcp': {'ecutwfc': 40.0},
                         'pw': {'ecutwfc': 40.0, 'nbnd': 10},
                         'w90_occ': {'conv_window': 5, },
                         'w90_emp': {'conv_window': 5, 'dis_froz_max': 10.6, 'dis_win_max': 16.9},
                         'ui': {'smooth_int_factor': 2},
                         'plot': {'Emin': -10, 'Emax': 4, 'degauss': 0.5}
                         }

# ozone
ozone = molecule('O3', vacuum=5.0, pbc=False)
ozone_master_calc_params = {'pw': {'ecutwfc': 20.0, 'nbnd': 10}}


def test_singlepoint_h2o_all_dscf(mockable, tmp_path, sys2file):
    with chdir(tmp_path):
        parameters = {'functional': 'all',
                      'n_max_sc_steps': 2,
                      'from_scratch': True,
                      'orbital_groups_self_hartree_tol': 100.0
                      }
        wf = SinglepointWorkflow(atoms=water,
                                 parameters=parameters,
                                 master_calc_params=water_master_calc_params)
        wf.run()


def test_singlepoint_si_ki_dscf(mockable, tmp_path, sys2file):
    with chdir(tmp_path):
        parameters = {'functional': 'ki',
                      'method': 'dscf',
                      'from_scratch': True,
                      'mp_correction': True,
                      'eps_inf': 13.02,
                      'init_orbitals': 'mlwfs',
                      'alpha_guess': 0.077,
                      'orbital_groups_self_hartree_tol': 100.0}

        wf = SinglepointWorkflow(atoms=si,
                                 parameters=parameters,
                                 kgrid=[2, 2, 2],
                                 kpath='GL',
                                 projections=si_projs,
                                 master_calc_params=si_master_calc_params)
        wf.run()


def test_singlepoint_si_ki_dfpt(mockable, tmp_path, sys2file):
    with chdir(tmp_path):

        parameters = {'functional': 'ki',
                      'method': 'dfpt',
                      'from_scratch': True,
                      'eps_inf': 13.02,
                      'init_orbitals': 'mlwfs',
                      'alpha_guess': 0.077,
                      'orbital_groups_self_hartree_tol': 100.0}

        wf = SinglepointWorkflow(atoms=si,
                                 parameters=parameters,
                                 kgrid=[2, 2, 2],
                                 kpath='GXG',
                                 projections=si_projs,
                                 master_calc_params=si_master_calc_params)
        wf.run()


def test_singlepoint_ozone_ki_dfpt(mockable, tmp_path, sys2file):
    with chdir(tmp_path):
        parameters = {'functional': 'ki',
                      'method': 'dfpt',
                      'from_scratch': True,
                      'init_orbitals': 'kohn-sham',
                      'npool': 1,
                      'orbital_groups_self_hartree_tol': 100.0,
                      'pseudo_library': 'sg15_v1.0'
                      }
        wf = SinglepointWorkflow(atoms=ozone,
                                 parameters=parameters,
                                 master_calc_params=ozone_master_calc_params)
        wf.run()
