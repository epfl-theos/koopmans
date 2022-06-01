import pytest
from koopmans.workflows import SinglepointWorkflow
from koopmans.projections import ProjectionBlocks
from koopmans.io import write
from koopmans.utils import chdir


@pytest.mark.espresso
def test_singlepoint_h2o_ki_dscf_explicit(water, espresso_patch, tmp_path, sys2file):
    '''
    Test of H2O calculation of water that explicitly calls Quantum ESPRESSO
    '''
    with chdir(tmp_path):
        parameters = {'functional': 'ki',
                      'n_max_sc_steps': 1,
                      'from_scratch': True,
                      'orbital_groups_self_hartree_tol': 100.0
                      }
        wf = SinglepointWorkflow(parameters=parameters, **water)
        wf.run()


def test_singlepoint_h2o_all_dscf(water, workflow_patch, tmp_path, sys2file):
    with chdir(tmp_path):
        parameters = {'functional': 'all',
                      'n_max_sc_steps': 2,
                      'from_scratch': True,
                      'orbital_groups_self_hartree_tol': 100.0
                      }
        wf = SinglepointWorkflow(parameters=parameters, **water)
        wf.run()


def test_singlepoint_si_ki_dscf(silicon, workflow_patch, tmp_path, sys2file):
    with chdir(tmp_path):
        parameters = {'functional': 'ki',
                      'method': 'dscf',
                      'from_scratch': True,
                      'mp_correction': True,
                      'eps_inf': 13.02,
                      'init_orbitals': 'mlwfs',
                      'alpha_guess': 0.077,
                      'orbital_groups_self_hartree_tol': 100.0}

        wf = SinglepointWorkflow(parameters=parameters, kgrid=[2, 2, 2], kpath='GL', kpath_density=50, **silicon)
        wf.run()


@pytest.mark.espresso
def test_singlepoint_si_ki_dfpt_explicit(silicon, espresso_patch, tmp_path, sys2file):
    with chdir(tmp_path):

        parameters = {'functional': 'ki',
                      'method': 'dfpt',
                      'from_scratch': True,
                      'eps_inf': 13.02,
                      'init_orbitals': 'mlwfs',
                      'alpha_guess': 0.077,
                      'orbital_groups_self_hartree_tol': 100.0}

        wf = SinglepointWorkflow(parameters=parameters, kgrid=[2, 2, 2], kpath='GXG', **silicon)
        wf.run()


def test_singlepoint_si_ki_dfpt(silicon, workflow_patch, tmp_path, sys2file):
    with chdir(tmp_path):

        parameters = {'functional': 'ki',
                      'method': 'dfpt',
                      'from_scratch': True,
                      'eps_inf': 13.02,
                      'init_orbitals': 'mlwfs',
                      'alpha_guess': 0.077,
                      'orbital_groups_self_hartree_tol': 100.0}

        wf = SinglepointWorkflow(parameters=parameters, kgrid=[2, 2, 2], kpath='GXG', **silicon)
        wf.run()


def test_singlepoint_ozone_ki_dfpt(ozone, workflow_patch, tmp_path, sys2file):
    with chdir(tmp_path):
        parameters = {'functional': 'ki',
                      'method': 'dfpt',
                      'from_scratch': True,
                      'init_orbitals': 'kohn-sham',
                      'npool': 1,
                      'orbital_groups_self_hartree_tol': 100.0,
                      'pseudo_library': 'sg15_v1.0'
                      }
        wf = SinglepointWorkflow(parameters=parameters, **ozone)
        wf.run()
