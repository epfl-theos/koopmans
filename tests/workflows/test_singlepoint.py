import copy

import pytest

from koopmans import workflows
from koopmans.io import write
from koopmans.kpoints import Kpoints
from koopmans.projections import ProjectionBlocks
from koopmans.utils import chdir


@pytest.mark.espresso
def test_singlepoint_h2o_ki_dscf_explicit(water, espresso_patch, tmp_path, sys2file):
    '''
    Test of H2O calculation of water that explicitly calls Quantum ESPRESSO
    '''
    with chdir(tmp_path):
        parameters = {'functional': 'ki',
                      'alpha_numsteps': 1,
                      'keep_tmpdirs': False,
                      'orbital_groups_self_hartree_tol': 100.0
                      }
        wf = workflows.SinglepointWorkflow(parameters=parameters, **water)
        wf.run()


def test_singlepoint_h2o_all_dscf(water, workflow_patch, tmp_path, sys2file):
    with chdir(tmp_path):
        parameters = {'functional': 'all',
                      'alpha_numsteps': 2,
                      'keep_tmpdirs': False,
                      'orbital_groups_self_hartree_tol': 100.0
                      }
        wf = workflows.SinglepointWorkflow(parameters=parameters, **water)
        wf.run()


@pytest.mark.parametrize('spin_polarized', [True, False])
def test_singlepoint_si_ki_dscf(spin_polarized, silicon, workflow_patch, tmp_path, sys2file):
    with chdir(tmp_path):
        parameters = {'functional': 'ki',
                      'method': 'dscf',
                      'spin_polarized': spin_polarized,
                      'keep_tmpdirs': False,
                      'mp_correction': True,
                      'eps_inf': 13.02,
                      'init_orbitals': 'mlwfs',
                      'alpha_guess': 0.077,
                      'orbital_groups_self_hartree_tol': 100.0,
                      'pseudo_library': 'pseudo_dojo_standard'}

        if spin_polarized:
            projs = silicon.pop('projections')
            spin_projs = []
            for spin in ['up', 'down']:
                for proj in projs:
                    spin_proj = copy.deepcopy(proj)
                    spin_proj.spin = spin
                    spin_projs.append(spin_proj)
            silicon['projections'] = ProjectionBlocks(spin_projs, silicon['atoms'])
        silicon['kpoints'] = Kpoints(grid=[2, 2, 2], path='GL', density=50.0, cell=silicon['atoms'].cell)
        wf = workflows.SinglepointWorkflow(parameters=parameters, **silicon)
        wf.run()


@pytest.mark.espresso
def test_singlepoint_si_ki_dfpt_explicit(silicon, espresso_patch, tmp_path, sys2file):
    with chdir(tmp_path):

        parameters = {'functional': 'ki',
                      'method': 'dfpt',
                      'keep_tmpdirs': False,
                      'eps_inf': 13.02,
                      'init_orbitals': 'mlwfs',
                      'alpha_guess': 0.077,
                      'orbital_groups_self_hartree_tol': 100.0}

        wf = workflows.SinglepointWorkflow(parameters=parameters, **silicon)
        wf.run()


def test_singlepoint_si_ki_dfpt(silicon, workflow_patch, tmp_path, sys2file):
    with chdir(tmp_path):

        parameters = {'functional': 'ki',
                      'method': 'dfpt',
                      'keep_tmpdirs': False,
                      'eps_inf': 13.02,
                      'init_orbitals': 'mlwfs',
                      'alpha_guess': 0.077,
                      'orbital_groups_self_hartree_tol': 100.0}

        wf = workflows.SinglepointWorkflow(parameters=parameters, **silicon)
        wf.run()


def test_singlepoint_ozone_ki_dfpt(ozone, workflow_patch, tmp_path, sys2file):
    with chdir(tmp_path):
        parameters = {'functional': 'ki',
                      'method': 'dfpt',
                      'keep_tmpdirs': False,
                      'init_orbitals': 'kohn-sham',
                      'npool': 1,
                      'orbital_groups_self_hartree_tol': 100.0,
                      'pseudo_library': 'sg15_v1.0'
                      }
        wf = workflows.SinglepointWorkflow(parameters=parameters, **ozone)
        wf.run()


@pytest.mark.espresso
def test_singlepoint_gaas_wan2odd(gaas, espresso_patch, tmp_path, sys2file):
    with chdir(tmp_path):
        parameters = {'functional': 'ki',
                      'method': 'dscf',
                      'keep_tmpdirs': False,
                      'calculate_alpha': False,
                      'init_orbitals': 'mlwfs',
                      'npool': 1,
                      'pseudo_library': 'sg15_v1.0',
                      'calculate_bands': True
                      }
        wf = workflows.SinglepointWorkflow(parameters=parameters, kgrid=[2, 2, 2], **gaas)
        wf.run()
