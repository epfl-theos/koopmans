import pytest
from ase.build import molecule
from koopmans.workflows import SinglepointWorkflow
from koopmans.io import write
from koopmans.utils import chdir


def test_singlepoint_h2o_all_dscf(tmp_path):
    with chdir(tmp_path):
        atoms = molecule('H2O', vacuum=5.0, pbc=False)
        wf = SinglepointWorkflow(atoms=atoms, parameters={'functional': 'all', 'n_max_sc_steps': 2, 'from_scratch': True, 'orbital_groups_self_hartree_tol': 100.0}, master_calc_params={
                                 'kcp': {'ecutwfc': 20.0, 'nbnd': 5}})
        write(wf, wf.name + '.json')
        wf.run()
