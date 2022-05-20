import pytest
import numpy as np
from pathlib import Path
from ase.spacegroup import crystal
from koopmans.utils import chdir, read_wannier_hr_file
from koopmans.io import write
from koopmans import workflows
from koopmans.projections import ProjectionBlocks


a = 4.6
c = 2.95

# Rutile TiO2:
tio2 = crystal(['Ti', 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
               spacegroup=136, cellpar=[a, a, c, 90, 90, 90])

tio2_projs = ProjectionBlocks.fromprojections([["Ti:l=0"], ["Ti:l=1"], ["O:l=0"], ["O:l=1"], ["Ti:l=0", "Ti:l=2"]],
                                              fillings=[True, True, True, True, False],
                                              spins=[None, None, None, None, None],
                                              atoms=tio2)

tio2_master_calc_params = {'pw': {'ecutwfc': 40, 'nbnd': 36}}


def test_wannierize_tio2(tmp_path, sys2file):
    with chdir(tmp_path):
        parameters = {
            "init_orbitals": "mlwfs",
            "init_empty_orbitals": "projwfs",
            "pseudo_library": "sg15_v1.0"}
        wf = workflows.WannierizeWorkflow(tio2, parameters, tio2_master_calc_params, kgrid=[
            2, 2, 2], kpath='GXG', projections=tio2_projs)
        wf.run()


def test_wannierize_merge_hr_files(tmp_path, datadir):
    with chdir(tmp_path):
        dirs_in = sorted((datadir / 'w90').glob('occ_block*'))
        workflows.WannierizeWorkflow.merge_wannier_hr_files(dirs_in, Path('test'), 'wann')
        ham, rvec, weights, num_wann = read_wannier_hr_file(Path('test/wann_hr.dat'))
        ham_ref, rvec_ref, weights_ref, num_wann_ref = read_wannier_hr_file(datadir / 'w90' / 'occ' / 'wann_hr.dat')

        assert np.allclose(ham, ham_ref)
        assert np.all(rvec == rvec_ref)
        assert weights == weights_ref
        assert num_wann == num_wann_ref
