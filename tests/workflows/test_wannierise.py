from pathlib import Path

import numpy as np
import pytest

from koopmans import workflows
from koopmans.io import write
from koopmans.utils import chdir, read_wannier_hr_file


def test_wannierize_tio2(tio2, tmp_path, sys2file, workflow_patch):
    with chdir(tmp_path):
        parameters = {
            "init_orbitals": "mlwfs",
            "init_empty_orbitals": "projwfs",
            "keep_tmpdirs": False,
            "pseudo_library": "pseudo_dojo_standard"}
        wf = workflows.WannieriseWorkflow(parameters=parameters, kgrid=[2, 2, 2], kpath='GXG', **tio2)
        wf.run()


def test_wannierize_merge_hr_files(tmp_path, datadir):
    with chdir(tmp_path):
        dirs_in = sorted((datadir / 'w90').glob('occ_block*'))
        workflows.WannieriseWorkflow.merge_wannier_hr_files(dirs_in, Path('test'), 'wann')
        ham, rvec, weights, num_wann = read_wannier_hr_file(Path('test/wann_hr.dat'))
        ham_ref, rvec_ref, weights_ref, num_wann_ref = read_wannier_hr_file(datadir / 'w90' / 'occ' / 'wann_hr.dat')

        assert np.allclose(ham, ham_ref)
        assert np.all(rvec == rvec_ref)
        assert weights == weights_ref
        assert num_wann == num_wann_ref