import pytest
import numpy as np
import itertools
from pathlib import Path
from koopmans import utils, workflows

wann_files_dir = Path(__file__).parent / 'w90_example_files'


def test_write_read_wannier_hr_file(tmp_path):
    '''
    Test function for utils.write_wannier_hr_file and utils.read_wannier_hr_file
    '''
    with utils.chdir(tmp_path):
        # Dummy data
        rvec_in = [x for x in itertools.product([-1, 0, 1], repeat=3)]
        n = 5
        ham_in = np.arange(len(rvec_in) * n**2, dtype=complex)
        weights_in = [1 for _ in rvec_in]

        # Write hr file
        utils.write_wannier_hr_file(Path('wann_hr.dat'), ham_in.reshape((len(rvec_in), n, n)), rvec_in, weights_in)

        # Read hr file
        ham_out, rvec_out, weights_out, num_wann_out = utils.read_wannier_hr_file(Path('wann_hr.dat'))

        # Compare the two
        assert np.allclose(ham_out, ham_in)
        assert np.allclose(rvec_out, rvec_in)
        assert weights_out == weights_in
        assert num_wann_out == len(rvec_in)


def test_write_read_wannier_u_file(tmp_path):
    '''
    Test function for utils.write_wannier_u_file and utils.read_wannier_u_file
    '''

    with utils.chdir(tmp_path):
        # Dummy data
        kpts = np.array([[0, 0, 0], [0.5, 0, 0]])
        n_k = len(kpts)
        n_w = 4
        umat = np.arange(n_k * n_w**2, dtype=complex).reshape((n_k, n_w, n_w))

        # Write u file
        utils.write_wannier_u_file('test_u.mat', umat, kpts)

        # Read u file
        umat_out, kpts_out, _ = utils.read_wannier_u_file('test_u.mat')

        # Compare the two
        assert np.allclose(umat, umat_out)
        assert np.allclose(kpts, kpts_out)
