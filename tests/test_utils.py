import itertools
from pathlib import Path

import numpy as np
import pytest

from koopmans import utils

wann_files_dir = Path(__file__).parent / 'w90_example_files'


def test_generate_and_parse_wannier_hr_file_contents():
    '''
    Test function for utils.generate_wannier_hr_file_contents and utils.parse_wannier_hr_file_contents
    '''
    # Dummy data
    rvec_in = [x for x in itertools.product([-1, 0, 1], repeat=3)]
    n = 5
    ham_in = np.arange(len(rvec_in) * n**2, dtype=complex)
    weights_in = [1 for _ in rvec_in]

    # Generate hr file contents
    filecontents = utils.generate_wannier_hr_file_contents(ham_in.reshape((len(rvec_in), n, n)), rvec_in, weights_in)

    # Parse hr file contents
    ham_out, rvec_out, weights_out, num_wann_out = utils.parse_wannier_hr_file_contents(filecontents)

    # Compare the two
    assert np.allclose(ham_out, ham_in)
    assert np.allclose(rvec_out, rvec_in)
    assert weights_out == weights_in
    assert num_wann_out == len(rvec_in)


def test_generate_and_parse_wannier_u_file_contents():
    '''
    Test function for utils.generate_wannier_u_file_contents and utils.parse_wannier_u_file_contents
    '''

    # Dummy data
    kpts = np.array([[0, 0, 0], [0.5, 0, 0]])
    n_k = len(kpts)
    n_w = 4
    umat = np.arange(n_k * n_w**2, dtype=complex).reshape((n_k, n_w, n_w))

    # Generate u file contents
    filecontents = utils.generate_wannier_u_file_contents(umat, kpts)

    # Read u file
    umat_out, kpts_out, _ = utils.parse_wannier_u_file_contents(filecontents)

    # Compare the two
    assert np.allclose(umat, umat_out)
    assert np.allclose(kpts, kpts_out)


def test_parse_wannier_centers_file_contents(tmp_path, datadir):
    '''
    Test function for utils.parse_wannier_centers_file_contents
    '''
    with utils.chdir(tmp_path):
        xyz_file = (datadir / 'w90' / 'example_centres.xyz').resolve()

        with open(xyz_file) as f:
            contents = f.readlines()

        centers, atoms = utils.parse_wannier_centers_file_contents(contents)

        assert all(atoms.symbols == 'C12H12')


def test_generate_and_parse_wannier_centers_file(silicon, tmp_path, datadir):
    '''
    Test function for utils.read_wannier_centers_file and utils.write_wannier_centers_file
    '''
    with utils.chdir(tmp_path):
        # Use silicon as the atoms
        atoms = silicon['atoms']

        # Generate some pseudo-random centers
        np.random.seed(100)
        centers = 10 * np.random.random_sample((3, 10))

        # Generate the centers file contents
        contents = utils.generate_wannier_centers_file_contents(centers, atoms)

        # Parse the centers file contents
        centers_in, atoms_in = utils.parse_wannier_centers_file_contents(contents)

        # Check that the centers and atoms are unchanged
        assert np.allclose(centers_in, centers)
        assert np.allclose(atoms.positions, atoms_in.positions)
        assert all(atoms.symbols == atoms_in.symbols)
