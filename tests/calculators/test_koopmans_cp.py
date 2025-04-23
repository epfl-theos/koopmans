"""Tests for the KoopmansCP calculator."""

from pathlib import Path

import numpy as np
import pytest  # noqa
from hypothesis import given
from hypothesis import strategies as st

from koopmans import settings, utils
from koopmans.calculators import KoopmansCPCalculator
from koopmans.calculators._koopmans_cp import (allowed,
                                               convert_flat_alphas_for_kcp,
                                               good_fft)
from koopmans.engines import LocalhostEngine
from koopmans.files import LocalFile
from koopmans.pseudopotentials import local_base_directory


def test_convert_flat_alphas_for_kcp():
    """Test the conversion of flat alphas to the correct shape for the KoopmansCP calculator."""
    nbnd = 10
    nspin = 2
    nelup = 5
    neldw = 4
    nelec = nelup + neldw
    flat_alphas = [float(x) for x in range(nbnd * nspin)]
    parameters = settings.KoopmansCPSettingsDict(nbnd=nbnd, nelup=nelup, neldw=neldw, nspin=nspin, nelec=nelec)

    alphas = convert_flat_alphas_for_kcp(flat_alphas, parameters)

    np_alphas = np.array(alphas)
    assert np_alphas.shape == (nspin, nbnd)
    assert np.allclose(np_alphas[0, :nelup], flat_alphas[:nelup])
    assert np.allclose(np_alphas[1, :neldw], flat_alphas[nelup:nelec])
    assert np.allclose(np_alphas[0, nelup:], flat_alphas[nelec:nbnd + neldw])
    assert np.allclose(np_alphas[1, neldw:], flat_alphas[nbnd + neldw:])


def test_read_write_ham_pkl(water, tmp_path):
    """Test reading and writing Hamiltonian pkl files."""
    with utils.chdir(tmp_path):
        # Create a kcp calculator
        pseudo_library = water.pop('pseudo_library')
        pseudo_dir = local_base_directory / pseudo_library
        calc = KoopmansCPCalculator(outdir='tmp', nspin=2, **water, pseudo_dir=pseudo_dir)
        calc.directory = Path()

        # generate a random array for our "Hamiltonian", making sure to set the random seed in order to always
        # generate the same random array
        np.random.seed(100)
        ham = np.random.random_sample((calc.parameters.nbnd, calc.parameters.nbnd))
        + 1.0j * np.random.random_sample((calc.parameters.nbnd, calc.parameters.nbnd))

        # Write the hamiltonian to file and then read it in
        calc.write_ham_pkl_files([ham])
        [ham_in] = calc.read_ham_pkl_files()

        # Check that the two are the same
        assert np.allclose(ham, ham_in)


def test_read_ham(water, datadir, tmp_path):
    """Test reading Hamiltonian files."""
    with utils.chdir(tmp_path):
        # Create a kcp calculator
        pseudo_library = water.pop('pseudo_library')
        pseudo_dir = local_base_directory / pseudo_library
        calc = KoopmansCPCalculator(outdir='tmp', nspin=2, nelec=8, ndw=50,
                                    prefix='test_read_ham', **water, pseudo_dir=pseudo_dir)
        calc.directory = Path()
        calc.engine = LocalhostEngine()

        # Copy over the XML Hamiltonian files
        destdir = calc.write_directory / 'K00001/'
        destdir.mkdir(parents=True)
        for f in LocalFile(datadir / 'kcp').glob('ham*'):
            calc.engine.copy_file(f, destdir)

        # Read the XML Hamiltonian files (and in so doing, write them in pkl format)
        screened_lambda = calc.read_ham_files()
        bare_lambda = calc.read_ham_files(bare=True)
        for ham in [screened_lambda, bare_lambda]:
            for sham in ham:
                assert sham.shape == (calc.parameters.nbnd, calc.parameters.nbnd)

        # Read the pkl files
        screened_lambda_pkl = calc.read_ham_files()
        bare_lambda_pkl = calc.read_ham_files(bare=True)

        assert np.allclose(screened_lambda, screened_lambda_pkl)
        assert np.allclose(bare_lambda, bare_lambda_pkl)


small_ints = st.integers(min_value=0, max_value=2)


@given(small_ints, small_ints, small_ints, small_ints, small_ints, small_ints)
def test_allowed_nr(p2: int, p3: int, p5: int, p7: int, p11: int, p13: int):
    """Return True if nr has no prime factors greater than 5."""
    nr = 2 ** p2 * 3 ** p3 * 5 ** p5 * 7 ** p7 * 11 ** p11 * 13 * p13

    allow = allowed(nr)
    assert allow == ((p7 == 0 and p11 == 0 and p13 == 0) and nr != 0)


small_pos_ints = st.integers(min_value=1, max_value=2)


@given(small_pos_ints, small_pos_ints, small_pos_ints, small_pos_ints, small_pos_ints, small_pos_ints)
def test_good_fft(p2: int, p3: int, p5: int, p7: int, p11: int, p13: int):
    """Return an allowed integer."""
    nr = 2 ** p2 * 3 ** p3 * 5 ** p5 * 7 ** p7 * 11 ** p11 * 13 * p13

    # This nr should not be allowed
    assert not allowed(nr)

    # Find a good nr
    nr_new = good_fft(nr)

    # Check that new_nr is indeed good
    if nr <= 2050:
        assert allowed(nr_new)
        assert nr_new <= 2050
    else:
        assert nr_new == nr
