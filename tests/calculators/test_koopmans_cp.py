import shutil

import numpy as np
import pytest

from koopmans import settings, utils, workflows
from koopmans.calculators._koopmans_cp import convert_flat_alphas_for_kcp
from koopmans.testing import benchmark_filename


def test_convert_flat_alphas_for_kcp():
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
    with utils.chdir(tmp_path):
        # Create a kcp calculator to match the one that was used to generate the pdos files
        wf = workflows.KoopmansDSCFWorkflow(**water)
        calc = wf.new_kcp_calculator('ki_final')

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
    with utils.chdir(tmp_path):
        # Create a kcp calculator to match the one that was used to generate the pdos files
        wf = workflows.KoopmansDSCFWorkflow(**water)
        calc = wf.new_kcp_calculator('ki_final')

        # Copy over the XML Hamiltonian files
        destdir = calc.parameters.outdir / f'{calc.parameters.prefix}_{calc.parameters.ndw}.save' / 'K00001'
        destdir.mkdir(parents=True)
        for f in (datadir / 'kcp').glob('ham*'):
            shutil.copy(f, destdir)

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
