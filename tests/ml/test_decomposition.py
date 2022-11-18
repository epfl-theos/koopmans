import os
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

from koopmans import ml, utils
from koopmans.bands import Bands
from koopmans.io import read_kwf as read_encoded_json
from koopmans.io import write_kwf as write_encoded_json

benchmark_dir = Path(__file__).parents[1] / 'benchmarks'


def test_decomposition(tmpdir, datadir, pytestconfig):
    with utils.chdir(tmpdir):

        # set up the test system
        n_max = 6
        l_max = 6
        r_min = 0.5
        r_max = 6.0

        water = Atoms('H2O', positions=[(3.2774, 4.01205, 3.61285),
                      (3.6068, 2.88085, 2.64155), (3.0000, 3.11915, 3.35845)])
        water.set_cell(np.diag([6.8929, 6.8929, 6.8929]))

        r_cut = min(water.get_cell_lengths_and_angles()[:3])

        # say which bands we want to extract
        n_bands = 1
        bands_to_extract = Bands(n_bands=n_bands, n_spin=1, spin_polarized=False)

        # provide the Wannier centers for the decomposition
        centers = np.array([[3.159166, -3.286943, -3.412441]])

        # set up the directories
        ml_dir = tmpdir / 'ml' / 'TMP'
        dir_suffix = '_'.join(str(x) for x in [n_max, l_max, r_min, r_max])
        dirs = {
            'ml': ml_dir,
            'xml': ml_dir / 'xml',
            'alphas': ml_dir / 'alphas',
            'betas': ml_dir / 'betas',
            'coeff': ml_dir / ('coefficients_' + dir_suffix),
            'coeff_orb': ml_dir / ('coefficients_' + dir_suffix) / 'coeff_orb',
            'coeff_tot': ml_dir / ('coefficients_' + dir_suffix) / 'coeff_tot',
        }
        for dir in dirs.values():
            os.makedirs(dir)

        # copy the required input files from the data directory
        filename_total_density = 'charge-density.xml'
        filename_orbital_density = 'orbital.occ.0.00001.xml'

        shutil.copy(datadir / 'ml' / filename_total_density, dirs['xml'] / filename_total_density)
        shutil.copy(datadir / 'ml' / filename_orbital_density, dirs['xml'] / filename_orbital_density)

        # compute the coefficient vectors
        ml.precompute_parameters_of_radial_basis(n_max, l_max, r_min, r_max, dirs)
        ml.compute_decomposition(n_max, l_max, r_min, r_max, r_cut, dirs, bands_to_extract, water, centers)

        # check if the coefficient vectors match the reference solution
        benchmark_file_orb = benchmark_dir / 'test_compute_decomposition_orb.json'
        benchmark_file_tot = benchmark_dir / 'test_compute_decomposition_tot.json'
        coeff_orb = np.loadtxt(dirs['coeff_orb'] / 'coff.orbital.occ.1.txt')
        coeff_tot = np.loadtxt(dirs['coeff_tot'] / 'coff.total.occ.1.txt')

        if pytestconfig.getoption('generate_benchmark'):
            with open(benchmark_file_orb, 'w') as fd:
                write_encoded_json(coeff_orb, fd)
            with open(benchmark_file_tot, 'w') as fd:
                write_encoded_json(coeff_tot, fd)
        else:
            with open(benchmark_file_orb, 'r') as fd:
                coeff_orb_ref = read_encoded_json(fd)
            with open(benchmark_file_tot, 'r') as fd:
                coeff_tot_ref = read_encoded_json(fd)

            assert np.allclose(coeff_orb, coeff_orb_ref)
            assert np.allclose(coeff_tot, coeff_tot_ref)


def test_compute_power(tmpdir, datadir, pytestconfig):
    with utils.chdir(tmpdir):

        # set up the test system
        n_max = 6
        l_max = 6
        r_min = 0.5
        r_max = 6.0

        # specify which bands we want to extract
        n_bands = 1
        bands_to_extract = Bands(n_bands=n_bands, n_spin=1, spin_polarized=False)

        # set up the directories
        ml_dir = tmpdir / 'ml' / 'TMP'
        dir_suffix = '_'.join(str(x) for x in [n_max, l_max, r_min, r_max])
        dirs = {
            'ml': ml_dir,
            'coeff_orb': ml_dir / ('coefficients_' + dir_suffix) / 'coeff_orb',
            'coeff_tot': ml_dir / ('coefficients_' + dir_suffix) / 'coeff_tot',
            'power': ml_dir / ('power_spectra_' + dir_suffix)
        }
        for dir in dirs.values():
            os.makedirs(dir)

        # copy the required input files from the data directory
        filename_coeff_tot = 'coff.total.occ.1.txt'
        filename_coeff_orb = 'coff.orbital.occ.1.txt'

        shutil.copy(datadir / 'ml' / filename_coeff_tot, dirs['coeff_tot'] / filename_coeff_tot)
        shutil.copy(datadir / 'ml' / filename_coeff_orb, dirs['coeff_orb'] / filename_coeff_orb)

        # compute the power spectra
        input_vectors_for_ml = {}
        ml.compute_power(n_max, l_max, dirs, bands_to_extract, input_vectors_for_ml)

        # check if the computed power spectra match the reference solution
        power = np.loadtxt(dirs['power'] / 'power_spectrum.orbital.occ.1.txt')

        benchmark_file = benchmark_dir / 'test_compute_power.json'

        if pytestconfig.getoption('generate_benchmark'):
            with open(benchmark_file, 'w') as fd:
                write_encoded_json(power, fd)
        else:
            with open(benchmark_file, 'r') as fd:
                power_ref = read_encoded_json(fd)
            assert np.allclose(power, power_ref)
