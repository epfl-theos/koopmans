import os
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms, units

from koopmans import calculators, ml_utils, utils
from koopmans.bands import Bands
from koopmans.io import read_kwf as read_encoded_json
from koopmans.io import write_kwf as write_encoded_json
from koopmans.ml_utils import (load_density_into_array,
                               load_grid_dimension_from_xml_file)

benchmark_dir = Path(__file__).parents[1] / 'benchmarks'


@pytest.mark.espresso
def test_bin2xml(tmpdir, datadir, pytestconfig):
    with utils.chdir(tmpdir):
        n_bands = 1
        num_bands_occ = 1

        # only check for one occupied orbital
        bands_to_extract = Bands(n_bands=n_bands, n_spin=1, spin_polarized=False)

        # set up the xml directory
        ml_dir = tmpdir / 'ml' / 'TMP'
        dirs_xml = ml_dir / 'xml'
        os.makedirs(dirs_xml)

        # write the band we want to solve into a file
        with open(dirs_xml / 'bands_to_solve.txt', "w") as f:
            f.write(f"{n_bands}\n")
            f.write(f"{bands_to_extract[0].index}, {int(bands_to_extract[0].filled)}, {bands_to_extract[0].spin}\n")

        # extract the xml file from the binary
        command = ' '.join(str(x) for x in ['bin2xml_real_space_density.x', datadir / 'ml', dirs_xml, num_bands_occ])
        utils.system_call(command)

        # check if the extracted xml-file matches the reference
        filename_total_density = 'charge-density.xml'
        filename_orbital_density = 'orbital.occ.0.00001.xml'

        if pytestconfig.getoption('generate_benchmark'):
            shutil.copy(dirs_xml / filename_total_density, benchmark_dir / filename_total_density)
            shutil.copy(dirs_xml / filename_orbital_density, benchmark_dir / filename_orbital_density)
        else:
            nr_xml = load_grid_dimension_from_xml_file(benchmark_dir / filename_total_density)
            norm_const = 1/(units.Bohr)**3

            result_orbital, _ = load_density_into_array(dirs_xml / filename_orbital_density, nr_xml, norm_const)
            ref_orbital, _ = load_density_into_array(benchmark_dir / filename_orbital_density, nr_xml, norm_const)
            result_total, _ = load_density_into_array(
                dirs_xml / filename_total_density, nr_xml, norm_const, 'CHARGE-DENSITY')
            ref_total, _ = load_density_into_array(
                benchmark_dir / filename_total_density, nr_xml, norm_const, 'CHARGE-DENSITY')

            assert np.allclose(result_orbital, ref_orbital)
            assert np.allclose(result_total, ref_total)


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
        ml_utils.precompute_parameters_of_radial_basis(n_max, l_max, r_min, r_max, dirs)
        ml_utils.compute_decomposition(n_max, l_max, r_min, r_max, r_cut, dirs, bands_to_extract, water, centers)

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
        ml_utils.compute_power(n_max, l_max, dirs, bands_to_extract, input_vectors_for_ml)

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
