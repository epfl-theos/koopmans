from fileinput import filename

import numpy as np

from ase import Atoms
from koopmans import base_directory, calculators, ml_utils, utils
from koopmans.bands import Bands
from koopmans.io import read_kwf as read_encoded_json
from koopmans.io import write_kwf as write_encoded_json


def test_compute_decomposition(tmpdir, datadir, pytestconfig):
    with utils.chdir(tmpdir):
        n_max = 6
        l_max = 6
        r_min = 0.5
        r_max = 6.0

        water = Atoms('H2O', positions=[(3.2774, 4.01205, 3.61285),
                      (3.6068, 2.88085, 2.64155), (3.0000, 3.11915, 3.35845)])
        water.set_cell(np.diag([6.8929, 6.8929, 6.8929]))

        # real values: num_bands_occ=4, n_bands=6
        num_bands_occ = 1
        n_bands = 1
        bands_to_extract = Bands(n_bands=n_bands, n_spin=1, spin_polarized=False)
        for i, band in enumerate(bands_to_extract):
            if i >= num_bands_occ:
                band.filled = False

        r_cut = min(water.get_cell_lengths_and_angles()[:3])
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
            'power': ml_dir / ('power_spectra_' + dir_suffix)
        }
        for dir in dirs.values():
            utils.system_call(f'mkdir -p {dir}')  # utils mkdir

        centers = np.array([[3.159166, -3.286943, -3.412441],
                            [3.332131,  3.005906,  2.967669],
                            [2.713486,  3.144959,  3.239604],
                            [3.073549,  2.939136, -3.295453],
                            [-3.333020, -2.131308, -3.034869],
                            [-2.822225,  2.595511,  2.101965]])

        with open(dirs['xml'] / 'bands_to_solve.txt', "w") as f:
            f.write(f"{n_bands}\n")
            for i in range(n_bands):
                f.write(
                    f"{bands_to_extract[i].index}, {int(bands_to_extract[i].filled)}, {bands_to_extract[i].spin}\n")

        command = str(calculators.bin_directory / 'bin2xml_real_space_density.x ') + ' '.join(str(x)
                                                                                              for x in [datadir / 'real_space_densities', dirs['xml'], num_bands_occ])
        utils.system_call(command)

        ml_utils.precompute_parameters_of_radial_basis(n_max, l_max, r_min, r_max, dirs)
        ml_utils.compute_decomposition(n_max, l_max, r_min, r_max, r_cut, dirs, bands_to_extract, water, centers)
        input_vectors_for_ml = {}
        ml_utils.compute_power(n_max, l_max, dirs, bands_to_extract, input_vectors_for_ml)

        power = np.loadtxt(dirs['power'] / 'power_spectrum.orbital.occ.1.txt')

        benchmark_dir = base_directory / 'tests' / 'benchmarks'
        benchmark_file = benchmark_dir / 'test_compute_decomposition.json'

        if pytestconfig.getoption('generate_benchmark'):
            # Write the power spectrum tensor to file
            with open(benchmark_file, 'w') as fd:
                write_encoded_json(power, fd)
        else:
            # Compare with the power spectrum on file
            with open(benchmark_file, 'r') as fd:
                power_ref = read_encoded_json(fd)
            assert np.allclose(power, power_ref)
