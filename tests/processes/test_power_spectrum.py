from pathlib import Path

import numpy as np
import pytest
from ase_koopmans import Atoms

from koopmans import utils
from koopmans.bands import Bands
from koopmans.engines.localhost import LocalhostEngine
from koopmans.files import LocalFile
from koopmans.processes.power_spectrum import (
    ComputePowerSpectrumProcess, ExtractCoefficientsFromXMLProcess)


def test_extract_coefficients_from_xml_process(tmpdir, datadir, check_patch):

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

        # copy the required input files from the data directory
        total_density_xml = LocalFile(datadir / 'power_spectrum' / 'charge-density.xml')
        orbital_density_xml = LocalFile(datadir / 'power_spectrum' / 'orbital.occ.0.00001.xml')

        process = ExtractCoefficientsFromXMLProcess(n_max=n_max, l_max=l_max, r_min=r_min, r_max=r_max, r_cut=r_cut,
                                                    wannier_centers=centers, total_density_xml=total_density_xml,
                                                    orbital_densities_xml=[orbital_density_xml], bands=bands_to_extract, cell=water.cell)

        process.directory = Path()
        process.engine = LocalhostEngine()

        process.run()


def test_compute_power_spectrum_process(tmpdir, datadir, check_patch):
    with utils.chdir(tmpdir):

        # set up the test system
        n_max = 6
        l_max = 6

        # copy the required input files from the data directory
        orbital_coefficients = LocalFile(datadir / 'power_spectrum' / 'coeff.orbital.occ.1.npy')
        total_coefficients = LocalFile(datadir / 'power_spectrum' / 'coeff.total.occ.1.npy')

        process = ComputePowerSpectrumProcess(n_max=n_max, l_max=l_max, orbital_coefficients=orbital_coefficients,
                                              total_coefficients=total_coefficients)

        process.directory = Path()
        process.engine = LocalhostEngine()

        process.run()
