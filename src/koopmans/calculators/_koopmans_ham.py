"""KCWHam calculator module for koopmans."""

from typing import List, Optional

import numpy as np
from ase_koopmans import Atoms
from ase_koopmans.calculators.espresso import KoopmansHam
from ase_koopmans.dft.kpoints import BandPath

from koopmans import settings, utils

from ._calculator import CalculatorABC, KCWannCalculator, ReturnsBandStructure


class KoopmansHamCalculator(KCWannCalculator, KoopmansHam, ReturnsBandStructure, CalculatorABC):
    """Subclass of KCWannCalculator for calculating the Koopmans Hamiltonian with kcw.x."""

    ext_in = '.khi'
    ext_out = '.kho'
    code = "kcw_ham"

    def __init__(self, atoms: Atoms, alphas: Optional[List[float]] = None, *args, **kwargs):
        # Define the valid settings
        self.parameters = settings.KoopmansHamSettingsDict()

        # Initialize using the ASE parent, and then CalculatorExt
        KoopmansHam.__init__(self, atoms=atoms)
        KCWannCalculator.__init__(self, *args, **kwargs)

        # Store the alphas
        self.alphas = alphas

    def write_alphas(self):
        """Write the screening parameters to a file.

        self.alphas is a list of alpha values indexed by band index. Meanwhile, kcw.x takes a
        single file for the alphas (rather than splitting between filled/empty)
        """
        fake_filling = [True for _ in self.alphas]
        utils.write_alpha_files(self, self.alphas, fake_filling)

    def _pre_calculate(self):
        super()._pre_calculate()
        self.write_alphas()

    def _post_calculate(self):
        super()._post_calculate()
        if isinstance(self.parameters.kpts, BandPath) and len(self.parameters.kpts.kpts) > 1:
            # Add the bandstructure to the results
            self.generate_band_structure()

    def get_k_point_weights(self):
        """Return the k-point weights for the calculation."""
        utils.warn('Need to properly define k-point weights')
        return np.ones(len(self.parameters['kpath'].kpts))

    def get_number_of_spins(self):
        """Return the number of spins in the calculation."""
        return 1

    def get_eigenvalues(self, kpt=None, spin=0):
        """Return the eigenvalues of the computed Hamiltonian."""
        if spin != 0:
            raise NotImplementedError('Koopmans Hamiltonian calculator is not implemented for spin-polarized systems')

        if 'band structure' not in self.results:
            raise ValueError('You must first calculate the band structure before you try to access the KS eigenvalues')

        if kpt is None:
            return self.results['band structure'].energies[spin, :]
        else:
            return self.results['band structure'].energies[spin, kpt]

    def get_fermi_level(self):
        """Return the Fermi level."""
        return 0

    def vbm_energy(self) -> float:
        """Calculate the energy of the valence band maximum."""
        eigenvalues_np = self.eigenvalues_from_results()
        return np.max(eigenvalues_np[:, :, self.parameters.num_wann_occ - 1])

    def eigenvalues_from_results(self):
        """Extract the eigenvalues from the results dictionary."""
        assert 'eigenvalues' in self.results, 'Please call {0}.calculate() prior to calling {0}.band_structure'.format(
            self.__class__.__name__)

        return np.array([self.results['eigenvalues']])

    def is_converged(self):
        """Return True if the calculation is converged."""
        raise NotImplementedError('TODO')

    def check_convergence(self) -> None:
        """Check if the calculation has converged. is_converged has not been implemented yet for this calculator."""
        return

    def read_input(self, **kwargs):
        """Read the input file and store the parameters.

        A .khi file doesn't have the requisite information to reconstruct the bandpath, so in the event that kpts
        are already provided in self.parameters, don't overwrite them."
        """
        kpts = self.parameters.kpts

        super().read_input(**kwargs)

        if kpts is not None:
            self.parameters.kpts = kpts
        return
