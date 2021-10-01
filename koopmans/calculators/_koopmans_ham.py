"""

kc_ham calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

import os
import numpy as np
from ase import Atoms
from ase.calculators.espresso import KoopmansHam
from koopmans import utils, settings
from ._utils import KCWannCalculator, qe_bin_directory
from koopmans.commands import ParallelCommand


class KoopmansHamCalculator(KCWannCalculator, KoopmansHam):
    # Subclass of KCWannCalculator for performing calculations with kc_wann.x

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = settings.KoopmansHamSettingsDict()

        # Define the appropriate file extensions
        self.ext_in = '.khi'
        self.ext_out = '.kho'

        # Initialise using the ASE parent, and then ExtendedCalculator
        KoopmansHam.__init__(self, atoms=atoms)
        KCWannCalculator.__init__(*args, **kwargs)

        self.command = ParallelCommand(
            f'{qe_bin_directory}{os.path.sep}kc_ham.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out}')

    def write_alphas(self):
        assert 'alphas' in self.results

        # kc_ham.x takes a single file for the alphas (rather than splitting between filled/empty) and does not have
        # duplicated results for spin up then spin down
        alphas = self.results['alphas']
        filling = [True for _ in range(len(alphas))]
        utils.write_alpha_file(self.directory, alphas, filling)

    def get_k_point_weights(self):
        utils.warn('Need to properly define k-point weights')
        return np.ones(len(self.parameters['kpath'].kpts))

    def get_number_of_spins(self):
        return 1

    def get_eigenvalues(self, kpt=None, spin=0):
        if spin != 0:
            raise NotImplementedError(
                f'Koopmans Hamiltonian calculator is not implemented for spin-polarised systems')

        if 'band structure' not in self.results:
            raise ValueError('You must first calculate the band structure before you try to access the KS eigenvalues')

        if kpt is None:
            return self.results['band structure'].energies[spin, :]
        else:
            return self.results['band structure'].energies[spin, kpt]

    def get_fermi_level(self):
        return 0
