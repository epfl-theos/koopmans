"""

kc_ham calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

import numpy as np
from ase.calculators.espresso import KoopmansHam
from ase.io.espresso import koopmans_ham as kch_io
from koopmans import io, utils
from koopmans.calculators.generic import KCWannCalc, qe_bin_directory
from koopmans.calculators.commands import ParallelCommand


class KoopmansHamCalc(KCWannCalc):
    # Subclass of KCWannCalc for performing calculations with kc_wann.x

    # Point to the appropriate ASE IO module
    _io = kch_io
    _ase_calc_class = KoopmansHam

    # Define the appropriate file extensions
    ext_in = '.khi'
    ext_out = '.kho'

    def __init__(self, *args, **kwargs):
        self.defaults.update({'do_bands': True, 'use_ws_distance': True,
                              'write_hr': True, 'l_alpha_corr': False, 'lrpa': False})
        super().__init__(*args, **kwargs)
        self.calc.command = ParallelCommand(
            f'{qe_bin_directory}kc_ham.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out}')

    def write_alphas(self):
        assert 'alphas' in self.results

        # kc_ham.x takes a single file for the alphas (rather than splitting between filled/empty) and does not have
        # duplicated results for spin up then spin down
        alphas = self.results['alphas']
        filling = [True for _ in range(len(alphas))]
        io.write_alpha_file(self.directory, alphas, filling)

    def get_k_point_weights(self):
        utils.warn('Need to properly define k-point weights')
        return np.ones(len(self.calc.parameters['kpath'].kpts))

    def get_number_of_spins(self):
        return 1

    def get_eigenvalues(self, kpt=None, spin=0):
        if spin != 0:
            raise NotImplementedError(
                f'Koopmans hamiltonian calculator is not implemented for spin-polarised systems')

        if 'band structure' not in self.results:
            raise ValueError('You must first calculate the band structure before you try to access the KS eigenvalues')

        if kpt is None:
            return self.results['band structure'].energies[spin, :]
        else:
            return self.results['band structure'].energies[spin, kpt]

    def get_fermi_level(self):
        return 0
