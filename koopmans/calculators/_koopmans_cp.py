"""

kcp calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os
import numpy as np
from pandas.core.series import Series
from ase.io.espresso import koopmans_cp as kcp_io
from koopmans import utils, settings
from ._utils import ExtendedCalculator, kcp_bin_directory
from koopmans.commands import ParallelCommand


class KoopmansCPCalculator(ExtendedCalculator, kcp_io.Espresso_kcp):
    # Subclass of ExtendedCalculator for performing calculations with kcp.x

    def __init__(self, calc=None, qe_files=[], skip_qc=False, alphas=None, filling=None, **kwargs):

        defaults = {'calculation': 'cp',
                    'outdir': './TMP-CP/',
                    'iprint': 1,
                    'prefix': 'kc',
                    'verbosity': 'low',
                    'disk_io': 'high',
                    'write_hr': False,
                    'do_wf_cmplx': True,
                    'do_ee': True,
                    'electron_dynamics': 'cg',
                    'nspin': 2,
                    'ortho_para': 1,
                    'passop': 2.0,
                    'ion_dynamics': 'none',
                    'ion_nstepe': 5,
                    'ion_radius(1)': 1.0,
                    'ion_radius(2)': 1.0,
                    'ion_radius(3)': 1.0,
                    'ion_radius(4)': 1.0,
                    'do_innerloop_cg': True,
                    'innerloop_cg_nreset': 20,
                    'innerloop_cg_nsd': 2,
                    'innerloop_init_n': 3,
                    'innerloop_nmax': 100,
                    'hartree_only_sic': False,
                    'empty_states_nbnd': 0,
                    'conv_thr': '1.0e-9*nelec',
                    'esic_conv_thr': '1.0e-9*nelec'}

        self.parameters = settings.SettingsDict(valid=[k for sublist in kcp_io.KEYS.values() for k in sublist],
                                                defaults=defaults,
                                                are_paths=['outdir', 'pseudo_dir'],
                                                to_not_parse=['assume_isolated'])

        # Define the appropriate file extensions
        self.ext_in = '.cpi'
        self.ext_out = '.cpo'

        self._ase_calc_class = kcp_io.Espresso_kcp

        super().__init__(calc, qe_files, skip_qc, **kwargs)

        self.results_for_qc = ['energy', 'homo_energy', 'lumo_energy']
        if not isinstance(self.command, ParallelCommand):
            self.command = ParallelCommand(os.environ.get(
                'ASE_ESPRESSO_KCP_COMMAND', kcp_bin_directory + self.command))

        self.alphas = alphas
        self.filling = filling

    def is_complete(self):
        return self.results['job_done']

    def is_converged(self):
        # Checks convergence of the calculation
        if 'conv_thr' not in self.parameters:
            raise ValueError(
                'Cannot check convergence when "conv_thr" is not set')
        return self._ase_is_converged()

    def _ase_is_converged(self):
        if 'convergence' not in self.results:
            raise ValueError(
                'Could not locate calculation details to check convergence')

        # Check convergence for both filled and empty, allowing for the possibility
        # of do_outerloop(_empty) = False meaning the calculation is immediately
        # 'converged'
        do_outers = [self.parameters.do_outerloop, self.parameters.do_outerloop_empty]
        convergence_data = self.results['convergence'].values()
        converged = []
        for do_outer, convergence in zip(do_outers, convergence_data):
            if not do_outer:
                converged.append(True)
            else:
                converged.append(
                    convergence[-1]['delta_E'] < self.parameters.conv_thr * utils.units.Hartree)
        return all(converged)

    @property
    def alphas(self):
        return self._alphas

    @alphas.setter
    def alphas(self, val):

        if val is None:
            return
        elif isinstance(val, Series):
            val = val.to_numpy()

        # alphas can be a 1D or 2D array. If it is 1D, convert it to 2D so that self.alphas
        # is indexed by [i_spin, i_orbital]
        if isinstance(val[0], float):
            val = [val for _ in range(self.parameters.nspin)]

        if len(val) == 1 and self.nspin == 2:
            val = [val[0] for _ in range(self.parameters.nspin)]

        self._alphas = val

    @property
    def filling(self):
        # Filling is indexed by [i_spin, i_orbital]
        # Filling is written in this way such that we can calculate it automatically,
        # but if it is explicitly set then we will always use that value instead
        if self._filling is None:
            n_filled_bands = self.parameters.nelec // 2
            n_empty_bands = self.parameters.empty_states_nbnd
            if n_empty_bands is None:
                n_empty_bands = 0
            filled_spin_channel = [True for _ in range(n_filled_bands)] + [False for _ in range(n_empty_bands)]
            return [filled_spin_channel for _ in range(self.parameters.nspin)]
        else:
            return self._filling

    @filling.setter
    def filling(self, val):

        if val is not None:
            # val can be a 1D or 2D array. If it is 1D, convert it to 2D so that filling
            # is indexed by [i_spin, i_orbital]
            if isinstance(val[0], bool):
                val = [val for _ in range(self.nspin)]

        self._filling = val

    def write_alphas(self):
        '''
        Generates file_alpharef.txt and file_alpharef_empty.txt

        '''

        if not self.do_orbdep or not self.odd_nkscalfact:
            return

        flat_alphas = [a for sublist in self.alphas for a in sublist]
        flat_filling = [f for sublist in self.filling for f in sublist]
        utils.write_alpha_file(self.directory, flat_alphas, flat_filling)

    def read_alphas(self):
        '''
        Reads in file_alpharef.txt and file_alpharef_empty.txt from this calculation's directory

        Output:
           alphas -- a list of alpha values (1 per orbital)
        '''

        if not self.do_orbdep or not self.odd_nkscalfact:
            return

        alphas = utils.read_alpha_file(self.directory)

        if self.nspin == 2:
            # Remove duplicates
            alphas = alphas[:self.nelec // 2] + alphas[self.nelec:self.nelec + self.empty_states_nbnd]
            alphas = [alphas, alphas]
        return alphas

    def transform_to_supercell(self, matrix, **kwargs):
        super().transform_to_supercell(matrix, **kwargs)

        # Also multiply all extensive properties by the appropriate prefactor
        prefactor = np.prod(np.diag(matrix))
        for attr in ['nelec', 'nelup', 'neldw', 'empty_states_nbnd', 'conv_thr', 'esic_conv_thr']:
            value = getattr(self, attr, None)
            if value is not None:
                setattr(self, attr, prefactor * value)

    # The following functions enable DOS generation via ase.dft.dos.DOS(<KoopmansCPCalculator object>)
    def get_k_point_weights(self):
        return [1]

    def get_number_of_spins(self):
        return 1

    def get_eigenvalues(self, kpt=None, spin=0):
        if 'eigenvalues' not in self.results:
            raise ValueError('You must first perform a calculation before you try to access the KS eigenvalues')

        if kpt is None:
            return [self.results['eigenvalues'][spin]]
        elif kpt == 0:
            return self.results['eigenvalues'][spin]
        else:
            raise ValueError(f'{self.__class__.__name__} does not have k-point-resolved KS eigenvalues')

    def get_fermi_level(self):
        return 0
