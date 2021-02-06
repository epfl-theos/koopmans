"""

kcp calculator module for python_KI

Written by Edward Linscott Sep 2020

"""

import numpy as np
from ase.io import espresso_kcp as kcp_io
from koopmans import io, utils
from koopmans.calculators.generic import QE_calc


class KCP_calc(QE_calc):
    # Subclass of QE_calc for performing calculations with kcp.x

    # Point to the appropriate ASE IO module
    _io = kcp_io

    # Define the appropriate file extensions
    ext_in = '.cpi'
    ext_out = '.cpo'

    # Adding all kcp.x keywords as decorated properties of the KCP_calc class.
    # This means one can set and get kcp.x keywords as self.<keyword> but
    # internally they are stored as self._settings['keyword'] rather than
    # self.<keyword>
    _recognised_keywords = []

    for keywords in _io.KEYS.values():
        for k in keywords:
            _recognised_keywords.append(k)

            # We need to use these make_get/set functions so that get/set_k are
            # evaluated immediately (otherwise we run into late binding and 'k'
            # is not defined when get/set_k are called)
            def make_get(key):
                def get_k(self):
                    # Return 'None' rather than an error if the keyword has not
                    # been defined
                    return self._settings.get(key, None)
                return get_k

            def make_set(key):
                def set_k(self, value):
                    self._settings[key] = value
                return set_k

            get_k = make_get(k)
            set_k = make_set(k)
            locals()[k] = property(get_k, set_k)

    def __init__(self, calc=None, qe_files=[], skip_qc=False, alphas=None, filling=None, **kwargs):
        self.settings_to_not_parse = ['pseudo_dir', 'assume_isolated']

        super().__init__(calc, qe_files, skip_qc, **kwargs)

        self.results_for_qc = ['energy', 'homo_energy', 'lumo_energy']

        self.alphas = alphas
        self.filling = filling

    @property
    def calc(self):
        # First, update the param block
        self._ase_calc.parameters['input_data'] = self.construct_namelist()

        return self._ase_calc

    @calc.setter
    def calc(self, value):
        self._ase_calc = value

    def is_complete(self):
        return self.results['job_done']

    def is_converged(self):
        # Checks convergence of the calculation
        if self.conv_thr is None:
            raise ValueError(
                'Cannot check convergence when "conv_thr" is not set')
        return self._ase_is_converged()

    def _ase_calculate(self):
        # Before running the calculation, update the keywords for the ASE calculator object
        self._ase_calc.parameters['input_data'] = self.construct_namelist()
        super()._ase_calculate()

    def construct_namelist(self):
        # Returns a namelist of settings, grouped by their Quantum Espresso headings
        return kcp_io.construct_namelist(**self._settings, warn=True)

    def _update_settings_dict(self):
        # Updates self._settings based on self._ase_calc
        self._settings = {}
        for namelist in self._ase_calc.parameters.get('input_data', {}).values():
            for key, val in namelist.items():
                self._settings[key] = val

    def _ase_is_converged(self):
        if 'convergence' not in self.results:
            raise ValueError(
                'Could not locate calculation details to check convergence')

        # Check convergence for both filled and empty, allowing for the possibility
        # of do_outerloop(_empty) = False meaning the calculation is immediately
        # 'converged'
        do_outers = [self.do_outerloop, self.do_outerloop_empty]
        convergence_data = self.results['convergence'].values()
        converged = []
        for do_outer, convergence in zip(do_outers, convergence_data):
            if not do_outer:
                converged.append(True)
            else:
                converged.append(
                    convergence[-1]['delta_E'] < self.conv_thr * utils.units.Hartree)
        return all(converged)

    @property
    def alphas(self):
        return self._alphas

    @alphas.setter
    def alphas(self, val):

        if val is None:
            return

        # convert to an array
        val = np.array(val)

        # alphas can be a 1D or 2D array. If it is 1D, convert it to 2D so that self.alphas
        # is indexed by [i_spin, i_orbital]
        if isinstance(val[0], float):
            val = [val for _ in range(self.nspin)]

        if len(val) == 1 and self.nspin == 2:
            val = [val[0] for _ in range(self.nspin)]

        self._alphas = val

    @property
    def filling(self):
        # Filling is indexed by [i_spin, i_orbital]
        # Filling is written in this way such that we can calculate it automatically,
        # but if it is explicitly set then we will always use that value instead
        if self._filling is None:
            n_filled_bands = self.nelec // 2
            n_empty_bands = self.empty_states_nbnd
            if n_empty_bands is None:
                n_empty_bands = 0
            filled_spin_channel = [True for _ in range(n_filled_bands)] + [False for _ in range(n_empty_bands)]
            return [filled_spin_channel for _ in range(self.nspin)]
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
        io.write_alpha_file(self.directory, flat_alphas, flat_filling)

    def read_alphas(self):
        '''
        Reads in file_alpharef.txt and file_alpharef_empty.txt from this calculation's directory

        Output:
           alphas -- a list of alpha values (1 per orbital)
        '''

        if not self.do_orbdep or not self.odd_nkscalfact:
            return

        alphas = io.read_alpha_file(self.directory)

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
