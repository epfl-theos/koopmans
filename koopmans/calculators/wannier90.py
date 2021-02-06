"""

wannier90 calculator module for python_KI

Written by Edward Linscott Sep 2020

"""

import numpy as np
from koopmans.utils import warn
from ase.io import wannier90 as w90_io
from ase.calculators.wannier90 import Wannier90
from ase.dft.kpoints import BandPath
from koopmans.calculators.generic import QE_calc


class W90_calc(QE_calc):
    # Link to relevant ase io module
    _io = w90_io

    # Define the appropriate file extensions
    ext_in = '.win'
    ext_out = '.wout'

    # Adding all wannier90 keywords as decorated properties of the W90_calc class.
    # This means one can set and get wannier90 keywords as self.<keyword> but
    # internally they are stored as self._settings['keyword'] rather than
    # self.<keyword>
    _recognised_keywords = ['num_bands', 'num_wann', 'exclude_bands',
                            'num_iter', 'conv_window', 'conv_tol', 'num_print_cycles',
                            'dis_froz_max', 'dis_num_iter', 'dis_win_max', 'guiding_centres',
                            'bands_plot', 'mp_grid', 'kpoint_path', 'projections', 'write_hr']

    for k in _recognised_keywords:
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

    def __init__(self, *args, **kwargs):
        self._ase_calc_class = Wannier90
        self.settings_to_not_parse = ['exclude_bands']
        super().__init__(*args, **kwargs)

    def calculate(self):
        if self.mp_grid is None:
            self.generate_kpoints()
        super().calculate()

    def generate_kpoints(self):
        self.mp_grid = self.calc.parameters['kpts']
        kpts = np.indices(self.calc.parameters['kpts']).transpose(1, 2, 3, 0).reshape(-1, 3)
        kpts = kpts / self.calc.parameters['kpts']
        kpts[kpts >= 0.5] -= 1
        kpts = BandPath(self.calc.atoms.cell, kpts)
        self.calc.parameters['kpts'] = kpts.kpts[:, :3]

    def is_converged(self):
        warn("is_converged is not properly implemented")
        return True

    def is_complete(self):
        return self.results['job done']
