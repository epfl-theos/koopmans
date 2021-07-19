"""

kc_screen calculator module for python_KI

Written by Edward Linscott Feb 2021

"""

import numpy as np
from ase.io.espresso import koopmans_screen as kcs_io
from ase.calculators.espresso import KoopmansScreen
from koopmans import io, utils
from koopmans.calculators.generic import KCWannCalc, qe_bin_directory
from koopmans.calculators.commands import ParallelCommandWithPostfix


class KoopmansScreenCalc(KCWannCalc):
    # Subclass of KCWannCalc for performing calculations with kc_screen.x

    # Point to the appropriate ASE IO module
    _io = kcs_io
    _ase_calc_class = KoopmansScreen

    # Define the appropriate file extensions
    ext_in = '.ksi'
    ext_out = '.kso'

    def __init__(self, *args, **kwargs):
        self.defaults.update({'tr2_ph': 1.0e-18, 'nmix_ph': 4, 'niter_ph': 33, 'lrpa': False, 'check_spread': True})
        super().__init__(*args, **kwargs)
        self.calc.command = ParallelCommandWithPostfix(
            f'{qe_bin_directory}kc_screen.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out}')

    def calculate(self):
        # Populate mp1-3
        kpoints = self.calc.parameters.get('kpts', [1, 1, 1])
        [self.mp1, self.mp2, self.mp3] = kpoints

        # Check eps infinity
        if np.max(kpoints) > 1 and self.eps_inf is None:
            utils.warn('You have not specified a value for eps_inf. This will mean that the screening parameters will '
                       'converge very slowly with respect to the k- and q-point grids')

        super().calculate()
