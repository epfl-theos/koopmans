"""

kc_screen calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

import numpy as np
from ase.io.espresso import koopmans_screen as kcs_io
from ase.calculators.espresso import KoopmansScreen
from koopmans import utils, settings
from ._utils import KCWannCalculator, qe_bin_directory, kc_wann_defaults
from koopmans.commands import ParallelCommandWithPostfix


class KoopmansScreenCalculator(KCWannCalculator):
    # Subclass of KCWannCalculator for performing calculations with kc_screen.x

    def __init__(self, *args, **kwargs):
        # Link to corresponding ASE Calculator
        self._ase_calc_class = KoopmansScreen

        # Define the valid settings
        self.parameters = settings.SettingsDict(valid=[k for block in kcs_io.KEYS for k in block],
                                                defaults={'tr2_ph': 1.0e-18, 'nmix_ph': 4, 'niter_ph': 33, 'lrpa': False,
                                                          'check_spread': True, **kc_wann_defaults},
                                                are_paths=['outdir', 'pseudo_dir'],
                                                to_not_parse=['assume_isolated'])

        # Define the appropriate file extensions
        self.ext_in = '.ksi'
        self.ext_out = '.kso'

        super().__init__(*args, **kwargs)
        self.command = ParallelCommandWithPostfix(
            f'{qe_bin_directory}kc_screen.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out}')

    def calculate(self):
        # Populate mp1-3
        kpoints = self.parameters.get('kpts', [1, 1, 1])
        [self.parameters.mp1, self.parameters.mp2, self.parameters.mp3] = kpoints

        # Check eps infinity
        if np.max(kpoints) > 1 and self.parameters.eps_inf is None:
            utils.warn('You have not specified a value for eps_inf. This will mean that the screening parameters will '
                       'converge very slowly with respect to the k- and q-point grids')

        super().calculate()
