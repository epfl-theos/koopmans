"""

kc_screen calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

import os
import numpy as np
from ase import Atoms
from ase.calculators.espresso import KoopmansScreen
from koopmans import utils, settings
from ._utils import KCWannCalculator, CalculatorABC, qe_bin_directory
from koopmans.commands import ParallelCommandWithPostfix


class KoopmansScreenCalculator(KCWannCalculator, KoopmansScreen, CalculatorABC):
    # Subclass of KCWannCalculator for performing calculations with kc_screen.x
    ext_in = '.ksi'
    ext_out = '.kso'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = settings.KoopmansScreenSettingsDict()

        # Initialise first using the ASE parent and then CalculatorExt
        KoopmansScreen.__init__(self, atoms=atoms)
        KCWannCalculator.__init__(self, *args, **kwargs)
        super().__init__(*args, **kwargs)

        self.command = ParallelCommandWithPostfix(
            f'{qe_bin_directory}{os.path.sep}kc_screen.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out}')

    def calculate(self):
        # Populate mp1-3
        kpoints = self.parameters.get('kpts', [1, 1, 1])
        [self.parameters.mp1, self.parameters.mp2, self.parameters.mp3] = kpoints

        # Check eps infinity
        if np.max(kpoints) > 1 and self.parameters.eps_inf is None:
            utils.warn('You have not specified a value for eps_inf. This will mean that the screening parameters will '
                       'converge very slowly with respect to the k- and q-point grids')

        super().calculate()
