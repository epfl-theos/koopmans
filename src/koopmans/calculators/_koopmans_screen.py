"""

KCWScreen calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

import os

import numpy as np
from ase import Atoms
from ase.calculators.espresso import KoopmansScreen

from koopmans import settings, utils
from koopmans.commands import ParallelCommandWithPostfix

from ._utils import CalculatorABC, KCWannCalculator


class KoopmansScreenCalculator(KCWannCalculator, KoopmansScreen, CalculatorABC):
    # Subclass of KCWannCalculator for calculating screening parameters with kcw.x
    ext_in = '.ksi'
    ext_out = '.kso'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = settings.KoopmansScreenSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        KoopmansScreen.__init__(self, atoms=atoms)
        KCWannCalculator.__init__(self, *args, **kwargs)
        super().__init__(*args, **kwargs)

        self.command = ParallelCommandWithPostfix(f'kcw.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out} 2>&1')

    def calculate(self):
        # Check eps infinity
        kpoints = [self.parameters.mp1, self.parameters.mp2, self.parameters.mp3]
        if np.max(kpoints) > 1 and self.parameters.eps_inf is None:
            utils.warn('You have not specified a value for eps_inf. This will mean that the screening parameters will '
                       'converge very slowly with respect to the k- and q-point grids')

        super().calculate()

    def is_converged(self):
        raise NotImplementedError('TODO')

    def check_convergence(self) -> None:
        # is_converged has not been implemented yet for this calculator
        return
