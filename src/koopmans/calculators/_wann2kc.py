"""

Wann2KCW calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

import os

from ase import Atoms
from ase.calculators.espresso import Wann2KC

from koopmans.commands import ParallelCommandWithPostfix
from koopmans.settings import Wann2KCSettingsDict

from ._utils import CalculatorABC, KCWannCalculator


class Wann2KCCalculator(KCWannCalculator, Wann2KC, CalculatorABC):
    # Subclass of KCWannCalculator for converting Wannier functions to a KCW format with kcw.x
    ext_in = '.w2ki'
    ext_out = '.w2ko'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = Wann2KCSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        Wann2KC.__init__(self, atoms=atoms)
        KCWannCalculator.__init__(self, *args, **kwargs)

        self.command = ParallelCommandWithPostfix(f'kcw.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out} 2>&1')

    def is_converged(self):
        return True

    def is_complete(self):
        return self.results['job_done']
