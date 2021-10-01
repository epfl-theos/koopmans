"""

wann2kc calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

import os
from ase import Atoms
from ase.calculators.espresso import Wann2KC
from ._utils import KCWannCalculator, qe_bin_directory
from koopmans.settings import Wann2KCSettingsDict
from koopmans.commands import ParallelCommandWithPostfix


class Wann2KCCalculator(KCWannCalculator, Wann2KC):
    # Subclass of KCWannCalculator for performing calculations with wann_to_kc.x

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = Wann2KCSettingsDict()

        # Define the appropriate file extensions
        self.ext_in = '.w2ki'
        self.ext_out = '.w2ko'

        # Initialise first using the ASE parent and then ExtendedCalculator
        Wann2KC.__init__(self, atoms=atoms)
        KCWannCalculator.__init__(self, *args, **kwargs)

        self.command = ParallelCommandWithPostfix(
            f'{qe_bin_directory}{os.path.sep}wann2kc.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out}')
