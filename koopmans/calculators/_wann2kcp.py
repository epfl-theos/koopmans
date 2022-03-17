"""

wann2kcp calculator module for koopmans

Written by Riccardo De Gennaro Mar 2022

"""

import os
from ase import Atoms
from ase.calculators.espresso import Wann2KCP
from koopmans.settings import Wann2KCPSettingsDict
from koopmans.commands import ParallelCommand
from ._utils import CalculatorExt, CalculatorABC, bin_directory


class Wann2KCPCalculator(CalculatorExt, Wann2KCP, CalculatorABC):

    ext_in = '.wki'
    ext_out = '.wko'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        self.parameters = Wann2KCPSettingsDict()

        # Initialise first using the ASE parent and then CalculatorExt
        Wann2KCP.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        self.command = ParallelCommand(os.environ.get('ASE_WANN2KCP_COMMAND',
                                                      str(bin_directory) + os.path.sep + self.command))

    def is_converged(self):
        return True

    def is_complete(self):
        return self.results['job done']
