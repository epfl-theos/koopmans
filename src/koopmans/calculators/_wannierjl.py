"""

wannierjl calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os

from ase_koopmans import Atoms
from ase_koopmans.calculators.wannierjl import WannierJL

from koopmans.commands import Command
from koopmans.settings import WannierJLSettingsDict
from koopmans.utils import CalculatorNotConvergedWarning, warn

from ._calculator import CalculatorABC, CalculatorExt


class WannierJLCalculator(CalculatorExt, WannierJL, CalculatorABC):
    ext_in = '.wjli'  # this a dummy value -- wjl does not take an input file
    ext_out = '.wjlo'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the list of parameters
        self.parameters = WannierJLSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        WannierJL.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        # Set up the command for running this calculator
        self.command = Command(os.environ.get('ASE_WANNIERJL_COMMAND', self.command))

    def is_converged(self):
        return True

    def is_complete(self):
        return self.results.get('job_done', False)

    def check_convergence(self) -> None:
        return
