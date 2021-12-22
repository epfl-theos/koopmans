"""

projwfc.x calculator module for koopmans

"""

import os
from ase import Atoms
from koopmans.commands import Command, ParallelCommand
from koopmans.settings import ProjwfcSettingsDict
from ase.calculators.espresso import Projwfc
from ._utils import CalculatorExt, CalculatorABC, qe_bin_directory


class ProjwfcCalculator(CalculatorExt, Projwfc, CalculatorABC):
    # Subclass of CalculatorExt for performing calculations with projwfc.x
    ext_in = '.pri'
    ext_out = '.pro'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = ProjwfcSettingsDict()

        # Initialise first using the ASE parent and then CalculatorExt
        Projwfc.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        self.results_for_qc = ['DOS']
        if not isinstance(self.command, Command):
            self.command = ParallelCommand(os.environ.get(
                'ASE_PROJWFC_COMMAND', str(qe_bin_directory) + os.path.sep + self.command))

        self.results_for_qc = []
    def is_complete(self):
        return True
    def is_converged(self):
        return True
