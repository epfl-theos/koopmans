"""

wannier90 calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os
import numpy as np
from ase import Atoms
from ase.calculators.wannier90 import Wannier90
from ase.dft.kpoints import BandPath
from koopmans.utils import warn
from koopmans.settings import Wannier90SettingsDict
from koopmans.commands import Command
from ._utils import CalculatorExt, CalculatorABC, qe_bin_directory


class Wannier90Calculator(CalculatorExt, Wannier90, CalculatorABC):
    ext_in = '.win'
    ext_out = '.wout'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the list of parameters
        self.parameters = Wannier90SettingsDict()

        # Initialise first using the ASE parent and then CalculatorExt
        Wannier90.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        # Set up the command for running this calculator
        self.command = Command(os.environ.get('ASE_WANNIER90_COMMAND',
                                              str(qe_bin_directory) + os.path.sep + self.command))

        self.results_for_qc = ['centers', 'spreads']

    def is_converged(self):
        return self.results['convergence']

    def is_complete(self):
        return self.results['job done']
