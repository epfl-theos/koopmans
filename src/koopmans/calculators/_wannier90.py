"""

wannier90 calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os

import numpy as np
from ase import Atoms
from ase.calculators.wannier90 import Wannier90
from ase.dft.kpoints import BandPath

from koopmans.commands import Command
from koopmans.settings import Wannier90SettingsDict
from koopmans.utils import CalculatorNotConvergedWarning, warn

from ._utils import CalculatorABC, CalculatorExt


class Wannier90Calculator(CalculatorExt, Wannier90, CalculatorABC):
    ext_in = '.win'
    ext_out = '.wout'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the list of parameters
        self.parameters = Wannier90SettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        Wannier90.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        # Set up the command for running this calculator
        self.command = Command(os.environ.get('ASE_WANNIER90_COMMAND', self.command))

    def is_converged(self):
        return self.results['convergence']

    def is_complete(self):
        return self.results['job done']

    def check_convergence(self) -> None:
        # For projwfs (num_iter=0) and preproc calculations the convergence check
        # cannot be applied; for mlwfs a warning is printed out in case the calculation
        # is not converged (but we allow the calculation to proceed)
        if self.parameters.num_iter == 0 or 'preproc' in self.prefix:
            pass
        elif not self.is_converged():
            warn(f'{self.directory}/{self.prefix} did not converge; proceed with caution',
                 CalculatorNotConvergedWarning)
