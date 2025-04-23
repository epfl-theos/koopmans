"""pw2wannier calculator module for koopmans."""

import os

from ase_koopmans import Atoms
from ase_koopmans.calculators.espresso import PW2Wannier

from koopmans.commands import ParallelCommand
from koopmans.settings import PW2WannierSettingsDict

from ._calculator import CalculatorABC, CalculatorExt


class PW2WannierCalculator(CalculatorExt, PW2Wannier, CalculatorABC):
    """A `PW2Wannier` calculator."""

    ext_in = '.p2wi'
    ext_out = '.p2wo'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        self.parameters = PW2WannierSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        PW2Wannier.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        self.command = ParallelCommand(os.environ.get('ASE_PW2WANNIER_COMMAND', self.command))

    def is_converged(self):
        """Return True; a PW2Wannier calculation is never not "converged"."""
        return True

    def is_complete(self):
        """Return True if the calculation is complete."""
        return self.results['job done']
