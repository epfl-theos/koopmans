"""Wann2KCW calculator module for koopmans."""

from ase_koopmans import Atoms
from ase_koopmans.calculators.espresso import Wann2KC

from koopmans.settings import Wann2KCSettingsDict

from ._calculator import CalculatorABC, KCWannCalculator


class Wann2KCCalculator(KCWannCalculator, Wann2KC, CalculatorABC):
    """Subclass of KCWannCalculator for converting Wannier functions to a KCW format with kcw.x."""

    ext_in = '.w2ki'
    ext_out = '.w2ko'
    code = "kcw_wannier"

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = Wann2KCSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        Wann2KC.__init__(self, atoms=atoms)
        KCWannCalculator.__init__(self, *args, **kwargs)

    def is_converged(self):
        """Return True; a Wann2KC calculation is never not "converged"."""
        return True

    def is_complete(self):
        """Return True if the calculation is complete."""
        return self.results['job_done']
