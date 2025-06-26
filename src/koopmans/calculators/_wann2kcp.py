"""wann2kcp calculator module for koopmans."""

from ase_koopmans import Atoms
from ase_koopmans.calculators.espresso import Wann2KCP

from koopmans.settings import Wann2KCPSettingsDict

from ._calculator import CalculatorABC, CalculatorExt


class Wann2KCPCalculator(CalculatorExt, Wann2KCP, CalculatorABC):
    """A `Wann2KCP` calculator."""

    ext_in = '.wki'
    ext_out = '.wko'
    code = "wann2kcp"

    def __init__(self, atoms: Atoms, *args, **kwargs):
        self.parameters = Wann2KCPSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        Wann2KCP.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

    def is_converged(self):
        """Return True; a Wann2KCP calculation is never not "converged"."""
        return True

    def is_complete(self):
        """Return True if the calculation is complete."""
        return self.results['job done']
