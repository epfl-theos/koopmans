"""WannierJL calculator module for koopmans."""

from ase_koopmans import Atoms
from ase_koopmans.calculators.wannierjl import WannierJL

from koopmans.settings import WannierJLSettingsDict

from ._calculator import CalculatorABC, CalculatorExt


class WannierJLCalculator(CalculatorExt, WannierJL, CalculatorABC):
    """Calculator for WannierJL."""

    ext_in = '.wjli'  # this a dummy value -- wjl does not take an input file
    ext_out = '.wjlo'
    code = "wannierjl"

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the list of parameters
        self.parameters = WannierJLSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        WannierJL.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

    def is_converged(self):
        """Check if the calculation has converged."""
        return True

    def is_complete(self):
        """Check if the calculation is complete."""
        return self.results.get('job_done', False)

    def check_convergence(self) -> None:
        """Check the convergence of the calculation."""
        return
