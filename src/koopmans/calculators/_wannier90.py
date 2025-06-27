"""wannier90 calculator module for koopmans."""

from ase_koopmans import Atoms
from ase_koopmans.calculators.wannier90 import Wannier90

from koopmans.settings import Wannier90SettingsDict
from koopmans.utils import CalculatorNotConvergedWarning, warn

from ._calculator import CalculatorABC, CalculatorExt


class Wannier90Calculator(CalculatorExt, Wannier90, CalculatorABC):
    """A `Wannier90` calculator."""

    ext_in = '.win'
    ext_out = '.wout'
    code = "wannier90"

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the list of parameters
        self.parameters = Wannier90SettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        Wannier90.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

    def is_converged(self):
        """Return True if the calculation is converged."""
        return self.results['convergence']

    def is_complete(self):
        """Return True if the calculation is complete."""
        return self.results['job done']

    def check_convergence(self) -> None:
        """Check the convergence of the Wannier90 calculation.

        For projwfs (num_iter=0) and preproc calculations the convergence check
        cannot be applied; for mlwfs a warning is printed out in case the calculation
        is not converged (but we allow the calculation to proceed)
        """
        if self.parameters.num_iter == 0 or ' -pp ' in self.command:
            pass
        elif not self.is_converged():
            warn(f'`{self.directory}/{self.prefix}` did not converge; proceed with caution',
                 CalculatorNotConvergedWarning)
