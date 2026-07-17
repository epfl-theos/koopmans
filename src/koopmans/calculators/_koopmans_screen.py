"""KCWScreen calculator module for koopmans."""

import numpy as np
from ase_koopmans import Atoms
from ase_koopmans.calculators.espresso import KoopmansScreen

from koopmans import settings
from koopmans.utils.warnings import warn

from ._calculator import CalculatorABC, KCWannCalculator


class KoopmansScreenCalculator(KCWannCalculator, KoopmansScreen, CalculatorABC):
    """Subclass of KCWannCalculator for calculating screening parameters with kcw.x."""

    ext_in = '.ksi'
    ext_out = '.kso'
    code = "kcw_screen"

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = settings.KoopmansScreenSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        KoopmansScreen.__init__(self, atoms=atoms)
        KCWannCalculator.__init__(self, *args, **kwargs)
        super().__init__(*args, **kwargs)

    def _pre_calculate(self):
        # Check eps infinity
        kpoints = [self.parameters.mp1, self.parameters.mp2, self.parameters.mp3]
        if np.max(kpoints) > 1 and self.parameters.eps_inf is None:
            warn('You have not specified a value for `eps_inf`. This will mean that the screening parameters '
                 'will converge very slowly with respect to the k- and q-point grids')

        super()._pre_calculate()

    def is_converged(self):
        """Check if the calculation has converged."""
        raise NotImplementedError('TODO')

    def check_convergence(self) -> None:
        """Check if the calculation has converged. Has not been implemented yet for this calculator."""
        return
