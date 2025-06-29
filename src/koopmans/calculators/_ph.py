"""PhCalculator calculator module for koopmans."""

import numpy as np
from ase_koopmans import Atoms
from ase_koopmans.calculators.espresso import EspressoPh

from koopmans.settings import PhSettingsDict

from ._calculator import CalculatorABC, CalculatorExt


class PhCalculator(CalculatorExt, EspressoPh, CalculatorABC):
    """Subclass of CalculatorExt for performing calculations with ph.x."""

    ext_in = '.phi'
    ext_out = '.pho'
    code = "ph"

    def __init__(self, atoms: Atoms, *args, **kwargs):
        self.parameters = PhSettingsDict()
        self.parent_process = None

        # Initialise first using the ASE parent and then CalculatorExt
        EspressoPh.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

    def is_converged(self):
        """Return True; a ph calculation is never not "converged"."""
        return True

    def is_complete(self):
        """Return True if the calculation is complete."""
        return self.results['job done']

    def _post_calculate(self):
        super()._post_calculate()
        if self.parameters.trans:
            self.read_dynG()

    def read_dynG(self):
        """Read the dynamical matrix from the output file."""
        with open(self.parameters.fildyn, 'r') as fd:
            flines = fd.readlines()

        i = [x.strip() for x in flines].index('Dielectric Tensor:')
        k = [x.strip() for x in flines].index('Effective Charges E-U: Z_{alpha}{s,beta}')
        epsilon = np.array([x.split() for x in flines[i + 2: k - 1]], dtype=float)
        self.results['dielectric tensor'] = epsilon
