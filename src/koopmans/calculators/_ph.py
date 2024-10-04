

"""

PhCalculator calculator module for koopmans

Written by Marija Stojkovic  May 2022

"""

import os

import numpy as np
from ase import Atoms
from ase.calculators.espresso import EspressoPh

from koopmans.commands import Command, ParallelCommand
from koopmans.settings import PhSettingsDict

from ._calculator import CalculatorABC, CalculatorExt


class PhCalculator(CalculatorExt, EspressoPh, CalculatorABC):

    ext_in = '.phi'
    ext_out = '.pho'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        self.parameters = PhSettingsDict()
        self.parent = None

        # Initialise first using the ASE parent and then CalculatorExt
        EspressoPh.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        self.command = ParallelCommand(f'ph.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out} 2>&1')

    def is_converged(self):
        return True

    def is_complete(self):
        return self.results['job done']

    def _post_calculate(self):
        super()._post_calculate()
        if self.parameters.trans:
            self.read_dynG()

    def read_dynG(self):
        with open(self.parameters.fildyn, 'r') as fd:
            flines = fd.readlines()

        i = [x.strip() for x in flines].index('Dielectric Tensor:')
        k = [x.strip() for x in flines].index('Effective Charges E-U: Z_{alpha}{s,beta}')
        epsilon = np.array([x.split() for x in flines[i + 2: k - 1]], dtype=float)
        self.results['dielectric tensor'] = epsilon
