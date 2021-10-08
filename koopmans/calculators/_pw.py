"""

pw calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os
from ase import Atoms
from ase.calculators.espresso import Espresso
from koopmans.settings import PWSettingsDict
from ._utils import CalculatorExt, CalculatorABC, qe_bin_directory
from koopmans.commands import ParallelCommandWithPostfix, Command


class PWCalculator(CalculatorExt, Espresso, CalculatorABC):
    # Subclass of CalculatorExt for performing calculations with pw.x
    ext_in = '.pwi'
    ext_out = '.pwo'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = PWSettingsDict()

        # Initialise first using the ASE parent and then CalculatorExt
        Espresso.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        self.results_for_qc = ['energy']
        if not isinstance(self.command, Command):
            self.command = ParallelCommandWithPostfix(os.environ.get(
                'ASE_ESPRESSO_COMMAND', str(qe_bin_directory) + os.path.sep + self.command))

        self.results_for_qc = ['energy']

    def calculate(self):
        if self.parameters.calculation == 'bands':
            if 'kpath' in self.parameters:
                self.parameters.kpts = self.parameters.pop('kpath')
            else:
                raise KeyError('You are running a calculation that requires a kpath; please provide it in the input '
                               'file')
        else:
            self.parameters.pop('kpath', None)
        super().calculate()

    def is_complete(self):
        return self.results.get('job done', False)

    def is_converged(self):
        return self.results.get('energy', None) is not None
