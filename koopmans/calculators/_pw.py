"""

pw calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os
from ase import Atoms
from ase.calculators.espresso import Espresso
from koopmans.settings import PWSettingsDict
from ._utils import ExtendedCalculator, qe_bin_directory
from koopmans.commands import ParallelCommandWithPostfix, Command


class PWCalculator(ExtendedCalculator, Espresso):
    # Subclass of ExtendedCalculator for performing calculations with pw.x

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = PWSettingsDict()

        # Define the appropriate file extensions
        self.ext_in = '.pwi'
        self.ext_out = '.pwo'

        # Initialise first using the ASE parent and then ExtendedCalculator
        Espresso.__init__(self, atoms=atoms)
        ExtendedCalculator.__init__(self, *args, **kwargs)

        self.results_for_qc = ['energy']
        if not isinstance(self.command, Command):
            self.command = ParallelCommandWithPostfix(os.environ.get(
                'ASE_ESPRESSO_COMMAND', str(qe_bin_directory) + os.path.sep + self.command))

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
