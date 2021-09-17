"""

pw calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os
from ase.io.espresso import pw as pw_io
from ase.calculators.espresso import Espresso
from koopmans.settings import SettingsDict
from ._utils import ExtendedCalculator, qe_bin_directory
from koopmans.commands import ParallelCommandWithPostfix, Command


class PWCalculator(ExtendedCalculator, Espresso):
    # Subclass of ExtendedCalculator for performing calculations with pw.x

    def __init__(self, *args, **kwargs):
        # Link to corresponding ASE Calculator
        self._ase_calc_class = Espresso

        # Define the valid settings
        self.parameters = SettingsDict(valid=[k for block in pw_io.KEYS for k in block],
                                       defaults={'calculation': 'scf', 'outdir': './TMP/', 'prefix': 'kc',
                                                 'conv_thr': '2.0e-9*nelec'},
                                       are_paths=['outdir', 'pseudo_dir'],
                                       to_not_parse=['assume_isolated'])

        # Define the appropriate file extensions
        self.ext_in = '.pwi'
        self.ext_out = '.pwo'

        super().__init__(*args, **kwargs)

        self.results_for_qc = ['energy']
        if not isinstance(self.command, Command):
            self.command = ParallelCommandWithPostfix(os.environ.get(
                'ASE_ESPRESSO_COMMAND', qe_bin_directory + self.command))

    def calculate(self):
        if self.calculation == 'bands':
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
