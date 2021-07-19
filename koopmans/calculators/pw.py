"""

pw calculator module for python_KI

Written by Edward Linscott Sep 2020

"""

import os
from ase.io.espresso import pw as pw_io
from ase.calculators.espresso import Espresso
from koopmans import io
from koopmans.calculators.generic import EspressoCalc, qe_bin_directory
from koopmans.calculators.commands import ParallelCommandWithPostfix, Command


class PW_calc(EspressoCalc):
    # Subclass of EspressoCalc for performing calculations with pw.x

    # Point to the appropriate ASE IO module
    _io = pw_io

    # Define the appropriate file extensions
    ext_in = '.pwi'
    ext_out = '.pwo'

    _settings_that_are_paths = ['outdir', 'pseudo_dir']

    def __init__(self, *args, **kwargs):
        self._ase_calc_class = Espresso
        self.settings_to_not_parse = ['assume_isolated']
        super().__init__(*args, **kwargs)
        self.results_for_qc = ['energy']
        if not isinstance(self.calc.command, Command):
            self.calc.command = ParallelCommandWithPostfix(os.environ.get(
                'ASE_ESPRESSO_COMMAND', qe_bin_directory + self.calc.command))

    def is_complete(self):
        return self.results.get('job done', False)

    def is_converged(self):
        return self.results.get('energy', None) is not None

    @property
    def nelec(self):
        # PW does not directly have an "nelec" keyword, but it is useful for setting
        # extensive properties such as energy convergence thresholds
        return io.nelec_from_pseudos(self.calc)

    @nelec.setter
    def nelec(self, value):
        raise ValueError('You tried to set PW_calc.nelec, which is invalid')

    defaults = {'calculation': 'scf',
                'outdir': './TMP/',
                'prefix': 'kc',
                'conv_thr': '2.0e-9*nelec'}
