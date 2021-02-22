"""

pw calculator module for python_KI

Written by Edward Linscott Sep 2020

"""

from ase.io import espresso as pw_io
from ase.calculators.espresso import Espresso
from koopmans import io
from koopmans.calculators.generic import EspressoCalc


class PW_calc(EspressoCalc):
    # Subclass of EspressoCalc for performing calculations with pw.x

    # Point to the appropriate ASE IO module
    _io = pw_io

    # Define the appropriate file extensions
    ext_in = '.pwi'
    ext_out = '.pwo'

    # Create a list of the valid settings
    _valid_settings = [k for sublist in _io.KEYS.values() for k in sublist]
    _settings_that_are_paths = ['outdir', 'pseudo_dir']

    def __init__(self, *args, **kwargs):
        self._ase_calc_class = Espresso
        self.settings_to_not_parse = ['assume_isolated']
        super().__init__(*args, **kwargs)
        self.results_for_qc = ['energy']

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
                'outdir': './TMP-PW/',
                'prefix': 'kc',
                'conv_thr': '2.0e-9*nelec'}
