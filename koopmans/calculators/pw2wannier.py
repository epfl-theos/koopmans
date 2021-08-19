"""

pw2wannier calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os
from koopmans.utils import warn
from ase.io.espresso import pw2wannier as p2w_io
from ase.calculators.espresso import PW2Wannier
from koopmans.calculators.generic import GenericCalc, qe_bin_directory
from koopmans.calculators.commands import ParallelCommand


class PW2Wannier_calc(GenericCalc):
    # Link to relevant ase io module
    _io = p2w_io

    # Define the appropriate file extensions
    ext_in = '.p2wi'
    ext_out = '.p2wo'

    # Adding all wannier90 keywords as decorated properties of the W90_calc class.
    # This means one can set and get wannier90 keywords as self.<keyword> but
    # internally they are stored as self._settings['keyword'] rather than
    # self.<keyword>
    _valid_settings = ['outdir', 'prefix', 'seedname', 'write_mmn', 'write_amn', 'write_uHu',
                       'uHu_formatted', 'write_unk', 'reduce_unk', 'wan_mode', 'wannier_plot',
                       'wannier_plot_list', 'split_evc_file', 'gamma_trick', 'spin_component']

    _settings_that_are_paths = ['outdir']

    def __init__(self, *args, **kwargs):
        self._ase_calc_class = PW2Wannier
        self.settings_to_not_parse = ['wan_mode']
        super().__init__(*args, **kwargs)
        self.calc.command = ParallelCommand(os.environ.get(
            'ASE_PW2WANNIER_COMMAND', qe_bin_directory + self.calc.command))

    def _update_settings_dict(self):
        if 'inputpp' not in self._ase_calc.parameters:
            self._ase_calc.parameters['inputpp'] = {}
        self._settings = self._ase_calc.parameters['inputpp']

    def is_converged(self):
        return True

    def is_complete(self):
        return self.results['job done']

    defaults = {'outdir': 'TMP',
                'prefix': 'kc',
                'seedname': 'wann',
                'wan_mode': 'standalone',
                'write_unk': True}
