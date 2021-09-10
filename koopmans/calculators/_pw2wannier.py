"""

pw2wannier calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

from koopmans.settings import Setting, SettingsDict
import os
from koopmans.utils import warn
from ase.io.espresso import pw2wannier as p2w_io
from ase.calculators.espresso import PW2Wannier
from ._utils import ExtendedCalculator, qe_bin_directory
from koopmans.commands import ParallelCommand


class PW2WannierCalculator(ExtendedCalculator):

    # Adding all wannier90 keywords as decorated properties of the Wannier90Calculator class.
    # This means one can set and get wannier90 keywords as self.<keyword> but
    # internally they are stored as self._settings['keyword'] rather than
    # self.<keyword>

    def __init__(self, *args, **kwargs):
        self.parameters = SettingsDict(valid=['outdir', 'prefix', 'seedname', 'write_mmn',
                                              'write_amn', 'write_uHu', 'uHu_formatted',
                                              'write_unk', 'reduce_unk', 'wan_mode', 'wannier_plot',
                                              'wannier_plot_list', 'split_evc_file', 'gamma_trick',
                                              'spin_component'],
                                       defaults={'outdir': 'TMP', 'prefix': 'kc', 'seedname': 'wann',
                                                 'wan_mode': 'standalone'},
                                       are_paths=['outdir'],
                                       to_not_parse=['wan_mode'])

        # Link to corresponding ASE IO module and Calculator
        self._io = p2w_io
        self._ase_calc_class = PW2Wannier

        # Define the appropriate file extensions
        self.ext_in = '.p2wi'
        self.ext_out = '.p2wo'

        super().__init__(*args, **kwargs)
        self.command = ParallelCommand(os.environ.get('ASE_PW2WANNIER_COMMAND', qe_bin_directory + self.command))

    def is_converged(self):
        return True

    def is_complete(self):
        return self.results['job done']
