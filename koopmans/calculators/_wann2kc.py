"""

wann2kc calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

from ase.calculators.espresso import Wann2KC
from ase.io.espresso import wann2kc as w2k_io
from ._utils import KCWannCalculator, qe_bin_directory, kc_wann_defaults
from koopmans.settings import SettingsDict
from koopmans.commands import ParallelCommandWithPostfix


class Wann2KCCalculator(KCWannCalculator):
    # Subclass of KCWannCalculator for performing calculations with wann_to_kc.x

    def __init__(self, *args, **kwargs):
        # Link to corresponding ASE Calculator
        self._ase_calc_class = Wann2KC

        # Define the valid settings
        self.parameters = SettingsDict(valid=[k for block in w2k_io.KEYS for k in block],
                                       defaults=kc_wann_defaults,
                                       are_paths=['outdir', 'pseudo_dir'],
                                       to_not_parse=['assume_isolated'])

        # Define the appropriate file extensions
        self.ext_in = '.w2ki'
        self.ext_out = '.w2ko'

        super().__init__(*args, **kwargs)
        self.command = ParallelCommandWithPostfix(
            f'{qe_bin_directory}wann2kc.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out}')
