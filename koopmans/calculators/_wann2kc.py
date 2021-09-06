"""

wann2kc calculator module for koopmans

Written by Edward Linscott Feb 2021

"""

from ase.calculators.espresso import Wann2KC
from ase.io.espresso import wann2kc as w2k_io
from ._utils import KCWannCalculator, qe_bin_directory
from koopmans.commands import ParallelCommandWithPostfix


class Wann2KCCalculator(KCWannCalculator):
    # Subclass of KCWannCalculator for performing calculations with wann_to_kc.x

    # Point to the appropriate ASE IO module
    _io = w2k_io
    _ase_calc_class = Wann2KC

    # Define the appropriate file extensions
    ext_in = '.w2ki'
    ext_out = '.w2ko'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.calc.command = ParallelCommandWithPostfix(
            f'{qe_bin_directory}wann2kc.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out}')
