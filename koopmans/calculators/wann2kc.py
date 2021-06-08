"""

wann2kc calculator module for python_KI

Written by Edward Linscott Feb 2021

"""

import os
import numpy as np
from ase.calculators.espresso import Wann2KC
from ase.io.espresso import wann2kc as w2k_io
from koopmans import io, utils
from koopmans.calculators.generic import KCWannCalc
from koopmans.calculators.commands import ParallelCommandWithPostfix


class Wann2KCCalc(KCWannCalc):
    # Subclass of KCWannCalc for performing calculations with wann_to_kc.x

    # Point to the appropriate ASE IO module
    _io = w2k_io
    _ase_calc_class = Wann2KC

    # Define the appropriate file extensions
    ext_in = '.w2ki'
    ext_out = '.w2ko'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.calc.command = ParallelCommandWithPostfix(f'wann2kc.x -in PREFIX{self.ext_in} > PREFIX{self.ext_out}')
