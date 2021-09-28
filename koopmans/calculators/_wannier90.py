"""

wannier90 calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os
from typing import List
from ase.calculators.calculator import CalculationFailed
import numpy as np
from koopmans.utils import warn
from ase.io import wannier90 as w90_io
from ase.calculators.wannier90 import Wannier90
from ase.dft.kpoints import BandPath
from ._utils import ExtendedCalculator, qe_bin_directory
from koopmans.settings import SettingsDict
from koopmans.commands import Command


class Wannier90Calculator(ExtendedCalculator):

    def __init__(self, *args, **kwargs):
        # Define the list of parameters
        self.parameters = SettingsDict(valid=['num_bands', 'num_wann', 'exclude_bands',
                                              'num_iter', 'conv_window', 'conv_tol', 'num_print_cycles',
                                              'dis_froz_max', 'dis_num_iter', 'dis_win_max', 'guiding_centres',
                                              'bands_plot', 'mp_grid', 'kpoint_path', 'projections', 'write_hr',
                                              'write_u_matrices', 'write_xyz', 'wannier_plot', 'gamma_only'],
                                       defaults={'num_iter': 10000, 'conv_tol': 1.e-10, 'conv_window': 5,
                                                 'write_hr': True, 'guiding_centres': True, 'gamma_only': False},
                                       to_not_parse=['exclude_bands'])

        # Define the appropraite file extensions
        self.ext_in = '.win'
        self.ext_out = '.wout'

        # Link the relevant ASE Calculator
        self._ase_calc_class = Wannier90

        super().__init__(*args, **kwargs)
        self.command = Command(os.environ.get('ASE_WANNIER90_COMMAND', str(qe_bin_directory) + self.command))

    def calculate(self):
        if self.parameters.mp_grid is None:
            self.generate_kpoints()
        super().calculate()

        # Afterwards, check the real vs imaginary component
        if self.parameters.wannier_plot and '-pp' not in self.command.flags:
            max_imre = np.max(self.results['Im/Re ratio'])
            if max_imre > 1e-6:
                warn(f'Im/Re ratio of {max_imre} detected during Wannierisation')

    def generate_kpoints(self):
        self.parameters.mp_grid = self.parameters.kpts
        kpts = np.indices(self.parameters.kpts).transpose(1, 2, 3, 0).reshape(-1, 3)
        kpts = kpts / self.parameters.kpts
        kpts[kpts >= 0.5] -= 1
        kpts = BandPath(self.atoms.cell, kpts)
        self.parameters.kpts = kpts.kpts[:, :3]

    def is_converged(self):
        return self.results['convergence']

    def is_complete(self):
        return self.results['job done']
