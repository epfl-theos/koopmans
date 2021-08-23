"""

wannier90 calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os
from ase.calculators.calculator import CalculationFailed
import numpy as np
from koopmans.utils import warn
from ase.io import wannier90 as w90_io
from ase.calculators.wannier90 import Wannier90
from ase.dft.kpoints import BandPath
from koopmans.calculators.generic import GenericCalc, qe_bin_directory
from koopmans.calculators.commands import Command


class W90_calc(GenericCalc):
    # Link to relevant ase io module
    _io = w90_io

    # Define the appropriate file extensions
    ext_in = '.win'
    ext_out = '.wout'

    # Create a record of the valid settings
    _valid_settings = ['num_bands', 'num_wann', 'exclude_bands',
                       'num_iter', 'conv_window', 'conv_tol', 'num_print_cycles',
                       'dis_froz_max', 'dis_num_iter', 'dis_win_max', 'guiding_centres',
                       'bands_plot', 'mp_grid', 'kpoint_path', 'projections', 'write_hr',
                       'write_u_matrices', 'write_xyz', 'wannier_plot']
    _settings_that_are_paths = []

    def __init__(self, *args, **kwargs):
        self._ase_calc_class = Wannier90
        self.settings_to_not_parse = ['exclude_bands']
        super().__init__(*args, **kwargs)
        self.calc.command = Command(os.environ.get('ASE_WANNIER90_COMMAND', qe_bin_directory + self.calc.command))

    def calculate(self):
        if self.mp_grid is None:
            self.generate_kpoints()
        super().calculate()

        # Afterwards, check the real vs imaginary component
        if self.wannier_plot and '-pp' not in self.calc.command.flags:
            max_imre = np.max(self.results['Im/Re ratio'])
            if max_imre > 1e-6:
                warn(f'Im/Re ratio of {max_imre} detected during Wannierisation')

    def generate_kpoints(self):
        self.mp_grid = self.calc.parameters['kpts']
        kpts = np.indices(self.calc.parameters['kpts']).transpose(1, 2, 3, 0).reshape(-1, 3)
        kpts = kpts / self.calc.parameters['kpts']
        kpts[kpts >= 0.5] -= 1
        kpts = BandPath(self.calc.atoms.cell, kpts)
        self.calc.parameters['kpts'] = kpts.kpts[:, :3]

    def is_converged(self):
        return self.results['convergence']

    def is_complete(self):
        return self.results['job done']

    @property
    def defaults(self):
        return {'num_iter': 10000,
                'conv_tol': 1.e-10,
                'conv_window': 5,
                'write_hr': True,
                'guiding_centres': True,
                'gamma_only': False}
