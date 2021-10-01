"""

wannier90 calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os
import numpy as np
from ase import Atoms
from ase.calculators.wannier90 import Wannier90
from ase.dft.kpoints import BandPath
from koopmans.utils import warn
from koopmans.settings import Wannier90SettingsDict
from koopmans.commands import Command
from ._utils import ExtendedCalculator, qe_bin_directory


class Wannier90Calculator(ExtendedCalculator, Wannier90):

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the list of parameters
        self.parameters = Wannier90SettingsDict()

        # Define the appropriate file extensions
        self.ext_in = '.win'
        self.ext_out = '.wout'

        # Remove kpts from the kwargs
        mp_grid = kwargs.pop('kpts', [1, 1, 1])

        # Initialise first using the ASE parent and then ExtendedCalculator
        Wannier90.__init__(self, atoms=atoms)
        ExtendedCalculator.__init__(self, *args, **kwargs)

        # Convert the kpts data into the format expected by W90
        self.parameters.mp_grid = np.array(mp_grid)
        kpts = np.indices(self.parameters.mp_grid).transpose(1, 2, 3, 0).reshape(-1, 3)
        kpts = kpts / self.parameters.mp_grid
        kpts[kpts >= 0.5] -= 1
        kpts = BandPath(atoms.cell, kpts)
        self.parameters.kpts = kpts.kpts[:, :3]

        # Set up the command for running this calculator
        self.command = Command(os.environ.get('ASE_WANNIER90_COMMAND',
                                              str(qe_bin_directory) + os.path.sep + self.command))

    def calculate(self):
        super().calculate()

        # Afterwards, check the real vs imaginary component
        if self.parameters.wannier_plot and '-pp' not in self.command.flags:
            max_imre = np.max(self.results['Im/Re ratio'])
            if max_imre > 1e-6:
                warn(f'Im/Re ratio of {max_imre} detected during Wannierisation')

    def is_converged(self):
        return self.results['convergence']

    def is_complete(self):
        return self.results['job done']
