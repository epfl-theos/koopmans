"""

pw calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os
import numpy as np
from ase import Atoms
from ase.calculators.espresso import Espresso
from ase.dft.kpoints import BandPath
from koopmans.settings import PWSettingsDict
from ._utils import CalculatorExt, CalculatorABC, bin_directory, ReturnsBandStructure
from koopmans.commands import ParallelCommandWithPostfix, Command


class PWCalculator(CalculatorExt, Espresso, ReturnsBandStructure, CalculatorABC):
    # Subclass of CalculatorExt for performing calculations with pw.x
    ext_in = '.pwi'
    ext_out = '.pwo'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = PWSettingsDict()

        # Initialise first using the ASE parent and then CalculatorExt
        Espresso.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        self.results_for_qc = ['energy']
        if not isinstance(self.command, Command):
            self.command = ParallelCommandWithPostfix(os.environ.get(
                'ASE_ESPRESSO_COMMAND', str(bin_directory) + os.path.sep + self.command))

        self.results_for_qc = ['energy']

    def calculate(self):
        if self.parameters.calculation == 'bands':
            if not isinstance(self.parameters.kpts, BandPath):
                raise KeyError('You are running a calculation that requires a kpoint path; please provide a BandPath '
                               'as the kpts parameter')
        super().calculate()

        if isinstance(self.parameters.kpts, BandPath):
            # Add the bandstructure to the results. This is very un-ASE-y and might eventually be replaced
            self.generate_band_structure()

    def is_complete(self):
        return self.results.get('job done', False)

    def is_converged(self):
        if self.parameters.calculation == 'scf':
            return self.results.get('energy', None) is not None
        else:
            return True

    def check_convergence(self) -> None:
        if self.parameters.calculation == 'scf':
            return super().check_convergence()

    def vbm_energy(self) -> float:
        return 0.0

    def eigenvalues_from_results(self):
        class_name = self.__class__.__name__
        assert getattr(self, 'kpts', None) is not None, f'Please call {class_name}.calculate() prior to calling ' \
            f'{class_name}.eigenvalues_from_results()'

        i_spins = [i for i in range(2) if i in [k.s for k in self.kpts]]
        return np.array([[k.eps_n for k in self.kpts if k.s == i_spin] for i_spin in i_spins])
