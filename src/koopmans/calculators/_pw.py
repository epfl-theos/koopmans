"""

pw calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os

import numpy as np
from ase import Atoms
from ase.calculators.espresso import Espresso
from ase.dft.kpoints import BandPath

from koopmans.cell import cell_follows_qe_conventions, cell_to_parameters
from koopmans.commands import Command, ParallelCommandWithPostfix
from koopmans.pseudopotentials import nelec_from_pseudos
from koopmans.settings import PWSettingsDict

from ._utils import CalculatorABC, CalculatorExt, ReturnsBandStructure


class PWCalculator(CalculatorExt, Espresso, ReturnsBandStructure, CalculatorABC):
    # Subclass of CalculatorExt for performing calculations with pw.x
    ext_in = '.pwi'
    ext_out = '.pwo'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = PWSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        Espresso.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        if not isinstance(self.command, Command):
            self.command = ParallelCommandWithPostfix(os.environ.get(
                'ASE_ESPRESSO_COMMAND', self.command))

    def calculate(self):
        # Update ibrav and celldms
        if cell_follows_qe_conventions(self.atoms.cell):
            self.parameters.update(**cell_to_parameters(self.atoms.cell))
        else:
            self.parameters.ibrav = 0
        super().calculate()

    def _calculate(self):
        if self.parameters.calculation == 'bands':
            if not isinstance(self.parameters.kpts, BandPath):
                raise KeyError('You are running a calculation that requires a kpoint path; please provide a BandPath '
                               'as the kpts parameter')

        super()._calculate()

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
        # Fetch the eigenvalues
        eigenvals = self.eigenvalues_from_results()

        # Fetch the total number of electrons in the system
        nelec = nelec_from_pseudos(self.atoms, self.parameters.pseudopotentials,
                                   self.parameters.pseudo_dir) + self.parameters.get('tot_charge', 0)

        # Determine the number of occupied bands in each spin channel
        if self.parameters.nspin == 1:
            n_occs = [nelec // 2]
        else:
            mag = self.parameters.get('tot_magnetization', nelec % 2)
            n_occs = [int(nelec / 2 + mag / 2), int(nelec / 2 - mag / 2)]

        # Return the energy of the highest occupied band
        return np.nanmax([eigs[:, :n_occ] for eigs, n_occ in zip(eigenvals, n_occs)])

    def eigenvalues_from_results(self):
        class_name = self.__class__.__name__
        assert getattr(self, 'kpts', None) is not None, f'Please call {class_name}.calculate() prior to calling ' \
            f'{class_name}.eigenvalues_from_results()'

        i_spins = [i for i in range(2) if i in [k.s for k in self.kpts]]
        return np.array([[k.eps_n for k in self.kpts if k.s == i_spin] for i_spin in i_spins])
