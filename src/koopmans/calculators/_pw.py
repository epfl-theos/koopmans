"""pw calculator module for koopmans."""

import numpy as np
from ase_koopmans import Atoms
from ase_koopmans.calculators.espresso import Espresso
from ase_koopmans.dft.kpoints import BandPath

from koopmans.cell import cell_follows_qe_conventions, cell_to_parameters
from koopmans.settings import PWSettingsDict

from ._calculator import CalculatorABC, CalculatorExt, ReturnsBandStructure


class PWCalculator(CalculatorExt, Espresso, ReturnsBandStructure, CalculatorABC):
    """Subclass of CalculatorExt for performing calculations with pw.x."""

    ext_in = '.pwi'
    ext_out = '.pwo'
    code = "pw"

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = PWSettingsDict()

        # Initialize first using the ASE parent and then CalculatorExt
        Espresso.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

    def _pre_calculate(self):
        # Update ibrav and celldms
        if cell_follows_qe_conventions(self.atoms.cell):
            self.parameters.update(**cell_to_parameters(self.atoms.cell))
        else:
            self.parameters.ibrav = 0

        # Make sure kpts has been correctly provided
        if self.parameters.calculation == 'bands':
            if not isinstance(self.parameters.kpts, BandPath):
                raise KeyError('You are running a calculation that requires a kpoint path; please provide a `BandPath` '
                               'as the `kpts` parameter')

        super()._pre_calculate()

        return

    def _post_calculate(self):

        super()._post_calculate()

        if isinstance(self.parameters.kpts, BandPath):
            # Add the bandstructure to the results. This is very un-ASE-y and might eventually be replaced
            self.generate_band_structure()

        return

    def is_complete(self):
        """Return True if the calculation is complete."""
        return self.results.get('job done', False)

    def is_converged(self):
        """Return True if the calculation is converged."""
        if self.parameters.calculation == 'scf':
            return self.results.get('energy', None) is not None
        else:
            return True

    def check_convergence(self) -> None:
        """Check if the calculation is converged."""
        if self.parameters.calculation == 'scf':
            return super().check_convergence()

    def vbm_energy(self) -> float:
        """Calculate the energy of the valence band maximum."""
        # Fetch the eigenvalues
        eigenvals = self.eigenvalues_from_results()

        # Fetch the total number of electrons in the system
        nelec = self.results['nelec']

        # Determine the number of occupied bands in each spin channel
        if self.parameters.nspin == 1:
            n_occs = [int(nelec // 2)]
        else:
            mag = self.parameters.get('tot_magnetization', nelec % 2)
            n_occs = [int(nelec / 2 + mag / 2), int(nelec / 2 - mag / 2)]

        # Return the energy of the highest occupied band
        return np.max([np.nanmax(eigs[:, :n_occ]) for eigs, n_occ in zip(eigenvals, n_occs)])

    def eigenvalues_from_results(self):
        """Extract the eigenvalues from the results dictionary."""
        class_name = self.__class__.__name__
        assert getattr(self, 'kpts', None) is not None, f'Please call {class_name}.calculate() prior to calling ' \
            f'{class_name}.eigenvalues_from_results()'

        i_spins = [i for i in range(2) if i in [k.s for k in self.kpts]]
        return np.array([[k.eps_n for k in self.kpts if k.s == i_spin] for i_spin in i_spins])
