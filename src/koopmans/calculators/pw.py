"""pw calculator module for koopmans."""

import numpy as np
from ase_koopmans import Atoms
from ase_koopmans.calculators.espresso import Espresso
from ase_koopmans.dft.kpoints import BandPath
from ase_koopmans.spectrum.band_structure import BandStructure

from koopmans.cell import cell_follows_qe_conventions, cell_to_parameters
from koopmans.parameters.pw import PWInputParameters

from .calculator import CalculatorInputs, CalculatorOutputs, CalculatorProcess


class PWInputs(CalculatorInputs):
    parameters: PWInputParameters
    pseudopotentials: dict[str, str]
    pseudo_family: str


class PWOutputs(CalculatorOutputs):
    total_energy: float | None = None
    band_structure: BandStructure | None = None


class PWCalculator(CalculatorProcess[PWInputs, PWOutputs]):
    """Calculator process for performing calculations with pw.x."""

    ext_in = '.pwi'
    ext_out = '.pwo'
    code = "pw"
    id = "pw"
    input_model = PWInputs
    output_model = PWOutputs

    def _pre_calculate(self):
        # Update ibrav and celldms
        if cell_follows_qe_conventions(self.inputs.atoms.cell):
            self.inputs.parameters.update(**cell_to_parameters(self.inputs.atoms.cell))
        else:
            self.inputs.parameters.ibrav = 0

        # Make sure kpts has been correctly provided
        if self.inputs.parameters.calculation == 'bands':
            if not isinstance(self.inputs.parameters.kpts, BandPath):
                raise KeyError('You are running a calculation that requires a kpoint path; please provide a `BandPath` '
                               'as the `kpts` parameter')

        super()._pre_calculate()

        return

    def _post_calculate(self):

        super()._post_calculate()

        if isinstance(self.inputs.parameters.kpts, BandPath):
            # Add the bandstructure to the results. This is very un-ASE-y and might eventually be replaced
            self.generate_band_structure()

        return

    def is_converged(self):
        """Return True if the calculation is converged."""
        if self.inputs.parameters.control.calculation == 'scf':
            return self.outputs.total_energy is not None
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
