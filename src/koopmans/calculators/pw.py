"""pw calculator module for koopmans."""

from typing import Annotated, Any

import numpy as np
import numpy.typing as npt
from ase.dft.kpoints import BandPath
from ase.spectrum.band_structure import BandStructure
from pydantic import BeforeValidator, field_validator, model_validator

from koopmans.cell import cell_follows_qe_conventions, cell_to_parameters
from koopmans.files import File, FileViaCopy
from koopmans.kpoints import GammaOnlyKpoints, Kpoints
from koopmans.parameters.pw import PWInputParameters
from koopmans.serialization import (deserialize_ase_bandstructure,
                                    deserialize_numpy)

from .calculator import CalculatorInputs, CalculatorOutputs, CalculatorProcess


class PWInputs(CalculatorInputs):
    """Inputs for calculations with pw.x."""

    parameters: PWInputParameters
    pseudopotentials: dict[str, str]
    kpoints: Kpoints
    pseudo_family: str
    outdir: FileViaCopy | None = None

    @model_validator(mode='after')
    def check_outdir_provided_if_required(self):
        """Check if outdir has been provided if using `restart_mode = "restart"`."""
        if self.parameters.control.restart_mode == 'restart' and self.outdir is None:
            raise ValueError('Please provide an `outdir` when using `restart_mode = "restart"`')
        return self


class PWOutputs(CalculatorOutputs):
    """Outputs of calculations with pw.x."""

    total_energy: float | None = None
    eigenvalues: Annotated[npt.NDArray[np.float64], BeforeValidator(deserialize_numpy)] | None = None
    fermi_level: list[float]
    band_structure: Annotated[BandStructure, BeforeValidator(deserialize_ase_bandstructure)] | None = None
    outdir: File


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

    def is_converged(self):
        """Return True if the calculation is converged."""
        raise NotImplementedError('Need to revisit')
        if self.inputs.parameters.control.calculation == 'scf':
            return self.outputs.total_energy is not None
        else:
            return True

    def check_convergence(self) -> None:
        """Check if the calculation is converged."""
        raise NotImplementedError('Need to revisit')
        if self.parameters.calculation == 'scf':
            return super().check_convergence()

    def vbm_energy(self) -> float:
        """Calculate the energy of the valence band maximum."""
        raise NotImplementedError('Need to revisit')
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
