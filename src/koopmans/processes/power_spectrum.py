"""
Processes used during the machine learning workflows
"""
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from ase_koopmans.cell import Cell
from pydantic import ConfigDict

from koopmans import ml, utils
from koopmans.bands import Band
from koopmans.files import FilePointer

from ._process import IOModel, Process


class ExtractCoefficientsFromXMLInput(IOModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    n_max: int
    l_max: int
    r_min: float
    r_max: float
    r_cut: float
    wannier_centers: List[List[float]]
    total_density_xml: FilePointer
    cell: Cell
    orbital_densities_xml: List[FilePointer]
    bands: List[Band]


class ExtractCoefficientsFromXMLOutput(IOModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    precomputed_alphas: FilePointer
    precomputed_betas: FilePointer
    total_coefficients: List[FilePointer]
    orbital_coefficients: List[FilePointer]


class ExtractCoefficientsFromXMLProcess(Process):

    input_model = ExtractCoefficientsFromXMLInput
    output_model = ExtractCoefficientsFromXMLOutput

    def _run(self):
        """
        Performs the decomposition into radial basis functions and spherical harmonics
        """

        # Precompute the parameters of the radial basis functions
        alphas, betas = ml.precompute_parameters_of_radial_basis(self.inputs.n_max, self.inputs.l_max,
                                                                 self.inputs.r_min, self.inputs.r_max)

        # Save the parameters to file
        suffix = '_'.join(str(x) for x in [self.inputs.n_max, self.inputs.l_max, self.inputs.r_min, self.inputs.r_max])
        assert self.directory is not None
        alpha_file = f'alphas_{suffix}.npy'
        beta_file = f'betas_{suffix}.npy'
        alpha_filepointer = FilePointer(self, alpha_file)
        beta_filepointer = FilePointer(self, beta_file)
        utils.write_binary_content(self.directory / alpha_file, alphas.tobytes())
        utils.write_binary_content(self.directory / beta_file, betas.tobytes())

        # Compute the decomposition
        orbital_files, total_files = ml.compute_decomposition(
            alpha_file=alpha_filepointer, beta_file=beta_filepointer, **self.inputs.dict())

        self.outputs = ExtractCoefficientsFromXMLOutput(precomputed_alphas=alpha_filepointer,
                                                        precomputed_betas=beta_filepointer,
                                                        total_coefficients=[FilePointer(self, f) for f in total_files],
                                                        orbital_coefficients=[FilePointer(self, f) for f in orbital_files])


class ComputePowerSpectrumInput(IOModel):
    n_max: int
    l_max: int
    orbital_coefficients: FilePointer
    total_coefficients: FilePointer

    class Config:
        arbitrary_types_allowed = True


class ComputePowerSpectrumOutput(IOModel):
    power_spectrum: FilePointer

    class Config:
        arbitrary_types_allowed = True


def read_coeff_matrix(coeff_orb: np.ndarray, coeff_tot: np.ndarray, n_max: int, l_max: int) -> np.ndarray:
    """
    Reads the flat coefficient vector into a matrix with the correct dimensions.
    """

    coeff_matrix = np.zeros((2, n_max, l_max + 1, 2 * l_max + 1), dtype=float)
    idx = 0
    for n in range(n_max):
        for l in range(l_max + 1):
            for m in range(2 * l + 1):
                coeff_matrix[0, n, l, m] = coeff_orb[idx]
                coeff_matrix[1, n, l, m] = coeff_tot[idx]
                idx += 1
    return coeff_matrix


def compute_power_mat(coeff_matrix: np.ndarray, n_max: int, l_max: int) -> np.ndarray:
    """
    Computes the power_spectrum from the coefficient matrices.
    """

    # Note that we only store the inequivalent entries and hence the second for-loops of each elements iterate only over
    # the indices that are equal or larger than the corresponding index from the first for-loop.
    power = []
    for i1, _ in enumerate(['orb', 'tot']):
        for i2 in range(i1, 2):
            for n1 in range(n_max):
                for n2 in range(n1, n_max):
                    for l in range(l_max+1):
                        sum_current = sum(coeff_matrix[i1, n1, l, m]*coeff_matrix[i2, n2, l, m] for m in range(2*l+1))
                        power.append(sum_current)
    return np.array(power)


class ComputePowerSpectrumProcess(Process):

    input_model = ComputePowerSpectrumInput
    output_model = ComputePowerSpectrumOutput

    def _run(self):
        """
        Loads the coefficient vectors corresponding to the orbital and to the total density, computes the corresponding
        power spectrum and saves it to file.
        """

        # Load the coefficients from file
        coeff_orb = np.frombuffer(utils.get_binary_content(*self.inputs.orbital_coefficients))
        coeff_tot = np.frombuffer(utils.get_binary_content(*self.inputs.total_coefficients))
        coeff_matrix = read_coeff_matrix(coeff_orb, coeff_tot, self.inputs.n_max, self.inputs.l_max)

        # Compute the power spectrum
        power_mat = compute_power_mat(coeff_matrix, self.inputs.n_max, self.inputs.l_max)

        # Write the power spectrum to file
        utils.write_binary_content(f'power_spectrum.npy', power_mat.tobytes())
        self.outputs = self.output_model(power_spectrum=FilePointer(self, 'power_spectrum.npy'))
