"""
Processes used during the machine learning workflows
"""
from pathlib import Path
from typing import List

import numpy as np
from ase_koopmans.cell import Cell
from pydantic import ConfigDict

from koopmans import ml
from koopmans.bands import Band
from koopmans.files import File

from ._process import IOModel, Process


class ExtractCoefficientsFromXMLInput(IOModel):
    n_max: int
    l_max: int
    r_min: float
    r_max: float
    r_cut: float
    wannier_centers: List[List[float]]
    total_density_xml: File
    cell: Cell
    orbital_densities_xml: List[File]
    bands: List[Band]
    model_config = ConfigDict(arbitrary_types_allowed=True)


class ExtractCoefficientsFromXMLOutput(IOModel):
    precomputed_alphas: File
    precomputed_betas: File
    total_coefficients: List[File]
    orbital_coefficients: List[File]
    model_config = ConfigDict(arbitrary_types_allowed=True)


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

        # Save the parameters to files
        suffix = '_'.join(str(x) for x in [self.inputs.n_max, self.inputs.l_max, self.inputs.r_min, self.inputs.r_max])
        assert self.directory is not None

        alpha_file = File(self, f'alphas_{suffix}.npy')
        alpha_file.write_bytes(alphas.tobytes())

        beta_file = File(self, f'betas_{suffix}.npy')
        beta_file.write_bytes(betas.tobytes())

        # Compute the decomposition
        orbital_files, total_files = ml.compute_decomposition(
            alpha_file=alpha_file, beta_file=beta_file,
            **self.inputs.dict())

        # Write the files
        for name, content in list(orbital_files.items()) + list(total_files.items()):
            dst = File(self, name)
            if isinstance(content, bytes):
                dst.write_bytes(content)
            else:
                dst.write_text(content)

        total_coefficients = [File(self, f) for f in total_files.keys()]
        orbital_coefficients = [File(self, f) for f in orbital_files.keys()]
        self.outputs = ExtractCoefficientsFromXMLOutput(precomputed_alphas=alpha_file,
                                                        precomputed_betas=beta_file,
                                                        total_coefficients=total_coefficients,
                                                        orbital_coefficients=orbital_coefficients)


class ComputePowerSpectrumInput(IOModel):
    n_max: int
    l_max: int
    orbital_coefficients: File
    total_coefficients: File
    model_config = ConfigDict(arbitrary_types_allowed=True)


class ComputePowerSpectrumOutput(IOModel):
    power_spectrum: File
    model_config = ConfigDict(arbitrary_types_allowed=True)


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
                        sum_current = sum(coeff_matrix[i1, n1, l, m]*coeff_matrix[i2, n2, l, m] for m in range(2*l + 1))
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
        coeff_orb = np.frombuffer(self.inputs.orbital_coefficients.read_bytes())
        coeff_tot = np.frombuffer(self.inputs.total_coefficients.read_bytes())
        coeff_matrix = read_coeff_matrix(coeff_orb, coeff_tot, self.inputs.n_max, self.inputs.l_max)

        # Compute the power spectrum
        power_mat = compute_power_mat(coeff_matrix, self.inputs.n_max, self.inputs.l_max)

        # Write the power spectrum to file
        power_spectrum_file = File(self, Path('power_spectrum.npy'))
        power_spectrum_file.write_bytes(power_mat.tobytes())

        # Populate self.outputs
        self.outputs = self.output_model(power_spectrum=power_spectrum_file)
