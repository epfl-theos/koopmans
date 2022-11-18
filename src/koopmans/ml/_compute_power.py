from pathlib import Path
from typing import Dict

import numpy as np

from koopmans.bands import Bands


def read_coeff_matrix(coff_orb: np.ndarray, coff_tot: np.ndarray, n_max: int, l_max: int) -> np.ndarray:
    """
    Reads the flat coefficient vector into a matrix with the correct dimensions.
    """

    coff_matrix = np.zeros((2, n_max, l_max + 1, 2 * l_max + 1), dtype=float)
    idx = 0
    for n in range(n_max):
        for l in range(l_max + 1):
            for m in range(2 * l + 1):
                coff_matrix[0, n, l, m] = coff_orb[idx]
                coff_matrix[1, n, l, m] = coff_tot[idx]
                idx += 1
    return coff_matrix


def compute_power_mat(coff_matrix: np.ndarray, n_max: int, l_max: int) -> np.ndarray:
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
                        sum_current = sum(coff_matrix[i1, n1, l, m]*coff_matrix[i2, n2, l, m] for m in range(2*l+1))
                        power.append(sum_current)
    return np.array(power)


def compute_power(n_max: int, l_max: int, dirs: Dict[str, Path], bands: Bands,
                  input_vectors_fo_ml: Dict[str, np.ndarray]):
    """
    Loads the coefficient vectors corresponding to the orbital and to the total density, computes the corresponding
    power spectrum and saves it to a file.
    """

    for band in bands:

        if band.filled:
            filled_str = 'occ'
        else:
            filled_str = 'emp'

        coff_orb = np.atleast_1d(np.loadtxt(dirs['coeff_orb'] / f'coff.orbital.{filled_str}.{band.index}.txt'))
        coff_tot = np.atleast_1d(np.loadtxt(dirs['coeff_tot'] / f'coff.total.{filled_str}.{band.index}.txt'))
        coff_matrix = read_coeff_matrix(coff_orb, coff_tot, n_max, l_max)
        power_mat = compute_power_mat(coff_matrix, n_max, l_max)
        np.savetxt(dirs['power'] / f"power_spectrum.orbital.{filled_str}.{band.index}.txt", power_mat)
        input_vectors_fo_ml[f"power_spectrum.orbital.{filled_str}.{band.index}"] = power_mat
