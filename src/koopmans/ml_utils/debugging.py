from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from ase import Atoms

from koopmans.bands import Band


def get_reconstructed_orbital_densities(total_basis_array: np.ndarray, coefficients: List[float]) -> np.ndarray:
    """
    Reconstruct the density with the truncated expansion coefficients multiplied with the corresponding basis functions.
    """

    rho_r_reconstruced = np.einsum('ijkl,l->ijk', total_basis_array, coefficients)
    return rho_r_reconstruced


def map_again_to_original_grid(f_new: np.ndarray, wfc_center_index: Tuple[int, int, int], nr_xml: Tuple[int, int, int], nr_new_integration_domain: Tuple[int, int, int]):
    """
    Maps the function f_new defined on the new integration domain back to the unit cell.
    """

    f_on_reg_grid = np.zeros((nr_xml[2]-1, nr_xml[1]-1, nr_xml[0]-1), dtype=float)

    for k_new, k in enumerate(range(wfc_center_index[0]-nr_new_integration_domain[2], wfc_center_index[0]+nr_new_integration_domain[2]+1)):
        for j_new, j in enumerate(range(wfc_center_index[1]-nr_new_integration_domain[1], wfc_center_index[1]+nr_new_integration_domain[1]+1)):
            for i_new, i in enumerate(range(wfc_center_index[2]-nr_new_integration_domain[0], wfc_center_index[2]+nr_new_integration_domain[0]+1)):
                f_on_reg_grid[k % (nr_xml[2]-1), j % (nr_xml[1]-1), i % (nr_xml[0]-1)] = f_new[k_new, j_new, i_new]

    return f_on_reg_grid
