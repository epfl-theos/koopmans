from typing import List, Tuple

import numpy as np


def reconstruct_orbital_densities(total_basis_array: np.ndarray, coefficients: List[float]) -> np.ndarray:
    """
    Reconstruct the density with the truncated expansion coefficients multiplied with the corresponding basis functions
    """

    return np.einsum('ijkl,l->ijk', total_basis_array, coefficients)


def map_to_original_grid(array: np.ndarray, wfc_center_index: Tuple[int, int, int], nr_xml: Tuple[int, int, int], nr_new: Tuple[int, int, int]):
    """
    Maps the array defined on the new integration domain back to the unit cell.
    """

    array_on_reg_grid = np.zeros((nr_xml[2]-1, nr_xml[1]-1, nr_xml[0]-1), dtype=float)

    for k_new, k in enumerate(range(wfc_center_index[0] - nr_new[2], wfc_center_index[0] + nr_new[2]+1)):
        for j_new, j in enumerate(range(wfc_center_index[1] - nr_new[1], wfc_center_index[1] + nr_new[1]+1)):
            for i_new, i in enumerate(range(wfc_center_index[2] - nr_new[0], wfc_center_index[2] + nr_new[0]+1)):
                array_on_reg_grid[k % (nr_xml[2] - 1), j % (nr_xml[1] - 1), i %
                                  (nr_xml[0] - 1)] = array[k_new, j_new, i_new]

    return array_on_reg_grid
