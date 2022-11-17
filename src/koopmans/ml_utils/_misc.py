from typing import List, Tuple

import numpy as np


def cart2sph_array(r_cartesian: np.ndarray) -> np.ndarray:
    """
    Converts an array of cartesian coordinates into the corresponding spherical ones.

    Note that cartesian is z, y, x; spherical is r, theta, phi
    """

    xy2 = r_cartesian[:, :, :, 2]**2 + r_cartesian[:, :, :, 1]**2
    r_spherical = np.zeros_like(r_cartesian)
    r_spherical[:, :, :, 0] = np.linalg.norm(r_cartesian, axis=-1)
    r_spherical[:, :, :, 1] = np.arctan2(r_cartesian[:, :, :, 0], np.sqrt(xy2)) + np.pi / 2.0
    r_spherical[:, :, :, 2] = np.arctan2(r_cartesian[:, :, :, 1], r_cartesian[:, :, :, 2]) + np.pi
    return r_spherical


def compute_3d_integral_naive(f: np.ndarray, r: np.ndarray) -> np.ndarray:
    """
    Computes the 3d-integral of an array of functions f on r with a simple trapezoidal rule.

    Note this assumes an orthorhombic cell!
    """

    z = r[:, 0, 0, 0]
    y = r[0, :, 0, 1]
    x = r[0, 0, :, 2]
    result_n = np.sum(f, axis=(0, 1, 2))*(x[-1]-x[0])*(y[-1]-y[0])*(z[-1]-z[0])/((len(x)-1)*(len(y)-1)*(len(z)-1))
    return result_n


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
