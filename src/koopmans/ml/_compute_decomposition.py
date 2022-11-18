from functools import partial
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import numpy as np
from ase import Atoms, units
from numpy.linalg import norm

from koopmans.bands import Bands
from koopmans.utils import read_xml_array, read_xml_nr

from ._basis_functions import g as g_basis
from ._basis_functions import \
    real_spherical_harmonics as real_spherical_harmonics_basis_functions
from ._misc import cart2sph_array, compute_3d_integral_naive

# Basis functions types
RadialBasisFunctions = Callable[[np.ndarray, int, int, int], np.ndarray]
SphericalBasisFunctions = Callable[[np.ndarray, np.ndarray, int, int], np.ndarray]


def precompute_basis_function(radial_basis_functions: RadialBasisFunctions,
                              spherical_basis_functions: SphericalBasisFunctions,
                              r_cartesian: np.ndarray, r_spherical: np.ndarray, n_max: int, l_max: int) -> np.ndarray:
    """
    Precomputes the total basis function (radial_basis_function*spherical_basis_function) for each point on the
    integration domain.
    """

    # Define the vector containing the values of the spherical basis function for each grid point for each pair of
    # (n,l).
    Y_array_all = np.zeros((np.shape(r_cartesian)[:3] + (l_max+1, 2*l_max+1)))
    for l in range(l_max+1):
        for i, m in enumerate(range(-l, l+1)):
            Y_array_all[:, :, :, l, i] = spherical_basis_functions(
                r_spherical[:, :, :, 1], r_spherical[:, :, :, 2], l, m)

    # Define the vector containing the values of the radial basis function for each grid point for each pair of (l,m).
    g_array_all = np.zeros((np.shape(r_cartesian)[:3] + (n_max, l_max+1)))
    for n in range(n_max):
        for l in range(l_max+1):
            g_array_all[:, :, :, n, l] = radial_basis_functions(r_spherical[:, :, :, 0], n, n_max, l)

    # Compute the vector containing the values of the total basis function for each grid point for each pair of (n,l,m)
    # All values corresponding to different values for m are stored in the last axis of total_basis_function_array
    number_of_l_elements = sum(2*l+1 for l in range(0, l_max+1))
    total_basis_function_array = np.zeros((np.shape(r_cartesian)[:3] + (n_max*number_of_l_elements,)))
    idx = 0
    for n in range(n_max):
        for l in range(l_max+1):
            total_basis_function_array[:, :, :, idx:(
                idx+2*l+1)] = np.expand_dims(g_array_all[:, :, :, n, l], axis=3)*Y_array_all[:, :, :, l, 0:2*l+1]
            idx += 2*l+1

    return total_basis_function_array


# functions to make sure that integration domain is chosen such that the center of the orbital density is the center
# of the integration domain

def get_index(r: np.ndarray, vec: np.ndarray) -> Tuple[int, int, int]:
    """
    Returns the index of the array r that is closest to vec.
    """

    norms = norm(r - vec, axis=3)
    idx_tmp = np.unravel_index(np.argmin(norms), np.shape(r)[:-1])
    idx = (int(idx_tmp[0]), int(idx_tmp[1]), int(idx_tmp[2]))
    return idx


def generate_integration_box(r: np.ndarray, r_cut: float) -> Tuple[Tuple[int, int, int], np.ndarray]:
    """
    Defines the cartesian coordinates of the new integration domain.

    This new integration domain is cubic, has the same grid spacing (dx,dy,dz) as the original grid
    but the mesh-size can be smaller (depending) on the cutoff value r_cut.
    """

    z = r[:, 0, 0, 0]
    y = r[0, :, 0, 1]
    x = r[0, 0, :, 2]
    dz = z[1]-z[0]
    dy = y[1]-y[0]
    dx = x[1]-x[0]

    nr_new_integration_domain: Tuple[int, int, int] = (
        min(int(r_cut/dx), len(x)//2-1), min(int(r_cut/dy), len(y)//2-1), min(int(r_cut/dz), len(z)//2-1))

    z_ = dz*np.arange(-nr_new_integration_domain[2], nr_new_integration_domain[2]+1)
    y_ = dy*np.arange(-nr_new_integration_domain[1], nr_new_integration_domain[1]+1)
    x_ = dx*np.arange(-nr_new_integration_domain[0], nr_new_integration_domain[0]+1)
    r_new = np.zeros((2 * nr_new_integration_domain[2] + 1,
                      2 * nr_new_integration_domain[1] + 1,
                      2 * nr_new_integration_domain[0] + 1, 3))
    z_, y_, x_ = np.meshgrid(z_, y_, x_, indexing='ij')
    r_new[:, :, :, 0] = z_
    r_new[:, :, :, 1] = y_
    r_new[:, :, :, 2] = x_

    return nr_new_integration_domain, r_new


def translate_to_new_integration_domain(f: np.ndarray, wfc_center_index: Tuple[int, int, int],
                                        nr_new_integration_domain: Tuple[int, int, int]) -> np.ndarray:
    """
    Rolls the array f, such that it is centered around wfc_center_index and brings it into the same shape as the new
    integration domain.
    """

    f_rolled = np.roll(f, (-(wfc_center_index[0] - nr_new_integration_domain[2]),
                           -(wfc_center_index[1] - nr_new_integration_domain[1]),
                           -(wfc_center_index[2] - nr_new_integration_domain[0])),
                       axis=(0, 1, 2))
    f_new = f_rolled[:2 * nr_new_integration_domain[2] + 1,
                     :2 * nr_new_integration_domain[1] + 1,
                     :2 * nr_new_integration_domain[0] + 1]
    return f_new


# functions to compute the expansion coefficients

def get_coefficients(rho: np.ndarray, rho_total: np.ndarray, r_cartesian: np.ndarray,
                     total_basis_function_array: np.ndarray) -> Tuple[List[float], List[float]]:
    """
    Computes the expansion coefficients of rho and rho_total wrt the basis defined in total_basis_function_array.
    """

    coefficients: List[float] = []
    coefficients_total: List[float] = []

    rho_tmp = np.expand_dims(rho, axis=3)
    rho_total_tmp = np.expand_dims(rho_total, axis=3)

    integrand_rho = rho_tmp*total_basis_function_array
    integrand_rho_total = rho_total_tmp*total_basis_function_array

    c = compute_3d_integral_naive(integrand_rho, r_cartesian).flatten()
    coefficients[len(coefficients):] = list(c)

    c_total = compute_3d_integral_naive(integrand_rho_total, r_cartesian).flatten()
    coefficients_total[len(coefficients_total):] = list(c_total)

    return coefficients, coefficients_total


def compute_decomposition(n_max: int, l_max: int, r_min: float, r_max: float, r_cut: float, dirs: Dict[str, Path],
                          bands: Bands, atoms: Atoms, centers: np.ndarray):
    """
    Computes the expansion coefficients of the total and orbital densities.
    """

    # Define the normalisation constant for densities
    norm_const = 1/(units.Bohr)**3

    # load the grid dimensions nr_xml from charge-density-file
    file_rho = dirs['xml'] / 'charge-density.xml'
    nr_xml = read_xml_nr(file_rho, 'CHARGE-DENSITY')

    # load the lattice parameters
    cell_parameters = atoms.get_cell()
    lat_vecs = np.array([cell_parameters[2, 2], cell_parameters[1, 1], cell_parameters[0, 0]])

    # Define the cartesian grid
    r_xsf = np.zeros((nr_xml[2], nr_xml[1], nr_xml[0], 3), dtype=float)
    r = np.zeros((nr_xml[2]-1, nr_xml[1]-1, nr_xml[0]-1, 3), dtype=float)
    for k in range(nr_xml[2]):
        for j in range(nr_xml[1]):
            for i in range(nr_xml[0]):
                r_xsf[k, j, i, :] = np.multiply(
                    np.array([float(k % (nr_xml[2] - 1))/(nr_xml[2] - 1),
                              float(j % (nr_xml[1] - 1))/(nr_xml[1] - 1),
                              float(i % (nr_xml[0] - 1))/(nr_xml[0] - 1)]), lat_vecs)
    r[:, :, :, :] = r_xsf[:-1, :-1, :-1, :]

    # Define an alternative grid which is used to perform the integrations. This can be identical or smaller as the
    # original grid but not larger
    nr_new_integration_domain, r_cartesian = generate_integration_box(r, r_cut)

    # Convert the cartesian coordinates to spherical coordinates for simplifying the integration with spherical
    # harmonics
    r_spherical = cart2sph_array(r_cartesian)

    # Define our radial basis functions, which are partially parametrised by precomputed vectors
    betas = np.fromfile(dirs['betas'] / ('betas_' + '_'.join(str(x)
                        for x in [n_max, l_max, r_min, r_max]) + '.dat')).reshape((n_max, n_max, l_max+1))
    alphas = np.fromfile(dirs['alphas'] / ('alphas_' + '_'.join(str(x) for x in [n_max, l_max, r_min, r_max])
                                           + '.dat')).reshape(n_max, l_max+1)
    radial_basis_functions: RadialBasisFunctions = partial(g_basis, betas=betas, alphas=alphas)

    # Compute R_nl Y_lm for each point on the integration domain
    total_basis_array = precompute_basis_function(
        radial_basis_functions, real_spherical_harmonics_basis_functions, r_cartesian, r_spherical, n_max, l_max)

    # load the total charge density
    total_density_r = read_xml_array(dirs['xml'] / 'charge-density.xml', norm_const, 'CHARGE-DENSITY')

    # Compute the decomposition for each band
    for band in bands:

        if band.filled:
            filled_str = 'occ'
        else:
            filled_str = 'emp'

        # load the orbital density
        rho_r = read_xml_array(dirs['xml'] / f'orbital.{filled_str}.{band.spin}.{band.index:05d}.xml', norm_const)

        # Bring the the density to the same integration domain as the precomputed basis, centered around the orbital's
        # center, making sure that the center is in within the unit cell
        wfc_center_tmp = centers[band.index-1]
        wfc_center = np.array([wfc_center_tmp[2] % lat_vecs[0], wfc_center_tmp[1] %
                               lat_vecs[1], wfc_center_tmp[0] % lat_vecs[2]])
        center_index = get_index(r, wfc_center)
        rho_r_new = translate_to_new_integration_domain(rho_r, center_index, nr_new_integration_domain)
        total_density_r_new = translate_to_new_integration_domain(
            total_density_r, center_index, nr_new_integration_domain)

        # compute the decomposition coefficients
        coefficients_orbital, coefficients_total = get_coefficients(
            rho_r_new, total_density_r_new, r_cartesian, total_basis_array)

        # save the decomposition coefficients in files
        np.savetxt(dirs['coeff_orb'] / f'coff.orbital.{filled_str}.{band.index}.txt', coefficients_orbital)
        np.savetxt(dirs['coeff_tot'] / f'coff.total.{filled_str}.{band.index}.txt', coefficients_total)
