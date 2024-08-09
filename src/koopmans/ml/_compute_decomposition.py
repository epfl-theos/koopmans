from functools import partial
from pathlib import Path
from typing import Callable, Dict, List, Tuple
from xml.etree import ElementTree as ET

import numpy as np
from ase import Atoms, units
from ase.cell import Cell
from numpy.linalg import norm

from koopmans.bands import Band
from koopmans.files import FilePointer
from koopmans.utils import (get_binary_content, get_content,
                            write_binary_content)

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
                     total_basis_function_array: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the expansion coefficients of rho and rho_total wrt the basis defined in total_basis_function_array.
    """

    rho_tmp = np.expand_dims(rho, axis=3)
    rho_total_tmp = np.expand_dims(rho_total, axis=3)

    integrand_rho = rho_tmp*total_basis_function_array
    integrand_rho_total = rho_total_tmp*total_basis_function_array

    coefficients = compute_3d_integral_naive(integrand_rho, r_cartesian).flatten()

    coefficients_total = compute_3d_integral_naive(integrand_rho_total, r_cartesian).flatten()

    return coefficients, coefficients_total


def compute_decomposition(n_max: int, l_max: int, r_min: float, r_max: float, r_cut: float, total_density_xml: FilePointer,
                          orbital_densities_xml: Dict[str, FilePointer],
                          bands: List[Band], cell: Cell, wannier_centers: np.ndarray, alphas: FilePointer, betas: FilePointer) -> Tuple[List[str], List[str]]:
    """
    Computes the expansion coefficients of the total and orbital densities.
    """

    # Define the normalisation constant for densities
    norm_const = 1/(units.Bohr)**3

    # Load the grid dimensions nr_xml from charge-density-file
    raw_filecontents = get_content(*total_density_xml)
    xml_root = ET.fromstringlist(raw_filecontents)
    xml_charge_density = xml_root.find('CHARGE-DENSITY')
    assert xml_charge_density is not None
    xml_info = xml_charge_density.find('INFO')
    assert xml_info is not None
    nr_xml = tuple([int(x) + 1 for x in [xml_info.get(f'nr{i+1}') for i in range(3)] if x is not None])
    assert len(nr_xml) == 3

    # load the lattice parameters
    lat_vecs = np.array([cell[2, 2], cell[1, 1], cell[0, 0]])

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

    # Define our radial basis functions, which are partially parameterized by precomputed vectors
    alphas = np.frombuffer(get_binary_content(*alphas)).reshape((n_max, l_max+1))
    betas = np.frombuffer(get_binary_content(*betas)).reshape((n_max, n_max, l_max+1))
    radial_basis_functions: RadialBasisFunctions = partial(g_basis, betas=betas, alphas=alphas)

    # Compute R_nl Y_lm for each point on the integration domain
    total_basis_array = precompute_basis_function(
        radial_basis_functions, real_spherical_harmonics_basis_functions, r_cartesian, r_spherical, n_max, l_max)

    # load the total charge density
    raw_filecontents = get_content(*total_density_xml)
    xml_root = ET.fromstringlist(raw_filecontents)
    assert xml_root is not None
    xml_charge_density = xml_root.find('CHARGE-DENSITY')
    assert xml_charge_density is not None
    total_density_r = parse_xml_array(xml_charge_density, nr_xml, norm_const)

    orbital_files = []
    total_files = []

    # Compute the decomposition for each band
    assert len(orbital_densities_xml) == len(bands)
    for orbital_density_xml, band in zip(orbital_densities_xml, bands):

        if band.filled:
            filled_str = 'occ'
        else:
            filled_str = 'emp'

        # load the orbital density
        raw_filecontents = get_content(*orbital_density_xml)
        xml_root = ET.fromstringlist(raw_filecontents)
        assert xml_root is not None
        xml_charge_density = xml_root.find('EFFECTIVE-POTENTIAL')
        assert xml_charge_density is not None
        rho_r = parse_xml_array(xml_charge_density, nr_xml, norm_const)

        # Bring the the density to the same integration domain as the precomputed basis, centered around the orbital's
        # center, making sure that the center is in within the unit cell
        assert band.index is not None
        wfc_center_tmp = wannier_centers[band.index-1]
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
        orbital_file = f'coeff.orbital.{filled_str}.{band.index}.txt'
        write_binary_content(orbital_file, coefficients_orbital.tobytes())
        total_file = f'coeff.total.{filled_str}.{band.index}.txt'
        write_binary_content(total_file, coefficients_total.tobytes())

        orbital_files.append(orbital_file)
        total_files.append(total_file)

    return orbital_files, total_files


def parse_xml_array(
    xml_root: ET.Element, nr: Tuple[int, int, int], norm_const: float, retain_final_element: bool = False
) -> np.ndarray:
    """
    Loads an array from an xml file.

    :param xml_root: The xml root containing the array
    :param norm_const: The normalization constant to multiply the array with (in our case 1/((Bohr radii)^3)
    :param string: The name of the field in the xml file that contains the array, in our case either
    'EFFECTIVE-POTENTIAL' or 'CHARGE-DENSITY'
    :param retain_final_element: If True, the array is returned in with periodic boundary conditions, i.e. the last
    element in each dimension is equal to the first element in each dimension. This is required for the xsf format.

    :return: The array
    """

    # Extract the array
    array_xml = np.zeros((nr[2], nr[1], nr[0]), dtype=float)

    for k in range(nr[2]):
        current_name = 'z.' + str(k % (nr[2] - 1) + 1)
        entry = xml_root.find(current_name)
        assert isinstance(entry, ET.Element)
        text = entry.text
        assert isinstance(text, str)
        rho_tmp = np.array(text.split(), dtype=float)
        for j in range(nr[1]):
            for i in range(nr[0]):
                array_xml[k, j, i] = rho_tmp[(j % (nr[1] - 1))*(nr[0] - 1) + (i % (nr[0] - 1))]
    array_xml *= norm_const

    if retain_final_element:
        # the xsf format requires an array where the last element is equal to the first element in each dimension
        return array_xml
    else:
        return array_xml[:-1, :-1, :-1]
