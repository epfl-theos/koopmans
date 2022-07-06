from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
from ase import Atoms
from numpy.linalg import norm
import xml.etree.ElementTree as ET
from numpy.linalg import norm
import math
from scipy.special import sph_harm
import sys
import koopmans.ml_utils._got_basis as gb

from koopmans.bands import Bands
from koopmans import utils

import koopmans.ml_utils.debugging as debugging


# Say what you want to compute
Debug = True
write_to_xsf = True


def real_spherical_harmonics(theta: float, phi: float, l: int, m: int):
    Y = sph_harm(abs(m), l, phi, theta)
    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real
    return Y.real


def cart2sph(x: float, y: float, z: float):
    XsqPlusYsq = x**2 + y**2
    r = math.sqrt(XsqPlusYsq + z**2)
    theta = math.atan2(z, math.sqrt(XsqPlusYsq))+np.pi/2.0
    phi = math.atan2(y, x) + np.pi
    return r, theta, phi


def cart2sph_array(r_cartesian: np.ndarray):
    (k_max, j_max, i_max, _) = np.shape(r_cartesian)
    r_spherical = np.zeros_like(r_cartesian)
    for k in range(k_max):
        for j in range(j_max):
            for i in range(i_max):
                r_spherical[k, j, i, :] = cart2sph(
                    r_cartesian[k, j, i, 2], r_cartesian[k, j, i, 1], r_cartesian[k, j, i, 0])
    return r_spherical


def precompute_basis_function(r_cartesian: np.ndarray, r_spherical: np.ndarray, n_max: int, l_max: int, betas: np.ndarray, alphas: np.ndarray):
    Y_array_all = np.zeros((np.shape(r_cartesian)[0], np.shape(r_cartesian)[
                           1], np.shape(r_cartesian)[2], l_max+1, 2*l_max+1))
    for l in range(l_max+1):
        for i, m in enumerate(range(-l, l+1)):
            Y_array_all[:, :, :, l, i] = real_spherical_harmonics(
                r_spherical[:, :, :, 1], r_spherical[:, :, :, 2], l, m)

    g_array_all = np.zeros((np.shape(r_cartesian)[0], np.shape(r_cartesian)[
                           1], np.shape(r_cartesian)[2], n_max, l_max+1))

    for n in range(n_max):
        for l in range(l_max+1):
            g_array_all[:, :, :, n, l] = gb.radial_basis_function(r_spherical[:, :, :, 0], n, n_max, l, betas, alphas)

    number_of_l_elements = sum(2*l+1 for l in range(0, l_max+1))
    total_basis_function_array = np.zeros((np.shape(r_cartesian)[0], np.shape(
        r_cartesian)[1], np.shape(r_cartesian)[2], n_max*number_of_l_elements))
    idx = 0
    for n in range(n_max):
        for l in range(l_max+1):
            total_basis_function_array[:, :, :, idx:(
                idx+2*l+1)] = np.expand_dims(g_array_all[:, :, :, n, l], axis=3)*Y_array_all[:, :, :, l, 0:2*l+1]
            idx += 2*l+1

    return total_basis_function_array


def get_coefficients(rho: np.ndarray, rho_total: np.ndarray, r_cartesian: np.ndarray, total_basis_function_array: np.ndarray):
    coefficients: List[np.ndarray] = []
    coefficients_total: List[np.ndarray] = []

    rho_tmp = np.expand_dims(rho, axis=3)
    rho_total_tmp = np.expand_dims(rho_total, axis=3)

    integrand_rho = rho_tmp*total_basis_function_array
    integrand_rho_total = rho_total_tmp*total_basis_function_array

    c = compute_3d_integral_naive_specialized(integrand_rho, r_cartesian).flatten()
    coefficients[len(coefficients):] = list(c)

    c_total = compute_3d_integral_naive_specialized(integrand_rho_total, r_cartesian).flatten()
    coefficients_total[len(coefficients_total):] = list(c_total)

    if Debug:
        print("np.norm(coefficients): ", np.linalg.norm(coefficients))
    return coefficients, coefficients_total


def get_index(r: np.ndarray, value: np.ndarray):
    norms = norm(r - value, axis=3)
    idx = np.unravel_index(np.argmin(norms), np.shape(r)[:-1])
    return idx


def compute_3d_integral_naive_specialized(f: np.ndarray, r: np.ndarray):
    z = r[:, 0, 0, 0]
    y = r[0, :, 0, 1]
    x = r[0, 0, :, 2]
    result_n = np.sum(f, axis=(0, 1, 2))*(x[-1]-x[0])*(y[-1]-y[0])*(z[-1]-z[0])/((len(x)-1)*(len(y)-1)*(len(z)-1))
    return result_n


def generate_integration_box(r: np.ndarray, r_cut: float):
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
    r_new = np.zeros((2*nr_new_integration_domain[2]+1, 2 *
                     nr_new_integration_domain[1]+1, 2*nr_new_integration_domain[0]+1, 3))
    z_, y_, x_ = np.meshgrid(z_, y_, x_, indexing='ij')
    r_new[:, :, :, 0] = z_
    r_new[:, :, :, 1] = y_
    r_new[:, :, :, 2] = x_

    return nr_new_integration_domain, r_new


def generate_integration_domain(f: np.ndarray, wfc_center_index: Tuple[int, int, int], nr_new_integration_domain: Tuple[int, int, int]):

    f_rolled = np.roll(f, (-(wfc_center_index[0]-nr_new_integration_domain[2]), -(wfc_center_index[1] -
                       nr_new_integration_domain[1]), -(wfc_center_index[2]-nr_new_integration_domain[0])), axis=(0, 1, 2))
    f_new = f_rolled[:2*nr_new_integration_domain[2]+1, :2 *
                     nr_new_integration_domain[1]+1, :2*nr_new_integration_domain[0]+1]

    if Debug:
        print("Max of f on original grid: ", np.max(f))
        print("Max of f on new grid:      ", np.max(f_new))

    return f_new


def load_grid_dimension_from_xml_file(file: Path, string: str = 'CHARGE-DENSITY') -> Tuple[int, int, int]:

    with open(file, 'r') as fd:
        tree = ET.parse(fd)
        rho_file = tree.getroot()
        rho_file_charge_density = rho_file.find(string)
        assert isinstance(rho_file_charge_density, ET.Element)
        info = rho_file_charge_density.find('INFO')
        assert isinstance(info, ET.Element)
        nr_xml_list = [0, 0, 0]
        for i in range(3):
            attr = 'nr' + str(i+1)
            info_i = info.get(attr)
            assert isinstance(info_i, str)
            nr_xml_list[i] = int(info_i) + 1
        nr_xml: Tuple[int, int, int] = (nr_xml_list[0], nr_xml_list[1], nr_xml_list[2])
    return nr_xml


def load_density_into_array(file_rho: Path, nr_xml: Tuple[int, int, int], norm_const: float, string: str = 'EFFECTIVE-POTENTIAL'):
    rho_r_xsf = np.zeros((nr_xml[2], nr_xml[1], nr_xml[0]), dtype=float)
    rho_r = np.zeros((nr_xml[2]-1, nr_xml[1]-1, nr_xml[0]-1), dtype=float)

    with open(file_rho, 'r') as fd:
        tree = ET.parse(fd)
    rho_file = tree.getroot()

    for k in range(nr_xml[2]):
        current_name = 'z.' + str(k % (nr_xml[2]-1)+1)
        rho_file_str = rho_file.find(string)
        assert isinstance(rho_file_str, ET.Element)
        rho_file_str_current_name = rho_file_str.find(current_name)
        assert isinstance(rho_file_str_current_name, ET.Element)
        rho_file_str_current_name_text = rho_file_str_current_name.text
        assert isinstance(rho_file_str_current_name_text, str)
        rho_tmp = np.array(rho_file_str_current_name_text.split('\n')[1:-1], dtype=float)
        for j in range(nr_xml[1]):
            for i in range(nr_xml[0]):
                rho_r_xsf[k, j, i] = rho_tmp[(j % (nr_xml[1]-1))*(nr_xml[0]-1)+(i %
                                                                                (nr_xml[0]-1))]  # np.dot(r[k,j,i, :], r[k,j,i, :])#
    rho_r_xsf *= norm_const
    rho_r[:, :, :] = rho_r_xsf[:-1, :-1, :-1]

    return rho_r, rho_r_xsf


def func_compute_decomposition(n_max: int, l_max: int, r_min: float, r_max: float, r_cut: float, dirs: Dict[str, Path], bands: Bands, atoms: Atoms, centers: np.ndarray):
    if write_to_xsf:
        dir_xsf = dirs['ml'] / 'xsf'
        utils.system_call(f'mkdir -p {dir_xsf}')
        dirs.update({'xsf': dir_xsf})

    # redirect print-statements away from std-output
    orig_stdout = sys.stdout
    with open(dirs['ml'] / 'orbitals_to_power_spectra.out', 'a') as f:
        sys.stdout = f

        # TODO: find the normalization constant 6.748334698446981
        norm_const = 6.748334698446981

        # load the grid dimensions nr_xml from charge-density-file
        file_rho = dirs['xml'] / 'charge-density.xml'
        nr_xml = load_grid_dimension_from_xml_file(file_rho)

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
                        np.array([float(k % (nr_xml[2]-1))/(nr_xml[2]-1), float(j % (nr_xml[1]-1))/(nr_xml[1]-1), float(i % (nr_xml[0]-1))/(nr_xml[0]-1)]), lat_vecs)
        r[:, :, :, :] = r_xsf[:-1, :-1, :-1, :]

        # Define an alternative grid which is used to perform the integrations. This can be identical or smaller as the original grid but not larger
        nr_new_integration_domain, r_cartesian = generate_integration_box(r, r_cut)

        # Convert the cartesian coordinates to spherical coordinates for simplifying the integration with spherical harmonics
        r_spherical = cart2sph_array(r_cartesian)

        if Debug:
            debugging.test_cart2sph(r_spherical)

        # load precomputed vectors defining the radial basis functions
        betas = np.fromfile(dirs['betas'] / ('betas_' + '_'.join(str(x)
                            for x in [n_max, l_max, r_min, r_max]) + '.dat')).reshape((n_max, n_max, l_max+1))
        alphas = np.fromfile(dirs['alphas'] / ('alphas_' + '_'.join(str(x)
                             for x in [n_max, l_max, r_min, r_max]) + '.dat')).reshape(n_max, l_max+1)

        # Compute R_nl Y_lm for each point on the integration domain
        total_basis_array = precompute_basis_function(r_cartesian, r_spherical, n_max, l_max, betas, alphas)

        # load the total charge density
        file_rho = dirs['xml'] / 'charge-density.xml'
        total_density_r, total_density_r_xsf = load_density_into_array(file_rho, nr_xml, norm_const, 'CHARGE-DENSITY')

        # for debugging print the total charge density to a xsf file that can be plotted with xcrysden
        if Debug:
            if write_to_xsf:
                filename_xsf = dirs['xsf'] / 'charge-density.xsf'
                debugging.print_to_xsf_file(filename_xsf, atoms, [total_density_r_xsf], nr_xml)

        # Compute the decomposition for each band
        for band in bands:

            if band.filled:
                filled_str = 'occ'
            else:
                filled_str = 'emp'

            # load the orbital density
            file_rho = dirs['xml'] / 'orbital.{}.{}.{:05d}.xml'.format(filled_str, band.spin, band.index)
            rho_r, rho_r_xsf = load_density_into_array(file_rho, nr_xml, norm_const)

            # for debugging print the orbital charge density to a xsf file that can be plotted with xcrysden
            if Debug:
                filename_xsf = dirs['xsf'] / 'orbital.{}.{:05d}.xsf'.format(filled_str, band.index)
                debugging.print_to_xsf_file(filename_xsf, atoms, [rho_r_xsf], nr_xml)

            # center the integration around the center of the orbital
            wfc_center_tmp = centers[band.index-1]
            wfc_center = np.array([wfc_center_tmp[2] % lat_vecs[0], wfc_center_tmp[1] %
                                   lat_vecs[1], wfc_center_tmp[0] % lat_vecs[2]])
            center_index = get_index(r, wfc_center)
            rho_r_new = generate_integration_domain(rho_r, center_index, nr_new_integration_domain)
            total_density_r_new = generate_integration_domain(total_density_r, center_index, nr_new_integration_domain)

            if Debug:
                print(wfc_center)

            # compute the decomposition coefficients
            if band.filled:
                coefficients_orbital, coefficients_total = get_coefficients(
                    rho_r_new, total_density_r_new, r_cartesian, total_basis_array)
            else:
                coefficients_orbital, coefficients_total = get_coefficients(
                    rho_r_new, total_density_r_new, r_cartesian, total_basis_array)

            # save the decomposition coefficients in files
            np.savetxt(dirs['coeff_orb'] / f'coff.orbital.{filled_str}.{band.index}.txt', coefficients_orbital)
            np.savetxt(dirs['coeff_tot'] / f'coff.total.{filled_str}.{band.index}.txt', coefficients_total)

            if Debug:
                debugging.test_decomposition(total_basis_array, rho_r, rho_r_xsf, total_density_r_xsf, coefficients_orbital,
                                             coefficients_total, nr_xml, nr_new_integration_domain, center_index, band, atoms, write_to_xsf, dirs)

    sys.stdout = orig_stdout
