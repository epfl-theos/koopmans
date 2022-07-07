from pathlib import Path
from typing import List, Tuple, Dict
from ase import Atoms
import numpy as np

from koopmans.bands import Band


def test_cart2sph(r_spherical: np.ndarray, debug_out: Path):
    """
    Write some quantities to file to test the function cart2sph.
    """

    with open(debug_out, 'a') as file:
        r_min = np.min(r_spherical[:, :, :, 0])
        r_max = np.max(r_spherical[:, :, :, 0])
        theta_min = np.min(r_spherical[:, :, :, 1])
        theta_max = np.max(r_spherical[:, :, :, 1])
        phi_min = np.min(r_spherical[:, :, :, 2])
        phi_max = np.max(r_spherical[:, :, :, 2])
        file.write(f"r_min     = {r_min}")
        file.write(f"r_max     = {r_max}")
        file.write(f"theta_min = {theta_min}")
        file.write(f"theta_max = {theta_max}")
        file.write(f"phi_min   = {phi_min}")
        file.write(f"phi_max   = {phi_max}")


def test_translate_to_new_integration_domain(f: np.ndarray, f_new: np.ndarray, wfc_center: np.ndarray, debug_out: Path):
    """
    Write some quantities to file to find problems in the function translate_to_new_integration_domain. 
    """

    with open(debug_out, 'a') as file:
        file.write(f"wfc_center                = ({wfc_center[0]}, {wfc_center[1]}, {wfc_center[2]})")
        file.write(f"Max of f on original grid = {np.max(f)}")
        file.write(f"Max of f on new grid:     = {np.max(f_new)}")


def test_decomposition(total_basis_array: np.ndarray, rho_r: np.ndarray, rho_r_xsf: np.ndarray, total_density_r_xsf: np.ndarray, coefficients_orbital: List[float], coefficients_total: List[float], nr_xml: Tuple[int, int, int], nr_new_integration_domain: Tuple[int, int, int], center_index: Tuple[int, int, int], band: Band, atoms: Atoms, write_to_xsf: bool, dirs: Dict[str, Path], debug_out: Path):
    """
    Write some quantities to outputfile to find problems in the compuation of the expansion coefficients.
    """

    with open(debug_out, 'a') as file:
        if band.filled:
            filled_str = 'occ'
        else:
            filled_str = 'emp'

        rho_r_reconstruced = get_reconstructed_orbital_densities(total_basis_array, coefficients_orbital)
        file.write(f"max rho_r_reconstructed                  = {np.max(rho_r_reconstruced)}")
        file.write(f"writing reconstructed orbital to to xsf file")
        rho_r_reconstruced = map_again_to_original_grid(
            rho_r_reconstruced, center_index, nr_xml, nr_new_integration_domain)
        file.write(f"max rho_r_reconstructed on original grid = {np.max(rho_r_reconstruced)}")
        difference = np.linalg.norm(rho_r_reconstruced-rho_r)
        file.write(f"Difference to original density           = {difference}")

        if write_to_xsf:
            rho_r_reconstructed_xsf = get_orbital_density_to_xsf_grid(rho_r_reconstruced, nr_xml)
            assert isinstance(band.index, int)
            filename_xsf = dirs['xsf'] / 'orbital.reconstructed.{}.{:05d}.xsf'.format(filled_str, band.index)
            print_to_xsf_file(filename_xsf, atoms, [rho_r_reconstructed_xsf], nr_xml)

        file.write(f"reconstruct total density")
        rho_r_reconstruced = get_reconstructed_orbital_densities(total_basis_array, coefficients_total)

        if write_to_xsf:
            file.write(f"writing reconstructed orbital to to xsf file")
            rho_r_reconstruced = map_again_to_original_grid(
                rho_r_reconstruced, center_index, nr_xml, nr_new_integration_domain)
            rho_r_reconstructed_xsf = get_orbital_density_to_xsf_grid(rho_r_reconstruced, nr_xml)
            assert isinstance(band.index, int)
            filename_xsf = dirs['xsf'] / 'total.reconstructed.{}.{:05d}.xsf'.format(filled_str, band.index)
            print_to_xsf_file(filename_xsf, atoms, [rho_r_reconstructed_xsf], nr_xml)

            file.write(f"writing total density minus orbital density to xsf file")
            assert isinstance(band.index, int)
            filename_xsf = dirs['xsf'] / 'total_minus.reconstructed.{}.{:05d}.xsf'.format(filled_str, band.index)
            print_to_xsf_file(filename_xsf, atoms, [total_density_r_xsf-rho_r_xsf], nr_xml)


def get_orbital_density_to_xsf_grid(rho_r_reconstructed: np.ndarray, nr_xml: Tuple[int, int, int]) -> np.ndarray:
    """
    Convert rho_r_reconstructed to the corresponding quantity in the xsf-format: (the last row at for each axis is identital to the first one).
    """

    rho_r_reconstruced_xsf = np.zeros((nr_xml[2], nr_xml[1], nr_xml[0]))
    for k in range(nr_xml[2]):
        for j in range(nr_xml[1]):
            for i in range(nr_xml[0]):
                rho_r_reconstruced_xsf[k, j, i] = rho_r_reconstructed[k %
                                                                      (nr_xml[2]-1), j % (nr_xml[1]-1), i % (nr_xml[0]-1)]
    return rho_r_reconstruced_xsf


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


def print_to_xsf_file(filename: Path, atoms: Atoms, list_vectors: List[np.ndarray], nr_xml: Tuple[int, int, int], wfc_centers: List[np.ndarray] = []):
    """
    Writes a quantity defined on the real space grid to a xsf file, which can be plotted with xcrysden    
    """

    cell_parameters = atoms.get_cell()
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    with open(filename, 'w') as out:
        out.write('# xsf file \n')
        out.write('CRYSTAL\n\n')
        out.write('PRIMVEC\n\n')
        for i in range(3):
            curr_string = ''
            for j in range(3):
                curr_string += "{:13.10f}".format(cell_parameters[i][j]) + " "
            out.write("\t" + curr_string + "\n")
        out.write('PRIMCOORD\n')
        out.write("\t" + str(len(symbols)+len(wfc_centers)) + '\t1\n')
        for i in range(len(symbols)):
            curr_string = symbols[i] + " "
            for j in range(3):
                curr_string += "{:13.10f}".format(positions[i][j]) + " "
            out.write("\t" + curr_string + "\n")
        # adding Nitrogen to the list of atoms to visualize my computed centers
        for i in range(len(wfc_centers)):
            curr_string = "N" + " "
            curr_string += "{:13.10f}".format(wfc_centers[i][2]) + " "
            curr_string += "{:13.10f}".format(wfc_centers[i][1]) + " "
            curr_string += "{:13.10f}".format(wfc_centers[i][0]) + " "
            out.write("\t" + curr_string + "\n")
        # end adding Nitrogen to list of atoms
        out.write('BEGIN_BLOCK_DATAGRID_3D\n')
        out.write("\t" + 'my_first_example_of_3D_datagrid\n')
        for i in range(len(list_vectors)):
            vector_r = list_vectors[i]
            out.write("\t" + 'BEGIN_DATAGRID_3D_this_is_3Dgrid#' + str(i+1) + '\n')
            out.write("\t" + "\t" + str(nr_xml[0]) + '\t' + str(nr_xml[1]) + '\t' + str(nr_xml[2]) + '\t\n')
            out.write("\t" + "\t" + str(0.0) + '\t' + str(0.0) + '\t' + str(0.0) + '\t\n')  # origin of the data grid
            # third spanning vector of the data grid
            out.write("\t" + "\t" + str(cell_parameters[0][0]) + '\t' + str(0.0) + '\t' + str(0.0) + '\t\n')
            # second spanning vector of the data grid
            out.write("\t" + "\t" + str(0.0) + '\t' + str(cell_parameters[1][1]) + '\t' + str(0.0) + '\t\n')
            out.write("\t" + "\t" + str(0.0) + '\t' + str(0.0) + '\t' +
                      str(cell_parameters[2][2]) + '\t\n')  # first spanning vector of the data grid
            for k in range(nr_xml[2]):
                for j in range(nr_xml[1]):
                    out.write("\t\t")
                    for i in range(nr_xml[0]):
                        out.write("{:.15E}\t".format(vector_r[k, j, i]))
                    out.write('\n')
                out.write("\n\n")
            out.write("\n\t" + 'END_DATAGRID_3D\n')
        out.write('END_BLOCK_DATAGRID_3D')
