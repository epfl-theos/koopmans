from pathlib import Path
from typing import Dict, List
import numpy as np
from ase import Atoms
from numpy.linalg import norm
import xml.etree.ElementTree as ET
from numpy.linalg import norm 
import math 
from scipy.special import sph_harm
import sys

from koopmans.bands import Bands


def get_orbital_density_to_xsf_grid(rho_r_reconstructed: np.ndarray):
    nr3m1, nr2m1, nr1m1 = np.shape(rho_r_reconstructed)
    nr3 = nr3m1 + 1 
    nr2 = nr2m1 + 1
    nr1 = nr1m1 + 1
    rho_r_reconstruced_xsf = np.zeros((nr3, nr2, nr1))
    for k in range(nr3):
        for j in range(nr2):
            for i in range(nr1):
                rho_r_reconstruced_xsf[k,j,i] = rho_r_reconstructed[k%(nr3-1),j%(nr2-1), i%(nr1-1)]
    return rho_r_reconstruced_xsf


def get_reconstructed_orbital_densities(total_basis_array: np.ndarray, coefficients: np.ndarray):
    rho_r_reconstruced = np.einsum('ijkl,l->ijk', total_basis_array, coefficients)
    return rho_r_reconstruced

def phi(r: np.ndarray,l: int, alpha: float):
    return r**l*np.exp(-alpha*r**2)


def g(r: np.ndarray, n: int,n_max: int,l,betas: np.ndarray,alphas: np.ndarray):
    return sum(betas[n_prime, n, l]*phi(r, l, alphas[n_prime]) for n_prime in range(n_max))


def radial_basis_function(r: np.ndarray,n:int ,n_max:int ,l:int,betas: np.ndarray,alphas: np.ndarray):
    return g(r,n,n_max,l,betas,alphas)


# from https://scipython.com/blog/visualizing-the-real-forms-of-the-spherical-harmonics/

def real_spherical_harmonics(theta:float, phi:float, l:int, m:int):
    Y = sph_harm(abs(m), l, phi, theta) #Yannick: not sure about order of theta and phi
    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real
    return Y.real #Yannick: added .real to only get real values

def cart2sph(x:float,y:float,z:float):
    XsqPlusYsq = x**2 + y**2
    r = math.sqrt(XsqPlusYsq + z**2)               # r
    theta = math.atan2(z,math.sqrt(XsqPlusYsq))+np.pi/2.0    # theta
    phi = math.atan2(y,x) + np.pi                          # phi
    return r, theta, phi

def cart2sph_array(r_cartesian: np.ndarray):
    (k_max, j_max, i_max, _) = np.shape(r_cartesian)
    r_spherical                = np.zeros_like(r_cartesian)
    for k in range(k_max):
        for j in range(j_max):
            for i in range(i_max):
                r_spherical[k,j,i,:] = cart2sph(r_cartesian[k,j,i,2], r_cartesian[k,j,i,1], r_cartesian[k,j,i,0])
    return r_spherical




def precompute_basis_function(r_cartesian: np.ndarray, r_spherical: np.ndarray, n_max: int, l_max: int, betas: np.ndarray, alphas: np.ndarray):
    Y_array_all = np.zeros((np.shape(r_cartesian)[0], np.shape(r_cartesian)[1], np.shape(r_cartesian)[2], l_max, 2*l_max+1))
    for l in range(l_max):
        for i, m in enumerate(range(-l, l+1)):
            Y_array_all[:,:,:,l,i] = real_spherical_harmonics(r_spherical[:,:,:,1], r_spherical[:,:,:,2], l, m)
    
    g_array_all = np.zeros((np.shape(r_cartesian)[0], np.shape(r_cartesian)[1], np.shape(r_cartesian)[2], n_max, l_max))

    for n in range(n_max):
        for l in range(l_max):
            g_array_all[:,:,:,n,l] = radial_basis_function(r_spherical[:,:,:,0],n,n_max,l,betas,alphas)

    number_of_l_elements       = sum(2*l+1 for l in range(0,l_max))
    total_basis_function_array = np.zeros((np.shape(r_cartesian)[0], np.shape(r_cartesian)[1], np.shape(r_cartesian)[2], n_max*number_of_l_elements))
    idx = 0
    for n in range(n_max):
        for l in range(l_max):
            total_basis_function_array[:,:,:,idx:(idx+2*l+1)]     = np.expand_dims(g_array_all[:,:,:,n,l], axis=3)*Y_array_all[:,:,:,l,0:2*l+1]
            idx += 2*l+1
    
    return total_basis_function_array


def get_coefficients(rho: np.ndarray, rho_total: np.ndarray, r_cartesian: np.ndarray, total_basis_function_array: np.ndarray):    
    print("starting with the loop to compute the coefficients")
    coefficients: List[np.ndarray] = []
    coefficients_total: List[np.ndarray] = []

    rho_tmp       = np.expand_dims(rho, axis=3)
    rho_total_tmp = np.expand_dims(rho_total, axis=3)

    integrand_rho         = rho_tmp*total_basis_function_array
    integrand_rho_total   = rho_total_tmp*total_basis_function_array

    c       = compute_3d_integral_naive_specialized(integrand_rho, r_cartesian).flatten()
    coefficients[len(coefficients):] = list(c)
            
    c_total = compute_3d_integral_naive_specialized(integrand_rho_total, r_cartesian).flatten()
    coefficients_total[len(coefficients_total):] = list(c_total)

    if Debug:
        print("np.norm(coefficients): ", np.linalg.norm(coefficients))
    return coefficients, coefficients_total


def shift_coordinates(r: np.ndarray, index:List):
    return r-r[index[0],index[1],index[2],:]


def get_index(r: np.ndarray, value: np.ndarray):
    norms = norm(r - value, axis=3)
    idx = np.unravel_index(np.argmin(norms), np.shape(r)[:-1])
    return idx


def compute_3d_integral_naive_specialized(f: np.ndarray, r: np.ndarray):
    z = r[:,0,0,0]
    y = r[0,:,0,1]
    x = r[0,0,:,2]
    result_n  = np.sum(f, axis=(0,1,2))*(x[-1]-x[0])*(y[-1]-y[0])*(z[-1]-z[0])/((len(x)-1)*(len(y)-1)*(len(z)-1))
    return result_n


def generate_integration_box(r: np.ndarray, r_cut:float):
    z = r[:,0,0,0]
    y = r[0,:,0,1]
    x = r[0,0,:,2]
    dz = z[1]-z[0]
    dy = y[1]-y[0]
    dx = x[1]-x[0]
    # number_of_x_grid_points = int(r_cut/dx)
    # number_of_y_grid_points = int(r_cut/dy)
    # number_of_z_grid_points = int(r_cut/dz)

    number_of_x_grid_points = min(int(r_cut/dx),len(x)//2-1)
    number_of_y_grid_points = min(int(r_cut/dy),len(y)//2-1)
    number_of_z_grid_points = min(int(r_cut/dz),len(z)//2-1)

    z_ = dz*np.arange(-number_of_z_grid_points, number_of_z_grid_points+1)
    y_ = dy*np.arange(-number_of_y_grid_points, number_of_y_grid_points+1)
    x_ = dx*np.arange(-number_of_x_grid_points, number_of_x_grid_points+1)
    r_new = np.zeros((2*number_of_z_grid_points+1, 2*number_of_y_grid_points+1, 2*number_of_x_grid_points+1,3))
    z_, y_, x_ = np.meshgrid(z_, y_, x_, indexing = 'ij')
    r_new[:,:,:,0] = z_
    r_new[:,:,:,1] = y_
    r_new[:,:,:,2] = x_

    return number_of_x_grid_points, number_of_y_grid_points, number_of_z_grid_points, r_new


def generate_integration_domain(f: np.ndarray, wfc_center_index: np.ndarray, number_of_x_grid_points:int, number_of_y_grid_points:int, number_of_z_grid_points:int):
    
    f_rolled = np.roll(f, (-(wfc_center_index[0]-number_of_z_grid_points), -(wfc_center_index[1]-number_of_y_grid_points), -(wfc_center_index[2]-number_of_x_grid_points)), axis=(0,1,2))    
    f_new = f_rolled[:2*number_of_z_grid_points+1, :2*number_of_y_grid_points+1, :2*number_of_x_grid_points+1]

    if Debug:
        print("Max of f on original grid: ", np.max(f))
        print("Max of f on new grid:      ", np.max(f_new))

    return f_new
    
    


def map_again_to_original_grid(f_new: np.ndarray, wfc_center_index: np.ndarray, number_of_x_grid_points:int, number_of_y_grid_points:int, number_of_z_grid_points:int, nr1:int, nr2:int, nr3:int):

    f_on_reg_grid = np.zeros((nr3-1, nr2-1, nr1-1), dtype=float)

    for k_new, k in enumerate(range(wfc_center_index[0]-number_of_z_grid_points, wfc_center_index[0]+number_of_z_grid_points+1)):
            for j_new, j in enumerate(range(wfc_center_index[1]-number_of_y_grid_points, wfc_center_index[1]+number_of_y_grid_points+1)):
                for i_new, i in enumerate(range(wfc_center_index[2]-number_of_x_grid_points, wfc_center_index[2]+number_of_x_grid_points+1)):
                    f_on_reg_grid[k%(nr3-1),j%(nr2-1),i%(nr1-1)] = f_new[k_new,j_new,i_new]
    
    return f_on_reg_grid



def get_coefficient_matrix(c: np.ndarray, n_max:int, l_max:int):
    c_matrix = np.zeros((n_max, l_max))
    idx = 0
    for n in range(n_max):
        for l in range(l_max):
            for m in range(-l, l+1):
                c_matrix[n,l] += np.abs(c[idx])
                idx += 1
    return c_matrix



def print_to_xsf_file(filename:Path, cell_parameters: np.ndarray, positions: np.ndarray, symbols:List[str], list_of_evcs:List[np.ndarray], nr1:int, nr2:int, nr3:int, wfc_centers:List[np.ndarray]=[]):
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
        for i in range(len(list_of_evcs)):
            evc_r = list_of_evcs[i]
            out.write("\t" + 'BEGIN_DATAGRID_3D_this_is_3Dgrid#' + str(i+1) + '\n')
            out.write("\t" + "\t" + str(nr1) + '\t' + str(nr2)+ '\t' + str(nr3) + '\t\n')
            out.write("\t" + "\t" + str(0.0) + '\t' + str(0.0)+ '\t' + str(0.0) + '\t\n') #origin of the data grid
            out.write("\t" + "\t" + str(cell_parameters[0][0]) + '\t' + str(0.0)+ '\t' + str(0.0) + '\t\n') #third spanning vector of the data grid
            out.write("\t" + "\t" + str(0.0) + '\t' + str(cell_parameters[1][1])+ '\t' + str(0.0) + '\t\n') #second spanning vector of the data grid
            out.write("\t" + "\t" + str(0.0) + '\t' + str(0.0)+ '\t' + str(cell_parameters[2][2]) + '\t\n') #first spanning vector of the data grid
            
            for k in range(nr3):
                for j in range(nr2):
                    out.write("\t\t")
                    for i in range(nr1):
                        out.write("{:.15E}\t".format(evc_r[k, j, i]))
                    out.write('\n')
                out.write("\n\n")
            out.write("\n\t" + 'END_DATAGRID_3D\n')                     
        out.write('END_BLOCK_DATAGRID_3D')


def load_density_into_array(file_rho:Path, nr1:int, nr2:int, nr3:int, norm_const:float, string:str='EFFECTIVE-POTENTIAL'):#string='EFFECTIVE-POTENTIAL'):
    rho_r_xsf = np.zeros((nr3, nr2, nr1), dtype=float)
    rho_r     = np.zeros((nr3-1, nr2-1, nr1-1), dtype=float)

    with open(file_rho, 'r') as fd:
        tree = ET.parse(fd)
    rho_file = tree.getroot()

    for k in range(nr3):
        current_name = 'z.' + str(k%(nr3-1)+1)
        rho_file_str = rho_file.find(string)
        assert isinstance(rho_file_str, ET.Element)
        rho_file_str_current_name = rho_file_str.find(current_name)
        assert isinstance(rho_file_str_current_name, ET.Element)
        rho_file_str_current_name_text = rho_file_str_current_name.text
        assert isinstance(rho_file_str_current_name_text, str)
        rho_tmp = np.array(rho_file_str_current_name_text.split('\n')[1:-1], dtype=float)
        for j in range(nr2):
            for i in range(nr1):
                rho_r_xsf[k, j, i] = rho_tmp[(j%(nr2-1))*(nr1-1)+(i%(nr1-1))]#np.dot(r[k,j,i, :], r[k,j,i, :])#
    rho_r_xsf   *= norm_const
    rho_r[:,:,:] = rho_r_xsf[:-1,:-1,:-1]

    return rho_r, rho_r_xsf

# Say what you want to compute
compute_decomposition                   = True
Debug                                   = False
write_to_xsf                            = False

def func_compute_decomposition(n_max:int, l_max:int, r_min:float, r_max:float, r_cut:float, dirs: Dict[str, Path], bands:Bands, atoms:Atoms, centers: np.ndarray):
    orig_stdout = sys.stdout
    with open(dirs['ML'] / 'orbitals_to_power_spectra.out', 'a') as f:
        sys.stdout = f       

        print("\ncompute decomposition\n")

        # TODO Don't hardcode these parameters
        norm_const      = 6.748334698446981
        # TODO end
        # TODO: find the normalization constant 6.748334698446981

        

        file_rho   = dirs['xml'] / 'charge-density.xml'
        with open(file_rho, 'r') as fd:
            tree = ET.parse(fd)
        rho_file = tree.getroot()
        rho_file_charge_density = rho_file.find('CHARGE-DENSITY')
        nrs = [0,0,0]
        assert isinstance(rho_file_charge_density, ET.Element)
        info = rho_file_charge_density.find('INFO')
        assert isinstance(info, ET.Element)
        for i in range(3):
            attr = 'nr' + str(i+1)
            info_i = info.get(attr)
            assert isinstance(info_i, str)
            nrs[i]  = int(info_i) + 1
        nr1 = nrs[0]
        nr2 = nrs[1]
        nr3 = nrs[2]
        
        positions       = atoms.get_positions()
        symbols         = atoms.get_chemical_symbols()
        cell_parameters = atoms.get_cell()

        # centers_occ = parse_w90(file_path_rho_base+'wann_occ.wout')


        # load precomputed vectors defining the radial basis functions 
        betas  = np.fromfile(dirs['betas']  / ('betas_'    + '_'.join(str(x) for x in [n_max, l_max, r_min, r_max]) + '.dat')).reshape((n_max, n_max, l_max))
        alphas = np.fromfile(dirs['alphas'] / ('alphas_'    + '_'.join(str(x) for x in [n_max, l_max, r_min, r_max]) + '.dat')).reshape(n_max)

        lat_vecs = np.array([cell_parameters[2,2], cell_parameters[1,1], cell_parameters[0,0]])


        print("writing cartesian grid into an array")
        r_xsf     = np.zeros((nr3, nr2, nr1, 3), dtype=float)
        r         = np.zeros((nr3-1, nr2-1, nr1-1, 3), dtype=float)
        for k in range(nr3):
            for j in range(nr2):
                for i in range(nr1):
                    r_xsf[k,j,i, :] = np.multiply(np.array([1.0*(k%(nr3-1))/(nr3-1), 1.0*(j%(nr2-1))/(nr2-1), 1.0*(i%(nr1-1))/(nr1-1)]), lat_vecs)

        r[:,:,:,:] = r_xsf[:-1,:-1,:-1,:]


        print("computing the alternative cartesian grid on which the actual integration is performed")
        number_of_x_grid_points, number_of_y_grid_points, number_of_z_grid_points, r_cartesian = generate_integration_box(r, r_cut)
        print("compute r_spherical from r_cartesian")
        r_spherical = cart2sph_array(r_cartesian)
        
        if Debug:
            r_min     = np.min(r_spherical[:,:,:,0])
            r_max     = np.max(r_spherical[:,:,:,0])
            theta_min = np.min(r_spherical[:,:,:,1])
            theta_max = np.max(r_spherical[:,:,:,1])
            phi_min   = np.min(r_spherical[:,:,:,2])
            phi_max   = np.max(r_spherical[:,:,:,2])
            print("r_min   = ", r_min) 
            print("r_max   = ", r_max)
            print("theta_min = ", theta_min) 
            print("theta_max = ", theta_max)
            print("phi_min   = ", phi_min) 
            print("phi_max   = ", phi_max)

        
        print("writing the total density into an array")
        file_rho = dirs['xml'] / 'charge-density.xml'
        total_density_r, total_density_r_xsf = load_density_into_array(file_rho, nr1, nr2, nr3, norm_const, 'CHARGE-DENSITY')

        if Debug:
            differences_to_original           = []

            if write_to_xsf:
                print("writing total density to xsf file")
                filename_xsf = dirs['ML'] / 'charge-density.xsf'
                print_to_xsf_file(filename_xsf, cell_parameters, positions, symbols, [total_density_r_xsf], nr1, nr2, nr3, [])


        total_basis_array = precompute_basis_function(r_cartesian, r_spherical, n_max, l_max, betas, alphas)

        for band in bands:
            print(f'calculating decomposition of band {band.index}') 

            if band.filled:
                filled_str = 'occ'
            else:
                filled_str = 'emp'
            


            file_rho = dirs['xml']  / 'orbital.{}.{}.{:05d}.xml'.format(filled_str,band.spin, band.index) 
            
            rho_r, rho_r_xsf = load_density_into_array(file_rho, nr1, nr2, nr3, norm_const)

            if Debug:
                print("writing orbital to xsf file")
                filename_xsf = dirs['ML'] /  'orbital.{}.{:05d}.xsf'.format(filled_str, band.index) 
                print_to_xsf_file(filename_xsf, cell_parameters, positions, symbols, [rho_r_xsf], nr1, nr2, nr3, [])
                # os.system(str('rm ' + file_rho))


            if compute_decomposition:
                print("Set up the new integration domain")
                wfc_center_tmp      = centers[band.index-1]
                wfc_center          = np.array([wfc_center_tmp[2]%lat_vecs[0], wfc_center_tmp[1]%lat_vecs[1], wfc_center_tmp[0]%lat_vecs[2]])
                print(wfc_center)
                center_index        = get_index(r, wfc_center)
                rho_r_new           = generate_integration_domain(rho_r,           center_index, number_of_x_grid_points, number_of_y_grid_points, number_of_z_grid_points)                
                total_density_r_new = generate_integration_domain(total_density_r, center_index, number_of_x_grid_points, number_of_y_grid_points, number_of_z_grid_points)
                print("computing coefficients of the orbital and total density")
                if band.filled:
                    coefficients_orbital, coefficients_total = get_coefficients(rho_r_new, total_density_r_new, r_cartesian, total_basis_array)
                else:
                    coefficients_orbital, coefficients_total = get_coefficients(rho_r_new, total_density_r_new, r_cartesian, total_basis_array)

                np.savetxt(dirs['coeff_orb'] / f'coff.orbital.{filled_str}.{band.index}.txt', coefficients_orbital)
                np.savetxt(dirs['coeff_tot'] / f'coff.total.{filled_str}.{band.index}.txt', coefficients_total)


                if Debug:
                    print("reconstruct orbital density")
                    rho_r_reconstruced      = get_reconstructed_orbital_densities(total_basis_array, coefficients_orbital)
                    print("max rho_r_reconstructed                  = ", np.max(rho_r_reconstruced))
                    print("writing reconstructed orbital to to xsf file")
                    rho_r_reconstruced      = map_again_to_original_grid(rho_r_reconstruced, center_index, number_of_x_grid_points, number_of_y_grid_points, number_of_z_grid_points, nr1, nr2, nr3)
                    print("max rho_r_reconstructed on original grid = ", np.max(rho_r_reconstruced))
                    difference              = np.linalg.norm(rho_r_reconstruced-rho_r)
                    print("Difference to original density           = ", difference)
                    differences_to_original.append(difference)
                    


                    # c_matrix                = get_coefficient_matrix(coefficients_orbital, n_max, l_max)
                    if write_to_xsf:
                        rho_r_reconstructed_xsf = get_orbital_density_to_xsf_grid(rho_r_reconstruced)
                        filename_xsf = dirs['ML'] /  'orbital.{}.{:05d}.reconstructed.xsf'.format(band.filled, band.index) 
                        print_to_xsf_file(filename_xsf, cell_parameters, positions, symbols, [rho_r_reconstructed_xsf], nr1, nr2, nr3, [])


                    print("reconstruct total density")
                    rho_r_reconstruced      = get_reconstructed_orbital_densities(total_basis_array, coefficients_total)
                    # c_matrix                = get_coefficient_matrix(coefficients_total, n_max, l_max)
                    
                    if write_to_xsf:
                        print("writing reconstructed orbital to to xsf file")
                        rho_r_reconstruced      = map_again_to_original_grid(rho_r_reconstruced, center_index, number_of_x_grid_points, number_of_y_grid_points, number_of_z_grid_points, nr1, nr2, nr3)
                        rho_r_reconstructed_xsf = get_orbital_density_to_xsf_grid(rho_r_reconstruced)
                        filename_xsf = dirs['ML'] /  'total.{}.{:05d}.reconstructed.xsf'.format(band.filled, band.index) 
                        print_to_xsf_file(filename_xsf, cell_parameters, positions, symbols, [rho_r_reconstructed_xsf], nr1, nr2, nr3, [])


                        print("writing total density minus orbital density to xsf file")
                        filename_xsf = dirs['ML'] /  'total_minus.{}.{:05d}.reconstructed.xsf'.format(band.filled, band.index) 
                        print_to_xsf_file(filename_xsf, cell_parameters, positions, symbols, [total_density_r_xsf-rho_r_xsf], nr1, nr2, nr3, [])

                print("length of coefficient vector = ", len(coefficients_orbital))
        if Debug:
            return np.mean(differences_to_original)
    sys.stdout = orig_stdout
  


