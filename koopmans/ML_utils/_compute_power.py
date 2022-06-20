import numpy as np
import os
import sys

def read_in(coff_orb, coff_tot, n_max, l_max):
    coff_matrix = np.zeros((2, n_max, l_max, 2*l_max+1), dtype=float)
    idx = 0
    for n in range(n_max):
        for l in range(l_max):
            for m in range(2*l+1):
                coff_matrix[0, n,l,m] = coff_orb[idx]
                idx += 1
    idx = 0
    for n in range(n_max):
        for l in range(l_max):
            for m in range(2*l+1):
                coff_matrix[1, n,l,m] = coff_tot[idx]
                idx += 1
    return coff_matrix

def compute_power(coff_matrix, n_max, l_max):
    power = []
    for i1, _ in enumerate(['orb', 'tot']):
        for i2 in range(i1, 2):
            for n1 in range(n_max):
                for n2 in range(n1, n_max):
                    for l in range(l_max):
                        sum_current = sum(coff_matrix[i1, n1, l, m]*coff_matrix[i2,n2,l,m] for m in range(2*l+1))
                        power.append(sum_current)
    return np.array(power)


def main_compute_power(n_max, l_max, r_min, r_max, ML_directory, power_dir, nbnd):
    dir_coeff = ML_directory + '/coefficients_'    + str(n_max) + '_' + str(l_max) + '_' + str(r_min) + '_' + str(r_max)
    dir_power = power_dir
    os.system('mkdir -p ' + dir_power)
    dir_orb   = dir_coeff + '/coff_orb/'
    dir_tot   = dir_coeff + '/coff_tot/'

    size_power_mat = sum(i for i in range(1,n_max+1))*l_max*3

    print(size_power_mat)


    if nbnd[1] > 0:

        num_orbitals = nbnd[1]#len([name for name in os.listdir(dir_orb)])# only emp orbital for now

        power_mat    = np.zeros((num_orbitals, size_power_mat))
        print(np.shape(power_mat))

        for i in range(num_orbitals):
            coff_orb        = np.atleast_1d(np.loadtxt(dir_orb + 'coff.orbital.occ.' + str(i+1) + '.txt'))
            coff_tot        = np.atleast_1d(np.loadtxt(dir_tot + 'coff.total.occ.' + str(i+1) + '.txt'))
            coff_matrix     = read_in(coff_orb, coff_tot, n_max, l_max)
            power_mat[i,:]  = compute_power(coff_matrix, n_max, l_max)
            np.savetxt(dir_power + "/power_spectrum.orbital.occ." + str(i) + ".txt", power_mat)
    
    if nbnd[0] > 0:

        num_orbitals = nbnd[0]#len([name for name in os.listdir(dir_orb)])# only emp orbital for now

        power_mat = np.zeros((num_orbitals, size_power_mat))
        print(np.shape(power_mat))

        for i in range(num_orbitals):
            coff_orb        = np.atleast_1d(np.loadtxt(dir_orb + 'coff.orbital.emp.' + str(i+1) + '.txt'))
            coff_tot        = np.atleast_1d(np.loadtxt(dir_tot + 'coff.total.emp.' + str(i+1) + '.txt'))
            coff_matrix     = read_in(coff_orb, coff_tot, n_max, l_max)
            power_mat[i,:]  = compute_power(coff_matrix, n_max, l_max)
            np.savetxt(dir_power+ "/power_spectrum.orbital.emp." + str(i) + ".txt", power_mat)
