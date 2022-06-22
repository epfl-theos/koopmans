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


def main_compute_power(n_max, l_max, r_min, r_max, ML_directory, dir_power, bands):
    dir_coeff = ML_directory / ('coefficients_' + '_'.join(str(x) for x in [n_max, l_max, r_min, r_max]))
    dir_power.mkdir(exist_ok=True)
    dir_orb   = dir_coeff / 'coff_orb'
    dir_tot   = dir_coeff / 'coff_tot'


    for band in bands:

        if band.filled:
            filled_str = 'occ'
        else:
            filled_str = 'emp'
        
        coff_orb        = np.atleast_1d(np.loadtxt(dir_orb / f'coff.orbital.{filled_str}.{band.index}.txt'))
        coff_tot        = np.atleast_1d(np.loadtxt(dir_tot / f'coff.total.{filled_str}.{band.index}.txt'))
        coff_matrix     = read_in(coff_orb, coff_tot, n_max, l_max)
        power_mat       = compute_power(coff_matrix, n_max, l_max)
        np.savetxt(dir_power / f"power_spectrum.orbital.{filled_str}.{band.index}.txt", power_mat)
        

