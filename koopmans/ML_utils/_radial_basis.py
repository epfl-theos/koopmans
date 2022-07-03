from pathlib import Path
from typing import Dict
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import os
import sys

def phi(r: np.ndarray,l: int,alpha: float):
    return r**l*np.exp(-alpha*r**2)

def compute_overlap(n: int,n_prime: int,l: int,alphas: np.ndarray):
    integrand = lambda r: r**2*phi(r,l,alphas[n])*phi(r,l,alphas[n_prime])
    return quad(integrand, 0.0, np.inf)

# https://github.com/pyscf/pyscf.github.io/blob/master/examples/dmrg/32-dmrg_casscf_nevpt2_for_FeS.py
def lowdin(s: np.ndarray):
    e, v = np.linalg.eigh(s)
    return np.dot(v/np.sqrt(e), v.T.conj())

def compute_s(n_max: int, l:int, alphas:np.ndarray):
    s = np.zeros((n_max, n_max))
    for n in range(n_max):
        for n_prime in range(n_max):
            tmp = compute_overlap(n,n_prime, l, alphas)[0]
            s[n,n_prime] = tmp
    return s

def compute_beta(n_max:int , l: int, alphas: np.ndarray):
    s = compute_s(n_max, l, alphas)
    beta = lowdin(s)
    return beta


def compute_alphas(n_max:int, r_thrs:np.ndarray, thr:float=10**(-3.0)):    
    alphas = np.zeros(n_max)
    for n in range(n_max):
        r_thr = r_thrs[n]
        f = lambda alpha: phi(r_thr,0,alpha)-thr
        root = fsolve(f, [0.01])
        alphas[n] = root
    return alphas

def g(r: np.ndarray,n:int,n_max:int,l:int,betas: np.ndarray,alphas: np.ndarray):
    return sum(betas[n_prime, n, l]*phi(r, l, alphas[n_prime]) for n_prime in range(n_max))


def precompute_radial_basis(n_max:int , l_max:int, r_min_thr:float, r_max_thr:float, dirs:Dict[str,Path]):
    orig_stdout = sys.stdout
    with open(dirs['ML'] / 'orbitals_to_power_spectra.out', 'a') as f:
        sys.stdout = f 

        print("\nprecompute radial basis\n") 
        

        thr   = 10**(-3)
        r_thrs = np.zeros(n_max)
        # for i in range(n_max):
        #     r_thrs[i] = r_min_thr + i*(r_max_thr-1)/n_max

        if n_max==1:
            r_thrs[0] = r_min_thr
        else:
            for i in range(n_max):
                r_thrs[i] = r_min_thr + i*(r_max_thr-r_min_thr)/(n_max-1.0)

        print('Threshold values = ', r_thrs)

        alphas = compute_alphas(n_max, r_thrs, thr)

        betas = np.zeros((n_max, n_max, l_max))

        for l in range(l_max):
            betas[:,:,l] = compute_beta(n_max, l, alphas)
        

        betas.tofile(dirs['betas'] / ('betas_'    + '_'.join(str(x) for x in [n_max, l_max, r_min_thr, r_max_thr]) + '.dat'))
        alphas.tofile(dirs['alphas'] / ('alphas_'    + '_'.join(str(x) for x in [n_max, l_max, r_min_thr, r_max_thr]) + '.dat'))
    sys.stdout = orig_stdout



    # r = np.linspace(0,4,1000)
    # for l in range(l_max):
    #     plt.figure()
    #     for n in range(n_max):
    #         f = lambda r: phi(r,l,alphas[n])
    #         plt.plot(r, f(r), label = "n = " + str(n+1))
    #     if(l==0):
    #         plt.vlines(r_thrs, ymin=0.0, ymax=1.0)
    #     plt.plot(r, thr*np.ones(len(r)), label = "threshold value")
    #     plt.legend()
    #     plt.savefig(curr_folder + '/plots/phi_n'+ str(l) + '.png')


    # r = np.linspace(0,3,1000)
    # for l in range(l_max):
    #     plt.figure()
    #     for n in range(n_max):
    #         f = lambda r: g(r,n,n_max,l,betas,alphas)
    #         plt.plot(r, f(r), label = "n = " + str(n+1))
    #     plt.legend()
    #     plt.savefig(curr_folder + '/plots/g_n'+ str(l) + '.png')
