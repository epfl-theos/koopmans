from pathlib import Path
from typing import Dict
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import koopmans.ml_utils._got_basis as gb
import os
import sys


def compute_alphas(n_max: int, l_max: int, r_thrs: np.ndarray, thr: float = 10**(-3.0)):
    alphas = np.zeros((n_max, l_max+1))
    for n in range(n_max):
        for l in range(l_max+1):
            alphas[n, l] = -1/r_thrs[n]**2*np.log(thr/(r_thrs[n]**l))
    return alphas


def compute_overlap(n: int, n_prime: int, l: int, alphas: np.ndarray):
    def integrand(r): return r**2*gb.phi(r, l, alphas[n, l])*gb.phi(r, l, alphas[n_prime, l])
    return quad(integrand, 0.0, np.inf)[0]


def lowdin(s: np.ndarray):
    e, v = np.linalg.eigh(s)
    return np.dot(v/np.sqrt(e), v.T.conj())


def compute_s(n_max: int, l: int, alphas: np.ndarray):
    s = np.zeros((n_max, n_max))
    for n in range(n_max):
        for n_prime in range(n_max):
            tmp = compute_overlap(n, n_prime, l, alphas)
            s[n, n_prime] = tmp
    return s


def compute_beta(n_max: int, l: int, alphas: np.ndarray):
    s = compute_s(n_max, l, alphas)
    beta = lowdin(s)
    return beta


def precompute_radial_basis(n_max: int, l_max: int, r_min_thr: float, r_max_thr: float, dirs: Dict[str, Path]):
    orig_stdout = sys.stdout
    with open(dirs['ml'] / 'orbitals_to_power_spectra.out', 'a') as f:
        sys.stdout = f

        print("\nprecompute radial basis\n")

        thr = 10**(-3)
        r_thrs = np.zeros(n_max)

        r_thrs = np.linspace(r_min_thr, r_max_thr, n_max)
        print('Threshold values = ', r_thrs)

        alphas = compute_alphas(n_max, l_max, r_thrs, thr)

        betas = np.zeros((n_max, n_max, l_max+1))

        for l in range(l_max+1):
            betas[:, :, l] = compute_beta(n_max, l, alphas)

        betas.tofile(dirs['betas'] / ('betas_' + '_'.join(str(x)
                                                          for x in [n_max, l_max, r_min_thr, r_max_thr]) + '.dat'))
        alphas.tofile(dirs['alphas'] / ('alphas_' + '_'.join(str(x)
                                                             for x in [n_max, l_max, r_min_thr, r_max_thr]) + '.dat'))
    sys.stdout = orig_stdout

    # r = np.linspace(0, 4, 1000)
    # for l in range(l_max):
    #     plt.figure()
    #     for n in range(n_max):
    #         def f(r): return gb.phi(r, l, alphas[n, l])
    #         plt.plot(r, f(r), label="n = " + str(n+1))
    #     if(l == 0):
    #         plt.vlines(r_thrs, ymin=0.0, ymax=1.0)
    #     plt.plot(r, thr*np.ones(len(r)), label="threshold value")
    #     plt.legend()
    #     plt.savefig(curr_folder + '/plots/phi_n' + str(l) + '.png')

    # r = np.linspace(0, 3, 1000)
    # for l in range(l_max):
    #     plt.figure()
    #     for n in range(n_max):
    #         def f(r): return gb.g(r, n, n_max, l, betas, alphas)
    #         plt.plot(r, f(r), label="n = " + str(n+1))
    #     plt.legend()
    #     plt.savefig(curr_folder + '/plots/g_n' + str(l) + '.png')

    # return alphas, betas
