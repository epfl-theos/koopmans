from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

import koopmans.ml_utils._basis_functions as basis

Debug = False


def compute_alphas(n_max: int, l_max: int, r_thrs: np.ndarray, thr: float):
    """
    Computes the decay-coefficients alpha_nl, by demanding that for each r_thrs[n],
    the corresponding phi_nl decays to threshold value thr at a cutoff radius of r_thr.
    """

    alphas = np.zeros((n_max, l_max+1))
    for n in range(n_max):
        for l in range(l_max+1):
            alphas[n, l] = -1/r_thrs[n]**2*np.log(thr/(r_thrs[n]**l))
    return alphas


def compute_overlap(n: int, n_prime: int, l: int, alphas: np.ndarray):
    """
    Computes the overelap between two radial basis functions phi_nl, phi_n'l.
    """

    def integrand(r): return r**2*basis.phi(r, l, alphas[n, l])*basis.phi(r, l, alphas[n_prime, l])
    return quad(integrand, 0.0, np.inf)[0]


def compute_s(n_max: int, l: int, alphas: np.ndarray) -> np.ndarray:
    """
    Computes the overlap matrix S (as in eq. (26) in Hilmanen et al 2020).
    """

    s = np.zeros((n_max, n_max))
    for n in range(n_max):
        for n_prime in range(n_max):
            tmp = compute_overlap(n, n_prime, l, alphas)
            s[n, n_prime] = tmp
    return s


def lowdin(s: np.ndarray):
    """
    Computes the Löwdin orthogonalization of the matrix s.
    """

    e, v = np.linalg.eigh(s)
    return np.dot(v/np.sqrt(e), v.T.conj())


def compute_beta(n_max: int, l: int, alphas: np.ndarray) -> np.ndarray:
    """
    Computes beta such that the corresponding basis is orthogonal (as in eq. (25) Hilmanen et al 2020).
    """

    s = compute_s(n_max, l, alphas)

    beta = lowdin(s)
    return beta


def precompute_parameters_of_radial_basis(n_max: int, l_max: int, r_min_thr: float, r_max_thr: float, dirs: Dict[str, Path]):
    """
    Precomputes the alphas and betas needed to define the basis functions (as in Hilmanen et al 2020).
    """

    thr = 10**(-3)
    r_thrs = np.zeros(n_max)

    r_thrs = np.linspace(r_min_thr, r_max_thr, n_max)

    if Debug:
        debug_out = dirs['ml'] / 'orbitals_to_power_spectra_debug.out'
        with open(debug_out, 'a') as file:
            file.write(f"r_thrs = {r_thrs}\n")

    alphas = compute_alphas(n_max, l_max, r_thrs, thr)

    betas = np.zeros((n_max, n_max, l_max+1))

    try:
        for l in range(l_max+1):
            betas[:, :, l] = compute_beta(n_max, l, alphas)
    except:
        raise ValueError(
            f"Failed to precompute the radial basis. You might want to try a larger r_min, e.g. r_min=0.5.")

    betas.tofile(dirs['betas'] / ('betas_' + '_'.join(str(x)
                                                      for x in [n_max, l_max, r_min_thr, r_max_thr]) + '.dat'))
    alphas.tofile(dirs['alphas'] / ('alphas_' + '_'.join(str(x)
                                                         for x in [n_max, l_max, r_min_thr, r_max_thr]) + '.dat'))

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
