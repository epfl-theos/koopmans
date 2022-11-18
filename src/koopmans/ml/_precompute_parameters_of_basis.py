from pathlib import Path
from typing import Dict, Optional

import numpy as np
from scipy.integrate import quad

from ._basis_functions import phi


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
    Computes the overelap between two radial basis functions phi_nl, phi_n'l
    """

    def integrand(r): return r**2*phi(r, l, alphas[n, l])*phi(r, l, alphas[n_prime, l])
    return quad(integrand, 0.0, np.inf)[0]


def compute_s(n_max: int, l: int, alphas: np.ndarray) -> np.ndarray:
    """
    Computes the overlap matrix S (as in eq. (26) in Himanen et al 2020)
    """

    s = np.zeros((n_max, n_max))
    for n in range(n_max):
        for n_prime in range(n_max):
            tmp = compute_overlap(n, n_prime, l, alphas)
            s[n, n_prime] = tmp
    return s


def lowdin(s: np.ndarray):
    """
    Computes the Löwdin orthogonalization of the matrix s
    """

    e, v = np.linalg.eigh(s)
    return np.dot(v/np.sqrt(e), v.T.conj())


def compute_beta(n_max: int, l: int, alphas: np.ndarray) -> np.ndarray:
    """
    Computes beta such that the corresponding basis is orthogonal (as in eq. (25) Himanen et al 2020)
    """

    s = compute_s(n_max, l, alphas)

    beta = lowdin(s)
    return beta


def precompute_parameters_of_radial_basis(n_max: int, l_max: int, r_min_thr: float, r_max_thr: float,
                                          dirs: Optional[Dict[str, Path]] = None):
    """
    Precomputes the alphas and betas needed to define the basis functions (as in Himanen et al 2020)
    """

    thr = 10**(-3)
    r_thrs = np.linspace(r_min_thr, r_max_thr, n_max)
    alphas = compute_alphas(n_max, l_max, r_thrs, thr)
    betas = np.zeros((n_max, n_max, l_max+1))

    try:
        for l in range(l_max+1):
            betas[:, :, l] = compute_beta(n_max, l, alphas)
    except:
        raise ValueError(
            f"Failed to precompute the radial basis. You might want to try a larger r_min, e.g. r_min=1.0")

    if np.isnan(betas).any() or np.isnan(alphas).any():
        raise ValueError(
            f"Failed to precompute the radial basis. You might want to try a larger r_min, e.g. r_min=1.0")

    if dirs is not None:
        betas.tofile(dirs['betas'] / ('betas_' + '_'.join(str(x)
                                                          for x in [n_max, l_max, r_min_thr, r_max_thr]) + '.dat'))
        alphas.tofile(dirs['alphas'] / ('alphas_' + '_'.join(str(x)
                                                             for x in [n_max, l_max, r_min_thr, r_max_thr]) + '.dat'))
