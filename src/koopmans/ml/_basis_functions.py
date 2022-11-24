from typing import Union

import numpy as np
from scipy.special import sph_harm


def phi(r: np.ndarray, l: int, alpha: float) -> np.ndarray:
    """Implementation of phi_nl from eq. (24) in Himanen et al 2020."""
    return r**l * np.exp(-alpha * r**2)


def g(r: np.ndarray, n: int, n_max: int, l, betas: np.ndarray, alphas: np.ndarray) -> np.ndarray:
    """Implementation of g from eq. (23) in Himanen et al 2020."""
    g_vec = np.zeros_like(r)
    g_vec[:, :, :] = sum(betas[n_prime, n, l] * phi(r, l, alphas[n_prime, l]) for n_prime in range(n_max))
    return g_vec


def real_spherical_harmonics(theta: np.ndarray, phi: np.ndarray, l: int, m: int) -> np.ndarray:
    """Implementation of Y_lm from eq. (20) in Himanen et al 2020."""
    Y = sph_harm(abs(m), l, phi, theta)
    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real
    return Y.real
