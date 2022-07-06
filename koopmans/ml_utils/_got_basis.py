import numpy as np


def phi(r: np.ndarray, l: int, alpha: float) -> np.ndarray:
    """Corresponding to phi_nl eq. (24) in Hilmanen et al 2020."""
    return r**l*np.exp(-alpha*r**2)


def g(r: np.ndarray, n: int, n_max: int, l, betas: np.ndarray, alphas: np.ndarray) -> np.ndarray:
    """Corresponding to """
    return sum(betas[n_prime, n, l]*phi(r, l, alphas[n_prime, l]) for n_prime in range(n_max))


def radial_basis_function(r: np.ndarray, n: int, n_max: int, l: int, betas: np.ndarray, alphas: np.ndarray):
    return g(r, n, n_max, l, betas, alphas)
