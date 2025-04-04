import numpy as np
from scipy.integrate import quad

from ._basis_functions import phi


def compute_alphas(n_max: int, l_max: int, r_thrs: np.ndarray, thr: float):
    """Compute the decay-coefficients alpha_nl.

    Does this by demanding that for each r_thrs[n], the corresponding phi_nl decays to threshold value
    thr at a cutoff radius of r_thr.
    """
    alphas = np.zeros((n_max, l_max + 1))
    for n in range(n_max):
        for l in range(l_max + 1):
            alphas[n, l] = -1 / r_thrs[n]**2 * np.log(thr / (r_thrs[n]**l))
    return alphas


def compute_overlap(n: int, n_prime: int, l: int, alphas: np.ndarray):
    """Compute the overelap between two radial basis functions phi_nl, phi_n'l."""
    def integrand(r):
        return r**2 * phi(r, l, alphas[n, l]) * phi(r, l, alphas[n_prime, l])
    return quad(integrand, 0.0, np.inf)[0]


def compute_s(n_max: int, l: int, alphas: np.ndarray) -> np.ndarray:
    """Compute the overlap matrix S (as in eq. (26) in Himanen et al 2020)."""
    s = np.zeros((n_max, n_max))
    for n in range(n_max):
        for n_prime in range(n_max):
            tmp = compute_overlap(n, n_prime, l, alphas)
            s[n, n_prime] = tmp
    return s


def lowdin(s: np.ndarray):
    """Compute the Löwdin orthogonalization of the matrix s."""
    e, v = np.linalg.eigh(s)
    return np.dot(v / np.sqrt(e), v.T.conj())


def compute_beta(n_max: int, l: int, alphas: np.ndarray) -> np.ndarray:
    """Compute beta such that the corresponding basis is orthogonal (as in eq. (25) Himanen et al 2020)."""
    s = compute_s(n_max, l, alphas)

    beta = lowdin(s)
    return beta


def precompute_parameters_of_radial_basis(n_max: int, l_max: int, r_min_thr: float, r_max_thr: float):
    """Precompute the alphas and betas needed to define the basis functions (as in Himanen et al 2020)."""
    thr = 10**(-3)
    r_thrs = np.linspace(r_min_thr, r_max_thr, n_max)
    alphas = compute_alphas(n_max, l_max, r_thrs, thr)
    betas = np.zeros((n_max, n_max, l_max + 1))

    for l in range(l_max + 1):
        betas[:, :, l] = compute_beta(n_max, l, alphas)

    if np.isnan(betas).any() or np.isnan(alphas).any():
        raise ValueError(
            "Failed to precompute the radial basis. You might want to try a larger `r_min`, e.g. `r_min = 1.0`")

    return alphas, betas
