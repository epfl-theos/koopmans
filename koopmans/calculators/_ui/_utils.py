"""
Utilities module for the UI calculator

Originally written by Riccardo De Gennaro as part of the standalone 'unfolding and interpolate' code
Integrated within koopmans by Edward Linscott Jan 2021

"""

from typing import List, TypeVar, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray


def crys_to_cart(vec: NDArray[np.float_], trmat: NDArray[np.float_], typ: int) -> NDArray[np.float_]:
    """
    Function to transform the numpy array vec (or a list/array of numpy arrays) from
    crystal to cartesian (in alat units), or viceversa, as it is done in QE: typ=+1
    for crystal-to-cartesian, typ=-1 for cartesian-to-crystal.
    For a real space vector trmat in input must be at if typ=+1 and bvec if typ=-1.
    For a k-space vector trmat in input must be bg if typ=+1 and avec if typ=-1.
    """

    if typ == +1:
        # crystal-to-cartesian conversion
        vec_tr = np.dot(vec, trmat)
    elif typ == -1:
        # cartesian-to-crystal conversion
        vec_tr = np.dot(vec, trmat.transpose())
    else:
        raise ValueError(f'typ = {typ} in crys_to_cart call must be either +1 or -1')

    return vec_tr


def extract_hr(hr: NDArray[np.complex_], rvect: NDArray[np.int_], nr1: int, nr2: int, nr3: int) -> NDArray[np.complex_]:
    """
    Function to select the Wannier Hamiltonian only on the primitive cell R-vectors.
    The Hamiltonian coming from a Wannier90 calculation with k-points is indeed
    defined on the Wigner-Seitz lattice vectors. In the case smooth=True, all the
    matrix elements corresponding to R-vectors exceeding the boundaries of the
    original supercell are ignored.
    """

    Rvec = latt_vect(nr1, nr2, nr3)
    rgrid = [nr1, nr2, nr3]
    hr_new = []
    ir = 0

    for R in Rvec:
        for ir, rvec in enumerate(rvect):
            if all(x < 1 for x in rvec / rgrid):
                rvec %= rgrid
                if all(rvec == R):
                    hr_new.append(hr[ir, :, :])
                    break

    assert len(hr_new) == np.prod(rgrid), f'Wrong number ({len(hr_new)}) of R-vectors in extract_hr'

    return np.array(hr_new, dtype=complex)


def latt_vect(nr1: int, nr2: int, nr3: int) -> NDArray[np.int_]:
    """
    Function for generating lattice vectors {R} of the primitive cell
    commensurate to the supercell. The R-vectors are given in crystal units.
    """

    Rvec = [[i, j, k] for i in range(nr1) for j in range(nr2) for k in range(nr3)]
    return np.array(Rvec)
