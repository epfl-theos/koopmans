"""Utilities for the UI calculator."""

import numpy as np
from numpy.typing import NDArray


def crys_to_cart(vec: NDArray[np.float64], trmat: NDArray[np.float64], typ: int) -> NDArray[np.float64]:
    """Transform the numpy array vec from crystal to cartesian (in alat units), or vice versa.

    Follows the transformation as done in Quantum ESPRESSO:
    typ=+1 for crystal-to-cartesian, typ=-1 for cartesian-to-crystal.
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
        raise ValueError(f'`typ = {typ}` in `crys_to_cart` must be either +1 or -1')

    return vec_tr


def extract_hr(hr: NDArray[np.complex128], rvect: NDArray[np.int_],
               nr1: int, nr2: int, nr3: int) -> NDArray[np.complex128]:
    """Select the Wannier Hamiltonian only on the primitive cell R-vectors.

    The Hamiltonian coming from a Wannier90 calculation with k-points is indeed
    defined on the Wigner-Seitz lattice vectors. In the case smooth=True, all the
    matrix elements corresponding to R-vectors exceeding the boundaries of the
    original supercell are ignored.
    """
    Rvec = latt_vect(nr1, nr2, nr3)
    rgrid = [nr1, nr2, nr3]
    hr_new: list[NDArray[np.complex128]] = []
    ir = 0

    for R in Rvec:
        for ir, rvec in enumerate(rvect):
            assert isinstance(rvec, np.ndarray)
            if all(x < 1 for x in rvec / rgrid):
                rvec %= rgrid
                if all(rvec == R):
                    hr_new.append(hr[ir, :, :])
                    break

    assert len(hr_new) == np.prod(rgrid), f'Wrong number ({len(hr_new)}) of R-vectors in `extract_hr`'

    return np.array(hr_new, dtype=np.complex128)


def latt_vect(nr1: int, nr2: int, nr3: int) -> NDArray[np.int_]:
    """Generate lattice vectors {R} of the primitive cell commensurate to the supercell.

    The R-vectors are given in crystal units.
    """
    Rvec = [[i, j, k] for i in range(nr1) for j in range(nr2) for k in range(nr3)]
    return np.array(Rvec)
