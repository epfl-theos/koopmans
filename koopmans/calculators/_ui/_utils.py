"""
Utilities module for the UI calculator

Originally written by Riccardo De Gennaro as part of the standalone 'unfolding and interpolate' code
Integrated within koopmans by Edward Linscott Jan 2021

"""

import numpy as np


def crys_to_cart(vec, trmat, typ):
    """
    Function to transform the vector coordinates from crystal to cartesian (in alat
    units), or viceversa, as it is done in QE: typ=+1 for crystal-to-cartesian, typ=-1
    for cartesian-to-crystal.
    For a real space vector trmat in input must be at if typ=+1 and bvec if typ=-1.
    For a k-space vector trmat in input must be bg if typ=+1 and avec if typ=-1.
    """
    if typ == +1:
        # crystal-to-cartesian conversion
        vec_tr = np.dot(trmat.transpose(), vec)
    elif typ == -1:
        # cartesian-to-crystal conversion
        vec_tr = np.dot(trmat, vec)
    else:
        raise ValueError(f'typ = {typ} in crys_to_cart call must be either +1 or -1')
    return vec_tr


def extract_hr(hr, rvect, nr1, nr2, nr3):
    """
    Function to select the Wannier Hamiltonian only on the primitive cell R-vectors.
    The Hamiltonian coming from a Wannier90 calculation with k-points is indeed
    defined on the Wigner-Seitz lattice vectors. In the case smooth=True, all the
    matrix elements corresponding to R-vectors exceeding the boundaries of the
    original supercell are ignored.
    """
    Rvecs = latt_vect(nr1, nr2, nr3)
    hr_new = []

    for R in Rvecs:
        for ir in range(len(rvect)):

            if (abs(float(rvect[ir][0]) / nr1) >= 1) or (abs(float(rvect[ir][1]) / nr2) >= 1) \
                    or (abs(float(rvect[ir][2]) / nr3) >= 1):
                continue

            rvect[ir][0] = rvect[ir][0] % nr1
            rvect[ir][1] = rvect[ir][1] % nr2
            rvect[ir][2] = rvect[ir][2] % nr3

            if ((rvect[ir] == R).all()):
                hr_new.append(hr[ir, :, :])
                break

    assert len(hr_new) == nr1 * nr2 * nr3, f'Wrong number ({len(hr_new)}) of R-vectors in extract_hr'

    hr_new = np.array(hr_new, dtype=complex)

    return hr_new


def latt_vect(nr1, nr2, nr3):
    """
    Function for generating lattice vectors {R} of the primitive cell
    commensurate to the supercell. The R-vectors are given in crystal units.
    """
    Rvec = []
    for i in range(nr1):
        for j in range(nr2):
            for k in range(nr3):
                Rvec.append(np.array([i, j, k]))
    return Rvec
