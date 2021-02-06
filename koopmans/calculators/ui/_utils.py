"""
Utilities module for the UI calculator

Originally written by Riccardo De Gennaro as part of the standalone 'unfolding and interpolate' code
Integrated within python_KI by Edward Linscott Jan 2021

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

    hr_new = np.array(hr_new)

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


def MP_mesh(mp1, mp2, mp3):
    """
    Function for generating a regular Monkhorst-Pack mesh of dimension mp1*mp2*mp3
    """
    kmesh = []
    for i in range(mp1):
        for j in range(mp2):
            for k in range(mp3):
                kmesh.append(np.array((i / mp1, j / mp2, k / mp3), dtype=float))
    return kmesh


def generate_path(k_path):
    """
    Function for generating the path for the band structure interpolation.
    The input variable here (taken from the JSON file) must be a list in the
    Quantum Espresso 'crystal_b' format, i.e.:

                   [[ X1_1, X1_2, X1_3, N1 ],
                    [ X2_1, X2_2, X2_3, N2 ],
                    [  ***   ***   ***,  * ],
                    [  ***   ***   ***,  * ],
                    [ Xn_1, Xn_2, Xn_3, Nn ]]

    where Xi_j is the coordinate j of the point i (in crystal units) and Ni
    is the number of points along the line connecting the points i and i+1.
    NB: Nn must be always set to 1.
    """

    msg = '\n\'k_path\' in generate_path call must be a LIST (in crystal units) as follows:\n\n' \
        '\t\t[[ X1_1, X1_2, X1_3, N1 ],\n'\
        '\t\t [ X2_1, X2_2, X2_3, N2 ],\n'\
        '\t\t [  ***   ***   ***   * ],\n'\
        '\t\t [  ***   ***   ***   * ],\n'\
        '\t\t [ Xn_1, Xn_2, Xn_3,  1 ]]'

    for n in range(len(k_path)):
        try:
            kpt = np.array(k_path[n], dtype=float)
        except ValueError:
            raise ValueError(msg)

        assert len(kpt) == 4, msg

    kvec = []
    for n in range(len(k_path)):
        if (n == len(k_path) - 1):
            assert int(k_path[n][3]) == 1, 'msg'
            kvec.append(np.array((k_path[n][:3]), dtype=float))
            break
        else:
            npts = int(k_path[n][3])
            if npts == 0:
                npts = 1
            dist = np.array(k_path[n + 1][:3]) - np.array(k_path[n][:3])
            for m in range(npts):
                kvec.append(np.array(k_path[n][:3]) + dist * m / npts)

    return kvec
