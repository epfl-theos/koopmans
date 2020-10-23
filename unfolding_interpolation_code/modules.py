import numpy as np
import sys


"""

Module with some functions used by the unfolding and interpolation code.
Written by Riccardo De Gennaro 2019 (EPFL)


It contains the following functions:
 - crys_to_cart
 - ang_to_bohr
 - latt_vect
 - MP_mesh
 - generate_path

"""


"""
Function to transform the vector coordinates from crystal to cartesian (in alat
units), or viceversa, as it is done in QE: typ=+1 for crystal-to-cartesian, typ=-1 
for cartesian-to-crystal.
For a real space vector trmat in input must be at if typ=+1 and bvec if typ=-1.
For a k-space vector trmat in input must be bg if typ=+1 and avec if typ=-1.
"""
def crys_to_cart(vec, trmat, typ):
    if ( typ == +1 ):                   # crystal-to-cartesian conversion
        vec_tr = np.dot(trmat.transpose(),vec)
    elif ( typ == -1 ):                 # cartesian-to-crystal conversion
        vec_tr = np.dot(trmat,vec)
    else:
        sys.exit('\ntyp in crys_to_cart call must be either +1 or -1 -> EXIT(%s)\n' %typ)
    return vec_tr


"""
Function to select the Wannier Hamiltonian only on the primitive cell R-vectors.
The Hamiltonian coming from a Wannier90 calculation with k-points is indeed
defined on the Wigner-Seitz lattice vectors.
"""
def order_hr(hr, rvect, nr1, nr2, nr3):
    Rvec = latt_vect(nr1, nr2, nr3)
    hr_new = np.zeros(( nr1*nr2*nr3, hr.shape[1], hr.shape[2] ))

    for ir in range(len(rvect)):
        rvect[ir][0] = rvect[ir][0]%nr1 
        rvect[ir][1] = rvect[ir][1]%nr2
        rvect[ir][2] = rvect[ir][2]%nr3

        for jr in range(nr1*nr2*nr3):
            if ( (rvect[ir] == Rvec[jr]).all() ):
                hr_new[jr,:,:] = hr[ir,:,:]

    return hr_new


"""
Function for Ang-Bohr conversions
typ = +1 conversion from Ang to Bohr
typ = -1 conversion from Bohr to Ang
"""
def ang_to_bohr(x, typ):
    if ( typ == +1 ):
        return x / 0.52917721067
    elif ( typ == -1 ):
        return x * 0.52917721067
    else:
        sys.exit('\ntyp in ang_to_bohr call must be either +1 or -1 -> EXIT(%s)\n' %typ)


"""
Function for generating lattice vectors {R} of the primitive cell
commensurate to the supercell. The R-vectors are given in crystal units.
"""
def latt_vect(nr1, nr2, nr3):
    Rvec = []
    for i in range(nr1):
        for j in range(nr2):
            for k in range(nr3):
                Rvec.append(np.array([i,j,k]))
    return Rvec
    

"""
Function for generating a regular Monkhorst-Pack mesh of dimension mp1*mp2*mp3
"""
def MP_mesh(mp1, mp2, mp3):
    kmesh = []
    for i in range(mp1):
        for j in range(mp2):
            for k in range(mp3):
                kmesh.append(np.array((i/mp1, j/mp2, k/mp3), dtype=float))
    return kmesh


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
Nn must be always set to 1.
"""
def generate_path(k_path):

    msg = '\n\'k_path\' in generate_path call must be a LIST (in crystal units) as follows:\n\n\
\t\t[[ X1_1, X1_2, X1_3, N1 ],\n\
\t\t [ X2_1, X2_2, X2_3, N2 ],\n\
\t\t [  ***   ***   ***   * ],\n\
\t\t [  ***   ***   ***   * ],\n\
\t\t [ Xn_1, Xn_2, Xn_3,  1 ]]\n'

    for n in range(len(k_path)):
        try:
            kpt = np.array(k_path[n], dtype=float)
        except ValueError:
            sys.exit(msg)

        if ( len(kpt) != 4 ):    sys.exit(msg)

    kvec = []
    for n in range(len(k_path)):
        if ( n == len(k_path)-1 ):
            if ( int(k_path[n][3]) != 1 ):    sys.exit(msg)
            kvec.append(np.array((k_path[n][:3]), dtype=float))
            break
        else:
            npts = int(k_path[n][3])
            dist = np.array(k_path[n+1][:3]) - np.array(k_path[n][:3])
            for m in range(k_path[n][3]):
                kvec.append(np.array(k_path[n][:3]) + dist * m / npts)

    return kvec
