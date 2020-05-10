import numpy as np
import sys
from time import time

from modules import latt_vect, crys_to_cart

"""

Module carrying out the main work of the unfolding and interpolation code.
Written by Riccardo De Gennaro 2019 (EPFL)

"""


"""
interpolate is the main program in this module and it calls consecutively
            the three independent functions:
            - map_wannier
            - calc_bands
            - calc_dos

"""
def interpolate(data,start_time):
    
    # Step 1: map the WFs
    map_wannier(data)
    print('\tBuilding the map |i> --> |Rn> in:\t%.3f sec' % (time()-start_time))
    reset = time()

    # Step 2: calculate the electronic bands along k_path
    calc_bands(data)
    print('\tCalculating bands in:\t\t\t%.3f sec' % (time()-reset))
    reset = time()

    # Step 3 (optional) : calculate the density-of-states
    if ( data.do_dos ):
        calc_dos(data)
        print('\tCalculating DOS in:\t\t\t%.3f sec' % (time()-reset))

    return


"""
map_wannier builds the map |i> --> |Rn> between the WFs in the SC and in the PC.
        
       IMPORTANT: when running this function we conceptually pass from the supercell
                  to a primitive cell description and so all the related quantities
                  are adapted (at, bg, alat).

"""
def map_wannier(data):

    # shift the WFs within the SC
    data.centers = data.centers - np.floor(data.centers)

    # generate R-vectors of the PC
    Rvec = latt_vect(data.nr1, data.nr2, data.nr3)
    
    # define number of WFs in the PC
    data.num_wann_pc = int(data.num_wann / len(Rvec))
    
    centers = []
    spreads = []
    index = []
    num_wann_pc = 0

    # converting at, bg and alat from the SC to the PC
    data.alat_sc = data.alat
    data.alat = data.alat / data.nr1
    #
    data.at_sc = np.zeros((3,3))
    data.at_sc[0] = data.at[0]    
    data.at_sc[1] = data.at[1]    
    data.at_sc[2] = data.at[2]
    data.at[0] = data.at[0]                          # a_1 --> a_1
    data.at[1] = data.at[1] * data.nr1 / data.nr2    # a_2 --> a_2 * a / b
    data.at[2] = data.at[2] * data.nr1 / data.nr3    # a_3 --> a_3 * a / c
    #
    data.bg_sc = np.zeros((3,3))
    data.bg_sc[0] = data.bg[0]    
    data.bg_sc[1] = data.bg[1]    
    data.bg_sc[2] = data.bg[2]
    data.bg[0] = data.bg[0]                          # b_1 --> b_1
    data.bg[1] = data.bg[1] * data.nr2 / data.nr1    # b_2 --> b_2 * b / a
    data.bg[2] = data.bg[2] * data.nr3 / data.nr1    # b_3 --> b_3 * c / a
 
    # here we identify the WFs within the R=0 cell
    for n in range(data.num_wann):
        # converting centers from crystal units of SC to crystal units of PC
        data.centers[n,0] = data.centers[n,0] * data.nr1      
        data.centers[n,1] = data.centers[n,1] * data.nr2
        data.centers[n,2] = data.centers[n,2] * data.nr3
        
        if ( (data.centers[n,0] - 1 < 1.e-3) and \
             (data.centers[n,1] - 1 < 1.e-3) and \
             (data.centers[n,2] - 1 < 1.e-3) ):
            centers.append(data.centers[n])
            spreads.append(data.spreads[n])
            index.append(n)
            num_wann_pc += 1

    # check on the WFs found in the R=0 cell
    if ( num_wann_pc != data.num_wann_pc ):
        sys.exit('\nDid not find the right number of WFs in the R=0 cell -> EXIT(%s)\n' \
                                                                    %num_wann_pc)
    
    # here we identify with |Rn> the WFs in the rest of the SC
    # the WFs are now ordered as (R0,1),(R0,2),...,(R0,n),(R1,1),...
    for rvect in Rvec[1:]:
        count = 0
        for m in range(data.num_wann_pc):
            for n in range(data.num_wann):
                wf_dist = data.centers[n] - centers[m]
                if ( (np.linalg.norm(wf_dist - rvect) < 1.e-3) and \
                     (abs(data.spreads[n] - spreads[m] < 1.e-3)) ):
                    centers.append(data.centers[n])
                    spreads.append(data.spreads[n])
                    index.append(n)
                    count += 1
#                    break
        
        if ( count != data.num_wann_pc ):
            sys.exit('\nDid not find the right number of WFs in the %s cell -> EXIT(%s)\n' \
                                                                        %(rvect,count))

    # redefine phases and hr in order to follow the new order of WFs
    phases = []
    hr = np.zeros((data.num_wann,data.num_wann), dtype=float)
    for n in range(data.num_wann):
        phases.append(data.phases[index[n]])
        for m in range(data.num_wann):
            hr[m,n] = data.hr[index[m],index[n]]

    data.Rvec = Rvec
    data.centers = centers
    data.spreads = spreads
    data.hr = hr
    data.phases = phases

    return


"""
calc_bands interpolates the k-space hamiltonian along the input path, by Fourier
           transforming the Wannier hamiltonian H(R). The function generates two
           new attirubtes for the instance data:
           - data.hk containing H(k) for any k-vector in the input path
           - data.bands containing the interpolated electronic energies

"""
def calc_bands(data):

    # renormalize H(R) on the WF phases
    hr = np.array(data.hr, dtype=complex)
    for m in range(data.num_wann):
        for n in range(data.num_wann):
            hr[m,n] = data.phases[m].conjugate() * hr[m,n] * data.phases[n]

    # here we build the interpolated H(k)
    hk = np.zeros((len(data.kvec),data.num_wann_pc,data.num_wann_pc), dtype=complex)
    bands = []
    for ik in range(len(data.kvec)):
        kpt = data.kvec[ik]
        for m in range(data.num_wann_pc):
            for n in range(data.num_wann_pc):
                for ir in range(len(data.Rvec)):

                    nn = ir * data.num_wann_pc + n
                    phase = correct_phase(data.centers[m], data.centers[nn], data.Rvec[ir], kpt, data)

                    hk[ik,m,n] = hk[ik,m,n] + np.exp(1j * 2 * np.pi * np.dot(kpt,data.Rvec[ir])) * \
                                              phase * hr[m,nn]

        bands.append(np.linalg.eigvalsh(hk[ik]))

    data.hk = hk
    data.bands = bands

    return


"""
calc_dos calculates the density of states using a gaussian smearing. The DOS is saved
         as a list [ [E1, DOS(E1)], [E2, DOS[E2]], ... , [En, DOS(En)] ]

"""
def calc_dos(data):

    dE = (data.Emax - data.Emin) / data.nstep
    ene = data.Emin
    eigvl = np.array(data.bands, dtype=float).reshape(data.num_wann_pc*len(data.kvec))
    dos = [ [ene, sum(np.exp( - ((ene - eigvl) / data.degauss)**2 )) / (data.degauss * np.pi**.5)] ]

    for n in range(data.nstep):
        ene = ene + dE
        dos.append([ ene, sum(np.exp( - ((ene - eigvl) / data.degauss)**2 )) / \
                                       (data.degauss * np.pi**.5)])
    
    data.dos = dos

    return


"""
correct_phase calculate the correct phase factor to put in the Fourier transform
              to get the interpolated k-space hamiltonian. The correction consists
              of finding the right distance, i.e. the right R-vector, considering
              also the BVK boundary conditions.
              if ws_distance=True, the function accounts also for the intracell 
              distance between Wannier functions, otherwise only the intercell
              distances are considered.
       
       IMPORTANT: the input vectors (center_ref, center, Rvec and kvec) must all 
                  be in cartesian units (normalized or not by alat does not make
                  any difference), otherwise the distances are not properly evaluated.      

"""
def correct_phase(center_ref, center, rvect, kvect, data):

    if ( data.ws_distance ):
        wf_dist = crys_to_cart(center - center_ref, data.at, +1)
    else:
        wf_dist = crys_to_cart(rvect, data.at, +1)

    dist_min = np.linalg.norm(wf_dist)
    Tlist = []

    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                tvect = np.array([ i*data.nr1, j*data.nr2, k*data.nr3 ])
                Tvec = crys_to_cart(tvect, data.at, +1)
                dist = np.linalg.norm(wf_dist + Tvec)

                if ( abs(dist - dist_min) < 1.e-3 ):
                    #
                    # an equidistant replica is found
                    #
                    Tlist.append(tvect)

                elif ( dist < dist_min ):
                    #
                    # a new nearest replica is found
                    # reset dist_min and reinitialize Tlist
                    #
                    dist_min = dist
                    Tlist = [tvect]

                else:
                    #
                    # this replica is rejected
                    #
                    continue

    phase = sum(np.exp(1j * 2 * np.pi * np.dot(Tlist,kvect))) / len(Tlist)

    return phase


