import numpy as np
from datetime import datetime as dt

from modules import crys_to_cart



"""

Module to write the output of the unfolding and interpolation code.
Written by Riccardo De Gennaro 2019 (EPFL)

"""


"""
write_results calls write_bands and write_dos if the DOS was calculated

"""
def write_results(data):

    kvec = []
    for n in range(len(data.kvec)):
        kvec.append(crys_to_cart(data.kvec[n], data.bg, +1))

    write_bands(kvec, data.bands)

    if ( data.do_dos ):
        write_dos(data.dos)

    return


"""
write_bands prints the interpolated bands, in the QE format, in a file called
            'bands_interpolated.dat'.
            (see PP/src/bands.f90 around line 574 for the linearized path)

"""
def write_bands(kvec, bands):

    kx = [0.]
    for ik in range(1,len(kvec)):
        dxmod = np.linalg.norm(kvec[ik] - kvec[ik-1])
        if ( ik == 1 ):
            dxmod_save = dxmod

        if ( dxmod > 5*dxmod_save ):
            kx.append(kx[ik-1])

        elif ( dxmod > 1.e-4 ):
            kx.append(kx[ik-1] + dxmod)
            dxmod_save = dxmod

        else:
            kx.append(kx[ik-1] + dxmod)

    ofile = open('bands_interpolated.dat', 'w')
    ofile.write('# Written on %d-%d-%d at %d:%d:%02d'  %(dt.now().day,
        dt.now().month,dt.now().year,dt.now().hour,dt.now().minute,dt.now().second))
    for n in range(len(bands[0])):
        for ik in range(len(kvec)):
            ofile.write('\n%10.4f%10.4f' %(kx[ik],bands[ik][n]))
        ofile.write('\n')
    ofile.close()

    return
        
           
"""
write_dos prints the DOS in a file called 'dos_interpolated.dat', in a format (E , DOS(E))

"""
def write_dos(dos):
    with open('dos_interpolated.dat','w') as ofile:
        ofile.write('# Written on %d-%d-%d at %d:%d:%02d'  %(dt.now().day,
            dt.now().month,dt.now().year,dt.now().hour,dt.now().minute,dt.now().second))
        for n in range(len(dos)):
            ofile.write('\n%10.4f%12.6f' %(dos[n][0],dos[n][1]))
        ofile.write('\n')

    return
    
