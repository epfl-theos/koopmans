"""

Program for the unfolding and interpolation of bands for the Koopmans CP-code.
Written by Riccardo De Gennaro 2019 (EPFL)

For the moment the code works only with cubic, tetragonal and orthorombic systems.

"""


import sys
from time import time

from parser import Parse_Data
from interpolate import interpolate
from write_results import write_results


try:
    file_hr = sys.argv[1]
except:
    sys.exit('\nMissing argument with hamiltonian file -> EXIT\n')


start = time()
reset = time()

print()
print('UNFOLDING & INTERPOLATION\n')


"""
 1) Parse data:
    - calc parameters from the JSON file
    - other parameters from W90 output file
    - H(R) from the file passed as sys.argv[1]
    - WFs phases read from file wf_phases.dat
"""

sys_data = Parse_Data()
sys_data.parse_data(file_hr)

print('\tParsing input in:\t\t\t%.3f sec' % (time()-reset))
reset = time()


"""
 2) Core of the unfolding and interpolation code:
    - build the map |i> ---> |Rn>
    - calc interpolated (if needed) bands
    - calc DOS (if needed)
"""

interpolate(sys_data,reset)

reset = time()


"""
 3) Print out the results:
    - bands into 'bands_interpolated.dat' file
    - DOS into 'dos_interpolated.dat' file
"""

write_results(sys_data)

print('\tPriting output in:\t\t\t%.3f sec' % (time()-reset))
print()
print('\tTotal time:\t\t\t\t%.3f sec' %(time()-start))

print('\nALL DONE.\n')

