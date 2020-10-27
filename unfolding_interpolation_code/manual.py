import argparse
import sys

"""

Module parsing the command line arguments and creating
the manual for the program 

"""


parser = argparse.ArgumentParser( formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''
#####################################################
### Manual for the UNFOLDING & INTERPOLATION code ###
#####################################################

The code must be run from the directory containing:
  * the file 'input_ui.json'
  * the Wannier90 output file 'seedname.wout'
  * (optional) the file 'wf_phases.dat'

Input keywords description
--------------------------
seedname        :        seedname as in the Wannier90 calculation
alat_sc         :        supercell lattice parameter (in BOHR)
nr1             :        units of repetition of the PC in the SC along a_1
nr2             :        units of repetition of the PC in the SC along a_2
nr3             :        units of repetition of the PC in the SC along a_3
w90_calc        :        type of Wannier90 calculation: 'pc' or 'sc'
do_map          :        if True, it maps the WFs from the SC to the PC
use_ws_distance :        if True, uses the Wigner-Seitz distance between WFs centers
k_path          :        k_path for bands interpolation (in crystal units)
smooth_int      :        if True, the smooth interpolation method is used
file_hr_coarse  :        
file_hr_smooth  :
do_dos          :        if True, the density-of-states is calculated
degauss         :        gaussian broadening (in eV)
nstep           :        number of steps for DOS plot
Emin            :        min energy (in eV) for DOS plot
Emax            :        max energy (in eV) for DOS plot
                                         ''',
                                  epilog=' ' )

list_of_keywords = [ 'seedname', 'alat_sc', 'nr1', 'nr2', 'nr3', 'w90_calc', 'do_map', \
                     'use_ws_distance', 'k_path', 'smooth_int', 'file_hr_coarse', 'file_hr_smooth', \
                     'do_dos', 'degauss', 'nstep', 'Emin', 'Emax', 'all' ]

# optional arguments
parser.add_argument( "-i", "--info",
                     help="followed by one (or multiple) input keyword provides a detailed " +
                          "description - type 'all' for a full description of all the input keywords",
                     nargs="+",
                     choices=list_of_keywords,
                     metavar='' )

# positional arguments
parser.add_argument( "file_hr",
                     help="file with the Koopmans or Wannier90 hamiltonian to interpolate", 
                     nargs="?" )


args = parser.parse_args()

keyword_descriptor = { 
       'seedname'        : 'the seedname must be equal to the seedname used in the previous ' +
                           'Wannier90 calculation. The code will look for a file called seedname.wout.',
       'alat_sc'         : 'the lattice parameter (in BOHR) of the supercell, as celldm(1) in QE. ' +
                           'NB: it is important to put the supercell (and not the primitive cell) ' +
                           'lattice parameter, otherwise the result will be wrong.',
       'nr1'             : 'units of repetition of the primitive cell within the supercell along ' +
                           'the direction of a_1. Equivalently this has to match the number of ' +
                           'k-points in the Monkhorst-Pack mesh along the direction of b_1.',
       'nr2'             : 'units of repetition of the primitive cell within the supercell along ' +
                           'the direction of a_2. Equivalently this has to match the number of ' +
                           'k-points in the Monkhorst-Pack mesh along the direction of b_2.',
       'nr3'             : 'units of repetition of the primitive cell within the supercell along ' +
                           'the direction of a_3. Equivalently this has to match the number of ' +
                           'k-points in the Monkhorst-Pack mesh along the direction of b_3.',
       'w90_calc'        : 'this keyword specifies the type of PW/Wannier90 calculation preceding ' +
                           'the koopmans calculation. If the latter is done in a supercell at Gamma ' +
                           'then w90_calc must be equal to \'sc\', otherwise if it comes from a ' +
                           'calculation with k-points it must be equal to \'pc\'.\n' +
                           'The default value is \'pc\'.',
       'do_map'          : 'if True, it realizes the map |m> --> |Rn>, that connects the Wannier ' +
                           'functions in the supercell to those in the primitive cell. This is ' +
                           'basically the unfolding procedure. It can be activated only if ' +
                           'w90_calc=\'sc\'.\n' +
                           'The default value is False.',
       'use_ws_distance' : 'if True, the real Wigner-Seitz distance between the Wannier functions ' +
                           'centers is considered as in the Wannier90 code. In particular, this ' +
                           'accounts for the periodic boundary conditions and it is crucial for ' +
                           'a good interpolation when using coarse MP meshes or, equivalently, ' +
                           'small supercells.\n' +
                           'The default value is True.',
       'k_path'          : 'path in the Brillouin zone for the band structure. The logic is the ' +
                           'crystal_b units in QE: the path must be defined by providing the ' +
                           'initial and final point crystal coordinates of each line, followed by ' +
                           'the number of points along the line.\n\n' +
                           'Example: path Gamma-X-Gamma with 10 points along each line\n\n' +
                           '\tk_path :\t[ [ 0.000, 0.000, 0.000, 10 ],\n' +
                           '\t        \t  [ 0.500, 0.000, 0.000, 10 ],\n' +
                           '\t        \t  [ 0.000, 0.000, 0.000,  1 ] ]\n' +
                           '\nThe band structure is then written into a file called ' +
                           '\'bands_interpolated.dat\'.',
       'smooth_int'      : 'if True, the smooth interpolation method is used. This consists of ' +
                           'removing the DFT part of the hamiltonian from the full Koopmans ' +
                           'hamiltonian and adding the DFT hamiltonian from a calculation with ' +
                           'a denser k-points mesh. This works only for a non self-consistent ' +
                           'Koopmans calculation using Wannier since, to be consistent, all the ' +
                           'hamiltonians must be in the same gauge, i.e. the Wannier gauge.\n' +
                           'The default value is False.',
       'file_hr_coarse'  : '',
       'file_hr_smooth'  : '',
       'do_dos'          : 'if True, the density-of-states is interpolated along the input k_path. ' +
                           'The DOS is written to a file called \'dos_interpolated.dat\'.\n' +
                           'The default value is False.',
       'degauss'         : 'gaussian broadening (in eV) for the DOS interpolation, as in QE.\n' +
                           'The default value is 0.05 eV.',
       'nstep'           : 'number of steps for the plot of the interpolated DOS.\n' +
                           'The default value is 1000.',
       'Emin'            : 'minimum energy for the plot of the interpolated DOS.\n' +
                           'The default value is the lowest calculated energy.',
       'Emax'            : 'maximum energy for the plot of the interpolated DOS.\n' +
                           'The default value is the highest calculated energy.'
                     }

if ( args.info ):
    #
    for arg in args.info:
        if ( arg == 'all' ):
            args.info = list_of_keywords
            args.info.remove('all')
    #
    for key in args.info:
        print( "\n%s :\n%s\n" %(key,keyword_descriptor[key]) )
    #
    sys.exit('\n')

if ( not args.file_hr ):
    sys.exit('\nMissing the positional argument file_hr -> EXIT\n')
