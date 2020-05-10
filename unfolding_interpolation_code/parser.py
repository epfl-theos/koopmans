import numpy as np
import sys,os
import json

from modules import *


"""

Module for parsing the input data used by the unfolding and interpolation program.
Written by Riccardo De Gennaro 2019 (EPFL)

"""


class Parse_Data():


    """
    parse_data simply gets all the parsable data:
               - JSON file 
               - Wannier90 output
               - hamiltonian H(R)
               - Wannier functions phases

    """
    def parse_data(self, file_hr):
        
        self.parse_json()
        self.parse_w90()
        self.parse_hr(file_hr)
        self.parse_phases()

        return


    """
    parse_json reads the input_ui.json and eventully assigns a default
               value (if present) to the missing parameters

    seedname        : as in W90 input
    alat            : lattice parameter (in BOHR units). Like celldm(1) in QE it refers to a_1
    nr1             : SC dimension along PC primitive vector a_1
    nr2             : SC dimension along PC primitive vector a_2
    nr3             : SC dimension along PC primitive vector a_3
    use_ws_distance : if True the Wigner-Seitz distance between WF centers is used
    k_path          : k_path for bands interpolation (in crystal units)
    smooth_int      : if True, smooth interpolation method is used
    file_hr_pw      : look at documentation (needed when smooth_int=True)
    file_hr_w90     : look at documentation (needed when smooth_int=True)
    do_dos          : if True, the density-of-states is calculated    
    degauss         : gaussian broadening in eV (as in QE except for the units)
    nstep           : number of steps for DOS plot
    Emin            : min energy (in eV) for DOS plot
    Emax            : max energy (in eV) for DOS plot
    
    """
    def parse_json(self):

        with open('input_ui.json','r') as json_file:
            json_data = json.load(json_file)

        # list of available input arguments
        self.list_of_args = ['seedname', 'alat', 'nr1', 'nr2', 'nr3', 'use_ws_distance',\
                             'k_path', 'smooth_int', 'file_hr_pw', 'file_hr_w90',\
                             'do_dos', 'degauss', 'nstep', 'Emin', 'Emax']
        
        if 'seedname' not in json_data.keys():    sys.exit('\nMissing \'seedname\' in input -> EXIT\n')
        if 'alat' not in json_data.keys():        sys.exit('\nMissing \'alat\' in input -> EXIT\n')
        if 'nr1' not in json_data.keys():         sys.exit('\nMissing \'nr1\' in input -> EXIT\n')
        if 'nr2' not in json_data.keys():         sys.exit('\nMissing \'nr2\' in input -> EXIT\n')
        if 'nr3' not in json_data.keys():         sys.exit('\nMissing \'nr3\' in input -> EXIT\n')
        
        ### Reading json dictionary and defining defaults ###
        self.seedname = json_data.get("seedname", None)
        self.alat = json_data.get("alat", None)
        self.nr1 = json_data.get("nr1", None)
        self.nr2 = json_data.get("nr2", None)
        self.nr3 = json_data.get("nr3", None)
        self.use_ws_distance = json_data.get("use_ws_distance", True)
        self.k_path = json_data.get("k_path", None)
        self.smooth_int = json_data.get("smooth_int", False)
        self.file_hr_pw = json_data.get("file_hr_pw", None)
        self.file_hr_w90 = json_data.get("file_hr_w90", None)
        self.do_dos = json_data.get("do_dos", False)
        self.degauss = json_data.get("degauss", 0.05)
        self.nstep = json_data.get("nstep", 1000)
        self.Emin = json_data.get("Emin", -10)
        self.Emax = json_data.get("Emax", 10)
        ####################################################

        # checks on the input arguments
        if ( type(self.seedname) is not str ):
            sys.exit('\n\'seedname\' must be a string -> EXIT(%s)' %type(self.seedname))
        if ( type(self.alat) is not float and type(self.alat) is not int ):
            sys.exit('\n\'alat\' must be a float or a int -> EXIT(%s)' %type(self.alat))
        if ( type(self.nr1) is not int ):
            sys.exit('\n\'nr1\' must be a int -> EXIT(%s)' %type(self.nr1))
        if ( type(self.nr2) is not int ):
            sys.exit('\n\'nr2\' must be a int -> EXIT(%s)' %type(self.nr2))
        if ( type(self.nr3) is not int ):
            sys.exit('\n\'nr3\' must be a int -> EXIT(%s)' %type(self.nr3))
        if ( type(self.use_ws_distance) is not bool ):
            sys.exit('\n\'use_ws_distance\' must be a bool -> EXIT(%s)' %type(self.use_ws_distance))
        if ( type(self.k_path) is not list ):
            sys.exit('\nk_path must be a list -> EXIT(%s)' %type(self.k_path))
        if ( type(self.smooth_int) is not bool ):
            sys.exit('\n\'smooth_int\' must be a bool -> EXIT(%s)' %type(self.smooth_int))
        if ( type(self.file_hr_pw) is not str ):
            sys.exit('\n\'file_hr_pw\' must be a string -> EXIT(%s)' %type(self.file_hr_pw))
        if ( type(self.file_hr_w90) is not str ):
            sys.exit('\n\'file_hr_w90\' must be a string -> EXIT(%s)' %type(self.file_hr_w90))
        if ( type(self.do_dos) is not bool ):
            sys.exit('\n\'do_dos\' must be a bool -> EXIT(%s)' %type(self.do_dos))
        if ( type(self.degauss) is not float and type(self.degauss) is not int ):
            sys.exit('\n\'degauss\' must be a float or a int -> EXIT(%s)' %type(degauss))
        if ( type(self.nstep) is not float and type(self.nstep) is not int ):
            sys.exit('\n\'nstep\' must be a float or a int -> EXIT(%s)' %type(self.nstep))
        if ( type(self.Emin) is not float and type(self.Emin) is not int ):
            sys.exit('\n\'Emin\' must be a float or a int -> EXIT(%s)' %type(self.Emin))
        if ( type(self.Emax) is not float and type(self.Emax) is not int ):
            sys.exit('\n\'Emax\' must be a float or a int -> EXIT(%s)' %type(self.Emax))
           
        for item in json_data.keys():
            if ( item not in self.list_of_args ):
                print('\nWARNING: argument \'%s\' is unknown!\n' %item)
       
        self.alat = ang_to_bohr(self.alat, -1)    # conversion to Ang
  
        if ( self.k_path == None ):
            self.kvec = MP_mesh(self.nr1,self.nr2,self.nr3)
            print('\nWARNING: \'k_path\' missing in input, the energies are calculated on a\
 commensurate Monkhorst-Pack mesh\n')
        else:
            self.kvec = generate_path(self.k_path)
        
        if ( self.smooth_int ):
            if ( 'file_hr_pw' not in json_data.keys() ):
                sys.exit('\nMissing \'file_hr_pw\' for smooth interpolation -> EXIT\n')
            if ( 'file_hr_w90' not in json_data.keys() ):
                sys.exit('\nMissing \'file_hr_w90\' for smooth interpolation -> EXIT\n')
        
        if ( self.do_dos ):
            if 'degauss' not in json_data.keys():
                print('\nWARNING: \'degauss\' missing in input, using default value\n')
            if 'nstep' not in json_data.keys():
                print('\nWARNING: \'nstep\' missing in input, using default value\n')
            if 'Emin' not in json_data.keys():
                print('\nWARNING: \'Emin\' missing in input, using default value\n')
            if 'Emax' not in json_data.keys():
                print('\nWARNING: \'Emax\' missing in input, using default value\n')

        return


    """
    parse_w90 gets from the W90 output the lattice vectors, and the centers and spreads
              of the Wannier functions. 

    at      : basis vectors of direct lattice (in alat units)
    bg      : basis vectors of reciprocal lattice (in 2pi/alat units)
    centers : centers of WFs (in crystal units)
    spreads : spreads of WFs (in Ang^2)
    
    """
    def parse_w90(self):

        with open(self.seedname+'.wout','r') as ifile:
            lines = ifile.readlines()

        self.at = []
        self.centers = []
        self.spreads = []
        count = 0

        for line in lines:
            if ( 'Number of Wannier Functions' in line ):
                self.num_wann = int(line.split()[6])
            if ( ' a_1 ' in line ):
                self.at.append(np.array(line.split()[1:], dtype=float))
            if ( ' a_2 ' in line ):
                self.at.append(np.array(line.split()[1:], dtype=float))
            if ( ' a_3 ' in line ):
                self.at.append(np.array(line.split()[1:], dtype=float))
            if ( count > 0 and count <= self.num_wann ):
                start = line.find('(')
                end = line.find(')')
                self.centers.append(np.array(line[start+1:end].replace(',',' ').split(), \
                                                                           dtype=float))
                #
                self.spreads.append(float(line.split()[-1]))
                #
                count += 1
            if ( 'Final State' in line ):
                count += 1

        self.at = np.array(self.at, dtype=float).reshape(3,3) / self.alat
        self.bg = np.linalg.inv(self.at).transpose()
        
        for n in range(self.num_wann):
            self.centers[n] = self.centers[n] / self.alat
            self.centers[n] = crys_to_cart(self.centers[n], self.bg, -1)
       
        return
    
    
    """
    parse_hr reads the hamiltonian file passed as sys.argv[1] and it sets it up
             as attribute self.hr
             
    There are 3 possible types of file:
      - W90 file normally called seedname_hr.dat
      - KC_occ file normally called hamiltonian1.xml
      - KC_emp file normally called hamiltonian_emp.dat

    NB: KC_emp must be called 'hamiltonian_emp.dat' otherwise the code may crash
        or misread the matrix elements. If the file name is different the code 
        should be updated.

    """
    def parse_hr(self, file_hr):

        with open(file_hr,'r') as ifile:
            lines = ifile.readlines()

        if ( 'written on' in lines[0] ):                        hr_type = 'w90'
        elif ( 'xml version' in lines[0]):                      hr_type = 'kc_occ'
        elif ( file_hr[-19:] == 'hamiltonian_emp.dat' ):    hr_type = 'kc_emp'
        else:		sys.exit('\nHamiltonian file not recognised -> EXIT\n')

        self.hr = []

        if ( hr_type == 'w90' ):
            for line in lines[4:]:
                if ( abs(float(line.split()[6])) > 1.e-6 ):
                    sys.exit('\nThe hamiltonian must be real, found a complex component -> EXIT\n')
                self.hr.append(line.split()[5])
        
        if ( hr_type == 'kc_occ' ):
            for line in lines[5:-1]:
                self.hr.append(line.split()[0])

        if ( hr_type == 'kc_emp' ):
            for line in lines:
                self.hr.append(line.split()[0])

        if ( len(self.hr) != self.num_wann**2 ):
            sys.exit('\nWrong number of matrix elements for the input hamiltonian -> EXIT(%s)\n' \
                                                                             %(len(self.hr)))

        self.hr = np.array(self.hr, dtype=float).reshape(self.num_wann,self.num_wann)

        # conversion to eV (hamiltonian from CP Koopmans code is in Hartree)
        if ( hr_type == 'kc_occ' or hr_type == 'kc_emp' ):
            self.hr = self.hr * 27.21138386

        # check the hermiticity of the hamiltonian
        for m in range (self.num_wann):
            for n in range(self.num_wann):
                if ( self.hr[m,n] - self.hr[n,m].conjugate() > 1.e-6 ):
                    sys.exit('\nHamiltonian matrix not hermitian -> EXIT(%s,%s)\n' %(m,n))

        return


    """
    parse_phases gets the phases of WFs from the file 'wf_phases.dat'. If the file 
                 is not found a warning is print out and the WFs phases are ignored.
                 The phases are set up as the attribute self.phases
    
    """
    def parse_phases(self):

        self.phases = []

        try:
            with open('wf_phases.dat','r') as ifile:
                lines = ifile.readlines()
        except FileNotFoundError:
            print('\nWARNING: file \'wf_phases.dat\' not found, phases are ignored.\n')
            self.phases = [1] * self.num_wann
            return

        for line in lines:
            self.phases.append(float(line.split()[0]) + \
                               float(line.split()[1])*1j)
        
        return


    """
    scell_centers prints out the list of centers as parsed from W90 output.
                  By passing (optional) argument one can choose between the
                  following units:
                  - 'ang' for Angstrom units (this is the default)
                  - 'bohr' for Bohr units
                  - 'alat' for units of lattice parameter 
                  - 'crys' for crystal units

    """
    def scell_centers(self, units='crys', cell='sc'):

        if ( units is not 'ang'  and units is not 'bohr' and \
             units is not 'alat' and units is not 'crys' ):
            units = 'crys'
            print('\nThe first argument of scell_centers can be \
\'ang\' or \'bohr\' or \'alat\' or \'crys\'. Using default \'crys\'\n')

        if ( cell is not 'sc' and cell is not 'pc' ):
            cell = 'sc'
            print('\nThe second argument of scell_centers can be \
\'sc\' for supercell or \'pc\' for primitive cell. Using default \'sc\'\n')

        if ( hasattr(self, 'alat_sc') ):
            if ( cell is 'sc' ):
                alat = self.alat_sc
                at = self.at_sc
            else:
                alat = self.alat
                at = self.at 
        else:
            alat = self.alat
            at = self.at
            if ( cell is 'pc' ):
                cell = 'sc'
                print('\nPC not available yet. Need to call map_wannier first. \
Using default \'sc\'\n')

        if ( units == 'crys' ):
            return self.centers
        else:
            c_tmp = []
            for n in range(self.num_wann):
                c_tmp.append(crys_to_cart(self.centers[n], at, +1))
                if ( units == 'alat' ):
                    continue
                if ( units == 'ang' ):
                    c_tmp[n] = c_tmp[n] * alat
                if ( units == 'bohr' ):
                    c_tmp[n] = c_tmp[n] * ang_to_bohr(alat, +1)

            return c_tmp


    """
    print_centers simply prints out the centers in the following Xcrysden-readable format:

                  X  0.000  0.000  0.000
                  X  0.000  0.000  0.000
                  *    *      *      *
                  *    *      *      *
                  *    *      *      *
                  X  0.000  0.000  0.000

    """
    def print_centers(self, centers=None):
        if ( centers == None ):
            centers = self.centers
        for n in range(self.num_wann):
            print(' X  %10.6f  %10.6f  %10.6f' %(centers[n][0],centers[n][1],centers[n][2]))
        return
        
  
