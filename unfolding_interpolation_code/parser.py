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
    alat_sc         : SC lattice parameter (in BOHR units). Like celldm(1) in QE it refers to a_1
    nr1             : SC dimension along PC primitive vector a_1
    nr2             : SC dimension along PC primitive vector a_2
    nr3             : SC dimension along PC primitive vector a_3
    w90_calc        : type of Wannier90 calculation: 'pc' for PC with k-points, 'sc' for SC Gamma-only
    do_map          : if True the algorithm to map the WFs from the SC to the PC is activated
    use_ws_distance : if True the Wigner-Seitz distance between WF centers is used
    k_path          : k_path for bands interpolation (in crystal units)
    smooth_int      : if True, smooth interpolation method is used
    file_hr_coarse  : look at documentation (needed when smooth_int=True)
    file_hr_smooth  : look at documentation (needed when smooth_int=True)
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
        self.list_of_args = [ 'seedname', 'alat_sc', 'nr1', 'nr2', 'nr3', 'w90_calc', 'do_map',\
                              'use_ws_distance', 'k_path', 'smooth_int', 'file_hr_coarse', \
                              'file_hr_smooth', 'do_dos', 'degauss', 'nstep', 'Emin', 'Emax' ]
        
        if 'seedname' not in json_data.keys():    sys.exit('\nMissing \'seedname\' in input -> EXIT\n')
        if 'alat_sc' not in json_data.keys():     sys.exit('\nMissing \'alat_sc\' in input -> EXIT\n')
        if 'nr1' not in json_data.keys():         sys.exit('\nMissing \'nr1\' in input -> EXIT\n')
        if 'nr2' not in json_data.keys():         sys.exit('\nMissing \'nr2\' in input -> EXIT\n')
        if 'nr3' not in json_data.keys():         sys.exit('\nMissing \'nr3\' in input -> EXIT\n')
        
        ### Reading json dictionary and defining defaults ###
        self.seedname = json_data.get("seedname", None)
        self.alat_sc = json_data.get("alat_sc", None)
        self.nr1 = json_data.get("nr1", None)
        self.nr2 = json_data.get("nr2", None)
        self.nr3 = json_data.get("nr3", None)
        self.w90_calc = json_data.get("w90_calc", "pc")
        self.do_map = json_data.get("do_map", False)
        self.use_ws_distance = json_data.get("use_ws_distance", True)
        self.k_path = json_data.get("k_path", None)
        self.smooth_int = json_data.get("smooth_int", False)
        self.file_hr_coarse = json_data.get("file_hr_coarse", None)
        self.file_hr_smooth = json_data.get("file_hr_smooth", None)
        self.do_dos = json_data.get("do_dos", False)
        self.degauss = json_data.get("degauss", 0.05)
        self.nstep = json_data.get("nstep", 1000)
        self.Emin = json_data.get("Emin")
        self.Emax = json_data.get("Emax")
        ####################################################

        # checks on the input arguments
        if ( type(self.seedname) is not str ):
            sys.exit('\n\'seedname\' must be a string -> EXIT(%s)' %type(self.seedname))
        if ( type(self.alat_sc) is not float and type(self.alat_sc) is not int ):
            sys.exit('\n\'alat_sc\' must be a float or a int -> EXIT(%s)' %type(self.alat_sc))
        if ( type(self.nr1) is not int ):
            sys.exit('\n\'nr1\' must be a int -> EXIT(%s)' %type(self.nr1))
        if ( type(self.nr2) is not int ):
            sys.exit('\n\'nr2\' must be a int -> EXIT(%s)' %type(self.nr2))
        if ( type(self.nr3) is not int ):
            sys.exit('\n\'nr3\' must be a int -> EXIT(%s)' %type(self.nr3))
        if ( type(self.w90_calc) is not str ):
            sys.exit('\n\'w90_calc\' must be a string -> EXIT(%s)' %type(self.w90_calc))
        if ( type(self.do_map) is not bool ):
            sys.exit('\n\'do_map\' must be a bool -> EXIT(%s)' %type(self.do_map))
        if ( type(self.use_ws_distance) is not bool ):
            sys.exit('\n\'use_ws_distance\' must be a bool -> EXIT(%s)' %type(self.use_ws_distance))
        if ( type(self.k_path) is not list ):
            sys.exit('\nk_path must be a list -> EXIT(%s)' %type(self.k_path))
        if ( type(self.smooth_int) is not bool ):
            sys.exit('\n\'smooth_int\' must be a bool -> EXIT(%s)' %type(self.smooth_int))
        if ( type(self.file_hr_coarse) is not str ):
            sys.exit('\n\'file_hr_coarse\' must be a string -> EXIT(%s)' %type(self.file_hr_coarse))
        if ( type(self.file_hr_smooth) is not str ):
            sys.exit('\n\'file_hr_smooth\' must be a string -> EXIT(%s)' %type(self.file_hr_smooth))
        if ( type(self.do_dos) is not bool ):
            sys.exit('\n\'do_dos\' must be a bool -> EXIT(%s)' %type(self.do_dos))
        if ( type(self.degauss) is not float and type(self.degauss) is not int ):
            sys.exit('\n\'degauss\' must be a float or a int -> EXIT(%s)' %type(degauss))
        if ( type(self.nstep) is not float and type(self.nstep) is not int ):
            sys.exit('\n\'nstep\' must be a float or a int -> EXIT(%s)' %type(self.nstep))
        if ( self.Emin is not None and type(self.Emin) is not float and type(self.Emin) is not int ):
            sys.exit('\n\'Emin\' must be a float or a int -> EXIT(%s)' %type(self.Emin))
        if ( self.Emax is not None and type(self.Emax) is not float and type(self.Emax) is not int ):
            sys.exit('\n\'Emax\' must be a float or a int -> EXIT(%s)' %type(self.Emax))
           
        for item in json_data.keys():
            if ( item not in self.list_of_args ):
                print('WARNING: argument \'%s\' is unknown!' %item)
       
        self.alat_sc = ang_to_bohr(self.alat_sc, -1)    # conversion to Ang

        self.w90_calc = self.w90_calc.lower()
        if ( self.w90_calc != 'pc' and self.w90_calc != 'sc' ):
            sys.exit('\n\'w90_calc\' must be \'pc\' or \'sc\' -> EXIT(%s)' %(w90_calc))

        if ( self.w90_calc == 'sc' ):
            self.w90_input_sc = True
        else:
            self.w90_input_sc = False
            if ( self.do_map ):
                sys.exit('\n\do_map=True incompatible with w90_calc=\'pc\' -> EXIT\n')
  
        if ( self.k_path == None ):
            self.kvec = MP_mesh(self.nr1,self.nr2,self.nr3)
            print('WARNING: \'k_path\' missing in input, the energies are calculated on a ' +
                                                    'commensurate Monkhorst-Pack mesh')
        else:
            self.kvec = generate_path(self.k_path)
        
        if ( self.smooth_int ):
            if ( 'file_hr_coarse' not in json_data.keys() ):
                sys.exit('\nMissing \'file_hr_coarse\' for smooth interpolation -> EXIT\n')
            if ( 'file_hr_smooth' not in json_data.keys() ):
                sys.exit('\nMissing \'file_hr_smooth\' for smooth interpolation -> EXIT\n')
        
        if ( self.do_dos ):
            if 'degauss' not in json_data.keys():
                print('WARNING: \'degauss\' missing in input, using default value')
            if 'nstep' not in json_data.keys():
                print('WARNING: \'nstep\' missing in input, using default value')
            if 'Emin' not in json_data.keys():
                print('WARNING: \'Emin\' missing in input, using default value')
            if 'Emax' not in json_data.keys():
                print('WARNING: \'Emax\' missing in input, using default value')

        return


    """
    parse_w90 gets from the W90 output the lattice vectors, and the centers and spreads
              of the Wannier functions. 

    at      : basis vectors of direct lattice (in PC alat units)
    bg      : basis vectors of reciprocal lattice (in PC 2pi/alat units)
    centers : centers of WFs (in PC crystal units)
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
                num_wann = int(line.split()[6])
            if ( ' a_1 ' in line ):
                self.at.append(np.array(line.split()[1:], dtype=float))
            if ( ' a_2 ' in line ):
                self.at.append(np.array(line.split()[1:], dtype=float))
            if ( ' a_3 ' in line ):
                self.at.append(np.array(line.split()[1:], dtype=float))
            if ( count > 0 and count <= num_wann ):
                start = line.find('(')
                end = line.find(')')
                self.centers.append(np.array(line[start+1:end].replace(',',' ').split(), \
                                                                           dtype=float))
                self.spreads.append(float(line.split()[-1]))

                count += 1
            if ( 'Final State' in line ):
                count += 1

        # primitive cell lattice parameter
        self.alat = self.alat_sc / self.nr1

        self.at = np.array(self.at, dtype=float).reshape(3,3) / self.alat
        self.bg = np.linalg.inv(self.at).transpose()
        self.at_sc = np.zeros((3,3))
        self.bg_sc = np.zeros((3,3))
        self.Rvec = latt_vect(self.nr1, self.nr2, self.nr3)
        
        if ( self.w90_input_sc ):
            self.num_wann_sc = num_wann
            self.num_wann = int(num_wann / (self.nr1*self.nr2*self.nr3))
        
            # converting at and bg from the SC to the PC
            self.at_sc[0] = self.at[0]
            self.at_sc[1] = self.at[1]
            self.at_sc[2] = self.at[2]
            self.at[0] = self.at[0] / self.nr1
            self.at[1] = self.at[1] / self.nr2
            self.at[2] = self.at[2] / self.nr3
            #
            self.bg_sc[0] = self.bg[0]
            self.bg_sc[1] = self.bg[1]
            self.bg_sc[2] = self.bg[2]
            self.bg[0] = self.bg[0] * self.nr1
            self.bg[1] = self.bg[1] * self.nr2
            self.bg[2] = self.bg[2] * self.nr3

        else:
            self.num_wann = num_wann
            self.num_wann_sc = num_wann * (self.nr1*self.nr2*self.nr3)

            self.at_sc[0] = self.at[0] * self.nr1
            self.at_sc[1] = self.at[1] * self.nr2
            self.at_sc[2] = self.at[2] * self.nr3
            #
            self.bg_sc[0] = self.bg[0] / self.nr1
            self.bg_sc[1] = self.bg[1] / self.nr2
            self.bg_sc[2] = self.bg[2] / self.nr3


        for n in range(num_wann):
            self.centers[n] = self.centers[n] / self.alat
            self.centers[n] = crys_to_cart(self.centers[n], self.bg, -1)
        
        # generate the centers and spreads of all the other (R/=0) WFs
        if ( not self.w90_input_sc ):
            for rvect in self.Rvec[1:]:
                for n in range(self.num_wann):
                    self.centers.append( self.centers[n] + rvect )
                    self.spreads.append( self.spreads[n] )

        return
    
    
    """
    parse_hr reads the hamiltonian file passed as sys.argv[1] and it sets it up
             as attribute self.hr
             
    there are 3 possible types of file:
      - w90 file normally called seedname_hr.dat
      - kc_occ file normally called hamiltonian1.xml
      - kc_emp file normally called hamiltonian_emp.dat

    nb: kc_emp must be called 'hamiltonian_emp.dat' otherwise the code may crash
        or misread the matrix elements. if the file name is different the code 
        should be updated.

    """
    def parse_hr(self, file_hr):

        with open(file_hr, 'r') as ifile:
            lines = ifile.readlines()

        if ( 'written on' in lines[0] ):                        hr_type = 'w90'
        elif ( 'xml version' in lines[0]):                      hr_type = 'kc_occ'
        elif ( file_hr[-19:] == 'hamiltonian_emp.dat' ):        hr_type = 'kc_emp'
        else:		sys.exit('\nHamiltonian file not recognised -> EXIT\n')

        self.hr = []

        if ( hr_type == 'w90' ):

            if ( self.w90_input_sc ):
                if ( int(lines[1].split()[0]) != self.num_wann_sc ):
                    sys.exit('\nIn parse_hr inconsistency in num_wann\n')
            else:
                if ( int(lines[1].split()[0]) != self.num_wann ):
                    sys.exit('\nIn parse_hr inconsistency in num_wann\n')
            
            if ( not self.w90_input_sc ): rvect = [] 
            nrpts = int(lines[2].split()[0])
            lines_to_skip = 3 + int(nrpts/15)
            if ( nrpts%15 > 0 ): lines_to_skip += 1
            counter = 0

            for line in lines[lines_to_skip:]:
                if ( abs(float(line.split()[6])) > 1.e-6 ):
                    sys.exit('\nThe hamiltonian must be real, found a complex component -> EXIT\n')

                self.hr.append(line.split()[5])

                counter += 1
                if ( not self.w90_input_sc and counter == self.num_wann**2 ):
                    rvect.append( np.array(line.split()[0:3], dtype=int) )
                    counter = 0

        if ( hr_type == 'kc_occ' ):
            for line in lines[5:-2]:
                self.hr.append(line.split()[0])

        if ( hr_type == 'kc_emp' ):
            for line in lines:
                self.hr.append(line.split()[0])

        if ( hr_type == 'w90' and not self.w90_input_sc ):
            if ( len(self.hr) != nrpts*self.num_wann**2 ):
                sys.exit('\nWrong number of matrix elements for the input hamiltonian -> EXIT(%s)\n' \
                                                                             %(len(self.hr)))
            self.hr = np.array(self.hr, dtype=float).reshape(nrpts,self.num_wann,self.num_wann)
            self.hr = extract_hr(self.hr, rvect, self.nr1, self.nr2, self.nr3)
            self.hr = self.hr.reshape(self.num_wann_sc,self.num_wann)
        else:
            if ( len(self.hr) != self.num_wann_sc**2 ):
                sys.exit('\nWrong number of matrix elements for the input hamiltonian -> EXIT(%s)\n' \
                                                                             %(len(self.hr)))
            self.hr = np.array(self.hr, dtype=float).reshape(self.num_wann_sc,self.num_wann_sc) 

        # conversion to eV (hamiltonian from CP Koopmans code is in Hartree)
        if ( hr_type == 'kc_occ' or hr_type == 'kc_emp' ):
            self.hr = self.hr * 27.21138386

        # check the hermiticity of the hamiltonian (except for H_R(m,n))
        if ( not (hr_type == 'w90' and not self.w90_input_sc) ):
            for m in range (self.num_wann):
                for n in range(self.num_wann):
                    if ( self.hr[m,n] - self.hr[n,m].conjugate() > 1.e-6 ):
                        sys.exit('\nHamiltonian matrix not hermitian -> EXIT(%s,%s)\n' %(m,n))


        # reading the 2 hamiltonians for the smooth interpolation method
        if ( self.smooth_int ):
            self.hr_coarse = []
            self.hr_smooth = []            

            # parsing hr_coarse
            with open(self.file_hr_coarse, 'r') as ifile:
                lines = ifile.readlines()

            if ( self.w90_input_sc ):
                if ( int(lines[1].split()[0]) != self.num_wann_sc ):
                    sys.exit('\nIn parse_hr inconsistency in num_wann in hr_coarse\n')
            else:
                if ( int(lines[1].split()[0]) != self.num_wann ):
                    sys.exit('\nIn parse_hr inconsistency in num_wann in hr_coarse\n')
            
            if ( not self.w90_input_sc ): rvect = []
            nrpts = int(lines[2].split()[0])
            lines_to_skip = 3 + int(nrpts/15)
            if ( nrpts%15 > 0 ): lines_to_skip += 1
            counter = 0

            for line in lines[lines_to_skip:]:
                if ( abs(float(line.split()[6])) > 1.e-6 ):
                    sys.exit('\nhr_coarse must be real, found a complex component -> EXIT\n')

                self.hr_coarse.append(line.split()[5])

                counter += 1
                if ( not self.w90_input_sc and counter == self.num_wann**2 ):
                    rvect.append( np.array(line.split()[0:3], dtype=int) )
                    counter = 0

            if ( self.w90_input_sc ):
                if ( len(self.hr_coarse) != self.num_wann_sc**2 ):
                    sys.exit('\nWrong number of matrix elements for hr_coarse -> EXIT(%s)\n' \
                                                                             %(len(self.hr_coarse)))
                self.hr_coarse = np.array(self.hr_coarse, dtype=float)
                self.hr_coarse = self.hr_coarse.reshape(self.num_wann_sc,self.num_wann_sc)
                self.hr_coarse = self.hr_coarse[:,:self.num_wann]
            else:
                if ( len(self.hr_coarse) != nrpts*self.num_wann**2 ):
                    sys.exit('\nWrong number of matrix elements for hr_coarse -> EXIT(%s)\n' \
                                                                             %(len(self.hr_coarse)))
                self.hr_coarse = np.array(self.hr_coarse, dtype=float)
                self.hr_coarse = self.hr_coarse.reshape(nrpts,self.num_wann,self.num_wann)
                self.hr_coarse = extract_hr(self.hr_coarse, rvect, self.nr1, self.nr2, self.nr3)
                self.hr_coarse = self.hr_coarse.reshape(self.num_wann_sc,self.num_wann)

            # parsing hr_smooth
            with open(self.file_hr_smooth, 'r') as ifile:
                lines = ifile.readlines()

            if ( int(lines[1].split()[0]) != self.num_wann ):
                sys.exit('\nIn parse_hr inconsistency in num_wann in hr_smooth\n')
            
            weights = []
            rvect = []
            nrpts = int(lines[2].split()[0])
            lines_to_skip = 3 + int(nrpts/15)
            if ( nrpts%15 > 0 ): lines_to_skip += 1
            counter = 0

            for line in lines[3:lines_to_skip]:
                for n in range(len(line.split())):
                    weights.append(int(line.split()[n]))

            for line in lines[lines_to_skip:]:
                if ( abs(float(line.split()[6])) > 1.e-6 ):
                    sys.exit('\nhr_smooth must be real, found a complex component -> EXIT\n')

                self.hr_smooth.append(line.split()[5])

                counter += 1
                if ( counter == self.num_wann**2 ):
                    rvect.append( np.array(line.split()[0:3], dtype=int) )
                    counter = 0

            if ( len(self.hr_smooth) != nrpts*self.num_wann**2 ):
                sys.exit('\nWrong number of matrix elements for hr_smooth -> EXIT(%s)\n' \
                                                                         %(len(self.hr_smooth)))
            self.wRs = weights
            self.Rsmooth = rvect
            self.hr_smooth = np.array(self.hr_smooth, dtype=float)
            self.hr_smooth = self.hr_smooth.reshape(nrpts,self.num_wann,self.num_wann)

            # check consistency between hr_coarse and hr_smooth
#            hr_smooth = np.copy(self.hr_smooth)
#            hr_smooth = extract_hr(hr_smooth, rvect, self.nr1, self.nr2, self.nr3)
#            hr_smooth = hr_smooth.reshape(self.num_wann_sc,self.num_wann)
#            if ( np.max(abs( self.hr_coarse - hr_smooth )) > 1.e-3 ):
#                print('WARNING: hr_coarse and hr_smooth differ. Be careful with the results')

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
            self.phases = [1] * self.num_wann_sc
            if ( self.w90_input_sc ):
                print('WARNING: file \'wf_phases.dat\' not found, phases are ignored.')
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
                print('\nPC not available yet. Using default \'sc\'\n')

        if ( units == 'crys' ):
            return self.centers
        else:
            c_tmp = []
            for n in range(self.num_wann_sc):
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
        for n in range(self.num_wann_sc):
            print(' X  %10.6f  %10.6f  %10.6f' %(centers[n][0],centers[n][1],centers[n][2]))
        return
