import numpy as np
import re
import sys


# Find num_wann, primitive lattice vectors and WFs centers from Koopmans output
def read_koopmans_output():
	ifile = open('313.out','r')
	lines = ifile.readlines()
	ifile.close()
	
	latt_vec = np.zeros((3,3))
	centers = []
	for line in lines:
		if re.search('Number of Electron',line):
			num_wann = int(line.split()[4])/2
			if num_wann%2!=0:
				sys.exit('\nThe number of wannier function must be even!\n')
		if re.search(' a1 ',line):
			latt_vec[0] = np.array(line.split()[2:],dtype=float)
		if re.search(' a2 ',line):
			latt_vec[1] = np.array(line.split()[2:],dtype=float)
		if re.search(' a3 ',line):
			latt_vec[2] = np.array(line.split()[2:],dtype=float)

	counter = 0
	for line in lines:
		if re.search('OCC',line) and counter<num_wann:
			centers.append(np.array(line.split()[5:8],dtype=float))
			counter = counter + 1

	centers = np.array(centers)
	return num_wann,latt_vec,centers


def read_wannier_scell(seedname):
        ifile = open(seedname+'.win','r')
        lines = ifile.readlines()
        ifile.close()

        for n in range(len(lines)):
                if re.search('num_wann',lines[n]):
                        num_wann = int(lines[n].split()[2])
                        break

        try:            num_wann
        except:         sys.exit('\nCannot find num_wann !\n')

        ifile = open(seedname+'.wout','r')
        lines = ifile.readlines()
        ifile.close()
        check_latt = 0
        check_wann = 0
        latt_vec = np.zeros((3,3))
        centers = np.zeros((num_wann,3))
	spreads = np.zeros(num_wann)
	signatures = np.zeros((num_wann,20))
	index = 0

        for n in range(len(lines)):
                if check_latt>0 and n<line_latt+3:
                        latt_vec[n-line_latt] = np.array(lines[n].split()[1:],dtype=float)
                if re.search('Lattice Vectors',lines[n]):
                        line_latt = n + 1
                        check_latt = check_latt + 1
                if check_wann>0 and n<line_wann+num_wann:
                        centers[n-line_wann][0] = float(lines[n].split()[6].replace(',',''))
                        centers[n-line_wann][1] = float(lines[n].split()[7].replace(',',''))
                        centers[n-line_wann][2] = float(lines[n].split()[8].replace(',',''))
			spreads[n-line_wann] = float(lines[n].split()[10])  
                if re.search('Final State',lines[n]):
                        line_wann = n + 1
                        check_wann = check_wann + 1
		if re.search('Wannier function:',lines[n]):
			for m in range(1,21):
				signatures[index,m-1] = float(lines[n+m].split()[1])
			index = index + 1

        if check_latt > 1:
                print '\nWARNING: multiple lattice vectors found: here the last set of lattice vectors is taken.\n'
        elif check_latt == 0:
                sys.exit('\nNo lattice vectors found !\n')

        if check_wann > 1:
                print '\nWARNING: multiple final states found: here the last final state is evaluated.\n'
        elif check_wann == 0:
                sys.exit('\nNo final state found !\n')

        return num_wann,latt_vec,centers,spreads,signatures


# Read Hamiltonian H(R) either from Koopmans output or from Wannier90 output
def read_hr(file_hr,hr_type,num_wann):

	if hr_type=='koopmans':
		ifile = open(file_hr,'r')
		lines = ifile.readlines()
		ifile.close()
		check = False
		hr = []
		for line in lines[:-2]:
			if check:
				hr.append(float(line))
			if re.search('size',line):
				check = True
		hr = np.array(hr).reshape(num_wann,num_wann)
		return hr

	if hr_type=='wannier':
		ifile = open(file_hr,'r')
	        lines = ifile.readlines()
        	ifile.close()
	        hr = np.zeros((num_wann,num_wann),dtype=complex)
        	for i in range(num_wann):
                	for j in range(num_wann):
                        	hr[i,j] = float(lines[i*num_wann+j+4].split()[5]) + 1j*float(lines[i*num_wann+j+4].split()[6])
	        return hr
