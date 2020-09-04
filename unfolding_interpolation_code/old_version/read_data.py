import numpy as np
import re
import sys


# Find num_wann, primitive lattice vectors and WFs centers from Koopmans output
def read_koopmans_output():
	ifile = open('313.out','r')
	lines = ifile.readlines()
	ifile.close()
	
	avec = np.zeros((3,3))
	centers = []
	for line in lines:
		if re.search('Number of Electron',line):
			num_wann = int(line.split()[4])/2
			if num_wann%2!=0:
				sys.exit('\nThe number of wannier function must be even!\n')
		if re.search(' a1 ',line):
			avec[0] = np.array(line.split()[2:],dtype=float)
		if re.search(' a2 ',line):
			avec[1] = np.array(line.split()[2:],dtype=float)
		if re.search(' a3 ',line):
			avec[2] = np.array(line.split()[2:],dtype=float)

	counter = 0
	for line in lines:
		if re.search('OCC',line) and counter<num_wann:
			centers.append(np.array(line.split()[5:8],dtype=float))
			counter = counter + 1

	centers = np.array(centers)
	return num_wann,avec,centers


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
        avec = np.zeros((3,3))
        centers = np.zeros((num_wann,3))
	spreads = np.zeros(num_wann)
	signatures = np.zeros((num_wann,20))
	index = 0

        for n in range(len(lines)):
                if check_latt>0 and n<line_latt+3:
                        avec[n-line_latt] = np.array(lines[n].split()[1:],dtype=float)
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

        return num_wann,avec,centers,spreads,signatures


# Read Hamiltonian H(R) either from Koopmans output or from Wannier90 output
def read_hr_scell(file_hr,hr_type,num_wann,emp_states):

	if emp_states:
		ifile = open(file_hr,'r')
		lines = ifile.readlines()
		ifile.close()
		hr = []
		for line in lines:
			hr.append(float(line))
		hr = np.array(hr).reshape(num_wann,num_wann)
		hr = hr * 27.2113845
		return hr

	else:
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
			hr = hr * 27.2113845    	# in Koopmans code energy is in Hartree units --> conversion to eV
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


# Read Hamiltonian H(R) from the calculation with k-points (used for smooth interpolation)
def read_hr_pcell(file_hr_dft):

	ifile = open(file_hr_dft,'r')
	lines = ifile.readlines()
	ifile.close()

	num_wann_pc = int(lines[1].split()[0])
	nrpts = int(lines[2].split()[0])
	weights = []
	hr_dft = []
	R_pcell = []
	
	nlines_to_skip = 3 + nrpts/15
	if nrpts%15 != 0:	nlines_to_skip += 1

	for n in range(3,nlines_to_skip):
		for m in range(len(lines[n].split())):
			weights.append(lines[n].split()[m])
	for n in range(nlines_to_skip,len(lines),num_wann_pc**2):
		R_pcell.append(lines[n].split()[:3])	
	for n in range(nlines_to_skip,len(lines)):
		hr_dft.append(float(lines[n].split()[5]) + 1.j*float(lines[n].split()[6]))
	weights = np.array(weights,dtype=int)
	R_pcell = np.array(R_pcell,dtype=int)
	hr_dft = np.array(hr_dft).reshape(nrpts,num_wann_pc,num_wann_pc)

	return weights,R_pcell,hr_dft



# Read file 'wf_phases.dat' with WFs phases from Wannier90
def read_phases():

	ifile = open('wf_phases.dat','r')
	lines = ifile.readlines()
	ifile.close()

	phase = []
	for line in lines:
		phase.append(float(line.split()[0]) + 1.j*float(line.split()[1]))

	phase = np.array(phase)

	return phase



