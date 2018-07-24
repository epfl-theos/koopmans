import numpy as np
import re
import sys


# Find the k-mesh from wannier input file of the PC
def read_wannier_pcell(dirname,seedname):
	ifile = open(dirname+seedname+'.win','r')
	lines = ifile.readlines()
	ifile.close()
	check_kpoints = False
	k_mesh = []

	for n in range(len(lines)):
		if re.search('mp_grid',lines[n]):
			nk1,nk2,nk3 = np.array(lines[n].split()[2:],dtype=int)
			nktot = nk1*nk2*nk3
		if check_kpoints and n<line+nktot:
			k_mesh.append(np.array(lines[n].split(),dtype=float))
 		if re.search('begin kpoints',lines[n]):
			line = n + 1
			check_kpoints = True
	
	try:		nk1,nk2,nk3
	except:		sys.exit('\nCannot find the k-mesh !\n')

	if not check_kpoints:	sys.exit('\nCannot find kpoints !\n')

	k_mesh = np.array(k_mesh)

	return nk1,nk2,nk3,k_mesh


# Find num_wann, primitive lattice vectors and WFs from wannier input/output files of the SC
def read_wannier_scell(dirname,seedname):
	ifile = open(dirname+seedname+'.win','r')
	lines = ifile.readlines()
	ifile.close()
	
	for n in range(len(lines)):
		if re.search('num_wann',lines[n]):
			num_wann = int(lines[n].split()[2])
			break
	
	try:		num_wann
	except:		sys.exit('\nCannot find num_wann !\n')
	
	ifile = open(dirname+seedname+'.wout','r')
	lines = ifile.readlines()
	ifile.close()
	check_latt = 0
	check_wann = 0
	latt_vec = np.zeros((3,3))
	centers = np.zeros((num_wann,3))

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
		if re.search('Final State',lines[n]):
			line_wann = n + 1
			check_wann = check_wann + 1

	if check_latt > 1:
		print '\nWARNING: multiple lattice vectors found: here the last set of lattice vectors is taken.\n'
	elif check_latt == 0:
		sys.exit('\nNo lattice vectors found !\n')

	if check_wann > 1:
		print '\nWARNING: multiple final states found: here the last final state is evaluated.\n'
	elif check_wann == 0:
		sys.exit('\nNo final state found !\n')

	return num_wann,latt_vec,centers


# Read wannier Hamiltonian 
def read_hr(dirname,seedname,num_wann):
	ifile = open(dirname+seedname+'_hr.dat','r')
	lines = ifile.readlines()
	ifile.close()
	hr = np.zeros((num_wann,num_wann),dtype=complex)

	for i in range(num_wann):
		for j in range(num_wann):
#			print lines[i*num_wann+j+4].split()[5]
		        hr[i,j] = float(lines[i*num_wann+j+4].split()[5]) + 1j*float(lines[i*num_wann+j+4].split()[6])

	return hr
