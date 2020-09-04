import numpy as np
import math
import re
import sys,time,os
from read_data import *
from write_data import *
from modules import *

sys.path.append(os.getcwd())
from input import *


try:
	hr_type=='koopmans' or hr_type=='wannier'
except:
	sys.exit("\nhr_type musts be equal to 'koopmans' or 'wannier'.\n")

if do_smooth_interpolation:
	try:
		hr_type_smooth=='koopmans' or hr_type_smooth=='wannier'
	except:
		sys.exit("\nhr_type_smooth musts be equal to 'koopmans' or 'wannier'.\n")

start = time.time()
reset = time.time()

print '\n'
print '\t############################################'
print '\t### UNFOLDING AND INTERPOLATION OF BANDS ###'
print '\t############################################\n'

if double_R:	print '\t ---> The translational symmetry of the Hamiltonian, H(r+R)=H(r), is NOT exploited.'
else:		print '\t ---> The translational symmetry of the Hamiltonian, H(r+R)=H(r), is exploited.'

if read_wf_phases:	print '\t ---> WF phases are taken from Wannier90 and exploited.'
else:			print '\n\tWARNING: WF phases are not exploited!\n'


k_mesh = MP_mesh(nr1,nr2,nr3)
nrtot = nr1*nr2*nr3

###################################
########## READING FILES ##########
###################################

num_wann,avec_sc,centers,spreads,signatures = read_wannier_scell(seedname)
#num_wann,avec_sc,centers = read_koopmans_output()
if do_smooth_interpolation:
	weights,R_pcell,hr_pcell = read_hr_pcell(file_hr_pcell)
if read_wf_phases:	
	wf_phase = read_phases()
	if wf_phase.shape[0] != num_wann:	sys.exit('\nwf_phases.dat not consistent with num_wann\n')

print '\n!!! WARNING: avec is forced to be positive...to check'
avec_sc = abs(avec_sc)

print '\nReading input in:\t\t\t\t\t%.3f sec' %(time.time()-reset)
reset = time.time()

#####################################
#####################################


avec_pc = (avec_sc.transpose()/np.array((nr1,nr2,nr3))).transpose()
bvec_pc = rec_latt_vec(avec_pc)

if num_wann%nrtot!=0:	sys.exit('\nSomething wrong !\nnum_wann must be a multiple of nrtot\n')

if do_smooth_interpolation:		# Check the consistency of the hamiltonian from k-points calculation and of the permutation --- Permute the matrix elements of the hamiltonian of the primitive cell to get the same WF representation of the supercell
	if hr_pcell.shape[1]!=num_wann/nrtot or hr_pcell.shape[2]!=num_wann/nrtot:
		sys.exit('\nSomething wrong !\nThe hamiltonian from the primitive cell calculation is not consistent with the number of WFs\n')
	if len(permutation)!=num_wann/nrtot:	sys.exit('\nSomething wrong !\nThe order of the permutation is not consistent with the number of WFs in the primitive cell\n')
	hr_pcell = permute_hr(hr_pcell,permutation)


# Shift the WFs inside the SC if they are outside
for n in range(num_wann):
	x = np.linalg.solve(avec_sc,centers[n])
	for i in range(3):
		if x[i] >= 0:					x[i] = x[i] - int(x[i])
		elif (x[i] < 0 and abs(x[i]) < cutoff):		x[i] = abs(x[i])
		else:						x[i] = x[i] - int(x[i]) + 1
	centers[n] = np.dot(x,avec_sc)


# Building the map |m> --> |Rn>.
# The variable wann_pc will contain all the WFs of the SC, but reordered following the PC indexing |Rn>.
# In the following loop we find the WFs in the PC only, that is |0m>. Further we find also the other WFs.
wann_pc = np.zeros((nrtot,num_wann/nrtot,25))
index = 0
signature = np.zeros(21)
for n in range(num_wann):
	x = np.linalg.solve(avec_pc,centers[n])
	if (x >= np.zeros(3)).all() and (x < np.ones(3)).all() and (x < (1-cutoff)*np.ones(3)).all():
		wann_pc[0,index,0] = n + 1
		wann_pc[0,index,1:4] = centers[n,:]
		wann_pc[0,index,4] = spreads[n]
		wann_pc[0,index,5:] = signatures[n]
		index = index + 1


# Here we find all the other WFs
R = Rvec(nr1,nr2,nr3)
for m in range(num_wann):
	for n in range(num_wann/nrtot):
		if wann_match(centers[m],spreads[m],signatures[m],wann_pc[0,n,1:],avec_pc,cutoff):
			dist = centers[m] - wann_pc[0,n,1:4]
			x = np.linalg.solve(avec_pc,dist)
			for i in range(nrtot):
				if (np.round(x) == R[i]).all():
					wann_pc[i,n,0] = m + 1
					wann_pc[i,n,1:4] = centers[m,:]
					wann_pc[i,n,4] = spreads[m]
					wann_pc[i,n,5:] = signatures[m]
					break
			break

 
print 'Building the map |m> --> |Rn> in:\t\t\t%.3f sec' %(time.time()-reset)
reset = time.time()


# Read wannier Hamiltonian and rewrite it in terms of the mapped WFs: <i|H|j> --> <Rm|H|R'n> (the translational symmetry is not exploited and all the matrix elements are calculated explicitly!)
hr_old = read_hr_scell(file_hr,hr_type,num_wann,emp_states)
hr_new = np.zeros((nrtot,nrtot,num_wann/nrtot,num_wann/nrtot),dtype=complex)
if do_smooth_interpolation:
	hr_old_smooth = read_hr_scell(file_hr_smooth,hr_type_smooth,num_wann,emp_states_smooth)
	hr_old = hr_old - hr_old_smooth			# Here we define hr_old as Delta(R) = H_KI(R) - H_PBE(R) and then we interpolate Delta(R)
#	print np.linalg.eigvalsh(hr_old)

#for i in range(num_wann):
#	for j in range(num_wann):
#		if i!=j and abs(hr_old[i,j])>1.E-6:	print hr_old[i,j]
#		if i==j:	print hr_old[i,j]

for r1 in range(nrtot):
	for r2 in range(nrtot):
		for m in range(num_wann/nrtot):
			for n in range(num_wann/nrtot):
				i = int(wann_pc[r1,m,0]) - 1
				j = int(wann_pc[r2,n,0]) - 1
				if read_wf_phases:	hr_new[r1,r2,m,n] = hr_old[i,j] * np.conj(wf_phase[i]) * wf_phase[j]
				else:			hr_new[r1,r2,m,n] = hr_old[i,j]

#auxfile = open('V_ki_matrix.dat','w')
#for r1 in range(nrtot):
#	for r2 in range(nrtot):
#		for m in range(num_wann/nrtot):
#			for n in range(num_wann/nrtot):
#				auxfile.write('(R%d,R%d,%d,%d) = %11.8f %11.8f\n' %(r1,r2,m,n,hr_new[r1,r2,m,n].real,hr_new[r1,r2,m,n].imag))
#				if (abs(hr_new[r1,r2,m,n].real) >= 1.E-6) or (abs(hr_new[r1,r2,m,n].imag) >= 1.E-6):
#					print '(R%d,R%d,%d,%d) = %11.8f %11.8f' %(r1,r2,m,n,hr_new[r1,r2,m,n].real,hr_new[r1,r2,m,n].imag)
#
#auxfile.close()

print 'Mapping the Hamiltonian <i|H|j> --> <Rm|H|R\'n> in:\t%.3f sec' %(time.time()-reset)
reset = time.time()


if interpolation:	
	if do_smooth_interpolation:	print '\n\t---> Performing interpolation of bands using smooth interpolation method.\n\n\t!!! WARNING: check the consistency between the WFs in the primitive cell and in the supercell.\n'
	else:				print '\n\t---> Performing interpolation of bands.\n'
else:			print '\n\t---> Calculation of the eigenvalues on the original k-mesh.\n'


# Calculate H(k) and its eigenvalues
eig_k = []
if interpolation:
	hk = np.zeros((k_path.shape[0],num_wann/nrtot,num_wann/nrtot),dtype=complex)
	for kn in range(k_path.shape[0]):
		k = k_path[kn]
		for m in range(num_wann/nrtot):
			for n in range(num_wann/nrtot):
				for i in range(nrtot):
					if not double_R:
						phase_factor,Tnn = ws_distance(nr1,nr2,nr3,avec_sc,k,cutoff,wann_pc[i,n,1:4],wann_pc[0,m,1:4],R[i],m,n)
						hk[kn,m,n] = hk[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k,R[i])) * phase_factor * hr_new[0,i,m,n]	# Hmn(k) = sum_R e^(-ikR) * ( sum_Tnn e^(-ikTnn) ) * <0m|H|Rn>
					else:
						for jj in range(nrtot):
                                                        phase_factor,Tnn = ws_distance(nr1,nr2,nr3,avec_sc,k,cutoff,wann_pc[i,n,1:4],wann_pc[jj,m,1:4])
	                                                hk[kn,m,n] = hk[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k,R[i]-R[jj])) * phase_factor * hr_new[jj,i,m,n] / nrtot     # Hmn(k) = 1/N * sum_{R,R'} e^-ik(R-R')] * <R'm|H|Rn>

				if do_smooth_interpolation:
					hk_pcell = np.zeros((k_path.shape[0],num_wann/nrtot,num_wann/nrtot),dtype=complex)
					for rn in range(R_pcell.shape[0]):
						hk_pcell[kn,m,n] = hk_pcell[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k,R_pcell[rn])) * hr_pcell[rn,m,n] / weights[rn]
					hk[kn,m,n] = hk[kn,m,n] + hk_pcell[kn,m,n]
					
		eig_k.append(np.linalg.eigvalsh(hk[kn,:,:]))

else:	# else the Hamiltonian is calculated over the k-mesh
	hk = np.zeros((k_mesh.shape[0],num_wann/nrtot,num_wann/nrtot),dtype=complex)
	for kn in range(k_mesh.shape[0]):
		k = k_mesh[kn]
		for m in range(num_wann/nrtot):
			for n in range(num_wann/nrtot):
				for i in range(nrtot):
					if not double_R:
						hk[kn,m,n] = hk[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k,R[i])) * hr_new[0,i,m,n]			# Hmn(k) = sum_R e^(-ikR) * <0m|H|Rn>
					else:
						hk[kn,m,n] = hk[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k,abs(R[i]-R[jj]))) * hr_new[jj,i,m,n] / nrtot               # Hmn(k) = 1/N * sum_{R,R'} e^-ik(R-R')] * <R'm|H|Rn>
#				if do_smooth_interpolation:
#                                        hk_pcell = np.zeros((k_mesh.shape[0],num_wann/nrtot,num_wann/nrtot),dtype=complex)
#                                        for rn in range(R_pcell.shape[0]):
#                                                hk_pcell[kn,m,n] = hk_pcell[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k,R_pcell[rn])) * hr_pcell[rn,m,n] / weights[rn]
#                                        hk[kn,m,n] = hk[kn,m,n] + hk_pcell[kn,m,n]

		eig_k.append(np.linalg.eigvalsh(hk[kn,:,:]))

eig_k = np.array(eig_k)


print 'Calculating the energies in:\t\t\t\t%.3f sec' %(time.time()-reset)
reset = time.time()


if DOS:
	if not interpolation:	print '\nWARNING: the DOS will be calculated on the original k-mesh !\n'
	
	deltaE = float(Emax-Emin)/nstep
	ene = Emin
	gaussian = []
	for n in range(nstep+1):
		gaussian.append(2 * np.sum(np.exp(-(ene-eig_k)**2/degauss)) * degauss**.5 / np.sqrt(np.pi) / len(eig_k))
		ene = ene + deltaE

	print 'Calculating the DOS in:\t\t\t\t\t%.3f sec' %(time.time()-reset)
	reset = time.time()



# Write the output
if interpolation:	
	k_path = path_to_plot(k_path,bvec_pc,eq_points,point)
	write_output_eigk(eig_k,'band',k_path)
else:
	write_output_eigk(eig_k,'kvec')	

if DOS:			write_dos(gaussian,Emin,deltaE)

print 'Writing the output in:\t\t\t\t\t%.3f sec' %(time.time()-reset)
print '\nTotal time:\t\t\t\t\t\t%.3f sec\n' %(time.time()-start)
print '\tALL DONE.\n'


if (os.path.isfile('input.pyc')):
	os.remove('input.pyc')

