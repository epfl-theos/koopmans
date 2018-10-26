import numpy as np
import math
import re
import sys,time
from inputfile import *
from read_data import *
from write_data import *
from modules import *



try:
	hr_type=='koopmans' or hr_type=='wannier'
except:
	sys.exit("\nhr_type musts be equal to 'koopmans' or 'wannier'.\n")


start = time.time()
reset = time.time()

k_mesh = MP_mesh(nk1,nk2,nk3)
nktot = nk1*nk2*nk3

num_wann,latt_vec_sc,centers,spreads,signatures = read_wannier_scell(dir_wann,seedname)
#num_wann,latt_vec_sc,centers = read_koopmans_output(dir_wann)

print '\nReading input in:\t\t\t\t\t%.3f sec' %(time.time()-reset)
reset = time.time()

latt_vec_pc = (latt_vec_sc.transpose()/np.array((nk1,nk2,nk3))).transpose()

if num_wann%nktot!=0:	sys.exit('\nSomething wrong !\nnum_wann must be a multiple of nktot\n')

#for n in range(num_wann):
#	print ' X  %.6f  %.6f  %.6f' %(centers[n,0]/0.52917721067,centers[n,1]/0.52917721067,centers[n,2]/0.52917721067)
#	print ' %d  %.6f  %.6f  %.6f' %(n,centers[n,0],centers[n,1],centers[n,2])


# Shift the WFs inside the SC if they are outside
for n in range(num_wann):
	x = np.linalg.solve(latt_vec_sc,centers[n])
	for i in range(3):
		if x[i] >= 0:					x[i] = x[i] - int(x[i])
		elif (x[i] < 0 and abs(x[i]) < cutoff):		x[i] = abs(x[i])
		else:						x[i] = x[i] - int(x[i]) + 1
	centers[n] = np.dot(x,latt_vec_sc)


# Building the map |m> --> |Rn>.
# The variable wann_pc will contain all the WFs of the SC, but reordered following the PC indexing |Rn>.
# In the following loop we find the WFs in the PC only, that is |0m>. Further we find also the other WFs.
wann_pc = np.zeros((nktot,num_wann/nktot,25))
index = 0
signature = np.zeros(21)
for n in range(num_wann):
	x = np.linalg.solve(latt_vec_pc,centers[n])
	if (x >= np.zeros(3)).all() and (x < np.ones(3)).all() and (x < (1-cutoff)*np.ones(3)).all():
		wann_pc[0,index,0] = n + 1
		wann_pc[0,index,1:4] = centers[n,:]
		wann_pc[0,index,4] = spreads[n]
		wann_pc[0,index,5:] = signatures[n]
		index = index + 1


# Here we find all the other WFs
R = Rvec(nk1,nk2,nk3)
for m in range(num_wann):
	for n in range(num_wann/nktot):
		if wann_match(centers[m],spreads[m],signatures[m],wann_pc[0,n,1:],latt_vec_pc,cutoff):
			dist = centers[m] - wann_pc[0,n,1:4]
			x = np.linalg.solve(latt_vec_pc,dist)
			for i in range(nktot):
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
hr_old = read_hr(file_hr,hr_type,num_wann)

if hr_type=='koopmans':		hr_new = np.zeros((nktot,nktot,num_wann/nktot,num_wann/nktot))
elif hr_type=='wannier':	hr_new = np.zeros((nktot,nktot,num_wann/nktot,num_wann/nktot),dtype=complex)
 
for r1 in range(nktot):
	for r2 in range(nktot):
		for m in range(num_wann/nktot):
			for n in range(num_wann/nktot):
				i = int(wann_pc[r1,m,0]) - 1
				j = int(wann_pc[r2,n,0]) - 1
				hr_new[r1,r2,m,n] = hr_old[i,j]

print 'Mapping the Hamiltonian <i|H|j> --> <Rm|H|R\'n> in:\t%.3f sec' %(time.time()-reset)
reset = time.time()


# Calculate H(k) and its eigenvalues
eig_k = []
if interpolation:
	hk = np.zeros((k_path.shape[0],num_wann/nktot,num_wann/nktot),dtype=complex)
	for kn in range(k_path.shape[0]):
		k = k_path[kn]
		for m in range(num_wann/nktot):
			for n in range(num_wann/nktot):
				for i in range(nktot):
					phase_factor,Tnn = ws_distance(nk1,nk2,nk3,latt_vec_sc,k,cutoff,wann_pc[i,n,1:4],wann_pc[0,m,1:4])
					hk[kn,m,n] = hk[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k,R[i])) * phase_factor * hr_new[0,i,m,n]	# Hmn(k) = sum_R e^(-ikR) * ( sum_Tnn e^(-ikTnn) ) * <0m|H|Rn>
		eig_k.append(np.linalg.eigvalsh(hk[kn,:,:]))

else:	# else the Hamiltonian is calculated over the k-mesh
	hk = np.zeros((k_mesh.shape[0],num_wann/nktot,num_wann/nktot),dtype=complex)
	for kn in range(k_mesh.shape[0]):
		k = k_mesh[kn]
		for m in range(num_wann/nktot):
			for n in range(num_wann/nktot):
				for i in range(nktot):
					hk[kn,m,n] = hk[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k,R[i])) * hr_new[0,i,m,n]			# Hmn(k) = sum_R e^(-ikR) * <0m|H|Rn>
		eig_k.append(np.linalg.eigvalsh(hk[kn,:,:]))

eig_k = np.array(eig_k)

if hr_type=='koopmans':		eig_k = eig_k * 27.2113845	# in Koopmans code energy is in Hartree units --> conversion to eV

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
if interpolation:	write_output_eigk(eig_k,'band',outdir,k_path)
else:			write_output_eigk(eig_k,'kvec',outdir)
if DOS:			write_dos(outdir,gaussian,Emin,deltaE)

print 'Writing the output in:\t\t\t\t\t%.3f sec' %(time.time()-reset)
print '\nTotal time:\t\t\t\t\t\t%.3f sec\n' %(time.time()-start)

