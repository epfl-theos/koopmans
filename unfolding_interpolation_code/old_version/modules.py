import numpy as np
import copy


# Definition of the Monkhorst-Pack k-mesh (in crystal coordinates) of the PC commensurate to the SC
def MP_mesh(nr1,nr2,nr3):
	k_mesh = []
	for i in range(nr1):
		for j in range(nr2):
			for k in range(nr3):
				k_mesh.append(np.array((i,j,k),dtype=float))
	k_mesh = np.array(k_mesh) / np.array((nr1,nr2,nr3))
	return k_mesh


# Calculate the nrtot R-vectors in crystal coordinates
def Rvec(nr1,nr2,nr3):
	Rvec = []
	for i in range(nr1):
		for j in range(nr2):
			for k in range(nr3):
				Rvec.append(np.array((i,j,k)))
	Rvec = np.array(Rvec)
	return Rvec


# Compare two WFs in the SC and return True if they differ just by a primitive lattice vector (of the primitive cell)
def wann_match(centers,spreads,signatures,wann_ref,avec,cutoff):
	dist = centers - wann_ref[:3]
	x = np.linalg.solve(avec,dist)
	if np.linalg.norm(abs(x-np.round(x))) < cutoff:
		if abs(wann_ref[3]-spreads) < cutoff:
			if (abs(wann_ref[4:]-signatures)<cutoff*np.ones(20)).all():
				return True
	else:
		return False


def ws_distance(nr1,nr2,nr3,avec,kvec,cutoff,wann,wann_ref,rvect,iband,jband):
	dist_min = np.linalg.norm(wann - wann_ref)
	Tnn = []		# List of degenerate T-vectors, corresponding to equidistant nn WFs
	for i in range(-1,2):
		for j in range(-1,2):
			for k in range(-1,2):
				T = i*avec[0] + j*avec[1] + k*avec[2]	# T is a lattice vector of the SC
				dist = np.linalg.norm(wann - wann_ref + T)		# T-translation of the WF (evaluation of the copies of the given WF)
				if abs(dist - dist_min) < cutoff:			# The position (in PC crystal coordinates) of the copy of the given WF is accepted	
					Tnn.append(np.array((i*nr1,j*nr2,k*nr3),dtype=int))	
				elif dist < dist_min:					# A new nearest neighbor is found	
					dist_min = dist
					Tnn = [np.array((i*nr1,j*nr2,k*nr3),dtype=int)]
				else:
					continue

	Tnn = np.array(Tnn)
	phase_factor = np.sum(np.exp(1j*2*np.pi*np.dot(Tnn,kvec))) / len(Tnn)

	return phase_factor,Tnn


def rec_latt_vec(avec):
	bvec = np.zeros((3,3))
	bvec[0] = np.linalg.solve(avec,np.array(([2*np.pi,0,0])))
	bvec[1] = np.linalg.solve(avec,np.array(([0,2*np.pi,0])))
	bvec[2] = np.linalg.solve(avec,np.array(([0,0,2*np.pi])))
	return bvec


def path_to_plot(k_path,bvec,eq_points,point):
	k_path_tmp = []

	if eq_points:
		point = np.dot(point,bvec)/abs(bvec[0,0])
		x = 0.
		check = False
                for kn in range(k_path.shape[0]):
			k_path[kn] = np.dot(k_path[kn],bvec)/abs(bvec[0,0])
			if (k_path[kn-1] == point).all():
				check = True
                        if kn > 0 and not (k_path[kn-1] == point).all():      
				x = x + np.linalg.norm(k_path[kn]-k_path[kn-1])
                        k_path_tmp.append(x)
		if not check:	sys.exit('\nDid not find the expected equivalent points.\n')
	else:
		x = 0.
		for kn in range(k_path.shape[0]):
                        k_path[kn] = np.dot(k_path[kn],bvec)/abs(bvec[0,0])
			if kn > 0:	x = x + np.linalg.norm(k_path[kn]-k_path[kn-1])
			k_path_tmp.append(x)

	k_path = np.array(k_path_tmp)
	return k_path


def permute_hr(hr_pcell,permutation):

	h_tmp = np.zeros((hr_pcell.shape[1],hr_pcell.shape[2]),dtype=complex)

	for rn in range(hr_pcell.shape[0]):
		for i in range(len(permutation)):
			pi = permutation[i]
			for j in range(len(permutation)):
				pj = permutation[j]
				h_tmp[pi,pj] = hr_pcell[rn,i,j]
		hr_pcell[rn,:,:] = copy.copy(h_tmp)

	return hr_pcell
