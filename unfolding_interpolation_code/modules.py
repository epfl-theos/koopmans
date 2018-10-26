import numpy as np


# Definition of the Monkhorst-Pack k-mesh (in crystal coordinates) of the PC commensurate to the SC
def MP_mesh(nk1,nk2,nk3):
	k_mesh = []
	for i in range(nk1):
		for j in range(nk2):
			for k in range(nk3):
				k_mesh.append(np.array((i,j,k),dtype=float))
	k_mesh = np.array(k_mesh) / np.array((nk1,nk2,nk3))
	return k_mesh


# Calculate the nktot R-vectors in crystal coordinates
def Rvec(nk1,nk2,nk3):
	Rvec = []
	for i in range(nk1):
		for j in range(nk2):
			for k in range(nk3):
				Rvec.append(np.array((i,j,k)))
	Rvec = np.array(Rvec)
	return Rvec


# Compare two WFs in the SC and return True if they differ just by a primitive lattice vector (of the primitive cell)
def wann_match(centers,spreads,signatures,wann_ref,latt_vec,cutoff):
	dist = centers - wann_ref[:3]
	x = np.linalg.solve(latt_vec,dist)
	if np.linalg.norm(abs(x-np.round(x))) < cutoff:
		if abs(wann_ref[3]-spreads) < cutoff:
			if (abs(wann_ref[4:]-signatures)<cutoff*np.ones(20)).all():
				return True
	else:
		return False


def check_R(R,alat_pc,wann,wann_ref,latt_vec_sc):
	dist = wann - wann_ref
	comp = np.linalg.solve(latt_vec_sc,dist)
	for i in range(3):
		if abs(comp[i])>.5:
			R = R - comp[i]*latt_vec_sc[i]/alat_pc
	return R


def ws_distance(nk1,nk2,nk3,latt_vec,kvec,cutoff,wann,wann_ref):
	dist_min = np.linalg.norm(wann - wann_ref)
	Tnn = []		# List of degenerate T-vectors, corresponding to equidistant nn WFs
	for i in range(-1,2):
		for j in range(-1,2):
			for k in range(-1,2):
				T = i*latt_vec[0] + j*latt_vec[1] + k*latt_vec[2]	# T is a lattice vector of the SC
				dist = np.linalg.norm(wann - wann_ref + T)		# T-translation of the WF (evaluation of the copies of the given WF)
				if abs(dist - dist_min) < cutoff:			# The position (in PC crystal coordinates) of the copy of the given WF is accepted	
					Tnn.append(np.array((i*nk1,j*nk2,k*nk3),dtype=float))	
				elif dist < dist_min:					# A new nearest neighbor is found	
					dist_min = dist
					Tnn = [np.array((i*nk1,j*nk2,k*nk3),dtype=float)]
				else:
					continue
	
	Tnn = np.array(Tnn)
	phase_factor = np.sum(np.exp(1j*2*np.pi*np.dot(Tnn,kvec))) / len(Tnn)
	return phase_factor,Tnn
