import numpy as np
import sys
from datetime import datetime as dt


################################################
############### INPUT PARAMETERS ###############
################################################
#
kmesh = 4,4,4                               # k-mesh from the KC calculation
#
file_hr_dft_occ = "Si_hr.dat"			    # file with DFT hamiltonian for occupied states (expected to be in the Wannier90 format)
file_hr_dft_emp = "Si_emp_hr.dat"		    # file with DFT hamiltonian for empty states (expected to be in the Wannier90 format)
file_hr_kc_occ = "Si_KC_ham_occ.dat"		# file with KC potential for occupied states (expected to be in the Koopmans code format)
file_hr_kc_emp = "Si_KC_ham_emp.dat"		# file with KC potential for empty states (expected to be in the Koopmans code format)
#
only_occ_states = False		# if True it will calculate only the occupied states and ignore file_hr_dft_emp and file_hr_kc_emp
only_emp_states = False		# if True it will calculate only the empty states and ignore file_hr_dft_occ and file_hr_kc_occ
#
avec = np.array(([-2.7155,0.,2.7155],[0.,2.7155,2.7155],[-2.7155,2.7155,0.]))	# primitive vectors in Ang
#
# K-POINTS TO EVALUATE
k_path = np.array(([0.500000,0.250000,0.750000],[0.500000,0.263889,0.736111],[0.500000,0.277778,0.722222],[0.500000,0.291667,0.708333],[0.500000,0.305556,0.694444],[0.500000,0.319444,0.680556],[0.500000,0.333333,0.666667],[0.500000,0.347222,0.652778],[0.500000,0.361111,0.638889],[0.500000,0.375000,0.625000],[0.500000,0.388889,0.611111],[0.500000,0.402778,0.597222],[0.500000,0.416667,0.583333],[0.500000,0.430556,0.569444],[0.500000,0.444444,0.555556],[0.500000,0.458333,0.541667],[0.500000,0.472222,0.527778],[0.500000,0.486111,0.513889],[0.500000,0.500000,0.500000],[0.488636,0.488636,0.488636],[0.477273,0.477273,0.477273],[0.465909,0.465909,0.465909],[0.454545,0.454545,0.454545],[0.443182,0.443182,0.443182],[0.431818,0.431818,0.431818],[0.420455,0.420455,0.420455],[0.409091,0.409091,0.409091],[0.397727,0.397727,0.397727],[0.386364,0.386364,0.386364],[0.375000,0.375000,0.375000],[0.363636,0.363636,0.363636],[0.352273,0.352273,0.352273],[0.340909,0.340909,0.340909],[0.329545,0.329545,0.329545],[0.318182,0.318182,0.318182],[0.306818,0.306818,0.306818],[0.295455,0.295455,0.295455],[0.284091,0.284091,0.284091],[0.272727,0.272727,0.272727],[0.261364,0.261364,0.261364],[0.250000,0.250000,0.250000],[0.238636,0.238636,0.238636],[0.227273,0.227273,0.227273],[0.215909,0.215909,0.215909],[0.204545,0.204545,0.204545],[0.193182,0.193182,0.193182],[0.181818,0.181818,0.181818],[0.170455,0.170455,0.170455],[0.159091,0.159091,0.159091],[0.147727,0.147727,0.147727],[0.136364,0.136364,0.136364],[0.125000,0.125000,0.125000],[0.113636,0.113636,0.113636],[0.102273,0.102273,0.102273],[0.090909,0.090909,0.090909],[0.079545,0.079545,0.079545],[0.068182,0.068182,0.068182],[0.056818,0.056818,0.056818],[0.045455,0.045455,0.045455],[0.034091,0.034091,0.034091],[0.022727,0.022727,0.022727],[0.011364,0.011364,0.011364],[0.000000,0.000000,0.000000],[0.033333,0.033333,0.000000],[0.066667,0.066667,0.000000],[0.100000,0.100000,0.000000],[0.133333,0.133333,0.000000],[0.166667,0.166667,0.000000],[0.200000,0.200000,0.000000],[0.233333,0.233333,0.000000],[0.266667,0.266667,0.000000],[0.300000,0.300000,0.000000],[0.333333,0.333333,0.000000],[0.366667,0.366667,0.000000],[0.400000,0.400000,0.000000],[0.433333,0.433333,0.000000],[0.466667,0.466667,0.000000],[0.500000,0.500000,0.000000],[0.500000,0.487500,0.037500],[0.500000,0.475000,0.075000],[0.500000,0.462500,0.112500],[0.500000,0.450000,0.150000],[0.500000,0.437500,0.187500],[0.500000,0.425000,0.225000],[0.500000,0.412500,0.262500],[0.500000,0.400000,0.300000],[0.500000,0.387500,0.337500],[0.500000,0.375000,0.375000],[0.500000,0.362500,0.412500],[0.500000,0.350000,0.450000],[0.500000,0.337500,0.487500],[0.500000,0.325000,0.525000],[0.500000,0.312500,0.562500],[0.500000,0.300000,0.600000],[0.500000,0.287500,0.637500],[0.500000,0.275000,0.675000],[0.500000,0.262500,0.712500],[0.500000,0.250000,0.750000],[0.468750,0.281250,0.750000],[0.437500,0.312500,0.750000],[0.406250,0.343750,0.750000]))
#k_path = MP_mesh(4,4,4)		# Define automatically a Monkhorst-Pack mesh
#
################################################
############## END INPUT SECTION ###############
################################################
#
#
#
#
#
#
#
#
#
################################################
################### MODULES ####################
################################################

def read_hr(file_hr,file_type):

        ifile = open(file_hr,'r')
        lines = ifile.readlines()
        ifile.close()

	if file_type=='wannier':

	        num_wann = int(lines[1].split()[0])
        	nrpts = int(lines[2].split()[0])
	        weights = []
	        hr = []
        	Rvec = []

        	nlines_to_skip = 3 + nrpts/15
        	if nrpts%15 != 0:       nlines_to_skip += 1

        	for n in range(3,nlines_to_skip):
                	for m in range(len(lines[n].split())):
                        	weights.append(lines[n].split()[m])

	        for n in range(nlines_to_skip,len(lines),num_wann**2):
        	        Rvec.append(lines[n].split()[:3])

	        for n in range(nlines_to_skip,len(lines)):
        	        hr.append(float(lines[n].split()[5]) + 1.j*float(lines[n].split()[6]))

	        weights = np.array(weights,dtype=int)
        	Rvec = np.array(Rvec,dtype=int)
        	hr = np.array(hr).reshape(nrpts,num_wann,num_wann)

	        return hr,Rvec,weights

	if file_type=='koopmans':

	        num_wann = int(lines[0].split()[0])
        	nrpts = int(lines[1].split()[0])
	        hr = []
        	Rvec = []

	        for n in range(2,len(lines),num_wann**2):
        	        Rvec.append(lines[n].split()[:3])

	        for n in range(2,len(lines)):
        	        hr.append(float(lines[n].split()[5]) + 1.j*float(lines[n].split()[6]))

        	Rvec = np.array(Rvec,dtype=int)
        	hr = np.array(hr).reshape(nrpts,num_wann,num_wann)

	        return hr,Rvec


#def ws_distance(Rvec,kmesh):

        
def path_to_plot(k_path,avec):
        k_path_tmp = []
        bvec = np.zeros((3,3))
        bvec[0] = np.linalg.solve(avec,np.array(([2*np.pi,0,0])))
        bvec[1] = np.linalg.solve(avec,np.array(([0,2*np.pi,0])))
        bvec[2] = np.linalg.solve(avec,np.array(([0,0,2*np.pi])))
        x = 0.
        for kn in range(k_path.shape[0]):
       		k_path[kn] = np.dot(k_path[kn],bvec)/abs(bvec[0,0])
        	if kn > 0:      x = x + np.linalg.norm(k_path[kn]-k_path[kn-1])
       		k_path_tmp.append(x)

        k_path = np.array(k_path_tmp)
        return k_path


def write_output_eigk(eig_k,k_path,output):
        ofile = open(output,'w')
        ofile.write('# Written on %d-%d-%d at %d:%d:%02d\n' %(dt.now().day,dt.now().month,dt.now().year,dt.now().hour,dt.now().minute,dt.now().second))
        for n in range(eig_k.shape[1]):                 # Loop over the band index
                for kn in range(eig_k.shape[0]):        # Loop over the k-points (of the path)
                        ofile.write('\n%.4f\t%.6f' %(k_path[kn],eig_k[kn,n]))
                ofile.write('\n')
        ofile.close()


def MP_mesh(nr1,nr2,nr3):
        k_mesh = []
        for i in range(nr1):
                for j in range(nr2):
                        for k in range(nr3):
                                k_mesh.append(np.array((i,j,k),dtype=float))
        k_mesh = np.array(k_mesh) / np.array((nr1,nr2,nr3))
        return k_mesh


################################################
################# MAIN PROGRAM #################
################################################

if only_occ_states and not only_emp_states:
	hr_dft_occ,Rvec_dft,weights = read_hr(file_hr_dft_occ,'wannier')
	hr_kc_occ,Rvec_kc = read_hr(file_hr_kc_occ,'koopmans')
#
elif not only_occ_states and only_emp_states:
	hr_dft_emp,Rvec_dft,weights = read_hr(file_hr_dft_emp,'wannier')
	hr_kc_emp,Rvec_kc = read_hr(file_hr_kc_emp,'koopmans')
#
elif not only_occ_states and not only_emp_states:
	hr_dft_occ,Rvec_dft1,weights1 = read_hr(file_hr_dft_occ,'wannier')
	hr_dft_emp,Rvec_dft2,weights2 = read_hr(file_hr_dft_emp,'wannier')
	hr_kc_occ,Rvec_kc1 = read_hr(file_hr_kc_occ,'koopmans')
	hr_kc_emp,Rvec_kc2 = read_hr(file_hr_kc_emp,'koopmans')
	#
	if (Rvec_dft1 != Rvec_dft2).any():	sys.exit("\nR-vectors from DFT H(R) for empty and occupied states do not match!\n")
	else:					Rvec_dft = Rvec_dft1
	#
	if (weights1 != weights2).any():	sys.exit("\nweights from DFT H(R) for empty and occupied states do not match!\n")
	else:					weights = weights1
	#
	if (Rvec_kc1 != Rvec_kc2).any():	sys.exit("\nR-vectors from KC H(R) for empty and occupied states do not match!\n")
	else:					Rvec_kc = Rvec_kc1
#
else:
	sys.exit("\nCannot be True both only_occ_states and only_emp_states at the same time!\n")

hk_dft_occ = np.zeros((k_path.shape[0],hr_dft_occ.shape[2],hr_dft_occ.shape[2]),dtype=complex)
hk_dft_emp = np.zeros((k_path.shape[0],hr_dft_occ.shape[2],hr_dft_occ.shape[2]),dtype=complex)
hk_kc_occ = np.zeros((k_path.shape[0],hr_dft_occ.shape[2],hr_dft_occ.shape[2]),dtype=complex)
hk_kc_emp = np.zeros((k_path.shape[0],hr_dft_occ.shape[2],hr_dft_occ.shape[2]),dtype=complex)
eig_k_occ = []
eig_k_emp = []
kmesh = np.array(kmesh)

for kn in range(k_path.shape[0]):
	for m in range(hr_dft_occ.shape[1]):
		for n in range(hr_dft_occ.shape[2]):
			for i in range(Rvec_dft.shape[0]):
				if not only_emp_states:		hk_dft_occ[kn,m,n] = hk_dft_occ[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k_path[kn],Rvec_dft[i])) * hr_dft_occ[i,m,n] / weights[i]
				if not only_occ_states:		hk_dft_emp[kn,m,n] = hk_dft_emp[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k_path[kn],Rvec_dft[i])) * hr_dft_emp[i,m,n] / weights[i]
			for j in range(Rvec_kc.shape[0]):
#				phase_factor = ws_distance(Rvec_kc[j],kmesh)
				if not only_emp_states:		hk_kc_occ[kn,m,n] = hk_kc_occ[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k_path[kn],Rvec_kc[j])) * hr_kc_occ[j,m,n]
				if not only_occ_states:		hk_kc_emp[kn,m,n] = hk_kc_emp[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k_path[kn],Rvec_kc[j])) * hr_kc_emp[j,m,n]
#				if not only_emp_states:		hk_kc_occ[kn,m,n] = hk_kc_occ[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k_path[kn],Rvec_kc[j])) * phase_factor * hr_kc_occ[j,m,n]
#				if not only_occ_states:		hk_kc_emp[kn,m,n] = hk_kc_emp[kn,m,n] + np.exp(1j*2*np.pi*np.dot(k_path[kn],Rvec_kc[j])) * phase_factor * hr_kc_emp[j,m,n]

			if not only_emp_states:		hk_occ = hk_dft_occ + hk_kc_occ
			if not only_occ_states:		hk_emp = hk_dft_emp + hk_kc_emp

	if not only_emp_states:		eig_k_occ.append(np.linalg.eigvalsh(hk_occ[kn,:,:]))
	if not only_occ_states:		eig_k_emp.append(np.linalg.eigvalsh(hk_emp[kn,:,:]))

if not only_emp_states:		eig_k_occ = np.array(eig_k_occ)
if not only_occ_states:		eig_k_emp = np.array(eig_k_emp)

k_path = path_to_plot(k_path,avec)
if not only_emp_states:		write_output_eigk(eig_k_occ,k_path,'occ_bands.dat')
if not only_occ_states:		write_output_eigk(eig_k_emp,k_path,'emp_bands.dat')

################################################
############### END OF PROGRAM #################
################################################
