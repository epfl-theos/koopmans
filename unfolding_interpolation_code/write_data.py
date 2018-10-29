import numpy as np
from datetime import datetime as dt


# Write eigenvalues E_k.
# In case of interpolation the output file is called bands_interpolated.dat
# With no interpolation the ouput file is called eigk.dat
def write_output_eigk(eig_k,order,*k_path):

	if order == 'band':		# Group eigenvalues with same band index
		ofile = open('bands_interpolated.dat','w')
		ofile.write('# Written on %d-%d-%d at %d:%d:%02d\n' %(dt.now().day,dt.now().month,dt.now().year,dt.now().hour,dt.now().minute,dt.now().second))
		k_path = np.array(k_path)
		k_path = k_path.reshape(k_path.shape[1],k_path.shape[2])	
		for n in range(eig_k.shape[1]):
			x = 0.
			for kn in range(eig_k.shape[0]):
				if kn > 0:
					x = x + np.linalg.norm(k_path[kn]-k_path[kn-1])
				ofile.write('\n%.4f\t%.6f' %(x,eig_k[kn,n]))
			ofile.write('\n')
		ofile.close()
	
	if order == 'kvec':		# Group eigenvalues with same k-vector
		ofile = open('eigk.dat','w')
		ofile.write('# Written on %d-%d-%d at %d:%d:%02d\n' %(dt.now().day,dt.now().month,dt.now().year,dt.now().hour,dt.now().minute,dt.now().second))
		for kn in range(eig_k.shape[0]):
			for n in range(eig_k.shape[1]):
				ofile.write('\n%.6f' %eig_k[kn,n])
			ofile.write('\n')
		ofile.close()


def write_dos(gaussian,Emin,deltaE):
	ofile = open('dos_interpolated.dat','w')
	ofile.write('# Written on %d-%d-%d at %d:%d:%02d\n' %(dt.now().day,dt.now().month,dt.now().year,dt.now().hour,dt.now().minute,dt.now().second))
	ene = Emin
	for n in range(len(gaussian)):
		ene = ene + deltaE
		ofile.write('\n%.4f\t%.6E' %(ene,gaussian[n]))
	ofile.close()
