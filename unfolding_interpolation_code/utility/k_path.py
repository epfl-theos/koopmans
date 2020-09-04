import numpy as np


numk = input('\nHow many k-points define the path? ')
k = np.zeros((numk,3))
length = np.zeros(numk-1)	# length of each portion of path
point_name = []

print '\nInsert the coordinates (separated by a comma) of each k-vector in crystal units:'
for n in range(numk):
	k[n] = input('k%d = ' %(n+1))
	point_name.append(str(raw_input('Point name (optional): ')))
	if n>0:	length[n-1] = np.linalg.norm(k[n]-k[n-1])

nktot = input('\nHow many k-points along the path? ')

Nij = np.zeros(numk-1)
for n in range(1,numk-1):
	Nij[n] = length[n]/length[n-1]

Nij[0] = float(nktot) / (1 + np.sum(Nij[1:]))
Nij[1:] = Nij[1:]*Nij[0]
N = np.array(np.round(Nij),dtype=int)
Ntot = np.sum(N)

kpath = []

counter = 0
pos_input_points = np.zeros(numk,dtype=int)
for p in range(numk-1):
	direction = (k[p+1] - k[p]) / N[p]
	for n in range(N[p]):
		kpath.append(k[p]+n*direction)
		counter += 1
	pos_input_points[p+1] = counter
kpath.append(k[numk-1])

kpath = np.array(kpath)
print kpath

ofile = open('kpath.txt','w')
for n in range(kpath.shape[0]):
	ofile.write('[%f,%f,%f],' %(kpath[n,0],kpath[n,1],kpath[n,2]))
ofile.close()


ofile = open('kpath_qe.txt','w')
ofile.write('K_POINTS crystal\n %d\n' %kpath.shape[0])
weight = 1./kpath.shape[0]
counter = 0
for n in range(kpath.shape[0]):
	if n==pos_input_points[counter]:
		ofile.write('  %f  %f  %f  %e\t! %s\n' %(kpath[n,0],kpath[n,1],kpath[n,2],weight,point_name[counter]))
		counter += 1
	else:
		ofile.write('  %f  %f  %f  %e\n' %(kpath[n,0],kpath[n,1],kpath[n,2],weight))
ofile.close()
