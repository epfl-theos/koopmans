import numpy as np
import sys,re
from warnings import warn



BOHR_RADIUS=0.52917721067

try:
	ifile = open(sys.argv[1],'r')
except:
	sys.exit('\nCannot find input file\n')


lines = ifile.readlines()
ifile.close()

line_apos = len(lines)
atom_type = []
apos = []
cell = np.zeros((3,3))
celldm = False
nspin = 1

for n in range(len(lines)):
    if re.search('ibrav',lines[n]):
        ibrav = int(lines[n].split()[2])
        if ibrav > 3:
            sys.exit('\nWe can deal only with ibrav=0 and cubic systems (ibrav=1,2,3).\n')
	
    if re.search('nat',lines[n]):
        nat = int(lines[n].split()[2])

    if re.search('ATOMIC_POSITIONS',lines[n]):
        apos_units = lines[n].split()[1].replace('{','').replace('}','')
        line_apos = n + 1

    if n>=line_apos and n<line_apos+nat:
        atom_type.append(lines[n].split()[0])
        apos.append(np.array(lines[n].split()[1:4],dtype=float))

    if re.search('celldm\(1\)',lines[n]):
        celldm = True
        alat = float(lines[n].split()[2])
	
    if re.search('CELL_PARAMETERS',lines[n]):
        cell_units = lines[n].split()[1].replace('{','').replace('}','')
        cell[0] = np.array(lines[n+1].split(),dtype=float)
        cell[1] = np.array(lines[n+2].split(),dtype=float)
        cell[2] = np.array(lines[n+3].split(),dtype=float)
	
    if re.search('K_POINTS',lines[n]):
        if lines[n].split()[1]!='automatic' and lines[n].split()[1]!='{automatic}':
            sys.exit('\nList of kpoints must be automatic.\n')
        kmesh = np.array(lines[n+1].split()[:3],dtype=int)
  
    if re.search('nbnd',lines[n]):
        nbnd = int(lines[n].split()[2])

    if re.search('nspin',lines[n]):
        nspin = int(lines[n].split()[2])
 
    if re.search('OCCUPATIONS',lines[n]):
        try:
            nbnd
        except NameError:
            warn('nbnd is needed to parse OCCUPATIONS. Skipping parsing of OCCUPATIONS')

        occ_line = n+1
        counter = 0
        occup = []
        while counter < nbnd*nspin:
            counter += len(lines[occ_line].split())
            for ii in range(len(lines[occ_line].split())):
                occup.append(float(lines[occ_line].split()[ii]))
            occ_line += 1

				
if not celldm and (cell_units=='alat' or apos_units=='alat'):
    sys.exit('\nalat coordinates need celldm(1).\n')


try:
    apos = np.array(apos)
except:
    sys.exit('\nCannot find atomic positions!')

try:
    ibrav
except:
    sys.exit('\nCannot find ibrav!\n')

try:
    nat
except:
    sys.exit('\nCannot find nat!\n')

try:
    kmesh
except:
    sys.exit('\nCannot find kpoints!')

if ibrav == 0:
    try:
        cell
    except:
        sys.exit('\nibrav=0 needs CELL_PARAMETERS card.\n')


# HERE WE BUILD THE SUPERCELL
if ibrav==1:
    cell_sc = np.identity(3)*kmesh
    cell_units = 'alat'
if ibrav==2:	# We use PW convention
    cell = np.array(([-0.5, 0.0, 0.5],
                     [ 0.0, 0.5, 0.5],
                     [-0.5, 0.5, 0.0]))
    cell_units = 'alat'
if ibrav==3:	# We use PW convention
    cell = np.array(([ 0.5, 0.5, 0.5],
                     [-0.5, 0.5, 0.5],
                     [-0.5,-0.5, 0.5]))
    cell_units = 'alat'

cell_sc = (cell.transpose()*kmesh).transpose()


# The new atomic positions are always calculated in crystal units
apos_sc = np.zeros((np.prod(kmesh)*nat,3))
if apos_units=='alat':
    if cell_units=='alat':
        apos = np.linalg.solve(cell_sc,apos.transpose()).transpose()
    if cell_units=='angstrom':
        apos = np.linalg.solve(cell_sc/BOHR_RADIUS,apos.transpose()*celldm).transpose()
    if cell_units=='bohr':
        apos = np.linalg.solve(cell_sc,apos.transpose()*celldm).transpose()

if apos_units=='bohr':
    if cell_units=='alat':
        apos = np.linalg.solve(cell_sc*celldm,apos.transpose()).transpose()
    if cell_units=='angstrom':
        apos = np.linalg.solve(cell_sc/BOHR_RADIUS,apos.transpose()).transpose()
    if cell_units=='bohr':
        apos = np.linalg.solve(cell_sc,apos.transpose()).transpose()

if apos_units=='angstrom':
    if cell_units=='alat':
        apos = np.linalg.solve(cell_sc*celldm*BOHR_RADIUS,apos.transpose()).transpose()
    if cell_units=='angstrom':
        apos = np.linalg.solve(cell_sc,apos.transpose()).transpose()
    if cell_units=='bohr':
        apos = np.linalg.solve(cell_sc*BOHR_RADIUS,apos.transpose()).transpose()

if apos_units=='crystal':
    apos[:,0] = apos[:,0] / kmesh[0]
    apos[:,1] = apos[:,1] / kmesh[1]
    apos[:,2] = apos[:,2] / kmesh[2]

for nk1 in range(kmesh[0]):
    for nk2 in range(kmesh[1]):
        for nk3 in range(kmesh[2]):
            counter = nk1*kmesh[1]*kmesh[2] + nk2*kmesh[2] + nk3
            apos_sc[counter*nat:(counter+1)*nat] = apos + \
            np.array((float(nk1)/kmesh[0],float(nk2)/kmesh[1],float(nk3)/kmesh[2]))

try:
    occup
    nbnd_sc = nbnd * np.prod(kmesh)
    occup_sc = [1.0] * (nbnd_sc-nbnd) + occup[:nbnd]
    if nspin>1:
        occup_sc += [1.0] * (nbnd_sc-nbnd) + occup[nbnd:]
except NameError:
    pass
    


# WRITING THE OUTPUT
ofile = open('supercell.txt','w')
ofile.write('nat = %d\n' %(nat*np.prod(kmesh)))

try:
    nbnd_sc
    ofile.write('nbnd = %d\n' %nbnd_sc)
except NameError:
    pass

ofile.write('\n----------------------- OPTION #1 -----------------------\n')
ofile.write('ibrav = 0\n')
if cell_units=='alat':
    ofile.write('celldm(1) = %f' %alat)

ofile.write('\nCELL_PARAMETERS %s\n' %cell_units)
for i in range(3):
    ofile.write(' %12.8f  %12.8f  %12.8f\n' %(cell_sc[i,0],cell_sc[i,1],cell_sc[i,2]))
ofile.write('---------------------------------------------------------\n')

if ibrav!=0 and kmesh[0]==kmesh[1]==kmesh[2]:
    ofile.write('\n----------------------- OPTION #2 -----------------------\n')
    ofile.write('ibrav = %d\n' %ibrav)
    ofile.write('celldm(1) = %f\n' %(alat*kmesh[0]))
    ofile.write('no CELL_PARAMETERS\n')
    ofile.write('---------------------------------------------------------\n')

ofile.write('\nATOMIC_POSITIONS crystal\n')
for n in range(apos_sc.shape[0]):
    ofile.write(' %s  %.8f  %.8f  %.8f\n' %(atom_type[n%nat],apos_sc[n,0],apos_sc[n,1],apos_sc[n,2]))

ofile.write('\nK_POINTS gamma\n')

try:
    occup_sc
    ofile.write('\nOCCUPATIONS\n')
    for s in range(nspin):
        occup_sc_ = occup_sc[nbnd_sc*s:nbnd_sc*(s+1)]
        remaining_occ = len(occup_sc_)
        while remaining_occ>=8:
            for occ in occup_sc_[nbnd_sc-remaining_occ:nbnd_sc-remaining_occ+8]:
                ofile.write(' %4.2f' %occ)
            ofile.write('\n')
            remaining_occ -= 8
        if remaining_occ>0:
            for occ in occup_sc_[nbnd_sc-remaining_occ:]:
                ofile.write(' %4.2f' %occ)
            ofile.write('\n')
except NameError:
    pass

ofile.close()
