import numpy as np
import sys,re


try:
	state = sys.argv[1]
except:
	sys.exit('\nYou have to indicate which kind of state: occ or emp ?\n')


# Read etot1 and lambda1
ifile = open('001.out','r')
lines = ifile.readlines()
ifile.close()

for line in lines:
	if re.search('total energy =',line):
		etot1 = float(line.split()[3])
	if re.search('fixed_lambda',line):
		lambda1 = float(line.split()[3])


# Read lambda_2 and alpha0
ifile = open('002.out','r')
lines = ifile.readlines()
ifile.close()

for line in lines:
	if re.search('fixed_lambda',line):
		lambda2 = float(line.split()[3])
	if re.search('NK scaling factor',line):
		alpha0 = float(line.split()[4])


# Read etot3
ifile = open('003.out','r')
lines = ifile.readlines()
ifile.close()

for line in lines:
	if re.search('total energy =',line):
		etot3 = float(line.split()[3])

if state=='occ':	alpha = alpha0 * ( etot1 - etot3 - lambda1 ) / ( lambda2 - lambda1 )
if state=='emp':	alpha = alpha0 * ( etot3 - etot1 - lambda1 ) / ( lambda2 - lambda1 )

print('\n alpha = %.8f\n' %alpha)
