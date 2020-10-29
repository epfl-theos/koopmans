import numpy as np
import argparse
import sys,re



parser = argparse.ArgumentParser()
parser.add_argument( "state", choices=['occ','emp'] )
state = parser.parse_args().state


# Read etot1 and lambda1
try:
    if ( state == 'occ' ): ifile = open( 'cp_pbe_N.out', 'r' )
    if ( state == 'emp' ): ifile = open( 'cp_pbe_N+1-1.out', 'r' )
except:
    ifile = open('001.out','r')

lines = ifile.readlines()
ifile.close()

for line in lines:
    if re.search('total energy =',line):
        etot1 = float(line.split()[3])
    if re.search('fixed_lambda',line):
        lambda1 = float(line.split()[3])


# Read lambda_2 and alpha0
try:
    if ( state == 'occ' ): ifile = open( 'cp_ki_N.out', 'r' )
    if ( state == 'emp' ): ifile = open( 'cp_ki_N+1-1.out', 'r' )
except:
    ifile = open('002.out','r')

lines = ifile.readlines()
ifile.close()

for line in lines:
    if re.search('fixed_lambda',line):
        lambda2 = float(line.split()[3])
    if re.search('NK scaling factor',line):
        alpha0 = float(line.split()[4])


# Read etot3
try:
    if ( state == 'occ' ): ifile = open( 'cp_pbe_N-1.out', 'r' )
    if ( state == 'emp' ): ifile = open( 'cp_pbe_N+1.out', 'r' )
except:
    ifile = open('003.out','r')

lines = ifile.readlines()
ifile.close()

for line in lines:
    if re.search('total energy =',line):
        etot3 = float(line.split()[3])

if ( state == 'occ' ):    alpha = alpha0 * ( etot1 - etot3 - lambda1 ) / ( lambda2 - lambda1 )
if ( state == 'emp' ):    alpha = alpha0 * ( etot3 - etot1 - lambda1 ) / ( lambda2 - lambda1 )

print('\n alpha = %.8f\n' %alpha)
