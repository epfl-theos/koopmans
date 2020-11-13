import numpy as np
import argparse
import sys,re



parser = argparse.ArgumentParser()
parser.add_argument( "state", choices=['occ','emp'] )
parser.add_argument( "orb", type=int, nargs='?', default=False )
args = parser.parse_args()


ki_calc = False
kipz_calc = False

# Read etot1 and lambda1
try:

    if ( args.state == 'occ' ):
        try:
            ifile = open( 'cp_pbe_N.out', 'r' )
            ki_calc = True
        except FileNotFoundError:
            ifile = open( 'cp_kipz_N.out', 'r' )
            kipz_calc = True

    if ( args.state == 'emp' ):
        try:
            ifile = open( 'cp_pbe_N+1-1.out', 'r' )
            ki_calc = True
        except FileNotFoundError:
            ifile = open( 'cp_kipz_N+1-1.out', 'r' )
            kipz_calc = True

except:
    ifile = open('001.out','r')

lines = ifile.readlines()
ifile.close()

for line in lines:
    if re.search('total energy =',line):
        etot1 = float(line.split()[3])
    if re.search('fixed_lambda',line):
        lambda1 = float(line.split()[3])


odd_alpha = False
if ( ki_calc ):
    if ( args.state == 'occ' ): ifile = open( 'cp_ki_N.in', 'r' )
    if ( args.state == 'emp' ): ifile = open( 'cp_ki_N+1-1.in', 'r' )
    lines = ifile.readlines()
    for line in lines:
        if ( re.search('odd_nkscalfact',line) ):
            if ( 'true' in line ): odd_alpha = True
        if ( re.search('fixed_band', line) ):
            orb = int(line.split()[2])
        if ( args.state == 'emp' and re.search('nelec',line) ):
            nelec = int(line.split()[2])

if ( kipz_calc ):
    if ( args.state == 'occ' ): ifile = open( 'cp_kipz_N.in', 'r' )
    if ( args.state == 'emp' ): ifile = open( 'cp_kipz_N+1-1.in', 'r' )
    lines = ifile.readlines()
    for line in lines:
        if ( re.search('odd_nkscalfact',line) ):
            if ( 'true' in line ): odd_alpha = True

if ( odd_alpha ):
    if ( args.state == 'occ' ): ifile = open( 'file_alpharef.txt', 'r' )
    if ( args.state == 'emp' ): ifile = open( 'file_alpharef_empty.txt', 'r' )
    ### TO FINISH ...




# Read lambda_2 and alpha0
try:

    if ( args.state == 'occ' ):
        try:
            ifile = open( 'cp_ki_N.out', 'r' )
        except FileNotFoundError:
            ifile = open( 'cp_pbe_N.out', 'r' )

    if ( args.state == 'emp' ):
        try:
            ifile = open( 'cp_ki_N+1-1.out', 'r' )
        except FileNotFoundError:
            ifile = open( 'cp_pbe_N+1-1.out', 'r' )

except:
    ifile = open('002.out','r')

lines = ifile.readlines()
ifile.close()

for line in lines:
    if ( re.search('fixed_lambda',line) ):
        lambda2 = float(line.split()[3])
    if ( not args.orb and re.search('NK scaling factor',line) ):
        alpha0 = float(line.split()[4])


# Read etot3
try:

    if ( args.state == 'occ' ):
        try:
            ifile = open( 'cp_pbe_N-1.out', 'r' )
        except FileNotFoundError:
            ifile = open( 'cp_kipz_N-1.out', 'r' )

    if ( args.state == 'emp' ):
        try:
            ifile = open( 'cp_pbe_N+1.out', 'r' )
        except FileNotFoundError:
            ifile = open( 'cp_kipz_N+1.out', 'r' )

except:
    ifile = open('003.out','r')

lines = ifile.readlines()
ifile.close()

for line in lines:
    if re.search('total energy =',line):
        etot3 = float(line.split()[3])

if ( args.state == 'occ' ):    alpha = alpha0 * ( etot1 - etot3 - lambda1 ) / ( lambda2 - lambda1 )
if ( args.state == 'emp' ):    alpha = alpha0 * ( etot3 - etot1 - lambda1 ) / ( lambda2 - lambda1 )

print('\n alpha = %.8f\n' %alpha)
