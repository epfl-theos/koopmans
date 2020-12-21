#! /bin/bash
#
# Benzene molecule
# 
#================================================================
#
# Input flags for this script (./run.sh FLAG): 
#
MANUAL=" Usage
   run.sh [FLAG]

 where FLAG is one of the following 
 (no FLAG will print this manual page) :
 
 init            generates the random wfcs using 1proc
 gram            performs Gram-Schmidt orthogonalization steps
 lda             LDA calculation
 nk0             NK0 calculation
 nk              NK calculation
 pz              Perdew-Zunger SIC calculation
 hf              Hartree-Fock calculation

 odd             performs the orbital dependent density (ODD) calculations
 all             perform all the above described steps

 check           check results with the reference outputs
 clean           delete all output files and the temporary directory
"
#
#================================================================
#

#
# source common enviroment
. ../environment.conf
#
# source low level macros for test
. ../script/libtest.sh

#
# macros
SUFFIX=


#
# evaluate the starting choice about what is to run 

INIT=
GRAM=
LDA=
NK0=
NK=
PZ=
HF=
CHECK=
CLEAN=

if [ $# = 0 ] ; then echo "$MANUAL" ; exit 0 ; fi
INPUT=`echo $1 | tr [:upper:] [:lower:]`

case $INPUT in 
   (init)           INIT=yes ;;
   (gram)           GRAM=yes ;;
   (lda)            LDA=yes ;;
   (nk0)            NK0=yes ;;
   (nk)             NK=yes ;;
   (pz)             PZ=yes ;;
   (hf)             HF=yes ;;

   (odd)            NK0=yes  ; NK=yes; PZ=yes; HF=yes ;;
   (all)            INIT=yes ; GRAM=yes ; LDA=yes ;
                    NK0=yes  ; NK=yes; PZ=yes; HF=yes ;;
   (check)          CHECK=yes ;;
   (clean)          CLEAN=yes ;;
   (*)              echo " Invalid input FLAG, type ./run.sh for help" ; exit 1 ;;
esac


#
# initialize
#
if [ -z "$CLEAN" ] ; then
   test_init
fi
#


#-----------------------------------------------------------------------------

#
# running calculations
#
run_cp  NAME=INIT   SUFFIX=$SUFFIX  RUN=$INIT  PARALLEL=no
#
run_cp  NAME=GRAM   SUFFIX=$SUFFIX  RUN=$GRAM
#
run_cp  NAME=LDA    SUFFIX=$SUFFIX  RUN=$LDA
#
run_cp  NAME=NK0    SUFFIX=$SUFFIX  RUN=$NK0
#
run_cp  NAME=NK     SUFFIX=$SUFFIX  RUN=$NK 
#
run_cp  NAME=PZ     SUFFIX=$SUFFIX  RUN=$PZ
#
run_cp  NAME=HF     SUFFIX=$SUFFIX  RUN=$HF


#
# running CHECK
#
if [ "$CHECK" = yes ] ; then  
   echo "running CHECK... " 
   #
   cd $TEST_HOME
   cd Reference
       list=`ls -1 *out | grep -v _US`
   cd ..
   #
   for file in $list
   do
      ../script/check.sh $file
   done
fi


#
# eventually clean
#
run_clean  RUN=$CLEAN


#
# exiting
exit 0

