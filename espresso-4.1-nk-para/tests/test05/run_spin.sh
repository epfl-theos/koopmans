#! /bin/bash
#
# Calculation of empty states (clorine mol, nspin=2)
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
 lda             LDA calculation (occupied states)
 lda_empty       LDA empty states
 lda_diag        LDA (occ+emp) by diagonalization
 nk0             NK0 calculation
 nk0_empty       NK0 empty states
 nk              NK calculation
 nk_empty        NK empty states
 pz              Perdew-Zunger SIC calculation
 pz_empty        Perdew-Zunger SIC empty states
 hf              Hartree-Fock calculation
 hf_empty        Hartree-Fock empty states

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
SUFFIX=_spin


#
# evaluate the starting choice about what is to run 

INIT=
GRAM=
LDA=
LDA_EMPTY=
LDA_DIAG=
NK0=
NK0_EMPTY=
NK=
NK_EMPTY=
PZ=
PZ_EMPTY=
HF=
HF_EMPTY=
CHECK=
CLEAN=

if [ $# = 0 ] ; then echo "$MANUAL" ; exit 0 ; fi
INPUT=`echo $1 | tr [:upper:] [:lower:]`

case $INPUT in 
   (init)           INIT=yes ;;
   (gram)           GRAM=yes ;;
   (lda)            LDA=yes ;;
   (lda_empty)      LDA_EMPTY=yes ;;
   (lda_diag)       LDA_DIAG=yes ;;
   (nk0)            NK0=yes ;;
   (nk0_empty)      NK0_EMPTY=yes ;;
   (nk)             NK=yes ;;
   (nk_empty)       NK_EMPTY=yes ;;
   (pz)             PZ=yes ;;
   (pz_empty)       PZ_EMPTY=yes ;;
   (hf)             HF=yes ;;
   (hf_empty)       HF_EMPTY=yes ;;

   (odd)            NK0=yes  ; NK=yes; PZ=yes; HF=yes ;
                    NK0_EMPTY=yes ; NK_EMPTY=yes; PZ_EMPTY=yes; HF_EMPTY=yes ;;
   (all)            INIT=yes ; GRAM=yes ; LDA=yes ; LDA_EMPTY=yes ; LDA_DIAG=yes ;
                    NK0=yes  ; NK=yes; PZ=yes; HF=yes ;
                    NK0_EMPTY=yes ; NK_EMPTY=yes; PZ_EMPTY=yes; HF_EMPTY=yes ;;
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
run_cp  NAME=INIT         SUFFIX=$SUFFIX  RUN=$INIT  PARALLEL=no
#
run_cp  NAME=GRAM         SUFFIX=$SUFFIX  RUN=$GRAM
#
run_cp  NAME=LDA          SUFFIX=$SUFFIX  RUN=$LDA
#
if [ "$LDA_EMPTY" = "yes" -a -e ./SCRATCH/clorine_spin.emp ] ; then
    rm -f ./SCRATCH/clorine_spin.emp 2> /dev/null
fi
#
run_cp  NAME=LDA_EMPTY    SUFFIX=$SUFFIX  RUN=$LDA_EMPTY
run_pw  NAME=LDA_DIAG     SUFFIX=$SUFFIX  RUN=$LDA_DIAG
#
run_cp  NAME=NK0          SUFFIX=$SUFFIX  RUN=$NK0
run_cp  NAME=NK0_EMPTY    SUFFIX=$SUFFIX  RUN=$NK0_EMPTY
#
run_cp  NAME=NK           SUFFIX=$SUFFIX  RUN=$NK 
run_cp  NAME=NK_EMPTY     SUFFIX=$SUFFIX  RUN=$NK_EMPTY 
#
run_cp  NAME=PZ           SUFFIX=$SUFFIX  RUN=$PZ
run_cp  NAME=PZ_EMPTY     SUFFIX=$SUFFIX  RUN=$PZ_EMPTY
#
run_cp  NAME=HF           SUFFIX=$SUFFIX  RUN=$HF
run_cp  NAME=HF_EMPTY     SUFFIX=$SUFFIX  RUN=$HF_EMPTY


#
# running CHECK
#
if [ "$CHECK" = yes ] ; then  
   echo "running CHECK... " 
   #
   cd $TEST_HOME
   cd Reference
       list=`ls -1 *_spin.out `
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

