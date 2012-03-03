#! /bin/bash
#
# Minimization pf PZ functional using CG and inner-loop techniques
# 
#================================================================
#
# Input flags for this script (./run.sh FLAG): 
#
MANUAL=" Usage
   run.sh [FLAG]

 where FLAG is one of the following 
 (no FLAG will print this manual page) :
 
 PZ_DD           PZ: Damped dynamics (outer loop only)
 PZ_CG           PZ: Conjugate gradient (outer loop only)
 PZ_CG_ICG       PZ: Conjugate gradient (outer loop) + Conjugate gradient inner-loop
 PZ_DD_ICG       PZ: Damped dynamics (outer loop) + Conjugate gradient inner-loop

 HF_DD           HF: Damped dynamics (outer loop only)
 HF_CG           HF: Conjugate gradient (outer loop only)
 PBE0_DD         PBE0: Damped dynamics (outer loop only)
 PBE0_CG         PBE0: Conjugate gradient (outer loop only)

 odd             performs all ODD (PZ) calculations
 all             performs all the above

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

PZ_DD=
PZ_CG=
HF_DD=
HF_CG=
PBE0_DD=
PBE0_CG=
PZ_CG_ICG=
PZ_DD_ICG=
CHECK=
CLEAN=

if [ $# = 0 ] ; then echo "$MANUAL" ; exit 0 ; fi
INPUT=`echo $1 | tr [:upper:] [:lower:]`

case $INPUT in 
   (pz_dd)          PZ_DD=yes ;;
   (pz_cg)          PZ_CG=yes ;;
   (pz_cg_icg)      PZ_CG_ICG=yes ;;
   (pz_dd_icg)      PZ_DD_ICG=yes ;;

   (hf_dd)          HF_DD=yes ;;
   (hf_cg)          HF_CG=yes ;;
   (pbe0_dd)        PBE0_DD=yes ;;
   (pbe0_cg)        PBE0_CG=yes ;;

   (odd)            PZ_DD=yes ; PZ_CG=yes ; PZ_CG_ICG=yes ; PZ_DD_ICG=yes ;;
   (all)            PZ_DD=yes ; PZ_CG=yes ; PZ_CG_ICG=yes ; PZ_DD_ICG=yes ;
                    HF_DD=yes ; HF_CG=yes ; PBE0_DD=yes   ; PBE0_CG=yes ;;
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
run_cp  NAME=PZ_DD        SUFFIX=$SUFFIX  RUN=$PZ_DD
run_cp  NAME=PZ_CG        SUFFIX=$SUFFIX  RUN=$PZ_CG
#
run_cp  NAME=PZ_DD_ICG    SUFFIX=$SUFFIX  RUN=$PZ_DD_ICG
run_cp  NAME=PZ_CG_ICG    SUFFIX=$SUFFIX  RUN=$PZ_CG_ICG
#
run_cp  NAME=HF_DD        SUFFIX=$SUFFIX  RUN=$HF_DD
run_cp  NAME=HF_CG        SUFFIX=$SUFFIX  RUN=$HF_CG
#
run_cp  NAME=PBE0_DD      SUFFIX=$SUFFIX  RUN=$PBE0_DD
run_cp  NAME=PBE0_CG      SUFFIX=$SUFFIX  RUN=$PBE0_CG


#
# running CHECK
#
if [ "$CHECK" = yes ] ; then  
   echo "running CHECK... " 
   #
   cd $TEST_HOME
   cd Reference
       list=`ls -1 *out | grep -v _occ.out | grep -v _spin.out `
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

