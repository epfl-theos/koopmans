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
 PZ_CG_SPIN      PZ: Conjugate gradient w/ spin (outer loop only)
 PZ_CG_INN       PZ: Conjugate gradient (outer loop) + simple inner-loop
 PZ_CG_INN_SPIN  PZ: Conjugate gradient (outer loop) + simple inner-loop, spin polariz
 PZ_CG_ICG       PZ: Conjugate gradient (outer loop) + Conjugate gradient inner-loop
 PZ_DD_ICG       PZ: Damped dynamics (outer loop) + Conjugate gradient inner-loop

 LDA_DD          LDA:  Damped dynamics (outer loop only)
 LDA_CG          LDA:  Conjugate gradient (outer loop only)
 HF_DD           HF:   Damped dynamics (outer loop only)
 HF_CG           HF:   Conjugate gradient (outer loop only)
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

LDA_DD=
LDA_CG=
PZ_DD=
PZ_CG=
PZ_CG_SPIN=
PZ_CG_INN=
PZ_CG_INN_SPIN=
PZ_CG_ICG=
PZ_DD_ICG=
HF_DD=
HF_CG=
PBE0_DD=
PBE0_CG=
CHECK=
CLEAN=

if [ $# = 0 ] ; then echo "$MANUAL" ; exit 0 ; fi
INPUT=`echo $1 | tr [:upper:] [:lower:]`

case $INPUT in 
   (lda_dd)         LDA_DD=yes ;;
   (lda_cg)         LDA_CG=yes ;;

   (pz_dd)          PZ_DD=yes ;;
   (pz_cg)          PZ_CG=yes ;;
   (pz_cg_spin)     PZ_CG_SPIN=yes ;;
   (pz_cg_inn)      PZ_CG_INN=yes ;;
   (pz_cg_inn_spin) PZ_CG_INN_SPIN=yes ;;
   (pz_cg_icg)      PZ_CG_ICG=yes ;;
   (pz_dd_icg)      PZ_DD_ICG=yes ;;

   (hf_dd)          HF_DD=yes ;;
   (hf_cg)          HF_CG=yes ;;
   (pbe0_dd)        PBE0_DD=yes ;;
   (pbe0_cg)        PBE0_CG=yes ;;

   (odd)            PZ_DD=yes ; PZ_CG=yes ; PZ_CG_SPIN=yes ; 
                    PZ_CG_INN=yes ; PZ_CG_INN_SPIN=yes ; PZ_CG_ICG=yes ; PZ_DD_ICG=yes ;;
   (all)            PZ_DD=yes ; PZ_CG=yes ; PZ_CG_SPIN=yes ; 
                    PZ_CG_INN=yes ; PZ_CG_INN_SPIN=yes ; PZ_CG_ICG=yes ; PZ_DD_ICG=yes ;
                    HF_DD=yes ; HF_CG=yes ; PBE0_DD=yes    ; PBE0_CG=yes ;
                    LDA_DD=yes ; LDA_CG=yes ;;
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
run_cp  NAME=PZ_CG_SPIN   SUFFIX=$SUFFIX  RUN=$PZ_CG_SPIN
#
run_cp  NAME=PZ_CG_INNER      SUFFIX=$SUFFIX  RUN=$PZ_CG_INN
run_cp  NAME=PZ_CG_INNER_SPIN SUFFIX=$SUFFIX  RUN=$PZ_CG_INN_SPIN
run_cp  NAME=PZ_CG_INNERCG    SUFFIX=$SUFFIX  RUN=$PZ_CG_ICG
run_cp  NAME=PZ_DD_INNERCG    SUFFIX=$SUFFIX  RUN=$PZ_DD_ICG
#
run_cp  NAME=LDA_DD       SUFFIX=$SUFFIX  RUN=$LDA_DD
run_cp  NAME=LDA_CG       SUFFIX=$SUFFIX  RUN=$LDA_CG
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
       list=`ls -1 *out `
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

