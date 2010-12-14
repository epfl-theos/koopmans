#! /bin/bash
#
# LiH molecule
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
 nk              NK calculation

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
SUFFIX=_LiH


#
# evaluate the starting choice about what is to run 

INIT=
GRAM=
NK=
CHECK=
CLEAN=

if [ $# = 0 ] ; then echo "$MANUAL" ; exit 0 ; fi
INPUT=`echo $1 | tr [:upper:] [:lower:]`

case $INPUT in 
   (init)           INIT=yes ;;
   (gram)           GRAM=yes ;;
   (nk)             NK=yes ;;

   (all)            INIT=yes ; GRAM=yes ; NK=yes ;;
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
run_cp  NAME=NK     SUFFIX=$SUFFIX  RUN=$NK 


#
# running CHECK
#
if [ "$CHECK" = yes ] ; then  
   echo "running CHECK... " 
   #
   cd $TEST_HOME
   cd Reference
       list=`ls *_LiH.out`
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

