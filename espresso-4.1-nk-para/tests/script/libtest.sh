#! /bin/bash
#
# script to manage the lowlevel launch of codes
# needed by the test suite
#

# get script homedir
LOCAL_DIR=`echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname

# source environment
source $LOCAL_DIR/../environment.conf

# source base def
source $LOCAL_DIR/../script/basedef.sh

# set redirection
case $INPUT_TYPE in
 ("from_stdin")       INPUT_REDIRECT="<" ;;
 ("from_file")        INPUT_REDIRECT="-input" ;;
 (*)                  INPUT_REDIRECT="$INPUT_TYPE" ;;
esac

# few basic definitions
TEST_HOME=$(pwd)
TEST_NAME=$(echo $TEST_HOME | awk -v FS=\/ '{print $NF}' )


#
#----------------------
test_init () {
#----------------------
#
   
   #
   # exit if TMPDIR dies not exist
   if [ ! -d $TMPDIR ] ; then 
       echo "TMPDIR = $TMPDIR   does not exist " ; exit 71 
   fi

   #
   # if the case, create local test dir
   test -d $TMPDIR/$TEST_NAME || mkdir $TMPDIR/$TEST_NAME
   #
   # create SCRATCH link
   test -e $TEST_HOME/SCRATCH && rm $TEST_HOME/SCRATCH
   cd $TEST_HOME
   ln -sf $TMPDIR/$TEST_NAME ./SCRATCH
   #
   # create HOME link
   test -e $TMPDIR/$TEST_NAME/HOME && rm $TMPDIR/$TEST_NAME/HOME
   cd $TMPDIR/$TEST_NAME
   ln -sf $TEST_HOME ./HOME
   #
   test -e $TMPDIR/$TEST_NAME/CRASH && rm $TMPDIR/$TEST_NAME/CRASH
   test -e $TEST_HOME/CRASH && rm $TEST_HOME/CRASH
   #
   cd $TEST_HOME

}

#
#----------------------
exit_if_no_etsf_support () {
#----------------------
#
   make_sys_file="$TEST_HOME/../../make.sys"
   #
   if [ ! -e "$make_sys_file" ] ; then 
      echo "ERROR: make.sys not present" ; exit 10 
   fi
   #
   str=`grep '__ETSF_IO' $make_sys_file`
   #
   if [ -z "$str" ] ; then
      echo "no ETSF-IO support... exit "
      exit 0
   fi
}


#
#----------------------
run_clean () {
#----------------------
#
   local RUN=

   for arg
   do
         [[ "$arg" == RUN=* ]]       && RUN="${arg#RUN=}"
   done

   [[ "$RUN" != "yes" ]]  && return

   #
   # actual clean up
   #
   cd $TEST_HOME
      rm -rf *.out *.dat 2> /dev/null
      test -e SCRATCH && rm SCRATCH
      test -e CRASH   && rm CRASH

   cd $TMPDIR
      test -d $TEST_NAME && rm -rf $TEST_NAME
   
}


#
#----------------------
run () {
#----------------------
#
# low level tool to launch generic executables
#
   local NAME=
   local EXEC=
   local INPUT=
   local OUTPUT=
   local PARALLEL=
   local INPUT_TYPE_LOC=$INPUT_TYPE

   for arg 
   do
         [[ "$arg" == NAME=* ]]        && NAME="${arg#NAME=}"
         [[ "$arg" == EXEC=* ]]        && EXEC="${arg#EXEC=}"
         [[ "$arg" == INPUT=* ]]       && INPUT="${arg#INPUT=}"
         [[ "$arg" == OUTPUT=* ]]      && OUTPUT="${arg#OUTPUT=}"
         [[ "$arg" == PARALLEL=* ]]    && PARALLEL="${arg#PARALLEL=}"
         [[ "$arg" == INPUT_TYPE=* ]]  && INPUT_TYPE_LOC="${arg#INPUT_TYPE=}"
   done

   if [ -z "$NAME" ]   ; then echo "empty NAME card"   ; exit 1 ; fi 
   if [ -z "$EXEC" ]   ; then echo "empty EXEC card"   ; exit 1 ; fi 
   if [ -z "$INPUT" ]  ; then echo "empty INPUT card"  ; exit 1 ; fi 
   if [ -z "$OUTPUT" ] ; then echo "empty OUTPUT card" ; exit 1 ; fi 
   
   if [ ! -x $EXEC ] ; then
      #
      echo "$EXEC not executable... exit "
      exit 0
      #
   fi

   if [ ! -z $NAME ] ; then
      #
      echo $ECHO_N "running $NAME calculation... $ECHO_C"
      #
   fi

   #
   if [ "$INPUT_TYPE_LOC" = "from_stdin" ] ; then
      #
      if [ "$PARALLEL" = "yes" ] ; then
         $PARA_PREFIX $EXEC $PARA_POSTFIX < $INPUT > $OUTPUT
      else
         $EXEC < $INPUT > $OUTPUT
      fi
   fi
   #
   if [ "$INPUT_TYPE_LOC" = "from_file" ] ; then
      #
      if [ "$PARALLEL" = "yes" ] ; then
         $PARA_PREFIX $EXEC $PARA_POSTFIX -input $INPUT > $OUTPUT
      else
         $EXEC -input $INPUT > $OUTPUT
      fi
   fi
   #
   if [ $? = 0 ] ; then
      if [ ! -z $NAME ] ; then echo "$ECHO_T done" ; fi
   else
      echo "$ECHO_T problems found" ; exit 1
   fi

}


#
#----------------------
run_cp () {
#----------------------
#
   local NAME=CP
   local EXEC=$QE_BIN/cp.x
   local RUN=yes
   local INPUT=
   local OUTPUT=
   local SUFFIX=
   local PARALLEL=yes
   local name_tmp
   
   for arg 
   do
         [[ "$arg" == NAME=* ]]      && NAME="${arg#NAME=}"
         [[ "$arg" == INPUT=* ]]     && INPUT="${arg#INPUT=}"
         [[ "$arg" == OUTPUT=* ]]    && OUTPUT="${arg#OUTPUT=}"
         [[ "$arg" == SUFFIX=* ]]    && SUFFIX="${arg#SUFFIX=}"
         [[ "$arg" == RUN=* ]]       && RUN="${arg#RUN=}"
         [[ "$arg" == PARALLEL=* ]]  && PARALLEL="${arg#PARALLEL=}"
   done

   [[ "$RUN" != "yes" ]]  && return
   
   name_tmp=`echo $NAME | tr [:upper:] [:lower:]`
   if [ -z "$INPUT" ]  ; then  INPUT=$TEST_HOME/cp_$name_tmp$SUFFIX.in  ; fi
   if [ -z "$OUTPUT" ] ; then OUTPUT=$TEST_HOME/cp_$name_tmp$SUFFIX.out ; fi

   run NAME=$NAME INPUT=$INPUT OUTPUT=$OUTPUT EXEC=$EXEC PARALLEL=$PARALLEL
}


