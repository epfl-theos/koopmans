#! /bin/awk -f

#
# parse_output.awk -- extract sensible information from the
#                     output files the CP-NK code
# 

BEGIN{ 
    status_error=0; 
    status_terminated=0; 
    program="" 

    iteration=0
    niter=-10
    etot=""
    ekinc=""

}

{ 
   #
   # first find the program which wrote the output
   #
   if ( NR < 15 ) 
   { 
      if ( match($0, "CP: variable-cell Car-Parrinello molecular dynamics") ) {
         program="cp";
      } else if ( match($0, "Program PWSCF") ) {
         program="pw";
      }
   }

   #
   # select the suitable parsing to do
   #
   if ( program == "cp" ) 
   {
       parse_cp();
   }
   else if ( program == "pw" )
   {
       parse_pw();
   }

}

END{ 
   if ( program = "cp") {
       #
       print "ITERATION@"iteration"@1e-3";
       print "ETOT@"etot"@1e-3";
       print "EKINC@"ekinc"@1e-3";
       #
   } else if ( program = "pw") {
       # 
   }
   #
   if ( status_error ) {
       print "STATUS@ERROR";
   } else if ( status_terminated ) { 
       print "STATUS@TERMINATED";
   } else {
       print "STATUS@UNKNOWN";
   }
}

function parse_cp()
{
   #
   # first, check whether 
   # the calculation is converged or not
   #
   if ( match($0, ": error #") ) {
      status_error=1;
   }
   if ( match($0, "This run was terminated") ) {
      status_terminated=1;
   }
   
   #
   # now, perform full check
   #
   if ( ! staus_error ) {
      check_line_cp();
   }
}

function check_line_cp()
{
   #
   # for each "sensible" value found, 
   # print KEYWORK @ value @ tollerance   (without blanks)
   #
   if ( match($0, "nfi    ekinc  temph  tempp        etot      enthal") ) 
      {
         nline=NR
      }
   else if ( match($0, "* Physical Quantities at step:") )
      {
         nline=NR
         iteration=$NF
      }
   else if ( NR == nline+1 && NF > 5 )
      {
         etot=$5
         ekinc=$2
         nline=-10
   }
}
   

function parse_pw()
{
   #
   # first, check whether
   # the calculation is converged or not
   #
   if ( match($0, ": error #") ) {
      status_error=1;
   }
   if ( match($0, "This run was terminated") ) {
      status_terminated=1;
   }

   #
   # now, perform full check
   #
   if ( ! staus_error ) {
      check_line_pw();
   }
}

function check_line_pw()
{

}

