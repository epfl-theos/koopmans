#!/bin/bash

nspin2dir=$1
nspin1dir=$2

if [ ! -d "$nspin1dir" ]; then
   echo "$nspin1dir does not exist"
   exit 1
fi

if [ ! -d "$nspin2dir" ]; then
   echo "$nspin2dir does not exist"
   exit 1
fi

for f in evc0.dat evc0_empty1.dat evcm.dat evc.dat evcm.dat hamiltonian.xml eigenval.xml evc_empty1.dat lambda01.dat lambdam1.dat
do
   if [[ $f == *1.* ]]; then
      prefix=${f%1.*}
      suffix=${f#*1.}
   else
      prefix=${f%.*}
      suffix=${f#*.}
   fi

   f_out=$nspin1dir/$f
   f_in=$nspin2dir/$prefix"1."$suffix

   if [ -f "$f_in" ]; then
      cp -r $f_in $f_out 2> /dev/null
      sed -i 's/nk="2"/nk="1"/g' $f_out
      sed -i 's/nspin="2"/nspin="1"/g' $f_out
   fi

   for i in 1 2
   do
      to_delete=$nspin2dir/$prefix$i"."$suffix
      if [ $to_delete != $f_out ]; then
         rm $to_delete 2> /dev/null
      fi
   done
done
