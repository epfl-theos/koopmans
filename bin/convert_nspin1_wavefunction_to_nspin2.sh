#!/bin/bash

nspin1dir=$1
nspin2dir=$2

if [ ! -d "$nspin1dir" ]; then
   echo "$nspin1dir does not exist"
   exit 1
fi

if [ ! -d "$nspin2dir" ]; then
   echo "$nspin2dir does not exist"
   exit 1
fi

for f in evc0.dat evc0_empty1.dat evcm.dat evc.dat evcm.dat hamiltonian.xml eigenval.xml evc_empty1.dat lambda01.dat
do
   if [[ $f == *1.* ]]; then
      prefix=${f%1.*}
      suffix=${f#*1.}
   else
      prefix=${f%.*}
      suffix=${f#*.}
   fi

   f_in=$nspin1dir/$f

   if [ -f "$f_in" ]; then
      for i in 1 2
      do
         f_out=$nspin2dir/$prefix$i"."$suffix
         cp -r $f_in $f_out
         sed -i 's/nk="1"/nk="2"/g' $f_out
         sed -i 's/nspin="1"/nspin="2"/g' $f_out
      done
      sed -i 's/ik="1"/ik="2"/g' $f_out
      sed -i 's/ispin="1"/ispin="2"/g' $f_out
   fi
done
