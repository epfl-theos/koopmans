#!/bin/bash
# Concatenates all of the .cpo files from an entire workflow

# Find all of the files, sorting by time, and excluding copies of init files in calc_alpha/orbital_? directories
files=$(find ./ \( -name '*.cp?' -o -name '*.pw?' -o -name '*.p2w?' -o -name '*.win' -o -name '*wout' -o -name '*.ui?' -o -name '*.w2k?' -o -name '*.ks?' -o -name '*.kh?' \) -exec stat -c '%.9Y %n' {} + | sort -n | grep -v "alpha.*init" | awk '{print $2}')

for f in $files
do 
   echo Contents of $f:
   cat $f
   echo ''
done
