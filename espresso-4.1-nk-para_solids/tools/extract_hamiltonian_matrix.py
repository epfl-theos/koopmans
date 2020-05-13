#!/usr/bin/env python

# Needed modules

import os
import sys
import time
import numpy
import optparse
import numpy

# Main function

def main():
  #
  # parsing command-line arguments
  parser = optparse.OptionParser()
  #
  parser.add_option('-p','--path',
                    dest="path",
                    default="./",
                    help="path to the hamiltonian file")
  parser.add_option('-n','--name',
                    dest="filename",
                    default="hamiltonian.xml",
                    help="name of the output matrix file")
  parser.add_option('-r','--rydberg',
                    dest="rydberg",
                    default=False,
                    action="store_true",
                    help="express the matrix elements in Rydberg")
  parser.add_option('--row_major',
                    dest="row_major",
                    default=False,
                    action="store_true",
                    help="reads the hamiltonian matrix assuming row major format")
  parser.add_option('-t', '--tolerance',
                    dest="tolerance",
                    default=0.001,
                    type="float",
                    help="tolerance on the norm of the anti-hermitian part of the matrix")
  #
  options, remainder = parser.parse_args()
  #
  if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)
  #
  # calling the main routine to extract the matrix
  # Default behavior:
  # - matrix elements are read from file as expressed in Hartree(s)
  # - matrix elements are written to file in eV
  # - matrix is read from file in Fortran format (column major)
  # - matrix is written to file in C format (row major)
  # - matrix is symmetrized to ensure hermiticity
  extract_and_print(options.path,
                    options.filename,
                    options.rydberg,
                    options.row_major,
                    options.tolerance)
  #
  return

def extract_and_print(path,filename,rydberg,row_major,tolerance):
  #
  # let us try to open the 'hamiltonian*.xml' file
  if not os.path.exists(path):
    print ""
    print " Error in extract_hamiltonian_matrix.py :"
    print " '%s'" %(path)
    print " This file does not exist"
    print ""
    sys.exit(0)
  if not os.path.isfile(path):
    print ""
    print " Error in extract_hamiltonian_matrix.py :"
    print " '%s'" %(path)
    print " This is not a regular file"
    print ""
    sys.exit(0)
  try:
    raw_file = open(path,'r')
    lines = raw_file.readlines()
    raw_file.close()
  except:
    print ""
    print " Error in extract_hamiltonian_matrix.py :"
    print " Could not open file :"
    print " '%s'" %(path)
    print ""
    sys.exit(0)
  #
  # now that we have the content of the file in memory
  # we will extract the data.
  # first the starting and finishing lines for the ham
  # matrix
  startline = 0
  for line in lines:
    if not line.startswith('<HAMILTONIAN'):
      startline += 1
    else:
      break
  stopline = 0
  for line in lines:
    if not line.startswith('</HAMILTONIAN'):
      stopline += 1
    else:
      break
  #
  # now the type of data and the size of the matrix
  strings = lines[startline].split()
  for element in strings:
    if 'type=' in element:
      data_type = element[6:-1]  # this will extract only the type from the string 'element'
    if 'size=' in element:
      size = element[6:-2]       # this will extract only the integer corresponding to the size
  #
  if data_type not in ['real']:
    print ""
    print " Error in extract_hamiltonian_matrix.py :"
    print " data type is : '%s'" %(data_type)
    print " Expecting 'real'"
    print ""
    sys.exit(0)
  else:
    type = 'double'
  try:
    matrix_size = int(size)
  except:
    print ""
    print " Error in extract_hamiltonian_matrix.py :"
    print " cound not extract matrix size."
    print " size = '%s'" %(size)
    print ""
    sys.exit(0)
  linear_size = int( round( float(matrix_size)**(0.5) ) )
  if linear_size**2!=matrix_size:
    print ""
    print " Error in extract_hamiltonian_matrix.py :"
    print " matrix size '%s' does not lead to a square matrix" %(size)
    print ""
    sys.exit(0)
  #
  # now reading the matrix in column-major format
  if not row_major:
    hamiltonian = numpy.zeros((linear_size,linear_size),dtype=type)
    iterator = 0
    for i in xrange(startline+1,stopline):
      strings = lines[i].split()
      if len(strings) > 0:
        for element in strings:
          coefficiant = numpy.double(element)
          col = iterator/linear_size
          row = iterator%linear_size
          hamiltonian[row,col] = coefficiant
          iterator += 1
  #
  # or in row-major format if requested
  else:
    hamiltonian = numpy.zeros((linear_size,linear_size),dtype=type)
    iterator = 0
    for i in xrange(startline+1,stopline):
      strings = lines[i].split()
      if len(strings) > 0:
        for element in strings:
          coefficiant = numpy.double(element)
          row = iterator/linear_size
          col = iterator%linear_size
          hamiltonian[row,col] = coefficiant
          iterator += 1
  #
  # decomposing the matrix into hermitian and anti-hermitian parts
  # here the matrix is real so hermitian=symmetric...
  hermitian_part      = 0.5 * ( hamiltonian + numpy.transpose(hamiltonian) )
  anti_hermitian_part = 0.5 * ( hamiltonian - numpy.transpose(hamiltonian) )
  #
  # computing the Frobenius norm of the anti-hermitian part
  norm = numpy.linalg.norm(anti_hermitian_part)
  if norm > tolerance:
    print ""
    print " Error in extract_hamiltonian_matrix.py :"
    print " the Frobenius norm of the anti-hermitian part"
    print " of the matrix is : %10.6f" %(norm)
    print " tolerance is     : %10.6f" %(tolerance)
    answer = "dummy"
    while answer not in ['yes', 'no']:
      answer = raw_input(' Do you wish to continue (yes/no) ?')
    if answer=="no":
      print ""
      sys.exit(0)
  #
  # computing the conversion factor giving the physical unit
  # of the matrix element (original unit is Hartree)
  if rydberg:
    conversion = numpy.double(0.5)
  else:
    conversion = numpy.double(27.211396132)
  hermitian_part = conversion * hermitian_part
  #
  # finally we create the hamiltonian matrix file
  outfile = open(filename,'w')
  outfile.write(" # file written by 'extract_hamiltonian_matrix.py' on %s at %s\n" %(time.strftime("%A, %d %B %Y",time.localtime()),
                                                                                     time.strftime("%H:%M:%S",time.localtime())))
  outfile.write(" %6i\n" %(linear_size))
  for row in xrange(hermitian_part.shape[0]):
    for col in xrange(hermitian_part.shape[1]):
      outfile.write(" %10.6f" %(hermitian_part[row,col]))
    outfile.write("\n")
  outfile.close()
  #
  # printing the norm of the anti-hermitian part for information
  print ""
  print " Frobenius norm of the anti-hermitian part"
  print " of the matrix : %10.6f" %(norm)
  print " tolerance is  : %10.6f" %(tolerance)
  print ""
  #
  return

if __name__=="__main__":
  main()

