#=======================================
# WANNIER90 	
#=======================================

TOPDIR = /home/yshubert/code/koopmans/quantum_espresso/cp_koopmans

F90=mpiifort
FCOPTS=$(FFLAGS) -nomodule 
LDOPTS= 

LIBS =    -lmkl_intel_lp64  -lmkl_sequential -lmkl_core 
