# Makefile for Koopmans utils

-include ../../qe_koopmans/make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../qe_koopmans/PW/src \
          $(MOD_FLAG)../../qe_koopmans/Modules/ \
          $(MOD_FLAG)../../qe_koopmans/dft-d3/


PWOBJS = ../../qe_koopmans/PW/src/libpw.a ../../qe_koopmans/KS_Solvers/libks_solvers.a ../../qe_koopmans/dft-d3/libdftd3qe.a ../../qe_koopmans/PP/src/libpp.a
QEMODS = $(BASEMODS)

WAN2KCP_OBJS = wannier.o read_wannier.o fft_supercell.o scell_wfc.o plot_wan2kcp.o \
               cp_files.o wannier2kcp.o

MODULES = $(PWOBJS) $(QEMODS)

all : pwlib pplib epsilon.x wann2kcp.x merge_evc.x

pwlib:
	cd ../../qe_koopmans/ ; $(MAKE) pwlibs || exit 1
pplib:
	cd ../../qe_koopmans/PP/src ; if [ ! -e libpp.a ]; then $(MAKE) libpp.a ; fi || exit 1 

merge_evc.x : merge_evc.o
	$(LD) $(LDFLAGS) -o $@ merge_evc.o
	- ( cd ../bin ; ln -fs ../src/$@ . )

wann2kcp.x : wann2kcp.o $(WAN2KCP_OBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		wann2kcp.o $(WAN2KCP_OBJS) $(MODULES) $(LIBOBJS) $(QELIBS)
	- ( cd ../bin ; ln -fs ../src/$@ . )

epsilon.x : epsilon.o $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ epsilon.o $(MODULES) \
	$(LIBOBJS) $(QELIBS)
	- ( cd ../bin ; ln -fs ../src/$@ . )

clean :
	- /bin/rm -f *.x *.o *~ *_tmp.f90 *.d *.mod *.i *.L *genmod*

include make.depend
