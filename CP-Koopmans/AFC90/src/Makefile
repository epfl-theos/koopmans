include ../../make.sys

#include make.rules

#LIBS=-llapack -lblas -lfftw3

# targets

OBJS=\
afcc.o \
besseli.o \
besselk.o \
cylharm.o \
cylharm0.o \
cylharmasympt.o \
cylharmlogslope.o \
cylharmrk.o \
cylharmseries.o \
cylharmslope.o \
ei.o \
eimlog.o \
fft1d.o \
fft2d.o \
fft3d.o \
gaussian.o \
gaussiank.o \
gaussianl.o \
interp.o \
main.o \
nfft.o \
phi0d.o \
phi1d.o \
phi2d.o \
phi3d.o \
reciprocal.o \
rtoaxis.o \
rungekutta.o \
steprk.o \
vectorproduct.o \
volume.o \
ydot.o

AFCFILES= $(OBJS)

all : $(OBJS) libafc90.a

libafc90.a: $(AFCFILES)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f cppp.x *.o *.mod

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

include make.depend
