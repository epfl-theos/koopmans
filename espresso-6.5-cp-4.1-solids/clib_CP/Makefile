# Makefile for clib

include ../make_cp.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS_CP)

OBJS = \
stack.o \
c_mkdir.o \
cptimer.o \
eval_infix.o \
fft_stick.o \
indici.o \
memstat.o \
qsort.o

all : clib.a

clib.a : $(OBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

source :
	co -l $(OBJS:.o=.c)

clean :
	- rm -f clib.a *.o *.mod *.i core*

include make.depend
