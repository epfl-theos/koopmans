# Makefile for PostProc
sinclude make.inc

default: all

all:
	if test -d src ; then \
	( cd src ; $(MAKE) || exit 1 ) ; fi

clean :
	if test -d src ; then \
	( cd src ; $(MAKE) clean ) ; fi
