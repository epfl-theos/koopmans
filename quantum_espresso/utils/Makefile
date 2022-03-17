# Makefile for PostProc

default : all

all : bindir
	if test -d src ; then \
	( cd src ; $(MAKE) || exit 1 ) ; fi

bindir :
	test -d bin || mkdir bin

links : bindir
	( cd bin/ ; \
	rm -f *.x ; \
	for exe in ../src/*.x ; do \
	    if test ! -L $$exe ; then ln -fs $$exe . ; fi \
	done ; )

clean :
	if test -d src ; then \
	( cd src ; $(MAKE) clean ) ; fi
	- /bin/rm -f bin/*.x