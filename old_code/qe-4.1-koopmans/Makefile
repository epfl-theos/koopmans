default :
	@echo 'to install, type at the shell prompt:'
	@echo '  ./configure'
	@echo '  make target'
	@echo 'where target is one of the following:'
	@echo '  pw           basic code for scf, structure optimization, MD'
	@echo '  cp           CP code: CP MD with ultrasoft pseudopotentials'
	@echo '  pp           postprocessing programs'
	@echo '  gamma        Gamma-only version of phonon code'
	@echo '  w90          Maximally localised Wannier Functions'
	@echo '  ld1          utilities for pseudopotential generation'
	@echo '  pwall        same as "make pw pp gamma"'
	@echo '  all          same as "make pwall cp ld1"'
	@echo '  clean        remove executables and objects'
	@echo '  veryclean    revert distribution to the original status'
	@echo '  tar          create a tarball of the source tree'
	@echo '  tar-gui      create a tarball of the GUI sources'
	@echo '  doc          build a documentation'
	@echo '  log          create ChangeLog and ChangeLog.html files'

pw : bindir mods libs libiotk eelib
	if test -d PW ; then \
	( cd PW ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

cp : bindir mods libs libiotk afclib
	if test -d CPV ; then \
	( cd CPV ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= cp ; \
	else $(MAKE) $(MFLAGS) TLDEPS= cp ; fi ) ; fi

afc90 :
	if test -d AFC90 ; then \
	( cd AFC90 ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

pp : bindir mods libs pw
	if test -d PP ; then \
	( cd PP ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

gamma : bindir mods libs pw
	if test -d Gamma ; then \
	( cd Gamma ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

w90   : bindir mods libs
	if test -d W90 ; then \
	( cd W90 ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= wannier ; \
	else $(MAKE) $(MFLAGS) TLDEPS= wannier ; fi ) ; fi

ld1 : bindir mods libs pw
	if test -d atomic ; then \
	( cd atomic ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

libiotk :
	if test -d iotk ; then \
	( cd iotk ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= lib+util ; \
	else $(MAKE) $(MFLAGS) TLDEPS= lib+util ; fi ) ; fi

pw_export : libiotk bindir mods libs pw
	if test -d PP ; then \
	( cd PP ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= pw_export.x ; \
	else $(MAKE) $(MFLAGS) TLDEPS= pw_export.x ; fi ) ; fi

pwall : pw pp gamma 
all   : pwall cp ld1 w90

mods : libiotk
	( cd Modules ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )

afclib:
	( if test -d AFC90/src; then cd AFC90/src ; \
        if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= libafc90.a ;\
        else $(MAKE) $(MFLAGS) TLDEPS= libafc90.a ; fi; fi )

libs : mods
	( cd clib ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )
	( cd flib ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )

eelib :
	( cd EE ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all  ; \
        else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )

bindir :
	test -d bin || mkdir bin

# remove object files and executables
clean :
	touch make.sys 
	for dir in \
		CPV Gamma Modules PP PW EE atomic clib flib \
		iotk W90 Doc \
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= clean ; \
		else $(MAKE) $(MFLAGS) TLDEPS= clean ; fi ) \
	    fi \
	done
	- /bin/rm -rf bin/*.x tmp
	#- cd tests; /bin/rm -rf CRASH *.out *.out2 

# remove configuration files too
distclean veryclean : clean
	- /bin/rm -rf make.sys \
		      config.log configure.msg config.status autom4te.cache \
		      espresso.tar.gz CPV/version.h ChangeLog* \
		      intel.pcl */intel.pcl
	- cd examples ; ./make_clean
	- cd atomic_doc ; ./make_clean
	- if test -d GUI ; then \
	( cd GUI ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= veryclean ; \
		else $(MAKE) $(MFLAGS) TLDEPS= veryclean ; fi ) \
	  fi

tar :
	@if test -f espresso.tar.gz ; then /bin/rm espresso.tar.gz ; fi
	# do not include unneeded stuff 
	find ./ -type f | grep -v -e /CVS/ -e /results/ -e'/\.' -e'\.o$$' \
             -e'\.mod$$' -e'\.a$$' -e'\.d$$' -e'\.i$$' -e'\.F90$$' -e'\.x$$' \
	     -e'~$$' -e'\./GUI' | xargs tar rvf espresso.tar
	gzip espresso.tar

# TAR-GUI works only if we have CVS-sources !!!
tar-gui :
	@if test -d GUI/PWgui ; then \
	    cd GUI/PWgui ; \
	    if test "$(MAKE)" = "" ; then \
		make $(MFLAGS) TLDEPS= clean cvsinit pwgui-source; \
	    else $(MAKE) $(MFLAGS) TLDEPS= clean cvsinit pwgui-source; fi; \
	    mv PWgui-*.tgz ../.. ; \
	else \
	    echo ; \
	    echo "  Sorry, tar-gui works only for CVS-sources !!!" ; \
	    echo ; \
	fi

doc :
	if test -d Doc ; then \
	( cd Doc ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi


log :
	-perl ./cvs2cl.pl
	-perl ./cvs2cl.pl --xml --header /dev/null --stdout \
		| perl ./cl2html.pl --entries 0 > ChangeLog.html

links : bindir
	( cd bin/ ; \
	for exe in \
	    ../CPV/cp.x \
	    ../CPV/cppp.x \
	    ../Gamma/phcg.x \
	    ../PP/average.x ../PP/bands.x ../PP/dos.x \
	      ../PP/plotband.x ../PP/plotrho.x ../PP/pmw.x \
	      ../PP/pp.x ../PP/projwfc.x ../PP/pw2casino.x ../PP/pw2wan.x \
	      ../PP/voronoy.x ../PP/pw_export.x \
	    ../PW/pw.x \
	    ../W90/wannier90.x \
	    ../W90/postw90.x \
	    ../W90/w90chk2chk.x \
	    ../atomic/ld1.x \
	; do \
	      if test -f $$exe ; then ln -fs $$exe . ; fi \
	done \
	)

depend:
	@echo 'Checking dependencies...'
	- ( if test -x ./makedeps.sh ; then ./makedeps.sh ; fi)
# DO NOT DELETE
