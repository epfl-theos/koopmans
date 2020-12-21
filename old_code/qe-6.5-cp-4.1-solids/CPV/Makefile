# Makefile for CP/FPMD

include	../make_cp.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS_CP)

FOBJS = \
atoms_type.o \
bessel.o \
berryion.o \
bforceion.o \
cg.o \
cg_sub.o \
cglib.o \
chargedensity.o \
chargemix.o \
chi2.o \
compute_fes_grads.o \
compute_scf.o \
cp_autopilot.o \
cp_emass.o \
cp_fpmd.o \
cp_interfaces.o \
cp_restart.o \
cplib.o \
cpr_mod.o \
cpr.o \
cprsub.o \
dealloc.o \
dforceb.o \
eelib.o \
efermi.o \
efield.o \
eigs0.o \
electrons.o \
emptystates.o \
ensemble_dft.o \
environment.o \
exch_corr.o \
fft.o \
forces.o \
fromscra.o \
gradrho.o \
gtable.o \
init.o \
init_run.o \
inner_loop.o \
inner_loop_cold.o \
inner_loop_smear.o \
input.o \
io_pot_sic_xml.o \
ions_positions.o \
ksstates.o \
main.o \
mainvar.o \
main_loops.o \
cplib_meta.o \
metaxc.o \
modules.o \
move_electrons.o \
nksiclib.o \
hflib.o \
nl_base.o \
nlcc.o \
ortho_base.o \
ortho.o \
para.o \
path_routines.o \
phasefactor.o \
polarization.o \
potentials.o \
pres_ai_mod.o \
print_out.o \
problem_size.o \
pseudo_base.o \
pseudopot.o \
pseudopot_sub.o \
qmatrixd.o \
qqberry.o \
read_pseudo.o \
restart.o \
restart_sub.o \
runcp.o \
spharmonic.o \
spline.o \
stop_run.o \
stress.o \
turbo.o \
util.o \
vanderwaals.o \
vol_clu.o \
cp_version.o \
wannier_base.o \
wannier.o \
waveinit.o \
wave.o \
wave_types.o \
wf.o \
writetofile.o \
gram_swap.o \
odd_alpha.o \
cg_empty_sub.o \
wave_init_wannier.o \
inner_loop_generalize.o \
makov_payne.o \
pc3nc_fixed.o \
empty_koopmans_pp.o 

# dgradcorr.o 
# setup_gga.o

LOBJS = \
adjef.o \
entropy.o \
forceconv.o \
geninv.o \
indices.o

MODULES = \
../Modules_CP/atom.o \
../Modules_CP/autopilot.o \
../Modules_CP/basic_algebra_routines.o \
../Modules_CP/berry_phase.o \
../Modules_CP/cell_base.o \
../Modules_CP/check_stop.o \
../Modules_CP/clocks.o \
../Modules_CP/constants.o \
../Modules_CP/constraints_module.o \
../Modules_CP/control_flags.o \
../Modules_CP/descriptors.o \
../Modules_CP/dspev_drv.o \
../Modules_CP/electrons_base.o \
../Modules_CP/energies.o \
../Modules_CP/error_handler.o \
../Modules_CP/fft_base.o \
../Modules_CP/fft_parallel.o \
../Modules_CP/fft_scalar.o \
../Modules_CP/fft_types.o \
../Modules_CP/functionals.o \
../Modules_CP/griddim.o \
../Modules_CP/input_parameters.o \
../Modules_CP/io_files.o \
../Modules_CP/io_global.o \
../Modules_CP/ions_base.o \
../Modules_CP/ions_nose.o \
../Modules_CP/kind.o \
../Modules_CP/mp_global.o \
../Modules_CP/mp_wave.o \
../Modules_CP/mp.o \
../Modules_CP/mp_base.o \
../Modules_CP/metagga.o \
../Modules_CP/metadyn_base.o \
../Modules_CP/metadyn_io.o \
../Modules_CP/metadyn_vars.o \
../Modules_CP/parallel_types.o \
../Modules_CP/path_base.o \
../Modules_CP/path_formats.o \
../Modules_CP/path_variables.o \
../Modules_CP/path_opt_routines.o \
../Modules_CP/path_io_routines.o \
../Modules_CP/path_reparametrisation.o \
../Modules_CP/parallel_include.o \
../Modules_CP/parameters.o \
../Modules_CP/parser.o \
../Modules_CP/printout_base.o \
../Modules_CP/pseudo_types.o \
../Modules_CP/ptoolkit.o \
../Modules_CP/radial_grids.o \
../Modules_CP/random_numbers.o \
../Modules_CP/read_upf_v1.o \
../Modules_CP/read_upf_v2.o \
../Modules_CP/read_cards.o \
../Modules_CP/read_namelists.o \
../Modules_CP/read_uspp.o \
../Modules_CP/recvec.o \
../Modules_CP/shmem_include.o \
../Modules_CP/sic.o \
../Modules_CP/smallbox.o \
../Modules_CP/splinelib.o \
../Modules_CP/stick_base.o \
../Modules_CP/task_groups.o \
../Modules_CP/timestep.o \
../Modules_CP/twin_types.o \
../Modules_CP/upf_to_internal.o \
../Modules_CP/upf.o \
../Modules_CP/uspp.o \
../Modules_CP/version.o \
../Modules_CP/vxc_t.o \
../Modules_CP/wavefunctions.o \
../Modules_CP/wave_base.o \
../Modules_CP/write_upf_v2.o \
../Modules_CP/xml_input.o \
../Modules_CP/xml_io_base.o \
../Modules_CP/zhpev_drv.o \
../Modules_CP/wannier_new.o \
../Modules_CP/wrappers.o \
../Modules_CP/compute_dipole.o

AFCLIB=../AFC90_CP/src/libafc90.a

TLDEPS= bindir mods_CP libs_CP libiotk_CP afclib_CP

all : cp
cp : tldeps libcp.a cp.x cppp.x

cp.x : cprstart.o libcp.a $(LIBOBJS_CP) $(AFCLIB)
	$(LD) $(LDFLAGS) -o cp.x \
		cprstart.o $(MODULES) libcp.a \
                $(LIBOBJS_CP) $(AFCLIB) $(LIBS)
	- ( cd ../bin ; ln -fs ../CPV/cp.x . )

cp_test : tldeps libcp.a cp_test.x cppp.x

cp_test.x : cprstart.o libcp.a $(LIBOBJS_CP)
	$(LD) $(LDFLAGS) -o cp_test.x \
		cprstart.o $(MODULES) libcp.a $(LIBOBJS_CP) $(LIBS)
	- ( cd ../bin ; ln -fs ../CPV/cp_test.x . )

cp_bad_sd : tldeps libcp.a cp_bad_sd.x cppp.x

cp_bad_sd.x : cprstart.o libcp.a $(LIBOBJS_CP)
	$(LD) $(LDFLAGS) -o cp_bad_sd.x \
		cprstart.o $(MODULES) libcp.a $(LIBOBJS_CP) $(LIBS)
	- ( cd ../bin ; ln -fs ../CPV/cp_bad_sd.x . )

libcp.a : $(FOBJS) $(LOBJS) 
	 $(AR) $(ARFLAGS) $@ $?
	 $(RANLIB) $@

cp_version.o : cpver.h

cpver.h :
	echo "CHARACTER(LEN=70), PARAMETER :: version_date = '"`date`"'" \
		> cpver.h

PPOBJS = \
../Modules_CP/xml_io_base.o \
../Modules_CP/mp.o \
../Modules_CP/mp_base.o \
../Modules_CP/mp_global.o \
../Modules_CP/io_global.o \
../Modules_CP/io_files.o \
../Modules_CP/mp_wave.o \
../Modules_CP/parser.o \
../Modules_CP/kind.o \
../Modules_CP/control_flags.o \
../Modules_CP/parameters.o \
../Modules_CP/parallel_include.o \
../Modules_CP/error_handler.o \
../Modules_CP/constants.o \
../Modules_CP/wrappers.o

cppp.x : fpmdpp.o $(PPOBJS) $(LIBOBJS_CP)
	$(LD) $(LDFLAGS) -o cppp.x fpmdpp.o $(PPOBJS) $(LIBOBJS_CP) $(LIBS)
	- (cd ../bin ; ln -fs ../CPV/cppp.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f cppp.x *.o *.mod cpver.h *.i core* *.F90 fort* \
		*.cpp *.d *.L *.a *.s cp.x

include make.depend
