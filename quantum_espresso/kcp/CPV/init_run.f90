!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE init_run()
   !----------------------------------------------------------------------------
   !
   ! ... this routine initialise the cp code and allocates (calling the
   ! ... appropriate routines) the memory
   !
   USE kinds, ONLY: DP
   USE control_flags, ONLY: nbeg, nomore, lwf, iprsta, iprint, &
                            tfor, tprnfor, tpre, &
                            newnfi, tnewnfi, ndw, non_ortho, &
                            iprint_manifold_overlap
   !above, added:giovanni non_ortho, iprint_manifold_overlap
   USE cp_electronic_mass, ONLY: emass_cutoff
   USE ions_base, ONLY: na, nax, nat, nsp, pmass, cdms
   USE ions_positions, ONLY: tau0, taum, taup, taus, tausm, tausp, &
                             vels, velsm, velsp, fion
   USE gvecw, ONLY: ngw, ngwt, ggp
   USE gvecb, ONLY: ngb
   USE gvecs, ONLY: ngs
   USE gvecp, ONLY: ngm
   USE reciprocal_vectors, ONLY: gzero
   USE grid_dimensions, ONLY: nnrx, nr1, nr2, nr3
   USE fft_base, ONLY: dfftp
   USE electrons_base, ONLY: nspin, nbsp, nbspx, nupdwn, f
   USE uspp, ONLY: nkb, vkb, deeq, becsum, nkbus
   USE smooth_grid_dimensions, ONLY: nnrsx
   USE wavefunctions_module, ONLY: c0, cm, cp, cdual, cmdual, cstart
   USE cdvan, ONLY: dbec, drhovan
   USE ensemble_dft, ONLY: tens, z0t, tsmear
   USE cg_module, ONLY: tcg
   USE electrons_base, ONLY: nudx
   USE efield_module, ONLY: tefield, tefield2
   USE uspp_param, ONLY: nhm
   USE ions_nose, ONLY: xnhp0, xnhpm, vnhp, nhpcl, nhpdim
   USE cell_base, ONLY: h, hold, hnew, velh, tpiba2, ibrav, &
                        alat, a1, a2, a3, b1, b2, b3
   USE cp_main_variables, ONLY: lambda, lambdam, ema0bg, &
                                acc, acc_this_run, wfill, hamilt, &
                                edft, nfi, ht0, iprint_stdout, &
                                gamma_only, do_wf_cmplx !added:giovanni gamma_only, do_wf_cmplx
   USE cp_main_variables, ONLY: allocate_mainvar, nlax, descla, nrlx, nlam
   USE energies, ONLY: ekincm
   USE time_step, ONLY: tps
   USE electrons_nose, ONLY: xnhe0, xnhem, vnhe
   USE cell_nose, ONLY: xnhh0, xnhhm, vnhh
   USE funct, ONLY: dft_is_meta
   USE metagga, ONLY: crosstaus, dkedtaus, gradwfc
   !
   USE efcalc, ONLY: clear_nbeg
   USE local_pseudo, ONLY: allocate_local_pseudo
   USE cp_electronic_mass, ONLY: emass_precond
   USE wannier_subroutines, ONLY: wannier_startup
   USE cp_interfaces, ONLY: readfile, symm_wannier
   USE ions_base, ONLY: ions_cofmass
   USE ensemble_dft, ONLY: id_matrix_init, allocate_ensemble_dft, h_matrix_init
   USE efield_module, ONLY: allocate_efield, allocate_efield2
   USE cg_module, ONLY: allocate_cg
   USE wannier_module, ONLY: allocate_wannier
   USE io_files, ONLY: outdir, prefix
   USE io_global, ONLY: stdout
   USE printout_base, ONLY: printout_base_init
   USE wave_types, ONLY: wave_descriptor_info
   USE xml_io_base, ONLY: restart_dir, create_directory
   USE orthogonalize_base, ONLY: mesure_diag_perf, mesure_mmul_perf
   USE step_constraint, ONLY: step_con
   USE ions_base, ONLY: ions_reference_positions
   USE ldau
   use eecp_mod, ONLY: do_comp, which_compensation, tcc_odd
   use efield_mod, ONLY: do_efield
   USE nksic, ONLY: do_orbdep
   use input_parameters, ONLY: odd_nkscalfact, restart_odd_nkscalfact, restart_mode, &
                               restart_from_wannier_cp, wannier_empty_only, &
                               restart_from_wannier_pwscf, impose_bloch_symm
   use wavefunctions_module, ONLY: c0_fixed
   USE twin_types !added:giovanni
   !
   IMPLICIT NONE
   !
   INTEGER            :: i, iss, ndim
   CHARACTER(LEN=256) :: dirname
   LOGICAL :: lgam
   !
   !
   CALL start_clock('initialize')
   !
   lgam = gamma_only .and. .not. do_wf_cmplx
   !
   ! ... initialize directories
   !
   write (6, *) "nbeg", nbeg
   IF (nbeg < 0) THEN
      CALL create_directory(outdir)
   END IF
   !
   CALL printout_base_init(outdir, prefix)
   !
   dirname = restart_dir(outdir, ndw)
   !
   ! ... Create main restart directory
   !
   CALL create_directory(dirname)
   !
   ! ... initialize g-vectors, fft grids
   ! ... The number of g-vectors are based on the input celldm!
   !
   CALL init_dimensions()
   !
   ! ... initialize atomic positions and cell
   !
   CALL init_geometry()
   !
   ! ... mesure performances of parallel routines
   !
   CALL mesure_mmul_perf(nudx)
   !
   CALL mesure_diag_perf(nudx)
   !
   IF (lwf) CALL clear_nbeg(nbeg)
   !
   !=======================================================================
   !     allocate and initialize nonlocal potentials
   !=======================================================================
   !
   CALL nlinit()
   !
   !=======================================================================
   !     allocation of all arrays not already allocated in init and nlinit
   !=======================================================================
   !
   CALL allocate_mainvar(ngw, ngwt, ngb, ngs, ngm, nr1, nr2, nr3, dfftp%nr1x, &
                         dfftp%nr2x, dfftp%npl, nnrx, nnrsx, nat, nax, nsp, &
                         nspin, nbsp, nbspx, nupdwn, nkb, gzero, nudx, &
                         tpre)
   !
   CALL allocate_local_pseudo(ngs, nsp)
   !
   !  initialize wave functions descriptors and allocate wf
   !
   ALLOCATE (c0(ngw, nbspx))
   !gvn23 ALLOCATE(c0i(ngw, nbspx))
   ALLOCATE (cm(ngw, nbspx))
   !gvn23 ALLOCATE(cmi(ngw, nbspx))
   ALLOCATE (cp(ngw, nbspx))
   !gvn23 ALLOCATE(cpi(ngw, nbspx))
   IF (iprint_manifold_overlap > 0) THEN
      ALLOCATE (cstart(ngw, nbspx))
   END IF
   !
   IF (odd_nkscalfact) ALLOCATE (c0_fixed(ngw, nbspx))
   !
   IF (non_ortho) THEN
      ALLOCATE (cdual(ngw, nbspx))
      ALLOCATE (cmdual(ngw, nbspx))
   END IF
   !
   IF (iprsta > 2) THEN
      !
      CALL wave_descriptor_info(wfill, 'wfill', stdout)
      !
   END IF
   !
   ! Depending on the verbosity set the frequency of
   ! verbose information to stdout
   !
   IF (iprsta < 1) iprint_stdout = 100*iprint
   IF (iprsta == 1) iprint_stdout = 10*iprint
   IF (iprsta > 1) iprint_stdout = iprint
   !
   acc = 0.D0
   acc_this_run = 0.D0
   !
   edft%ent = 0.D0
   edft%esr = 0.D0
   edft%evdw = 0.D0
   edft%ekin = 0.D0
   edft%enl = 0.D0
   edft%etot = 0.D0
   !
   ALLOCATE (becsum(nhm*(nhm + 1)/2, nat, nspin))
   ALLOCATE (deeq(nhm, nhm, nat, nspin))
   IF (tpre) THEN
      ALLOCATE (dbec(nkb, 2*nlam, 3, 3))
      ALLOCATE (drhovan(nhm*(nhm + 1)/2, nat, nspin, 3, 3))
   END IF
   !
   ALLOCATE (vkb(ngw, nkb))
   !
   IF (dft_is_meta() .AND. tens) &
      CALL errore('cprmain ', 'ensemble_dft not implimented for metaGGA', 1)
   !
   IF (dft_is_meta() .AND. tpre) THEN
      !
      ALLOCATE (crosstaus(nnrsx, 6, nspin))
      ALLOCATE (dkedtaus(nnrsx, 3, 3, nspin))
      ALLOCATE (gradwfc(nnrsx, 3))
      !
   END IF
   !
   IF (lwf) CALL allocate_wannier(nbsp, nnrsx, nspin, ngm)
   !
   IF (tens .OR. tcg .OR. tsmear) &
      CALL allocate_ensemble_dft(nkb, nbsp, ngw, nudx, nspin, nbspx, nnrsx, nat, nlax, nrlx, lgam)
   !
   IF (tcg) CALL allocate_cg(ngw, nbspx, nkbus)
   !
   IF (tefield) CALL allocate_efield(ngw, ngwt, nbspx, nhm, nax, nsp)
   IF (tefield2) CALL allocate_efield2(ngw, nbspx, nhm, nax, nsp)
   !
   IF (ALLOCATED(deeq)) deeq(:, :, :, :) = 0.D0
   !
   IF (ALLOCATED(lambda)) THEN
      DO iss = 1, size(lambda)
         IF (lambda(iss)%isalloc) THEN
            IF (.not. lambda(iss)%iscmplx) THEN
               lambda(iss)%rvec = 0.D0
            ELSE
               lambda(iss)%cvec = CMPLX(0.D0, 0.D0)
            END IF
         END IF
      END DO
   END IF

   IF (ALLOCATED(lambdam)) THEN
      DO iss = 1, size(lambdam)
         IF (lambdam(iss)%isalloc) THEN
            IF (.not. lambdam(iss)%iscmplx) THEN
               lambdam(iss)%rvec = 0.D0
            ELSE
               lambdam(iss)%cvec = CMPLX(0.D0, 0.D0)
            END IF
         END IF
      END DO
   END IF
   !
   taum = tau0
   taup = 0.D0
   tausm = taus
   tausp = 0.D0
   vels = 0.D0
   velsm = 0.D0
   velsp = 0.D0
   !
   hnew = h
   !
   cm = (0.D0, 0.D0)
   !gvn23 cmi = (0.D0, 0.D0)
   c0 = (0.D0, 0.D0)
   !gvn23 c0i = (0.D0, 0.D0)
   cp = (0.D0, 0.D0)
   !gvn23 cpi = (0.D0, 0.D0)
   !
   IF (odd_nkscalfact) c0_fixed = (0.D0, 0.D0)
   !
   IF (tens .OR. tsmear) then
      !
      CALL id_matrix_init(descla, nspin)
      CALL h_matrix_init(descla, nspin)
      !
   END IF
   !
   IF (lwf) CALL wannier_startup(ibrav, alat, a1, a2, a3, b1, b2, b3)
   !
   ! ... Calculate: ema0bg = ecutmass /  MAX( 1.0d0, (2pi/alat)^2 * |G|^2 )
   !
   CALL emass_precond(ema0bg, ggp, ngw, tpiba2, emass_cutoff)
   !
   CALL print_legend()
   !
   step_con = .FALSE.
   !
   CALL ldau_init()
   !
   CALL nksic_init()
   !
   CALL hf_init()
   !
   CALL ee_init()
   !
   CALL efield_init()
   !
   IF (do_comp) THEN
      !
      write (stdout, *) "USING TCC FOR ODD", tcc_odd
      !
      IF (trim(which_compensation) == 'tcc1d') THEN
         CALL ee_green_1d_init(ht0)
         IF (tcc_odd) THEN
            CALL ee_green_0d_init(ht0)
         END IF
      ELSE IF (trim(which_compensation) == 'tcc2d') THEN
         CALL ee_green_2d_init(ht0)
         IF (tcc_odd) THEN
            CALL ee_green_0d_init(ht0)
         END IF
      ELSE
         CALL ee_green_0d_init(ht0)
      END IF
      !
   END IF
   !
   IF (do_orbdep .AND. iprsta > 1) THEN
      !
      ndim = MAXVAL(nupdwn(:))

      ALLOCATE (hamilt(nspin))
      DO iss = 1, nspin
         call init_twin(hamilt(iss), lgam)
         call allocate_twin(hamilt(iss), ndim, ndim, lgam)
      END DO
!       ALLOCATE( hamilt( ndim, ndim, nspin) )
      !
   ELSE
      ALLOCATE (hamilt(1))
      call init_twin(hamilt(1), lgam)
      call allocate_twin(hamilt(1), 1, 1, lgam)
   END IF
   !
   DO iss = 1, size(hamilt)
      if (.not. hamilt(iss)%iscmplx) then
         hamilt(iss)%rvec = 0.0d0
      else
         hamilt(iss)%cvec = CMPLX(0.0d0, 0.d0)
      end if
   END DO

   IF (do_efield) CALL ee_efieldpot_init(ht0)

   IF (nbeg < 0) THEN
      !
      !======================================================================
      !     Initialize from scratch nbeg = -1
      !======================================================================
      !
      nfi = 0
      !
      CALL from_scratch()
      !
   ELSE
      !
      !======================================================================
      !     nbeg = 0, nbeg = 1
      !======================================================================
      !
      i = 1
      !gvn22 !change readfile?
      CALL readfile(i, h, hold, nfi, c0, cm, taus, &
                    tausm, vels, velsm, acc, lambda, lambdam, xnhe0, xnhem, &
                    vnhe, xnhp0, xnhpm, vnhp, nhpcl, nhpdim, ekincm, xnhh0, xnhhm, &
                    vnhh, velh, fion, tps, z0t, f)
      !
      CALL from_restart()
      !
   END IF
   !
   !=======================================================================
   !     restart with new averages and nfi=0
   !=======================================================================
   !
   ! ... reset some variables if nbeg < 0
   ! ... ( new simulation or step counter reset to 0 )
   !
   IF (nbeg <= 0) THEN
      !
      acc = 0.D0
      nfi = 0
      !
   END IF
   !
   IF (.NOT. tfor .AND. .NOT. tprnfor) fion(:, :) = 0.D0
   !
   IF (tnewnfi) nfi = newnfi
   !
   nomore = nomore + nfi
   !
   !  Set center of mass for scaled coordinates
   !
   CALL ions_cofmass(taus, pmass, na, nsp, cdms)
   !
   IF (nbeg <= 0 .OR. lwf) THEN
      !
      CALL ions_reference_positions(tau0)
      !
   END IF
   !
   ! here we provide an option to restart wfc from wannier orbitals
   ! for occupied many-folds
   !
   IF (restart_from_wannier_cp .or. restart_from_wannier_pwscf .and. .not. wannier_empty_only) THEN
      !
      write (stdout, *) "in init_run from wannier start Linh"
      !
      IF (TRIM(restart_mode) == "from_scratch") THEN
         !
         CALL errore('init_run ', 'A restart from wannier orbitals needs restart_mode = restart', 1)
         !
      END IF
      !
      IF (restart_from_wannier_cp .and. restart_from_wannier_pwscf) THEN
         !
         CALL errore('init_run ', 'choose either restart_from_wannier_pwscf or restart_from_wannier_cp == true', 1)
         !
      END IF
      !
      IF (restart_from_wannier_pwscf) CALL wave_init_wannier_pwscf(c0, nbspx)
      !
      IF (restart_from_wannier_cp) CALL wave_init_wannier_cp(c0, ngw, nbspx, .True.)
      !
      write (stdout, *) "in init_run from wannier end Linh"
      !
   END IF
   !
   IF (odd_nkscalfact) THEN
      !
      IF (.not. restart_odd_nkscalfact) then
         c0_fixed(:, :) = c0(:, :)
      END IF
      !
   END IF
   !
   CALL stop_clock('initialize')
   !
   IF (impose_bloch_symm) CALL symm_wannier(c0, nbspx, .false.)
   !
   RETURN
   !
END SUBROUTINE init_run
