!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE iosys()
  !-----------------------------------------------------------------------------
  !
  ! ...  this subroutine reads input data from standard input ( unit 5 )
  ! ...  Use "-input filename" to read input from file "filename":
  ! ...  may be useful if you have trouble reading from standard input
  ! ...  ---------------------------------------------------------------
  !
  ! ...  access the modules renaming the variables that have the same name
  ! ...  as the input parameters, this is required in order to use a code
  ! ...  independent input parser
  !
  !
  USE kinds,         ONLY : DP
  USE funct,         ONLY : enforce_input_dft, dft_has_finite_size_correction, set_finite_size_volume
  USE constants,     ONLY : autoev, eV_to_kelvin, pi, rytoev, &
                            uakbar, amconv, bohr_radius_angs
  USE mp_global,     ONLY : npool, nproc_pool
  !
  USE io_global,     ONLY : stdout, ionode
  !
  USE bp,            ONLY : nppstr_    => nppstr, &
                            gdir_      => gdir, &
                            lberry_    => lberry, &
                            lelfield_  => lelfield, &
                            efield_    => efield, &
                            nberrycyc_ => nberrycyc, &
                            efield_cart_ => efield_cart
  !
  USE cell_base,     ONLY : at, bg, alat, omega, &
                            celldm_ => celldm, &
                            ibrav_  => ibrav, &
                            init_dofree
  !
  USE ions_base,     ONLY : if_pos, ityp, tau, extfor, &
                            ntyp_ => nsp, &
                            nat_  => nat, &
                            amass
  !
  USE basis,         ONLY : atomic_positions, startingconfig, &
                            starting_wfc, starting_pot
  !
  USE printout_base, ONLY : title_ => title
  !
  USE cellmd,        ONLY : cmass, omega_old, at_old, ntcheck, &
                            cell_factor_ => cell_factor , &
                            press_       => press, &
                            calc, lmovecell
  !
  USE dynamics_module, ONLY : control_temp, temperature, thermostat, &
                              dt_         => dt, &
                              delta_t_    => delta_t, &
                              nraise_     => nraise, &
                              refold_pos_ => refold_pos
  !
  USE extfield,      ONLY : tefield_  => tefield, &
                            dipfield_ => dipfield, &
                            edir_     => edir, &
                            emaxpos_  => emaxpos, &
                            eopreg_   => eopreg, &
                            eamp_     => eamp, &
                            forcefield
  !
  USE io_files,      ONLY : input_drho, output_drho, trimcheck, &
                            psfile, tmp_dir, wfc_dir, &
                            prefix_     => prefix, &
                            pseudo_dir_ => pseudo_dir
  !
  USE force_mod,     ONLY : lforce, lstres, force
  !
  USE gvect,         ONLY : dual, &
                            nr1_     => nr1, &
                            nr2_     => nr2, &
                            nr3_     => nr3,  &
                            ecutwfc_ => ecutwfc, &
                            ecfixed_ => ecfixed, &
                            qcutz_   => qcutz, &
                            q2sigma_ => q2sigma
  !
  USE gsmooth,       ONLY : nr1s_ => nr1s, &
                            nr2s_ => nr2s, &
                            nr3s_ => nr3s
  !
  USE klist,         ONLY : ngauss, two_fermi_energies, &
                            degauss_           => degauss, &
                            nelec_             => nelec, &
                            nelup_             => nelup, &
                            neldw_             => neldw, &
                            tot_charge_        => tot_charge, &
                            tot_magnetization_ => tot_magnetization, &
                            multiplicity_      => multiplicity
  !
  USE ktetra,        ONLY : ltetra, nk1, nk2, nk3
  !
  USE ldaU,          ONLY : Hubbard_U_     => hubbard_u, &
                            Hubbard_alpha_ => hubbard_alpha, &
                            lda_plus_u_    => lda_plus_u, &
                            niter_with_fixed_ns, starting_ns, U_projection
  !
! DCC
  USE ee_mod,        ONLY :           &
                            do_comp_ => do_comp, &
                            do_coarse_ => do_coarse, &
                            do_mltgrid_ => do_mltgrid, &
                            n_self_interaction_ => n_self_interaction, &
                            mixing_charge_compensation_ => mixing_charge_compensation, &
                            mr1_ => mr1, &
                            mr2_ => mr2, &
                            mr3_ => mr3, &
                            ecutcoarse_ => ecutcoarse, &
                            whichbc_ => whichbc, &
                            errtol_ => errtol, &
                            itmax_ => itmax, &
                            which_compensation_ => which_compensation, &
                            n_charge_compensation_ => n_charge_compensation, & 
                            ncompx_ => ncompx, &
                            ncompy_ => ncompy, &
                            ncompz_ => ncompz, &
                            rhocut_ => rhocut, &
                            rhocution_ => rhocution, &
                            cellmin_ => cellmin, &
                            cellmax_ => cellmax, &
                            n_smoothing_ => n_smoothing, &
                            comp_thr_ => comp_thr, &
                            smoothspr_ => smoothspr, &
                            nlev_ =>  nlev, &
                            deltapot_ => deltapot, &
                            rhoionmax_ => rhoionmax, &
                            which_smoothing_ => which_smoothing
  !
  USE a2F,           ONLY : la2F_ => la2F
  !
  USE exx,           ONLY : x_gamma_extrapolation_ => x_gamma_extrapolation, &
                            nqx1_ => nq1, &
                            nqx2_ => nq2, &
                            nqx3_ => nq3
  !
  !
  USE lsda_mod,      ONLY : nspin_                  => nspin, &
                            starting_magnetization_ => starting_magnetization, &
                            lsda
  !
  USE relax,         ONLY : epse, epsf, epsp, starting_scf_threshold
  !
  USE control_flags, ONLY : isolve, max_cg_iter, david, tr2, imix, &
                            nmix, iverbosity, niter, pot_order, wfc_order, &
                            remove_rigid_rot_ => remove_rigid_rot, &
                            diago_full_acc_   => diago_full_acc, &
                            tolp_             => tolp, &
                            upscale_          => upscale, &
                            mixing_beta_      => mixing_beta, &
                            nstep_            => nstep, &
                            iprint_           => iprint, &
                            nosym_            => nosym, &
                            nosym_evc_        => nosym_evc, &
                            noinv_            => noinv, &
                            modenum_          => modenum, &
                            lkpoint_dir_      => lkpoint_dir, &
                            tqr_              => tqr, &
                            io_level, ethr, lscf, lbfgs, lmd, lpath, lneb,   &
                            lsmd, lphonon, ldamped, lbands, lmetadyn, llang, &
                            lconstrain, lcoarsegrained, restart, twfcollect, &
                            use_para_diag, llondon, nofrac
  USE control_flags, ONLY : ortho_para_ => ortho_para
  !
  USE wvfct,         ONLY : nbnd_ => nbnd
  !
  USE fixed_occ,     ONLY : tfixed_occ, f_inp
  !
  USE path_variables, ONLY : nstep_path, lsteep_des, lquick_min, lbroyden, &
                             llangevin, &
                             ds_              => ds, &
                             use_masses_      => use_masses, &
                             CI_scheme_       => CI_scheme, &
                             fixed_tan_       => fixed_tan, &
                             use_freezing_    => use_freezing, &
                             k_max_           => k_max, &
                             k_min_           => k_min, &
                             num_of_images_   => num_of_images, &
                             first_last_opt_  => first_last_opt, &
                             temp_req_        => temp_req, &
                             path_thr_        => path_thr
  !
  USE noncollin_module, ONLY : i_cons, mcons, &
                               noncolin_  => noncolin, &
                               lambda_    => lambda, &
                               angle1_    => angle1, &
                               angle2_    => angle2, &
                               report_    => report
  !
  USE spin_orb, ONLY : lspinorb_ => lspinorb
  !
  USE bfgs_module,   ONLY : bfgs_ndim_        => bfgs_ndim, &
                            trust_radius_max_ => trust_radius_max, &
                            trust_radius_min_ => trust_radius_min, &
                            trust_radius_ini_ => trust_radius_ini, &
                            w_1_              => w_1, &
                            w_2_              => w_2
  USE wannier_new, ONLY :   use_wannier_      => use_wannier, &
                            use_energy_int_   => use_energy_int, &
                            nwan_             => nwan, &
                            print_wannier_coeff_    => print_wannier_coeff

  USE realus,                ONLY : real_space_ => real_space

  !
  ! ... CONTROL namelist
  !
  USE input_parameters, ONLY : title, calculation, verbosity, restart_mode, &
                               nstep, iprint, tstress, tprnfor, dt, outdir, &
                               wfcdir, prefix, etot_conv_thr, forc_conv_thr, &
                               pseudo_dir, disk_io, tefield, dipfield, lberry, &
                               gdir, nppstr, wf_collect,lelfield, efield, &
                               nberrycyc, lkpoint_dir, efield_cart
  !
  ! ... SYSTEM namelist
  !
  USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                               nat, ntyp, nbnd, nelec, nelup, neldw,        &
                               tot_charge, tot_magnetization, multiplicity, &
                               ecutwfc, ecutrho, nr1, nr2, nr3, nr1s, nr2s, &
                               nr3s, noinv, nosym, nosym_evc,               &
                               force_symmorphic, starting_magnetization,    &
                               occupations, degauss, smearing, nspin,       &
                               ecfixed, qcutz, q2sigma, lda_plus_U,         &
                               Hubbard_U, Hubbard_alpha, input_dft, la2F,   &
                               starting_ns_eigenvalue, U_projection_type,   &
#if defined (EXX)
                               x_gamma_extrapolation, nqx1, nqx2, nqx3, &
#endif
                               edir, emaxpos, eopreg, eamp, noncolin, lambda, &
                               angle1, angle2, constrained_magnetization,     &
                               B_field, fixed_magnetization, report, lspinorb,&
                               assume_isolated, spline_ps, london, london_s6, &
                               london_rcut,                                   &
! DCC
                               do_ee
  !
  ! ... EE namelist
  !
  USE input_parameters, ONLY : which_compensation,                          &
                               n_charge_compensation,                       &
                               mixing_charge_compensation,                  &
                               mr1, mr2, mr3, ecutcoarse,comp_thr,          &
                               errtol,itmax,whichbc,                        &
                               ncompx, ncompy, ncompz,                      &
                               cellmin, cellmax,                            &
                               nlev 
  !
  ! ... ELECTRONS namelist
  !
  USE input_parameters, ONLY : electron_maxstep, mixing_mode, mixing_beta, &
                               mixing_ndim, mixing_fixed_ns, conv_thr,     &
                               tqr, diago_thr_init, diago_cg_maxiter,      &
                               diago_david_ndim, diagonalization,          &
                               diago_full_acc, startingwfc, startingpot,   &
                               real_space
  USE input_parameters, ONLY : ortho_para
  !
  ! ... IONS namelist
  !
  USE input_parameters, ONLY : phase_space, ion_dynamics, ion_positions, tolp, &
                               tempw, delta_t, nraise, ion_temperature,        &
                               refold_pos, remove_rigid_rot, upscale,          &
                               pot_extrapolation,  wfc_extrapolation,          &
                               num_of_images, path_thr, CI_scheme, opt_scheme, &
                               use_masses, first_last_opt, temp_req, k_max,    &
                               k_min, ds, use_freezing, fixed_tan,             &
                               w_1, w_2, trust_radius_max, trust_radius_min,   &
                               trust_radius_ini, bfgs_ndim
  !
  ! ... CELL namelist
  !
  USE input_parameters, ONLY : cell_parameters, cell_dynamics, press, wmass,  &
                               cell_temperature, cell_factor, press_conv_thr, &
                               cell_dofree
  !
  ! ... WANNIER_NEW namelist
  !
  USE input_parameters, ONLY : use_wannier, nwan, &
                               constrain_pot, use_energy_int, print_wannier_coeff

  !
  ! ... "path" specific
  !
  USE input_parameters, ONLY : pos, full_phs_path_flag
  !
  USE input_parameters, ONLY : nconstr_inp, ncolvar_inp
  !
  USE constraints_module,    ONLY : init_constraint
  USE metadyn_vars,          ONLY : init_metadyn_vars
  USE read_namelists_module, ONLY : read_namelists, sm_not_set
  USE london_module,         ONLY : init_london, lon_rcut, scal6
  USE us, ONLY : spline_ps_ => spline_ps
  !
  IMPLICIT NONE
  !
  INTEGER  :: ia, image, nt
  REAL(DP) :: theta, phi
  !
  !
  IF ( ionode ) THEN
    CALL input_from_file()
    WRITE(stdout, '(5x,a)') "Waiting for input..."
  ENDIF
  !
  ! ... all namelists are read
  !
  CALL read_namelists( 'PW' )

  !
  ! ... various initializations of control variables
  !
  lscf      = .FALSE.
  lmd       = .FALSE.
  lmetadyn  = .FALSE.
  lpath     = .FALSE.
  lneb      = .FALSE.
  lsmd      = .FALSE.
  lmovecell = .FALSE.
  lphonon   = .FALSE.
  lbands    = .FALSE.
  lbfgs     = .FALSE.
  ldamped   = .FALSE.
  lforce    = tprnfor
  calc      = ' '
  !
  SELECT CASE( TRIM( calculation ) )
  CASE( 'scf' )
     !
     lscf  = .TRUE.
     nstep = 1
     !
  CASE( 'nscf' )
     !
     lforce = .FALSE.
     nstep  = 1
     !
  CASE( 'bands' )
     !
     lforce = .FALSE.
     lbands = .TRUE.
     nstep  = 1
     !
  CASE( 'phonon' )
     !
     CALL errore ('iosys','calculation="phonon" no longer available',1)
     lforce  = .FALSE.
     lphonon = .TRUE.
     !
     nstep = 1
     !
  CASE( 'relax' )
     !
     lscf   = .TRUE.
     lforce = .TRUE.
     !
     epse = etot_conv_thr
     epsf = forc_conv_thr
     !
     SELECT CASE( TRIM( ion_dynamics ) )
     CASE( 'bfgs' )
        !
        lbfgs = .TRUE.
        !
        IF ( epse <= 20.D0 * ( tr2 / upscale ) ) &
           CALL errore( 'iosys', 'required etot_conv_thr is too small:' // &
                      & ' conv_thr must be reduced', 1 )
        !
     CASE ( 'damp' )
        !
        lmd     = .TRUE.
        ldamped = .TRUE.
        !
        ntcheck = nstep + 1
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
  CASE( 'md' )
     !
     lscf   = .TRUE.
     lmd    = .TRUE.
     lforce = .TRUE.
     !
     SELECT CASE( TRIM( ion_dynamics ) )
     CASE( 'verlet' )
        !
        CONTINUE
        !
     CASE( 'langevin' )
        !
        llang       = .TRUE.
        temperature = tempw
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
     END SELECT
     !
  CASE( 'vc-relax' )
     !
     lscf      = .TRUE.
     lmd       = .TRUE.
     lmovecell = .TRUE.
     lforce    = .TRUE.
     ldamped   = .TRUE.
     !
     epse =  etot_conv_thr
     epsf =  forc_conv_thr
     epsp = press_conv_thr
     !
     SELECT CASE( TRIM( cell_dynamics ) )
     CASE( 'none' )
        !
        calc    = 'mm'
        ntcheck = nstep + 1
        !
     CASE( 'damp-pr' )
        !
        calc    = 'cm'
        ntcheck = nstep + 1
        !
     CASE( 'damp-w' )
        !
        calc    = 'nm'
        ntcheck = nstep + 1
        !
     CASE( 'bfgs' )
        !
        lbfgs = .TRUE.
        lmd   = .FALSE.
        ldamped = .FALSE.
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys', 'calculation=' // TRIM( calculation ) // &
                   & ': cell_dynamics=' // TRIM( cell_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
     IF ( .NOT. ldamped .AND. .NOT. lbfgs) &
        CALL errore( 'iosys', 'calculation='// TRIM( calculation ) // &
                   & ': incompatible ion (' // TRIM( ion_dynamics )// &
                   & ') and cell dynamics ('// TRIM(cell_dynamics )// ')', 1 )
     !
  CASE( 'vc-md' )
     !
     lscf      = .TRUE.
     lmd       = .TRUE.
     lmovecell = .TRUE.
     lforce    = .TRUE.
     !
     ntcheck = nstep + 1
     !
     SELECT CASE( TRIM( cell_dynamics ) )
     CASE( 'none' )
        !
        calc = 'md'
        !
     CASE( 'pr' )
        !
        calc = 'cd'
        !
     CASE( 'w' )
        !
        calc = 'nd'
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
     IF ( TRIM( ion_dynamics ) /= 'beeman' ) &
        CALL errore( 'iosys', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
     !
  CASE( 'neb' )
     !
     lscf  = .TRUE.
     lpath = .TRUE.
     lneb  = .TRUE.
     !
  CASE( 'smd' )
     !
     lscf  = .TRUE.
     lpath = .TRUE.
     lsmd  = .TRUE.
     !
  CASE( 'metadyn' )
     !
     lscf    = .TRUE.
     lmetadyn= .TRUE.
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', 'calculation ' // &
                & TRIM( calculation ) // ' not implemented', 1 )
     !
  END SELECT
  !
  ! ... translate from input to internals of PWscf, various checks
  !
  if (input_dft /='none') call enforce_input_dft (input_dft)
  !
  IF ( tefield .AND. ( .NOT. nosym ) ) THEN
     nosym = .TRUE.
     WRITE( stdout, &
            '(5x,"Presently no symmetry can be used with electric field",/)' )
  END IF
  IF ( tefield .AND. tstress ) THEN
     tstress = .FALSE.
     WRITE( stdout, &
            '(5x,"Presently stress not available with electric field",/)' )
  END IF
  IF ( tefield .AND. ( nspin > 2 ) ) THEN
     CALL errore( 'iosys', 'LSDA not available with electric field' , 1 )
  END IF
  !
  twfcollect = wf_collect
  !
  ! ... Set Values for electron and bands
  !
  tfixed_occ = .FALSE.
  ltetra     = .FALSE.
  !
  SELECT CASE( TRIM( occupations ) )
  CASE( 'fixed' )
     !
     ngauss = 0
     !
     IF ( degauss /= 0.D0 ) THEN
        CALL errore( ' iosys ', &
                   & ' fixed occupations, gauss. broadening ignored', -1 )
        degauss = 0.D0
     END IF
     !
  CASE( 'smearing' )
     !
     IF ( degauss == 0.D0 ) &
        CALL errore( ' iosys ', &
                   & ' smearing requires gaussian broadening', 1 )
     !
     SELECT CASE ( TRIM( smearing ) )
     CASE ( 'gaussian', 'gauss' )
        ngauss = 0
     CASE ( 'methfessel-paxton', 'm-p', 'mp' )
        ngauss = 1
     CASE ( 'marzari-vanderbilt', 'cold', 'm-v', 'mv' )
        ngauss = -1
     CASE ( 'fermi-dirac', 'f-d', 'fd' )
        ngauss = -99
     END SELECT
     !
  CASE( 'tetrahedra' )
     !
     ngauss = 0
     ltetra = .TRUE.
     !
  CASE( 'from_input' )
     !
     ngauss     = 0
     tfixed_occ = .TRUE.
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys','occupations ' // TRIM( occupations ) // &
                & 'not implemented', 1 )
     !
  END SELECT
  !
  IF( nbnd < 1 ) &
     CALL errore( 'iosys', 'nbnd less than 1', nbnd )
  !
  IF( nelec < 0_DP ) &
     CALL errore( 'iosys', 'nelec less than 0', 1 )
  !
  IF ( nelup < 0_DP ) &
     CALL errore( 'iosys', 'nelup less than 0', 1 )
  !
  IF ( neldw < 0_DP ) &
     CALL errore( 'iosys', 'neldw less than 0', 1 )
  !
  SELECT CASE( nspin )
  CASE( 1 )
     !
     lsda = .FALSE.
     IF ( noncolin ) nspin = 4
     !
  CASE( 2 )
     !
     lsda = .TRUE.
     IF ( noncolin ) &
        CALL errore( 'iosys', &
                     'noncolin .and. nspin==2 are conflicting flags', 1 )
     !
  CASE( 4 )
     !
     lsda = .FALSE.
     noncolin = .TRUE.
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', 'wrong input value for nspin', 1 )
     !
  END SELECT
  !
  IF ( nelup == 0.D0 .AND. neldw == 0.D0 .AND. &
       tot_magnetization < 0 .AND. multiplicity == 0) THEN
     !
     two_fermi_energies = .FALSE.
     !
  ELSE
     !
     two_fermi_energies = .TRUE.
     !
     IF ( .NOT. lsda ) &
        CALL errore( 'iosys', 'fixed nelup/neldw requires nspin=2', 1 )
     !
     IF ( tot_magnetization < 0 .AND. multiplicity == 0 .AND. &
          ABS( nelup + neldw - nelec ) > 1.D-10 ) &
        CALL errore( 'iosys', 'nelup + neldw must be equal to nelec', 1 )
     !
  END IF
  !
  ! ... starting_magnetization(nt) = sm_not_set means "not set" 
  ! ... stop if starting_magnetization is not set for at least
  ! ... one atomic type and occupations are not set in any other way
  !
  IF ( lscf .AND. nspin == 2 .AND. .NOT. tfixed_occ .AND. &
          nelup == 0.d0 .AND. neldw == 0.d0 .AND. &
          multiplicity == 0 .AND. tot_magnetization == -1    .AND. &
          ALL(starting_magnetization == sm_not_set) ) THEN
      CALL errore('iosys','some starting_magnetization MUST be set', 1 )
  END IF
  !
  DO nt = 1, ntyp
     !
     IF ( starting_magnetization(nt) == sm_not_set ) THEN
        starting_magnetization(nt) = 0.0_dp
     ELSE IF ( starting_magnetization(nt) > 1.0_dp ) THEN
        starting_magnetization(nt) = 1.0_dp
     ELSE IF ( starting_magnetization(nt) <-1.0_dp ) THEN
        starting_magnetization(nt) =-1.0_dp
     END IF
     !
  END DO
  !
  IF (noncolin) THEN
     DO nt = 1, ntyp
        !
        angle1(nt) = pi * angle1(nt) / 180.D0
        angle2(nt) = pi * angle2(nt) / 180.D0
        !
     END DO
  ELSE
     angle1=0.d0
     angle2=0.d0
  ENDIF
  !
  SELECT CASE( TRIM( constrained_magnetization ) )
  CASE( 'none' )
     !
     i_cons = 0
     !
  CASE( 'total' )
     !
     IF ( nspin == 4 ) THEN
        !
        i_cons = 3
        !
        mcons(1,1) = fixed_magnetization(1)
        mcons(2,1) = fixed_magnetization(2)
        mcons(3,1) = fixed_magnetization(3)
        !
     ELSE IF ( nspin == 2 ) THEN
        !
        i_cons = 5
        !
        two_fermi_energies = .TRUE.
        !
        mcons(3,1) = fixed_magnetization(3)
        !
        IF ( fixed_magnetization(1) /= 0.D0 .OR. &
             fixed_magnetization(2) /= 0.D0 ) &
           CALL errore( 'iosys', 'only fixed_magnetization(3)' // &
                      & ' can be specified with nspin=2 ', 1 )
        !
     ELSE
        !
        CALL errore( 'iosys','constrained total magnetization ' // &
                   & 'requires nspin=2 or 4 ', 1 )
        !
     END IF
     !
  CASE( 'atomic' )
     !
     IF ( nspin == 1 ) &
        CALL errore( 'iosys','constrained atomic magnetizations ' // &
                   & 'require nspin=2 or 4 ', 1 )
     !
     i_cons = 1
     !
     if (nspin == 4) then
        ! non-collinear case
        DO nt = 1, ntyp
           !
           theta = angle1(nt)
           phi   = angle2(nt)
           !
           mcons(1,nt) = starting_magnetization(nt) * SIN( theta ) * COS( phi )
           mcons(2,nt) = starting_magnetization(nt) * SIN( theta ) * SIN( phi )
           mcons(3,nt) = starting_magnetization(nt) * COS( theta )
           !
        END DO
     else
        ! collinear case
        DO nt = 1, ntyp
           !
           mcons(1,nt) = starting_magnetization(nt) 
           !
        END DO
     end if
     !
  CASE( 'total direction' )
     i_cons = 6
     mcons(3,1) = fixed_magnetization(3)
     !
  CASE( 'atomic direction' )
     !
     IF ( nspin == 1 ) &
        CALL errore( 'iosys','constrained atomic magnetization ' // &
                   & 'directions require nspin=2 or 4 ', 1 )
     !
     i_cons = 2
     !
     DO nt = 1, ntyp
        !
        theta = angle1(nt)
        !
        mcons(3,nt) = cos(theta)
        !
     END DO
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys','constrained magnetization ' // &
                & TRIM( constrained_magnetization ) // 'not implemented', 1 )
     !
  END SELECT
  !

  IF ( B_field(1) /= 0.D0 .OR. &
       B_field(2) /= 0.D0 .OR. &
       B_field(3) /= 0.D0 ) THEN
     CALL errore( 'iosys', &
                & 'B field currently not implemented', 1 )
     !
     !IF ( nspin == 1 ) &
     !   CALL errore( 'iosys', &
     !              & 'non-zero external B_field requires nspin=2 or 4', 1 )
     !
     !IF ( TRIM( constrained_magnetization ) /= 'none' ) &
     !   CALL errore( 'iosys', 'constrained_magnetization and ' // &
     !              & 'non-zero external B_field are conflicting flags', 1 )
     !
     !IF ( nspin == 2 .AND. &
     !     ( B_field(1) /= 0.D0 .OR. B_field(2) /= 0.D0 ) ) &
     !   CALL errore( 'iosys', &
     !              & 'only B_field(3) can be specified with nspin=2', 1 )
     !
  END IF
  !
  IF ( occupations == 'fixed' .AND. nspin == 2  .AND. lscf ) THEN
     !
     IF ( two_fermi_energies ) THEN
        !
        IF ( ABS( NINT( nelup ) - nelup ) > 1.D-10 ) &
           CALL errore( 'iosys', &
                      & 'fixed occupations requires integer nelup', 1 )
        IF ( ABS( NINT( neldw ) - neldw ) > 1.D-10 ) &
           CALL errore( 'iosys', &
                      & 'fixed occupations requires integer neldw', 1 )
        !
     ELSE
        !
        CALL errore( 'iosys', &
                   & 'fixed occupations and lsda need nelup and neldw', 1 )
        !
     END IF
     !
  END IF
  !
  IF ( ecutrho <= 0.D0 ) THEN
     !
     dual = 4.D0
     !
  ELSE
     !
     dual = ecutrho / ecutwfc
     !
     IF ( dual <= 1.D0 ) &
        CALL errore( 'iosys', 'invalid dual?', 1 )
     !
  END IF
  !
  SELECT CASE( TRIM( restart_mode ) )
  CASE( 'from_scratch' )
     !
     restart        = .FALSE.
     IF ( lscf ) THEN
        startingconfig = 'input'
     ELSE
        startingconfig = 'file'
     END IF
     !
  CASE( 'restart' )
     !
     IF ( lneb .OR. lsmd ) THEN
        !
        ! ... "path" specific
        !
        restart = .FALSE.
        !
     ELSE
        !
        restart = .TRUE.
        !
        IF ( TRIM( ion_positions ) == 'from_input' ) THEN
           !
           startingconfig = 'input'
           !
        ELSE
           !
           startingconfig = 'file'
           !
        END IF
        !
     END IF
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', &
                & 'unknown restart_mode ' // TRIM( restart_mode ), 1 )
     !
  END SELECT
  !
  SELECT CASE( TRIM( disk_io ) )
  CASE( 'high' )
     !
     io_level = 2
     !
  CASE ( 'low' )
     !
     io_level = 0
     restart  = .FALSE.
     !
  CASE ( 'none' )
     !
     io_level = -1
     restart  = .FALSE.
     IF ( twfcollect ) THEN
        CALL infomsg('iosys', 'minimal I/O required, wf_collect reset to FALSE')
        twfcollect= .FALSE.
     END IF
     !
  CASE DEFAULT
     !
     io_level = 1
     !
     IF ( lscf ) restart  = .FALSE.
     !
  END SELECT
  !
  Hubbard_U(:)    = Hubbard_U(:) / rytoev
  Hubbard_alpha(:)= Hubbard_alpha(:) / rytoev
  !
  ethr = diago_thr_init
  !
  IF ( startingpot /= 'atomic' .AND. startingpot /= 'file' ) THEN
     !
     CALL infomsg( 'iosys', 'wrong startingpot: use default (1)' )
     !
     IF (       lscf ) startingpot = 'atomic'
     IF ( .NOT. lscf ) startingpot = 'file'
     !
  END IF
  !
  IF ( .NOT. lscf .AND. startingpot /= 'file' ) THEN
     !
     CALL infomsg( 'iosys', 'wrong startingpot: use default (2)' )
     !
     startingpot = 'file'
     !
  END IF
  !
  IF (      startingwfc /= 'atomic' .AND. &
            startingwfc /= 'random' .AND. &
            startingwfc /= 'atomic+random' .AND. &
            startingwfc /= 'file' ) THEN
     !
     CALL infomsg( 'iosys', 'wrong startingwfc: use default' )
     !
     startingwfc = 'atomic'
     !
  END IF
  !
  ! parallel distributed diagonalization is the default
  ! (with more than 4 procs)
  !
  use_para_diag = .TRUE.
  IF( ortho_para_ < 1 ) ortho_para_ = ortho_para
  !
  SELECT CASE( TRIM( diagonalization ) )
  CASE ( 'cg' )
     !
     isolve = 1
     max_cg_iter = diago_cg_maxiter
     !
  CASE ( 'cg-serial', 'cg+serial' )
     !
     isolve = 1
     max_cg_iter = diago_cg_maxiter
     use_para_diag = .FALSE.
     !
  CASE ( 'david', 'davidson', 'david+para', 'david+distpara' )
     !
     isolve = 0
     david = diago_david_ndim
     !
  CASE ( 'david+serial', 'david-serial' )
     !
     isolve = 0
     david = diago_david_ndim
     use_para_diag = .FALSE.
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', 'diagonalization ' // &
                & TRIM( diagonalization ) // ' not implemented', 1 )
     !
  END SELECT
  !
  tr2   = conv_thr
  niter = electron_maxstep
  !
  SELECT CASE( TRIM( pot_extrapolation ) )
  CASE( 'from_wfcs', 'from-wfcs' )
     !
     pot_order =-1
     !
  CASE( 'none' )
     !
     pot_order = 0
     !
  CASE( 'atomic' )
     !
     pot_order = 1
     !
  CASE( 'first_order', 'first-order', 'first order' )
     !
     pot_order = 2
     !
  CASE( 'second_order', 'second-order', 'second order' )
     !
     pot_order = 3
     !
  CASE DEFAULT
     !
     pot_order = 1
     !
  END SELECT
  !
  SELECT CASE( TRIM( wfc_extrapolation ) )
  CASE( 'none' )
     !
     wfc_order = 0
     !
  CASE( 'first_order', 'first-order', 'first order' )
     !
     wfc_order = 2
     !
  CASE( 'second_order', 'second-order', 'second order' )
     !
     wfc_order = 3
     !
  CASE DEFAULT
     !
     wfc_order = 0
     !
  END SELECT
  !
  lcoarsegrained  = lmetadyn .OR. ( TRIM( phase_space ) == 'coarse-grained' )
  !
  IF ( lcoarsegrained ) THEN
     !
     lmd = .TRUE.
     !
     SELECT CASE( TRIM( ion_dynamics ) )
     CASE( 'verlet' )
        !
        CONTINUE
        !
     CASE( 'langevin' )
        !
        llang = .TRUE.
        !
     CASE( 'damp' )
        !
        ldamped = .TRUE.
        !
        epse = etot_conv_thr
        epsf = forc_conv_thr
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
  END IF
  !
  ! ... "path" specific initialization of control variables
  !
  IF ( lpath ) THEN
     !
     nstep_path = nstep
     !
     IF ( num_of_images < 2 ) &
        CALL errore( 'iosys', 'calculation=' // TRIM( calculation ) // &
                   & ': num_of_images must be at least 2', 1 )
     !
     IF ( ( CI_scheme /= "no-CI"  ) .AND. &
          ( CI_scheme /= "auto"   ) .AND. &
          ( CI_scheme /= "manual" ) ) THEN
        !
        CALL errore( 'iosys', 'calculation=' // TRIM( calculation ) // &
                   & ': unknown CI_scheme', 1 )
        !
     END IF
     !
     ! ... initialization of logical variables
     !
     lsteep_des  = .FALSE.
     lquick_min  = .FALSE.
     lbroyden    = .FALSE.
     !
     SELECT CASE( opt_scheme )
     CASE( "sd" )
        !
        lsteep_des = .TRUE.
        !
     CASE( "quick-min" )
        !
        lquick_min = .TRUE.
        !
     CASE( "broyden" )
        !
        lbroyden     = .TRUE.
        !
     CASE( "langevin" )
        !
        llangevin = .TRUE.
        !
        IF ( lneb ) &
           CALL errore( 'iosys','calculation=' // TRIM( calculation ) // &
                      & ': langevin dynamics not implemented', 1 )
        !
        temp_req = temp_req / ( eV_to_kelvin * autoev )
        !
        IF ( temp_req <= 0.D0 ) &
           CALL errore( 'iosys','calculation=' // TRIM( calculation ) // &
                      & ': tepm_req has not been set', 1 )
        !
        IF ( use_freezing ) &
           WRITE( UNIT = stdout, &
                  FMT = '(5X,"warning: freezing cannot be used in langevin")' )
        !
        use_freezing = .FALSE.
        !
     CASE default
        !
        CALL errore( 'iosys','calculation=' // TRIM( calculation ) // &
                   & ': unknown opt_scheme', 1 )
        !
     END SELECT
     !
  END IF
  !
  SELECT CASE( TRIM( ion_temperature ) )
  CASE( 'not_controlled', 'not-controlled', 'not controlled' )
     !
     control_temp = .FALSE.
     !
  CASE( 'rescaling' )
     !
     control_temp = .TRUE.
     thermostat   = TRIM( ion_temperature )
     temperature  = tempw
     tolp_        = tolp
     !
  CASE( 'rescale-v', 'rescale-V', 'rescale_v', 'rescale_V' )
     !
     control_temp = .TRUE.
     thermostat   = TRIM( ion_temperature )
     temperature  = tempw
     nraise_      = nraise
     !
  CASE( 'reduce-T', 'reduce-t', 'reduce_T', 'reduce_t' )
     !
     control_temp = .TRUE.
     thermostat   = TRIM( ion_temperature )
     temperature  = tempw
     delta_t_     = delta_t
     nraise_      = nraise
     !
  CASE( 'rescale-T', 'rescale-t', 'rescale_T', 'rescale_t' )
     !
     control_temp = .TRUE.
     thermostat   = TRIM( ion_temperature )
     temperature  = tempw
     delta_t_     = delta_t
     !
  CASE( 'berendsen', ' Berendsen' )
     !
     control_temp = .TRUE.
     thermostat   = TRIM( ion_temperature )
     temperature  = tempw
     nraise_      = nraise
     !
  CASE( 'andersen', 'Andersen' )
     !
     control_temp = .TRUE.
     thermostat   = TRIM( ion_temperature )
     temperature  = tempw
     nraise_      = nraise
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', &
                & 'unknown ion_temperature ' // TRIM( ion_temperature ), 1 )
     !
  END SELECT
  !
  SELECT CASE( TRIM( mixing_mode ) )
  CASE( 'plain' )
     !
     imix = 0
     !
  CASE( 'TF' )
     !
     imix = 1
     !
  CASE( 'local-TF' )
     !
     imix = 2
     !
  CASE( 'potential' )
     !
     CALL errore( 'iosys', 'potential mixing no longer implemented', 1 )
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', 'unknown mixing ' // TRIM( mixing_mode ), 1 )
     !
  END SELECT
  !
  starting_scf_threshold = tr2
  nmix = mixing_ndim
  niter_with_fixed_ns = mixing_fixed_ns
  !
  SELECT CASE( TRIM( verbosity ) )
  CASE( 'high' )
     !
     iverbosity = 1
     !
  CASE DEFAULT
     !
     iverbosity = 0
     !
  END SELECT
  !
  tmp_dir = trimcheck ( outdir )
  lstres = ( tstress .AND. lscf )
  !
  IF ( lberry .AND. npool > 1 ) &
     CALL errore( 'iosys', 'Berry Phase not implemented with pools', 1 )
  IF ( lberry .AND. noncolin ) &
     CALL errore( 'iosys', 'Noncolinear Berry Phase not implemented', 1 )
  !
  ! ... Copy values from input module to PW internals
  !
  nppstr_     = nppstr
  gdir_       = gdir
  lberry_     = lberry
  lelfield_   = lelfield
  efield_     = efield
  nberrycyc_  = nberrycyc
  efield_cart_ = efield_cart
  tqr_        = tqr
  real_space_ = real_space
  !
  title_      = title
  lkpoint_dir_=lkpoint_dir
  dt_         = dt
  tefield_    = tefield
  dipfield_   = dipfield
  prefix_     = TRIM( prefix )
  pseudo_dir_ = TRIM( pseudo_dir )
  nstep_      = nstep
  iprint_     = iprint
  !
  celldm_  = celldm
  ibrav_   = ibrav
  nat_     = nat
  ntyp_    = ntyp
  edir_    = edir
  emaxpos_ = emaxpos
  eopreg_  = eopreg
  eamp_    = eamp
  nr1_     = nr1
  nr2_     = nr2
  nr3_     = nr3
  ecutwfc_ = ecutwfc
  ecfixed_ = ecfixed
  qcutz_   = qcutz
  q2sigma_ = q2sigma
  nr1s_    = nr1s
  nr2s_    = nr2s
  nr3s_    = nr3s
  degauss_ = degauss
  nelec_   = nelec
  nelup_   = nelup
  neldw_   = neldw
  !
  tot_charge_        = tot_charge
  tot_magnetization_ = tot_magnetization
  multiplicity_      = multiplicity
  !
  lspinorb_ = lspinorb
  noncolin_ = noncolin
  angle1_   = angle1
  angle2_   = angle2
  report_   = report
  lambda_   = lambda
  !
  spline_ps_ = spline_ps
  !
  Hubbard_U_(1:ntyp)      = hubbard_u(1:ntyp)
  Hubbard_alpha_(1:ntyp)  = hubbard_alpha(1:ntyp)
  lda_plus_u_             = lda_plus_u
  la2F_                   = la2F
  nspin_                  = nspin
  starting_magnetization_ = starting_magnetization
  starting_ns             = starting_ns_eigenvalue
  U_projection            = U_projection_type
  noinv_                  = noinv
  nosym_                  = nosym
  nosym_evc_              = nosym_evc
  nofrac                  = force_symmorphic
  nbnd_                   = nbnd
  !
#if defined (EXX)
  !
  x_gamma_extrapolation_ = x_gamma_extrapolation
  !
  nqx1_ = nqx1
  nqx2_ = nqx2
  nqx3_ = nqx3
  !
#endif
  !
  diago_full_acc_ = diago_full_acc
  starting_wfc    = startingwfc
  starting_pot    = startingpot
  mixing_beta_    = mixing_beta
  !
  remove_rigid_rot_ = remove_rigid_rot
  upscale_          = upscale
  refold_pos_       = refold_pos
  press_            = press
  cell_factor_      = cell_factor
  !
  ! ... "path"-optimization variables
  !
  ds_             = ds
  num_of_images_  = num_of_images
  first_last_opt_ = first_last_opt
  use_masses_     = use_masses
  use_freezing_   = use_freezing
  temp_req_       = temp_req
  path_thr_       = path_thr
  CI_scheme_      = CI_scheme
  k_max_          = k_max
  k_min_          = k_min
  fixed_tan_      = fixed_tan
  !
  ! ... for WANNIER_AC
  use_wannier_ = use_wannier
  use_energy_int_ = use_energy_int
  nwan_ = nwan
  print_wannier_coeff_ = print_wannier_coeff
  !
  !
  ! ... BFGS specific
  !
  bfgs_ndim_        = bfgs_ndim
  trust_radius_max_ = trust_radius_max
  trust_radius_min_ = trust_radius_min
  trust_radius_ini_ = trust_radius_ini
  w_1_              = w_1
  w_2_              = w_2
! DCC
  !
  ! ...  Charge Compensation
  !
  which_compensation_ = which_compensation
  ecutcoarse_ = ecutcoarse
  mixing_charge_compensation_ = mixing_charge_compensation 
  n_charge_compensation_ = n_charge_compensation 
  comp_thr_ = comp_thr
  nlev_ = nlev
  ! more DCC variables
  mr1_ = mr1
  mr2_ = mr2
  mr3_ = mr3
  whichbc_ = whichbc
  errtol_ = errtol
  itmax_ = itmax
  ncompx_ = ncompx
  ncompy_ = ncompy
  ncompz_ = ncompz
  cellmin_ = cellmin
  cellmax_ = cellmax
  !
  llondon     = london
  lon_rcut    = london_rcut
  scal6       = london_s6
  !
  SELECT CASE( TRIM( which_compensation ) )
      !
    CASE( 'dcc' )
      !
      which_compensation_ = 'dcc'
      do_comp_ = .TRUE.
      do_coarse_ = .FALSE.
      do_mltgrid_ = .TRUE.
      !
    CASE( 'dcc+' )
      !
      which_compensation_ = 'dcc+'
      do_comp_ = .TRUE.
      do_coarse_ = .FALSE.
      do_mltgrid_ = .FALSE.
      !
    CASE DEFAULT
      !
      do_comp_ = .FALSE.
      do_coarse_ = .FALSE.
      do_mltgrid_ = .FALSE.
      !
  END SELECT
  !
  ! ... read following cards
  !
  ALLOCATE( ityp( nat_ ) )
  ALLOCATE( tau(    3, nat_ ) )
  ALLOCATE( force(  3, nat_ ) )
  ALLOCATE( if_pos( 3, nat_ ) )
  ALLOCATE( extfor( 3, nat_ ) )
  IF ( tfixed_occ ) THEN
     IF ( nspin_ == 4 ) THEN
        ALLOCATE( f_inp( nbnd_, 1 ) )
     ELSE
        ALLOCATE( f_inp( nbnd_, nspin_ ) )
     END IF
  END IF
  !
  IF ( tefield ) ALLOCATE( forcefield( 3, nat_ ) )
  !
  CALL read_cards( psfile, atomic_positions )
  !
  ! ... set up atomic positions and crystal lattice
  !
  IF ( celldm_(1) == 0.D0 .AND. a /= 0.D0 ) THEN
     !
     celldm_(1) = a / bohr_radius_angs
     celldm_(2) = b / a
     celldm_(3) = c / a
     !
     IF ( ibrav_ /= 14 ) THEN
        !
        ! ... trigonal and monoclinic lattices
        !
        celldm_(4) = cosab
     ELSE
        !
        ! ... triclinic lattice
        !
        celldm_(4) = cosbc
        celldm_(5) = cosac
        celldm_(6) = cosab
     END IF
     !
  ELSE IF ( celldm_(1) /= 0.D0 .AND. a /= 0.D0 ) THEN
     !
     CALL errore( 'input', 'do not specify both celldm and a,b,c!', 1 )
     !
  END IF
  !
  ! ... generate at (in atomic units) from ibrav and celldm
  !
  CALL latgen( ibrav_, celldm_, at(1,1), at(1,2), at(1,3), omega )
  !
  ! ... define alat
  !
  alat = celldm_(1)
  !
  ! ... convert at to unit of alat
  !
  at = at / alat
  !
  ! ... Generate the reciprocal lattice vectors
  !
  CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
  !
  IF ( full_phs_path_flag ) THEN
     !
     ! ... "path" optimizations specific
     !
     DO image = 1, num_of_images_
        !
        tau = RESHAPE( pos(1:3*nat_,image), (/ 3 , nat_ /) )
        !
        CALL convert_tau ( atomic_positions, nat_, tau)
        !
        ! ... note that this positions array is in Bohr
        !
        pos(1:3*nat_,image) = RESHAPE( tau, (/ 3 * nat_ /) ) * alat
        !
     END DO
     !
  ELSE
     !
     CALL convert_tau ( atomic_positions, nat_, tau)
     !
  END IF
  !
  IF ( wmass == 0.D0 ) THEN
     !
     ! ... set default value of wmass
     !
#if defined __PGI
     DO ia = 1, nat_
        wmass = wmass + amass( ityp(ia) )
     END DO
#else
     wmass = SUM( amass(ityp(:)) )
#endif
     !
     wmass = wmass * amconv
     IF ( calc == 'nd' .OR. calc == 'nm' ) THEN
        wmass = 0.75D0 * wmass / pi / pi / omega**( 2.D0 / 3.D0 )
     ELSE IF ( calc == 'cd' .OR. calc == 'cm' ) THEN
        wmass = 0.75D0 * wmass / pi / pi
     END IF
     !
     cmass  = wmass
     !
  ELSE
     !
     ! ... wmass is given in amu, Renata's dynamics uses masses in atomic units
     !
     cmass  = wmass * amconv
     !
  END IF
  !
  ! ... unit conversion for pressure
  !
  press_ = press_ / uakbar
  !
  ! ... set constraints for cell dynamics/optimization
  !
  CALL init_dofree ( cell_dofree ) 
  !
  ! ... read pseudopotentials
  !
  CALL readpp()
  !
  ! ... if DFT finite size corrections are needed, define the appropriate volume
  !
  if (dft_has_finite_size_correction()) &
      call set_finite_size_volume(REAL(omega*nk1*nk2*nk3))
  !
  ! ... In the case of variable cell dynamics save old cell variables
  ! ... and initialize a few other variables
  !
  IF ( lmovecell ) THEN
     !
     at_old    = at
     omega_old = omega
     lstres    = .TRUE.
     !
     IF ( cell_factor_ <= 0.D0 ) cell_factor_ = 1.2D0
     !
     IF ( cmass <= 0.D0 ) &
        CALL errore( 'iosys', &
                   & 'vcsmd: a positive value for cell mass is required', 1 )
     !
  ELSE
     !
     cell_factor_ = 1.D0
     !
  END IF
  !
  !
  ! ... allocate arrays for dispersion correction
  !
  IF ( llondon) CALL init_london ( )
  !
  ! ... variables for constrained dynamics are set here
  !
  lconstrain = ( ncolvar_inp + nconstr_inp > 0 )
  !
  IF ( lconstrain ) CALL init_constraint( nat, tau, ityp, alat )
  !
  ! ... set variables for metadynamics
  !
  IF ( lcoarsegrained ) CALL init_metadyn_vars()
  !
  ! ... read atomic positions and unit cell from data file
  ! ... must be done before "verify_tmpdir" because the latter
  ! ... removes the data file in a run from scratch
  !
  IF ( startingconfig == 'file' ) CALL read_config_from_file()
  !
  CALL verify_tmpdir( tmp_dir )
  !
  IF ( .NOT. TRIM( wfcdir ) == 'undefined' ) THEN
     !
     wfc_dir = trimcheck ( wfcdir )
     !
     CALL verify_tmpdir( wfc_dir )
     !
  ENDIF
  !
  CALL restart_from_file()
  !
  ! ... Files
  !
  input_drho  = ' '
  output_drho = ' '
  ! 
  if (real_space ) then
   WRITE( stdout, '(5x,"Real space treatment of Beta functions, &
         &V.1 (BE SURE TO CHECK MANUAL!)",/)' )
  endif
  RETURN
  !
END SUBROUTINE iosys
!
!----------------------------------------------------------------------------
SUBROUTINE read_cards( psfile, atomic_positions_ )
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE input_parameters,   ONLY : atom_label, atom_pfile, atom_mass, taspc, &
                                 tapos, rd_pos, atomic_positions, if_pos,  &
                                 sp_pos, k_points, xk, wk, nk1, nk2, nk3,  &
                                 k1, k2, k3, nkstot, cell_symmetry, rd_ht, &
                                 trd_ht, f_inp, rd_for
  USE cell_base,          ONLY : at, ibrav, symm_type
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp, tau, atm, extfor
  USE klist,              ONLY : nkstot_ => nkstot
  USE ktetra,             ONLY : nk1_   => nk1, &
                                 nk2_   => nk2, &
                                 nk3_   => nk3, &
                                 k1_    => k1,  &
                                 k2_    => k2,  &
                                 k3_    => k3
  USE klist,              ONLY : lxkcry, &
                                 xk_    => xk, &
                                 wk_    => wk
  USE fixed_occ,          ONLY : tfixed_occ, &
                                 f_inp_ => f_inp
  USE ions_base,          ONLY : fixatom, &
                                 if_pos_ =>  if_pos
  USE ions_base,          ONLY : amass
  USE control_flags,      ONLY : lfixatom, gamma_only, textfor
  USE read_cards_module,  ONLY : read_cards_base => read_cards
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=256) :: psfile(ntyp)
  CHARACTER (LEN=30)  :: atomic_positions_
  INTEGER, EXTERNAL :: atomic_number
  REAL(DP), EXTERNAL :: atom_weight
  !
  LOGICAL :: tcell = .FALSE.
  INTEGER :: is, ia
  !
  !
  amass = 0
  !
  CALL read_cards_base( 'PW' )
  !
  IF ( .NOT. taspc ) &
     CALL errore( 'read_cards', 'atomic species info missing', 1 )
  IF ( .NOT. tapos ) &
     CALL errore( 'read_cards', 'atomic position info missing', 1 )
  !
  DO is = 1, ntyp
     !
     amass(is)  = atom_mass(is)
     psfile(is) = atom_pfile(is)
     atm(is)    = atom_label(is)
     !
     IF ( amass(is) <= 0.0_DP ) amass(is)= &
              atom_weight(atomic_number(TRIM(atm(is))))

     IF ( amass(is) <= 0.D0 ) CALL errore( 'read_cards', 'invalid  mass', is )
     !
  END DO
  !
  textfor = .FALSE.
  IF( ANY( rd_for /= 0.0_DP ) ) textfor = .TRUE.
  !
  DO ia = 1, nat
     !
     tau(:,ia) = rd_pos(:,ia)
     ityp(ia)  = sp_pos(ia)
     extfor(:,ia) = rd_for(:,ia)
     !
  END DO
  !
  ! ... calculate fixatom
  !
  fixatom = 0
  !
  IF ( ANY( if_pos(:,1:nat) == 0 ) ) lfixatom = .TRUE.
  !
  DO ia = 1, nat
     !
     IF ( if_pos(1,ia) /= 0 .OR. &
          if_pos(2,ia) /= 0 .OR. &
          if_pos(3,ia) /= 0 ) CYCLE
     !
     fixatom = fixatom + 1
     !
  END DO
  !
  ! ... The constrain on fixed coordinates is implemented using the array
  ! ... if_pos whose value is 0 when the coordinate is to be kept fixed, 1
  ! ... otherwise. fixatom is maintained for compatibility. ( C.S. 15/10/2003 )
  !
  if_pos_(:,:) = if_pos(:,1:nat)
  !
  atomic_positions_ = TRIM( atomic_positions )
  !
  CALL reset_k_points ( )
  gamma_only = ( k_points == 'gamma' )
  !
  IF ( tfixed_occ ) THEN
     !
     IF ( nkstot_ > 1 .OR. ( nk1 * nk2 * nk3 ) > 1 ) &
        CALL errore( 'read_cards', &
                   & 'only one k point with fixed occupations', 1 )
     !
     f_inp_ = f_inp
     !
     DEALLOCATE ( f_inp )
     !
  END IF
  !
  IF ( trd_ht ) THEN
    !
    symm_type = cell_symmetry
    at        = TRANSPOSE( rd_ht )
    tcell     = .TRUE.
    !
  END IF
  !
  IF ( ibrav == 0 .AND. .NOT. tcell ) &
     CALL errore( 'read_cards', 'ibrav=0: must read cell parameters', 1 )
  IF ( ibrav /= 0 .AND. tcell ) &
     CALL errore( 'read_cards', 'redundant data for cell parameters', 2 )
  !
  RETURN
  !
END SUBROUTINE read_cards
!
!-----------------------------------------------------------------------
SUBROUTINE convert_tau (atomic_positions, nat_, tau)
!-----------------------------------------------------------------------
  !
  ! ... convert input atomic positions to internally used format:
  ! ... tau in a0 units
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : bohr_radius_angs
  USE cell_base,     ONLY : at, alat
  IMPLICIT NONE
  CHARACTER (LEN=*), INTENT(IN)  :: atomic_positions
  INTEGER, INTENT(IN)  :: nat_
  REAL (DP), INTENT(INOUT) :: tau(3,nat_)
  ! 
  SELECT CASE( atomic_positions )
  CASE( 'alat' )
     !
     ! ... input atomic positions are divided by a0: do nothing
     !
  CASE( 'bohr' )
     !
     ! ... input atomic positions are in a.u.: divide by alat
     !
     tau = tau / alat
     !
  CASE( 'crystal' )
     !
     ! ... input atomic positions are in crystal axis
     !
     CALL cryst_to_cart( nat_, tau, at, 1 )
     !
  CASE( 'angstrom' )
     !
     ! ... atomic positions in A: convert to a.u. and divide by alat
     !
     tau = tau / bohr_radius_angs / alat
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys','atomic_positions=' // &
                & TRIM( atomic_positions ) // ' not implemented', 1 )
     !
  END SELECT
  !
END SUBROUTINE convert_tau
!-----------------------------------------------------------------------
SUBROUTINE verify_tmpdir( tmp_dir )
  !-----------------------------------------------------------------------
  !
  USE wrappers,         ONLY : f_mkdir
  USE input_parameters, ONLY : restart_mode
  USE control_flags,    ONLY : lpath, lbands
  USE io_files,         ONLY : prefix, xmlpun, delete_if_present
  USE pw_restart,       ONLY : pw_readfile
  USE path_variables,   ONLY : num_of_images
  USE mp_global,        ONLY : mpime, nproc, nimage
  USE io_global,        ONLY : ionode
  USE mp,               ONLY : mp_barrier
  USE xml_io_base,      ONLY : copy_file
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(INOUT) :: tmp_dir
  !
  INTEGER             :: ios, image, proc, nofi
  LOGICAL             :: exst
  CHARACTER (LEN=256) :: file_path, tmp_dir_saved, filename
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  ios = 0
  !
  file_path = TRIM( tmp_dir ) // 'pwscf'
  !
  OPEN( UNIT = 4, FILE = TRIM( file_path ) // TRIM( int_to_char( mpime ) ), &
      & STATUS = 'UNKNOWN',  FORM = 'UNFORMATTED', IOSTAT = ios )
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !
  IF ( ios /= 0 ) CALL errore( 'outdir: ', TRIM( tmp_dir ) // &
                             & ' non existent or non writable', 1 )
  !
  ! ... if starting from scratch all temporary files are removed
  ! ... from tmp_dir ( only by the master node )
  !
  file_path = TRIM( tmp_dir ) // TRIM( prefix )
  !
  IF ( restart_mode == 'from_scratch' ) THEN
     !
     ! ... xml data file in save directory is removed
     !     but, header is read anyway to store qexml version
     !
     CALL pw_readfile( 'header', ios )
     !
     IF ( ionode ) THEN
        !
        IF ( .NOT. lbands ) THEN
            ! 
            ! save a bck copy of datafile.xml (AF)
            ! 
            filename = TRIM( file_path ) // '.save/' // TRIM( xmlpun )
            INQUIRE( FILE = filename, EXIST = exst )
            !
            IF ( exst ) CALL copy_file( TRIM(filename), TRIM(filename) // '.bck' )
            !
            CALL delete_if_present( TRIM(filename) )
            !
        ENDIF
        !
        ! ... extrapolation file is removed
        !
        CALL delete_if_present( TRIM( file_path ) // '.update' )
        !
        ! ... MD restart file is removed
        !
        CALL delete_if_present( TRIM( file_path ) // '.md' )
        !
        ! ... BFGS restart file is removed
        !
        CALL delete_if_present( TRIM( file_path ) // '.bfgs' )
        !
     END IF
     !
  END IF
  !
  ! ... "path" optimisation specific :   in the scratch directory the tree of
  ! ... subdirectories needed by "path" calculations are created
  !
#if defined(EXX)

     IF ( lpath .or. nimage > 1 ) THEN
#else
     IF ( lpath ) THEN
#endif

     IF ( ionode ) THEN
        !
        ! ... files needed by parallelization among images are removed
        !
        CALL delete_if_present( TRIM( file_path ) // '.BLOCK' )
        CALL delete_if_present( TRIM( file_path ) // '.newimage' )
        !
        ! ... file containing the broyden's history
        !
        IF ( restart_mode == 'from_scratch' ) THEN
           !
           CALL delete_if_present( TRIM( tmp_dir ) // &
                                 & TRIM( prefix ) // '.broyden' )
           !
        END IF
        !
     END IF
     !
     tmp_dir_saved = tmp_dir
     !
     nofi = num_of_images
#if defined(EXX)
     if( nimage > 1 ) nofi = nimage
#endif
     !

     DO image = 1, nofi
        !
        tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) //"_" // &
                  TRIM( int_to_char( image ) ) // '/'
        !
        DO proc = 0, nproc - 1
           !
           ! ... a scratch directory for this image is created sequentially
           ! ... by all the cpus
           !
           IF ( proc == mpime ) &
              ios = f_mkdir( TRIM( tmp_dir ) )
           !
           CALL mp_barrier()
           !
        END DO
        !
        ! ... each job checks whether the scratch directory is accessible
        ! ... or not
        !
        OPEN( UNIT = 4, FILE = TRIM( tmp_dir ) // TRIM( prefix ) // &
            & TRIM( int_to_char( mpime ) ), &
            & STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', IOSTAT = ios )
        CLOSE( UNIT = 4, STATUS = 'DELETE' )
        !
        IF ( ios /= 0 ) &
           CALL errore( 'outdir: ', TRIM( tmp_dir ) // &
                      & ' non existent or non writable', 1 )
        !
        ! ... if starting from scratch all temporary files are removed
        ! ... from tmp_dir ( by all the cpus in sequence )
        !
        IF ( restart_mode == 'from_scratch' ) THEN
           !
           DO proc = 0, nproc - 1
              !
              IF ( proc == mpime ) THEN
                 !
                 ! ... extrapolation file is removed
                 !
                 CALL delete_if_present( TRIM( tmp_dir ) // &
                                       & TRIM( prefix ) // '.update' )
                 !
                 ! ... standard output of the self-consistency is removed
                 !
                 CALL delete_if_present( TRIM( tmp_dir ) // 'PW.out' )
                 !
              END IF
              !
              CALL mp_barrier()
              !
           END DO
           !
        END IF
        !
     END DO
     !
     tmp_dir = tmp_dir_saved
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE verify_tmpdir


#if defined (AAAEXX)
SUBROUTINE iosys_image( )
  !
  CALL verify_tmpdir( tmp_dir )
  !
  IF ( .NOT. TRIM( wfcdir ) == 'undefined' ) THEN
     !
     wfc_dir = trimcheck ( wfcdir )
     !
     CALL verify_tmpdir( wfc_dir )
     !
  ENDIF
  !
  CALL restart_from_file()
  !
  IF ( startingconfig == 'file' ) CALL read_config_from_file()
  !
  RETURN
  !
END SUBROUTINE iosys_image
#endif
