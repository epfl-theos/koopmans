!
! Copyright (C) 2002-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE cprmain( tau_out, fion_out, etot_out )
  !----------------------------------------------------------------------------
  !
  USE kinds,                    ONLY : DP
  USE constants,                ONLY : bohr_radius_angs, amu_au, autoev
  USE control_flags,            ONLY : iprint, isave, thdyn, tpre, iprsta,     &
                                       tfor, remove_rigid_rot, taurdr,         &
                                       tprnfor, tsdc, lconstrain, lwf, lneb,   &
                                       lcoarsegrained, ndr, ndw, nomore, tsde, &
                                       tortho, tnosee, tnosep, trane, tranp,   &
                                       tsdp, tcp, tcap, ampre, amprp, tnoseh,  &
                                       tolp, ortho_eps, ortho_max, printwfc,   &
                                       tprojwfc, textfor, non_ortho
  USE core,                     ONLY : nlcc_any, rhoc
  USE uspp_param,               ONLY : nhm, nh
  USE cvan,                     ONLY : nvb, ish
  USE uspp,                     ONLY : nkb, vkb, becsum, deeq, okvan
  USE energies,                 ONLY : eht, epseu, exc, etot, eself, enl, &
                                       ekin, atot, entropy, egrand, enthal, &
                                       ekincm, print_energies, debug_energies
  USE electrons_base,           ONLY : nbspx, nbsp, ispin, f, nspin
  USE electrons_base,           ONLY : nel, nelt, iupdwn, nupdwn, nudx
  USE electrons_module,         ONLY : ei, sort_spreads, wfc_spreads, wfc_centers
  USE efield_module,            ONLY : efield, epol, tefield, allocate_efield, &
                                       efield_update, ipolp, qmat, gqq, evalue,&
                                       berry_energy, pberryel, pberryion,      &
                                       efield2, epol2, tefield2,               &
                                       allocate_efield2, efield_update2,       &
                                       ipolp2, qmat2, gqq2, evalue2,           &
                                       berry_energy2, pberryel2, pberryion2
  USE ensemble_dft,             ONLY : tens, e0, z0t, fmat0, fmat0_diag, fmat0_diag_set, gibbsfe, &
                                       tsmear, ef, ismear, degauss => etemp, &
                                       c0diag, becdiag, id_matrix_init, psihpsi, nfroz_occ
  USE cg_module,                ONLY : tcg,  cg_update, c0old
  USE gvecp,                    ONLY : ngm
  USE gvecs,                    ONLY : ngs
  USE gvecb,                    ONLY : ngb
  USE gvecw,                    ONLY : ngw
  USE reciprocal_vectors,       ONLY : gstart, mill_l
  USE ions_base,                ONLY : na, nat, pmass, nax, nsp, rcmax
  USE ions_base,                ONLY : ind_srt, ions_cofmass, ions_kinene, &
                                       ions_temp, ions_thermal_stress, if_pos, extfor
  USE ions_base,                ONLY : ions_vrescal, fricp, greasp, &
                                       iforce, ndfrz, ions_shiftvar, ityp, &
                                       atm, ind_bck, cdm, cdms, ions_cofmsub
  USE cell_base,                ONLY : a1, a2, a3, b1, b2, b3, ainv, frich, &
                                       greash, tpiba2, omega, alat, ibrav,  &
                                       celldm, h, hold, hnew, velh, deth,   &
                                       wmass, press, iforceh, cell_force,   &
                                       thdiag
  USE grid_dimensions,          ONLY : nnrx, nr1, nr2, nr3
  USE smooth_grid_dimensions,   ONLY : nnrsx, nr1s, nr2s, nr3s
  USE smallbox_grid_dimensions, ONLY : nr1b, nr2b, nr3b
  USE local_pseudo,             ONLY : allocate_local_pseudo
  USE io_global,                ONLY : io_global_start, &
                                       stdout, ionode, ionode_id
  USE dener,                    ONLY : detot, denl, dekin6
  USE cdvan,                    ONLY : dbec, drhovan
  USE gvecw,                    ONLY : ggp
  USE constants,                ONLY : pi, k_boltzmann_au, au_ps
  USE io_files,                 ONLY : psfile, pseudo_dir
  USE wave_base,                ONLY : wave_steepest, wave_verlet
  USE wave_base,                ONLY : wave_speed2, frice, grease
  USE control_flags,            ONLY : conv_elec, tconvthrs, gamma_only, do_wf_cmplx !added:giovanni gamma_only, do_wf_cmplx
  USE check_stop,               ONLY : check_stop_now
  USE efcalc,                   ONLY : clear_nbeg, ef_force
  USE ions_base,                ONLY : zv, ions_vel
  USE cp_electronic_mass,       ONLY : emass, emass_cutoff, emass_precond
  USE ions_positions,           ONLY : tau0, taum, taup, taus, tausm, tausp, &
                                       vels, velsm, velsp, ions_hmove,       &
                                       ions_move, fion, fionm
  USE ions_nose,                ONLY : gkbt, kbt, qnp, ndega, nhpcl, nhpdim, &
                                       nhpbeg, nhpend,               &
                                       vnhp, xnhp0, xnhpm, xnhpp,    &
                                       atm2nhp, ions_nosevel, ions_noseupd,  &
                                       tempw, ions_nose_nrg, gkbt2nhp,       &
                                       ekin2nhp, anum2nhp
  USE electrons_nose,           ONLY : qne, ekincw, xnhe0, xnhep, xnhem,  &
                                       vnhe, electrons_nose_nrg,    &
                                       electrons_nose_shiftvar,           &
                                       electrons_nosevel, electrons_noseupd
  USE pres_ai_mod,              ONLY : P_ext, P_in, P_fin, pvar, volclu, &
                                       surfclu, Surf_t, abivol, abisur
  USE wavefunctions_module,     ONLY : c0, cm, phi => cp
  USE wannier_module,           ONLY : allocate_wannier
  USE cp_interfaces,            ONLY : printout_new, move_electrons, rhoofr, new_ns
  USE printout_base,            ONLY : printout_base_open, &
                                       printout_base_close, &
                                       printout_pos, printout_cell, &
                                       printout_stress
  USE cell_nose,                ONLY : xnhh0, xnhhm, xnhhp, vnhh, temph, &
                                       qnh, cell_nosevel, cell_noseupd,  &
                                       cell_nose_nrg, cell_nose_shiftvar
  USE cell_base,                ONLY : cell_kinene, cell_gamma, &
                                       cell_move, cell_hmove
  USE gvecw,                    ONLY : ecutw
  USE gvecp,                    ONLY : ecutp
  USE time_step,                ONLY : delt, tps, dt2, dt2by2, twodelt
  USE cp_interfaces,            ONLY : cp_print_rho, nlfh, print_lambda
  USE cp_main_variables,        ONLY : acc, bec, lambda, lambdam, lambdap, &
                                       ema0bg, sfac, eigr, ei1, ei2, ei3,  &
                                       irb, becdr, taub, eigrb, rhog, rhos, &
                                       rhor, bephi, becp, nfi, descla, iprint_stdout, &
                                       drhor, drhog, nlax, hamilt, collect_zmat
  USE cp_main_variables,        ONLY : vpot
  USE autopilot,                ONLY : event_step, event_index, &
                                       max_event_step, restart_p
  USE cell_base,                ONLY : s_to_r, r_to_s
  USE wannier_subroutines,      ONLY : wannier_startup, wf_closing_options, &
                                       ef_enthalpy
  USE cp_interfaces,            ONLY : readfile, writefile, eigs, strucf, phfacs, eigs_non_ortho
  USE cp_interfaces,            ONLY : empty_cp, ortho, elec_fakekine, print_projwfc
  USE constraints_module,       ONLY : check_constraint, remove_constr_force
  USE metadyn_base,             ONLY : set_target, mean_force
  USE cp_autopilot,             ONLY : pilot
  USE ions_nose,                ONLY : ions_nose_allocate, ions_nose_shiftvar
  USE orthogonalize_base,       ONLY : updatc
  USE control_flags,            ONLY : force_pairing
  USE mp,                       ONLY : mp_bcast
  USE mp_global,                ONLY : root_image, intra_image_comm, np_ortho, me_ortho, ortho_comm, &
                                       me_image
  USE ldaU,                     ONLY : lda_plus_u, vupsi
  USE nksic,                    ONLY : do_orbdep, pink, complexification_index
  USE step_constraint
  USE small_box,                ONLY : ainvb
  USE descriptors,          ONLY: descla_siz_
  USE twin_types
  !
  IMPLICIT NONE
  !
  ! ... input/output variables
  !
  REAL(DP), INTENT(OUT) :: tau_out(3,nat)
  REAL(DP), INTENT(OUT) :: fion_out(3,nat)
  REAL(DP), INTENT(OUT) :: etot_out
  !
  ! ... control variables
  !
  LOGICAL :: tfirst, tlast, tstop, tconv
  LOGICAL :: ttprint, tfile, tstdout
  LOGICAL :: tprint_ham
    !  logical variable used to control printout
  !
  ! ... forces on ions
  !
  REAL(DP) :: maxfion, fion_tot(3)
  !
  ! ... work variables
  !
  REAL(DP) :: tempp, savee, saveh, savep, epot, epre, &
              enow, econs, econt, fccc, ccc, bigr, dt2bye
  REAL(DP) :: ekinc0, ekinp, ekinpr, ekinc
  REAL(DP) :: temps(nat)
  REAL(DP) :: ekinh, temphc, randy
  REAL(DP) :: delta_etot
  REAL(DP) :: ftmp, enb, enbi
  INTEGER  :: is, nacc, ia, j, iter, i, isa, ipos, iat
  INTEGER  :: k, ii, l, m, iss
  REAL(DP) :: hgamma(3,3), temphh(3,3)
  REAL(DP) :: fcell(3,3)
  REAL(DP) :: deltaP, ekincf
  REAL(DP) :: stress_gpa(3,3), thstress(3,3), stress(3,3)
  !
  REAL(DP), ALLOCATABLE :: usrt_tau0(:,:), usrt_taup(:,:), usrt_fion(:,:)
    ! temporary array used to store unsorted positions and forces for
    ! constrained dynamics
  REAL(DP), ALLOCATABLE :: tauw(:,:)  
    ! temporary array used to printout positions
  CHARACTER(LEN=3) :: labelw( nat )
    ! for force_pairing
  INTEGER   :: nspin_sub , i1, i2
  !
  REAL(DP),  ALLOCATABLE :: forceh(:,:)
  REAL(DP),  ALLOCATABLE :: fmat0_repl(:,:)
  COMPLEX(DP),  ALLOCATABLE :: fmat0_repl_c(:,:)
  LOGICAL :: lgam
  COMPLEX(DP), PARAMETER :: c_zero = CMPLX(0.d0,0.d0)

!$$
  LOGICAL :: ttest
!$$
  iter=0
  lgam=gamma_only.and..not.do_wf_cmplx
  !
  !
  dt2bye   = dt2 / emass
  etot_out = 0.D0
  enow     = 1.D9
  !
  tfirst = .TRUE.
  tlast  = .FALSE.
  nacc   = 5
  !
  nspin_sub = nspin
  IF( force_pairing ) nspin_sub = 1
  !
  ! ... Check for restart_p from Autopilot Feature Suite
  !
  IF ( restart_p ) THEN
     !
     ! ... do not add past nfi
     !
     nomore = nomore
     !
  END IF
  !
  IF ( lda_plus_u ) ALLOCATE( forceh( 3, nat ) )

  !
  ! setup data to printout herm and anti-herm
  ! measures of the hamiltonian  H_ij = < i | h_j | j >
  !
  ! matrix elements stored in  hamilt
  !
  tprint_ham = do_orbdep .AND. ( iprsta > 1 ) 
  !
  DO iss=1,size(hamilt)
    call set_twin(hamilt(iss), c_zero)
  END DO
  !
  !======================================================================
  !
  !           basic loop for molecular dynamics starts here
  !
  !======================================================================
  !
  main_loop: DO
     !
     CALL start_clock( 'total_time' )
!$$ For CG calculation, one minimization is enough
     if(tcg) tlast = .true.
!$$
     !
     nfi     = nfi + 1
     tlast   = ( nfi == nomore ) .OR. tlast
     ttprint = ( MOD( nfi, iprint ) == 0 ) .OR. tlast 
     tfile   = ( MOD( nfi, iprint ) == 0 )
     tstdout = ( MOD( nfi, iprint_stdout ) == 0 ) .OR. tlast
     !
     IF ( abivol ) THEN
        IF ( pvar ) THEN
           IF ( nfi .EQ. 1 ) THEN
              deltaP = (P_fin - P_in) / DBLE(nomore)
              P_ext = P_in
           ELSE
              P_ext = P_ext + deltaP
           END IF
        END IF
     END IF
     !
     IF ( ionode .AND. tstdout ) &
        WRITE( stdout, '(/," * Physical Quantities at step:",I6)' ) nfi
     !
     IF ( tnosee ) THEN
        fccc = 1.D0 / ( 1.D0 + 0.5D0 * delt * vnhe )
     ELSE IF ( tsde ) THEN
        fccc = 1.D0 
     ELSE
        fccc = 1.D0 / ( 1.D0 + frice )
     END IF
     !
     ! ... calculation of velocity of nose-hoover variables
     !
     IF ( tnosep ) THEN
        !
        CALL ions_nosevel( vnhp, xnhp0, xnhpm, delt, nhpcl, nhpdim )
        !
     END IF
     !
     IF ( tnosee ) THEN
        !
        CALL electrons_nosevel( vnhe, xnhe0, xnhem, delt )
        !
     END IF
     !
     IF ( tnoseh ) THEN
        !
        CALL cell_nosevel( vnhh, xnhh0, xnhhm, delt )
        !
        velh(:,:) = 2.D0 * ( h(:,:) - hold(:,:) ) / delt - velh(:,:)
        !
     END IF
     ! 
     IF ( (okvan .or. nlcc_any ) .AND. (tfor .OR. thdyn .OR. tfirst) ) THEN
        !
        CALL initbox( tau0, taub, irb, ainv, a1, a2, a3 )
        !
        CALL phbox( taub, eigrb, ainvb )
        !
     END IF
     !
     IF ( tfor .OR. thdyn ) THEN
        !
        CALL phfacs( ei1, ei2, ei3, eigr, mill_l, taus, nr1, nr2, nr3, nat )
        !
        ! ... strucf calculates the structure factor sfac
        !
        CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
        !
     END IF
     !
     IF ( thdyn ) THEN
        !
        CALL formf( tfirst, eself )
        !
     END IF
     !
     ! ... why this call ??? from Paolo Umari
     !
     IF ( tefield .or. tefield2 ) THEN
        !
        CALL calbec( 1, nsp, eigr, c0, bec ) ! ATTENZIONE  
        !
     END IF
     !
     ! Autopilot (Dynamic Rules) Implimentation    
     !
     call pilot(nfi)
     !
     IF ( ( tfor .OR. tfirst ) .AND. tefield ) CALL efield_update( eigr )
     IF ( ( tfor .OR. tfirst ) .AND. tefield2 ) CALL efield_update2( eigr )
     !
     IF ( lda_plus_u ) then
        ! forceh    ! Forces on ions due to Hubbard U 
        forceh=0.0d0
        ! vupsi     ! potentials on electrons due to Hubbard U
        vupsi=(0.0d0,0.0d0)
        ! vpsi_con  ! potentials on electrons due to occupation constraints ...not yet implemented...
        vpsi_con=(0.0d0,0.0d0)
        !
        CALL new_ns(c0,eigr,vkb,vupsi,vpsi_con,forceh)
        if ( mod(nfi,iprint).eq.0 ) call write_ns
     endif
     
     !=======================================================================
     !
     !    electronic degrees of freedom are updated here
     !
     !=======================================================================
     !
     IF( force_pairing ) THEN
          c0(:,iupdwn(2):nbsp)       =     c0(:,1:nupdwn(2))
          cm(:,iupdwn(2):nbsp)       =     cm(:,1:nupdwn(2))
         phi(:,iupdwn(2):nbsp)       =    phi(:,1:nupdwn(2))
!begin_modified:giovanni warning: is nspin=2?
         IF(.not.lambda(1)%iscmplx) THEN
              lambda(2)%rvec(:,:) = lambda(1)%rvec(:,:)
         ELSE
              lambda(2)%cvec(:,:) = lambda(1)%cvec(:,:)
         ENDIF
!end_modified:giovanni
     ENDIF
     !
     ! ... fake electronic kinetic energy
     !
     IF ( .NOT. tcg ) THEN
        !
        ekincf = 0.0d0

        CALL elec_fakekine( ekincf, ema0bg, emass, cm, c0, ngw, nbsp, 1, delt )
        !
     END IF
     !
     !
     CALL move_electrons( nfi, tfirst, tlast, b1, b2, b3, fion, &
                          enthal, enb, enbi, fccc, ccc, dt2bye, stress, &
                          tprint_ham = tprint_ham )
     !
     IF (lda_plus_u) fion = fion + forceh
     !
     IF ( tpre ) THEN
        !
        CALL nlfh( stress, bec, dbec, lambda )
        !
        CALL ions_thermal_stress( stress, pmass, omega, h, vels, nsp, na )
        !
     END IF
     !
     !=======================================================================
     !
     !              verlet algorithm
     !
     !     loop which updates cell parameters and ionic degrees of freedom
     !     hnew=h(t+dt) is obtained from hold=h(t-dt) and h=h(t)
     !     tausp=pos(t+dt) from tausm=pos(t-dt) taus=pos(t) h=h(t)
     !
     !           guessed displacement of ions
     !=======================================================================
     !
     hgamma(:,:) = 0.D0
     !
     IF ( thdyn ) THEN
        !
        CALL cell_force( fcell, ainv, stress, omega, press, wmass )
        !
        CALL cell_move( hnew, h, hold, delt, iforceh, &
                        fcell, frich, tnoseh, vnhh, velh, tsdc )
        !
        velh(:,:) = ( hnew(:,:) - hold(:,:) ) / twodelt
        !
        CALL cell_gamma( hgamma, ainv, h, velh )
        !
     END IF
     !
     !======================================================================
     !
     IF ( tfor ) THEN
        !
        IF ( lwf ) CALL ef_force( fion, na, nsp, zv )
        !
        IF( textfor ) THEN 
           !
           FORALL( ia = 1:nat ) fion(:,ia) = fion(:,ia) + extfor(:,ia)
           !
           fion_tot(:) = SUM( fion(:,:), DIM = 2 ) / DBLE( nat )
           !
           FORALL( ia = 1:nat ) fion(:,ia) = fion(:,ia) - fion_tot(:)
           !
        END IF
        !
        IF ( remove_rigid_rot ) &
           CALL remove_tot_torque( nat, tau0, pmass(ityp(ind_srt(:))), fion )
        !
        IF ( lconstrain ) THEN
           !
           IF ( ionode ) THEN
              !
              ALLOCATE( usrt_tau0( 3, nat ) )
              ALLOCATE( usrt_taup( 3, nat ) )
              ALLOCATE( usrt_fion( 3, nat ) )
              !
              usrt_tau0(:,:) = tau0(:,ind_bck(:))
              usrt_fion(:,:) = fion(:,ind_bck(:))
              !
              IF ( lcoarsegrained ) CALL set_target()
              !
              ! ... we first remove the component of the force along the 
              ! ... constrain gradient (this constitutes the initial guess 
              ! ... for the lagrange multiplier)
              !
              CALL remove_constr_force( nat, usrt_tau0, if_pos, ityp, 1.D0, usrt_fion )
              !
              fion(:,:) = usrt_fion(:,ind_srt(:))
              !
           END IF
           !
           CALL mp_bcast( fion, ionode_id, intra_image_comm )
           !
        END IF
        !
        CALL ions_move( tausp, taus, tausm, iforce, pmass, fion, ainv, &
                        delt, na, nsp, fricp, hgamma, vels, tsdp, tnosep, &
                        fionm, vnhp, velsp, velsm, nhpcl, nhpdim, atm2nhp )
        !
        IF ( lconstrain ) THEN
           !
           ! ... constraints are imposed here
           !
           IF ( ionode ) THEN
              !
              CALL s_to_r( tausp, taup, na, nsp, hnew )
              !
              usrt_taup(:,:) = taup(:,ind_bck(:))
              !
              CALL check_constraint( nat, usrt_taup, usrt_tau0, usrt_fion, &
                                     if_pos, ityp, 1.D0, delt, amu_au )
              !
              taup(:,:) = usrt_taup(:,ind_srt(:))
              fion(:,:) = usrt_fion(:,ind_srt(:))
              !
              ! ... average value of the lagrange multipliers
              !
              IF ( lcoarsegrained ) CALL mean_force( nfi, etot, 1.D0 )
              !
              DEALLOCATE( usrt_tau0, usrt_taup, usrt_fion )
              !
           END IF
           !
           CALL mp_bcast( taup, ionode_id, intra_image_comm )
           CALL mp_bcast( fion, ionode_id, intra_image_comm )
           !
           CALL r_to_s( taup, tausp, na, nsp, ainv )
           !
        END IF
        !
        CALL ions_cofmass( tausp, pmass, na, nsp, cdm )
        !
        IF ( ndfrz == 0 ) &
           CALL ions_cofmsub( tausp, iforce, nat, cdm, cdms )
        !
        CALL s_to_r( tausp, taup, na, nsp, hnew )
        !
     END IF
     !     
     !--------------------------------------------------------------------------
     !              initialization with guessed positions of ions
     !--------------------------------------------------------------------------
     !
     ! ... if thdyn=true g vectors and pseudopotentials are recalculated for 
     ! ... the new cell parameters
     !
     IF ( tfor .OR. thdyn ) THEN
        !
        IF ( thdyn ) THEN
           !
           hold = h
           h    = hnew
           !
           CALL newinit( h )
           !
           CALL newnlinit()
           !
        ELSE
           !
           hold = h
           !
        END IF
        !
        ! ... phfac calculates eigr
        !
        CALL phfacs( ei1, ei2, ei3, eigr, mill_l, tausp, nr1, nr2, nr3, nat ) 
        !
        ! ... prefor calculates vkb
        !
        CALL prefor( eigr, vkb )
        !
     END IF
     !
     !--------------------------------------------------------------------------
     !                    imposing the orthogonality
     !--------------------------------------------------------------------------
     !
     !
     IF ( .NOT. tcg ) THEN
        !
        IF ( tortho ) THEN
           !
!            write(6,*) "calling ortho_cp_twin", phi !added:giovanni:debug
           CALL ortho_cp_twin( eigr(1:ngw,1:nat), cm(1:ngw,1:nbsp), phi(1:ngw,1:nbsp), ngw, &
                                                 lambda(1:nspin), descla(1:descla_siz_ , 1:nspin) &
                                                 , bigr, iter, ccc, bephi, becp, nbsp, nspin, nupdwn, iupdwn )
           !
        ELSE IF(.not.non_ortho) THEN
           !
           CALL gram( vkb, bec, nkb, cm, ngw, nbsp )
           !
           IF ( iprsta > 4 ) CALL dotcsc( eigr, cm, ngw, nbsp, lgam )!added:giovanni lgam
           !
        END IF
        !
        !  correction to displacement of ions
        !
        !IF ( iprsta >= 3 ) CALL print_lambda( lambda, nbsp, 9, 1.D0 )
        !
        IF ( tortho ) THEN
           DO iss = 1, nspin_sub
              i1 = (iss-1)*nlax+1
              i2 = iss*nlax
              IF(.not.lambda(iss)%iscmplx) THEN
		  CALL updatc( ccc, nbsp, lambda(iss)%rvec, SIZE(lambda(iss)%rvec,1),  &
			    phi, SIZE(phi,1),bephi%rvec(:,i1:i2), SIZE(bephi%rvec,1), becp%rvec,  &
			    bec%rvec, cm, nupdwn(iss), iupdwn(iss), descla(:,iss) )
              ELSE
		  CALL updatc( ccc, nbsp, lambda(iss)%cvec, SIZE(lambda(iss)%cvec,1), &
			    phi, SIZE(phi,1), bephi%cvec(:,i1:i2), SIZE(bephi%cvec,1), becp%cvec, &
			     bec%cvec, cm, nupdwn(iss), iupdwn(iss), descla(:,iss) )
              ENDIF
           END DO
        END IF
        !
        IF( force_pairing ) THEN
              c0(:,iupdwn(2):nbsp)       =     c0(:,1:nupdwn(2))
              cm(:,iupdwn(2):nbsp)       =     cm(:,1:nupdwn(2))
             phi(:,iupdwn(2):nbsp)       =    phi(:,1:nupdwn(2)) 
             IF(.not.lambda(1)%iscmplx) THEN
		  lambda(2)%rvec(:,:) = lambda(1)%rvec(:,:)
             ELSE
		  lambda(2)%cvec(:,:) = lambda(1)%cvec(:,:)
             ENDIF
        ENDIF
        !
        CALL calbec( nvb+1, nsp, eigr, cm, bec )
        !
        IF ( tpre ) THEN
           CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, cm, dbec )
        END IF
        !
        IF ( iprsta >= 3 ) CALL dotcsc( eigr, cm, ngw, nbsp, lgam )!added:giovanni lgam
        !
     END IF
     !
     !--------------------------------------------------------------------------
     !                  temperature monitored and controlled
     !--------------------------------------------------------------------------
     !
     ekinp  = 0.D0
     ekinpr = 0.D0
     tempp  = 0.D0
     temps  = 0.D0
     ekinc0 = 0.0d0
     ekinc = 0.0d0
     !
     !
     ! ... ionic kinetic energy and temperature
     !
     IF ( tfor ) THEN
        !
        CALL ions_vel( vels, tausp, tausm, na, nsp, delt )
        !
        CALL ions_kinene( ekinp, vels, na, nsp, hold, pmass )
        !
        CALL ions_temp( tempp, temps, ekinpr, vels, na, nsp, &
                        hold, pmass, ndega, nhpdim, atm2nhp, ekin2nhp )
        !
     END IF
     !
     ! ... fake electronic kinetic energy
     !
     IF ( .NOT. tcg ) THEN
        !
        CALL elec_fakekine( ekinc0, ema0bg, emass, c0, cm, ngw, nbsp, 1, delt )
        !
        ekinc0 = (ekinc0 + ekincf)*0.5d0
        !
        ekinc = ekinc0
        !
     END IF
     !
     ! ... fake cell-parameters kinetic energy
     !
     ekinh = 0.D0
     !
     IF ( thdyn ) THEN
        !
        CALL cell_kinene( ekinh, temphh, velh )
        !
     END IF
     !
     IF ( COUNT( iforceh == 1 ) > 0 ) THEN
        !
        temphc = 2.D0 / k_boltzmann_au * ekinh / DBLE( COUNT( iforceh == 1 ) )
        !
     ELSE
        !
        temphc = 0.D0
        !
     END IF
     !
     ! ... udating nose-hoover friction variables
     !
     IF ( tnosep ) CALL ions_noseupd( xnhpp, xnhp0, xnhpm, delt, qnp, &
                                      ekin2nhp, gkbt2nhp, vnhp, kbt,  &
                                      nhpcl, nhpdim, nhpbeg, nhpend )
     !
     IF ( tnosee ) CALL electrons_noseupd( xnhep, xnhe0, xnhem, &
                                           delt, qne, ekinc, ekincw, vnhe )
     !
     IF ( tnoseh ) CALL cell_noseupd( xnhhp, xnhh0, xnhhm, &
                                      delt, qnh, temphh, temph, vnhh )
     !
     ! ... warning:  thdyn and tcp/tcap are not compatible yet!!!
     !
     IF ( tcp .OR. tcap .AND. tfor .AND. .NOT.thdyn ) THEN
        !
        IF ( tempp > (tempw+tolp) .OR. &
             tempp < (tempw-tolp) .AND. tempp /= 0.D0 ) THEN
           !
           CALL  ions_vrescal( tcap, tempw, tempp, taup, &
                               tau0, taum, na, nsp, fion, iforce, pmass, delt )
           !
        END IF
        !
     END IF
     !
     IF ( MOD( nfi, iprint ) == 0 .OR. tlast ) THEN
        !
!$$        IF ( tortho ) THEN
!$$ In order to calculate the eigenvalues for CG case
        IF ( tortho .or. tcg ) THEN
!$$

!$$ test orthonormality of wavefunctions here
!$$
            !
            IF( force_pairing )  THEN
             IF(.not.lambda(1)%iscmplx) THEN
		  lambda(2)%rvec(:,:) = lambda(1)%rvec(:,:)
             ELSE
		  lambda(2)%cvec(:,:) = lambda(1)%cvec(:,:)
             ENDIF
             IF(.not.lambdap(1)%iscmplx) THEN
		  lambdap(2)%rvec(:,:) = lambdap(1)%rvec(:,:)
             ELSE
		  lambdap(2)%cvec(:,:) = lambdap(1)%cvec(:,:)
             ENDIF
                WRITE( stdout, '("Occupations in CPR:")' )
                WRITE( stdout, '(10F9.6)' ) ( f(i), i = 1, nbspx )  
            END IF
            !
!             write(0,*) lambdap(1)%iscmplx, lambdap(1)%rvec(:,:)
            IF(non_ortho) THEN
               CALL eigs_non_ortho( nfi, lambdap, lambda )
            ELSE
               CALL eigs( nfi, lambdap, lambda )
            ENDIF
!$$
!            do iss=1,nspin
!              if(ionode) write(101,*) 'eigenvalues for spin',iss,'are',(13.6056923 * ei(i,iss),i=1,nupdwn(iss))
!            enddo
!$$
            !
            ! ... Compute empty states
            !
            CALL empty_cp ( nfi, c0, rhos )
            !
        ELSE
            !   
            WRITE( stdout, '(6x,"NOTE: eigenvalues are not computed without ortho")' )
            !
            CALL empty_cp ( nfi, c0, rhos )
            ! 
        ENDIF
        !
     END IF
     !
     IF ( lwf ) CALL ef_enthalpy( enthal, tau0 )
     !
     ! poorman eDFT 
     ! occupations are re-computed according 
     ! to the eigenvalues (recomputed as well)
     !
     IF ( tsmear .AND. tortho .AND. ( MOD(nfi, nfroz_occ)==0 .OR. .NOT. fmat0_diag_set ) ) THEN
         !
         !CALL inner_loop_smear( c0, bec, rhos, psihpsi )
         DO iss=1,nspin
            call copy_twin(psihpsi(iss),lambda(iss))
         END DO
         !
         CALL inner_loop_diag( c0, bec%rvec, psihpsi, z0t, e0 )
         !
         CALL efermi( nelt, nbsp, degauss, 1, f, ef, e0, entropy, ismear, nspin )
         !
         fmat0_diag_set = .TRUE.

         !!
         !! recompute rho
         !!
         !CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, &
         !             becsum, rhor, rhog, rhos, enl, denl, ekin, dekin6 )
         !!
         !IF(nlcc_any) CALL set_cc( irb, eigrb, rhoc )
         
         !!
         !! calculates the SCF potential, the total energy
         !! and the ionic forces
         !!
         !vpot = rhor
         !!
         !CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, &
         !            tlast, ei1, ei2, ei3, irb, eigrb, sfac, &
         !            tau0, fion )         
         
         !
         ! recompute the proper density matrix, once z0t is given
         ! and store its diagonal components
         !
         call calcmt_twin( f, z0t, fmat0, .false. )
         !
         DO iss = 1, nspin
             !
             IF(.not.fmat0(iss)%iscmplx) THEN
	      ALLOCATE( fmat0_repl(nupdwn(iss),nupdwn(iss)) )
	      !
	      CALL collect_zmat( fmat0_repl, fmat0(iss)%rvec(:,:), descla(:,iss) )
	      !
	      DO  i = 1, nupdwn(iss)
		  !
		  fmat0_diag( i+iupdwn(iss)-1 ) = fmat0_repl(i,i)
		  !
	      ENDDO
	      !
	      DEALLOCATE( fmat0_repl )
             ELSE
	      ALLOCATE( fmat0_repl_c(nupdwn(iss),nupdwn(iss)) )
	      !
	      CALL collect_zmat( fmat0_repl_c, fmat0(iss)%cvec(:,:), descla(:,iss) )
	      !
	      DO  i = 1, nupdwn(iss)
		  !
		  fmat0_diag( i+iupdwn(iss)-1 ) = fmat0_repl_c(i,i)
		  !
	      ENDDO
	      !
	      DEALLOCATE( fmat0_repl_c )
             ENDIF
             !
         ENDDO
         !
     ENDIF
     !
     IF ( tens ) THEN
        !
        IF ( MOD( nfi, iprint ) == 0 .OR. tlast ) THEN
            !
            WRITE( stdout, '("Occupations  :")' )
            WRITE( stdout, '(10F9.6)' ) ( f(i), i = 1, nbsp )
            !
        ENDIF
        !
     END IF
     !
     epot = eht + epseu + exc
     !
     IF ( .NOT. tcg ) THEN
        !
        econs = ekinp + ekinh + enthal
        econt = econs + ekinc
        !
     ELSE
        !
        IF ( .NOT. tens ) THEN
           !
           econs = ekinp + etot
           atot  = etot
           econt = econs
           !
        ELSE
           !
           gibbsfe = atot
           econs   = ekinp + atot
           econt   = econs
           !
        END IF
        !
     END IF
     !
     ! ... add energies of thermostats
     !
     IF ( tnosep ) &
        econt = econt + ions_nose_nrg( xnhp0, vnhp, qnp, &
                                       gkbt2nhp, kbt, nhpcl, nhpdim )
     IF ( tnosee ) &
        econt = econt + electrons_nose_nrg( xnhe0, vnhe, qne, ekincw )
     IF ( tnoseh ) &
        econt = econt + cell_nose_nrg( qnh, xnhh0, vnhh, temph, iforceh )
     ! 
     tps = tps + delt * au_ps
     !
     if (abivol) etot = etot - P_ext*volclu
     if (abisur) etot = etot - Surf_t*surfclu
     !
     ! CALL debug_energies()
     !
     IF(do_orbdep.and.(tstdout.or.(MOD(nfi,iprint_stdout)==0))) THEN
        ! 
        !Sort wavefunctions with respect to spread !!added:giovanni
        !
        IF(allocated(wfc_spreads)) THEN
           !
           call spread_sort(ngw, nspin, nbsp, nudx, nupdwn, iupdwn, &
                wfc_spreads, wfc_centers, sort_spreads)
           !
        ENDIF
        !
        IF(do_wf_cmplx) THEN
           !
           ! Compute complexification index
           !
             call compute_complexification_index(ngw, nnrsx, nbsp, nbspx, nspin, ispin, &
              iupdwn, nupdwn, c0, bec, complexification_index)
           !
        ENDIF
        !
     ENDIF

     CALL printout_new( nfi, tfirst, tfile, tstdout, ttprint, tps, hold, stress, &
                        tau0, vels, fion, ekinc, temphc, tempp, temps, etot, &
                        enthal, econs, econt, vnhh, xnhh0, vnhp, xnhp0, atot, &
                        ekin, epot, tprnfor, tpre, hamilt, tprint_ham, lgam )
     !
     if (abivol) etot = etot + P_ext*volclu
     if (abisur) etot = etot + Surf_t*surfclu
     !
     IF( tfor ) THEN
        !
        ! ... new variables for next step
        !
        CALL ions_shiftvar( taup,  tau0, taum  )   !  real positions
        CALL ions_shiftvar( tausp, taus, tausm )   !  scaled positions         
        CALL ions_shiftvar( velsp, vels, velsm )   !  scaled velocities
        !
        IF ( tnosep ) CALL ions_nose_shiftvar( xnhpp, xnhp0, xnhpm )
        IF ( tnosee ) CALL electrons_nose_shiftvar( xnhep, xnhe0, xnhem )
        IF ( tnoseh ) CALL cell_nose_shiftvar( xnhhp, xnhh0, xnhhm )
        !
     END IF
     !
     write(6,*) "tnosep=",tnosep,"thdyn=",thdyn
     !
     IF ( thdyn ) CALL emass_precond( ema0bg, ggp, ngw, tpiba2, emass_cutoff )
     !
     ekincm = ekinc0
     !  
     ! ... cm=c(t+dt) c0=c(t)
     !
     IF( .NOT. tcg ) THEN
        !
        CALL dswap( 2*ngw*nbsp, c0, 1, cm, 1 )
        !
     ELSE
        !
        CALL cg_update( tfirst, nfi, c0 )
        !
        IF ( tfor .AND. .NOT. ( tens ) .AND. &
             ( ( MOD( nfi, isave ) == 0 ) .OR. tlast ) ) THEN
           !
           ! ... in this case optimize c0 and lambda for smooth
           ! ... restart with CP
           !
           IF ( okvan .or. nlcc_any ) THEN
              CALL initbox( tau0, taub, irb, ainv, a1, a2, a3 )
              CALL phbox( taub, eigrb, ainvb )
           END IF
           CALL r_to_s( tau0, taus, na, nsp, ainv )
           CALL phfacs( ei1, ei2, ei3, eigr, mill_l, taus, nr1, nr2, nr3, nat )
           CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
           !
           IF ( thdyn )    CALL formf( tfirst, eself )
           IF ( tefield )  CALL efield_update( eigr )
           IF ( tefield2 ) CALL efield_update2( eigr )
           !
           lambdam = lambda
           !
           CALL move_electrons( nfi, tfirst, tlast, b1, b2, b3, &
                                fion, enthal, enb, enbi, fccc, ccc, dt2bye, stress )
           !
        END IF
        !
     END IF
     !
     ! ... now:  cm=c(t) c0=c(t+dt)
     !
     tfirst = .FALSE.
     !
     ! ... write on file ndw each isave
     !
     IF ( ( MOD( nfi, isave ) == 0 ) .AND. ( nfi < nomore ) ) THEN
        !
        IF ( tcg ) THEN
          !
          CALL writefile( h, hold ,nfi, c0, c0old, taus, tausm,  &
                          vels, velsm, acc, lambda, lambdam, xnhe0, xnhem,     &
                          vnhe, xnhp0, xnhpm, vnhp, nhpcl,nhpdim,ekincm, xnhh0,&
                          xnhhm, vnhh, velh, fion, tps, z0t, f, rhor )
           !
        ELSE
           !
           CALL writefile( h, hold, nfi, c0, cm, taus,  &
                           tausm, vels, velsm, acc,  lambda, lambdam, xnhe0,   &
                           xnhem, vnhe, xnhp0, xnhpm, vnhp, nhpcl, nhpdim, ekincm,&
                           xnhh0, xnhhm, vnhh, velh, fion, tps, z0t, f, rhor )
           !
        END IF
        !
     END IF
     !
     epre = enow
     enow = etot
     !
     frice = frice * grease
     fricp = fricp * greasp
     frich = frich * greash
     !
     !======================================================================
     !
     CALL stop_clock( 'total_time' )
     !
     delta_etot = ABS( epre - enow )
     !
     tstop = check_stop_now() .OR. tlast
!$$
     ttest = check_stop_now()
     tstop = ttest .OR. tlast
!     if(ionode) write(137,*) 'check_stop_now tlast', ttest,tlast
!     if(ionode) write(137,*) 'nfi nomore', nfi,nomore
!$$
     !
     tconv = .FALSE.
     !
     IF ( tconvthrs%active ) THEN
        !
        ! ... electrons
        !
        tconv = ( ekinc < tconvthrs%ekin .AND. delta_etot < tconvthrs%derho )
        !
        IF ( tfor ) THEN
           !
           ! ... ions
           !
           maxfion = MAXVAL( ABS( fion(:,1:nat) ) )
           !
           tconv = tconv .AND. ( maxfion < tconvthrs%force )
           !
        END IF
        !
     END IF
!$$
!     if(ionode) write(1111,'("etot=",F20.13," delta_etot=",F20.13," conv=",F20.13)') enow,delta_etot,tconvthrs%derho
     !if(ionode) write(1111,*) 'tconvthrs%active,tconv,ttest,tlast,tstop,tfor'
     !if(ionode) write(1111,*) tconvthrs%active,tconv,ttest,tlast,tstop,tfor
!$$$$
!     if(nfi.lt.10000) tconv=.false.
!$$$$
!$$
     !
     ! ... in the case cp-wf the check on convergence is done starting
     ! ... from the second step 
     !
     IF ( lwf .AND. tfirst ) tconv = .FALSE.
     !
     IF ( tconv ) THEN
        !
        tlast = .TRUE.
        !
        IF ( ionode ) THEN
           !
           WRITE( stdout, &
                & "(/,3X,'MAIN:',10X,'EKINC   (thr)', &
                & 10X,'DETOT   (thr)',7X,'MAXFORCE   (thr)')" )
           WRITE( stdout, "(3X,'MAIN: ',3(D14.6,1X,D8.1))" ) &
               ekinc, tconvthrs%ekin, delta_etot,                  &
               tconvthrs%derho, 0.D0, tconvthrs%force
           WRITE( stdout, &
                  "(3X,'MAIN: convergence achieved for system relaxation')" )
           !
        END IF
        !
     END IF
     !
     IF ( lwf ) &
        CALL wf_closing_options( nfi, c0, cm, bec, eigr, eigrb, taub, &
                                 irb, ibrav, b1, b2, b3, taus, tausm, vels, &
                                 velsm, acc, lambda, lambdam, xnhe0, xnhem, &
                                 vnhe, xnhp0, xnhpm, vnhp, nhpcl, nhpdim, &
                                 ekincm, xnhh0, xnhhm, vnhh, velh, ecutp, &
                                 ecutw, delt, celldm, fion, tps, z0t, f, rhor )
!$$
!     if(ionode) write(137,*) 'tconv', tconv
!$$
     !
!$$
!$$     exit main_loop
!$$
     IF ( tstop ) EXIT main_loop
     !
  END DO main_loop
  !
  !===================== end of main loop of molecular dynamics ===============
  !
  ! ... Here copy relevant physical quantities into the output arrays/variables
  !
  etot_out = etot
  !
  isa = 0
  !
  DO is = 1, nsp
     !
     DO ia = 1, na(is)
        !
        isa = isa + 1
        !
        ipos = ind_srt( isa )
        !
        tau_out(:,ipos) = tau0(:,isa)
        !
        fion_out(:,ipos) = fion(:,isa)
        !
     END DO
     !
  END DO
  !
  IF ( lneb ) fion_out(:,1:nat) = fion(:,1:nat) * DBLE( if_pos(:,1:nat) )
  !
  conv_elec = .TRUE.
  !
  IF ( tcg ) cm = c0old
  !
  CALL writefile( h, hold, nfi, c0, cm, taus, tausm, &
                  vels, velsm, acc, lambda, lambdam, xnhe0, xnhem, vnhe,    &
                  xnhp0, xnhpm, vnhp, nhpcl,nhpdim,ekincm, xnhh0, xnhhm,    &
                  vnhh, velh, fion, tps, z0t, f, rhor )
  !
  IF( tprojwfc ) CALL print_projwfc( c0, lambda, eigr, vkb )
  !
  IF( iprsta > 1 ) CALL print_lambda( lambda, nbsp, nbsp, 1.D0 )
  !
  IF (lda_plus_u) DEALLOCATE( forceh )
  !
  RETURN
  !
END SUBROUTINE cprmain
!
!----------------------------------------------------------------------------
SUBROUTINE terminate_run()
  !----------------------------------------------------------------------------
  !
  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout, ionode
  USE cp_main_variables, ONLY : acc
  USE cg_module,         ONLY : tcg, print_clock_tcg
  USE mp,                ONLY : mp_report
  USE control_flags,     ONLY : use_task_groups
  USE nksic,             ONLY : do_orbdep
  USE hfmod,             ONLY : do_hf
  !
  IMPLICIT NONE
  !
  ! ...  print statistics
  !
  CALL printacc()
  !
  CALL print_clock( 'initialize' )
  CALL print_clock( 'total_time' )
  CALL print_clock( 'main_loop' )
  CALL print_clock( 'formf' )
  CALL print_clock( 'rhoofr' )
  CALL print_clock( 'vofrho' )
  CALL print_clock( 'vofrho1' )
  CALL print_clock( 'vofrho2' )
  CALL print_clock( 'vofrho3' )
  CALL print_clock( 'vofrho4' )
  CALL print_clock( 'vofrho5' )
  CALL print_clock( 'dforce' )
  CALL print_clock( 'dforce1' )
  CALL print_clock( 'dforce2' )
  CALL print_clock( 'calphi' )
  CALL print_clock( 'ortho' )
  CALL print_clock( 'ortho_iter' )
  CALL print_clock( 'rsg' )
  CALL print_clock( 'rhoset' )
  CALL print_clock( 'updatc' )
  CALL print_clock( 'gram' )
  CALL print_clock( 'newd' )
  CALL print_clock( 'calbec' )
  CALL print_clock( 'prefor' )
  CALL print_clock( 'strucf' )
  CALL print_clock( 'nlfl' )
  CALL print_clock( 'nlfq' )
  CALL print_clock( 'set_cc' )
  CALL print_clock( 'rhov' )
  CALL print_clock( 'nlsm1' )
  CALL print_clock( 'nlsm2' )
  CALL print_clock( 'forcecc' )
  CALL print_clock( 'fft' )
  CALL print_clock( 'ffts' )
  CALL print_clock( 'fftw' )
  CALL print_clock( 'fftb' )
  CALL print_clock( 'cft3s' )
  CALL print_clock( 'fft_scatter' )
  IF ( ionode ) WRITE(stdout,"()")
  !
  IF ( do_orbdep ) THEN
     !
     CALL print_clock( 'nk_drv' )
        CALL print_clock( 'nk_orbrho' )
        CALL print_clock( 'nk_rhoref' )
        CALL print_clock( 'nk_eforce' )
     CALL print_clock( 'nk_corr' )
        CALL print_clock( 'nk_corr_h' )
        CALL print_clock( 'nk_corr_vxc' )
        CALL print_clock( 'nk_corr_fxc' )
     !
     IF ( ionode ) WRITE(stdout,"()")
     CALL print_clock( 'nk_rot_emin' )
     CALL print_clock( 'nk_innerloop' )
     CALL print_clock( 'nk_rot_cm' )
     CALL print_clock( 'nk_get_vsicah' )
     CALL print_clock( 'nk_getOmattot' )
     CALL print_clock( 'nk_getOmat1' )
     CALL print_clock( 'nk_rotwfn' )
     IF ( ionode ) WRITE(stdout,"()")
     !
  ENDIF
  !
  IF ( do_hf ) THEN
     !
     CALL print_clock( 'hf_potential' )
     CALL print_clock( 'hf_corr' )
     IF ( ionode ) WRITE(stdout,"()")
     !
  ENDIF
  !
  IF ( tcg ) THEN
     !
     CALL print_clock( 'outer_loop' )
     !CALL print_clock( 'inner_loop' )
     call print_clock_tcg()
     !
  ENDIF
  !
  IF( use_task_groups ) THEN
     !
     CALL print_clock( 'ALLTOALL' )
     !
  END IF
  !
  CALL mp_report()
  !
END SUBROUTINE terminate_run
