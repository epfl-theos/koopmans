!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
SUBROUTINE from_scratch( )
    !
    USE kinds,                ONLY : DP
    USE control_flags,        ONLY : tranp, trane, iprsta, tpre, tcarpar,  &
                                     tzeroc, tzerop, tzeroe, tfor, thdyn, &
                                     lwf, tprnfor, non_ortho, tortho, amprp, ampre,  &
                                     tsde, ortho_eps, ortho_max, program_name, &
                                     force_pairing, use_task_groups, gamma_only, &
                                     do_wf_cmplx, iprint_manifold_overlap !added:giovanni gamma_only, do_wf_cmplx
    USE ions_positions,       ONLY : taus, tau0, tausm, vels, fion, fionm, atoms0
    USE ions_base,            ONLY : na, nsp, randpos, zv, ions_vel, pmass
    USE ions_base,            ONLY : taui, cdmi, nat, iforce
    USE ions_nose,            ONLY : xnhp0, xnhpm, vnhp
    USE cell_base,            ONLY : ainv, h, s_to_r, ibrav, omega, press, &
                                     hold, r_to_s, deth, wmass, iforceh,   &
                                     cell_force, boxdimensions, velh, a1,  &
                                     a2, a3, b1, b2, b3
    USE cell_nose,            ONLY : xnhh0, xnhhm, vnhh
    USE electrons_nose,       ONLY : xnhe0, xnhem, vnhe
    use electrons_base,       ONLY : nbsp, f, nspin, nupdwn, iupdwn, nbspx
    USE electrons_module,     ONLY : occn_info
    USE energies,             ONLY : entropy, eself, enl, ekin, enthal, etot, ekincm
    USE energies,             ONLY : dft_energy_type, debug_energies
    USE dener,                ONLY : denl, denl6, dekin6, detot
    USE uspp,                 ONLY : vkb, becsum, deeq, nkb, okvan
    USE io_global,            ONLY : stdout, ionode, meta_ionode
    USE core,                 ONLY : nlcc_any, rhoc
    USE gvecw,                ONLY : ngw
    USE gvecs,                ONLY : ngs
    USE gvecp,                ONLY : ngm
    USE reciprocal_vectors,   ONLY : gstart, mill_l, gx, ig_l2g!added:giovanni ig_l2g
    USE cvan,                 ONLY : nvb
    USE cp_electronic_mass,   ONLY : emass
    USE efield_module,        ONLY : tefield, efield_berry_setup, berry_energy, &
                                     tefield2, efield_berry_setup2, berry_energy2
    USE cg_module,            ONLY : tcg
    USE ensemble_dft,         ONLY : tens, compute_entropy
    USE cp_interfaces,        ONLY : runcp_uspp, runcp_uspp_force_pairing, &
                                     strucf, phfacs, nlfh, nlfl
    USE cp_interfaces,        ONLY : rhoofr, ortho, wave_rand_init, elec_fakekine, &
                                     wave_sine_init, wave_atom_init
    USE cp_interfaces,        ONLY : vofrhos, compute_stress
    USE cp_interfaces,        ONLY : printout, print_lambda
    USE printout_base,        ONLY : printout_pos
    USE orthogonalize_base,   ONLY : updatc, calphi 
    USE atoms_type_module,    ONLY : atoms_type
    USE wave_base,            ONLY : wave_steepest
    USE wavefunctions_module, ONLY : c0, cm, phi => cp, cdual
    USE grid_dimensions,      ONLY : nr1, nr2, nr3
    USE time_step,            ONLY : delt
    USE cp_main_variables,    ONLY : setval_lambda, descla, bephi, becp, becdr, nfi, &
                                     sfac, eigr, ei1, ei2, ei3, bec, taub, irb, eigrb, &
                                     lambda, lambdam, lambdap, ema0bg, rhog, rhor, rhos, &
                                     vpot, ht0, edft, nlax, becdual
    USE mp_global,            ONLY : np_ortho, me_ortho, ortho_comm
    USE small_box,            ONLY : ainvb
    USE cdvan,                ONLY : dbec
    USE nksic,                ONLY : do_spinsym
    USE mp_global,      ONLY : mpime !added:giovanni:debug
    USE control_flags,        ONLY : tatomicwfc, trane
    USE descriptors,          ONLY: descla_siz_

    !
    IMPLICIT NONE
    !
    REAL(DP),    ALLOCATABLE :: emadt2(:), emaver(:)
    REAL(DP)                 :: verl1, verl2
    REAL(DP)                 :: bigr, dum
    INTEGER                  :: i, j, iter, iss, ierr, nspin_wfc
    LOGICAL                  :: tlast = .FALSE.
    REAL(DP)                 :: gam(1,1,1)
    REAL(DP)                 :: fcell(3,3), ccc, enb, enbi, fccc
    LOGICAL                  :: ttforce
    LOGICAL                  :: tstress
    LOGICAL, PARAMETER       :: ttprint = .TRUE.
    REAL(DP)                 :: ei_unp  
    REAL(DP)                 :: dt2bye
    INTEGER                  :: n_spin_start 
    LOGICAL                  :: tfirst = .TRUE.
    REAL(DP)                 :: stress(3,3)
    INTEGER                  :: i1, i2 
    LOGICAL                  :: lgam !added:giovanni
    COMPLEX(DP), parameter :: c_zero = CMPLX(0.d0,0.d0) !added:giovanni
    !
    ! ... Subroutine body
    !
    lgam=gamma_only.and..not.do_wf_cmplx !added:giovanni
    nfi = 0
    !
    ttforce = tfor  .or. tprnfor
    tstress = thdyn .or. tpre
    !
    stress = 0.0d0
    !
    IF( tsde ) THEN
       fccc = 1.0d0
    ELSE
       fccc = 0.5d0
    END IF
    !
    dt2bye = delt * delt / emass
    !
    IF( ANY( tranp( 1:nsp ) ) ) THEN
       !
       CALL invmat( 3, h, ainv, deth )
       !
       CALL randpos( taus, na, nsp, tranp, amprp, ainv, iforce )
       !
       CALL s_to_r( taus, tau0, na, nsp, h )
       !
    END IF
    !
    CALL phfacs( ei1, ei2, ei3, eigr, mill_l, atoms0%taus, nr1, nr2, nr3, atoms0%nat )
    !
    CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
    !     
    IF ( okvan .OR. nlcc_any ) THEN
       CALL initbox ( tau0, taub, irb, ainv, a1, a2, a3 )
       CALL phbox( taub, eigrb, ainvb )
    END IF
    !
    !     wfc initialization with random numbers
    !     
    IF ( ionode ) &
       WRITE( stdout, fmt = '(//,3X, "Wave Initialization: random initial wave-functions" )' )
    !
    IF ( .NOT. do_spinsym .OR. nspin == 1 ) then
       !
       IF(tatomicwfc) THEN
          !
          call wave_atom_init( cm, nbsp, 1 )
          !
       ELSE
          !
          CALL wave_rand_init( cm, nbsp, 1 )!modified:giovanni
          !
       ENDIF
! !begin_added:giovanni:debug
!         DO i=1,size(cm(:,1))
!             write(201+mpime,'(5((F18.12)5x))') gx(1:3,i), cm(i,1) 
!         ENDDO
! !end_added:giovanni:debug
        !
    ELSE 
        !
        IF ( nupdwn(1) < nupdwn(2) ) CALL errore('from_scratch','unexpec nupdwn(1) < nupdwn(2)',10)
        !
        CALL wave_rand_init( cm, nupdwn(1) , 1 )
        !
        DO i = 1, MIN(nupdwn(1),nupdwn(2))
            !
            j=i+iupdwn(2)-1
            cm(:,j)=cm(:,i)
            !
        ENDDO
        !
        IF( meta_ionode ) write(stdout, "(24x, 'spin symmetry applied to init wave')" )
        !
    ENDIF

    !
    ! ... prefor calculates vkb (used by gram)
    !
    CALL prefor( eigr, vkb )
    !
    nspin_wfc = nspin
    IF( force_pairing ) nspin_wfc = 1
    !

    DO iss = 1, nspin_wfc
       !
       CALL gram( vkb, bec, nkb, cm(1,iupdwn(iss)), ngw, nupdwn(iss) )
       !
    END DO
! begin_added:giovanni:debug
!         DO i=1,size(cm(:,1))
!             write(201+mpime,'(5((F18.12)5x))') gx(1:3,i), cm(i,1) 
!         ENDDO
!        stop
! end_added:giovanni:debug
    !
    IF( force_pairing ) cm(:,iupdwn(2):iupdwn(2)+nupdwn(2)-1) = cm(:,1:nupdwn(2))
    !
    if( iprsta >= 3 ) CALL dotcsc( eigr, cm, ngw, nbsp, lgam )!added:giovanni lgam
    !
    ! ... initialize bands
    !
    CALL occn_info( f )
    !
    atoms0%for  = 0.D0
    atoms0%vels = 0.D0
    hold = h
    velh = 0.0d0
    fion = 0.0d0
    tausm = taus
    !
    ! ... compute local form factors
    !
    CALL formf( tfirst, eself )
    !
    edft%eself = eself

    IF( tefield ) THEN
      CALL efield_berry_setup( eigr, tau0 )
    END IF
    IF( tefield2 ) THEN
      CALL efield_berry_setup2( eigr, tau0 )
    END IF
    !
    IF( .NOT. tcg ) THEN
       !
! ! begin_added:giovanni:debug
!         DO i=1,size(cm(:,1))
!              write(202,'(5((F18.12)5x))') gx(1:3,i), eigr( i,1) 
!         ENDDO
! ! end_added:giovanni:debug
       CALL calbec ( 1, nsp, eigr, cm, bec )
       !
       if ( tstress ) CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, cm, dbec ) !warning:giovanni still to be modified
       !
       IF(non_ortho) THEN
          call compute_duals(cm,cdual,nbsp,1)
          call calbec(1,nsp,eigr,cdual,becdual)
       ENDIF
       !
       CALL rhoofr ( nfi, cm(:,:), irb, eigrb, bec, becsum, rhor, rhog, rhos, enl, denl, ekin, dekin6 )
       !
       edft%enl  = enl
       edft%ekin = ekin
       !
    END IF

    !
    !     put core charge (if present) in rhoc(r)
    !
    if ( nlcc_any ) CALL set_cc( irb, eigrb, rhoc ) !warning:giovanni not yet implemented
    !
    IF( program_name == 'CP90' ) THEN

       IF( .NOT. tcg ) THEN
   
         IF( tens ) THEN
           CALL compute_entropy( entropy, f(1), nspin )
           entropy = entropy * nbsp
         END IF
         !
         vpot = rhor
         !
         CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
        &  ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion )

         IF( tefield ) THEN !warning:giovanni modified but not checked
           CALL berry_energy( enb, enbi, bec, cm(:,:), fion )
           etot = etot + enb + enbi
         END IF
         IF( tefield2 ) THEN !warning:giovanni modified but not checked
           CALL berry_energy2( enb, enbi, bec, cm(:,:), fion )
           etot = etot + enb + enbi
         END IF

         CALL compute_stress( stress, detot, h, omega )

         if(iprsta.gt.2) &
             CALL printout_pos( stdout, fion, nat, head = ' fion ' )

         CALL newd( vpot, irb, eigrb, becsum, fion )
         !
         IF( force_pairing ) THEN
            !
            CALL runcp_uspp_force_pairing( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec%rvec, cm, &
        &                 c0, ei_unp, fromscra = .TRUE. ) !warning:giovanni not yet modified
            !
            IF(.not.lambda(2)%iscmplx) THEN
		CALL setval_lambda( lambda(2)%rvec, nupdwn(1), nupdwn(1), 0.d0, descla(:,1) )
            ELSE
                CALL setval_lambda( lambda(2)%cvec, nupdwn(1), nupdwn(1), c_zero, descla(:,1) )
            ENDIF
            !
         ELSE
            !
! !begin_added:giovanni:debug ---- RUNCP WAVEFUNCTION
!             write(6,*) "wavefunctions before runcp"
!             write(6,*) "wavefunction_0", c0(1,2), c0(2,2), c0(3,2)
!             write(6,*) "wavefunction_m", cm(1,2), cm(2,2), cm(3,2)
! !end_added:giovanni:debug ---- RUNCP WAVEFUNCTION
            CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, cm, c0, fromscra = .TRUE. )
! !begin_added:giovanni:debug ---- RUNCP WAVEFUNCTION
!             write(6,*) "wavefunctions after runcp"
!             write(6,*) "check_allocated_5", associated(bec%rvec), ubound(bec%rvec)!added:giovanni:debug
!             write(6,*) "wavefunction_0", c0(1,2), c0(2,2), c0(3,2)
!             write(6,*) "wavefunction_m", cm(1,2), cm(2,2), cm(3,2)
! !end_added:giovanni:debug ---- RUNCP WAVEFUNCTION
            !
         ENDIF
         !
         !     nlfq needs deeq bec
         !
         if( ttforce ) CALL nlfq( cm, eigr, bec, becdr, fion, lgam )
         !
         !     calphi calculates phi
         !     the electron mass rises with g**2
         !
         CALL calphi( cm, ngw, bec, nkb, vkb, phi, nbsp, lgam, ema0bg)
         !
         IF( force_pairing ) &
         &   phi( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) =    phi( :, 1:nupdwn(2))


         IF( tortho ) THEN
            CALL ortho_cp_twin( eigr(1:ngw,1:nat), c0(1:ngw,1:nbsp), & 
                                phi(1:ngw,1:nbsp), ngw, lambda, descla(1:descla_siz_ , 1:nspin), &
                                bigr, iter, ccc, bephi, becp, nbsp, nspin, nupdwn, iupdwn )
         ELSE IF(.not.non_ortho) THEN
            !
            CALL gram( vkb, bec, nkb, c0, ngw, nbsp )
            !
         ENDIF
         !
         !
         if ( ttforce ) CALL nlfl_twin( bec, becdr, lambda, fion, lgam ) !warning:giovanni may not work

         if ( iprsta >= 3 ) CALL print_lambda( lambda, nbsp, 9, ccc )

         if ( tstress ) CALL nlfh( stress, bec, dbec, lambda )
         !
         IF ( tortho ) THEN
            DO iss = 1, nspin_wfc
               i1 = (iss-1)*nlax+1
               i2 = iss*nlax
               IF(.not.bec%iscmplx) THEN
		    CALL updatc( ccc, nbsp, lambda(iss)%rvec, SIZE(lambda(iss)%rvec,1), phi, SIZE(phi,1), &
				  bephi%rvec(:,i1:i2), SIZE(bephi%rvec,1), becp%rvec, bec%rvec, c0, nupdwn(iss), iupdwn(iss), &
				  descla(:,iss) )
               ELSE
		    CALL updatc( ccc, nbsp, lambda(iss)%cvec, SIZE(lambda(iss)%cvec,1), phi, SIZE(phi,1), &
				  bephi%cvec(:,i1:i2), SIZE(bephi%cvec,1), becp%cvec, bec%cvec, c0, nupdwn(iss), iupdwn(iss), &
				  descla(:,iss) )
               ENDIF
            END DO
         END IF
         !
         IF( force_pairing ) THEN
            !
            c0 ( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) =  c0( :, 1:nupdwn(2))
            phi( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) = phi( :, 1:nupdwn(2))

            IF(.not.lambda(1)%iscmplx) THEN
		lambda(2)%rvec(:,:) = lambda(1)%rvec(:,:)
            ELSE
		lambda(2)%cvec(:,:) = lambda(1)%cvec(:,:)
            ENDIF
            !
         ENDIF
         !
         CALL calbec ( nvb+1, nsp, eigr, c0, bec)

         if ( tstress ) CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, cm, dbec )

         if ( iprsta >= 3 ) CALL dotcsc( eigr, c0, ngw, nbsp , lgam)
         !
         xnhp0 = 0.0d0
         xnhpm = 0.0d0
         vnhp  = 0.0d0
         fionm = 0.0d0
         !
         CALL ions_vel( vels, taus, tausm, na, nsp, delt )
         !
         xnhh0(:,:) = 0.0d0
         xnhhm(:,:) = 0.0d0
         vnhh (:,:) = 0.0d0
         velh (:,:) = ( h(:,:) - hold(:,:) ) / delt
         !
         CALL elec_fakekine( ekincm, ema0bg, emass, c0, cm, ngw, nbsp, 1, delt )

         xnhe0 = 0.0d0
         xnhem = 0.0d0
         vnhe  = 0.0d0

         DO iss=1,nspin
            IF(.not.lambda(iss)%iscmplx) THEN
                lambdam(iss)%rvec = lambda(iss)%rvec
            ELSE
                lambdam(iss)%cvec = lambda(iss)%cvec
            ENDIF
         ENDDO
         !
       ELSE 
          !
          c0 = cm
          !
       END IF

    ELSE
       !
       IF( ttforce ) call nlfq( cm, eigr, bec, becdr, atoms0%for, lgam )
       !
       CALL vofrhos( ttprint, ttforce, tstress, rhor, rhog, atoms0, &
                  vpot, bec%rvec, cm, f, eigr, ei1, ei2, ei3, sfac, ht0, edft ) 
       !
       IF( iprsta > 1 ) CALL debug_energies( edft )
       !
       IF ( .NOT. force_pairing ) THEN
          !
          CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, vpot, bec, cm, c0, fromscra = .TRUE. )
          !
       ELSE
          !
          c0 = cm
          !
       END IF
       !
       IF ( tortho .AND. ( .NOT. force_pairing ) ) THEN
          !
          ccc = fccc * dt2bye
          !
          CALL ortho( cm, c0, lambda, descla, ccc, nupdwn, iupdwn, nspin )
          !
       ELSE IF(.not. non_ortho) THEN
          !
          DO iss = 1, nspin_wfc
            !
            CALL gram( vkb, bec, nkb, c0(1,iupdwn(iss)), SIZE(c0,1), nupdwn( iss ) )
            !
          END DO
          !
          IF( force_pairing ) c0(1,iupdwn(2):iupdwn(2)+nupdwn(2)-1) = c0(1,1:nupdwn(2))
          !
       END IF
       !
       !
    END IF

    IF( iprsta > 1 ) THEN
       !
       !  Printout values at step 0, useful for debugging
       !
       CALL printout( nfi, atoms0, 0.0d0, 0.0d0, ttprint, ht0, edft )
       !
    END IF
    !
    RETURN
    !
END SUBROUTINE from_scratch
