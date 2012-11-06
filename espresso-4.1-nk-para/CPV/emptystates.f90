!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"


!-----------------------------------------------------------------------
   SUBROUTINE empty_cp_twin_x ( nfi, c0, v, tcg )
!-----------------------------------------------------------------------
!
! Performs the minimization on the empty state subspace keeping the
! occupied manyfold fixed. A proper orthogonalization of the two
! manyfolds is performed. 
!
      USE kinds,                ONLY : DP
      USE control_flags,        ONLY : iprsta, tsde, program_name, gamma_only, do_wf_cmplx, &
                                           tortho
      USE io_global,            ONLY : ionode, stdout
      USE cp_main_variables,    ONLY : eigr, ema0bg, collect_lambda, &
                                       rhor, rhog, rhos, eigr, eigrb, irb, bec
      USE descriptors,          ONLY : descla_siz_ , descla_init, nlax_, lambda_node_
      USE cell_base,            ONLY : omega
      USE uspp,                 ONLY : vkb, nkb, okvan
      USE uspp_param,           ONLY : nhm
      USE grid_dimensions,      ONLY : nnrx
      USE electrons_base,       ONLY : nbsp, nbspx, ispin, nspin, f, nudx, iupdwn, nupdwn
      USE electrons_module,     ONLY : iupdwn_emp, nupdwn_emp, n_emp, ei_emp,  &
                                       max_emp, ethr_emp, etot_emp, eodd_emp
      USE ions_base,            ONLY : nat, nsp
      USE gvecw,                ONLY : ngw
      USE orthogonalize_base,   ONLY : calphi, updatc
      USE reciprocal_vectors,   ONLY : gzero, gstart
      USE wave_base,            ONLY : wave_steepest, wave_verlet, frice
      USE cvan,                 ONLY : nvb
      USE cp_electronic_mass,   ONLY : emass
      USE time_step,            ONLY : delt
      USE check_stop,           ONLY : check_stop_now
      USE cp_interfaces,        ONLY : writeempty, readempty, gram_empty, ortho, &
                                       wave_rand_init, wave_atom_init, elec_fakekine, &
                                       crot, dforce, nlsm1, grabec, &
                                       bec_csv
      USE mp,                   ONLY : mp_comm_split, mp_comm_free, mp_sum
      USE mp_global,            ONLY : intra_image_comm, me_image
      USE nksic,                ONLY : do_orbdep, do_pz, do_wxd, vsicpsi, wtot, sizwtot
      USE nksic,                ONLY : do_spinsym, pink_emp, allocate_nksic_empty
      USE hfmod,                ONLY : do_hf, vxxpsi
      USE twin_types !added:giovanni
      USE control_flags,        ONLY : tatomicwfc, trane
      USE electrons_module,     ONLY : wfc_centers_emp, wfc_spreads_emp, icompute_spread
      !
      IMPLICIT NONE
      !
      INTEGER,    INTENT(IN) :: nfi
      COMPLEX(DP)            :: c0(:,:)
      REAL(DP)               :: v(:,:)
      logical, optional, intent(IN) :: tcg
      !
      INTEGER  :: i, iss, j, in, in_emp, iter, iter_ortho
      INTEGER  :: n_occs, n_emps, n_empx, nudx_emp, issw, n
      INTEGER  :: nlax_emp, nlam_emp
      LOGICAL  :: exst, do_wxd_, tcg_
      !
      REAL(DP) :: fccc, ccc, csv, dt2bye, bigr
      REAL(DP) :: verl1, verl2, verl3
      REAL(DP) :: dek, ekinc, ekinc_old, detothf
      !
      REAL(DP),    ALLOCATABLE :: emadt2(:)
      REAL(DP),    ALLOCATABLE :: emaver(:)
      COMPLEX(DP), ALLOCATABLE :: c2(:), c3(:)
      COMPLEX(DP), ALLOCATABLE :: c0_emp(:,:), cm_emp(:,:), phi_emp(:,:)
!       REAL(DP),    ALLOCATABLE :: bec_emp(:,:)
      type(twin_matrix) :: bec_emp !modified:giovanni
      REAL(DP),    ALLOCATABLE :: becsum_emp(:,:,:)
      type(twin_matrix) :: bephi_emp! !modified:giovanni
      type(twin_matrix) :: becp_emp !modified:giovanni
      type(twin_matrix) :: bec_occ !(:,:) !modified:giovanni
      type(twin_matrix), dimension(:),  ALLOCATABLE :: lambda_emp !(:,:,:) !, 
      REAL(DP),    ALLOCATABLE :: f_emp(:)
      REAL(DP),    ALLOCATABLE :: lambda_rep(:,:)
      COMPLEX(DP),    ALLOCATABLE :: lambda_rep_c(:,:)
      INTEGER,     ALLOCATABLE :: ispin_emp(:)
      REAL(DP),    ALLOCATABLE :: fsic_emp(:)
      REAL(DP),    ALLOCATABLE :: vsic_emp(:,:)
      REAL(DP),    ALLOCATABLE :: wxd_emp(:,:)
      REAL(DP),    ALLOCATABLE :: deeq_sic_emp(:,:,:,:)
      COMPLEX(DP), ALLOCATABLE :: vxxpsi_emp(:,:)
      REAL(DP),    ALLOCATABLE :: exx_emp(:)
!       REAL(DP),    ALLOCATABLE :: pink_emp(:)
      !
      INTEGER, SAVE :: np_emp(2), me_emp(2), emp_comm, color
      INTEGER, SAVE :: desc_emp( descla_siz_ , 2 )
      LOGICAL, SAVE :: first = .true.
      LOGICAL :: lgam !added:giovanni
      LOGICAL :: done_extra !added:giovanni
      COMPLEX(DP), PARAMETER :: c_zero=CMPLX(0.d0,0.d0)

      lgam=gamma_only.and..not.do_wf_cmplx
      !
      if(present(tcg)) THEN
         !
         tcg_=tcg
         !
      ELSE
         !
         tcg_=.false.
         !
      ENDIF
      !
      ! ...  quick exit if empty states have not to be computed
      !
      IF( n_emp < 1 ) RETURN
      !
      !  Here set the group of processors for empty states
      !
      IF( .NOT. first ) THEN
         CALL mp_comm_free( emp_comm )
      END IF 
      !
      np_emp = 1
      IF( me_image < np_emp(1) * np_emp(2) ) THEN
         color = 1
      ELSE
         color = 0
      END IF
      CALL mp_comm_split( intra_image_comm, color, me_image, emp_comm )
      
      if( me_image <  np_emp(1) * np_emp(2) ) then
          me_emp(1) = me_image / np_emp(1)
          me_emp(2) = MOD( me_image, np_emp(1) ) 
      else
          me_emp(1) = me_image
          me_emp(2) = me_image
      endif
      !
      first = .FALSE.
      !
      !  Done with the group
      !
      ! n_occs    == nbsp
      ! n_emps    => nbsp   (corresponds to)
      ! n_empx    => nbspx  
      ! nudx_emp  => nudx
      !
      n_occs = nupdwn( 1 )
      IF( nspin == 2 ) n_occs = n_occs + nupdwn( 2 )
      !
      n_emps = nupdwn_emp( 1 )
      IF( nspin == 2 ) n_emps = n_emps + nupdwn_emp( 2 )
      !
      nudx_emp = nupdwn_emp( 1 ) !+ MOD( nupdwn_emp( 1 ), 2)
      IF( nspin == 2 ) nudx_emp = MAX( nudx_emp, nupdwn_emp( 2 ) )
      !
      n_empx = nupdwn_emp( 1 )
      IF( nspin == 2 ) n_empx = n_empx + nupdwn_emp( 2 )
      n_empx = n_empx + MOD( n_empx, 2)
      !
      DO iss = 1, nspin
         CALL descla_init( desc_emp( :, iss ), nupdwn_emp( iss ), nudx_emp, np_emp, me_emp, emp_comm, color )
      END DO
      !
      nlax_emp = MAXVAL( desc_emp( nlax_, 1:nspin ) )
      nlam_emp = 1
      IF ( ANY( desc_emp( lambda_node_, : ) > 0 ) ) nlam_emp = nlax_emp
      !
      ALLOCATE( c0_emp( ngw, n_empx ) )
      ALLOCATE( cm_emp( ngw, n_empx ) )
      ALLOCATE( phi_emp( ngw, n_empx ) )

      call init_twin(bec_emp, lgam)
      call allocate_twin(bec_emp, nkb, n_emps, lgam)
      call init_twin(becp_emp, lgam)
      call allocate_twin(becp_emp, nkb, n_emps, lgam)
      call init_twin(bec_occ, lgam)
      call allocate_twin(bec_occ, nkb, n_occs, lgam)
      call init_twin(bephi_emp, lgam)
      call allocate_twin(bephi_emp, nkb, n_emps, lgam)

      ALLOCATE(lambda_emp(nspin))
! 
      DO iss=1,nspin
	call init_twin(lambda_emp(iss), lgam)
	call allocate_twin(lambda_emp(iss), nlam_emp, nlam_emp, lgam)
      ENDDO

      ALLOCATE( f_emp( n_empx ) )
      ALLOCATE( ispin_emp( n_empx ) )
      !
      c0_emp     = 0.0
      cm_emp     = 0.0
      !
      phi_emp    = 0.0d0

      call set_twin(bec_emp,c_zero)
      call set_twin(bec_occ,c_zero)
      call set_twin(bephi_emp,c_zero)
      call set_twin(becp_emp,c_zero)
      DO iss=1,nspin
         call set_twin(lambda_emp(iss),c_zero)
      ENDDO

      f_emp      = 2.0d0 / DBLE(nspin)
      !
      ispin_emp = 0
      ispin_emp( 1:nupdwn_emp( 1 ) ) = 1
      IF( nspin == 2 ) ispin_emp( iupdwn_emp(2) : ) = 2
      !
      !
      IF ( do_orbdep ) THEN
          !
          ALLOCATE( fsic_emp( n_empx ) )
!           n_empx_odd=n_empx
          ALLOCATE( vsic_emp(nnrx, n_empx) )
          ALLOCATE( wxd_emp (nnrx, 2) )
          ALLOCATE( deeq_sic_emp (nhm,nhm,nat,n_empx) )
          ALLOCATE( becsum_emp(nhm*(nhm+1)/2,nat,nspin))
          CALL allocate_nksic_empty(n_empx)
          !
          fsic_emp = 0.0d0
          vsic_emp = 0.0d0
          wxd_emp  = 0.0d0
          pink_emp = 0.0d0
          ! 
      ELSE
          !
          ALLOCATE( fsic_emp( 1 ) )
!           n_empx_odd=1
          ALLOCATE( vsic_emp(1, 1) )
          ALLOCATE( wxd_emp (1, 2) )
          ALLOCATE( deeq_sic_emp (1,1,1,1) )
          ALLOCATE( becsum_emp( 1,1, 1 ) )
          call allocate_nksic_empty(1)
          !
      ENDIF
      !
      IF ( do_hf ) THEN
          !
          ALLOCATE( fsic_emp(n_empx ) )
          ALLOCATE( vxxpsi_emp( ngw, n_empx) )
          ALLOCATE( exx_emp( n_empx ) )
          !
          fsic_emp   = 0.0d0
          vxxpsi_emp = 0.0d0
          exx_emp    = 0.0d0
          !
      ENDIF
      !
      ! init wfcs
      !
      exst = readempty( c0_emp, n_empx )
      !
      CALL prefor( eigr, vkb )
      !
      DO iss = 1, nspin
          issw = iupdwn( iss )
          CALL nlsm1 ( nupdwn( iss ), 1, nvb, eigr, c0( 1:1, issw:issw ), bec_occ, iupdwn(iss), lgam ) !warning:giovanni:definition
      ENDDO
      !
      IF( .NOT. exst ) THEN
         !
         ! ...  initial random states orthogonal to filled ones
         !
         IF ( .NOT. do_spinsym .OR. nspin == 1 ) THEN
             !
             IF(tatomicwfc) THEN
                !
                call wave_atom_init( c0_emp, n_emps, 1 )
                !
             ELSE
                !
                CALL wave_rand_init( c0_emp, n_emps, 1 )
                !
             ENDIF
             !
         ELSE
             !  
             IF ( nupdwn_emp(1) < nupdwn_emp(2) ) CALL errore('empty_cp','unexpec emp nupdwn(1) < nupdwn(2)',10)
                !
                IF(tatomicwfc) THEN
                   !
                   call wave_atom_init( c0_emp, nupdwn_emp(1), 1 )
                   !
                ELSE
                   !
                   CALL wave_rand_init( c0_emp, nupdwn_emp(1) , 1 )
                   !
                ENDIF
                !
             DO i = 1, MIN(nupdwn_emp(1),nupdwn_emp(2))
                 !
                 j=i+iupdwn_emp(2)-1
                 c0_emp(:,j)=c0_emp(:,i)
                 !
             ENDDO
             !
             IF( ionode ) write(stdout, "(24x, 'spin symmetry applied to init wave')" )
             !
         ENDIF
         !
         IF ( gzero ) THEN
            c0_emp( 1, : ) = (0.0d0, 0.0d0)
         END IF
         !
         CALL nlsm1 ( n_emps, 1, nvb, eigr, c0_emp, bec_emp, 1, lgam )
         !
         DO iss = 1, nspin
            !
            in     = iupdwn(iss)
            in_emp = iupdwn_emp(iss)
            !
            issw   = iupdwn( iss )
            !
            if(nupdwn(iss)>0.and.nupdwn_emp(iss)>0) THEN
               !
               CALL gram_empty( .false. , eigr, vkb, bec_emp, bec_occ, nkb, &
                c0_emp( :, in_emp: ), c0( :, issw: ), ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, in )
               !
            ENDIF
            !
         END DO
         !
      END IF

       
      !
      !WRITE(stdout,*) '------- filled filled -------'
      !DO i = 1, nupdwn( 1 )
      !   DO j = 1, nupdwn( 1 )
      !      CALL dotcsv( csv, eigr, c0(1,i), c0(1,j), ngw )
      !      WRITE(stdout,"(2i5,f15.9)") i, j, csv
      !   END DO
      !END DO
      !WRITE(stdout,*) '------- empty filled -------'
      !DO i = 1, n_emps
      !  DO j = 1, nupdwn( 1 )
      !      CALL dotcsv( csv, eigr, c0_emp(1,i), c0(1,j), ngw )
      !      WRITE(stdout,"(2i5,f15.9)") i, j, csv
      !   END DO
      !END DO
      !IF( nspin == 2 ) THEN
      !DO i = iupdwn_emp(2), iupdwn_emp(2) + nupdwn_emp( 2 ) - 1
      !  DO j = iupdwn(2), iupdwn(2) + nupdwn( 2 ) - 1
      !      CALL dotcsv( csv, eigr, c0_emp(1,i), c0(1,j), ngw )
      !      WRITE(stdout,"(2i5,f15.9)") i, j, csv
      !   END DO
      !END DO
      !END IF
      !WRITE(stdout,*) '------- empty empty -------'
      !DO i = 1, n_emps
      !   DO j = 1, n_emps
      !      CALL dotcsv( csv, eigr, c0_emp(1,i), c0_emp(1,j), ngw )
      !      WRITE(stdout,"(2i5,f15.9)") i, j, csv
      !   END DO
      !END DO
      

      CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
      !
      ! ...  set verlet variables
      !
      IF( tsde ) THEN
          fccc = 1.0d0
      ELSE   
          fccc = 1.0d0 / ( 1.0d0 + frice )
      END IF
      !
      verl1 = 2.0d0 * fccc
      verl2 = 1.0d0 - verl1
      verl3 = 1.0d0 * fccc
      !
      ALLOCATE( c2( ngw ) )
      ALLOCATE( c3( ngw ) )
      ALLOCATE( emadt2( ngw ) )
      ALLOCATE( emaver( ngw ) )

      dt2bye = delt * delt / emass

      ccc    = fccc   * dt2bye
      emadt2 = dt2bye * ema0bg
      emaver = emadt2 * verl3
 
      cm_emp = c0_emp

      ekinc_old = 0.0
      ekinc     = 0.0

      !
      ! init xd potential
      !
      ! we need to use wtot from previous calls with occupied states
      ! we save here wtot in wxd_emp
      !
      IF ( do_orbdep ) THEN
          !
          wxd_emp(:,:) = 0.0_DP
          !
          IF ( do_wxd .AND. .NOT. do_pz ) THEN
              !
              wxd_emp(:,:) = wtot(:,:)
              !
          ENDIF
      ENDIF
      !
      !
      IF( ionode ) THEN
          WRITE( stdout, "(/,3X,'Empty states minimization starting ', &
                        & /,3x,'nfi         dekinc         ekinc' )")
      ENDIF
      
      IF(tcg_) THEN ! compute empty states with conjugate gradient
      
         call runcg_uspp_emp(c0_emp, cm_emp, bec_emp, f_emp, fsic_emp, n_empx,&
                          n_emps, ispin_emp, iupdwn_emp, nupdwn_emp, phi_emp, lambda_emp, &
                          max_emp, wxd_emp, vsic_emp, pink_emp, nnrx, becsum_emp, &
                          deeq_sic_emp, nudx_emp, eodd_emp, etot_emp, v, &
                          nfi, .true., .true., eigr, bec, irb, eigrb, &
                          rhor, rhog, rhos, ema0bg, desc_emp)
      
      ELSE ! compute empty states with damped dynamics
      !
         done_extra=.false.
         !
         ITERATIONS: DO iter = 1, max_emp

            IF ( do_orbdep ) THEN 
                !
                ! In the nksic case, we do not need to compute wxd here, 
                ! because no contribution comes from empty states.
                !
                ! Instead, wxd from all occupied states is already computed
                ! by the previous calls to nksic_potentials, and stored wxe_emp
                !
                fsic_emp(:) = 0.0
                !
                ! the two lines below were removed by Giovanni, passing do_wxd as input to nksic_potential
                !do_wxd_ = do_wxd
                !do_wxd  = .FALSE.
                !
                IF(done_extra.or.iter==max_emp) THEN
                   !
                   icompute_spread=.true.
                   !
                ENDIF
!                 !
!                 write(6,*) "checkbounds", ubound(wfc_centers_emp), ubound(wfc_spreads_emp), nudx_emp, nspin
                !
                call nksic_potential( n_emps, n_empx, c0_emp, fsic_emp, &
                                      bec_emp, becsum_emp, deeq_sic_emp, &
                                      ispin_emp, iupdwn_emp, nupdwn_emp, rhor, rhog, &
                                      wtot, sizwtot, vsic_emp, .false., pink_emp, nudx_emp, &
                                      wfc_centers_emp, wfc_spreads_emp, &
                                      icompute_spread, .false.)
!                                                       write(6,*) "checkbounds", ubound(wfc_centers_emp), ubound(wfc_spreads_emp), nudx_emp, nspin

                ! line below removed by Giovanni, introduced do_wxd=.false. into call to nksic_potential
                !do_wxd = do_wxd_
                !
                DO i = 1, n_emps
                    vsic_emp(:,i) = vsic_emp(:,i) + wxd_emp(:, ispin_emp(i) )
                ENDDO
                !
                !
            ENDIF
            !
            ! HF contribution
            !
            IF ( do_hf ) THEN
                !
                vxxpsi_emp = 0.0d0
                !
                CALL hf_potential( nbsp,   nbspx,  c0,     f, ispin, iupdwn, nupdwn, &
                                   n_emps, n_empx, c0_emp, fsic_emp, ispin_emp, &
                                   iupdwn_emp, nupdwn_emp, rhor, rhog, vxxpsi_emp, exx_emp )
                !
            ENDIF

            !
            ! std terms
            !
            DO i = 1, n_emps, 2
                !
                CALL dforce( i, bec_emp, vkb, c0_emp, c2, c3, v, SIZE(v,1), ispin_emp, f_emp, n_emps, nspin )

                !
                ! ODD terms
                !
                IF ( do_orbdep ) THEN
                    !
                    CALL nksic_eforce( i, n_emps, n_empx, vsic_emp, deeq_sic_emp, bec_emp, ngw, &
                                       c0_emp(:,i), c0_emp(:,i+1), vsicpsi, lgam )
                    !
                    c2(:) = c2(:) - vsicpsi(:,1) * f_emp(i)
                    c3(:) = c3(:) - vsicpsi(:,2) * f_emp(i+1)
                    !
                ENDIF
                !
                ! HF terms
                !
                IF ( do_hf ) THEN
                    !
                    c2(:) = c2(:) - vxxpsi_emp(:,i)   * f_emp(i)
                    c3(:) = c3(:) - vxxpsi_emp(:,i+1) * f_emp(i+1)
                    !
                ENDIF


                IF( tsde ) THEN
                   CALL wave_steepest( cm_emp(:, i  ), c0_emp(:, i  ), emaver, c2 )
                   CALL wave_steepest( cm_emp(:, i+1), c0_emp(:, i+1), emaver, c3 )
                ELSE
                  CALL wave_verlet( cm_emp(:, i  ), c0_emp(:, i  ), verl1, verl2, emaver, c2 )
                  CALL wave_verlet( cm_emp(:, i+1), c0_emp(:, i+1), verl1, verl2, emaver, c3 )
                END IF
                !
                IF(lgam) THEN
                   IF ( gstart == 2) THEN
                       cm_emp(1,  i)=CMPLX(DBLE(cm_emp(1,  i)),0.d0)
                       cm_emp(1,i+1)=CMPLX(DBLE(cm_emp(1,i+1)),0.d0)
                   ENDIF
                ENDIF
                !
            ENDDO

            DO iss = 1, nspin
                !
                in     = iupdwn(iss)
                in_emp = iupdwn_emp(iss)
                !
                issw   = iupdwn( iss )
                !
                IF(nupdwn(iss)>0.and.nupdwn_emp(iss)>0) THEN
                   !
                   CALL gram_empty( .true. , eigr, vkb, bec_emp, bec_occ, nkb, &
                            cm_emp( :, in_emp: ), c0( :, issw: ), ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, in )
                   !
                ENDIF
                !
            ENDDO

            ! ... calphi calculates phi
            ! ... the electron mass rises with g**2
            !
            CALL calphi( c0_emp, ngw, bec_emp, nkb, vkb, phi_emp, n_emps, lgam, ema0bg )
            !
            !
            IF( tortho ) THEN
               !
!                write(6,*) "n_emps", n_emps
               CALL ortho_cp_twin( eigr(1:ngw,1:nat), cm_emp(1:ngw,1:n_emps), phi_emp(1:ngw,1:n_emps), ngw, lambda_emp(1:nspin), desc_emp(1:descla_siz_,1:nspin), &
                           bigr, iter_ortho, ccc, bephi_emp, becp_emp, n_emps, nspin, &
                           nupdwn_emp, iupdwn_emp )
               !
            ELSE
               !
               CALL gram( vkb, bec_emp, nkb, cm_emp, ngw, n_emps )
               !
            ENDIF
            !
            DO iss = 1, nspin
                !
                IF(.not.lambda_emp(iss)%iscmplx) THEN
                   CALL updatc( ccc, n_emps, lambda_emp(iss)%rvec(:,:), SIZE(lambda_emp(iss)%rvec,1), &
                                 phi_emp, ngw, bephi_emp%rvec, nkb, becp_emp%rvec, bec_emp%rvec, &
                                 cm_emp, nupdwn_emp(iss), iupdwn_emp(iss), desc_emp(:,iss) )
                ELSE
                   CALL updatc( ccc, n_emps, lambda_emp(iss)%cvec(:,:), SIZE(lambda_emp(iss)%cvec,1), & 
                                 phi_emp, ngw, bephi_emp%cvec(:,:), nkb, becp_emp%cvec(:,:), bec_emp%cvec(:,:), &
                                 cm_emp, nupdwn_emp(iss), iupdwn_emp(iss), desc_emp(:,iss) )
                ENDIF
                !
            ENDDO
            !
            CALL nlsm1 ( n_emps, 1, nsp, eigr, cm_emp, bec_emp, 1, lgam )
            !
            CALL elec_fakekine( ekinc, ema0bg, emass, c0_emp, cm_emp, ngw, n_emps, 1, delt )
            !
            CALL dswap( 2*SIZE( c0_emp ), c0_emp, 1, cm_emp, 1 )

            dek = ekinc - ekinc_old
            IF( ionode ) WRITE( stdout,113) ITER, dek, ekinc
         
            ! ...   check for exit
            !
            IF ( check_stop_now() ) THEN
                EXIT ITERATIONS
            ENDIF
         
            ! ...   check for convergence
            !     
            IF( ( ekinc <  ethr_emp ) .AND. ( iter > 3 ) ) THEN
                IF(done_extra) THEN
                   IF( ionode ) WRITE( stdout,112)
                   EXIT ITERATIONS
                ELSE
                   done_extra=.true.
                ENDIF
            ENDIF
         
            ekinc_old = ekinc


         ENDDO ITERATIONS
         
      ENDIF !if clause to choose between cg and damped dynamics
      !
      IF ( ionode ) WRITE( stdout, "()")

      ! ...  Compute eigenvalues and bring wave functions on Kohn-Sham orbitals

      IF(.not.lambda_emp(1)%iscmplx) THEN
	  ALLOCATE( lambda_rep( nudx_emp, nudx_emp ) )
      ELSE
	  ALLOCATE( lambda_rep_c( nudx_emp, nudx_emp ) )
      ENDIF
      !
      DO iss = 1, nspin
          !
          i = iupdwn_emp(iss)
          n = nupdwn_emp(iss)
          IF(.not.lambda_emp(iss)%iscmplx) THEN
	      CALL collect_lambda( lambda_rep, lambda_emp(iss)%rvec(:,:), desc_emp( :, iss ) )
	      CALL crot( cm_emp, c0_emp, ngw, n, i, i, lambda_rep, nudx_emp, ei_emp(:,iss) )
          ELSE
	      CALL collect_lambda( lambda_rep_c, lambda_emp(iss)%cvec(:,:), desc_emp( :, iss ) )
	      CALL crot( cm_emp, c0_emp, ngw, n, i, i, lambda_rep_c, nudx_emp, ei_emp(:,iss) )
          ENDIF
          ei_emp( 1:n, iss ) = ei_emp( 1:n, iss ) / f_emp( i : i + n - 1 )
          !
      ENDDO
      !
      IF(.not.lambda_emp(1)%iscmplx) THEN
	  DEALLOCATE( lambda_rep)
      ELSE
	  DEALLOCATE( lambda_rep_c)
      ENDIF


      ! ...   Save emptystates to disk

      CALL writeempty( cm_emp, n_empx )
      ! 
      DEALLOCATE( ispin_emp )
      DEALLOCATE( f_emp )
      DEALLOCATE( emadt2 )
      DEALLOCATE( emaver )
      DEALLOCATE( c2 )
      DEALLOCATE( c3 )
      DEALLOCATE( c0_emp )
      DEALLOCATE( cm_emp )
      DEALLOCATE( phi_emp )
      CALL deallocate_twin(bec_emp)
      CALL deallocate_twin(bec_occ)
      CALL deallocate_twin(bephi_emp)
      CALL deallocate_twin(becp_emp)
      DO iss=1,nspin
         CALL deallocate_twin(lambda_emp(iss))
      ENDDO

      !
      IF ( do_orbdep ) THEN
          DEALLOCATE( fsic_emp ) 
          DEALLOCATE( vsic_emp ) 
          DEALLOCATE( wxd_emp ) 
          DEALLOCATE( deeq_sic_emp )
          DEALLOCATE( becsum_emp )
!           DEALLOCATE( pink_emp )
      ENDIF
      !
      IF ( do_hf ) THEN
          DEALLOCATE( fsic_emp ) 
          DEALLOCATE( vxxpsi_emp )
          DEALLOCATE( exx_emp )
      ENDIF
  
112   FORMAT(/,3X,'Empty states: convergence achieved')
113   FORMAT(I5,2X,2D14.6)

      RETURN
   END SUBROUTINE empty_cp_twin_x


!-------------------------------------------------------------------------
   SUBROUTINE gram_empty_real_x  &
      ( tortho, eigr, betae, bec_emp, bec_occ, nkbx, c_emp, c_occ, ngwx, n_emp, n_occ )
!-----------------------------------------------------------------------
!     gram-schmidt orthogonalization of the empty states ( c_emp ) 
!     c_emp are orthogonalized among themself and to the occupied states c_occ
!
      USE uspp,           ONLY : nkb, nkbus
      USE cvan,           ONLY : nvb
      USE gvecw,          ONLY : ngw
      USE kinds,          ONLY : DP
      USE mp,             ONLY : mp_sum
      USE mp_global,      ONLY : intra_image_comm
      USE ions_base,      ONLY : nat
      USE cp_interfaces, ONLY : bec_csv, smooth_csv, grabec
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: nkbx, ngwx, n_emp, n_occ
      COMPLEX(DP)   :: eigr(ngwx,nat)
      REAL(DP)      :: bec_emp( nkbx, n_emp )
      REAL(DP)      :: bec_occ( nkbx, n_occ )
      COMPLEX(DP)   :: betae( ngwx, nkb )
      COMPLEX(DP)   :: c_emp( ngwx, n_emp )
      COMPLEX(DP)   :: c_occ( ngwx, n_occ )
      LOGICAL, INTENT(IN) :: tortho
!
      REAL(DP) :: anorm, cscnorm
      REAL(DP), ALLOCATABLE :: csc_emp( : )
      REAL(DP), ALLOCATABLE :: csc_occ( : )
      INTEGER :: i, k, inl
      EXTERNAL cscnorm
!
      ALLOCATE( csc_emp( n_emp ) )
      ALLOCATE( csc_occ( n_occ ) )
      !
      ! orthogonalize empty states to the occupied one and among each other
      !
      DO i = 1, n_emp
         !
         csc_emp = 0.0d0
         csc_occ = 0.0d0
         !
         ! compute scalar product csc_occ(k) = <c_emp(i)|c_occ(k)>
         !
         CALL smooth_csv( c_emp(1:,i), c_occ(1:,1:), ngwx, csc_occ, n_occ )
         !
         IF( .NOT. tortho ) THEN
           !
           ! compute scalar product csc_emp(k) = <c_emp(i)|c_emp(k)>
           !
           CALL smooth_csv( c_emp(1:,i), c_emp(1:,1:), ngwx, csc_emp, i-1 )
           !
           CALL mp_sum( csc_emp, intra_image_comm )
           !
         END IF
         !
         CALL mp_sum( csc_occ, intra_image_comm )
         !
         IF( nvb > 1 ) THEN
            !
            CALL grabec( bec_emp(1:,i), nkbx, betae, c_emp(1:,i:), ngwx )
            !
            CALL mp_sum( bec_emp(1:nkbus,i), intra_image_comm )
            !
            CALL bec_csv( bec_emp(1:,i), bec_occ, nkbx, csc_occ, n_occ )
            !
            IF( .NOT. tortho ) THEN
               CALL bec_csv( bec_emp(1:,i), bec_emp, nkbx, csc_emp, i-1 )
            END IF
            !
            DO k = 1, n_occ
               DO inl = 1, nkbx
                  bec_emp( inl, i ) = bec_emp( inl, i ) - csc_occ(k) * bec_occ( inl, k )
               END DO
            END DO
            !
            IF( .NOT. tortho ) THEN
               DO k = 1, i-1
                  DO inl = 1, nkbx
                     bec_emp( inl, i ) = bec_emp( inl, i ) - csc_emp(k) * bec_emp( inl, k )
                  END DO
               END DO
            END IF
            !
         END IF
         !
         ! calculate orthogonalized c_emp(i) : |c_emp(i)> = |c_emp(i)> - SUM_k    csv(k)|c_occ(k)>
         !                          c_emp(i) : |c_emp(i)> = |c_emp(i)> - SUM_k<i  csv(k)|c_emp(k)>
         !
         DO k = 1, n_occ
            CALL DAXPY( 2*ngw, -csc_occ(k), c_occ(1,k), 1, c_emp(1,i), 1 )
         END DO
         IF( .NOT. tortho ) THEN
            DO k = 1, i - 1
               CALL DAXPY( 2*ngw, -csc_emp(k), c_emp(1,k), 1, c_emp(1,i), 1 )
            END DO
         END IF
         !
         !
         IF( .NOT. tortho ) THEN
            anorm = cscnorm( bec_emp, nkbx, c_emp, ngwx, i, n_emp )
            !
            CALL DSCAL( 2*ngw, 1.0d0/anorm, c_emp(1,i), 1 )
            !
            IF( nvb > 1 ) THEN
               CALL DSCAL( nkbx, 1.0d0/anorm, bec_emp(1,i), 1 )
            END IF
         END IF
         !
      END DO

      DEALLOCATE( csc_emp )
      DEALLOCATE( csc_occ )
      !
      RETURN
   END SUBROUTINE gram_empty_real_x

!-------------------------------------------------------------------------
   SUBROUTINE gram_empty_twin_x  & 
      ( tortho, eigr, betae, bec_emp, bec_occ, nkbx, c_emp, c_occ, ngwx, n_emp, n_occ, l_emp, l_occ )
!-----------------------------------------------------------------------
!
!     gram-schmidt orthogonalization of the empty states ( c_emp ) 
!     c_emp are orthogonalized among themself and to the occupied states c_occ
!
      USE uspp,           ONLY : nkb, nkbus
      USE cvan,           ONLY : nvb
      USE gvecw,          ONLY : ngw
      USE kinds,          ONLY : DP
      USE mp,             ONLY : mp_sum
      USE mp_global,      ONLY : intra_image_comm
      USE ions_base,      ONLY : nat
      USE control_flags,  ONLY : gamma_only, do_wf_cmplx !added:giovanni
      USE cp_interfaces,  ONLY : grabec, smooth_csv,bec_csv
      USE twin_types
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: nkbx, ngwx, n_emp, n_occ, l_emp, l_occ
      COMPLEX(DP)   :: eigr(ngwx,nat)
      type(twin_matrix)      :: bec_emp !( nkbx, n_emp )
      type(twin_matrix)      :: bec_occ !( nkbx, n_occ )
      COMPLEX(DP)   :: betae( ngwx, nkb )
      COMPLEX(DP)   :: c_emp( ngwx, n_emp )
      COMPLEX(DP)   :: c_occ( ngwx, n_occ )
      LOGICAL, INTENT(IN) :: tortho
!
      REAL(DP) :: anorm, cscnorm
      type(twin_matrix) :: csc_emp, csc_occ
!       COMPLEX(DP), ALLOCATABLE :: csc_emp!( : )
!       COMPLEX(DP), ALLOCATABLE :: csc_occ !( : )
      INTEGER :: i, k, inl
      LOGICAL :: lgam
      EXTERNAL cscnorm

      lgam=gamma_only.and..not.do_wf_cmplx

      !Quick return if there are either no filled or no empty states with this spin (no need to orthogonalize them)
      IF(n_emp.le.0.or.n_occ.le.0) THEN
         !
         return
         !
      ENDIF
      
      call init_twin(csc_emp, lgam)
      call allocate_twin(csc_emp, n_emp, 1, lgam)

      call init_twin(csc_occ, lgam)
      call allocate_twin(csc_occ, n_occ, 1, lgam)

      !
      ! orthogonalize empty states to the occupied one and among each other
      !
      DO i = 1, n_emp
         !
         call set_twin(csc_emp, CMPLX(0.d0,0.d0))
         call set_twin(csc_occ, CMPLX(0.d0,0.d0))

         !
         ! compute scalar product csc_occ(k) = <c_emp(i)|c_occ(k)> .. is it <k,i>? Yes! watch out!
         !
         CALL smooth_csv( c_emp(1:ngwx,i), c_occ(1:,1:), ngwx, csc_occ, n_occ )
         !
         IF( .NOT. tortho ) THEN
           !
           ! compute scalar product csc_emp(k) = <c_emp(i)|c_emp(k)>
           !
           CALL smooth_csv( c_emp(1:,i), c_emp(1:,1:), ngwx, csc_emp, i-1 )
           !
           IF(.not.csc_emp%iscmplx) THEN
	      CALL mp_sum( csc_emp%rvec, intra_image_comm )
           ELSE
	      CALL mp_sum( csc_emp%cvec, intra_image_comm )
           ENDIF
           !
         END IF
         !
	IF(.not.csc_occ%iscmplx) THEN
	  CALL mp_sum( csc_occ%rvec, intra_image_comm )
	ELSE
	  CALL mp_sum( csc_occ%cvec, intra_image_comm )
	ENDIF
         !
         IF( nvb > 1 ) THEN
            !
            CALL grabec( bec_emp, nkbx, betae, c_emp(1:,i:), ngwx,i )
            !
            IF(.not. bec_emp%iscmplx) THEN
		CALL mp_sum( bec_emp%rvec(1:nkbus,i), intra_image_comm )
            ELSE
		CALL mp_sum( bec_emp%cvec(1:nkbus,i), intra_image_comm )
            ENDIF
            !
            CALL bec_csv( bec_emp, bec_occ, nkbx, csc_occ, n_occ, i  )
            !
            IF( .NOT. tortho ) THEN
               CALL bec_csv( bec_emp, bec_emp, nkbx, csc_emp, i-1, i )
            END IF
            !
            IF(.not.bec_emp%iscmplx) THEN
	      DO k = 1, n_occ
		DO inl = 1, nkbx
		    bec_emp%rvec( inl, i ) = bec_emp%rvec( inl, i ) - csc_occ%rvec(k,1) * bec_occ%rvec( inl, k )
		END DO
	      END DO
            ELSE
	      DO k = 1, n_occ
		DO inl = 1, nkbx
		    bec_emp%cvec( inl, i ) = bec_emp%cvec( inl, i ) - CONJG(csc_occ%cvec(k,1)) * bec_occ%cvec( inl, k )
		END DO
	      END DO
            ENDIF
            !
            IF( .NOT. tortho ) THEN
               IF(.not.bec_emp%iscmplx) THEN
		  DO k = 1, i-1
		      DO inl = 1, nkbx
			bec_emp%rvec( inl, i ) = bec_emp%rvec( inl, i ) - csc_emp%rvec(k,1) * bec_emp%rvec( inl, k )
		      END DO
		  END DO
               ELSE
		  DO k = 1, i-1
		      DO inl = 1, nkbx
			bec_emp%cvec( inl, i ) = bec_emp%cvec( inl, i ) - CONJG(csc_emp%cvec(k,1)) * bec_emp%cvec( inl, k )
		      END DO
		  END DO
               ENDIF
            END IF
            !
         END IF
         !
         ! calculate orthogonalized c_emp(i) : |c_emp(i)> = |c_emp(i)> - SUM_k    csv(k)|c_occ(k)>
         !                          c_emp(i) : |c_emp(i)> = |c_emp(i)> - SUM_k<i  csv(k)|c_emp(k)>
         !
         IF(.not.csc_occ%iscmplx) THEN
	    DO k = 1, n_occ
		CALL DAXPY( 2*ngw, -csc_occ%rvec(k,1), c_occ(1,k), 1, c_emp(1,i), 1 )!warning:giovanni tochange
	    END DO
	    IF( .NOT. tortho ) THEN
		DO k = 1, i - 1
		  CALL DAXPY( 2*ngw, -csc_emp%rvec(k,1), c_emp(1,k), 1, c_emp(1,i), 1 )!warning:giovanni tochange
		END DO
	    END IF
         ELSE
	    DO k = 1, n_occ
		CALL ZAXPY( ngw, -csc_occ%cvec(k,1), c_occ(1,k), 1, c_emp(1,i), 1 )
	    END DO
	    IF( .NOT. tortho ) THEN
		DO k = 1, i - 1
		  CALL ZAXPY( ngw, -csc_emp%cvec(k,1), c_emp(1,k), 1, c_emp(1,i), 1 )
		END DO
	    END IF
         ENDIF
         !
         !
         IF( .NOT. tortho ) THEN
	    anorm = cscnorm( bec_emp, nkbx, c_emp, ngwx, i, n_emp, lgam)
            IF(.not.bec_emp%iscmplx) THEN
	      !
	      CALL DSCAL( 2*ngw, 1.0d0/anorm, c_emp(1,i), 1 )
	      !
	      IF( nvb > 1 ) THEN
		CALL DSCAL( nkbx, 1.0d0/anorm, bec_emp%rvec(1,i), 1 )
	      END IF
            ELSE
! 	      anorm = cscnorm( bec_emp, nkbx, c_emp, ngwx, i, n_emp, lgam)
	      !
	      CALL ZSCAL( ngw, CMPLX(1.0d0/anorm,0.d0), c_emp(1,i), 1 )
	      !
	      IF( nvb > 1 ) THEN
		CALL ZSCAL( nkbx, CMPLX(1.0d0/anorm,0.d0), bec_emp%cvec(1,i), 1 )
	      END IF
            ENDIF
         END IF
         !
      END DO

      call deallocate_twin( csc_emp )
      call deallocate_twin( csc_occ )
      !
      RETURN
   END SUBROUTINE gram_empty_twin_x

!-----------------------------------------------------------------------
   LOGICAL FUNCTION readempty_x( c_emp, ne )
!-----------------------------------------------------------------------

        ! ...   This subroutine reads empty states from unit emptyunit

        USE kinds,              ONLY: DP
        USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
        USE io_global,          ONLY: stdout, ionode, ionode_id
        USE mp,                 ONLY: mp_bcast, mp_sum
        USE mp_wave,            ONLY: splitwf
        USE io_files,           ONLY: outdir, prefix
        USE io_files,           ONLY: empty_file, emptyunit
        USE reciprocal_vectors, ONLY: ig_l2g
        USE gvecw,              ONLY: ngw

        IMPLICIT none

        COMPLEX(DP), INTENT(OUT) :: c_emp(:,:)
        INTEGER,     INTENT(IN)    :: ne

        LOGICAL :: exst
        INTEGER :: ierr, ig, i, iss
        INTEGER :: ngw_rd, ne_rd, ngw_l
        INTEGER :: ngw_g

        CHARACTER(LEN=256) :: fileempty

        COMPLEX(DP), ALLOCATABLE :: ctmp(:)
        !
        ! ... Subroutine Body
        !
        ngw_g    = ngw
        ngw_l    = ngw
        !
        CALL mp_sum( ngw_g, intra_image_comm )
        !
        ALLOCATE( ctmp(ngw_g) )
        !
        empty_file = TRIM(prefix) // '.emp'
        !
        IF( LEN( TRIM( outdir ) ) > 1 ) THEN
          fileempty = TRIM( outdir ) // '/' // TRIM( empty_file )
        ELSE
          fileempty = TRIM( empty_file )
        END IF

        IF ( ionode ) THEN

          INQUIRE( FILE = TRIM(fileempty), EXIST = EXST )

          IF ( EXST ) THEN
            !
            OPEN( UNIT=emptyunit, FILE=TRIM(fileempty), STATUS='OLD', FORM='UNFORMATTED' )
            !
            READ(emptyunit) ngw_rd, ne_rd
            !
            IF( ne > ne_rd ) THEN
              EXST = .false.
              WRITE( stdout,10)  TRIM(fileempty) 
              WRITE( stdout,20)  ngw_rd, ne_rd
              WRITE( stdout,20)  ngw_g, ne
            END IF
            !
          END IF

        END IF

 10     FORMAT('*** EMPTY STATES : wavefunctions dimensions changed  ', A )
 20     FORMAT('*** NGW = ', I8, ' NE = ', I4)

        CALL mp_bcast(exst,   ionode_id, intra_image_comm)
        CALL mp_bcast(ne_rd,  ionode_id, intra_image_comm)
        CALL mp_bcast(ngw_rd, ionode_id, intra_image_comm)

        IF ( exst ) THEN

           DO i = 1, MIN( ne, ne_rd )
              IF ( ionode ) THEN
                  READ(emptyunit) ( ctmp(ig), ig = 1, MIN( SIZE(ctmp), ngw_rd ) )
              END IF
              IF( i <= ne ) THEN
                  CALL splitwf(c_emp(:,i), ctmp, ngw_l, ig_l2g, me_image, &
                       nproc_image, ionode_id, intra_image_comm)
              END IF
           END DO

        END IF

        IF ( ionode .AND. EXST ) THEN
          CLOSE(emptyunit)
        END IF

        readempty_x = exst
        DEALLOCATE(ctmp)

        RETURN
   END FUNCTION readempty_x


!-----------------------------------------------------------------------
   SUBROUTINE writeempty_x( c_emp, ne )
!-----------------------------------------------------------------------

        ! ...   This subroutine writes empty states to unit emptyunit

        USE xml_io_base
        USE kinds,              ONLY: DP
        USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
        USE mp_wave,            ONLY: mergewf
        USE mp,                 ONLY: mp_sum
        USE io_files,           ONLY: empty_file, emptyunit, outdir, prefix
        USE io_global,          ONLY: ionode, ionode_id, stdout
        USE reciprocal_vectors, ONLY: ig_l2g
        USE gvecw,              ONLY: ngw

        COMPLEX(DP), INTENT(IN) :: c_emp(:,:)
        INTEGER,     INTENT(IN) :: ne

        INTEGER :: ig, i, ngw_g, iss, ngw_l
        LOGICAL :: exst
        COMPLEX(DP), ALLOCATABLE :: ctmp(:)
        CHARACTER(LEN=256) :: fileempty
        !
        ! ... Subroutine Body
        !
        ngw_g    = ngw
        ngw_l    = ngw
        !
        CALL mp_sum( ngw_g, intra_image_comm )
        !
        ALLOCATE( ctmp( ngw_g ) )
        !
        empty_file = TRIM(prefix) // '.emp'
        !
        IF( LEN( TRIM( outdir ) ) > 1  ) THEN
          fileempty = TRIM( outdir ) // '/' // TRIM( empty_file )
        ELSE
          fileempty = TRIM( empty_file )
        END IF

        IF( ionode ) THEN
          OPEN( UNIT = emptyunit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'UNFORMATTED' )
          REWIND( emptyunit )
          WRITE (emptyunit)  ngw_g, ne
         ! WRITE( stdout,10)  TRIM(fileempty)
         ! WRITE( stdout,20)  ngw_g, ne
        END IF

 10     FORMAT('*** EMPTY STATES : writing wavefunctions  ', A )
 20     FORMAT('*** NGW = ', I8, ' NE = ', I4)

        DO i = 1, ne
           ctmp = 0.0d0
           CALL MERGEWF( c_emp(:,i), ctmp(:), ngw_l, ig_l2g, me_image, &
                         nproc_image, ionode_id, intra_image_comm )
           IF( ionode ) THEN
              WRITE (emptyunit) ( ctmp(ig), ig=1, ngw_g )
           END IF
        END DO

        IF( ionode ) THEN
          CLOSE (emptyunit)
        END IF

        DEALLOCATE(ctmp)

        RETURN
   END SUBROUTINE writeempty_x


