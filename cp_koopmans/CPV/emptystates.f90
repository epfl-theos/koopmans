        
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
      USE constants,            ONLY : autoev
      USE control_flags,        ONLY : iprsta, tsde, program_name, gamma_only, do_wf_cmplx, &
                                           tortho
      USE io_global,            ONLY : ionode, stdout
      USE cp_main_variables,    ONLY : eigr, ema0bg, collect_lambda, &
                                       rhor, rhog, rhos, eigr, eigrb, irb, bec, bec_emp
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
                                       bec_csv, readempty_twin, writeempty_twin, &
                                       write_hamiltonian
      USE mp,                   ONLY : mp_comm_split, mp_comm_free, mp_sum
      USE mp_global,            ONLY : intra_image_comm, me_image
      USE nksic,                ONLY : do_orbdep, do_pz, do_wxd, vsicpsi, wtot, sizwtot, &
                                       odd_alpha, valpsi, nkscalfact, odd_alpha_emp
      USE nksic,                ONLY : do_spinsym, pink_emp, allocate_nksic_empty
      USE hfmod,                ONLY : do_hf, vxxpsi
      USE twin_types !added:giovanni
      USE control_flags,        ONLY : tatomicwfc, trane, ndr, ndw
      USE electrons_module,     ONLY : wfc_centers_emp, wfc_spreads_emp, icompute_spread
      USE core,                 ONLY : nlcc_any, rhoc
      USE input_parameters,     ONLY : odd_nkscalfact_empty, &
                                       restart_from_wannier_cp, wannier_empty_only, &
                                       fixed_band, print_wfc_anion, wo_odd_in_empty_run, &
                                       odd_nkscalfact, index_empty_to_save, write_hr
      USE wavefunctions_module, ONLY : c0fixed_emp
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
      REAL(DP),    ALLOCATABLE :: becsum_emp(:,:,:)
      type(twin_matrix) :: bephi_emp! !modified:giovanni
      type(twin_matrix) :: becp_emp !modified:giovanni
      type(twin_matrix) :: bec_occ !(:,:) !modified:giovanni
      type(twin_matrix), dimension(:),  ALLOCATABLE :: lambda_emp !(:,:,:) !, 
      REAL(DP),    ALLOCATABLE :: f_emp(:)
      REAL(DP),    ALLOCATABLE :: f_aux(:)
      REAL(DP),    ALLOCATABLE :: lambda_rep(:,:)
      COMPLEX(DP), ALLOCATABLE :: lambda_rep_c(:,:)
      INTEGER,     ALLOCATABLE :: ispin_emp(:)
      REAL(DP),    ALLOCATABLE :: fsic_emp(:)
      REAL(DP),    ALLOCATABLE :: vsic_emp(:,:)
      REAL(DP),    ALLOCATABLE :: wxd_emp(:,:)
      REAL(DP),    ALLOCATABLE :: deeq_sic_emp(:,:,:,:)
      COMPLEX(DP), ALLOCATABLE :: vxxpsi_emp(:,:)
      REAL(DP),    ALLOCATABLE :: exx_emp(:)
      REAL(DP),    ALLOCATABLE :: old_odd_alpha(:)
      !
      INTEGER, SAVE :: np_emp(2), me_emp(2), emp_comm, color
      INTEGER, SAVE :: desc_emp( descla_siz_ , 2 )
      LOGICAL, SAVE :: first = .true.
      LOGICAL :: lgam !added:giovanni
      LOGICAL :: done_extra !added:giovanni
      COMPLEX(DP), PARAMETER :: c_zero=CMPLX(0.d0,0.d0)
      INTEGER :: sizvsic_emp
      INTEGER :: ndr_loc, ndw_loc
      !
      LOGICAL :: odd_nkscalfact_old
      INTEGER :: nbnd_, ib, start_is
      COMPLEX(DP), ALLOCATABLE :: c0_anion(:,:)
      !
      lgam=gamma_only.and..not.do_wf_cmplx
      !
      odd_nkscalfact_old = odd_nkscalfact
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
      ! restart directories
      !
      IF ( first ) THEN
          ndr_loc = ndr
          ndw_loc = ndw
      ELSE
          ndr_loc = ndw
          ndw_loc = ndw
      ENDIF
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
      !
      IF ( wannier_empty_only .and. odd_nkscalfact_empty) THEN
         !
         ALLOCATE( c0fixed_emp( ngw, n_empx ) )
         !
      ENDIF
      !
      call init_twin(bec_emp, lgam)
      call allocate_twin(bec_emp, nkb, n_emps, lgam)
      call init_twin(becp_emp, lgam)
      call allocate_twin(becp_emp, nkb, n_emps, lgam)
      call init_twin(bec_occ, lgam)
      call allocate_twin(bec_occ, nkb, n_occs, lgam)
      call init_twin(bephi_emp, lgam)
      call allocate_twin(bephi_emp, nkb, n_emps, lgam)
      !
      ALLOCATE(lambda_emp(nspin))
      ! 
      DO iss=1,nspin
          CALL init_twin(lambda_emp(iss), lgam)
          CALL allocate_twin(lambda_emp(iss), nlam_emp, nlam_emp, lgam)
      ENDDO
      !
      ALLOCATE( f_emp( n_empx ) )
      ALLOCATE( f_aux( n_empx ) )
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
      !
      DO iss=1,nspin
         call set_twin(lambda_emp(iss),c_zero)
      ENDDO
      !
      f_emp      = 2.0d0 / DBLE(nspin)
      f_aux      = 1.0d0
      !
      ispin_emp = 0
      ispin_emp( 1:nupdwn_emp( 1 ) ) = 1
      IF( nspin == 2 ) ispin_emp( iupdwn_emp(2) : ) = 2
      !
      IF ( do_orbdep ) THEN
          !
          ALLOCATE( fsic_emp( n_empx ) )
          ! n_empx_odd=n_empx
          ALLOCATE( vsic_emp(nnrx, n_empx) )
          ALLOCATE( wxd_emp (nnrx, 2) )
          ALLOCATE( deeq_sic_emp (nhm,nhm,nat,n_empx) )
          ALLOCATE( becsum_emp(nhm*(nhm+1)/2,nat,nspin))
          CALL allocate_nksic_empty(n_empx)
          sizvsic_emp=nnrx
          !
          fsic_emp = 0.0d0
          vsic_emp = 0.0d0
          wxd_emp  = 0.0d0
          ! 
      ELSE
          !
          ALLOCATE( fsic_emp( n_empx ) )
          ! n_empx_odd=1
          ALLOCATE( vsic_emp(1, n_empx) )
          ALLOCATE( wxd_emp (1, 2) )
          ALLOCATE( deeq_sic_emp (nhm,nhm,nat,n_empx) )
          ALLOCATE( becsum_emp(nhm*(nhm+1)/2,nat,nspin) )
          !
          call allocate_nksic_empty(n_empx)
          sizvsic_emp=1
          !
      ENDIF
      !
      IF ( do_hf ) THEN
          !
          !ALLOCATE( fsic_emp(n_empx ) )
          ALLOCATE( vxxpsi_emp( ngw, n_empx) )
          ALLOCATE( exx_emp( n_empx ) )
          !
          !fsic_emp   = 0.0d0
          vxxpsi_emp = 0.0d0
          exx_emp    = 0.0d0
          !
      ENDIF
      !
      CALL prefor( eigr, vkb )
      !
      CALL nlsm1 ( n_occs, 1, nvb, eigr, c0, bec_occ, 1, lgam )
      !
      ! here is initialize wfcs
      !
      IF ( wannier_empty_only .and. (.not. odd_nkscalfact_empty) ) THEN
         !
         ! (1) read canonical orbital evctot_empty
         !
         exst = readempty( c0_emp, n_empx, ndr_loc )
         !
         IF (exst) THEN
            !
            ! (2) tranform evc to wannier evc0
            ! 
            CALL wave_init_wannier_cp (c0_emp, ngw, n_empx, .false.)
            ! 
         ELSE
            !
            CALL errore( 'empty_run ', 'A restart from wannier orbitals needs evctot_empty in ./our dir', 1 )
            !
         ENDIF
         !
      ELSEIF ( wannier_empty_only .and. odd_nkscalfact_empty) THEN 
         !
         ! (3) read wannier evc0_empty  
         !
         !exst = readempty_twin( c0_emp, n_empx, ndr_loc, .true., .false.)
         !
         !exst = readempty_twin( c0fixed_emp, n_empx, ndr_loc, .false., .false.)
         !
      ELSE
         !
         ! here is restart from the minimizing empty orbitals (evctot_empty)
         ! it can be the wannierized orbitals
         !  
         IF (.not. wo_odd_in_empty_run)THEN
            exst = readempty_twin( c0_emp, n_empx, ndr_loc )
         ELSE
            exst = readempty( c0_emp, n_empx, ndr_loc )
         ENDIF  
         !
         IF ( .NOT. exst ) THEN
            !
            write(stdout, * ) 'Linh: oopp restart from minimizing orbital does not work for emptystate'
            write(stdout, * ) 'Linh: initialize random states and orthogonalize to filled ones'
            !
            ! ...  initial random states orthogonal to filled ones
            !
            IF ( .NOT. do_spinsym .OR. nspin == 1 ) THEN
               !
               IF (tatomicwfc) THEN
                  !
                  CALL wave_atom_init( c0_emp, n_emps, 1 )
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
               IF (tatomicwfc) THEN
                  !
                  CALL wave_atom_init( c0_emp, nupdwn_emp(1), 1 )
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
               !
               c0_emp( 1, : ) = (0.0d0, 0.0d0)
               ! 
            ENDIF
            !
            CALL nlsm1 ( n_emps, 1, nvb, eigr, c0_emp, bec_emp, 1, lgam )
            !
            DO iss = 1, nspin
               !
               in_emp = iupdwn_emp(iss)
               !
               issw   = iupdwn(iss)
               !
               IF (nupdwn(iss)>0.and.nupdwn_emp(iss)>0) THEN
                  !
                  CALL gram_empty( .false., eigr, vkb, bec_emp, bec_occ, nkb, &
                                  c0_emp( :, in_emp: ), c0( :, issw: ), ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, issw )
                  !
               ENDIF
               !
            ENDDO
            !
         ELSE
            !
            write(stdout, * ) 'Linh: the code restarts not random wfc'
            !
         ENDIF
         !
      ENDIF
      !
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
      ! 
      cm_emp = c0_emp

      ekinc_old = 0.0
      ekinc     = 0.0
      !
      ! init xd potential
      !
      ! we need to use wtot from previous calls with occupied states
      ! we save here wtot in wxd_emp
      !
      IF ( do_orbdep .and. (.not.wo_odd_in_empty_run) ) THEN
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
      IF ( do_orbdep .and. (.not. wo_odd_in_empty_run)  ) THEN
         !
         IF (odd_nkscalfact_empty) THEN
            !
            allocate (old_odd_alpha(nbspx) )
            old_odd_alpha(:) = odd_alpha(:)
            ! here, deallocate the memory of odd_alpha for occupied states
            if(allocated(odd_alpha)) deallocate(odd_alpha)
            if(allocated(valpsi)) deallocate(valpsi)
            !
            ! reallocate the memory of odd_alpha for empty states
            allocate (odd_alpha(n_empx))
            allocate (valpsi(n_empx, ngw))
            !
         ENDIF
         !
      ENDIF
      !
      IF( ionode ) THEN
          WRITE( stdout, "(/,3X,'Empty states minimization starting ', &
                        & /,3x,'nfi         dekinc         ekinc' )")
      ENDIF
      !
      IF (tcg_) THEN ! compute empty states with conjugate gradient
         ! 
         call runcg_uspp_emp(c0_emp, cm_emp, bec_emp, f_emp, fsic_emp, n_empx,&
                             n_emps, ispin_emp, iupdwn_emp, nupdwn_emp, phi_emp, lambda_emp, &
                             max_emp, wxd_emp, vsic_emp, sizvsic_emp, pink_emp, nnrx, becsum_emp, &
                             deeq_sic_emp, nudx_emp, eodd_emp, etot_emp, v, &
                             nfi, .true., .true., eigr, bec, irb, eigrb, &
                             rhor, rhog, rhos, rhoc, ema0bg, desc_emp)     !!! Added rhoc NICOLA 
         !
      ELSE ! compute empty states with damped dynamics
         !
         done_extra=.false.
         !
         ITERATIONS: DO iter = 1, max_emp
            !
            IF ( do_orbdep .and. (.not. wo_odd_in_empty_run)  ) THEN 
                !
                IF (odd_nkscalfact_empty) THEN
                   !
                   valpsi(:,:)  = (0.0_DP, 0.0_DP)
                   odd_alpha(:) = 0.0_DP
                   !
                   CALL odd_alpha_routine(c0_emp, n_emps, n_empx, lgam, .true.)
                   !
                   
                ELSE
                   !
                   ! here, we want to use only one value alpha for all empty states,
                   ! that value alpha is defined from in input file. 
                   ! This require to deactive the odd_nkscalfact here so that 
                   ! it does not do odd_alpha in nksic_potential.
                   !  
                   odd_nkscalfact = .false. 
                   ! 
                ENDIF
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
                !IF(done_extra.or.iter==max_emp) THEN
                   !
                   icompute_spread=.true.
                   !
                !ENDIF
                !
                call nksic_potential( n_emps, n_empx, c0_emp, fsic_emp, &
                                      bec_emp, becsum_emp, deeq_sic_emp, &
                                      ispin_emp, iupdwn_emp, nupdwn_emp, rhor, rhoc, &
                                      wtot, sizwtot, vsic_emp, .false., pink_emp, nudx_emp, &
                                      wfc_centers_emp, wfc_spreads_emp, &
                                      icompute_spread, .false.)
                !
                ! line below removed by Giovanni, introduced do_wxd=.false. into call to nksic_potential
                !do_wxd = do_wxd_
                !
                ! Print spreads infor
                !
                WRITE( stdout, *) "sum spreads:1", sum(wfc_spreads_emp(1:nupdwn_emp(1), 1, 1)), &
                                                   sum(wfc_spreads_emp(1:nupdwn_emp(2), 2, 1))
                WRITE( stdout, *) "sum spreads:2", sum(wfc_spreads_emp(1:nupdwn_emp(1), 1, 2)), &
                                                   sum(wfc_spreads_emp(1:nupdwn_emp(2), 2, 2))
                !
                DO i = 1, n_emps
                    !  
                    ! Here wxd_emp <-> wtot that computed from nksic_potential of occupied states.
                    ! wtot is scaled with nkscalfact constant, we thus need to rescaled it here with
                    ! odd_alpha
                    !
                    IF(odd_nkscalfact_empty) wxd_emp(:,:) = wxd_emp(:,:)*odd_alpha(i)/nkscalfact 
                    !  
                    vsic_emp(:,i) = vsic_emp(:,i) + wxd_emp(:, ispin_emp(i))
                    !
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
            ! standard terms
            !
            DO i = 1, n_emps, 2
                !
                CALL dforce( i, bec_emp, vkb, c0_emp, c2, c3, v, SIZE(v,1), ispin_emp, f_aux, n_emps, nspin )
                !
                ! ODD terms
                !
                IF ( do_orbdep .and. (.not.wo_odd_in_empty_run) ) THEN
                    !
                    IF ( odd_nkscalfact_empty ) THEN
                       !
                       c2(:) = c2(:) - valpsi(i, :)   * f_aux(i)
                       c3(:) = c3(:) - valpsi(i+1, :) * f_aux(i+1)
                       !
                    ENDIF
                    !   
                    CALL nksic_eforce( i, n_emps, n_empx, vsic_emp, deeq_sic_emp, bec_emp, ngw, &
                                       c0_emp(:,i), c0_emp(:,i+1), vsicpsi, lgam )
                    !
                    c2(:) = c2(:) - vsicpsi(:,1) * f_aux(i)
                    c3(:) = c3(:) - vsicpsi(:,2) * f_aux(i+1)
                    !
                ENDIF
                !
                ! HF terms
                !
                IF ( do_hf ) THEN
                    !
                    c2(:) = c2(:) - vxxpsi_emp(:,i)   * f_aux(i)
                    c3(:) = c3(:) - vxxpsi_emp(:,i+1) * f_aux(i+1)
                    !
                ENDIF
                !
                IF( tsde ) THEN
                    CALL wave_steepest( cm_emp(:, i  ), c0_emp(:, i  ), emaver, c2 )
                    CALL wave_steepest( cm_emp(:, i+1), c0_emp(:, i+1), emaver, c3 )
                ELSE
                    CALL wave_verlet( cm_emp(:, i  ), c0_emp(:, i  ), verl1, verl2, emaver, c2 )
                    CALL wave_verlet( cm_emp(:, i+1), c0_emp(:, i+1), verl1, verl2, emaver, c3 )
                ENDIF
                !
                IF(lgam) THEN
                   IF ( gstart == 2) THEN
                       cm_emp(1,  i)=CMPLX(DBLE(cm_emp(1,  i)),0.d0)
                       cm_emp(1,i+1)=CMPLX(DBLE(cm_emp(1,i+1)),0.d0)
                   ENDIF
                ENDIF
                !
            ENDDO
            ! 
            ! ortho cm_emp with c0
            !     
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
                            cm_emp( :, in_emp: ), c0( :, issw: ), ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, issw )
                   !
                ENDIF
                !
            ENDDO
            !
            ! ... calphi calculates phi
            ! ... the electron mass rises with g**2
            !
            CALL calphi( c0_emp, ngw, bec_emp, nkb, vkb, phi_emp, n_emps, lgam, ema0bg )
            !
            IF ( tortho ) THEN
               !
               CALL ortho_cp_twin( eigr(1:ngw,1:nat), cm_emp(1:ngw,1:n_emps), phi_emp(1:ngw,1:n_emps), &
                                   ngw, lambda_emp(1:nspin), desc_emp(1:descla_siz_,1:nspin), &
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
            !
            dek = ekinc - ekinc_old
            !
            IF( ionode ) WRITE( stdout,113) ITER, dek, ekinc
            !    
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
            ! 
            ekinc_old = ekinc
            !
         ENDDO ITERATIONS
         !  
      ENDIF !if clause to choose between cg and damped dynamics
      !
      IF ( ionode ) WRITE( stdout, "()")
      ! 
      ! ...  Compute eigenvalues and bring wave functions on Kohn-Sham orbitals
      !
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
          !
          IF(.not.lambda_emp(iss)%iscmplx) THEN
              CALL collect_lambda( lambda_rep, lambda_emp(iss)%rvec(:,:), desc_emp( :, iss ) )
!              IF (iss==1) THEN
!                OPEN(27,FILE='hamiltonian_emp.dat',FORM='formatted',status='UNKNOWN')
!                WRITE(27,'(1E16.10)') lambda_rep
!                CLOSE(27)
!              ENDIF
              CALL crot( cm_emp, c0_emp, ngw, n, i, i, lambda_rep, nudx_emp, ei_emp(:,iss) )
              IF ( write_hr ) CALL write_hamiltonian( lambda_rep, n, iss, .true. )
          ELSE
              CALL collect_lambda( lambda_rep_c, lambda_emp(iss)%cvec(:,:), desc_emp( :, iss ) )
              CALL crot( cm_emp, c0_emp, ngw, n, i, i, lambda_rep_c, nudx_emp, ei_emp(:,iss) )
              IF ( write_hr ) CALL write_hamiltonian( lambda_rep_c, n, iss, .true. )
          ENDIF
          !   
          ei_emp( 1:n, iss ) = ei_emp( 1:n, iss ) / f_aux( i : i + n - 1 )
          !
      ENDDO
      !
      IF(.not.lambda_emp(1)%iscmplx) THEN
          DEALLOCATE( lambda_rep)
      ELSE
          DEALLOCATE( lambda_rep_c)
      ENDIF
      !
!write(stdout,*) "Empty Eigenvalues"
!write(stdout,*) ei_emp
      !
      IF (do_orbdep .and. wo_odd_in_empty_run ) THEN
         !
         CALL empty_koopmans_pp ( n_emps, ispin_emp, cm_emp)
         !
      ENDIF
      !
      ! ...   Save canonical empty orbitals to disk
      !
      CALL writeempty( cm_emp, n_empx, ndw_loc )
      ! 
      ! ...   Save minimizing empty orbitals to disk
      !
      CALL writeempty_twin( c0_emp, n_empx, ndw_loc, .true. )
      !
      IF (print_wfc_anion) THEN
         !
         write(stdout,*) "\n Writing on file the anion WFC \n" 
         ! ...   Save N+1 orbitals to disk to be used in the future anion calculation
         !
         IF ( nupdwn(1) == nupdwn(2) ) THEN !! NsC

           nbnd_    = nupdwn(1)! Why??? only works if spin_up = spin_dw... NsC
           !
           write(stdout,*) " Case nup = ndw", nupdwn(1), nupdwn(2), (nbnd_+1) * nspin
           ! 
           allocate(c0_anion(ngw, (nbnd_+1) * nspin))
           !
           do iss = 1, nspin
              !
              ib = iupdwn( iss )
              !
              start_is = iupdwn( iss )
              !  
              if (iss == 2) start_is = start_is + 1
              !
              c0_anion(:, start_is:start_is + nbnd_ - 1 ) = c0(:, ib:ib + nbnd_ - 1)
              !
              if (iss == 1 ) then
                 !
                 c0_anion(:, start_is + nbnd_) =  c0_emp(:, index_empty_to_save)
                 !
              endif
              !
           enddo
           !
           call writeempty_twin( c0_anion, ((nbnd_+1) * nspin), ndw_loc, .false.)
           !
           deallocate(c0_anion)
           !
         ELSE ! here the case for nup>ndw
           !
           nbnd_    = nupdwn(1)+nupdwn(2) ! Not sure it works for REAL wcf when nup is ODD
           ! 
           write(stdout,*) " Case nup > ndw", nupdwn(1), nupdwn(2), nbnd_
           write(stdout,*) " Warning: NOT TESTED!"
           !
           allocate(c0_anion(ngw, (nbnd_+1) ))
           !
           do iss = 1, nspin
              !
              ib = iupdwn( iss )
              !
              start_is = iupdwn( iss )
              !  
              write(stdout,*) " iss, start_is ", iss, start_is
              !if (iss == 2) start_is = start_is + 1
              !
              c0_anion(:, start_is:start_is + nupdwn(1) - 1 ) = c0(:, ib:ib + nupdwn(1) - 1)
              !
              if (iss == 2 ) then
                 !
                 c0_anion(:, start_is + nupdwn(2)) =  c0_emp(:, index_empty_to_save)
                 !
              endif
              !
           enddo
           !
           call writeempty_twin( c0_anion, (nbnd_+1), ndw_loc,.false.)
           !
           deallocate(c0_anion)
           !
         ENDIF  !! NsC
         !
      ENDIF
      !
      odd_nkscalfact = odd_nkscalfact_old
      !
      IF ( do_orbdep .and. (.not. wo_odd_in_empty_run)  ) THEN
         !
         IF (odd_nkscalfact_empty) THEN
            !
            odd_alpha_emp(:) = odd_alpha(:)
            ! here, deallocate the memory of odd_alpha for empty states
            if(allocated(odd_alpha)) deallocate(odd_alpha)
            ! reallocate the memory of odd_alpha for occupied states
            allocate (odd_alpha(nbspx))
            !
            odd_alpha (:) = old_odd_alpha(:)
            ! 
            deallocate (old_odd_alpha)
            ! 
         ENDIF
         !
      ENDIF 
      ! 
      DEALLOCATE( ispin_emp )
      DEALLOCATE( f_emp )
      DEALLOCATE( f_aux )
      DEALLOCATE( emadt2 )
      DEALLOCATE( emaver )
      DEALLOCATE( c2 )
      DEALLOCATE( c3 )
      DEALLOCATE( c0_emp )
      DEALLOCATE( cm_emp )
      DEALLOCATE( phi_emp )
      IF(ALLOCATED (c0fixed_emp))DEALLOCATE( c0fixed_emp )
      !
      CALL deallocate_twin(bec_emp)
      CALL deallocate_twin(bec_occ)
      CALL deallocate_twin(bephi_emp)
      CALL deallocate_twin(becp_emp)
      !
      DO iss=1,nspin
         CALL deallocate_twin(lambda_emp(iss))
      ENDDO
      !
      DEALLOCATE( fsic_emp ) 
      !
      DEALLOCATE( vsic_emp ) 
      DEALLOCATE( wxd_emp ) 
      DEALLOCATE( deeq_sic_emp )
      DEALLOCATE( becsum_emp )
      !
      IF ( do_hf ) THEN
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
      INTEGER :: i, k, inl
      LOGICAL :: lgam
      EXTERNAL cscnorm
      ! 
      lgam=gamma_only.and..not.do_wf_cmplx
      !
      ! Quick return if there are either no filled or no empty states with this spin (no need to orthogonalize them)
      !
      IF(n_emp.le.0.or.n_occ.le.0) THEN
         !
         return
         !
      ENDIF
      ! 
      call init_twin(csc_emp, lgam)
      call allocate_twin(csc_emp, n_emp, 1, lgam)
      !   
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
         CALL smooth_csv( c_emp(1:ngwx,i), c_occ(1:ngwx,1:n_occ), ngwx, csc_occ, n_occ )
         !
         IF( .NOT. tortho ) THEN
             !
             ! compute scalar product csc_emp(k) = <c_emp(i)|c_emp(k)>
             !
             CALL smooth_csv( c_emp(1:,i), c_emp(1:,1:), ngwx, csc_emp, i-1 )
             !
             IF(.not.csc_emp%iscmplx) THEN
                 CALL mp_sum( csc_emp%rvec(1:n_emp,1:1), intra_image_comm )
             ELSE
                 CALL mp_sum( csc_emp%cvec(1:n_emp,1:1), intra_image_comm )
             ENDIF
             !
         ENDIF
         !
         IF (.not.csc_occ%iscmplx) THEN
            CALL mp_sum( csc_occ%rvec(1:n_occ,1:1), intra_image_comm )
         ELSE
             CALL mp_sum( csc_occ%cvec(1:n_occ,1:1), intra_image_comm )
         ENDIF
         !
         IF ( nvb > 1 ) THEN
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
            IF (.not.bec_emp%iscmplx) THEN
               !
               DO k = 1, n_occ
                  DO inl = 1, nkbx
                     bec_emp%rvec( inl, i ) = bec_emp%rvec( inl, i ) - csc_occ%rvec(k,1) * bec_occ%rvec( inl, k )
                  ENDDO
               ENDDO
               !
            ELSE
               !
               DO k = 1, n_occ
                  DO inl = 1, nkbx
                     bec_emp%cvec( inl, i ) = bec_emp%cvec( inl, i ) - CONJG(csc_occ%cvec(k,1)) * bec_occ%cvec( inl, k )
                  ENDDO
               ENDDO
               !
            ENDIF
            !
            IF( .NOT. tortho ) THEN
               IF(.not.bec_emp%iscmplx) THEN
                  DO k = 1, i-1
                     DO inl = 1, nkbx
                        bec_emp%rvec( inl, i ) = bec_emp%rvec( inl, i ) - csc_emp%rvec(k,1) * bec_emp%rvec( inl, k )
                     ENDDO
                  ENDDO
               ELSE
                  DO k = 1, i-1
                     DO inl = 1, nkbx
                        bec_emp%cvec( inl, i ) = bec_emp%cvec( inl, i ) - CONJG(csc_emp%cvec(k,1)) * bec_emp%cvec( inl, k )
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
            !
         ENDIF
         !
         ! calculate orthogonalized c_emp(i) : |c_emp(i)> = |c_emp(i)> - SUM_k    csv(k)|c_occ(k)>
         !                          c_emp(i) : |c_emp(i)> = |c_emp(i)> - SUM_k<i  csv(k)|c_emp(k)>
         !
         IF (.not.csc_occ%iscmplx) THEN
            ! 
            DO k = 1, n_occ
               CALL DAXPY( 2*ngw, -csc_occ%rvec(k,1), c_occ(:,k), 1, c_emp(:,i), 1 )!warning:giovanni tochange
            ENDDO
            ! 
            IF (.NOT. tortho ) THEN
               DO k = 1, i-1
                  CALL DAXPY( 2*ngw, -csc_emp%rvec(k,1), c_emp(:,k), 1, c_emp(:,i), 1 )!warning:giovanni tochange
               END DO
            ENDIF
            !
         ELSE
            !
            DO k = 1, n_occ
               CALL ZAXPY( ngw, -csc_occ%cvec(k,1), c_occ(:,k), 1, c_emp(:,i), 1 )
            ENDDO
            !
            IF (.NOT. tortho ) THEN
               DO k = 1, i-1
                  CALL ZAXPY( ngw, -csc_emp%cvec(k,1), c_emp(:,k), 1, c_emp(:,i), 1 )
               ENDDO
            ENDIF
            !
         ENDIF
         !
         IF ( .NOT. tortho ) THEN
            !   
            anorm = cscnorm( bec_emp, nkbx, c_emp, ngwx, i, n_emp, lgam)
            !   
            IF (.not.bec_emp%iscmplx) THEN
               !
               CALL DSCAL( 2*ngw, 1.0d0/anorm, c_emp(:,i), 1 )
               !
               IF ( nvb > 1 ) THEN
                  CALL DSCAL( nkbx, 1.0d0/anorm, bec_emp%rvec(:,i), 1 )
               ENDIF
               !
            ELSE
               ! 
               CALL ZSCAL( ngw, CMPLX(1.0d0/anorm,0.d0), c_emp(:,i), 1 )
               !
               IF ( nvb > 1 ) THEN
                  CALL ZSCAL( nkbx, CMPLX(1.0d0/anorm,0.d0), bec_emp%cvec(:,i), 1 )
               ENDIF
               !
            ENDIF
            !
         ENDIF
         !
      ENDDO
      !
      call deallocate_twin( csc_emp )
      call deallocate_twin( csc_occ )
      !
      RETURN
      !
   END SUBROUTINE gram_empty_twin_x

!-----------------------------------------------------------------------
   LOGICAL FUNCTION readempty_x( c_emp, ne, ndi )
!-----------------------------------------------------------------------
        ! 
        ! ...   This subroutine reads canonical empty states from unit emptyunit
        !
        USE kinds,              ONLY: DP
        USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
        USE io_global,          ONLY: stdout, ionode, ionode_id
        USE mp,                 ONLY: mp_bcast, mp_sum
        USE mp_wave,            ONLY: splitwf
        USE io_files,           ONLY: outdir, prefix
        USE io_files,           ONLY: empty_file, emptyunit
        USE reciprocal_vectors, ONLY: ig_l2g
        USE gvecw,              ONLY: ngw
        USE xml_io_base,        ONLY: restart_dir, wfc_filename
        USE electrons_base,     ONLY: nspin
        USE electrons_module,   ONLY: iupdwn_emp, nupdwn_emp

        IMPLICIT none

        COMPLEX(DP), INTENT(OUT) :: c_emp(:,:)
        INTEGER,     INTENT(IN)  :: ne
        INTEGER,     INTENT(IN)  :: ndi

        LOGICAL :: exst
        INTEGER :: ierr, ig, i, iss
        INTEGER :: ngw_rd, ne_rd, ngw_l
        INTEGER :: ngw_g

        CHARACTER(LEN=256) :: fileempty, dirname

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
        dirname   = restart_dir( outdir, ndi )
        !
        DO iss = 1, nspin
          IF ( ionode ) THEN
            fileempty = TRIM( wfc_filename( dirname, 'evc_empty', 1, iss ) )
            INQUIRE( FILE = TRIM(fileempty), EXIST = EXST )

            IF ( EXST ) THEN
              !
              OPEN( UNIT=emptyunit, FILE=TRIM(fileempty), STATUS='OLD', FORM='UNFORMATTED' )
              !
              READ(emptyunit) ngw_rd, ne_rd
              !
              IF( nupdwn_emp(iss) > ne_rd ) THEN
                EXST = .false.
                WRITE( stdout,10)  TRIM(fileempty) 
                WRITE( stdout,20)  ngw_rd, ne_rd
                WRITE( stdout,20)  ngw_g, nupdwn_emp(iss)
              END IF
              !
            END IF
            !
          END IF

 10       FORMAT('*** EMPTY STATES : wavefunctions dimensions changed  ', A )
 20       FORMAT('*** NGW = ', I8, ' NE = ', I4)

          CALL mp_bcast(exst,   ionode_id, intra_image_comm)
          CALL mp_bcast(ne_rd,  ionode_id, intra_image_comm)
          CALL mp_bcast(ngw_rd, ionode_id, intra_image_comm)

          IF ( exst ) THEN

             DO i = 1, nupdwn_emp(iss)
                IF ( ionode ) THEN
                    READ(emptyunit) ( ctmp(ig), ig = 1, MIN( SIZE(ctmp), ngw_rd ) )
                END IF
                CALL splitwf(c_emp(:, i + iupdwn_emp(iss) - 1), ctmp, ngw_l, ig_l2g, me_image, &
                     nproc_image, ionode_id, intra_image_comm)
             END DO

          END IF

          IF ( ionode .AND. EXST ) THEN
            CLOSE(emptyunit)
          END IF

          readempty_x = exst

          IF (.NOT. readempty_x) RETURN
        END DO

        DEALLOCATE(ctmp)

        RETURN
   END FUNCTION readempty_x


!-----------------------------------------------------------------------
   SUBROUTINE writeempty_x( c_emp, ne, ndi )
!-----------------------------------------------------------------------
        !
        ! ...   This subroutine write canonical empty states from unit emptyunit
        !
        USE kinds,              ONLY: DP
        USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
        USE mp_wave,            ONLY: mergewf
        USE mp,                 ONLY: mp_sum
        USE io_files,           ONLY: empty_file, emptyunit, outdir, prefix
        USE io_global,          ONLY: ionode, ionode_id, stdout
        USE reciprocal_vectors, ONLY: ig_l2g
        USE gvecw,              ONLY: ngw
        USE xml_io_base,        ONLY: restart_dir, wfc_filename
        USE wrappers,           ONLY: f_mkdir
        USE electrons_base,     ONLY: nspin
        USE electrons_module,   ONLY: iupdwn_emp, nupdwn_emp
        !
        IMPLICIT NONE

        COMPLEX(DP), INTENT(IN) :: c_emp(:,:)
        INTEGER,     INTENT(IN) :: ne
        INTEGER,     INTENT(IN) :: ndi

        INTEGER :: ig, i, ngw_g, iss, ngw_l
        LOGICAL :: exst, ierr
        COMPLEX(DP), ALLOCATABLE :: ctmp(:)
        CHARACTER(LEN=256) :: fileempty, dirname
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
        dirname   = restart_dir( outdir, ndi )
        !
        ierr = f_mkdir( TRIM(dirname)//"/K00001" )
        !
        DO iss = 1, nspin
          fileempty = TRIM( wfc_filename( dirname, 'evc_empty', 1, iss ) )
          IF( ionode ) THEN
            OPEN( UNIT = emptyunit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'UNFORMATTED' )
            REWIND( emptyunit )
            WRITE (emptyunit)  ngw_g, ne
          END IF

 10       FORMAT('*** EMPTY STATES : writing wavefunctions  ', A )
 20       FORMAT('*** NGW = ', I8, ' NE = ', I4)

          DO i = 1, nupdwn_emp(iss)
             ctmp = 0.0d0
             CALL MERGEWF( c_emp(:,i + iupdwn_emp(iss) - 1), ctmp(:), ngw_l, ig_l2g, me_image, &
                           nproc_image, ionode_id, intra_image_comm )
             IF( ionode ) THEN
                WRITE (emptyunit) ( ctmp(ig), ig=1, ngw_g )
             END IF
          END DO

          IF( ionode ) THEN
            CLOSE (emptyunit)
          END IF
        END DO

        DEALLOCATE(ctmp)

        RETURN
   END SUBROUTINE writeempty_x

!-----------------------------------------------------------------------
LOGICAL FUNCTION reademptytwin_x( c_emp, ne, ndi )
!-----------------------------------------------------------------------
        !
        ! ...   This subroutine reads canonical, or mininimzing empty states
        !
        USE kinds,              ONLY: DP
        USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
        USE io_global,          ONLY: stdout, ionode, ionode_id
        USE mp,                 ONLY: mp_bcast, mp_sum
        USE mp_wave,            ONLY: splitwf
        USE io_files,           ONLY: outdir, prefix
        USE io_files,           ONLY: empty_file, emptyunitc0, emptyunit, emptyunitc0fixed
        USE reciprocal_vectors, ONLY: ig_l2g
        USE gvecw,              ONLY: ngw
        USE xml_io_base,        ONLY: restart_dir, wfc_filename
        USE electrons_base,     ONLY: nspin
        USE electrons_module,   ONLY: iupdwn_emp, nupdwn_emp
        !
        IMPLICIT none
        ! 
        COMPLEX(DP), INTENT(OUT) :: c_emp(:,:)
        INTEGER,     INTENT(IN)  :: ne
        INTEGER,     INTENT(IN)  :: ndi
        ! 
        LOGICAL :: exst
        INTEGER :: ierr, ig, i, iss
        INTEGER :: ngw_rd, ne_rd, ngw_l
        INTEGER :: ngw_g
        ! 
        CHARACTER(LEN=256) :: fileempty, dirname
        !
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
        dirname   = restart_dir( outdir, ndi )
        !
        DO iss = 1, nspin
          fileempty = TRIM( wfc_filename( dirname, 'evc0_empty', 1, iss ) )
          !
          IF ( ionode ) THEN
             !
             INQUIRE( FILE = TRIM(fileempty), EXIST = EXST )
             !
             IF ( exst ) THEN
                !
                OPEN( UNIT=emptyunitc0, FILE=TRIM(fileempty), STATUS='OLD', FORM='UNFORMATTED' )
                !
                READ(emptyunitc0) ngw_rd, ne_rd
                !
                IF ( nupdwn_emp(iss) > ne_rd ) THEN
                   !
                   exst = .false.
                   WRITE( stdout,10)  TRIM(fileempty) 
                   WRITE( stdout,20)  ngw_rd, ne_rd
                   WRITE( stdout,20)  ngw_g, nupdwn_emp(iss)
                   !
                ENDIF
                ! 
             ENDIF
             !
          ENDIF
          !
 10       FORMAT('*** EMPTY STATES : wavefunctions dimensions changed  ', A )
 20       FORMAT('*** NGW = ', I8, ' NE = ', I4)
          !  
          CALL mp_bcast(exst,   ionode_id, intra_image_comm)
          CALL mp_bcast(ne_rd,  ionode_id, intra_image_comm)
          CALL mp_bcast(ngw_rd, ionode_id, intra_image_comm)
          !
          IF ( exst ) THEN
             ! 
             DO i = 1, nupdwn_emp(iss)
                !
                IF ( ionode ) THEN
                   ! 
                   READ(emptyunitc0) ( ctmp(ig), ig = 1, MIN( SIZE(ctmp), ngw_rd ) )
                   !
                ENDIF
                ! 
                CALL splitwf(c_emp(:,i + iupdwn_emp(iss) - 1), ctmp, ngw_l, ig_l2g, me_image, &
                      nproc_image, ionode_id, intra_image_comm)
                !
             ENDDO
             !
          ENDIF
          ! 
          IF ( ionode .AND. exst ) THEN
             !
             CLOSE(emptyunitc0) 
             !
          ENDIF
          !
          reademptytwin_x = exst
          !
          IF (.NOT. reademptytwin_x) RETURN
        END DO
        ! 
        DEALLOCATE(ctmp)
        !
        RETURN
        !
END FUNCTION reademptytwin_x


!-----------------------------------------------------------------------
SUBROUTINE writeemptytwin_x( c_emp, ne, ndi, write_evc0)
!-----------------------------------------------------------------------
        !
        ! ...   This subroutine writes empty states to unit emptyunitc0
        ! 
        USE kinds,              ONLY: DP
        USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
        USE mp_wave,            ONLY: mergewf
        USE mp,                 ONLY: mp_sum
        USE io_files,           ONLY: empty_file, outdir, prefix, &
                                      emptyunitc0, emptyunitc0fixed
        USE io_global,          ONLY: ionode, ionode_id, stdout
        USE reciprocal_vectors, ONLY: ig_l2g
        USE gvecw,              ONLY: ngw
        USE xml_io_base,        ONLY: restart_dir, wfc_filename
        USE wrappers,           ONLY: f_mkdir
        USE electrons_base,     ONLY: nspin
        USE electrons_module,   ONLY: iupdwn_emp, nupdwn_emp
        !
        IMPLICIT NONE
        ! 
        COMPLEX(DP), INTENT(IN) :: c_emp(:,:)
        INTEGER,     INTENT(IN) :: ne
        INTEGER,     INTENT(IN) :: ndi
        LOGICAL,     INTENT(IN) :: write_evc0
        !
        INTEGER :: ig, i, ngw_g, iss, ngw_l, funit, ne_loc, i_start
        LOGICAL :: exst, ierr
        COMPLEX(DP), ALLOCATABLE :: ctmp(:)
        CHARACTER(LEN=256) :: fileempty, dirname
        !
        ! ... Subroutine Body
        !
        ngw_g    = ngw
        ngw_l    = ngw
        !
        IF (write_evc0) THEN
          funit = emptyunitc0
        ELSE
          funit = emptyunitc0fixed
        END IF
        !
        CALL mp_sum( ngw_g, intra_image_comm )
        !
        ALLOCATE( ctmp( ngw_g ) )
        !
        dirname   = restart_dir( outdir, ndi )
        !
        ierr = f_mkdir( TRIM(dirname)//"/K00001" )
        !
        DO iss = 1, nspin
          ! We only split the files by spin if evc0_empty = .true.. If .false.,
          ! we ignore the iss index, do everything at once, and exit after
          ! completing the loop once
          !
          IF (write_evc0) THEN
             !    
             fileempty = TRIM( wfc_filename( dirname, 'evc0_empty', 1, iss ) )
             ne_loc = nupdwn_emp(iss)
             i_start = iupdwn_emp(iss)
             !
          ELSE
             !
             fileempty = TRIM( wfc_filename( dirname, 'evcfixed_empty', 1 ) )
             ne_loc = ne
             i_start = 1
             !
          ENDIF   
          !
          IF ( ionode ) THEN
             OPEN( UNIT = funit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'UNFORMATTED' )
             !
             REWIND( funit )
             !
             WRITE (funit)  ngw_g, ne_loc
          ENDIF
          !
 10       FORMAT('*** EMPTY STATES : writing wavefunctions  ', A )
 20       FORMAT('*** NGW = ', I8, ' NE = ', I4)
          !
          DO i = 1, ne_loc
             !
             ctmp = 0.0d0
             !
             CALL MERGEWF( c_emp(:,i + i_start - 1), ctmp(:), ngw_l, ig_l2g, me_image, &
                           nproc_image, ionode_id, intra_image_comm )
             !
             IF (ionode ) THEN
                !
                WRITE (funit) ( ctmp(ig), ig=1, ngw_g )
                ! 
             ENDIF
             !
          ENDDO
          ! 
          IF ( ionode ) CLOSE(funit)

          ! For evc0_fixed, exit this 'spin' loop
          IF (.NOT. write_evc0) EXIT
        END DO

        DEALLOCATE(ctmp)
        
        !
        RETURN
        ! 
   END SUBROUTINE writeemptytwin_x
