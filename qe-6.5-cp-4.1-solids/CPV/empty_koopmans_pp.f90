        
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE empty_koopmans_pp (n_emps_evc, ispin_evc, evc)
!-----------------------------------------------------------------------
      !
      ! Performs the minimization on the empty state subspace keeping the
      ! occupied manyfold fixed. A proper orthogonalization of the two
      ! manyfolds is performed. 
      !
      USE kinds,                ONLY : DP
      USE constants,            ONLY : autoev
      USE control_flags,        ONLY : gamma_only, do_wf_cmplx,tortho, &
                                       ndr, ndw
      USE io_global,            ONLY : ionode, stdout
      USE cp_main_variables,    ONLY : eigr, rhor
      USE core,                 ONLY : nlcc_any, rhoc
      USE electrons_base,       ONLY : nspin, nbspx
      USE uspp,                 ONLY : nkb
      USE uspp_param,           ONLY : nhm
      USE ions_base,            ONLY : nat, nsp
      USE grid_dimensions,      ONLY : nnrx
      USE gvecw,                ONLY : ngw
      USE reciprocal_vectors,   ONLY : ng0 => gstart
      USE cp_interfaces,        ONLY : readempty_twin, writeempty_twin, nlsm1, readempty
      USE mp,                   ONLY : mp_comm_split, mp_comm_free, mp_sum
      USE mp_global,            ONLY : intra_image_comm, me_image
      USE nksic,                ONLY : do_orbdep, do_pz, do_wxd, vsicpsi, wtot, sizwtot, &
                                       odd_alpha, valpsi, nkscalfact, odd_alpha_emp
      USE nksic,                ONLY : allocate_nksic_empty
      USE input_parameters,     ONLY : odd_nkscalfact_empty, odd_nkscalfact, aux_empty_nbnd 
      USE electrons_module,     ONLY : ei_emp 
      USE twin_types 
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)      :: n_emps_evc
      INTEGER, INTENT(IN)      :: ispin_evc(n_emps_evc)
      COMPLEX(DP), INTENT (INOUT) :: evc(ngw, n_emps_evc)
      !
      INTEGER :: i, m, ig, iss, nbnd 
      !
      INTEGER :: n_emp   ! number of empty states from input
      LOGICAL :: exst
      !
      INTEGER :: n_emps, n_empx, nudx_emp
      INTEGER :: nupdwn_emp(nspin) 
      INTEGER :: iupdwn_emp(nspin) 
      !
      COMPLEX(DP), ALLOCATABLE :: c2(:), c3(:)
      COMPLEX(DP), ALLOCATABLE :: c0_emp(:,:), gi(:,:)
      REAL(DP),    ALLOCATABLE :: becsum_emp(:,:,:)
      !
      INTEGER,     ALLOCATABLE :: ispin_emp(:)
      REAL(DP),    ALLOCATABLE :: fsic_emp(:)
      REAL(DP),    ALLOCATABLE :: vsic_emp(:,:)
      REAL(DP),    ALLOCATABLE :: wxd_emp(:,:)
      REAL(DP),    ALLOCATABLE :: deeq_sic_emp(:,:,:,:)
      REAL(DP),    ALLOCATABLE :: wfc_centers_emp(:,:,:)
      REAL(DP),    ALLOCATABLE :: wfc_spreads_emp(:,:,:)
      REAL(DP),    ALLOCATABLE :: old_odd_alpha(:)
      LOGICAL :: icompute_spread
      !
      LOGICAL :: lgam !added:giovanni
      INTEGER :: ndr_loc, ndw_loc
      COMPLEX(DP), PARAMETER :: c_zero=CMPLX(0.d0,0.d0)
      COMPLEX(DP) :: delta_eig(n_emps_evc), scar1, scar2 
      TYPE(twin_matrix)     :: bec_emp
      REAL(dp), allocatable :: pink_emp(:)
      !
      LOGICAL :: odd_nkscalfact_old
      !
      IF (aux_empty_nbnd == 0) THEN
         !
         WRITE( stdout, "(/,3X,' aux_empty_nbnd == 0, no ODD correction ' )")
         RETURN
         ! 
      ENDIF
      !
      lgam = gamma_only.and..not.do_wf_cmplx
      !
      odd_nkscalfact_old = odd_nkscalfact
      !
      ! restart directories
      !
      ndr_loc = ndr
      ndw_loc = ndw
      !
      ! Setting number electrons
      !
      n_emp = aux_empty_nbnd 
      !
      nupdwn_emp(1) = n_emp
      iupdwn_emp(1) = 1
      IF ( nspin == 2 ) THEN
         nupdwn_emp(2) = n_emp
         iupdwn_emp(2) = 1 + n_emp
         ! 
      ENDIF
      !
      n_emps = nupdwn_emp( 1 )
      IF( nspin == 2 ) n_emps = n_emps + nupdwn_emp( 2 )
      !
      nudx_emp = nupdwn_emp( 1 )
      IF( nspin == 2 ) nudx_emp = MAX( nudx_emp, nupdwn_emp( 2 ) )
      !
      n_empx = nupdwn_emp( 1 )
      IF( nspin == 2 ) n_empx = n_empx + nupdwn_emp( 2 )
      n_empx = n_empx + MOD( n_empx, 2)
      !
      ALLOCATE( c0_emp( ngw, n_empx ) )
      ALLOCATE( ispin_emp( n_empx ) )
      !
      call init_twin(bec_emp, lgam)
      call allocate_twin(bec_emp, nkb, n_emps, lgam)
      call set_twin(bec_emp, c_zero)
      !
      ispin_emp (:) = 0
      ispin_emp ( 1:nupdwn_emp( 1 ) ) = 1
      IF ( nspin == 2 ) ispin_emp( iupdwn_emp(2) : ) = 2
      !
      ALLOCATE( fsic_emp( n_empx ) )
      ALLOCATE( vsic_emp(nnrx, n_empx) )
      ALLOCATE( wxd_emp (nnrx, 2) )
      ALLOCATE( deeq_sic_emp (nhm,nhm,nat,n_empx))
      ALLOCATE( becsum_emp (nhm*(nhm+1)/2,nat,nspin))
      ALLOCATE( wfc_centers_emp(4, nudx_emp, nspin )) 
      ALLOCATE( wfc_spreads_emp(nudx_emp, nspin, 2 ))
      ALLOCATE( pink_emp(n_emps))
      ! 
      CALL allocate_nksic_empty(n_empx)
      !
      fsic_emp = 0.0d0
      vsic_emp = 0.0d0
      wxd_emp  = 0.0d0
      !
      ! read auxilary orbitals
      ! 
      exst = readempty_twin( c0_emp, n_empx, ndr_loc )
      !
      IF ( .NOT. exst ) THEN
         !
         CALL errore( 'empty_koopmans_pp', 'there is no auxilary orbital', 1 )
         !
      ENDIF
      !
      CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
      !
      IF ( ionode ) THEN
         ! 
         WRITE( stdout, "(/,3X,'Compute empty ODD correction for KS eigenvalues ' )")
         ! 
      ENDIF
      !  
      ! init xd potential
      !
      ! we need to use wtot from previous calls with occupied states
      ! we save here wtot in wxd_emp
      !
      wxd_emp(:,:) = 0.0_DP
      !
      IF ( do_wxd .AND. .NOT. do_pz ) THEN
         !
         wxd_emp(:,:) = wtot(:,:)
         !
      ENDIF
      !
      IF (odd_nkscalfact_empty) THEN
         !
         allocate(old_odd_alpha (nbspx))
         old_odd_alpha(:) = odd_alpha(:)
         ! here, deallocate the memory of odd_alpha for occupied states
         if(allocated(odd_alpha)) deallocate(odd_alpha)
         if(allocated(valpsi)) deallocate(valpsi)
         !
         ! reallocate the memory of odd_alpha for empty states
         allocate (odd_alpha(n_emps))
         allocate (valpsi(n_emps, ngw))
         !
      ENDIF
      ! 
      IF (odd_nkscalfact_empty) THEN
         !
         valpsi(:,:)  = (0.0_DP, 0.0_DP)
         odd_alpha(:) =  0.0_DP
         !
         CALL odd_alpha_routine(c0_emp, n_emps, n_empx, lgam, .true.)
         !
         odd_nkscalfact = .true. 
         !
      ELSE
         !
         ! here, we want to use only one value alpha for all empty states,
         ! that value alpha is defined from in input file. 
         ! This require to deactive here the odd_nkscalfact so that 
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
      icompute_spread=.true.
      !
      CALL nksic_potential( n_emps, n_empx, c0_emp, fsic_emp, &
                            bec_emp, becsum_emp, deeq_sic_emp, &
                            ispin_emp, iupdwn_emp, nupdwn_emp, rhor, rhoc, &
                            wtot, sizwtot, vsic_emp, .false., pink_emp, nudx_emp, &
                            wfc_centers_emp, wfc_spreads_emp, &
                            icompute_spread, .true.)
      !
      ! Print spreads infor
      !
      WRITE( stdout, *) "sum spreads:1", sum(wfc_spreads_emp(1:nupdwn_emp(1), 1, 1)), &
                                         sum(wfc_spreads_emp(1:nupdwn_emp(2), 2, 1))
      WRITE( stdout, *) "sum spreads:2", sum(wfc_spreads_emp(1:nupdwn_emp(1), 1, 2)), &
                                         sum(wfc_spreads_emp(1:nupdwn_emp(2), 2, 2))

                write(stdout, *) "Localization of orbitals from PZS localization"
                do i = 1, nupdwn_emp(2)
                   write(stdout, *) i, wfc_spreads_emp(i, 2, 2), pink_emp (nupdwn_emp(1)+i)
                enddo
      !
      DO i = 1, n_emps
         !  
         ! Here wxd_emp <-> wtot that computed from nksic_potential of occupied states.
         ! wtot is scaled with nkscalfact constant, we thus need to rescaled it here with
         ! odd_alpha
         !
         IF(odd_nkscalfact_empty) wxd_emp(:, ispin_emp(i)) = wxd_emp(:, ispin_emp(i))*odd_alpha(i)/nkscalfact 
         !  
         vsic_emp(:,i) = vsic_emp(:,i) + wxd_emp(:, ispin_emp(i))
         !
      ENDDO
      ! 
      ! Compute hpsi terms
      !
      allocate( c2(ngw), c3(ngw), gi(ngw, n_emps) )
      gi(:,:) = cmplx(0.0d0, 0.0d0)
      ! 
      DO i = 1, n_emps, 2
         ! 
         c2(:) = cmplx(0.0d0, 0.0d0)
         c3(:) = cmplx(0.0d0, 0.0d0)
         !
         IF ( odd_nkscalfact_empty ) THEN
            !
            c2(:) = c2(:) + valpsi(i, :)
            c3(:) = c3(:) + valpsi(i+1, :)
            !
         ENDIF
         !   
         CALL nksic_eforce( i, n_emps, n_empx, vsic_emp, deeq_sic_emp, bec_emp, ngw, &
                            c0_emp(:,i), c0_emp(:,i+1), vsicpsi, lgam )
         !
         c2(:) = c2(:) + vsicpsi(:,1) 
         c3(:) = c3(:) + vsicpsi(:,2) 
         !
         gi(:, i) = c2(:)
         ! 
         IF (i+1 <= n_emps) gi(:, i+1) = c3(:)
         ! 
         IF (lgam .and. ng0.eq.2) THEN
            ! 
            gi(1, i) = CMPLX(DBLE(gi(1, i)),0.d0)
            !
            IF (i+1 <= n_emps) gi(1, i+1) = CMPLX(DBLE(gi(1, i+1)),0.d0)
            !    
         ENDIF
         !
      ENDDO
      !
      ! compute delta_eig(m) \sum_i <evc(m)|gi(i)><c0_emp(i)|evc(m)>
      !
      delta_eig(:) = cmplx(0.0d0, 0.0d0)
      DO m = 1, n_emps_evc
         !
         DO i = 1, n_emps 
            !
            scar1 = cmplx(0.0d0, 0.0d0)
            scar2 = cmplx(0.0d0, 0.0d0)
            !
            IF (ispin_emp(i) == ispin_evc(m)) THEN
               ! 
               IF (lgam) THEN
                  !
                  DO ig = 1, ngw
                     !
                     scar1 = scar1 + cmplx(dble(conjg(evc(ig, m)) * gi(ig, i)),0.0d0)
                     scar2 = scar2 + cmplx(dble(conjg(c0_emp(ig, i)) * evc(ig, m)), 0.0d0) 
                     !
                  ENDDO
                  !
                  scar1 = 2.d0*scar1
                  scar2 = 2.d0*scar2
                  !
                  IF (ng0.eq.2) THEN
                     !
                     scar1 = cmplx(dble(scar1), 0.0d0) &
                           - cmplx(dble(conjg(evc(1, m)) * gi(1, i)), 0.0d0)
                     ! 
                     scar2 = cmplx(dble(scar2), 0.0d0) &
                           - cmplx(dble(conjg(c0_emp(1, i)) * evc(1, m)), 0.0d0)
                     ! 
                  ELSE
                     !
                     scar1 = cmplx(dble(scar1), 0.0d0)
                     scar2 = cmplx(dble(scar2), 0.0d0)
                     ! 
                  ENDIF
                  ! 
               ELSE
                  !
                  DO ig = 1, ngw
                     !
                     scar1 = scar1 + conjg(evc(ig, m)) * gi(ig, i)
                     scar2 = scar2 + conjg(c0_emp(ig, i)) * evc(ig, m)
                     !
                  ENDDO
                  !  
               ENDIF
               !
               CALL mp_sum( scar1, intra_image_comm )
               CALL mp_sum( scar2, intra_image_comm )
               !
               delta_eig(m) = delta_eig(m) + scar1 * scar2
               !  
            ENDIF
            !
if (m==1) then
write(stdout, *)  'Linh test Amn:', i, scar1, scar2  
endif
            !
         ENDDO
         !
      ENDDO
      !
      DO iss = 1, nspin
         i = 0 
         DO m = 1, n_emps_evc
            IF ( ispin_evc(m) == iss) THEN
               i = i + 1
               WRITE(stdout, *) i, m, ei_emp( i, iss ) * autoev, (ei_emp( i, iss ) + dble(delta_eig(m)))*autoev
               ei_emp(i, iss) = ei_emp(i, iss) + dble(delta_eig(m))
            ENDIF
         ENDDO
      ENDDO
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
         deallocate(old_odd_alpha)
         ! 
      ENDIF
      !
      DEALLOCATE( gi )
      DEALLOCATE( c2 )
      DEALLOCATE( c3 )
      DEALLOCATE( c0_emp )
      !
      DEALLOCATE( ispin_emp )
      DEALLOCATE( fsic_emp ) 
      !
      CALL deallocate_twin(bec_emp)
      !
      DEALLOCATE( vsic_emp ) 
      DEALLOCATE( wxd_emp ) 
      DEALLOCATE( deeq_sic_emp )
      DEALLOCATE( becsum_emp )
      DEALLOCATE( wfc_centers_emp )
      DEALLOCATE( wfc_spreads_emp )
      DEALLOCATE( pink_emp )
      ! 
      RETURN
      !  
   END SUBROUTINE empty_koopmans_pp
