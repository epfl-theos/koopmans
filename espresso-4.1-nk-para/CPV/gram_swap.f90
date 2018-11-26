!-------------------------------------------------------------------------
SUBROUTINE gram_swap( betae, bec, nkbx, cp, ngwx, n, fixed_index )
!-----------------------------------------------------------------------
! 
!     This is a modified gram-schmidt orthogonalization
!     to adapt the case 1 wfc is fixed. 
!     Using the fact that gram-schmidt does not change the first vector,
!     we swap the fixed wfc with the first one, then do as normal,
!     then in the end, we do another swap to return to input order.
!
      USE uspp,           ONLY : nkb, nhsavb=> nkbus
      USE gvecw,          ONLY : ngw
      USE kinds,          ONLY : DP
      USE control_flags,  ONLY : gamma_only, do_wf_cmplx !added:giovanni
      USE electrons_base,  ONLY : ispin
      USE twin_types !added:giovanni
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nkbx, ngwx, n, fixed_index
      type(twin_matrix)   :: bec!( nkbx, n )!modified:giovanni
      COMPLEX(DP)   :: cp( ngwx, n ), betae( ngwx, nkb )
      !
      REAL(DP) :: anorm, cscnorm
      COMPLEX(DP), ALLOCATABLE :: csc( : ) !modified:giovanni
      INTEGER :: i,k, ispin_aux
      LOGICAL :: lgam !i added:giovanni
      EXTERNAL cscnorm
      !
      lgam=gamma_only.and..not.do_wf_cmplx !added:giovanni
      !
      CALL start_clock( 'gram_swap' )
      !
      ALLOCATE( csc( n ) )
      ! 
      csc=CMPLX(0.d0,0.d0)
      !
      WRITE(*,*) "NICOLA in gram_swap fixed_inex=", fixed_index
      CALL dswap( 2*ngw*1, cp(:,fixed_index), 1, cp(:,1), 1 ) 
      ! NsC Exchange the spin indeces as well
      ispin_aux = ispin(fixed_index) ! Store the original spin 
      ispin(1) = ispin(fixed_index)
      ispin(fixed_index) = 1
      !
      DO i = 1, n
         !
         CALL gracsc( bec, nkbx, betae, cp, ngwx, i, csc, n, lgam )!added:giovanni lgam
         !
         ! calculate orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
         !
         DO k = 1, i - 1
            CALL ZAXPY(ngw,-csc(k),cp(1,k),1,cp(1,i),1)!modified:giovanni
         END DO
         !
         anorm = cscnorm( bec, nkbx, cp, ngwx, i, n, lgam)
         !
         ! below is the same CALL ZSCAL( ngw, CMPLX(1.0d0/anorm, 0.d0) , cp(1,i), 1 )
         !  
         DO k = 1, ngw
              cp(k,i) = cp(k,i)*CMPLX(1.0d0/anorm, 0.d0)
         ENDDO 
         !
         ! these are the final bec's
         !
         IF (nkbx > 0 ) THEN
            IF(.not.bec%iscmplx) THEN
              CALL DSCAL( nkbx, 1.0d0/anorm, bec%rvec(1:nkbx,i), 1 )!modified:giovanni
            ELSE
              CALL ZSCAL( nkbx, CMPLX(1.0d0/anorm,0.d0) , bec%cvec(1:nkbx,i), 1 )!added:giovanni
            ENDIF
         ENDIF
         !
      END DO
      !
      DEALLOCATE( csc )
      !
      CALL dswap( 2*ngw*1, cp(:,1), 1, cp(:,fixed_index), 1 ) 
      ! Nsc
      ispin(1) = 1
      ispin(fixed_index) = ispin_aux
      !
      CALL stop_clock( 'gram_swap' )
      !
      RETURN
      !
END SUBROUTINE gram_swap
