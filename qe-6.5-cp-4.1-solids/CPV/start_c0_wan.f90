!=----------------------------------------------------------------------------=!
   SUBROUTINE start_c0_wan ( ctot, c0, ngw, n, u_matrix, nx )
!=----------------------------------------------------------------------------=!
      ! 
      !  this routine rotates the KS wave functions to very localize minimizing 
      !  wfc. That will be used as starting of KS
      !  it works with a block-like distributed matrix
      !  of the Lagrange multipliers ( lambda ).
      !
      ! ... declare modules
      !
      USE kinds,            ONLY: DP
      ! 
      IMPLICIT NONE
      !  
      ! ... declare subroutine arguments
      ! 
      INTEGER,     INTENT(IN)    :: ngw, n, nx
      COMPLEX,     INTENT(IN)    :: u_matrix(:,:,:)
      COMPLEX(DP), INTENT(INOUT) :: ctot(:,:)
      COMPLEX(DP), INTENT(IN)    :: c0(:,:)
      !       
      ! ... declare other variables
      !
      INTEGER               :: i, j, k
      !
      IF( nx < 1 ) THEN
        RETURN
      END IF
      ! 
      c0(:,:) = (0.0, 0.0)
      !
      DO ispin = 1, nspin
         !
         istart = iupdwn (ispin) 
         !
         DO j = 1, nbnd
            ! 
            DO i = 1, nbnd
               ! 
               CALL DAXPY( 2*ngw, u_matrix(j,i,ispin), c0(1, istart:istart+nbnd-1), 1, ctot(1,istart:istart+nbnd-1), 1 )
               !
            ENDDO
            ! 
         ENDDO
         !  
      ENDDO
      !
      RETURN
      ! 
   END SUBROUTINE start_c0_wan 
