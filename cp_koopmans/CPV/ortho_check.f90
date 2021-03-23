!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
SUBROUTINE ortho_check_cmplx( c0_emp, lgam )
  !-----------------------------------------------------------------------
  !
  ! ...  This subroutine checks the projectability of each empty
  ! ...  wavefunction on the occupied manifold (for complex wfc)
  !
  USE io_global,              ONLY : stdout
  USE kinds,                  ONLY : DP
  USE wavefunctions_module,   ONLY : c0
  USE gvecw,                  ONLY : ngw
  USE electrons_base,         ONLY : nbspx, nspin
  USE electrons_module,       ONLY : nupdwn_emp
  USE mp_global,              ONLY : intra_pool_comm
  USE mp,                     ONLY : mp_sum
  USE reciprocal_vectors,     ONLY : gstart
  !
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: c0_emp(:,:)
  LOGICAL, INTENT(IN) :: lgam
  !
  INTEGER :: n_empx
  INTEGER :: m, n, ig
  REAL(DP) :: sqmod
  REAL(DP) :: proj, proj_tot
  !
  !
  WRITE( stdout, '(/, A)' ) " -----------------------------------------"
  WRITE( stdout, '(A)'    ) " Projectability EMP states on OCC manifold"
  WRITE( stdout, '(A, /)' ) " -----------------------------------------"
  !
  !
  n_empx = nupdwn_emp( 1 )
  IF ( nspin == 2 ) n_empx = n_empx + nupdwn_emp( 2 )
  n_empx = n_empx + MOD( n_empx, 2 )
  !
  proj_tot = 0.D0
  !
  !
  DO m = 1, n_empx
    !
    proj = 0.D0
    !
    DO n = 1, nbspx
      DO ig = 1, ngw
        !
        sqmod = c0_emp(ig,m) * CONJG( c0(ig,n) )
        !
        IF ( lgam .AND. ( gstart .NE. 2 .OR. ig .GT. 1 ) ) THEN
          sqmod = sqmod * 2
        ENDIF
        !
        proj = proj + sqmod
        !
      ENDDO
    ENDDO
    !
    CALL mp_sum( proj, intra_pool_comm )
    WRITE( stdout, 100 ) m, proj
    proj_tot = proj_tot + proj
    !
  ENDDO
  !
  WRITE( stdout, 101 ) proj_tot
  !
  !
100 FORMAT( 4X, "orbital # ", I4, " : ", F12.8 )
101 FORMAT( /, 6X, "Total projectability = ", F12.8, / )
  !
  !
END SUBROUTINE ortho_check_cmplx
