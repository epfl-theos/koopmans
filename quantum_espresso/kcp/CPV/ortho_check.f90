!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
SUBROUTINE ortho_check_cmplx(c0_emp, lgam)
   !-----------------------------------------------------------------------
   !
   ! ...  This subroutine checks the projectability of each empty
   ! ...  wavefunction on the occupied manifold (for complex wfc)
   !
   USE io_global, ONLY: stdout
   USE kinds, ONLY: DP
   USE wavefunctions_module, ONLY: c0
   USE gvecw, ONLY: ngw
   USE electrons_base, ONLY: nbspx, nspin
   USE electrons_module, ONLY: nbsp_emp
   USE mp_global, ONLY: intra_pool_comm
   USE mp, ONLY: mp_sum
   USE reciprocal_vectors, ONLY: gstart
   !
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), INTENT(IN) :: c0_emp(:, :)
   LOGICAL, INTENT(IN) :: lgam
   !
   INTEGER :: m, n
   REAL(DP) :: proj, proj_tot
   COMPLEX(DP) :: g0comp(nbspx)
   !
   COMPLEX(DP), EXTERNAL :: ZDOTC
   !
   !
   WRITE (stdout, '(/, A)') " -----------------------------------------"
   WRITE (stdout, '(A)') " Projectability EMP states on OCC manifold"
   WRITE (stdout, '(A, /)') " -----------------------------------------"
   !
   !
   proj_tot = 0.D0
   !
   !
   DO m = 1, nbsp_emp
      !
      proj = 0.D0
      !
      DO n = 1, nbspx
         !
         proj = proj + ZDOTC(ngw, c0(:, n), 1, c0_emp(:, m), 1)
         !
      END DO
      !
      ! when gamma-trick is used ...
      !
      IF (lgam) THEN
         !
         ! ... account for G<0 vectors
         !
         proj = proj*2
         !
         IF (gstart == 2) THEN
            !
            ! ... and remove double counting for G=0 component
            !
            g0comp(:) = CONJG(c0(1, :))*c0_emp(1, m)
            proj = proj - SUM(g0comp(:))
            !
         END IF
         !
      END IF
      !
      CALL mp_sum(proj, intra_pool_comm)
      WRITE (stdout, 100) m, proj
      proj_tot = proj_tot + proj
      !
   END DO
   !
   !
   ! ... normalize proj_tot to 1
   !
   proj_tot = proj_tot/nbsp_emp
   WRITE (stdout, 101) proj_tot
   !
   !
100 FORMAT(4X, "orbital # ", I4, " : ", F12.8)
101 FORMAT(/, 6X, "Total projectability = ", F12.8,/)
   !
   !
END SUBROUTINE ortho_check_cmplx

