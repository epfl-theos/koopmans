!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
FUNCTION exx_divergence()
  !-----------------------------------------------------------------------
  !
  ! ...  This function calculates the G=0 term of the Coulomb potential
  ! ...  in presence  of a unitary charge (q=1). Used for counter-charge
  ! ...  corrections.
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : fpi, e2
  USE cell_base,          ONLY : bg, at, alat, omega
  USE reciprocal_vectors, ONLY : gx
  USE gvecp,              ONLY : ngm
  USE gvecw,              ONLY : ecutw
  USE io_global,          ONLY : stdout
  USE control_flags,      ONLY : gamma_only
  USE mp_global,          ONLY : intra_pool_comm
  USE mp,                 ONLY : mp_sum
  !
  !
  IMPLICIT NONE
  !
  REAL(DP) :: exx_divergence
  !
  REAL(DP) :: yukawa = 0.d0
  INTEGER :: nq1=1, nq2=1, nq3=1   ! integers defining the X integration mesh
  INTEGER :: nqs                   ! number of points in the q-grid
  LOGICAL :: x_gamma_extrapolation=.false.
  LOGICAL :: on_double_grid=.false.
  REAL(DP) :: grid_factor = 8.d0/7.d0
  REAL(DP) :: eps=1.d-6
  INTEGER :: iq1,iq2,iq3, ig
  REAL(DP) :: div, dq1, dq2, dq3, xq(3), q_, qq, tpiba2, alpha, x, q(3)
  INTEGER :: nqq, iq
  REAL(DP) :: aa, dq
  !
  !
  CALL start_clock( 'exx_div' )
  !
  WRITE( stdout, '(/,A,3I4)' ) " EXX : q-grid dimensions are ", nq1, nq2, nq3
  WRITE( stdout, '(A,5X,L)' )  " EXX : Gamma Extrapolation", x_gamma_extrapolation
  !
  IF ( x_gamma_extrapolation ) THEN
     WRITE( stdout, '(A)' ) " EXX : q->0 dealt with 8/7 -1/7 trick"
     grid_factor = 8.d0 / 7.d0
  ELSE
     WRITE( stdout, '(A)' ) " EXX : q->0 term not estimated"
     grid_factor = 1.d0
  ENDIF
  !
  nqs = nq1 * nq2 * nq3
  tpiba2 = ( fpi / 2.d0 / alat ) ** 2
  alpha  = 10.d0 * tpiba2 / ecutw
  !
  dq1= 1.d0/DBLE(nq1)
  dq2= 1.d0/DBLE(nq2) 
  dq3= 1.d0/DBLE(nq3)
  !
  div = 0.d0
  !
  DO iq1 = 1, nq1
    DO iq2 = 1, nq2
      DO iq3 = 1, nq3
        !
        xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                bg(:,2) * (iq2-1) * dq2 + &
                bg(:,3) * (iq3-1) * dq3
        !
        DO ig = 1, ngm
          !
          q(1) = xq(1) + gx(1,ig)
          q(2) = xq(2) + gx(2,ig)
          q(3) = xq(3) + gx(3,ig)
          qq = q(1)*q(1) + q(2)*q(2) + q(3)*q(3)
          !
          IF ( x_gamma_extrapolation ) THEN
            !
            on_double_grid = .true.
            x = 0.5d0*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
            on_double_grid = on_double_grid .AND. (abs(x-nint(x))<eps)
            x = 0.5d0*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
            on_double_grid = on_double_grid .AND. (abs(x-nint(x))<eps)
            x = 0.5d0*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
            on_double_grid = on_double_grid .AND. (abs(x-nint(x))<eps)
            !
          ENDIF
          !
          IF ( .NOT. on_double_grid ) THEN
            IF ( qq > 1.d-8 ) THEN
              div = div + EXP( - alpha * qq ) / ( qq + yukawa / tpiba2 ) &
                                                  * grid_factor
            ENDIF
          ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  CALL mp_sum( div, intra_pool_comm )
  !
  IF ( gamma_only ) div = 2.d0 * div
  !
  IF ( .NOT. x_gamma_extrapolation ) THEN
    IF ( yukawa < 1.d-8) THEN
      div = div - alpha
    ELSE
      div = div + tpiba2 / yukawa
    ENDIF
  ENDIF
  !
  div = div * e2 * fpi / tpiba2 / nqs
  !
  alpha = alpha / tpiba2
  !
  nqq = 100000
  dq = 5.0d0 / SQRT( alpha ) / nqq
  aa = 0.d0
  !
  DO iq = 0, nqq
    !
    q_ = dq * ( iq + 0.5d0 )
    qq = q_ * q_
    aa = aa - EXP( - alpha * qq ) * yukawa / ( qq + yukawa ) * dq
    !
  ENDDO
  !
  aa = aa * 8.d0 / fpi
  aa = aa + 1.d0 / SQRT( alpha * 0.25d0 * fpi )
  !  
  div = div - e2 * omega * aa
  !
  ! RdG: not clear what is this (absent in recent QE)
  !div = div - e2 * omega / SQRT( alpha * 0.25d0 * fpi )
  !
  exx_divergence = div * nqs
  !
  WRITE( stdout, '(A,1ES15.5,/)' ) " EXX : Coulomb G0 ", exx_divergence
  !
  CALL stop_clock( 'exx_div' )
  !CALL print_clock( 'exx_div' )
  !
  RETURN
  !
  !
END FUNCTION exx_divergence 
