!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc &
            ( npwx, npw, nstart, gstart, nbnd, psi, npol, overlap, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Driver routine (maybe it should be an interface) for
  ! ... Hamiltonian diagonalization in the subspace spanned
  ! ... by nstart states psi ( atomic or random wavefunctions ).
  ! ... Produces on output nbnd eigenvectors ( nbnd <= nstart ) in evc.
  ! ... Calls h_psi, s_psi to calculate H|psi> ans S|psi>
  ! ... It only uses an auxiliary array of the same size as psi.
  !
  USE kinds,         ONLY : DP
  USE control_flags, ONLY : use_para_diag, gamma_only
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd, gstart, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! first G with nonzero norm
    ! number of spin polarizations
  LOGICAL, INTENT(IN) :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nstart), evc(npwx*npol,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP), INTENT(OUT) :: e(nbnd)
    ! eigenvalues
  !
  CALL start_clock( 'wfcrot' )
  !
  IF( use_para_diag ) THEN
     !
     ! use data distributed subroutine
     !
     IF ( gamma_only ) THEN
        !
        CALL protate_wfc_gamma &
            ( npwx, npw, nstart, gstart, nbnd, psi, overlap, evc, e )
        !
     ELSE
        !
        CALL protate_wfc_k &
            ( npwx, npw, nstart, nbnd, npol, psi, overlap, evc, e )
        !
     END IF
     !
  ELSE
     !
     ! use serial subroutines
     !
     IF ( gamma_only ) THEN
        !
        CALL rotate_wfc_gamma &
            ( npwx, npw, nstart, gstart, nbnd, psi, overlap, evc, e )
        !
     ELSE
        !
        CALL rotate_wfc_k &
            ( npwx, npw, nstart, nbnd, npol, psi, overlap, evc, e )
        !
     END IF
     !
  END IF
  !
  CALL stop_clock( 'wfcrot' )
  !
END SUBROUTINE rotate_wfc
