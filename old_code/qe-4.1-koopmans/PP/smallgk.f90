!
! Copyright (C) 2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
SUBROUTINE smallgk (xk, at, bg, s, ftau, t_rev, sname, nsym, sk, ftauk, gk, &
                    t_revk, snamek, nsymk)
!-----------------------------------------------------------------------
!
! This routine selects, among the symmetry matrices of the point group
! of a crystal, the symmetry operations which leave k unchanged.
!
!
USE kinds, ONLY : DP
IMPLICIT NONE
CHARACTER(LEN=45) :: snamek(48), sname(48)

REAL(DP) :: bg (3, 3), at (3, 3), xk (3)
! input: the reciprocal lattice vectors
! input: the direct lattice vectors
! input: the k point of the crystal

INTEGER :: s (3, 3, 48), ftau(3,48), t_rev(48), nsym, sk (3, 3, 48), &
           ftauk(3,48), t_revk(48), gk(3,48), nsymk
! input: the symmetry matrices
! input: fractional translation associated to each rotation
! input: possible time reversal associated to the rotation
! input: dimension of the point group
! output: the symmetry matrices of the small group of k
! output: the fract. trans. associated to the operations of the small group of k
! output: the time reversal associated to the operations of the small group of k
! output: the G vector which connects k and the rotated k.

  REAL(DP) :: ak (3), rak (3), zero (3)
  ! k vector in crystal basis
  ! the rotated of the k vector
  ! the zero vector

  INTEGER :: isym, ipol, jpol
  ! counter on symmetry operations
  ! counter on polarizations
  ! counter on polarizations

  LOGICAL :: eqvect
  ! logical function, check if two vectors are equal
  !
  !  Set to zero some variables and transform xq to the crystal basis
  !
  zero = 0.d0
  ak = xk
  CALL cryst_to_cart (1, ak, at, - 1)
  !
  !   test all symmetries to see if the operation S sends k in k+G ...
  !
  nsymk = 0
  DO isym = 1, nsym
     rak = 0.d0
     DO ipol = 1, 3
        DO jpol = 1, 3
           rak (ipol) = rak (ipol) + DBLE (s (ipol, jpol, isym) ) * &
                ak (jpol)
        ENDDO
     ENDDO
     IF ((t_rev(isym)==0 .AND. eqvect(rak, ak, zero)) .OR. &
         (t_rev(isym)==1 .AND. eqvect(rak, -ak, zero)) ) THEN
        nsymk=nsymk+1
        sk(:,:,nsymk)=s(:,:,isym)
        ftauk(:,nsymk)=ftau(:,isym)
        snamek(nsymk)=sname(isym)
        t_revk(nsymk)=t_rev(isym)
        IF (t_rev(isym)==0) THEN
           gk(:,nsymk)=NINT(rak(:)-ak(:))
        ELSEIF (t_rev(isym)==1) THEN
           gk(:,nsymk)=NINT(rak(:)+ak(:))
        ELSE
           CALL errore('smallgk','wrong t_rev',1)
        ENDIF
     END IF
  ENDDO
  !
  RETURN
END SUBROUTINE smallgk

