!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine symz (phi, nsym, s, nat, irt)
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds, only : DP
  implicit none

  integer :: nsym, s (3, 3, 48), nat, irt (48, nat)
  ! input: the number of symmetries
  ! input: the rotation matrix
  ! input: the number of atoms
  ! input: correspondence between rotated atoms

  real(DP) :: phi (3, 3, nat)
  ! matrix to symmetrize

  integer :: isym, i, j, k, l, na, sna
  ! counter on symmetries
  ! counter on points
  ! counter on atoms
  ! the rotated atom

  real(DP) :: work (3, 3, nat)
  ! auxiliary space
  !
  if (nsym == 1) return
  work = 0.d0      
  !
  do na = 1, nat
     do isym = 1, nsym
        sna = irt (isym, na)
        do i = 1, 3
           do j = 1, 3
              do k = 1, 3
                 do l = 1, 3
                    work (i, j, na) = work (i, j, na) + &
                       s (i, k, isym) * s (j, l, isym) * phi (k, l, sna)
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  phi(:, :, :) = work (:, :, :) / DBLE (nsym)
  !
  return
end subroutine symz
