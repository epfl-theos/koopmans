!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine calc_desic(desic)
  !---------------------------------------------------------------
  !
  ! calculate self-interaction energy (allows fractional occupations)
  !
  use kinds, only: dp
  use radial_grids, only: ndmx
  use ld1inc, only: grid, psi, vsic, nwf,ll
 
  implicit none
  integer:: n, i
  real(DP), intent(out) :: desic(nwf)
  real(DP) :: work(ndmx)
  real(DP), external :: int_0_inf_dr
  !
  do n = 1, nwf
     work=0.0_dp
     do i=1,grid%mesh
        work(i)=vsic(i,n)*psi(i,1,n)**2
     enddo
     desic(n) = -int_0_inf_dr(work,grid,grid%mesh,2*(ll(n)+1))
  enddo
  !
  return
end subroutine calc_desic
