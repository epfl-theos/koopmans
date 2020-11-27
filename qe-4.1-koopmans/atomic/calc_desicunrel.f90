!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine calc_desicunrel(desicunrel,nunrel,vsicunrel)
  !---------------------------------------------------------------
  !
  ! calculate self-interaction energy (allows fractional occupations)
  !
  use kinds, only: dp
  use radial_grids, only: ndmx
  use ld1inc, only: grid, psi, nwf,ll
 
  implicit none
  integer:: i
  integer,  intent(in) :: nunrel
  real(DP), intent(out) :: desicunrel
  real(DP) :: work(ndmx), vsicunrel(ndmx)
  real(DP), external :: int_0_inf_dr
  !
  work=0.0_dp
  do i=1,grid%mesh
    work(i)=vsicunrel(i)*psi(i,1,nunrel)**2
  enddo
  desicunrel = -int_0_inf_dr(work,grid,grid%mesh,2*(ll(nunrel)+1))
  !
  return
end subroutine calc_desicunrel
