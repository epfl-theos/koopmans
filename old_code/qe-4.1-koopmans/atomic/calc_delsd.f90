!
! Copyright (C) 2004-2007 QUANTUM-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine calc_delsd(zed,grid,vxt,vh,vxc,nwf,nspin,delsd)   
  !---------------------------------------------------------------
  !
  use kinds, only : DP
  use radial_grids, only: ndmx, radial_grid_type
  use funct, only: get_iexch
  use ld1inc, only : isw,psi
  implicit none
  integer, intent(in) :: nwf, nspin
  type(radial_grid_type),intent(in)::grid
  real(DP), intent(in)  :: zed
  real(DP), intent(in)  :: vxt(ndmx), vh(ndmx), vxc(ndmx,2)
  real(DP), intent(out) :: delsd(nwf)
  real(DP),allocatable :: f1(:), f2(:), f3(:), f4(:), f5(:)
  real(DP), external :: int_0_inf_dr
  real(DP):: rhoele
  integer:: i,n,is,ierr

  allocate(f1(grid%mesh),stat=ierr)
  allocate(f2(grid%mesh),stat=ierr)
  allocate(f3(grid%mesh),stat=ierr)
  allocate(f4(grid%mesh),stat=ierr)
  allocate(f5(grid%mesh),stat=ierr)

  do n=1,nwf
    do i=1,grid%mesh
!    
     rhoele=psi(i,1,n)**2
!
!   The integral for the energy due to the interaction with nuclei
!
     f1(i)=-2.0_DP*zed/grid%r(i) * rhoele
!
!   The integral for the Hartree energy
!
     f2(i)= vh (i) * rhoele
!
!   The integral for the interaction with an external potential
!
     f4(i)= vxt(i) * rhoele
!
!   The integral to be subtracted to the sum of the eigenvalues to
!   get the kinetic energy
!
     f5(i) =vxc(i,isw(n))*rhoele+f1(i)+f2(i)+f4(i)
    enddo
!  The kinetic energy is the sum of the eigenvalues plus the f5 integral
!
   delsd(n) = int_0_inf_dr(f5,grid,grid%mesh,1)
  enddo

  deallocate(f5)
  deallocate(f4)
  deallocate(f3)
  deallocate(f2)
  deallocate(f1)

  return
end subroutine calc_delsd
