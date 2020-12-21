!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
subroutine nscf
  !---------------------------------------------------------------
  !
  use kinds, only : dp
  use radial_grids, only : ndmx
  use constants, only: e2
  use ld1inc, only : grid, zed, psi, isic, vpot, vh, vxt, rho, iter, &
                     lsd, rel, latt, enne, beta, nspin, tr2, eps0, &
                     nwf, nn, ll, jj, enl, oc, isw, core_state, frozen_core, &
                     vsic, vsicnew, vhn1, egc, relpert, wsic, wsictot, &
                     w2sic
  implicit none
  integer:: n, i, is
  real(DP) ::  rhoc1(ndmx), ze2, f(nwf)
  !
  ! 
  ze2 = - zed * e2
  rhoc1=0.0_dp
     !
     ! calculate charge density (spherical approximation)
     !
     rho=0.0_dp
     do n=1,nwf
        do i=1,grid%mesh
           rho(i,isw(n))=rho(i,isw(n))+oc(n)*(psi(i,1,n)**2+psi(i,2,n)**2)
        enddo
     enddo
     !
     ! calculate new potential
     !
     call new_potential(ndmx,grid%mesh,grid,zed,vxt,&
          lsd,.false.,latt,enne,rhoc1,rho,vh,vpot,1)
     !
end subroutine nscf
