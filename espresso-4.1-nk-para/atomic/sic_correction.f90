!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Non-Koopmans method developed and implemented by I. Dabo (Ecole des Pont, ParisTech)
!
!---------------------------------------------------------------
subroutine sic_correction(n,f,vhn1,vhn2,egc) 
  !---------------------------------------------------------------
  !   set up the orbital-dependent selfconsistent potential generated
  !   by the n-th wavefunction - for self-interaction correction
  !
  use kinds, only : dp
  use radial_grids, only : ndmx
  use constants, only: e2, fpi
  use ld1inc, only : nspin, lsd, rel, nlcc, rhoc, grid, psi, rho, fref, &
                     isw, rhobarfact, do_nkmix, nkmixfact, nkscalfact
  use funct, only: dft_is_gradient
  use radial_grids, only: hartree
  implicit none
  integer :: n
  real(DP):: f,vhn1(ndmx),vhn2(ndmx), egc(ndmx)
  !
  integer :: i, is
  real(DP):: rh(2), rhc, exc_t, vxc(2), vxcref(2), vxc0(2)
  real(DP):: vgc(ndmx,2),  egc0(ndmx), rhoele(ndmx), rhobar(ndmx,2)
  logical :: gga
  vhn1=0.0_dp
  vhn2=0.0_dp
  rhc=0.0_dp
  egc=0.0_dp
  gga=dft_is_gradient()
  nspin=1
  if (lsd.eq.1) nspin=2
  !
  !   compute hartree potential with the charge of orbital n
  !
  rhoele=0.0_dp
  do i=1,grid%mesh
    rhoele(i)=psi(i,1,n)**2
  enddo
  call hartree(0,2,grid%mesh,grid,rhoele,vhn1)
  !
  !
  rhobar=0.0_dp
  do i=1,grid%mesh
    do is=1,nspin
      rhobar(i,is)=rho(i,is)
      if(is==isw(n)) rhobar(i,is)=rhobar(i,is)-f*rhoele(i)
    enddo
  enddo
  !
  !    add exchange and correlation contribution: LDA or LSDA terms
  !
  do i=1,grid%mesh
     !
     vhn1(i) = e2*vhn1(i)
     !
     rh=0.0_dp
     do is=1,nspin
       rh(is)=rhobar(i,is)*rhobarfact
       if(is==isw(n)) rh(is)=rh(is)+fref*rhoele(i)
       rh(is)=rh(is)/grid%r2(i)/fpi
     enddo
     vxcref=0.0_dp
     call vxc_t(rh,rhc,lsd,vxcref)
     egc(i)=-vxcref(isw(n))*rhoele(i)*f
     !
     rh=0.0_dp
     do is=1,nspin
       rh(is)=rhobar(i,is)*rhobarfact
       if(is==isw(n)) rh(is)=rh(is)+f*rhoele(i)
       rh(is)=rh(is)/grid%r2(i)/fpi
     enddo
     vxc=0.0_dp
     call vxc_t(rh,rhc,lsd,vxc)
     egc(i)=egc(i)+exc_t(rh,rhc,lsd)*(rh(1)+rh(2))*grid%r2(i)*fpi
     !
     rh=0.0_dp
     do is=1,nspin
       rh(is)=rhobar(i,is)*rhobarfact
       rh(is)=rh(is)/grid%r2(i)/fpi
     enddo
     egc(i)=egc(i)-exc_t(rh,rhc,lsd)*(rh(1)+rh(2))*grid%r2(i)*fpi
     !
     vhn2(i)= vhn1(i)*(f-fref)+vxc(isw(n))-vxcref(isw(n))
     !
     if( do_nkmix ) then
       rh=0.0_dp
       do is=1,nspin
         if(is==isw(n)) rh(is)=f*rhoele(i)
         rh(is)=rh(is)/grid%r2(i)/fpi
       enddo
       vxc=0.0_dp
       call vxc_t(rh,rhc,lsd,vxc)
       egc(i)=(1.0-nkmixfact)*egc(i)+nkmixfact*exc_t(rh,rhc,lsd)*(rh(1)+rh(2))*grid%r2(i)*fpi
       vhn2(i)=(1.0-nkmixfact)*vhn2(i)+nkmixfact*(vhn1(i)*f+vxc(isw(n)))
     end if
     !
     egc(i)=egc(i)*nkscalfact
     vhn2(i)=vhn2(i)*nkscalfact
     !
  enddo
  !
  return
  !
end subroutine sic_correction
