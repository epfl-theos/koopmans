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
subroutine sic_wkercorr(n,f,wsic,w2sic) 
  !---------------------------------------------------------------
  !   set up the orbital-dependent sic potential related to the kernel f_xc
  !
  use kinds, only : dp
  use radial_grids, only : ndmx
  use constants, only: e2, fpi
  use ld1inc, only : nspin, lsd, rel, nlcc, rhoc, grid, psi, rho, isw, fref, &
  ll, rhobarfact,nkscalfact, do_nkpz, do_wref, do_wxd
  use radial_grids, only: hartree
  use funct, only: dmxc_spin
  implicit none
  integer :: n
  real(DP):: f,wsic(ndmx,2),w2sic(ndmx)
  real(DP):: int_0_inf_dr
  !
  integer :: i
  integer :: is
  real(DP):: rh(2), rhc, vxc(2), vxc0(2), dvxc(2), dmuxc(2,2)
  real(DP):: rhoele(ndmx), rhobar(ndmx,2), w2cst, waux(ndmx), vh(ndmx)
  real(DP),parameter:: vanishing_rh=1.0e-12_dp
  logical:: do_w2
  !
  do_w2 = .true.
  nspin=1
  rhoc=0.0_dp
  if (lsd.eq.1) nspin=2
  !
  !    calculate electronic densities and electrostatic potential
  ! 
  rhoele=0.0_dp
  rhobar=0.0_dp
  vh=0.0_dp
  do i=1,grid%mesh
    rhoele(i)=psi(i,1,n)**2
    do is=1,nspin
      rhobar(i,is)=rho(i,is)
      if(is==isw(n)) rhobar(i,is)=rhobar(i,is)-f*rhoele(i)
    enddo
  enddo
  call hartree(0,2,grid%mesh,grid,rhoele,vh)
  vh(1:ndmx)=e2*vh(1:ndmx)
  !
  !    calculate the potential contribution related to the kernel  
  !
  rhc=0.0_dp
  rh=0.0_dp
  w2sic=0.0_dp
  wsic=0.0_dp
  waux=0.0_dp
  do i=1,grid%mesh
     do is=1,nspin
       rh(is)=rhobar(i,is)/grid%r2(i)/fpi*rhobarfact
       if(is==isw(n)) rh(is)=rh(is)+fref*rhoele(i)/grid%r2(i)/fpi
     enddo
     dmuxc=0.0_dp
     call dmxc_spin_ld1(rh(1),rh(2),dmuxc(1,1),dmuxc(1,2), &
                        dmuxc(2,1),dmuxc(2,2),vanishing_rh)
     dvxc(1)=dmuxc(1,isw(n))*rhoele(i)/grid%r2(i)/fpi*f
     dvxc(2)=dmuxc(2,isw(n))*rhoele(i)/grid%r2(i)/fpi*f
     w2sic(i)=dmuxc(isw(n),isw(n))*rhoele(i)/grid%r2(i)/fpi+vh(i)
     waux(i)=w2sic(i)*psi(i,1,n)**2
     do is=1,nspin
       rh(is)=rhobar(i,is)/grid%r2(i)/fpi*rhobarfact
     enddo
     vxc0=0.0_dp
     call vxc_t(rh,rhc,lsd,vxc0)
     do is=1,nspin
       rh(is)=rhobar(i,is)/grid%r2(i)/fpi*rhobarfact
       if(is==isw(n)) rh(is)=rh(is)+f*rhoele(i)/grid%r2(i)/fpi
     enddo
     vxc=0.0_dp
     call vxc_t(rh,rhc,lsd,vxc)
     do is=1,nspin
       if(dvxc(is).ne.0.0_dp) &
         wsic(i,is) = rhobarfact*(vxc0(is)+dvxc(is)-vxc(is))
     enddo
  enddo
  !
  !     add contribution for NK on top of PZ
  !
  if(do_nkpz) then
    do i=1,grid%mesh
      rh=0.0_dp
      do is=1,nspin
        if(is==isw(n)) rh(is)=fref*rhoele(i)/grid%r2(i)/fpi
      enddo
      dmuxc=0.0_dp
      call dmxc_spin_ld1(rh(1),rh(2),dmuxc(1,1),dmuxc(1,2), &
                         dmuxc(2,1),dmuxc(2,2),vanishing_rh)
      w2sic(i)=w2sic(i)-dmuxc(isw(n),isw(n))*rhoele(i)/grid%r2(i)/fpi-vh(i)
      waux(i)=w2sic(i)*psi(i,1,n)**2
    enddo
  endif
  if( do_w2 ) then
    w2cst = int_0_inf_dr(waux,grid,grid%mesh,2*ll(n)+2)
    do i=1,grid%mesh
       w2sic(i) = fref*(w2sic(i)-w2cst)
    enddo
  endif
  !
  wsic=nkscalfact*wsic
  w2sic=nkscalfact*w2sic
  !
  if(.not.do_wref) w2sic=0.0_dp
  if(.not.do_wxd) wsic=0.0_dp
  !
  return
  !
end subroutine sic_wkercorr
!---------------------------------------------------------------
      subroutine dmxc_spin_ld1(rhoup, rhodw,dmuxc_uu, &
                              dmuxc_ud, dmuxc_du, dmuxc_dd, small)
!---------------------------------------------------------------
!
! ... derivative of the xc potential with respect to the local density
!
      USE kinds,                ONLY : DP
      USE funct,                ONLY : xc_spin, get_iexch, get_icorr
      implicit none
      !
      real(dp), intent(in) :: rhoup, rhodw, small
      real(dp), intent(out) :: dmuxc_uu, dmuxc_ud, dmuxc_du, dmuxc_dd
      !
      real(dp) :: rhotot, rs, zeta, fz, fz1, fz2, ex, vx, ecu, ecp, vcu, &
           vcp, dmcu, dmcp, aa, bb, cc, dr, dz, ec, vxupm, vxdwm, vcupm, &
           vcdwm, rho, vxupp, vxdwp, vcupp, vcdwp, dzm, dzp
      real(dp), external :: dpz, dpz_polarized
      integer :: iflg
      !
      real(dp), parameter :: e2 = 2.0_dp, &
           pi34=0.75_dp/3.141592653589793_dp, third=1.0_dp/3.0_dp, &
           p43=4.0_dp/3.0_dp, p49=4.0_dp/ 9.0_dp, m23=-2.0_dp/3.0_dp
      real(dp) :: f, alpha
      parameter (f=-1.10783814957303361d0, alpha=2.0d0/3.0d0)
      !
      ! ... initialize variable
      !
      dmuxc_uu=0.0_dp
      dmuxc_du=0.0_dp
      dmuxc_ud=0.0_dp
      dmuxc_dd=0.0_dp
      !
      rhotot=rhoup+rhodw
      if(rhotot<small) return
      !
      zeta=(rhoup-rhodw)/rhotot
      if(abs(zeta)>1.0_dp) zeta=sign(1.0_dp,zeta)
      !
      ! ... calculate exchange contribution (analytical)
      !
      if (get_iexch().ne.1) goto 100
      !
      if(rhoup>small) then
        rs=(pi34/(2.0_dp*rhoup))**third
        call slater(rs,ex,vx)
        dmuxc_uu=vx/(3.0_dp*rhoup)
      endif
      if(rhodw>small) then
        rs=(pi34/(2.0_dp*rhodw))**third
        call slater(rs,ex,vx)
        dmuxc_dd=vx/(3.0_dp*rhodw)
      endif
      !
100   continue
      !
      ! ... calculate correlation contribution (numerical)
      !
      if(get_icorr().ne.1) goto 200
      !
      dr=min(1.e-6_dp,1.e-4_dp*rhotot)
      call xc_spin(rhotot-dr,zeta,ex,ec,vxupm,vxdwm,vcupm,vcdwm)
      call xc_spin(rhotot+dr,zeta,ex,ec,vxupp,vxdwp,vcupp,vcdwp)
      dmuxc_uu=dmuxc_uu+(vcupp-vcupm)/(2.0_dp*dr)
      dmuxc_ud=dmuxc_ud+(vcupp-vcupm)/(2.0_dp*dr)
      dmuxc_dd=dmuxc_dd+(vcdwp-vcdwm)/(2.0_dp*dr)
      dmuxc_du=dmuxc_du+(vcdwp-vcdwm)/(2.0_dp*dr)
      !
      dz=1.e-6_dp
      dzp=min(1.0,zeta+dz)-zeta
      dzm=-max(-1.0,zeta-dz)+zeta
      call xc_spin(rhotot,zeta-dzm,ex,ec,vxupm,vxdwm,vcupm,vcdwm)
      call xc_spin(rhotot,zeta+dzp,ex,ec,vxupp,vxdwp,vcupp,vcdwp)
      dmuxc_uu=dmuxc_uu+(vcupp-vcupm)*(1.0_dp-zeta)/rhotot/(dzp+dzm)
      dmuxc_ud=dmuxc_ud-(vcupp-vcupm)*(1.0_dp+zeta)/rhotot/(dzp+dzm)
      dmuxc_du=dmuxc_du+(vcdwp-vcdwm)*(1.0_dp-zeta)/rhotot/(dzp+dzm)
      dmuxc_dd=dmuxc_dd-(vcdwp-vcdwm)*(1.0_dp+zeta)/rhotot/(dzp+dzm)
      !
200   continue
      !
      ! ... convert to rydberg units
      !
      dmuxc_uu=dmuxc_uu*e2
      dmuxc_ud=dmuxc_ud*e2
      dmuxc_du=dmuxc_du*e2
      dmuxc_dd=dmuxc_dd*e2
      !
      return
      !
!---------------------------------------------------------------
      end subroutine dmxc_spin_ld1
!---------------------------------------------------------------
