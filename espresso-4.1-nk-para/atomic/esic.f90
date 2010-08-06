!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Non-Koopmans method 
! Developed and implemented by I. Dabo (Universite Paris-Est, Ecole des Ponts, ParisTech)
!
!---------------------------------------------------------------
subroutine esic
  !---------------------------------------------------------------
  !
  ! calculate self-interaction energy (allows fractional occupations)
  !
  use kinds, only: dp
  use radial_grids, only: ndmx
  use ld1inc, only: grid, etot, ehrt, ecxc, ekin, rel, dhrsic, dxcsic, &
                    nwf, ll,  oc, psi, fref, enne, isw, vsic, wsic, w2sic, &
                    wsictot,nkscalfact,do_nkpz
 
  implicit none
  ! output
  ! local
  integer:: n, i
  real(DP) :: int_0_inf_dr,deksic  
  real(DP) :: work1(ndmx),work(ndmx),v(ndmx),f(nwf)
  external int_0_inf_dr
  !
  deksic = 0.0_DP
  dhrsic = 0.0_DP
  dxcsic = 0.0_DP
  wsictot=0.0_dp
  vsic=0.0_dp
  wsic=0.0_dp
  w2sic=0.0_dp
  f=0.0_dp
  do n=1,nwf
     f(n)=min(oc(n),1.0_dp)
     call sic_wkercorr(n,f(n),wsic(1,1,n),w2sic(1,n))
     if(oc(n)>1.0_dp) then
       wsictot(:,:) = wsictot(:,:)+oc(n)*wsic(:,:,n)
     else
       wsictot(:,:) = wsictot(:,:)+wsic(:,:,n)
     endif
  enddo
  do n=1,nwf
     call sic_correction(n,f(n),v,vsic(1,n),work1)
     vsic(:,n) = vsic(:,n)-(wsictot(:,isw(n))-wsic(:,isw(n),n))-w2sic(:,n)
     !
     work=0.0_dp
     do i=1,grid%mesh
        v(i)=v(i)*psi(i,1,n)**2
        work(i)=vsic(i,n)*psi(i,1,n)**2
     enddo
     if(oc(n)>1.0_dp) then
       deksic = deksic + oc(n)*int_0_inf_dr(work,grid,grid%mesh,2*(ll(n)+1))
       dhrsic = dhrsic + oc(n)*f(n)*(2.0_dp*fref-f(n)) &
                       * 0.5_dp*int_0_inf_dr(v,grid,grid%mesh,2)*nkscalfact
       if(do_nkpz) dhrsic=dhrsic-oc(n)*f(n)*fref &
                         *int_0_inf_dr(v,grid,grid%mesh,2)*nkscalfact
       dxcsic = dxcsic - oc(n)*int_0_inf_dr(work1,grid,grid%mesh,2)
     else
       deksic = deksic + f(n)*int_0_inf_dr(work,grid,grid%mesh,2*(ll(n)+1))
       dhrsic = dhrsic + f(n)*(2.0_dp*fref-f(n)) &
                       * 0.5_dp*int_0_inf_dr(v,grid,grid%mesh,2)*nkscalfact
       if(do_nkpz) dhrsic=dhrsic-f(n)*fref &
                         *int_0_inf_dr(v,grid,grid%mesh,2)*nkscalfact
       dxcsic = dxcsic - int_0_inf_dr(work1,grid,grid%mesh,2)
     endif
  enddo
  !
  ekin=ekin+deksic
  ehrt=ehrt
  etot=etot+dhrsic+dxcsic+deksic
  !
  return
end subroutine esic
