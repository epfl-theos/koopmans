!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Non-Koopmans method 
! developed and implemented by I. Dabo (Universite Paris-Est, Ecole des Ponts, ParisTech)
!
!---------------------------------------------------------------
subroutine esic_unrel(enlunrel,nunrel,funrel,ekinnl)
  !---------------------------------------------------------------
  !
  use kinds, only: dp
  use radial_grids, only: ndmx
  use ld1inc, only: grid, etot, ehrt, ecxc, ekin, rel, dhrsic, dxcsic, &
                    nwf, ll,  oc, psi, wsic, w2sic, wsictot, fref, &
                    enne, isw,nkscalfact,do_nkpz
 
  implicit none
  ! output
  ! local
  integer:: n, i
  real(DP) :: int_0_inf_dr
  real(DP) :: work1(ndmx),vsic(ndmx),vsicunrel(ndmx),wsicunrel(ndmx,2), &
              w2sicunrel(ndmx),v(ndmx),f(nwf),nfill(nwf),desicunrel
  real(DP) :: ekinnl(nwf)
  real(DP) :: work1fr(ndmx)
  real(DP), intent(in) :: funrel
  real(DP), intent(out) :: enlunrel
  integer, intent(in) :: nunrel
  external int_0_inf_dr
  !
  dhrsic = 0.0_DP
  dxcsic = 0.0_DP
  wsictot=0.0_dp
  wsic=0.0_dp
  w2sic=0.0_dp
  wsicunrel=0.0_dp
  w2sicunrel=0.0_dp
  oc(nunrel)=oc(nunrel)-funrel
  f=0.0_dp
  !
  do n=1,nwf
     if(oc(n)>0.0_dp) f(n)=1.0_dp
     if(f(n)>0.0_dp) then
        call sic_wkercorr(n,f(n),wsic(1,1,n),w2sic(1,n))
        wsictot(:,:) = wsictot(:,:)+oc(n)*wsic(:,:,n)
        ! the contribution from the orbital nunrel is not included
     endif
  enddo
  call sic_wkercorr(nunrel,funrel,wsicunrel,w2sicunrel)
  !
  do n=1,nwf
     v=0.0_dp
     work1=0.0_dp
     vsic=0.0_dp
     call sic_correction(n,f(n),v,vsic,work1)
     do i=1,grid%mesh
        v(i)=v(i)*psi(i,1,n)**2
     enddo
     if(oc(n)>1.0_dp) then
       dhrsic = dhrsic + oc(n)*f(n)*(2.0_dp*fref-f(n)) &
                       * 0.5_dp*int_0_inf_dr(v,grid,grid%mesh,2)*nkscalfact
       if(do_nkpz) dhrsic = dhrsic - oc(n)*f(n)*fref &
                       *int_0_inf_dr(v,grid,grid%mesh,2)*nkscalfact
       dxcsic = dxcsic - oc(n)*int_0_inf_dr(work1,grid,grid%mesh,2)
     else
       dhrsic = dhrsic + f(n)*(2.0_dp*fref-f(n)) &
                       * 0.5_dp*int_0_inf_dr(v,grid,grid%mesh,2)*nkscalfact
       if(do_nkpz) dhrsic = dhrsic - f(n)*fref &
                       *int_0_inf_dr(v,grid,grid%mesh,2)*nkscalfact
       dxcsic = dxcsic - int_0_inf_dr(work1,grid,grid%mesh,2)
     endif
  enddo
  !
  call sic_correction(nunrel,funrel,v,vsicunrel,work1)
  vsicunrel(:) = vsicunrel(:)-wsictot(:,isw(nunrel))-w2sicunrel(:)
  !
  do i=1,grid%mesh
    v(i)=v(i)*psi(i,1,nunrel)**2
  enddo
  dhrsic = dhrsic + funrel*(2.0_dp*fref-funrel) &
                  * 0.5_dp*int_0_inf_dr(v,grid,grid%mesh,2)*nkscalfact
  dxcsic = dxcsic - int_0_inf_dr(work1,grid,grid%mesh,2)
  !
  etot=etot+dhrsic+dxcsic
  !
  call calc_desicunrel(desicunrel,nunrel,vsicunrel)
  enlunrel=enlunrel+desicunrel
  !
  oc(nunrel)=oc(nunrel)+funrel
  !
  return
  !
end subroutine esic_unrel
