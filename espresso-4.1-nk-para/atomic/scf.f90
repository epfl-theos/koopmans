!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
subroutine scf(ic)
  !---------------------------------------------------------------
  !
  !   this routine performs the atomic self-consistent procedure
  !   self-interaction-correction allowed
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

  integer, intent(in) :: ic

  logical:: conv
  integer:: nerr, nstop, n, i, id, nin, mch
  real(DP) ::  vnew(ndmx,2), rhoc1(ndmx), ze2, f(nwf)
  integer, parameter :: maxter=200
  real(DP), parameter :: thresh=1.0e-10_dp
  !
  ! 
  ze2 = - zed * e2
  rhoc1=0.0_dp
  id=3
  if (.not.frozen_core.or.ic==1) psi=0.0_dp
  !!!
  if (isic /= 0 .and. relpert)  id=1 ! 
  !!!
  do iter=1,maxter
     nerr=0
     vnew=vpot
     do n=1,nwf
        !if (oc(n) >= 0.0_dp) then
           if (ic==1.or..not.frozen_core.or..not.core_state(n)) then
              if (isic /= 0 .and. iter > 1) vnew(:,isw(n))=vpot(:,isw(n))-vsic(:,n)
              if (rel == 0) then
                 call ascheq (nn(n),ll(n),enl(n),grid%mesh,grid,vnew(1,isw(n)),&
                      ze2,thresh,psi(1,1,n),nstop)
              elseif (rel == 1) then
                 call lschps (1,zed,grid,nin,mch,nn(n),ll(n),enl(n),&
                             psi(1,1,n),vnew(1,isw(n)),nstop)
                 if (nstop>0.and.oc(n)<1.e-10_DP) nstop=0
              elseif (rel == 2) then
                 call dirsol (ndmx,grid%mesh,nn(n),ll(n),jj(n),iter,enl(n), &
                      thresh,grid,psi(1,1,n),vnew(1,isw(n)),nstop)
              else
                 call errore('scf','relativistic not programmed',1)
              endif
              !      write(6,*) el(n),enl(n)
              ! if (nstop /= 0) write(6,'(4i6)') iter,nn(n),ll(n),nstop
              nerr=nerr+nstop
          !endif
        else
           enl(n)=0.0_dp
           psi(:,:,n)=0.0_dp
        endif
     enddo
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
          lsd,.false.,latt,enne,rhoc1,rho,vh,vnew,1)
     !
     ! calculate SIC correction potential (if present)
     !
     if (isic /= 0) then
        vsicnew=0.0_dp
        wsictot=0.0_dp
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
            call sic_correction(n,f(n),vhn1,vsicnew,egc)
            vsicnew(:) = vsicnew(:)-(wsictot(:,isw(n))-wsic(:,isw(n),n))-w2sic(:,n)
            !
            ! use simple mixing for SIC correction
            !
            vsic(:,n) = (1.0_dp-beta)*vsic(:,n)+beta*vsicnew(:)
        enddo
     endif
     !
     ! mix old and new potential
     !
     call vpack(grid%mesh,ndmx,nspin,vnew,vpot,1)
     call dmixp(grid%mesh*nspin,vnew,vpot,beta,tr2,iter,id,eps0,conv,maxter)
     call vpack(grid%mesh,ndmx,nspin,vnew,vpot,-1)
!        write(6,*) iter, eps0
     !
     if (conv) then
        if (nerr /= 0) call infomsg ('scf','errors in KS equations')
        goto 45
     endif
  enddo
  call infomsg('scf','warning: convergence not achieved')
45 continue
   return

end subroutine scf
