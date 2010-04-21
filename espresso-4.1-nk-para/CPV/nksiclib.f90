!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Non-Koopmans method
! Developed and implemented by I. Dabo (Universite Paris-Est, Ecole des Ponts, ParisTech)
! Parallelized by Andrea Ferretti (MIT)
!
!-----------------------------------------------------------------------
      subroutine calc_sic_potential(c,orbitalrhor,vsic,wxdsic,wrefsic,  &
                                    rhorefsic,vsicpsi,pink,tfirst)
!-----------------------------------------------------------------------
!
! ... calculate orbital densities and non-Koopmans potentials
!
      use kinds,              only: dp
      use gvecp,              only: ngm
      use gvecs,              only: ngs, nps, nms
      use gvecw,              only: ngw
      use recvecs_indexes,    only: np, nm
      use reciprocal_vectors, only: gstart
      use grid_dimensions,    only: nr1, nr2, nr3, nr1x, nr2x, nr3x, &
                                    nnrx
      use cell_base,          only: omega, a1, a2, a3
      use smooth_grid_dimensions, only: nr1s, nr2s, nr3s, &
                                        nr1sx, nr2sx, nr3sx, nnrsx
      use electrons_base,     only: nx => nbspx, n => nbsp, f, ispin, &
                                    nspin, nelt, nupdwn, iupdwn
      use constants,          only: pi, fpi
      use mp,                 only: mp_sum
      use io_global,          only: stdout, ionode
      use cp_interfaces,      only: fwfft, invfft
      use fft_base,           only: dffts, dfftp
!
      implicit none
      complex(dp) c(ngw,nx)
      complex(dp) vsicpsi(ngw,nx)
      real(dp) vsic(nnrx,nx)
      real(dp) pink(nx)
      real(dp) orbitalrhor(nnrx,2,nx)
      real(dp) wxdsic(nnrx,2,nx)
      real(dp) wrefsic(nnrx,nx)
      real(dp) rhorefsic(nnrx,2,nx)
      logical tfirst

      ! local variables

      integer  :: iss, isup, isdw, iss1, iss2, ios, &
                  i, ir, ig, k, j, ibnd
      real(dp) :: sa1, sa2,delta1,delta2,delta3
      real(dp) :: delta1s,delta2s,delta3s
      complex(dp) :: ci,fp,fm
      character(len=6), external :: int_to_char
      complex(dp), allocatable :: psi(:),psis(:),rhog(:,:),vsics(:), &
                                  vsicg(:),vsicpsis(:)
      real(dp), allocatable :: aux(:,:),rhor(:,:),rhos(:,:),wtot(:,:)
      complex(dp), allocatable :: orbitalrhog(:,:)
      real(dp), allocatable :: orbitalrhos(:,:)

      ci = (0.0d0,1.0d0)
      allocate(rhos(nnrsx,2))
      allocate(rhor(nnrx,2))
      allocate(rhog(ngm,2))
      rhor=0.d0
      rhos=0.d0
      rhog=(0.d0,0.d0)
      !
      ! ... calculate total charge density
      ! ... (this step can be removed depending on where the 
      ! ... subroutine is called, ID)
      !
      allocate(psis(nnrsx)) 
      !
      do i=1,n,2
        !
        call c2psi(psis,nnrsx,c(1,i),c(1,i+1),ngw,2)
        call invfft('Wave',psis,dffts)
        !
        iss1=ispin(i)
        sa1=f(i)/omega
        if(i.ne.n) then
          iss2=ispin(i+1)
          sa2=f(i+1)/omega
        else
          iss2=iss1
          sa2=0.0d0
        end if
        !
        do ir=1,nnrsx
          rhos(ir,iss1)=rhos(ir,iss1)+sa1*(dble(psis(ir)))**2
          rhos(ir,iss2)=rhos(ir,iss2)+sa2*(aimag(psis(ir)))**2
        end do
        !
      end do
      !
      deallocate(psis) 
      !
      allocate(psis(nnrsx)) 
      !
      isup=1
      isdw=2
      !
      do ir=1,nnrsx
        psis(ir)=cmplx(rhos(ir,isup),rhos(ir,isdw))
      end do
      !
      call fwfft('Smooth',psis,dffts)
      !
      do ig=1,ngs
         fp=psis(nps(ig))+psis(nms(ig))
         fm=psis(nps(ig))-psis(nms(ig))
         rhog(ig,isup)=0.5d0*cmplx( dble(fp),aimag(fm))
         rhog(ig,isdw)=0.5d0*cmplx(aimag(fp),-dble(fm))
      end do
      !
      allocate(psi(nnrx))
      !
      isup=1
      isdw=2
      psi(:)=(0.d0,0.d0)
      do ig=1,ngs
        psi(nm(ig))=conjg(rhog(ig,isup))+ci*conjg(rhog(ig,isdw))
        psi(np(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
      end do
      !
      call invfft('Dense',psi,dfftp)
      !
      do ir=1,nnrx
        rhor(ir,isup)=dble(psi(ir))
        rhor(ir,isdw)=aimag(psi(ir))
      end do
      !
      deallocate(psi) 
      deallocate(psis) 
      !
      ! ... calculate the density of each orbital
      !
      allocate(orbitalrhog(ngm,2))
      allocate(orbitalrhos(nnrsx,2))
      orbitalrhor=0.d0
      !
      do i=1,n
        !
        orbitalrhos=0.d0
        orbitalrhog=(0.d0,0.d0)
        !
        allocate(psis(nnrsx)) 
        !
        psis=0.d0
        do ig=1,ngw
          psis(nms(ig))=conjg(c(ig,i))
          psis(nps(ig))=c(ig,i)
        end do
        call invfft('Wave',psis,dffts)
        !
        iss1=ispin(i)
        sa1=1.d0/omega
        do ir=1,nnrsx
          orbitalrhos(ir,iss1)=sa1*(dble(psis(ir)))**2
        end do
        !
        deallocate(psis) 
        !
        allocate(psis(nnrsx)) 
        !
        isup=1
        isdw=2
        do ir=1,nnrsx
          psis(ir)=cmplx(orbitalrhos(ir,isup),orbitalrhos(ir,isdw))
        end do
        call fwfft('Smooth',psis, dffts )
        do ig=1,ngs
           fp=psis(nps(ig))+psis(nms(ig))
           fm=psis(nps(ig))-psis(nms(ig))
           orbitalrhog(ig,isup)=0.5d0*cmplx(dble(fp),aimag(fm))
           orbitalrhog(ig,isdw)=0.5d0*cmplx(aimag(fp),-dble(fm))
        end do
        !
        allocate(psi(nnrx))
        !
        isup=1
        isdw=2
        psi (:) = (0.d0, 0.d0)
        do ig=1,ngs
          psi(nm(ig))=conjg(orbitalrhog(ig,isup)) &
                     +ci*conjg(orbitalrhog(ig,isdw))
          psi(np(ig))=orbitalrhog(ig,isup)+ci*orbitalrhog(ig,isdw)
        end do
        !
        call invfft('Dense',psi,dfftp)
        !
        do ir=1,nnrx
          orbitalrhor(ir,isup,i)= dble(psi(ir))
          orbitalrhor(ir,isdw,i)=aimag(psi(ir))
        end do
        !
        deallocate(psi) 
        deallocate(psis) 
        !
      end do
      !
      allocate(aux(nnrx,2))
      !
      delta1=a1(1)/dble(nr1)
      delta2=a2(2)/dble(nr2)
      delta3=a3(3)/dble(nr3)
      !
      delta1s=a1(1)/dble(nr1s)
      delta2s=a2(2)/dble(nr2s)
      delta3s=a3(3)/dble(nr3s)
      !
      aux=0.0_dp
      do i=1,n
        aux(:,1)=aux(:,1)+orbitalrhor(:,1,i)
        aux(:,2)=aux(:,2)+orbitalrhor(:,2,i)
        call writetofile(orbitalrhor(:,1,i),nnrx, &
                        'rhoup'//trim(int_to_char(i))//'z.dat',dfftp, &
                        delta1,delta2,delta3,'az')
        call writetofile(orbitalrhor(:,2,i),nnrx, &
                        'rhodw'//trim(int_to_char(i))//'z.dat',dfftp, &
                        delta1,delta2,delta3,'az')
        call writetofile(orbitalrhor(:,1,i),nnrx, &
                        'rhoup'//trim(int_to_char(i))//'x.dat',dfftp, &
                        delta1,delta2,delta3,'ax')
        call writetofile(orbitalrhor(:,2,i),nnrx, &
                        'rhodw'//TRIM(int_to_char(i))//'x.dat',dfftp, &
                        delta1,delta2,delta3,'ax')
      end do
      !
      call writetofile(rhor(:,1),nnrx,'rhoupz.dat',dfftp, &
                       delta1,delta2,delta3,'az')
      call writetofile(rhor(:,2),nnrx,'rhodwz.dat',dfftp, &
                       delta1,delta2,delta3,'az')
      call writetofile(rhor(:,1),nnrx,'rhoupx.dat',dfftp, &
                       delta1,delta2,delta3,'ax')
      call writetofile(rhor(:,2),nnrx,'rhodwx.dat',dfftp, &
                       delta1,delta2,delta3,'ax' )
      !
      allocate(vsicg(nnrx))
      allocate(wtot(nnrx,2))
      allocate(vsics(nnrsx))
      allocate(psis(nnrsx))
      allocate(vsicpsis(nnrsx))
      vsicpsi=0.0_dp
      wtot=0.0_dp
      pink=0.0_dp
      !
      do i=1,n
        ibnd=i
        if(i.ge.iupdwn(2)) ibnd=i-iupdwn(2)+1
        call sic_correction(f(i),ispin(i),rhor(1,1),&
                            orbitalrhor(1,1,i),vsic(1,i),wxdsic(1,1,i), &
                            wrefsic(1,i),rhorefsic(1,1,i),pink(i),ibnd,tfirst)
        wtot(1:nnrx,1:2)=wtot(1:nnrx,1:2)+wxdsic(1:nnrx,1:2,i)
      end do
      !
      do i=1,n
        !
        call writetofile(vsic(:,i),nnrx, &
                         'vsic'//trim(int_to_char(i))//'.dat',dfftp, &
                         delta1,delta2,delta3,'az')
        vsic(1:nnrx,i)=vsic(1:nnrx,i)+wrefsic(1:nnrx,i)
        vsic(1:nnrx,i)=vsic(1:nnrx,i)+wtot(1:nnrx,ispin(i)) &
                      -wxdsic(1:nnrx,ispin(i),i)
        call writetofile(vsic(1,i),nnrx, &
                         'v+wxdsic'//trim(int_to_char(i))//'.dat',dfftp, &
                         delta1,delta2,delta3,'az')
        !
        psis=0.d0
        do ig=1,ngw
          psis(nms(ig))=conjg(c(ig,i))
          psis(nps(ig))=c(ig,i)
        end do
        call invfft('Wave',psis,dffts)
        !
        vsicg(1:nnrx)=vsic(1:nnrx,i)
        call fwfft('Dense',vsicg,dfftp)
        vsics=0.0_dp
        do ig=1,ngs
          vsics(nps(ig))=vsicg(np(ig))
          vsics(nms(ig))=conjg(vsicg(np(ig)))
        end do
        call invfft('Smooth',vsics,dffts)
        !
        vsicpsis=0.0_dp
        do ir = 1, nnrsx
          vsicpsis(ir)=cmplx(dble(vsics(ir))*dble(psis(ir)),0.0_dp)
        end do
        call fwfft('Wave',vsicpsis,dffts)
        !
        do ig=1,ngw
          vsicpsi(ig,i)=vsicpsis(nps(ig))
        end do
        !
      end do
      !
      deallocate(vsicg)
      deallocate(wtot)
      deallocate(vsics)
      deallocate(psis)
      deallocate(vsicpsis)
      !
      deallocate(rhos)
      deallocate(rhor)
      deallocate(rhog)
      deallocate(orbitalrhog)
      deallocate(orbitalrhos)
      deallocate(aux)
      !
      return
      !
!-----------------------------------------------------------------------
      end subroutine calc_sic_potential
!-----------------------------------------------------------------------
!---------------------------------------------------------------
      subroutine sic_correction(f,ispin,rho,rhoele,vsic,wxdsic,wrefsic, &
                                rhorefsic,pink,ibnd,tfirst) 
!---------------------------------------------------------------
!
! ... calculate the non-Koopmans potential from the orbital density
!
      use kinds,                only : dp
      use constants,            only : e2, fpi
      use cell_base,            only : tpiba2,omega
      use nksic,                only : fref, rhobarfact, nkmixfact, &
                                       do_nkmix, nknmax, nkscalfact, &
                                   vanishing_rho => vanishing_rho_w, &
                                       do_wref, do_wxd, update_rhoref
      use grid_dimensions,      only : nnrx, nr1, nr2, nr3
      use gvecp,                only : ngm
      use recvecs_indexes,      only : np, nm
      use reciprocal_vectors,   only : gstart, g
      use eecp_mod,             only : do_comp
      use cp_interfaces,        only : fwfft, invfft
      use fft_base,             only : dfftp
      use funct,                only : dmxc_spin
      use mp,                   only : mp_sum
      use mp_global,            only : intra_image_comm
      use io_global,            only : stdout
      use control_flags,        only : iprsta
      !
      implicit none
      real(dp),    intent(in)  :: f,rho(nnrx,2),rhoele(nnrx,2)
      real(dp),    intent(out) :: vsic(nnrx),wrefsic(nnrx), &
                                  wxdsic(nnrx,2),rhorefsic(nnrx,2)
      real(dp),    intent(out) :: pink
      integer,     intent(in)  :: ispin, ibnd
      logical,     intent(in)  :: tfirst
      !
      integer     :: i, is, ig, ir
      real(dp)    :: rh(2), rhc(2), exc_t, etxc, fact, pink_pz, ehele
      real(dp)    :: etxcref, etxc0, w2cst, dvxc(2), dmuxc(2,2)
      !
      real(dp),    allocatable :: rhobar(:,:)
      real(dp),    allocatable :: grhoraux(:,:,:)
      real(dp),    allocatable :: haux(:,:,:)
      real(dp),    allocatable :: rhoraux(:,:)
      real(dp),    allocatable :: vxc(:,:)
      real(dp),    allocatable :: vxc0(:,:)
      real(dp),    allocatable :: vxcref(:,:)
      complex(dp), allocatable :: aux(:)
      complex(dp), allocatable :: rhotmp(:)
      complex(dp), allocatable :: vtemp(:)
      complex(dp), allocatable :: vcorrtmp(:)
      integer,       parameter :: lsd=1
      !
      fact=omega/dble(nr1*nr2*nr3)
      !
      allocate(rhobar(nnrx,2))
      allocate(rhotmp(ngm))
      allocate(vtemp(ngm))
      allocate(vcorrtmp(ngm))
      allocate(aux(nnrx))
      !
      vsic=0.0_dp
      wrefsic=0.0_dp
      wxdsic=0.0_dp
      pink=0.0_dp
      !
      if(ibnd>nknmax.and.nknmax>0) return
      !
      rhc=0.D0
      rhobar=rho-f*rhoele
      !
      if(tfirst.and..not.update_rhoref) then
        rhorefsic=rhobar*rhobarfact+fref*rhoele
        if(iprsta>1) write(stdout,2000) 
      endif
      !
      !   compute self-hartree contributions
      !
      if(.not.update_rhoref) then
        !
        rhotmp=0.0_dp
        do is=1,2
          aux(:)=rhoele(:,is)
          call fwfft('Dense',aux,dfftp )
          do ig=1,ngm
            rhotmp(ig)=rhotmp(ig)+aux(np(ig))
          end do
        end do
        !
        if(gstart==2) vtemp(1)=(0.d0,0.d0)
        do ig=gstart,ngm
          vtemp(ig)=rhotmp(ig)*fpi/(tpiba2*g(ig))
        enddo
        !
        if(do_comp) then
          call calc_tcc_potential(vcorrtmp,rhotmp)
          vtemp=vtemp+vcorrtmp
        endif
        !
        aux=0.0_dp
        do ig=1,ngm
          aux(np(ig))=vtemp(ig)
          aux(nm(ig))=conjg(vtemp(ig))
        enddo
        call invfft('Dense',aux,dfftp)
        ehele=0.5_dp*f*f*sum(dble(aux(1:nnrx))*rhoele(1:nnrx,ispin))
        !
      endif 
      !
      rhotmp=0.0_dp
      do is=1,2
        if(update_rhoref) then
          aux(:)=rhoele(:,is) ! (fref-f) is included afterwards
        else
          aux(:)=rhorefsic(:,is)-rhobar(:,is)*rhobarfact-f*rhoele(:,is)
        endif
        call fwfft('Dense',aux,dfftp )
        do ig=1,ngm
          rhotmp(ig)=rhotmp(ig)+aux(np(ig))
        end do
      end do
      !
      if(gstart==2) vtemp(1)=(0.d0,0.d0)
      do ig=gstart,ngm
        vtemp(ig)=rhotmp(ig)*fpi/(tpiba2*g(ig))
      enddo
      !
      if(do_comp) then
        call calc_tcc_potential(vcorrtmp,rhotmp)
        vtemp=vtemp+vcorrtmp
      endif
      !
      aux=0.0_dp
      do ig=1,ngm
        aux(np(ig))=vtemp(ig)
        aux(nm(ig))=conjg(vtemp(ig))
      enddo
      call invfft('Dense',aux,dfftp)
      !
      if(update_rhoref) then
        ehele=0.5_dp*f*f*sum(dble(aux(1:nnrx))*rhoele(1:nnrx,ispin))
        aux=(fref-f)*aux
      end if
      !
      !   add self-xc contributions
      !
      allocate(grhoraux(nnrx,3,2))
      allocate(rhoraux(nnrx,2))
      allocate(vxc(nnrx,2))
      allocate(vxc0(nnrx,2))
      allocate(vxcref(nnrx,2))
      allocate(haux(nnrx,2,2))
      !
      grhoraux=0.0_dp
      haux=0.0_dp
      etxcref=0.0_dp
      vxcref=0.0_dp
      !
      if(update_rhoref) then
        rhoraux=rhobar*rhobarfact+fref*rhoele
      else
        rhoraux=rhorefsic
      endif
      call exch_corr_wrapper(nnrx,2,grhoraux,rhoraux,etxcref,vxcref,haux)
      !
      grhoraux=0.0_dp
      haux=0.0_dp
      etxc=0.0_dp
      vxc=0.0_dp
      !
      rhoraux=rhobar*rhobarfact+f*rhoele
      call exch_corr_wrapper(nnrx,2,grhoraux,rhoraux,etxc,vxc,haux)
      !
      grhoraux=0.0_dp
      haux=0.0_dp
      etxc0=0.0_dp
      vxc0=0.0_dp
      !
      rhoraux=rhobar*rhobarfact
      call exch_corr_wrapper(nnrx,2,grhoraux,rhoraux,etxc0,vxc0,haux)
      !
      vsic(1:nnrx)=dble(aux(1:nnrx)) &
                  +vxcref(1:nnrx,ispin)-vxc(1:nnrx,ispin)
      !
      pink=(etxc0-etxc) &
          +f*sum(vxcref(1:nnrx,ispin)*rhoele(1:nnrx,ispin)) &
          +ehele+f*sum(dble(aux(1:nnrx))*rhoele(1:nnrx,ispin))
      pink=pink*fact
      !
      call mp_sum(pink,intra_image_comm)
      !
      !   calculate the potential contributions related to the xc kernel  
      !
      if(do_wref.or.do_wxd) then
        !  
        if(update_rhoref) then
          rhoraux=rhobar*rhobarfact+fref*rhoele
        else
          rhoraux=rhorefsic
          if(iprsta>1) write(stdout,2020) 
        endif
        !
        do ir=1,nnrx
          dmuxc=0.0_dp
          rh(1:2)=rhoraux(ir,1:2)
          call dmxc_spin_cp(rh(1),rh(2),dmuxc(1,1),dmuxc(1,2), &
                            dmuxc(2,1),dmuxc(2,2),vanishing_rho)
          dvxc(1)=dmuxc(1,ispin)*rhoele(ir,ispin)*f
          dvxc(2)=dmuxc(2,ispin)*rhoele(ir,ispin)*f
          wrefsic(ir)=dmuxc(ispin,ispin)*rhoele(ir,ispin)
          do is=1,2
            if(dvxc(is).ne.0.0_dp) & 
              wxdsic(ir,is)=rhobarfact*(vxc0(ir,is)+dvxc(is)-vxc(ir,is))
          enddo
        enddo
        !
        wrefsic(1:nnrx)=wrefsic(1:nnrx)+dble(aux(1:nnrx))
        w2cst=sum(wrefsic(1:nnrx)*rhoele(1:nnrx,ispin))*fact
        call mp_sum(w2cst,intra_image_comm)
        do ir=1,nnrx
          wrefsic(ir)=fref*(wrefsic(ir)-w2cst)
        end do
        !
      endif
      !
      !   rescale contributions with the nkscalfact parameter
      !
      pink=pink*nkscalfact
      vsic=vsic*nkscalfact
      wxdsic=wxdsic*nkscalfact
      wrefsic=wrefsic*nkscalfact
      !
      !   functional admixture
      !
      if(do_nkmix.and.iprsta>1) write(stdout,2010) 
      !
      !   non-variational formulations
      !
      if(.not.do_wxd) wxdsic=0.d0
      if(.not.do_wref) wrefsic=0.d0
      !
      deallocate(grhoraux)
      deallocate(rhoraux)
      deallocate(vxc)
      deallocate(vxc0)
      deallocate(vxcref)
      deallocate(haux)
      deallocate(rhobar)
      deallocate(rhotmp)
      deallocate(vtemp)
      deallocate(vcorrtmp)
      deallocate(aux)
      !
2000  format(3X,'update reference density')
2010  format(3X,'do_nkmix is now obsolete, nothing done')
2020  format(3X,'warning: do_w used without update_rhoref')
      !
      return
      !
!---------------------------------------------------------------
      end subroutine sic_correction
!---------------------------------------------------------------
!---------------------------------------------------------------
      subroutine dmxc_spin_cp(rhoup,rhodw,dmuxc_uu, &
                              dmuxc_ud,dmuxc_du,dmuxc_dd,small)
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
      if(get_icorr().ne.1) return
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
      return
      !
!---------------------------------------------------------------
end subroutine dmxc_spin_cp
!---------------------------------------------------------------
