!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Non-Koopmans method
! Developed and implemented by I. Dabo 
!      (Universite Paris-Est, Ecole des Ponts, ParisTech)
! Parallelized and optimized by Andrea Ferretti 
!      (MIT, University of Oxford)
!
!-----------------------------------------------------------------------
      subroutine nksic_potential( c, rhor, tfirst )
!-----------------------------------------------------------------------
!
! ... calculate orbital densities and non-Koopmans potentials
!
      use kinds,                      only: dp
      use constants,                  only: pi, fpi
      use gvecp,                      only: ngm
      use gvecs,                      only: ngs, nps, nms
      use gvecw,                      only: ngw
      use recvecs_indexes,            only: np, nm
      use reciprocal_vectors,         only: gstart
      use grid_dimensions,            only: nr1, nr2, nr3, &
                                            nr1x, nr2x, nr3x, nnrx
      use cell_base,                  only: omega, a1, a2, a3
      use smooth_grid_dimensions,     only: nr1s, nr2s, nr3s, &
                                            nr1sx, nr2sx, nr3sx, nnrsx
      use electrons_base,             only: nx => nbspx, n => nbsp, f, &
                                            ispin, nspin, nelt, nupdwn, iupdwn
      use mp,                         only: mp_sum
      use io_global,                  only: stdout, ionode
      use cp_interfaces,              only: fwfft, invfft
      use fft_base,                   only: dffts, dfftp
      !
      use nksic,                      only: orb_rhor, vsic, wxdsic, &
                                            wrefsic, rhorefsic, vsicpsi, pink
      !
      implicit none
      !
      ! in/out vars
      !
      complex(dp), intent(in) :: c(ngw,nx)
      real(dp),    intent(in) :: rhor(nnrx,nspin)
      logical,     intent(in) :: tfirst

      !
      ! local variables
      !
      integer     :: i, ir, ig, k, j, ibnd
      !
      complex(dp),   allocatable :: psi(:),psis(:),vsics(:), &
                                    vsicg(:),vsicpsis(:)
      real(dp),      allocatable :: wtot(:,:)
      !
      character(len=6), external :: int_to_char

      !
      ! main body
      !
      CALL start_clock( 'nksic_drv' )

      !
      ! computing orbital densities
      !
      CALL nksic_get_orbitalrho( ngw, n, nnrx, c, orb_rhor )


      !
      ! compute the potentials
      !
      CALL start_clock( 'nksic_drv_2' )
      !
      allocate(vsicg(nnrx))
      allocate(wtot(nnrx,2))
      allocate(vsics(nnrsx))
      allocate(psis(nnrsx))
      allocate(vsicpsis(nnrsx))
      !
      vsicpsi=0.0_dp
      wtot=0.0_dp
      pink=0.0_dp
      !
      do i=1,n
          !
          ibnd=i
          if( i >= iupdwn(2) ) ibnd = i-iupdwn(2)+1
          !
          call nksic_correction( f(i), ispin(i), rhor,&
                       orb_rhor(:,i), vsic(:,i), wxdsic(:,:,i), &
                       wrefsic(:,i), rhorefsic(:,:,i), pink(i), ibnd, tfirst)
          !
          wtot(1:nnrx,1:2) = wtot(1:nnrx,1:2) + wxdsic(1:nnrx,1:2,i)
          !
      enddo
      !
      CALL stop_clock( 'nksic_drv_2' )
      CALL start_clock( 'nksic_drv_3' )
      !
      do i=1,n
          !
          call writetofile(vsic(:,i),nnrx, &
                        'vsic'//trim(int_to_char(i))//'.dat',dfftp, 'az')

          vsic(1:nnrx,i) = vsic(1:nnrx,i) + wrefsic(1:nnrx,i)
          vsic(1:nnrx,i) = vsic(1:nnrx,i) + wtot(1:nnrx,ispin(i)) &
                           -wxdsic(1:nnrx,ispin(i),i)

          call writetofile(vsic(1,i),nnrx, &
                        'v+wxdsic'//trim(int_to_char(i))//'.dat',dfftp, 'az')
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
          !
          vsics=0.0_dp
          do ig=1,ngs
              vsics(nps(ig))=vsicg(np(ig))
              vsics(nms(ig))=conjg(vsicg(np(ig)))
          end do
          !
          call invfft('Smooth',vsics,dffts)
          !
          vsicpsis=0.0_dp
          do ir = 1, nnrsx
              vsicpsis(ir)=cmplx(dble(vsics(ir))*dble(psis(ir)),0.0_dp)
          enddo
          !
          call fwfft('Wave',vsicpsis,dffts)
          !
          do ig=1,ngw
             vsicpsi(ig,i)=vsicpsis(nps(ig))
          enddo
          !
      enddo
      !
      CALL stop_clock( 'nksic_drv_3' )
      !
      deallocate(vsicg)
      deallocate(wtot)
      deallocate(vsics)
      deallocate(psis)
      deallocate(vsicpsis)
      !
      CALL stop_clock( 'nksic_drv' )
      return
      !
!-----------------------------------------------------------------------
      end subroutine nksic_potential
!-----------------------------------------------------------------------
!---------------------------------------------------------------
      subroutine nksic_correction( f, ispin, rho, orb_rho, vsic, wxdsic, &
                                   wrefsic, rhorefsic, pink, ibnd, tfirst) 
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
      real(dp),    intent(in)  :: f, rho(nnrx,2), orb_rho(nnrx)
      integer,     intent(in)  :: ispin, ibnd
      logical,     intent(in)  :: tfirst
      real(dp),    intent(out) :: vsic(nnrx),wrefsic(nnrx)
      real(dp),    intent(out) :: wxdsic(nnrx,2)
      real(dp),  intent(inout) :: rhorefsic(nnrx,2)
      real(dp),    intent(out) :: pink
      !
      integer     :: i, is, ig, ir
      real(dp)    :: rh(2), rhc(2), exc_t, etxc, fact, pink_pz, ehele
      real(dp)    :: etxcref, etxc0, w2cst, dvxc(2), dmuxc(2,2)
      !
      real(dp),    allocatable :: rhobar(:,:)
      real(dp),    allocatable :: rhoele(:,:)
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
      !
      CALL start_clock( 'nksic_corr' )
      CALL start_clock( 'nksic_corr_1' )

!
! XXX
! self-hartree has no contribution as a potential
! if fref = 1.2
!
      !
      fact=omega/dble(nr1*nr2*nr3)
      !
      allocate(rhobar(nnrx,2))
      allocate(rhoele(nnrx,2))
      allocate(rhotmp(ngm))
      allocate(vtemp(ngm))
      allocate(vcorrtmp(ngm))
      allocate(aux(nnrx))
      !
      rhoele=0.0d0
      rhoele(:,ispin) = orb_rho(:)
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
      if( tfirst .and. .not. update_rhoref ) then
          !
          rhorefsic = rhobar*rhobarfact + fref*rhoele
          if(iprsta>1) write(stdout,2000) 
          !
      endif
      !
      ! Compute self-hartree contributions
      !
      if( .not. update_rhoref ) then
          !
          rhotmp=0.0_dp
          do is=1,2
              !
              aux(:)=rhoele(:,is)
              call fwfft('Dense',aux,dfftp )
              !
              do ig=1,ngm
                  rhotmp(ig)=rhotmp(ig)+aux(np(ig))
              enddo
              !
          enddo
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
          !
          if( update_rhoref ) then
              aux(:)=rhoele(:,is) ! (fref-f) is included afterwards
          else
              aux(:)=rhorefsic(:,is)-rhobar(:,is)*rhobarfact-f*rhoele(:,is)
          endif
          !
          call fwfft('Dense',aux,dfftp )
          !
          do ig=1,ngm
              rhotmp(ig)=rhotmp(ig)+aux(np(ig))
          enddo
          !
      enddo
      !
      if(gstart==2) vtemp(1)=(0.d0,0.d0)
      do ig=gstart,ngm
          vtemp(ig)=rhotmp(ig)*fpi/(tpiba2*g(ig))
      enddo
      !
      if(do_comp) then
          !
          call calc_tcc_potential(vcorrtmp,rhotmp)
          vtemp=vtemp+vcorrtmp
          !
      endif
      !
      aux=0.0_dp
      do ig=1,ngm
          !
          aux(np(ig))=vtemp(ig)
          aux(nm(ig))=conjg(vtemp(ig))
          !
      enddo
      call invfft('Dense',aux,dfftp)
      !
      if( update_rhoref ) then
          !
          ehele=0.5_dp*f*f*sum(dble(aux(1:nnrx))*rhoele(1:nnrx,ispin))
          !
      endif
      !
      CALL stop_clock( 'nksic_corr_1' )
      CALL start_clock( 'nksic_corr_2' )
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
      if ( update_rhoref ) then
          !
          rhoraux = rhobar*rhobarfact + fref*rhoele
      else
          rhoraux = rhorefsic
          !
      endif
      call exch_corr_wrapper(nnrx,2,grhoraux,rhoraux,etxcref,vxcref,haux)
      !
      grhoraux=0.0_dp
      haux=0.0_dp
      etxc=0.0_dp
      vxc=0.0_dp
      !
!
! XXX maybe this call could be eliminated
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
! 
! XXX
! using a blas here ?
! 
      if( update_rhoref) then
          !
          vsic(1:nnrx)=(fref-f)*dble(aux(1:nnrx)) &
                       +vxcref(1:nnrx,ispin)-vxc(1:nnrx,ispin)
          !
          pink=(etxc0-etxc) &
              +f*sum(vxcref(1:nnrx,ispin)*rhoele(1:nnrx,ispin)) &
              +ehele+f*(fref-f)*sum(dble(aux(1:nnrx))*rhoele(1:nnrx,ispin))
          !
      else
          !
          vsic(1:nnrx)=dble(aux(1:nnrx)) &
                      +vxcref(1:nnrx,ispin)-vxc(1:nnrx,ispin)
          !
          pink=(etxc0-etxc) &
              +f*sum(vxcref(1:nnrx,ispin)*rhoele(1:nnrx,ispin)) &
              +ehele+f*sum(dble(aux(1:nnrx))*rhoele(1:nnrx,ispin))
          !
      endif        
      !
      pink=pink*fact
      !
      call mp_sum(pink,intra_image_comm)
      !
      CALL stop_clock( 'nksic_corr_2' )
      !
      !   calculate the potential contributions related to the xc kernel  
      !
      CALL start_clock( 'nksic_corr_3' )
      !
      if( do_wref .or. do_wxd ) then
        !  
        if(update_rhoref) then
          rhoraux=rhobar*rhobarfact+fref*rhoele
        else
          rhoraux=rhorefsic
          if(tfirst) write(stdout,2020) 
        endif
        !
        do ir=1,nnrx
          !
          dmuxc=0.0_dp
          rh(1:2)=rhoraux(ir,1:2)
          !
! XXX
! is it really needed to call this routine for 
! spin up and down at the same time ?
!
          call nksic_dmxc_spin_cp(rh(1),rh(2),dmuxc(1,1),dmuxc(1,2), &
                                  dmuxc(2,1),dmuxc(2,2),vanishing_rho)
          !
          dvxc(1)=dmuxc(1,ispin)*rhoele(ir,ispin)*f
          dvxc(2)=dmuxc(2,ispin)*rhoele(ir,ispin)*f
          wrefsic(ir)=dmuxc(ispin,ispin)*rhoele(ir,ispin)
          !
! XXX
! loop to be optimized
!
          do is=1,2
            if( dvxc(is).ne.0.0_dp ) & 
              wxdsic(ir,is)=rhobarfact*(vxc0(ir,is)+dvxc(is)-vxc(ir,is))
          enddo
        enddo
        !
        wrefsic(1:nnrx)=wrefsic(1:nnrx)+dble(aux(1:nnrx))
        w2cst=sum(wrefsic(1:nnrx)*rhoele(1:nnrx,ispin))*fact
        !
        call mp_sum(w2cst,intra_image_comm)
        !
!
! XXX
! loop to be avoided
!
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
      CALL stop_clock( 'nksic_corr_3' )
      !
      !   functional admixture
      !
      if( do_nkmix .and. iprsta>1 ) write(stdout,2010) 
      !
      !   non-variational formulations
      !
      if( .not. do_wxd )  wxdsic=0.d0
      if( .not. do_wref ) wrefsic=0.d0
      !
      deallocate(grhoraux)
      deallocate(rhoraux)
      deallocate(vxc)
      deallocate(vxc0)
      deallocate(vxcref)
      deallocate(haux)
      deallocate(rhobar)
      deallocate(rhoele)
      deallocate(rhotmp)
      deallocate(vtemp)
      deallocate(vcorrtmp)
      deallocate(aux)
      !
2000  format(3X,'update reference density')
2010  format(3X,'do_nkmix is now obsolete, nothing done')
2020  format(3X,'warning: do_w used without update_rhoref')
      !
      CALL stop_clock( 'nksic_corr' )
      return
      !
!---------------------------------------------------------------
      end subroutine nksic_correction
!---------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine nksic_get_orbitalrho( ngw, n, nnrx, c, orb_rhor )
!-----------------------------------------------------------------------
!
! Computes orbital densities on the real (not smooth) grid
!
      use kinds,                      only: dp
      use cp_interfaces,              only: fwfft, invfft
      use fft_base,                   only: dffts, dfftp
      use cell_base,                  only: omega
      use gvecp,                      only: ngm
      use gvecs,                      only: ngs, nps, nms
      use recvecs_indexes,            only: np, nm
      use smooth_grid_dimensions,     only: nnrsx

      !
      implicit none

      !
      ! input/output vars
      !
      integer,       intent(in)  :: ngw, n, nnrx
      complex(dp),   intent(in)  :: c(ngw,n)
      real(dp),      intent(out) :: orb_rhor(nnrx,n) 

      !
      ! local vars
      !
      character(20) :: subname='nksic_get_orbitalrho'
      integer       :: i, ir, ig, ierr
      real(dp)      :: sa1
      complex(dp)   :: ci, fm, fp
      complex(dp),  allocatable :: psis(:), psi(:)
      complex(dp),  allocatable :: orb_rhog(:,:)
      real(dp),     allocatable :: orb_rhos(:)

      !
      !====================
      ! main body
      !====================
      !
      call start_clock( 'nksic_orbrho' )

      ci = (0.0d0,1.0d0)
      !
      allocate( psis(nnrsx), stat=ierr )
      if ( ierr/=0 ) call errore(subname,'allocating psis',abs(ierr))
      allocate( psi(nnrx), stat=ierr )
      if ( ierr/=0 ) call errore(subname,'allocating psi',abs(ierr))
      allocate( orb_rhos(2), orb_rhog(ngm,2), stat=ierr )
      if ( ierr/=0 ) call errore(subname,'allocating orb_rhos, orb_rhog',abs(ierr))

      sa1 = 1.0d0 / omega 
      !
#define __SHORTCUT
#ifdef __NO_SHORTCUT
      do i = 1, n, 2
          !
          CALL c2psi( psis, nnrsx, c( 1, i ), c( 1, i+1 ), ngw, 2 )
          !
          CALL invfft('Wave',psis, dffts )
          !
          ! computing the orbital charge 
          ! in real space on the smooth grid
          !
          do ir = 1, nnrsx
              !
              orb_rhos(1) = sa1 * ( DBLE(psis(ir)) )**2 
              orb_rhos(2) = sa1 * ( AIMAG(psis(ir)) )**2 
              !
              psis( ir )  = cmplx( orb_rhos(1), orb_rhos(2) ) 
          enddo
          !
          ! orbital charges are taken to the G space
          !
          CALL fwfft('Smooth',psis, dffts )
          !
          do ig = 1, ngs
              !
              fp=psis(nps(ig))+psis(nms(ig))
              fm=psis(nps(ig))-psis(nms(ig))
              orb_rhog(ig,1)=0.5d0*cmplx(dble(fp),aimag(fm))
              orb_rhog(ig,2)=0.5d0*cmplx(aimag(fp),-dble(fm))
              !
          enddo
          !
          psi (:) = (0.d0, 0.d0)
          do ig=1,ngs
              !
              psi(nm(ig)) = conjg( orb_rhog(ig,1) ) &
                            +ci*conjg( orb_rhog(ig,2) )
              psi(np(ig)) = orb_rhog(ig,1) +ci*orb_rhog(ig,2)
              !
          enddo
          !
          call invfft('Dense',psi,dfftp)
          !
          do ir=1,nnrx
              !
              orb_rhor(ir,i)   =  dble(psi(ir))
              orb_rhor(ir,i+1) = aimag(psi(ir))
          enddo
          !
      enddo

#else

      do i = 1, n, 2
          !
          CALL c2psi( psi, nnrx, c( 1, i ), c( 1, i+1 ), ngw, 2 )
          !
          CALL invfft('Dense', psi, dfftp )
          !
          ! computing the orbital charge 
          ! in real space on the full grid
          !
          do ir = 1, nnrx
              !
              orb_rhor(ir,i)   = sa1 * ( DBLE(psi(ir)) )**2 
              orb_rhor(ir,i+1) = sa1 * ( AIMAG(psi(ir)) )**2 
              !
          enddo
          !
      enddo
#endif

      !
      !
      deallocate( psis )
      deallocate( psi )
      deallocate( orb_rhog )
      deallocate( orb_rhos )


      call stop_clock( 'nksic_orbrho' )
      return
      !
!---------------------------------------------------------------
end subroutine nksic_get_orbitalrho
!---------------------------------------------------------------

!---------------------------------------------------------------
      subroutine nksic_dmxc_spin_cp(rhoup,rhodw,dmuxc_uu, &
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
      !CALL start_clock( 'nksic_dmxc_spin_cp' )
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
      !CALL stop_clock( 'nksic_dmxc_spin_cp' )
      return
      !
!---------------------------------------------------------------
end subroutine nksic_dmxc_spin_cp
!---------------------------------------------------------------
