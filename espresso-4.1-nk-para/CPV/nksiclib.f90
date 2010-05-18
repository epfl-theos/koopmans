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
                                            wrefsic, rhoref, rhobar, pink
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
      integer     :: i, j, jj, ibnd
      real(dp)    :: focc
      real(dp),   allocatable :: wtot(:,:)

      !
      ! main body
      !
      CALL start_clock( 'nksic_drv' )

      !
      !
      ! compute the potentials
      !
      allocate(wtot(nnrx,2))
      !
      wtot=0.0_dp
      pink=0.0_dp

      !
      ! loop over bands
      ! (perform 2 ffts at the same time)
      !
      do j = 1, n, 2

          !
          ! computing orbital densities
          ! to be checked: nspin = 1, n odd
          !
          call nksic_get_orbitalrho( ngw, nnrx, c(:,j), c(:,j+1), orb_rhor )

          !
          ! computing orbital potantials
          !
          do jj = 1, 2
              ! 
              i = j + jj -1
              !
              ibnd=i
              if( nspin == 2 .and. i >= iupdwn(2) ) ibnd = i-iupdwn(2)+1
              ! NOTE: iupdwn(2) is set to zero if nspin = 1
              !
              focc = f(i) * dble( nspin ) / 2.0d0
 
              !
              ! define rhoref and rhobar
              !
              call nksic_get_rhoref( i, nnrx, ispin(i), nspin, &
                                     focc, rhor, orb_rhor(:,jj),&
                                     rhoref, rhobar, tfirst )

              !
              ! compute all the piecies to build the potentials
              !
              call nksic_correction( focc, ispin(i), orb_rhor(:,jj), &
                                     rhor, rhoref, rhobar, vsic(:,i), &
                                     wxdsic(:,:), wrefsic(:), pink(i), ibnd )

              !
              ! here information is accumulated over states
              ! This is way, wtot is added in the next loop
              !
              wtot(1:nnrx,1:2) = wtot(1:nnrx,1:2) + wxdsic(1:nnrx,1:2)

              !
              ! take care of spin symmetry
              !
              if ( nspin == 1 ) then
                  !
                  pink(i) = 2.0_dp * pink(i) 
                  wtot(1:nnrx,1) = wtot(1:nnrx,1) + wxdsic(1:nnrx,2)
                  wtot(1:nnrx,2) = wtot(1:nnrx,2) + wxdsic(1:nnrx,1)
                  !
              endif

              !
              ! ths sic potential is partly updated here 
              ! to save some memory
              !
              vsic(1:nnrx,i) = vsic(1:nnrx,i) + wrefsic(1:nnrx) &
                             -wxdsic( 1:nnrx,ispin(i) )
              ! 
          enddo
          !
      enddo

      !
      ! now wtot is completely built
      ! and can be added to vsic
      !
      do i = 1, n
          !
          vsic(1:nnrx,i) = vsic(1:nnrx,i) + wtot( 1:nnrx, ispin(i) ) 
          !
      enddo
      !
      deallocate(wtot)
      !
      CALL stop_clock( 'nksic_drv' )
      return
      !
!-----------------------------------------------------------------------
      end subroutine nksic_potential
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine nksic_get_orbitalrho( ngw, nnrx, c1, c2, orb_rhor )
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
      integer,       intent(in)  :: ngw, nnrx
      complex(dp),   intent(in)  :: c1(ngw), c2(ngw)
      real(dp),      intent(out) :: orb_rhor(nnrx,2) 

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
      allocate( psi(nnrx), stat=ierr )
      if ( ierr/=0 ) call errore(subname,'allocating psi',abs(ierr))

      sa1 = 1.0d0 / omega 
      !
#define __SHORTCUT
#ifdef __NO_SHORTCUT

      allocate( psis(nnrsx), stat=ierr )
      if ( ierr/=0 ) call errore(subname,'allocating psis',abs(ierr))
      allocate( orb_rhos(2), orb_rhog(ngm,2), stat=ierr )
      if ( ierr/=0 ) call errore(subname,'allocating orb_rhos, orb_rhog',abs(ierr))
      !
      CALL c2psi( psis, nnrsx, c1, c2, ngw, 2 )
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
          orb_rhor(ir,1) =  dble(psi(ir))
          orb_rhor(ir,2) = aimag(psi(ir))
      enddo

      deallocate( psis )
      deallocate( orb_rhog )
      deallocate( orb_rhos )

#else
      !
      CALL c2psi( psi, nnrx, c1, c2, ngw, 2 )
      !
      CALL invfft('Dense', psi, dfftp )
      !
      ! computing the orbital charge 
      ! in real space on the full grid
      !
      do ir = 1, nnrx
          !
          orb_rhor(ir,1) = sa1 * ( DBLE(psi(ir)) )**2 
          orb_rhor(ir,2) = sa1 * ( AIMAG(psi(ir)) )**2 
          !
      enddo
      !
#endif
      !
      deallocate( psi )


      call stop_clock( 'nksic_orbrho' )
      return
      !
!---------------------------------------------------------------
end subroutine nksic_get_orbitalrho
!---------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine nksic_get_rhoref( i, nnrx, ispin, nspin, f, &
                                   rhor, orb_rhor, rhoref_, rhobar_, tfirst )
!-----------------------------------------------------------------------
!
! Computes rhoref and rhobar
!
      use kinds,         only: dp
      use nksic,         only: rhoref0, update_rhoref, &
                               fref, rhobarfact
      !
      implicit none

      !
      ! input/output vars
      !
      integer,       intent(in)  :: i, nnrx
      integer,       intent(in)  :: ispin, nspin
      real(dp),      intent(in)  :: f
      real(dp),      intent(in)  :: rhor(nnrx,nspin)
      real(dp),      intent(in)  :: orb_rhor(nnrx)
      real(dp),      intent(out) :: rhoref_(nnrx,2)
      real(dp),      intent(out) :: rhobar_(nnrx,2)
      logical,       intent(in)  :: tfirst
      !
      integer :: ir, isp

      !
      ! main body
      !
      call start_clock( 'nksic_rhoref' )

      !
      ! define rhobar_i = rho - f_i * rho_i
      !
      if ( nspin == 1 ) then
          rhobar_(:,1) = rhor(:,1) * 0.5_dp
          rhobar_(:,2) = rhor(:,1) * 0.5_dp
      else
          rhobar_(:,1:2) = rhor(:,1:2)
      endif
      !
      rhobar_(:,ispin) = rhobar_(:,ispin) -f * orb_rhor(:)
      !
      ! probably obsolete
      if ( rhobarfact < 1.0d0 ) then 
          rhobar_ = rhobar_ * rhobarfact
      endif

      !
      ! define rhoref = rho + (f_ref -f_i) rho_i = rhobar_i + f_ref * rho_i
      !
      if ( .not. update_rhoref .and. .not. tfirst ) then
          !
          ! rhoref is not updated and it is taken
          ! from rhoref0
          !
          rhoref_(:,:) = rhoref0(:,:,i)
          !
      else
          !
          ! build rhoref from scratch
          !
          rhoref_(:,1:2)   = rhobar_(:,1:2)
          rhoref_(:,ispin) = rhoref_(:,ispin) + fref * orb_rhor(:)
          !
      endif
      !
      ! store rhoref if the case
      if ( tfirst .and. .not. update_rhoref ) then
          !
          rhoref0(:,:,i) = rhoref_(:,:)
          !
      endif
      !
      call stop_clock( 'nksic_rhoref' )
      return
      !
!---------------------------------------------------------------
end subroutine nksic_get_rhoref
!---------------------------------------------------------------

!---------------------------------------------------------------
      subroutine nksic_correction( f, ispin, orb_rhor, rhor, rhoref, rhobar, &
                                   vsic, wxdsic, wrefsic, pink, ibnd ) 
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
                                       do_wref, do_wxd, update_rhoref, &
                                       etxc, vxc => vxcsic
      use grid_dimensions,      only : nnrx, nr1, nr2, nr3
      use gvecp,                only : ngm
      use recvecs_indexes,      only : np, nm
      use reciprocal_vectors,   only : gstart, g
      use eecp_mod,             only : do_comp
      use cp_interfaces,        only : fwfft, invfft
      use fft_base,             only : dfftp
      use funct,                only : dmxc_spin, dft_is_gradient
      use mp,                   only : mp_sum
      use mp_global,            only : intra_image_comm
      use io_global,            only : stdout, ionode
      use control_flags,        only : iprsta
      use electrons_base,       only : nspin
      !
      implicit none
      integer,     intent(in)  :: ispin, ibnd
      real(dp),    intent(in)  :: f, orb_rhor(nnrx)
      real(dp),    intent(in)  :: rhor(nnrx,nspin)
      real(dp),    intent(in)  :: rhoref(nnrx,2)
      real(dp),    intent(in)  :: rhobar(nnrx,2)
      real(dp),    intent(out) :: vsic(nnrx), wrefsic(nnrx)
      real(dp),    intent(out) :: wxdsic(nnrx,2)
      real(dp),    intent(out) :: pink
      !
      integer     :: i, is, ig, ir
      real(dp)    :: rh(2), rhc(2), exc_t, fact, pink_pz, ehele, ehele2
      real(dp)    :: etxcref, etxc0, w2cst, dvxc(2), dmuxc(2,2)
      !
      real(dp),    allocatable :: rhoele(:,:)
      real(dp),    allocatable :: rhoraux(:,:)
      real(dp),    allocatable :: vxc0(:,:)
      real(dp),    allocatable :: vxcref(:,:)
      complex(dp), allocatable :: vhaux(:)
      complex(dp), allocatable :: vcorr(:)
      complex(dp), allocatable :: rhogaux(:)
      complex(dp), allocatable :: vtmp(:)
      !
      real(dp),    allocatable :: grhoraux(:,:,:)
      real(dp),    allocatable :: haux(:,:,:)

      !
      !==================
      ! main body
      !==================
      !
      if( ibnd > nknmax .and. nknmax > 0 ) return
      !
      CALL start_clock( 'nksic_corr' )
      CALL start_clock( 'nksic_corr_h' )

      !
      fact=omega/dble(nr1*nr2*nr3)
      !
      allocate(rhoele(nnrx,2))
      allocate(rhogaux(ngm))
      allocate(vtmp(ngm))
      allocate(vcorr(ngm))
      allocate(vhaux(nnrx))
      !
      rhoele=0.0d0
      rhoele(:,ispin) = orb_rhor(:)
      !
      vsic=0.0_dp
      wrefsic=0.0_dp
      wxdsic=0.0_dp
      pink=0.0_dp
      !
      rhc=0.D0

      !
      ! Compute self-hartree contributions
      !
      if( .not. update_rhoref ) then
          !
          ! when rhoref is not updated, we need to compute
          ! the hartree term for more then one density.
          ! The extra term is computed here.
          !
          rhogaux=0.0_dp
          do is=1,2
              !
              vhaux(:)=rhoele(:,is)
              call fwfft('Dense',vhaux,dfftp )
              !
              do ig=1,ngm
                  rhogaux(ig)=rhogaux(ig)+vhaux(np(ig))
              enddo
              !
          enddo
          !
          if(gstart==2) vtmp(1)=(0.d0,0.d0)
          do ig=gstart,ngm
              vtmp(ig)=rhogaux(ig)*fpi/(tpiba2*g(ig))
          enddo
          !
          if(do_comp) then
              call calc_tcc_potential(vcorr,rhogaux)
              vtmp=vtmp+vcorr
          endif
          !
          vhaux=0.0_dp
          do ig=1,ngm
              vhaux(np(ig))=vtmp(ig)
              vhaux(nm(ig))=conjg(vtmp(ig))
          enddo
          call invfft('Dense',vhaux,dfftp)
          !
          ehele=0.5_dp*f*f*sum(dble(vhaux(1:nnrx))*rhoele(1:nnrx,ispin))
          !
      endif 
      

      rhogaux=0.0_dp
      if ( update_rhoref ) then
          !
          vhaux(:) = rhoele(:,ispin) ! (fref-f) is included afterwards
          !
          call fwfft('Dense',vhaux,dfftp )
          !
          do ig=1,ngm
              rhogaux(ig) = vhaux( np(ig) )
          enddo
          !
      else
          !
          do is = 1, 2
              !
              vhaux(:) = rhoref(:,is)-rhobar(:,is)-f*rhoele(:,is)
              !
              call fwfft('Dense',vhaux,dfftp )
              !
              do ig=1,ngm
                  rhogaux(ig)=rhogaux(ig)+vhaux(np(ig))
              enddo
              !
          enddo
      endif

      !    
      ! compute hartree-like potential
      !
      if( gstart == 2 ) vtmp(1)=(0.d0,0.d0)
      do ig=gstart,ngm
          vtmp(ig) = rhogaux(ig) * fpi/( tpiba2*g(ig) )
      enddo
      !
      ! compute periodic corrections
      !
      if( do_comp ) then
          !
          call calc_tcc_potential( vcorr, rhogaux)
          vtmp(:) = vtmp(:) + vcorr(:)
          !
      endif
      
      vhaux=0.0_dp
      do ig=1,ngm
          !
          vhaux(np(ig)) = vtmp(ig)
          vhaux(nm(ig)) = conjg(vtmp(ig))
          !
      enddo
      call invfft('Dense',vhaux,dfftp)
      !
      ! init here wref sic to save some memrory
      !
      wrefsic(1:nnrx) = dble( vhaux(1:nnrx) )
      

      !
      ! self-hartree contrib to pink
      ! and init vsic
      !
      if( update_rhoref ) then
          !
          !ehele=0.5_dp * sum(dble(vhaux(1:nnrx))*rhoele(1:nnrx,ispin))
          !
          ehele = 2.0_dp * DBLE ( DOT_PRODUCT( vtmp(1:ngm), rhogaux(1:ngm)))
          if ( gstart == 2 ) ehele = ehele -DBLE ( CONJG( vtmp(1) ) * rhogaux(1) )
          !
          ehele = 0.5_dp * ehele * omega / fact

          !
          ! fref-f has to be included explicitly in rhoele
          !
          vsic(1:nnrx)=(fref-f)*dble(vhaux(1:nnrx)) 
          !
      else
          ! ehele has been computed at the beginning in this case
          ! compute ehele2
          ehele2 = f * sum(dble(vhaux(1:nnrx))*rhoele(1:nnrx,ispin))

          !
          ! fref-f has already included
          !
          vsic(1:nnrx)=dble(vhaux(1:nnrx)) 
          !
      endif
      !
      deallocate(rhogaux)
      deallocate(vtmp)
      deallocate(vcorr)
      deallocate(vhaux)
      !
      CALL stop_clock( 'nksic_corr_h' )
      CALL start_clock( 'nksic_corr_vxc' )
      !
      !   add self-xc contributions
      !
      if ( dft_is_gradient() ) then
         allocate(grhoraux(nnrx,3,2))
         allocate(haux(nnrx,2,2))
      else
         allocate(grhoraux(1,1,1))
         allocate(haux(1,1,1))
         grhoraux=0.0_dp
      endif
         
      if ( nspin == 1 .or. .not. update_rhoref ) then
          allocate(rhoraux(nnrx,2))
      endif
      allocate(vxc0(nnrx,2))
      allocate(vxcref(nnrx,2))


      etxcref=0.0_dp
      vxcref=0.0_dp
      !
      !rhoraux = rhoref
      call exch_corr_wrapper(nnrx,2,grhoraux,rhoref,etxcref,vxcref,haux)
      

      !
      ! this term is always computed if rhoref is not updated.
      ! otherwise is computed only for ibnd, ispin == 1 and stored
      !

      if ( .not. update_rhoref .or. ( ibnd == 1 .and. ispin == 1 ) ) then
          !
          etxc=0.0_dp
          vxc=0.0_dp
          !
          ! this is just a complication to save some memory 
          if ( nspin == 2 ) then 
              call exch_corr_wrapper(nnrx,2,grhoraux,rhor,etxc,vxc,haux)
          else
              !
              rhoraux = rhobar + f*rhoele
              call exch_corr_wrapper(nnrx,2,grhoraux,rhoraux,etxc,vxc,haux)
              !
          endif
          !
      endif


      etxc0=0.0_dp
      vxc0=0.0_dp
      !
      !rhoraux = rhobar
      call exch_corr_wrapper(nnrx,2,grhoraux,rhobar,etxc0,vxc0,haux)
      !
! 
! XXX
! using a blas here ?
! 
      if( update_rhoref) then
          !
          vsic(1:nnrx) = vsic(1:nnrx) &
                       + vxcref(1:nnrx,ispin)-vxc(1:nnrx,ispin)
          !
          pink=(etxc0-etxc) &
              +f*sum(vxcref(1:nnrx,ispin)*rhoele(1:nnrx,ispin)) &
              +(2.0_dp*fref-f) * ehele
          !
      else
          !
          vsic(1:nnrx) = vsic(1:nnrx) &
                       + vxcref(1:nnrx,ispin)-vxc(1:nnrx,ispin)
          !
          pink=(etxc0-etxc) &
              +f*sum(vxcref(1:nnrx,ispin)*rhoele(1:nnrx,ispin)) &
              +ehele +ehele2
          !
      endif        
      !
      pink=pink*fact
      !
      call mp_sum(pink,intra_image_comm)
      !
      CALL stop_clock( 'nksic_corr_vxc' )

      !
      !   calculate wref and wxd
      !
      CALL start_clock( 'nksic_corr_fxc' )
      !
      if( do_wref .or. do_wxd ) then
          !  
          ! note that vxd and wref are updated 
          ! (and not overwritten) by the next call
          wxdsic(:,:) = 0.0d0
          !
          call nksic_dmxc_spin_cp( nnrx, rhoref, f, ispin, rhoele, &
                                   vanishing_rho, & wrefsic, wxdsic )

!        do ir=1,nnrx
!            !
!            dmuxc=0.0_dp
!            rh(1:2) = rhoref( ir, 1:2)
!            !
!            call nksic_dmxc_spin_cp(rh(1),rh(2),dmuxc(1,1),dmuxc(1,2), &
!                                    dmuxc(2,1),dmuxc(2,2),vanishing_rho)
!            !
!            dvxc(1)=dmuxc(1,ispin)*rhoele(ir,ispin)*f
!            dvxc(2)=dmuxc(2,ispin)*rhoele(ir,ispin)*f
!            !
!            wrefsic(ir) = wrefsic(ir) + dmuxc(ispin,ispin)*rhoele(ir,ispin)
!            !
!            wxdsic(ir,1)=rhobarfact*(vxc0(ir,1)+dvxc(1)-vxc(ir,1))
!            wxdsic(ir,2)=rhobarfact*(vxc0(ir,2)+dvxc(2)-vxc(ir,2))     
!            !
!            !do is=1,2
!            !   if( dvxc(is).ne.0.0_dp ) & 
!            !     wxdsic(ir,is)=rhobarfact*(vxc0(ir,is)+dvxc(is)-vxc(ir,is))
!            !enddo
!            !
!        enddo
          !
          !
          if ( do_wref ) then
              !
              w2cst = sum( wrefsic(1:nnrx)*rhoele(1:nnrx,ispin) )*fact
              !
              call mp_sum(w2cst,intra_image_comm)
              !
              do ir=1,nnrx
                  wrefsic(ir)=fref*(wrefsic(ir)-w2cst)
              enddo
              !
          endif
          !
          if ( do_wxd ) then
              !
              wxdsic(:,1:2)= rhobarfact *( wxdsic(:,1:2) &
                                         + vxc0(:,1:2) -vxc(:,1:2) )
              !
          endif
          !
      endif
      !
      CALL stop_clock( 'nksic_corr_fxc' )

      !
      !   rescale contributions with the nkscalfact parameter
      !   take care of non-variational formulations
      !
      pink = pink * nkscalfact
      vsic = vsic * nkscalfact
      !
      if( do_wxd ) then 
          wxdsic = wxdsic * nkscalfact
      else
          wxdsic = 0.d0
      endif
      !
      if( do_wref ) then
          wrefsic = wrefsic * nkscalfact
      else
          wrefsic = 0.d0
      endif


      !
      !   functional admixture
      !
      if( do_nkmix .and. iprsta>1 .and. ionode ) &
         write(stdout,"(3X,'do_nkmix is now obsolete, nothing done')") 
      !
      if ( allocated(rhoraux) ) deallocate(rhoraux)
      deallocate(vxc0)
      deallocate(vxcref)
      deallocate(rhoele)
      !
      deallocate(grhoraux)
      deallocate(haux)
      !
      CALL stop_clock( 'nksic_corr' )
      return
      !
!---------------------------------------------------------------
      end subroutine nksic_correction
!---------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine nksic_eforce( i, ngw, c1, c2, vsicpsi )
!-----------------------------------------------------------------------
!
! Compute vsic potential for orbitals i and i+1 (c1 and c2) 
!
      use kinds,                only: dp
      use cp_interfaces,        only: fwfft, invfft
      use fft_base,             only: dffts, dfftp
      use gvecs,                only: ngs, nps, nms
      use nksic,                only: vsic
      use grid_dimensions,      only: nnrx

      !
      implicit none

      !
      ! input/output vars
      !
      integer,       intent(in)  :: i, ngw
      complex(dp),   intent(in)  :: c1(ngw), c2(ngw)
      complex(dp),   intent(out) :: vsicpsi(ngw, 2)

      !
      ! local vars
      !
      character(12) :: subname='nksic_eforce'
      integer       :: ir, ig, ierr, j
      real(dp)      :: wfc(2)
      complex(dp)   :: fm, fp
      complex(dp),  allocatable :: psi(:)

      !
      !====================
      ! main body
      !====================
      !
      call start_clock( 'nksic_eforce' )
      !
      allocate( psi(nnrx), stat=ierr )
      if ( ierr/=0 ) call errore(subname,'allocating psi',abs(ierr))
      !
      !

#ifdef __I_AM_SLOW
      ! this part is to be eliminated

      CALL nksic_eforce_std()

#else
      !
      CALL c2psi( psi, nnrx, c1, c2, ngw, 2 )
      !
      CALL invfft('Dense', psi, dfftp )

      !
      ! computing the orbital wfcs
      ! and the potentials in real space on the full grid
      !
      do ir = 1, nnrx
          !
          wfc(1)    =  DBLE( psi(ir) )
          wfc(2)    = AIMAG( psi(ir) )
          !
          psi( ir ) = CMPLX( wfc(1) * vsic(ir,i), &
                             wfc(2) * vsic(ir,i+1), DP ) 
          !
      enddo
      !
      CALL fwfft('Dense', psi, dfftp )
      !
      vsicpsi(:,:)=0.0_dp
      !
      do ig=1,ngw
          !
          fp = psi(nps(ig))+psi(nms(ig))
          fm = psi(nps(ig))-psi(nms(ig))
          !
          vsicpsi(ig,1)=0.5d0*cmplx(dble(fp),aimag(fm))
          vsicpsi(ig,2)=0.5d0*cmplx(aimag(fp),-dble(fm))
          !
      enddo
#endif
      !
      !
      deallocate( psi )

      call stop_clock( 'nksic_eforce' )
      return

#ifdef __I_AM_SLOW
      ! this part is to be eliminated
     
      !
      ! std way
      !
      CONTAINS
      !
      subroutine nksic_eforce_std()
      !
      use smooth_grid_dimensions,     only: nnrsx
      use recvecs_indexes,            only: np, nm
      implicit none
      !     
      complex(dp) :: c(ngw,2)
      complex(dp) :: psis(nnrsx)
      complex(dp) :: vsicg(ngs)
      complex(dp) :: vsics(nnrsx)
      complex(dp) :: vsicpsis(nnrsx)

      c(:,1) = c1
      c(:,2) = c2

      do j = 1, 2
          !
          psis=0.d0
          do ig=1,ngw
              psis(nms(ig))=conjg(c(ig,j))
              psis(nps(ig))=c(ig,j)
          end do
          call invfft('Wave',psis,dffts)
          !
          vsicg(1:nnrx)=vsic(1:nnrx,i+j-1)
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
             vsicpsi(ig,j)=vsicpsis(nps(ig))
          enddo
          !
      enddo
      !
      end subroutine nksic_eforce_std
#endif
      !
!---------------------------------------------------------------
end subroutine nksic_eforce
!---------------------------------------------------------------


!---------------------------------------------------------------
      subroutine nksic_dmxc_spin_cp( nnrx, rhoref, f, ispin, rhoele, &
                                     small, wref, wxd )
!---------------------------------------------------------------
!
! the derivative of the xc potential with respect to the local density
! is computed. 
! In order to save time, the loop over space coordinates is performed 
! inside this routine (inlining). 
!
! NOTE: wref and wsic are UPDATED and NOT OVERWRITTEN by this subroutine
!
      USE kinds,                ONLY : DP
      USE funct,                ONLY : xc_spin, get_iexch, get_icorr
      implicit none
      !
      integer,  intent(in)    :: nnrx, ispin
      real(dp), intent(in)    :: rhoref(nnrx,2), rhoele(nnrx,2)
      real(dp), intent(in)    :: f, small
      real(dp), intent(inout) :: wref(nnrx), wxd(nnrx,2)
      !
      real(dp) :: rhoup, rhodw, rhotot, zeta
      real(dp) :: dmuxc(2,2)
      real(dp) :: rs, fz, fz1, fz2, ex, vx, ecu, ecp, vcu, &
           vcp, dmcu, dmcp, aa, bb, cc, dr, dz, ec, vxupm, vxdwm, vcupm, &
           vcdwm, rho, vxupp, vxdwp, vcupp, vcdwp, dzm, dzp, fact
      real(dp), external :: dpz, dpz_polarized
      integer :: ir
      logical :: do_exch, do_corr
      !
      real(dp), parameter :: e2 = 2.0_dp, &
           pi34=0.75_dp/3.141592653589793_dp, third=1.0_dp/3.0_dp, &
           p43=4.0_dp/3.0_dp, p49=4.0_dp/ 9.0_dp, m23=-2.0_dp/3.0_dp

      !
      ! mian body
      !
      !CALL start_clock( 'nksic_dmxc_spin_cp' )
      !
      do_exch = ( get_iexch() == 1 )
      do_corr = ( get_icorr() == 1 )
      !
      ! main loop
      !
      do ir = 1, nnrx
          ! 
          dmuxc(:,:)=0.0_dp
          !
          rhoup  = rhoref(ir,1)
          rhodw  = rhoref(ir,2)
          rhotot = rhoup + rhodw
          !
          if( rhotot < small) cycle
          !
          zeta = (rhoup-rhodw)/rhotot
          if(abs(zeta)>1.0_dp) zeta=sign(1.0_dp,zeta)

          !
          ! calculate exchange contribution (analytical)
          !
          if ( do_exch ) then
              !
              if ( rhoup > small) then
                  rs=(pi34/(2.0_dp*rhoup))**third
                  call slater(rs,ex,vx)
                  dmuxc(1,1)=vx/(3.0_dp*rhoup)
              endif
              !
              if( rhodw > small) then
                  rs=(pi34/(2.0_dp*rhodw))**third
                  call slater(rs,ex,vx)
                  dmuxc(2,2)=vx/(3.0_dp*rhodw)
              endif
              !
          endif
          !

          !
          ! calculate correlation contribution (numerical)
          !
          if( do_corr ) then
              !
              dr   = min(1.e-6_dp,1.e-4_dp*rhotot)
              fact = 0.5d0 / dr
              !
              call xc_spin(rhotot-dr,zeta,ex,ec,vxupm,vxdwm,vcupm,vcdwm)
              call xc_spin(rhotot+dr,zeta,ex,ec,vxupp,vxdwp,vcupp,vcdwp)
              !
              dmuxc(1,1) = dmuxc(1,1) +(vcupp-vcupm) * fact
              dmuxc(1,2) = dmuxc(1,2) +(vcupp-vcupm) * fact
              dmuxc(2,1) = dmuxc(2,1) +(vcdwp-vcdwm) * fact
              dmuxc(2,2) = dmuxc(2,2) +(vcdwp-vcdwm) * fact
      
              dz=1.e-6_dp
              dzp=min(1.0,zeta+dz)-zeta
              dzm=-max(-1.0,zeta-dz)+zeta
              !
              fact = 1.0d0 / ( rhotot * (dzp+dzm) )
              ! 
              call xc_spin(rhotot,zeta-dzm,ex,ec,vxupm,vxdwm,vcupm,vcdwm)
              call xc_spin(rhotot,zeta+dzp,ex,ec,vxupp,vxdwp,vcupp,vcdwp)
              !
              dmuxc(1,1) = dmuxc(1,1) +(vcupp-vcupm)*(1.0_dp-zeta)*fact
              dmuxc(1,2) = dmuxc(1,2) -(vcupp-vcupm)*(1.0_dp+zeta)*fact
              dmuxc(2,1) = dmuxc(2,1) +(vcdwp-vcdwm)*(1.0_dp-zeta)*fact
              dmuxc(2,2) = dmuxc(2,2) -(vcdwp-vcdwm)*(1.0_dp+zeta)*fact
              !
          endif
          !
          wxd(ir,1) = wxd(ir,1) + dmuxc(1,ispin) * rhoele(ir,ispin)*f
          wxd(ir,2) = wxd(ir,2) + dmuxc(2,ispin) * rhoele(ir,ispin)*f
          !   
          wref(ir)  = wref(ir)  + dmuxc(ispin,ispin)*rhoele(ir,ispin)
          !   
      enddo

      return
      !
!---------------------------------------------------------------
end subroutine nksic_dmxc_spin_cp
!---------------------------------------------------------------
