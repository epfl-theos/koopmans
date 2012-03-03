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
! Further developed and optimized by Andrea Ferretti 
!      (MIT, University of Oxford)
!
!-----------------------------------------------------------------------
      subroutine nksic_potential( nbsp, nx, c, f_diag, bec, becsum, &
                                  deeq_sic, ispin, iupdwn, nupdwn, &
                                  rhor, rhog, wtot, vsic, pink)
!-----------------------------------------------------------------------
!
! ... calculate orbital dependent potentials, 
!     following the Non-Koopmans' (NK) scheme, 
!     but also Perdew-Zunger (PZ),
!     Non-Koopmans' integral definition (NKI),
!     Non-Joopmans on Perdew Zunger (PZNK)
!
      use kinds,                      only: dp
      use gvecp,                      only: ngm
      use gvecw,                      only: ngw
      use grid_dimensions,            only: nnrx
      use electrons_base,             only: nspin
      use nksic,                      only: orb_rhor, wxdsic, &
                                            wrefsic, rhoref, rhobar, &
                                            do_nk, do_nki, do_pz, do_nkpz, &
                                            grhobar, fion_sic
      use ions_base,                  only: nat
      use uspp,                       only: nkb
      use uspp_param,                 only: nhm
      !
      implicit none
      !
      ! in/out vars
      !
      integer,     intent(in)  :: nbsp, nx
      complex(dp), intent(in)  :: c(ngw,nx)
      real(dp),    intent(in)  :: bec(nkb,nbsp)
      real(dp),    intent(in)  :: becsum( nhm*(nhm+1)/2, nat, nspin)
      integer,     intent(in)  :: ispin(nx)
      integer,     intent(in)  :: iupdwn(nspin), nupdwn(nspin)
      real(dp),    intent(in)  :: f_diag(nx)
      real(dp),    intent(in)  :: rhor(nnrx,nspin)
      complex(dp), intent(in)  :: rhog(ngm,nspin)
      real(dp),    intent(out) :: vsic(nnrx,nx), wtot(nnrx,2)
      real(dp),    intent(out) :: deeq_sic(nhm,nhm,nat,nx)
      real(dp),    intent(out) :: pink(nx)

      !
      ! local variables
      !
      integer  :: i,j,jj,ibnd
      real(dp) :: focc,pinkpz
      real(dp), allocatable :: vsicpz(:)
      !
      ! main body
      !
      CALL start_clock( 'nksic_drv' )

      !
      ! compute potentials
      !
      if ( do_nk .or. do_nkpz .or. do_nki ) then
          wtot=0.0_dp
      endif
      !
      if ( do_nkpz ) then
          allocate(vsicpz(nnrx))
          vsicpz=0.0_dp
      endif
      !
      pink=0.0_dp
      vsic=0.0_dp


      !
      ! loop over bands (2 ffts at the same time)
      !
      do j=1,nbsp,2
        !
        ! compute orbital densities
        ! n odd => c(:,n+1) is already set to zero
        !
        call nksic_get_orbitalrho( ngw, nnrx, bec, ispin, nbsp, &
                                   c(:,j), c(:,j+1), orb_rhor, j, j+1)

        !
        ! compute orbital potentials
        !
        inner_loop: do jj=1,2
          ! 
          i=j+jj-1
          !
          ! this condition is important when n is odd
          !
          if ( i > nbsp ) exit inner_loop
          !
          ibnd=i
          if( nspin==2 ) then
              if ( i >= iupdwn(2) ) ibnd=i-iupdwn(2)+1
          endif
          !
          ! note: iupdwn(2) is set to zero if nspin = 1
          !
          focc=f_diag(i)*dble(nspin)/2.0d0
          !
          ! define rhoref and rhobar
          !
          call nksic_get_rhoref( i, nnrx, ispin(i), nspin,  &
                                 focc, rhor, orb_rhor(:,jj), &
                                 rhoref, rhobar, grhobar )

          !
          ! compute nk pieces to build the potentials and the energy
          !
          if ( do_nk .or. do_nkpz ) then
              !
              call nksic_correction_nk( focc, ispin(i), orb_rhor(:,jj), &
                                        rhor, rhoref, rhobar, grhobar, &
                                        vsic(:,i), wxdsic, wrefsic, &
                                        pink(i), ibnd)

              !
              ! here information is accumulated over states
              ! (wtot is added in the next loop)
              !
              wtot(1:nnrx,1:2) = wtot(1:nnrx,1:2) + wxdsic(1:nnrx,1:2)
              !
              ! ths sic potential is partly updated here to save some memory
              !
              vsic(1:nnrx,i) = vsic(1:nnrx,i) + wrefsic(1:nnrx) &
                             - wxdsic( 1:nnrx, ispin(i) )
              !
          endif

          ! 
          ! compute nkpz pieces to build the potential and the energy
          !
          if( do_nkpz ) then
              !
              call nksic_correction_nkpz( focc, orb_rhor(:,jj), vsicpz, &
                                          wrefsic, pinkpz, ibnd, ispin)
              !
              vsic(1:nnrx,i) = vsic(1:nnrx,i) + vsicpz(1:nnrx) &
                             + wrefsic(1:nnrx)
              !
              pink(i) = pink(i) +pinkpz
              !
          endif

          !
          ! compute pz potentials and energy
          !
          if ( do_pz ) then 
              !
              call nksic_correction_pz ( focc, ispin(i), orb_rhor(:,jj), &
                                         vsic(:,i), pink(i), ibnd )
              !
          endif

          !
          ! compute nki pieces to build the potentials and the energy
          !
          if ( do_nki ) then
              !
              call nksic_correction_nki( focc, ispin(i), orb_rhor(:,jj), &
                                         rhor, rhoref, rhobar, grhobar, &
                                         vsic(:,i), wxdsic, pink(i), ibnd)
              !
              ! here information is accumulated over states
              ! (wtot is added in the next loop)
              !
              wtot(1:nnrx,1:2) = wtot(1:nnrx,1:2) + wxdsic(1:nnrx,1:2)
              !
              ! ths sic potential is partly updated here to save some memory
              !
              vsic(1:nnrx,i) = vsic(1:nnrx,i) - wxdsic( 1:nnrx, ispin(i) )
              !
          endif

          !
          ! take care of spin symmetry
          !
          pink(i) = f_diag(i) * pink(i)
          !
          if ( do_nk .or. do_nkpz .or. do_nki ) then
              !
              if( nspin== 1 ) then
                  !
                  wtot(1:nnrx,1) = wtot(1:nnrx,1) + wxdsic(1:nnrx,2)
                  wtot(1:nnrx,2) = wtot(1:nnrx,2) + wxdsic(1:nnrx,1)
                  !
              endif
              !
          endif
          !
        enddo inner_loop
        !
      enddo

      !
      ! now wtot is completely built and can be added to vsic
      !
      if ( do_nk .or. do_nkpz .or. do_nki ) then
          !
          do i = 1, nbsp
              !
              vsic(1:nnrx,i) = vsic(1:nnrx,i) + wtot( 1:nnrx, ispin(i) ) 
              ! 
          enddo
          !
      endif
      !
      if ( allocated(vsicpz) ) deallocate(vsicpz)
      
      !
      ! USPP:
      ! compute corrections to the D coefficients of the pseudopots
      ! due to vsic(r, i) in the case of orbital dependent functionals.
      ! The corresponding contributions to the forces are computed.
      !
      ! IMPORTANT: the following call makes use of newd. 
      !            It must be done before we call newd for the
      !            total potentials, because deeq is overwritten at every call
      !
      fion_sic(:,:)     = 0.0d0
      !
      IF ( nhm > 0 ) then
          !
          deeq_sic(:,:,:,:) = 0.0d0
          !
          DO i = 1, nbsp
              !
              CALL nksic_newd( i, nnrx, ispin, nspin, vsic(:,i), nat, nhm, &
                               becsum, fion_sic, deeq_sic(:,:,:,i) )
              !
          ENDDO
          !
      ENDIF
      ! 
      ! 
      CALL stop_clock( 'nksic_drv' )
      return
      !
!-----------------------------------------------------------------------
      end subroutine nksic_potential
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine nksic_get_orbitalrho( ngw, nnrx, bec, ispin, nbsp, &
                                       c1, c2, orb_rhor, i1, i2 )
!-----------------------------------------------------------------------
!
! Computes orbital densities on the real (not smooth) grid
!
      use kinds,                      only: dp
      use constants,                  only: ci
      use cp_interfaces,              only: fwfft, invfft
      use fft_base,                   only: dffts, dfftp
      use cell_base,                  only: omega
      use gvecp,                      only: ngm
      use gvecs,                      only: ngs, nps, nms
      use recvecs_indexes,            only: np, nm
      use smooth_grid_dimensions,     only: nnrsx
      use cp_main_variables,          only: irb,eigrb
      use uspp_param,                 only: nhm
      use electrons_base,             only: nspin
      use ions_base,                  only: nat
      use uspp,                       only: okvan, nkb
      !
      implicit none

      !
      ! input/output vars
      !
      integer,     intent(in) :: ngw,nnrx,i1,i2
      integer,     intent(in) :: nbsp, ispin(nbsp)
      real(dp),    intent(in) :: bec(nkb, nbsp)
      complex(dp), intent(in) :: c1(ngw),c2(ngw)
      real(dp),   intent(out) :: orb_rhor(nnrx,2) 

      !
      ! local vars
      !
      character(20) :: subname='nksic_get_orbitalrho'
      integer       :: ir, ig, ierr
      real(dp)      :: sa1
      complex(dp)   :: fm, fp
      complex(dp), allocatable :: psis(:), psi(:)
      complex(dp), allocatable :: orb_rhog(:,:)
      real(dp),    allocatable :: orb_rhos(:)
      real(dp),    allocatable :: rhovan(:,:,:)
      real(dp),    allocatable :: rhovanaux(:,:,:)
      !
      !====================
      ! main body
      !====================
      !
      call start_clock( 'nksic_orbrho' )

      !
      if ( okvan ) then
          !
          allocate(rhovan(nhm*(nhm+1)/2,nat,nspin), stat=ierr )
          if ( ierr/=0 ) call errore(subname,'allocating rhovan',abs(ierr))
          allocate(rhovanaux(nhm*(nhm+1)/2,nat,nspin), stat=ierr)
          if ( ierr/=0 ) call errore(subname,'allocating rhovanaux',abs(ierr))
          !
      endif
      !
      allocate(psi(nnrx),stat=ierr)
      if ( ierr/=0 ) call errore(subname,'allocating psi',abs(ierr))
      !
      allocate(orb_rhog(ngm,2),stat=ierr)
      if ( ierr/=0 ) call errore(subname,'allocating orb_rhog',abs(ierr))

      sa1 = 1.0d0 / omega 


      !
      ! check whether it is necessary to 
      ! deal with the smooth and dense grids separately
      !
      if ( nnrsx == nnrx ) then
          !
          ! This case should be the one when using NCPP
          !
          CALL c2psi( psi, nnrx, c1, c2, ngw, 2 )
          !
          CALL invfft('Dense', psi, dfftp )
          !
          ! computing the orbital charge in real space on the full grid
          !
          do ir = 1, nnrx
              !
              orb_rhor(ir,1) = sa1 * ( DBLE(psi(ir)) )**2 
              orb_rhor(ir,2) = sa1 * ( AIMAG(psi(ir)) )**2 
              !
          enddo
          !
      else
          !
          ! this is the general case, 
          ! normally used with USPP
          !

          allocate( psis(nnrsx), stat=ierr )
          if ( ierr/=0 ) call errore(subname,'allocating psis',abs(ierr))
          allocate( orb_rhos(2), stat=ierr )
          if ( ierr/=0 ) call errore(subname,'allocating orb_rhos',abs(ierr))
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
          deallocate( orb_rhos )

      endif

      !
      ! add Vanderbilt contribution to orbital density
      !
      if( okvan ) then
        !
        rhovan(:,:,:) = 0.0d0
        !
        if ( nspin == 2 ) then 
            !
            if ( i2 <= nbsp ) then
                call calrhovan(rhovanaux,bec,i1)
                rhovan(:,:,1)=rhovanaux(:,:,ispin(i1))
            endif
            !
            if ( i2 <= nbsp ) then
                call calrhovan(rhovanaux,bec,i2)
                rhovan(:,:,2)=rhovanaux(:,:,ispin(i2))
            endif
            !
            call rhov(irb,eigrb,rhovan,orb_rhog,orb_rhor)
        else
            !
            if ( i2 <= nbsp ) then
                call calrhovan(rhovanaux,bec,i1)
                rhovan(:,:,1)=rhovanaux(:,:,ispin(i1))
                !
                call rhov(irb,eigrb,rhovan,orb_rhog(:,1),orb_rhor(:,1))
            endif
            !
            if ( i2 <= nbsp ) then
                call calrhovan(rhovanaux,bec,i2)
                rhovan(:,:,1)=rhovanaux(:,:,ispin(i2))
                !
                call rhov(irb,eigrb,rhovan,orb_rhog(:,2),orb_rhor(:,2))
            endif
            !
        endif
        !
      endif
      !
      deallocate(psi)
      deallocate(orb_rhog)
      !
      if ( okvan ) then
          deallocate(rhovan)
          deallocate(rhovanaux)
      endif
      !
      call stop_clock('nksic_orbrho')
      !
      return
      !
!---------------------------------------------------------------
end subroutine nksic_get_orbitalrho
!---------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine nksic_get_rhoref( i, nnrx, ispin, nspin, f, &
                                   rhor, orb_rhor, &
                                   rhoref_, rhobar_, grhobar_)
!-----------------------------------------------------------------------
!
! Computes rhoref and rhobar
!
      use kinds,                      only : dp
      use gvecp,                      only : ngm
      use funct,                      only : dft_is_gradient
      use cp_interfaces,              only : fwfft, invfft, fillgrad
      use fft_base,                   only : dfftp
      use recvecs_indexes,            only : np, nm
      use nksic,                      only : fref, rhobarfact
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
      real(dp),      intent(out) :: grhobar_(nnrx,3,2)
      !
      integer      :: ig
      complex(dp)  :: fp, fm
      complex(dp),   allocatable :: psi(:)
      complex(dp),   allocatable :: rhobarg(:,:)

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
      ! build rhoref from scratch
      !
      rhoref_(:,1:2)   = rhobar_(:,1:2)
      rhoref_(:,ispin) = rhoref_(:,ispin) + fref * orb_rhor(:)
      !
 
      !
      ! compute the gradient of rhobar if needed
      !
      if ( dft_is_gradient() ) then
          !
          allocate( rhobarg(ngm,2) )
          allocate( psi(nnrx) )
          !
          psi(:) = cmplx ( rhobar_(:,1), rhobar_(:,2) )
          !
          call fwfft('Dense',psi,dfftp )
          !
          do ig=1,ngm
              fp = psi( np(ig) ) +psi( nm(ig) )
              fm = psi( np(ig) ) -psi( nm(ig) )
              !
              rhobarg(ig,1) = 0.5d0 *cmplx( dble(fp),aimag(fm))
              rhobarg(ig,2) = 0.5d0 *cmplx(aimag(fp),-dble(fm))
          enddo
          !
          call fillgrad( 2, rhobarg, grhobar_ )
          !
          deallocate( psi )
          deallocate( rhobarg )
          !
      endif
      !
      call stop_clock( 'nksic_rhoref' )
      return
      !
!---------------------------------------------------------------
end subroutine nksic_get_rhoref
!---------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine nksic_newd( i, nnrx, ispin, nspin, vsic, nat, nhm, &
                             becsum, fion, deeq_sic )
!-----------------------------------------------------------------------
!
! computes the deeq coefficients (contributions to the D coeff of USPP)
! for the given orbital i. Coefficients are sotred in deeq_sic
!
      use kinds,                      only : dp
      use uspp,                       only : okvan, deeq
      use cp_main_variables,          only : irb, eigrb
      !
      implicit none

      !
      ! input/output vars
      !
      integer,       intent(in)    :: i, nnrx, nat, nhm
      integer,       intent(in)    :: ispin, nspin
      real(dp),      intent(in)    :: vsic(nnrx)
      real(dp),      intent(in)    :: becsum(nhm*(nhm+1)/2,nat,nspin) 
      real(dp),      intent(inout) :: fion(3,nat)
      real(dp),      intent(out)   :: deeq_sic(nhm,nhm,nat) 
      !
      ! local vars
      !
      real(dp),      allocatable   :: vsic_aux(:,:)

      !
      ! main body
      !
      if ( .not. okvan ) then 
          deeq_sic(:,:,:) = 0.0d0
          return
      endif
      !
      call start_clock( 'nksic_newd' )
      !
      allocate( vsic_aux(nnrx,nspin) )

      !
      ! fion are updated
      ! deeq coefficients are overwritten
      !
      vsic_aux = 0.0d0
      vsic_aux(:, ispin ) = vsic(:)
      !
      call newd( vsic_aux, irb, eigrb, becsum, fion )
      !
      deeq_sic(:,:,:) = deeq(:,:,:,ispin)

      deallocate( vsic_aux )
      !
      call stop_clock( 'nksic_newd' )
      return
      !
!---------------------------------------------------------------
end subroutine nksic_newd
!---------------------------------------------------------------

!---------------------------------------------------------------
      subroutine nksic_correction_nk( f, ispin, orb_rhor, rhor, &
                                      rhoref, rhobar, grhobar,  &
                                      vsic, wxdsic, wrefsic, pink, ibnd ) 
!---------------------------------------------------------------
!
! ... calculate the non-Koopmans potential from the orbital density
!
      use kinds,                only : dp
      use constants,            only : e2, fpi
      use cell_base,            only : tpiba2,omega
      use nksic,                only : fref, rhobarfact, nknmax, &
                                       vanishing_rho_w, &
                                       nkscalfact, do_wref, do_wxd, &
                                       etxc => etxc_sic, vxc => vxc_sic
      use grid_dimensions,      only : nnrx, nr1, nr2, nr3
      use gvecp,                only : ngm
      use recvecs_indexes,      only : np, nm
      use reciprocal_vectors,   only : gstart, g
      use eecp_mod,             only : do_comp
      use cp_interfaces,        only : fwfft, invfft, fillgrad
      use fft_base,             only : dfftp
      use funct,                only : dmxc_spin, dft_is_gradient
      use mp,                   only : mp_sum
      use mp_global,            only : intra_image_comm
      use electrons_base,       only : nspin
      !
      implicit none
      integer,     intent(in)  :: ispin, ibnd
      real(dp),    intent(in)  :: f, orb_rhor(nnrx)
      real(dp),    intent(in)  :: rhor(nnrx,nspin)
      real(dp),    intent(in)  :: rhoref(nnrx,2)
      real(dp),    intent(in)  :: rhobar(nnrx,2)
      real(dp),    intent(in)  :: grhobar(nnrx,3,2)
      real(dp),    intent(out) :: vsic(nnrx), wrefsic(nnrx)
      real(dp),    intent(out) :: wxdsic(nnrx,2)
      real(dp),    intent(out) :: pink
      !
      !character(19) :: subname='nksic_correction_nk'
      integer       :: ig, ir
      real(dp)      :: fact, ehele, etmp
      real(dp)      :: etxcref, etxc0, w2cst
      !
      real(dp),    allocatable :: rhoele(:,:)
      real(dp),    allocatable :: rhoraux(:,:)
      real(dp),    allocatable :: vxc0(:,:)
      real(dp),    allocatable :: vxcref(:,:)
      complex(dp), allocatable :: vhaux(:)
      complex(dp), allocatable :: vcorr(:)
      complex(dp), allocatable :: rhogaux(:,:)
      complex(dp), allocatable :: vtmp(:)
      !
      real(dp),    allocatable :: grhoraux(:,:,:)
      real(dp),    allocatable :: orb_grhor(:,:,:)
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
      allocate(rhogaux(ngm,1))
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
      ! Compute self-hartree contributions
      !
      rhogaux=0.0_dp
      !
      ! rhoele has no occupation
      !
      ! f-fref is NOT included here in vhaux
      ! (will be added afterwards)
      !
      vhaux(:) = rhoele(:,ispin)
      !
      call fwfft('Dense',vhaux,dfftp )
      !
      do ig=1,ngm
          rhogaux(ig,1) = vhaux( np(ig) )
      enddo

      !    
      ! compute hartree-like potential
      !
      if( gstart == 2 ) vtmp(1)=(0.d0,0.d0)
      do ig=gstart,ngm
          vtmp(ig) = rhogaux(ig,1) * fpi/( tpiba2*g(ig) )
      enddo
      !
      ! compute periodic corrections
      !
      if( do_comp ) then
          !
          call calc_tcc_potential( vcorr, rhogaux(:,1) )
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
      ! init here wref sic to save some memory
      !
      ! this is just the self-hartree potential 
      ! (to be multiplied by fref later on)
      !
      wrefsic(1:nnrx) = dble( vhaux(1:nnrx) )
      

      !
      ! self-hartree contrib to pink
      ! and init vsic
      !
      !ehele=0.5_dp * sum(dble(vhaux(1:nnrx))*rhoele(1:nnrx,ispin))
      !
      ehele = 2.0_dp * DBLE ( DOT_PRODUCT( vtmp(1:ngm), rhogaux(1:ngm,1)))
      if ( gstart == 2 ) ehele = ehele -DBLE ( CONJG( vtmp(1) ) * rhogaux(1,1) )
      !
      ! the f * (2.0d0 * fref-f) term is added here
      ehele = 0.5_dp * f * (2.0_dp * fref-f) * ehele * omega / fact

      !
      ! fref-f has to be included explicitly in rhoele
      !
      vsic(1:nnrx)=(fref-f)*dble(vhaux(1:nnrx)) 

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
           !
           allocate(grhoraux(nnrx,3,2))
           allocate(orb_grhor(nnrx,3,1))
           allocate(haux(nnrx,2,2))
           !
           ! compute the gradient of n_i(r)
           call fillgrad( 1, rhogaux, orb_grhor(:,:,1:1) )
           !
      else
           allocate(grhoraux(1,1,1))
           allocate(haux(1,1,1))
           grhoraux=0.0_dp
           !
      endif
      !
      deallocate(rhogaux)
      !   
      allocate(vxc0(nnrx,2))
      allocate(vxcref(nnrx,2))
      !
      etxcref=0.0_dp
      vxcref=0.0_dp

      !
      !rhoraux = rhoref
      !
      if ( dft_is_gradient() ) then
          !
          grhoraux(:,:,1:2)   = grhobar(:,:,1:2)
          grhoraux(:,:,ispin) = grhobar(:,:,ispin) &
                              + fref * orb_grhor(:,:,1)
      endif    
      !
      call exch_corr_wrapper(nnrx,2,grhoraux,rhoref,etxcref,vxcref,haux)
      

      !
      ! this term is computed for ibnd, ispin == 1 and stored
      ! or if rhobarfact < 1
      !
      if ( ( ibnd == 1 .and. ispin == 1) .OR. rhobarfact < 1.0_dp ) then
          !
          etxc=0.0_dp
          vxc=0.0_dp
          !
          ! some meory can be same in the nspin-2 case, 
          ! considering that rhobar + f*rhoele is identical to rho
          ! when rhobarfact == 1 
          !
          ! call exch_corr_wrapper(nnrx,2,grhoraux,rhor,etxc,vxc,haux)
          !
          allocate( rhoraux(nnrx, 2) )
          !
          rhoraux = rhobar + f*rhoele
          !
          if ( dft_is_gradient() ) then
              !
              grhoraux(:,:,1:2)   = grhobar(:,:,1:2)
              grhoraux(:,:,ispin) = grhobar(:,:,ispin) &
                                  + f * orb_grhor(:,:,1)
          endif
          !
          call exch_corr_wrapper(nnrx,2,grhoraux,rhoraux,etxc,vxc,haux)
          !
          deallocate( rhoraux )
          !
      endif

      etxc0=0.0_dp
      vxc0=0.0_dp
      !
      !rhoraux = rhobar
      !
      call exch_corr_wrapper(nnrx,2,grhobar,rhobar,etxc0,vxc0,haux)

      !
      ! update vsic pot 
      ! 
      vsic(1:nnrx) = vsic(1:nnrx) &
                   + vxcref(1:nnrx,ispin)-vxc(1:nnrx,ispin)
      !
      ! define pink
      !
      etmp = f*sum( vxcref(1:nnrx,ispin) * rhoele(1:nnrx,ispin) )
      !
      pink = ( etxc0-etxc ) + etmp + ehele
      pink = pink*fact
      !
      call mp_sum(pink,intra_image_comm)
      !
      call stop_clock( 'nksic_corr_vxc' )

      !
      !   calculate wref and wxd
      !
      CALL start_clock( 'nksic_corr_fxc' )
      !
      wxdsic(:,:) = 0.0d0
      !
      if( do_wref .or. do_wxd ) then
          !  
          ! note that vxd and wref are updated 
          ! (and not overwritten) by the next call
          !
          call nksic_dmxc_spin_cp( nnrx, rhoref, f, ispin, rhoele, &
                                   vanishing_rho_w, wrefsic, wxdsic )
          !
          !
          if ( do_wref ) then
              !
              w2cst = sum( wrefsic(1:nnrx) * rhoele(1:nnrx,ispin) ) * fact
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
      deallocate(vxc0)
      deallocate(vxcref)
      deallocate(rhoele)
      !
      deallocate(grhoraux)
      deallocate(haux)
      !
      if ( allocated(orb_grhor) ) deallocate(orb_grhor)
      !
      CALL stop_clock( 'nksic_corr' )
      return
      !
!---------------------------------------------------------------
      end subroutine nksic_correction_nk
!---------------------------------------------------------------

!---------------------------------------------------------------
      subroutine nksic_correction_pz( f, ispin, orb_rhor, &
                                      vsic, pink, ibnd ) 
!---------------------------------------------------------------
!
! ... calculate the non-Koopmans potential from the orbital density
!
      use kinds,                only : dp
      use constants,            only : e2, fpi
      use cell_base,            only : tpiba2,omega
      use nksic,                only : etxc => etxc_sic, vxc => vxc_sic, nknmax
!$$
      use nksic,                only : nkscalfact
!$$
      use grid_dimensions,      only : nnrx, nr1, nr2, nr3
      use gvecp,                only : ngm
      use recvecs_indexes,      only : np, nm
      use reciprocal_vectors,   only : gstart, g
      use eecp_mod,             only : do_comp
      use cp_interfaces,        only : fwfft, invfft, fillgrad
      use fft_base,             only : dfftp
      use funct,                only : dft_is_gradient
      use mp,                   only : mp_sum
      use mp_global,            only : intra_image_comm
      !
      implicit none
      integer,     intent(in)  :: ispin, ibnd
      real(dp),    intent(in)  :: f, orb_rhor(nnrx)
      real(dp),    intent(out) :: vsic(nnrx)
      real(dp),    intent(out) :: pink
      !
      !character(19) :: subname='nksic_correction_pz'
      integer       :: ig
      real(dp)      :: ehele, fact
      !
      real(dp),    allocatable :: rhoelef(:,:)
      complex(dp), allocatable :: rhogaux(:,:)
      complex(dp), allocatable :: vhaux(:)
      complex(dp), allocatable :: vcorr(:)
      complex(dp), allocatable :: vtmp(:)
      !
      real(dp),    allocatable :: grhoraux(:,:,:)
      real(dp),    allocatable :: haux(:,:,:)

      !
      !==================
      ! main body
      !==================
      !
      vsic=0.0_dp
      pink=0.0_dp
      !
      if ( ibnd > nknmax .and. nknmax > 0 ) return
      if ( f < 1.0d-6 ) return
      !
      CALL start_clock( 'nksic_corr' )
      CALL start_clock( 'nksic_corr_h' )

      !
      fact=omega/dble(nr1*nr2*nr3)
      !
      allocate(rhoelef(nnrx,2))
      allocate(rhogaux(ngm,1))
      allocate(vtmp(ngm))
      allocate(vcorr(ngm))
      allocate(vhaux(nnrx))
      !
      rhoelef=0.0d0
      rhoelef(:,ispin) = f * orb_rhor(:)

      !
      ! Compute self-hartree contributions
      !
      rhogaux=0.0_dp
      !
      ! rhoelef contains occupations
      !
      vhaux(:) = rhoelef(:,ispin)
      !
      call fwfft('Dense',vhaux,dfftp )
      !
      do ig=1,ngm
          rhogaux(ig,1) = vhaux( np(ig) )
      enddo

      !    
      ! compute hartree-like potential
      !
      if( gstart == 2 ) vtmp(1)=(0.d0,0.d0)
      do ig=gstart,ngm
          vtmp(ig) = rhogaux(ig,1) * fpi/( tpiba2*g(ig) )
      enddo
      !
      ! compute periodic corrections
      !
      if( do_comp ) then
          !
          call calc_tcc_potential( vcorr, rhogaux(:,1) )
          vtmp(:) = vtmp(:) + vcorr(:)
          !
      endif
      ! 
      vhaux=0.0_dp
      do ig=1,ngm
          !
          vhaux(np(ig)) = vtmp(ig)
          vhaux(nm(ig)) = conjg(vtmp(ig))
          !
      enddo
      call invfft('Dense',vhaux,dfftp)
      !
      ! init vsic
      !
      vsic(1:nnrx) =  -dble( vhaux(1:nnrx) )
      ehele        =   0.5_dp * sum( dble( vhaux(1:nnrx) ) &
                              * rhoelef(1:nnrx,ispin) )
      
      !
      ! partial cleanup
      !
      deallocate( vtmp )
      deallocate( vcorr )
      deallocate( vhaux )
      !
      CALL stop_clock( 'nksic_corr_h' )


      !
      ! Compute xc-contributions
      !
      if ( dft_is_gradient() ) then
          allocate(grhoraux(nnrx,3,2))
          allocate(haux(nnrx,2,2))
          !
          ! note: rhogaux contains the occupation
          !
          grhoraux=0.0_dp
          call fillgrad( 1, rhogaux, grhoraux(:,:,ispin:ispin) ) 
      else
          allocate(grhoraux(1,1,1))
          allocate(haux(1,1,1))
          !
          grhoraux=0.0_dp
      endif
      !
      deallocate( rhogaux )
      !
      vxc=0.0_dp
      haux=0.0_dp
      etxc=0.0_dp
      !
      call exch_corr_wrapper(nnrx,2,grhoraux,rhoelef,etxc,vxc,haux)
      !
!$$
      vsic(1:nnrx) =  vsic(1:nnrx) -vxc(1:nnrx,ispin)
!      vsic(1:nnrx) = -vxc(1:nnrx,ispin)
!$$
      !
      ! energy correction terms
      !
!$$
      pink = fact * ( -etxc -ehele )
!$$
!      pink = fact * ( -ehele )
!      pink = fact * ( -etxc )
!$$

!$$ This is for screened pz functional; apparently, I should have used a different variable name.
      !
      !   rescale contributions with the nkscalfact parameter
      !   take care of non-variational formulations
      !
      pink = pink * nkscalfact
      vsic = vsic * nkscalfact
      !
!$$

      call mp_sum(pink,intra_image_comm)
      
      !
      ! cleanup
      !
      deallocate( rhoelef )
      deallocate( grhoraux )
      deallocate( haux )
      !
      CALL stop_clock( 'nksic_corr' )
      !
      return
      !
!---------------------------------------------------------------
end subroutine nksic_correction_pz
!---------------------------------------------------------------


!---------------------------------------------------------------
      subroutine nksic_correction_nkpz( f, orb_rhor, vsic, wrefsic, pink, ibnd, ispin ) 
!---------------------------------------------------------------
!
! ... calculate the non-Koopmans potential on top of Perdew-Zunger, 
!     from the orbital densities
!
      use kinds,                only : dp
      use constants,            only : e2, fpi
      use cell_base,            only : tpiba2,omega
      use nksic,                only : fref, nkscalfact, &
                                       do_wref, vanishing_rho_w
      use grid_dimensions,      only : nnrx, nr1, nr2, nr3
      use gvecp,                only : ngm
      use recvecs_indexes,      only : np, nm
      use reciprocal_vectors,   only : gstart, g
      use eecp_mod,             only : do_comp
      use cp_interfaces,        only : fwfft, invfft, fillgrad
      use fft_base,             only : dfftp
      use funct,                only : dmxc_spin, dft_is_gradient
      use mp_global,            only : intra_image_comm
      use mp,                   only : mp_sum
      !
      implicit none
      real(dp),    intent(in)  :: f, orb_rhor(nnrx)
      integer,     intent(in)  :: ispin, ibnd
      real(dp),    intent(out) :: vsic(nnrx), wrefsic(nnrx)
      real(dp),    intent(out) :: pink
      !
      integer     :: ig, ir
      real(dp)    :: fact, etxcref
      real(dp)    :: w2cst
      !
      real(dp),    allocatable :: rhoele(:,:)
      real(dp),    allocatable :: rhoref(:,:)
      real(dp),    allocatable :: vxcref(:,:)
      real(dp),    allocatable :: wxdsic(:,:)
      real(dp),    allocatable :: grhoraux(:,:,:)
      real(dp),    allocatable :: haux(:,:,:)
      complex(dp), allocatable :: vhaux(:)
      complex(dp), allocatable :: vcorr(:)
      complex(dp), allocatable :: rhogaux(:,:)
      complex(dp), allocatable :: vtmp(:)
      !
      CALL start_clock( 'nksic_corr' )
      CALL start_clock( 'nksic_corr_h' )
      !
      fact=omega/dble(nr1*nr2*nr3)
      !
      allocate(wxdsic(nnrx,2))
      allocate(rhoele(nnrx,2))
      allocate(rhoref(nnrx,2))
      allocate(rhogaux(ngm,1))
      allocate(vtmp(ngm))
      allocate(vcorr(ngm))
      allocate(vhaux(nnrx))
      !
      rhoele=0.0d0
      rhoele(:,ispin)=orb_rhor(:)
      !
      vsic=0.0_dp
      wrefsic=0.0_dp
      wxdsic=0.0_dp
      pink=0.0_dp
      !
      ! compute self-hartree contributions
      !
      rhogaux=0.0_dp
      !
      ! rhoele has no occupation
      !
      vhaux(:) = rhoele(:,ispin)
      !
      call fwfft('Dense',vhaux,dfftp )
      !
      do ig=1,ngm
        rhogaux(ig,1) = vhaux( np(ig) )
      enddo
      !    
      ! compute hartree-like potential
      !
      if( gstart == 2 ) vtmp(1)=(0.d0,0.d0)
      do ig=gstart,ngm
        vtmp(ig)=rhogaux(ig,1)*fpi/(tpiba2*g(ig))
      enddo
      !
      ! compute periodic corrections
      !
      if( do_comp ) then
        !
        call calc_tcc_potential( vcorr, rhogaux(:,1))
        vtmp(:) = vtmp(:) + vcorr(:)
        !
      endif
      !
      vhaux=0.0_dp
      do ig=1,ngm
        !
        vhaux(np(ig)) = vtmp(ig)
        vhaux(nm(ig)) = conjg(vtmp(ig))
        !
      enddo
      !
      call invfft('Dense',vhaux,dfftp)
      !
      ! init here wref sic to save some memory
      !
      ! this is just the self-hartree potential 
      ! (to be multiplied by fref later on)
      !
      wrefsic(1:nnrx)=dble(vhaux(1:nnrx))
      !
      ! the term - fref has to be included explicitly in rhoele
      !
      vsic(1:nnrx)=-fref*dble(vhaux(1:nnrx)) 
      !
      deallocate(vtmp)
      deallocate(vcorr)
      deallocate(vhaux)
      !
      call stop_clock( 'nksic_corr_h' )
      call start_clock( 'nksic_corr_vxc' )
      !
      !   add self-xc contributions
      !
      rhoref=fref*rhoele
      !
      if ( dft_is_gradient() ) then
         allocate(grhoraux(nnrx,3,2))
         allocate(haux(nnrx,2,2))
         !
         grhoraux=0.0_dp
         call fillgrad( 1, rhogaux, grhoraux(:,:,ispin:ispin) ) 
         !
         grhoraux(:,:,ispin) = grhoraux(:,:,ispin) * fref
      else
         allocate(grhoraux(1,1,1))
         allocate(haux(1,1,1))
         grhoraux=0.0_dp
      endif
      !   
      deallocate(rhogaux)


      allocate(vxcref(nnrx,2))
      !
      etxcref=0.0_dp
      vxcref=0.0_dp
      !
      call exch_corr_wrapper(nnrx,2,grhoraux,rhoref,etxcref,vxcref,haux)
      !
      ! update vsic pot 
      ! 
      vsic(1:nnrx)=vsic(1:nnrx)-vxcref(1:nnrx,ispin)
      !
      ! define pink
      !
      pink=f*sum(vsic(1:nnrx)*rhoele(1:nnrx,ispin))*fact
      call mp_sum(pink,intra_image_comm)
      !
      call stop_clock( 'nksic_corr_vxc' )
      !
      !   calculate wref
      !
      CALL start_clock( 'nksic_corr_fxc' )
      !
      if( do_wref ) then
          !  
          ! note that wxd and wref are updated 
          ! (and not overwritten) by the next call
          !
          call nksic_dmxc_spin_cp(nnrx,rhoref,f,ispin,rhoele, &
                                  vanishing_rho_w,wrefsic,wxdsic)
          !
          w2cst=sum(wrefsic(1:nnrx)*rhoele(1:nnrx,ispin))*fact
          !
          call mp_sum(w2cst,intra_image_comm)
          !
          do ir=1,nnrx
            wrefsic(ir)=-fref*(wrefsic(ir)-w2cst)
          enddo
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
      if( do_wref ) then
          wrefsic = wrefsic * nkscalfact
      else
          wrefsic = 0.d0
      endif
      !
      deallocate(wxdsic)
      deallocate(vxcref)
      deallocate(rhoele)
      deallocate(rhoref)
      deallocate(grhoraux)
      deallocate(haux)
      !
      CALL stop_clock( 'nksic_corr' )
      return
      !
!---------------------------------------------------------------
      end subroutine nksic_correction_nkpz
!---------------------------------------------------------------


!---------------------------------------------------------------
      subroutine nksic_correction_nki( f, ispin, orb_rhor, rhor, &
                                       rhoref, rhobar, grhobar,  &
                                       vsic, wxdsic, pink, ibnd ) 
!---------------------------------------------------------------
!
! ... calculate the non-Koopmans (integrated, NKI) 
!     potential from the orbital density
!
!     note that fref=1.0 when performing NKI (i.e. it has a diff meaning)
!     then  rho_ref = rho - rho_i + n_i
!           rho_bar = rho - rho_i
!
      use kinds,                only : dp
      use constants,            only : e2, fpi
      use cell_base,            only : tpiba2,omega
      use nksic,                only : fref, rhobarfact, nknmax, &
                                       nkscalfact, do_wxd, &
                                       etxc => etxc_sic, vxc => vxc_sic
      use grid_dimensions,      only : nnrx, nr1, nr2, nr3
      use gvecp,                only : ngm
      use recvecs_indexes,      only : np, nm
      use reciprocal_vectors,   only : gstart, g
      use eecp_mod,             only : do_comp
      use cp_interfaces,        only : fwfft, invfft, fillgrad
      use fft_base,             only : dfftp
      use funct,                only : dmxc_spin, dft_is_gradient
      use mp,                   only : mp_sum
      use mp_global,            only : intra_image_comm
      use electrons_base,       only : nspin
      !
      implicit none
      integer,     intent(in)  :: ispin, ibnd
      real(dp),    intent(in)  :: f, orb_rhor(nnrx)
      real(dp),    intent(in)  :: rhor(nnrx,nspin)
      real(dp),    intent(in)  :: rhoref(nnrx,2)
      real(dp),    intent(in)  :: rhobar(nnrx,2)
      real(dp),    intent(in)  :: grhobar(nnrx,3,2)
      real(dp),    intent(out) :: vsic(nnrx)
      real(dp),    intent(out) :: wxdsic(nnrx,2)
      real(dp),    intent(out) :: pink
      !
      !character(20) :: subname='nksic_correction_nki'
      integer       :: ig
      real(dp)      :: fact, ehele, etmp
      real(dp)      :: etxcref, etxc0, w2cst
      !
      real(dp),    allocatable :: rhoele(:,:)
      real(dp),    allocatable :: rhoraux(:,:)
      real(dp),    allocatable :: vxc0(:,:)
      real(dp),    allocatable :: vxcref(:,:)
      complex(dp), allocatable :: vhaux(:)
      complex(dp), allocatable :: vcorr(:)
      complex(dp), allocatable :: rhogaux(:,:)
      complex(dp), allocatable :: vtmp(:)
      !
      real(dp),    allocatable :: grhoraux(:,:,:)
      real(dp),    allocatable :: orb_grhor(:,:,:)
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
      allocate(rhogaux(ngm,1))
      allocate(vtmp(ngm))
      allocate(vcorr(ngm))
      allocate(vhaux(nnrx))
      !
      rhoele=0.0d0
      rhoele(:,ispin) = orb_rhor(:)
      !
      vsic=0.0_dp
      wxdsic=0.0_dp
      pink=0.0_dp

      !
      ! Compute self-hartree contributions
      !
      rhogaux=0.0_dp
      !
      ! rhoele has no occupation
      !
      vhaux(:) = rhoele(:,ispin)
      !
      call fwfft('Dense',vhaux,dfftp )
      !
      do ig=1,ngm
          rhogaux(ig,1) = vhaux( np(ig) )
      enddo

      !    
      ! compute hartree-like potential
      !
      if( gstart == 2 ) vtmp(1)=(0.d0,0.d0)
      do ig=gstart,ngm
          vtmp(ig) = rhogaux(ig,1) * fpi/( tpiba2*g(ig) )
      enddo
      !
      ! compute periodic corrections
      !
      if( do_comp ) then
          !
          call calc_tcc_potential( vcorr, rhogaux(:,1))
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
      ! init here vsic to save some memory
      !
      ! this is just the self-hartree potential 
      !
      vsic(1:nnrx) = (1.0_dp-f) * dble( vhaux(1:nnrx) )
      

      !
      ! self-hartree contrib to pink
      ! and w2cst for vsic
      !
      !ehele=0.5_dp * sum(dble(vhaux(1:nnrx))*rhoele(1:nnrx,ispin))
      !
      ehele = 2.0_dp * DBLE ( DOT_PRODUCT( vtmp(1:ngm), rhogaux(1:ngm,1)))
      if ( gstart == 2 ) ehele = ehele -DBLE ( CONJG( vtmp(1) ) * rhogaux(1,1) )
      !
      ! -self-hartree energy to be added to the vsic potential
      w2cst = -0.5_dp * ehele * omega 
      call mp_sum(w2cst,intra_image_comm)
      !
      vsic  = vsic +w2cst
      !
      ! the f * (1-f) term is added here
      ehele = 0.5_dp * f * (1.0_dp-f) * ehele * omega / fact


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
           !
           allocate(grhoraux(nnrx,3,2))
           allocate(orb_grhor(nnrx,3,1))
           allocate(haux(nnrx,2,2))
           !
           ! compute the gradient of n_i(r)
           call fillgrad( 1, rhogaux, orb_grhor(:,:,1:1) )
           !
      else
           allocate(grhoraux(1,1,1))
           allocate(haux(1,1,1))
           grhoraux=0.0_dp
           !
      endif
      !
      deallocate(rhogaux)


      allocate(vxc0(nnrx,2))
      allocate(vxcref(nnrx,2))

      !
      ! this term is computed for ibnd, ispin == 1 and stored
      ! or if rhobarfact < 1
      !
      if ( ( ibnd == 1 .and. ispin == 1) .OR. rhobarfact < 1.0_dp ) then
          !
          etxc=0.0_dp
          vxc=0.0_dp
          !
          ! some meory can be same in the nspin-2 case, 
          ! considering that rhobar + f*rhoele is identical to rho
          ! when rhobarfact == 1 
          !
          ! call exch_corr_wrapper(nnrx,2,grhoraux,rhor,etxc,vxc,haux)
          !
          if ( dft_is_gradient() ) then
              !
              grhoraux(:,:,1:2)   = grhobar(:,:,1:2)
              grhoraux(:,:,ispin) = grhobar(:,:,ispin) &
                                  + f * orb_grhor(:,:,1)
          endif
          !
          allocate( rhoraux(nnrx, 2) )
          !
          rhoraux = rhobar + f*rhoele
          call exch_corr_wrapper(nnrx,2,grhoraux,rhoraux,etxc,vxc,haux)
          !
          deallocate( rhoraux )
          !
      endif

      !
      !rhoraux = rhoref
      !
      etxcref=0.0_dp
      vxcref=0.0_dp
      !
      if ( f == 1.0_dp ) then
          !
          vxcref=vxc
          etxcref=etxc
          !
      else
          !
          if ( dft_is_gradient() ) then
              !
              grhoraux(:,:,1:2)   = grhobar(:,:,1:2)
              grhoraux(:,:,ispin) = grhobar(:,:,ispin) &
                                  + fref * orb_grhor(:,:,1)
          endif
          !
          call exch_corr_wrapper(nnrx,2,grhoraux,rhoref,etxcref,vxcref,haux)
          !
      endif


      !
      !rhoraux = rhobar
      !
      etxc0=0.0_dp
      vxc0=0.0_dp
      !
      call exch_corr_wrapper(nnrx,2,grhobar,rhobar,etxc0,vxc0,haux)

      !
      ! update potential (including other constant terms)
      ! and define pink
      !
      etmp  = sum( vxcref(1:nnrx,ispin) * rhoele(1:nnrx,ispin) )
      w2cst = ( etxcref-etxc0 ) -etmp
      w2cst = w2cst * fact
      !
      call mp_sum(w2cst,intra_image_comm)
      !
      pink = (1.0_dp-f)*etxc0 -etxc + f*etxcref + ehele
      pink = pink*fact
      !
      call mp_sum(pink,intra_image_comm)
      !
      !
      ! update vsic pot 
      ! 
      vsic(1:nnrx) = vsic(1:nnrx) &
                   + vxcref(1:nnrx,ispin) -vxc(1:nnrx,ispin) + w2cst

      !
      !   calculate wxd
      !
      wxdsic(:,:) = 0.0d0
      !
      if( do_wxd ) then
          !
          wxdsic(:,1:2)= (1.0_dp-f)*vxc0(:,1:2) -vxc(:,1:2) +f*vxcref(:,1:2)
          !
      endif
      !
      call stop_clock( 'nksic_corr_vxc' )

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
      deallocate(vxc0)
      deallocate(vxcref)
      deallocate(rhoele)
      !
      deallocate(grhoraux)
      deallocate(haux)
      !
      if ( allocated(orb_grhor) ) deallocate(orb_grhor)
      !
      CALL stop_clock( 'nksic_corr' )
      return
      !
!---------------------------------------------------------------
      end subroutine nksic_correction_nki
!---------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine nksic_eforce( i, nbsp, nx, vsic, deeq_sic, bec, ngw, c1, c2, vsicpsi )
!-----------------------------------------------------------------------
!
! Compute vsic potential for orbitals i and i+1 (c1 and c2) 
!
      use kinds,                    only : dp
      use cp_interfaces,            only : fwfft, invfft
      use fft_base,                 only : dffts, dfftp
      use gvecs,                    only : ngs, nps, nms
      use grid_dimensions,          only : nnrx
      use smooth_grid_dimensions,   only : nnrsx
      use uspp,                     only : nkb, vkb
      use uspp_param,               only : nhm, nh
      use cvan,                     only : ish
      use ions_base,                only : nsp, na, nat

      !
      implicit none

      !
      ! input/output vars
      !
      integer,       intent(in)  :: i, nbsp, nx, ngw
      real(dp),      intent(in)  :: vsic(nnrx,nx)
      real(dp),      intent(in)  :: deeq_sic(nhm,nhm,nat,nx)
      real(dp),      intent(in)  :: bec(nkb,nbsp)
      complex(dp),   intent(in)  :: c1(ngw), c2(ngw)
      complex(dp),   intent(out) :: vsicpsi(ngw, 2)

      !
      ! local vars
      !
      character(12) :: subname='nksic_eforce'
      integer       :: ir, ig, ierr, j
      integer       :: is, iv, jv, isa, ism
      integer       :: ivoff, jvoff, ia, inl, jnl
      real(dp)      :: wfc(2), dd
      complex(dp)   :: fm, fp
      complex(dp),  allocatable :: psi(:)
      real(dp),     allocatable :: aa(:,:)

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
      ! init
      !
      vsicpsi(:,:) = 0.0d0
      !
      ! take advantage of the smooth and the dense grids 
      ! being equal (NCPP case)
      !
      if ( nnrsx == nnrx ) then
         !
         ! no need to take care of the double grid.
         ! typically, NCPP case

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

      else
          !
          ! here we take properly into account the 
          ! smooth and the dense grids
          ! typically, USPP case
          !
          CALL nksic_eforce_std()
          !
      endif
      !
      deallocate( psi )

      !
      ! add USPP non-local contribution 
      ! to the potantial 
      ! (this comes from the orbital-dependent piece of
      ! the potential)
      !
      if( nkb > 0 ) then
          !
          !     aa_i,i,n = sum_j d_i,ij <beta_i,j|c_n>
          ! 
          allocate( aa( nkb, 2 ) )
          !
          aa = 0.0d0
          !
          !
          do is = 1, nsp
              !
              do iv = 1, nh(is)
              do jv = 1, nh(is)
                  !
                  isa = 0
                  do ism = 1, is-1
                      isa = isa + na( ism )
                  enddo
                  !
                  ivoff = ish(is)+(iv-1)*na(is)
                  jvoff = ish(is)+(jv-1)*na(is)
                  ! 
                  if( i /= nbsp ) then
                      ! 
                      do ia=1,na(is)
                          inl = ivoff + ia
                          jnl = jvoff + ia
                          isa = isa + 1
                          !
                          dd  = deeq_sic(iv,jv,isa,i) 
                          aa(inl,1) = aa(inl,1) +dd * bec(jnl,i)
                          !
                          dd  = deeq_sic(iv,jv,isa,i+1) 
                          aa(inl,2) = aa(inl,2) +dd * bec(jnl,i+1)
                          !
                      enddo
                      ! 
                  else
                      ! 
                      do ia=1,na(is)
                          inl = ivoff + ia
                          jnl = jvoff + ia
                          isa = isa + 1
                          !
                          dd  = deeq_sic(iv,jv,isa,i) 
                          aa(inl,1) = aa(inl,1) +dd * bec(jnl,i)
                          !
                      enddo
                      ! 
                  endif
                  ! 
              enddo
              enddo
              !
          enddo
          ! 
          call DGEMM ( 'N', 'N', 2*ngw, 2, nkb, 1.0d0, &
                       vkb, 2*ngw, aa, nkb, 1.0d0, vsicpsi(:,:), 2*ngw)
          !
          deallocate( aa )
          !
      endif


      call stop_clock( 'nksic_eforce' )
      return
     
!
! implementation to deal with both 
! the smooth and the dense grids
!
CONTAINS
      !
      subroutine nksic_eforce_std()
      !
      use smooth_grid_dimensions,     only : nnrsx
      use recvecs_indexes,            only : np
      implicit none
      !     
      complex(dp) :: c(ngw,2)
      complex(dp) :: psis(nnrsx)
      complex(dp) :: vsicg(nnrx)
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
      USE kinds,                ONLY : dp
      USE funct,                ONLY : xc_spin, get_iexch, get_icorr
      implicit none
      !
      integer,  intent(in)    :: nnrx, ispin
      real(dp), intent(in)    :: rhoref(nnrx,2), rhoele(nnrx,2)
      real(dp), intent(in)    :: f, small
      real(dp), intent(inout) :: wref(nnrx), wxd(nnrx,2)
      !
      character(18) :: subname='nksic_dmxc_spin_cp'
      real(dp) :: rhoup, rhodw, rhotot, zeta
      real(dp) :: dmuxc(2,2)
      real(dp) :: rs, ex, vx, dr, dz, ec, &
                  vcupm, vcdwm, vcupp, vcdwp, dzm, dzp, fact
      !real(dp) :: vxupp, vxdwp, vxupm, vxdwm
      real(dp), external :: dpz, dpz_polarized
      integer :: ir
      !logical :: do_exch, do_corr
      !
      real(dp), parameter :: e2 = 2.0_dp, &
           pi34    = 0.6203504908994_DP,  & ! redefined to pi34=(3/4pi)^(1/3)
           pi34_old= 0.75_dp/3.141592653589793_dp, third=1.0_dp/3.0_dp, &
           p43=4.0_dp/3.0_dp, p49=4.0_dp/ 9.0_dp, m23=-2.0_dp/3.0_dp

      !
      ! mian body
      !
      !CALL start_clock( 'nksic_dmxc_spin_cp' )
      !
      ! the current implementation works only on top
      ! of LSD and LDA. Other functionals have to
      ! be implemented explicitly. To do that, we need to
      ! call the proper xc-routine (at the moment we call
      ! slater and pz_corr)
      !
      if ( get_iexch() /= 1 .or. get_icorr() /= 1 ) &
         call errore(subname,'only LDA/LSD PZ functionals implemented',10)
      !
      !do_exch = ( get_iexch() == 1 )
      !do_corr = ( get_icorr() == 1 )
      !
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
          if ( rhoup > small) then
              rs = pi34 / (2.0_dp*rhoup)**third
              call slater(rs,ex,vx)
              dmuxc(1,1)=vx/(3.0_dp*rhoup)
          endif
          !
          if( rhodw > small) then
              rs = pi34 / (2.0_dp*rhodw)**third
              call slater(rs,ex,vx)
              dmuxc(2,2)=vx/(3.0_dp*rhodw)
          endif

          !
          ! calculate correlation contribution (numerical)
          !
          dr   = min(1.e-6_dp,1.e-4_dp*rhotot)
          fact = 0.5d0 / dr
          !
          ! the explicit call to the correlation part only
          ! are performed instead of calling xc_spin.
          ! this saves some CPU time.
          ! unfortunately, different functionals have then
          ! to be treated explicitly
          !
          !call xc_spin(rhotot-dr,zeta,ex,ec,vxupm,vxdwm,vcupm,vcdwm)
          !call xc_spin(rhotot+dr,zeta,ex,ec,vxupp,vxdwp,vcupp,vcdwp)
          !
          rs = pi34 / (rhotot-dr)**third
          call pz_spin (rs, zeta, ec, vcupm, vcdwm)
          rs = pi34 / (rhotot+dr)**third
          call pz_spin (rs, zeta, ec, vcupp, vcdwp)
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
          !call xc_spin(rhotot,zeta-dzm,ex,ec,vxupm,vxdwm,vcupm,vcdwm)
          !call xc_spin(rhotot,zeta+dzp,ex,ec,vxupp,vxdwp,vcupp,vcdwp)
          !
          rs = pi34 / (rhotot)**third
          call pz_spin (rs, zeta-dzm, ec, vcupm, vcdwm)
          call pz_spin (rs, zeta+dzp, ec, vcupp, vcdwp)
          
          dmuxc(1,1) = dmuxc(1,1) +(vcupp-vcupm)*(1.0_dp-zeta)*fact
          dmuxc(1,2) = dmuxc(1,2) -(vcupp-vcupm)*(1.0_dp+zeta)*fact
          dmuxc(2,1) = dmuxc(2,1) +(vcdwp-vcdwm)*(1.0_dp-zeta)*fact
          dmuxc(2,2) = dmuxc(2,2) -(vcdwp-vcdwm)*(1.0_dp+zeta)*fact

          !
          ! add corrections to the nksic potentials
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

!$$
!-----------------------------------------------------------------------
      subroutine nksic_rot_emin(nouter,ninner,etot,Omattot)
!-----------------------------------------------------------------------
!
! ... Finds the orthogonal rotation matrix Omattot that minimizes
!     the orbital-dependent and hence the total energy, and then
!     rotate the wavefunction c0 accordingly.
!     We may need Omattot for further rotation of the gradient for outer loop CG.
!     Right now we do not do that because we set resetcg=.true. after inner loop
!     minimization routine, i.e., setting the search direction to be gradient direction.
!     (Ultrasoft pseudopotential case is not implemented.)
!
      use kinds,                      only : dp
      use grid_dimensions,            only : nnrx
      use gvecw,                      only : ngw
      use io_global,                  only : ionode
      use electrons_base,             only : nbsp, nbspx, nspin, &
                                             iupdwn,nupdwn
      use cp_interfaces,              only: invfft
      use fft_base,                   only: dfftp
      use ions_base,                  only: nsp, nat
      use uspp_param,                 only: nhm
      use nksic,                      only: vsic, pink,&
                                            do_nk, do_wref, do_wxd,&
                                            innerloop_nmax
      use uspp,                       only : nkb
      use cp_main_variables,          only : bec
      use wavefunctions_module,       only : c0, cm
      use control_flags,              only : esic_conv_thr
      !
      implicit none
      !
      ! in/out vars
      !
      integer                  :: ninner
      integer,     intent(in)  :: nouter
      real(dp),    intent(in)  :: etot
      real(dp)                 :: Omattot(nbspx,nbspx)


      !
      ! local variables
      !
      real(dp) :: esic,esic_old
      real(dp)                 :: bec1(nkb,nbsp)
      integer                  :: nbnd1,nbnd2
      real(dp)                 :: Omat1tot(nbspx,nbspx)
      real(dp)                 :: vsic1(nnrx,nbspx)
      complex(dp), allocatable :: Umat(:,:)
      complex(dp)              :: Umatbig(nbspx,nbspx)
      real(dp), allocatable    :: Heig(:)
      real(dp)                 :: Heigbig(nbspx)
      complex(dp)              :: wfn_ctmp(ngw,nbspx)
      integer                  :: npassofailmax
      integer, save            :: npassofail=0
      real(dp), save           :: passoprod=0.3d0
      real(dp)                 :: pink1(nbspx),dtmp,dalpha
      integer                  :: isp
      real(dp), allocatable    :: vsicah(:,:)
      real(dp)                 :: vsicah2sum,deigrms,dmaxeig
      integer :: nfile
      logical                  :: do_nonvar,lstopinner

      !
      ! variables for test calculations - along gradient line direction
      !
      logical :: ldotest
      
      !
      npassofailmax = 5 ! when to stop dividing passoprod by 2
      ldotest=.false.
      esic_old=0.d0

      Umatbig(:,:)=(0.d0,0.d0)
      Heigbig(:)=0.d0
      deigrms = 0.d0

      ! main body
      !
      CALL start_clock( 'nksic_rot_emin' )

      Omattot(:,:)=0.d0
      do nbnd1=1,nbspx
        Omattot(nbnd1,nbnd1)=1.d0
      enddo
           
      ninner = 0

      ldotest=.false.

      do while (.true.)

        ninner = ninner + 1

        if(ninner.gt.innerloop_nmax) then
          if(ionode) write(1031,*) '# innerloop_nmax reached.'
          if(ionode) write(1031,*)
          exit
        endif
        
!        if(mod(ninner,10).eq.1) ldotest=.true.
        if(mod(ninner,10).eq.1.or.ninner.le.5) ldotest=.true.
!        ldotest=.true.

!$$ Now do the test
        if(ldotest) then
          dtmp = 4.d0*3.141592d0
!          call nksic_rot_test(dtmp,201,nouter,ninner,etot)
          ldotest=.false.
        endif

        Omat1tot(:,:)=0.d0
        do nbnd1=1,nbspx
          Omat1tot(nbnd1,nbnd1)=1.d0
        enddo

!$$ This part calculates the anti-hermitian part of the hamiltonian
!$$ vsicah and see whether a convergence has been achieved

        wfn_ctmp(:,:) = (0.d0,0.d0)
        deigrms = 0.d0

        do isp=1,nspin

          allocate(Umat(nupdwn(isp),nupdwn(isp)))
          allocate(Heig(nupdwn(isp)))
          allocate(vsicah(nupdwn(isp),nupdwn(isp)))

          call nksic_getvsicah(isp,vsicah,vsicah2sum)
          call nksic_getHeigU(isp,vsicah,Heig,Umat)

          Umatbig(iupdwn(isp):iupdwn(isp)-1+nupdwn(isp),iupdwn(isp):iupdwn(isp)-1+nupdwn(isp)) = Umat(:,:)
          Heigbig(iupdwn(isp):iupdwn(isp)-1+nupdwn(isp)) = Heig(:)

!$$ The following file prints out the eigenvalues of the force matrix: for debugging
!$$
!$$          if(ionode) then
!$$            nfile=10000+isp
!$$            write(nfile,'(2I10,100F10.6)') ninner,nouter,sum(Heig(:)**2),Heig(:)
!$$          endif
!$$

          deigrms = deigrms + sum(Heig(:)**2)

          deallocate(Umat)
          deallocate(Heig)
          deallocate(vsicah)

        enddo ! do isp=1,nspin

        dmaxeig = max( abs(Heigbig(iupdwn(1))), abs(Heigbig(iupdwn(1)+nupdwn(1)-1)) )
        do isp=2,nspin
          dmaxeig = max(dmaxeig,abs(Heigbig(iupdwn(isp))))
          dmaxeig = max(dmaxeig,abs(Heigbig(iupdwn(isp)+nupdwn(isp)-1)))
        enddo

        ! how severe the transform is
        deigrms = sqrt(deigrms/nbsp)

!$$ print out ESIC part & other total energy
        esic=sum(pink(:))
        if(ionode) write(1031,'(2I10,3F24.13)') ninner, nouter,etot,esic,deigrms

        dalpha = passoprod/dmaxeig
        call nksic_getOmattot(dalpha,Heigbig,Umatbig,c0,wfn_ctmp,Omat1tot,bec1,vsic1,pink1,dtmp)

!$$ for nk0 which is not derived from an energy functional variation
        do_nonvar = (do_nk.and.(.not.do_wref.or..not.do_wxd))

        if(.not.do_nonvar) then
          lstopinner = dtmp.ge.esic
        else
          lstopinner = (ninner.ge.2.and.(esic-dtmp)*(esic-esic_old).gt.0.d0)
        endif
        lstopinner = lstopinner.or.(abs(esic-dtmp).lt.esic_conv_thr)

        if(lstopinner) then

          npassofail = npassofail+1

          if(ionode) then
            write(1031,'("# procedure  ",I4," / ",I4," is finished.")') npassofail,npassofailmax
            write(1031,*)
          endif

          if(npassofail.ge.npassofailmax) then
!$$ if we reach at the maximum allowed npassofail number, we exit without further update
            ninner = ninner + 1
            exit
          endif
          passoprod = passoprod * 0.5d0

!          ldotest=.true.

          cycle
        endif

!$$ we keep track of all the rotations to rotate cm later
        Omattot = MATMUL(Omattot,Omat1tot)

        pink(:) = pink1(:)
        vsic(:,:) = vsic1(:,:)
        c0(:,:) = wfn_ctmp(:,:)
        bec(:,:) = bec1(:,:)
        esic_old = esic

      enddo  !$$ do while (.true.)

!$$ Wavefunction cm rotation according to Omattot
! XXX We can put a if ( .not. tcg ) here (cm is relevant only for damped dynamics)
      if(ninner.ge.2) then
        wfn_ctmp(:,:) = (0.d0,0.d0)
        do nbnd1=1,nbspx
          do nbnd2=1,nbspx
            wfn_ctmp(:,nbnd1)=wfn_ctmp(:,nbnd1) + cm(:,nbnd2) * Omattot(nbnd2,nbnd1)
! XXX (we can think to use a blas, here)
          enddo
        enddo
        cm(:,1:nbspx) = wfn_ctmp(:,1:nbspx)
      endif

      CALL stop_clock( 'nksic_rot_emin' )
      return
      !
!---------------------------------------------------------------
end subroutine nksic_rot_emin
!---------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine nksic_rot_test(passoprod,nsteps,nouter,ninner,etot)
!-----------------------------------------------------------------------
!
! ... prints out esic by varying the wavefunction along a search direction.
!     (Ultrasoft pseudopotential case is not implemented.)
!
      use kinds,                      only : dp
      use grid_dimensions,            only : nnrx
      use gvecw,                      only : ngw
      use io_global,                  only : ionode
      use electrons_base,             only : nbsp, nbspx, nspin, &
                                             iupdwn, nupdwn
      use uspp,                       only : nkb
      use wavefunctions_module,       only : c0
      !
      implicit none
      !
      ! in/out vars
      !
      real(dp),    intent(in)  :: passoprod
      integer,     intent(in)  :: nsteps, ninner, nouter
      real(dp),    intent(in)  :: etot


      !
      ! local variables
      !
      real(dp) :: esic
      real(dp)                 :: bec1(nkb,nbsp)
      real(dp)                 :: Omat1tot(nbspx,nbspx)
      real(dp)                 :: vsic1(nnrx,nbspx)
      complex(dp), allocatable :: Umat(:,:)
      complex(dp)              :: Umatbig(nbspx,nbspx)
      real(dp), allocatable    :: Heig(:)
      real(dp)                 :: Heigbig(nbspx)
      complex(dp)              :: wfn_ctmp(ngw,nbspx)
      real(dp)                 :: dalpha,dmaxeig
      real(dp)                 :: pink1(nbspx)
      integer                  :: isp,istep
      real(dp), allocatable    :: vsicah(:,:)
      real(dp)                 :: vsicah2sum,deigrms
      integer                  :: nfile

      !
      ! variables for test calculations - along gradient line direction
      !

      !
      ! main body
      !
      CALL start_clock( 'nksic_rot_test' )


      Umatbig(:,:)=(0.d0,0.d0)
      Heigbig(:)=0.d0
      deigrms = 0.d0

      do isp=1,nspin

        allocate(Umat(nupdwn(isp),nupdwn(isp)))
        allocate(Heig(nupdwn(isp)))
        allocate(vsicah(nupdwn(isp),nupdwn(isp)))

        call nksic_getvsicah(isp,vsicah,vsicah2sum)
        call nksic_getHeigU(isp,vsicah,Heig,Umat)

        Umatbig(iupdwn(isp):iupdwn(isp)-1+nupdwn(isp),iupdwn(isp):iupdwn(isp)-1+nupdwn(isp)) = Umat(:,:)      
        Heigbig(iupdwn(isp):iupdwn(isp)-1+nupdwn(isp)) = Heig(:)      

        deigrms = deigrms + sum(Heig(:)**2)

        deallocate(Umat)
        deallocate(Heig)
        deallocate(vsicah)

      enddo ! do isp=1,nspin

      ! how severe the transform is
      deigrms = sqrt(deigrms/nbsp)

      dmaxeig = max( abs(Heigbig(iupdwn(1))), abs(Heigbig(iupdwn(1)+nupdwn(1)-1)) )
      do isp=2,nspin
        dmaxeig = max(dmaxeig,abs(Heigbig(iupdwn(isp))))
        dmaxeig = max(dmaxeig,abs(Heigbig(iupdwn(isp)+nupdwn(isp)-1)))
      enddo

      nfile = 10000+100*nouter+ninner
      if(ionode) write(nfile,*) '# passoprod',passoprod

      do istep=1,nsteps
        if(nsteps.ne.1) then
          dalpha = passoprod*(2.d0*istep-nsteps-1.d0)/(nsteps-1.d0) / dmaxeig
        else
          dalpha = 0.d0
        endif

        call nksic_getOmattot(dalpha,Heigbig,Umatbig,c0,wfn_ctmp,Omat1tot,bec1,vsic1,pink1,esic)

        if(ionode) write(nfile,'(5F24.13,2I10)') dalpha/3.141592*dmaxeig, dmaxeig, etot, esic, deigrms,ninner, nouter

      enddo  !$$ do istep=1,nsteps

      if(ionode) write(nfile,*)

      CALL stop_clock( 'nksic_rot_test' )
      return
      !
!---------------------------------------------------------------
end subroutine nksic_rot_test
!---------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine nksic_rot_emin_cg(nouter,ninner,etot,Omattot)
!-----------------------------------------------------------------------
!
! ... Finds the orthogonal rotation matrix Omattot that minimizes
!     the orbital-dependent and hence the total energy, and then
!     rotate the wavefunction c0 accordingly using cg minimization.
!     We may need Omattot for further rotation of the gradient for outer loop CG.
!     Right now we do not do that because we set resetcg=.true. after inner loop
!     minimization routine, i.e., setting the search direction to be gradient direction.
!     (Ultrasoft pseudopotential case is not implemented.)
!
      use kinds,                      only : dp
      use grid_dimensions,            only : nnrx
      use gvecw,                      only : ngw
      use io_global,                  only : ionode
      use electrons_base,             only : nbsp, nbspx, nspin, &
                                             iupdwn,nupdwn
      use cp_interfaces,              only: invfft
      use fft_base,                   only: dfftp
      use ions_base,                  only: nsp, nat
      use uspp_param,                 only: nhm
      use nksic,                      only: vsic, pink, &
                                            innerloop_cg_nsd, innerloop_cg_nreset,&
                                            innerloop_nmax
      use uspp,                       only : nkb
      use cp_main_variables,          only : bec
      use wavefunctions_module,       only : c0, cm
      use control_flags,              only : esic_conv_thr
      !
      implicit none
      !
      ! in/out vars
      ! 
      integer                  :: ninner
      integer,     intent(in)  :: nouter
      real(dp),    intent(in)  :: etot
      real(dp)                 :: Omattot(nbspx,nbspx)
        
      ! 
      ! local variables for cg routine
      ! 
      real(dp)                 :: bec1(nkb,nbsp),bec2(nkb,nbsp)
      integer                  :: nbnd1,nbnd2
      real(dp)                 :: Omat1tot(nbspx,nbspx),Omat2tot(nbspx,nbspx)
      real(dp)                 :: vsic1(nnrx,nbspx),vsic2(nnrx,nbspx)
      complex(dp), allocatable :: Umat(:,:)
      real(dp), allocatable    :: Heig(:)
      complex(dp)              :: Umatbig(nbspx,nbspx)
      real(dp)                 :: Heigbig(nbspx)
      complex(dp)              :: wfn_ctmp(ngw,nbspx),wfn_ctmp2(ngw,nbspx)
      real(dp)                 :: pink1(nbspx),pink2(nbspx),dtmp
      integer                  :: isp
      logical :: ldotest
      real(dp)                 :: ene0,ene1,enesti,enever,dene0
      real(dp)                 :: passo,passov,passof,passomax,spasso
      real(dp)                 :: gi(nbsp,nbsp) ! gradient
      real(dp)                 :: hi(nbsp,nbsp) ! search direction
      real(dp), allocatable    :: vsicah(:,:)
      real(dp)                 :: vsicah2sum,vsicah2sum_prev
      integer                  :: nidx1,nidx2
      real(dp)                 :: dPI,dalpha,dmaxeig,deigrms
      real(dp)                 :: pinksumprev,passoprod
      integer :: nfile
!

      !
      ! main body
      !
      ldotest=.false.

      pinksumprev=1.d8
      dPI = 2.0*asin(1.0)
      passoprod = (0.3d0/dPI)*dPI

      CALL start_clock( 'nksic_rot_emin' )

      Omattot(:,:)=0.d0
      do nbnd1=1,nbspx
        Omattot(nbnd1,nbnd1)=1.d0
      enddo

      ninner = 0

      ldotest=.false.

      hi(:,:) = 0.d0

      do while (.true.)

        ninner = ninner + 1

        if(ninner.gt.innerloop_nmax) then
          if(ionode) write(1031,*) '# innerloop_nmax reached.'
          if(ionode) write(1031,*)
          exit
        endif

!        call nksic_printoverlap(ninner,nouter)

!        if(mod(ninner,10).eq.1.or.ninner.le.5) ldotest=.true.
        if(ninner.eq.31.or.ninner.eq.61.or.ninner.eq.91) ldotest=.true.
!        if(ninner.le.10.and.nouter.eq.1) ldotest=.true.
!         ldotest=.true.
!        if(ninner.ge.25) ldotest=.true.
!$$ Now do the test
        if(ldotest) then
!          dtmp = 1.0d0*3.141592d0
          dtmp = 4.d0*3.141592d0
!          call nksic_rot_test(dtmp,201,nouter,ninner,etot)
          ldotest=.false.
        endif

!$$ print out ESIC part & other total energy
        ene0 = sum(pink(:))

!$$ test convergence
        if(abs(ene0-pinksumprev).lt.esic_conv_thr) then
!        if(ninner.ge.20) then
          if(ionode) then
            write(1037,*) '# inner-loop converged.'
            write(1037,*)

            write(1031,*) '# inner-loop converged.'
            write(1031,*)
          endif
          exit
        endif

        pinksumprev=ene0


        Omat1tot(:,:)=0.d0
        Omat2tot(:,:)=0.d0
        do nbnd1=1,nbspx
          Omat1tot(nbnd1,nbnd1)=1.d0
          Omat2tot(nbnd1,nbnd1)=1.d0
        enddo


!$$ This part calculates the anti-hermitian part of the Hamiltonian vsicah
!$$ and see whether a convergence has been achieved

        wfn_ctmp(:,:) = (0.d0,0.d0)
        wfn_ctmp2(:,:) = (0.d0,0.d0)

        gi(:,:) = 0.d0
        vsicah2sum = 0.d0
        Heigbig(:)=0.d0
        Umatbig(:,:)=(0.d0,0.d0)
        deigrms = 0.d0

!$$ For this run, we obtain the gradient

        do isp=1,nspin
          allocate(vsicah(nupdwn(isp),nupdwn(isp)))
          call nksic_getvsicah(isp,vsicah,dtmp)
          gi(iupdwn(isp):iupdwn(isp)-1+nupdwn(isp),iupdwn(isp):iupdwn(isp)-1+nupdwn(isp)) = vsicah(:,:)
          vsicah2sum = vsicah2sum + dtmp
          deallocate(vsicah)
        enddo ! do isp=1,nspin

!        call nksic_printoverlap(ninner,nouter)

        if(ninner.ne.1) dtmp = vsicah2sum/vsicah2sum_prev

        if(ninner.le.innerloop_cg_nsd.or.mod(ninner,innerloop_cg_nreset).eq.0) then
          hi(:,:) = gi(:,:)
        else
          hi(:,:) = gi(:,:) + dtmp*hi(:,:)
        endif

        do isp=1,nspin

          allocate(vsicah(nupdwn(isp),nupdwn(isp)))
          allocate(Umat(nupdwn(isp),nupdwn(isp)))
          allocate(Heig(nupdwn(isp)))

          vsicah(:,:) = hi(iupdwn(isp):iupdwn(isp)-1+nupdwn(isp),iupdwn(isp):iupdwn(isp)-1+nupdwn(isp))

          Heig(:)=0.d0
          Umat(:,:)=(0.d0,0.d0)
          call nksic_getHeigU(isp,vsicah,Heig,Umat)


          deigrms = deigrms + sum(Heig(:)**2)

          Umatbig(iupdwn(isp):iupdwn(isp)-1+nupdwn(isp),iupdwn(isp):iupdwn(isp)-1+nupdwn(isp)) = Umat(:,:)
          Heigbig(iupdwn(isp):iupdwn(isp)-1+nupdwn(isp)) = Heig(:)

          deallocate(vsicah)
          deallocate(Umat)
          deallocate(Heig)
        enddo ! do isp=1,nspin

        ! how severe the transform is
        deigrms = sqrt(deigrms/nbsp)

        if(ionode) write(1031,'(2I10,3F24.13)') ninner, nouter,etot,ene0,deigrms

        dmaxeig = max( abs(Heigbig(iupdwn(1))), abs(Heigbig(iupdwn(1)+nupdwn(1)-1)) )
        do isp=2,nspin
          dmaxeig = max(dmaxeig,abs(Heigbig(iupdwn(isp))))
          dmaxeig = max(dmaxeig,abs(Heigbig(iupdwn(isp)+nupdwn(isp)-1)))
        enddo

        passomax=passoprod/dmaxeig

        if(ninner.eq.1) then
          passof = passomax
          if(ionode) write(1031,*) '# passof set to passomax'
        endif

!$$$$        if(passof .gt. passomax*2.d0) then
!$$$$          passof = passomax*2.d0
!$$$$          if(ionode) write(1031,*) '# passof > twice passomax'
!$$$$        endif

!        if(ionode) then
!          write(1037,*)'# deigrms = ',deigrms
!          write(1037,*)'# vsicah2sum = ',vsicah2sum
!          if(ninner.ne.1) write(1037,*)'# vsicah2sum/vsicah2sum_prev = ',dtmp
!        endif

        vsicah2sum_prev = vsicah2sum

        dene0 = 0.d0
        do isp=1,nspin
          do nbnd1=1,nupdwn(isp)
            do nbnd2=1,nupdwn(isp)
              nidx1=nbnd1-1+iupdwn(isp)
              nidx2=nbnd2-1+iupdwn(isp)
              dene0 = dene0 - gi(nidx1,nidx2)*hi(nidx1,nidx2)
            enddo
          enddo
        enddo

!$$
!$$        dene0 = dene0 * 2.d0/nspin
!$$
!$$  Be careful, the following is correct because A_ji = - A_ij, i.e., the number of
!$$  linearly independent variables is half the number of total variables!
!$$
        dene0 = dene0 * 1.d0/nspin

        spasso = 1.d0
        if(dene0 .gt. 0.d0) spasso=-1.d0

        dalpha = spasso*passof
        call nksic_getOmattot(dalpha,Heigbig,Umatbig,c0,wfn_ctmp,Omat1tot,bec1,vsic1,pink1,ene1)
        call minparabola(ene0,spasso*dene0,ene1,passof,passo,enesti)

!$$$$ We neglect this step for paper writing purposes
!$$$$
        if(passo .gt. passomax) then
          passo = passomax
          if(ionode) write(1031,*) '# passo > passomax'
        endif

        passov = passof
        passof = 2.d0*passo

        dalpha = spasso*passo
!$$ The following line is for dene0 test
!        if(ninner.ge.15) dalpha = spasso*passo*0.00001
!$$
        call nksic_getOmattot(dalpha,Heigbig,Umatbig,c0,wfn_ctmp2,Omat2tot,bec2,vsic2,pink2,enever)

        if(ionode) then
          write(1037,*) ninner, nouter
          write(1037,'("ene0,ene1,enesti,enever")')
          write(1037,'(a3,4f20.10)') 'CG1',ene0,ene1,enesti,enever
          write(1037,'("spasso,passov,passo,passomax,dene0,&
                       & (enever-ene0)/passo/dene0")')
          write(1037,'(a3,4f12.7,e20.10,f12.7)')  &
            'CG2',spasso,passov,passo,passomax,dene0,(enever-ene0)/passo/dene0
          write(1037,*)
        endif

        if(ene0.lt.ene1.and.ene0.lt.enever) then
          if(ionode) then
            write(1037,'("# ene0<ene1 and ene0<enever!!  we exit!")')
            write(1037,*)
            write(1031,*) '# inner-loop NOT converged. we exit!'
            write(1031,*)
          endif
          ninner = ninner + 1
          exit
        endif

        if(ene1.ge.enever) then
          pink(:) = pink2(:)
          vsic(:,:) = vsic2(:,:)
          c0(:,:) = wfn_ctmp2(:,:)
          bec(:,:) = bec2(:,:)
          Omattot = MATMUL(Omattot,Omat2tot)
        else
          pink(:) = pink1(:)
          vsic(:,:) = vsic1(:,:)
          c0(:,:) = wfn_ctmp(:,:)
          bec(:,:) = bec1(:,:)
          Omattot = MATMUL(Omattot,Omat1tot)
          if(ionode) then
            write(1037,'("# It happened that ene1 < enever!!")')
            write(1037,*)
          endif
        endif

      enddo  !$$ do while (.true.)


!$$ Wavefunction cm rotation according to Omattot
      if(ninner.ge.2) then
        wfn_ctmp(:,:) = (0.d0,0.d0)
        do nbnd1=1,nbspx
          do nbnd2=1,nbspx
            wfn_ctmp(:,nbnd1)=wfn_ctmp(:,nbnd1) + cm(:,nbnd2) * Omattot(nbnd2,nbnd1)
          enddo
        enddo
        cm(:,1:nbspx) = wfn_ctmp(:,1:nbspx)
      endif
!$$ We need this because outer loop could be damped dynamics.

      CALL stop_clock( 'nksic_rot_emin' )
      return
      !
!---------------------------------------------------------------
end subroutine nksic_rot_emin_cg
!---------------------------------------------------------------

!---------------------------------------------------------------
      subroutine nksic_getOmattot(dalpha,Heigbig,Umatbig,wfn0,wfn1,Omat1tot,bec1,vsic1,pink1,ene1)
!---------------------------------------------------------------
!
! ... This routine rotates the wavefunction wfn0 into wfn1 according to
!     the force matrix (Heigbig, Umatbig) and the step of size dalpha.
!     Other quantities such as bec, vsic, pink are also calculated for wfn1.
!

      use kinds,                      only : dp
      use grid_dimensions,            only : nnrx
      use gvecw,                      only : ngw
      use electrons_base,             only : nbsp, nbspx, nspin, ispin, &
                                             iupdwn, nupdwn
      use ions_base,                  only : nsp
      use uspp,                       only : becsum,nkb
      use cp_main_variables,          only : eigr, rhor, rhog
      use nksic,                      only : deeq_sic, wtot, fsic
      !
      implicit none
      !
      ! in/out vars
      !
      real(dp),           intent(in) :: dalpha
      complex(dp),        intent(in) :: Umatbig(nbspx,nbspx)
      real(dp),           intent(in) :: Heigbig(nbspx)
      complex(dp),        intent(in) :: wfn0(ngw,nbspx)
      !
      ! XXX AF: variables to be properly allocated
      !
      complex(dp)  :: wfn1(ngw,nbspx)
      real(dp)     :: Omat1tot(nbspx,nbspx)
      real(dp)     :: bec1(nkb,nbsp)
      real(dp)     :: vsic1(nnrx,nbspx)
      real(dp)     :: pink1(nbspx)
      real(dp)     :: ene1

      !
      ! local variables for cg routine
      !
      integer                  :: isp,nbnd1
      real(dp), allocatable    :: Omat1(:,:)
      complex(dp), allocatable :: Umat(:,:)
      real(dp), allocatable    :: Heig(:)
      real(dp) :: dmaxeig
!

      Omat1tot(:,:) = 0.d0
      do nbnd1=1,nbspx
        Omat1tot(nbnd1,nbnd1)=1.d0
      enddo

      wfn1(:,:) = (0.d0,0.d0)

      dmaxeig = max( abs(Heigbig(iupdwn(1))), abs(Heigbig(iupdwn(1)+nupdwn(1)-1)) )
      do isp=2,nspin
        dmaxeig = max(dmaxeig,abs(Heigbig(iupdwn(isp))))
        dmaxeig = max(dmaxeig,abs(Heigbig(iupdwn(isp)+nupdwn(isp)-1)))
      enddo

!$$ For this run, we calculate ...
      do isp=1,nspin

        allocate(Umat(nupdwn(isp),nupdwn(isp)))
        allocate(Heig(nupdwn(isp)))
        allocate(Omat1(nupdwn(isp),nupdwn(isp)))

        Umat(:,:) = Umatbig(iupdwn(isp):iupdwn(isp)-1+nupdwn(isp),iupdwn(isp):iupdwn(isp)-1+nupdwn(isp))
        Heig(:) = Heigbig(iupdwn(isp):iupdwn(isp)-1+nupdwn(isp))

        call nksic_getOmat1(isp,Heig,Umat,dalpha,Omat1)

!$$ Wavefunction wfn0 is rotated into wfn0 using Omat1
        call nksic_rotwfn(isp,Omat1,wfn0,wfn1)

! Assigning the rotation matrix for a specific spin isp
        Omat1tot(iupdwn(isp):iupdwn(isp)-1+nupdwn(isp),iupdwn(isp):iupdwn(isp)-1+nupdwn(isp)) = Omat1(:,:)

        deallocate(Umat)
        deallocate(Heig)
        deallocate(Omat1)

      enddo ! do isp=1,nspin

!$$ recalculate bec & vsic according to the new wavefunction
      call calbec(1,nsp,eigr,wfn1,bec1)

      vsic1(:,:) = 0.d0
      pink1(:) = 0.d0
      call nksic_potential( nbsp, nbspx, wfn1, fsic, bec1, becsum, deeq_sic, &
                 ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic1, pink1 )

      ene1=sum(pink1(:))

      return
      !
!---------------------------------------------------------------
end subroutine nksic_getOmattot
!---------------------------------------------------------------



!-----------------------------------------------------------------------
      subroutine nksic_rotwfn(isp,Omat1,wfn1,wfn2)
!-----------------------------------------------------------------------
!
! ... Simple rotation of wfn1 into wfn2 by Omat1.
!     wfn2(n) = sum_m wfn1(m) Omat1(m,n)
!
      use electrons_base,             only : iupdwn,nupdwn,nbspx
      use gvecw,                      only : ngw
      use kinds,                      only : dp
      !
      implicit none
      !
      ! in/out vars
      !
      integer,     intent(in)  :: isp
      real(dp),    intent(in)  :: Omat1(nupdwn(isp),nupdwn(isp))
      complex(dp), intent(in)  :: wfn1(ngw,nbspx)
      complex(dp)              :: wfn2(ngw,nbspx)

      !
      ! local variables for cg routine
      !
      integer                  :: nbnd1,nbnd2

      CALL start_clock('nksic_rotwfn')
      !
      wfn2(:,iupdwn(isp):iupdwn(isp)-1+nupdwn(isp))=(0.d0,0.d0)

      !
      ! a blas could be used here XXX
      !
      do nbnd1=1,nupdwn(isp)
        do nbnd2=1,nupdwn(isp)
          wfn2(:,iupdwn(isp)-1 + nbnd1)=wfn2(:,iupdwn(isp)-1 + nbnd1) &
              + wfn1(:,iupdwn(isp)-1 + nbnd2) * Omat1(nbnd2,nbnd1)
        enddo
      enddo

      CALL stop_clock('nksic_rotwfn')
      !
      return
      !
!---------------------------------------------------------------
end subroutine nksic_rotwfn
!---------------------------------------------------------------




!-----------------------------------------------------------------------
      subroutine nksic_getHeigU(isp,vsicah,Heig,Umat)
!-----------------------------------------------------------------------
!
! ... solves the eigenvalues (Heig) and eigenvectors (Umat) of the force
!     matrix vsicah.
!     (Ultrasoft pseudopotential case is not implemented.)
!
      use kinds,                      only : dp
      use mp,                         only : mp_bcast
      use mp_global,                  only : intra_image_comm
      use io_global,                  only : ionode, ionode_id
      use electrons_base,             only : nupdwn
      !
      implicit none
      !
      ! in/out vars
      !
      integer,     intent(in)  :: isp
      real(dp)     :: Heig(nupdwn(isp))
      complex(dp)  :: Umat(nupdwn(isp),nupdwn(isp))
      real(dp)     :: vsicah(nupdwn(isp),nupdwn(isp))


      !
      ! local variables
      !
      complex(dp)              :: Hmat(nupdwn(isp),nupdwn(isp))
      complex(dp)              :: ci

      ci = (0.d0,1.d0)

!$$ Now this part diagonalizes Hmat = iWmat
      Hmat(:,:) = ci * vsicah(:,:)
!$$ diagonalize Hmat
!      if(ionode) then
      CALL zdiag(nupdwn(isp),nupdwn(isp),Hmat(1,1),Heig(1),Umat(1,1),1)
!      endif

!      CALL mp_bcast(Umat, ionode_id, intra_image_comm)
!      CALL mp_bcast(Heig, ionode_id, intra_image_comm)

      return
      !
!---------------------------------------------------------------
end subroutine nksic_getHeigU
!---------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine nksic_printoverlap(ninner,nouter)
!-----------------------------------------------------------------------
!
! ... Calculates the anti-hermitian part of the SIC hamiltonian, vsicah.
!
      use kinds,                      only : dp
      use grid_dimensions,            only : nr1x, nr2x, nr3x, nnrx
      use gvecw,                      only : ngw
      use mp,                         only : mp_sum
      use mp_global,                  only : intra_image_comm
      use io_global,                  only : ionode
      use electrons_base,             only : nbspx
      use cp_interfaces,              only : invfft
      use fft_base,                   only : dfftp
      use nksic,                      only : vsic
      use wavefunctions_module,       only : c0
      !
      implicit none
      !
      ! in/out vars
      !
      integer   :: ninner, nouter
      real(dp)  :: overlap(nbspx,nbspx),vsicah(nbspx,nbspx)

      !
      ! local variables
      !
      complex(dp)              :: psi1(nnrx), psi2(nnrx)
      real(dp)                 :: overlaptmp,vsicahtmp
      integer                  :: i,nbnd1,nbnd2
      real(dp)                 :: dwfnnorm

      dwfnnorm = 1.0/(dble(nr1x)*dble(nr2x)*dble(nr3x))

      vsicah(:,:) = 0.d0
      overlap(:,:) = 0.d0

      do nbnd1=1,nbspx
        CALL c2psi( psi1, nnrx, c0(:,nbnd1), (0.d0,0.d0), ngw, 1)
        CALL invfft('Dense', psi1, dfftp )

        do nbnd2=1,nbspx
          if(nbnd2.lt.nbnd1) then
            vsicahtmp = -vsicah(nbnd2,nbnd1)
            overlaptmp = overlap(nbnd2,nbnd1)
          else
            CALL c2psi( psi2, nnrx, c0(:,nbnd2), (0.d0,0.d0), ngw, 1)
            CALL invfft('Dense', psi2, dfftp )

            vsicahtmp = 0.d0
            overlaptmp = 0.d0

            do i=1,nnrx
!$$ Imposing Pederson condition
              vsicahtmp = vsicahtmp &
                  + 2.d0 * dble( conjg(psi1(i)) * (vsic(i,nbnd2)  &
                  - vsic(i,nbnd1) ) * psi2(i) ) * dwfnnorm
!$$ The following two lines give exactly the same results: checked
              overlaptmp = overlaptmp + dble( conjg(psi1(i)) * psi2(i) ) * dwfnnorm
!              overlaptmp = overlaptmp + dble(psi1(i)) * dble(psi2(i)) * dwfnnorm
            enddo

            CALL mp_sum(vsicahtmp,intra_image_comm)
            CALL mp_sum(overlaptmp,intra_image_comm)
          endif ! if(nbnd2.lt.nbnd1)

          vsicah(nbnd1,nbnd2) = vsicahtmp
          overlap(nbnd1,nbnd2) = overlaptmp

        enddo ! nbnd2=1,nbspx

      enddo ! nbspx

      if(ionode) then
        write(1021,*) ninner,nouter
        write(1022,*) ninner,nouter
        do nbnd1=1,nbspx
          write(1021,'(100F12.7)') (overlap(nbnd1,nbnd2),nbnd2=1,nbspx)
          write(1022,'(100F12.7)') (vsicah(nbnd1,nbnd2),nbnd2=1,nbspx)
        enddo
        write(1021,*)
        write(1022,*)
      endif

      return
      !
!---------------------------------------------------------------
end subroutine nksic_printoverlap
!---------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine nksic_getvsicah(isp,vsicah,vsicah2sum)
!-----------------------------------------------------------------------
!
! ... Calculates the anti-hermitian part of the SIC hamiltonian, vsicah.
!
      use kinds,                      only : dp
      use grid_dimensions,            only : nr1x, nr2x, nr3x, nnrx
      use gvecw,                      only : ngw
      use mp,                         only : mp_sum
      use mp_global,                  only : intra_image_comm
      use electrons_base,             only : nspin, iupdwn,nupdwn
      use cp_interfaces,              only : invfft
      use fft_base,                   only : dfftp
      use nksic,                      only : vsic,fsic
      use wavefunctions_module,       only : c0
      !
      implicit none
      !
      ! in/out vars
      ! 
      integer,     intent(in)  :: isp
      real(dp)  :: vsicah(nupdwn(isp),nupdwn(isp))
      !real(dp)  :: overlap(nupdwn(isp),nupdwn(isp))
      real(dp)  :: vsicah2sum

      !
      ! local variables
      !     
      complex(dp)              :: psi1(nnrx), psi2(nnrx)
      real(dp)                 :: vsicahtmp !,overlaptmp
      integer                  :: i,nbnd1,nbnd2
      real(dp)                 :: dwfnnorm
!      complex(dp)              :: Htest(2,2), Utest(2,2) 
!      real(dp)                 :: Eigtest(2)

    
      CALL start_clock('get_vsicah')
      !
      dwfnnorm = 1.0/(dble(nr1x)*dble(nr2x)*dble(nr3x))


      vsicah(:,:) = 0.d0
!      overlap(:,:) = 0.d0

      vsicah2sum = 0.d0

      do nbnd1=1,nupdwn(isp)
        CALL c2psi( psi1, nnrx, c0(:,iupdwn(isp)-1 + nbnd1), (0.d0,0.d0), ngw, 1)
        CALL invfft('Dense', psi1, dfftp )

        do nbnd2=1,nupdwn(isp)
          if(nbnd2.lt.nbnd1) then
            vsicahtmp = -vsicah(nbnd2,nbnd1)
!            overlaptmp = overlap(nbnd2,nbnd1)
          else
            CALL c2psi( psi2, nnrx, c0(:,iupdwn(isp)-1 + nbnd2), (0.d0,0.d0), ngw, 1)
            CALL invfft('Dense', psi2, dfftp )

            vsicahtmp = 0.d0
!            overlaptmp = 0.d0

            do i=1,nnrx
!$$ Imposing Pederson condition
              vsicahtmp = vsicahtmp &
                 + 2.d0 * dble( conjg(psi1(i)) * nspin * 0.5d0 &
                  * (vsic(i,iupdwn(isp)+nbnd2-1) * fsic(iupdwn(isp)+nbnd2-1) &
                  -  vsic(i,iupdwn(isp)+nbnd1-1) * fsic(iupdwn(isp)+nbnd1-1)) * psi2(i) ) * dwfnnorm
!                  + 2.d0 * dble( conjg(psi1(i)) * (vsic(i,iupdwn(isp)+nbnd2-1)  &
!                  - vsic(i,iupdwn(isp)+nbnd1-1) ) * psi2(i) ) * dwfnnorm
!              overlaptmp = overlaptmp + dble( conjg(psi1(i)) * psi2(i) ) * dwfnnorm
            enddo

! XXXX Use reciprocal space

            CALL mp_sum(vsicahtmp,intra_image_comm)
!            CALL mp_sum(overlaptmp,intra_image_comm)
          endif ! if(nbnd2.lt.nbnd1)

          vsicah(nbnd1,nbnd2) = vsicahtmp
!          overlap(nbnd1,nbnd2) = overlaptmp

          vsicah2sum = vsicah2sum + vsicahtmp*vsicahtmp
        enddo ! nbnd2=1,nupdwn(isp)

      enddo ! nbnd1=1,nupdwn(isp)

      CALL stop_clock('get_vsicah')
      !
      return
      !
!---------------------------------------------------------------
end subroutine nksic_getvsicah
!---------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine nksic_getOmat1(isp,Heig,Umat,passof,Omat1)
!-----------------------------------------------------------------------
!
! ... Obtains the rotation matrix from the force-related matrices Heig and Umat
!     and also from the step size (passof).
!
      use kinds,                      only : dp
      use constants,                  only : ci 
      use electrons_base,             only : nupdwn
      !
      implicit none
      !
      ! in/out vars
      !
      integer,      intent(in) :: isp
      real(dp),     intent(in) :: Heig(nupdwn(isp))
      complex(dp),  intent(in) :: Umat(nupdwn(isp),nupdwn(isp))
      real(dp),     intent(in) :: passof
      real(dp)                 :: Omat1(nupdwn(isp),nupdwn(isp))

      !
      ! local variables
      !
      complex(dp) :: Cmattmp(nupdwn(isp),nupdwn(isp))
      complex(dp) :: exp_iHeig(nupdwn(isp))

      integer     :: nbnd1
      real(dp)    :: dtmp


!$$ We set the step size in such a way that the phase change
!$$ of the wavevector with largest eigenvalue upon rotation is fixed
!          passof = passoprod/max(abs(Heig(1)),abs(Heig(nupdwn(isp))))
!$$ Now the above step is done outside.
            
          do nbnd1=1,nupdwn(isp)
            dtmp =  passof * Heig(nbnd1)
            exp_iHeig(nbnd1) = cos(dtmp) + ci*sin(dtmp)
          enddo
            
!$$ Cmattmp = exp(i * passof * Heig) * Umat^dagger   ; Omat = Umat * Cmattmp
          do nbnd1=1,nupdwn(isp)
            Cmattmp(:,nbnd1) = Umat(:,nbnd1) * exp_iHeig(nbnd1)
          enddo

          Omat1 = dble ( MATMUL(Cmattmp, conjg(transpose(Umat)) ) )

      return
!---------------------------------------------------------------
end subroutine nksic_getOmat1
!---------------------------------------------------------------

!$$
