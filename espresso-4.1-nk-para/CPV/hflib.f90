!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Hartree-Fock method
! Implemented by I. Dabo (Universite Paris-Est, Ecole des Ponts, ParisTech)
! Optimized and Parallelized by Andrea Ferretti (MIT)
!
!-----------------------------------------------------------------------
      subroutine calc_hf_potential(c, rhor, rhog, vxxpsi, exx, detothf)
!-----------------------------------------------------------------------
!
! ... calculate orbital densities and Hartree-Fock potentials
! ... (the calculation of the exchange potential uses
! ... periodic-image corrections)
!
      use kinds,                   only: dp
      use gvecp,                   only: ngm
      use gvecs,                   only: ngs, nps, nms
      use gvecw,                   only: ngw
      use recvecs_indexes,         only: np, nm
      use reciprocal_vectors,      only: gstart
      use grid_dimensions,         only: nr1, nr2, nr3, &
                                         nr1x, nr2x, nr3x, nnrx
      use cell_base,               only: omega, a1, a2, a3
      use smooth_grid_dimensions,  only: nr1s, nr2s, nr3s, &
                                         nr1sx, nr2sx, nr3sx, nnrsx
      use electrons_base,          only: nx => nbspx, nbsp, &
                                         f, ispin, nspin, nelt, &
                                         nupdwn, iupdwn
      use constants,               only: pi, fpi
      use mp,                      only: mp_sum
      use mp_global,               only: intra_image_comm
      use io_global,               only: stdout, ionode
      use funct,                   only: dft_is_gradient
      use cp_interfaces,           only: fwfft, invfft, fillgrad
      use fft_base,                only: dffts, dfftp
      use hfmod,                   only: hfscalfact
      !
      implicit none
      !
      ! I/O vars
      !
      complex(dp), intent(in)  :: c(ngw,nx)
      real(dp),    intent(in)  :: rhor(nnrx,nspin)
      complex(dp), intent(out) :: rhog(ngm,nspin)
      complex(dp), intent(out) :: vxxpsi(ngw,nx)
      real(dp),    intent(out) :: exx(nx)
      real(dp),    intent(out) :: detothf
      !
      ! local vars
      !
      real(dp), parameter :: vanishing_f = 1.e-6_dp
      !
      integer     :: iss, isup, isdw, iss1, iss2 
      integer     :: ios, i, ir, ig, k, j
      integer     :: istart, iend
      real(dp)    :: sa1, sa2, eaux(2), etxc, fact
      real(dp)    :: faux_i, faux_j0, faux_jp
      complex(dp) :: ci,fp,fm
      !
      character(6),   external :: int_to_char
      complex(dp), allocatable :: psi(:), psis1(:), psis2(:), &
                                  vxxs(:), vxxd(:), &
                                  vxxpsis(:), psis(:)
      complex(dp), allocatable :: orbitalrhog(:,:)
      real(dp),    allocatable :: orbitalrhos(:,:)
      real(dp),    allocatable :: orbitalrhor(:,:)
      real(dp),    allocatable :: aux(:,:)
      real(dp),    allocatable :: vxc(:,:)
      real(dp),    allocatable :: grhor(:,:,:)
      real(dp),    allocatable :: h(:,:,:)
      !
      ci=(0.0d0,1.0d0)
      !
      ! ... calculate the density of each individual orbital
      !
      allocate(psis1(nnrsx)) 
      allocate(psis2(nnrsx)) 
      allocate(vxc(nnrx,nspin))
      !
      if ( dft_is_gradient() ) then
          !
          allocate(grhor(nnrx,3,nspin))
          allocate(h(nnrx,nspin,nspin))
          !
          call fillgrad( nspin, rhog, grhor)     
      else
          !
          allocate(grhor(1,1,1))
          allocate(h(1,1,1))
          !
          grhor=0.0_dp
          h=0.0_dp
          !
      endif
      !
      fact=omega/dble(nr1*nr2*nr3)
      !
      vxxpsi=0.0_dp
      exx=0.0_dp

      !
      ! the calculation of the total charge density 
      ! has been removed, now it comes from input
      !
      vxc=0.0_dp
      etxc=0.0_dp
      !
      call exch_corr_wrapper(nnrx,nspin,grhor,rhor,etxc,vxc,h)
      detothf=-etxc*fact
      !
      deallocate(grhor)
      deallocate(h)
      !
      call mp_sum(detothf,intra_image_comm)
      !
  
      !
      ! main loop over states
      !
      outer_loop: &
      do i = 1, nbsp

          !
          ! psi_i on the smooth grid
          !
          psis1=0.d0
          do ig=1,ngw
              psis1(nms(ig))=conjg(c(ig,i))
              psis1(nps(ig))=c(ig,i)
          enddo
          call invfft('Wave',psis1,dffts)


          !
          ! inner loop
          !
          if ( nspin == 1 ) then
             istart = 1
             iend   = nbsp
          else
             istart = iupdwn( ispin(i) )
             iend   = iupdwn( ispin(i) ) + nupdwn( ispin(i) ) -1
          endif
          !          
          inner_loop: &
          do j  = istart, iend, 2
              !
              ! if( ispin(i) /= ispin(j) ) cycle
              !
              ! take into account spin multiplicity
              !
              faux_i  = f(i)   * dble( nspin ) / 2.0_dp 
              faux_j0 = f(j)   * dble( nspin ) / 2.0_dp 
              !
              if ( j+1 <= iend ) then
                 faux_jp = f(j+1) * dble( nspin ) / 2.0_dp 
              else
                 faux_jp = 0
              endif
              !
              !if ( faux_j < 1.0e-6 ) cycle
              !if ( faux_i < 1.0e-6 ) cycle

              !
              ! allocate mem
              !
              allocate(orbitalrhog(ngm,2))
              allocate(orbitalrhos(nnrsx,2))
              allocate(orbitalrhor(nnrx,2))
              !

              orbitalrhog(1:ngw,1) = c(1:ngw, j)
              !
              if ( j+1 <= iend ) then
                  orbitalrhog(1:ngw,2) = c(1:ngw, j+1)
              else
                  orbitalrhog(1:ngw,2) = 0.0_dp
              endif
              
              !
              ! psi_j's on the smooth grid
              !
              call c2psi( psis2, nnrsx, orbitalrhog(:,1), &
                                        orbitalrhog(:,2), ngw, 2)
              !
              call invfft('Wave', psis2, dffts)
          

              !
              ! psi_i * psi_j on the smooth grid
              !
              sa1=1.d0/omega
              orbitalrhos=0.d0
              !
              do ir=1,nnrsx
                  orbitalrhos(ir,1) = sa1*dble(psis1(ir))*dble(psis2(ir))
                  orbitalrhos(ir,2) = sa1*dble(psis1(ir))*aimag(psis2(ir))
              enddo


              !
              ! move to the dense grid
              !
              if ( nnrsx == nnrx ) then
                  !
                  orbitalrhor(:,:) = orbitalrhos(:,:)
                  !
              else
                  !
                  ! transform to reciprocal space
                  !
                  allocate(psis(nnrsx)) 
                  !
                  psis=0.0_dp
                  do ir=1,nnrsx
                      psis(ir) = cmplx( orbitalrhos(ir,1), &
                                        orbitalrhos(ir,2)  )
                  enddo
                  !
                  call fwfft('Smooth',psis,dffts)
                  orbitalrhog=(0.d0,0.d0)
                  do ig=1,ngs
                      fp=psis(nps(ig))+psis(nms(ig))
                      fm=psis(nps(ig))-psis(nms(ig))
                      orbitalrhog(ig,1)=0.5d0*cmplx(dble(fp),aimag(fm))
                      orbitalrhog(ig,2)=0.5d0*cmplx(aimag(fp),-dble(fm))
                  enddo
                  !
                  ! switch to dense grid in real space
                  !
                  allocate(psi(nnrx))
                  !
                  psi(:)=(0.d0,0.d0)
                  do ig=1,ngs
                      psi(nm(ig)) = conjg(orbitalrhog(ig,1)) &
                                  + ci*conjg(orbitalrhog(ig,2))
                      psi(np(ig)) = orbitalrhog(ig,1)&
                                  + ci*orbitalrhog(ig,2)
                  enddo
                  !
                  call invfft('Dense',psi,dfftp)
                  do ir=1,nnrx
                      orbitalrhor(ir,1) =  dble(psi(ir))
                      orbitalrhor(ir,2) = aimag(psi(ir))
                  enddo
                  !
                  deallocate(psi) 
                  deallocate(psis)
                  !
              endif
              !
              deallocate(orbitalrhog)
              deallocate(orbitalrhos)
              !
              !!
              !! The following is left here for DEBUG purposes only,
              !! it is definitely too bad for the performance in real calculations
              !!
              !if(i.eq.j) then
              !  call writetofile(orbitalrhor(:,1),nnrx, &
              !              'rhoup'//TRIM(int_to_char(i))//'.dat',dfftp,'az')
              !  call writetofile(orbitalrhor(:,2),nnrx, &
              !              'rhodw'//TRIM(int_to_char(i))//'.dat',dfftp,'az')
              !end if
    
              !
              ! calculate exchange contribution 
              !
              allocate( aux( nnrx,2) )
              !
              call hf_correction( orbitalrhor(:,1), aux(:,1), eaux(1) )
              call hf_correction( orbitalrhor(:,2), aux(:,2), eaux(2) )
              !
              deallocate(orbitalrhor)
              allocate(vxxd(nnrx))
              !
              vxxd(1:nnrx) = cmplx( faux_j0 * aux(1:nnrx,1), &
                                    faux_jp * aux(1:nnrx,2)  )
              !
              exx(i)=exx(i) - 0.5_dp*faux_i*faux_j0*eaux(1)  &
                            - 0.5_dp*faux_i*faux_jp*eaux(2)
              !
              deallocate(aux)
              !
              !
              if ( i == j ) then
                  vxxd(1:nnrx) = vxxd(1:nnrx) + vxc(1:nnrx,ispin(i))
              elseif ( i == j+1 ) then
                  vxxd(1:nnrx) = vxxd(1:nnrx) + ci * vxc(1:nnrx,ispin(i))
              endif

              !
              ! change grid if the case
              !
              allocate(vxxs(nnrsx))
              !
              if ( nnrsx == nnrx ) then
                  !
                  vxxs(1:nnrsx) = vxxd(1:nnrx)
                  !
              else
                  !
                  call fwfft('Dense',vxxd,dfftp )
                  !
                  vxxs=0.0_dp
                  do ig=1,ngs
                      vxxs(ig) = vxxd(np(ig))
                  enddo
                  !
                  call invfft('Smooth',vxxs,dffts)
                  !
              endif
              !
              deallocate(vxxd)

                 
              allocate(vxxpsis(nnrsx))
              !
              vxxpsis=0.0_dp
              do ir=1,nnrsx
                  vxxpsis(ir)= cmplx( dble(vxxs(ir)) * dble(psis2(ir)), & 
                                     aimag(vxxs(ir)) *aimag(psis2(ir))  )
              enddo
              call fwfft('Wave',vxxpsis,dffts)
              !
              deallocate(vxxs)
              !
              do ig=1,ngw
                  !
                  fp = vxxpsis(nps(ig)) +vxxpsis(nms(ig))
                  fm = vxxpsis(nps(ig)) -vxxpsis(nms(ig))
                  !
                  vxxpsi(ig,i) = vxxpsi(ig,i) &
                               - 0.5_dp * cmplx( dble(fp),aimag(fm) ) &
                               - 0.5_dp * cmplx( aimag(fp),-dble(fm) )
              enddo
              !
              deallocate(vxxpsis)
              !
          enddo inner_loop
          !
          if ( nspin == 1 ) exx(i) = 2.0_dp * exx(i)
          !
      enddo outer_loop
      !
      vxxpsi  = hfscalfact *vxxpsi
      detothf = hfscalfact *detothf


      deallocate(psis1) 
      deallocate(psis2) 
      !
      deallocate(vxc)
      return
      !
!-----------------------------------------------------------------------
      end subroutine calc_hf_potential
!-----------------------------------------------------------------------

!---------------------------------------------------------------
      SUBROUTINE hf_correction(rhoele,vxx,exx) 
!---------------------------------------------------------------
      !
      ! ... calculate Hartree-Fock potential from orbital density
      !
      use kinds,                only : dp
      use constants,            only : e2, fpi
      use cell_base,            only : tpiba2,omega
      use grid_dimensions,      only : nnrx, nr1, nr2, nr3
      use gvecp,                only : ngm
      use recvecs_indexes,      only : np, nm
      use reciprocal_vectors,   only : gstart, g
      use eecp_mod,             only : do_comp
      use cp_interfaces,        only : fwfft, invfft
      use fft_base,             only : dfftp
      use mp,                   only : mp_sum
      use mp_global,            only : intra_image_comm
      !
      implicit none
      !
      real(dp),    intent(in) :: rhoele(nnrx)
      real(dp),   intent(out) :: vxx(nnrx)
      real(dp),   intent(out) :: exx
      !
      integer  :: i, is, ig, ir
      real(dp) :: fact
      complex(dp),allocatable :: aux(:)
      complex(dp),allocatable :: rhotmp(:)
      complex(dp),allocatable :: vtemp(:)
      complex(dp),allocatable :: vcorrtmp(:)
      !
      allocate(rhotmp(ngm))
      allocate(vtemp(ngm))
      allocate(vcorrtmp(ngm))
      allocate(aux(nnrx))
      !
      ! ... compute self-hartree potential
      !
      aux(:)=rhoele(:)
      call fwfft('Dense',aux,dfftp)
      !
      rhotmp=0.0_dp
      do ig=1,ngm
          rhotmp(ig) = aux(np(ig))
      enddo
      !
      if(gstart==2) vtemp(1)=(0.d0,0.d0)
      DO ig=gstart,ngm
        vtemp(ig)=rhotmp(ig)*fpi/(tpiba2*g(ig))
      END DO
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
      vxx=dble(aux)
      !
      fact=omega/dble(nr1*nr2*nr3)
      exx=sum(vxx(1:nnrx)*rhoele(1:nnrx))*fact
      !
      CALL mp_sum(exx,intra_image_comm)
      !
      deallocate(rhotmp)
      deallocate(vtemp)
      deallocate(vcorrtmp)
      deallocate(aux)
      !
      return
      !
!---------------------------------------------------------------
      end subroutine hf_correction
!---------------------------------------------------------------


