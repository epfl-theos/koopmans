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
      subroutine hf_potential(nbsp1, nx1, c1, f1, ispin1, iupdwn1, nupdwn1, &
                              nbsp2, nx2, c2, f2, ispin2, iupdwn2, nupdwn2, &
                              rhor, rhog, vxxpsi, exx )
!-----------------------------------------------------------------------
!
! Calculate orbital densities and Hartree-Fock potentials
! (the calculation of the exchange potential uses
! periodic-image corrections)
!
! 1 -> inner quantities (those used to build the density matrix)
! 2 -> outer quantities (those referring to the wfcs to which the
!              HF potential is applied)
!
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
      use electrons_base,          only: nspin
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
      integer,     intent(in)  :: nbsp1, nx1
      complex(dp), intent(in)  :: c1(ngw,nx1)
      integer,     intent(in)  :: ispin1(nx1)
      integer,     intent(in)  :: iupdwn1(nspin), nupdwn1(nspin)
      real(dp),    intent(in)  :: f1(nx1)
      !
      integer,     intent(in)  :: nbsp2, nx2
      complex(dp), intent(in)  :: c2(ngw,nx2)
      integer,     intent(in)  :: ispin2(nx2)
      integer,     intent(in)  :: iupdwn2(nspin), nupdwn2(nspin)
      real(dp),    intent(in)  :: f2(nx2)
      !
      real(dp),    intent(in)  :: rhor(nnrx,nspin)
      complex(dp), intent(in)  :: rhog(ngm,nspin)
      complex(dp), intent(out) :: vxxpsi(ngw,nx2)
      real(dp),    intent(out) :: exx(nx2)
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
      CALL start_clock( 'hf_potential' )
      !
      ! ... local workspace
      !
      allocate(psis1(nnrsx)) 
      allocate(psis2(nnrsx)) 
      !
      vxxpsi=0.0_dp
      exx=0.0_dp

  
      !
      ! main loop over states
      !
      outer_loop: &
      do i = 1, nbsp2

          !
          ! psi_i on the smooth grid
          !
          psis2=0.d0
          do ig=1,ngw
              psis2(nms(ig))=conjg(c2(ig,i))
              psis2(nps(ig))=c2(ig,i)
          enddo
          call invfft('Wave',psis2,dffts)


          !
          ! inner loop
          !
          if ( nspin == 1 ) then
             istart = 1
             iend   = nbsp1
          else
             istart = iupdwn1( ispin2(i) )
             iend   = iupdwn1( ispin2(i) ) + nupdwn1( ispin2(i) ) -1
          endif
          !          
          inner_loop: &
          do j  = istart, iend, 2
              !
              ! if( ispin(i) /= ispin(j) ) cycle
              !
              ! take into account spin multiplicity
              !
              faux_i  = f2(i)   * dble( nspin ) / 2.0_dp 
              faux_j0 = f1(j)   * dble( nspin ) / 2.0_dp 
              !
              if ( j+1 <= iend ) then
                 faux_jp = f1(j+1) * dble( nspin ) / 2.0_dp 
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

              orbitalrhog(1:ngw,1) = c1(1:ngw, j)
              !
              if ( j+1 <= iend ) then
                  orbitalrhog(1:ngw,2) = c1(1:ngw, j+1)
              else
                  orbitalrhog(1:ngw,2) = 0.0_dp
              endif
              
              !
              ! psi_j's on the smooth grid
              !
              call c2psi( psis1, nnrsx, orbitalrhog(:,1), &
                                        orbitalrhog(:,2), ngw, 2)
              !
              call invfft('Wave', psis1, dffts)
          

              !
              ! psi_i * psi_j on the smooth grid
              ! psis2 <=  psi_i
              ! psis1 <=  psi_j, psi_j+1
              !
              sa1=1.d0/omega
              orbitalrhos=0.d0
              !
              do ir=1,nnrsx
                  orbitalrhos(ir,1) = sa1*dble(psis2(ir))*dble(psis1(ir))
                  orbitalrhos(ir,2) = sa1*dble(psis2(ir))*aimag(psis1(ir))
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
                  vxxpsis(ir)= cmplx( dble(vxxs(ir)) * dble(psis1(ir)), & 
                                     aimag(vxxs(ir)) *aimag(psis1(ir))  )
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
      exx     = hfscalfact *exx

      deallocate(psis2) 
      deallocate(psis1) 
      !
      CALL stop_clock( 'hf_potential' )
      !
      return
      !
!-----------------------------------------------------------------------
      end subroutine hf_potential
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
      

      CALL start_clock( 'hf_corr' )
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
      CALL stop_clock( 'hf_corr' )
      return
      !
!---------------------------------------------------------------
      end subroutine hf_correction
!---------------------------------------------------------------


