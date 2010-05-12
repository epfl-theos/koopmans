!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Hartree-Fock method
! Implemented by I. Dabo (Universite Paris-Est, Ecole des Ponts, ParisTech)
! Parallelized by Andrea Ferretti (MIT)
!
!-----------------------------------------------------------------------
      subroutine calc_hf_potential(c,vxxpsi,exx,detothf)
!-----------------------------------------------------------------------
!
! ... calculate orbital densities and Hartree-Fock potentials
! ... (the calculation of the exchange potential uses
! ... periodic-image corrections)
!
      use kinds,              only: dp
      use gvecp,              only: ngm
      use gvecs,              only: ngs, nps, nms
      use gvecw,              only: ngw
      use recvecs_indexes,    only: np, nm
      use reciprocal_vectors, only: gstart
      use grid_dimensions,    only: nr1, nr2, nr3, &
                                    nr1x, nr2x, nr3x, nnrx
      use cell_base,          only: omega, a1, a2, a3
      use smooth_grid_dimensions, only: nr1s, nr2s, nr3s, &
                                        nr1sx, nr2sx, nr3sx, nnrsx
      use electrons_base,     only: nx => nbspx, n => nbsp, f, ispin, &
                                    nspin, nelt
      use constants,          only: pi, fpi
      use mp,                 only: mp_sum
      use mp_global,          only : intra_image_comm
      use io_global,          only: stdout, ionode
      use cp_interfaces,      only: fwfft, invfft
      use fft_base,           only: dffts, dfftp
      use hfmod,              only: hfscalfact
      !
      implicit none
      !
      complex(dp) :: c(ngw,nx)
      complex(dp) :: vxxpsi(ngw,nx)
      real(dp)    :: exx(nx)
      real(dp)    :: detothf
      !
      integer  :: iss, isup, isdw, iss1, iss2, ios, i, ir, ig, k, j
      real(dp) :: sa1, sa2, eaux, etxc, fact
      real(dp), parameter :: vanishing_f = 1.e-6_dp
      complex(dp) :: ci,fp,fm
      character(len=6), external :: int_to_char
      complex(dp), allocatable :: psi(:), psis1(:), psis2(:), &
                                  rhog(:,:), &
                                  vxxs(:), vxxg(:), vxxpsis(:), psis(:)
      complex(dp), allocatable :: orbitalrhog(:,:)
      real(dp), allocatable :: orbitalrhos(:,:)
      real(dp), allocatable :: orbitalrhor(:,:)
      real(dp), allocatable :: aux(:),rhor(:,:),rhos(:,:)
      real(dp), allocatable :: vxc(:,:)
      real(dp), allocatable :: grhor(:,:,:)
      real(dp), allocatable :: h(:,:,:)
      !
      ci=(0.0d0,1.0d0)
      !
      ! ... calculate the density of each individual orbital
      !
      allocate(orbitalrhog(ngm,2))
      allocate(orbitalrhos(nnrsx,2))
      allocate(orbitalrhor(nnrx,2))
      allocate(aux(nnrx))
      allocate(psis1(nnrsx)) 
      allocate(psis2(nnrsx)) 
      allocate(vxxg(nnrx))
      allocate(vxxs(nnrsx))
      allocate(vxxpsis(nnrsx))
      allocate(vxc(nnrx,2))
      allocate(grhor(nnrx,3,2))
      allocate(h(nnrx,2,2))
      allocate(rhos(nnrsx,2))
      allocate(rhor(nnrx,2))
      allocate(rhog(ngm,2))
      !
      fact=omega/dble(nr1*nr2*nr3)
      !
      vxxpsi=0.0_dp
      exx=0.0_dp
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
      grhor=0.0_dp
      h=0.0_dp
      vxc=0.0_dp
      etxc=0.0_dp
      !
      call exch_corr_wrapper(nnrx,2,grhor,rhor,etxc,vxc,h)
      detothf=-etxc*fact
      !
      call mp_sum(detothf,intra_image_comm)
      !
      do i=1,n
        !
        ! psi_i on the smooth grid
        !
        psis1=0.d0
        do ig=1,ngw
          psis1(nms(ig))=conjg(c(ig,i))
          psis1(nps(ig))=c(ig,i)
        end do
        call invfft('Wave',psis1,dffts)
        !
        do j=1,n
          !
          if(ispin(i).ne.ispin(j)) cycle
          !
          ! psi_j on the smooth grid
          !
          psis2=0.d0
          do ig=1,ngw
            psis2(nms(ig))=conjg(c(ig,j))
            psis2(nps(ig))=c(ig,j)
          end do
          call invfft('Wave',psis2,dffts)
          !
          ! psi_i * psi_j on the smooth grid
          !
          iss1=ispin(i)
          sa1=1.d0/omega
          orbitalrhos=0.d0
          do ir=1,nnrsx
            orbitalrhos(ir,iss1)=sa1*dble(psis1(ir))*dble(psis2(ir))
          end do
          !
          ! transform to reciprocal space
          !
          allocate(psis(nnrsx)) 
          isup=1
          isdw=2
          psis=0.0_dp
          do ir=1,nnrsx
            psis(ir)=cmplx(orbitalrhos(ir,isup),orbitalrhos(ir,isdw))
          end do
          call fwfft('Smooth',psis,dffts)
          orbitalrhog=(0.d0,0.d0)
          do ig=1,ngs
            fp=psis(nps(ig))+psis(nms(ig))
            fm=psis(nps(ig))-psis(nms(ig))
            orbitalrhog(ig,isup)=0.5d0*cmplx(dble(fp),aimag(fm))
            orbitalrhog(ig,isdw)=0.5d0*cmplx(aimag(fp),-dble(fm))
          end do
          !
          ! switch to dense grid in real space
          !
          allocate(psi(nnrx))
          isup=1
          isdw=2
          psi(:)=(0.d0,0.d0)
          do ig=1,ngs
            psi(nm(ig))=conjg(orbitalrhog(ig,isup)) &
                       +ci*conjg(orbitalrhog(ig,isdw))
            psi(np(ig))=orbitalrhog(ig,isup)+ci*orbitalrhog(ig,isdw)
          end do
          call invfft('Dense',psi,dfftp)
          do ir=1,nnrx
            orbitalrhor(ir,isup)= dble(psi(ir))
            orbitalrhor(ir,isdw)=aimag(psi(ir))
          end do
          !
          deallocate(psi) 
          deallocate(psis)
          !
          if(i.eq.j) then
            call writetofile(orbitalrhor(:,1),nnrx, &
                        'rhoup'//TRIM(int_to_char(i))//'.dat',dfftp,'az')
            call writetofile(orbitalrhor(:,2),nnrx, &
                        'rhodw'//TRIM(int_to_char(i))//'.dat',dfftp,'az')
          end if
          !
          ! calculate exchange contribution 
          !
          call hf_correction(ispin(i),orbitalrhor,aux,eaux)
          !
          vxxg(1:nnrx)=f(j)*aux(1:nnrx)
          if(i.eq.j) vxxg(1:nnrx)=vxxg(1:nnrx)+vxc(1:nnrx,ispin(i))
          call fwfft('Dense',vxxg,dfftp )
          vxxs=0.0_dp
          do ig=1,ngs
            vxxs(nps(ig))=vxxg(np(ig))
            vxxs(nms(ig))=conjg(vxxg(np(ig)))
          end do
          call invfft('Smooth',vxxs,dffts)
          !
          vxxpsis=0.0_dp
          do ir=1,nnrsx
            vxxpsis(ir)=cmplx(dble(vxxs(ir))*dble(psis2(ir)),0.0_dp)
          end do
          call fwfft('Wave',vxxpsis,dffts)
          !
          do ig=1,ngw
            vxxpsi(ig,i)=vxxpsi(ig,i)-vxxpsis(nps(ig))
          end do
          !
          exx(i)=exx(i)-0.5_dp*f(i)*f(j)*eaux
          !
        end do
        !
      end do
      !
      vxxpsi=hfscalfact*vxxpsi
      detothf=hfscalfact*detothf
      !
      deallocate(orbitalrhog)
      deallocate(orbitalrhos)
      deallocate(orbitalrhor)
      deallocate(aux)
      deallocate(psis1) 
      deallocate(psis2) 
      deallocate(vxxg)
      deallocate(vxxs)
      deallocate(vxxpsis)
      deallocate(vxc)
      deallocate(grhor)
      deallocate(h)
      deallocate(rhos)
      deallocate(rhor)
      deallocate(rhog)
      !
      return
      !
!-----------------------------------------------------------------------
      end subroutine calc_hf_potential
!-----------------------------------------------------------------------
!---------------------------------------------------------------
      SUBROUTINE hf_correction(ispin,rhoele,vxx,exx) 
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
      real(dp),    intent(in) :: rhoele(nnrx,2)
      integer,     intent(in) :: ispin
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
      rhotmp=0.0_dp
      do is=1,2
        aux(:)=rhoele(:,is)
        call fwfft('Dense',aux,dfftp)
        do ig=1,ngm
          rhotmp(ig)=rhotmp(ig)+aux(np(ig))
        end do
      end do
      !
      if(gstart==2) vtemp(1)=(0.d0,0.d0)
      DO ig=gstart,ngm
        vtemp(ig)=rhotmp(ig)*fpi/(tpiba2*g(ig))
      END DO
      !
      if(do_comp) then
        call calc_tcc_potential(vcorrtmp,rhotmp)
        vtemp=vtemp+vcorrtmp
      end if
      !
      aux=0.0_dp
      do ig=1,ngm
        aux(np(ig))=vtemp(ig)
        aux(nm(ig))=conjg(vtemp(ig))
      end do
      call invfft('Dense',aux,dfftp)
      !
      vxx=dble(aux)
      !
      fact=omega/dble(nr1*nr2*nr3)
      exx=sum(vxx(1:nnrx)*rhoele(1:nnrx,ispin))*fact
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


