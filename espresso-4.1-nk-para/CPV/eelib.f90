!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Electrostatic embedding methods 
! Developed and implemented by I. Dabo (Universite Paris-Est, Ecole des Ponts, ParisTech)
! Parallelized by Andrea Ferretti (MIT)
!
!-----------------------------------------------------------------------
      subroutine ee_green_0d_init(box)
!-----------------------------------------------------------------------
!
! ... initialize Green's functions for periodic-image correction
! ... in 0D setttings (e.g., isolated molecule, cluster)
!
      use kinds,              only : dp
      use cell_base,          only : a1, a2, a3, omega, tpiba2, s_to_r, &
                                     boxdimensions
      use constants,          only : fpi, pi
      use io_global,          only : stdout
      use grid_dimensions,    only : nnrx, nr1,  nr2,  nr3, nr1x, nr2x, &
                                     nr3x, nr1l, nr2l, nr3l
      use recvecs_indexes,    only : np, nm
      use reciprocal_vectors, only : gstart, g
      use gvecp,              only : ngm
      use fft_base,           only : dfftp
      use cp_interfaces,      only : fwfft, invfft
      use eecp_mod,           only : gcorr,gcorr_fft
      use mp_global,          only : me_image
      !
      implicit none
      !
      type(boxdimensions), intent(in) :: box
      !
      real(dp),      parameter :: sigma=2.0_dp !afcmodified:giovanni 2.d0
      real(dp),      parameter :: vanishing_dist=1.0e-3_dp
      !
      complex(dp), allocatable :: vtemp(:)
      real(dp),    allocatable :: vtempr(:)
      real(dp) :: aux(dfftp%nr1,dfftp%nr2,dfftp%nr3)
      !
      integer             :: ig, ir1, ir2, ir3, ir, i, j, k
      real(dp)            :: sv(3), lv(3) ,dist
      real(dp)            :: a(3,3)
      integer             :: npt(3)
      logical             :: tperiodic(3)
      real(dp),  external :: qe_erf
      !
      interface 
        !
        function afc(a,npt,tperiodic,spreadopt)
          !
          real(8), intent(in), optional :: spreadopt
          real(8), intent(in), dimension(3,3) :: a
          integer, intent(in), dimension(3) :: npt
          logical, intent(in), dimension(3) :: tperiodic
          real(8) :: afc(npt(1),npt(2),npt(3))
          !
       end function
       !
      end interface
      !
      ! main body
      !
      allocate(vtemp(nnrx))
      allocate(vtempr(nnrx))
      ! 
      vtemp=0.0_dp
      vtempr=0.0_dp
      gcorr=0.0_dp
      !
      ir1=1
      ir2=1
      ir3=1
      do k=1,me_image
        ir3=ir3+dfftp%npp(k)
      enddo      
      !
      npt(1)=dfftp%nr1
      npt(2)=dfftp%nr2
      npt(3)=dfftp%nr3
      tperiodic=.false.
      a(1:3,1)=a1(1:3)
      a(1:3,2)=a2(1:3)
      a(1:3,3)=a3(1:3)
      aux=afc(a,npt,tperiodic,sigma)
      !
      do k=1,nr3l
        do j=1,nr2l
          do i=1,nr1l
            !
            ir=i+(j-1)*dfftp%nr1x+(k-1)*dfftp%nr1x*dfftp%nr2x
            !
            gcorr(ir)=aux(i+ir1-1,j+ir2-1,k+ir3-1)
            !
          end do
        end do
      end do
      !call writetofile(gcorr,nnrx,'afc0d.dat',dfftp,'az')
      vtemp=gcorr
      call fwfft('Dense',vtemp,dfftp)
      do ig=1,ngm
        gcorr_fft(ig)=vtemp(np(ig))
      enddo
      !
      deallocate(vtempr)
      deallocate(vtemp)
      !
      return
      !
!-----------------------------------------------------------------------
      end subroutine ee_green_0d_init
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine ee_green_1d_init(box)
!-----------------------------------------------------------------------
!
! ... initialize Green's functions for periodic-image correction
! ... for 1D setttings (e.g., nanotube, polymer chain)
!
      use kinds,              only : dp
      use cell_base,          only : a1, a2, a3, omega, tpiba2, s_to_r, &
                                     boxdimensions
      use constants,          only : fpi, pi
      use io_global,          only : stdout
      use grid_dimensions,    only : nnrx, nr1, nr2, nr3, nr1x, nr2x, nr3x, &
                                           nr1l, nr2l, nr3l
      use recvecs_indexes,    only : np, nm
      use reciprocal_vectors, only : gstart, g, gx
      use gvecp,              only : ngm
      use fft_base,           only : dfftp
      use cp_interfaces,      only : fwfft, invfft
      use eecp_mod,           only : gcorr1d,gcorr1d_fft
      use mp_global,          only : me_image
      !
      implicit none
      !
      type(boxdimensions), intent(in) :: box
      !
      real(dp),      parameter :: sigma=2.0_dp
      real(dp),      parameter :: vanishing_dist=1.0e-3_dp
      real(dp),      parameter :: vanishing_g=1.0e-3_dp
      real(dp),      parameter :: euler_gamma=0.57721566490153286061d0
      !
      complex(dp), allocatable :: vtemp(:)
      real(dp),    allocatable :: vtempr(:)
      real(dp) :: aux(dfftp%nr1,dfftp%nr2,dfftp%nr3)
      !
      integer                  :: ig, ir1, ir2, ir3, ir, i, j, k
      real(dp)                 :: sv(3), lv(3), dist
      real(dp)            :: a(3,3)
      integer             :: npt(3)
      logical             :: tperiodic(3)
      !
      real(dp),       external :: qe_erf
      real(dp),       external :: eimlmg
      !
      interface 
        !
        function afc(a,npt,tperiodic,spreadopt)
          !
          real(8), intent(in), optional :: spreadopt
          real(8), intent(in), dimension(3,3) :: a
          integer, intent(in), dimension(3) :: npt
          logical, intent(in), dimension(3) :: tperiodic
          real(8) :: afc(npt(1),npt(2),npt(3))
          !
        end function
        !
      end interface
      !
      ! main body
      !
      allocate(vtemp(nnrx))
      allocate(vtempr(nnrx))
      ! 
      vtemp=0.0_dp
      vtempr=0.0_dp
      gcorr1d=0.0_dp
      !
      ir1=1
      ir2=1
      ir3=1
      do k=1,me_image
        ir3=ir3+dfftp%npp(k)
      enddo
      !      
      npt(1)=dfftp%nr1
      npt(2)=dfftp%nr2
      npt(3)=dfftp%nr3
      tperiodic(1)=.false.
      tperiodic(2)=.false.
      tperiodic(3)=.true.
      a(1:3,1)=a1(1:3)
      a(1:3,2)=a2(1:3)
      a(1:3,3)=a3(1:3)
      aux=afc(a,npt,tperiodic,sigma)
      do k=1,nr3l
        do j=1,nr2l
          do i=1,nr1l
            !
            ir=i+(j-1)*dfftp%nr1x+(k-1)*dfftp%nr1x*dfftp%nr2x
            gcorr1d(ir)=aux(i+ir1-1,j+ir2-1,k+ir3-1)
            !
          end do
        end do
      end do
      call writetofile(gcorr1d,nnrx,'afc1d.dat',dfftp, 'ax')
      !
      vtemp=gcorr1d
      call fwfft('Dense',vtemp,dfftp)
      do ig=1,ngm
        gcorr1d_fft(ig)=vtemp(np(ig))
      enddo
      !
      deallocate(vtempr)
      deallocate(vtemp)
      !
      return
      !
!-----------------------------------------------------------------------
      end subroutine ee_green_1d_init
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine ee_green_2d_init(box)
!-----------------------------------------------------------------------
!
! ... initialize Green's functions for periodic-image correction
! ... for 2D setttings (e.g., surface, thin film)
!
      use kinds,              only : dp
      use cell_base,          only : a1, a2, a3, omega, tpiba2, s_to_r, &
                                     boxdimensions
      use constants,          only : fpi, pi
      use io_global,          only : stdout
      use grid_dimensions,    only : nnrx, nr1, nr2, nr3, nr1x, nr2x, nr3x, &
                                           nr1l, nr2l, nr3l
      use recvecs_indexes,    only : np, nm
      use reciprocal_vectors, only : gstart, g, gx
      use gvecp,              only : ngm
      use fft_base,           only : dfftp
      use cp_interfaces,      only : fwfft, invfft
      use eecp_mod,           only : gcorr2d,gcorr2d_fft
      use mp_global,          only : me_image
      !
      implicit none
      !
      type(boxdimensions), intent(in) :: box
      !
      real(dp),      parameter :: sigma=2.0_dp
      real(dp),      parameter :: vanishing_dist=1.0e-3_dp
      real(dp),      parameter :: vanishing_g=1.0e-3_dp
      real(dp),      parameter :: euler_gamma=0.57721566490153286061d0
      !
      complex(dp), allocatable :: vtemp(:)
      real(dp),    allocatable :: vtempr(:)
      real(dp) :: aux(dfftp%nr1,dfftp%nr2,dfftp%nr3)
      !
      integer                  :: ig, ir1, ir2, ir3, ir, i, j, k
      real(dp)                 :: sv(3), lv(3), dist
      real(dp)            :: a(3,3)
      integer             :: npt(3)
      logical             :: tperiodic(3)
      !
      real(dp),       external :: qe_erf
      real(dp),       external :: eimlmg
      !
      interface 
        !
        function afc(a,npt,tperiodic,spreadopt)
          !
          real(8), intent(in), optional :: spreadopt
          real(8), intent(in), dimension(3,3) :: a
          integer, intent(in), dimension(3) :: npt
          logical, intent(in), dimension(3) :: tperiodic
          real(8) :: afc(npt(1),npt(2),npt(3))
          !
        end function
        !
      end interface
      !
      ! main body
      !
      allocate(vtemp(nnrx))
      allocate(vtempr(nnrx))
      ! 
      vtemp=0.0_dp
      vtempr=0.0_dp
      gcorr2d=0.0_dp
      !
      ir1=1
      ir2=1
      ir3=1
      do k=1,me_image
        ir3=ir3+dfftp%npp(k)
      enddo
      !   
      npt(1)=dfftp%nr1
      npt(2)=dfftp%nr2
      npt(3)=dfftp%nr3
      tperiodic(1)=.true.
      tperiodic(2)=.true.
      tperiodic(3)=.false.
      a(1:3,1)=a1(1:3)
      a(1:3,2)=a2(1:3)
      a(1:3,3)=a3(1:3)
      aux=afc(a,npt,tperiodic,sigma)
      do k=1,nr3l
        do j=1,nr2l
          do i=1,nr1l
            !
            ir=i+(j-1)*dfftp%nr1x+(k-1)*dfftp%nr1x*dfftp%nr2x
            gcorr2d(ir)=aux(i+ir1-1,j+ir2-1,k+ir3-1)
            !
          end do
        end do
      end do
      call writetofile(gcorr2d,nnrx,'afc2d.dat',dfftp, 'ax')
      !
      vtemp=gcorr2d
      call fwfft('Dense',vtemp,dfftp)
      do ig=1,ngm
        gcorr2d_fft(ig)=vtemp(np(ig))
      enddo
      !
      deallocate(vtempr)
      deallocate(vtemp)
      !
      return
      !
!-----------------------------------------------------------------------
      end subroutine ee_green_2d_init
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine calc_compensation_potential(vcorr_fft,rho_fft,odd_flag)
!-----------------------------------------------------------------------
!
! ... driver for calculating the truncated countercharge (TCC) 
! ... for 0,1,2d periodicity
!
      use kinds,              only: dp
      use gvecp,              only: ngm
      use eecp_mod,           only: gcorr_fft, which_compensation, tcc_odd
      use cell_base,          only: omega

      implicit none
      complex(dp) :: rho_fft(ngm)
      complex(dp) :: vcorr_fft(ngm)
      logical :: odd_flag !true if compensation is being computed
                          !for self-interaction correction

      select case(trim(which_compensation))
          !
        case('tcc')
          !
          call calc_tcc_potential(vcorr_fft,rho_fft)
          !
        case('tcc1d')
          !
          IF((.not.odd_flag).or.(.not.tcc_odd)) THEN
            call calc_tcc1d_potential(vcorr_fft,rho_fft)
          ELSE IF(odd_flag.and.tcc_odd) THEN
            call calc_tcc_potential(vcorr_fft,rho_fft)
          ENDIF
          !
        case('tcc2d')
          !
          IF((.not.odd_flag).or.(.not.tcc_odd)) THEN
            call calc_tcc2d_potential(vcorr_fft,rho_fft)
          ELSE IF(odd_flag.and.tcc_odd) THEN
            call calc_tcc_potential(vcorr_fft,rho_fft)
          ENDIF
          !
        case('none')
          !
          continue
          !
        case default
          !
          call errore('vofrho','Invalid correction: '//TRIM(which_compensation), 10)
          !
        end select
!-----------------------------------------------------------------------
      end subroutine calc_compensation_potential
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      subroutine calc_tcc_potential(vcorr_fft,rho_fft)
!-----------------------------------------------------------------------
!
! ... calculate the truncated countercharge (TCC) 
! ... periodic-image correction potential in
! ... reciprocal space for 0D settings
!
      use kinds,              only: dp
      use gvecp,              only: ngm
      use eecp_mod,           only: gcorr_fft
      use cell_base,          only: omega
      !
      implicit none
      complex(dp) :: rho_fft(ngm)
      complex(dp) :: vcorr_fft(ngm)
      integer :: ig
      !
      do ig=1,ngm
        vcorr_fft(ig)=omega*gcorr_fft(ig)*rho_fft(ig)
      end do
      !
      return
!
!-----------------------------------------------------------------------
      end subroutine calc_tcc_potential
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine calc_tcc1d_potential(vcorr_fft,rho_fft)
!-----------------------------------------------------------------------
!
! ... calculate the truncated countercharge (TCC) 
! ... periodic-image correction potential in
! ... reciprocal space for 1D settings
!
      use kinds,              only: dp
      use gvecp,              only: ngm
      use eecp_mod,           only: gcorr1d_fft
      use cell_base,          only: omega
      !
      implicit none
      complex(dp) :: rho_fft(ngm)
      complex(dp) :: vcorr_fft(ngm)
      integer :: ig
      !
      do ig=1,ngm
        vcorr_fft(ig)=omega*gcorr1d_fft(ig)*rho_fft(ig)
      end do
      !
      return
!
!-----------------------------------------------------------------------
      end subroutine calc_tcc1d_potential
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine calc_tcc2d_potential(vcorr_fft,rho_fft)
!-----------------------------------------------------------------------
!
! ... calculate the truncated countercharge (TCC) 
! ... periodic-image correction potential in
! ... reciprocal space for 2D settings
!
      use kinds,              only: dp
      use gvecp,              only: ngm
      use eecp_mod,           only: gcorr2d_fft
      use cell_base,          only: omega
      !
      implicit none
      complex(dp) :: rho_fft(ngm)
      complex(dp) :: vcorr_fft(ngm)
      integer :: ig
      !
      do ig=1,ngm
        vcorr_fft(ig)=omega*gcorr2d_fft(ig)*rho_fft(ig)
      end do
      !
      return
!
!-----------------------------------------------------------------------
      end subroutine calc_tcc2d_potential
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine calc_tcc_energy(ecomp,vcorr_fft,rho_fft)
!-----------------------------------------------------------------------
!
! ... calculate the truncated countercharge (TCC) 
! ... periodic-image corrective energy in for 0D settings
!
      use kinds,              only : dp
      use gvecp,              only : ngm
      use eecp_mod,           only : gcorr_fft
      use cell_base,          only : omega
      use reciprocal_vectors, only : gstart
      use mp,                 only : mp_sum
      use mp_global,          only : intra_image_comm
      use control_flags,      only : gamma_only, do_wf_cmplx!added:giovanni
      !
      implicit none
      !
      real(dp),    intent(out) :: ecomp
      complex(dp), intent(in)  :: rho_fft(ngm)
      complex(dp), intent(in)  :: vcorr_fft(ngm)
      !
      complex(dp), allocatable :: aux(:)
      integer      :: ig
      complex(dp)  :: zh
      real(dp) :: wz
      logical      :: lgam
      !
      lgam = gamma_only.and..not.do_wf_cmplx
      !
      IF(lgam) THEN
         wz=2.d0
      ELSE
         wz=1.d0
      ENDIF
      !
      allocate(aux(ngm))
      !
      aux=0.0_dp
      !
      if(gstart.ne.1) then
        aux(1)=0.5d0*omega*vcorr_fft(1)*conjg(rho_fft(1))
      end if
      !
      do ig=gstart,ngm
        aux(ig)=0.5d0*wz*omega*vcorr_fft(ig)*conjg(rho_fft(ig))
      end do
      !
      zh=0.0_dp
      do ig=1,ngm
        zh=zh+aux(ig)
      enddo
      ecomp=dble(zh)
      !
      call mp_sum(ecomp,intra_image_comm)
      !
      deallocate(aux)
      !
      return
      !
!-----------------------------------------------------------------------
      end subroutine calc_tcc_energy
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine ee_efieldpot_init(box)
!-----------------------------------------------------------------------
!
! ... add electric field with a sawtooth potential in real space
! ... (the system is assumed to be centered at the origin of the cell)
!
      use kinds,              only : dp
      use cell_base,          only : a1, a2, a3, omega, tpiba2, s_to_r, &
                                     boxdimensions
      use constants,          only : fpi, pi
      use io_global,          only : stdout
      use grid_dimensions,    only : nr1,  nr2,  nr3,  nr1x, nr2x, nr3x, &
                                     nr1l, nr2l, nr3l, nnrx
      use recvecs_indexes,    only : np, nm
      use reciprocal_vectors, only : gstart, g
      use gvecp,              only : ngm
      use fft_base,           only : dfftp
      use cp_interfaces,      only : fwfft, invfft
      use efield_mod,         only : efieldpot, efieldpotg, ampfield
      use mp_global,          only : me_image
      !
      implicit none
      !
      type(boxdimensions), intent(in) :: box
      !
      integer    :: ig, ir1, ir2, ir3, ir, i, j, k
      real(dp)   :: sv(3), lv(3)
      complex(dp), allocatable :: vtemp(:)
      !
      allocate(vtemp(nnrx))
      ! 
      efieldpot=0.0_dp
      efieldpotg=0.0_dp
      !
      ir1=1
      ir2=1
      ir3=1
      do k=1,me_image
        ir3=ir3+dfftp%npp(k)
      enddo      
      !
      do k=1,nr3l
        do j=1,nr2l
          do i=1,nr1l
            !
            sv(1)=dble((i-1)+(ir1-1))/nr1 
            sv(2)=dble((j-1)+(ir2-1))/nr2
            sv(3)=dble((k-1)+(ir3-1))/nr3
            !
            if(sv(1)>0.5_dp) sv(1)=sv(1)-1.0_dp
            if(sv(2)>0.5_dp) sv(2)=sv(2)-1.0_dp
            if(sv(3)>0.5_dp) sv(3)=sv(3)-1.0_dp
            !
            call s_to_r(sv,lv,box%hmat)
            !
            ir=i+(j-1)*dfftp%nr1x+(k-1)*dfftp%nr1x*dfftp%nr2x
            !
            efieldpot(ir)=dot_product(lv,ampfield)
            !
          end do
        end do
      end do
      !
      vtemp=efieldpot
      call fwfft('Dense',vtemp,dfftp)
      !
      do ig=1,ngm
        efieldpotg(ig)=vtemp(np(ig))
      end do
      !
      call writetofile(efieldpot,nnrx,'efieldpot.dat',dfftp, 'az')
      !
      deallocate( vtemp )
      !
      return
      !
!-----------------------------------------------------------------------
      end subroutine ee_efieldpot_init
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine calc_dipole(wfc,box)
!-----------------------------------------------------------------------
!
! ... calculate the electronic dipole in real space
!
      use kinds,              only : dp
      use cell_base,          only : a1, a2, a3, omega, tpiba2, s_to_r, &
                                     boxdimensions
      use constants,          only : fpi, pi
      use io_global,          only : stdout, meta_ionode
      use gvecw,              only : ngw
      use grid_dimensions,    only : nr1l, nr2l, nr3l, nnrx, &
                                     nr1,  nr2, nr3, nr1x, nr2x, nr3x
      use smooth_grid_dimensions, &
                              only : nr1s, nr2s, nr3s, &
                                     nr1sx, nr2sx, nr3sx, nnrsx
      use recvecs_indexes,    only : np, nm
      use reciprocal_vectors, only : gstart, g
      use gvecs,              only : ngs, nps, nms
      use gvecp,              only : ngm
      use fft_base,           only : dfftp, dffts
      use cp_interfaces,      only : fwfft, invfft
      use electrons_base,     only : nbspx, nbsp, &
                                     f, ispin, nspin, nelt
      use mp_global,          only : me_image, intra_image_comm
      use mp,                 only : mp_sum
      !
      implicit none
      !
      type(boxdimensions),  intent(in) :: box
      complex(dp),       intent(inout) :: wfc(ngw,nbspx)
      !
      integer       :: ig, ir1, ir2, ir3, ir, i, j, k
      integer       :: iss, isup, isdw, iss1, iss2, ios
      real(dp)      :: sa1, sa2, dipole(3), sv(3), lv(3) 
      complex(dp)   :: ci, fp, fm
      !
      complex(dp), allocatable :: vtemp(:)
      complex(dp), allocatable :: psi(:), psis(:), rhog(:,:)
      real(dp),    allocatable :: rhor(:,:), rhos(:,:)
      !
      ci=(0.0d0,1.0d0)
      !
      allocate(rhos(nnrsx,2))
      allocate(rhor(nnrx,2))
      allocate(rhog(ngm,2))
      !
      rhor=0.0_dp
      rhos=0.0_dp
      rhog=(0.0_dp,0.0_dp)
      !
      ! ... calculate total charge density 
      ! ... (this step can be removed depending on where the 
      ! ... subroutine is called, ID)
      !
      allocate(psis(nnrsx)) 
      !
      do i=1,nbsp,2
        !
        call c2psi(psis,nnrsx,wfc(1,i),wfc(1,i+1),ngw,2)
        call invfft('wave',psis,dffts)
        !
        iss1=ispin(i)
        sa1=f(i)/omega
        !
        if(i/=nbsp) then
          iss2=ispin(i+1)
          sa2=f(i+1)/omega
        else
          iss2=iss1
          sa2=0.0_dp
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
         fp= psis(nps(ig))+psis(nms(ig))
         fm= psis(nps(ig))-psis(nms(ig))
         rhog(ig,isup)=0.5d0*cmplx(dble(fp),aimag(fm))
         rhog(ig,isdw)=0.5d0*cmplx(aimag(fp),-dble(fm))
      end do
      !
      allocate(psi(nnrx))
      !
      isup=1
      isdw=2
      psi(:)=(0.0_dp,0.0_dp)
      !
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
      call writetofile(rhor,nnrx,'rhodipolez.dat',dfftp, 'az')
      call writetofile(rhor,nnrx,'rhodipolex.dat',dfftp, 'ax')
      !
      deallocate(psi) 
      deallocate(psis) 
      ! 
      dipole(1:3)=0.0_dp
      !
      ir1=1
      ir2=1
      ir3=1
      do k=1,me_image
        ir3=ir3+dfftp%npp(k)
      end do      
      !
      do k=1,nr3l
        do j=1,nr2l
          do i=1,nr1l
            !
            sv(1)=dble((i-1)+(ir1-1))/nr1 
            sv(2)=dble((j-1)+(ir2-1))/nr2
            sv(3)=dble((k-1)+(ir3-1))/nr3
            !
            if(sv(1)>0.5_dp) sv(1)=sv(1)-1.0_dp
            if(sv(2)>0.5_dp) sv(2)=sv(2)-1.0_dp
            if(sv(3)>0.5_dp) sv(3)=sv(3)-1.0_dp
            !
            call s_to_r(sv,lv,box%hmat)
            !
            dipole(:)=dipole(:)+lv(:)*sum(rhor(ir,1:nspin))
            !
          end do
        end do
      end do
      !
      call mp_sum(dipole,intra_image_comm)
      !
      if(meta_ionode) write(stdout,2015) dipole
2015  format( 3x,'electronic dipole moment = ',3f12.6)
      !
      return
      !
!--------------------------------------------------------------------
      end subroutine calc_dipole
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      function eimlmg( xx )
!--------------------------------------------------------------------
      !
      ! ... exponential integral minus logarithm minus gamma function
      ! ...
      ! ... calculate eimlmg( x ) = ei( x ) - ln | x | - euler_gamma
      ! ...                       = x / 1 / 1! + x ^ 2 / 2 / 2! + ...
      ! ...                       + x ^ n / n / n! + ...
      !
      use kinds,         only : dp
      !
      implicit none
      !
      real( dp )             :: eimlmg
      real( dp )             :: xx
      !
      real( dp )             :: fact
      real( dp )             :: term
      real( dp )             :: total1
      real( dp )             :: total2
      !
      integer                :: k
      !
      integer, parameter     :: maxit=100000
      !
      real( dp ), parameter  :: eps=1.d-20
      real( dp ), parameter  :: euler_gamma=0.57721566490153286061d0
      real( dp ), parameter  :: xxlim=-20.d0
      !
      if ( xx > xxlim ) then
        eimlmg = 0.d0
        fact = 1.d0
        summation : do k = 1, maxit
          fact = fact * xx / dble( k )
          term = fact / dble( k )
          eimlmg = eimlmg + term
          if( abs( term ) .lt. eps ) exit summation
        end do summation
      else
        eimlmg = - log( abs( xx ) ) - euler_gamma
      end if
      !
      return
      !
!--------------------------------------------------------------------
      end function eimlmg
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      function pinterp( x, side, bound )
!--------------------------------------------------------------------
      !
      ! ... calculate the interpolation polynomial satifying
      ! ... P'(0)=0,P'(1)=0,P(0)=1,P(1)=0 for side=0,bound=0
      ! ... P'(0)=0,P'(1)=0,P(0)=0,P(1)=1 for side=1,bound=0
      ! ... P'(1)=0,P(0)=0,P(1)=1 for side=1,bound=-1
      ! ... P'(1)=0,P(0)=1,P(1)=0 for side=0,bound=-1
      ! ... P'(0)=0,P(0)=0,P(1)=1 for side=1,bound=1
      ! ... P'(0)=0,P(0)=1,P(1)=0 for side=0,bound=1
      !
      use kinds,         only : dp
      !
      implicit none
      !
      real( dp ) :: pinterp
      real( dp ) :: x
      integer    :: side
      integer    :: bound
      !
      if( bound == 0 .and. side == 1 ) then
        pinterp = 3.d0 * x * x - 2.d0 * x * x * x
      else if( bound == 0 .and. side == 0 ) then
        pinterp = 1.d0 - 3.d0 * x * x + 2.d0 * x * x * x
      else if( bound == - 1 .and. side == 0 ) then
        pinterp = 1.d0 - 2.d0 * x + x * x
      else if( bound == - 1 .and. side ==  1 ) then
        pinterp = 2.d0 * x - x * x
      else if( bound == 1 .and. side == 1 ) then
        pinterp = x * x
      else if( bound == 1 .and. side == 0 ) then
        pinterp =  1 - x * x
      end if
      !
      return
      !
!--------------------------------------------------------------------
      end function pinterp
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      function dpinterp( x, side, bound )
!--------------------------------------------------------------------
      !
      ! ... calculate the derivative of P(x)
      !
      use kinds,         only : dp
      !
      implicit none
      !
      real( dp ) :: dpinterp
      real( dp ) :: x
      integer    :: side
      integer    :: bound
      !
      if( bound == 0 .and. side == 1 ) then
        dpinterp = 6.d0 * x  - 6.d0 * x * x
      else if( bound == 0 .and. side == 0 ) then
        dpinterp = - 6.d0 * x + 6.d0 * x * x
      else if( bound == - 1 .and. side == 0 ) then
        dpinterp = - 2.d0  + 2.d0 * x
      else if( bound == - 1 .and. side ==  1 ) then
        dpinterp = 2.d0 - 2.d0 * x
      else if( bound == 1 .and. side == 1 ) then
        dpinterp = 2.d0 * x
      else if( bound == 1 .and. side == 0 ) then
        dpinterp =  - 2.d0 * x
      end if
      !
      return
      !
!--------------------------------------------------------------------
      end function dpinterp
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      function qinterp( x, side, bound )
!--------------------------------------------------------------------
      !
      ! ... calculate the interpolation polynomial satifying
      ! ... Q'(0)=1,Q'(1)=0,Q(0)=0,Q(1)=0 for side=0,bound=0
      ! ... Q'(0)=0,Q'(1)=1,Q(0)=0,Q(1)=0 for side=1,bound=0
      ! ... Q'(1)=1,Q(0)=0,Q(1)=0 for side=1,bound=-1
      ! ... Q'(0)=1,Q(0)=0,Q(1)=0 for side=0,bound=-1
      ! ... Q'(1)=1,Q(0)=0,Q(1)=0 for side=1,bound=1
      ! ... Q'(0)=1,Q(0)=0,Q(1)=0 for side=0,bound=1
      !
      use kinds,         only : dp
      !
      implicit none
      !
      real( dp ) :: qinterp
      real( dp ) :: x
      integer    :: side
      integer    :: bound
      !
      if( bound == 0 .and. side == 1 ) then
        qinterp = - x * x + x * x * x
      else if( bound == 0 .and. side == 0 ) then
        qinterp = x - 2.d0 * x * x + x * x * x
      else if( bound == - 1 .and. side == 0 ) then
        qinterp = 0.d0
      else if( bound == 1 .and. side == 1 ) then
        qinterp =  0.d0
      else if( bound == - 1 .and. side == 1 ) then
        qinterp = - x + x * x
      else if( bound == 1 .and. side == 0 ) then
        qinterp = x - x * x
      end if
      !
      return
      !
!--------------------------------------------------------------------
      end function qinterp
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      function dqinterp( x, side, bound )
!--------------------------------------------------------------------
      !
      ! ... calculate the derivative of Q(x)
      !
      use kinds,         only : dp
      !
      implicit none
      !
      real( dp ) :: dqinterp
      real( dp ) :: x
      integer    :: side
      integer    :: bound
      !
      if( bound == 0 .and. side == 1 ) then
        dqinterp = - 2.d0 * x  + 3.d0 * x * x
      else if( bound == 0 .and. side == 0 ) then
        dqinterp = 1 - 4.d0 * x  + 3 * x * x
      else if( bound == - 1 .and. side == 0 ) then
        dqinterp = 0.d0
      else if( bound == 1 .and. side == 1 ) then
        dqinterp =  0.d0
      else if( bound == - 1 .and. side == 1 ) then
        dqinterp = - 1.d0 + 2.d0 * x
      else if( bound == 1 .and. side == 0 ) then
        dqinterp = 1.d0 - 2.d0 * x
      end if
      !
      return
      !
!--------------------------------------------------------------------
      end function dqinterp
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      function compmod( ir, nr )
!--------------------------------------------------------------------
      !
      ! ... calculate the composite grid index corresponding to ir
      !
      implicit none
      !
      integer :: compmod
      integer, intent(in) :: ir
      integer, intent(in) :: nr
      !
      compmod = modulo( ir - 1, nr ) + 1
      !
      return
      !
!--------------------------------------------------------------------
      end function compmod
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      function bound( i, n )
!--------------------------------------------------------------------
      !
      ! ... return -1 if i = 0, 1 if i = n - 1, and 0 otherwise
      !
      implicit none
      !
      integer :: i
      integer :: n
      integer :: bound
      !
      if( i == 1 ) then
        bound = -1
      else if( i == n - 1 ) then
        bound = 1
      else
        bound = 0
      end if
      !
      return
      !
!--------------------------------------------------------------------
      end function bound
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      subroutine calc_fcorr(fion,vcorr,taus,nat,na,nsp,box,dfft)
!--------------------------------------------------------------------
      !
      ! ... calculate the interatomic force contribution due to
      ! ... the corrective potential
      !
      use kinds,              only : dp
      use ions_base,          only : zv
      use cell_base,          only : a1, a2, a3, omega, tpiba2, s_to_r, &
                                     boxdimensions
      use io_global,          only : meta_ionode, meta_ionode_id
      use mp,                 only : mp_sum, mp_barrier, mp_gather
      use mp_global,          only : intra_image_comm, me_image, &
                                     nproc_image
      use grid_dimensions,    only : nr1, nr2, nr3, nr1x, nr2x, nr3x, &
                                     nnrx
      use fft_types,          only : fft_dlay_descriptor
      !
      implicit none
      !
      real(dp)               :: fion(3,nat)
      real(dp)               :: taus(3,nat)
      real(dp)               :: vcorr(nnrx)
      integer                :: nat,na(nsp),nsp
      type(boxdimensions)    :: box
      type(fft_dlay_descriptor) :: dfft
      !
      integer                :: ir,ir1,ir2,ir3,ia,a,b,c,&
                                bound1,bound2,bound3,na_loc,ia_e,ia_s, &
                                proc,is,isa
      real(dp)               :: taur(3,nat)
      real(dp)               :: delta1,delta2,delta3, &
                                t1,t2,t3,df1,df2,df3,f,g1,g2,g3
      !
      integer,  allocatable :: displs(:),recvcount(:)
      real(dp), allocatable    :: vcomp(:)
      real(dp), allocatable    :: fcorr(:,:)
      complex(dp), allocatable :: aux(:)
      !
      integer, external      :: compindex
      integer, external      :: compmod
      integer, external      :: bound
      integer, external      :: ldim_block
      integer, external      :: gind_block
      real(dp), external     :: pinterp
      real(dp), external     :: qinterp
      real(dp), external     :: dpinterp
      real(dp), external     :: dqinterp
      !
      ! ... initializes variables
      !
      delta1=a1(1)/dble(nr1)
      delta2=a2(2)/dble(nr2)
      delta3=a3(3)/dble(nr3)
      !
      df1=0.0_dp
      df2=0.0_dp
      df3=0.0_dp
      !
      bound1=0
      bound2=0
      bound3=0
      !
      allocate(fcorr(3,nat))
      fcorr(:,:)=0.0_dp
      !
      ! ... collect vcomp and scatter across nodes (cf. old_write_rho)
      !
      nr1=dfft%nr1
      nr2=dfft%nr2
      nr3=dfft%nr3
      !
      nr1x=dfft%nr1x
      nr2x=dfft%nr2x
      nr3x=dfft%nr3x
      !
      allocate(displs(nproc_image),recvcount(nproc_image))
      allocate(vcomp(nr1x*nr2x*nr3x))
      !
      vcomp=0.0_dp
      !
      if(nproc_image>1) then
        !
        do proc=1,nproc_image
          !
          recvcount(proc)=dfft%nnp*dfft%npp(proc)
          !
          if(proc==1) then
            displs(proc)=0
          else
            displs(proc)=displs(proc-1)+recvcount(proc-1)
          endif
          !
        enddo
        !
        call mp_barrier()
        call mp_gather(vcorr,vcomp,recvcount,displs, &
                       meta_ionode_id,intra_image_comm)
        call mp_sum(vcomp,intra_image_comm)
        !
      else
        !  
        if (nr1/=nr1x.or.nr2/=nr2x.or.nr3/=nr3x) &
          call errore('calc_fcorr','dimension mistmatch',10)
        !
        vcomp(1:nr1x*nr2x*nr3x)=vcorr(1:nnrx)
        !
      endif
      !
      ! ... distributes atoms over processors (cf. vofesr)
      !
      na_loc=ldim_block(nat,nproc_image,me_image)
      ia_s=gind_block(1,nat,nproc_image,me_image )
      ia_e=ia_s+na_loc-1
      !
      do ia=ia_s,ia_e
        !
        call s_to_r(taus(1:3,ia),taur(1:3,ia),box%hmat)
        !
        t1=taur(1,ia)/delta1
        t2=taur(2,ia)/delta2
        t3=taur(3,ia)/delta3
        !
        ir1=int(t1)+1
        ir2=int(t2)+1
        ir3=int(t3)+1
        !
        t1=t1-dble(ir1-1)
        t2=t2-dble(ir2-1)
        t3=t3-dble(ir3-1)
        !
        ir1=compmod(ir1,nr1)
        ir2=compmod(ir2,nr2)
        ir3=compmod(ir3,nr3)
        !
        ! ... with TCC, we use a periodic interpolation
        !
        !bound1=bound(ir1,nr1)
        !bound2=bound(ir2,nr2)
        !bound3=bound(ir3,nr3)
        !
        f=0 
        g1=0
        g2=0
        g3=0
        !
        do a=0,1
         do b=0,1
          do c=0,1
           !
           f=vcomp(compindex(ir1+a,ir2+b,ir3+c,nr1x,nr2x,nr3x))
           !
           g1=f*dpinterp(t1,a,bound1)*pinterp(t2,b,bound2) &
             *pinterp(t3,c,bound3)
           g2=f*pinterp(t1,a,bound1)*dpinterp(t2,b,bound2) &
             *pinterp(t3,c,bound3)
           g3=f*pinterp(t1,a,bound1)*pinterp(t2,b,bound2) &
             *dpinterp(t3,c,bound3)
           !
           df1=0.5_dp*(vcomp(compindex(ir1+a+1,ir2+b,ir3+c,nr1x,nr2x,nr3x)) &
                      -vcomp(compindex(ir1+a-1,ir2+b,ir3+c,nr1x,nr2x,nr3x)))
           df2=0.5_dp*(vcomp(compindex(ir1+a,ir2+b+1,ir3+c,nr1x,nr2x,nr3x)) &
                      -vcomp(compindex(ir1+a,ir2+b-1,ir3+c,nr1x,nr2x,nr3x)))
           df3=0.5_dp*(vcomp(compindex(ir1+a,ir2+b,ir3+c+1,nr1x,nr2x,nr3x)) &
                      -vcomp(compindex(ir1+a,ir2+b,ir3+c-1,nr1x,nr2x,nr3x)))
           !     
           fcorr(1,ia)=fcorr(1,ia)+g1 &
                      +df1*dqinterp(t1,a,bound1)*pinterp(t2,b,bound2) &
                      *pinterp(t3,c,bound3)                       
           fcorr(2,ia)=fcorr(2,ia)+g2 &
                      +df2*pinterp(t1,a,bound1)*dqinterp(t2,b,bound2) &
                      *pinterp(t3,c,bound3)                      
           fcorr(3,ia)=fcorr(3,ia)+g3 &
                      +df3*pinterp(t1,a,bound1)*pinterp(t2,b,bound2) &
                      *dqinterp(t3,c,bound3)
           !
          end do
         end do
        end do
        !
        fcorr(1,ia)=fcorr(1,ia)/delta1
        fcorr(2,ia)=fcorr(2,ia)/delta2
        fcorr(3,ia)=fcorr(3,ia)/delta3
        !
      end do
      !
      call mp_sum(fcorr,intra_image_comm)
      !
      isa=0
      do is=1,nsp
        do ia=1,na(is)
          isa=isa+1
          fion(1,isa)=fion(1,isa)+fcorr(1,isa)*zv(is)
          fion(2,isa)=fion(2,isa)+fcorr(2,isa)*zv(is)
          fion(3,isa)=fion(3,isa)+fcorr(3,isa)*zv(is)
        end do
      end do
      !
      deallocate(displs,recvcount)
      deallocate(vcomp)
      deallocate(fcorr)
      !
!--------------------------------------------------------------------
      end subroutine calc_fcorr
!--------------------------------------------------------------------


