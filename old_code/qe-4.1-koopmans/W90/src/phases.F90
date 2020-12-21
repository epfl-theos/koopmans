!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!

!! This module calculates the phases of the WFs
subroutine compute_phases

  use w90_constants,  only : dp,cmplx_0,cmplx_i,twopi,cmplx_1
  use w90_io,         only : io_error,stdout,io_file_unit,seedname, &
                             io_date,io_stopwatch
  use w90_parameters, only : num_wann,num_bands,num_kpts,u_matrix,spin, &
       ngs=>wannier_plot_supercell,kpt_latt,real_lattice,have_disentangled, &
       ndimwin,lwindow,u_matrix_opt,wvfn_formatted,timing_level,wannier_plot_format

  implicit none

  real(kind=dp) :: scalfac,tmax,tmaxx
  real(kind=dp) :: w_real,w_imag,ratmax,ratio
  complex(kind=dp), allocatable :: wann_func(:,:,:,:)
  complex(kind=dp), allocatable :: r_wvfn(:,:)
  complex(kind=dp), allocatable :: r_wvfn_tmp(:,:)
  complex(kind=dp) :: catmp,wmod

  logical :: have_file
  integer :: i,j,nbnd,counter,ierr
  integer :: loop_kpt,ik,ix,iy,iz,nk,ngx,ngy,ngz,nxx,nyy,nzz
  integer :: loop_b,nx,ny,nz,npoint,file_unit,loop_w,num_inc
  character(len=11) :: wfnname
  character(len=9)  :: cdate, ctime
  logical           :: inc_band(num_bands)  
  !
  if (timing_level>1) call io_stopwatch('plot: wannier',1)
  !
  write(wfnname,200 ) 1,spin
  inquire(file=wfnname,exist=have_file)
  if(.not.have_file) call io_error('plot_wannier: file '//wfnname//' not found') 

  file_unit=io_file_unit()
  if(wvfn_formatted) then
     open(unit=file_unit,file=wfnname,form='formatted')
     read(file_unit,*) ngx,ngy,ngz,nk,nbnd
  else
     open(unit=file_unit,file=wfnname,form='unformatted')
     read(file_unit) ngx,ngy,ngz,nk,nbnd
  end if
  close(file_unit)

200 format ('UNK',i5.5,'.',i1)
  allocate(wann_func(-((ngs(1))/2)*ngx:((ngs(1)+1)/2)*ngx-1,&
       -((ngs(2))/2)*ngy:((ngs(2)+1)/2)*ngy-1,&
       -((ngs(3))/2)*ngz:((ngs(3)+1)/2)*ngz-1,num_wann),stat=ierr ) 
  if (ierr/=0) call io_error('Error in allocating wann_func in plot_wannier')
  if(have_disentangled) then
     allocate(r_wvfn_tmp(ngx*ngy*ngz,maxval(ndimwin)),stat=ierr )
     if (ierr/=0) call io_error('Error in allocating r_wvfn_tmp in plot_wannier')
  end if
  allocate(r_wvfn(ngx*ngy*ngz,num_wann),stat=ierr )
!  write(stdout,*) 'RICCARDO: ', -((ngs(1))/2)*ngx,':',((ngs(1)+1)/2)*ngx-1,-((ngs(2))/2)*ngy,':',((ngs(2)+1)/2)*ngy-1,-((ngs(3))/2)*ngz,':',((ngs(3)+1)/2)*ngz-1,num_wann
!  write(stdout,*) 'RICCARDO: ',ngx*ngy*ngz,num_wann
  if (ierr/=0) call io_error('Error in allocating r_wvfn in plot_wannier')
  wann_func=cmplx_0


  call io_date(cdate, ctime)
  do loop_kpt=1,num_kpts

     inc_band=.true.
     num_inc=num_wann
     if(have_disentangled) then
        inc_band(:)=lwindow(:,loop_kpt)
        num_inc=ndimwin(loop_kpt)
     end if

     write(wfnname,200 ) loop_kpt,spin
     file_unit=io_file_unit()
     if(wvfn_formatted) then
        open(unit=file_unit,file=wfnname,form='formatted')
        read(file_unit,*) ix,iy,iz,ik,nbnd
     else
        open(unit=file_unit,file=wfnname,form='unformatted')
        read(file_unit) ix,iy,iz,ik,nbnd
     end if

     if ( (ix/=ngx) .or. (iy/=ngy) .or. (iz/=ngz) .or. (ik/=loop_kpt) ) then
        write(stdout,'(1x,a,a)') 'WARNING: mismatch in file',trim(wfnname)
        write(stdout,'(1x,5(a6,I5))')  '   ix=',ix ,'   iy=',iy ,'   iz=',iz ,'   ik=',ik      ,' nbnd=',nbnd
        write(stdout,'(1x,5(a6,I5))')  '  ngx=',ngx,'  ngy=',ngy,'  ngz=',ngz,'  kpt=',loop_kpt,'bands=',num_bands
        call io_error('plot_wannier')
     end if

     if(have_disentangled) then
        counter=1
        do loop_b=1,num_bands
           if(counter>num_inc) exit
           if(wvfn_formatted) then
              do nx=1,ngx*ngy*ngz
                 read(file_unit,*) w_real, w_imag
                 r_wvfn_tmp(nx,counter) = cmplx(w_real,w_imag,kind=dp)
              end do
           else
              read(file_unit) (r_wvfn_tmp(nx,counter),nx=1,ngx*ngy*ngz)
           end if
           if(inc_band(loop_b)) counter=counter+1
        end do
     else
        do loop_b=1,num_bands 
           if(wvfn_formatted) then
              do nx=1,ngx*ngy*ngz
                 read(file_unit,*) w_real, w_imag
                 r_wvfn(nx,loop_b) = cmplx(w_real,w_imag,kind=dp)
              end do
           else
              read(file_unit) (r_wvfn(nx,loop_b),nx=1,ngx*ngy*ngz)
           end if
        end do
     end if

     close(file_unit)

     if(have_disentangled) then
        r_wvfn=cmplx_0
        do loop_w=1,num_wann
           do loop_b=1,num_inc
              r_wvfn(:,loop_w)=r_wvfn(:,loop_w)+ &
                   u_matrix_opt(loop_b,loop_w,loop_kpt)*r_wvfn_tmp(:,loop_b)
           end do
        end do
     end if


     ! nxx, nyy, nzz span a parallelogram in the real space mesh, of side
     ! 2*nphir, and centered around the maximum of phi_i, nphimx(i, 1 2 3)
     !
     ! nx ny nz are the nxx nyy nzz brought back to the unit cell in
     ! which u_nk(r)=cptwrb(r,n)  is represented
     !
     ! There is a big performance improvement in looping over num_wann
     ! in the inner loop. This is poor memory access for wann_func and
     ! but the reduced number of operations wins out. 

     do nzz =-((ngs(3))/2)*ngz,((ngs(3)+1)/2)*ngz-1
        nz=mod(nzz,ngz)
        if(nz.lt.1) nz=nz+ngz
        do nyy=-((ngs(2))/2)*ngy,((ngs(2)+1)/2)*ngy-1
           ny=mod(nyy,ngy)
           if(ny.lt.1) ny=ny+ngy
           do nxx=-((ngs(1))/2)*ngx,((ngs(1)+1)/2)*ngx-1
              nx=mod(nxx,ngx)
              if(nx.lt.1) nx=nx+ngx

              scalfac=kpt_latt(1,loop_kpt)*real(nxx-1,dp)/real(ngx,dp)+ &
                   kpt_latt(2,loop_kpt)*real(nyy-1,dp)/real(ngy,dp)+ &
                   kpt_latt(3,loop_kpt)*real(nzz-1,dp)/real(ngz,dp)
              npoint=nx+(ny-1)*ngx+(nz-1)*ngy*ngx
              catmp=exp(twopi*cmplx_i*scalfac)
              do loop_b=1,num_wann
                 do loop_w=1,num_wann            
                    wann_func(nxx,nyy,nzz,loop_w)=wann_func(nxx,nyy,nzz,loop_w)+ &
                         u_matrix(loop_b,loop_w,loop_kpt)*r_wvfn(npoint,loop_b)*catmp
                 end do
              end do
           end do
        end do

     end do

  end do !loop over kpoints

  ! fix the global phase by setting the wannier to
  ! be real at the point where it has max. modulus


  ! Write the phases in 'wf_phases.dat'
  file_unit=io_file_unit()
  open(unit=file_unit,file='wf_phases.dat',form='formatted',status='unknown')

  do loop_w=1,num_wann
     tmaxx=0.0
     wmod=cmplx_1
     do nzz=-((ngs(3))/2)*ngz,((ngs(3)+1)/2)*ngz-1
        do nyy=-((ngs(2))/2)*ngy,((ngs(2)+1)/2)*ngy-1
           do nxx=-((ngs(1))/2)*ngx,((ngs(1)+1)/2)*ngx-1
              wann_func(nxx,nyy,nzz,loop_w)= wann_func(nxx,nyy,nzz,loop_w)/ real(num_kpts,dp)
              tmax=real(wann_func(nxx,nyy,nzz,loop_w)* & 
                   conjg(wann_func(nxx,nyy,nzz,loop_w)),dp)
              if (tmax>tmaxx) then
                 tmaxx=tmax
                 wmod=wann_func(nxx,nyy,nzz,loop_w)
              end if
           end do
        end do
     end do
     wmod=wmod/sqrt(real(wmod)**2+aimag(wmod)**2)
     wann_func(:,:,:,loop_w)=wann_func(:,:,:,loop_w)/wmod
     !
     write(file_unit,'(2E13.5)') real(wmod),aimag(wmod)
     !
  end do

  close(file_unit)
  !
  ! Check the 'reality' of the WF
  !
  do loop_w=1,num_wann
     ratmax=0.0_dp
     do nzz=-((ngs(3))/2)*ngz,((ngs(3)+1)/2)*ngz-1
        do nyy=-((ngs(2))/2)*ngy,((ngs(2)+1)/2)*ngy-1
           do nxx=-((ngs(1))/2)*ngx,((ngs(1)+1)/2)*ngx-1
              if (abs(real(wann_func(nxx,nyy,nzz,loop_w),dp))>=0.01_dp) then
                 ratio=abs(aimag(wann_func(nxx,nyy,nzz,loop_w)))/ &
                      abs(real(wann_func(nxx,nyy,nzz,loop_w),dp))
                 ratmax=max(ratmax,ratio)
              end if
           end do
        end do
     end do
     write(stdout,'(6x,a,i4,7x,a,f11.6)') 'Wannier Function Num: ',loop_w,&
          'Maximum Im/Re Ratio = ',ratmax
  end do
  write(stdout,*) ' '
  
  return

end subroutine compute_phases
