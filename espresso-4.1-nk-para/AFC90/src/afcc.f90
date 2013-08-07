!
! LibAFCC - Library for auxiliary-function countercharge correction 
! Copyright (c) 2010-2011 I. Dabo
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. See GPL/gpl-3.0.txt. 
! If not, see <http://www.gnu.org/licenses/>.
!
function afc(a,npt,tperiodic,spreadopt)
  !
  implicit none
  !
  real(8), intent(in), optional :: spreadopt
  real(8), intent(in), dimension(3,3) :: a
  integer, intent(in), dimension(3) :: npt
  logical, intent(in), dimension(3) :: tperiodic
  real(8), dimension(npt(1),npt(2),npt(3)) :: afc
  !
  real(8), dimension(3,3) :: e
  real(8), dimension(3,3) :: aaux
  real(8), allocatable, dimension(:) :: r
  real(8), allocatable, dimension(:,:) :: phi
  real(8), allocatable, dimension(:,:,:) :: phi0
  real(8), allocatable, dimension(:,:,:) :: phi1
  real(8), allocatable, dimension(:,:,:) :: phi2
  real(8), allocatable, dimension(:,:,:) :: phi3
  integer, dimension(3) :: nptaux
  real(8), dimension(3) :: l
  real(8) :: spread
  integer, dimension(3) :: i,j
  integer :: nperiodic,n,m,p
  real(8), parameter :: pi=3.141592653589793d0
  !
  integer :: ierr
  CHARACTER(len=5) :: subname="afcc"
  !
  interface 
    !
    function fft1d(in,sign)
      complex(8), intent(in), dimension(:) :: in
      integer, intent(in) :: sign
      complex(8), dimension(size(in)) :: fft1d
    end function
    !
    function fft2d(in,sign)
      complex(8), intent(in), dimension(:,:) :: in
      integer, intent(in) :: sign
      complex(8), dimension(size(in,1),size(in,2)) :: fft2d
    end function
    !
    function fft3d(in,sign)
      complex(8), intent(in), dimension(:,:,:) :: in
      integer, intent(in) :: sign
      complex(8), dimension(size(in,1),size(in,2),size(in,3)) :: fft3d
    end function
    !
    function gaussianl(sigma,g,z)
      real(8), intent(in) :: sigma
      real(8), intent(in) :: z,g
      real(8) :: gaussianl
    end function
    !
    function phi0d(s,a,npt)
      real(8), intent(in) :: s
      integer, intent(in), dimension(3) ::  npt
      real(8), intent(in), dimension(3,3) :: a
      real(8), dimension(npt(1),npt(2),npt(3)) :: phi0d
    end function
    !
    function phi1d(s,a,npt)
      real(8), intent(in) :: s
      integer, intent(in), dimension(3) ::  npt
      real(8), intent(in), dimension(3,3) :: a
      real(8), dimension(npt(1),npt(2),npt(3)) :: phi1d
    end function
    !
    function phi2d(s,a,npt)
      real(8), intent(in) :: s
      integer, intent(in), dimension(3) ::  npt
      real(8), intent(in), dimension(3,3) :: a
      real(8), dimension(npt(1),npt(2),npt(3)) :: phi2d
    end function
    !
    function phi3d(s,a,npt)
      real(8), intent(in) :: s
      integer, intent(in), dimension(3) ::  npt
      real(8), intent(in), dimension(3,3) :: a
      real(8), dimension(npt(1),npt(2),npt(3)) :: phi3d
    end function
    !
    function volume1(a)
      real(8), intent(in), dimension(3,3) :: a
      real(8) :: volume1
    end function
    !
  end interface
  !
  l(1)=sqrt(sum(a(1:3,1)**2))
  l(2)=sqrt(sum(a(1:3,2)**2))
  l(3)=sqrt(sum(a(1:3,3)**2))
  e(1:3,1)=a(1:3,1)/l(1)
  e(1:3,2)=a(1:3,2)/l(2)
  e(1:3,3)=a(1:3,3)/l(3)
  !
  if (present(spreadopt)) then
    spread=spreadopt
  else
    spread=minval(l(1:3)/npt(1:3))*5.d0
  endif
  !
  nperiodic=count(tperiodic)
  !
#ifdef __AFC90_DEBUG
  print *, '#lattice vectors'
  print *, '#',a(1:3,1)
  print *, '#',a(1:3,2)
  print *, '#',a(1:3,3)
  print *, '#volume1'
  write(6,*) volume1(a), "volume1", a, lbound(a), ubound(a)
  !print *, '#',volume1(a)
  print *, '#periodicity'
  print *, '#',tperiodic(1:3)
  print *, '#grid'
  print *, '#',npt(1:3)
  print *, '#periodicity dimension'
  print *, '#',nperiodic
  print *, '#Gaussian spread'
  print *, '#',spread
  print *, '#lattice unit vectors'
  print *, '#',e(1:3,1)
  print *, '#',e(1:3,2)
  print *, '#',e(1:3,3)
#endif
  !
  if (nperiodic.eq.0) then
    !
    allocate(phi0(npt(1),npt(2),npt(3)))
    allocate(phi3(npt(1),npt(2),npt(3)))
    !write(6,*) "allocated phi0phi3"
    phi0=phi0d(spread,a,npt)
    !write(6,*) "called phi0"
    phi3=phi3d(spread,a,npt)
    !write(6,*) "called phi3"
    afc=phi0-phi3+pi/volume1(a)*spread*spread
#ifdef __AFC90_DEBUG
    print *, '#phi3(0)', phi3(1,1,1)
    print *, '#phi0(0)', phi0(1,1,1)
    print *, '#afc(0)', afc(1,1,1)
    print *, 'npt', npt
#endif
    deallocate(phi0,phi3)
    !
  elseif (nperiodic.eq.1) then
    !
    do n=1,3
      if (tperiodic(n)) i(3)=n
    enddo
    i(1)=mod(i(3),3)+1
    i(2)=mod(i(1),3)+1
#ifdef __AFC90_DEBUG
    print *, '#re-indexing'
    print *, '#',i
#endif
    allocate(phi1(npt(1),npt(2),npt(3)))
    allocate(phi3(npt(1),npt(2),npt(3)))
    do n=1,3
      aaux(1:3,n)=a(1:3,i(n))
      !write(6,*) "done", n
      nptaux(n)=npt(i(n))
    enddo
    phi3=phi1d(spread,aaux,nptaux)
    !write(6,*) "done phi3", n
    do m=1,npt(i(1))
      do n=1,npt(i(2))
        do p=1,npt(i(3))
           j(i(1))=m
           j(i(2))=n
           j(i(3))=p
           phi1(j(1),j(2),j(3))=phi3(m,n,p)
        enddo
      enddo
    enddo
    phi3=phi3d(spread,a,npt)
    afc=phi1-phi3+pi/volume1(a)*spread*spread
#ifdef __AFC90_DEBUG
    print *, '#phi1(0)', phi1(1,1,1)
    print *, '#phi3(0)', phi3(1,1,1)
    print *, '#afc(0)', afc(1,1,1)+2.d0*log(l(i(1)))-2.d0*log(spread)
#endif
    deallocate(phi1,phi3, STAT=ierr)
    IF(ierr/=0) call errore(subname,"deallocating phi1,phi3", abs(ierr))
    !
  elseif (nperiodic.eq.2) then
    !
    do n=1,3
      if (.not.tperiodic(n)) i(3)=n
    enddo
    i(1)=mod(i(3),3)+1
    i(2)=mod(i(1),3)+1
#ifdef __AFC90_DEBUG
    print *, '#re-indexing'
    print *, '#',i
#endif
    allocate(phi2(npt(1),npt(2),npt(3)))
    allocate(phi3(npt(1),npt(2),npt(3)))
    do n=1,3
      aaux(1:3,n)=a(1:3,i(n))
      nptaux(n)=npt(i(n))
    enddo
    phi3=phi2d(spread,aaux,nptaux)
    do m=1,npt(i(1))
      do n=1,npt(i(2))
        do p=1,npt(i(3))
           j(i(1))=m
           j(i(2))=n
           j(i(3))=p
           phi2(j(1),j(2),j(3))=phi3(m,n,p)
        enddo
      enddo
    enddo
    phi3=phi3d(spread,a,npt)
    afc=phi2-phi3+pi/volume1(a)*spread*spread
#ifdef __AFC90_DEBUG
    print *, '#phi2(0)', phi2(1,1,1)
    print *, '#phi3(0)', phi3(1,1,1)
    print *, '#afc(0)', afc(1,1,1)
#endif
    deallocate(phi2,phi3)
    !
  endif
  !
!   do m=1,npt(1)
!     do n=1,npt(2)
!       do p=1,npt(3)
!         !
!         print *, afc(m,n,p)
!         !
!       enddo
!     enddo
!   enddo
  !
end function afc
