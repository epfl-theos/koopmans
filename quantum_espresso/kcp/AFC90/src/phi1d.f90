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
function phi1d(s,a,npt)
  !
  implicit none
  !
  real(8), intent(in) :: s
  integer, intent(in), dimension(3) :: npt
  real(8), intent(in), dimension(3,3) :: a
  real(8), dimension(npt(1),npt(2),npt(3)) :: phi1d
  real(8), parameter :: pi=3.141592653589793d0
  real(8), parameter :: sthr=1.d-50
  real(8) :: sigma,g,l,phase,z
  real(8), dimension(npt(1)*npt(2)) :: x
  real(8), dimension(npt(1)*npt(2)) :: r
  real(8), dimension(3) :: epara
  complex(8), dimension(npt(1)*npt(2),npt(3)) :: phi
  complex(8), dimension(npt(3)) :: aux1
  integer :: i,j,k,n,m,p,nr
  !
  interface 
    !
    function fft1d(in,sign)
      complex(8), intent(in), dimension(:) :: in
      integer, intent(in) :: sign
      complex(8), dimension(size(in)) :: fft1d
    end function
    !
    function gaussiank(sigma,x)
      real(8), intent(in) :: sigma
      real(8), intent(in), dimension(:) :: x
      real(8), dimension(size(x)) :: gaussiank
    end function
    !
    function eimlog(x)
      real(8), intent(in) :: x
      real(8) :: eimlog
    end function
    !
    function nfft(k,n)
      integer, intent(in) :: k,n
      integer :: nfft
    end function
    !
    function rtoaxis(a1,a2,n1,n2,e3)
      real(8), intent(in), dimension(3) :: a1,a2,e3
      integer, intent(in) :: n1,n2
      real(8), dimension(n1*n2) :: rtoaxis
    end function
    !
  end interface
  !
  nr=npt(1)*npt(2)
  l=sqrt(sum(a(1:3,3)**2))
  epara=a(1:3,3)/l
  r=rtoaxis(a(1:3,1),a(1:3,2),npt(1),npt(2),epara)
  !
  phi1d=0.d0
  do i=1,npt(3)
    g=abs(2.d0*pi/l*nfft(i,npt(3)))
    if (g.eq.0.d0) then
      if (s>sthr) then
        do j=1,nr
          phi(j,i)=eimlog(-r(j)*r(j)/s/s)
        enddo
      else
        do j=1,nr
          phi(j,i)=-2.d0*log(r(j))
        enddo
      endif
    else
      x=g*r
      sigma=g*s
      phi(1:nr,i)=2.d0*gaussiank(sigma,x)
    endif
    phi(1:nr,i)=phi(1:nr,i)*exp(-g*g*s*s/4.d0)/l
  enddo
  !
  do i=1,npt(1)
    m=nfft(i,npt(1))
    do j=1,npt(2)
      n=nfft(j,npt(2))
      do k=1,npt(3)
        p=nfft(k,npt(3))
        g=2.d0*pi/l*p
        z=sum(epara(1:3)*(m*a(1:3,1)/npt(1)+n*a(1:3,2)/npt(2)))
        phase=g*z
        aux1(k)=phi(i+npt(1)*(j-1),k)*cmplx(cos(phase),sin(phase))
      enddo
      aux1=fft1d(aux1,1)
      !aux1=0.d0
      phi1d(i,j,1:npt(3))=aux1(1:npt(3))
    enddo
  enddo
  !
end function phi1d
