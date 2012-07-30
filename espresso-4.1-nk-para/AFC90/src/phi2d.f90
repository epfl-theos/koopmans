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
function phi2d(s,a,npt)
  !
  implicit none
  !
  real(8), intent(in) :: s
  integer, intent(in), dimension(3) :: npt
  real(8), intent(in), dimension(3,3) :: a
  real(8), dimension(npt(1),npt(2),npt(3)) :: phi2d
  real(8), dimension(3,3) :: b
  real(8), dimension(3) :: eperp,apara
  real(8) :: z,surface,g,phase
  complex(8), dimension(npt(1),npt(2)) :: aux
  complex(8), dimension(npt(1),npt(2),npt(3)) :: phi
  integer :: i,j,k,m,n,p
  !
  interface 
    !
    function fft2d(in,sign)
      complex(8), intent(in), dimension(:,:) :: in
      integer, intent(in) :: sign
      complex(8), dimension(size(in,1),size(in,2)) :: fft2d
    end function
    !
    function gaussianl(sigma,g,z)
      real(8), intent(in) :: sigma
      real(8), intent(in) :: z,g
      real(8) :: gaussianl
    end function
    !
    function nfft(k,n)
      integer, intent(in) :: k,n
      integer :: nfft
    end function
    !
    function reciprocal(a)
      real(8), intent(in), dimension(3,3) :: a
      real(8), dimension(3,3) :: reciprocal
    end function
    !
    function vectorproduct(u,v)
      real(8), intent(in), dimension(3) :: u,v
      real(8), dimension(3) :: vectorproduct
    end function
    !
  end interface
  !
  eperp=vectorproduct(a(1:3,1),a(1:3,2))
  surface=sqrt(abs(sum(eperp**2)))
  if (sum(a(1:3,3)*eperp(1:3))>0.d0) then
    eperp=eperp/surface
  else
    eperp=-eperp/surface
  endif
  apara(1:3)=a(1:3,3)-sum(a(1:3,3)*eperp(1:3))*eperp(1:3)
  !
  b=a
  b(1:3,3)=a(1:3,3)-apara(1:3)
  b=reciprocal(b)
  !
  do i=1,npt(1)
    m=nfft(i,npt(1))
    do j=1,npt(2)
      n=nfft(j,npt(2))
      g=sqrt(sum((m*b(1:3,1)+n*b(1:3,2))**2))
      phase=sum((m*b(1:3,1)+n*b(1:3,2))*apara(1:3))
      do k=1,npt(3)
        p=nfft(k,npt(3))
        z=p*sum(eperp(1:3)*a(1:3,3))/npt(3)
        phi(i,j,k)=gaussianl(s,g,z)*cmplx(cos(p*phase/npt(3)),sin(p*phase/npt(3)))
      enddo
    enddo
  enddo
  !
  do k=1,npt(3)
    aux=phi(:,:,k)
    aux=fft2d(aux,1)
    phi2d(:,:,k)=aux
  enddo
  !
  phi2d=phi2d/surface
  !
end function phi2d
