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
function phi3d(s,a,npt)
  !
  implicit none
  !
  real(8), intent(in) :: s
  integer, intent(in), dimension(3) :: npt
  real(8), intent(in), dimension(3,3) :: a
  real(8), dimension(npt(1),npt(2),npt(3)) :: phi3d
  complex(8), dimension(npt(1),npt(2),npt(3)) :: aux
  real(8), parameter :: pi=3.141592653589793d0
  real(8), parameter :: gthr=1.d-10
  real(8) :: omega, g
  real(8), dimension(3,3) :: b
  integer :: i,j,k,n,m,p
  !
  interface 
    !
    function fft3d(in,sign)
      complex(8), intent(in), dimension(:,:,:) :: in
      integer, intent(in) :: sign
      complex(8), dimension(size(in,1),size(in,2),size(in,3)) :: fft3d
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
    function volume1(a)
      real(8), intent(in), dimension(3,3) :: a
      real(8) :: volume1
    end function
    !
  end interface
  !
  b=reciprocal(a)
  omega=volume1(a)
  !
  do i=1,npt(1)
    m=nfft(i,npt(1))
    do j=1,npt(2)
      n=nfft(j,npt(2))
      do k=1,npt(3)
        p=nfft(k,npt(3))
        g=sqrt(sum((m*b(1:3,1)+n*b(1:3,2)+p*b(1:3,3))**2))
        if (g.ge.gthr) then
          phi3d(i,j,k)=4.d0*pi/omega/g/g*exp(-s*s*g*g/4.d0)
        else
          phi3d(i,j,k)=0.d0
        endif
      enddo
    enddo
  enddo
  !
  aux=phi3d
  aux=fft3d(aux,1)
  phi3d=aux
  !
end function phi3d
