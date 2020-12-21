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
function phi0d(s,a,npt)
  !
  implicit none
  !
  real(8), intent(in) :: s
  integer, intent(in), dimension(3) :: npt
  real(8), intent(in), dimension(3,3) :: a
  real(8), dimension(npt(1),npt(2),npt(3)) :: phi0d
  real(8), parameter :: pi=3.141592653589793d0
  real(8), parameter :: rthr=1.d-10
  real(8) :: r
  integer :: i,j,k,n,m,p
  !
  interface 
    !
    function nfft(k,n)
      integer, intent(in) :: k,n
      integer :: nfft
    end function
    !
    function volume1(a)
      real(8), intent(in), dimension(3,3) :: a
      real(8) :: volume1
    end function
    !
  end interface
  !
  do i=1,npt(1)
    m=nfft(i,npt(1))
    do j=1,npt(2)
      n=nfft(j,npt(2))
      do k=1,npt(3)
        p=nfft(k,npt(3))
        r=sqrt(sum((m*a(1:3,1)/npt(1)+n*a(1:3,2)/npt(2)+p*a(1:3,3)/npt(3))**2))
        if (r>rthr) then
          phi0d(i,j,k)=erf(r/s)/r
        else
          phi0d(i,j,k)=2.d0/sqrt(pi)/s
        endif
      enddo
    enddo
  enddo
  !
end function phi0d
