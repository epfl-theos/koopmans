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
function rtoaxis(a1,a2,n1,n2,e3)
  !
  implicit none
  !
  real(8), intent(in), dimension(3) :: a1,a2,e3
  integer, intent(in) :: n1,n2
  real(8), dimension(n1*n2) :: rtoaxis
  integer :: i1,i2
  real(8), dimension(3) :: u
  !
  interface 
    !
    function nfft(k,n)
      integer, intent(in) :: k,n
      integer :: nfft
    end function
    !
  end interface
  !
  do i1=1,n1
    do i2=1,n2
      u=nfft(i1,n1)*a1/n1+nfft(i2,n2)*a2/n2
      u=u-sum(e3(1:3)*u(1:3))*u
      rtoaxis(i1+n1*(i2-1))=sqrt(sum(u(1:3)**2))
    enddo
  enddo
  !
end function rtoaxis


  
