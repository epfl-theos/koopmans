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
function ei(x)
  !
  implicit none
  !
  real(8), intent(in) :: x
  real(8) :: ei
  real(8), parameter :: tol=1d-20
  real(8), parameter :: eulergamma=0.5772156649015328606d0
  integer, parameter :: nmax=100
  integer :: n
  real(8) :: y,dy
  !
  n=1
  dy=x
  ei=log(abs(x))+eulergamma+dy
  !
  do while ((abs(dy)>tol*abs(ei)).and.(n<nmax))
    n=n+1
    dy=dy*x/n
    ei=ei+dy/n
  enddo
  !
end function ei
