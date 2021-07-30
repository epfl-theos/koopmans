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
function gaussian(sigma,x)
  !
  implicit none
  !
  real(8), intent(in) :: sigma
  real(8), intent(in) :: x
  real(8) :: gaussian
  real(8), parameter :: pi=3.141592653589793d0
  real(8), parameter :: sigmathr=1.d-50
  !
  if (sigma>sigmathr) then
    gaussian=1.d0/pi/sigma/sigma*exp(-x*x/sigma/sigma)
  else
    gaussian=0.d0
  endif
  !
end function gaussian
