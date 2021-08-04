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
function ydot(sigma,t,x)
  !
  implicit none
  !
  real(8), intent(in) :: sigma,t
  real(8), intent(in), dimension(2) :: x
  real(8), dimension(2) :: ydot
  real(8), parameter :: pi=3.141592653589793d0
  !
  interface 
    !
    function gaussian(sigma,x)
      real(8), intent(in) :: sigma, x
      real(8) :: gaussian
    end function
    !
  end interface
  !
  ydot(1)=x(2)
  ydot(2)=(x(1)-2.d0*pi*gaussian(sigma,exp(t)))*exp(2*t)
  !
end function ydot
