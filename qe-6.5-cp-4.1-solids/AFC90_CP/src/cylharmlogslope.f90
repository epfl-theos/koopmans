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
function cylharmlogslope(dt,sigma,x)
  !
  implicit none
  !
  real(8), intent(in) :: dt,sigma,x
  real(8) :: cylharmlogslope
  real(8), dimension(2) :: y
  !
  interface 
    !
    function cylharmrk(dt,sigma,x)
      real(8), intent(in) :: dt,sigma,x
      real(8), dimension(2) :: cylharmrk
    end function
    !
  end interface
  !
  y=cylharmrk(dt,sigma,x)
  !
  cylharmlogslope=y(2)
  !
end function cylharmlogslope
