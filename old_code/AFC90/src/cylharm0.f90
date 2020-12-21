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
function cylharm0(sigma)
  !
  implicit none
  !
  real(8), intent(in) :: sigma
  real(8) :: cylharm0
  !
  interface 
    !
    function ei(x)
      real(8), intent(in) :: x
      real(8) :: ei
    end function
    !
  end interface
  !
  cylharm0=-0.5d0*exp(sigma*sigma/4.d0)*ei(-sigma*sigma/4.d0)
  !
end function cylharm0
