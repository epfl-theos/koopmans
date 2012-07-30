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
function volume(a)
  !
  implicit none
  !
  real(8), intent(in), dimension(3,3) :: a
  real(8) :: volume
  !
  interface
    !
    function vectorproduct(u,v)
      real(8), intent(in), dimension(3) :: u,v
      real(8), dimension(3) :: vectorproduct
    end function
    !
  end interface 
  !
  volume=abs(sum(vectorproduct(a(1:3,1),a(1:3,2))*a(1:3,3)))
  !
end function volume
