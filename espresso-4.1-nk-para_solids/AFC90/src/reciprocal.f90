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
function reciprocal(a)
  !
  implicit none
  !
  real(8), intent(in), dimension(3,3) :: a
  real(8), dimension(3,3) :: reciprocal
  real(8), parameter :: pi=3.141592653589793d0
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
  reciprocal(1:3,1)=vectorproduct(a(1:3,2),a(1:3,3))
  reciprocal(1:3,1)=reciprocal(1:3,1)/sum(reciprocal(1:3,1)*a(1:3,1))
  reciprocal(1:3,2)=vectorproduct(a(1:3,3),a(1:3,1))
  reciprocal(1:3,2)=reciprocal(1:3,2)/sum(reciprocal(1:3,2)*a(1:3,2))
  reciprocal(1:3,3)=vectorproduct(a(1:3,1),a(1:3,2))
  reciprocal(1:3,3)=reciprocal(1:3,3)/sum(reciprocal(1:3,3)*a(1:3,3))
  !
  reciprocal=reciprocal*2.d0*pi
  !
end function reciprocal
