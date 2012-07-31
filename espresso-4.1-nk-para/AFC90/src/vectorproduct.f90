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
function vectorproduct(u,v)
  !
  implicit none
  !
  real(8), intent(in), dimension(3) :: u,v
  real(8), dimension(3) :: vectorproduct
  !
  vectorproduct(1)=u(2)*v(3)-u(3)*v(2)
  vectorproduct(2)=u(3)*v(1)-u(1)*v(3)
  vectorproduct(3)=u(1)*v(2)-u(2)*v(1)
  !
  return
end function vectorproduct
