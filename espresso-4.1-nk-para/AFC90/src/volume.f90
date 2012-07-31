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
real(8) function volume1(a)
   !
   implicit none
   !
   real(8), intent(IN) :: a(3,3)
   real(8) :: tmp(3)=0.d0

 interface
   function vectorproduct(u,v)
     real(8), intent(in), dimension(3) :: u,v
     real(8), dimension(3) :: vectorproduct
   end function
 end interface

   tmp=vectorproduct(a(1:3,1),a(1:3,2))
   !volume1=1.000
   volume1=abs(dot_product(tmp,a(1:3,3)))
   return
   !
end function volume1
