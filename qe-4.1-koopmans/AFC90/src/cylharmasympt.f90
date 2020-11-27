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
function cylharmasympt(x)
  !
  implicit none
  !
  real(8), intent(in) :: x
  real(8) :: cylharmasympt
  real(8), parameter :: pi=3.141592653589793d0
  !
  cylharmasympt=sqrt(0.5d0*pi)*exp(-exp(x))*exp(-0.5d0*x)
  !
end function cylharmasympt
function cylharmasymptdot(x)
  !
  implicit none
  !
  real(8), intent(in) :: x
  real(8) :: cylharmasymptdot
  real(8), parameter :: pi=3.141592653589793d0
  !
  cylharmasymptdot=sqrt(0.5d0*pi)*exp(-exp(x))*exp(-0.5d0*x)*(-exp(x)-0.5d0)
  !
end function cylharmasymptdot
function cylharmasymptdotdot(x)
  !
  implicit none
  !
  real(8), intent(in) :: x
  real(8) :: cylharmasymptdotdot
  real(8), parameter :: pi=3.141592653589793d0
  !
  cylharmasymptdotdot=sqrt(0.5d0*pi)*exp(-exp(x))*exp(-0.5d0*x)*((-exp(x)-0.5d0)**2-exp(x))
  !
end function cylharmasymptdotdot
