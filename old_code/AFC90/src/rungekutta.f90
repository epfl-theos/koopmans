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
function rungekutta(sigma,t,dt,x)
  !
  implicit none
  !
  real(8), intent(in), dimension(2) :: x
  real(8), intent(in) :: sigma,t,dt
  real(8), dimension(2) :: rungekutta
  real(8), dimension(2) :: k1,k2,k3,k4
  !
  interface 
    !
    function ydot(sigma,t,x)
      real(8), intent(in) :: sigma,t
      real(8), intent(in), dimension(2) :: x
      real(8), dimension(2) :: ydot
    end function
    !
  end interface
  !
  k1=ydot(sigma,t,x)
  k2=ydot(sigma,t+0.5d0*dt,x+k1*0.5d0*dt)
  k3=ydot(sigma,t+0.5d0*dt,x+k2*0.5d0*dt)
  k4=ydot(sigma,t+dt,x+k3*dt)
  rungekutta=x+(k1+2.d0*k2+2.d0*k3+k4)*dt/6.d0
  !
end function rungekutta
