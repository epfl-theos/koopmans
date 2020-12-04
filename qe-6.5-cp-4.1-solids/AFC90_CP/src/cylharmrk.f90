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
function cylharmrk(dt,sigma,x)
  !
  implicit none
  !
  real(8), intent(in) :: dt,sigma,x
  real(8), dimension(2) :: cylharmrk
  real(8), parameter :: tol=1.d-10
  real(8), parameter :: tmax=20.d0
  real(8), dimension(2) :: y,yd
  real(8) :: t0,error,t,tx
  integer :: n
  !
  interface 
    !
    function cylharmasympt(x)
      real(8), intent(in) :: x
      real(8) :: cylharmasympt
    end function
    !
    function cylharmasymptdot(x)
      real(8), intent(in) :: x
      real(8) :: cylharmasymptdot
    end function
    !
    function cylharmasymptdotdot(x)
      real(8), intent(in) :: x
      real(8) :: cylharmasymptdotdot
    end function
    !
    function ydot(sigma,t,x)
      real(8), intent(in) :: sigma,t
      real(8), intent(in), dimension(2) :: x
      real(8), dimension(2) :: ydot
    end function
    !
    function rungekutta(sigma,t,dt,x)
      real(8), intent(in), dimension(2) :: x
      real(8), intent(in) :: sigma,t,dt
      real(8), dimension(2) :: rungekutta
    end function
    !
  end interface
  !
  tx=log(abs(x))
  t0=tx
  y=(/cylharmasympt(t0),cylharmasymptdot(t0)/)
  yd=ydot(sigma,t0,y)
  error=abs(yd(2)-cylharmasymptdotdot(t0))
  searcht0: do while (error>tol)
    if (t0+dt>tmax) then
      print *, 'warning: maximum t0 reached in subroutine cylharmrk'
      exit searcht0
    endif
    t0=t0+dt
    y=(/cylharmasympt(t0),cylharmasymptdot(t0)/)
    yd=ydot(sigma,t0,y)
    error=abs(yd(2)-cylharmasymptdotdot(t0))
  enddo searcht0
  y=(/cylharmasympt(t0),cylharmasymptdot(t0)/)
  !
  t=t0
  rk: do while (t>tx)
    y=rungekutta(sigma,t,-dt,y)
    t=t-dt
  enddo rk
  !
  cylharmrk=y
  !
end function cylharmrk
