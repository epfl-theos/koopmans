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
function steprk(sigma)
  !
  implicit none
  !
  real(8), intent(in) :: sigma
  real(8) :: steprk
  real(8), parameter :: x=exp(-20.d0)
  real(8), parameter :: tol=1.d-3 !afcmodified:giovanni 1.d-3
  real(8), parameter :: dt0=1.d-2 !afcmodified:giovanni 1.d-2
  integer, parameter :: nmax=5  !afcmodified:giovanni 5
  real(8) :: logslope,dy,dt,y
  integer :: n
  !
  interface 
    !
    function cylharm(dt,sigma,x)
      real(8), intent(in) :: dt,sigma,x
      real(8) :: cylharm
    end function
    !
    function cylharm0(sigma)
      real(8), intent(in) :: sigma
      real(8) :: cylharm0
    end function
    !
    function besselk(x)
      real(8), intent(in) :: x
      real(8) :: besselk
    end function
    !
    function cylharmslope(tol,dt,sigma)
      real(8), intent(in) :: tol,dt,sigma
      real(8) :: cylharmslope
    end function
    !
    function cylharmseries(a0,b0,sigma,x)
      real(8), intent(in) :: a0,b0,sigma,x
      real(8) :: cylharmseries
    end function
    !
  end interface
  !
  dt=dt0
  logslope=cylharmslope(1.d-10,dt,sigma)
  y=cylharmseries(cylharm0(sigma),0.d0,sigma,x)
  dy=abs(cylharm(dt,sigma,x)+logslope*besselk(x)-y)
  n=0
  step: do while (dy>tol*abs(y))
    n=n+1
    if (n>nmax) then
      print *, 'warning: nmax exceeded in steprk'
      print *, 'rk:', cylharm(dt,sigma,x)+logslope*besselk(x)
      print *, 'series:',cylharmseries(cylharm0(sigma),0.d0,sigma,x)
      exit step
    endif
    dt=dt/2.d0
    logslope=cylharmslope(tol,dt,sigma)
    y=cylharmseries(cylharm0(sigma),0.d0,sigma,x)
    dy=abs(cylharm(dt,sigma,x)+logslope*besselk(x)-y)
  enddo step
  !
  steprk=dt
  !
end function steprk
