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
function cylharmslope(tol,dt,sigma)
  !
  implicit none
  !
  real(8), intent(in) :: tol,dt,sigma
  real(8) :: cylharmslope
  integer, parameter :: nmax=100 !afcmodified:giovanni 100
  real(8) :: t,x,dlogslope
  integer :: n
  !
  interface 
    !
    function cylharmlogslope(dt,sigma,x)
      real(8), intent(in) :: dt,sigma,x
      real(8) :: cylharmlogslope
    end function
    !
  end interface
  !
  t=0.d0
  x=exp(t)
  cylharmslope=cylharmlogslope(dt,sigma,x)
  t=t-1.d0
  x=exp(t)
  dlogslope=cylharmslope
  cylharmslope=cylharmlogslope(dt,sigma,x)
  dlogslope=abs(cylharmslope-dlogslope)
  n=0
  logslope: do while (dlogslope>tol*abs(cylharmslope))
    n=n+1
    if (n>nmax) then
      print *, 'warning: nmax exceeded in subroutine cylharmslope'
      exit logslope
    endif
    t=t-1.d0
    x=exp(t)
    dlogslope=cylharmslope
    cylharmslope=cylharmlogslope(dt,sigma,x)
    dlogslope=abs(cylharmslope-dlogslope)
  enddo logslope
  !
end function cylharmslope
