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
function eimlog(x)
  !
  implicit none
  !
  real(8), intent(in) :: x
  real(8) :: eimlog
  real(8), parameter :: tol=1d-20
  real(8), parameter :: xthr=25.d0
  real(8), parameter :: eulergamma=0.5772156649015328606d0
  integer, parameter :: nmax=100
  integer :: n
  real(8) :: dy
  !
  if (x>-xthr) then
    n=1
    dy=x
    eimlog=eulergamma+dy
    loop: do while ((abs(dy)>tol*abs(eimlog)).and.(n<nmax))
      n=n+1
      if (n>nmax) then
        print *, 'warning: nmax exceeded in eimlog'
        exit loop
      endif 
      dy=dy*x/n
      eimlog=eimlog+dy/n
    enddo loop
  else
    eimlog=-log(abs(x))
  endif
  !
end function eimlog
