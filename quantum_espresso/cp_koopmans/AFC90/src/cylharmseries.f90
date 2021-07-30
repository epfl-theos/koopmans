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
function cylharmseries(a0,b0,sigma,x)
  !
  implicit none
  !
  real(8), intent(in) :: a0,b0,sigma,x
  real(8) :: cylharmseries
  real(8), parameter :: tol=1d-20
  real(8), parameter :: pi=3.141592653589793d0
  real(8), parameter :: sigmathr=1.d-50
  integer, parameter :: nmax=50
  integer :: n
  real(8) :: an,bn,cn,dy
  !
  n=0
  an=a0
  bn=b0
  cn=0.d0
  if (sigma>sigmathr) cn=-2.d0/sigma/sigma
  dy=a0+b0*log(abs(x))
  cylharmseries=dy
  !
  series: do while (n<nmax)
    !
    n=n+1
    bn=bn/4.d0/n/n
    an=(an+cn)/4.d0/n/n-bn/n
    if (sigma>sigmathr) cn=-cn/n/sigma/sigma
    dy=an+bn*log(abs(x))
    cylharmseries=cylharmseries/x/x+dy
    if ((abs(dy)<tol*abs(cylharmseries)).or.(cylharmseries.eq.0.d0)) exit series
    !
  enddo series
  !if (n.ge.nmax) then
  !  print *, 'nmax exceeded in cylharmseries'
  !endif
  !
  cylharmseries=cylharmseries*x**n*x**n
  !
end function cylharmseries
