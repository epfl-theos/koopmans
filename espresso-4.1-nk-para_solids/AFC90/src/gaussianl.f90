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
function gaussianl(sigma,g,z)
  !
  implicit none
  !
  real(8), intent(in) :: sigma
  real(8), intent(in) :: z,g
  real(8) :: gaussianl
  real(8), parameter :: pi=3.141592653589793d0
  real(8), parameter :: sigmathr=1.d-50
  real(8), parameter :: gthr=1.d-50
  real(8), parameter :: expthr=500.d0
  !
  if (g>gthr) then
    if (g*abs(z)<expthr) then
      if (sigma>sigmathr) then
        gaussianl=pi/g*(exp(g*z)*(1.d0-erf(z/sigma+sigma*g/2.d0)) &
                       +exp(-g*z)*(1.d0+erf(z/sigma-sigma*g/2.d0))) 
      else
        gaussianl=2.d0*pi/g*exp(-g*abs(z))
      endif
      !
    else
      print *, 'warning: expthr exceeded in gaussianl'
      gaussianl=0.d0
    endif
  else
    if (sigma>sigmathr) then
      gaussianl=-2.d0*pi*(z*erf(z/sigma)+sigma/sqrt(pi)*exp(-z*z/sigma/sigma))
    else
      gaussianl=-2.d0*pi*abs(z)
    endif
  endif
  !
end function gaussianl
