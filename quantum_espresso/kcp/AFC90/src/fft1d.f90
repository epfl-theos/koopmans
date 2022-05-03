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
function fft1d(in,sign)
  !
  implicit none
  !
  complex(8), intent(in), dimension(:) :: in
  integer, intent(in) :: sign
  complex(8), dimension(size(in)) :: fft1d
  complex(8), dimension(size(in)) :: aux
  integer(8) :: plan
  integer :: n
  !
  include "fftw3.f"
  !
  n=size(in)
  !
  if (sign>0) then
    !
    call dfftw_plan_dft_1d(plan,n,in,aux,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)
    !
  else
    !
    call dfftw_plan_dft_1d(plan,n,in,aux,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)
    aux=aux/n
    !
  endif
  !
  fft1d=aux
  !
end function fft1d
