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
function fft2d(in,sign)
  !
  implicit none
  !
  complex(8), intent(in), dimension(:,:) :: in
  integer, intent(in) :: sign
  complex(8), dimension(size(in,1),size(in,2)) :: fft2d
  complex(8), dimension(size(in,1),size(in,2)) :: aux
  integer(8) :: plan
  integer :: nx,ny
  !
  include "fftw3.f"
  !
  nx=size(in,1)
  ny=size(in,2)
  !
  if (sign>0) then
    !
    call dfftw_plan_dft_2d(plan,nx,ny,in,aux,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)
    !
  else
    !
    call dfftw_plan_dft_2d(plan,nx,ny,in,aux,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)
    aux=aux/nx/ny
    !
  endif
  !
  fft2d=aux
  !
end function fft2d
