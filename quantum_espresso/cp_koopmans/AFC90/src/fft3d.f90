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
function fft3d(in,sign)
  !
  implicit none
  !
  complex(8), intent(in), dimension(:,:,:) :: in
  integer, intent(in) :: sign
  complex(8), dimension(size(in,1),size(in,2),size(in,3)) :: fft3d
  complex(8), dimension(size(in,3)) :: aux1
  complex(8), dimension(size(in,1),size(in,2)) :: aux2
  integer :: i,j,nx,ny,nz
  !
  interface 
    !
    function fft1d(in,sign)
      complex(8), intent(in), dimension(:) :: in
      integer, intent(in) :: sign
      complex(8), dimension(size(in)) :: fft1d
    end function
    !
    function fft2d(in,sign)
      complex(8), intent(in), dimension(:,:) :: in
      integer, intent(in) :: sign
      complex(8), dimension(size(in,1),size(in,2)) :: fft2d
    end function
    !
  end interface
  !
  nx=size(in,1)
  ny=size(in,2)
  nz=size(in,3)
  !
  fft3d=in
  !
  do i=1,nz
    aux2(1:nx,1:ny)=fft3d(1:nx,1:ny,i)
    aux2=fft2d(aux2,sign)
    fft3d(1:nx,1:ny,i)=aux2(1:nx,1:ny)
  enddo
  !
  do i=1,nx
    do j=1,ny
      aux1(1:nz)=fft3d(i,j,1:nz)
      aux1=fft1d(aux1,sign)
      fft3d(i,j,1:nz)=aux1(1:nz)
    enddo
  enddo
  !
end function fft3d
