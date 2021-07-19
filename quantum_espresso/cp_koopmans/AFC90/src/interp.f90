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
function pinterp(x,side,bound)
  !
  ! interpolation polynomial satifying
  ! p'(0)=0,p'(1)=0,p(0)=1,p(1)=0 for side=0,bound=0
  ! p'(0)=0,p'(1)=0,p(0)=0,p(1)=1 for side=1,bound=0
  ! p'(1)=0,p(0)=0,p(1)=1 for side=1,bound=-1
  ! p'(1)=0,p(0)=1,p(1)=0 for side=0,bound=-1
  ! p'(0)=0,p(0)=0,p(1)=1 for side=1,bound=1
  ! p'(0)=0,p(0)=1,p(1)=0 for side=0,bound=1
  !
  implicit none
  !
  real(8) :: pinterp
  real(8) :: x
  integer :: side
  integer :: bound
  !
  if( bound == 0 .and. side == 1 ) then
    pinterp = 3.d0 * x * x - 2.d0 * x * x * x
  else if( bound == 0 .and. side == 0 ) then
    pinterp = 1.d0 - 3.d0 * x * x + 2.d0 * x * x * x
  else if( bound == - 1 .and. side == 0 ) then
    pinterp = 1.d0 - 2.d0 * x + x * x
  else if( bound == - 1 .and. side ==  1 ) then
    pinterp = 2.d0 * x - x * x
  else if( bound == 1 .and. side == 1 ) then
    pinterp = x * x
  else if( bound == 1 .and. side == 0 ) then
    pinterp =  1 - x * x
  end if
  !
  return
  !
end function pinterp

function dpinterp(x,side,bound)
  !
  ! derivative of pinterp
  !
  implicit none
  !
  real(8) :: dpinterp
  real(8) :: x
  integer    :: side
  integer    :: bound
  !
  if( bound == 0 .and. side == 1 ) then
    dpinterp = 6.d0 * x  - 6.d0 * x * x
  else if( bound == 0 .and. side == 0 ) then
    dpinterp = - 6.d0 * x + 6.d0 * x * x
  else if( bound == - 1 .and. side == 0 ) then
    dpinterp = - 2.d0  + 2 * x
  else if( bound == - 1 .and. side ==  1 ) then
    dpinterp = 2.d0 - 2 * x
  else if( bound == 1 .and. side == 1 ) then
    dpinterp = 2 * x
  else if( bound == 1 .and. side == 0 ) then
    dpinterp =  - 2 * x
  end if
  !
  return
  !
end function dpinterp

function qinterp(x,side,bound)
  !
  ! interpolation polynomial satifying
  ! q'(0)=1,q'(1)=0,q(0)=0,q(1)=0 for side=0,bound=0
  ! q'(0)=0,q'(1)=1,q(0)=0,q(1)=0 for side=1,bound=0
  ! q'(1)=1,q(0)=0,q(1)=0 for side=1,bound=-1
  ! q'(0)=1,q(0)=0,q(1)=0 for side=0,bound=-1
  ! q'(1)=1,q(0)=0,q(1)=0 for side=1,bound=1
  ! q'(0)=1,q(0)=0,q(1)=0 for side=0,bound=1
  !
  implicit none
  !
  real(8) :: qinterp
  real(8) :: x
  integer :: side
  integer :: bound
  !
  if( bound == 0 .and. side == 1 ) then
    qinterp = - x * x + x * x * x
  else if( bound == 0 .and. side == 0 ) then
    qinterp = x - 2.d0 * x * x + x * x * x
  else if( bound == - 1 .and. side == 0 ) then
    qinterp = 0.d0
  else if( bound == 1 .and. side == 1 ) then
    qinterp =  0.d0
  else if( bound == - 1 .and. side == 1 ) then
    qinterp = - x + x * x
  else if( bound == 1 .and. side == 0 ) then
    qinterp = x - x * x
  end if
  !
  return
  !
end function qinterp

function dqinterp(x,side,bound)
  !
  ! derivative of qinterp
  !
  implicit none
  !
  real(8) :: dqinterp
  real(8) :: x
  integer :: side
  integer :: bound
  !
  if( bound == 0 .and. side == 1 ) then
    dqinterp = - 2 * x  + 3 * x * x
  else if( bound == 0 .and. side == 0 ) then
    dqinterp = 1 - 4.d0 * x  + 3 * x * x
  else if( bound == - 1 .and. side == 0 ) then
    dqinterp = 0.d0
  else if( bound == 1 .and. side == 1 ) then
    dqinterp =  0.d0
  else if( bound == - 1 .and. side == 1 ) then
    dqinterp = - 1 + 2 * x
  else if( bound == 1 .and. side == 0 ) then
    dqinterp = 1 - 2 * x
  end if
  !
  return
  !
end function dqinterp
