# 1 "iotk_unit_list.spp"
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004-2006 Giovanni Bussi
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 29 "iotk_unit_list.spp"
#include "iotk_auxmacros.h"
# 31 "iotk_unit_list.spp"

module iotk_unit_list_module
use iotk_base
implicit none

# 35 "iotk_unit_list.spp"
type iotk_unit_list
# 35 "iotk_unit_list.spp"
  type (iotk_unit_list), pointer :: next
# 35 "iotk_unit_list.spp"
  type (iotk_unit),            pointer :: ptr
# 35 "iotk_unit_list.spp"
end type iotk_unit_list

logical,              save :: iotk_units_init = .false.
type(iotk_unit_list), save :: iotk_units

contains


# 42 "iotk_unit_list.spp"
subroutine iotk_unit_list_init(list)
# 42 "iotk_unit_list.spp"
  type (iotk_unit_list), intent(out) :: list
# 42 "iotk_unit_list.spp"
  nullify(list%ptr)
# 42 "iotk_unit_list.spp"
  nullify(list%next)
# 42 "iotk_unit_list.spp"
end subroutine iotk_unit_list_init
# 42 "iotk_unit_list.spp"

# 42 "iotk_unit_list.spp"
subroutine iotk_unit_list_destroy(list)
# 42 "iotk_unit_list.spp"
  type (iotk_unit_list), intent(inout) :: list
# 42 "iotk_unit_list.spp"
  type (iotk_unit_list), pointer       :: this,next
# 42 "iotk_unit_list.spp"
  if(.not.associated(list%next)) return
# 42 "iotk_unit_list.spp"
  this=>list%next
# 42 "iotk_unit_list.spp"
  do
# 42 "iotk_unit_list.spp"
    if(associated(this%ptr))deallocate(this%ptr)
# 42 "iotk_unit_list.spp"
    next=>this%next
# 42 "iotk_unit_list.spp"
    deallocate(this)
# 42 "iotk_unit_list.spp"
    if(.not.associated(next)) exit
# 42 "iotk_unit_list.spp"
    this=>next
# 42 "iotk_unit_list.spp"
  end do
# 42 "iotk_unit_list.spp"
end subroutine iotk_unit_list_destroy
# 42 "iotk_unit_list.spp"

# 42 "iotk_unit_list.spp"
subroutine iotk_unit_list_add(list,ptr)
# 42 "iotk_unit_list.spp"
  type (iotk_unit_list), intent(inout) :: list
# 42 "iotk_unit_list.spp"
  type (iotk_unit),            pointer       :: ptr
# 42 "iotk_unit_list.spp"
  type (iotk_unit_list), pointer       :: this
# 42 "iotk_unit_list.spp"
  allocate(this)
# 42 "iotk_unit_list.spp"
  this%next => list%next
# 42 "iotk_unit_list.spp"
  list%next => this
# 42 "iotk_unit_list.spp"
  allocate(this%ptr)
# 42 "iotk_unit_list.spp"
  ptr => this%ptr
# 42 "iotk_unit_list.spp"
end subroutine iotk_unit_list_add
# 42 "iotk_unit_list.spp"

# 42 "iotk_unit_list.spp"
subroutine iotk_unit_list_del(list,ptr)
# 42 "iotk_unit_list.spp"
  type (iotk_unit_list), intent(inout) :: list
# 42 "iotk_unit_list.spp"
  type (iotk_unit),            pointer       :: ptr
# 42 "iotk_unit_list.spp"
  type (iotk_unit_list), pointer       :: this,next_save
# 42 "iotk_unit_list.spp"
  if(.not.associated(list%next)) return
# 42 "iotk_unit_list.spp"
  if(associated(list%next%ptr,ptr)) then
# 42 "iotk_unit_list.spp"
    deallocate(list%next%ptr)
# 42 "iotk_unit_list.spp"
    next_save => list%next%next
# 42 "iotk_unit_list.spp"
    deallocate(list%next)
# 42 "iotk_unit_list.spp"
    list%next => next_save
# 42 "iotk_unit_list.spp"
    nullify(ptr)
# 42 "iotk_unit_list.spp"
    return
# 42 "iotk_unit_list.spp"
  end if
# 42 "iotk_unit_list.spp"
  this => list%next
# 42 "iotk_unit_list.spp"
  do
# 42 "iotk_unit_list.spp"
    if(.not.associated(this%next)) return
# 42 "iotk_unit_list.spp"
    if(associated(this%next%ptr,ptr)) exit
# 42 "iotk_unit_list.spp"
    this => this%next
# 42 "iotk_unit_list.spp"
  end do
# 42 "iotk_unit_list.spp"
  deallocate(this%next%ptr)
# 42 "iotk_unit_list.spp"
  next_save => this%next%next
# 42 "iotk_unit_list.spp"
  deallocate(this%next)
# 42 "iotk_unit_list.spp"
  this%next => next_save
# 42 "iotk_unit_list.spp"
  nullify(ptr)
# 42 "iotk_unit_list.spp"
end subroutine iotk_unit_list_del
# 42 "iotk_unit_list.spp"

# 42 "iotk_unit_list.spp"

# 42 "iotk_unit_list.spp"
 subroutine iotk_unit_list_search(list,ptr,unit)
# 42 "iotk_unit_list.spp"
   type (iotk_unit_list), intent(in) :: list
# 42 "iotk_unit_list.spp"
   type (iotk_unit),            pointer    :: ptr
# 42 "iotk_unit_list.spp"
 
# 42 "iotk_unit_list.spp"
  integer, optional,intent(in) :: unit
# 42 "iotk_unit_list.spp"

# 42 "iotk_unit_list.spp"
   type (iotk_unit_list), pointer    :: this
# 42 "iotk_unit_list.spp"
   nullify(ptr)
# 42 "iotk_unit_list.spp"
   this => list%next
# 42 "iotk_unit_list.spp"
   if(.not.associated(this)) return
# 42 "iotk_unit_list.spp"
   do
# 42 "iotk_unit_list.spp"
     if(.not.associated(this%ptr)) goto 1000
# 42 "iotk_unit_list.spp"

# 42 "iotk_unit_list.spp"

# 42 "iotk_unit_list.spp"
     if(present(unit)) then
# 42 "iotk_unit_list.spp"
       if(this%ptr%unit /= unit) goto 1000
# 42 "iotk_unit_list.spp"
     end if
# 42 "iotk_unit_list.spp"

# 42 "iotk_unit_list.spp"

# 42 "iotk_unit_list.spp"
     ptr => this%ptr
# 42 "iotk_unit_list.spp"
     exit
# 42 "iotk_unit_list.spp"
1000 continue
# 42 "iotk_unit_list.spp"
     if(.not.associated(this%next)) exit
# 42 "iotk_unit_list.spp"
     this => this%next
# 42 "iotk_unit_list.spp"
   end do
# 42 "iotk_unit_list.spp"
 end subroutine iotk_unit_list_search
# 42 "iotk_unit_list.spp"
 

end module iotk_unit_list_module
