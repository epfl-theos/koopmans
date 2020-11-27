# 48 "iotk_attr.spp"
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
# 65 "iotk_attr.spp"


!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 74 "iotk_attr.spp"
#include "iotk_auxmacros.h"
# 76 "iotk_attr.spp"

# 78 "iotk_attr.spp"

#ifdef __IOTK_COMPLEX1
#if 0 <= __IOTK_MAXRANK

# 83 "iotk_attr.spp"
! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_COMPLEX1(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
# 92 "iotk_attr.spp"
    COMPLEX (kind=iotk_COMPLEX1), intent(out) :: out(n)
    COMPLEX (kind=iotk_COMPLEX1), intent(in)  :: in(n)
# 95 "iotk_attr.spp"
    out = in
end subroutine iotk_private_pack_COMPLEX1

# 100 "iotk_attr.spp"
subroutine iotk_write_COMPLEX1(val,string,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  COMPLEX(kind=iotk_COMPLEX1), intent(in) :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)              :: string
#else
  character(len=*), intent(out) :: string
#endif
  integer, intent(out) :: ierr
# 116 "iotk_attr.spp"
  character(len=100) :: tmpval
# 118 "iotk_attr.spp"
  integer :: index,iostat
  ierr = 0
  iostat = 0 
  string(1:1) = iotk_eos
  if(size(val)==0) return
  if(len(string)==0) then
    call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 124 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.7 ")
    return
  end if
  do index=1,size(val)
# 141 "iotk_attr.spp"
    write(tmpval,trim(iotk_wfmt("COMPLEX",kind(val),1,-1," ")),iostat=iostat) val(index)
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 143 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.7 ")
# 143 "iotk_attr.spp"
call iotk_error_msg(ierr,' ')
# 143 "iotk_attr.spp"
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
    call iotk_strcat(string,trim(adjustl(tmpval))//" ",ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 148 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.7 ")
      return
    end if
# 152 "iotk_attr.spp"
  end do
! the last blank is deleted
  string(iotk_strlen(string):iotk_strlen(string)) = iotk_eos
end subroutine iotk_write_COMPLEX1
# 158 "iotk_attr.spp"

# 162 "iotk_attr.spp"
subroutine iotk_read_COMPLEX1(val,string,index,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  COMPLEX(kind=iotk_COMPLEX1), intent(inout) :: val(:)
  character(len=*), intent(in) :: string
  integer, intent(inout) :: index
  integer, intent(out) :: ierr
# 175 "iotk_attr.spp"
  integer :: pos,pos1,iostat
  integer :: maxindex
# 178 "iotk_attr.spp"
  real(kind=iotk_COMPLEX1) :: tmpreal
  complex(kind=iotk_COMPLEX1) :: tmpcomplex
# 181 "iotk_attr.spp"
#ifdef __IOTK_WORKAROUND9
  character(len=100) :: tmpstr ! debug
#endif
  pos = 0
  pos1= 0
  ierr = 0
  iostat = 0
# 189 "iotk_attr.spp"
   maxindex = 2 * size(val)
# 193 "iotk_attr.spp"
! for the moment, commas are considered as blanks
  do
    pos = verify(string(pos1+1:)," ,")+pos1
    if(pos==pos1) exit
    pos = pos - 1
    pos1 = scan(string(pos+1:)," ,")+pos
    if(pos1==pos) pos1 = len(string) + 1
!READ string(pos+1:pos1-1)
    index = index+1
    if(index>maxindex) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 203 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.7 ")
# 203 "iotk_attr.spp"
call iotk_error_msg(ierr,'Too many data')
    end if
# 217 "iotk_attr.spp"
#ifdef __IOTK_WORKAROUND9
    tmpstr = TRIM( string(pos+1:pos1-1) )
    read( tmpstr,"(G100.95)",iostat=iostat) tmpreal
#else
    read(string(pos+1:pos1-1),"(G100.95)",iostat=iostat) tmpreal
#endif
    if(modulo(index,2)==1) then
      tmpcomplex = cmplx(tmpreal,aimag((val((index+1)/2))),kind=iotk_COMPLEX1)
    else
      tmpcomplex = cmplx(real(val((index+1)/2)),tmpreal,kind=iotk_COMPLEX1)
    end if
    val((index+1)/2) = tmpcomplex
# 236 "iotk_attr.spp"
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.7 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierr,'Error reading a COMPLEX number from string')
# 237 "iotk_attr.spp"
call iotk_error_write(ierr,"string",string(pos+1:pos1-1))
# 237 "iotk_attr.spp"
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
# 241 "iotk_attr.spp"
    if(pos1>=len(string)) exit
  end do
end subroutine iotk_read_COMPLEX1
# 246 "iotk_attr.spp"

# 249 "iotk_attr.spp"
subroutine iotk_write_attr_COMPLEX1_0(attr,name,val,dummy,first,newline,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=iotk_COMPLEX1), intent(in)  :: val 
  type(iotk_dummytype), optional :: dummy
  logical, optional, intent(in)  :: first
  logical, optional, intent(in)  :: newline
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  character :: delim
# 271 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  logical :: nl
  if(present(newline)) then
    nl = newline
  else
    nl = .false.
  endif
  ierrl = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  attlen = iotk_strlen_trim(attr)
  namlen = iotk_strlen_trim(name)
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 285 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 285 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 285 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
# 302 "iotk_attr.spp"
  delim = '"'
# 304 "iotk_attr.spp"
  call iotk_write((/val/),tmpval,ierrl)
# 308 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 309 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
    goto 1
  end if
# 313 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 315 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 315 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute dummy argument is too short')
    goto 1
  end if
  if(.not. nl) then
    attr(attlen+1:attlen+vallen+namlen+5) = " "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  else
    attr(attlen+1:attlen+vallen+namlen+len(iotk_newline)+5) &
       = iotk_newline//" "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_attr_COMPLEX1_0

# 333 "iotk_attr.spp"
subroutine iotk_scan_attr_COMPLEX1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=iotk_COMPLEX1)                        :: val 
#else
  COMPLEX(kind=iotk_COMPLEX1), intent(out)           :: val 
#endif
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  COMPLEX(kind=iotk_COMPLEX1), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 364 "iotk_attr.spp"
  integer :: index
  COMPLEX(kind=iotk_COMPLEX1), allocatable :: tmpval (:)
# 367 "iotk_attr.spp"
  ierrl = 0
  attlen=iotk_strlen(attr)
  namlen=iotk_strlen_trim(name)
  foundl = .false.
  equal = 0
  do
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) exit
    equal = equal + pos
    pos = scan(attr(equal+1:attlen),"=")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 378 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 378 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 378 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 385 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 391 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 405 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
  else
    goto 1
  end if
# 426 "iotk_attr.spp"
  allocate(tmpval(1))
  index = 0
  call iotk_str_clean(valc)
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 431 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
    goto 1
  end if
# 435 "iotk_attr.spp"
  if(index/=2*1) then
# 439 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 439 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 439 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 439 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc)
# 439 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 443 "iotk_attr.spp"
  val = tmpval(1)
# 447 "iotk_attr.spp"
  deallocate(tmpval)
# 449 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 453 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 453 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 453 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 466 "iotk_attr.spp"
    val = default
# 468 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX1_0
# 476 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_COMPLEX1_0
  write(0,*)
end subroutine iotk_attr_dummy_COMPLEX1_0

# 45 "iotk_attr.spp"

# 65 "iotk_attr.spp"


!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 74 "iotk_attr.spp"
#include "iotk_auxmacros.h"
# 76 "iotk_attr.spp"

# 78 "iotk_attr.spp"

#ifdef __IOTK_COMPLEX1
#if 1 <= __IOTK_MAXRANK

# 158 "iotk_attr.spp"

# 246 "iotk_attr.spp"

# 249 "iotk_attr.spp"
subroutine iotk_write_attr_COMPLEX1_1(attr,name,val,dummy,first,newline,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=iotk_COMPLEX1), intent(in)  :: val (:)
  type(iotk_dummytype), optional :: dummy
  logical, optional, intent(in)  :: first
  logical, optional, intent(in)  :: newline
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  character :: delim
# 271 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  logical :: nl
  if(present(newline)) then
    nl = newline
  else
    nl = .false.
  endif
  ierrl = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  attlen = iotk_strlen_trim(attr)
  namlen = iotk_strlen_trim(name)
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 285 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 285 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 285 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
# 302 "iotk_attr.spp"
  delim = '"'
# 306 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 308 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 309 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
    goto 1
  end if
# 313 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 315 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 315 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute dummy argument is too short')
    goto 1
  end if
  if(.not. nl) then
    attr(attlen+1:attlen+vallen+namlen+5) = " "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  else
    attr(attlen+1:attlen+vallen+namlen+len(iotk_newline)+5) &
       = iotk_newline//" "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_attr_COMPLEX1_1

# 333 "iotk_attr.spp"
subroutine iotk_scan_attr_COMPLEX1_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=iotk_COMPLEX1)                        :: val (:)
#else
  COMPLEX(kind=iotk_COMPLEX1), intent(out)           :: val (:)
#endif
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  COMPLEX(kind=iotk_COMPLEX1), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 364 "iotk_attr.spp"
  integer :: index
  COMPLEX(kind=iotk_COMPLEX1), allocatable :: tmpval (:)
# 367 "iotk_attr.spp"
  ierrl = 0
  attlen=iotk_strlen(attr)
  namlen=iotk_strlen_trim(name)
  foundl = .false.
  equal = 0
  do
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) exit
    equal = equal + pos
    pos = scan(attr(equal+1:attlen),"=")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 378 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 378 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 378 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 385 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 391 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 405 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
  else
    goto 1
  end if
# 426 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_str_clean(valc)
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 431 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
    goto 1
  end if
# 435 "iotk_attr.spp"
  if(index/=2*size(val)) then
# 439 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 439 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 439 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 439 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc)
# 439 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 445 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 447 "iotk_attr.spp"
  deallocate(tmpval)
# 449 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 453 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 453 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 453 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 466 "iotk_attr.spp"
    val = default
# 468 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX1_1
# 476 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_COMPLEX1_1
  write(0,*)
end subroutine iotk_attr_dummy_COMPLEX1_1

# 45 "iotk_attr.spp"

# 65 "iotk_attr.spp"


!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 74 "iotk_attr.spp"
#include "iotk_auxmacros.h"
# 76 "iotk_attr.spp"

# 78 "iotk_attr.spp"

#ifdef __IOTK_COMPLEX1
#if 2 <= __IOTK_MAXRANK

# 158 "iotk_attr.spp"

# 246 "iotk_attr.spp"

# 249 "iotk_attr.spp"
subroutine iotk_write_attr_COMPLEX1_2(attr,name,val,dummy,first,newline,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=iotk_COMPLEX1), intent(in)  :: val (:,:)
  type(iotk_dummytype), optional :: dummy
  logical, optional, intent(in)  :: first
  logical, optional, intent(in)  :: newline
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  character :: delim
# 271 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  logical :: nl
  if(present(newline)) then
    nl = newline
  else
    nl = .false.
  endif
  ierrl = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  attlen = iotk_strlen_trim(attr)
  namlen = iotk_strlen_trim(name)
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 285 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 285 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 285 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
# 302 "iotk_attr.spp"
  delim = '"'
# 306 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 308 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 309 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
    goto 1
  end if
# 313 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 315 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 315 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute dummy argument is too short')
    goto 1
  end if
  if(.not. nl) then
    attr(attlen+1:attlen+vallen+namlen+5) = " "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  else
    attr(attlen+1:attlen+vallen+namlen+len(iotk_newline)+5) &
       = iotk_newline//" "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_attr_COMPLEX1_2

# 333 "iotk_attr.spp"
subroutine iotk_scan_attr_COMPLEX1_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=iotk_COMPLEX1)                        :: val (:,:)
#else
  COMPLEX(kind=iotk_COMPLEX1), intent(out)           :: val (:,:)
#endif
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  COMPLEX(kind=iotk_COMPLEX1), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 364 "iotk_attr.spp"
  integer :: index
  COMPLEX(kind=iotk_COMPLEX1), allocatable :: tmpval (:)
# 367 "iotk_attr.spp"
  ierrl = 0
  attlen=iotk_strlen(attr)
  namlen=iotk_strlen_trim(name)
  foundl = .false.
  equal = 0
  do
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) exit
    equal = equal + pos
    pos = scan(attr(equal+1:attlen),"=")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 378 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 378 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 378 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 385 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 391 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 405 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
  else
    goto 1
  end if
# 426 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_str_clean(valc)
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 431 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
    goto 1
  end if
# 435 "iotk_attr.spp"
  if(index/=2*size(val)) then
# 439 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 439 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 439 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 439 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc)
# 439 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 445 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 447 "iotk_attr.spp"
  deallocate(tmpval)
# 449 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 453 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 453 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 453 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 466 "iotk_attr.spp"
    val = default
# 468 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX1_2
# 476 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_COMPLEX1_2
  write(0,*)
end subroutine iotk_attr_dummy_COMPLEX1_2

# 45 "iotk_attr.spp"

