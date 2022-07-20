# 1 "iotk_stream.spp"
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

# 28 "iotk_stream.spp"
#include "iotk_auxmacros.h"
# 30 "iotk_stream.spp"

# 33 "iotk_stream.spp"

subroutine iotk_stream_read_x(unit,header,setpos,getpos,ierr)
  use iotk_base
  use iotk_stream_interf
  implicit none
  integer,                                    intent(in)  :: unit
  integer(iotk_header_kind),                  intent(out) :: header
  integer,                          optional, intent(in)  :: setpos
  integer,                          optional, intent(out) :: getpos
  integer,                          optional, intent(out) :: ierr
  integer :: aa(1)
  call iotk_stream_read(unit,header,aa,setpos,getpos,.true.,ierr)
end subroutine iotk_stream_read_x


# 51 "iotk_stream.spp"

#ifdef __IOTK_LOGICAL1
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_LOGICAL1(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  LOGICAL(kind=iotk_LOGICAL1),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_LOGICAL1
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_LOGICAL2
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_LOGICAL2(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  LOGICAL(kind=iotk_LOGICAL2),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_LOGICAL2
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_LOGICAL3
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_LOGICAL3(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  LOGICAL(kind=iotk_LOGICAL3),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_LOGICAL3
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_LOGICAL4
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_LOGICAL4(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  LOGICAL(kind=iotk_LOGICAL4),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_LOGICAL4
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_INTEGER1
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_INTEGER1(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  INTEGER(kind=iotk_INTEGER1),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_INTEGER1
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_INTEGER2
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_INTEGER2(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  INTEGER(kind=iotk_INTEGER2),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_INTEGER2
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_INTEGER3
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_INTEGER3(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  INTEGER(kind=iotk_INTEGER3),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_INTEGER3
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_INTEGER4
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_INTEGER4(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  INTEGER(kind=iotk_INTEGER4),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_INTEGER4
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_REAL1
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_REAL1(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  REAL(kind=iotk_REAL1),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_REAL1
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_REAL2
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_REAL2(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  REAL(kind=iotk_REAL2),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_REAL2
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_REAL3
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_REAL3(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  REAL(kind=iotk_REAL3),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_REAL3
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_REAL4
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_REAL4(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  REAL(kind=iotk_REAL4),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_REAL4
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_COMPLEX1
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_COMPLEX1(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  COMPLEX(kind=iotk_COMPLEX1),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_COMPLEX1
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_COMPLEX2
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_COMPLEX2(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  COMPLEX(kind=iotk_COMPLEX2),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_COMPLEX2
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_COMPLEX3
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_COMPLEX3(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  COMPLEX(kind=iotk_COMPLEX3),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_COMPLEX3
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_COMPLEX4
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_COMPLEX4(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  COMPLEX(kind=iotk_COMPLEX4),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_COMPLEX4
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_CHARACTER1
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_CHARACTER1(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  CHARACTER(kind=iotk_CHARACTER1,len=*),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_CHARACTER1
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_CHARACTER2
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_CHARACTER2(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  CHARACTER(kind=iotk_CHARACTER2,len=*),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_CHARACTER2
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_CHARACTER3
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_CHARACTER3(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  CHARACTER(kind=iotk_CHARACTER3,len=*),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_CHARACTER3
#endif
# 51 "iotk_stream.spp"

#ifdef __IOTK_CHARACTER4
# 54 "iotk_stream.spp"
subroutine iotk_stream_read_CHARACTER4(unit,header,val,setpos,getpos,noval,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,                                      intent(in)  :: unit
  integer(iotk_header_kind),                    intent(out) :: header
  CHARACTER(kind=iotk_CHARACTER4,len=*),      intent(out) :: val(:)
  integer,                            optional, intent(in)  :: setpos
  integer,                            optional, intent(out) :: getpos
  logical,                            optional, intent(in)  :: noval
  integer,                            optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec,rec1
  integer :: iostat,lpos
#endif
  logical :: lnoval
  integer :: ierrl
  lnoval = .false.
  if(present(noval)) lnoval = noval
  ierrl = 0
#ifdef __IOTK_STREAMS
  if(present(setpos)) then
    lpos=setpos
  else
    inquire(unit,pos=lpos,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 80 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 80 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  if(present(getpos)) getpos = lpos
  read(unit,pos=lpos,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 87 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 87 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,iostat=iostat) header
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 92 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 92 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(.not.lnoval) then
    read(unit,iostat=iostat) val
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 98 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 98 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  read(unit,pos=lpos+iotk_record_length+rec,iostat=iostat) rec1
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 104 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 104 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  if(rec1/=rec) then
    call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 108 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
    goto 1
  end if
#else
  header = 0
  getpos = 0
  call iotk_error_issue(ierrl,"iotk_stream_read",__FILE__,__LINE__)
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 114 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_read_CHARACTER4
#endif
# 126 "iotk_stream.spp"

# 128 "iotk_stream.spp"
subroutine iotk_stream_backspace_x(unit,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
  integer,           intent(in)  :: unit
  integer, optional, intent(out) :: ierr
#ifdef __IOTK_STREAMS
  integer(iotk_record_kind) :: rec
  integer :: pos,iostat
#endif
  integer :: ierrl
  ierrl=0
#ifdef __IOTK_STREAMS
  inquire(unit,pos=pos,iostat=iostat)
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_backspace",__FILE__,__LINE__)
# 143 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 143 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 143 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,pos=pos-iotk_record_length,iostat=iostat) rec
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_backspace",__FILE__,__LINE__)
# 148 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 148 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 148 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
  read(unit,pos=pos-2*iotk_record_length-rec,iostat=iostat)
  if(iostat/=0) then
    call iotk_error_issue(ierrl,"iotk_stream_backspace",__FILE__,__LINE__)
# 153 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 153 "iotk_stream.spp"
call iotk_error_msg(ierrl,'""')
# 153 "iotk_stream.spp"
call iotk_error_write(ierrl,"iostat",iostat)
    goto 1
  end if
#else
  call iotk_error_issue(ierrl,"iotk_stream_backspace",__FILE__,__LINE__)
# 157 "iotk_stream.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 157 "iotk_stream.spp"
call iotk_error_msg(ierrl,'Streams are not implemented')
#endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_stream_backspace_x

