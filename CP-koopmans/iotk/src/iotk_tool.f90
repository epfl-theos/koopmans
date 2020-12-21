# 1 "iotk_tool.spp"
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2006 Giovanni Bussi
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

# 28 "iotk_tool.spp"
#include "iotk_auxmacros.h"
# 33 "iotk_tool.spp"

# 36 "iotk_tool.spp"

# 38 "iotk_tool.spp"
subroutine iotk_tool_x(args)
  use iotk_base
  use iotk_error_interf
  use iotk_str_interf
  use iotk_tool_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  character(len=*), intent(in) :: args(:)
  integer :: iarg,ierrl
  character(iotk_linlenx) :: arg
  logical :: print_help_options,print_help_commands,print_help_basic,print_copyright,print_version
  logical :: check
  integer :: linlen,indent,maxindent
  ierrl = 0
  iarg = 1

  print_version = .false.
  print_help_options  = .false.
  print_help_commands = .false.
  print_help_basic = .false.
  print_copyright = .false.

  if(size(args)==0) then
    print_help_basic = .true.
  end if

  do iarg = 1 , size(args)
    arg = args(iarg)
    if(iotk_strcomp(arg(1:1),"-")) then
! options here
      if(iotk_strcomp(arg,"--help") .or. iotk_strcomp(arg,"-H")) then
        print_help_basic = .true.
        exit
      else if(iotk_strcomp(arg,"--version")) then
        print_version = .true.
        exit
      else if(iotk_strcomp(arg,"--do-nothing")) then
        exit
      else if(iotk_strcomp(arg,"--copyright")) then
        print_copyright = .true.
        exit
      else if(iotk_strcomp(arg,"--help-options")) then
        print_help_options = .true.
        exit
      else if(iotk_strcomp(arg,"--help-commands")) then
        print_help_commands = .true.
        exit
      else if(arg(1:13)=="--set-linlen=") then
        call iotk_atoi(linlen,arg(14:iotk_strlen(arg)),check=check)
        if(.not.check) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
# 89 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 89 "iotk_tool.spp"
call iotk_error_msg(ierrl,'')
          goto 1
        end if
        call iotk_set(linlen=linlen,ierr=ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
# 94 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 94 "iotk_tool.spp"
call iotk_error_msg(ierrl,'')
          goto 1
        end if
      else if(arg(1:13)=="--set-indent=") then
        call iotk_atoi(indent,arg(14:iotk_strlen(arg)),check=check)
        if(.not.check) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
# 100 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 100 "iotk_tool.spp"
call iotk_error_msg(ierrl,'')
          goto 1
        end if
        call iotk_set(indent=indent,ierr=ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
# 105 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 105 "iotk_tool.spp"
call iotk_error_msg(ierrl,'')
          goto 1
        end if
      else if(arg(1:16)=="--set-maxindent=") then
        call iotk_atoi(maxindent,arg(17:iotk_strlen(arg)),check=check)
        if(.not.check) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
# 111 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 111 "iotk_tool.spp"
call iotk_error_msg(ierrl,'')
          goto 1
        end if
        call iotk_set(maxindent=maxindent,ierr=ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
# 116 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 116 "iotk_tool.spp"
call iotk_error_msg(ierrl,'')
          goto 1
        end if
      else
        write(iotk_error_unit,"(a)") "unrecognized option `"//arg(1:iotk_strlen(arg))//"'"
        print_help_basic = .true.
        exit
      end if
    else
! commands here
      if(iotk_strcomp(arg,"convert")) then
        call iotk_tool_convert(args(iarg+1:),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
# 129 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 129 "iotk_tool.spp"
call iotk_error_msg(ierrl,'Error converting file')
          goto 1
        end if
      else if(iotk_strcomp(arg,"dump")) then
        call iotk_tool_dump(args(iarg+1:),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
# 135 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 135 "iotk_tool.spp"
call iotk_error_msg(ierrl,'Error dumping file')
          goto 1
        end if
      else if(iotk_strcomp(arg,"info")) then
        call iotk_tool_info(args(iarg+1:),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
# 141 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 141 "iotk_tool.spp"
call iotk_error_msg(ierrl,'Error')
          goto 1
        end if
      else if(iotk_strcomp(arg,"man")) then
        call iotk_tool_man(args(iarg+1:),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool",__FILE__,__LINE__)
# 147 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 147 "iotk_tool.spp"
call iotk_error_msg(ierrl,'Error')
          goto 1
        end if
      else
        write(iotk_error_unit,"(a)") "Unknown command: `"//arg(1:iotk_strlen(arg))//"'"
        write(iotk_error_unit,"(a)") ""
        print_help_commands = .true.
      end if
      exit
    end if
  end do

  if(print_help_basic) then
    write(iotk_error_unit,"(a)") "Usage: iotk [iotk-options] command [command-options-and-arguments]"
    write(iotk_error_unit,"(a)") "  where iotk-options are ..."
    write(iotk_error_unit,"(a)") "    (specify --help-options for a list of options)"
    write(iotk_error_unit,"(a)") "  where command is convert, dump, etc."
    write(iotk_error_unit,"(a)") "    (specify --help-commands for a list of commands)"
    write(iotk_error_unit,"(a)") "  where command-options-and-arguments depend on the specific command"
    write(iotk_error_unit,"(a)") "    (specify a command followed by --help for command-specific help)"
    write(iotk_error_unit,"(a)") "  Specify --help to receive this message"
  end if

  if(print_help_commands) then
    write(iotk_error_unit,"(a)") "IOTK commands are:"
    write(iotk_error_unit,"(a)") "  convert    to convert a file"
    write(iotk_error_unit,"(a)") "  dump       to dump a file"
    write(iotk_error_unit,"(a)") "  info       to obtain informations about how iotk was compiled"
    write(iotk_error_unit,"(a)") "  man        to print manual pages"
  end if

  if(print_help_options) then
    write(iotk_error_unit,"(a)") "IOTK options are:"
    write(iotk_error_unit,"(a)") "  --copyright        print copyright informations"
    write(iotk_error_unit,"(a)") "  --version          print version informations"
    write(iotk_error_unit,"(a)") "  --help             print a short, generic help"
    write(iotk_error_unit,"(a)") "  --help-options     print a list of options (this list)"
    write(iotk_error_unit,"(a)") "  --help-commands    print a list of commands"
    write(iotk_error_unit,"(a)") "  --set-linlen=N     to set the length of an output line"
    write(iotk_error_unit,"(a)") "  --set-indent=N     to set the number of spaces for an indent level"
    write(iotk_error_unit,"(a)") "  --set-maxindent=N  to set the maximum number of spaces when indenting"
  end if

  if(print_version) then
    write(*,"(a)") "Input/Output Tool Kit (IOTK) version: "//trim(iotk_version)
  end if

  if(print_copyright) then
    write(iotk_error_unit,"(a)") "Input/Output Tool Kit (IOTK)"
    write(iotk_error_unit,"(a)") "Copyright (C) 2004-2006 Giovanni Bussi"
    write(iotk_error_unit,"(a)") ""
    write(iotk_error_unit,"(a)") "This library is free software; you can redistribute it and/or"
    write(iotk_error_unit,"(a)") "modify it under the terms of the GNU Lesser General Public"
    write(iotk_error_unit,"(a)") "License as published by the Free Software Foundation; either"
    write(iotk_error_unit,"(a)") "version 2.1 of the License, or (at your option) any later version."
    write(iotk_error_unit,"(a)") ""
    write(iotk_error_unit,"(a)") "This library is distributed in the hope that it will be useful,"
    write(iotk_error_unit,"(a)") "but WITHOUT ANY WARRANTY; without even the implied warranty of"
    write(iotk_error_unit,"(a)") "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU"
    write(iotk_error_unit,"(a)") "Lesser General Public License for more details."
    write(iotk_error_unit,"(a)") ""
    write(iotk_error_unit,"(a)") "You should have received a copy of the GNU Lesser General Public"
    write(iotk_error_unit,"(a)") "License along with this library; if not, write to the Free Software"
    write(iotk_error_unit,"(a)") "Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA"
  end if

1 continue
  if(ierrl/=0) call iotk_error_handler(ierrl)

end subroutine iotk_tool_x

# 219 "iotk_tool.spp"
subroutine iotk_tool_convert_x(args,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_str_interf
  use iotk_misc_interf
  use iotk_files_interf
  implicit none
  character(len=*),           intent(in)  :: args(:)
  integer,          optional, intent(out) :: ierr
  integer :: iarg,ierrl,outfile_len
  character(len=iotk_fillenx) :: infile,outfile
  logical :: binary
  character(len=iotk_attlenx) :: attr
  character(len=iotk_taglenx) :: root
  integer :: maxsize
  logical :: autofmt
  infile=""
  outfile=""
  binary=.true.
  maxsize=-1
  ierrl = 0
  autofmt = .true.
  do iarg = 1 , size(args)
    if(iotk_strcomp(args(iarg)(1:1),"-")) then
      if(iotk_strcomp(args(iarg),"--help")) then
        write(iotk_error_unit,"(a)") "Usage: iotk convert [OPTIONS] infile outfile"
        write(iotk_error_unit,"(a)") "Options:"
        write(iotk_error_unit,"(a)") "  --mode=X  set the output file to be X, where X can be"
        write(iotk_error_unit,"(a)") "            'textual', 'binary' or 'auto'."
        write(iotk_error_unit,"(a)") "  -b        equivalent to --mode=binary"
        write(iotk_error_unit,"(a)") "  -t        equivalent to --mode=textual"
        write(iotk_error_unit,"(a)") "This command converts a iotk data file into another iotk data file."
        write(iotk_error_unit,"(a)") "The infile can be textual or binary, and its format is automatically detected."
        write(iotk_error_unit,"(a)") "The outfile can be textual or binary depending on the --mode option."
        write(iotk_error_unit,"(a)") "If the mode is 'auto', the decision is driven by outfile extension,"
        write(iotk_error_unit,"(a)") "i.e. a file matching *.txt of *.xml will be considered textual, otherwise binary"
        goto 1
      else if(iotk_strcomp(args(iarg),"-b") .or. iotk_strcomp(args(iarg),"--mode=binary")) then
        binary = .true.
        autofmt = .false.
      else if(iotk_strcomp(args(iarg),"-t") .or. iotk_strcomp(args(iarg),"--mode=textual")) then
        binary = .false.
        autofmt = .false.
      else if(iotk_strcomp(args(iarg),"--mode=auto")) then
        binary = .true.
        autofmt = .true.
      else
        call iotk_error_issue(ierrl,"iotk_tool_convert",__FILE__,__LINE__)
# 266 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 266 "iotk_tool.spp"
call iotk_error_msg(ierrl,'Unknown option')
        goto 1
      end if
    else
      if(infile=="") then
        call iotk_strcpy(infile,args(iarg),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool_convert",__FILE__,__LINE__)
# 273 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 273 "iotk_tool.spp"
call iotk_error_msg(ierrl,'File name too long')
          goto 1
        end if
      else if(outfile=="") then
        call iotk_strcpy(outfile,args(iarg),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool_convert",__FILE__,__LINE__)
# 279 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 279 "iotk_tool.spp"
call iotk_error_msg(ierrl,'File name too long')
          goto 1
        end if
      else
        call iotk_error_issue(ierrl,"iotk_tool_convert",__FILE__,__LINE__)
# 283 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 283 "iotk_tool.spp"
call iotk_error_msg(ierrl,'Three files. What do you mean?')
        goto 1
      end if
    end if
  end do
  if(outfile=="") then
    call iotk_error_issue(ierrl,"iotk_tool_convert",__FILE__,__LINE__)
# 289 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 289 "iotk_tool.spp"
call iotk_error_msg(ierrl,'Convert: bad usage')
    goto 1
  end if

  outfile_len = iotk_strlen(outfile)
  if(outfile_len>3) then
    select case(outfile(outfile_len-3:outfile_len))
    case(".xml")
      binary = .false.
    case(".txt")
      binary = .false.
    case default
      binary = .true.
    end select
  end if

  call iotk_open_read(60,infile,root=root,attr=attr)
  call iotk_open_write(61,outfile,binary=binary,root=root,attr=attr)
  call iotk_copy_tag(60,61,maxsize=-1)
  call iotk_close_write(61)
  call iotk_close_read(60)

1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_tool_convert_x


# 321 "iotk_tool.spp"
subroutine iotk_tool_dump_x(args,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_str_interf
  use iotk_misc_interf
  use iotk_files_interf
  implicit none
  character(len=*),           intent(in)  :: args(:)
  integer,          optional, intent(out) :: ierr
  integer :: iarg,ierrl
  character(len=iotk_fillenx) :: infile
  character(len=iotk_attlenx) :: attr
  character(len=iotk_taglenx) :: root
  integer :: maxsize
  infile=""
  maxsize=-1
  ierrl = 0
  do iarg = 1 , size(args)
    if(iotk_strcomp(args(iarg)(1:1),"-")) then
      if(iotk_strcomp(args(iarg),"--help")) then
        write(iotk_error_unit,"(a)") "Usage: iotk dump file"
        write(iotk_error_unit,"(a)") "This command dumps a iotk data file on standard out."
        write(iotk_error_unit,"(a)") "The file can be textual or binary, and its format is automatically detected."
        goto 1
      else
        call iotk_error_issue(ierrl,"iotk_tool_dump",__FILE__,__LINE__)
# 346 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 346 "iotk_tool.spp"
call iotk_error_msg(ierrl,'Unknown option')
        goto 1
      end if
    else
      if(infile=="") then
        call iotk_strcpy(infile,args(iarg),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool_dump",__FILE__,__LINE__)
# 353 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 353 "iotk_tool.spp"
call iotk_error_msg(ierrl,'File name too long')
          goto 1
        end if
      else
        call iotk_error_issue(ierrl,"iotk_tool_dump",__FILE__,__LINE__)
# 357 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 357 "iotk_tool.spp"
call iotk_error_msg(ierrl,'Two files. What do you mean?')
        goto 1
      end if
    end if
  end do

  call iotk_open_read(60, trim(infile),root=root,attr=attr)
  call iotk_open_write(iotk_output_unit,root=root,attr=attr)
  call iotk_copy_tag(60,iotk_output_unit,maxsize=-1)
  call iotk_close_write(iotk_output_unit)
  call iotk_close_read(60)

1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_tool_dump_x

subroutine iotk_tool_info_x(args,ierr)
  use iotk_base
  use iotk_misc_interf
  use iotk_xtox_interf
  use iotk_error_interf
  implicit none
  character(len=*),           intent(in)  :: args(:)
  integer,          optional, intent(out) :: ierr
  integer :: ierrl
  ierrl = 0
  write(*,"(a)") "IOTK (Input/Output Tool Kit) version: "//trim(iotk_version)
  write(*,"(a)") "Limits:"
  write(*,"(a)") "  maximum rank (soft limit): "//trim(iotk_itoa(iotk_maxrank))
  write(*,"(a)") "  maximum rank (hard limit): "//trim(iotk_itoa(iotk_maxrank_hard))
  write(*,"(a)") "Special kinds:"
  write(*,"(a)") "  headers in binary files are integer(kind="//trim(iotk_itoa(iotk_header_kind))//")"
  write(*,"(a)") "  default integers are integer(kind="//trim(iotk_itoa(iotk_integer_defkind))//")"
  write(*,"(a)") "  default logicals are logical(kind="//trim(iotk_itoa(iotk_logical_defkind))//")"
  write(*,"(a)") "  default characters are character(kind="//trim(iotk_itoa(iotk_character_defkind))//")"
  write(*,"(a)") "Kinds configured for i/o operations:"
# 398 "iotk_tool.spp"
#ifdef __IOTK_LOGICAL1
  write(*,"(a)") "  logical(kind="//trim(iotk_itoa(iotk_logical1))//")"
#endif
# 398 "iotk_tool.spp"
#ifdef __IOTK_LOGICAL2
  write(*,"(a)") "  logical(kind="//trim(iotk_itoa(iotk_logical2))//")"
#endif
# 398 "iotk_tool.spp"
#ifdef __IOTK_LOGICAL3
  write(*,"(a)") "  logical(kind="//trim(iotk_itoa(iotk_logical3))//")"
#endif
# 398 "iotk_tool.spp"
#ifdef __IOTK_LOGICAL4
  write(*,"(a)") "  logical(kind="//trim(iotk_itoa(iotk_logical4))//")"
#endif
# 403 "iotk_tool.spp"
#ifdef __IOTK_INTEGER1
  write(*,"(a)") "  integer(kind="//trim(iotk_itoa(iotk_integer1))//")"
#endif
# 403 "iotk_tool.spp"
#ifdef __IOTK_INTEGER2
  write(*,"(a)") "  integer(kind="//trim(iotk_itoa(iotk_integer2))//")"
#endif
# 403 "iotk_tool.spp"
#ifdef __IOTK_INTEGER3
  write(*,"(a)") "  integer(kind="//trim(iotk_itoa(iotk_integer3))//")"
#endif
# 403 "iotk_tool.spp"
#ifdef __IOTK_INTEGER4
  write(*,"(a)") "  integer(kind="//trim(iotk_itoa(iotk_integer4))//")"
#endif
# 408 "iotk_tool.spp"
#ifdef __IOTK_REAL1
  write(*,"(a)") "  real(kind="//trim(iotk_itoa(iotk_real1))//")"
#endif
# 408 "iotk_tool.spp"
#ifdef __IOTK_REAL2
  write(*,"(a)") "  real(kind="//trim(iotk_itoa(iotk_real2))//")"
#endif
# 408 "iotk_tool.spp"
#ifdef __IOTK_REAL3
  write(*,"(a)") "  real(kind="//trim(iotk_itoa(iotk_real3))//")"
#endif
# 408 "iotk_tool.spp"
#ifdef __IOTK_REAL4
  write(*,"(a)") "  real(kind="//trim(iotk_itoa(iotk_real4))//")"
#endif
# 413 "iotk_tool.spp"
#ifdef __IOTK_REAL1
  write(*,"(a)") "  complex(kind="//trim(iotk_itoa(iotk_real1))//")"
#endif
# 413 "iotk_tool.spp"
#ifdef __IOTK_REAL2
  write(*,"(a)") "  complex(kind="//trim(iotk_itoa(iotk_real2))//")"
#endif
# 413 "iotk_tool.spp"
#ifdef __IOTK_REAL3
  write(*,"(a)") "  complex(kind="//trim(iotk_itoa(iotk_real3))//")"
#endif
# 413 "iotk_tool.spp"
#ifdef __IOTK_REAL4
  write(*,"(a)") "  complex(kind="//trim(iotk_itoa(iotk_real4))//")"
#endif
# 417 "iotk_tool.spp"
  write(*,"(a)") "  character(kind="//trim(iotk_itoa(iotk_character1))//")"

1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_tool_info_x

subroutine iotk_tool_man_x(args,ierr)
  use iotk_base
  use iotk_misc_interf
  use iotk_xtox_interf
  use iotk_error_interf
  use iotk_str_interf
  implicit none
  character(len=*),           intent(in)  :: args(:)
  integer,          optional, intent(out) :: ierr
  character(len=iotk_linlenx) :: keyword
  integer :: ierrl,iarg
  logical :: printme,printlist

  ierrl = 0
  printme = .false.
  printlist = .false.
  keyword(1:1) = iotk_eos

  do iarg = 1 , size(args)
    if(iotk_strcomp(args(iarg)(1:1),"-")) then
      if(iotk_strcomp(args(iarg),"--help")) then
        write(iotk_error_unit,"(a)") "Usage: iotk man [keyword]"
        write(iotk_error_unit,"(a)") "This command prints on stdout the page of the built-in manual associated with the keyword."
        write(iotk_error_unit,"(a)") "If the keyword is not given a list of all the available keywords will be printed."
        goto 1
      else
        call iotk_error_issue(ierrl,"iotk_tool_dump",__FILE__,__LINE__)
# 453 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 453 "iotk_tool.spp"
call iotk_error_msg(ierrl,'Unknown option')
        goto 1
      end if
    else
      if(iotk_strcomp(keyword,"")) then
        call iotk_strcpy(keyword,args(iarg),ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_tool_dump",__FILE__,__LINE__)
# 460 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
          goto 1
        end if
      else
        call iotk_error_issue(ierrl,"iotk_tool_dump",__FILE__,__LINE__)
# 464 "iotk_tool.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 464 "iotk_tool.spp"
call iotk_error_msg(ierrl,'Two keywords. What do you mean?')
        goto 1
      end if
    end if
  end do

  if(iotk_strcomp(keyword,"")) then
    write(iotk_output_unit,"(a)") "List of available pages:"
    printlist = .true.
  end if
#ifndef __IOTK_WORKAROUND8
  if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" intro"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'intro')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: INTRODUCTION"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"The input/output tool kit (IOTK) is a FORTRAN90 library intended"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"to provide an easy access to tagged files formatted using some specific rule."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"In this context, a tagged file is a file containing tags and data."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Tagged files can be textual, in which case a XML-like format is used,"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"or binary, in which case a special format is used."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Note that IOTK is not an XML parser, but it can be used as a writer/parser"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"for a limited subset of the XML language."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"To use the IOTK library from a FORTRAN90 source, the user should"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"use the module 'iotk_module'."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"To minimize the possibility of name clashes, all public names exported"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"from this module has the "//'"'//&
# 475 "iotk_tool.spp"
"iotk_"//'"'//&
# 475 "iotk_tool.spp"
" prefix."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Communication between user and library is based on"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integers, characters and logicals of the default kind (notice that"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"these kinds can be changed using proper compiler options, so that"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"the actual kinds depend on how the library was compiled on your machine)."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"However, the library can handle formatted input/output for"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"all intrinsic datatypes, kinds and ranks if properly configured."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"This is obtained interfacing procedures which acts on all kinds,"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"types and (in almost all cases) ranks. Thus, a single generic"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"name has to be remembered for each subroutine, and the compiler will"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"link the correct one dependening on type, kind and rank of the arguments."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Backward API compatibility will be mantained (as long as it is possible)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"in future versions."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Backward file compatibility will be mantained (as long as it is possible) in"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"future versions."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"The library writes on files informations about the version of the library."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"It also writes informations about the version of the file format (file_version)."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"The later has to be older or equal to the format supported in the actual library."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" error_handling"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'error_handling')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: ERROR HANDLING"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"The way iotk handles error is sophisticated and allows for a trace back"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"of the error condition inside the library."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Every iotk routines which possibly leads to an error condition has an optional"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"intent(out) integer argument ierr. The returned value is conventionally"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"0 when the routine returns correctly, and different from 0 when the routines"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"raise an error. The value is effectively a handler for a more complex"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"object containing the error message. When an error is raised in a low-level"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk routine, a message is written on the error object. Any intermediate routine"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"can add other messages to the error object, at least the number of the line in"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"the source file. In this way, the error message contains a complete trace of"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"the error plus some additional information."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"At any point in the chain the messages can be exctracted from the error object."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"At some point in the chain the error is really handled, usually by writing the"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"message on the appropriate unit and aborting the execution."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Scanning routines (iotk_scan_*) have an optional logical argument "//'"'//&
# 475 "iotk_tool.spp"
"found"//'"'//&
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"which returns true or false. When scanning for data, also a "//'"'//&
# 475 "iotk_tool.spp"
"default"//'"'//&
# 475 "iotk_tool.spp"
" argument"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"can be used. If one of these two argument is present, the searched"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"object is considered as an optional object. Otherwise, it is considered as a needed object."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If the ierr optional argument is absent, the error handling is leaved to the iotk library."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"In this case, if a needed object is not present, the library handles the error with a"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"forced stop."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If the ierr optional argument is present, it returns an error code."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"ierr = 0 means that no error has occurred"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"ierr > 0 means that an error has occurred probably related to file corruption"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"ierr < 0 means that the item that was searched for has not been found"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"(it is possible only for scanning routines and only if the"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"found and the default keywords are both missing, i.e. only for no-optional objetcs)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"In scanning routines, if the argument "//'"'//&
# 475 "iotk_tool.spp"
"found"//'"'//&
# 475 "iotk_tool.spp"
" is present it returns .true."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"if the item has been found, .false. otherwise."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If a library routine returns an ierr /= 0 it is STRONGLY RECOMMENDED to"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"clear that error with "//'"'//&
# 475 "iotk_tool.spp"
"call iotk_error_clear(ierr)"//'"'//&
# 475 "iotk_tool.spp"
" before proceeding."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Thus, the final recipe is:"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"* if you want to handle errors, always use the 'ierr' optional argument."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"looking at the sign, you will discern between lacking data and file corruption."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"with iotk_error_print you can obtain a description of the error."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"* if you want to leave the error handling to the library, don't use"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"the 'ierr' optional argument."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"- if the object you are searching is optional, use 'found' or 'default' optional arguments."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"- if the object you are searching is non-optional, don't use 'found' nor 'default' optional arguments."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" binary_and_textual_files"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'binary_and_textual_files')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: BINARY AND TEXTUAL FILES"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Units can be opened on textual or binary files."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"The word 'binary' is used instead of the fortran 'unformatted' since"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"using this libray also binary files have a degree of formattation."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"After a unit has been opened, the library automatically detects"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"its format through an INQUIRE and acts consequently."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Note that the iotk routines check for necessary properties of an opened unit"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"access="//'"'//&
# 475 "iotk_tool.spp"
"sequential"//'"'//&
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"blank ="//'"'//&
# 475 "iotk_tool.spp"
"null"//'"'//&
# 475 "iotk_tool.spp"
" (only textual i/o)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"pad   ="//'"'//&
# 475 "iotk_tool.spp"
"yes"//'"'//&
# 475 "iotk_tool.spp"
"  (only textual i/o)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Moreover, a textual or binary unit can be designed as raw."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"In that case, no tags are placed on the file and everything"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"has to be read and written in the same order."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"This feature is provided for compatibility reasons but it should be"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"used as few as possible."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" optional_arguments"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'optional_arguments')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: OPTIONAL ARGUMENTS"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Most iotk routines accept optional arguments."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"The calling routine will not compile if the names of the"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"arguments are not indicated.  For instance, use"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"call iotk_scan_dat(10,"//'"'//&
# 475 "iotk_tool.spp"
"pippo"//'"'//&
# 475 "iotk_tool.spp"
",aa(:),ierr=ii)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"and NOT:"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"call iotk_scan_dat(10,"//'"'//&
# 475 "iotk_tool.spp"
"pippo"//'"'//&
# 475 "iotk_tool.spp"
",aa(:),ii)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"The only exeption is the attr argument, for which the name can be"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"omitted if it is placed as the first of the optional arguments."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"In any case, it is better to always explicitly label optional arguments."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" basic_writing_routines iotk_write_begin iotk_write_end iotk_write_empty iotk_write_pi iotk_write_comment"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'basic_writing_routines')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_write_begin')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_write_end')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_write_empty')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_write_pi')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_write_comment')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: BASIC WRITING ROUTINES"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_write_begin  (unit,name[,attr][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_write_end    (unit,name[,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_write_empty  (unit,name[,attr][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_write_pi     (unit,name[,attr][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_write_comment(unit,text[,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(in) :: unit"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in) :: name"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in) :: text"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in) :: attr"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(out):: ierr ! see error_handling page"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"These routines write a tag named 'name' on fortran unit 'unit'."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"The type of the tag is determined from the name of the routine:"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_write_begin   => <name attr>"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_write_end     => </name>"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_write_empty   => <name attr/>"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_write_pi      => <?name attr?>"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_write_comment => <!--text-->"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"An optional attribute string can be supplied in 'attr'"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"In end tags, no attribute is allowed."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"To build the attribute string, use iotk_write_attr."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"DON'T TRY TO MANIPULATE THE ATTRIBUTE STRING DIRECTLY!"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" basic_scanning_routines iotk_scan_begin iotk_scan_end iotk_scan_empty iotk_scan_pi"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'basic_scanning_routines')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_scan_begin')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_scan_end')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_scan_empty')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_scan_pi')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: BASIC SCANNING ROUTINES"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_scan_begin(unit,name[,attr][,found][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_scan_end  (unit,name[,found][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_scan_empty(unit,name[,attr][,found][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_scan_pi   (unit,name[,attr][,found][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(in) :: unit"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in) :: name  ! len less or equal iotk_namlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(out):: attr  ! len possibily equal iotk_attlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(out):: found ! see error_handling page"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(out):: ierr  ! see error_handling page"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"These routines scan for a tag named 'name' on fortran unit 'unit'."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"The type of the tag is determined from the name of the routine:"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_scan_begin => <name attr>"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_scan_end   => </name>"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_scan_empty => <name attr/>"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_scan_pi    => <?name attr?>"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"These routines (except for iotk_scan_end) also fills the"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"attr string, which can be subsequently decoded with iotk_scan_attr."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"DON'T TRY TO MANIPULATE THE ATTRIBUTE STRING DIRECTLY!"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" writing_attributes iotk_write_attr"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'writing_attributes')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_write_attr')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: WRITING ATTRIBUTES"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_write_attr (attr,name,val[,first][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(out):: attr  ! len less or equal iotk_namlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in) :: name  ! len less or equal iotk_attlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"TYPE(KIND),       intent(in) :: val   ! any type, any kind, any rank [but only scalars for character]"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(in) :: first"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(out):: ierr"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"This routine adds one attribute to the 'attr' string."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"To clean the string (for the first attribute) use first=.true."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Example:"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"call iotk_write_attr(attr,"//'"'//&
# 475 "iotk_tool.spp"
"pippo"//'"'//&
# 475 "iotk_tool.spp"
",1,first=.true.)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"call iotk_write_attr(attr,"//'"'//&
# 475 "iotk_tool.spp"
"paperino"//'"'//&
# 475 "iotk_tool.spp"
",2)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"call iotk_write_attr(attr,"//'"'//&
# 475 "iotk_tool.spp"
"pluto"//'"'//&
# 475 "iotk_tool.spp"
",3)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"This is equivalent to attr="//'"'//&
# 475 "iotk_tool.spp"
""//'"'//&
# 475 "iotk_tool.spp"
" before the call, but more efficient."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"The attribute is added in the form name="//'"'//&
# 475 "iotk_tool.spp"
"value"//'"'//&
# 475 "iotk_tool.spp"
","
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"where "//'"'//&
# 475 "iotk_tool.spp"
"value"//'"'//&
# 475 "iotk_tool.spp"
" is a string containing a textual representation"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"of the val variable."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If one of <>&"//'"'//&
# 475 "iotk_tool.spp"
"' appears in val, it is automatically escaped."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" scanning_attributes iotk_scan_attr"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'scanning_attributes')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_scan_attr')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: SCANNING ATTRIBUTES"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_scan_attr  (attr,name,val[,found][,default][,eos][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in) :: attr    ! len possibily equal iotk_attlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in) :: name    ! len less or equal iotk_namlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"TYPE(KIND),       intent(out):: val     ! any type, any kind, any rank [but only scalars for character]"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(out):: found   ! see error_handling page"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"TYPE(KIND),       intent(in) :: default ! same type, kind and rank as val"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(in) :: eos"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(out):: ierr    ! see error_handling page"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"This routine scans for one attribute named 'name' from the 'attr' string."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If the attribute is found, it is read to variable 'val'."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If it is not found and default is present, default is copied onto val."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If TYPE is character and eos is present and true,"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"an end-of-string terminator will be attached at the end of the read string,"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"and the following bytes will not be touched. This is faster, but requires"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"the user to take care directly of the end-of-string. Thus, it is discouraged."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"The attribute can be delimited with "//'"'//&
# 475 "iotk_tool.spp"
""//'"'//&
# 475 "iotk_tool.spp"
" or with ''"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" writing_data iotk_write_dat"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'writing_data')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_write_dat')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: WRITING DATA"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_write_dat  (unit,name,dat[,fmt][,columns][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(in) :: unit"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in) :: name    ! len less or equal iotk_namlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"TYPE(KIND),       intent(in) :: dat     ! any type, any kind, any rank"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in) :: fmt"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(in) :: columns"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(out):: ierr    ! see error_handling page"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"This routines write a data object, that is a self-described"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"object containg fortran data."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"A single data object has the following form"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"<name type="//'"'//&
# 475 "iotk_tool.spp"
"TYPE"//'"'//&
# 475 "iotk_tool.spp"
" kind="//'"'//&
# 475 "iotk_tool.spp"
"KIND"//'"'//&
# 475 "iotk_tool.spp"
" size="//'"'//&
# 475 "iotk_tool.spp"
"SIZE"//'"'//&
# 475 "iotk_tool.spp"
" columns="//'"'//&
# 475 "iotk_tool.spp"
"COLUMNS"//'"'//&
# 475 "iotk_tool.spp"
" len="//'"'//&
# 475 "iotk_tool.spp"
"LEN"//'"'//&
# 475 "iotk_tool.spp"
" fmt="//'"'//&
# 475 "iotk_tool.spp"
"FMT"//'"'//&
# 475 "iotk_tool.spp"
">"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
".. DATA ..."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"</name>"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"where"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"TYPE    is the intrinsic type (logical,integer,real,complex or character),"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"KIND    is the data kind (stored in binary files only)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"SIZE    is the array size (shape informations are not stored)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"COLUMNS is the number of data per line"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"LEN     is the string length"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"FMT     is a fortran format string used to write data"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If the optional 'fmt' is not passed, default format ('columns' element per line)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"is used and the fmt attribute is not written. Otherwise, the string"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"fmt is used as a FORTRAN format specifierfor the write statement. In this"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"case it is also written on the file (and used for reading the data back)."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"fmt="//'"'//&
# 475 "iotk_tool.spp"
"*"//'"'//&
# 475 "iotk_tool.spp"
" can be used and correspond to the "//'"'//&
# 475 "iotk_tool.spp"
"write(unit,*)"//'"'//&
# 475 "iotk_tool.spp"
" statement."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If the optional 'columns' is not passed, it is assumed to be 1 and"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"the columns attribute is not written. Note that this attribute is completely"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"ininfluent when reading."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"columns and fmt arguments are incompatible."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"For complex data, one element is a couple of comma separated real numbers."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If one of <>& appears in dat, it is escaped."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" scanning_data iotk_scan_dat"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'scanning_data')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_scan_dat')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: SCANNING DATA"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_scan_dat  (unit,name,dat[,found][,default][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(in) :: unit"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in) :: name    ! len less or equal iotk_namlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"TYPE(KIND),       intent(out):: dat     ! any type, any kind, any rank"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(out):: found   ! see error_handling page"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"TYPE(KIND),       intent(in) :: default ! same type, kind and rank as dat"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(out):: ierr    ! see error_handling page"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"A data object written with iotk_write_dat is read."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If it is not found and default is present, default is copied onto dat."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If a keyword is absent in the file, the value is deduced from the"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"dat argument and no check is performed. This allows to write"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"rapidly by hand data objects. For instance"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"<datum> 1.0 </datum>"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"can be read with"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"real :: val"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"call iotk_scan_dat(unit,"//'"'//&
# 475 "iotk_tool.spp"
"datum"//'"'//&
# 475 "iotk_tool.spp"
",val)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If fmt is not present on file, the default format is used."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Types and sizes are checked."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Different kinds (for binary i/o) are automatically converted."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Length (for characters) are not checked. If strings on files"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"are longer then len(dat), only the first characters are read; if strings"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"on files are shorter, dat is padded with blanks."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" opening_files iotk_open_write iotk_open_read"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'opening_files')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_open_write')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_open_read')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: OPENING FILES"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_open_write(unit[,file][,attr][,binary][,raw][,new][,root][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(in)  :: unit"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in)  :: file"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in)  :: attr"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(in)  :: binary"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(in)  :: new"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(in)  :: raw"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in)  :: root   ! len less or equal iotk_namlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(out) :: ierr   ! see error_handling page"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If file is present, this routines opens file 'file' on"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"unit 'unit' with the proper options."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If binary is present and true, the file is binary."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If new is present and true, the file must not exist already."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If raw is present and true, the file is considered as a raw data file"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"(use of raw data files is discouraged)."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If file is not present, unit is assumed to be already connected."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If root is present, it is used as the name of the root begin/end tag."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If it is absent, the default "//'"'//&
# 475 "iotk_tool.spp"
"Root"//'"'//&
# 475 "iotk_tool.spp"
" is used."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"An optional attribute string can be supplied in 'attr', and will be used"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"as an attribute list for the begin root tag."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Also informations about iotk version and binary format are written as"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"pi informations."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_open_read(unit[,file][,attr][,binary][,raw][,root][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(in)  :: unit"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in)  :: file"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(out) :: attr"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(in)  :: binary"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(in)  :: raw"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(out) :: root   ! len possibly equal iotk_namlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(out) :: ierr   ! see error_handling page"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If file is present, this routines opens file 'file' on"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"unit 'unit' with the proper options."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If binary is present and true, the file is binary."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If raw is present and true, the file is considered as a raw data file"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"(use of raw data files is discouraged)."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If file is not present, unit is assumed to be already connected."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If root is present, the name of root in file is read onto that variable."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If attr is present, the attributes of the root tag are read onto that variable,"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"which can be subsequently decoded with iotk_scan_attr."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"DON'T TRY TO MANIPULATE THE ATTRIBUTE STRING DIRECTLY!"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" closing_files iotk_close_write iotk_close_read"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'closing_files')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_close_write')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_close_read')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: CLOSING FILES"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_close_write(unit[,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_close_read(unit[,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,      intent(in)  :: unit"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,      intent(out) :: ierr ! see error_handling page"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"This routines close a file opened with iotk_open_*"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Note that if the units were already connected before iotk_open_*, they"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"are left connected here."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" multiple_files iotk_link"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'multiple_files')) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'iotk_link')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: MULTIPLE FILES"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"When reading, if a begin tag with an attribute iotk_link="//'"'//&
# 475 "iotk_tool.spp"
"FILENAME"//'"'//&
# 475 "iotk_tool.spp"
" is found,"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"file FILENAME is mounted in its place"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"If FILENAME begins with a "//'"'//&
# 475 "iotk_tool.spp"
"/"//'"'//&
# 475 "iotk_tool.spp"
", the path is absolute, otherwise it is relative"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"to the original file."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Note that the mounting is completely transparent for users, which can access"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"the new file using the old unit. However, if the user wants to access"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"directly the new file, iotk_physical_unit should be used."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"When writing, the user can switch a logical unit to a different file using"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"the following routine"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_link(unit,name,file,dummy[,binary][,raw][,create][,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(in)  :: unit"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in)  :: name"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*), intent(in)  :: file"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(in)  :: binary"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(in)  :: raw"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"logical,          intent(in)  :: create"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer,          intent(out) :: ierr"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"name is the name of the tag which represents the link."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"file is the name of the new file"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"if binary is present and true, the new file will be binary"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"if raw is present and true, the new file will be raw"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"if create is present and true, the new file is actually created"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"and the next write statement will act on this new file automatically."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Otherwise, only the symbolic link is created."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printlist) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
" utilities"
# 475 "iotk_tool.spp"
printme=.false.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,"all")) printme=.true.
# 475 "iotk_tool.spp"
if(iotk_strcomp(keyword,'utilities')) printme=.true.
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"IOTK: OTHER UTILITIES"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Here a number of additional routines/parameters available"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"from the iotk_module is listed"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*) :: iotk_index (index)"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer, intent(in) :: index ! scalar or rank 1"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Returns a string representing the index in an array."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"Example: index = (/1,2,3/) => iotk_index = "//'"'//&
# 475 "iotk_tool.spp"
".1.2.3"//'"'//&
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"The correct way for writing an array of derived types is"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"to build the names as follows"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"! ONE-DIMENSIONAL ARRAY"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"do i = 1 , n"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"call iotk_write_begin(unit,"//'"'//&
# 475 "iotk_tool.spp"
"dummy"//'"'//&
# 475 "iotk_tool.spp"
"//iotk_index(i))"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"! WRITE THE OBJECT HERE"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"call iotk_write_end  (unit,"//'"'//&
# 475 "iotk_tool.spp"
"dummy"//'"'//&
# 475 "iotk_tool.spp"
"//iotk_index(i))"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"end do"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"do i = 1 , n"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"do j = 1 , m"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"! NOTE THE ORDER OF INDEXES, THE FASTER IS THE LAST"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"call iotk_write_begin(unit,"//'"'//&
# 475 "iotk_tool.spp"
"dummy"//'"'//&
# 475 "iotk_tool.spp"
"//iotk_index((/i,j/)))"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"! WRITE THE OBJECT HERE"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"call iotk_write_end  (unit,"//'"'//&
# 475 "iotk_tool.spp"
"dummy"//'"'//&
# 475 "iotk_tool.spp"
"//iotk_index((/i,j/)))"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"end do"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"end do"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"iotk_free_unit(unit[,ierr])"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer, intent(out) :: unit"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer, intent(out) :: ierr"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"This routine returns the number of a free FORTRAN unit."
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character(len=*) :: iotk_version"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"version string of iotk"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character :: iotk_newline"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"newline sequence"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"character :: iotk_eos"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"end-of-string character"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer :: iotk_taglenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"max length of a tag"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer :: iotk_namlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"max length of a tag or attribute name"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer :: iotk_attlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"max length of the attribute string"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer :: iotk_vallenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"max length of the value of an attribute"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer :: iotk_linlenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"max length of a line in textual files"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer :: iotk_fillenx"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"max length of a file name"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer :: iotk_header_kind"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
"integer kind of headers in binary files"
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
# 475 "iotk_tool.spp"
if(printme) write(iotk_output_unit,"(a)") &
# 475 "iotk_tool.spp"
""
#endif
  1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_tool_man_x


