module bin2xml
    use iotk_module
    contains 
    subroutine write_bin2xml(source, dest)
        character(len=256) :: source
        character(len=256) :: dest
        integer            :: maxsize
        integer            :: ierr
        integer            :: iunout
        maxsize = -1
        call iotk_open_read(60,  source, binary= .True.)
        call iotk_open_write(61, dest,   binary= .False.)
        call iotk_copy_tag(60, 61, maxsize=maxsize)
        call iotk_close_write(61)
        call iotk_close_read(60)
    end subroutine write_bin2xml
end module bin2xml

program bin2xml_real_space_density
    use bin2xml
    character(len=256) :: source_dir
    character(len=256) :: dest_dir
    character(len=5)   :: nbsp_occ_char
    character(len=5)   :: nbsp_emp_char

    integer            :: nbsp_occ
    integer            :: nbsp_emp
    character(LEN=256) :: source_filename
    character(LEN=256) :: dest_filename
    character(LEN=5)   :: file_number
    character(LEN=5)   :: file_number_emp
    

    if(command_argument_count().ne.4) then 
        call errore('bin2xml_real_space_density', 'Wrong value number of input arguments', 1 )
    end if

    call get_command_argument(1, source_dir)
    call get_command_argument(2, dest_dir)
    call get_command_argument(3, nbsp_occ_char)
    call get_command_argument(4, nbsp_emp_char)

    read(nbsp_occ_char,'(i)') nbsp_occ
    read(nbsp_emp_char,'(i)') nbsp_emp


    ! First write total density to XML
    source_filename =  TRIM(source_dir)//'/charge-density.dat'
    dest_filename   =  TRIM(dest_dir)//'/charge-density.xml'
    call write_bin2xml(source_filename, dest_filename)

    ! Then write all orbital densities to XML
    do i = 1, nbsp_occ
        write(file_number, "(I5.5)") i
        source_filename =  TRIM(source_dir)//'/sic_potential.occ.'//TRIM(file_number)//'.dat'
        dest_filename   =  TRIM(dest_dir)//'/orbital.occ.'//TRIM(file_number)//'.xml'
        call write_bin2xml(source_filename, dest_filename)
    end do

    ! Then write all emp densities
    do i = 1, nbsp_emp
        write(file_number, "(I5.5)") (i + 2*nbsp_occ)
        write(file_number_emp, "(I5.5)") (i + nbsp_occ)
        source_filename =  TRIM(source_dir)//'/sic_potential.emp.'//TRIM(file_number)//'.dat'
        dest_filename   =  TRIM(dest_dir)//'/orbital.emp.'//TRIM(file_number_emp)//'.xml'
        call write_bin2xml(source_filename, dest_filename)
    end do
end program 