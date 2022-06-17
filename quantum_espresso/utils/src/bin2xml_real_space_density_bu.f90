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

    if(command_argument_count().ne.2) then 
        call errore('bin2xml_real_space_density', 'Wrong value number of input arguments', 1 )
    end if

    call get_command_argument(1, source_filename)
    call get_command_argument(2, dest_filename)


    if (compute_charge_density==.True.) then
        source_filename =  TRIM(source_dir)//'charge-density.dat'
        dest_filename   =  TRIM(dest_dir)//'charge-density.xml'
        call write_bin2xml(source_filename, dest_filename)
    end if 

    do i = 1, nbsp_occ_start+nbsp_occ
        if(i<10) then
            write(file_number, "(I1)") i
        else if(i<100) then
            write(file_number, "(I2)") i
        else if
            write(file_number, "(I3)") i
        end if(i<1000) then
            write(file_number, "(I4)") i
        else if(i<10000)
            write(file_number, "(I5)") i
        else 
            call errore('bin2xml_real_space_density', 'Too many orbitals for this program', 1 )
        source_filename =  TRIM(source_dir)//'sic_potential.occ.'//TRIM(file_number)//'.dat'
        dest_filename   =  TRIM(dest_dir)//'orbital.occ.'//TRIM(file_number)//'.xml'
        call copy_file_yannick(source_filename, dest_filename)
    end do

    do i = nbsp_empty_start, nbsp_empty
        if(i<10) then
            write(file_number, "(I1)") i
        else if(i<100) then
            write(file_number, "(I2)") i
        else
            write(file_number, "(I3)") i
        end if
        source_filename =  TRIM(source_dir)//'sic_potential.empty.'//TRIM(file_number)//'.dat'
        dest_filename   =  TRIM(dest_dir)//'orbital.empty.'//TRIM(file_number)//'.xml'
        call copy_file_yannick(source_filename, dest_filename)
    end do
end program 