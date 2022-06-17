module bin2xml
    use iotk_module
    contains 
    subroutine write_bin2xml(source, dest)
        character(LEN=256) :: source
        character(LEN=256) :: dest
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
    character(LEN=256) :: source_filename
    character(LEN=256) :: dest_filename

    if(command_argument_count().ne.2) then 
        call errore('bin2xml_real_space_density', 'Wrong value number of input arguments', 1 )
    end if

    call get_command_argument(1, source_filename)
    call get_command_argument(2, dest_filename)


    ! write(*,*) source_filename
    ! write(*,*) dest_filename


    ! source_filename     = '/scratch/yshubert/All_Water/path-integral-nqe_0/snapshot_1/kc_70.save/sic_potential.occ.1.dat'
    ! dest_filename       = '/home/yshubert/Documents/master_project/final_ML_complete/extract_descriptor/All_Water/path-integral-nqe_0/snapshot_1/orbital.occ.1.xml'

    call write_bin2xml(source_filename, dest_filename)

end program 