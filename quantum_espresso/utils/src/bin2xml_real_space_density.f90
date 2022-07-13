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
    character(len=1)   :: nspin_char

    integer            :: nbsp_occ
    integer            :: nbsp_emp
    integer            :: nspin
    integer            :: spin
    character(LEN=256) :: source_filename
    character(LEN=256) :: dest_filename
    character(LEN=5)   :: orbital_identifier, orbital_identifier_xml
    character(LEN=1)   :: spin_identifier

    integer              :: num_bands
    integer, allocatable :: bands(:)

    call get_command_argument(1, source_dir)
    call get_command_argument(2, dest_dir)
    call get_command_argument(3, nbsp_occ_char)
    call get_command_argument(4, nbsp_emp_char)
    call get_command_argument(5, nspin_char)


    read(nbsp_occ_char,'(i)') nbsp_occ
    read(nbsp_emp_char,'(i)') nbsp_emp
    read(nspin_char,'(i)')    nspin
    

    open(unit=1, file=TRIM(dest_dir)//'/bands_to_solve.txt', status='old', action='read')
    read(1, *) num_bands
    allocate(bands(num_bands))
    write(*,*) "Beginn read"
    
    do i=1,num_bands
        read(1, *) bands(i)
    enddo
    CLOSE (UNIT=1)
    write(*,*) "Beginn write"

    do i=1,num_bands
        write(*,*), bands(i)
    enddo

    

    if(command_argument_count().ne.5) then 
        call errore('bin2xml_real_space_density', 'Wrong value number of input arguments', 1 )
    end if

    


    ! First write total density to XML
    source_filename =  TRIM(source_dir)//'/charge-density.dat'
    dest_filename   =  TRIM(dest_dir)//'/charge-density.xml'
    call write_bin2xml(source_filename, dest_filename)

    do spin = 0, nspin-1
        ! Then write all orbital densities to XML
        do i = 1, nbsp_occ
            write(spin_identifier, "(I1.1)")       spin
            write(orbital_identifier, "(I5.5)")    i

            source_filename =  TRIM(source_dir)//'/sic_potential.occ.' // TRIM(spin_identifier) // '.' //TRIM(orbital_identifier)//'.dat'
            dest_filename   =  TRIM(dest_dir)//'/orbital.occ.' // TRIM(spin_identifier) // '.' //TRIM(orbital_identifier)//'.xml'
            call write_bin2xml(source_filename, dest_filename)
        end do

        ! Then write all emp densities
        do i = 1, nbsp_emp
            write(spin_identifier, "(I1.1)")       spin
            write(orbital_identifier, "(I5.5)")    i
            write(orbital_identifier_xml, "(I5.5)")    i+nbsp_occ

            source_filename =  TRIM(source_dir)//'/sic_potential.emp.' // TRIM(spin_identifier) // '.' //TRIM(orbital_identifier)//'.dat'
            dest_filename   =  TRIM(dest_dir)//'/orbital.emp.' // TRIM(spin_identifier) // '.' //TRIM(orbital_identifier_xml)//'.xml'
            call write_bin2xml(source_filename, dest_filename)
        end do
    end do
end program 