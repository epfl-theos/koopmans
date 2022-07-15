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
    character(len=256) :: source_filename
    character(len=256) :: dest_filename
    character(len=5)   :: orbital_identifier, orbital_identifier_dat
    character(len=1)   :: spin_identifier
    character(len=3)   :: occpuation_identifier 

    integer              :: num_bands
    integer, allocatable :: bands(:,:)

    call get_command_argument(1, source_dir)
    call get_command_argument(2, dest_dir)
    call get_command_argument(3, nbsp_occ_char)

    read(nbsp_occ_char,'(i)') nbsp_occ
    
    
    open(unit=1, file=TRIM(dest_dir)//'/bands_to_solve.txt', status='old', action='read')
    read(1, *) num_bands
    allocate(bands(num_bands, 3))
    do i=1,num_bands
        read(1, *) bands(i,:)
    enddo
    close (unit=1)

    if(command_argument_count().ne.3) then 
        call errore('bin2xml_real_space_density', 'Wrong value number of input arguments', 1 )
    end if

    ! First write total density to XML
    source_filename =  TRIM(source_dir)//'/charge-density.dat'
    dest_filename   =  TRIM(dest_dir)//'/charge-density.xml'
    call write_bin2xml(source_filename, dest_filename)

    ! Then the orbital densities
    do i = 1, num_bands
        write(spin_identifier, "(I1.1)")       bands(i,3)
        write(orbital_identifier, "(I5.5)")    bands(i,1)

        if (bands(i,2)==1) then 
            occpuation_identifier = 'occ'
        else 
            occpuation_identifier = 'emp'
        end if

        source_filename =  TRIM(source_dir)//'/real_space_orb_density.' // TRIM(occpuation_identifier) // '.' // TRIM(spin_identifier) // '.' //TRIM(orbital_identifier)//'.dat'
        dest_filename   =  TRIM(dest_dir)  //'/orbital.'                // TRIM(occpuation_identifier) // '.' // TRIM(spin_identifier) // '.' //TRIM(orbital_identifier)//'.xml'


        !     source_filename =  TRIM(source_dir)//'/real_space_orb_density.occ.' // TRIM(spin_identifier) // '.' //TRIM(orbital_identifier)//'.dat'
        !     dest_filename   =  TRIM(dest_dir)//'/orbital.occ.' // TRIM(spin_identifier) // '.' //TRIM(orbital_identifier)//'.xml'
        ! else 
        !     source_filename =  TRIM(source_dir)//'/real_space_orb_density.emp.' // TRIM(spin_identifier) // '.' //TRIM(orbital_identifier_dat)//'.dat'
        !     dest_filename   =  TRIM(dest_dir)//'/orbital.emp.' // TRIM(spin_identifier) // '.' //TRIM(orbital_identifier)//'.xml'
        ! end if 
        ! write(*,*) "Yannick Debug: orbital identifier"
        ! write(*,*) orbital_identifier

        call write_bin2xml(source_filename, dest_filename)
    end do 

end program 