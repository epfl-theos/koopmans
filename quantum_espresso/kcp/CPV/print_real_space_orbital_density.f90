module print_real_space_orbital_density
    integer, dimension(2)      :: spin_dependent_idx=(/1, 1/)
    contains
    subroutine print_orbr(bec, nbsp_filling, ispin, lgam, is_empty, c)
        ! simple routine to print real space orbital densities to output files
        use grid_dimensions,          only: nnrx
        use gvecw,                    only: ngw
        use twin_types
        use cp_interfaces,            only: nksic_get_orbitalrho ! computes the real space density
        use io_pot_sic_xml,           only: write_pot_sic        ! writes a scalar field to a binary file
        use nksic,                    only: orb_rhor             ! array to store the real space density in 
        
        type(twin_matrix), intent(in) :: bec 
        logical, intent(in)           :: lgam 
        logical, intent(in)           :: is_empty            ! occupation
        integer, intent(in)           :: nbsp_filling        ! number of occupied/empty orbitals
        integer, intent(in)           :: ispin(nbsp_filling) ! array containing the information of the spin channels
        complex(dp), intent(in)       :: c(ngw,nbsp_filling) ! coefficient vector

        ! local variables
        integer             :: j
        integer             :: jj
        integer             :: i
        character(len=256)  :: filename
        character(len=2)    :: spin_identifier
        character(len=5)    :: orbital_identifier

        
    
        ! looping of each orbital
        do j=1,nbsp_filling,2
            
            ! computing the real space density from the coefficient vectors. This function is implemented in a way, that it
            ! always takes two adjacent coeffient vectors and computes the real space density of both
            call nksic_get_orbitalrho( ngw, nnrx, bec, ispin, nbsp_filling, c(:,j), c(:,j+1), orb_rhor, j, j+1, lgam) 
            
            ! write both densities to output file
            do jj = 1, 2

                ! determine the name of the output file. It should contain information about the index of the orbital, its spin and its occupation
                i=j+jj-1
                write(spin_identifier, "(I1.1)")       (ispin(i)-1) ! spin 
                write(orbital_identifier, "(I5.5)")    spin_dependent_idx(ispin(i)) ! orbital

                if(is_empty) then    ! occupation
                    filename = 'emp.' // TRIM(spin_identifier) // '.' // TRIM(orbital_identifier) 
                else                    
                    filename = 'occ.' // TRIM(spin_identifier) // '.' // TRIM(orbital_identifier)
                end if

                ! write the density to an output file
                call write_pot_sic ( orb_rhor(:, jj), filename, field_specifier='real_space_orb_density')

                spin_dependent_idx(ispin(i)) = spin_dependent_idx(ispin(i)) + 1
            end do  
        end do
    end subroutine
end module