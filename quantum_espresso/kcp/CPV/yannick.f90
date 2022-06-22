module yannick_print_orbr
    contains
    subroutine print_orbr(bec, nbsp_filling, ispin, nbspx, lgam, is_empty)
        use grid_dimensions,          only: nnrx
        use nksic,                    only: orb_rhor
        use io_pot_sic_xml,           only: write_pot_sic
        use cp_interfaces,            only: nksic_get_orbitalrho
        use gvecw,                    only: ngw
        use wavefunctions_module,     only: c0
        use electrons_base,           only: nbsp      ! Yannick Debug
        use twin_types
    
        type(twin_matrix), intent(in) :: bec 
        integer, intent(in)           :: nbspx, nbsp_filling, ispin(nbspx)
        logical, intent(in)           :: lgam
        logical, intent(in)           :: is_empty

        ! local variables
        integer             :: j
        integer             :: jj
        integer             :: i
        character(len=256)  :: filename
        character(len=2)    :: spin_identifier
        character(len=5)    :: orbital_identifier
    
        do j=1,nbsp_filling,2
            call nksic_get_orbitalrho( ngw, nnrx, bec, ispin, nbsp_filling, &
                        c0(:,j), c0(:,j+1), orb_rhor, j, j+1, lgam) 
            
            do jj = 1, 2
                i=j+jj-1
                
                write(spin_identifier, "(I1.1)")       jj-1
                write(orbital_identifier, "(I5.5)")    (j+1)/2

                if(is_empty) then      
                    filename = 'emp.' // TRIM(spin_identifier) // '.' // TRIM(orbital_identifier)
                else
                    filename = 'occ.' // TRIM(spin_identifier) // '.' // TRIM(orbital_identifier)
                end if
                call write_pot_sic ( orb_rhor(:, jj), filename)
            end do  
        end do
    end subroutine
end module