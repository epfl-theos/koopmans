module yannick_print_orbr
    contains
    subroutine print_orbr(bec, nbsp_filling, ispin, nbspx, lgam, is_empty)
        use grid_dimensions,          only: nnrx
        use nksic,                    only: orb_rhor
        use io_pot_sic_xml,           only: write_pot_sic
        use cp_interfaces,            only: nksic_get_orbitalrho
        use gvecw,                    only: ngw
        use wavefunctions_module,     only: c0
        use electrons_base,           only:  nbsp      ! Yannick Debug
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
        character(len=256)  :: filename_complete
    
        do j=1,nbsp_filling,2
            call nksic_get_orbitalrho( ngw, nnrx, bec, ispin, nbsp_filling, &
                        c0(:,j), c0(:,j+1), orb_rhor, j, j+1, lgam) 
            do jj = 1, 2
                i=j+jj-1

                if(is_empty) then 
                    write(filename, "(I5.5)") i+nbsp
                    filename_complete = 'emp.' // TRIM(filename)
                else
                    write(filename, "(I5.5)") i
                    filename_complete = 'occ.' // TRIM(filename)
                end if
                call write_pot_sic ( orb_rhor(:, jj), filename_complete)
            end do  
        end do
    end subroutine
end module