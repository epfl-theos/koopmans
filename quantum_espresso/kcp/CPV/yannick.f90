module yannick_print_orbr
    contains
    subroutine print_orbr(bec, nbsp, ispin, nbspx, lgam, is_empty)
        use grid_dimensions,          only: nnrx
        use nksic,                    only: orb_rhor
        use io_pot_sic_xml,           only: write_pot_sic
        use io_global,                only: ionode
        use cp_interfaces,            only: nksic_get_orbitalrho
        use gvecw,                    only: ngw
        use wavefunctions_module,     only: c0
        use twin_types
    
        type(twin_matrix), intent(in) :: bec 
        integer, intent(in)           :: nbspx, nbsp, ispin(nbspx)
        logical, intent(in)           :: lgam
        logical, intent(in)           :: is_empty

        ! local variables
        integer              :: j
        integer              :: jj
        integer              :: i
        character(len=1024)  :: filename
        character(len=1024)  :: filename_complete
    
        do j=1,nbsp,2
            call nksic_get_orbitalrho( ngw, nnrx, bec, ispin, nbsp, &
                        c0(:,j), c0(:,j+1), orb_rhor, j, j+1, lgam) 
            do jj = 1, 2
                i=j+jj-1
                if(ionode) then
                    write(*,*) "Yannick Debug own_module: Printing the density of orbital ", i
                end if
                if(i<10) then
                    write(filename, "(I1)") i
                else if(i<100) then
                    write(filename, "(I2)") i
                else
                    write(filename, "(I3)") i
                end if
                if(is_empty) then 
                    filename_complete = 'empty.' // TRIM(filename)
                else
                    filename_complete = 'occ.' // TRIM(filename)
                end if
                call write_pot_sic ( orb_rhor(:, jj), filename_complete)
            end do  
        end do
    end subroutine
end module