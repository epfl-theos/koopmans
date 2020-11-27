subroutine pw2cp_empty (seedname, iknum, ikstart, empty_start_index, nbnd_empty_to_save, spin_component, force_spin_symmetry)
  !
  USE io_global,       only : stdout
  USE kinds,           only : DP
  use io_files,        only : iunwfc, iunatsicwfc, nwordwfc, nwordwann
  USE cell_base,       only : omega, tpiba2
  use gvect,           only : g, ngm, ecutwfc
  use gsmooth,         only: nls, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
  use wavefunctions_module, only : evc, psic
  use wvfct,           only : nbnd, npwx, npw, igk, g2kin
  use klist,           only : nkstot, xk, wk
  use lsda_mod,        only : nspin, isk
  use read_wannier_files
  !
  ! local variables
  !
  character(len=50), intent (in) :: seedname
  integer, intent (in) :: iknum, ikstart, empty_start_index, nbnd_empty_to_save
  !
  integer :: i, j, nn, ik, ibnd, ikevc, counter, counter_tot
  complex (DP), allocatable:: evc_tmp(:,:)
  complex (DP), allocatable:: c0_tmp_empty(:,:)
  ! 
  character(len=4) :: spin_component
  CHARACTER(LEN=256) :: filename
  logical :: force_spin_symmetry
  !
  ! main program
  !
  do ik=1, iknum
     !
     ikevc = ik + ikstart - 1
     call davcio (evc, nwordwfc, iunwfc, ikevc, -1)
     call gk_sort (xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
     !
     allocate(evc_tmp (npw, (nbnd-empty_start_index+1)))
     counter = 0
     do ibnd = 1, nbnd
        !
        if (ibnd .ge. empty_start_index) then
           counter = counter + 1
           evc_tmp(:,counter) = evc(:, ibnd)
        endif
        ! 
     enddo 
     !
  enddo ! k-points
  !
  counter_tot = counter
  !
  if (counter_tot .ne. (nbnd-empty_start_index+1)) then
      call errore('pw2cp',' counter_tot \= nbnd-emp_start_index+1 ', 1) 
  endif
  !
  if (nspin == 2) then
     ! 
     if (force_spin_symmetry) then
        allocate (c0_tmp_empty (npw, nbnd_empty_to_save*2))
     else
        allocate (c0_tmp_empty (npw, nbnd_empty_to_save))
     endif
  else
     allocate (c0_tmp_empty (npw, nbnd_empty_to_save))
  endif
  !
  c0_tmp_empty(:,:) = (0.0, 0.0)
  !
  do ibnd=1, nbnd_empty_to_save
     !
     c0_tmp_empty (:, ibnd) = evc_tmp (:, ibnd)
     !
  enddo
  !
  if (nspin == 2) then 
     !
     if (force_spin_symmetry) then
        counter = 0
        !
        do ibnd=nbnd_empty_to_save+1, nbnd_empty_to_save*2
           !
           counter = counter + 1
           c0_tmp_empty (:, ibnd) = c0_tmp_empty (:, counter)
           ! 
        enddo
        !
        call writeempty_state(c0_tmp_empty, nbnd_empty_to_save*2, 'evc_empty.dat')
        !
     else
        !
        filename='evc_empty_'//TRIM(spin_component)//'.dat'
        call writeempty_state(c0_tmp_empty, nbnd_empty_to_save, filename)
        !
     endif
     !
  else
     !  
     call writeempty_state(c0_tmp_empty, nbnd_empty_to_save, 'evc_empty.dat')
     !
  endif 
  !
  write(stdout, *) nspin, counter_tot, nbnd_empty_to_save 
  write(stdout, *) c0_tmp_empty(1:3, 1 )
  write(stdout, *) c0_tmp_empty(1:3, nbnd_empty_to_save+1 )
  !
  deallocate ( c0_tmp_empty, evc_tmp)
  !
  return
  !
end subroutine pw2cp_empty
