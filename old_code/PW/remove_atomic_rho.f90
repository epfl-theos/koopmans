!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine remove_atomic_rho
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE io_global, ONLY: stdout
  USE io_files, ONLY: output_drho
  USE kinds, ONLY: DP
  USE gvect, ONLY: nrxx
  USE lsda_mod, ONLY: nspin
  USE scf, ONLY: rho
  USE io_rho_xml, ONLY : write_rho
  implicit none

  real(DP), allocatable :: work (:,:)
  ! workspace, is the difference between the charge density
  ! and the superposition of atomic charges

  allocate ( work( nrxx, 1 ) )
  work = 0.d0
  !
  IF ( nspin > 1 ) CALL errore &
       ( 'remove_atomic_rho', 'spin polarization not allowed in drho', 1 )

  WRITE( stdout, '(/5x,"remove atomic charge density from scf rho")')
  !
  !     subtract the old atomic charge density
  !
  call atomic_rho (work, nspin)
  !
  work = rho%of_r - work
  !
  call write_rho ( work, 1, output_drho )
  !
  deallocate(work)
  return

end subroutine remove_atomic_rho

