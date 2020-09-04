!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine save_in_electrons (iter, dr2)
  !-----------------------------------------------------------------------
  USE kinds,         ONLY: DP
  USE io_files,      ONLY: iunres, prefix
  USE ener,          ONLY: etot
  USE klist,         ONLY: nks
  USE control_flags, ONLY: io_level, conv_elec, tr2, ethr
  USE wvfct,         ONLY: nbnd, et
  USE scf,           ONLY: vnew
#ifdef EXX
  USE funct,         ONLY : exx_is_active
  USE exx,           ONLY : fock0, fock1, fock2, dexx, x_occupation
#endif
  implicit none
  character :: where * 20
  ! are we in the right place?
  integer :: ik, ibnd, ik_, iter
  ! counters
  ! last completed kpoint
  ! last completed iteration
  logical :: exst

  real(DP) :: dr2
  if ( io_level < 2 ) return
  !
  ! open recover file
  !
  call seqopn (iunres, 'restart', 'unformatted', exst)
  !
  ! save restart information
  !
  if (conv_elec) then
     !
     ! step on electrons has been completed. restart from ions
     !
     where = 'IONS'
     write (iunres) where
     write (iunres) ( (et (ibnd, ik), ibnd = 1, nbnd), ik = 1, nks)
     write (iunres) etot, tr2
     ! vnew = V(in)-V(out) is needed in the scf correction term to forces
     write (iunres) vnew%of_r
  else
     !
     ! save iteration number
     !
     ! iteration iter has been completed
     ik_ = 0
     where = 'ELECTRONS'
     write (iunres) where
     write (iunres) ( (et (ibnd, ik), ibnd = 1, nbnd), ik = 1, nks)
     write (iunres) iter, ik_, dr2, tr2, ethr

  endif
#ifdef EXX
  write (iunres)  exx_is_active(), fock0, fock1, fock2, dexx
  write (iunres) ( (x_occupation (ibnd, ik), ibnd = 1, nbnd), ik = 1, nks)
#endif
  !

  close (unit = iunres, status = 'keep')
  return

end subroutine save_in_electrons
