!
! Copyright (C) 2001-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE read_config_from_file()
  !-----------------------------------------------------------------------
  !
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : stdout
  USE ions_base,      ONLY : nat, ityp, tau
  USE basis,          ONLY : startingconfig
  USE cell_base,      ONLY : at, bg, omega
  USE cellmd,         ONLY : at_old, omega_old, lmovecell
  USE io_files,       ONLY : tmp_dir, prefix
  USE pw_restart,     ONLY : pw_readfile
  !
  IMPLICIT NONE
  !
  INTEGER :: ierr
  !
  !
  IF ( TRIM( startingconfig ) /= 'file' ) RETURN
  !
  WRITE( stdout, '(/5X,"Atomic positions and unit cell read from directory:"/5X,A)') &
      TRIM( tmp_dir ) // TRIM( prefix ) // ".save/"
  !
  ! ... check if restart file is present, if yes read config parameters
  !
  CALL pw_readfile( 'config', ierr )
  !
  IF ( ierr > 0 ) THEN
     !
     WRITE( stdout, '(5X,"Failed to open directory or to read data file! " &
                       & "Using input atomic positions and unit cell",/)' )
     !
     RETURN
     !
  END IF
  !
  WRITE( stdout, * )
  !
  IF ( lmovecell ) THEN
     !
     ! ... input value of at and omega (currently stored in xxx_old variables)
     ! ... must be used to initialize G vectors and other things
     ! ... swap xxx and xxx_old variables and scale the atomic position to the
     ! ... input cell shape in order to check the symmetry.
     !
     CALL cryst_to_cart( nat, tau, bg, - 1 )
     !
     CALL dswap( 9, at, 1, at_old,1  )
     CALL dswap( 1, omega, 1, omega_old, 1 )
     !
     CALL cryst_to_cart( nat, tau, at, + 1 )
     !
     CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE read_config_from_file
