!
! Copyright (C) 2002-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE centers_and_spreads
  !---------------------------------------------------------------------
  !
  ! ... This module contains all the routines and variables important for
  ! ... the calculation of centers and spreads of the variational orbitals
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE reciprocal_vectors,   ONLY : g2_g, g, gx
  USE electrons_base,       ONLY : nspin
  USE constants,            ONLY : BOHR_RADIUS_ANGS
  !
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: read_wannier_centers, read_wannier_spreads
  !
  ! ...  end of module-scope declarations
  !
  !---------------------------------------------------------------------
  !
CONTAINS
  !
  !
  SUBROUTINE read_wannier_centers( centers, num_wann, emp )
    !---------------------------------------------------------------------
    !
    ! ...  This routine reads the centers of Wannier functions from .xyz
    ! ...  file print out by Wannier90, fold them into the R=0 primitive
    ! ...  cell and gives them in output (in crystal units)
    !
    ! ...  emp = .true. when reading empty states
    !
    USE kinds,                ONLY : DP
    USE cell_base,            ONLY : bg, alat
    USE constants,            ONLY : BOHR_RADIUS_ANGS
    USE io_files,             ONLY : prefix
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: num_wann        ! number of Wannier functions 
    LOGICAL, INTENT(IN) :: emp             
    !
    REAL(DP), INTENT(OUT) :: centers(3,num_wann)
    !
    LOGICAL :: exst
    INTEGER :: n 
    CHARACTER(LEN=268) :: filename
    CHARACTER(LEN=256) :: input_line
    !
    !
    IF ( emp ) THEN
      filename = trim(prefix)//'_emp_centres.xyz'
    ELSE
      filename = trim(prefix)//'_centres.xyz'
    ENDIF
    !
    INQUIRE( file=filename, exist=exst )
    !
    IF ( .not. exst ) CALL errore( 'read_wannier_centers', 'File not found', 1 )
    !
    OPEN( 100, file=filename, form='formatted', status='old' )
    !
    READ( 100, *, end=10, err=20 )    ! skip 1st line
    READ( 100, *, end=10, err=20 )    ! skip 2nd line
         !
    DO n = 1, num_wann
      !
      READ( 100, '(a256)', end=10, err=20 ) input_line
      !
      IF ( input_line(1:1) .ne. 'X' ) CALL errore( 'read_wannier_centers', &
              'X must precede each Wannier center line', 1 )
      !
      READ( input_line(2:), *, end=10, err=20 ) centers(:,n)
      !
    ENDDO
    !
    READ( 100, * ) input_line
    IF ( input_line(1:1) == 'X' ) CALL errore( 'read_wannier_centers', &
            'Missing some center!', 1 )
    !
    CLOSE( 100 )
    !
    !
    centers = centers / ( alat * BOHR_RADIUS_ANGS )
    !
    DO n = 1, num_wann
      !
      CALL cryst_to_cart( 1, centers(:,n), bg, -1 )
      !
    ENDDO
    !
    RETURN
    !
  10  CALL errore ( 'read_wannier_centers', 'end of file while reading', 1 )
  20  CALL errore ( 'read_wannier_centers', 'error while reading', 1 )
    !
    !
  END SUBROUTINE read_wannier_centers
  !
  !
  SUBROUTINE read_wannier_spreads( spreads, num_wann, emp )
    !---------------------------------------------------------------------
    !
    ! ...  This routine reads the spreads of Wannier functions from .wout
    ! ...  file print out by Wannier90, gives them in output (in Ang^2)
    !
    ! ...  emp = .true. when reading empty states
    !
    USE io_files,             ONLY : prefix
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: num_wann        ! number of Wannier functions 
    LOGICAL, INTENT(IN) :: emp
    !
    REAL(DP), INTENT(OUT) :: spreads(num_wann)
    !
    LOGICAL :: exst
    INTEGER :: n 
    CHARACTER(LEN=268) :: filename
    CHARACTER(LEN=256) :: input_line
    !
    !
    IF ( emp ) THEN
      filename = trim(prefix)//'_emp.wout'
    ELSE
      filename = trim(prefix)//'.wout'
    ENDIF
    !
    INQUIRE( file=filename, exist=exst )
    !
    IF ( .not. exst ) CALL errore( 'read_wannier_spreads', 'File not found', 1 )
    !
    OPEN( 200, file=filename, form='formatted', status='old' )
    !
    READ( 200, '(a256)', end=10, err=20 ) input_line
    DO WHILE ( input_line .ne. ' Final State' )
      READ( 200, '(a256)', end=10, err=20 ) input_line
    ENDDO
    !
    DO n = 1, num_wann
      !
      READ( 200, '(a256)', end=10, err=20 ) input_line
      READ( input_line(65:), * ) spreads(n)
      !
    ENDDO
    !
    CLOSE( 200 )
    !
    !
    RETURN
    !
  10  CALL errore ( 'read_wannier_spreads', 'end of file while reading', 1 )
  20  CALL errore ( 'read_wannier_spreads', 'error while reading', 1 )
    !
    !
  END SUBROUTINE read_wannier_spreads
  !
  !
END MODULE centers_and_spreads
