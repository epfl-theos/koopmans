!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE write_hamiltonian_real( ham, nbnd, ispin, empty )
  !---------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  USE io_global,           ONLY : ionode
  USE constants,           ONLY : AUTOEV
  !
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: ham(:,:)
  INTEGER, INTENT(IN) :: nbnd
  INTEGER, INTENT(IN) :: ispin
  LOGICAL, INTENT(IN) :: empty
  !
  INTEGER :: iunhr=24
  LOGICAL :: opnd
  CHARACTER(LEN=20) :: filename
  CHARACTER(LEN=33) :: header
  CHARACTER(LEN=9) :: cdate, ctime
  INTEGER :: i, j
  COMPLEX(DP) :: ham_(nbnd,nbnd)
  !
  !
  IF ( ionode ) THEN
    !
    INQUIRE( UNIT=iunhr, OPENED=opnd )
    IF ( opnd ) CALL errore( 'write_hamiltonian_real', 'file unit already opened', iunhr )
    !
    IF ( empty ) THEN
      WRITE( filename, FMT='( "ham_emp_", I1, ".dat" )' ) ispin
    ELSE
      WRITE( filename, FMT='( "ham_occ_", I1, ".dat" )' ) ispin
    ENDIF
    !
    OPEN( UNIT=iunhr, FILE=TRIM(filename), FORM='formatted', STATUS='unknown', ERR=101 )
    !
    ham_(:,:) = CMPLX( ham(1:nbnd,1:nbnd) ) * AUTOEV     ! Hartree to eV conversion
    !
    CALL date_and_tim( cdate, ctime )
    header = 'Written on '//cdate//' at '//ctime
    !
    WRITE( iunhr, * ) header
    WRITE( iunhr, * ) nbnd
    WRITE( iunhr, * ) 1
    WRITE( iunhr, '(I5)' ) 1
    !
    DO i = 1, nbnd
      DO j = 1, nbnd
        !
        WRITE( iunhr, FMT='( 5I5, 2F12.6 )' ) 0, 0, 0, i, j, ham_(j,i)
        !
      ENDDO
    ENDDO
    !
    CLOSE( iunhr )
    !
  ENDIF
  !
  RETURN
  !
101 CALL errore( 'write_hamiltonian_real', 'problem opening hamiltonian file', 1 ) 
  !
  !
END SUBROUTINE write_hamiltonian_real
!
!
!-----------------------------------------------------------------------
SUBROUTINE write_hamiltonian_cmplx( ham, nbnd, ispin, empty )
  !---------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  USE io_global,           ONLY : ionode
  USE constants,           ONLY : AUTOEV
  !
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: ham(:,:)
  INTEGER, INTENT(IN) :: nbnd
  INTEGER, INTENT(IN) :: ispin
  LOGICAL, INTENT(IN) :: empty
  !
  INTEGER :: iunhr=24
  LOGICAL :: opnd
  CHARACTER(LEN=20) :: filename
  CHARACTER(LEN=33) :: header
  CHARACTER(LEN=9) :: cdate, ctime
  INTEGER :: i, j
  COMPLEX(DP) :: ham_(nbnd,nbnd)
  !
  !
  IF ( ionode ) THEN
    !
    INQUIRE( UNIT=iunhr, OPENED=opnd )
    IF ( opnd ) CALL errore( 'write_hamiltonian_cmplx', 'file unit already opened', iunhr )
    !
    IF ( empty ) THEN
      WRITE( filename, FMT='( "ham_emp_", I1, ".dat" )' ) ispin
    ELSE
      WRITE( filename, FMT='( "ham_occ_", I1, ".dat" )' ) ispin
    ENDIF
    !
    OPEN( UNIT=iunhr, FILE=TRIM(filename), FORM='formatted', STATUS='unknown', ERR=101 )
    !
    ham_(:,:) = ham(1:nbnd,1:nbnd) * AUTOEV     ! Hartree to eV conversion
    !
    CALL date_and_tim( cdate, ctime )
    header = 'Written on '//cdate//' at '//ctime
    !
    WRITE( iunhr, * ) header
    WRITE( iunhr, * ) nbnd
    WRITE( iunhr, * ) 1
    WRITE( iunhr, '(I5)' ) 1
    !
    DO i = 1, nbnd
      DO j = 1, nbnd
        !
        WRITE( iunhr, FMT='( 5I5, 2F12.6 )' ) 0, 0, 0, i, j, ham_(j,i)
        !
      ENDDO
    ENDDO
    !
    CLOSE( iunhr )
    !
  ENDIF
  !
  RETURN
  !
101 CALL errore( 'write_hamiltonian_cmplx', 'problem opening hamiltonian file', 1 ) 
  !
  !
END SUBROUTINE write_hamiltonian_cmplx
