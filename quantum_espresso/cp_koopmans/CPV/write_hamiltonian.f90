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


!-----------------------------------------------------------------------
   SUBROUTINE write_ham_emp_xml (nspin, nudx_emp, lambda_emp, desc_emp, fname)
!-----------------------------------------------------------------------
!   
    USE iotk_module
    USE io_global,            ONLY : ionode, stdout, ionode_id
    USE io_files,             ONLY : prefix, iunpun, outdir, xmlpun
    USE xml_io_base,          ONLY : wfc_filename, restart_dir, create_directory, kpoint_dir
    USE twin_types !added:giovanni 
    USE cp_main_variables,    ONLY : collect_lambda, descla
    USE descriptors,          ONLY : descla_siz_ 
    USE control_flags,        ONLY : ndw
    USE mp_global,            ONLY : intra_image_comm
    USE mp,                   ONLY : mp_bcast
    ! 
    IMPLICIT NONE 
    !
    INTEGER, INTENT(IN) :: nspin
    INTEGER, INTENT(IN) :: nudx_emp
    CHARACTER(256), INTENT(IN) :: fname
    TYPE(twin_matrix), INTENT(IN) :: lambda_emp(nspin)
    !
    INTEGER :: ik = 1 ! There only 1 kpoint in CP 
    CHARACTER(LEN=256)    :: dirname, filename
    CHARACTER(LEN=4)      :: cspin
    INTEGER :: iss, ierr
    REAL(DP), ALLOCATABLE :: mrepl(:,:)
    COMPLEX(DP), ALLOCATABLE :: mrepl_c(:,:) !added:giovanni
    INTEGER, INTENT (IN) :: desc_emp( descla_siz_ , 2 )
    !
    dirname = restart_dir( outdir, ndw )
    CALL create_directory( kpoint_dir( dirname, ik ) )
    !
    IF (ionode) THEN 
       CALL iotk_open_write( iunpun, FILE = TRIM( dirname ) // '/' // &
                           & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
    ENDIF
    !
    CALL mp_bcast( ierr, ionode_id, intra_image_comm )
    CALL errore( 'write_ham_emp_xml ', 'cannot open restart file for writing', ierr )
    !
    DO iss = 1, nspin 
    !
    IF(.not. lambda_emp(1)%iscmplx) THEN
        ALLOCATE( mrepl( nudx_emp, nudx_emp ) )    
        CALL collect_lambda( mrepl, lambda_emp(iss)%rvec(:,:), desc_emp(:,iss) )
    ELSE
        ALLOCATE( mrepl_c( nudx_emp, nudx_emp) )
        CALL collect_lambda( mrepl_c, lambda_emp(iss)%cvec(:,:), desc_emp(:,iss))
    ENDIF
    !
    cspin = iotk_index( iss )
    ! 
    IF ( ionode ) THEN
       !
       !
       ! Changes by Nicolas Poilvert, Sep. 2010 for printing the lambda
       ! matrix at current time step into a formatted file.
       ! This matrix corresponds to the Hamiltonian matrix  in the case
       ! of Self-Interaction. Only in the basis  of minimizing orbitals
       ! do this matrix has an interpretation.
       !
       IF ( nspin == 1 ) THEN
           !
           filename = TRIM( wfc_filename( ".", fname, ik, EXTENSION='xml' ) )
           !
       ELSE
           !
           filename = TRIM( wfc_filename( ".", fname, ik, iss, EXTENSION='xml' ) )
           !
       ENDIF
       !
       CALL iotk_link( iunpun, "HAMILTONIAN" // TRIM( cspin ), &
                       filename, CREATE = .TRUE., BINARY = .FALSE. )
       !
       !
       IF(allocated(mrepl)) THEN
           CALL iotk_write_dat( iunpun, &
                                "HAMILTONIAN" // TRIM( cspin ), mrepl )
           !IF ( write_hr ) CALL write_hamiltonian( mrepl, nupdwn(iss), iss, .false. )
           !
           !
       ELSE IF(allocated(mrepl_c)) THEN
           CALL iotk_write_dat( iunpun, &
                                "HAMILTONIAN" // TRIM( cspin ), mrepl_c )
           !IF ( write_hr ) CALL write_hamiltonian( mrepl_c, nupdwn(iss), iss, .false. )
       ENDIF
       !
    ENDIF
    !
    IF(allocated(mrepl)) THEN
      DEALLOCATE( mrepl )
    ENDIF
    IF(allocated(mrepl_c)) THEN
      DEALLOCATE( mrepl_c )
    ENDIF
    !
    ENDDO
    !
    IF ( ionode ) CALL iotk_close_write( iunpun )
    !
END SUBROUTINE write_ham_emp_xml

