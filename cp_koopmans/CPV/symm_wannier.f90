!
! Copyright (C) 2002-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE symm_wannier( wfc, num_states, emp) 
  !---------------------------------------------------------------------
  !
  ! ...  This routine imposes the Bloch symmetry on the Wannier
  ! ...  functions in two steps:
  ! ...
  ! ...  1) finds the WFs inside a reference primitive cell
  ! ...  2) builds the other WFs by imposing w_R(r) = w_0(r-R), where
  ! ...     R are the lattice vectors of the corresponding primitive cell
  !
  !
  USE kinds,                ONLY : DP
  USE electrons_base,       ONLY : nspin
  USE gvecw,                ONLY : ngw
  USE input_parameters,     ONLY : mp1, mp2, mp3, &
                                   offset_centers_occ, offset_centers_emp 
  USE reciprocal_vectors,   ONLY : gx
  USE cell_base,            ONLY : at, alat
  USE constants,            ONLY : tpi
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : intra_image_comm
  USE constants,            ONLY : BOHR_RADIUS_ANGS
  USE centers_and_spreads,  ONLY : centers_occ, centers_emp, &
                                   spreads_occ, spreads_emp
  !
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: num_states
  LOGICAL, INTENT(IN) :: emp
  COMPLEX(DP), INTENT(INOUT) :: wfc(ngw,num_states)
  !
  INTEGER :: norb
  INTEGER :: norb_pc     ! number of ref WFs
  REAL(DP), ALLOCATABLE :: centers(:,:), spreads(:)
  REAL(DP), ALLOCATABLE :: centers_(:,:), spreads_(:)
  COMPLEX(DP), ALLOCATABLE :: wfc_aux(:,:)
  INTEGER :: n, ig
  INTEGER :: i, j, k
  INTEGER :: counter
  REAL(DP) :: offset(3)
  REAL(DP) :: rvect(3), Rvec(3)
  CHARACTER(3) :: typ='occ'
  COMPLEX(DP) :: imag = (0.D0,1.D0)
  !
  !
  norb = num_states / nspin
  norb_pc = norb / (mp1*mp2*mp3)
  ALLOCATE( centers(3,norb) )
  ALLOCATE( centers_(3,norb) )
  ALLOCATE( spreads(norb), spreads_(norb) )
  ALLOCATE( wfc_aux(ngw,norb_pc) )
  !
  IF ( emp ) typ = 'emp'
  !
  WRITE( stdout, 101 ) typ
  101 FORMAT( //, 3x, 'Forcing Bloch symmetry on the ',a3,' Wannier functions', / &
                  3x, '---------------------------------------------------' )
  !
  IF ( ionode ) CALL read_wannier_centers(centers, norb, emp)
  IF ( ionode ) CALL read_wannier_spreads(spreads_, norb, emp)
  CALL mp_bcast( centers, ionode_id, intra_image_comm )
  CALL mp_bcast( spreads_, ionode_id, intra_image_comm )
  spreads(:) = 0.D0
  !
  WRITE( stdout, * )
  WRITE( stdout, '( 3x, "Reference WFs found (crystal units):" )' )
  !
  offset = (/0.,0.,0./)
  IF ( emp ) THEN
    IF ( offset_centers_emp ) offset = (/0.1,0.1,0.1/)
  ELSE
    IF ( offset_centers_occ ) offset = (/0.1,0.1,0.1/)
  ENDIF
  !
  counter = 0
  DO n = 1, norb
    !
    ! shift the WFs inside the supercell (0,1)x(0,1)x(0,1)
    !
    centers(1,n) = centers(1,n) + offset(1) - floor(centers(1,n))
    centers(2,n) = centers(2,n) + offset(2) - floor(centers(2,n))
    centers(3,n) = centers(3,n) + offset(3) - floor(centers(3,n))
    !
    centers(:,n) = centers(:,n) - offset(:)
    !
    ! identify the WFs inside a ref primitive cell
    !
    IF ( (centers(1,n) + offset(1)) * mp1 - 1 < 1.e-3 .and. &
         (centers(2,n) + offset(2)) * mp2 - 1 < 1.e-3 .and. &
         (centers(3,n) + offset(3)) * mp3 - 1 < 1.e-3 ) THEN
      !
      counter = counter + 1
      wfc_aux(:,counter) = wfc(:,n)
      centers_(:,counter) = centers(:,n)
      spreads(counter) = spreads_(n)
      !
      WRITE ( stdout, 201 ) n, centers(:,n)
      !
    ENDIF
    !
  ENDDO
  !
  201 FORMAT( 5x, 'Wannier function # ', i3, ' :', 4x, '(', 3f10.6, ' )' )  
  !
  IF ( counter .ne. norb_pc ) THEN
    CALL errore( 'symm_wannier', 'Wrong number of ref Wannier functions', counter )
  ENDIF
  !
  centers(:,:) = 0.D0
  centers(:,:norb_pc) = centers_(:,:norb_pc)
  DEALLOCATE( centers_ )
  !
  wfc(:,1:norb_pc) = wfc_aux(:,:)
  DEALLOCATE( wfc_aux )
  !
  !
  ! Now we build all the other WFs by imposing w_R(g) = e^(igR) w_0(g)
  !
  WRITE( stdout, 301 )
  301 FORMAT( /, 3x, 'Building the other WFs  --->  w_Rn(r) = w_0n(r-R)', / )
  !
  DO i = 1, mp1
    DO j = 1, mp2
      DO k = 1, mp3
        !
        IF ( i == 1 .and. j == 1 .and. k == 1 ) CYCLE
        !
        rvect = (/ float(i-1) / mp1, &
                   float(j-1) / mp2, &
                   float(k-1) / mp3 /)
        Rvec(:) = rvect(:)
        CALL cryst_to_cart( 1, Rvec, at, 1 )
        !
        DO n = 1, norb_pc
          ! 
          counter = counter + 1
          centers(:,counter) = centers(:,n) + rvect(:)
          spreads(counter) = spreads_(n)
          !
          DO ig = 1, ngw
            !
            wfc(ig,counter) = wfc(ig,n) * EXP( imag * tpi * DOT_PRODUCT( gx(:,ig), Rvec ) )
            !
          ENDDO
        ENDDO
        !
      ENDDO
    ENDDO
  ENDDO
  !
  IF ( nspin == 2 ) wfc(:,norb+1:) = wfc(:,1:norb)
  !
  !
  CALL cryst_to_cart( norb, centers, at, 1 )
  centers = centers * alat * BOHR_RADIUS_ANGS 
  !
  IF ( emp ) THEN
    centers_emp = centers
    spreads_emp = spreads
  ELSE
    centers_occ = centers
    spreads_occ = spreads
  ENDIF
  !
  !  
END SUBROUTINE symm_wannier
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
  USE kinds,                ONLY : DP
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
SUBROUTINE write_centers_and_spreads( num_wann, centers, spreads, emp )
  !---------------------------------------------------------------------
  !
  ! ...  This routine prints the centers and spreads of Wannier functions
  ! ...  (occupied or empty) to file. The units are Ang.
  !
  ! ...  emp = .true. when printing empty states
  !
  USE kinds,                ONLY : DP
  USE io_files,             ONLY : prefix
  USE io_global,            ONLY : ionode, stdout
  !
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: num_wann        ! number of Wannier functions
  REAL(DP), INTENT(IN) :: centers(3,num_wann) 
  REAL(DP), INTENT(IN) :: spreads(num_wann) 
  LOGICAL, INTENT(IN) :: emp
  !
  CHARACTER(LEN=268) :: filename
  CHARACTER(LEN=9)  :: cdate, ctime
  CHARACTER(LEN=100) :: header
  INTEGER :: n
  !
  !
  IF ( emp ) THEN
    filename = TRIM(prefix)//'_cp_centers_emp.xyz'
  ELSE
    filename = TRIM(prefix)//'_cp_centers.xyz'
  ENDIF
  !
  CALL date_and_tim( cdate, ctime )
  header = 'Wannier centers and spreads, written by CP on '//cdate//' at '//ctime
  !
  IF ( ionode ) THEN
    !
    OPEN( 300, file=filename, status='unknown' )
    !
    WRITE( 300, '(i6)' ) num_wann
    WRITE( 300, * ) header
    !
    WRITE( 300, 55 ) ( centers(:,n), spreads(n), n = 1, num_wann )
    !
    CLOSE( 300 )
    !
  ENDIF
  !
  !
  55 FORMAT( 'X', 5x, 3f17.8, 2x, f17.8 )
  !
  !
END SUBROUTINE write_centers_and_spreads
