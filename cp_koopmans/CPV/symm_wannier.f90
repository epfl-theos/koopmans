!
! Copyright (C) 2002-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE symm_wannier_x( wfc, num_states, emp ) 
  !---------------------------------------------------------------------
  !
  ! ...  This routine imposes the Bloch symmetry on the variational
  ! ...  orbitals in two steps:
  ! ...
  ! ...  1) finds the orbitals inside a reference primitive cell; if
  ! ...     read_centers = .false. then the first N orbitals (where N is
  ! ...     the number of PC orbitals) are taken as reference 
  ! ...  2) builds the other orbitals by imposing w_R(r) = w_0(r-R), where
  ! ...     R are the lattice vectors of the corresponding primitive cell
  !
  !
  USE kinds,                ONLY : DP
  USE electrons_base,       ONLY : nspin
  USE gvecw,                ONLY : ngw
  USE input_parameters,     ONLY : mp1, mp2, mp3, read_centers, &
                                   offset_centers_occ, offset_centers_emp 
  USE reciprocal_vectors,   ONLY : gx
  USE cell_base,            ONLY : at, alat
  USE constants,            ONLY : tpi
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : intra_image_comm
  USE constants,            ONLY : BOHR_RADIUS_ANGS
  USE centers_and_spreads,  ONLY : read_wannier_centers, read_wannier_spreads
  !
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: wfc(:,:)
  INTEGER, INTENT(IN) :: num_states
  LOGICAL, INTENT(IN) :: emp
  !
  INTEGER :: norb
  INTEGER :: norb_pc     ! number of ref orbitals
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
  !
  IF ( emp ) typ = 'emp'
  !
  WRITE( stdout, 101 ) typ
  !
  IF ( read_centers ) THEN
    !
    ALLOCATE( centers(3,norb) )
    ALLOCATE( centers_(3,norb) )
    ALLOCATE( spreads(norb), spreads_(norb) )
    ALLOCATE( wfc_aux(ngw,norb_pc) )
    !
    IF ( ionode ) CALL read_wannier_centers(centers, norb, emp)
    IF ( ionode ) CALL read_wannier_spreads(spreads_, norb, emp)
    CALL mp_bcast( centers, ionode_id, intra_image_comm )
    CALL mp_bcast( spreads_, ionode_id, intra_image_comm )
    spreads(:) = 0.D0
    !
    WRITE( stdout, * )
    WRITE( stdout, '( 3x, "Reference orbitals found (crystal units):" )' )
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
      ! shift the orbitals inside the supercell (0,1)x(0,1)x(0,1)
      !
      centers(1,n) = centers(1,n) + offset(1) - floor(centers(1,n))
      centers(2,n) = centers(2,n) + offset(2) - floor(centers(2,n))
      centers(3,n) = centers(3,n) + offset(3) - floor(centers(3,n))
      !
      centers(:,n) = centers(:,n) - offset(:)
      !
      ! identify the orbitals inside a ref primitive cell
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
    !
    IF ( counter .ne. norb_pc ) THEN
      CALL errore( 'symm_wannier', 'Wrong number of ref orbitals', counter )
    ENDIF
    !
    centers(:,:) = 0.D0
    centers(:,:norb_pc) = centers_(:,:norb_pc)
    DEALLOCATE( centers_ )
    !
    wfc(:,1:norb_pc) = wfc_aux(:,:)
    DEALLOCATE( wfc_aux )
    !
  ELSE
    !
    WRITE( stdout, * )
    WRITE( stdout, 202 ) norb_pc
    ! 
  ENDIF
  !
  !
  ! Now we build all the other orbitals by imposing w_R(g) = e^(igR) w_0(g)
  !
  WRITE( stdout, 301 )
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
          IF ( read_centers ) THEN
            !
            counter = counter + 1
            centers(:,counter) = centers(:,n) + rvect(:)
            spreads(counter) = spreads_(n)
            !
          ENDIF
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
  IF ( read_centers ) THEN
    !
    CALL cryst_to_cart( norb, centers, at, 1 )
    centers = centers * alat * BOHR_RADIUS_ANGS 
    !
  ENDIF
  !
  !
101 FORMAT( //, 3x, 'Forcing Bloch symmetry on the ', a3, ' orbitals', / &
                3x, '------------------------------------------' )
201 FORMAT( 5x, 'orbital # ', i3, ' :', 4x, '(', 3f10.6, ' )' )  
202 FORMAT( 3x, 'Taking the first ', i4, 'orbitals as reference' )
301 FORMAT( /, 3x, 'Building the other orbitals  --->  w_Rn(r) = w_0n(r-R)', / )
  !
  !  
END SUBROUTINE symm_wannier_x
