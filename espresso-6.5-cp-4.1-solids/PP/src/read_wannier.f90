!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Written by Riccardo De Gennaro, EPFL (Sept 2020).
!
!
!-----------------------------------------------------------------------
MODULE read_wannier
  !---------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  !
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: read_wannier_chk
  !
  ! Remember to mp_bcast possible additions
  !
  INTEGER, PUBLIC :: num_bands
  INTEGER, PUBLIC :: num_wann
  INTEGER, PUBLIC :: num_kpts
  INTEGER, PUBLIC :: kgrid(3)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: u_mat(:,:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: u_mat_opt(:,:,:)
  !
  CONTAINS
  !
  !---------------------------------------------------------------------
  SUBROUTINE read_wannier_chk( seedname )
    !-------------------------------------------------------------------
    !
    ! ...  parse the Wannier90 chk file
    !
    USE io_global,           ONLY : ionode, ionode_id
    USE mp_global,           ONLY : intra_image_comm
    USE mp,                  ONLY : mp_bcast
    USE cell_base,           ONLY : bg
    USE constants,           ONLY : eps8
    USE klist,               ONLY : nkstot, xk
    USE lsda_mod,            ONLY : nspin
    !
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), INTENT(IN) :: seedname
    !
    CHARACTER(LEN=33) :: header
    INTEGER :: i, j, nkp, nn
    INTEGER :: chk_unit=124
    INTEGER :: num_exclude_bands, nntot
    LOGICAL :: exst
    LOGICAL :: checkpoint, have_disentangled
    LOGICAL, ALLOCATABLE :: exclude_bands(:)
    REAL(DP) :: at_(3,3), bg_(3,3)
    REAL(DP), ALLOCATABLE :: kpt_latt(:,:)
    REAL(DP), ALLOCATABLE :: kaux(:,:)
    REAL(DP), ALLOCATABLE :: centers(:,:)
    REAL(DP), ALLOCATABLE :: spreads(:)
    COMPLEX(DP), ALLOCATABLE :: uu_prod(:,:)
    COMPLEX(DP), ALLOCATABLE :: m_mat(:,:,:,:)
    !
    REAL(DP) :: omega_invariant           ! disentanglement parameters
    LOGICAL, ALLOCATABLE :: lwindow(:,:)  ! disentanglement parameters
    INTEGER, ALLOCATABLE :: ndimwin(:)    ! disentanglement parameters
    !
    !
    IF ( ionode ) THEN
      !
      INQUIRE( FILE=trim(seedname)//'.chk', EXIST=exst )
      IF ( .not. exst ) CALL errore( 'read_wannier_chk', 'chk file not found', 1 )
      !
      OPEN( UNIT=chk_unit, FILE=trim(seedname)//'.chk', STATUS='old', FORM='unformatted' )
      !
      READ( chk_unit ) header                   ! date and time
      READ( chk_unit ) num_bands                ! number of bands
      READ( chk_unit ) num_exclude_bands        ! number of excluded bands
      !
      IF ( num_exclude_bands .lt. 0 ) &
        CALL  errore( 'read_wannier_chk', 'Invalid value for num_exclude_bands', &
                      num_exclude_bands )
      !
      ALLOCATE( exclude_bands(num_exclude_bands) )
      exclude_bands(:) = .FALSE.
      !
      READ( chk_unit ) ( exclude_bands(i), i=1,num_exclude_bands )   ! excluded bands
      READ( chk_unit ) (( at_(i,j), i=1,3 ), j=1,3 )                 ! prim real latt vectors
      READ( chk_unit ) (( bg_(i,j), i=1,3 ), j=1,3 )                 ! prim recip latt vectors
      READ( chk_unit ) num_kpts                                      ! num of k-points
      !
      IF ( nspin == 2 ) THEN
        IF ( num_kpts .ne. nkstot/2 ) &
        CALL errore( 'read_wannier_chk', 'Invalid value for num_kpts', num_kpts )
      ELSE
        IF ( num_kpts .ne. nkstot ) &
        CALL errore( 'read_wannier_chk', 'Invalid value for num_kpts', num_kpts )
      ENDIF
      !
      READ( chk_unit ) ( kgrid(i), i=1,3 )                           ! MP grid
      !
      ALLOCATE( kpt_latt(3,num_kpts) )
      ALLOCATE( kaux(3,num_kpts) )
      !
      READ( chk_unit ) (( kpt_latt(i,nkp), i=1,3 ), nkp=1,num_kpts )
      !
      kaux(:,:) = kpt_latt(:,:)
      CALL cryst_to_cart( num_kpts, kaux, bg, 1 )
      IF ( ANY( kaux .ne. xk ) ) &
        CALL errore( 'read_wannier_chk', 'Invalid value for kpt_latt', 1 )
      !
      READ( chk_unit ) nntot                                         ! nntot
      READ( chk_unit ) num_wann                                      ! num of WFs
      READ( chk_unit ) checkpoint                                    ! checkpoint
      READ( chk_unit ) have_disentangled     ! .true. if disentanglement has been performed
      !
      IF ( have_disentangled ) THEN 
        !
        READ( chk_unit ) omega_invariant                             ! omega invariant
        !
        ALLOCATE( lwindow(num_bands,num_kpts) )
        ALLOCATE( ndimwin(num_kpts) )
        !
        READ( chk_unit ) (( lwindow(i,nkp), i=1,num_bands ), nkp=1,num_kpts )
        READ( chk_unit ) ( ndimwin(nkp), nkp=1,num_kpts )
        !
        ALLOCATE( u_mat_opt(num_bands,num_wann,num_kpts) )
        READ ( chk_unit ) ((( u_mat_opt(i,j,nkp), i=1,num_wann ), &  ! optimal U-matrix
                                                  j=1,num_wann ), &
                                                  nkp=1,num_kpts )
        !
      ENDIF
      !
      ALLOCATE( u_mat(num_wann,num_wann,num_kpts) )
      READ ( chk_unit ) ((( u_mat(i,j,nkp), i=1,num_wann ), &        ! U-matrix
                                            j=1,num_wann ), &
                                            nkp=1,num_kpts )
      !
      ! check that U(k) is unitary for each k
      ALLOCATE( uu_prod(num_wann,num_wann) )
      DO nkp = 1, num_kpts
        uu_prod = MATMUL( u_mat(:,:,nkp), CONJG(TRANSPOSE( u_mat(:,:,nkp) )) )
        DO i = 1, num_wann
          DO j = 1, num_wann
            IF ( i == j ) THEN
              IF ( ( ABS(DBLE( uu_prod(i,j) - 1 )) .ge. eps8 ) .or. &
                  ( ABS(AIMAG( uu_prod(i,j) )) .ge. eps8 ) ) &
                CALL errore( 'read_wannier_chk', 'u_mat is not unitary', nkp )
            ELSE
              IF ( ( ABS(DBLE( uu_prod(i,j) )) .ge. eps8 ) .or. &
                  ( ABS(AIMAG( uu_prod(i,j) )) .ge. eps8 ) ) &
                CALL errore( 'read_wannier_chk', 'u_mat is not unitary', nkp )
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      !
      ALLOCATE( m_mat(num_wann,num_wann,nntot,num_kpts) )            ! M-matrix
      !
      READ( chk_unit ) (((( m_mat(i,j,nn,nkp), i=1,num_wann ), &
                                               j=1,num_wann ), &
                                               nn=1,nntot ), &
                                               nkp=1,num_kpts )
      !
      ALLOCATE( centers(3,num_wann) )                                ! Wannier centers
      ALLOCATE( spreads(num_wann) )                                  ! Wannier spreads
      !
      READ( chk_unit ) (( centers(i,j), i=1,3 ), j=1,num_wann )
      READ( chk_unit ) ( spreads(i), i=1,num_wann )
      !
      CLOSE( chk_unit )
      !
    ENDIF
    !
    !
    CALL mp_bcast( num_bands, ionode_id, intra_image_comm )
    CALL mp_bcast( num_wann, ionode_id, intra_image_comm )
    CALL mp_bcast( num_kpts, ionode_id, intra_image_comm )
    CALL mp_bcast( kgrid, ionode_id, intra_image_comm )
    CALL mp_bcast( u_mat, ionode_id, intra_image_comm )
    CALL mp_bcast( u_mat_opt, ionode_id, intra_image_comm )
    !
    !
  END SUBROUTINE read_wannier_chk
  !
  !
END MODULE read_wannier
