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
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: read_u_matrices
  !
  CONTAINS
  !
  !---------------------------------------------------------------------
  SUBROUTINE read_u_matrices( u_mat, u_mat_opt )
    !-------------------------------------------------------------------
    !
    ! ...  Read from the Wannier90 chk file the matrices U(k) and
    ! ...  U_opt(k)
    !
    USE kinds,               ONLY : DP
    USE io_global,           ONLY : ionode, ionode_id
    USE mp_global,           ONLY : intra_image_comm
    USE mp,                  ONLY : mp_bcast
    USE wannier,             ONLY : iknum, seedname, mp_grid, &
                                    n_wannier, num_bands
    !
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(OUT) :: u_mat(:,:,:)
    COMPLEX(DP), INTENT(OUT) :: u_mat_opt(:,:,:)
    !
    INTEGER :: chk_unit=124
    INTEGER :: i, j, nkp, nn
    LOGICAL :: exst
    CHARACTER(LEN=33) :: header
    INTEGER :: nbands, num_exclude_bands
    LOGICAL, ALLOCATABLE :: exclude_bands(:)
    REAL(DP) :: at(3,3), bg(3,3)
    INTEGER :: num_kpts
    INTEGER :: kgrid(3)
    REAL(DP), ALLOCATABLE :: kpt_latt(:,:)
    INTEGER :: nntot, num_wann
    LOGICAL :: checkpoint, have_disentangled
    REAL(DP) :: omega_invariant           ! disentanglement parameters
    LOGICAL, ALLOCATABLE :: lwindow(:,:)  ! disentanglement parameters
    INTEGER, ALLOCATABLE :: ndimwin(:)    ! disentanglement parameters
    COMPLEX(DP), ALLOCATABLE :: m_mat(:,:,:,:)
    REAL(DP), ALLOCATABLE :: centers(:,:)
    REAL(DP), ALLOCATABLE :: spreads(:)
    !
    !
    IF ( ionode ) THEN
      !
      INQUIRE( FILE=trim(seedname)//'.chk', EXIST=exst )
      IF ( .not. exst ) CALL errore( 'read_u_matrix', 'chk file not found', 1 )
      !
      OPEN( UNIT=chk_unit, FILE=trim(seedname)//'.chk', STATUS='old', FORM='unformatted' )
      !
      READ( chk_unit ) header                   ! date and time
      READ( chk_unit ) nbands                   ! number of bands
      READ( chk_unit ) num_exclude_bands        ! number of excluded bands
      !
      IF ( nbands .ne. num_bands ) &
        CALL  errore( 'read_u_matrix', 'Invalid value for nbands', nbands )
      ! 
      IF ( num_exclude_bands .lt. 0 ) &
        CALL  errore( 'read_u_matrix', 'Invalid value for num_exclude_bands', &
                      num_exclude_bands )
      !
      ALLOCATE( exclude_bands(num_exclude_bands) )
      exclude_bands(:) = .FALSE.
      !
      READ( chk_unit ) ( exclude_bands(i), i=1,num_exclude_bands )   ! excluded bands
      READ( chk_unit ) (( at(i,j), i=1,3 ), j=1,3 )                  ! prim real latt vectors
      READ( chk_unit ) (( bg(i,j), i=1,3 ), j=1,3 )                  ! prim recip latt vectors
      READ( chk_unit ) num_kpts                                      ! num of k-points
      !
      IF ( num_kpts .ne. iknum ) &
        CALL errore( 'read_u_matrix', 'Invalid value for num_kpts', num_kpts )
      !
      READ( chk_unit ) ( kgrid(i), i=1,3 )                           ! MP grid
      !
      IF ( ANY( kgrid .ne. mp_grid ) ) &
        CALL errore( 'read_u_matrix', 'Invalid value for kgrid', 1 )
      !
      ALLOCATE( kpt_latt(3,num_kpts) )
      !
      READ( chk_unit ) (( kpt_latt(i,nkp), i=1,3 ), nkp=1,num_kpts )
      READ( chk_unit ) nntot                                         ! nntot
      READ( chk_unit ) num_wann                                      ! num of WFs
      READ( chk_unit ) checkpoint                                    ! checkpoint
      READ( chk_unit ) have_disentangled     ! .true. if disentanglement has been performed
      !
      IF ( num_wann .ne. n_wannier ) &
        CALL errore( 'read_u_matrix', 'Invalid value for num_wann', num_wann )
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
        READ ( chk_unit ) ((( u_mat_opt(i,j,nkp), i=1,num_wann ), &  ! optimal U-matrix
                                                  j=1,num_wann ), &
                                                  nkp=1,num_kpts )
        !
      ENDIF
      !
      READ ( chk_unit ) ((( u_mat(i,j,nkp), i=1,num_wann ), &        ! U-matrix
                                            j=1,num_wann ), &
                                            nkp=1,num_kpts )
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
    CALL mp_bcast( u_mat, ionode_id, intra_image_comm )
    CALL mp_bcast( u_mat_opt, ionode_id, intra_image_comm )
    !
    !
  END SUBROUTINE read_u_matrices
  !
  !
END MODULE read_wannier
