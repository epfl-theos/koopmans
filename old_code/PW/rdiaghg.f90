!
! Copyright (C) 2003-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE rdiaghg( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H symmetric matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... LAPACK version - uses both DSYGV and DSYGVX
  !
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : me_pool, root_pool, intra_pool_comm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculated
    ! leading dimension of h, as declared in the calling pgm unit
  REAL(DP), INTENT(INOUT) :: h(ldh,n), s(ldh,n)
    ! matrix to be diagonalized
    ! overlap matrix
  !
  REAL(DP), INTENT(OUT) :: e(n)
    ! eigenvalues
  REAL(DP), INTENT(OUT) :: v(ldh,m)
    ! eigenvectors (column-wise)
  !
  INTEGER               :: i, j, lwork, nb, mm, info
    ! mm = number of calculated eigenvectors
  REAL(DP)              :: abstol
  REAL(DP), PARAMETER   :: one = 1_DP
  REAL(DP), PARAMETER   :: zero = 0_DP
  INTEGER,  ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP), ALLOCATABLE :: work(:), sdiag(:), hdiag(:)
  LOGICAL               :: all_eigenvalues
  INTEGER,  EXTERNAL    :: ILAENV
    ! ILAENV returns optimal block size "nb"
  !
  CALL start_clock( 'rdiaghg' )
  !
  ! ... only the first processor diagonalize the matrix
  !
  IF ( me_pool == root_pool ) THEN
     !
     ! ... save the diagonal of input S (it will be overwritten)
     !
     ALLOCATE( sdiag( n ) )
     DO i = 1, n
        sdiag(i) = s(i,i)
     END DO
     !
     all_eigenvalues = ( m == n )
     !
     ! ... check for optimal block size
     !
     nb = ILAENV( 1, 'DSYTRD', 'U', n, -1, -1, -1 )
     !
     IF ( nb < 1 .OR. nb >= n ) THEN
        !
        lwork = 8*n
        !
     ELSE
        !
        lwork = ( nb + 3 )*n
        !
     END IF
     !
     ALLOCATE( work( lwork ) )
     !
     IF ( all_eigenvalues ) THEN
        !
        ! ... calculate all eigenvalues
        !
        v(:,:) = h(:,:)
        !
#if defined (__ESSL)
        !
        ! ... there is a name conflict between essl and lapack ...
        !
        CALL DSYGV( 1, v, ldh, s, ldh, e, v, ldh, n, work, lwork )
        !
        info = 0
#else
        CALL DSYGV( 1, 'V', 'U', n, v, ldh, s, ldh, e, work, lwork, info )
#endif
        !
     ELSE
        !
        ! ... calculate only m lowest eigenvalues
        !
        ALLOCATE( iwork( 5*n ) )
        ALLOCATE( ifail( n ) )
        !
        ! ... save the diagonal of input H (it will be overwritten)
        !
        ALLOCATE( hdiag( n ) )
        DO i = 1, n
           hdiag(i) = h(i,i)
        END DO
        !
        abstol = 0.D0
       ! abstol = 2.D0*DLAMCH( 'S' )
        !
        CALL DSYGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                     0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                     work, lwork, iwork, ifail, info )
        !
        DEALLOCATE( ifail )
        DEALLOCATE( iwork )
        !
        ! ... restore input H matrix from saved diagonal and lower triangle
        !
        DO i = 1, n
           h(i,i) = hdiag(i)
           DO j = i + 1, n
              h(i,j) = h(j,i)
           END DO
           DO j = n + 1, ldh
              h(j,i) = 0.0_DP
           END DO
        END DO
        !
        DEALLOCATE( hdiag )
        !
     END IF
     !
     DEALLOCATE( work )
     !
     CALL errore( 'rdiaghg', 'diagonalization (DSYGV*) failed', ABS( info ) )
     
     ! ... restore input S matrix from saved diagonal and lower triangle
     !
     DO i = 1, n
        s(i,i) = sdiag(i)
        DO j = i + 1, n
           s(i,j) = s(j,i)
        END DO
        DO j = n + 1, ldh
           s(j,i) = 0.0_DP
        END DO
     END DO
     !
     DEALLOCATE( sdiag )
     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
  CALL mp_bcast( e, root_pool, intra_pool_comm )
  CALL mp_bcast( v, root_pool, intra_pool_comm )
  !
  CALL stop_clock( 'rdiaghg' )
  !
  RETURN
  !
END SUBROUTINE rdiaghg
!
!----------------------------------------------------------------------------
SUBROUTINE prdiaghg( n, h, s, ldh, e, v, desc )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H symmetric matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... Parallel version with ful data distribution
  !
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : root_pool, intra_pool_comm
  USE dspev_module,     ONLY : pdspev_drv
  USE descriptors,      ONLY : descla_siz_ , lambda_node_ , nlax_ , &
                               la_npc_ , la_npr_ , la_me_ , la_comm_ , &
                               nlar_ , la_myc_ , la_myr_
#if defined __SCALAPACK
  USE mp_global,        ONLY : ortho_cntx, me_blacs, np_ortho, me_ortho
  USE dspev_module,     ONLY : pdsyevd_drv
#endif
  !
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n, ldh
    ! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
    ! leading dimension of h, as declared in the calling pgm unit
  REAL(DP), INTENT(INOUT) :: h(ldh,ldh), s(ldh,ldh)
    ! matrix to be diagonalized
    ! overlap matrix
  !
  REAL(DP), INTENT(OUT) :: e(n)
    ! eigenvalues
  REAL(DP), INTENT(OUT) :: v(ldh,ldh)
    ! eigenvectors (column-wise)
  INTEGER, INTENT(IN)   :: desc( descla_siz_ )
  !
  INTEGER               :: nx
    ! local block size
  REAL(DP), PARAMETER   :: one = 1_DP
  REAL(DP), PARAMETER   :: zero = 0_DP
  REAL(DP), ALLOCATABLE :: hh(:,:)
  REAL(DP), ALLOCATABLE :: ss(:,:)
#ifdef __SCALAPACK
  INTEGER     :: desch( 16 ), info
#endif
  !
  CALL start_clock( 'rdiaghg' )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     nx   = desc( nlax_ )
     !
     IF( nx /= ldh ) &
        CALL errore(" prdiaghg ", " inconsistent leading dimension ", ldh )
     !
     ALLOCATE( hh( nx, nx ) )
     ALLOCATE( ss( nx, nx ) )
     !
     hh(1:nx,1:nx) = h(1:nx,1:nx)
     ss(1:nx,1:nx) = s(1:nx,1:nx)
     !
  END IF
  !
  CALL start_clock( 'rdiaghg:choldc' )
  !
  ! ... Cholesky decomposition of s ( L is stored in s )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
#ifdef __SCALAPACK
     CALL descinit( desch, n, n, desc( nlax_ ), desc( nlax_ ), 0, 0, ortho_cntx, SIZE( hh, 1 ) , info )
  
     IF( info /= 0 ) CALL errore( ' cdiaghg ', ' descinit ', ABS( info ) )
#endif
     !
#ifdef __SCALAPACK
     CALL PDPOTRF( 'L', n, ss, 1, 1, desch, info )
     IF( info /= 0 ) CALL errore( ' rdiaghg ', ' problems computing cholesky ', ABS( info ) )
#else
     CALL qe_pdpotrf( ss, nx, n, desc )
#endif
     !
  END IF
  !
  CALL stop_clock( 'rdiaghg:choldc' )
  !
  ! ... L is inverted ( s = L^-1 )
  !
  CALL start_clock( 'rdiaghg:inversion' )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
#ifdef __SCALAPACK
     ! 
     CALL sqr_dsetmat( 'U', n, zero, ss, size(ss,1), desc )

     CALL PDTRTRI( 'L', 'N', n, ss, 1, 1, desch, info )
     !
     IF( info /= 0 ) CALL errore( ' rdiaghg ', ' problems computing inverse ', ABS( info ) )
#else
     CALL qe_pdtrtri ( ss, nx, n, desc )
#endif
     !
  END IF
  !
  CALL stop_clock( 'rdiaghg:inversion' )
  !
  ! ... v = L^-1*H
  !
  CALL start_clock( 'rdiaghg:paragemm' )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'N', 'N', n, ONE, ss, nx, hh, nx, ZERO, v, nx, desc )
     !
  END IF
  !
  ! ... h = ( L^-1*H )*(L^-1)^T
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'N', 'T', n, ONE, v, nx, ss, nx, ZERO, hh, nx, desc )
     !
  END IF
  !
  CALL stop_clock( 'rdiaghg:paragemm' )
  !
  IF ( desc( lambda_node_ ) > 0 ) THEN
     ! 
     !  Compute local dimension of the cyclically distributed matrix
     !
#ifdef __SCALAPACK
     CALL pdsyevd_drv( .true., n, desc( nlax_ ), hh, SIZE(hh,1), e, ortho_cntx )
#else
     CALL qe_pdsyevd( .true., n, desc, hh, SIZE(hh,1), e )
#endif
     !
  END IF
  !
  ! ... v = (L^T)^-1 v
  !
  CALL start_clock( 'rdiaghg:paragemm' )
  !
  IF ( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'T', 'N', n, ONE, ss, nx, hh, nx, ZERO, v, nx, desc )
     !
     DEALLOCATE( ss )
     DEALLOCATE( hh )
     !
  END IF
  !
  CALL mp_bcast( e, root_pool, intra_pool_comm )
  !
  CALL stop_clock( 'rdiaghg:paragemm' )
  !
  CALL stop_clock( 'rdiaghg' )
  !
  RETURN
  !
END SUBROUTINE prdiaghg
!
!----------------------------------------------------------------------------
SUBROUTINE prdiaghg_nodist( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H symmetric matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... Parallel version, matrices are NOT distributed
  !
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : use_para_diag
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : me_pool, root_pool, intra_pool_comm, &
                               ortho_comm, np_ortho, me_ortho, ortho_comm_id
  USE parallel_toolkit, ONLY : dsqmdst, dsqmcll
  USE dspev_module,     ONLY : pdspev_drv
  USE descriptors,      ONLY : descla_siz_ , descla_init , lambda_node_ , &
                               nlax_ , la_nrl_ , la_npc_ , la_npr_ , la_me_ ,&
                               la_comm_ , la_nrlx_
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculated
    ! leading dimension of h, as declared in the calling pgm unit
  REAL(DP), INTENT(INOUT) :: h(ldh,n), s(ldh,n)
    ! matrix to be diagonalized
    ! overlap matrix
  !
  REAL(DP), INTENT(OUT) :: e(n)
    ! eigenvalues
  REAL(DP), INTENT(OUT) :: v(ldh,m)
    ! eigenvectors (column-wise)
  !
  INTEGER               :: i, j, lwork, nb, mm, info
    ! mm = number of calculated eigenvectors
  INTEGER               :: nx, nrl, nrlx
    ! local block size
  REAL(DP)              :: abstol
  REAL(DP), PARAMETER   :: one = 1_DP
  REAL(DP), PARAMETER   :: zero = 0_DP
  INTEGER,  ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP), ALLOCATABLE :: sl(:,:), hl(:,:), vl(:,:), vv(:,:), diag(:,:)
  REAL(DP), ALLOCATABLE :: work(:), sdiag(:), hdiag(:)
  LOGICAL               :: all_eigenvalues
  INTEGER,  EXTERNAL    :: ILAENV
    ! ILAENV returns optimal block size "nb"
  INTEGER               :: desc( descla_siz_ )
  !
  CALL start_clock( 'rdiaghg' )
  !
  CALL descla_init( desc, n, n, np_ortho, me_ortho, ortho_comm, ortho_comm_id )
  !
  nx  = desc( nlax_ )
  nrl  = desc( la_nrl_ )
  nrlx = desc( la_nrlx_ )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     ALLOCATE( sl( nx , nx ) )
     ALLOCATE( vl( nx , nx ) )
     ALLOCATE( hl( nx , nx ) )
  END IF

  !  distribute matrixes s and h
  !
  CALL dsqmdst( n, s, ldh, sl, nx, desc )
  CALL dsqmdst( n, h, ldh, hl, nx, desc )
  !
  !  Call block parallel algorithm
  !
  CALL start_clock( 'rdiaghg:choldc' )
  !
  ! ... Cholesky decomposition of s ( L is stored in s )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL qe_pdpotrf( sl, nx, n, desc )
     !
  END IF
  !
  CALL stop_clock( 'rdiaghg:choldc' )
  !
  ! ... L is inverted ( s = L^-1 )
  !
  CALL start_clock( 'rdiaghg:inversion' )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL qe_pdtrtri ( sl, nx, n, desc )
     !
  END IF
  !
  CALL stop_clock( 'rdiaghg:inversion' )
  !
  ! ... v = L^-1*H
  !
  CALL start_clock( 'rdiaghg:paragemm' )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'N', 'N', n, ONE, sl, nx, hl, nx, ZERO, vl, nx, desc )
     !
  END IF
  !
  ! ... h = ( L^-1*H )*(L^-1)^T
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'N', 'T', n, ONE, vl, nx, sl, nx, ZERO, hl, nx, desc )
     !
  END IF
  !
  CALL stop_clock( 'rdiaghg:paragemm' )
  !
  IF ( desc( lambda_node_ ) > 0 ) THEN
     ! 
     CALL qe_pdsyevd( .true., n, desc, hl, SIZE(hl,1), e )
     !
  END IF
  !
  ! ... v = (L^T)^-1 v
  !
  CALL start_clock( 'rdiaghg:paragemm' )
  !
  IF ( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'T', 'N', n, ONE, sl, nx, hl, nx, ZERO, vl, nx, desc )
     !
  END IF
  !
  CALL mp_bcast( e, root_pool, intra_pool_comm )
  !
  CALL stop_clock( 'rdiaghg:paragemm' )
  !
  ! CALL prdiaghg( n, hl, sl, nx, e, vl, desc )
  !
  !  collect distributed matrix vl into replicated matrix v
  !
  CALL dsqmcll( n, vl, nx, v, ldh, desc, intra_pool_comm )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     DEALLOCATE( hl, sl, vl )
  END IF
  !
  RETURN
  !
END SUBROUTINE prdiaghg_nodist
