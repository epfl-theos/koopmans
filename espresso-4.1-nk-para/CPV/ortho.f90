!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_m( c0, cp, lambda, descla, ccc, nupdwn, iupdwn, nspin )
      !
      USE kinds,              ONLY: DP
      USE control_flags,      ONLY: force_pairing
      USE cp_main_variables,  ONLY: ema0bg
      USE descriptors,        ONLY: lambda_node_ , nlar_ , nlac_
      USE control_flags,      ONLY: ortho_eps, ortho_max
      USE orthogonalize_base, ONLY: calphi, updatc
      USE cp_interfaces,      ONLY: ortho_gamma
      !
      IMPLICIT NONE

      INTEGER,     INTENT(IN)    :: descla(:,:)
      INTEGER,     INTENT(IN)    :: nupdwn(:), iupdwn(:), nspin
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:), cp(:,:)
      REAL(DP),    INTENT(INOUT) :: lambda(:,:,:)
      REAL(DP),    INTENT(IN)    :: ccc
      !
      COMPLEX(DP), ALLOCATABLE :: phi(:,:)
      INTEGER                  :: iss, nsc, iwfc, nwfc, info
      INTEGER                  :: iter, i, j
      INTEGER                  :: ngwx, n, nr, nc, nx
      REAL(DP)                 :: diff
      REAL(DP),    ALLOCATABLE :: dum(:,:)
      REAL(DP),    ALLOCATABLE :: ddum(:,:)
      COMPLEX(DP), ALLOCATABLE :: cdum(:,:)
      !
      CALL start_clock( 'ortho' )  

      n    = SIZE( c0, 2 )
      ngwx = SIZE( c0, 1 )
      nx   = SIZE( lambda, 1 )

      ALLOCATE( dum( 1, n ) )
      ALLOCATE( ddum( 1, nx ) )
      ALLOCATE( cdum( ngwx, 1 ) )

      ALLOCATE( phi( ngwx, n ), STAT = info )
      IF( info /= 0 ) CALL errore( ' ortho ', ' allocating phi ', 3 )

      CALL calphi( c0, ngwx, dum, 1, cdum, phi, n, ema0bg )
      !
      nsc = nspin
      IF( force_pairing ) nsc = 1
      !
      DO iss = 1, nsc
          !
          nwfc = nupdwn(iss)
          iwfc = iupdwn(iss)
          !
          CALL ortho_gamma( 1, cp, ngwx, phi, dum, ddum, 1, dum, ddum, lambda(:,:,iss), nx, &
               descla(:,iss), diff, iter, n, nwfc, iwfc )
          !
          IF ( iter > ortho_max ) THEN
             call errore(' ortho ','  itermax ',iter)
          END IF
          !
          CALL updatc( 1.0d0, n, lambda(:,:,iss), nx, phi, ngwx, dum, 1, dum, dum, cp, &
                       nwfc, iwfc, descla(:,iss) )
          !     
          !     lagrange multipliers
          !
          IF( descla( lambda_node_ , iss ) > 0 ) THEN
             DO j = 1, descla( nlac_ , iss )
                DO i = 1, descla( nlar_ , iss )
                   lambda( i, j, iss ) = lambda( i, j, iss ) / ccc
                END DO
             END DO
          END IF
          !
      END DO
      !
      IF( force_pairing ) cp(:, iupdwn(2):iupdwn(2)+nupdwn(2)-1 ) = cp(:,1:nupdwn(2))
      !
      DEALLOCATE( phi )
      DEALLOCATE( dum )
      DEALLOCATE( ddum )
      DEALLOCATE( cdum )
      !
      CALL stop_clock( 'ortho' )
      !
      RETURN
   END SUBROUTINE ortho_m

SUBROUTINE ortho_m_twin(c0, cp, lambda, descla, ccc, nupdwn, iupdwn, nspin)
      !
      USE kinds,              ONLY: DP
      USE control_flags,      ONLY: force_pairing
      USE cp_main_variables,  ONLY: ema0bg
      USE descriptors,        ONLY: lambda_node_ , nlar_ , nlac_
      USE control_flags,      ONLY: ortho_eps, ortho_max
      USE orthogonalize_base, ONLY: calphi, updatc
      USE cp_interfaces,      ONLY: ortho_gamma
      USE twin_types
      USE mp_global,          ONLY: mpime
      !
      IMPLICIT NONE

      INTEGER,     INTENT(IN)    :: descla(:,:)
      INTEGER,     INTENT(IN)    :: nupdwn(:), iupdwn(:), nspin
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:), cp(:,:)
      TYPE(twin_matrix), DIMENSION(:), INTENT(INOUT) :: lambda
      REAL(DP),    INTENT(IN)    :: ccc
      !
      COMPLEX(DP), ALLOCATABLE :: phi(:,:)
      INTEGER                  :: iss, nsc, iwfc, nwfc, info
      INTEGER                  :: iter, i, j
      INTEGER                  :: ngwx, n, nr, nc, nx
      REAL(DP)                 :: diff
      REAL(DP),    ALLOCATABLE :: dum(:,:)
      REAL(DP),    ALLOCATABLE :: ddum(:,:)
      COMPLEX(DP), ALLOCATABLE :: cdum(:,:)
      COMPLEX(DP), ALLOCATABLE :: dum_c(:,:)
      COMPLEX(DP), ALLOCATABLE :: ddum_c(:,:)
      !
      CALL start_clock( 'ortho' )  

      n    = SIZE( c0, 2 )
      ngwx = SIZE( c0, 1 )
      nx   = SIZE( lambda, 1 )

      ALLOCATE( dum( 1, n ) )
      ALLOCATE( ddum( 1, nx ) )
      ALLOCATE( cdum( ngwx, 1 ) )
      ALLOCATE( dum_c( 1, n ) )
      ALLOCATE( ddum_c( 1, nx ) )

      ALLOCATE( phi( ngwx, n ), STAT = info )
      IF( info /= 0 ) CALL errore( ' ortho ', ' allocating phi ', 3 )

      CALL calphi( c0, ngwx, dum, 1, cdum, phi, n, ema0bg )
      !
      nsc = nspin
      IF( force_pairing ) nsc = 1
      !
      DO iss = 1, nsc
          !
          nwfc = nupdwn(iss)
          iwfc = iupdwn(iss)
          !
          IF(.not.lambda(iss)%iscmplx) THEN
	      CALL ortho_gamma( 1, cp, ngwx, phi, dum, ddum, 1, dum, ddum, lambda(iss)%rvec(:,:), nx, &
		  descla(:,iss), diff, iter, n, nwfc, iwfc )
          ELSE
              !write(313+mpime,*) "lambda_ref", lambda(iss)%cvec
	      CALL ortho_gamma( 1, cp, ngwx, phi, dum_c, ddum_c, 1, dum_c, ddum_c, lambda(iss)%cvec(:,:), nx, &
		  descla(:,iss), diff, iter, n, nwfc, iwfc )
          ENDIF
          !
          IF ( iter > ortho_max ) THEN
             call errore(' ortho ','  itermax ',iter)
          END IF
          !
          IF(.not.lambda(iss)%iscmplx) THEN
	      CALL updatc( 1.0d0, n, lambda(iss)%rvec(:,:), nx, phi, ngwx, dum, 1, dum, dum, cp, &
			  nwfc, iwfc, descla(:,iss) )
          ELSE
	      CALL updatc( 1.0d0, n, lambda(iss)%cvec(:,:), nx, phi, ngwx, cdum, 1, cdum, cdum, cp, &
			  nwfc, iwfc, descla(:,iss))
          ENDIF
          !     
          !     lagrange multipliers
          !
          IF( descla( lambda_node_ , iss ) > 0 ) THEN
             IF(.not.lambda(iss)%iscmplx) THEN
		DO j = 1, descla( nlac_ , iss )
		    DO i = 1, descla( nlar_ , iss )
		      lambda(iss)%rvec( i, j ) = lambda(iss)%rvec( i, j ) / ccc
		    END DO
		END DO
             ELSE
		DO j = 1, descla( nlac_ , iss )
		    DO i = 1, descla( nlar_ , iss )
		      lambda(iss)%cvec( i, j) = lambda(iss)%cvec( i, j) / ccc
		    END DO
		END DO
             ENDIF
          END IF
          !
      END DO
      !
      IF( force_pairing ) cp(:, iupdwn(2):iupdwn(2)+nupdwn(2)-1 ) = cp(:,1:nupdwn(2))
      !
      DEALLOCATE( phi )
      DEALLOCATE( dum )
      DEALLOCATE( ddum )
      DEALLOCATE( cdum )
      !
      CALL stop_clock( 'ortho' )
      !
      RETURN
   END SUBROUTINE ortho_m_twin

!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_gamma_real_x( iopt, cp, ngwx, phi, becp, qbecp, nkbx, bephi, qbephi, &
                           x0, nx0, descla, diff, iter, n, nss, istart )
!=----------------------------------------------------------------------------=!
      !
      USE kinds,              ONLY: DP
      USE orthogonalize_base, ONLY: rhoset, sigset, tauset, ortho_iterate,   &
                                    ortho_alt_iterate, diagonalize_serial,   &
                                    use_parallel_diag, diagonalize_parallel
      USE descriptors,        ONLY: lambda_node_ , nlar_ , nlac_ , ilar_ , ilac_ , &
                                    nlax_ , la_comm_ , descla_siz_
      USE mp_global,          ONLY: nproc_image, me_image, intra_image_comm
      USE mp,                 ONLY: mp_sum
      USE parallel_toolkit,            ONLY: sqr_tr_cannon

      IMPLICIT  NONE

      ! ... Arguments

      INTEGER,  INTENT(IN)  :: iopt
      INTEGER,  INTENT(IN)  :: ngwx, nkbx, nx0
      INTEGER,  INTENT(IN)  :: n, nss, istart
      COMPLEX(DP) :: phi( ngwx, n ), cp( ngwx, n )
      REAL(DP)    :: bephi( :, : ), becp( :, : )
      REAL(DP)    :: qbephi( :, : ), qbecp( :, : )
      REAL(DP)    :: x0( nx0, nx0 )
      INTEGER,  INTENT(IN)  :: descla( descla_siz_ )
      INTEGER,  INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff

      ! ... Locals

      REAL(DP),   ALLOCATABLE :: s(:,:), sig(:,:), tau(:,:), rhot(:,:)
      REAL(DP),   ALLOCATABLE :: wrk(:,:), rhoa(:,:), rhos(:,:), rhod(:)
      INTEGER  :: i, j, info, nr, nc, ir, ic
      INTEGER  :: nlam, nlax
      LOGICAL  :: la_proc

      !
      ! ...   Subroutine body
      !
      nlax    = descla( nlax_ )
      la_proc = ( descla( lambda_node_ ) > 0 )
      nlam    = 1
      if ( la_proc ) nlam = nlax

      IF( la_proc ) THEN
         !
         IF( nx0 /= descla( nlax_ ) ) &
            CALL errore( ' ortho_gamma ', ' inconsistent dimensions nx0 ' , nx0 )
         IF( nlam /= descla( nlax_ ) ) &
            CALL errore( ' ortho_gamma ', ' inconsistent dimensions nlam ' , nlam )
         !
         nr = descla( nlar_ )
         nc = descla( nlac_ )
         !
         ir = descla( ilar_ )
         ic = descla( ilac_ )
         !
      ELSE
         !
         nr = 1
         nc = 1
         !
         IF( nlam /= 1 ) &
            CALL errore( ' ortho_gamma ', ' inconsistent dimensions nlam, should be 1 ' , nlam )
         IF( nx0 /= 1 ) &
            CALL errore( ' ortho_gamma ', ' inconsistent dimensions nx0, should be 1 ' , nx0 )
         !
      END IF
      !
      ALLOCATE( rhos( nlam, nlam ) )
      ALLOCATE( rhoa( nlam, nlam ) )   !   antisymmetric part of rho
      ALLOCATE( s( nlam, nlam ) ) 
      ALLOCATE( sig( nlam, nlam ) ) 
      ALLOCATE( tau( nlam, nlam ) ) 
      !
      ALLOCATE( rhod( nss ) )
      !
      !     rho = <s'c0|s|cp>
      !
      CALL start_clock( 'rhoset' )
      !
      CALL rhoset( cp, ngwx, phi, bephi, nkbx, qbecp, n, nss, istart, rhos, nlam, descla )
      !
      IF( la_proc ) THEN
         !
         ALLOCATE( rhot( nlam, nlam ) )   !   transpose of rho
         !
         !    distributed array rhos contains "rho", 
         !    now transpose rhos and store the result in distributed array rhot
         !
         CALL sqr_tr_cannon( nss, rhos, nlam, rhot, nlam, descla ) !modified:giovanni
         !
         !  Compute the symmetric part of rho
         !
         DO j = 1, nc
            DO i = 1, nr
               rhos( i, j ) = 0.5d0 * ( rhos( i, j ) + rhot( i, j ) )
            END DO
         END DO
         !
         !    distributed array rhos now contains symmetric part of "rho", 
         !
         CALL consistency_check( rhos )
         !
         !  Antisymmetric part of rho, alredy distributed across ortho procs.
         !
         DO j = 1, nc
            DO i = 1, nr
               rhoa( i, j ) = rhos( i, j ) - rhot( i, j )
            END DO
         END DO
         !
         DEALLOCATE( rhot )
         !
      END IF

      CALL stop_clock( 'rhoset' )


      CALL start_clock( 'rsg' )
      !
      ! ...   Diagonalize symmetric part of rho (rhos)
      ! ...   "s" is the matrix of eigenvectors, "rhod" is the array of eigenvalues
      !
      IF( use_parallel_diag ) THEN
         !
         CALL diagonalize_parallel( nss, rhos, rhod, s, descla )
         !
      ELSE
         !
         IF( la_proc ) THEN
            !
            ALLOCATE( wrk( nss, nss ), STAT = info )
            IF( info /= 0 ) CALL errore( ' ortho ', ' allocating matrixes ', 1 )
            !
            CALL collect_matrix( wrk, rhos )
            !
            CALL diagonalize_serial( nss, wrk, rhod )
            !
            CALL distribute_matrix( wrk, s )
            !
            DEALLOCATE( wrk )
            !
         END IF
         !
      END IF
      !
      CALL stop_clock( 'rsg' )
      !
      !     sig = 1-<cp|s|cp>
      !
      CALL sigset( cp, ngwx, becp, nkbx, qbecp, n, nss, istart, sig, nlam, descla )
      !
      !     tau = <s'c0|s|s'c0>
      !
      CALL tauset( phi, ngwx, bephi, nkbx, qbephi, n, nss, istart, tau, nlam, descla )
      !
      CALL start_clock( 'ortho_iter' )
      !
      IF( iopt == 0 ) THEN
         !
         CALL ortho_iterate( iter, diff, s, nlam, rhod, x0, nx0, sig, rhoa, rhos, tau, nss, descla)
         !
      ELSE
         !
         CALL ortho_alt_iterate( iter, diff, s, nlam, rhod, x0, nx0, sig, rhoa, tau, nss, descla)
         !
      END IF
      !
      CALL stop_clock( 'ortho_iter' )
      !
      DEALLOCATE( rhoa, rhos, rhod, s, sig, tau )
      !
      IF( la_proc )  CALL consistency_check( x0 )

      RETURN

   CONTAINS

      SUBROUTINE distribute_matrix( a, b )
         REAL(DP) :: a(:,:), b(:,:)
         INTEGER :: i, j
         IF( la_proc ) THEN
            DO j = 1, nc
               DO i = 1, nr
                  b( i, j ) = a( i + ir - 1, j + ic - 1 )
               END DO
            END DO
         END IF
         RETURN
      END SUBROUTINE

      SUBROUTINE collect_matrix( a, b )
         REAL(DP) :: a(:,:), b(:,:)
         INTEGER :: i, j
         a = 0.0d0
         IF( la_proc ) THEN
            DO j = 1, nc
               DO i = 1, nr
                  a( ir + i - 1, ic + j - 1 ) = b( i, j )
               END DO
            END DO
         END IF
         CALL mp_sum( a, descla( la_comm_ ) )
         RETURN
      END SUBROUTINE

      SUBROUTINE consistency_check( a )
         REAL(DP) :: a(:,:)
         INTEGER :: i, j
         !
         ! on some machines (IBM RS/6000 for instance) the following test allows
         ! to distinguish between Numbers and Sodium Nitride (NaN, Not a Number).
         ! If a matrix of Not-Numbers is passed to rs, the most likely outcome is
         ! that the program goes on forever doing nothing and writing nothing.
         !
         DO j = 1, nc
            DO i = 1, nr
               IF (a(i,j) /= a(i,j)) CALL errore(' ortho ',' ortho went bananas ',1)
            END DO
         END DO
         RETURN
      END SUBROUTINE

   END SUBROUTINE ortho_gamma_real_x

!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_gamma_cmplx_x( iopt, cp, ngwx, phi, becp, qbecp, nkbx, bephi, qbephi, &
                           x0, nx0, descla, diff, iter, n, nss, istart )
!=----------------------------------------------------------------------------=!
      !
      USE kinds,              ONLY: DP
      USE orthogonalize_base, ONLY: rhoset, sigset, tauset, ortho_iterate,   &
                                    ortho_alt_iterate, diagonalize_serial,   &
                                    use_parallel_diag, diagonalize_parallel
      USE descriptors,        ONLY: lambda_node_ , nlar_ , nlac_ , ilar_ , ilac_ , &
                                    nlax_ , la_comm_ , descla_siz_
      USE mp_global,          ONLY: nproc_image, me_image, intra_image_comm
      USE mp,                 ONLY: mp_sum
      USE parallel_toolkit,      ONLY: sqr_tr_cannon

      IMPLICIT  NONE

      ! ... Arguments

      INTEGER,  INTENT(IN)  :: iopt
      INTEGER,  INTENT(IN)  :: ngwx, nkbx, nx0
      INTEGER,  INTENT(IN)  :: n, nss, istart
      COMPLEX(DP) :: phi( ngwx, n ), cp( ngwx, n )
      COMPLEX(DP)    :: bephi( :, : ), becp( :, : )
      COMPLEX(DP)    :: qbephi( :, : ), qbecp( :, : )
      COMPLEX(DP)    :: x0( nx0, nx0 )
      INTEGER,  INTENT(IN)  :: descla( descla_siz_ )
      INTEGER,  INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff

      ! ... Locals

      COMPLEX(DP),   ALLOCATABLE :: s(:,:), sig(:,:), tau(:,:), rhot(:,:)
      COMPLEX(DP),   ALLOCATABLE :: wrk(:,:), rhoa(:,:), rhos(:,:)
      REAL(DP), ALLOCATABLE :: rhod(:)
      INTEGER  :: i, j, info, nr, nc, ir, ic
      INTEGER  :: nlam, nlax
      LOGICAL  :: la_proc

!       write(6,*) "calling ortho_gamma"!added:giovanni:debug
!       do i=1, size(x0,1)
! 	do j=1, size(x0,2)
!            write(6,*) i,j,x0(i,j)
! 	enddo
!       enddo
      !
      ! ...   Subroutine body
      !
      nlax    = descla( nlax_ )
      la_proc = ( descla( lambda_node_ ) > 0 )
      nlam    = 1
      if ( la_proc ) nlam = nlax

      IF( la_proc ) THEN
         !
         IF( nx0 /= descla( nlax_ ) ) &
            CALL errore( ' ortho_gamma ', ' inconsistent dimensions nx0 ' , nx0 )
         IF( nlam /= descla( nlax_ ) ) &
            CALL errore( ' ortho_gamma ', ' inconsistent dimensions nlam ' , nlam )
         !
         nr = descla( nlar_ )
         nc = descla( nlac_ )
         !
         ir = descla( ilar_ )
         ic = descla( ilac_ )
         !
      ELSE
         !
         nr = 1
         nc = 1
         !
         IF( nlam /= 1 ) &
            CALL errore( ' ortho_gamma ', ' inconsistent dimensions nlam, should be 1 ' , nlam )
         IF( nx0 /= 1 ) &
            CALL errore( ' ortho_gamma ', ' inconsistent dimensions nx0, should be 1 ' , nx0 )
         !
      END IF
      !
      ALLOCATE( rhos( nlam, nlam ) )
      ALLOCATE( rhoa( nlam, nlam ) )   !   antisymmetric part of rho
      ALLOCATE( s( nlam, nlam ) ) 
      ALLOCATE( sig( nlam, nlam ) ) 
      ALLOCATE( tau( nlam, nlam ) ) 
      !
      ALLOCATE( rhod( nss ) )
      !
      !     rho = <s'c0|s|cp>
      !
      CALL start_clock( 'rhoset' )
      !
!       write(6,*) "becp", becp !added:giovanni:debug
!       write(6,*) "bephi", bephi, "qbecp", qbecp !added:giovanni:debug
      CALL rhoset( cp, ngwx, phi, bephi, nkbx, qbecp, n, nss, istart, rhos, nlam, descla )
!       write(6,*) "rhos", rhos !added:giovanni:debug
      !
      IF( la_proc ) THEN
         !
         ALLOCATE( rhot( nlam, nlam ) )   !   transpose of rho
         !
         !    distributed array rhos contains "rho", 
         !    now transpose rhos and store the result in distributed array rhot
         !
         CALL sqr_tr_cannon( nss, rhos, nlam, rhot, nlam, descla )
         !
         !  Compute the symmetric part of rho
         !
         DO j = 1, nc
            DO i = 1, nr
               rhos( i, j ) = 0.5d0 * ( rhos( i, j ) + rhot( i, j ) )
            END DO
         END DO
         !
         !    distributed array rhos now contains symmetric part of "rho", 
         !
         CALL consistency_check( rhos )
         !
         !  Antisymmetric part of rho, alredy distributed across ortho procs.
         !
         DO j = 1, nc
            DO i = 1, nr
               rhoa( i, j ) = rhos( i, j ) - rhot( i, j )
            END DO
         END DO
         !
         DEALLOCATE( rhot )
         !
      END IF
!      write(6,*) "rhoss", rhos !added:giovanni:debug
!      write(6,*) "rhoa", rhoa !added:giovanni:debug

      CALL stop_clock( 'rhoset' )


      CALL start_clock( 'rsg' )
      !
      ! ...   Diagonalize symmetric part of rho (rhos)
      ! ...   "s" is the matrix of eigenvectors, "rhod" is the array of eigenvalues
      !
      IF( use_parallel_diag ) THEN
         !
!          write(6,*) "calling diagonalize_parallel" ! added:giovanni:debug
         CALL diagonalize_parallel( nss, rhos, rhod, s, descla )
!      write(6,*) "rhod", rhod !added:giovanni:debug
         !
      ELSE
         !
         IF( la_proc ) THEN
            !
            ALLOCATE( wrk( nss, nss ), STAT = info )
            IF( info /= 0 ) CALL errore( ' ortho ', ' allocating matrixes ', 1 )
            !
            CALL collect_matrix( wrk, rhos )
            !
            CALL diagonalize_serial( nss, wrk, rhod )
            !
            CALL distribute_matrix( wrk, s )
            !
            DEALLOCATE( wrk )
            !
         END IF
         !
      END IF
      !
      CALL stop_clock( 'rsg' )
      !
      !     sig = 1-<cp|s|cp>
      !
      CALL sigset( cp, ngwx, becp, nkbx, qbecp, n, nss, istart, sig, nlam, descla )
!      write(6,*) "sig", sig !added:giovanni:debug
      !
      !     tau = <s'c0|s|s'c0>
      !
      CALL tauset( phi, ngwx, bephi, nkbx, qbephi, n, nss, istart, tau, nlam, descla )
!      write(6,*)"tau", tau!, "qbephi2", qbephi !added:giovanni:debug
      !
!       write(6,*) "calling ortho_iterate", x0!added:giovanni:debug
      !stop
      CALL start_clock( 'ortho_iter' )
      !
      IF( iopt == 0 ) THEN
         !
         CALL ortho_iterate( iter, diff, s, nlam, rhod, x0, nx0, sig, rhoa, rhos, tau, nss, descla)
         !
      ELSE
         !
         CALL ortho_alt_iterate( iter, diff, s, nlam, rhod, x0, nx0, sig, rhoa, tau, nss, descla)
         !
      END IF
      !
      CALL stop_clock( 'ortho_iter' )
      !
      DEALLOCATE( rhoa, rhos, rhod, s, sig, tau )
      !
      IF( la_proc )  CALL consistency_check( x0 )

      RETURN

   CONTAINS

      SUBROUTINE distribute_matrix( a, b )
         COMPLEX(DP) :: a(:,:), b(:,:)
         INTEGER :: i, j
         IF( la_proc ) THEN
            DO j = 1, nc
               DO i = 1, nr
                  b( i, j ) = a( i + ir - 1, j + ic - 1 )
               END DO
            END DO
         END IF
         RETURN
      END SUBROUTINE

      SUBROUTINE collect_matrix( a, b )
         COMPLEX(DP) :: a(:,:), b(:,:)
         INTEGER :: i, j
         a = CMPLX(0.0d0, 0.d0)
         IF( la_proc ) THEN
            DO j = 1, nc
               DO i = 1, nr
                  a( ir + i - 1, ic + j - 1 ) = b( i, j )
               END DO
            END DO
         END IF
         CALL mp_sum( a, descla( la_comm_ ) )
         RETURN
      END SUBROUTINE

      SUBROUTINE consistency_check( a )
         COMPLEX(DP) :: a(:,:)
         INTEGER :: i, j
         !
         ! on some machines (IBM RS/6000 for instance) the following test allows
         ! to distinguish between Numbers and Sodium Nitride (NaN, Not a Number).
         ! If a matrix of Not-Numbers is passed to rs, the most likely outcome is
         ! that the program goes on forever doing nothing and writing nothing.
         !
         DO j = 1, nc
            DO i = 1, nr
               IF (a(i,j) /= a(i,j)) CALL errore(' ortho ',' ortho went bananas ',1)
            END DO
         END DO
         RETURN
      END SUBROUTINE

   END SUBROUTINE ortho_gamma_cmplx_x


!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_cp_real( eigr, cp, phi, ngwx, x0, descla, diff, iter, ccc, &
                        bephi, becp, nbsp, nspin, nupdwn, iupdwn )
!=----------------------------------------------------------------------------=!
      !
      !     input = cp (non-orthonormal), beta
      !     input = phi |phi>=s'|c0>
      !     output= cp (orthonormal with s( r(t+dt) ) )
      !     output= bephi, becp
      !     the method used is similar to the version in les houches 1988
      !     'simple molecular systems at..'  p. 462-463  (18-22)
      !      xcx + b x + b^t x^t + a = 1
      !     where c = <s'c0|s|s'c0>   b = <s'c0|s cp>   a = <cp|s|cp>
      !     where s=s(r(t+dt)) and s'=s(r(t))  
      !     for vanderbilt pseudo pot - kl & ap
      !
      USE kinds,          ONLY: DP
      USE ions_base,      ONLY: na, nat
      USE cvan,           ONLY: ish, nvb
      USE uspp,           ONLY: nkb, qq
      USE uspp_param,     ONLY: nh
      USE electrons_base, ONLY: f
      USE gvecw,          ONLY: ngw
      USE control_flags,  ONLY: iprint, iprsta, ortho_max
      USE control_flags,  ONLY: force_pairing
      USE io_global,      ONLY: stdout, ionode
      USE cp_interfaces,  ONLY: ortho_gamma, nlsm1, nlsm1_dist
      USE descriptors,    ONLY: nlac_ , ilac_ , descla_siz_ , nlar_ , ilar_, lambda_node_, nlax_
      USE mp_global,          ONLY: nproc_image, me_image, intra_image_comm  ! DEBUG
      !
      IMPLICIT NONE
      !
      INTEGER,    INTENT(IN) :: ngwx, nbsp, nspin
      INTEGER,    INTENT(IN) :: nupdwn( nspin ), iupdwn( nspin )
      INTEGER,    INTENT(IN) :: descla(descla_siz_,nspin)
      COMPLEX(DP) :: cp(ngwx,nbsp), phi(ngwx,nbsp), eigr(ngwx,nat)
      REAL(DP)    :: x0(:,:,:), diff, ccc
      INTEGER     :: iter
      REAL(DP)    :: bephi(:,:), becp(:,:)
      !
      REAL(DP), ALLOCATABLE :: xloc(:,:)
      REAL(DP), ALLOCATABLE :: qbephi(:,:,:), qbecp(:,:,:), bephi_c(:,:)
      !
      INTEGER :: nlam, nlax
      LOGICAL :: la_proc
      !
      INTEGER :: nkbx
      INTEGER :: istart, nss, ifail, i, j, iss, iv, jv, ia, is, inl, jnl
      INTEGER :: nspin_sub, nx0, nc, ic, icc, nr, ir
      REAL(DP) :: qqf
      !
!       write(6,*) "inside ortho-cp-real" !added:giovanni:debug
      nkbx = nkb
      !
      nx0 = SIZE( x0, 1 )
      !
      nlax    = descla( nlax_ , 1)
      la_proc = ( descla( lambda_node_ , 1) > 0 )
      nlam    = 1
      if ( la_proc ) nlam = nlax

      IF( nx0 /= nlam ) &
         CALL errore( " ortho_cp_real ", " inconsistent dimensions for x0 ", nx0 )

      !
      !
      !     calculation of becp and bephi
      !
      CALL start_clock( 'ortho' )

      CALL nlsm1( nbsp, 1, nvb, eigr,  cp,  becp)

      CALL nlsm1_dist ( nbsp, 1, nvb, eigr, phi, bephi, nlax, nspin, descla )
      !
      !     calculation of qbephi and qbecp
      !
      ALLOCATE( qbephi( nkbx, nx0, nspin ) )
      !
      IF( nvb > 0 ) THEN
         ALLOCATE( bephi_c ( nkbx, nlax*nspin ) )
         CALL redist_row2col_real(nupdwn(1), bephi, bephi_c, nkbx, nlax, descla(1,1) ) !modified:giovanni
         IF( nspin == 2 ) THEN
            CALL redist_row2col_cmplx( nupdwn(2), bephi(1,nlax+1), bephi_c(1,nlax+1), nkbx, nlax, descla(1,2) ) !modified:giovanni
         END IF
      END IF
      !
      qbephi = 0.d0
      !
      DO is=1,nvb
         DO iv=1,nh(is)
            inl = ish(is)+(iv-1)*na(is)
            DO jv=1,nh(is)
               jnl = ish(is)+(jv-1)*na(is)
               qqf = qq(iv,jv,is)
               IF( ABS( qqf ) > 1.D-5 ) THEN
                  DO iss = 1, nspin
                     istart = iupdwn(iss)
                     nc     = descla( nlac_ , iss )
                     ic     = descla( ilac_ , iss ) + istart - 1
                     IF( la_proc ) THEN
                        DO i = 1, nc
                           icc=i+ic-1
                           CALL daxpy( na(is), qqf, bephi_c(jnl+1,i+(iss-1)*nlax),1,qbephi(inl+1,i,iss), 1 ) 
                        END DO
                     END IF
                  END DO
               ENDIF
            END DO
         END DO
      END DO

      IF( nvb > 0 ) DEALLOCATE( bephi_c )
      !
      ALLOCATE( qbecp ( nkbx, nx0, nspin ) )

      qbecp  = 0.d0

      DO is=1,nvb
         DO iv=1,nh(is)
            inl = ish(is)+(iv-1)*na(is)
            DO jv=1,nh(is)
               jnl = ish(is)+(jv-1)*na(is)
               qqf = qq(iv,jv,is)
               IF( ABS( qqf ) > 1.D-5 ) THEN
                  DO iss = 1, nspin
                     istart = iupdwn(iss)
                     nc     = descla( nlac_ , iss )
                     ic     = descla( ilac_ , iss ) + istart - 1
                     IF( la_proc ) THEN
                        DO i = 1, nc
                           icc=i+ic-1
                           CALL daxpy( na(is), qqf, becp (jnl+1,icc),1, qbecp(inl+1,i,iss), 1 )
                        END DO
                     END IF
                  END DO
               ENDIF
            END DO
         END DO
      END DO
      !
      ALLOCATE( xloc( nx0, nx0 ) )
      !
      nspin_sub = nspin 
      if( force_pairing ) nspin_sub = 1
      !
      DO iss = 1, nspin_sub

         nss    = nupdwn(iss)
         istart = iupdwn(iss)

         IF( la_proc ) xloc = x0(:,:,iss) * ccc

         CALL ortho_gamma( 0, cp, ngwx, phi, becp, qbecp(:,:,iss), nkbx, bephi(:,((iss-1)*nlax+1):iss*nlax), &
                           qbephi(:,:,iss), xloc, nx0, descla(:,iss), diff, iter, nbsp, nss, istart )

         IF( iter > ortho_max ) THEN
            WRITE( stdout, 100 ) diff, iter
            CALL errore('ortho','max number of iterations exceeded',iter)
         END IF

         IF( iprsta > 2 ) THEN
            WRITE( stdout, 100 ) diff, iter
         ENDIF
         !     
         IF( la_proc ) x0( :, :, iss ) = xloc / ccc
         !
      END DO
      !
      IF( force_pairing ) cp(:, iupdwn(2):iupdwn(2)+nupdwn(2)-1 ) = cp(:,1:nupdwn(2))
      !
      DEALLOCATE( xloc )
      DEALLOCATE( qbecp )
      DEALLOCATE( qbephi )
      !
      CALL stop_clock( 'ortho' )
      !
      RETURN
      !
100   FORMAT(3X,'diff = ',D18.10,' iter = ', I5 )
      !
   END SUBROUTINE ortho_cp_real

!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_cp_twin( eigr, cp, phi, ngwx, x0, descla, diff, iter, ccc, &
                        bephi, becp, nbsp, nspin, nupdwn, iupdwn )
!=----------------------------------------------------------------------------=!
      !
      !     input = cp (non-orthonormal), beta
      !     input = phi |phi>=s'|c0>
      !     output= cp (orthonormal with s( r(t+dt) ) )
      !     output= bephi, becp
      !     the method used is similar to the version in les houches 1988
      !     'simple molecular systems at..'  p. 462-463  (18-22)
      !      xcx + b x + b^t x^t + a = 1
      !     where c = <s'c0|s|s'c0>   b = <s'c0|s cp>   a = <cp|s|cp>
      !     where s=s(r(t+dt)) and s'=s(r(t))  
      !     for vanderbilt pseudo pot - kl & ap
      !
      USE kinds,          ONLY: DP
      USE ions_base,      ONLY: na, nat
      USE cvan,           ONLY: ish, nvb
      USE uspp,           ONLY: nkb, qq, vkb !added:giovanni:vkb
      USE uspp_param,     ONLY: nh
      USE electrons_base, ONLY: f
      USE gvecw,          ONLY: ngw
      USE control_flags,  ONLY: iprint, iprsta, ortho_max
      USE control_flags,  ONLY: force_pairing, gamma_only, do_wf_cmplx
      USE io_global,      ONLY: stdout, ionode
      USE cp_interfaces,  ONLY: ortho_gamma, nlsm1, nlsm1_dist
      USE descriptors,    ONLY: nlac_ , ilac_ , descla_siz_ , nlar_ , ilar_, lambda_node_, nlax_
      USE mp_global,          ONLY: nproc_image, me_image, intra_image_comm, mpime   ! DEBUG !added:giovanni mpime
      USE twin_types !added:giovanni
      USE parallel_toolkit,  ONLY: redist_row2col
      !
      IMPLICIT NONE
      !
      INTEGER,    INTENT(IN) :: ngwx, nbsp, nspin
      INTEGER,    INTENT(IN) :: nupdwn( nspin ), iupdwn( nspin )
      INTEGER,    INTENT(IN) :: descla(descla_siz_, nspin)
      COMPLEX(DP) :: cp(ngwx,nbsp), phi(ngwx,nbsp), eigr(ngwx,nat)
      type(twin_matrix), dimension(nspin) :: x0
      REAL(DP) :: diff, ccc
      INTEGER :: iter
      type(twin_matrix) :: bephi, becp !modified:giovanni
      !
      REAL(DP), ALLOCATABLE :: xloc(:,:)
      COMPLEX(DP), ALLOCATABLE :: xloc_c(:,:)
      REAL(DP), ALLOCATABLE :: qbephi(:,:,:), qbecp(:,:,:), bephi_c(:,:)
      COMPLEX(DP), ALLOCATABLE :: qbephi_c(:,:,:), qbecp_c(:,:,:), bephi_c_c(:,:)
      !
      INTEGER :: nlam, nlax
      LOGICAL :: la_proc
      !
      INTEGER :: nkbx
      INTEGER :: istart, nss, ifail, i, j, iss, iv, jv, ia, is, inl, jnl
      INTEGER :: nspin_sub, nx0, nc, ic, icc, nr, ir
      REAL(DP) :: qqf
      COMPLEX(DP) :: qqf_c
      LOGICAL :: lgam !added:giovanni
      COMPLEX(DP), PARAMETER :: c_one = CMPLX(1.d0, 0.d0)
      COMPLEX(DP), DIMENSION(nbsp) :: csc
      
      !
!       write(6,*) "inside ortho-cp-twin", bephi%cvec !added:giovanni:debug
      lgam=gamma_only.and..not.do_wf_cmplx
      nkbx = nkb
      !
      IF(.not.x0(1)%iscmplx) THEN
        nx0 = SIZE( x0(1)%rvec, 1 )
      ELSE
        nx0 = SIZE( x0(1)%cvec, 1 )
      ENDIF
      !
      nlax    = descla( nlax_ , 1)
      la_proc = ( descla( lambda_node_ , 1) > 0 )
      nlam    = 1
      if ( la_proc ) nlam = nlax
      write(6,*) "nx0", nx0, nlam, x0(1)%iscmplx
      IF( nx0 /= nlam ) THEN
         CALL errore( " ortho_cp_twin ", " inconsistent dimensions for x0 ", nx0 )
      ENDIF
      !
      !     calculation of becp and bephi
      !
      CALL start_clock( 'ortho' )
!       write(6,*) "starting nlsm1", cp !added:giovanni:debug
      CALL nlsm1( nbsp, 1, nvb, eigr,  cp,  becp, 1, lgam) !modified:giovanni
!       write(6,*) "finished becp", phi !added:giovanni:debug
      CALL nlsm1_dist( nbsp, 1, nvb, eigr, phi, bephi, nlax, nspin, descla, lgam )
!   write(6,*) "finished bephi", bephi%cvec !modified:giovanni
      !
      !     calculation of qbephi and qbecp
      !
      IF(.not.becp%iscmplx) THEN
         ALLOCATE( qbephi( nkbx, nx0, nspin ) )
      ELSE
         ALLOCATE( qbephi_c( nkbx, nx0, nspin ) )
      ENDIF
      !
      IF( nvb > 0 ) THEN
      !begin_modified:giovanni
         ! 
         IF(.not.becp%iscmplx) THEN
            ALLOCATE( bephi_c ( nkbx, nlax*nspin ) )
            CALL redist_row2col_real( nupdwn(1), bephi%rvec(1:nkbx,1:nlax), bephi_c, nkbx, nlax, descla(1,1) )
         ELSE 
            ALLOCATE( bephi_c_c ( nkbx, nlax*nspin ) )
            CALL redist_row2col_cmplx( nupdwn(1), bephi%cvec(1:nkbx,1:nlax), bephi_c_c, nkbx, nlax, descla(1,1) )
         ENDIF
         !
         IF( nspin == 2 ) THEN
	    IF(.not.becp%iscmplx) THEN
		CALL redist_row2col_real( nupdwn(2), bephi%rvec(1:nkbx, &!warning:giovanni ... need to conform to interface
                 nlax+1:nlax+nlax), bephi_c(1:nkbx,  &
                 nlax+1:nlax+nlax), nkbx, nlax, descla(1,2) )
! 		CALL redist_row2col( nupdwn(2), bephi%rvec(1:ubound(bephi%rvec,1), &
!                  nlax+1:ubound(bephi%rvec,2)), bephi_c(1:ubound(bephi_c,1),  &
!                  nlax+1:ubound(bephi_c,2)), nkbx, nlax, descla(1,2) )
	    ELSE
                CALL redist_row2col_cmplx(nupdwn(2), bephi%cvec(1, & !warning:giovanni ... need to conform to interface
                 nlax+1), bephi_c_c(1,  &
                 nlax+1), nkbx, nlax, descla(1,2))
!                 CALL redist_row2col( nupdwn(2), bephi%cvec(1:ubound(bephi%cvec,1), &
!                  nlax+1:ubound(bephi%cvec,2)), bephi_c_c(1:ubound(bephi_c_c,1),  &
!                  nlax+1:ubound(bephi_c_c,2)), nkbx, nlax, descla(1,2) )
	    ENDIF
         END IF
         !
      END IF
!end_modified:giovanni
     ! write(6,*) "bephi_cc", bephi_c_c !added:giovanni:debug
      !
      IF(.not.becp%iscmplx) THEN
	  qbephi = 0.d0
	  !
	  DO is=1,nvb
	    DO iv=1,nh(is)
		inl = ish(is)+(iv-1)*na(is)
		DO jv=1,nh(is)
		  jnl = ish(is)+(jv-1)*na(is)
		  qqf = qq(iv,jv,is)
		  IF( ABS( qqf ) > 1.D-5 ) THEN
		      DO iss = 1, nspin
			istart = iupdwn(iss)
			nc     = descla( nlac_ , iss )
			ic     = descla( ilac_ , iss ) + istart - 1
			IF( la_proc ) THEN
			    DO i = 1, nc
			      icc=i+ic-1
			      CALL daxpy( na(is), qqf, bephi_c(jnl+1,i+(iss-1)*nlax),1,qbephi(inl+1,i,iss), 1 ) 
			    END DO
			END IF
		      END DO
		  ENDIF
		END DO
	    END DO
	  END DO

	  IF( nvb > 0 ) DEALLOCATE( bephi_c )
      ELSE
	  qbephi_c = CMPLX(0.d0, 0.d0)
	  !
	  DO is=1,nvb
	    DO iv=1,nh(is)
		inl = ish(is)+(iv-1)*na(is)
		DO jv=1,nh(is)
		  jnl = ish(is)+(jv-1)*na(is)
		  qqf_c = CMPLX(qq(iv,jv,is),0.d0)
		  IF( ABS( qqf_c ) > 1.D-5 ) THEN
		      DO iss = 1, nspin
			istart = iupdwn(iss)
			nc     = descla( nlac_ , iss )
			ic     = descla( ilac_ , iss ) + istart - 1
			IF( la_proc ) THEN
			    DO i = 1, nc
			      icc=i+ic-1
			      CALL ZAXPY ( na(is), qqf_c, bephi_c_c(jnl+1,i+(iss-1)*nlax), 1, qbephi_c(inl+1,i,iss), 1 ) 
			    END DO
			END IF
		      END DO
		  ENDIF
		END DO
	    END DO
	  END DO

	  IF( nvb > 0 ) DEALLOCATE( bephi_c_c )
      ENDIF
      !
      IF(.not.becp%iscmplx) THEN
	ALLOCATE( qbecp ( nkbx, nx0, nspin ) )

	qbecp  = 0.d0

	DO is=1,nvb
	  DO iv=1,nh(is)
	      inl = ish(is)+(iv-1)*na(is)
	      DO jv=1,nh(is)
		jnl = ish(is)+(jv-1)*na(is)
		qqf = qq(iv,jv,is)
		IF( ABS( qqf ) > 1.D-5 ) THEN
		    DO iss = 1, nspin
		      istart = iupdwn(iss)
		      nc     = descla( nlac_ , iss )
		      ic     = descla( ilac_ , iss ) + istart - 1
		      IF( la_proc ) THEN
			  DO i = 1, nc
			    icc=i+ic-1
			    CALL daxpy( na(is), qqf, becp%rvec (jnl+1,icc),1, qbecp(inl+1,i,iss), 1 )
			  END DO
		      END IF
		    END DO
		ENDIF
	      END DO
	  END DO
	END DO
      ELSE

	ALLOCATE( qbecp_c ( nkbx, nx0, nspin ) )

	qbecp_c(:,:,:)  = CMPLX(0.d0,0.d0)

	DO is=1,nvb
	  DO iv=1,nh(is)
	      inl = ish(is)+(iv-1)*na(is)
	      DO jv=1,nh(is)
		jnl = ish(is)+(jv-1)*na(is)
		qqf_c = CMPLX(qq(iv,jv,is), 0.d0)
		IF( ABS( qqf_c ) > 1.D-5 ) THEN
		    DO iss = 1, nspin
		      istart = iupdwn(iss)
		      nc     = descla( nlac_ , iss )
		      ic     = descla( ilac_ , iss ) + istart - 1
		      IF( la_proc ) THEN
			  DO i = 1, nc
			    icc=i+ic-1
			    CALL ZAXPY(na(is), qqf_c, becp%cvec(jnl+1,icc),1, qbecp_c(inl+1,i,iss), 1)
			  END DO
		      END IF
		    END DO
		ENDIF
	      END DO
	  END DO
	END DO
      ENDIF
      !
      IF(.not.becp%iscmplx) THEN
         ALLOCATE( xloc( nx0, nx0 ) )
      ELSE
         ALLOCATE( xloc_c( nx0, nx0 ) )
      ENDIF
      !
      nspin_sub = nspin 
      if( force_pairing ) nspin_sub = 1
      !
!         csc=CMPLX(0.d0,0.d0) !added:giovanni:debug
!         call scalar_us(becp, nkbx, vkb, cp, ngwx, 15, csc, nbsp, lgam)
!         write(6,*) nbsp, "csc-ortho-before", csc
      DO iss = 1, nspin_sub
!       write(6,*) "calling ortho_gamma in ortho_cp_twin", phi, "qbecp_c", qbecp_c(:,:,iss), "bephi", bephi%cvec!, "qbephi", qbephi_c !added:giovanni:debug
         nss    = nupdwn(iss)
         istart = iupdwn(iss)
         IF(.not.x0(iss)%iscmplx) THEN
	      IF( la_proc ) xloc = x0(iss)%rvec(:,:) * ccc
              CALL ortho_gamma( 0, cp, ngwx, phi, becp%rvec, qbecp(:,:,iss), nkbx, bephi%rvec(:,((iss-1)*nlax+1):iss*nlax), &
	      qbephi(:,:,iss), xloc, nx0, descla(:,iss), diff, iter, nbsp, nss, istart )
         ELSE
              IF( la_proc ) xloc_c = x0(iss)%cvec(:,:) * ccc
              CALL ortho_gamma( 0, cp, ngwx, phi, becp%cvec, (qbecp_c(:,:,iss)), nkbx, bephi%cvec(:,((iss-1)*nlax+1):iss*nlax), &
			  (qbephi_c(:,:,iss)), xloc_c, nx0, descla(:,iss), diff, iter, nbsp, nss, istart )
         ENDIF
        !
!         csc=CMPLX(0.d0,0.d0) !added:giovanni:debug
!         call scalar_us(becp, nkbx, vkb, cp, ngwx, 15, csc, nbsp, lgam)
!         write(6,*) nbsp, "csc-ortho-after", csc
        ! 
	IF( iter > ortho_max ) THEN
	    WRITE( stdout, 100 ) diff, iter
	    CALL errore('ortho','max number of iterations exceeded',iter)
	END IF

	IF( iprsta > 2 ) THEN
	    WRITE( stdout, 100 ) diff, iter
	ENDIF
	!
        IF(.not.x0(iss)%iscmplx) THEN
	      IF( la_proc ) x0(iss)%rvec(:,:) = xloc/ccc
         ELSE
              IF( la_proc ) x0(iss)%cvec(:,:) = xloc_c/ccc
         ENDIF
        !
        
      END DO
      !
      IF( force_pairing ) cp(:, iupdwn(2):iupdwn(2)+nupdwn(2)-1 ) = cp(:,1:nupdwn(2))
      !
      IF(.not.becp%iscmplx) THEN
	DEALLOCATE( xloc )
	DEALLOCATE( qbecp )
	DEALLOCATE( qbephi )
      ELSE
	DEALLOCATE( xloc_c )
	DEALLOCATE( qbecp_c )
	DEALLOCATE( qbephi_c )
      ENDIF

      !
      CALL stop_clock( 'ortho' )
      !
      RETURN
      !
100   FORMAT(3X,'diff = ',D18.10,' iter = ', I5 )
      !
   END SUBROUTINE ortho_cp_twin
