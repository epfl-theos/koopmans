!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
      SUBROUTINE writetofile( f, nnrx, filename, dfft, which_print )
!--------------------------------------------------------------------
      !
      USE kinds,             ONLY : DP
      USE io_global,         ONLY : meta_ionode, meta_ionode_id
      USE mp_global,         ONLY : nproc_image, me_image, intra_image_comm
      USE mp,                ONLY : mp_barrier, mp_gather
      USE fft_types,         ONLY : fft_dlay_descriptor
      USE control_flags,     ONLY : iprsta
      USE cell_base,         ONLY : a1, a2, a3
      !
      IMPLICIT NONE
      !
      INTEGER,   INTENT(IN) :: nnrx
      REAL(DP),  INTENT(IN) :: f(nnrx)
      TYPE(fft_dlay_descriptor), INTENT(IN)  :: dfft
      !
      CHARACTER( LEN=* ), INTENT(IN) :: filename
      CHARACTER( LEN=* ), INTENT(IN) :: which_print
      !
      INTEGER               :: ir, ir1, ir2, ir3, num
      INTEGER               :: nr1, nr2, nr3, nr1x, nr2x, nr3x
      INTEGER               :: proc, ierr
      REAL(DP)              :: total, delta1m, delta2m, delta3m
      REAL(DP), ALLOCATABLE :: fdist(:)
      INTEGER,  ALLOCATABLE :: displs(:), recvcount(:)
      !
      INTEGER, EXTERNAL     :: compindex
      

      !
      ! write only if verbosity is "high"
      !
      IF ( iprsta < 3 ) RETURN
      !
      ! first collect f across the processors
      !
      nr1=dfft%nr1
      nr2=dfft%nr2
      nr3=dfft%nr3
      !
      nr1x=dfft%nr1x
      nr2x=dfft%nr2x
      nr3x=dfft%nr3x
      !
      delta1m = a1(1)/DBLE( nr1 ) 
      delta2m = a2(2)/DBLE( nr2 ) 
      delta3m = a3(3)/DBLE( nr3 ) 
      !
      ALLOCATE( displs( nproc_image ), recvcount( nproc_image ) )
      !
      IF ( meta_ionode ) THEN
          ALLOCATE( fdist(nr1x*nr2x*nr3x), STAT=ierr )
      ELSE
          ALLOCATE( fdist(1), STAT=ierr )
      ENDIF
      IF ( ierr/=0 ) CALL errore('writetofile','allocating fdist', ABS(ierr))
      !
      IF ( nproc_image > 1 ) THEN
          !
          DO proc=1,nproc_image
              !
              recvcount(proc) =  dfft%nnp * ( dfft%npp(proc) )
              !
              IF (proc == 1) THEN
                 displs(proc)=0
              ELSE
                 displs(proc)=displs(proc-1) + recvcount(proc-1)
              ENDIF
              !
          ENDDO
          !
          ! gather the charge density on the first node
          !
          call mp_barrier()
          call mp_gather( f, fdist, recvcount, displs, meta_ionode_id, intra_image_comm )
          !
      ELSE
          ! 
          ! one processor per image  
          !  
          IF ( nr1 /= nr1x .OR. nr2 /= nr2x .OR. nr3 /= nr3x ) &
             CALL errore('writetofile','dimension mistmatch',10)
          !
          fdist(1:nr1x*nr2x*nr3x) = f(1:nnrx)
          !
      ENDIF
 
      !
      ! write data to file
      !
      IF ( meta_ionode ) THEN
         ! 
         OPEN( 300, file = TRIM( filename ), status = 'unknown' )
         !
         SELECT CASE( TRIM( which_print ) )
           !
         CASE( 'all' )
           !
           WRITE( 300, * ) fdist(:)
           !
         CASE( 'x' )
           !
           DO ir1 = 1, nr1
             total = 0.D0
             num = 0
             DO ir2 = nr2 / 2, nr2 / 2
               DO ir3 = nr3 / 2, nr3 / 2
                 ir = compindex( ir1, ir2, ir3, nr1, nr2, nr3 )
                 total = total + fdist( ir )
                 num = num + 1
               END DO
             END DO
             total = total / DBLE( num )
             WRITE( 300, '(2E30.10)' ) DBLE( ir1 - 1 )   &
                  * delta1m, total
           END DO
           !
         CASE( 'y' )
           !
           DO ir2 = 1, nr2
             total = 0.D0
             num = 0
             DO ir1 = nr1 / 2, nr1 / 2
               DO ir3 = nr3 / 2, nr3 / 2
                 ir = compindex( ir1, ir2, ir3, nr1, nr2, nr3 )
                 total = total + fdist( ir )
                 num = num + 1
               END DO
             END DO
             total = total / DBLE( num )
             WRITE( 300, '(2E30.10)' ) DBLE( ir2 - 1 )   &
                  * delta2m, total
           END DO
           !
         CASE( 'z' )
           !
           DO ir3 = 1, nr3
             total = 0.D0
             num = 0
             DO ir2 = 0, 0 ! nr2 / 2, nr2 / 2
               DO ir1 = 0, 0 !nr1 / 2, nr1 / 2
                 ir = compindex( ir1, ir2, ir3, nr1, nr2, nr3 )
                 total = total + fdist( ir )
                 num = num + 1
               END DO
             END DO
             total = total / DBLE( num )
             WRITE( 300, '(2E30.10)' ) DBLE( ir3 - 1 )   &
                  * delta3m, total
           END DO
           !
         CASE( 'ax' )
           !
           DO ir1 = 1, nr1
             total = 0.D0
             num = 0
             DO ir2 = 1, nr2
               DO ir3 = 1, nr3
                 ir = compindex( ir1, ir2, ir3, nr1, nr2, nr3 )
                 total = total + fdist( ir )
                 num = num + 1
               END DO
             END DO
             total = total / DBLE( num )
             WRITE( 300, '(2E30.10)' ) DBLE( ir1 - 1 )   &
                  * delta1m, total
           END DO
           !
          CASE( 'ay' )
           !
           DO ir2 = 1, nr2
             total = 0.D0
             num = 0
             DO ir1 = 1, nr1
               DO ir3 = 1, nr3
                 ir = compindex( ir1, ir2, ir3, nr1, nr2, nr3 )
                 total = total + fdist( ir )
                 num = num + 1
               END DO
             END DO
             total = total / DBLE( num )
             WRITE( 300, '(2E30.10)' ) DBLE( ir2 - 1 )   &
                  * delta2m, total
           END DO
           !
         CASE( 'az' )
           !
           DO ir3 = 1, nr3
             total = 0.D0
             num = 0
             DO ir2 = 1, nr2
               DO ir1 = 1, nr1
                 ir = compindex( ir1, ir2, ir3, nr1, nr2, nr3 )
                 total = total + fdist( ir )
                 num = num + 1
               END DO
             END DO
             total = total / DBLE( num )
             WRITE( 300, '(2E30.10)' ) DBLE( ir3 - 1 )   &
                  * delta3m, total
           END DO
           !
         END SELECT
         !
         CLOSE( 300 )
         !
      ENDIF
    
      ! 
      ! cleanup local memory   
      !   
      DEALLOCATE( displs, recvcount )
      !
      DEALLOCATE( fdist, STAT=ierr )
      IF ( ierr/=0 ) CALL errore('writetofile','deallocating fdist', ABS(ierr))
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE writetofile
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION compindex( ir1, ir2, ir3, nr1, nr2, nr3 )
!--------------------------------------------------------------------
      ! ... Calculates the composite grid index corresponding
      ! ... to ir1, ir2, ir3
      !
      IMPLICIT NONE
      !
      INTEGER :: compindex
      INTEGER, INTENT(IN) :: ir1, ir2, ir3
      INTEGER, INTENT(IN) :: nr1, nr2, nr3
      INTEGER             :: jr1, jr2, jr3
      !
      jr1 = MODULO( ir1 - 1 , nr1 ) + 1
      jr2 = MODULO( ir2 - 1 , nr2 ) + 1
      jr3 = MODULO( ir3 - 1 , nr3 ) + 1
      !
      compindex = jr1 + ( jr2 -1 ) * nr1 + ( jr3 - 1 ) * nr1 * nr2
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION compindex
!--------------------------------------------------------------------
