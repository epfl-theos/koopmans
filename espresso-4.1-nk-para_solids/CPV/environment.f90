!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!==-----------------------------------------------------------------------==!
   MODULE environment
!==-----------------------------------------------------------------------==!

        USE kinds
        USE io_files, ONLY: crash_file, crashunit, &
                            stop_file, stopunit

        IMPLICIT NONE
        SAVE

        PRIVATE

        REAL(DP)          :: start_cclock_val

        PUBLIC :: environment_start
        PUBLIC :: environment_end
        PUBLIC :: opening_date_and_time, closing_date_and_time
        PUBLIC :: start_cclock_val
      
!==-----------------------------------------------------------------------==!
   CONTAINS
!==-----------------------------------------------------------------------==!
    

        SUBROUTINE environment_start( )

          USE io_global, ONLY: stdout, meta_ionode
          USE mp_global, ONLY: mpime, nproc, me_image, &
               my_image_id, root_image
          USE cp_version

          LOGICAL           :: texst
          INTEGER           :: nchar
          CHARACTER(LEN=80) :: uname
          CHARACTER(LEN=80) :: version_str
          REAL(DP),         EXTERNAL :: cclock
          CHARACTER(LEN=6), EXTERNAL :: int_to_char

#if defined __OPENMP
          INTEGER, EXTERNAL :: omp_get_max_threads
#endif

          CALL init_clocks( .TRUE. )
          CALL start_clock( 'CP' )

          start_cclock_val = cclock( )

          version_str = TRIM (version_number) // " - " // TRIM (version_date)

          ! ...  search for file CRASH and delete it

          IF( meta_ionode ) THEN
            INQUIRE( FILE=TRIM(crash_file), EXIST=texst )
            IF( texst ) THEN
              OPEN(  UNIT=crashunit, FILE=TRIM(crash_file), STATUS='OLD' )
              CLOSE( UNIT=crashunit, STATUS='DELETE' )
            END IF
          END IF

          ! ...       each processor other than me=1 opens its own standard output file,
          ! ...       this is mainly for debugging, usually only ionode writes to stdout

          IF( .NOT. meta_ionode ) THEN

             uname = 'out.' // trim(int_to_char( my_image_id )) // '_' // &
                  trim(int_to_char( me_image))
             nchar = INDEX(uname,' ') - 1

             !
             ! useful for debugging purposes 
             if (me_image == root_image) then
                open( unit = stdout, file = uname(1:nchar),status='unknown')
             else
                open( unit = stdout, file='/dev/null', status='unknown' )
             endif
             !open( unit = stdout, file = uname(1:nchar),status='unknown')
             !
             !open( unit = stdout, file='/dev/null', status='unknown' )

          END IF
!
          CALL opening_date_and_time( version_str )

#if defined __MPI

          WRITE( stdout,100)  nproc, mpime
100       FORMAT(3X,'MPI Parallel Build',/,3X,'Tasks =',I5,'  This task id =',I5)

#else

          WRITE( stdout,100)  
100       FORMAT(3X,'Serial Build')

#endif

#if defined __OPENMP
          WRITE( stdout,110) omp_get_max_threads()
110       FORMAT(3X,'Using OpenMP with',I5,' threads')
#endif

          RETURN
        END SUBROUTINE environment_start

!==-----------------------------------------------------------------------==!

        SUBROUTINE environment_end( )

          USE io_global, ONLY: stdout, meta_ionode

          IF ( meta_ionode ) WRITE( stdout, * )

          CALL stop_clock(  'CP' )
          CALL print_clock( 'CP' )

          CALL closing_date_and_time( )

          IF( meta_ionode ) THEN
            WRITE( stdout,'(A)')      '   JOB DONE.'
            WRITE( stdout,3335)
          END IF
 3335     FORMAT('=',78('-'),'=')

          RETURN
        END SUBROUTINE environment_end

!==-----------------------------------------------------------------------==!

        SUBROUTINE opening_date_and_time( version_str )

          USE io_global, ONLY: stdout, meta_ionode

          CHARACTER(LEN=*), INTENT(IN) :: version_str
          CHARACTER(LEN=9)  :: cdate, ctime
          CHARACTER(LEN=80) :: time_str

          CALL date_and_tim( cdate, ctime )
          time_str = 'This run was started on:  ' // ctime // ' ' // cdate

! ...     write program heading


          IF( meta_ionode ) THEN
            WRITE( stdout,3331) 
            WRITE( stdout,3332) version_str
            WRITE( stdout,3331) 
            WRITE( stdout,3334) time_str
          END IF

 3331     FORMAT('=',78('-'),'=')
 3332     FORMAT( /, 5X,'CP: variable-cell Car-Parrinello molecular dynamics',/&
        & ,5X,'using norm-conserving and ultrasoft Vanderbilt pseudopotentials',//&
        & ,5X,'Version: ',A60,/&
        & ,5X,'Authors: Alfredo Pasquarello, Kari Laasonen, Andrea Trave, Roberto Car,',/&
        & ,5X,'  Paolo Giannozzi, Nicola Marzari, Carlo Cavazzoni, Guido Chiarotti,',/&
        & ,5X,'  Sandro Scandolo, Paolo Focher, Gerardo Ballabio, and others',/)

 3334     FORMAT(/,3X,A60,/)
          RETURN
        END SUBROUTINE opening_date_and_time

!==-----------------------------------------------------------------------==!

        SUBROUTINE closing_date_and_time( )

          USE io_global, ONLY: stdout, meta_ionode

          CHARACTER(LEN=9)  :: cdate, ctime
          CHARACTER(LEN=80) :: time_str

          CALL date_and_tim( cdate, ctime )

          time_str = 'This run was terminated on:  ' // ctime // ' ' // cdate

          IF( meta_ionode ) THEN
            WRITE( stdout,*)
            WRITE( stdout,3334) time_str
            WRITE( stdout,3335)
          END IF

 3334     FORMAT(3X,A60,/)
 3335     FORMAT('=',78('-'),'=')

          RETURN
        END SUBROUTINE closing_date_and_time

!==-----------------------------------------------------------------------==!
   END MODULE environment
!==-----------------------------------------------------------------------==!
