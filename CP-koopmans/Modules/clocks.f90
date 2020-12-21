!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mytime
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER,  PARAMETER :: maxclock = 100
  REAL(DP), PARAMETER :: notrunning = - 1.0_DP
  ! 
  REAL(DP)          :: cputime(maxclock), t0cpu(maxclock)
  REAL(DP)          :: walltime, t0wall
  CHARACTER(LEN=12) :: clock_label(maxclock)
  INTEGER           :: called(maxclock)
  !
  INTEGER :: nclock = 0
  LOGICAL :: no
  INTEGER :: trace_depth = 0
  !
END MODULE mytime
!
!#define __TRACE
!----------------------------------------------------------------------------
SUBROUTINE init_clocks( go )
  !----------------------------------------------------------------------------
  !
  ! ... flag = .TRUE.  : clocks will run
  ! ... flag = .FALSE. : only clock #1 will run
  !
  USE kinds,  ONLY : DP
  USE mytime, ONLY : called, t0cpu, cputime, no, notrunning, maxclock, &
       clock_label, walltime, t0wall
  !
  IMPLICIT NONE
  !
  LOGICAL :: go
  INTEGER :: n
  !
  !
  no = .NOT. go
  !
  DO n = 1, maxclock
     !
     called(n)      = 0
     cputime(n)     = 0.0_DP
     t0cpu(n)       = notrunning
     clock_label(n) = ' '
     !
  END DO
  !
  t0wall = 0.0_DP
  walltime = 0.0_DP
  !
  RETURN
  !
END SUBROUTINE init_clocks
!
!----------------------------------------------------------------------------
SUBROUTINE start_clock( label )
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE mp_global, ONLY : mpime
  USE mytime,    ONLY : nclock, clock_label, notrunning, no, maxclock, &
                        t0cpu, t0wall, trace_depth
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*) :: label
  !
  CHARACTER(LEN=12) :: label_
  INTEGER          :: n
  REAL(DP), EXTERNAL :: scnds, cclock
  !
#if defined (__TRACE)
  WRITE( stdout, '("mpime = ",I2,", TRACE (depth=",I2,") Start: ",A12)') mpime, trace_depth, label
  trace_depth = trace_depth + 1
#endif
  !
  IF ( no .AND. ( nclock == 1 ) ) RETURN
  !
  ! ... prevent trouble if label is longer than 12 characters 
  !
  label_ = TRIM ( label ) 
  !
  DO n = 1, nclock
     !
     IF ( clock_label(n) == label_ ) THEN
        !
        ! ... found previously defined clock: check if not already started,
        ! ... store in t0cpu the starting time
        !
        IF ( t0cpu(n) /= notrunning ) THEN
!            WRITE( stdout, '("start_clock: clock # ",I2," for ",A12, &
!                           & " already started")' ) n, label_
        ELSE
#ifdef __WALLTIME
           t0cpu(n) = cclock()
#else
           t0cpu(n) = scnds()
#endif
           IF ( n == 1 ) t0wall = cclock()
        END IF
        !
        RETURN
        !
     END IF
     !
  END DO
  !
  ! ... clock not found : add new clock for given label
  !
  IF ( nclock == maxclock ) THEN
     !
     WRITE( stdout, '("start_clock: Too many clocks! call ignored")' )
     !
  ELSE
     !
     nclock              = nclock + 1
     clock_label(nclock) = label_
#ifdef __WALLTIME
     t0cpu(nclock)          = cclock()
#else
     t0cpu(nclock)          = scnds()
#endif
     IF ( nclock == 1 ) t0wall = cclock()
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE start_clock
!
!----------------------------------------------------------------------------
SUBROUTINE stop_clock( label )
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE mp_global, ONLY : mpime
  USE mytime,    ONLY : no, nclock, clock_label, cputime, walltime, &
                        notrunning, called, t0cpu, t0wall, trace_depth
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*) :: label
  !
  CHARACTER(LEN=12) :: label_
  INTEGER          :: n
  REAL(DP), EXTERNAL :: scnds, cclock
  !
#if defined (__TRACE)
  trace_depth = trace_depth - 1
  WRITE( *, '("mpime = ",I2,", TRACE (depth=",I2,") End: ",A12)') mpime, trace_depth, label
#endif
  !
  IF ( no ) RETURN
  !
  ! ... prevent trouble if label is longer than 12 characters 
  !
  label_ = TRIM ( label ) 
  ! 
  DO n = 1, nclock
     !
     IF ( clock_label(n) == label_ ) THEN
        !
        ! ... found previously defined clock : check if properly initialised,
        ! ... add elapsed time, increase the counter of calls
        !
        IF ( t0cpu(n) == notrunning ) THEN
           !
!            WRITE( stdout, '("stop_clock: clock # ",I2," for ",A12, &
!                           & " not running")' ) n, label
           !
        ELSE
           !
#ifdef __WALLTIME
           cputime(n) = cputime(n) + cclock() - t0cpu(n)
#else
           cputime(n) = cputime(n) + scnds() - t0cpu(n)
#endif
           IF ( n == 1 ) walltime = walltime + cclock() - t0wall
           t0cpu(n)      = notrunning
           called(n)  = called(n) + 1
           !
        END IF
        !
        RETURN
        !
     END IF
     !
  END DO
  !
  ! ... clock not found
  !
  WRITE( stdout, '("stop_clock: no clock for ",A12," found !")' ) label
  !
  RETURN
  !
END SUBROUTINE stop_clock
!
!----------------------------------------------------------------------------
SUBROUTINE print_clock( label )
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE mytime,    ONLY : nclock, clock_label
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*) :: label
  !
  CHARACTER(LEN=12) :: label_
  INTEGER          :: n
  !
  IF ( label == ' ' ) THEN
     !
     WRITE( stdout, * )
     !
     DO n = 1, nclock
        !
        CALL print_this_clock( n )
        !
     END DO
     !
  ELSE
     !
     ! ... prevent trouble if label is longer than 12 characters 
     !
     label_ = TRIM ( label ) 
     !
     DO n = 1, nclock
        !
        IF ( clock_label(n) == label_ ) THEN
           !
           CALL print_this_clock( n )
           !
           EXIT
           !
        END IF
        !
     END DO
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE print_clock
!
! !----------------------------------------------------------------------------
! SUBROUTINE print_clock_to_myfile( label )
!   !----------------------------------------------------------------------------
!   !
!   USE kinds,     ONLY : DP
!   USE io_global, ONLY : stdout
!   USE mytime,    ONLY : nclock, clock_label
!   !
!   IMPLICIT NONE
!   !
!   CHARACTER(LEN=*) :: label
!   !
!   CHARACTER(LEN=12) :: label_
!   INTEGER          :: n
!   !
!   open(unit=40,file='aida.times',access = 'append')
!   IF ( label == ' ' ) THEN
!      !
!      WRITE( stdout, * )
!      !
!      DO n = 1, nclock
!         !
!         CALL print_this_clock( n )
!         !
!      END DO
!      !
!   ELSE
!      !
!      ! ... prevent trouble if label is longer than 12 characters 
!      !
!      label_ = TRIM ( label ) 
!      !
!      DO n = 1, nclock
!         !
!         IF ( clock_label(n) == label_ ) THEN
!            !
!            CALL print_this_clock( n )
!            !
!            EXIT
!            !
!         END IF
!         !
!      END DO
!      !
!   END IF
!   !
!   RETURN
!   !
! END SUBROUTINE print_clock_to_myfile
!----------------------------------------------------------------------------
SUBROUTINE print_this_clock( n )
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE mytime,    ONLY : no, nclock, clock_label, cputime, walltime, &
                        notrunning, called, t0cpu, t0wall
  USE mp,        ONLY : mp_max
  USE mp_global, ONLY : intra_image_comm, my_image_id
  !
  IMPLICIT NONE
  !
  INTEGER  :: n
  REAL(DP) :: elapsed_cpu_time, elapsed_wall_time, nsec, msec
  INTEGER  :: nday, nhour, nmin, nmax, mday, mhour, mmin
  !
  REAL(DP), EXTERNAL :: scnds, cclock
  !
  !
  IF ( t0cpu(n) == notrunning ) THEN
     !
     ! ... clock stopped, print the stored value for the cpu time
     !
     elapsed_cpu_time = cputime(n)
     elapsed_wall_time= walltime
     !
  ELSE
     !
     ! ... clock not stopped, print the current value of the cpu time
     !
#ifdef __WALLTIME
     elapsed_cpu_time = cputime(n) + cclock() - t0cpu(n)
#else
     elapsed_cpu_time = cputime(n) + scnds() - t0cpu(n)
#endif
     elapsed_wall_time = walltime + cclock() - t0wall
     !
  END If
  !
  nmax = called(n)
  !
  ! ... In the parallel case there are several possible approaches
  ! ... The safest one is to leave each clock independent from the others
  ! ... Another possibility is to print the maximum across all processors
  ! ... This is done by uncommenting the following lines
  !
  ! CALL mp_max( elapsed_cpu_time, intra_image_comm )
  ! CALL mp_max( elapsed_wall_time, intra_image_comm )
  ! CALL mp_max( nmax, intra_image_comm )
  !
  ! ... In the last line we assume that the maximum cpu time
  ! ... is associated to the maximum number of calls
  ! ... NOTA BENE: by uncommenting the above lines you may run into
  ! ... serious trouble if clocks are not started on all nodes
  !
  IF ( n == 1 ) THEN
     !
     ! ... The first clock is written as days/hour/min/sec
     !
     nday  = elapsed_cpu_time / 86400
     nsec  = elapsed_cpu_time - 86400 * nday
     nhour = nsec / 3600 
     nsec  = nsec - 3600 * nhour
     nmin  = nsec / 60
     nsec  = nsec - 60 * nmin
     !
     ! ... The first clock writes elapsed (wall) time as well
     !
     mday  = elapsed_wall_time / 86400
     msec  = elapsed_wall_time - 86400 * mday
     mhour = msec / 3600 
     msec  = msec - 3600 * mhour
     mmin  = msec / 60
     msec  = msec - 60 * mmin
     !
     IF ( nday > 0 .OR. mday > 0 ) THEN
        !    
#ifdef __WALLTIME
        WRITE( stdout, &
               '(5X,A12," : ",3X,I2,"d",3X,I2,"h",I2, "m wall time"/)' ) &
             clock_label(n), mday, mhour, mmin
#else
        WRITE( stdout, &
               '(5X,A12," : ",3X,I2,"d",3X,I2,"h",I2, "m CPU time, ", &
           &            "   ",3X,I2,"d",3X,I2,"h",I2, "m wall time"/)' ) &
             clock_label(n), nday, nhour, nmin, mday, mhour, mmin
#endif
        !
     ELSE IF ( nhour > 0 .OR. mhour > 0 ) THEN
        !
#ifdef __WALLTIME
        WRITE( stdout, &
               '(5X,A12," : ",3X,I2,"h",I2,"m wall time"/)' ) &
             clock_label(n), mhour, mmin
#else
        WRITE( stdout, &
               '(5X,A12," : ",3X,I2,"h",I2,"m CPU time, ", &
           &            "   ",3X,I2,"h",I2,"m wall time"/)' ) &
             clock_label(n), nhour, nmin, mhour, mmin
#endif
        !
     ELSE IF ( nmin > 0 .OR. mmin > 0 ) THEN
        !
#ifdef __WALLTIME
        WRITE( stdout, &
               '(5X,A12," : ",I2,"m",F5.2,"s wall time"/)' ) &
             clock_label(n), mmin, msec
#else
        WRITE( stdout, &
               '(5X,A12," : ",I2,"m",F5.2,"s CPU time, ", &
               &        "   ",I2,"m",F5.2,"s wall time"/)' ) &
             clock_label(n), nmin, nsec, mmin, msec
#endif
        !
     ELSE
        !
#ifdef __WALLTIME
        WRITE( stdout, &
               '(5X,A12," : ",3X,F5.2,"s wall time"/)' )&
             clock_label(n), msec
#else
        WRITE( stdout, &
               '(5X,A12," : ",3X,F5.2,"s CPU time,",3X,F5.2,"s wall time"/)' )&
             clock_label(n), nsec, msec
#endif
        !
     END IF
     !
  ELSE IF ( nmax == 1 .OR. t0cpu(n) /= notrunning ) THEN
     !
     ! ... for clocks that have been called only once
     !
     WRITE( stdout, &
            '(5X,A12," :",F9.2,"s CPU")') clock_label(n), elapsed_cpu_time
     !
  ELSE IF ( nmax == 0 ) THEN
     !
     ! ... for clocks that have never been called
     !
     WRITE( stdout, &
            '("print_this: clock # ",I2," for ",A12," never called !")' ) &
         n, clock_label(n)
     !
  ELSE
     !
     ! ... for all other clocks
     !
     WRITE( stdout, &
            '(5X,A12," :",F9.2,"s CPU (",I8," calls,",F8.3," s avg)")' ) &
         clock_label(n), elapsed_cpu_time, &
         nmax, ( elapsed_cpu_time / nmax )
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE print_this_clock
!
!----------------------------------------------------------------------------
FUNCTION get_clock( label )
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE mytime,    ONLY : no, nclock, clock_label, cputime, &
                        notrunning, called, t0cpu
  USE mp,        ONLY : mp_max
  USE mp_global, ONLY : intra_image_comm 
  !
  IMPLICIT NONE
  !
  REAL(DP)         :: get_clock
  CHARACTER(LEN=*) :: label
  INTEGER          :: n
  !
  REAL(DP), EXTERNAL :: scnds, cclock
  !
  !
  IF ( no ) THEN
     !
     IF ( label == clock_label(1) ) THEN
        !
#ifdef __WALLTIME
        get_clock = cclock()
#else
        get_clock = scnds()
#endif
        !
     ELSE
        !
        get_clock = notrunning
        !
     END IF
     !
     RETURN
     !
  END IF
  !
  DO n = 1, nclock
     !
     IF ( label == clock_label(n) ) THEN
        !
        IF ( t0cpu(n) == notrunning ) THEN
           !
           get_clock = cputime(n)
           !
        ELSE
           !
#ifdef __WALLTIME
           get_clock = cputime(n) + cclock() - t0cpu(n)
#else
           get_clock = cputime(n) + scnds() - t0cpu(n)
#endif
           !
        END IF
        !
        ! ... See comments in subroutine print_this_clock
        !
        ! CALL mp_max( get_clock, intra_image_comm )
        !
        RETURN
        !
     END IF
     !
  END DO
  !
  ! ... clock not found
  !
  get_clock = notrunning
  !
  WRITE( stdout, '("get_clock: no clock for ",A12," found !")') label
  !
  RETURN
  !
END FUNCTION get_clock
