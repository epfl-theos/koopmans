!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"


MODULE kohn_sham_states


   USE io_files, ONLY: ksunit, ks_file, ks_emp_file

   IMPLICIT NONE
   SAVE

   PRIVATE

   ! ...   print KS states to file KS.indx_ksout if ksout true
   LOGICAL :: tksout                        
   LOGICAL :: tksout_emp 

   INTEGER, ALLOCATABLE :: indx_ksout(:,:)  ! (state inds, spin indxs)
   INTEGER, ALLOCATABLE :: n_ksout(:)       ! (spin indxs)
   INTEGER, ALLOCATABLE :: indx_ksout_emp(:,:)  ! (state inds, spin indxs)
   INTEGER, ALLOCATABLE :: n_ksout_emp(:)       ! (spin indxs)

   PUBLIC :: ks_states_init, ks_states_closeup
   PUBLIC :: n_ksout, indx_ksout, tksout, print_all_states

!  ----------------------------------------------
CONTAINS
!  ----------------------------------------------


   SUBROUTINE ks_states_init( nspin, nprnks, iprnks, nprnks_emp, iprnks_emp )

      INTEGER, INTENT(IN) :: nspin, nprnks(:), nprnks_emp(:)
      INTEGER, INTENT(IN) :: iprnks(:,:)
      INTEGER, INTENT(IN) :: iprnks_emp(:,:)

      INTEGER :: i, ip, k, nstates

      ! ...   Tell the code which Kohn-Sham state should be printed to file
      !
      IF( ALLOCATED( n_ksout    ) ) DEALLOCATE( n_ksout )
      IF( ALLOCATED( indx_ksout ) ) DEALLOCATE( indx_ksout )
      !
      tksout = ANY( nprnks > 0 )
      !
      IF( tksout ) THEN
         nstates = MAXVAL( nprnks )
         ALLOCATE( n_ksout( nspin ) )
         ALLOCATE( indx_ksout( nstates, nspin) )
         n_ksout( 1:nspin ) = nprnks( 1:nspin )
         DO i = 1, nspin
           DO k = 1, nprnks( i )
              indx_ksout( k, i ) = iprnks( k, i )
           END DO
         END DO
      END IF

      IF( ALLOCATED( n_ksout_emp    ) ) DEALLOCATE( n_ksout_emp )
      IF( ALLOCATED( indx_ksout_emp ) ) DEALLOCATE( indx_ksout_emp )
      !
      tksout_emp = ANY( nprnks_emp > 0 )
      !
      IF( tksout_emp ) THEN
         nstates = MAXVAL( nprnks_emp )
         ALLOCATE( n_ksout_emp( nspin ) )
         ALLOCATE( indx_ksout_emp( nstates, nspin ) )
         n_ksout_emp( 1:nspin ) = nprnks_emp( 1:nspin )
         DO i = 1, nspin
            DO k = 1, n_ksout_emp( i )
               indx_ksout_emp( k, i ) = iprnks_emp( k, i )
            END DO
         END DO
      END IF

      RETURN
   END SUBROUTINE ks_states_init

!  ----------------------------------------------

   SUBROUTINE ks_states_closeup()
      IF( ALLOCATED( indx_ksout ) ) DEALLOCATE( indx_ksout )
      IF( ALLOCATED( n_ksout ) ) DEALLOCATE( n_ksout )
      tksout = .FALSE.
      IF( ALLOCATED( indx_ksout_emp ) ) DEALLOCATE( indx_ksout_emp )
      IF( ALLOCATED( n_ksout_emp ) ) DEALLOCATE( n_ksout_emp )
      tksout_emp = .FALSE.
      RETURN
   END SUBROUTINE ks_states_closeup

!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE print_all_states( ctot, iupdwn_tot, nupdwn_tot )

        USE kinds,            ONLY : DP
        USE mp_global,        ONLY : intra_image_comm
        USE io_global,        ONLY : ionode
        USE io_global,        ONLY : stdout
        USE electrons_module, ONLY : iupdwn_emp, nupdwn_emp
        USE electrons_base,   ONLY : nupdwn, iupdwn, nspin

        IMPLICIT NONE

        ! ...   declare subroutine arguments
        COMPLEX(DP), INTENT(IN) :: ctot(:,:)
        INTEGER,     INTENT(IN) :: iupdwn_tot(2)
        INTEGER,     INTENT(IN) :: nupdwn_tot(2)

        ! ...   declare other variables
        INTEGER ::  i, iss, iks, itot

        CHARACTER(LEN=256) :: file_name
        CHARACTER(LEN=10), DIMENSION(2) :: spin_name
        CHARACTER (LEN=6), EXTERNAL :: int_to_char

        IF( tksout .OR. tksout_emp ) THEN

          IF (ionode) THEN
            WRITE( stdout,*) 
            WRITE( stdout,'( "   Khon Sham state")') 
            WRITE( stdout,'( "   ---------------")') 
          END IF

          IF( nspin == 2 ) THEN
            spin_name(1) = '_UP_'
            spin_name(2) = '_DW_'
          ELSE
            spin_name(1) = '_'
            spin_name(2) = '_'
          END IF

          DO iss = 1, nspin
            IF( tksout ) THEN
              DO i = 1, n_ksout(iss)
                iks = indx_ksout(i, iss)
                IF( ( iks > 0 ) .AND. ( iks <= nupdwn( iss ) ) ) THEN
                  itot = iks + iupdwn_tot(iss) - 1 
                  file_name = TRIM( ks_file ) // &
                            & trim(spin_name(iss)) // trim( int_to_char( iks ) )
                  CALL print_ks_states( ctot( :, itot ), file_name )
                END IF
              END DO
            END IF
            IF( tksout_emp ) THEN
              DO i = 1, n_ksout_emp(iss)
                iks = indx_ksout_emp(i, iss)
                IF( ( iks > 0 ) .AND. ( iks <= nupdwn_emp( iss ) ) ) THEN
                  itot = iks + iupdwn_tot(iss) + nupdwn( iss ) - 1 
                  file_name = TRIM( ks_emp_file ) // &
                            & trim(spin_name(iss)) // trim( int_to_char( iks ) )
                  CALL print_ks_states( ctot( :, itot ), file_name )
                END IF
              END DO
            END IF
          END DO

        END IF

        RETURN
        ! ...
      END SUBROUTINE print_all_states


!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE print_ks_states( c, file_name )

        USE kinds
        USE mp, ONLY: mp_sum
        USE io_global, ONLY: ionode, ionode_id
        USE io_global, ONLY: stdout
        USE gvecw, ONLY: ngw
        USE fft_base, ONLY: dfftp, dffts
        USE grid_dimensions, ONLY: nr1, nr2, nr3, nr1x, nr2x, nr3x, nnrx
        USE cp_interfaces, ONLY: invfft
        USE xml_io_base, ONLY: write_rho_xml
        USE mp_global,       ONLY: nproc_image, me_image, intra_image_comm

        IMPLICIT NONE

        COMPLEX(DP),      INTENT(IN) :: c(:)
        CHARACTER(LEN=*), INTENT(IN) :: file_name
        REAL(DP),    ALLOCATABLE :: rpsi2(:)
        COMPLEX(DP), ALLOCATABLE :: psi(:)
        INTEGER   ::  i
        REAL(DP) :: charge

        ALLOCATE( psi( nnrx ) )
        ALLOCATE( rpsi2( nnrx ) )

        CALL c2psi( psi, dffts%nnr, c, c, ngw, 1 )
        CALL invfft( 'Wave', psi, dffts )

        DO i = 1, nnrx
           rpsi2( i ) = DBLE( psi( i ) )**2
        END DO
        charge = SUM( rpsi2 )

        CALL write_rho_xml( file_name, rpsi2, &
                            nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
        
        CALL mp_sum( charge, intra_image_comm )

        IF ( ionode ) THEN
          WRITE( stdout,'(3X,A15," integrated charge : ",F14.5)')  &
     &      TRIM(file_name), charge / DBLE(nr1*nr2*nr3)
        END IF

        DEALLOCATE( rpsi2, psi )
        ! ...
        RETURN
        ! ...
      END SUBROUTINE print_ks_states

!  ----------------------------------------------
!
END MODULE kohn_sham_states
