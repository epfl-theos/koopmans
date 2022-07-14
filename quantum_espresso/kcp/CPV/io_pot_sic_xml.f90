!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE io_pot_sic_xml
  !----------------------------------------------------------------------------
  !
  USE kinds,       ONLY : DP
  USE xml_io_base, ONLY : create_directory, write_pot_xml, read_pot_xml, &
                          restart_dir
  !
  PRIVATE
  !
  PUBLIC :: write_pot_sic, read_pot_sic
  !
  ! {read|write}_pot_only:    read or write the real space charge density
  ! {read|write}_pot_general: as above, plus read or write ldaU ns coeffs
  !                           and PAW becsum coeffs.

  INTERFACE write_pot_sic
        MODULE PROCEDURE write_pot_sic_only
  END INTERFACE

  INTERFACE read_pot_sic
        MODULE PROCEDURE read_pot_sic_only
  END INTERFACE

  CONTAINS

    !------------------------------------------------------------------------
    SUBROUTINE write_pot_sic_only( pot, extension, field_specifier)
      !------------------------------------------------------------------------
      !
      ! ... this routine writes the charge-density in xml format into the
      ! ... '.save' directory
      ! ... the '.save' directory is created if not already present
      !
      USE io_files, ONLY : outdir, prefix
      USE fft_base, ONLY : dfftp
      USE io_global, ONLY : ionode
      USE mp_global, ONLY : intra_pool_comm, inter_pool_comm
      USE control_flags, ONLY : ndw
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(IN)           :: pot(dfftp%nnr)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: field_specifier
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension
      !
      CHARACTER(LEN=256)    :: dirname, file_base
      CHARACTER(LEN=256)    :: ext
      REAL(DP), ALLOCATABLE :: potaux(:)
      !
      !
      ext = ' '
      !
      !dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      dirname = restart_dir( outdir, ndw )
      !
      CALL create_directory( dirname )
      !
      IF ( PRESENT( extension ) ) ext = '.' // TRIM( extension )
      !
      IF ( PRESENT( extension ) ) THEN 
        file_base = TRIM( dirname ) // '/' // TRIM(field_specifier) // TRIM( ext )
      ELSE 
        file_base = TRIM( dirname ) // '/sic_potential' // TRIM( ext )
      END IF
      !
      ! 
      CALL write_pot_xml( file_base, pot(:), dfftp%nr1, dfftp%nr2, &
              dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, &
              ionode, intra_pool_comm, inter_pool_comm )
      RETURN
      !
    END SUBROUTINE write_pot_sic_only
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_pot_sic_only( pot, extension )
      !------------------------------------------------------------------------
      !
      ! ... this routine reads the effective potential in xml format from the
      ! ... files saved into the '.save' directory
      !
      USE io_files, ONLY : tmp_dir, prefix
      USE fft_base, ONLY : dfftp
      USE io_global, ONLY : ionode
      USE mp_global, ONLY : intra_pool_comm, inter_pool_comm
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(OUT)          :: pot(dfftp%nnr)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension
      !
      CHARACTER(LEN=256)    :: dirname, file_base
      CHARACTER(LEN=256)    :: ext
      !
      ext = ' '
      !
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      !
      IF ( PRESENT( extension ) ) ext = '.' // TRIM( extension )
      !
      file_base = TRIM( dirname ) // '/sic-potential' // TRIM( ext )
      !
      !
      CALL read_pot_xml( file_base, pot(:), dfftp%nr1, dfftp%nr2, &
               dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, &
               ionode, intra_pool_comm, inter_pool_comm ) 
      !
      RETURN
      !
    END SUBROUTINE read_pot_sic_only
    !
END MODULE io_pot_sic_xml
