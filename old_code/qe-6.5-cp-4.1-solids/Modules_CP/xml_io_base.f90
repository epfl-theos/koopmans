!
! Copyright (C) 2005-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE xml_io_base
  !----------------------------------------------------------------------------
  !
  ! ... this module contains some common subroutines used to read and write
  ! ... in XML format the data produced by Quantum-ESPRESSO package
  !
  ! ... written by Carlo Sbraccia (2005)
  ! ... modified by Andrea Ferretti (2006-08)
  !
  USE iotk_module
  !
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : tmp_dir, prefix, iunpun, xmlpun, &
                        current_fmt_version => qexml_version
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE parser,    ONLY : version_compare
  !
  IMPLICIT NONE
  PRIVATE
  !
  CHARACTER(5),  PARAMETER :: fmt_name = "QEXML"
  CHARACTER(5),  PARAMETER :: fmt_version = "1.4.0"
  !
  LOGICAL,       PARAMETER :: rho_binary = .TRUE.
  !
  LOGICAL,       PARAMETER :: pot_binary = .TRUE.
  !
  CHARACTER(iotk_attlenx)  :: attr
  !
  !
  PUBLIC :: fmt_name, fmt_version
  PUBLIC :: current_fmt_version
  !
  PUBLIC :: rho_binary
  PUBLIC :: pot_binary
  !
  PUBLIC :: attr
  !
  PUBLIC :: create_directory, kpoint_dir, wfc_filename, copy_file,       &
            restart_dir, check_restartfile, check_file_exst,             &
            pp_check_file, save_history, save_print_counter,             &
            read_print_counter, set_kpoints_vars,                        &
            write_header, write_control, write_control_ph,               &
            write_status_ph, write_q,                                    &
            write_cell, write_ions, write_symmetry, write_planewaves,    &
            write_efield, write_spin, write_magnetization, write_xc,     &
            write_occ, write_bz,     &
            write_phonon, write_rho_xml, write_wfc, write_wfc_cmplx, write_eig, &
            read_wfc, read_rho_xml,  write_pot_xml, read_pot_xml
  !
       INTERFACE write_wfc
	    module procedure write_wfc_real, write_wfc_cmplx
       END INTERFACE       

       INTERFACE write_planewaves
	    module procedure write_planewaves_real, write_planewaves_cmplx
       END INTERFACE       
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE create_directory( dirname )
      !------------------------------------------------------------------------
      !
      USE wrappers,  ONLY : f_mkdir
      USE mp,        ONLY : mp_barrier
      USE mp_global, ONLY : me_image, intra_image_comm
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      !
      INTEGER                    :: ierr

      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      IF ( ionode ) ierr = f_mkdir( TRIM( dirname ) )
      !
      CALL mp_bcast ( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'create_directory', &
                   'unable to create directory ' // TRIM( dirname ), ierr )
      !
      ! ... syncronize all jobs (not sure it is really useful)
      !
      CALL mp_barrier( intra_image_comm )
      !
      ! ... check whether the scratch directory is writable
      !
      IF ( ionode ) THEN
         !
         OPEN( UNIT = 4, FILE = TRIM( dirname ) // '/test' // &
               TRIM( int_to_char( me_image ) ), IOSTAT = ierr )
         CLOSE( UNIT = 4, STATUS = 'DELETE' )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'create_directory:', &
                   TRIM( dirname ) // ' non existent or non writable', ierr )
      !
      RETURN
      !
    END SUBROUTINE create_directory
    !
    !------------------------------------------------------------------------
    FUNCTION kpoint_dir( basedir, ik )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=256)           :: kpoint_dir
      CHARACTER(LEN=*), INTENT(IN) :: basedir
      INTEGER,          INTENT(IN) :: ik
      !
      CHARACTER(LEN=256) :: kdirname
      CHARACTER(LEN=5)   :: kindex
      CHARACTER(LEN=6)   :: kindex1
      !
      IF (ik<99999) THEN
         WRITE( kindex, FMT = '( I5.5 )' ) ik     
         kdirname = TRIM( basedir ) // '/K' // kindex
      ELSEIF (ik<999999) THEN
         WRITE( kindex1, FMT = '( I6.6 )' ) ik     
         kdirname = TRIM( basedir ) // '/K' // kindex1
      ELSE
         call errore('kpoint_dir','ik too large, increase format',1)
      ENDIF
      !
      kpoint_dir = TRIM( kdirname )
      !
      RETURN
      !
    END FUNCTION kpoint_dir
    !
    !------------------------------------------------------------------------
    FUNCTION wfc_filename( basedir, name, ik, ipol, tag, extension, dir )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=256)                 :: wfc_filename
      CHARACTER(LEN=*),       INTENT(IN) :: basedir
      CHARACTER(LEN=*),       INTENT(IN) :: name
      INTEGER,                INTENT(IN) :: ik
      INTEGER,      OPTIONAL, INTENT(IN) :: ipol
      CHARACTER(*), OPTIONAL, INTENT(IN) :: tag
      CHARACTER(*), OPTIONAL, INTENT(IN) :: extension
      LOGICAL,      OPTIONAL, INTENT(IN) :: dir
      !    
      CHARACTER(LEN=256) :: filename, tag_, ext_
      LOGICAL :: dir_true
      !
      !
      filename = ''
      tag_     = ''
      ext_     = '.dat'
      dir_true = .true.
      !
      IF ( PRESENT( tag ) )         tag_ = '_'//TRIM(tag)
      IF ( PRESENT( extension ) )   ext_ = '.'//TRIM(extension)
      !
      IF ( PRESENT( ipol ) ) THEN
         !      
         WRITE( filename, FMT = '( I1 )' ) ipol
         !
      END IF
      IF ( PRESENT( dir )) dir_true=dir
      !
      IF (dir_true) THEN
         filename = TRIM( kpoint_dir( basedir, ik ) ) // '/' // &
                 & TRIM( name ) // TRIM( filename ) // TRIM( tag_ ) // TRIM( ext_)
      ELSE
         filename = TRIM( kpoint_dir( basedir, ik ) ) // '_' // &
                 & TRIM( name ) // TRIM( filename ) // TRIM( tag_ ) // TRIM( ext_)
      ENDIF
      !
      wfc_filename = TRIM( filename )
      !
      RETURN
      !
    END FUNCTION
    !
    !------------------------------------------------------------------------
    SUBROUTINE copy_file( file_in, file_out )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*), INTENT(IN) :: file_in, file_out
      !
      CHARACTER(LEN=256) :: string
      INTEGER            :: iun_in, iun_out, ierr
      !
      !
      IF ( .NOT. ionode ) RETURN
      !
      CALL iotk_free_unit( iun_in,  ierr )
      CALL iotk_free_unit( iun_out, ierr )
      !
      CALL errore( 'copy_file', 'no free units available', ierr )
      !
      OPEN( UNIT = iun_in,  FILE = file_in,  STATUS = "OLD" )
      OPEN( UNIT = iun_out, FILE = file_out, STATUS = "UNKNOWN" )         
      !
      copy_loop: DO
         !
         READ( UNIT = iun_in, FMT = '(A256)', IOSTAT = ierr ) string
         !
         IF ( ierr < 0 ) EXIT copy_loop
         !
         WRITE( UNIT = iun_out, FMT = '(A)' ) TRIM( string )
         !
      END DO copy_loop
      !
      CLOSE( UNIT = iun_in )
      CLOSE( UNIT = iun_out )
      !
      RETURN
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    FUNCTION restart_dir( outdir, runit )
      !------------------------------------------------------------------------
      !
      ! KNK_nimage
      ! USE mp_global, ONLY:  my_image_id
      CHARACTER(LEN=256)           :: restart_dir
      CHARACTER(LEN=*), INTENT(IN) :: outdir
      INTEGER,          INTENT(IN) :: runit
      !
      CHARACTER(LEN=256)         :: dirname
      INTEGER                    :: strlen
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      ! ... main restart directory
      !
      ! ... keep the line below ( this is the old style RESTARTXX ) !!!
      !
      ! dirname = 'RESTART' // int_to_char( runit )
      ! the next line is to have seperate RESTART for each image
      ! KNK_nimage
      ! if (my_image_id > 0) dirname = trim(dirname) // '_' // trim(int_to_char( my_image_id ))
      !
      dirname = TRIM( prefix ) // '_' // TRIM( int_to_char( runit ) )// '.save'
      !
      IF ( LEN( outdir ) > 1 ) THEN
         !
         strlen = INDEX( outdir, ' ' ) - 1
         !
         dirname = outdir(1:strlen) // '/' // dirname
         !
      END IF
      !
      restart_dir = TRIM( dirname )
      !
      RETURN
      !
    END FUNCTION restart_dir
    !
    !------------------------------------------------------------------------
    FUNCTION check_restartfile( outdir, ndr )
      !------------------------------------------------------------------------
      !
      USE io_global, ONLY : ionode, ionode_id
      USE mp_global, ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      LOGICAL                      :: check_restartfile
      INTEGER,          INTENT(IN) :: ndr
      CHARACTER(LEN=*), INTENT(IN) :: outdir
      CHARACTER(LEN=256)           :: filename
      LOGICAL                      :: lval
      !
      !
      filename = restart_dir( outdir, ndr )
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( filename ) // '/' // TRIM( xmlpun )
         !
         INQUIRE( FILE = TRIM( filename ), EXIST = lval )
         !
      END IF
      !
      CALL mp_bcast( lval, ionode_id, intra_image_comm )
      !
      check_restartfile = lval
      !
      RETURN
      !
    END FUNCTION check_restartfile
    !
    !------------------------------------------------------------------------
    FUNCTION check_file_exst( filename )
      !------------------------------------------------------------------------
      !
      USE io_global, ONLY : ionode, ionode_id
      USE mp_global, ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      LOGICAL          :: check_file_exst
      CHARACTER(LEN=*) :: filename
      !
      LOGICAL :: lexists
      !
      IF ( ionode ) THEN 
         !
         INQUIRE( FILE = TRIM( filename ), EXIST = lexists )
         !
      ENDIF
      !
      CALL mp_bcast ( lexists, ionode_id, intra_image_comm )
      !
      check_file_exst = lexists
      RETURN
      !
    END FUNCTION check_file_exst
    !
    !------------------------------------------------------------------------
    FUNCTION pp_check_file()
      !------------------------------------------------------------------------
      !
      USE io_global,         ONLY : ionode, ionode_id
      USE mp_global,         ONLY : intra_image_comm
      USE control_flags,     ONLY : lkpoint_dir, tqr
      !
      IMPLICIT NONE
      !
      LOGICAL            :: pp_check_file
      CHARACTER(LEN=256) :: dirname, filename
      INTEGER            :: ierr
      LOGICAL            :: lval, found, back_compat
      !
      !
      dirname  = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      filename = TRIM( dirname ) // '/' // TRIM( xmlpun )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = filename, IERR = ierr )
      !
      CALL mp_bcast ( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'pp_check_file', 'file ' // &
                 & TRIM( dirname ) // ' not found', ierr )

      !
      ! set a flag for back compatibility (before fmt v1.4.0)
      !
      back_compat = .FALSE.
      !
      IF ( TRIM( version_compare( current_fmt_version, "1.4.0" )) == "older") &
         back_compat = .TRUE.
      !
      IF ( ionode ) THEN
         !
         IF ( .NOT. back_compat ) THEN
             !
             CALL iotk_scan_begin( iunpun, "CONTROL" ) 
             !
         ENDIF
         !
         CALL iotk_scan_dat( iunpun, "PP_CHECK_FLAG", lval, FOUND = found)
         !
         IF ( .NOT. found ) lval = .FALSE. 
         !
         CALL iotk_scan_dat( iunpun, "LKPOINT_DIR", lkpoint_dir, FOUND = found)
         !
         IF ( .NOT. found ) lkpoint_dir = .TRUE. 
         !
         CALL iotk_scan_dat( iunpun, "Q_REAL_SPACE", tqr, FOUND = found)
         !
         IF ( .NOT. found ) tqr = .FALSE. 
         !
         !
         IF ( .NOT. back_compat ) THEN
             !
             CALL iotk_scan_end( iunpun, "CONTROL" ) 
             !
         ENDIF
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( lval, ionode_id, intra_image_comm )
      !
      CALL mp_bcast( lkpoint_dir, ionode_id, intra_image_comm )
      !
      CALL mp_bcast( tqr, ionode_id, intra_image_comm )
      !
      pp_check_file = lval
      !
      RETURN
      !
    END FUNCTION pp_check_file
    !
    !------------------------------------------------------------------------
    SUBROUTINE save_history( dirname, iter )
      !------------------------------------------------------------------------
      !
      ! ... a copy of the xml descriptor (data-file.xml) is saved in the 
      ! ... history subdir
      !
      USE io_files, ONLY : xmlpun_base
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      INTEGER,          INTENT(IN) :: iter
      !
#if defined (__VERBOSE_SAVE)
      !
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=6)   :: hindex
      !
      CALL create_directory( TRIM( dirname ) // '/history' )
      !
      WRITE( hindex, FMT = '(I6.6)' ) iter
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // '/history/' // &
                  & TRIM( xmlpun_base ) // hindex // '.xml'
         !
         CALL copy_file( TRIM( dirname ) // "/" // TRIM( xmlpun ), &
                         TRIM( filename ) )
         !
      END IF
      !
#endif
      !
      RETURN
      !
    END SUBROUTINE save_history
    !
    !------------------------------------------------------------------------
    SUBROUTINE save_print_counter( iter, outdir, wunit )
      !------------------------------------------------------------------------
      !
      ! ... a counter indicating the last successful printout iteration is saved
      !
      USE io_global, ONLY : ionode, ionode_id
      USE mp_global, ONLY : intra_image_comm
      USE mp,        ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN) :: iter
      CHARACTER(LEN=*), INTENT(IN) :: outdir
      INTEGER,          INTENT(IN) :: wunit
      !
      INTEGER            :: ierr
      CHARACTER(LEN=256) :: filename, dirname
      !
      !
      dirname = restart_dir( outdir, wunit )
      !
      CALL create_directory( TRIM( dirname ) )
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // '/print_counter.xml'
         !
         CALL iotk_open_write( iunpun, FILE = filename, &
                             & ROOT = "PRINT_COUNTER",  IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'save_print_counter', &
                   'cannot open restart file for writing', ierr )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_begin( iunpun, "LAST_SUCCESSFUL_PRINTOUT" )
         CALL iotk_write_dat(   iunpun, "STEP", iter )
         CALL iotk_write_end(   iunpun, "LAST_SUCCESSFUL_PRINTOUT" )
         !
         CALL iotk_close_write( iunpun )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE save_print_counter
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_print_counter( nprint_nfi, outdir, runit )
      !------------------------------------------------------------------------
      !
      ! ... the counter indicating the last successful printout iteration 
      ! ... is read here
      !
      USE io_global, ONLY : ionode, ionode_id
      USE mp_global, ONLY : intra_image_comm
      USE mp,        ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: nprint_nfi
      CHARACTER(LEN=*), INTENT(IN)  :: outdir
      INTEGER,          INTENT(IN)  :: runit
      !
      INTEGER            :: ierr
      CHARACTER(LEN=256) :: filename, dirname
      !
      !
      dirname = restart_dir( outdir, runit )
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // '/print_counter.xml'
         !
         CALL iotk_open_read( iunpun, FILE = filename, IERR = ierr )
         !
         IF ( ierr > 0 ) THEN
            !
            nprint_nfi = -1
            !
         ELSE
            !
            CALL iotk_scan_begin( iunpun, "LAST_SUCCESSFUL_PRINTOUT" )
            CALL iotk_scan_dat(   iunpun, "STEP", nprint_nfi )
            CALL iotk_scan_end(   iunpun, "LAST_SUCCESSFUL_PRINTOUT" )
            !
            CALL iotk_close_read( iunpun )
            !
         END IF
         !
      END IF
      !
      CALL mp_bcast( nprint_nfi, ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_print_counter   
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                                 npool, ikt, iks, ike, igwx, ipmask, ipsour )
      !------------------------------------------------------------------------
      !
      ! ... set working variables for k-point index (ikt) and 
      ! ... k-points number (nkt)
      !
      USE mp,         ONLY : mp_sum, mp_get, mp_max
      USE mp_global,  ONLY : me_image, nproc_image, me_pool, my_pool_id, &
                             nproc_pool, intra_pool_comm, root_pool, &
                             intra_image_comm
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)  :: ik, nk, kunit
      INTEGER, INTENT(IN)  :: ngwl, igl(:)
      INTEGER, INTENT(OUT) :: npool
      INTEGER, INTENT(OUT) :: ikt, iks, ike, igwx
      INTEGER, INTENT(OUT) :: ipmask(:), ipsour
      !
      INTEGER :: ierr, i
      INTEGER :: nkl, nkr, nkbl, nkt
      !
      !
      ikt = ik
      nkt = nk
      !
      ! ... find out the number of pools
      !
      npool = nproc_image / nproc_pool 
      !
      ! ... find out number of k points blocks
      !
      nkbl = nkt / kunit  
      !
      ! ... k points per pool
      !
      nkl = kunit * ( nkbl / npool )
      !
      ! ... find out the reminder
      !
      nkr = ( nkt - nkl * npool ) / kunit
      !
      ! ... Assign the reminder to the first nkr pools
      !
      IF ( my_pool_id < nkr ) nkl = nkl + kunit
      !
      ! ... find out the index of the first k point in this pool
      !
      iks = nkl * my_pool_id + 1
      !
      IF ( my_pool_id >= nkr ) iks = iks + nkr * kunit
      !
      ! ... find out the index of the last k point in this pool
      !
      ike = iks + nkl - 1
      !
      ipmask = 0
      ipsour = ionode_id
      !
      ! ... find out the index of the processor which collect the data 
      ! ... in the pool of ik
      !
      IF ( npool > 1 ) THEN
         !
         IF ( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
            !
            IF ( me_pool == root_pool ) ipmask( me_image + 1 ) = 1
            !
         END IF
         !
         ! ... Collect the mask for all proc in the image
         !
         CALL mp_sum( ipmask, intra_image_comm )
         !
         DO i = 1, nproc_image
            !
            IF( ipmask(i) == 1 ) ipsour = ( i - 1 )
            !
         END DO
         !
      END IF
      !
      igwx = 0
      ierr = 0
      !
      IF ( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
         !
         IF ( ngwl > SIZE( igl ) ) THEN
            !
            ierr = 1
            !
         ELSE
            !
            igwx = MAXVAL( igl(1:ngwl) )
            !
         END IF
         !
      END IF
      !
      ! ... get the maximum index within the pool
      !
      CALL mp_max( igwx, intra_pool_comm )
      !
      ! ... now notify all procs if an error has been found 
      !
      CALL mp_max( ierr, intra_image_comm )
      !
      CALL errore( 'set_kpoint_vars ', 'wrong size ngl', ierr )
      !
      IF ( ipsour /= ionode_id ) &
         CALL mp_get( igwx, igwx, me_image, ionode_id, ipsour, 1, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE set_kpoints_vars
    !
    !
    ! ... writing subroutines
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_header( creator_name, creator_version ) 
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: creator_name, creator_version


      CALL iotk_write_begin( iunpun, "HEADER" )
      !
      CALL iotk_write_attr(attr, "NAME",TRIM(fmt_name), FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "VERSION",TRIM(fmt_version) )
      CALL iotk_write_empty( iunpun, "FORMAT", ATTR=attr )
      !
      CALL iotk_write_attr(attr, "NAME",TRIM(creator_name), FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "VERSION",TRIM(creator_version) )
      CALL iotk_write_empty( iunpun, "CREATOR", ATTR=attr )
      !
      CALL iotk_write_end( iunpun, "HEADER" )
      !
    END SUBROUTINE write_header
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_control( pp_check_flag, lkpoint_dir, q_real_space, beta_real_space) 
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      LOGICAL, OPTIONAL, INTENT(IN) :: pp_check_flag, lkpoint_dir, q_real_space, beta_real_space


      CALL iotk_write_begin( iunpun, "CONTROL" )
      !
      !  This flag is used to check if the file can be used for post-processing
      IF ( PRESENT( pp_check_flag ) ) &
         CALL iotk_write_dat( iunpun, "PP_CHECK_FLAG", pp_check_flag )
      !
      !  This flag says how eigenvalues are saved
      IF ( PRESENT( lkpoint_dir ) ) &
         CALL iotk_write_dat( iunpun, "LKPOINT_DIR", lkpoint_dir )
      !
      !  This flag says if Q in real space has to be used
      IF ( PRESENT( q_real_space ) ) &
         CALL iotk_write_dat( iunpun, "Q_REAL_SPACE", q_real_space )
      ! This flag says if Beta functions were treated in real space
      IF ( PRESENT( beta_real_space ) ) &
         CALL iotk_write_dat( iunpun, "BETA_REAL_SPACE", beta_real_space )
      !
      CALL iotk_write_end( iunpun, "CONTROL" )
      !
    END SUBROUTINE write_control
    !

    SUBROUTINE write_control_ph( ldisp, epsil, trans, elph, zue, &
                      lraman, elop ) 
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: ldisp, epsil, trans, elph, zue, &
                      lraman, elop


      CALL iotk_write_begin( iunpun, "CONTROL" )
      !
      CALL iotk_write_dat( iunpun, "DISPERSION_RUN", ldisp )
      CALL iotk_write_dat( iunpun, "ELECTRIC_FIELD", epsil )
      CALL iotk_write_dat( iunpun, "PHONON_RUN", trans )
      CALL iotk_write_dat( iunpun, "ELECTRON_PHONON", elph )
      CALL iotk_write_dat( iunpun, "EFFECTIVE_CHARGE_PH", zue )
      CALL iotk_write_dat( iunpun, "RAMAN_TENSOR", lraman )
      CALL iotk_write_dat( iunpun, "ELECTRO_OPTIC", elop )
      !
      CALL iotk_write_end( iunpun, "CONTROL" )
      !
      RETURN
    END SUBROUTINE write_control_ph

    SUBROUTINE write_status_ph(current_iq, done_bands)
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: current_iq
      LOGICAL, INTENT(IN) :: done_bands

      CALL iotk_write_begin( iunpun, "STATUS_PH" )
      !
      CALL iotk_write_dat( iunpun, "DONE_BANDS", done_bands )
      CALL iotk_write_dat( iunpun, "CURRENT_Q", current_iq )
      !
      CALL iotk_write_end( iunpun, "STATUS_PH" )
      !
      RETURN
    END SUBROUTINE write_status_ph
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_cell( ibrav, symm_type, &
                           celldm, alat, a1, a2, a3, b1, b2, b3 )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(IN) :: ibrav
      CHARACTER(LEN=*), INTENT(IN) :: symm_type
      REAL(DP),         INTENT(IN) :: celldm(6), alat
      REAL(DP),         INTENT(IN) :: a1(3), a2(3), a3(3)
      REAL(DP),         INTENT(IN) :: b1(3), b2(3), b3(3)
      !
      CHARACTER(LEN=256) :: bravais_lattice
      !
      CALL iotk_write_begin( iunpun, "CELL" )
      !
      SELECT CASE ( ibrav )
        CASE(  0 )
           bravais_lattice = "free"
        CASE(  1 )
           bravais_lattice = "cubic P (sc)"
        CASE(  2 )
           bravais_lattice = "cubic F (fcc)"
        CASE(  3 )
           bravais_lattice = "cubic I (bcc)"
        CASE(  4 )
           bravais_lattice = "Hexagonal and Trigonal P"
        CASE(  5 )
           bravais_lattice = "Trigonal R"
        CASE(  6 )
           bravais_lattice = "Tetragonal P (st)"
        CASE(  7 )
           bravais_lattice = "Tetragonal I (bct)"
        CASE(  8 )
           bravais_lattice = "Orthorhombic P"
        CASE(  9 )
           bravais_lattice = "Orthorhombic base-centered(bco)"
        CASE( 10 )
           bravais_lattice = "Orthorhombic face-centered"
        CASE( 11 )
           bravais_lattice = "Orthorhombic body-centered"
        CASE( 12 )
           bravais_lattice = "Monoclinic P"
        CASE( 13 )
           bravais_lattice = "Monoclinic base-centered"
        CASE( 14 )
           bravais_lattice = "Triclinic P"
      END SELECT
      !
      CALL iotk_write_dat( iunpun, &
                           "BRAVAIS_LATTICE", TRIM( bravais_lattice ) )
      !
      CALL iotk_write_dat( iunpun, "CELL_SYMMETRY", symm_type )
      !
      CALL iotk_write_attr( attr, "UNITS", "Bohr", FIRST = .TRUE. )
      CALL iotk_write_dat( iunpun, "LATTICE_PARAMETER", alat, ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "CELL_DIMENSIONS", celldm(1:6) )
      !
      CALL iotk_write_begin( iunpun, "DIRECT_LATTICE_VECTORS" )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_DIRECT_LATTICE_VECTORS", ATTR=attr )
      CALL iotk_write_dat(   iunpun, "a1", a1(:) * alat, COLUMNS=3 )
      CALL iotk_write_dat(   iunpun, "a2", a2(:) * alat, COLUMNS=3 )
      CALL iotk_write_dat(   iunpun, "a3", a3(:) * alat, COLUMNS=3 )
      CALL iotk_write_end(   iunpun, "DIRECT_LATTICE_VECTORS" )
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      CALL iotk_write_begin( iunpun, "RECIPROCAL_LATTICE_VECTORS" )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_RECIPROCAL_LATTICE_VECTORS", ATTR=attr )
      CALL iotk_write_dat(   iunpun, "b1", b1(:), COLUMNS=3 )
      CALL iotk_write_dat(   iunpun, "b2", b2(:), COLUMNS=3 )
      CALL iotk_write_dat(   iunpun, "b3", b3(:), COLUMNS=3 )
      CALL iotk_write_end(   iunpun, "RECIPROCAL_LATTICE_VECTORS" )
      !
      CALL iotk_write_end( iunpun, "CELL" )
      !
    END SUBROUTINE write_cell
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_ions( nsp, nat, atm, ityp, psfile, &
                           pseudo_dir, amass, tau, if_pos, dirname, pos_unit )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(IN) :: nsp, nat
      INTEGER,          INTENT(IN) :: ityp(:)
      CHARACTER(LEN=*), INTENT(IN) :: atm(:)
      CHARACTER(LEN=*), INTENT(IN) :: psfile(:)
      CHARACTER(LEN=*), INTENT(IN) :: pseudo_dir
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      REAL(DP),         INTENT(IN) :: amass(:)
      REAL(DP),         INTENT(IN) :: tau(:,:)
      INTEGER,          INTENT(IN) :: if_pos(:,:)
      REAL(DP),         INTENT(IN) :: pos_unit
      !
      INTEGER            :: i, flen
      CHARACTER(LEN=256) :: file_pseudo
      !
      !
      CALL iotk_write_begin( iunpun, "IONS" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_ATOMS", nat )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_SPECIES", nsp )
      !
      flen = LEN_TRIM( pseudo_dir )
      !
      CALL iotk_write_attr ( attr, "UNITS", "a.m.u.", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_ATOMIC_MASSES", ATTR = attr )
      !
      DO i = 1, nsp
         !
         CALL iotk_write_begin( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
         !
         CALL iotk_write_dat( iunpun, "ATOM_TYPE", atm(i) )
         !
         IF ( pseudo_dir(flen:flen) /= '/' ) THEN
            !
            file_pseudo = pseudo_dir(1:flen) // '/' // psfile(i)
            !
         ELSE
            !
            file_pseudo = pseudo_dir(1:flen) // psfile(i)
            !
         END IF
         !
         IF (TRIM( file_pseudo ).ne. TRIM( dirname ) // "/" // &
                           TRIM(psfile(i))) &
         CALL copy_file( TRIM( file_pseudo ), &
                            TRIM( dirname ) // "/" // TRIM( psfile(i) ) )
         !
         CALL iotk_write_dat( iunpun, "MASS", amass(i) )
         !
         CALL iotk_write_dat( iunpun, "PSEUDO", TRIM( psfile(i) ) )
         !
         CALL iotk_write_end( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
         !
      ENDDO
      !
      ! BEWARE: the following instruction is part of a ugly hack to allow
      !         restarting in parallel execution in machines without a
      !         parallel file system - See read_ions in pw_restart.f90
      !
      CALL iotk_write_dat( iunpun, "PSEUDO_DIR", TRIM( pseudo_dir) )
      !
      CALL iotk_write_attr( attr, "UNITS", "Bohr", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_ATOMIC_POSITIONS", ATTR = attr )
      !
      DO i = 1, nat
         !
         CALL iotk_write_attr( attr, "SPECIES", &
                             & atm( ityp(i) ), FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "INDEX",  ityp(i) )                     
         CALL iotk_write_attr( attr, "tau",    tau(:,i)*pos_unit )
         CALL iotk_write_attr( attr, "if_pos", if_pos(:,i) )
         CALL iotk_write_empty( iunpun, &
                              & "ATOM" // TRIM( iotk_index( i ) ), attr )
         !
      END DO
      !
      CALL iotk_write_end( iunpun, "IONS" )
      !
    END SUBROUTINE write_ions
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_symmetry( ibrav, symm_type, nrot, nsym, invsym, noinv, &
                               nr1, nr2, nr3, ftau, s, sname, irt, nat, t_rev )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(IN) :: ibrav, nrot, nsym,  nr1, nr2, nr3
      CHARACTER(LEN=*), INTENT(IN) :: symm_type
      LOGICAL,          INTENT(IN) :: invsym, noinv
      INTEGER,          INTENT(IN) :: s(:,:,:), ftau(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: sname(:)
      INTEGER,          INTENT(IN) :: irt(:,:), nat, t_rev(:)
      !
      INTEGER  :: i
      REAL(DP) :: tmp(3)
      !
      !
      CALL iotk_write_begin( iunpun, "SYMMETRIES" )
      !
      IF ( ibrav == 0 ) &
         CALL iotk_write_dat( iunpun, "CELL_SYMMETRY", symm_type )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_SYMMETRIES", nsym )
      CALL iotk_write_dat( iunpun, "NUMBER_OF_BRAVAIS_SYMMETRIES", nrot )
      !
      CALL iotk_write_dat( iunpun, "INVERSION_SYMMETRY", invsym )
      !
      CALL iotk_write_dat( iunpun, "DO_NOT_USE_TIME_REVERSAL", noinv )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_ATOMS", nat )
      !
      CALL iotk_write_attr( attr, "UNITS", "Crystal", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_SYMMETRIES", ATTR = attr )
      !
      DO i = 1, nsym
         !
         CALL iotk_write_begin( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
         !
         CALL iotk_write_attr ( attr, "NAME", TRIM( sname(i) ), FIRST=.TRUE. )
         CALL iotk_write_attr ( attr, "T_REV", t_rev(i) )
         CALL iotk_write_empty( iunpun, "INFO", ATTR = attr )
         !
         tmp(1) = ftau(1,i) / DBLE( nr1 )
         tmp(2) = ftau(2,i) / DBLE( nr2 )
         tmp(3) = ftau(3,i) / DBLE( nr3 )
         !
         CALL iotk_write_dat( iunpun, "ROTATION", s(:,:,i), COLUMNS=3 )
         CALL iotk_write_dat( iunpun, "FRACTIONAL_TRANSLATION", tmp(1:3), COLUMNS=3 )
         CALL iotk_write_dat( iunpun, "EQUIVALENT_IONS", irt(i,1:nat), COLUMNS=8 )
         !
         CALL iotk_write_end( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
         !
      ENDDO
      !
      ! ... the following are the symmetries of the Bravais lattice alone
      ! ... (they may be more than crystal, i.e. basis+lattice, symmetries)
      !
      DO i = nsym+1, nrot
         !
         CALL iotk_write_begin( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
         !
         CALL iotk_write_attr ( attr, "NAME", TRIM( sname(i) ), FIRST=.TRUE. )
         CALL iotk_write_empty( iunpun, "INFO", ATTR = attr )
         CALL iotk_write_dat( iunpun, "ROTATION", s(:,:,i), COLUMNS=3 )
         !
         CALL iotk_write_end( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
         !
      ENDDO
      !
      CALL iotk_write_end( iunpun, "SYMMETRIES" )
      !
    END SUBROUTINE write_symmetry
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_efield( tefield, dipfield, edir, emaxpos, eopreg, eamp )
      !------------------------------------------------------------------------
      !
      LOGICAL, INTENT(IN) :: &
           tefield,      &! if .TRUE. a finite electric field is added to the
                          ! local potential
           dipfield       ! if .TRUE. the dipole field is subtracted
      INTEGER, INTENT(IN) :: &
           edir           ! direction of the field
      REAL(DP), INTENT(IN) :: &
           emaxpos,  &! position of the maximum of the field (0<emaxpos<1)
           eopreg,   &! amplitude of the inverse region (0<eopreg<1)
           eamp       ! field amplitude (in a.u.) (1 a.u. = 51.44 10^11 V/m)
      !
      !
      CALL iotk_write_begin( iunpun, "ELECTRIC_FIELD" )
      !
      CALL iotk_write_dat( iunpun, "HAS_ELECTRIC_FIELD", tefield )
      !
      CALL iotk_write_dat( iunpun, "HAS_DIPOLE_CORRECTION", dipfield )
      !
      CALL iotk_write_dat( iunpun, "FIELD_DIRECTION", edir )
      !
      CALL iotk_write_dat( iunpun, "MAXIMUM_POSITION", emaxpos )
      !
      CALL iotk_write_dat( iunpun, "INVERSE_REGION", eopreg )
      !
      CALL iotk_write_dat( iunpun, "FIELD_AMPLITUDE", eamp )
      !
      CALL iotk_write_end( iunpun, "ELECTRIC_FIELD" )
      !
    END SUBROUTINE write_efield
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_planewaves_real( ecutwfc, dual, npwx, gamma_only, nr1, nr2, &
                                 nr3, ngm_g, nr1s, nr2s, nr3s, ngms_g, nr1b, &
                                 nr2b, nr3b, itmp, lgvec )
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : e2
      !
      INTEGER,  INTENT(IN) :: npwx, nr1, nr2, nr3, ngm_g, &
                              nr1s, nr2s, nr3s, ngms_g, nr1b, nr2b, nr3b
      INTEGER,  INTENT(IN) :: itmp(:,:)
      REAL(DP), INTENT(IN) :: ecutwfc, dual
      LOGICAL,  INTENT(IN) :: gamma_only, lgvec
      !
      !
      CALL iotk_write_begin( iunpun, "PLANE_WAVES" )
      !
      CALL iotk_write_attr ( attr, "UNITS", "Hartree", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_CUTOFF", ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "WFC_CUTOFF", ecutwfc / e2 )
      !
      CALL iotk_write_dat( iunpun, "RHO_CUTOFF", ecutwfc * dual / e2 )
      !
      CALL iotk_write_dat( iunpun, "MAX_NUMBER_OF_GK-VECTORS", npwx )
      !
      CALL iotk_write_dat( iunpun, "GAMMA_ONLY", gamma_only )
      !
      CALL iotk_write_attr( attr, "nr1", nr1, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2", nr2 )
      CALL iotk_write_attr( attr, "nr3", nr3 )
      CALL iotk_write_empty( iunpun, "FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "GVECT_NUMBER", ngm_g )
      !
      CALL iotk_write_attr( attr, "nr1s", nr1s, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2s", nr2s )
      CALL iotk_write_attr( attr, "nr3s", nr3s )
      CALL iotk_write_empty( iunpun, "SMOOTH_FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "SMOOTH_GVECT_NUMBER", ngms_g )
      !
      IF ( lgvec ) THEN
         !
         ! ... write the G-vectors
         !
         CALL iotk_link( iunpun, "G-VECTORS", &
                         "./gvectors.dat", CREATE = .TRUE., BINARY = .TRUE. )
         !
         CALL iotk_write_begin( iunpun, "G-VECTORS" )
         !
         CALL iotk_write_attr( attr, "nr1s", nr1s, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "nr2s", nr2s )
         CALL iotk_write_attr( attr, "nr3s", nr3s )
         CALL iotk_write_attr( attr, "gvect_number", ngm_g )
         CALL iotk_write_attr( attr, "gamma_only", gamma_only )
         CALL iotk_write_attr( attr, "units", "crystal" )
         CALL iotk_write_empty( iunpun, "INFO", ATTR = attr )
         !
         CALL iotk_write_dat  ( iunpun, "g", itmp(1:3,1:ngm_g), COLUMNS = 3 )
         CALL iotk_write_end  ( iunpun, "G-VECTORS" )
         !
      END IF
      !
      CALL iotk_write_attr( attr, "nr1b", nr1b , FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2b", nr2b )
      CALL iotk_write_attr( attr, "nr3b", nr3b )
      CALL iotk_write_empty( iunpun, "SMALLBOX_FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_end( iunpun, "PLANE_WAVES" )
      !
    END SUBROUTINE write_planewaves_real
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_planewaves_cmplx( ecutwfc, dual, npwx, do_wf_cmplx, gamma_only, nr1, nr2, & !added:giovanni do_wf_cmplx
                                 nr3, ngm_g, nr1s, nr2s, nr3s, ngms_g, nr1b,            &
                                 nr2b, nr3b, itmp, lgvec )
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : e2
      !
      INTEGER,  INTENT(IN) :: npwx, nr1, nr2, nr3, ngm_g, &
                              nr1s, nr2s, nr3s, ngms_g, nr1b, nr2b, nr3b
      INTEGER,  INTENT(IN) :: itmp(:,:)
      REAL(DP), INTENT(IN) :: ecutwfc, dual
      LOGICAL,  INTENT(IN) :: do_wf_cmplx, gamma_only, lgvec !added:giovanni do_wf_cmplx
      !
      !
      CALL iotk_write_begin( iunpun, "PLANE_WAVES" )
      !
      CALL iotk_write_attr ( attr, "UNITS", "Hartree", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_CUTOFF", ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "WFC_CUTOFF", ecutwfc / e2 )
      !
      CALL iotk_write_dat( iunpun, "RHO_CUTOFF", ecutwfc * dual / e2 )
      !
      CALL iotk_write_dat( iunpun, "MAX_NUMBER_OF_GK-VECTORS", npwx )
      !
      CALL iotk_write_dat( iunpun, "DO_WF_CMPLX", do_wf_cmplx ) !added:giovanni
      CALL iotk_write_dat( iunpun, "GAMMA_ONLY", gamma_only.and..not.do_wf_cmplx )
      !
      CALL iotk_write_attr( attr, "nr1", nr1, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2", nr2 )
      CALL iotk_write_attr( attr, "nr3", nr3 )
      CALL iotk_write_empty( iunpun, "FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "GVECT_NUMBER", ngm_g )
      !
      CALL iotk_write_attr( attr, "nr1s", nr1s, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2s", nr2s )
      CALL iotk_write_attr( attr, "nr3s", nr3s )
      CALL iotk_write_empty( iunpun, "SMOOTH_FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "SMOOTH_GVECT_NUMBER", ngms_g )
      !
      IF ( lgvec ) THEN
         !
         ! ... write the G-vectors
         !
         CALL iotk_link( iunpun, "G-VECTORS", &
                         "./gvectors.dat", CREATE = .TRUE., BINARY = .TRUE. )
         !
         CALL iotk_write_begin( iunpun, "G-VECTORS" )
         !
         CALL iotk_write_attr( attr, "nr1s", nr1s, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "nr2s", nr2s )
         CALL iotk_write_attr( attr, "nr3s", nr3s )
         CALL iotk_write_attr( attr, "gvect_number", ngm_g )
         CALL iotk_write_attr( attr, "do_wf_cmplx", do_wf_cmplx ) !added:giovanni
         CALL iotk_write_attr( attr, "gamma_only", gamma_only.and..not.do_wf_cmplx ) !modified:giovanni
         CALL iotk_write_attr( attr, "units", "crystal" )
         CALL iotk_write_empty( iunpun, "INFO", ATTR = attr )
         !
         CALL iotk_write_dat  ( iunpun, "g", itmp(1:3,1:ngm_g), COLUMNS = 3 )
         CALL iotk_write_end  ( iunpun, "G-VECTORS" )
         !
      END IF
      !
      CALL iotk_write_attr( attr, "nr1b", nr1b , FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2b", nr2b )
      CALL iotk_write_attr( attr, "nr3b", nr3b )
      CALL iotk_write_empty( iunpun, "SMALLBOX_FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_end( iunpun, "PLANE_WAVES" )
      !
    END SUBROUTINE write_planewaves_cmplx
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_spin( lsda, noncolin, npol, lspinorb, domag )
      !------------------------------------------------------------------------
      !
      LOGICAL, INTENT(IN) :: lsda, noncolin, lspinorb, domag
      INTEGER, INTENT(IN) :: npol
      !
      !
      CALL iotk_write_begin( iunpun, "SPIN" )
      !
      CALL iotk_write_dat( iunpun, "LSDA", lsda )
      !
      CALL iotk_write_dat( iunpun, "NON-COLINEAR_CALCULATION", noncolin )
      !
      IF ( noncolin ) &
         CALL iotk_write_dat( iunpun, "SPINOR_DIM", npol )
      !
      CALL iotk_write_dat( iunpun, "SPIN-ORBIT_CALCULATION", lspinorb )
      CALL iotk_write_dat( iunpun, "SPIN-ORBIT_DOMAG", domag )
      !
      CALL iotk_write_end( iunpun, "SPIN" )
      !
    END SUBROUTINE write_spin
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_magnetization(starting_magnetization, angle1, angle2, &
                                   nsp, two_fermi_energies, i_cons, mcons, bfield, &
                                   ef_up, ef_dw, nelup, neldw, lambda)
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : PI
      !
      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: nsp, i_cons
      REAL(DP), INTENT(IN) :: starting_magnetization(nsp), &
                              angle1(nsp), angle2(nsp), mcons(3,nsp), &
                              bfield(3), ef_up, ef_dw, nelup, neldw, lambda
      LOGICAL,  INTENT(IN) :: two_fermi_energies
      !
      INTEGER :: i
      !
      CALL iotk_write_begin( iunpun, "MAGNETIZATION_INIT" )

      CALL iotk_write_dat(iunpun,"CONSTRAINT_MAG", i_cons)

      CALL iotk_write_dat( iunpun, "NUMBER_OF_SPECIES", nsp ) 

      DO i = 1, nsp
         !
         CALL iotk_write_begin( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
         !
         CALL iotk_write_dat( iunpun, "STARTING_MAGNETIZATION",  &
                                       starting_magnetization(i) )
         CALL iotk_write_dat( iunpun, "ANGLE1", &
                                       angle1(i)*180.0_DP/PI )
         CALL iotk_write_dat( iunpun, "ANGLE2", &
                                       angle2(i)*180.0_DP/PI )
         IF (i_cons==1.OR.i_cons==2) THEN
            CALL iotk_write_dat( iunpun, "CONSTRANT_1", mcons(1,i) )
            CALL iotk_write_dat( iunpun, "CONSTRANT_2", mcons(2,i) )
            CALL iotk_write_dat( iunpun, "CONSTRANT_3", mcons(3,i) )
         ENDIF
         !
         CALL iotk_write_end( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
         !
      ENDDO

      IF (i_cons==3) THEN
         !
         CALL iotk_write_dat( iunpun, "FIXED_MAGNETIZATION_1", mcons(1,1) )
         CALL iotk_write_dat( iunpun, "FIXED_MAGNETIZATION_2", mcons(2,1) )
         CALL iotk_write_dat( iunpun, "FIXED_MAGNETIZATION_3", mcons(3,1) )
         !
      ELSE IF (i_cons==4) THEN
         !
         CALL iotk_write_dat( iunpun, "MAGNETIC_FIELD_1", bfield(1) )
         CALL iotk_write_dat( iunpun, "MAGNETIC_FIELD_2", bfield(2) )
         CALL iotk_write_dat( iunpun, "MAGNETIC_FIELD_3", bfield(3) )
         !
      ENDIF
      !
      CALL iotk_write_dat(iunpun,"TWO_FERMI_ENERGIES",two_fermi_energies)
      !
      IF (two_fermi_energies) THEN
         !
         CALL iotk_write_dat( iunpun, "FIXED_MAGNETIZATION", mcons(3,1) )
         CALL iotk_write_dat( iunpun, "ELECTRONS_UP", nelup )
         CALL iotk_write_dat( iunpun, "ELECTRONS_DOWN", neldw )
         CALL iotk_write_dat( iunpun, "FERMI_ENERGY_UP", ef_up )
         CALL iotk_write_dat( iunpun, "FERMI_ENERGY_DOWN", ef_dw )
         !
      ENDIF
      !
      IF (i_cons>0) CALL iotk_write_dat(iunpun,"LAMBDA",lambda)
      !
      CALL iotk_write_end( iunpun, "MAGNETIZATION_INIT" )
      !
    RETURN
    !
    END SUBROUTINE write_magnetization
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_xc( dft, nsp, lda_plus_u, &
                         Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_alpha )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*),   INTENT(IN) :: dft
      LOGICAL,            INTENT(IN) :: lda_plus_u
      INTEGER,  OPTIONAL, INTENT(IN) :: nsp
      INTEGER,  OPTIONAL, INTENT(IN) :: Hubbard_lmax
      INTEGER,  OPTIONAL, INTENT(IN) :: Hubbard_l(:)
      REAL(DP), OPTIONAL, INTENT(IN) :: Hubbard_U(:), Hubbard_alpha(:)
      !
      !
      CALL iotk_write_begin( iunpun, "EXCHANGE_CORRELATION" )
      !
      CALL iotk_write_dat( iunpun, "DFT", dft )
      !
      CALL iotk_write_dat( iunpun, "LDA_PLUS_U_CALCULATION", lda_plus_u )
      !
      IF ( lda_plus_u ) THEN
         !
         IF ( .NOT. PRESENT( Hubbard_lmax ) .OR. &
              .NOT. PRESENT( Hubbard_l )    .OR. & 
              .NOT. PRESENT( Hubbard_U )    .OR. &
              .NOT. PRESENT( nsp )          .OR. &
              .NOT. PRESENT( Hubbard_alpha ) ) &
            CALL errore( 'write_exchange_correlation', &
                         ' variables for LDA+U not present', 1 )
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_SPECIES", nsp )
         !
         CALL iotk_write_dat( iunpun, "HUBBARD_LMAX", Hubbard_lmax )
         !
         CALL iotk_write_dat( iunpun, "HUBBARD_L", &
                              Hubbard_l(1:nsp) )
         !
         CALL iotk_write_dat( iunpun, "HUBBARD_U", Hubbard_U(1:nsp) )
         !
         CALL iotk_write_dat( iunpun, "HUBBARD_ALPHA", Hubbard_alpha(1:nsp) )
         !
      END IF
      !
      CALL iotk_write_end( iunpun, "EXCHANGE_CORRELATION" )
      !
    END SUBROUTINE write_xc
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_occ( lgauss, ngauss, degauss, ltetra, ntetra, &
                          tetra, tfixed_occ, lsda, nstates_up, nstates_down, f_inp )
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : e2
      !
      LOGICAL,            INTENT(IN) :: lgauss, ltetra, tfixed_occ, lsda
      INTEGER,  OPTIONAL, INTENT(IN) :: ngauss, ntetra, nstates_up, nstates_down
      INTEGER,  OPTIONAL, INTENT(IN) :: tetra(:,:)
      REAL(DP), OPTIONAL, INTENT(IN) :: degauss, f_inp(:,:)      
      !
      INTEGER :: i
      !
      !
      CALL iotk_write_begin( iunpun, "OCCUPATIONS" )
      !
      CALL iotk_write_dat( iunpun, "SMEARING_METHOD", lgauss )
      !
      IF ( lgauss ) THEN
         !
         CALL iotk_write_dat( iunpun, "SMEARING_TYPE", ngauss )
         !
         CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
         !
         CALL iotk_write_dat( iunpun, "SMEARING_PARAMETER", &
                              degauss / e2, ATTR = attr )
         !
      END IF
      !
      CALL iotk_write_dat( iunpun, "TETRAHEDRON_METHOD", ltetra )
      !
      IF ( ltetra ) THEN
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_TETRAHEDRA", ntetra )
         !
         DO i = 1, ntetra
            !
            CALL iotk_write_dat( iunpun, "TETRAHEDRON" // &
                               & iotk_index( i ), tetra(1:4,i) )
            !
         END DO
         !
      END IF
      !
      CALL iotk_write_dat( iunpun, "FIXED_OCCUPATIONS", tfixed_occ )
      !
      IF ( tfixed_occ ) THEN
         !
         CALL iotk_write_attr( attr, "lsda" , lsda, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "nstates_up", nstates_up )
         CALL iotk_write_attr( attr, "nstates_down", nstates_down )
         !
         CALL iotk_write_empty( iunpun, 'INFO', ATTR = attr )
         !
         CALL iotk_write_dat( iunpun, "INPUT_OCC_UP", f_inp(1:nstates_up,1) )
         !
         IF ( lsda ) &
            CALL iotk_write_dat( iunpun, "INPUT_OCC_DOWN", f_inp(1:nstates_down,2) )
         !
      END IF
      !
      CALL iotk_write_end( iunpun, "OCCUPATIONS" )
      !
    END SUBROUTINE write_occ
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_bz( num_k_points, xk, wk, k1, k2, k3, nk1, nk2, nk3, &
                         qnorm, nks_start, xk_start, wk_start )
      !------------------------------------------------------------------------
      !
      INTEGER,  INTENT(IN) :: num_k_points, k1, k2, k3, nk1, nk2, nk3
      REAL(DP), INTENT(IN) :: xk(:,:), wk(:)
      REAL(DP), INTENT(IN) :: qnorm
      INTEGER,  INTENT(IN), OPTIONAL ::  nks_start
      REAL(DP), INTENT(IN), OPTIONAL :: xk_start(:,:), wk_start(:)
      !
      INTEGER :: ik, i
      !
      !
      CALL iotk_write_begin( iunpun, "BRILLOUIN_ZONE" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_K-POINTS", num_k_points )
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_K-POINTS", attr )
      !
      CALL iotk_write_attr( attr, "nk1", nk1, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nk2", nk2 )
      CALL iotk_write_attr( attr, "nk3", nk3 )
      CALL iotk_write_empty( iunpun, "MONKHORST_PACK_GRID", attr )
      CALL iotk_write_attr( attr, "k1", k1, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "k2", k2 )
      CALL iotk_write_attr( attr, "k3", k3 )
      CALL iotk_write_empty( iunpun, "MONKHORST_PACK_OFFSET", attr )
      !
      DO ik = 1, num_k_points
         !
         CALL iotk_write_attr( attr, "XYZ", xk(:,ik), FIRST = .TRUE. )
         !            
         CALL iotk_write_attr( attr, "WEIGHT", wk(ik) )
         !
         CALL iotk_write_empty( iunpun, "K-POINT" // &
                              & TRIM( iotk_index(ik) ), attr )
         !
      END DO
      !
      ! ... these are k-points and weights in the Irreducible BZ
      !
      IF (present(nks_start).and.present(xk_start).and.present(wk_start)) THEN
         CALL iotk_write_dat( iunpun, "STARTING_K-POINTS", nks_start )
         !
         DO ik = 1, nks_start
            !
            CALL iotk_write_attr( attr, "XYZ", xk_start(:,ik), FIRST = .TRUE. )
            !            
            CALL iotk_write_attr( attr, "WEIGHT", wk_start(ik) )
            !
            CALL iotk_write_empty( iunpun, "K-POINT_START" // &
                              & TRIM( iotk_index(ik) ), attr )
            !
         END DO
      ENDIF
      !
      CALL iotk_write_dat( iunpun, "NORM-OF-Q", qnorm )
      !
      CALL iotk_write_end( iunpun, "BRILLOUIN_ZONE" )
      !
    END SUBROUTINE write_bz
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_phonon( modenum, xqq )
      !------------------------------------------------------------------------
      !
      INTEGER,  INTENT(IN) :: modenum
      REAL(DP), INTENT(IN) :: xqq(:)
      !
      !
      CALL iotk_write_begin( iunpun, "PHONON" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_MODES", modenum )
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_Q-POINT", attr )
      !
      CALL iotk_write_dat( iunpun, "Q-POINT", xqq(:), COLUMNS=3 )
      !
      CALL iotk_write_end( iunpun, "PHONON" )
      !
    END SUBROUTINE write_phonon

    SUBROUTINE write_q( nqs, x_q, done_iq )
      !------------------------------------------------------------------------
      !
      INTEGER, INTENT(IN) :: nqs
      REAL(DP), INTENT(IN) :: x_q(3,nqs)
      INTEGER, INTENT(IN) :: done_iq(nqs)
      !
      CALL iotk_write_begin( iunpun, "Q_POINTS" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_Q_POINTS", nqs  )
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      !
      CALL iotk_write_empty( iunpun, "UNITS_FOR_Q-POINT", attr )
      !
      CALL iotk_write_dat( iunpun, "Q-POINT_COORDINATES", x_q(:,:), COLUMNS=3 )
      !
      CALL iotk_write_dat( iunpun, "Q-POINT_DONE", done_iq(:) )
      !
      CALL iotk_write_end( iunpun, "Q_POINTS" )
      !
      RETURN
    END SUBROUTINE write_q
    !
    ! ... methods to write and read effective_potential
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_pot_xml( pot_file_base, pot, &
                              nr1, nr2, nr3, nr1x, nr2x, ipp, npp, &
                              ionode, intra_group_comm, inter_group_comm )
      !------------------------------------------------------------------------
      !
      ! ... Writes effective-potential pot, one plane at a time.
      ! ... If ipp and npp are specified, planes are collected one by one from
      ! ... all processors, avoiding an overall collect of the effective-potential
      ! ... on a single proc.
      !
      USE mp,        ONLY : mp_get, mp_sum, mp_rank, mp_size
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*),  INTENT(IN) :: pot_file_base
      REAL(DP),          INTENT(IN) :: pot(:)
      INTEGER,           INTENT(IN) :: nr1, nr2, nr3
      INTEGER,           INTENT(IN) :: nr1x, nr2x
      INTEGER,           INTENT(IN) :: ipp(:)
      INTEGER,           INTENT(IN) :: npp(:)
      LOGICAL,           INTENT(IN) :: ionode
      INTEGER,           INTENT(IN) :: intra_group_comm, inter_group_comm
      !
      INTEGER               :: ierr, i, j, k, kk, ldr, ip
      CHARACTER(LEN=256)    :: pot_file
      CHARACTER(LEN=10)     :: pot_extension
      REAL(DP), ALLOCATABLE :: pot_plane(:)
      INTEGER,  ALLOCATABLE :: kowner(:)
      INTEGER               :: my_group_id, me_group, nproc_group, io_group_id, io_group
      !
      INTEGER,    PARAMETER :: potunit = 19
      !
      me_group    = mp_rank( intra_group_comm )
      nproc_group = mp_size( intra_group_comm )
      my_group_id = mp_rank( inter_group_comm )
      !
      pot_extension = '.dat'
      IF ( .NOT. pot_binary ) pot_extension = '.xml'
      !
      pot_file = TRIM( pot_file_base ) // TRIM( pot_extension )
      !
      IF ( ionode ) THEN 
         CALL iotk_open_write( potunit, FILE = pot_file,  BINARY = pot_binary, IERR = ierr )
         CALL errore( 'write_pot_xml', 'cannot open' // TRIM( pot_file ) // ' file for writing', ierr )
      END IF 
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_begin( potunit, "EFFECTIVE-POTENTIAL" )
         !
         CALL iotk_write_attr( attr, "nr1", nr1, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "nr2", nr2 )
         CALL iotk_write_attr( attr, "nr3", nr3 )
         !
         CALL iotk_write_empty( potunit, "INFO", attr )
         !
      END IF
      !
      ALLOCATE( pot_plane( nr1*nr2 ) )
      ALLOCATE( kowner( nr3 ) )
      !
      ! ... find the index of the group (pool) that will write potential
      !
      io_group_id = 0
      !
      IF ( ionode ) io_group_id = my_group_id
      !
      CALL mp_sum( io_group_id, intra_group_comm )
      CALL mp_sum( io_group_id, inter_group_comm )
      !
      ! ... find the index of the ionode within its own group (pool)
      !
      io_group = 0
      !
      IF ( ionode ) io_group = me_group
      !
      CALL mp_sum( io_group, intra_group_comm )
      !
      ! ... find out the owner of each "z" plane
      !
      DO ip = 1, nproc_group
         !
         kowner( (ipp(ip)+1):(ipp(ip)+npp(ip)) ) = ip - 1
         !
      END DO
      !
      ldr = nr1x*nr2x
      !
      DO k = 1, nr3
         !
         !  Only one subgroup write the effective-potential
         !
         IF( ( kowner(k) == me_group ) .AND. ( my_group_id == io_group_id ) ) THEN
            !
            kk = k - ipp( me_group + 1 )
            ! 
            DO j = 1, nr2
               !
               DO i = 1, nr1
                  !
                  pot_plane(i+(j-1)*nr1) = pot(i+(j-1)*nr1x+(kk-1)*ldr)
                  !
               END DO
               !
            END DO
            !
         END IF
         !
         IF ( kowner(k) /= io_group .AND. my_group_id == io_group_id ) &
            CALL mp_get( pot_plane, pot_plane, me_group, io_group, kowner(k), k, intra_group_comm )
         !
         IF ( ionode ) &
            CALL iotk_write_dat( potunit, "z" // iotk_index( k ), pot_plane )
         !
      END DO
      !
      DEALLOCATE( pot_plane )
      DEALLOCATE( kowner )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_end( potunit, "EFFECTIVE-POTENTIAL" )
         !
         CALL iotk_close_write( potunit )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE write_pot_xml
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_pot_xml( pot_file_base, pot, &
                             nr1, nr2, nr3, nr1x, nr2x, ipp, npp, &
                             ionode, intra_group_comm, inter_group_comm )
      !------------------------------------------------------------------------
      !
      ! ... Reads effective-potential pot, one plane at a time.
      ! ... If ipp and npp are specified, planes are collected one by one from
      ! ... all processors, avoiding an overall collect of the effective-potential
      ! ... on a single proc.
      !
      USE mp,        ONLY : mp_put, mp_sum, mp_rank, mp_size
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*),  INTENT(IN)  :: pot_file_base
      INTEGER,           INTENT(IN)  :: nr1, nr2, nr3
      INTEGER,           INTENT(IN)  :: nr1x, nr2x
      REAL(DP),          INTENT(OUT) :: pot(:)
      INTEGER, OPTIONAL, INTENT(IN)  :: ipp(:)
      INTEGER, OPTIONAL, INTENT(IN)  :: npp(:)
      LOGICAL,           INTENT(IN)  :: ionode
      INTEGER,           INTENT(IN)  :: intra_group_comm, inter_group_comm
      !
      INTEGER               :: ierr, i, j, k, kk, ldr, ip
      INTEGER               :: nr( 3 )
      CHARACTER(LEN=256)    :: pot_file
      REAL(DP), ALLOCATABLE :: pot_plane(:)
      INTEGER,  ALLOCATABLE :: kowner(:)
      LOGICAL               :: exst
      INTEGER               :: ngroup, my_group_id, me_group, nproc_group, io_group_id, io_group
      !
      INTEGER,   PARAMETER  :: potunit = 19      
      !
      me_group    = mp_rank( intra_group_comm )
      nproc_group = mp_size( intra_group_comm )
      my_group_id = mp_rank( inter_group_comm )
      ngroup      = mp_size( inter_group_comm )
      !
      pot_file = TRIM( pot_file_base ) // ".dat"
      exst = check_file_exst( TRIM(pot_file) ) 
      !
      IF ( .NOT. exst ) THEN
          !
          pot_file = TRIM( pot_file_base ) // ".xml"
          exst = check_file_exst( TRIM(pot_file) ) 
          !
      ENDIF
      ! 
      IF ( .NOT. exst ) CALL errore('read_pot_xml', 'searching for '//TRIM(pot_file), 10)
      !
      IF ( ionode ) THEN
         CALL iotk_open_read( potunit, FILE = pot_file, IERR = ierr )
         CALL errore( 'read_pot_xml', 'cannot open ' // TRIM( pot_file ) // ' file for reading', ierr )
      END IF
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( potunit, "EFFECTIVE-POTENTIAL" )
         !
         CALL iotk_scan_empty( potunit, "INFO", attr )
         !
         CALL iotk_scan_attr( attr, "nr1", nr(1) )
         CALL iotk_scan_attr( attr, "nr2", nr(2) )
         CALL iotk_scan_attr( attr, "nr3", nr(3) )
         !
         IF ( nr1 /= nr(1) .OR. nr2 /= nr(2) .OR. nr3 /= nr(3) ) &
            CALL errore( 'read_pot_xml', 'dimensions do not match', 1 )
         !
      END IF
      !
      ALLOCATE( pot_plane( nr1*nr2 ) )
      ALLOCATE( kowner( nr3 ) )
      !
      ! ... find the index of the pool that will write pot
      !
      io_group_id = 0
      !
      IF ( ionode ) io_group_id = my_group_id
      !
      CALL mp_sum( io_group_id, intra_group_comm )
      CALL mp_sum( io_group_id, inter_group_comm )
      !
      ! ... find the index of the ionode within its own pool
      !
      io_group = 0
      !
      IF ( ionode ) io_group = me_group
      !
      CALL mp_sum( io_group, intra_group_comm )
      CALL mp_sum( io_group, inter_group_comm )
      !
      ! ... find out the owner of each "z" plane
      !
      DO ip = 1, nproc_group
         !
         kowner((ipp(ip)+1):(ipp(ip)+npp(ip))) = ip - 1
         !
      END DO
      !
      ldr = nr1x*nr2x
      !
      ! ... explicit initialization to zero is needed because the physical
      ! ... dimensions pot may exceed the true size of the FFT grid 
      !
      pot(:) = 0.0_DP
      !
      DO k = 1, nr3
         !
         ! ... only ionode reads the potential planes
         !
         IF ( ionode ) &
            CALL iotk_scan_dat( potunit, "z" // iotk_index( k ), pot_plane )
         !
         ! ... planes are sent to the destination processor
         !
         IF( ngroup > 1 ) THEN
            !
            !  send to all proc/pools
            !
            IF( io_group_id == my_group_id ) THEN
               CALL mp_bcast( pot_plane, io_group, intra_group_comm )
            END IF
            CALL mp_bcast( pot_plane, io_group_id, inter_group_comm )
            !
         ELSE
            !
            !  send to the destination proc
            !
            IF ( kowner(k) /= io_group ) &
               CALL mp_put( pot_plane, pot_plane, me_group, io_group, kowner(k), k, intra_group_comm )
            !
         END IF
         !
         IF( kowner(k) == me_group ) THEN
            !
            kk = k - ipp( me_group + 1 )
            ! 
            DO j = 1, nr2
               !
               DO i = 1, nr1
                  !
                  pot(i+(j-1)*nr1x+(kk-1)*ldr) = pot_plane(i+(j-1)*nr1)
                  !
               END DO
               !
            END DO
            !
         END IF
         !
      END DO
      !
      DEALLOCATE( pot_plane )
      DEALLOCATE( kowner )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_end( potunit, "EFFECTIVE-POTENTIAL" )
         !
         CALL iotk_close_read( potunit )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE read_pot_xml   
    !
    ! ... methods to write and read charge_density
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_rho_xml( rho_file_base, rho, &
                              nr1, nr2, nr3, nr1x, nr2x, ipp, npp )
      !------------------------------------------------------------------------
      !
      ! ... Writes charge density rho, one plane at a time.
      ! ... If ipp and npp are specified, planes are collected one by one from
      ! ... all processors, avoiding an overall collect of the charge density
      ! ... on a single proc.
      !
      USE io_files,  ONLY : rhounit
      USE io_global, ONLY : ionode
      USE mp_global, ONLY : me_image, intra_image_comm, me_pool, nproc_pool, &
                            intra_pool_comm, my_pool_id
      USE mp,        ONLY : mp_get
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*),  INTENT(IN) :: rho_file_base
      INTEGER,           INTENT(IN) :: nr1, nr2, nr3
      INTEGER,           INTENT(IN) :: nr1x, nr2x
      REAL(DP),          INTENT(IN) :: rho(:)
      INTEGER, OPTIONAL, INTENT(IN) :: ipp(:)
      INTEGER, OPTIONAL, INTENT(IN) :: npp(:)
      !
      INTEGER               :: ierr, i, j, k, kk, ldr, ip
      CHARACTER(LEN=256)    :: rho_file
      CHARACTER(LEN=10)     :: rho_extension
      REAL(DP), ALLOCATABLE :: rho_plane(:)
      INTEGER,  ALLOCATABLE :: kowner(:)
      INTEGER               :: iopool_id, ionode_pool
      !
      !
      rho_extension = '.dat'
      IF ( .NOT. rho_binary ) rho_extension = '.xml'
      !
      rho_file = TRIM( rho_file_base ) // TRIM( rho_extension )
      !
      IF ( ionode ) &
         CALL iotk_open_write( rhounit, FILE = rho_file, &
                               BINARY = rho_binary, IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'write_rho_xml', 'cannot open' // &
                 & TRIM( rho_file ) // ' file for writing', ierr )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_begin( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_write_attr( attr, "nr1", nr1, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "nr2", nr2 )
         CALL iotk_write_attr( attr, "nr3", nr3 )
         !
         CALL iotk_write_empty( rhounit, "INFO", attr )
         !
      END IF
      !
      ALLOCATE( rho_plane( nr1*nr2 ) )
      ALLOCATE( kowner( nr3 ) )
      !
      ! ... find the index of the pool that will write rho
      !
      IF ( ionode ) iopool_id = my_pool_id
      !
      CALL mp_bcast( iopool_id, ionode_id, intra_image_comm )
      !
      ! ... find the index of the ionode within its own pool
      !
      IF ( ionode ) ionode_pool = me_pool
      !
      CALL mp_bcast( ionode_pool, ionode_id, intra_image_comm )
      !
      ! ... find out the owner of each "z" plane
      !
      IF ( PRESENT( ipp ) .AND. PRESENT( npp ) ) THEN
         !
         DO ip = 1, nproc_pool
            !
            kowner( (ipp(ip)+1):(ipp(ip)+npp(ip)) ) = ip - 1
            !
         END DO
         !
      ELSE
         !
         kowner = ionode_id
         !
      END IF
      !
      ldr = nr1x*nr2x
      !
      DO k = 1, nr3
         !
         IF( kowner(k) == me_pool ) THEN
            !
            kk = k
            !
            IF ( PRESENT( ipp ) ) kk = k - ipp(me_pool+1)
            ! 
            DO j = 1, nr2
               !
               DO i = 1, nr1
                  !
                  rho_plane(i+(j-1)*nr1) = rho(i+(j-1)*nr1x+(kk-1)*ldr)
                  !
               END DO
               !
            END DO
            !
         END IF
         !
         IF ( kowner(k) /= ionode_pool .AND. my_pool_id == iopool_id ) &
            CALL mp_get( rho_plane, rho_plane, &
                         me_pool, ionode_pool, kowner(k), k, intra_pool_comm )
         !
         IF ( ionode ) &
            CALL iotk_write_dat( rhounit, "z" // iotk_index( k ), rho_plane )
         !
      END DO
      !
      DEALLOCATE( rho_plane )
      DEALLOCATE( kowner )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_end( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_close_write( rhounit )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE write_rho_xml
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_rho_xml( rho_file_base, rho, &
                             nr1, nr2, nr3, nr1x, nr2x, ipp, npp )
      !------------------------------------------------------------------------
      !
      ! ... Reads charge density rho, one plane at a time.
      ! ... If ipp and npp are specified, planes are collected one by one from
      ! ... all processors, avoiding an overall collect of the charge density
      ! ... on a single proc.
      !
      USE io_files,  ONLY : rhounit
      USE io_global, ONLY : ionode, ionode_id
      USE mp_global, ONLY : me_image, intra_image_comm, me_pool, nproc_pool, &
                            intra_pool_comm, my_pool_id, npool
      USE mp,        ONLY : mp_put
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*),  INTENT(IN)  :: rho_file_base
      INTEGER,           INTENT(IN)  :: nr1, nr2, nr3
      INTEGER,           INTENT(IN)  :: nr1x, nr2x
      REAL(DP),          INTENT(OUT) :: rho(:)
      INTEGER, OPTIONAL, INTENT(IN)  :: ipp(:)
      INTEGER, OPTIONAL, INTENT(IN)  :: npp(:)
      !
      INTEGER               :: ierr, i, j, k, kk, ldr, ip
      INTEGER               :: nr( 3 )
      CHARACTER(LEN=256)    :: rho_file
      REAL(DP), ALLOCATABLE :: rho_plane(:)
      INTEGER,  ALLOCATABLE :: kowner(:)
      INTEGER               :: iopool_id, ionode_pool
      LOGICAL               :: exst
      !
      !
      rho_file = TRIM( rho_file_base ) // ".dat"
      exst = check_file_exst( TRIM(rho_file) ) 
      !
      IF ( .NOT. exst ) THEN
          !
          rho_file = TRIM( rho_file_base ) // ".xml"
          exst = check_file_exst( TRIM(rho_file) ) 
          !
      ENDIF
      !
      IF ( .NOT. exst ) CALL errore('read_rho_xml', 'searching for '//TRIM(rho_file), 10)
      !
      IF ( ionode ) &
         CALL iotk_open_read( rhounit, FILE = rho_file, IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'read_rho_xml', 'cannot open ' // &
                 & TRIM( rho_file ) // ' file for reading', ierr )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_scan_empty( rhounit, "INFO", attr )
         !
         CALL iotk_scan_attr( attr, "nr1", nr(1) )
         CALL iotk_scan_attr( attr, "nr2", nr(2) )
         CALL iotk_scan_attr( attr, "nr3", nr(3) )
         !
      END IF
      !
      CALL mp_bcast( nr, ionode_id, intra_image_comm )
      !
      IF ( nr1 /= nr(1) .OR. nr2 /= nr(2) .OR. nr3 /= nr(3) ) &
         CALL errore( 'read_rho_xml', 'dimensions do not match', 1 )
      !
      ALLOCATE( rho_plane( nr1*nr2 ) )
      ALLOCATE( kowner( nr3 ) )
      !
      ! ... find the index of the pool that will write rho
      !
      IF ( ionode ) iopool_id = my_pool_id
      !
      CALL mp_bcast( iopool_id, ionode_id, intra_image_comm )
      !
      ! ... find the index of the ionode within its own pool
      !
      IF ( ionode ) ionode_pool = me_pool
      !
      CALL mp_bcast( ionode_pool, ionode_id, intra_image_comm )
      !
      ! ... find out the owner of each "z" plane
      !
      IF ( PRESENT( ipp ) .AND. PRESENT( npp ) ) THEN
         !
         DO ip = 1, nproc_pool
            !
            kowner((ipp(ip)+1):(ipp(ip)+npp(ip))) = ip - 1
            !
         END DO
         !
      ELSE
         !
         kowner = ionode_id
         !
      END IF
      !
      ldr = nr1x*nr2x
      !
      ! ... explicit initialization to zero is needed because the physical
      ! ... dimensions rho may exceed the true size of the FFT grid 
      !
      rho(:) = 0.0_DP
      !
      DO k = 1, nr3
         !
         ! ... only ionode reads the charge planes
         !
         IF ( ionode ) &
            CALL iotk_scan_dat( rhounit, "z" // iotk_index( k ), rho_plane )
         !
         ! ... planes are sent to the destination processor
         !
         IF( npool > 1 ) THEN
            !
            !  send to all proc/pools
            !
            CALL mp_bcast( rho_plane, ionode_id, intra_image_comm )
            !
         ELSE
            !
            !  send to the destination proc
            !
            IF ( kowner(k) /= ionode_id ) &
               CALL mp_put( rho_plane, rho_plane, me_image, &
                            ionode_id, kowner(k), k, intra_image_comm )
            !
         END IF
         !
         IF( kowner(k) == me_pool ) THEN
            !
            kk = k
            !
            IF ( PRESENT( ipp ) ) kk = k - ipp(me_pool+1)
            ! 
            DO j = 1, nr2
               !
               DO i = 1, nr1
                  !
                  rho(i+(j-1)*nr1x+(kk-1)*ldr) = rho_plane(i+(j-1)*nr1)
                  !
               END DO
               !
            END DO
            !
         END IF
         !
      END DO
      !
      DEALLOCATE( rho_plane )
      DEALLOCATE( kowner )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_end( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_close_read( rhounit )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE read_rho_xml
    !
    ! ... methods to write and read wavefunctions
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_wfc_real( iuni, ik, nk, kunit, ispin, &
                          nspin, wf0, ngw, gamma_only, nbnd, igl, ngwl, filename, scalef )
      !------------------------------------------------------------------------
      !
      USE mp_wave,    ONLY : mergewf
      USE mp,         ONLY : mp_get
      USE mp_global,  ONLY : me_pool, nproc_image, nproc_pool, &
                             root_pool, intra_pool_comm, me_image, &
                             intra_image_comm
      !
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN) :: iuni
      INTEGER,            INTENT(IN) :: ik, nk, kunit, ispin, nspin
      COMPLEX(DP),        INTENT(IN) :: wf0(:,:)
      INTEGER,            INTENT(IN) :: ngw
      LOGICAL,            INTENT(IN) :: gamma_only
      INTEGER,            INTENT(IN) :: nbnd
      INTEGER,            INTENT(IN) :: ngwl
      INTEGER,            INTENT(IN) :: igl(:)
      CHARACTER(LEN=256), INTENT(IN) :: filename
      REAL(DP),           INTENT(IN) :: scalef    
        ! scale factor, usually 1.0 for pw and 1/SQRT( omega ) for CP
      !
      INTEGER                  :: j
      INTEGER                  :: iks, ike, ikt, igwx
      INTEGER                  :: npool, ipmask(nproc_image), ipsour
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      !
      !
      CALL set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                             npool, ikt, iks, ike, igwx, ipmask, ipsour )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_open_write( iuni, FILE = TRIM( filename ), ROOT="WFC", BINARY = .TRUE. )
         !
         CALL iotk_write_attr( attr, "ngw",          ngw, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "igwx",         igwx )
         CALL iotk_write_attr( attr, "gamma_only",   gamma_only )
         CALL iotk_write_attr( attr, "nbnd",         nbnd )
         CALL iotk_write_attr( attr, "ik",           ik )
         CALL iotk_write_attr( attr, "nk",           nk )
         CALL iotk_write_attr( attr, "ispin",        ispin )
         CALL iotk_write_attr( attr, "nspin",        nspin )
         CALL iotk_write_attr( attr, "scale_factor", scalef )
         !
         CALL iotk_write_empty( iuni, "INFO", attr )
         !
      END IF
      !
      ALLOCATE( wtmp( MAX( igwx, 1 ) ) )
      !
      wtmp = 0.0_DP
      !
      DO j = 1, nbnd
         !
         IF ( npool > 1 ) THEN
            !
            IF ( ikt >= iks .AND. ikt <= ike ) &      
               CALL mergewf( wf0(:,j), wtmp, ngwl, igl, me_pool, &
                             nproc_pool, root_pool, intra_pool_comm )
            !
            IF ( ipsour /= ionode_id ) &
               CALL mp_get( wtmp, wtmp, me_image, &
                            ionode_id, ipsour, j, intra_image_comm )
            !
         ELSE
            !
            CALL mergewf( wf0(:,j), wtmp, ngwl, igl, &
                          me_image, nproc_image, ionode_id, intra_image_comm )
            !
         END IF
         !
         IF ( ionode ) &
            CALL iotk_write_dat( iuni, "evc" // iotk_index( j ), wtmp(1:igwx) )
         !
      END DO
      !
      IF ( ionode ) CALL iotk_close_write( iuni )
      !
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE write_wfc_real
    !
!------------------------------------------------------------------------
    SUBROUTINE write_wfc_cmplx( iuni, ik, nk, kunit, ispin, &
                          nspin, wf0, ngw, do_wf_cmplx, gamma_only, nbnd, igl, ngwl, filename, scalef ) !added:giovanni do_wf_cmplx
      !------------------------------------------------------------------------
      !
      USE mp_wave,    ONLY : mergewf
      USE mp,         ONLY : mp_get
      USE mp_global,  ONLY : me_pool, nproc_image, nproc_pool, &
                             root_pool, intra_pool_comm, me_image, &
                             intra_image_comm
      !
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN) :: iuni
      INTEGER,            INTENT(IN) :: ik, nk, kunit, ispin, nspin
      COMPLEX(DP),        INTENT(IN) :: wf0(:,:)
      INTEGER,            INTENT(IN) :: ngw
      LOGICAL,            INTENT(IN) :: do_wf_cmplx ! added:giovanni
      LOGICAL,            INTENT(IN) :: gamma_only
      INTEGER,            INTENT(IN) :: nbnd
      INTEGER,            INTENT(IN) :: ngwl
      INTEGER,            INTENT(IN) :: igl(:)
      CHARACTER(LEN=256), INTENT(IN) :: filename
      REAL(DP),           INTENT(IN) :: scalef    
        ! scale factor, usually 1.0 for pw and 1/SQRT( omega ) for CP
      !
      INTEGER                  :: j
      INTEGER                  :: iks, ike, ikt, igwx
      INTEGER                  :: npool, ipmask(nproc_image), ipsour
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      !
      !
      CALL set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                             npool, ikt, iks, ike, igwx, ipmask, ipsour )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_open_write( iuni, FILE = TRIM( filename ), ROOT="WFC", BINARY = .TRUE. )
         !
         CALL iotk_write_attr( attr, "ngw",          ngw, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "igwx",         igwx )
         CALL iotk_write_attr( attr, "do_wf_cmplx",   do_wf_cmplx ) !added:giovanni
         CALL iotk_write_attr( attr, "gamma_only",   gamma_only.and..not.do_wf_cmplx )
         CALL iotk_write_attr( attr, "nbnd",         nbnd )
         CALL iotk_write_attr( attr, "ik",           ik )
         CALL iotk_write_attr( attr, "nk",           nk )
         CALL iotk_write_attr( attr, "ispin",        ispin )
         CALL iotk_write_attr( attr, "nspin",        nspin )
         CALL iotk_write_attr( attr, "scale_factor", scalef )
         !
         CALL iotk_write_empty( iuni, "INFO", attr )
         !
      END IF
      !
      ALLOCATE( wtmp( MAX( igwx, 1 ) ) )
      !
      wtmp = 0.0_DP
      !
      DO j = 1, nbnd
         !
         IF ( npool > 1 ) THEN
            !
            IF ( ikt >= iks .AND. ikt <= ike ) &      
               CALL mergewf( wf0(:,j), wtmp, ngwl, igl, me_pool, &
                             nproc_pool, root_pool, intra_pool_comm )
            !
            IF ( ipsour /= ionode_id ) &
               CALL mp_get( wtmp, wtmp, me_image, &
                            ionode_id, ipsour, j, intra_image_comm )
            !
         ELSE
            !
            CALL mergewf( wf0(:,j), wtmp, ngwl, igl, &
                          me_image, nproc_image, ionode_id, intra_image_comm )
            !
         END IF
         !
         IF ( ionode ) &
            CALL iotk_write_dat( iuni, "evc" // iotk_index( j ), wtmp(1:igwx) )
         !
      END DO
      !
      IF ( ionode ) CALL iotk_close_write( iuni )
      !
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE write_wfc_cmplx

    !------------------------------------------------------------------------
    SUBROUTINE read_wfc( iuni, ik, nk, kunit, ispin, &
                         nspin, wf, ngw, nbnd, igl, ngwl, filename, scalef, &
                         flink )
      !------------------------------------------------------------------------
      !
      USE mp_wave,   ONLY : splitwf
      USE mp,        ONLY : mp_put
      USE mp_global, ONLY : me_image, nproc_image, root_image, me_pool, my_pool_id, &
                            nproc_pool, intra_pool_comm, root_pool, my_image_id, &
                            intra_image_comm
      !
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN)    :: iuni
      COMPLEX(DP),        INTENT(OUT)   :: wf(:,:)
      INTEGER,            INTENT(IN)    :: ik, nk
      INTEGER,            INTENT(IN)    :: kunit
      INTEGER,            INTENT(INOUT) :: ngw, nbnd, ispin, nspin
      INTEGER,            INTENT(IN)    :: ngwl
      INTEGER,            INTENT(IN)    :: igl(:)
      CHARACTER(LEN=256), INTENT(IN)    :: filename
      REAL(DP),           INTENT(OUT)   :: scalef
      LOGICAL, OPTIONAL,  INTENT(IN)    :: flink
      !
      INTEGER                  :: j
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      INTEGER                  :: ierr
      INTEGER                  :: iks, ike, ikt
      INTEGER                  :: igwx, igwx_, ik_, nk_
      INTEGER                  :: npool, ipmask(nproc_image), ipdest
      LOGICAL                  :: flink_
      !
      flink_ = .FALSE.
      IF( PRESENT( flink ) ) flink_ = flink
      !
      CALL set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                             npool, ikt, iks, ike, igwx, ipmask, ipdest )
      !
      !  if flink = .true. we are following a link and the file is
      !  already opened for read
      !
      ierr = 0
      !
      IF ( ionode .AND. .NOT. flink_ ) &
         CALL iotk_open_read( iuni, FILE = filename, &
                              BINARY = .TRUE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'read_wfc ', &
                   'cannot open restart file for reading', ierr )
      !
      IF ( ionode ) THEN
          !
          CALL iotk_scan_empty( iuni, "INFO", attr )
          !
          CALL iotk_scan_attr( attr, "ngw",          ngw )
          CALL iotk_scan_attr( attr, "nbnd",         nbnd )
          CALL iotk_scan_attr( attr, "ik",           ik_ )
          CALL iotk_scan_attr( attr, "nk",           nk_ )
          CALL iotk_scan_attr( attr, "ispin",        ispin )
          CALL iotk_scan_attr( attr, "nspin",        nspin )
          CALL iotk_scan_attr( attr, "igwx",         igwx_ )
          CALL iotk_scan_attr( attr, "scale_factor", scalef )
          !
      END IF
      !
      CALL mp_bcast( ngw,    ionode_id, intra_image_comm )
      CALL mp_bcast( nbnd,   ionode_id, intra_image_comm )
      CALL mp_bcast( ik_,    ionode_id, intra_image_comm )
      CALL mp_bcast( nk_,    ionode_id, intra_image_comm )
      CALL mp_bcast( ispin,  ionode_id, intra_image_comm )
      CALL mp_bcast( nspin,  ionode_id, intra_image_comm )
      CALL mp_bcast( igwx_,  ionode_id, intra_image_comm )
      CALL mp_bcast( scalef, ionode_id, intra_image_comm )
      !
      ALLOCATE( wtmp( MAX( igwx_, igwx ) ) )
      !
      DO j = 1, nbnd
         !
         IF ( j <= SIZE( wf, 2 ) ) THEN
            !
            IF ( ionode ) THEN 
               !
               CALL iotk_scan_dat( iuni, &
                                   "evc" // iotk_index( j ), wtmp(1:igwx_) )
               !
               IF ( igwx > igwx_ ) wtmp((igwx_+1):igwx) = 0.0_DP
               !
            END IF
            !
            IF ( npool > 1 ) THEN
               !
               IF ( ipdest /= ionode_id ) &
                  CALL mp_put( wtmp, wtmp, me_image, &
                               ionode_id, ipdest, j, intra_image_comm )
               !
               IF ( ( ikt >= iks ) .AND. ( ikt <= ike ) ) &
                  CALL splitwf( wf(:,j), wtmp, ngwl, igl, me_pool, &
                                nproc_pool, root_pool, intra_pool_comm )
               !
            ELSE
               !
               CALL splitwf( wf(:,j), wtmp, ngwl, igl, &
                             me_image, nproc_image, ionode_id, intra_image_comm )
               !
            END IF
            !
         END IF
         !
      END DO
      !
      IF ( ionode .AND. .NOT. flink_ ) CALL iotk_close_read( iuni )
      !
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE read_wfc
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_eig( iuni, filename, nbnd, eig, energy_units, &
                          occ, ik, ispin, lkpoint_dir )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN) :: iuni
      INTEGER,            INTENT(IN) :: nbnd
      REAL(DP),           INTENT(IN) :: eig(:)
      CHARACTER(*),       INTENT(IN) :: energy_units
      REAL(DP), OPTIONAL, INTENT(IN) :: occ(:)
      INTEGER,  OPTIONAL, INTENT(IN) :: ik, ispin
      LOGICAL,  OPTIONAL, INTENT(IN) :: lkpoint_dir
      CHARACTER(LEN=256), INTENT(IN) :: filename
      LOGICAL :: lkpoint_dir0
      !
      lkpoint_dir0=.TRUE.
      IF (present(lkpoint_dir)) lkpoint_dir0=lkpoint_dir
      IF ( ionode ) THEN
         !
         if (lkpoint_dir0) CALL iotk_open_write ( iuni, &
                           FILE = TRIM( filename ), BINARY = .FALSE. )
         !
         CALL iotk_write_attr ( attr, "nbnd", nbnd, FIRST=.TRUE. )
         IF ( PRESENT( ik) )    CALL iotk_write_attr ( attr, "ik", ik )
         IF ( PRESENT( ispin) ) CALL iotk_write_attr ( attr, "ispin", ispin )
         CALL iotk_write_empty( iuni, "INFO", ATTR = attr )
         !
         CALL iotk_write_attr ( attr, "UNITS", TRIM(energy_units), FIRST = .TRUE. )
         CALL iotk_write_empty( iuni, "UNITS_FOR_ENERGIES", ATTR=attr)
         !
         CALL iotk_write_dat( iuni, "EIGENVALUES", eig(:) )
         !
         IF ( PRESENT( occ ) ) THEN
            !
            CALL iotk_write_dat( iuni, "OCCUPATIONS", occ(:) )
            !
         ENDIF
         !
         IF (lkpoint_dir0) CALL iotk_close_write ( iuni )
         !
      ENDIF
      !
    END SUBROUTINE write_eig
    ! 
END MODULE xml_io_base
