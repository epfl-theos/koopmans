! Copyright (C) 2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
MODULE cp_restart
  !-----------------------------------------------------------------------------
  ! 
  ! ... This module contains subroutines to write and read data required to
  ! ... restart a calculation from the disk  
  !
  USE iotk_module
  USE xml_io_base, ONLY : default_fmt_version => fmt_version
  USE xml_io_base
  !
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : ionode, ionode_id, stdout
  USE io_files,    ONLY : prefix, iunpun, xmlpun, qexml_version, qexml_version_init
  USE mp,          ONLY : mp_bcast
  USE parser,      ONLY : version_compare
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE :: read_cell
  !
  INTEGER, PRIVATE :: iunout
  !
  !
  ! variables to describe qexml current version
  ! and back compatibility
  !
  LOGICAL, PRIVATE :: qexml_version_before_1_4_0 = .FALSE.
  !
  INTERFACE cp_writefile
     module procedure cp_writefile_twin, cp_writefile_real
  END INTERFACE cp_writefile

  INTERFACE cp_readfile
     module procedure cp_readfile_twin, cp_readfile_real
  END INTERFACE cp_readfile
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE cp_writefile_twin( ndw, outdir, ascii, nfi, simtime, acc, nk, xk, &
                             wk, ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh,   &
                             taui, cdmi, stau0, svel0, staum, svelm, force,  &
                             vnhp, xnhp0, xnhpm, nhpcl, nhpdim, occ0, occm,  &
                             lambda0,lambdam,lambda_bare, xnhe0, xnhem, vnhe,&
                             ekincm, et, rho, c02, cm2, ctot, iupdwn, nupdwn,&
                             iupdwn_tot, nupdwn_tot, mat_z )
      !------------------------------------------------------------------------
      !
      USE control_flags,            ONLY : do_wf_cmplx, gamma_only, force_pairing, trhow, tksw !added:giovanni do_wf_cmplx!
      USE control_flags,            ONLY : evc_restart !added:giovanni evc_restart
      USE io_files,                 ONLY : psfile, pseudo_dir
      USE mp_global,                ONLY : intra_image_comm, me_image, nproc_image
      USE printout_base,            ONLY : title
      USE grid_dimensions,          ONLY : nr1, nr2, nr3, nr1x, nr2x, nr3l
      USE smooth_grid_dimensions,   ONLY : nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, ONLY : nr1b, nr2b, nr3b
      USE gvecp,                    ONLY : ngm, ngmt, ecutp, gcutp
      USE gvecs,                    ONLY : ngs, ngst, ecuts, gcuts, dual
      USE gvecw,                    ONLY : ngw, ngwt, ecutw, gcutw
      USE reciprocal_vectors,       ONLY : ig_l2g, mill_l
      USE electrons_base,           ONLY : nspin, nelt, nel, nudx
      USE cell_base,                ONLY : ibrav, alat, celldm, &
                                           symm_type, s_to_r
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, &
                                           pmass, amass, iforce, ind_bck
      USE funct,                    ONLY : get_dft_name
      USE energies,                 ONLY : enthal, ekin, eht, esr, eself, &
                                           epseu, enl, exc, vave
      USE mp_global,                ONLY : nproc, mpime
      USE mp,                       ONLY : mp_sum
      USE fft_base,                 ONLY : dfftp
      USE constants,                ONLY : pi
      USE cp_interfaces,            ONLY : n_atom_wfc
      USE global_version,           ONLY : version_number
      USE cp_main_variables,        ONLY : collect_lambda, descla, collect_zmat
      USE twin_types !added:giovanni
      USE electrons_base,       ONLY : nbsp !added:giovanni
      USE nksic,                ONLY : do_bare_eigs !added:giovanni
      USE input_parameters,     ONLY : odd_nkscalfact, restart_odd_nkscalfact, print_wfc_empty
      USE wavefunctions_module, ONLY : c0_fixed, c0_emp_aux, cm_emp_aux
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN) :: ndw          !
      CHARACTER(LEN=*),      INTENT(IN) :: outdir       !  directory used to store output and restart files
      LOGICAL,               INTENT(IN) :: ascii        !
      INTEGER,               INTENT(IN) :: nfi          ! index of the current step
      REAL(DP),              INTENT(IN) :: simtime      ! simulated time
      REAL(DP),              INTENT(IN) :: acc(:)       !  
      INTEGER,               INTENT(IN) :: nk           ! number of kpoints
      REAL(DP),              INTENT(IN) :: xk(:,:)      ! k-points coordinates 
      REAL(DP),              INTENT(IN) :: wk(:)        ! k-points weights
      REAL(DP),              INTENT(IN) :: ht(3,3)      ! 
      REAL(DP),              INTENT(IN) :: htm(3,3)     ! 
      REAL(DP),              INTENT(IN) :: htvel(3,3)   ! 
      REAL(DP),              INTENT(IN) :: gvel(3,3)    ! 
      REAL(DP),              INTENT(IN) :: xnhh0(3,3)   ! 
      REAL(DP),              INTENT(IN) :: xnhhm(3,3)   ! 
      REAL(DP),              INTENT(IN) :: vnhh(3,3)    ! 
      REAL(DP),              INTENT(IN) :: taui(:,:)    ! 
      REAL(DP),              INTENT(IN) :: cdmi(:)      ! 
      REAL(DP),              INTENT(IN) :: stau0(:,:)   ! 
      REAL(DP),              INTENT(IN) :: svel0(:,:)   ! 
      REAL(DP),              INTENT(IN) :: staum(:,:)   ! 
      REAL(DP),              INTENT(IN) :: svelm(:,:)   ! 
      REAL(DP),              INTENT(IN) :: force(:,:)   ! 
      REAL(DP),              INTENT(IN) :: xnhp0(:)     ! 
      REAL(DP),              INTENT(IN) :: xnhpm(:)     ! 
      REAL(DP),              INTENT(IN) :: vnhp(:)      ! 
      INTEGER,               INTENT(IN) :: nhpcl        ! 
      INTEGER,               INTENT(IN) :: nhpdim       ! 
      REAL(DP),              INTENT(IN) :: occ0(:)      !  occupations of electronic states
      REAL(DP),              INTENT(IN) :: occm(:)      ! 
      TYPE(twin_matrix), DIMENSION(:), INTENT(IN) :: lambda0
      TYPE(twin_matrix), DIMENSION(:), INTENT(IN) :: lambdam
      TYPE(twin_matrix), DIMENSION(:), INTENT(IN) :: lambda_bare
!       REAL(DP),              INTENT(IN) :: lambda0(:,:,:) ! !removed:giovanni
!       REAL(DP),              INTENT(IN) :: lambdam(:,:,:) ! !removed:giovanni
      REAL(DP),              INTENT(IN) :: xnhe0        ! 
      REAL(DP),              INTENT(IN) :: xnhem        ! 
      REAL(DP),              INTENT(IN) :: vnhe         ! 
      REAL(DP),              INTENT(IN) :: ekincm       ! 
      REAL(DP),              INTENT(IN) :: et(:,:)      !  eigenvalues
      REAL(DP),              INTENT(IN) :: rho(:,:)     ! 
      COMPLEX(DP),           INTENT(IN) :: c02(:,:)     ! 
      COMPLEX(DP),           INTENT(IN) :: cm2(:,:)     ! 
      COMPLEX(DP),           INTENT(IN) :: ctot(:,:)    ! 
      INTEGER,               INTENT(IN) :: iupdwn(:)    ! 
      INTEGER,               INTENT(IN) :: nupdwn(:)    ! 
      INTEGER,               INTENT(IN) :: iupdwn_tot(:)! 
      INTEGER,               INTENT(IN) :: nupdwn_tot(:)! 
!       REAL(DP),    OPTIONAL, INTENT(IN) :: mat_z(:,:,:) ! 
      TYPE(twin_matrix), DIMENSION(:), OPTIONAL, INTENT(IN) :: mat_z
      !
      LOGICAL               :: write_charge_density
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=256)    :: dirname, filename, rho_file_base
      CHARACTER(LEN=4)      :: cspin
      INTEGER               :: kunit, ib, ik_eff
      INTEGER               :: k1, k2, k3
      INTEGER               :: nk1, nk2, nk3
      INTEGER               :: j, i, iss, ig, nspin_wfc, iss_wfc
      INTEGER               :: is, ia, isa, ik, ierr
      INTEGER,  ALLOCATABLE :: mill(:,:)
      INTEGER,  ALLOCATABLE :: ftmp(:,:)
      INTEGER,  ALLOCATABLE :: ityp(:)
      REAL(DP), ALLOCATABLE :: tau(:,:)
      REAL(DP), ALLOCATABLE :: dtmp(:)
      REAL(DP), ALLOCATABLE :: rhoaux(:)
      REAL(DP)              :: omega, htm1(3,3), h(3,3)
      REAL(DP)              :: a1(3), a2(3), a3(3)
      REAL(DP)              :: b1(3), b2(3), b3(3)
      REAL(DP)              :: nelec
      REAL(DP)              :: scalef
      LOGICAL               :: lsda
      REAL(DP)              :: s0, s1, cclock
      INTEGER               :: nbnd_tot
      INTEGER               :: nbnd_emp
      INTEGER               :: nbnd_
      REAL(DP), ALLOCATABLE :: mrepl(:,:)
      COMPLEX(DP), ALLOCATABLE :: mrepl_c(:,:) !added:giovanni
      COMPLEX(DP) :: csc(nbsp, nbsp)
      INTEGER:: icsc, icsc2
      !
      write_charge_density = trhow
      !
      IF( nspin > 1 .AND. .NOT. force_pairing ) THEN
         !
         !  check if the array storing wave functions is large enought
         !
         IF( SIZE( c02, 2 ) < ( iupdwn( 2 ) + nupdwn(1) - 1 ) ) &
            CALL errore('cp_writefile',' wrong wave functions dimension ', 1 )
         !
      END IF
      !
      IF(  nupdwn_tot(1) < nupdwn(1) ) &
         CALL errore( " writefile ", " wrong number of states ", 1 )
      !
      nbnd_    = nupdwn(1) 
      nbnd_tot = MAX( nupdwn(1), nupdwn_tot(1) )
      nbnd_emp = MAX( 0, nupdwn_tot(1) - nupdwn(1) )
      !
      IF ( ionode ) THEN
         !
         ! ... look for an empty unit (only ionode needs it)
         !
         CALL iotk_free_unit( iunout, ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_writefile', &
                   'no free units to write wavefunctions', ierr )
      !
      dirname = restart_dir( outdir, ndw )
      !
      ! ... Create main restart directory
      !
      CALL create_directory( dirname )
      !
      ! ... Create k-points subdirectories
      ! ... note: in FPMD and CP k-points are not distributed to processors
      !
      DO i = 1, nk
         !
         CALL create_directory( kpoint_dir( dirname, i ) )
         !
      END DO
      !
      ! ... Some ( CP/FPMD ) default values
      !
      IF ( nspin == 2 ) THEN
         !
         kunit = 2
         !
      ELSE
         !
         kunit = 1
         !
      END IF
      !
      k1  = 0
      k2  = 0
      k3  = 0
      nk1 = 0
      nk2 = 0
      nk3 = 0
      !
      ! ... Compute Cell related variables
      !
      h = TRANSPOSE( ht )
      !
      CALL invmat( 3, ht, htm1, omega )
      !
      a1 = ht(1,:)
      a2 = ht(2,:)
      a3 = ht(3,:)
      !
      ! ... Beware: omega may be negative if axis are left-handed!
      !
      scalef = 1.D0 / SQRT( ABS (omega) )
      !
      ! ... Compute array ityp, and tau
      !
      ALLOCATE( ityp( nat ) )
      ALLOCATE( tau( 3, nat ) )
      !
      isa = 0
      !
      DO is = 1, nsp
         !
         DO ia = 1, na(is)
            !
            isa = isa + 1
            ityp(isa) = is
            !
         END DO
         !
      END DO
      !
      CALL s_to_r( stau0, tau, na, nsp, h )
      !   
      ! ... Collect G vectors
      !   
      ALLOCATE( mill( 3, ngmt ) )
      !
      mill = 0
      !
      mill(:,ig_l2g(1:ngm)) = mill_l(:,1:ngm)
      !
      CALL mp_sum( mill, intra_image_comm )
      !
      lsda = ( nspin == 2 )
      !
      ALLOCATE( ftmp( nbnd_tot , nspin ) )
      !
      ftmp = 0.0d0
      !
      DO iss = 1, nspin
         !
         ftmp( 1:nupdwn(iss), iss ) = occ0( iupdwn(iss) : iupdwn(iss) + nupdwn(iss) - 1 )
         !
      END DO
      !
      IF ( ionode ) THEN
         !
         ! ... Open XML descriptor
         !
         WRITE( stdout, '(/,3X,"writing restart file: ",A)' ) TRIM( dirname )
         !
         CALL iotk_open_write( iunpun, FILE = TRIM( dirname ) // '/' // &
                             & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_writefile ', 'cannot open restart file for writing', ierr )
      !
      s0 = cclock() 
      !
      IF ( ionode ) THEN

!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
         CALL write_header( "CP", TRIM(version_number) )
         !
!-------------------------------------------------------------------------------
! ... this flag is used to check if the file can be used for post-processing
!-------------------------------------------------------------------------------
         !
         CALL write_control( PP_CHECK_FLAG=.TRUE. )
         !
!-------------------------------------------------------------------------------
! ... STATUS
!-------------------------------------------------------------------------------
         !
         CALL iotk_write_begin( iunpun, "STATUS" )
         !
         CALL iotk_write_attr( attr, "ITERATION", nfi, FIRST = .TRUE. )
         CALL iotk_write_empty(iunpun, "STEP", attr )
         !
         CALL iotk_write_attr( attr, "UNITS", "pico-seconds", FIRST = .TRUE. ) 
         CALL iotk_write_dat( iunpun, "TIME", simtime, ATTR = attr )
         !
         CALL iotk_write_dat( iunpun, "TITLE", TRIM( title ) )
         !
         CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
         CALL iotk_write_dat( iunpun, "KINETIC_ENERGY", ekin,   ATTR = attr )
         CALL iotk_write_dat( iunpun, "HARTREE_ENERGY", eht,    ATTR = attr )
         CALL iotk_write_dat( iunpun, "EWALD_TERM",     esr,    ATTR = attr )
         CALL iotk_write_dat( iunpun, "GAUSS_SELFINT",  eself,  ATTR = attr )
         CALL iotk_write_dat( iunpun, "LPSP_ENERGY",    epseu,  ATTR = attr )
         CALL iotk_write_dat( iunpun, "NLPSP_ENERGY",   enl,    ATTR = attr )
         CALL iotk_write_dat( iunpun, "EXC_ENERGY",     exc,    ATTR = attr )
         CALL iotk_write_dat( iunpun, "AVERAGE_POT",    vave,   ATTR = attr )
         CALL iotk_write_dat( iunpun, "ENTHALPY",       enthal, ATTR = attr )
         !
         CALL iotk_write_end( iunpun, "STATUS" )      
         !
!-------------------------------------------------------------------------------
! ... CELL
!-------------------------------------------------------------------------------
         !
         a1 = a1 / alat
         a2 = a2 / alat
         a3 = a3 / alat
         !
         CALL recips( a1, a2, a3, b1, b2, b3 )
         !
         CALL write_cell( ibrav, symm_type, &
                          celldm, alat, a1, a2, a3, b1, b2, b3 )
         !
!-------------------------------------------------------------------------------
! ... IONS
!-------------------------------------------------------------------------------
         !
         CALL write_ions( nsp, nat, atm, ityp(ind_bck(:)), &
                          psfile, pseudo_dir, amass, tau(:,ind_bck(:)), &
                          iforce(:,ind_bck(:)), dirname, 1.D0 )
         !
!-------------------------------------------------------------------------------
! ... PLANE_WAVES
!-------------------------------------------------------------------------------
         !
         ! change to .TRUE. to write gvectors.dat for rho
         !
         CALL write_planewaves( ecutw, dual, ngwt, do_wf_cmplx, gamma_only, & !added:giovanni do_wf_cmplx
                                nr1, nr2, nr3, ngmt, nr1s, nr2s, nr3s,     &
                                ngst, nr1b, nr2b, nr3b, mill, .TRUE. )
         !
!-------------------------------------------------------------------------------
! ... SPIN
!-------------------------------------------------------------------------------
         !
         CALL write_spin( lsda, .FALSE., 1, .FALSE., .TRUE. )
         !
!-------------------------------------------------------------------------------
! ... EXCHANGE_CORRELATION
!-------------------------------------------------------------------------------
         !
         dft_name = get_dft_name()
         CALL write_xc( DFT = dft_name, NSP = nsp, LDA_PLUS_U = .FALSE. )
         !
!-------------------------------------------------------------------------------
! ... OCCUPATIONS
!-------------------------------------------------------------------------------
         !
         CALL write_occ( LGAUSS = .FALSE., LTETRA = .FALSE., &
                         TFIXED_OCC = .TRUE., LSDA = lsda, NSTATES_UP = nupdwn_tot(1), &
                         NSTATES_DOWN = nupdwn_tot(2), F_INP = DBLE( ftmp ) )
         !
!-------------------------------------------------------------------------------
! ... BRILLOUIN_ZONE
!-------------------------------------------------------------------------------
         !
         CALL write_bz( nk, xk, wk, k1, k2, k3, nk1, nk2, nk3, 0.0_DP )
         !
!-------------------------------------------------------------------------------
! ... PARALLELISM
!-------------------------------------------------------------------------------
         !
         CALL iotk_write_begin( iunpun, "PARALLELISM" )
         !
         CALL iotk_write_dat( iunpun, &
                              "GRANULARITY_OF_K-POINTS_DISTRIBUTION", kunit )
         !
         CALL iotk_write_end( iunpun, "PARALLELISM" )
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... CHARGE-DENSITY
!-------------------------------------------------------------------------------
      !
      IF (write_charge_density) then
         !
         rho_file_base = 'charge-density'
         !
         IF ( ionode )&
              CALL iotk_link( iunpun, "CHARGE-DENSITY", rho_file_base, &
              CREATE = .FALSE., BINARY = .TRUE. )
         !
         rho_file_base = TRIM( dirname ) // '/' // TRIM( rho_file_base )
         !
         IF ( nspin == 1 ) THEN
            !
            CALL write_rho_xml( rho_file_base, rho(:,1), &
                                nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
            !
         ELSE IF ( nspin == 2 ) THEN
            !
            ALLOCATE( rhoaux( SIZE( rho, 1 ) ) )
            !
            rhoaux = rho(:,1) + rho(:,2) 
            !
            CALL write_rho_xml( rho_file_base, rhoaux, &
                                nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
            !
            rho_file_base = 'spin-polarization'
            !
            IF ( ionode ) &
                 CALL iotk_link( iunpun, "SPIN-POLARIZATION", rho_file_base, &
                 CREATE = .FALSE., BINARY = .TRUE. )
            !
            rho_file_base = TRIM( dirname ) // '/' // TRIM( rho_file_base )
            !
            rhoaux = rho(:,1) - rho(:,2) 
            !
            CALL write_rho_xml( rho_file_base, rhoaux, &
                                nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
            !
            DEALLOCATE( rhoaux )
            !
         END IF
         !
      END IF ! write_charge_density
      !
!-------------------------------------------------------------------------------
! ... TIMESTEPS
!-------------------------------------------------------------------------------
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_attr( attr, "nt", 2, FIRST = .TRUE. )
         !
         CALL iotk_write_begin( iunpun, "TIMESTEPS", attr )
         !
         ! ... STEP0
         !
         CALL iotk_write_begin( iunpun, "STEP0" )
         !
         CALL iotk_write_dat( iunpun, "ACCUMULATORS", acc )
         !
         CALL iotk_write_begin( iunpun, "IONS_POSITIONS" )
         CALL iotk_write_dat(   iunpun, "stau",  stau0(1:3,1:nat),   COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "svel",  svel0(1:3,1:nat),   COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "taui",  taui(1:3,1:nat),    COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "cdmi",  cdmi(1:3),          COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "force", force(1:3,1:nat),   COLUMNS=3 )
         CALL iotk_write_end(   iunpun, "IONS_POSITIONS" )
         !
         CALL iotk_write_begin( iunpun, "IONS_NOSE" )
         CALL iotk_write_dat(   iunpun, "nhpcl", nhpcl )
         CALL iotk_write_dat(   iunpun, "nhpdim", nhpdim )
         CALL iotk_write_dat(   iunpun, "xnhp",  xnhp0(1:nhpcl*nhpdim) )
         CALL iotk_write_dat(   iunpun, "vnhp",  vnhp(1:nhpcl*nhpdim) )
         CALL iotk_write_end(   iunpun, "IONS_NOSE" )
         !
         CALL iotk_write_dat( iunpun, "ekincm", ekincm )
         !
         CALL iotk_write_begin( iunpun, "ELECTRONS_NOSE" )
         CALL iotk_write_dat(   iunpun, "xnhe", xnhe0 )
         CALL iotk_write_dat(   iunpun, "vnhe", vnhe )
         CALL iotk_write_end(   iunpun, "ELECTRONS_NOSE" )
         !
         CALL iotk_write_begin( iunpun, "CELL_PARAMETERS" )
         CALL iotk_write_dat(   iunpun, "ht",    ht )
         CALL iotk_write_dat(   iunpun, "htvel", htvel )
         CALL iotk_write_dat(   iunpun, "gvel",  gvel )
         CALL iotk_write_end(   iunpun, "CELL_PARAMETERS" )
         !
         CALL iotk_write_begin( iunpun, "CELL_NOSE" )
         CALL iotk_write_dat(   iunpun, "xnhh", xnhh0 )
         CALL iotk_write_dat(   iunpun, "vnhh", vnhh )
         CALL iotk_write_end(   iunpun, "CELL_NOSE" )
         !
         CALL iotk_write_end( iunpun, "STEP0" )
         !
         ! ... STEPM
         !
         CALL iotk_write_begin( iunpun, "STEPM" )
         !
         CALL iotk_write_begin( iunpun, "IONS_POSITIONS" )
         CALL iotk_write_dat(   iunpun, "stau", staum(1:3,1:nat),  COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "svel", svelm(1:3,1:nat),  COLUMNS=3 )
         CALL iotk_write_end(   iunpun, "IONS_POSITIONS" )
         !
         CALL iotk_write_begin( iunpun, "IONS_NOSE" )
         CALL iotk_write_dat(   iunpun, "nhpcl", nhpcl )
         CALL iotk_write_dat(   iunpun, "nhpdim", nhpdim )
         CALL iotk_write_dat(   iunpun, "xnhp",  xnhpm(1:nhpcl*nhpdim) )
         CALL iotk_write_end(   iunpun, "IONS_NOSE" )
         !
         CALL iotk_write_begin( iunpun, "ELECTRONS_NOSE" )
         CALL iotk_write_dat(   iunpun, "xnhe", xnhem )
         CALL iotk_write_end(   iunpun, "ELECTRONS_NOSE" )
         !
         CALL iotk_write_begin( iunpun, "CELL_PARAMETERS" )
         CALL iotk_write_dat(   iunpun, "ht",    htm )
         CALL iotk_write_end(   iunpun, "CELL_PARAMETERS" )
         !
         CALL iotk_write_begin( iunpun, "CELL_NOSE" )
         CALL iotk_write_dat(   iunpun, "xnhh", xnhhm )
         CALL iotk_write_end(   iunpun, "CELL_NOSE" )
         !
         CALL iotk_write_end( iunpun, "STEPM" )
         !
         CALL iotk_write_end( iunpun, "TIMESTEPS" )
         !
      END IF

!-------------------------------------------------------------------------------
! ... BAND_STRUCTURE_INFO
!-------------------------------------------------------------------------------

      IF ( ionode ) THEN

         ! 
         CALL iotk_write_begin( iunpun, "BAND_STRUCTURE_INFO" )
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_ATOMIC_WFC", n_atom_wfc() )
         !
         nelec = nelt
         !
         IF ( nspin == 2 ) THEN
            !
            CALL iotk_write_attr( attr, "UP", nel(1), FIRST = .TRUE. )
            CALL iotk_write_attr( attr, "DW", nel(2) )
            CALL iotk_write_dat( iunpun, &
                                 "NUMBER_OF_ELECTRONS", nelec, ATTR = attr )
            !
            CALL iotk_write_attr( attr, "UP", nupdwn_tot(1), FIRST = .TRUE. )
            CALL iotk_write_attr( attr, "DW", nupdwn_tot(2) )
            CALL iotk_write_dat( iunpun, &
                                 "NUMBER_OF_BANDS", nbnd_tot , ATTR = attr )
            !
         ELSE
            !
            CALL iotk_write_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec )
            !
            CALL iotk_write_dat( iunpun, "NUMBER_OF_BANDS", nbnd_tot )
            !
         END IF
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_EMPTY_STATES", nbnd_emp )
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_SPIN_COMPONENTS", nspin )
         !
         CALL iotk_write_end( iunpun, "BAND_STRUCTURE_INFO" )
         !
         CALL iotk_write_begin( iunpun, "EIGENVALUES" )
         !
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... EIGENVALUES
!-------------------------------------------------------------------------------
      !
      k_points_loop1: DO ik = 1, nk
         !
         IF ( ionode ) THEN
            !
            CALL iotk_write_begin( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
            !
            CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
            CALL iotk_write_dat( iunpun, &
                                 "K-POINT_COORDS", xk(:,ik), ATTR = attr )
            !
            CALL iotk_write_dat( iunpun, "WEIGHT", wk(ik) )
            !
            ALLOCATE( dtmp ( nbnd_tot ) )
            !
            DO iss = 1, nspin
               !
               cspin = iotk_index( iss )
               !
               dtmp = 0.0d0
               !
               IF( tksw ) THEN
                  ! 
                  !  writes data required by postproc and PW
                  !
                  IF( nspin == 2 ) THEN
                     IF( iss == 1 ) filename = wfc_filename( ".", 'eigenval1', ik, EXTENSION='xml' )
                     IF( iss == 2 ) filename = wfc_filename( ".", 'eigenval2', ik, EXTENSION='xml' )
                     !
                     IF( iss == 1 ) CALL iotk_link( iunpun, "DATAFILE.1", &
                                                    filename, CREATE = .FALSE., BINARY = .FALSE. )
                     IF( iss == 2 ) CALL iotk_link( iunpun, "DATAFILE.2", &
                                                    filename, CREATE = .FALSE., BINARY = .FALSE. )
   
                     IF( iss == 1 ) filename = wfc_filename( dirname, 'eigenval1', ik, EXTENSION='xml' )
                     IF( iss == 2 ) filename = wfc_filename( dirname, 'eigenval2', ik, EXTENSION='xml' )
                  ELSE
                     filename = wfc_filename( ".", 'eigenval', ik, EXTENSION='xml' )
                     CALL iotk_link( iunpun, "DATAFILE", filename, CREATE = .FALSE., BINARY = .FALSE. )
                     filename = wfc_filename( dirname, 'eigenval', ik, EXTENSION='xml' )
                  END IF

                  dtmp ( 1:nupdwn( iss ) ) = occ0( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) / wk(ik)
                  !
                  CALL write_eig( iunout, filename, nbnd_tot, et( 1:nbnd_tot, iss) , "Hartree", &
                               OCC = dtmp(:), IK=ik, ISPIN=iss )
                  !
               END IF
               !
               CALL iotk_write_dat( iunpun, "OCC0"  // TRIM( cspin ), &
                                    occ0( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) )
               !
               CALL iotk_write_dat( iunpun, "OCCM" // TRIM( cspin ), &
                                    occm( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) )
               !
            END DO
            !
            DEALLOCATE( dtmp )
            !
            CALL iotk_write_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )

         END IF
         !
      END DO k_points_loop1
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_end( iunpun, "EIGENVALUES" )
         !
         CALL iotk_write_begin( iunpun, "EIGENVECTORS" )
         !
         CALL iotk_write_dat  ( iunpun, "MAX_NUMBER_OF_GK-VECTORS", ngwt )
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... EIGENVECTORS
!-------------------------------------------------------------------------------
      !
      k_points_loop2: DO ik = 1, nk

         IF( ionode ) THEN

            CALL iotk_write_begin( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
            ! ... G+K vectors
            !
            CALL iotk_write_dat( iunpun, "NUMBER_OF_GK-VECTORS", ngwt )
            !
            !
            filename = TRIM( wfc_filename( ".", 'gkvectors', ik ) )
            !
            CALL iotk_link( iunpun, "GK-VECTORS", filename, CREATE = .FALSE., BINARY = .TRUE. )
            !
            filename = TRIM( wfc_filename( dirname, 'gkvectors', ik ) )
            !
         END IF
         !
         CALL write_gk( iunout, ik, mill, filename )
         !
         DO iss = 1, nspin
            ! 
            ik_eff = ik + ( iss - 1 ) * nk
            ! 
            iss_wfc = iss
            if( force_pairing ) iss_wfc = 1   ! only the WF for the first spin is allocated
            !
            IF( tksw ) THEN 
               ! 
               !   Save additional WF, 
               !   orthogonal KS states to be used for post processing and PW
               ! 
               IF ( ionode ) THEN
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( ".", 'evc', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( ".", 'evc', ik, iss ) )
                     !
                  END IF
                  !
                  IF( nspin == 2 ) THEN
                     CALL iotk_link( iunpun, "WFC" // TRIM( iotk_index (iss) ), &
                                     filename, CREATE = .FALSE., BINARY = .TRUE. )
                  ELSE
                     CALL iotk_link( iunpun, "WFC", filename, CREATE = .FALSE., BINARY = .TRUE. )
                  END IF
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc', ik, iss ) )
                  !
                  END IF
                  !
               END IF
               !
               ib = iupdwn_tot( iss_wfc )
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin,         & !added_giovanni do_wf_cmplx
                               ctot( :, ib : ib + nbnd_tot - 1 ), ngwt, do_wf_cmplx, &
                               gamma_only, nbnd_tot, ig_l2g, ngw, filename, scalef)
               !
            END IF
            !
            !  Save wave function at time t
            !
            IF ( ionode ) THEN
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( ".", 'evc0', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( ".", 'evc0', ik, iss ) )
                  !
               END IF
               !
               CALL iotk_link( iunpun, "WFC0" // TRIM( iotk_index (iss) ), &
                               filename, CREATE = .FALSE., BINARY = .TRUE. )
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( dirname, 'evc0', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( dirname, 'evc0', ik, iss ) )
                  !
               END IF
               !
            END IF
            !
            ib = iupdwn(iss_wfc)
            !
            IF(.not. evc_restart) then
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin, &
                            c02( :, ib : ib + nbnd_ - 1 ), ngwt, do_wf_cmplx, & !added:giovanni do_wf_cmplx
                            gamma_only, nbnd_, ig_l2g, ngw, filename, scalef)
            ELSE
               WRITE(*,*) "Careful: I am writing Kohn-Sham eigenstates as restart wavefunctions"
               WRITE(*,*) "Errors may happen"
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin, &
                            ctot( :, ib : ib + nbnd_ - 1 ), ngwt, do_wf_cmplx, & !added:giovanni do_wf_cmplx
                            gamma_only, nbnd_, ig_l2g, ngw, filename, scalef)
            ENDIF
            !
            !  Save wave function at time t - dt
            !
            IF ( ionode ) THEN
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( ".", 'evcm', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( ".", 'evcm', ik, iss ) )
                  !
               END IF
               !
               CALL iotk_link( iunpun, "WFCM" // TRIM( iotk_index (iss) ), &
                               filename, CREATE = .FALSE., BINARY = .TRUE. )
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( dirname, 'evcm', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( dirname, 'evcm', ik, iss ) )
                  !
               END IF
               !
            END IF
            !
            ib = iupdwn(iss_wfc)
            !
            CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin,     &
                            cm2( :, ib : ib + nbnd_ - 1 ), ngwt, do_wf_cmplx, & !added:giovanni do_wf_cmplx
                            gamma_only, nbnd_, ig_l2g, ngw, filename, scalef)
            !
            !  Save fixed wave function
            !
            IF (odd_nkscalfact) THEN
               !   
               IF ( ionode ) THEN
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( ".", 'evc0fixed', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( ".", 'evc0fixed', ik, iss ) )
                     !
                  END IF
                  !
                  CALL iotk_link( iunpun, "WFC0FIXED" // TRIM( iotk_index (iss) ), &
                                  filename, CREATE = .FALSE., BINARY = .TRUE. )
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc0fixed', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc0fixed', ik, iss ) )
                     !
                  END IF
                  !
               END IF
               !
               ib = iupdwn(iss_wfc)
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin, &
                               c0_fixed( :, ib : ib + nbnd_ - 1 ), ngwt, do_wf_cmplx, & !added:giovanni do_wf_cmplx
                               gamma_only, nbnd_, ig_l2g, ngw, filename, scalef)
               !
            ENDIF
            !
            !
            IF (print_wfc_empty .and. (nbnd_emp>0) ) THEN
               !
               IF ( ionode ) THEN
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( ".", 'evc0empty', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( ".", 'evc0empty', ik, iss ) )
                     !
                  END IF
                  !
                  CALL iotk_link( iunpun, "WFC0EMPTY" // TRIM( iotk_index (iss) ), &
                                  filename, CREATE = .FALSE., BINARY = .TRUE. )
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc0empty', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc0empty', ik, iss ) )
                     !
                  END IF
                  !
               END IF
               !
               ib = iupdwn(iss_wfc)
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin, &
                               c0_emp_aux( :, ib : ib + nbnd_emp - 1 ), ngwt, do_wf_cmplx, & !added:giovanni do_wf_cmplx
                               gamma_only, nbnd_emp, ig_l2g, ngw, filename, scalef)
               !
               IF ( ionode ) THEN
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( ".", 'evcmempty', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( ".", 'evcmempty', ik, iss ) )
                     !
                  END IF
                  !
                  CALL iotk_link( iunpun, "WFCMEMPTY" // TRIM( iotk_index (iss) ), &
                                  filename, CREATE = .FALSE., BINARY = .TRUE. )
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( dirname, 'evcmempty', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( dirname, 'evcmempty', ik, iss ) )
                     !
                  END IF
                  !
               END IF
               !
               ib = iupdwn(iss_wfc)
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin, &
                               cm_emp_aux( :, ib : ib + nbnd_emp - 1 ), ngwt, do_wf_cmplx, & !added:giovanni do_wf_cmplx
                               gamma_only, nbnd_emp, ig_l2g, ngw, filename, scalef)
               !
            ENDIF
            !
            cspin = iotk_index( iss )
            !
            ! ... write matrix lambda to file
            !
            IF(.not. lambda0(1)%iscmplx) THEN
                ALLOCATE( mrepl( nudx, nudx ) )    
                CALL collect_lambda( mrepl, lambda0(iss)%rvec(:,:), descla(:,iss) )
            ELSE
                ALLOCATE( mrepl_c( nudx, nudx) )
                CALL collect_lambda( mrepl_c, lambda0(iss)%cvec(:,:), descla(:,iss))
            ENDIF
            !
            !
            IF ( ionode ) THEN
               !
               filename = TRIM( wfc_filename( ".", 'lambda0', ik, iss ) )
               !
               CALL iotk_link( iunpun, "LAMBDA0" // TRIM( cspin ), &
                               filename, CREATE = .TRUE., BINARY = .TRUE. )
               !
               IF(.not. lambda0(1)%iscmplx) THEN
		  CALL iotk_write_dat( iunpun, &
					"LAMBDA0" // TRIM( cspin ), mrepl )
               ELSE
                  CALL iotk_write_dat( iunpun, &
					"LAMBDA0" // TRIM( cspin ), mrepl_c )
               ENDIF
               !
               ! Changes by Nicolas Poilvert, Sep. 2010 for printing the lambda
               ! matrix at current time step into a formatted file.
               ! This matrix corresponds to the Hamiltonian matrix  in the case
               ! of Self-Interaction. Only in the basis  of minimizing orbitals
               ! do this matrix has an interpretation.
               !
               IF ( nspin == 1 ) THEN
                   !
                   filename = TRIM( wfc_filename( ".", 'hamiltonian', ik, EXTENSION='xml' ) )
                   !
               ELSE
                   !
                   filename = TRIM( wfc_filename( ".", 'hamiltonian', ik, iss, EXTENSION='xml' ) )
                   !
               ENDIF
               !
               CALL iotk_link( iunpun, "HAMILTONIAN" // TRIM( cspin ), &
                               filename, CREATE = .TRUE., BINARY = .FALSE. )
               !
               IF(allocated(mrepl)) THEN
		  CALL iotk_write_dat( iunpun, &
					"HAMILTONIAN" // TRIM( cspin ), mrepl )
               ELSE IF(allocated(mrepl_c)) THEN
                  CALL iotk_write_dat( iunpun, &
					"HAMILTONIAN" // TRIM( cspin ), mrepl_c )
               ENDIF
               !
            ENDIF
 
            IF(do_bare_eigs) THEN
               !
               IF(.not. lambda_bare(1)%iscmplx) THEN
                  mrepl=0.d0
                  CALL collect_lambda( mrepl, lambda_bare(iss)%rvec(:,:), descla(:,iss) )
               ELSE
                  mrepl_c=0.d0
                  CALL collect_lambda( mrepl_c, lambda_bare(iss)%cvec(:,:), descla(:,iss))
               ENDIF
               !
               IF(ionode) THEN
                  !
                  IF ( nspin == 1 ) THEN
                      !
                      filename = TRIM( wfc_filename( ".", 'hamiltonian0', ik, EXTENSION='xml' ) )
                      !
                  ELSE
                      !
                      filename = TRIM( wfc_filename( ".", 'hamiltonian0', ik, iss, EXTENSION='xml' ) )
                      !
                  ENDIF
                  !
                  CALL iotk_link( iunpun, "HAMILTONIAN0" // TRIM( cspin ), &
                               filename, CREATE = .TRUE., BINARY = .FALSE. )
                  !
                  IF(allocated(mrepl)) THEN
		     CALL iotk_write_dat( iunpun, &
					"HAMILTONIAN0" // TRIM( cspin ), mrepl )
                  ELSE IF(allocated(mrepl_c)) THEN
                     CALL iotk_write_dat( iunpun, &
					"HAMILTONIAN0" // TRIM( cspin ), mrepl_c )
                  ENDIF
                  !
               ENDIF
               !
            END IF
            !
            IF(.not. lambdam(1)%iscmplx) THEN
                CALL collect_lambda( mrepl, lambdam(iss)%rvec(:,:), descla(:,iss) )
            ELSE
                CALL collect_lambda( mrepl_c, lambdam(iss)%cvec(:,:), descla(:,iss))
            ENDIF
            !
            IF ( ionode ) THEN
               !
               filename = TRIM( wfc_filename( ".", 'lambdam', ik, iss ) )
               !
               CALL iotk_link( iunpun, "LAMBDAM" // TRIM( cspin ), &
                               filename, CREATE = .TRUE., BINARY = .TRUE. )
               !
	      IF(allocated(mrepl)) THEN
		CALL iotk_write_dat( iunpun, &
				      "LAMBDAM" // TRIM( cspin ), mrepl )
	      ELSE IF(allocated(mrepl_c)) THEN
		CALL iotk_write_dat( iunpun, &
				      "LAMBDAM" // TRIM( cspin ), mrepl_c )
	      ENDIF

               !
            END IF
            !
	    ! 
            IF( PRESENT( mat_z ) ) THEN
               !
               IF(.not.mat_z(iss)%iscmplx) THEN
		  IF(.not.allocated(mrepl)) THEN
		    ALLOCATE( mrepl( nudx, nudx ) )
		  ENDIF
		  CALL collect_zmat( mrepl, mat_z(iss)%rvec(:,:), descla(:,iss) )
               ELSE
		  IF(.not.allocated(mrepl_c)) THEN
		    ALLOCATE( mrepl_c( nudx, nudx ) )
		  ENDIF
		  CALL collect_zmat( mrepl_c, mat_z(iss)%cvec(:,:), descla(:,iss) )
               ENDIF
               !
               IF ( ionode ) THEN
                  !
                  filename = TRIM( wfc_filename( ".", 'mat_z', ik, iss ) )
                  !
                  CALL iotk_link( iunpun, "MAT_Z" // TRIM( cspin ), &
                                  filename, CREATE = .TRUE., BINARY = .TRUE. )
                  !
                  IF(.not.mat_z(iss)%iscmplx) THEN
                    CALL iotk_write_dat( iunpun, "MAT_Z" // TRIM( cspin ), mrepl )
                  ELSE
		    CALL iotk_write_dat( iunpun, "MAT_Z" // TRIM( cspin ), mrepl_c )
                  ENDIF
                  !
               END IF
               !

            END IF
            !
	    IF(allocated(mrepl)) THEN
	      DEALLOCATE( mrepl )
	    ENDIF
	    IF(allocated(mrepl_c)) THEN
	      DEALLOCATE( mrepl_c )
	    ENDIF
            !
         END DO
         !
         IF ( ionode ) &
            CALL iotk_write_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         !
      END DO k_points_loop2
      !
      IF ( ionode ) CALL iotk_write_end( iunpun, "EIGENVECTORS" )
      !
      IF ( ionode ) CALL iotk_close_write( iunpun )
      !
!-------------------------------------------------------------------------------
! ... END RESTART SECTIONS
!-------------------------------------------------------------------------------
      !
      DEALLOCATE( ftmp )
      DEALLOCATE( tau  )
      DEALLOCATE( ityp )
      DEALLOCATE( mill )
      !
      CALL save_history( dirname, nfi )
      !
      s1 = cclock() 
      !
      IF ( ionode ) THEN
         !
         WRITE( stdout, &
                '(3X,"restart file written in ",F8.3," sec.",/)' ) ( s1 - s0 )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE cp_writefile_twin
    !
    SUBROUTINE cp_writefile_real( ndw, outdir, ascii, nfi, simtime, acc, nk, xk, &
                             wk, ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh,   &
                             taui, cdmi, stau0, svel0, staum, svelm, force,  &
                             vnhp, xnhp0, xnhpm, nhpcl, nhpdim, occ0, occm,  &
                             lambda0,lambdam, xnhe0, xnhem, vnhe, ekincm,    &
                             et, rho, c02, cm2, ctot, iupdwn, nupdwn,        &
                             iupdwn_tot, nupdwn_tot, mat_z )
      !------------------------------------------------------------------------
      !
      USE control_flags,            ONLY : gamma_only, do_wf_cmplx, force_pairing, trhow, tksw !added:giovanni do_wf_cmplx
      USE io_files,                 ONLY : psfile, pseudo_dir
      USE mp_global,                ONLY : intra_image_comm, me_image, nproc_image
      USE printout_base,            ONLY : title
      USE grid_dimensions,          ONLY : nr1, nr2, nr3, nr1x, nr2x, nr3l
      USE smooth_grid_dimensions,   ONLY : nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, ONLY : nr1b, nr2b, nr3b
      USE gvecp,                    ONLY : ngm, ngmt, ecutp, gcutp
      USE gvecs,                    ONLY : ngs, ngst, ecuts, gcuts, dual
      USE gvecw,                    ONLY : ngw, ngwt, ecutw, gcutw
      USE reciprocal_vectors,       ONLY : ig_l2g, mill_l
      USE electrons_base,           ONLY : nspin, nelt, nel, nudx
      USE cell_base,                ONLY : ibrav, alat, celldm, &
                                           symm_type, s_to_r
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, &
                                           pmass, amass, iforce, ind_bck
      USE funct,                    ONLY : get_dft_name
      USE energies,                 ONLY : enthal, ekin, eht, esr, eself, &
                                           epseu, enl, exc, vave
      USE mp_global,                ONLY : nproc, mpime
      USE mp,                       ONLY : mp_sum
      USE fft_base,                 ONLY : dfftp
      USE constants,                ONLY : pi
      USE cp_interfaces,            ONLY : n_atom_wfc
      USE global_version,           ONLY : version_number
      USE cp_main_variables,        ONLY : collect_lambda, descla, collect_zmat
      USE input_parameters,     ONLY : odd_nkscalfact, restart_odd_nkscalfact, print_wfc_empty
      USE wavefunctions_module, ONLY : c0_fixed, c0_emp_aux, cm_emp_aux
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN) :: ndw          !
      CHARACTER(LEN=*),      INTENT(IN) :: outdir       !  directory used to store output and restart files
      LOGICAL,               INTENT(IN) :: ascii        !
      INTEGER,               INTENT(IN) :: nfi          ! index of the current step
      REAL(DP),              INTENT(IN) :: simtime      ! simulated time
      REAL(DP),              INTENT(IN) :: acc(:)       !  
      INTEGER,               INTENT(IN) :: nk           ! number of kpoints
      REAL(DP),              INTENT(IN) :: xk(:,:)      ! k-points coordinates 
      REAL(DP),              INTENT(IN) :: wk(:)        ! k-points weights
      REAL(DP),              INTENT(IN) :: ht(3,3)      ! 
      REAL(DP),              INTENT(IN) :: htm(3,3)     ! 
      REAL(DP),              INTENT(IN) :: htvel(3,3)   ! 
      REAL(DP),              INTENT(IN) :: gvel(3,3)    ! 
      REAL(DP),              INTENT(IN) :: xnhh0(3,3)   ! 
      REAL(DP),              INTENT(IN) :: xnhhm(3,3)   ! 
      REAL(DP),              INTENT(IN) :: vnhh(3,3)    ! 
      REAL(DP),              INTENT(IN) :: taui(:,:)    ! 
      REAL(DP),              INTENT(IN) :: cdmi(:)      ! 
      REAL(DP),              INTENT(IN) :: stau0(:,:)   ! 
      REAL(DP),              INTENT(IN) :: svel0(:,:)   ! 
      REAL(DP),              INTENT(IN) :: staum(:,:)   ! 
      REAL(DP),              INTENT(IN) :: svelm(:,:)   ! 
      REAL(DP),              INTENT(IN) :: force(:,:)   ! 
      REAL(DP),              INTENT(IN) :: xnhp0(:)     ! 
      REAL(DP),              INTENT(IN) :: xnhpm(:)     ! 
      REAL(DP),              INTENT(IN) :: vnhp(:)      ! 
      INTEGER,               INTENT(IN) :: nhpcl        ! 
      INTEGER,               INTENT(IN) :: nhpdim       ! 
      REAL(DP),              INTENT(IN) :: occ0(:)      !  occupations of electronic states
      REAL(DP),              INTENT(IN) :: occm(:)      ! 
      REAL(DP),              INTENT(IN) :: lambda0(:,:,:) ! 
      REAL(DP),              INTENT(IN) :: lambdam(:,:,:) ! 
      REAL(DP),              INTENT(IN) :: xnhe0        ! 
      REAL(DP),              INTENT(IN) :: xnhem        ! 
      REAL(DP),              INTENT(IN) :: vnhe         ! 
      REAL(DP),              INTENT(IN) :: ekincm       ! 
      REAL(DP),              INTENT(IN) :: et(:,:)      !  eigenvalues
      REAL(DP),              INTENT(IN) :: rho(:,:)     ! 
      COMPLEX(DP),           INTENT(IN) :: c02(:,:)     ! 
      COMPLEX(DP),           INTENT(IN) :: cm2(:,:)     ! 
      COMPLEX(DP),           INTENT(IN) :: ctot(:,:)    ! 
      INTEGER,               INTENT(IN) :: iupdwn(:)    ! 
      INTEGER,               INTENT(IN) :: nupdwn(:)    ! 
      INTEGER,               INTENT(IN) :: iupdwn_tot(:)! 
      INTEGER,               INTENT(IN) :: nupdwn_tot(:)! 
      REAL(DP),    OPTIONAL, INTENT(IN) :: mat_z(:,:,:) ! 
      !
      LOGICAL               :: write_charge_density
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=256)    :: dirname, filename, rho_file_base
      CHARACTER(LEN=4)      :: cspin
      INTEGER               :: kunit, ib, ik_eff
      INTEGER               :: k1, k2, k3
      INTEGER               :: nk1, nk2, nk3
      INTEGER               :: j, i, iss, ig, nspin_wfc, iss_wfc
      INTEGER               :: is, ia, isa, ik, ierr
      INTEGER,  ALLOCATABLE :: mill(:,:)
      INTEGER,  ALLOCATABLE :: ftmp(:,:)
      INTEGER,  ALLOCATABLE :: ityp(:)
      REAL(DP), ALLOCATABLE :: tau(:,:)
      REAL(DP), ALLOCATABLE :: dtmp(:)
      REAL(DP), ALLOCATABLE :: rhoaux(:)
      REAL(DP)              :: omega, htm1(3,3), h(3,3)
      REAL(DP)              :: a1(3), a2(3), a3(3)
      REAL(DP)              :: b1(3), b2(3), b3(3)
      REAL(DP)              :: nelec
      REAL(DP)              :: scalef
      LOGICAL               :: lsda
      REAL(DP)              :: s0, s1, cclock
      INTEGER               :: nbnd_tot
      INTEGER               :: nbnd_emp
      INTEGER               :: nbnd_
      REAL(DP), ALLOCATABLE :: mrepl(:,:)
      !
      write_charge_density = trhow
      !
      IF( nspin > 1 .AND. .NOT. force_pairing ) THEN
         !
         !  check if the array storing wave functions is large enought
         !
         IF( SIZE( c02, 2 ) < ( iupdwn( 2 ) + nupdwn(1) - 1 ) ) &
            CALL errore('cp_writefile',' wrong wave functions dimension ', 1 )
         !
      END IF
      !
      IF(  nupdwn_tot(1) < nupdwn(1) ) &
         CALL errore( " writefile ", " wrong number of states ", 1 )
      !
      nbnd_    = nupdwn(1) 
      nbnd_tot = MAX( nupdwn(1), nupdwn_tot(1) )
      nbnd_emp = MAX( 0, nupdwn_tot(1) - nupdwn(1) )
      !
      IF ( ionode ) THEN
         !
         ! ... look for an empty unit (only ionode needs it)
         !
         CALL iotk_free_unit( iunout, ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_writefile', &
                   'no free units to write wavefunctions', ierr )
      !
      dirname = restart_dir( outdir, ndw )
      !
      ! ... Create main restart directory
      !
      CALL create_directory( dirname )
      !
      ! ... Create k-points subdirectories
      ! ... note: in FPMD and CP k-points are not distributed to processors
      !
      DO i = 1, nk
         !
         CALL create_directory( kpoint_dir( dirname, i ) )
         !
      END DO
      !
      ! ... Some ( CP/FPMD ) default values
      !
      IF ( nspin == 2 ) THEN
         !
         kunit = 2
         !
      ELSE
         !
         kunit = 1
         !
      END IF
      !
      k1  = 0
      k2  = 0
      k3  = 0
      nk1 = 0
      nk2 = 0
      nk3 = 0
      !
      ! ... Compute Cell related variables
      !
      h = TRANSPOSE( ht )
      !
      CALL invmat( 3, ht, htm1, omega )
      !
      a1 = ht(1,:)
      a2 = ht(2,:)
      a3 = ht(3,:)
      !
      ! ... Beware: omega may be negative if axis are left-handed!
      !
      scalef = 1.D0 / SQRT( ABS (omega) )
      !
      ! ... Compute array ityp, and tau
      !
      ALLOCATE( ityp( nat ) )
      ALLOCATE( tau( 3, nat ) )
      !
      isa = 0
      !
      DO is = 1, nsp
         !
         DO ia = 1, na(is)
            !
            isa = isa + 1
            ityp(isa) = is
            !
         END DO
         !
      END DO
      !
      CALL s_to_r( stau0, tau, na, nsp, h )
      !   
      ! ... Collect G vectors
      !   
      ALLOCATE( mill( 3, ngmt ) )
      !
      mill = 0
      !
      mill(:,ig_l2g(1:ngm)) = mill_l(:,1:ngm)
      !
      CALL mp_sum( mill, intra_image_comm )
      !
      lsda = ( nspin == 2 )
      !
      ALLOCATE( ftmp( nbnd_tot , nspin ) )
      !
      ftmp = 0.0d0
      !
      DO iss = 1, nspin
         !
         ftmp( 1:nupdwn(iss), iss ) = occ0( iupdwn(iss) : iupdwn(iss) + nupdwn(iss) - 1 )
         !
      END DO
      !
      IF ( ionode ) THEN
         !
         ! ... Open XML descriptor
         !
         WRITE( stdout, '(/,3X,"writing restart file: ",A)' ) TRIM( dirname )
         !
         CALL iotk_open_write( iunpun, FILE = TRIM( dirname ) // '/' // &
                             & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_writefile ', 'cannot open restart file for writing', ierr )
      !
      s0 = cclock() 
      !
      IF ( ionode ) THEN

!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
         CALL write_header( "CP", TRIM(version_number) )
         !
!-------------------------------------------------------------------------------
! ... this flag is used to check if the file can be used for post-processing
!-------------------------------------------------------------------------------
         !
         CALL write_control( PP_CHECK_FLAG=.TRUE. )
         !
!-------------------------------------------------------------------------------
! ... STATUS
!-------------------------------------------------------------------------------
         !
         CALL iotk_write_begin( iunpun, "STATUS" )
         !
         CALL iotk_write_attr( attr, "ITERATION", nfi, FIRST = .TRUE. )
         CALL iotk_write_empty(iunpun, "STEP", attr )
         !
         CALL iotk_write_attr( attr, "UNITS", "pico-seconds", FIRST = .TRUE. ) 
         CALL iotk_write_dat( iunpun, "TIME", simtime, ATTR = attr )
         !
         CALL iotk_write_dat( iunpun, "TITLE", TRIM( title ) )
         !
         CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
         CALL iotk_write_dat( iunpun, "KINETIC_ENERGY", ekin,   ATTR = attr )
         CALL iotk_write_dat( iunpun, "HARTREE_ENERGY", eht,    ATTR = attr )
         CALL iotk_write_dat( iunpun, "EWALD_TERM",     esr,    ATTR = attr )
         CALL iotk_write_dat( iunpun, "GAUSS_SELFINT",  eself,  ATTR = attr )
         CALL iotk_write_dat( iunpun, "LPSP_ENERGY",    epseu,  ATTR = attr )
         CALL iotk_write_dat( iunpun, "NLPSP_ENERGY",   enl,    ATTR = attr )
         CALL iotk_write_dat( iunpun, "EXC_ENERGY",     exc,    ATTR = attr )
         CALL iotk_write_dat( iunpun, "AVERAGE_POT",    vave,   ATTR = attr )
         CALL iotk_write_dat( iunpun, "ENTHALPY",       enthal, ATTR = attr )
         !
         CALL iotk_write_end( iunpun, "STATUS" )      
         !
!-------------------------------------------------------------------------------
! ... CELL
!-------------------------------------------------------------------------------
         !
         a1 = a1 / alat
         a2 = a2 / alat
         a3 = a3 / alat
         !
         CALL recips( a1, a2, a3, b1, b2, b3 )
         !
         CALL write_cell( ibrav, symm_type, &
                          celldm, alat, a1, a2, a3, b1, b2, b3 )
         !
!-------------------------------------------------------------------------------
! ... IONS
!-------------------------------------------------------------------------------
         !
         CALL write_ions( nsp, nat, atm, ityp(ind_bck(:)), &
                          psfile, pseudo_dir, amass, tau(:,ind_bck(:)), &
                          iforce(:,ind_bck(:)), dirname, 1.D0 )
         !
!-------------------------------------------------------------------------------
! ... PLANE_WAVES
!-------------------------------------------------------------------------------
         !
         ! change to .TRUE. to write gvectors.dat for rho
         !
         CALL write_planewaves( ecutw, dual, ngwt, do_wf_cmplx, gamma_only, & 
         !added:giovanni do_wf_cmplx
                                nr1, nr2, nr3, ngmt, nr1s, nr2s, nr3s, ngst, nr1b, &
                                nr2b, nr3b, mill, .FALSE. )
         !
!-------------------------------------------------------------------------------
! ... SPIN
!-------------------------------------------------------------------------------
         !
         CALL write_spin( lsda, .FALSE., 1, .FALSE., .TRUE. )
         !
!-------------------------------------------------------------------------------
! ... EXCHANGE_CORRELATION
!-------------------------------------------------------------------------------
         !
         dft_name = get_dft_name()
         CALL write_xc( DFT = dft_name, NSP = nsp, LDA_PLUS_U = .FALSE. )
         !
!-------------------------------------------------------------------------------
! ... OCCUPATIONS
!-------------------------------------------------------------------------------
         !
         CALL write_occ( LGAUSS = .FALSE., LTETRA = .FALSE., &
                         TFIXED_OCC = .TRUE., LSDA = lsda, NSTATES_UP = nupdwn_tot(1), &
                         NSTATES_DOWN = nupdwn_tot(2), F_INP = DBLE( ftmp ) )
         !
!-------------------------------------------------------------------------------
! ... BRILLOUIN_ZONE
!-------------------------------------------------------------------------------
         !
         CALL write_bz( nk, xk, wk, k1, k2, k3, nk1, nk2, nk3, 0.0_DP )
         !
!-------------------------------------------------------------------------------
! ... PARALLELISM
!-------------------------------------------------------------------------------
         !
         CALL iotk_write_begin( iunpun, "PARALLELISM" )
         !
         CALL iotk_write_dat( iunpun, &
                              "GRANULARITY_OF_K-POINTS_DISTRIBUTION", kunit )
         !
         CALL iotk_write_end( iunpun, "PARALLELISM" )
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... CHARGE-DENSITY
!-------------------------------------------------------------------------------
      !
      IF (write_charge_density) then
         !
         rho_file_base = 'charge-density'
         !
         IF ( ionode )&
              CALL iotk_link( iunpun, "CHARGE-DENSITY", rho_file_base, &
              CREATE = .FALSE., BINARY = .TRUE. )
         !
         rho_file_base = TRIM( dirname ) // '/' // TRIM( rho_file_base )
         !
         IF ( nspin == 1 ) THEN
            !
            CALL write_rho_xml( rho_file_base, rho(:,1), &
                                nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
            !
         ELSE IF ( nspin == 2 ) THEN
            !
            ALLOCATE( rhoaux( SIZE( rho, 1 ) ) )
            !
            rhoaux = rho(:,1) + rho(:,2) 
            !
            CALL write_rho_xml( rho_file_base, rhoaux, &
                                nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
            !
            rho_file_base = 'spin-polarization'
            !
            IF ( ionode ) &
                 CALL iotk_link( iunpun, "SPIN-POLARIZATION", rho_file_base, &
                 CREATE = .FALSE., BINARY = .TRUE. )
            !
            rho_file_base = TRIM( dirname ) // '/' // TRIM( rho_file_base )
            !
            rhoaux = rho(:,1) - rho(:,2) 
            !
            CALL write_rho_xml( rho_file_base, rhoaux, &
                                nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
            !
            DEALLOCATE( rhoaux )
            !
         END IF
         !
      END IF ! write_charge_density
      !
!-------------------------------------------------------------------------------
! ... TIMESTEPS
!-------------------------------------------------------------------------------
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_attr( attr, "nt", 2, FIRST = .TRUE. )
         !
         CALL iotk_write_begin( iunpun, "TIMESTEPS", attr )
         !
         ! ... STEP0
         !
         CALL iotk_write_begin( iunpun, "STEP0" )
         !
         CALL iotk_write_dat( iunpun, "ACCUMULATORS", acc )
         !
         CALL iotk_write_begin( iunpun, "IONS_POSITIONS" )
         CALL iotk_write_dat(   iunpun, "stau",  stau0(1:3,1:nat),   COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "svel",  svel0(1:3,1:nat),   COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "taui",  taui(1:3,1:nat),    COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "cdmi",  cdmi(1:3),          COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "force", force(1:3,1:nat),   COLUMNS=3 )
         CALL iotk_write_end(   iunpun, "IONS_POSITIONS" )
         !
         CALL iotk_write_begin( iunpun, "IONS_NOSE" )
         CALL iotk_write_dat(   iunpun, "nhpcl", nhpcl )
         CALL iotk_write_dat(   iunpun, "nhpdim", nhpdim )
         CALL iotk_write_dat(   iunpun, "xnhp",  xnhp0(1:nhpcl*nhpdim) )
         CALL iotk_write_dat(   iunpun, "vnhp",  vnhp(1:nhpcl*nhpdim) )
         CALL iotk_write_end(   iunpun, "IONS_NOSE" )
         !
         CALL iotk_write_dat( iunpun, "ekincm", ekincm )
         !
         CALL iotk_write_begin( iunpun, "ELECTRONS_NOSE" )
         CALL iotk_write_dat(   iunpun, "xnhe", xnhe0 )
         CALL iotk_write_dat(   iunpun, "vnhe", vnhe )
         CALL iotk_write_end(   iunpun, "ELECTRONS_NOSE" )
         !
         CALL iotk_write_begin( iunpun, "CELL_PARAMETERS" )
         CALL iotk_write_dat(   iunpun, "ht",    ht )
         CALL iotk_write_dat(   iunpun, "htvel", htvel )
         CALL iotk_write_dat(   iunpun, "gvel",  gvel )
         CALL iotk_write_end(   iunpun, "CELL_PARAMETERS" )
         !
         CALL iotk_write_begin( iunpun, "CELL_NOSE" )
         CALL iotk_write_dat(   iunpun, "xnhh", xnhh0 )
         CALL iotk_write_dat(   iunpun, "vnhh", vnhh )
         CALL iotk_write_end(   iunpun, "CELL_NOSE" )
         !
         CALL iotk_write_end( iunpun, "STEP0" )
         !
         ! ... STEPM
         !
         CALL iotk_write_begin( iunpun, "STEPM" )
         !
         CALL iotk_write_begin( iunpun, "IONS_POSITIONS" )
         CALL iotk_write_dat(   iunpun, "stau", staum(1:3,1:nat),  COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "svel", svelm(1:3,1:nat),  COLUMNS=3 )
         CALL iotk_write_end(   iunpun, "IONS_POSITIONS" )
         !
         CALL iotk_write_begin( iunpun, "IONS_NOSE" )
         CALL iotk_write_dat(   iunpun, "nhpcl", nhpcl )
         CALL iotk_write_dat(   iunpun, "nhpdim", nhpdim )
         CALL iotk_write_dat(   iunpun, "xnhp",  xnhpm(1:nhpcl*nhpdim) )
         CALL iotk_write_end(   iunpun, "IONS_NOSE" )
         !
         CALL iotk_write_begin( iunpun, "ELECTRONS_NOSE" )
         CALL iotk_write_dat(   iunpun, "xnhe", xnhem )
         CALL iotk_write_end(   iunpun, "ELECTRONS_NOSE" )
         !
         CALL iotk_write_begin( iunpun, "CELL_PARAMETERS" )
         CALL iotk_write_dat(   iunpun, "ht",    htm )
         CALL iotk_write_end(   iunpun, "CELL_PARAMETERS" )
         !
         CALL iotk_write_begin( iunpun, "CELL_NOSE" )
         CALL iotk_write_dat(   iunpun, "xnhh", xnhhm )
         CALL iotk_write_end(   iunpun, "CELL_NOSE" )
         !
         CALL iotk_write_end( iunpun, "STEPM" )
         !
         CALL iotk_write_end( iunpun, "TIMESTEPS" )
         !
      END IF

!-------------------------------------------------------------------------------
! ... BAND_STRUCTURE_INFO
!-------------------------------------------------------------------------------

      IF ( ionode ) THEN

         ! 
         CALL iotk_write_begin( iunpun, "BAND_STRUCTURE_INFO" )
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_ATOMIC_WFC", n_atom_wfc() )
         !
         nelec = nelt
         !
         IF ( nspin == 2 ) THEN
            !
            CALL iotk_write_attr( attr, "UP", nel(1), FIRST = .TRUE. )
            CALL iotk_write_attr( attr, "DW", nel(2) )
            CALL iotk_write_dat( iunpun, &
                                 "NUMBER_OF_ELECTRONS", nelec, ATTR = attr )
            !
            CALL iotk_write_attr( attr, "UP", nupdwn_tot(1), FIRST = .TRUE. )
            CALL iotk_write_attr( attr, "DW", nupdwn_tot(2) )
            CALL iotk_write_dat( iunpun, &
                                 "NUMBER_OF_BANDS", nbnd_tot , ATTR = attr )
            !
         ELSE
            !
            CALL iotk_write_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec )
            !
            CALL iotk_write_dat( iunpun, "NUMBER_OF_BANDS", nbnd_tot )
            !
         END IF
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_EMPTY_STATES", nbnd_emp )
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_SPIN_COMPONENTS", nspin )
         !
         CALL iotk_write_end( iunpun, "BAND_STRUCTURE_INFO" )
         !
         CALL iotk_write_begin( iunpun, "EIGENVALUES" )
         !
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... EIGENVALUES
!-------------------------------------------------------------------------------
      !
      k_points_loop1: DO ik = 1, nk
         !
         IF ( ionode ) THEN
            !
            CALL iotk_write_begin( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
            !
            CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
            CALL iotk_write_dat( iunpun, &
                                 "K-POINT_COORDS", xk(:,ik), ATTR = attr )
            !
            CALL iotk_write_dat( iunpun, "WEIGHT", wk(ik) )
            !
            ALLOCATE( dtmp ( nbnd_tot ) )
            !
            DO iss = 1, nspin
               !
               cspin = iotk_index( iss )
               !
               dtmp = 0.0d0
               !
               IF( tksw ) THEN
                  ! 
                  !  writes data required by postproc and PW
                  !
                  IF( nspin == 2 ) THEN
                     IF( iss == 1 ) filename = wfc_filename( ".", 'eigenval1', ik, EXTENSION='xml' )
                     IF( iss == 2 ) filename = wfc_filename( ".", 'eigenval2', ik, EXTENSION='xml' )
                     !
                     IF( iss == 1 ) CALL iotk_link( iunpun, "DATAFILE.1", &
                                                    filename, CREATE = .FALSE., BINARY = .FALSE. )
                     IF( iss == 2 ) CALL iotk_link( iunpun, "DATAFILE.2", &
                                                    filename, CREATE = .FALSE., BINARY = .FALSE. )
   
                     IF( iss == 1 ) filename = wfc_filename( dirname, 'eigenval1', ik, EXTENSION='xml' )
                     IF( iss == 2 ) filename = wfc_filename( dirname, 'eigenval2', ik, EXTENSION='xml' )
                  ELSE
                     filename = wfc_filename( ".", 'eigenval', ik, EXTENSION='xml' )
                     CALL iotk_link( iunpun, "DATAFILE", filename, CREATE = .FALSE., BINARY = .FALSE. )
                     filename = wfc_filename( dirname, 'eigenval', ik, EXTENSION='xml' )
                  END IF

                  dtmp ( 1:nupdwn( iss ) ) = occ0( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) / wk(ik)
                  !
                  CALL write_eig( iunout, filename, nbnd_tot, et( 1:nbnd_tot, iss) , "Hartree", &
                               OCC = dtmp(:), IK=ik, ISPIN=iss )
                  !
               END IF
               !
               CALL iotk_write_dat( iunpun, "OCC0"  // TRIM( cspin ), &
                                    occ0( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) )
               !
               CALL iotk_write_dat( iunpun, "OCCM" // TRIM( cspin ), &
                                    occm( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) )
               !
            END DO
            !
            DEALLOCATE( dtmp )
            !
            CALL iotk_write_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )

         END IF
         !
      END DO k_points_loop1
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_end( iunpun, "EIGENVALUES" )
         !
         CALL iotk_write_begin( iunpun, "EIGENVECTORS" )
         !
         CALL iotk_write_dat  ( iunpun, "MAX_NUMBER_OF_GK-VECTORS", ngwt )
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... EIGENVECTORS
!-------------------------------------------------------------------------------
      !
      k_points_loop2: DO ik = 1, nk

         IF( ionode ) THEN

            CALL iotk_write_begin( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
            ! ... G+K vectors
            !
            CALL iotk_write_dat( iunpun, "NUMBER_OF_GK-VECTORS", ngwt )
            !
            !
            filename = TRIM( wfc_filename( ".", 'gkvectors', ik ) )
            !
            CALL iotk_link( iunpun, "GK-VECTORS", filename, CREATE = .FALSE., BINARY = .TRUE. )
            !
            filename = TRIM( wfc_filename( dirname, 'gkvectors', ik ) )
            !
         END IF
         !
         CALL write_gk( iunout, ik, mill, filename )
         !
         DO iss = 1, nspin
            ! 
            ik_eff = ik + ( iss - 1 ) * nk
            ! 
            iss_wfc = iss
            if( force_pairing ) iss_wfc = 1   ! only the WF for the first spin is allocated
            !
            IF( tksw ) THEN 
               ! 
               !   Save additional WF, 
               !   orthogonal KS states to be used for post processing and PW
               ! 
               IF ( ionode ) THEN
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( ".", 'evc', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( ".", 'evc', ik, iss ) )
                     !
                  END IF
                  !
                  IF( nspin == 2 ) THEN
                     CALL iotk_link( iunpun, "WFC" // TRIM( iotk_index (iss) ), &
                                     filename, CREATE = .FALSE., BINARY = .TRUE. )
                  ELSE
                     CALL iotk_link( iunpun, "WFC", filename, CREATE = .FALSE., BINARY = .TRUE. )
                  END IF
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc', ik, iss ) )
                     !
                  END IF
                  !
               END IF
               !
               ib = iupdwn_tot( iss_wfc )
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin,        &
                               ctot( :, ib : ib + nbnd_tot - 1 ), ngwt, do_wf_cmplx, gamma_only,& !added:giovanni do_wf_cmplx
                               nbnd_tot, ig_l2g, ngw, filename, scalef )
               !
            END IF
            !
            !  Save wave function at time t
            !
            IF ( ionode ) THEN
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( ".", 'evc0', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( ".", 'evc0', ik, iss ) )
                  !
               END IF
               !
               CALL iotk_link( iunpun, "WFC0" // TRIM( iotk_index (iss) ), &
                               filename, CREATE = .FALSE., BINARY = .TRUE. )
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( dirname, 'evc0', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( dirname, 'evc0', ik, iss ) )
                  !
               END IF
               !
            END IF
            !
            ib = iupdwn(iss_wfc)
            !
            CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin,     &
                            c02( :, ib : ib + nbnd_ - 1 ), ngwt, do_wf_cmplx, gamma_only, & !added:giovanni do_wf_cmplx
                            nbnd_, ig_l2g, ngw, filename, scalef )
            !
            !  Save wave function at time t - dt
            !
            IF ( ionode ) THEN
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( ".", 'evcm', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( ".", 'evcm', ik, iss ) )
                  !
               END IF
               !
               CALL iotk_link( iunpun, "WFCM" // TRIM( iotk_index (iss) ), &
                               filename, CREATE = .FALSE., BINARY = .TRUE. )
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( dirname, 'evcm', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( dirname, 'evcm', ik, iss ) )
                  !
               END IF
               !
            END IF
            !
            ib = iupdwn(iss_wfc)
            !
            CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin,     &
                            cm2( :, ib : ib + nbnd_ - 1 ), ngwt, do_wf_cmplx, gamma_only, &
                            nbnd_, ig_l2g, ngw, filename, scalef )
            !
            !  Save fixed wave function
            !
            IF (odd_nkscalfact) THEN
               !   
               IF ( ionode ) THEN
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( ".", 'evc0fixed', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( ".", 'evc0fixed', ik, iss ) )
                     !
                  END IF
                  !
                  CALL iotk_link( iunpun, "WFC0FIXED" // TRIM( iotk_index (iss) ), &
                                  filename, CREATE = .FALSE., BINARY = .TRUE. )
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc0fixed', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc0fixed', ik, iss ) )
                     !
                  END IF
                  !
               END IF
               !
               ib = iupdwn(iss_wfc)
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin, &
                               c0_fixed( :, ib : ib + nbnd_ - 1 ), ngwt, do_wf_cmplx, & !added:giovanni do_wf_cmplx
                               gamma_only, nbnd_, ig_l2g, ngw, filename, scalef)
               !
            ENDIF
            !
            IF (print_wfc_empty .and. (nbnd_emp>0) ) THEN
               !
               IF ( ionode ) THEN
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( ".", 'evc0empty', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( ".", 'evc0empty', ik, iss ) )
                     !
                  END IF
                  !
                  CALL iotk_link( iunpun, "WFC0EMPTY" // TRIM( iotk_index (iss) ), &
                                  filename, CREATE = .FALSE., BINARY = .TRUE. )
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc0empty', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc0empty', ik, iss ) )
                     !
                  END IF
                  !
               END IF
               !
               ib = iupdwn(iss_wfc)
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin, &
                               c0_emp_aux( :, ib : ib + nbnd_emp - 1 ), ngwt, do_wf_cmplx, & !added:giovanni do_wf_cmplx
                               gamma_only, nbnd_emp, ig_l2g, ngw, filename, scalef)
               !
               IF ( ionode ) THEN
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( ".", 'evcmempty', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( ".", 'evcmempty', ik, iss ) )
                     !
                  END IF
                  !
                  CALL iotk_link( iunpun, "WFCMEMPTY" // TRIM( iotk_index (iss) ), &
                                  filename, CREATE = .FALSE., BINARY = .TRUE. )
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( dirname, 'evcmempty', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( dirname, 'evcmempty', ik, iss ) )
                     !
                  END IF
                  !
               END IF
               !
               ib = iupdwn(iss_wfc)
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin, &
                               cm_emp_aux( :, ib : ib + nbnd_emp - 1 ), ngwt, do_wf_cmplx, & !added:giovanni do_wf_cmplx
                               gamma_only, nbnd_emp, ig_l2g, ngw, filename, scalef)
               !
            ENDIF
            !
            cspin = iotk_index( iss )
            !
            ! ... write matrix lambda to file
            !
            ALLOCATE( mrepl( nudx, nudx ) )
            !
            CALL collect_lambda( mrepl, lambda0(:,:,iss), descla(:,iss) )
            !
            IF ( ionode ) THEN
               !
               filename = TRIM( wfc_filename( ".", 'lambda0', ik, iss ) )
               !
               CALL iotk_link( iunpun, "LAMBDA0" // TRIM( cspin ), &
                               filename, CREATE = .TRUE., BINARY = .TRUE. )
               !
               CALL iotk_write_dat( iunpun, &
                                    "LAMBDA0" // TRIM( cspin ), mrepl )
               !
               ! Changes by Nicolas Poilvert, Sep. 2010 for printing the lambda
               ! matrix at current time step into a formatted file.
               ! This matrix corresponds to the Hamiltonian matrix  in the case
               ! of Self-Interaction. Only in the basis  of minimizing orbitals
               ! do this matrix has an interpretation.
               !
               IF ( nspin == 1 ) THEN
                   !
                   filename = TRIM( wfc_filename( ".", 'hamiltonian', ik, EXTENSION='xml' ) )
                   !
               ELSE
                   !
                   filename = TRIM( wfc_filename( ".", 'hamiltonian', ik, iss, EXTENSION='xml' ) )
                   !
               ENDIF
               !
               CALL iotk_link( iunpun, "HAMILTONIAN" // TRIM( cspin ), &
                               filename, CREATE = .TRUE., BINARY = .FALSE. )
               !
               CALL iotk_write_dat( iunpun, &
                                    "HAMILTONIAN" // TRIM( cspin ), mrepl )
               !
            END IF
            !
            CALL collect_lambda( mrepl, lambdam(:,:,iss), descla(:,iss) )
            !
            IF ( ionode ) THEN
               !
               filename = TRIM( wfc_filename( ".", 'lambdam', ik, iss ) )
               !
               CALL iotk_link( iunpun, "LAMBDAM" // TRIM( cspin ), &
                               filename, CREATE = .TRUE., BINARY = .TRUE. )
               !
               CALL iotk_write_dat( iunpun, &
                                    "LAMBDAM" // TRIM( cspin ), mrepl )
               !
            END IF
            !
            IF( PRESENT( mat_z ) ) THEN
               !
               CALL collect_zmat( mrepl, mat_z(:,:,iss), descla(:,iss) )
               !
               IF ( ionode ) THEN
                  !
                  filename = TRIM( wfc_filename( ".", 'mat_z', ik, iss ) )
                  !
                  CALL iotk_link( iunpun, "MAT_Z" // TRIM( cspin ), &
                                  filename, CREATE = .TRUE., BINARY = .TRUE. )
                  !
                  CALL iotk_write_dat( iunpun, "MAT_Z" // TRIM( cspin ), mrepl )
                  !
               END IF
               !
            END IF
            !
            DEALLOCATE( mrepl )
            !
         END DO
         !
         IF ( ionode ) &
            CALL iotk_write_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         !
      END DO k_points_loop2
      !
      IF ( ionode ) CALL iotk_write_end( iunpun, "EIGENVECTORS" )
      !
      IF ( ionode ) CALL iotk_close_write( iunpun )
      !
!-------------------------------------------------------------------------------
! ... END RESTART SECTIONS
!-------------------------------------------------------------------------------
      !
      DEALLOCATE( ftmp )
      DEALLOCATE( tau  )
      DEALLOCATE( ityp )
      DEALLOCATE( mill )
      !
      CALL save_history( dirname, nfi )
      !
      s1 = cclock() 
      !
      IF ( ionode ) THEN
         !
         WRITE( stdout, &
                '(3X,"restart file written in ",F8.3," sec.",/)' ) ( s1 - s0 )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE cp_writefile_real

    !------------------------------------------------------------------------
    SUBROUTINE cp_readfile_real( ndr, outdir, ascii, nfi, simtime, acc, nk, xk,   &
                            wk, ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh,     &
                            taui, cdmi, stau0, svel0, staum, svelm, force,    &
                            vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ0, occm,      &
                            lambda0, lambdam, b1, b2, b3, xnhe0, xnhem, vnhe, &
                            ekincm, c02, cm2, mat_z )
      !------------------------------------------------------------------------
      !
      USE control_flags,            ONLY : do_wf_cmplx, gamma_only, force_pairing !added:giovanni do_wf_cmplx
      USE io_files,                 ONLY : iunpun, xmlpun
      USE printout_base,            ONLY : title
      USE grid_dimensions,          ONLY : nr1, nr2, nr3
      USE smooth_grid_dimensions,   ONLY : nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, ONLY : nr1b, nr2b, nr3b
      USE gvecp,                    ONLY : ngm, ngmt, ecutp
      USE gvecs,                    ONLY : ngs, ngst
      USE gvecw,                    ONLY : ngw, ngwt, ecutw
      USE electrons_base,           ONLY : nspin, nbnd, nelt, nel, &
                                           nupdwn, iupdwn, nudx
      USE cell_base,                ONLY : ibrav, alat, celldm, symm_type, &
                                           s_to_r, r_to_s
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, pmass, &
                                           sort_tau, ityp, ions_cofmass
      USE reciprocal_vectors,       ONLY : ig_l2g, mill_l
      USE cp_main_variables,        ONLY : nprint_nfi, distribute_lambda, descla, distribute_zmat
      USE mp,                       ONLY : mp_sum
      USE mp_global,                ONLY : intra_image_comm
      USE parameters,               ONLY : ntypx
      USE constants,                ONLY : eps8, angstrom_au, pi
      USE input_parameters,     ONLY : odd_nkscalfact, restart_odd_nkscalfact
      USE wavefunctions_module, ONLY : c0_fixed
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN)    :: ndr          !  I/O unit number
      CHARACTER(LEN=*),      INTENT(IN)    :: outdir       !
      LOGICAL,               INTENT(IN)    :: ascii        !
      INTEGER,               INTENT(INOUT) :: nfi          ! index of the current step
      REAL(DP),              INTENT(INOUT) :: simtime      ! simulated time
      REAL(DP),              INTENT(INOUT) :: acc(:)       !
      INTEGER,               INTENT(IN)    :: nk           ! number of kpoints
      REAL(DP),              INTENT(INOUT) :: xk(:,:)      ! k-points coordinates
      REAL(DP),              INTENT(INOUT) :: wk(:)        ! k-points weights
      REAL(DP),              INTENT(INOUT) :: ht(3,3)      !
      REAL(DP),              INTENT(INOUT) :: htm(3,3)     !
      REAL(DP),              INTENT(INOUT) :: htvel(3,3)   !
      REAL(DP),              INTENT(INOUT) :: gvel(3,3)    !
      REAL(DP),              INTENT(INOUT) :: xnhh0(3,3)   !
      REAL(DP),              INTENT(INOUT) :: xnhhm(3,3)   !
      REAL(DP),              INTENT(INOUT) :: vnhh(3,3)    !
      REAL(DP),              INTENT(INOUT) :: taui(:,:)    !
      REAL(DP),              INTENT(INOUT) :: cdmi(:)      !
      REAL(DP),              INTENT(INOUT) :: stau0(:,:)   !
      REAL(DP),              INTENT(INOUT) :: svel0(:,:)   !
      REAL(DP),              INTENT(INOUT) :: staum(:,:)   !
      REAL(DP),              INTENT(INOUT) :: svelm(:,:)   !
      REAL(DP),              INTENT(INOUT) :: force(:,:)   ! 
      REAL(DP),              INTENT(INOUT) :: xnhp0(:)     !      
      REAL(DP),              INTENT(INOUT) :: xnhpm(:)     ! 
      REAL(DP),              INTENT(INOUT) :: vnhp(:)      !  
      INTEGER,               INTENT(INOUT) :: nhpcl        !  
      INTEGER,               INTENT(INOUT) :: nhpdim       !  
      REAL(DP),              INTENT(INOUT) :: occ0(:)      ! occupations
      REAL(DP),              INTENT(INOUT) :: occm(:)      !
      REAL(DP),              INTENT(INOUT) :: lambda0(:,:,:) !
      REAL(DP),              INTENT(INOUT) :: lambdam(:,:,:) !
      REAL(DP),              INTENT(INOUT) :: b1(3)        !
      REAL(DP),              INTENT(INOUT) :: b2(3)        !
      REAL(DP),              INTENT(INOUT) :: b3(3)        !
      REAL(DP),              INTENT(INOUT) :: xnhe0        !
      REAL(DP),              INTENT(INOUT) :: xnhem        !
      REAL(DP),              INTENT(INOUT) :: vnhe         !  
      REAL(DP),              INTENT(INOUT) :: ekincm       !  
      COMPLEX(DP),           INTENT(INOUT) :: c02(:,:)     ! 
      COMPLEX(DP),           INTENT(INOUT) :: cm2(:,:)     ! 
      REAL(DP),    OPTIONAL, INTENT(INOUT) :: mat_z(:,:,:) ! 
      !
      CHARACTER(LEN=256)   :: dirname, kdirname, filename
      CHARACTER(LEN=5)     :: kindex
      CHARACTER(LEN=4)     :: cspin
      INTEGER              :: strlen
      INTEGER              :: kunit
      INTEGER              :: k1, k2, k3
      INTEGER              :: nk1, nk2, nk3
      INTEGER              :: i, j, iss, ig, nspin_wfc, ierr, ik
      REAL(DP)             :: omega, htm1(3,3), hinv(3,3), scalef
      LOGICAL              :: found
      INTEGER, ALLOCATABLE :: mill(:,:)
      !
      ! ... variables read for testing pourposes
      !
      INTEGER               :: ibrav_
      CHARACTER(LEN=9)      :: symm_type_
      CHARACTER(LEN=3)      :: atm_(ntypx)
      INTEGER               :: nat_, nsp_, na_
      INTEGER               :: nk_, ik_, nt_
      LOGICAL               :: do_wf_cmplx_, gamma_only_ , lsda_ !added:giovanni do_wf_cmplx
      REAL(DP)              :: alat_, a1_(3), a2_(3), a3_(3)
      REAL(DP)              :: pmass_, zv_ 
      REAL(DP)              :: celldm_(6)
      INTEGER               :: iss_, nspin_, ngwt_, nbnd_ , n_emp_ , nbnd_tot
      INTEGER               :: nstates_up_ , nstates_dw_ , ntmp, nel_(2)
      REAL(DP)              :: nelec_ 
      REAL(DP)              :: scalef_
      REAL(DP)              :: wk_
      INTEGER               :: nhpcl_, nhpdim_ 
      INTEGER               :: ib, nb
      INTEGER               :: ik_eff
      REAL(DP)              :: amass_(ntypx)
      INTEGER,  ALLOCATABLE :: ityp_(:) 
      INTEGER,  ALLOCATABLE :: isrt_(:) 
      REAL(DP), ALLOCATABLE :: tau_(:,:) 
      REAL(DP), ALLOCATABLE :: occ_(:) 
      INTEGER,  ALLOCATABLE :: if_pos_(:,:) 
      CHARACTER(LEN=256)    :: psfile_(ntypx)
      CHARACTER(LEN=80)     :: pos_unit
      REAL(DP)              :: s1, s0, cclock
      REAL(DP), ALLOCATABLE :: mrepl(:,:) 
      !
      ! ... look for an empty unit
      !
      CALL iotk_free_unit( iunout, ierr )
      !
      CALL errore( 'cp_readfile', &
                   'no free units to read wavefunctions', ierr )
      !
      kunit = 1
      found = .FALSE.
      !
      dirname = restart_dir( outdir, ndr )
      !
      ! ... Open XML descriptor
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // '/' // TRIM( xmlpun )
         !
         WRITE( stdout, '(/,3X,"reading restart file: ",A)' ) TRIM( dirname )
         !
         CALL iotk_open_read( iunpun, FILE = TRIM( filename ), IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_readfile', &
                   'cannot open restart file for reading', ierr )
      !
      s0 = cclock()
      !
      IF ( ionode ) THEN
         !
         qexml_version = " "
         !
         CALL iotk_scan_begin( iunpun, "HEADER", FOUND=found )
         !
         IF ( found ) THEN
            !
            CALL iotk_scan_empty( iunpun, "FORMAT", ATTR=attr )
            CALL iotk_scan_attr( attr, "VERSION", qexml_version )
            CALL iotk_scan_end( iunpun, "HEADER" )
            !
         ELSE
            !
            qexml_version = TRIM( default_fmt_version )
            !
         ENDIF
         !
         qexml_version_init = .TRUE.
        
         !
         ! init logical variables for versioning
         !
         qexml_version_before_1_4_0 = .FALSE.
         !
         IF ( TRIM( version_compare( qexml_version, "1.4.0" )) == "older" ) &
            qexml_version_before_1_4_0 = .TRUE.
         !
      ENDIF
      !
      CALL mp_bcast( qexml_version,               ionode_id, intra_image_comm )
      CALL mp_bcast( qexml_version_init,          ionode_id, intra_image_comm )
      CALL mp_bcast( qexml_version_before_1_4_0 , ionode_id, intra_image_comm )
      !
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "STATUS", FOUND = found )
         !
         IF ( found ) THEN
            !
            CALL iotk_scan_empty( iunpun, "STEP", attr )
            CALL iotk_scan_attr( attr, "ITERATION", nfi )
            CALL iotk_scan_dat( iunpun, "TIME", simtime )
            CALL iotk_scan_dat( iunpun, "TITLE", title )
            CALL iotk_scan_end( iunpun, "STATUS" )
            !
         END IF
         !
      END IF
      !
      ! ... Read cell and positions
      !
      ALLOCATE( tau_( 3, nat ) )
      ALLOCATE( if_pos_( 3, nat ) )
      ALLOCATE( ityp_( nat ) )
      !
      IF ( ionode ) THEN
         !
         CALL read_cell( ibrav_, symm_type_, celldm_, &
                         alat_, a1_, a2_, a3_, b1, b2, b3 )
         !
         CALL recips( a1_, a2_, a3_, b1, b2, b3 )
         !
      END IF
      !
      IF ( ionode ) THEN
         !
         CALL read_ions( nsp_, nat_, atm_, ityp_, &
                         psfile_, amass_, tau_, if_pos_, pos_unit, ierr )
         !
         IF ( ierr == 0 ) THEN
            !
            IF( nsp_ /= nsp .OR. nat_ /= nat ) ierr = 2
            !
            DO i = 1, nat
               !
               IF ( ityp_(i) /= ityp(i) ) ierr = 3
               !
            END DO
            !
         END IF
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_readfile', &
                   'cannot read positions from restart file', ierr )
      !
      !  Read SPIN infos
      !
      lsda_ = ( nspin == 2 )
      !
      IF( ionode ) THEN
         CALL iotk_scan_begin( iunpun, "SPIN", FOUND = found )
         IF( found ) THEN
            CALL iotk_scan_dat( iunpun, "LSDA", lsda_ )
            CALL iotk_scan_end( iunpun, "SPIN" )
         END IF
      END IF
      !
      CALL mp_bcast( lsda_ , ionode_id, intra_image_comm )
      !
      IF( lsda_ .AND. nspin == 1 ) &
         CALL errore( 'cp_readfile', 'LSDA restart file with a spinless run', ierr )

      !
      !  Read Occupations infos
      !
      nstates_up_ = nupdwn( 1 )
      nstates_dw_ = nupdwn( 2 )

      IF( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "OCCUPATIONS", FOUND = found )
         IF( found ) THEN
            ! 
            CALL iotk_scan_empty( iunpun, "INFO", attr, FOUND = found )
            !
            IF( lsda_ .AND. found ) THEN
               !
               IF ( qexml_version_before_1_4_0 ) THEN
                  !
                  CALL iotk_scan_attr( attr, "nelup",  nstates_up_ )
                  CALL iotk_scan_attr( attr, "neldw",  nstates_dw_ )
                  !
               ELSE
                  !
                  ! current version
                  !
                  CALL iotk_scan_attr( attr, "nstates_up",  nstates_up_ )
                  CALL iotk_scan_attr( attr, "nstates_down",  nstates_dw_ )
                  !
               ENDIF 
                !
            ENDIF
            !
            CALL iotk_scan_end( iunpun, "OCCUPATIONS" )
            !
         ENDIF
      ENDIF
      !
      CALL mp_bcast( nstates_up_ , ionode_id, intra_image_comm )
      CALL mp_bcast( nstates_dw_ , ionode_id, intra_image_comm )
      !
      IF( lsda_ ) THEN
         IF( ( nstates_up_ /=  nupdwn( 1 ) ) .OR. ( nstates_dw_ /=  nupdwn( 2 ) ) ) &
            CALL errore( 'cp_readfile', 'inconsistent number of spin states', ierr )
      END IF

      ! ... read MD timesteps variables
      !
      IF ( ionode ) &
         CALL iotk_scan_begin( iunpun, "TIMESTEPS", attr, FOUND = found )
      ! 
      ierr = 0
      ! 
      IF ( ionode .AND. found ) THEN
         !
         CALL iotk_scan_attr( attr, "nt", nt_ )
         !
         IF ( nt_ > 0 ) THEN
            !
            CALL iotk_scan_begin( iunpun, "STEP0" )
            !
            CALL iotk_scan_dat( iunpun, "ACCUMULATORS", acc )
            !
            CALL iotk_scan_begin( iunpun,"IONS_POSITIONS" )
            CALL iotk_scan_dat(   iunpun, "stau",  stau0(1:3,1:nat) )
            CALL iotk_scan_dat(   iunpun, "svel",  svel0(1:3,1:nat) )
            CALL iotk_scan_dat(   iunpun, "taui",  taui(1:3,1:nat) )
            CALL iotk_scan_dat(   iunpun, "cdmi",  cdmi(1:3) )
            CALL iotk_scan_dat(   iunpun, "force", force(1:3,1:nat) )
            CALL iotk_scan_end(   iunpun, "IONS_POSITIONS" )
            !
            CALL iotk_scan_begin( iunpun, "IONS_NOSE" )
            CALL iotk_scan_dat(   iunpun, "nhpcl", nhpcl_ )
            CALL iotk_scan_dat(   iunpun, "nhpdim", nhpdim_ )
            !
            IF ( nhpcl_ == nhpcl .AND. nhpdim_ == nhpdim ) THEN
               !
               CALL iotk_scan_dat( iunpun, "xnhp", xnhp0(1:nhpcl*nhpdim) )
               CALL iotk_scan_dat( iunpun, "vnhp", vnhp(1:nhpcl*nhpdim) )
               !
            ELSE
               !
               xnhp0(1:nhpcl*nhpdim) = 0.D0
               vnhp(1:nhpcl*nhpdim)  = 0.D0
               !
            END IF
            !
            CALL iotk_scan_end(   iunpun, "IONS_NOSE" )
            !
            CALL iotk_scan_dat( iunpun, "ekincm", ekincm )
            !
            CALL iotk_scan_begin( iunpun, "ELECTRONS_NOSE" )
            CALL iotk_scan_dat(   iunpun, "xnhe", xnhe0 )
            CALL iotk_scan_dat(   iunpun, "vnhe", vnhe )
            CALL iotk_scan_end(   iunpun, "ELECTRONS_NOSE" )
            !
            CALL iotk_scan_begin( iunpun, "CELL_PARAMETERS" )
            CALL iotk_scan_dat(   iunpun, "ht",    ht )
            CALL iotk_scan_dat(   iunpun, "htvel", htvel )
            CALL iotk_scan_dat(   iunpun, "gvel",  gvel )
            CALL iotk_scan_end(   iunpun, "CELL_PARAMETERS" )
            !
            CALL iotk_scan_begin( iunpun, "CELL_NOSE" )
            CALL iotk_scan_dat(   iunpun, "xnhh", xnhh0 )
            CALL iotk_scan_dat(   iunpun, "vnhh", vnhh )
            CALL iotk_scan_end(   iunpun, "CELL_NOSE" )
            !
            CALL iotk_scan_end( iunpun, "STEP0" )
            !
         ELSE
            !
            ierr = 40
            !
            GOTO 100
            !
         END IF
         !
         IF ( nt_ > 1 ) THEN
            !
            CALL iotk_scan_begin( iunpun, "STEPM" )
            !
            CALL iotk_scan_begin( iunpun, "IONS_POSITIONS" )
            CALL iotk_scan_dat(   iunpun, "stau", staum(1:3,1:nat) )
            CALL iotk_scan_dat(   iunpun, "svel", svelm(1:3,1:nat) )
            CALL iotk_scan_end(   iunpun, "IONS_POSITIONS" )
            !
            CALL iotk_scan_begin( iunpun, "IONS_NOSE" )
            CALL iotk_scan_dat(   iunpun, "nhpcl", nhpcl_ )
            CALL iotk_scan_dat(   iunpun, "nhpdim", nhpdim_ )
            !
            IF ( nhpcl_ == nhpcl .AND. nhpdim_ == nhpdim ) THEN
               !
               CALL iotk_scan_dat( iunpun, "xnhp",  xnhpm(1:nhpcl*nhpdim) )
               !
            ELSE
               !
               xnhpm(1:nhpcl*nhpdim) = 0.D0
               !
            END IF
            !
            CALL iotk_scan_end(   iunpun,"IONS_NOSE" )
            !
            CALL iotk_scan_begin( iunpun, "ELECTRONS_NOSE" )
            CALL iotk_scan_dat(   iunpun, "xnhe", xnhem )
            CALL iotk_scan_end(   iunpun, "ELECTRONS_NOSE" )
            !
            CALL iotk_scan_begin( iunpun, "CELL_PARAMETERS" )
            CALL iotk_scan_dat(   iunpun, "ht", htm )
            CALL iotk_scan_end(   iunpun, "CELL_PARAMETERS" )
            !
            CALL iotk_scan_begin( iunpun, "CELL_NOSE" )
            CALL iotk_scan_dat(   iunpun, "xnhh", xnhhm )
            CALL iotk_scan_end(   iunpun, "CELL_NOSE" )
            !
            CALL iotk_scan_end( iunpun, "STEPM" )
            !
         END IF
         !
         CALL iotk_scan_end( iunpun, "TIMESTEPS" )
         !
      ELSE IF ( ionode ) THEN
         !
         ! ... MD time steps not found, try to recover from CELL and POSITIONS
         ! 
         acc = 0.D0
         ! 
         ALLOCATE( isrt_( nat ) )
         !
         SELECT CASE( TRIM( pos_unit ) )
         CASE( "alat" )
            !
            tau_ = tau_ * alat_
            !
         CASE( "Angstrom" )
            !
            tau_ = tau_ * angstrom_au
            !
         CASE DEFAULT
            !
         END SELECT
         !
         CALL sort_tau( taui, isrt_ , tau_ , ityp_ , nat_ , nsp_ )
         ! 
         ht(1,:) = a1_
         ht(2,:) = a2_
         ht(3,:) = a3_
         !
         CALL invmat( 3, ht, htm1, omega )
         !
         hinv = TRANSPOSE( htm1 )
         !
         CALL r_to_s( taui, stau0, na, nsp, hinv )
         !
         CALL ions_cofmass( taui, amass_ , na, nsp, cdmi )
         !
         staum = stau0
         svel0 = 0.D0
         svelm = 0.D0
         force = 0.D0
         !
         htm   = ht
         htvel = 0.D0
         gvel  = 0.D0
         xnhh0 = 0.D0
         vnhh  = 0.D0
         xnhhm = 0.D0
         !
         xnhe0 = 0.D0
         xnhem = 0.D0
         vnhe  = 0.D0
         !
         ekincm = 0.D0
         !
         xnhp0 = 0.D0
         xnhpm = 0.D0
         vnhp  = 0.D0
         !
         DEALLOCATE( isrt_  )
         !
      END IF
      !
 100  CONTINUE
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF( ierr /= 0 ) THEN
         CALL mp_bcast( attr, ionode_id, intra_image_comm )
         CALL errore( 'cp_readfile ', TRIM( attr ), ierr )
      END IF
      !
      DEALLOCATE( tau_  )
      DEALLOCATE( if_pos_ )
      DEALLOCATE( ityp_ )
      !
      ! ... compute the scale factor
      !
      IF ( ionode ) CALL invmat( 3, ht, htm1, omega )
      !
      CALL mp_bcast( omega, ionode_id, intra_image_comm )
      !
      ! ... Beware: omega may be negative if axis are left-handed!
      !
      scalef = 1.D0 / SQRT( ABS( omega ) )
      !
      ! ... band Structure
      !
      IF ( ionode ) THEN
         !
         ierr = 0
         !
         CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE_INFO" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_SPIN_COMPONENTS", nspin_ )
         ! 
         IF ( nspin_ /= nspin ) THEN
            attr = "spin do not match"
            ierr = 31
            GOTO 90
         END IF
         !
         IF ( nspin == 2 ) THEN
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec_, ATTR = attr )
            CALL iotk_scan_attr( attr, "UP", nel_(1) )
            CALL iotk_scan_attr( attr, "DW", nel_(2) )
            !
            IF ( ( nel(1) /= nel_(1) ) .OR. ( nel(2) /= nel_(2) ) .OR. ( NINT( nelec_ ) /= nelt ) ) THEN
               attr = "electrons do not match"
               write(0,*) "from cp_readfile warning: electrons do not match"
               write(6,*) "from cp_readfile warning: electrons do not match"
               !ierr = 33
               GOTO 90
            END IF
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_BANDS", nbnd_tot , ATTR = attr )
            !
         ELSE
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec_ )
            !
            IF ( NINT( nelec_ ) /= nelt ) THEN
               attr = "electrons do not match"
               ierr = 33
               GOTO 90
            END IF
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_BANDS", nbnd_tot )
            !
         END IF
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_EMPTY_STATES", n_emp_, FOUND = found )
         !
         IF( .NOT. found ) n_emp_ = 0
         !
         nbnd_ = nbnd_tot - n_emp_
         !
         IF ( nbnd_ < nupdwn(1) ) THEN
            attr = "nbnd do not match"
            ierr = 32
            GOTO 90
         END IF
         !
         CALL iotk_scan_end( iunpun, "BAND_STRUCTURE_INFO" )
         !
      END IF
      !
 90   CONTINUE
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      IF( ierr /= 0 ) THEN
         CALL mp_bcast( attr, ionode_id, intra_image_comm )
         CALL errore( 'cp_readfile ', TRIM( attr ), ierr )
      END IF
      !
      IF( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "EIGENVALUES" )
         !
      END IF
      !
      k_points_loop1: DO ik = 1, nk
         !
         IF ( ionode ) THEN
            !
            CALL iotk_scan_begin( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
            !
            CALL iotk_scan_dat( iunpun, "WEIGHT", wk_ )
            !
         END IF
         !
         DO iss = 1, nspin
            !
            cspin = iotk_index( iss )
            !
            ik_eff = ik + ( iss - 1 ) * nk
            !
            IF ( ionode ) THEN
               !
               ALLOCATE( occ_ ( MAX( nudx , nbnd_tot ) ) )
               !
               occ_ = 0.0d0
               !
               CALL iotk_scan_dat( iunpun, "OCC0" // TRIM( cspin ), occ_ ( 1 :  nupdwn( iss ) ), FOUND = found )
               !
               IF( .NOT. found ) THEN
                  !
                  IF( nspin == 1 ) THEN
                     CALL iotk_scan_begin( iunpun, "DATAFILE", FOUND = found )
                  ELSE
                     CALL iotk_scan_begin( iunpun, "DATAFILE//TRIM(cspin)", FOUND = found )
                  END IF
                  !
                  CALL iotk_scan_dat  ( iunpun, "OCCUPATIONS", occ_( 1:nbnd_tot ) ) 
                  !
                  IF( nspin == 1 ) THEN
                     CALL iotk_scan_end( iunpun, "DATAFILE" )
                  ELSE
                     CALL iotk_scan_end( iunpun, "DATAFILE//TRIM(cspin)" )
                  END IF
                  !
                  IF( found ) THEN
                     occ0( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) = occ_ ( 1:nupdwn( iss ) ) * wk_
                     occm( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) = occ_ ( 1:nupdwn( iss ) ) * wk_
                  END IF
                  !
               ELSE
                  !
                  occ0( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) = occ_ ( 1:nupdwn( iss ) )
                  !
                  CALL iotk_scan_dat( iunpun, "OCCM" // TRIM( cspin ), occ_ ( 1 :  nupdwn( iss ) ), FOUND = found )
                  !
                  IF( found ) THEN
                     occm( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) = occ_ ( 1:nupdwn( iss ) )
                  END IF
                  !
               END IF
               !
               DEALLOCATE ( occ_ )
               !
            END IF
            !
            CALL mp_bcast( found, ionode_id, intra_image_comm )
            !
            IF( .NOT. found ) &
               CALL errore( " readfile ", " occupation numbers not found! ", 1 )
            !
         END DO

         IF ( ionode ) CALL iotk_scan_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         !
      END DO k_points_loop1

      IF ( ionode ) THEN
         CALL iotk_scan_end  ( iunpun, "EIGENVALUES" )
         CALL iotk_scan_begin( iunpun, "EIGENVECTORS" )
      END IF
            !
      k_points_loop2: DO ik = 1, nk
         !
         IF ( ionode ) THEN
            CALL iotk_scan_begin( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         END IF
         !
         DO iss = 1, nspin
            IF ( ionode ) THEN
               !
               CALL iotk_scan_begin( iunpun, "WFC0" // TRIM( iotk_index (iss) ), FOUND = found )
               !
               filename = "WFC0" // TRIM( iotk_index (iss) )
               !
               IF( .NOT. found ) THEN
                  !
                  IF( nspin == 2 ) THEN
                     CALL iotk_scan_begin( iunpun, "WFC" // TRIM( iotk_index (iss) ), FOUND = found )
                     filename = "WFC" // TRIM( iotk_index (iss) )
                  ELSE
                     CALL iotk_scan_begin( iunpun, "WFC", FOUND = found )
                     filename = "WFC"
                  END IF
                  !
               END IF
               !
            END IF
            !
            CALL mp_bcast( found, ionode_id, intra_image_comm )
            !
            IF( .NOT. found ) &
               CALL errore( " readfile ", " wave functions not found! ", 1 )
            !
            IF( .NOT. ( iss > 1 .AND. force_pairing ) ) THEN
               !
               ! Only WF with spin 1 are needed when force_pairing is active
               !
               ib = iupdwn(iss)
               nb = nupdwn(iss)
               !
               ! filename is not needed we are following the link!
               !
               CALL read_wfc( iunpun, ik_eff , nk, kunit, iss_, nspin_, &
                              c02( :, ib:ib+nb-1 ), ngwt_, nbnd_, ig_l2g, ngw, &
                              filename, scalef_, .TRUE. )
               !
            END IF
            !
            IF ( ionode ) &
               CALL iotk_scan_end( iunpun, TRIM(filename) )
            !
            IF ( ionode ) THEN 
               !
               CALL iotk_scan_begin( iunpun, "WFCM" // TRIM( iotk_index (iss) ), FOUND = found )
               !
               filename = "WFCM" // TRIM( iotk_index (iss) )
               !
            END IF
            !
            CALL mp_bcast( found, ionode_id, intra_image_comm )
            !
            IF( found ) THEN
               !
               IF( .NOT. ( iss > 1 .AND. force_pairing ) ) THEN
                  !
                  ! Only WF with spin 1 are needed when force_pairing is active
                  !
                  ib = iupdwn(iss)
                  nb = nupdwn(iss)
                  !
                  CALL read_wfc( iunpun, ik_eff, nk, kunit, iss_, nspin_, &
                                 cm2( :, ib:ib+nb-1 ), ngwt_, nbnd_, ig_l2g, ngw, &
                                 filename, scalef_ , .TRUE. )
                  !
               END IF
               !
               IF ( ionode ) &
                  CALL iotk_scan_end( iunpun, TRIM( filename ) )
               !
            ELSE
               !
               cm2 = c02
               !
            END IF
            !
            IF (odd_nkscalfact .and. restart_odd_nkscalfact) THEN
               ! 
               IF ( ionode ) THEN
                  !
                  CALL iotk_scan_begin( iunpun, "WFC0FIXED" // TRIM( iotk_index (iss) ), FOUND = found )
                  !
                  filename = "WFC0FIXED" // TRIM( iotk_index (iss) )
                  !
               END IF
               !
               CALL mp_bcast( found, ionode_id, intra_image_comm )
               !
               IF( .NOT. found ) &
                  CALL errore( " readfile ", " wave functions evc0fixed not found! ", 1 )
               !
               IF( .NOT. ( iss > 1 .AND. force_pairing ) ) THEN
                  !
                  ! Only WF with spin 1 are needed when force_pairing is active
                  !
                  ib = iupdwn(iss)
                  nb = nupdwn(iss)
                  !
                  ! filename is not needed we are following the link!
                  !
                  CALL read_wfc( iunpun, ik_eff , nk, kunit, iss_, nspin_, &
                                 c0_fixed( :, ib:ib+nb-1 ), ngwt_, nbnd_, ig_l2g, ngw, &
                                 filename, scalef_, .TRUE. )
                  !
               END IF
               !
               !
               IF ( ionode ) &
                  CALL iotk_scan_end( iunpun, TRIM(filename) )
               ! 
            ENDIF
            !
         END DO
         !
         DO iss = 1, nspin
            !
            ! ... read matrix lambda to file
            !
            ALLOCATE( mrepl( nudx, nudx ) )
            !
            IF( ionode ) THEN
               CALL iotk_scan_dat( iunpun, "LAMBDA0" // TRIM( cspin ), mrepl, FOUND = found )
               IF( .NOT. found ) THEN
                  WRITE( stdout, * ) 'WARNING lambda0 not read from restart file'
                  mrepl = 0.0d0
               END IF
            END IF

            CALL mp_bcast( mrepl, ionode_id, intra_image_comm )

            CALL distribute_lambda( mrepl, lambda0(:,:,iss), descla(:,iss) )

            IF( ionode ) THEN
               CALL iotk_scan_dat( iunpun, "LAMBDAM" // TRIM( cspin ), mrepl, FOUND = found )
               IF( .NOT. found ) THEN
                  WRITE( stdout, * ) 'WARNING lambdam not read from restart file'
                  mrepl = 0.0d0
               END IF
            END IF
            ! 
            CALL mp_bcast( mrepl, ionode_id, intra_image_comm )

            CALL distribute_lambda( mrepl, lambdam(:,:,iss), descla(:,iss) )
            !
            IF ( PRESENT( mat_z ) ) THEN
               !
               IF( ionode ) THEN
                  CALL iotk_scan_dat( iunpun, "MAT_Z" // TRIM( iotk_index( iss ) ), mrepl, FOUND = found )
                  IF( .NOT. found ) THEN
                     WRITE( stdout, * ) 'WARNING mat_z not read from restart file'
                     mrepl = 0.0d0
                  END IF
               END IF

               CALL mp_bcast( mrepl, ionode_id, intra_image_comm )

               CALL distribute_zmat( mrepl, mat_z(:,:,iss), descla(:,iss) )
               !
            END IF
            !
            DEALLOCATE( mrepl )
            !
         END DO
         !
         IF ( ionode ) CALL iotk_scan_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         !
      END DO k_points_loop2
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_end( iunpun, "EIGENVECTORS" )
         !
      END IF
      !
      CALL mp_bcast( qexml_version,      ionode_id, intra_image_comm )
      CALL mp_bcast( qexml_version_init, ionode_id, intra_image_comm )
      !
      CALL mp_bcast( nfi,     ionode_id, intra_image_comm )
      CALL mp_bcast( simtime, ionode_id, intra_image_comm )
      CALL mp_bcast( title,   ionode_id, intra_image_comm )
      CALL mp_bcast( acc,     ionode_id, intra_image_comm )
      !
      CALL mp_bcast( ht,    ionode_id, intra_image_comm )
      CALL mp_bcast( htm,   ionode_id, intra_image_comm )
      CALL mp_bcast( htvel, ionode_id, intra_image_comm )
      CALL mp_bcast( gvel,  ionode_id, intra_image_comm )
      CALL mp_bcast( xnhh0, ionode_id, intra_image_comm )
      CALL mp_bcast( xnhhm, ionode_id, intra_image_comm )
      CALL mp_bcast( vnhh,  ionode_id, intra_image_comm )
      CALL mp_bcast( b1,    ionode_id, intra_image_comm )
      CALL mp_bcast( b2,    ionode_id, intra_image_comm )
      CALL mp_bcast( b3,    ionode_id, intra_image_comm )
      !
      CALL mp_bcast( stau0, ionode_id, intra_image_comm )
      CALL mp_bcast( svel0, ionode_id, intra_image_comm )
      CALL mp_bcast( staum, ionode_id, intra_image_comm )
      CALL mp_bcast( svelm, ionode_id, intra_image_comm )
      CALL mp_bcast( taui,  ionode_id, intra_image_comm )
      CALL mp_bcast( force, ionode_id, intra_image_comm )
      CALL mp_bcast( cdmi,  ionode_id, intra_image_comm )
      CALL mp_bcast( xnhp0, ionode_id, intra_image_comm )
      CALL mp_bcast( xnhpm, ionode_id, intra_image_comm ) 
      CALL mp_bcast( vnhp,  ionode_id, intra_image_comm )
      !
      CALL mp_bcast( xnhe0, ionode_id, intra_image_comm )
      CALL mp_bcast( xnhem, ionode_id, intra_image_comm )
      CALL mp_bcast( vnhe,  ionode_id, intra_image_comm )
      !
      CALL mp_bcast( kunit, ionode_id, intra_image_comm )

      CALL mp_bcast( occ0, ionode_id, intra_image_comm )
      CALL mp_bcast( occm, ionode_id, intra_image_comm )
      !
      IF ( PRESENT( mat_z ) ) &
         CALL mp_bcast( mat_z(:,:,:), ionode_id, intra_image_comm )
      !
      IF ( ionode ) &
         CALL iotk_close_read( iunpun )

      !
      s1 = cclock()
      !
      IF ( ionode ) THEN
         !
         WRITE( stdout, &
                '(3X,"restart file read in ",F8.3," sec.",/)' )  ( s1 - s0 )
         !
      END IF
      !
      if (nprint_nfi.eq.-2) then
         write( stdout,*) 'nprint_nfi= ',nprint_nfi
         CALL read_print_counter( nprint_nfi, outdir, ndr )
         write( stdout,*) 'nprint_nfi= ',nprint_nfi
      endif
      !
      RETURN
      !
    END SUBROUTINE cp_readfile_real
    !------------------------------------------------------------------------
    SUBROUTINE cp_readfile_twin( ndr, outdir, ascii, nfi, simtime, acc, nk, xk,   &
                            wk, ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh,     &
                            taui, cdmi, stau0, svel0, staum, svelm, force,    &
                            vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ0, occm,      &
                            lambda0, lambdam, b1, b2, b3, xnhe0, xnhem, vnhe, &
                            ekincm, c02, cm2, mat_z )
      !------------------------------------------------------------------------
      !
      USE control_flags,            ONLY : do_wf_cmplx, gamma_only, force_pairing !added:giovanni do_wf_cmplx
      USE io_files,                 ONLY : iunpun, xmlpun
      USE printout_base,            ONLY : title
      USE grid_dimensions,          ONLY : nr1, nr2, nr3
      USE smooth_grid_dimensions,   ONLY : nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, ONLY : nr1b, nr2b, nr3b
      USE gvecp,                    ONLY : ngm, ngmt, ecutp
      USE gvecs,                    ONLY : ngs, ngst
      USE gvecw,                    ONLY : ngw, ngwt, ecutw
      USE electrons_base,           ONLY : nspin, nbnd, nelt, nel, &
                                           nupdwn, iupdwn, nudx
      USE cell_base,                ONLY : ibrav, alat, celldm, symm_type, &
                                           s_to_r, r_to_s
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, pmass, &
                                           sort_tau, ityp, ions_cofmass
      USE reciprocal_vectors,       ONLY : ig_l2g, mill_l
      USE cp_main_variables,        ONLY : nprint_nfi, distribute_lambda, descla, distribute_zmat
      USE mp,                       ONLY : mp_sum
      USE mp_global,                ONLY : intra_image_comm
      USE parameters,               ONLY : ntypx
      USE constants,                ONLY : eps8, angstrom_au, pi
      USE twin_types
      USE input_parameters,     ONLY : odd_nkscalfact, restart_odd_nkscalfact
      USE wavefunctions_module, ONLY : c0_fixed
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN)    :: ndr          !  I/O unit number
      CHARACTER(LEN=*),      INTENT(IN)    :: outdir       !
      LOGICAL,               INTENT(IN)    :: ascii        !
      INTEGER,               INTENT(INOUT) :: nfi          ! index of the current step
      REAL(DP),              INTENT(INOUT) :: simtime      ! simulated time
      REAL(DP),              INTENT(INOUT) :: acc(:)       !
      INTEGER,               INTENT(IN)    :: nk           ! number of kpoints
      REAL(DP),              INTENT(INOUT) :: xk(:,:)      ! k-points coordinates
      REAL(DP),              INTENT(INOUT) :: wk(:)        ! k-points weights
      REAL(DP),              INTENT(INOUT) :: ht(3,3)      !
      REAL(DP),              INTENT(INOUT) :: htm(3,3)     !
      REAL(DP),              INTENT(INOUT) :: htvel(3,3)   !
      REAL(DP),              INTENT(INOUT) :: gvel(3,3)    !
      REAL(DP),              INTENT(INOUT) :: xnhh0(3,3)   !
      REAL(DP),              INTENT(INOUT) :: xnhhm(3,3)   !
      REAL(DP),              INTENT(INOUT) :: vnhh(3,3)    !
      REAL(DP),              INTENT(INOUT) :: taui(:,:)    !
      REAL(DP),              INTENT(INOUT) :: cdmi(:)      !
      REAL(DP),              INTENT(INOUT) :: stau0(:,:)   !
      REAL(DP),              INTENT(INOUT) :: svel0(:,:)   !
      REAL(DP),              INTENT(INOUT) :: staum(:,:)   !
      REAL(DP),              INTENT(INOUT) :: svelm(:,:)   !
      REAL(DP),              INTENT(INOUT) :: force(:,:)   ! 
      REAL(DP),              INTENT(INOUT) :: xnhp0(:)     !      
      REAL(DP),              INTENT(INOUT) :: xnhpm(:)     ! 
      REAL(DP),              INTENT(INOUT) :: vnhp(:)      !  
      INTEGER,               INTENT(INOUT) :: nhpcl        !  
      INTEGER,               INTENT(INOUT) :: nhpdim       !  
      REAL(DP),              INTENT(INOUT) :: occ0(:)      ! occupations
      REAL(DP),              INTENT(INOUT) :: occm(:)      !
!       REAL(DP),              INTENT(INOUT) :: lambda0(:,:,:) !
!       REAL(DP),              INTENT(INOUT) :: lambdam(:,:,:) !
      TYPE(twin_matrix), dimension(:), INTENT(INOUT) :: lambda0
      TYPE(twin_matrix), dimension(:), INTENT(INOUT) :: lambdam
      REAL(DP),              INTENT(INOUT) :: b1(3)        !
      REAL(DP),              INTENT(INOUT) :: b2(3)        !
      REAL(DP),              INTENT(INOUT) :: b3(3)        !
      REAL(DP),              INTENT(INOUT) :: xnhe0        !
      REAL(DP),              INTENT(INOUT) :: xnhem        !
      REAL(DP),              INTENT(INOUT) :: vnhe         !  
      REAL(DP),              INTENT(INOUT) :: ekincm       !  
      COMPLEX(DP),           INTENT(INOUT) :: c02(:,:)     ! 
      COMPLEX(DP),           INTENT(INOUT) :: cm2(:,:)     ! 
      TYPE(twin_matrix), dimension(:), INTENT(INOUT), optional :: mat_z
      !
      CHARACTER(LEN=256)   :: dirname, kdirname, filename
      CHARACTER(LEN=5)     :: kindex
      CHARACTER(LEN=4)     :: cspin
      INTEGER              :: strlen
      INTEGER              :: kunit
      INTEGER              :: k1, k2, k3
      INTEGER              :: nk1, nk2, nk3
      INTEGER              :: i, j, iss, ig, nspin_wfc, ierr, ik
      REAL(DP)             :: omega, htm1(3,3), hinv(3,3), scalef
      LOGICAL              :: found
      INTEGER, ALLOCATABLE :: mill(:,:)
      !
      ! ... variables read for testing pourposes
      !
      INTEGER               :: ibrav_
      CHARACTER(LEN=9)      :: symm_type_
      CHARACTER(LEN=3)      :: atm_(ntypx)
      INTEGER               :: nat_, nsp_, na_
      INTEGER               :: nk_, ik_, nt_
      LOGICAL               :: do_wf_cmplx_, gamma_only_ , lsda_ !added:giovanni do_wf_cmplx
      REAL(DP)              :: alat_, a1_(3), a2_(3), a3_(3)
      REAL(DP)              :: pmass_, zv_ 
      REAL(DP)              :: celldm_(6)
      INTEGER               :: iss_, nspin_, ngwt_, nbnd_ , n_emp_ , nbnd_tot
      INTEGER               :: nstates_up_ , nstates_dw_ , ntmp, nel_(2)
      REAL(DP)              :: nelec_ 
      REAL(DP)              :: scalef_
      REAL(DP)              :: wk_
      INTEGER               :: nhpcl_, nhpdim_ 
      INTEGER               :: ib, nb
      INTEGER               :: ik_eff
      REAL(DP)              :: amass_(ntypx)
      INTEGER,  ALLOCATABLE :: ityp_(:) 
      INTEGER,  ALLOCATABLE :: isrt_(:) 
      REAL(DP), ALLOCATABLE :: tau_(:,:) 
      REAL(DP), ALLOCATABLE :: occ_(:) 
      INTEGER,  ALLOCATABLE :: if_pos_(:,:) 
      CHARACTER(LEN=256)    :: psfile_(ntypx)
      CHARACTER(LEN=80)     :: pos_unit
      REAL(DP)              :: s1, s0, cclock
      REAL(DP), ALLOCATABLE :: mrepl(:,:) 
      COMPLEX(DP), ALLOCATABLE :: mrepl_c(:,:) 
      !
      ! ... look for an empty unit
      !
      CALL iotk_free_unit( iunout, ierr )
      !
      CALL errore( 'cp_readfile', &
                   'no free units to read wavefunctions', ierr )
      !
      kunit = 1
      found = .FALSE.
      !
      dirname = restart_dir( outdir, ndr )
      !
      ! ... Open XML descriptor
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // '/' // TRIM( xmlpun )
         !
         WRITE( stdout, '(/,3X,"reading restart file: ",A)' ) TRIM( dirname )
         !
         CALL iotk_open_read( iunpun, FILE = TRIM( filename ), IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_readfile', &
                   'cannot open restart file for reading', ierr )
      !
      s0 = cclock()
      !
      IF ( ionode ) THEN
         !
         qexml_version = " "
         !
         CALL iotk_scan_begin( iunpun, "HEADER", FOUND=found )
         !
         IF ( found ) THEN
            !
            CALL iotk_scan_empty( iunpun, "FORMAT", ATTR=attr )
            CALL iotk_scan_attr( attr, "VERSION", qexml_version )
            CALL iotk_scan_end( iunpun, "HEADER" )
            !
         ELSE
            !
            qexml_version = TRIM( default_fmt_version )
            !
         ENDIF
         !
         qexml_version_init = .TRUE.
        
         !
         ! init logical variables for versioning
         !
         qexml_version_before_1_4_0 = .FALSE.
         !
         IF ( TRIM( version_compare( qexml_version, "1.4.0" )) == "older" ) &
            qexml_version_before_1_4_0 = .TRUE.
         !
      ENDIF
      !
      CALL mp_bcast( qexml_version,               ionode_id, intra_image_comm )
      CALL mp_bcast( qexml_version_init,          ionode_id, intra_image_comm )
      CALL mp_bcast( qexml_version_before_1_4_0 , ionode_id, intra_image_comm )
      !
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "STATUS", FOUND = found )
         !
         IF ( found ) THEN
            !
            CALL iotk_scan_empty( iunpun, "STEP", attr )
            CALL iotk_scan_attr( attr, "ITERATION", nfi )
            CALL iotk_scan_dat( iunpun, "TIME", simtime )
            CALL iotk_scan_dat( iunpun, "TITLE", title )
            CALL iotk_scan_end( iunpun, "STATUS" )
            !
         END IF
         !
      END IF
      !
      ! ... Read cell and positions
      !
      ALLOCATE( tau_( 3, nat ) )
      ALLOCATE( if_pos_( 3, nat ) )
      ALLOCATE( ityp_( nat ) )
      !
      IF ( ionode ) THEN
         !
         CALL read_cell( ibrav_, symm_type_, celldm_, &
                         alat_, a1_, a2_, a3_, b1, b2, b3 )
         !
         CALL recips( a1_, a2_, a3_, b1, b2, b3 )
         !
      END IF
      !
      IF ( ionode ) THEN
         !
         CALL read_ions( nsp_, nat_, atm_, ityp_, &
                         psfile_, amass_, tau_, if_pos_, pos_unit, ierr )
         !
         IF ( ierr == 0 ) THEN
            !
            IF( nsp_ /= nsp .OR. nat_ /= nat ) ierr = 2
            !
            DO i = 1, nat
               !
               IF ( ityp_(i) /= ityp(i) ) ierr = 3
               !
            END DO
            !
         END IF
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_readfile', &
                   'cannot read positions from restart file', ierr )
      !
      !  Read SPIN infos
      !
      lsda_ = ( nspin == 2 )
      !
      IF( ionode ) THEN
         CALL iotk_scan_begin( iunpun, "SPIN", FOUND = found )
         IF( found ) THEN
            CALL iotk_scan_dat( iunpun, "LSDA", lsda_ )
            CALL iotk_scan_end( iunpun, "SPIN" )
         END IF
      END IF
      !
      CALL mp_bcast( lsda_ , ionode_id, intra_image_comm )
      !
      IF( lsda_ .AND. nspin == 1 ) &
         CALL errore( 'cp_readfile', 'LSDA restart file with a spinless run', ierr )

      !
      !  Read Occupations infos
      !
      nstates_up_ = nupdwn( 1 )
      nstates_dw_ = nupdwn( 2 )

      IF( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "OCCUPATIONS", FOUND = found )
         IF( found ) THEN
            ! 
            CALL iotk_scan_empty( iunpun, "INFO", attr, FOUND = found )
            !
            IF( lsda_ .AND. found ) THEN
               !
               IF ( qexml_version_before_1_4_0 ) THEN
                  !
                  CALL iotk_scan_attr( attr, "nelup",  nstates_up_ )
                  CALL iotk_scan_attr( attr, "neldw",  nstates_dw_ )
                  !
               ELSE
                  !
                  ! current version
                  !
                  CALL iotk_scan_attr( attr, "nstates_up",  nstates_up_ )
                  CALL iotk_scan_attr( attr, "nstates_down",  nstates_dw_ )
                  !
               ENDIF 
                !
            ENDIF
            !
            CALL iotk_scan_end( iunpun, "OCCUPATIONS" )
            !
         ENDIF
      ENDIF
      !
      CALL mp_bcast( nstates_up_ , ionode_id, intra_image_comm )
      CALL mp_bcast( nstates_dw_ , ionode_id, intra_image_comm )
      !
      IF( lsda_ ) THEN
         IF( ( nstates_up_ /=  nupdwn( 1 ) ) .OR. ( nstates_dw_ /=  nupdwn( 2 ) ) ) &
            CALL errore( 'cp_readfile', 'inconsistent number of spin states', ierr )
      END IF

      ! ... read MD timesteps variables
      !
      IF ( ionode ) &
         CALL iotk_scan_begin( iunpun, "TIMESTEPS", attr, FOUND = found )
      ! 
      ierr = 0
      ! 
      IF ( ionode .AND. found ) THEN
         !
         CALL iotk_scan_attr( attr, "nt", nt_ )
         !
         IF ( nt_ > 0 ) THEN
            !
            CALL iotk_scan_begin( iunpun, "STEP0" )
            !
            CALL iotk_scan_dat( iunpun, "ACCUMULATORS", acc )
            !
            CALL iotk_scan_begin( iunpun,"IONS_POSITIONS" )
            CALL iotk_scan_dat(   iunpun, "stau",  stau0(1:3,1:nat) )
            CALL iotk_scan_dat(   iunpun, "svel",  svel0(1:3,1:nat) )
            CALL iotk_scan_dat(   iunpun, "taui",  taui(1:3,1:nat) )
            CALL iotk_scan_dat(   iunpun, "cdmi",  cdmi(1:3) )
            CALL iotk_scan_dat(   iunpun, "force", force(1:3,1:nat) )
            CALL iotk_scan_end(   iunpun, "IONS_POSITIONS" )
            !
            CALL iotk_scan_begin( iunpun, "IONS_NOSE" )
            CALL iotk_scan_dat(   iunpun, "nhpcl", nhpcl_ )
            CALL iotk_scan_dat(   iunpun, "nhpdim", nhpdim_ )
            !
            IF ( nhpcl_ == nhpcl .AND. nhpdim_ == nhpdim ) THEN
               !
               CALL iotk_scan_dat( iunpun, "xnhp", xnhp0(1:nhpcl*nhpdim) )
               CALL iotk_scan_dat( iunpun, "vnhp", vnhp(1:nhpcl*nhpdim) )
               !
            ELSE
               !
               xnhp0(1:nhpcl*nhpdim) = 0.D0
               vnhp(1:nhpcl*nhpdim)  = 0.D0
               !
            END IF
            !
            CALL iotk_scan_end(   iunpun, "IONS_NOSE" )
            !
            CALL iotk_scan_dat( iunpun, "ekincm", ekincm )
            !
            CALL iotk_scan_begin( iunpun, "ELECTRONS_NOSE" )
            CALL iotk_scan_dat(   iunpun, "xnhe", xnhe0 )
            CALL iotk_scan_dat(   iunpun, "vnhe", vnhe )
            CALL iotk_scan_end(   iunpun, "ELECTRONS_NOSE" )
            !
            CALL iotk_scan_begin( iunpun, "CELL_PARAMETERS" )
            CALL iotk_scan_dat(   iunpun, "ht",    ht )
            CALL iotk_scan_dat(   iunpun, "htvel", htvel )
            CALL iotk_scan_dat(   iunpun, "gvel",  gvel )
            CALL iotk_scan_end(   iunpun, "CELL_PARAMETERS" )
            !
            CALL iotk_scan_begin( iunpun, "CELL_NOSE" )
            CALL iotk_scan_dat(   iunpun, "xnhh", xnhh0 )
            CALL iotk_scan_dat(   iunpun, "vnhh", vnhh )
            CALL iotk_scan_end(   iunpun, "CELL_NOSE" )
            !
            CALL iotk_scan_end( iunpun, "STEP0" )
            !
         ELSE
            !
            ierr = 40
            !
            GOTO 100
            !
         END IF
         !
         IF ( nt_ > 1 ) THEN
            !
            CALL iotk_scan_begin( iunpun, "STEPM" )
            !
            CALL iotk_scan_begin( iunpun, "IONS_POSITIONS" )
            CALL iotk_scan_dat(   iunpun, "stau", staum(1:3,1:nat) )
            CALL iotk_scan_dat(   iunpun, "svel", svelm(1:3,1:nat) )
            CALL iotk_scan_end(   iunpun, "IONS_POSITIONS" )
            !
            CALL iotk_scan_begin( iunpun, "IONS_NOSE" )
            CALL iotk_scan_dat(   iunpun, "nhpcl", nhpcl_ )
            CALL iotk_scan_dat(   iunpun, "nhpdim", nhpdim_ )
            !
            IF ( nhpcl_ == nhpcl .AND. nhpdim_ == nhpdim ) THEN
               !
               CALL iotk_scan_dat( iunpun, "xnhp",  xnhpm(1:nhpcl*nhpdim) )
               !
            ELSE
               !
               xnhpm(1:nhpcl*nhpdim) = 0.D0
               !
            END IF
            !
            CALL iotk_scan_end(   iunpun,"IONS_NOSE" )
            !
            CALL iotk_scan_begin( iunpun, "ELECTRONS_NOSE" )
            CALL iotk_scan_dat(   iunpun, "xnhe", xnhem )
            CALL iotk_scan_end(   iunpun, "ELECTRONS_NOSE" )
            !
            CALL iotk_scan_begin( iunpun, "CELL_PARAMETERS" )
            CALL iotk_scan_dat(   iunpun, "ht", htm )
            CALL iotk_scan_end(   iunpun, "CELL_PARAMETERS" )
            !
            CALL iotk_scan_begin( iunpun, "CELL_NOSE" )
            CALL iotk_scan_dat(   iunpun, "xnhh", xnhhm )
            CALL iotk_scan_end(   iunpun, "CELL_NOSE" )
            !
            CALL iotk_scan_end( iunpun, "STEPM" )
            !
         END IF
         !
         CALL iotk_scan_end( iunpun, "TIMESTEPS" )
         !
      ELSE IF ( ionode ) THEN
         !
         ! ... MD time steps not found, try to recover from CELL and POSITIONS
         ! 
         acc = 0.D0
         ! 
         ALLOCATE( isrt_( nat ) )
         !
         SELECT CASE( TRIM( pos_unit ) )
         CASE( "alat" )
            !
            tau_ = tau_ * alat_
            !
         CASE( "Angstrom" )
            !
            tau_ = tau_ * angstrom_au
            !
         CASE DEFAULT
            !
         END SELECT
         !
         CALL sort_tau( taui, isrt_ , tau_ , ityp_ , nat_ , nsp_ )
         ! 
         ht(1,:) = a1_
         ht(2,:) = a2_
         ht(3,:) = a3_
         !
         CALL invmat( 3, ht, htm1, omega )
         !
         hinv = TRANSPOSE( htm1 )
         !
         CALL r_to_s( taui, stau0, na, nsp, hinv )
         !
         CALL ions_cofmass( taui, amass_ , na, nsp, cdmi )
         !
         staum = stau0
         svel0 = 0.D0
         svelm = 0.D0
         force = 0.D0
         !
         htm   = ht
         htvel = 0.D0
         gvel  = 0.D0
         xnhh0 = 0.D0
         vnhh  = 0.D0
         xnhhm = 0.D0
         !
         xnhe0 = 0.D0
         xnhem = 0.D0
         vnhe  = 0.D0
         !
         ekincm = 0.D0
         !
         xnhp0 = 0.D0
         xnhpm = 0.D0
         vnhp  = 0.D0
         !
         DEALLOCATE( isrt_  )
         !
      END IF
      !
 100  CONTINUE
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF( ierr /= 0 ) THEN
         CALL mp_bcast( attr, ionode_id, intra_image_comm )
         CALL errore( 'cp_readfile ', TRIM( attr ), ierr )
      END IF
      !
      DEALLOCATE( tau_  )
      DEALLOCATE( if_pos_ )
      DEALLOCATE( ityp_ )
      !
      ! ... compute the scale factor
      !
      IF ( ionode ) CALL invmat( 3, ht, htm1, omega )
      !
      CALL mp_bcast( omega, ionode_id, intra_image_comm )
      !
      ! ... Beware: omega may be negative if axis are left-handed!
      !
      scalef = 1.D0 / SQRT( ABS( omega ) )
      !
      ! ... band Structure
      !
      IF ( ionode ) THEN
         !
         ierr = 0
         !
         CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE_INFO" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_SPIN_COMPONENTS", nspin_ )
         ! 
         IF ( nspin_ /= nspin ) THEN
            attr = "spin do not match"
            ierr = 31
            GOTO 90
         END IF
         !
         IF ( nspin == 2 ) THEN
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec_, ATTR = attr )
            CALL iotk_scan_attr( attr, "UP", nel_(1) )
            CALL iotk_scan_attr( attr, "DW", nel_(2) )
            !
            IF ( ( nel(1) /= nel_(1) ) .OR. ( nel(2) /= nel_(2) ) .OR. ( NINT( nelec_ ) /= nelt ) ) THEN
               attr = "electrons do not match"
               write(0,*) "from cp_readfile warning: electrons do not match"
               write(6,*) "from cp_readfile warning: electrons do not match"
               !ierr = 33
               GOTO 90
            END IF
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_BANDS", nbnd_tot , ATTR = attr )
            !
         ELSE
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec_ )
            !
            IF ( NINT( nelec_ ) /= nelt ) THEN
               attr = "electrons do not match"
               ierr = 33
               GOTO 90
            END IF
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_BANDS", nbnd_tot )
            !
         END IF
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_EMPTY_STATES", n_emp_, FOUND = found )
         !
         IF( .NOT. found ) n_emp_ = 0
         !
         nbnd_ = nbnd_tot - n_emp_
         !
         IF ( nbnd_ < nupdwn(1) ) THEN
            attr = "nbnd do not match"
            ierr = 32
            GOTO 90
         END IF
         !
         CALL iotk_scan_end( iunpun, "BAND_STRUCTURE_INFO" )
         !
      END IF
      !
 90   CONTINUE
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      IF( ierr /= 0 ) THEN
         CALL mp_bcast( attr, ionode_id, intra_image_comm )
         CALL errore( 'cp_readfile ', TRIM( attr ), ierr )
      END IF
      !
      IF( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "EIGENVALUES" )
         !
      END IF
      !
      k_points_loop1: DO ik = 1, nk
         !
         IF ( ionode ) THEN
            !
            CALL iotk_scan_begin( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
            !
            CALL iotk_scan_dat( iunpun, "WEIGHT", wk_ )
            !
         END IF
         !
         DO iss = 1, nspin
            !
            cspin = iotk_index( iss )
            !
            ik_eff = ik + ( iss - 1 ) * nk
            !
            IF ( ionode ) THEN
               !
               ALLOCATE( occ_ ( MAX( nudx , nbnd_tot ) ) )
               !
               occ_ = 0.0d0
               !
               CALL iotk_scan_dat( iunpun, "OCC0" // TRIM( cspin ), occ_ ( 1 :  nupdwn( iss ) ), FOUND = found )
               !
               IF( .NOT. found ) THEN
                  !
                  IF( nspin == 1 ) THEN
                     CALL iotk_scan_begin( iunpun, "DATAFILE", FOUND = found )
                  ELSE
                     CALL iotk_scan_begin( iunpun, "DATAFILE//TRIM(cspin)", FOUND = found )
                  END IF
                  !
                  CALL iotk_scan_dat  ( iunpun, "OCCUPATIONS", occ_( 1:nbnd_tot ) ) 
                  !
                  IF( nspin == 1 ) THEN
                     CALL iotk_scan_end( iunpun, "DATAFILE" )
                  ELSE
                     CALL iotk_scan_end( iunpun, "DATAFILE//TRIM(cspin)" )
                  END IF
                  !
                  IF( found ) THEN
                     occ0( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) = occ_ ( 1:nupdwn( iss ) ) * wk_
                     occm( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) = occ_ ( 1:nupdwn( iss ) ) * wk_
                  END IF
                  !
               ELSE
                  !
                  occ0( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) = occ_ ( 1:nupdwn( iss ) )
                  !
                  CALL iotk_scan_dat( iunpun, "OCCM" // TRIM( cspin ), occ_ ( 1 :  nupdwn( iss ) ), FOUND = found )
                  !
                  IF( found ) THEN
                     occm( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) = occ_ ( 1:nupdwn( iss ) )
                  END IF
                  !
               END IF
               !
               DEALLOCATE ( occ_ )
               !
            END IF
            !
            CALL mp_bcast( found, ionode_id, intra_image_comm )
            !
            IF( .NOT. found ) &
               CALL errore( " readfile ", " occupation numbers not found! ", 1 )
            !
         END DO

         IF ( ionode ) CALL iotk_scan_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         !
      END DO k_points_loop1

      IF ( ionode ) THEN
         CALL iotk_scan_end  ( iunpun, "EIGENVALUES" )
         CALL iotk_scan_begin( iunpun, "EIGENVECTORS" )
      END IF
            !
      k_points_loop2: DO ik = 1, nk
         !
         IF ( ionode ) THEN
            CALL iotk_scan_begin( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         END IF
         !
         DO iss = 1, nspin
            IF ( ionode ) THEN
               !
               CALL iotk_scan_begin( iunpun, "WFC0" // TRIM( iotk_index (iss) ), FOUND = found )
               !
               filename = "WFC0" // TRIM( iotk_index (iss) )
               !
               IF( .NOT. found ) THEN
                  !
                  IF( nspin == 2 ) THEN
                     CALL iotk_scan_begin( iunpun, "WFC" // TRIM( iotk_index (iss) ), FOUND = found )
                     filename = "WFC" // TRIM( iotk_index (iss) )
                  ELSE
                     CALL iotk_scan_begin( iunpun, "WFC", FOUND = found )
                     filename = "WFC"
                  END IF
                  !
               END IF
               !
            END IF
            !
            CALL mp_bcast( found, ionode_id, intra_image_comm )
            !
            IF( .NOT. found ) &
               CALL errore( " readfile ", " wave functions not found! ", 1 )
            !
            IF( .NOT. ( iss > 1 .AND. force_pairing ) ) THEN
               !
               ! Only WF with spin 1 are needed when force_pairing is active
               !
               ib = iupdwn(iss)
               nb = nupdwn(iss)
               !
               ! filename is not needed we are following the link!
               !
               CALL read_wfc( iunpun, ik_eff , nk, kunit, iss_, nspin_, &
                              c02( :, ib:ib+nb-1 ), ngwt_, nbnd_, ig_l2g, ngw, &
                              filename, scalef_, .TRUE. )
               !
            END IF
            !
            IF ( ionode ) &
               CALL iotk_scan_end( iunpun, TRIM(filename) )
            !
            IF ( ionode ) THEN 
               !
               CALL iotk_scan_begin( iunpun, "WFCM" // TRIM( iotk_index (iss) ), FOUND = found )
               !
               filename = "WFCM" // TRIM( iotk_index (iss) )
               !
            END IF
            !
            CALL mp_bcast( found, ionode_id, intra_image_comm )
            !
            IF( found ) THEN
               !
               IF( .NOT. ( iss > 1 .AND. force_pairing ) ) THEN
                  !
                  ! Only WF with spin 1 are needed when force_pairing is active
                  !
                  ib = iupdwn(iss)
                  nb = nupdwn(iss)
                  !
                  CALL read_wfc( iunpun, ik_eff, nk, kunit, iss_, nspin_, &
                                 cm2( :, ib:ib+nb-1 ), ngwt_, nbnd_, ig_l2g, ngw, &
                                 filename, scalef_ , .TRUE. )
                  !
               END IF
               !
               IF ( ionode ) &
                  CALL iotk_scan_end( iunpun, TRIM( filename ) )
               !
            ELSE
               !
               cm2 = c02
               !
            END IF
            !
            IF (odd_nkscalfact .and. restart_odd_nkscalfact) THEN
               !
               IF ( ionode ) THEN
                  !
                  CALL iotk_scan_begin( iunpun, "WFC0FIXED" // TRIM( iotk_index (iss) ), FOUND = found )
                  !
                  filename = "WFC0FIXED" // TRIM( iotk_index (iss) )
                  !
               END IF
               !
               CALL mp_bcast( found, ionode_id, intra_image_comm )
               !
               IF( .NOT. found ) &
                  CALL errore( " readfile ", " wave functions evc0fixed not found! ", 1 )
               !
               IF( .NOT. ( iss > 1 .AND. force_pairing ) ) THEN
                  !
                  ! Only WF with spin 1 are needed when force_pairing is active
                  !
                  ib = iupdwn(iss)
                  nb = nupdwn(iss)
                  !
                  ! filename is not needed we are following the link!
                  !
                  CALL read_wfc( iunpun, ik_eff , nk, kunit, iss_, nspin_, &
                                 c0_fixed( :, ib:ib+nb-1 ), ngwt_, nbnd_, ig_l2g, ngw, &
                                 filename, scalef_, .TRUE. )
                  !
               END IF
               !
               IF ( ionode ) &
                  CALL iotk_scan_end( iunpun, TRIM(filename) )
               !
            ENDIF
            !
         END DO
         !
         DO iss = 1, nspin
            !
            ! ... read matrix lambda to file
            !
            IF(.not. lambda0(1)%iscmplx) THEN
                ALLOCATE(mrepl(nudx, nudx ))            
		!
		IF( ionode ) THEN
		  CALL iotk_scan_dat( iunpun, "LAMBDA0" // TRIM( cspin ), mrepl, FOUND = found )
		  IF( .NOT. found ) THEN
		      WRITE( stdout, * ) 'WARNING lambda0 not read from restart file'
		      mrepl = 0.0d0
		  END IF
		END IF

		CALL mp_bcast( mrepl, ionode_id, intra_image_comm )
		CALL distribute_lambda( mrepl, lambda0(iss)%rvec(:,:), descla(:,iss) )

		IF( ionode ) THEN
		  CALL iotk_scan_dat( iunpun, "LAMBDAM" // TRIM( cspin ), mrepl, FOUND = found )
		  IF( .NOT. found ) THEN
		      WRITE( stdout, * ) 'WARNING lambdam not read from restart file'
		      mrepl = 0.0d0
		  END IF
		END IF
		! 
		CALL mp_bcast( mrepl, ionode_id, intra_image_comm )

		CALL distribute_lambda( mrepl, lambdam(iss)%rvec(:,:), descla(:,iss) )
               DEALLOCATE(mrepl)
		!
            ELSE
                  WRITE( stdout, * ) 'here should be iotk first error'
                ALLOCATE(mrepl_c(nudx, nudx))
		IF( ionode ) THEN
		  CALL iotk_scan_dat( iunpun, "LAMBDA0" // TRIM( cspin ), mrepl_c, FOUND = found )
		  IF( .NOT. found ) THEN
		      WRITE( stdout, * ) 'WARNING lambda0 not read from restart file'
		      mrepl_c = CMPLX(0.0d0, 0.d0)
		  END IF
		END IF

		CALL mp_bcast( mrepl_c, ionode_id, intra_image_comm )

		CALL distribute_lambda(mrepl_c, lambda0(iss)%cvec(:,:), descla(:,iss) )

		IF( ionode ) THEN
		  CALL iotk_scan_dat( iunpun, "LAMBDAM" // TRIM( cspin ), mrepl_c, FOUND = found )
		  IF( .NOT. found ) THEN
		      WRITE( stdout, * ) 'WARNING lambdam not read from restart file'
		      mrepl_c = CMPLX(0.0d0,0.d0)
		  END IF
		END IF
		! 
		CALL mp_bcast( mrepl_c, ionode_id, intra_image_comm )
		CALL distribute_lambda( mrepl_c, lambdam(iss)%cvec(:,:), descla(:,iss) )
               DEALLOCATE(mrepl_c)

            ENDIF

            IF ( PRESENT( mat_z ) ) THEN
               !

               IF( ionode ) THEN
                  IF(.not.mat_z(iss)%iscmplx) THEN
		    IF(.not.allocated(mrepl)) THEN
		      ALLOCATE(mrepl(nudx,nudx))
		    ENDIF
		    CALL iotk_scan_dat( iunpun, "MAT_Z" // TRIM( iotk_index( iss ) ), mrepl, FOUND = found )
		    IF( .NOT. found ) THEN
		      WRITE( stdout, * ) 'WARNING mat_z not read from restart file'
		      mrepl = 0.0d0
		    END IF
		    CALL mp_bcast( mrepl, ionode_id, intra_image_comm )
		    CALL distribute_zmat( mrepl, mat_z(iss)%rvec(:,:), descla(:,iss) )
		    DEALLOCATE(mrepl)
                  ELSE
                   IF(.not.allocated(mrepl_c)) THEN
                    ALLOCATE(mrepl_c(nudx,nudx))
                   ENDIF
		    CALL iotk_scan_dat( iunpun, "MAT_Z" // TRIM( iotk_index( iss ) ), mrepl_c, FOUND = found )
		    IF( .NOT. found ) THEN
		      WRITE( stdout, * ) 'WARNING mat_z not read from restart file'
		      mrepl_c = CMPLX(0.0d0,0.d0)
		    END IF
		    CALL mp_bcast( mrepl_c, ionode_id, intra_image_comm )
		    CALL distribute_zmat( mrepl_c, mat_z(iss)%cvec(:,:), descla(:,iss) )
		    DEALLOCATE(mrepl_c)
                  ENDIF
               END IF
               !
            END IF
            !
            !
         END DO
         !
         IF ( ionode ) CALL iotk_scan_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         !
      END DO k_points_loop2
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_end( iunpun, "EIGENVECTORS" )
         !
      END IF
      !
      CALL mp_bcast( qexml_version,      ionode_id, intra_image_comm )
      CALL mp_bcast( qexml_version_init, ionode_id, intra_image_comm )
      !
      CALL mp_bcast( nfi,     ionode_id, intra_image_comm )
      CALL mp_bcast( simtime, ionode_id, intra_image_comm )
      CALL mp_bcast( title,   ionode_id, intra_image_comm )
      CALL mp_bcast( acc,     ionode_id, intra_image_comm )
      !
      CALL mp_bcast( ht,    ionode_id, intra_image_comm )
      CALL mp_bcast( htm,   ionode_id, intra_image_comm )
      CALL mp_bcast( htvel, ionode_id, intra_image_comm )
      CALL mp_bcast( gvel,  ionode_id, intra_image_comm )
      CALL mp_bcast( xnhh0, ionode_id, intra_image_comm )
      CALL mp_bcast( xnhhm, ionode_id, intra_image_comm )
      CALL mp_bcast( vnhh,  ionode_id, intra_image_comm )
      CALL mp_bcast( b1,    ionode_id, intra_image_comm )
      CALL mp_bcast( b2,    ionode_id, intra_image_comm )
      CALL mp_bcast( b3,    ionode_id, intra_image_comm )
      !
      CALL mp_bcast( stau0, ionode_id, intra_image_comm )
      CALL mp_bcast( svel0, ionode_id, intra_image_comm )
      CALL mp_bcast( staum, ionode_id, intra_image_comm )
      CALL mp_bcast( svelm, ionode_id, intra_image_comm )
      CALL mp_bcast( taui,  ionode_id, intra_image_comm )
      CALL mp_bcast( force, ionode_id, intra_image_comm )
      CALL mp_bcast( cdmi,  ionode_id, intra_image_comm )
      CALL mp_bcast( xnhp0, ionode_id, intra_image_comm )
      CALL mp_bcast( xnhpm, ionode_id, intra_image_comm ) 
      CALL mp_bcast( vnhp,  ionode_id, intra_image_comm )
      !
      CALL mp_bcast( xnhe0, ionode_id, intra_image_comm )
      CALL mp_bcast( xnhem, ionode_id, intra_image_comm )
      CALL mp_bcast( vnhe,  ionode_id, intra_image_comm )
      !
      CALL mp_bcast( kunit, ionode_id, intra_image_comm )

      CALL mp_bcast( occ0, ionode_id, intra_image_comm )
      CALL mp_bcast( occm, ionode_id, intra_image_comm )
      !
!       IF ( PRESENT( mat_z ) ) THEN !warning:giovanni this part is a bug???
! 	DO iss=1,nspin
! 	  IF(.not.mat_z(iss)%iscmplx) THEN  
! 	    CALL mp_bcast( mat_z(iss)%rvec(:,:), ionode_id, intra_image_comm )
!           ELSE
! 	    CALL mp_bcast( mat_z(iss)%cvec(:,:), ionode_id, intra_image_comm )
!           ENDIF
! 	ENDDO
!       ENDIF
      !
      IF ( ionode ) &
         CALL iotk_close_read( iunpun )

      !
      s1 = cclock()
      !
      IF ( ionode ) THEN
         !
         WRITE( stdout, &
                '(3X,"restart file read in ",F8.3," sec.",/)' )  ( s1 - s0 )
         !
      END IF
      !
      if (nprint_nfi.eq.-2) then
         write( stdout,*) 'nprint_nfi= ',nprint_nfi
         CALL read_print_counter( nprint_nfi, outdir, ndr )
         write( stdout,*) 'nprint_nfi= ',nprint_nfi
      endif
      !
      RETURN
      !
    END SUBROUTINE cp_readfile_twin
    ! 
    !------------------------------------------------------------------------
    SUBROUTINE cp_read_wfc( ndr, outdir, ik, nk, iss, nspin, c2, tag )
      !------------------------------------------------------------------------
      !
      USE electrons_base,     ONLY : iupdwn, nupdwn
      USE reciprocal_vectors, ONLY : ngwt, ngw, ig_l2g
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN)  :: ndr
      CHARACTER(LEN=*),      INTENT(IN)  :: outdir
      INTEGER,               INTENT(IN)  :: ik, iss, nk, nspin
      CHARACTER,             INTENT(IN)  :: tag
      COMPLEX(DP), OPTIONAL, INTENT(OUT) :: c2(:,:)
      !
      CHARACTER(LEN=256) :: dirname, filename
      INTEGER            :: ik_eff, ib, nb, kunit, iss_, nspin_, ngwt_, nbnd_
      REAL(DP)           :: scalef
      !
      kunit = 1
      !
      ik_eff = ik + ( iss - 1 ) * nk
      !
      dirname = restart_dir( outdir, ndr )
      !
      IF ( tag /= 'm' ) THEN
         !
         IF ( nspin == 1 ) THEN
            !
            filename = TRIM( wfc_filename( dirname, 'evc0', ik ) )
            !
         ELSE
            !
            filename = TRIM( wfc_filename( dirname, 'evc0', ik, iss ) )
            !
         END IF
         !
      ELSE
         !
         IF ( nspin == 1 ) THEN
            !
            filename = TRIM( wfc_filename( dirname, 'evcm', ik ) )
            !
         ELSE
            !
            filename = TRIM( wfc_filename( dirname, 'evcm', ik, iss ) )
            !
         END IF
         !
      END IF
      !
      ib = iupdwn(iss)
      nb = nupdwn(iss)
      !
      CALL read_wfc( iunout, ik_eff, nk, kunit, iss_, nspin_, &
                     c2(:,ib:ib+nb-1), ngwt_, nbnd_, ig_l2g, ngw,  &
                     filename, scalef )
      !
      RETURN
      !
    END SUBROUTINE cp_read_wfc
    !
    !------------------------------------------------------------------------
    SUBROUTINE cp_read_cell( ndr, outdir, ascii, ht, &
                             htm, htvel, gvel, xnhh0, xnhhm, vnhh )
      !------------------------------------------------------------------------
      !
      USE io_files,  ONLY : iunpun, xmlpun
      USE mp_global, ONLY : intra_image_comm
      USE mp,        ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN)    :: ndr
      CHARACTER(LEN=*), INTENT(IN)    :: outdir
      LOGICAL,          INTENT(IN)    :: ascii
      REAL(DP),         INTENT(INOUT) :: ht(3,3)
      REAL(DP),         INTENT(INOUT) :: htm(3,3)
      REAL(DP),         INTENT(INOUT) :: htvel(3,3)
      REAL(DP),         INTENT(INOUT) :: gvel(3,3)
      REAL(DP),         INTENT(INOUT) :: xnhh0(3,3)
      REAL(DP),         INTENT(INOUT) :: xnhhm(3,3)
      REAL(DP),         INTENT(INOUT) :: vnhh(3,3)
      !
      CHARACTER(LEN=256) :: dirname, filename
      INTEGER            :: strlen
      INTEGER            :: i, ierr, nt_
      LOGICAL            :: found
      !
      ! ... variables read for testing pourposes
      !
      INTEGER          :: ibrav_
      REAL(DP)         :: alat_
      REAL(DP)         :: celldm_(6)
      REAL(DP)         :: a1_(3), a2_(3), a3_(3)
      REAL(DP)         :: b1_(3), b2_(3), b3_(3)
      CHARACTER(LEN=9) :: symm_type_
      !
      !
      dirname = restart_dir( outdir, ndr ) 
      !
      filename = TRIM( dirname ) // '/' // TRIM( xmlpun )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( filename ), &
                              BINARY = .FALSE., ROOT = attr, IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_read_cell', &
                   'cannot open restart file for reading: ' // TRIM(filename), &
                   ierr )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "TIMESTEPS", attr, FOUND = found )
         !
         IF ( found ) THEN
            !
            CALL iotk_scan_attr( attr, "nt", nt_ )
            !
            IF ( nt_ > 0 ) THEN
               !
               CALL iotk_scan_begin( iunpun, "STEP0" )
               !
               CALL iotk_scan_begin( iunpun, "CELL_PARAMETERS" )
               CALL iotk_scan_dat(   iunpun, "ht",    ht )
               CALL iotk_scan_dat(   iunpun, "htvel", htvel )
               CALL iotk_scan_dat(   iunpun, "gvel",  gvel, &
                                     FOUND = found, IERR = ierr )
               !
               IF ( .NOT. found ) gvel = 0.D0
               !
               CALL iotk_scan_end( iunpun, "CELL_PARAMETERS" )
               !
               CALL iotk_scan_begin( iunpun, "CELL_NOSE" )
               CALL iotk_scan_dat(   iunpun, "xnhh", xnhh0 )
               CALL iotk_scan_dat(   iunpun, "vnhh", vnhh )
               CALL iotk_scan_end(   iunpun, "CELL_NOSE" )
               !
               CALL iotk_scan_end( iunpun, "STEP0" )
               !
            ELSE
               !
               ierr = 40
               !
               GOTO 100
               !
            END IF
            !
            IF( nt_ > 1 ) THEN
               !
               CALL iotk_scan_begin(iunpun,"STEPM")
               !
               CALL iotk_scan_begin( iunpun, "CELL_PARAMETERS" )
               CALL iotk_scan_dat(   iunpun, "ht", htm)
               CALL iotk_scan_end(   iunpun, "CELL_PARAMETERS" )
               !
               CALL iotk_scan_begin( iunpun, "CELL_NOSE" )
               CALL iotk_scan_dat(   iunpun, "xnhh", xnhhm )
               CALL iotk_scan_end(   iunpun, "CELL_NOSE" )
               !
               CALL iotk_scan_end( iunpun, "STEPM" )
               !
            END IF
            !
            CALL iotk_scan_end( iunpun, "TIMESTEPS" )
            !
         ELSE
            !
            ! ... MD steps have not been found, try to restart from cell data
            !
            CALL read_cell( ibrav_, symm_type_, celldm_, &
                            alat_, a1_, a2_, a3_, b1_, b2_, b3_ )
            !
            ht(1,:) = a1_
            ht(2,:) = a2_
            ht(3,:) = a3_
            !
            htm   = ht
            htvel = 0.D0
            gvel  = 0.D0
            xnhh0 = 0.D0
            vnhh  = 0.D0
            xnhhm = 0.D0
            !
         END IF
         !
      END IF
      !
 100  CONTINUE
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      CALL mp_bcast( attr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_read_cell ', attr, ierr )
      !
      CALL mp_bcast( ht,    ionode_id, intra_image_comm )
      CALL mp_bcast( htm,   ionode_id, intra_image_comm )
      CALL mp_bcast( htvel, ionode_id, intra_image_comm )
      CALL mp_bcast( gvel,  ionode_id, intra_image_comm )
      CALL mp_bcast( xnhh0, ionode_id, intra_image_comm )
      CALL mp_bcast( xnhhm, ionode_id, intra_image_comm )
      CALL mp_bcast( vnhh,  ionode_id, intra_image_comm )
      !
      IF ( ionode ) CALL iotk_close_read( iunpun )
      !
      RETURN
      !
    END SUBROUTINE cp_read_cell
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_cell( ibrav, symm_type, &
                          celldm, alat, a1, a2, a3, b1, b2, b3 )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(OUT) :: ibrav
      CHARACTER(LEN=*), INTENT(OUT) :: symm_type
      REAL(DP),         INTENT(OUT) :: celldm(6), alat
      REAL(DP),         INTENT(OUT) :: a1(3), a2(3), a3(3)
      REAL(DP),         INTENT(OUT) :: b1(3), b2(3), b3(3)
      !
      CHARACTER(LEN=256) :: bravais_lattice
      !
      !
      CALL iotk_scan_begin( iunpun, "CELL" )
      !
      CALL iotk_scan_dat( iunpun, "BRAVAIS_LATTICE", bravais_lattice )
      !
      SELECT CASE ( TRIM( bravais_lattice ) )
        CASE( "free" )
           ibrav = 0
        CASE( "cubic P (sc)" )
           ibrav = 1
        CASE( "cubic F (fcc)" )
           ibrav = 2
        CASE( "cubic I (bcc)" )
           ibrav = 3
        CASE( "Hexagonal and Trigonal P" )
           ibrav = 4
        CASE( "Trigonal R" )
           ibrav = 5
        CASE( "Tetragonal P (st)" )
           ibrav = 6
        CASE( "Tetragonal I (bct)" )
           ibrav = 7
        CASE( "Orthorhombic P" )
           ibrav = 8
        CASE( "Orthorhombic base-centered(bco)" )
           ibrav = 9
        CASE( "Orthorhombic face-centered" )
           ibrav = 10
        CASE( "Orthorhombic body-centered" )
           ibrav = 11
        CASE( "Monoclinic P" )
           ibrav = 12
        CASE( "Monoclinic base-centered" )
           ibrav = 13
        CASE( "Triclinic P" )
           ibrav = 14
      END SELECT
      !
      IF ( ibrav == 0 ) &
         CALL iotk_scan_dat( iunpun, "CELL_SYMMETRY", symm_type )
      !
      CALL iotk_scan_dat( iunpun, "LATTICE_PARAMETER", alat )
      CALL iotk_scan_dat( iunpun, "CELL_DIMENSIONS", celldm(1:6) )
      !
      CALL iotk_scan_begin( iunpun, "DIRECT_LATTICE_VECTORS" )
      CALL iotk_scan_dat(   iunpun, "a1", a1 )
      CALL iotk_scan_dat(   iunpun, "a2", a2 )
      CALL iotk_scan_dat(   iunpun, "a3", a3 )
      CALL iotk_scan_end(   iunpun, "DIRECT_LATTICE_VECTORS" )
      !
      CALL iotk_scan_begin( iunpun, "RECIPROCAL_LATTICE_VECTORS" )
      CALL iotk_scan_dat(   iunpun, "b1", b1 )
      CALL iotk_scan_dat(   iunpun, "b2", b2 )
      CALL iotk_scan_dat(   iunpun, "b3", b3 )
      CALL iotk_scan_end(   iunpun, "RECIPROCAL_LATTICE_VECTORS" )
      !
      CALL iotk_scan_end( iunpun, "CELL" )
      !
      RETURN
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_ions( nsp, nat, atm, ityp, psfile, &
                          amass, tau, if_pos, pos_unit, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,            INTENT(OUT) :: nsp, nat
      CHARACTER(LEN=3),   INTENT(OUT) :: atm(:)
      INTEGER,            INTENT(OUT) :: ityp(:)
      CHARACTER(LEN=256), INTENT(OUT) :: psfile(:)
      REAL(DP),           INTENT(OUT) :: amass(:)
      REAL(DP),           INTENT(OUT) :: tau(:,:)
      INTEGER,            INTENT(OUT) :: if_pos(:,:)
      INTEGER,            INTENT(OUT) :: ierr
      CHARACTER(LEN=*),   INTENT(OUT) :: pos_unit
      !
      LOGICAL          :: found, back_compat
      INTEGER          :: i
      CHARACTER(LEN=3) :: lab
      !
      ierr = 0
      !
      CALL iotk_scan_begin( iunpun, "IONS", FOUND = found )
      !
      IF ( .NOT. found ) THEN
         !
         ierr = 1
         !
         RETURN
         !
      END IF
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_ATOMS",   nat )
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_SPECIES", nsp )
      !
      IF ( nsp > SIZE( atm ) .OR. nat > SIZE( ityp ) ) THEN
         !
         ierr = 10
         !
         CALL iotk_scan_end( iunpun, "IONS" )
         !
         RETURN
         !
      END IF
      !
      !
      DO i = 1, nsp
         !
         IF ( qexml_version_before_1_4_0 ) THEN
            !
            CALL iotk_scan_dat( iunpun, "ATOM_TYPE", atm(i) )
            CALL iotk_scan_dat( iunpun, TRIM( atm(i) )//"_MASS", amass(i) )
            CALL iotk_scan_dat( iunpun, "PSEUDO_FOR_" // TRIM( atm(i) ), psfile(i) )
            !
         ELSE
            !
            ! current format
            !
            CALL iotk_scan_begin( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
            !
            CALL iotk_scan_dat( iunpun, "ATOM_TYPE", atm(i) )
            CALL iotk_scan_dat( iunpun, "MASS", amass(i) )
            CALL iotk_scan_dat( iunpun, "PSEUDO", psfile(i) )
            !
            CALL iotk_scan_end( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
            !
         ENDIF
         !
      ENDDO
      !
      CALL iotk_scan_empty( iunpun, "UNITS_FOR_ATOMIC_POSITIONS", attr )
      CALL iotk_scan_attr( attr, "UNITS", pos_unit  )
      !
      DO i = 1, nat
         !
         CALL iotk_scan_empty( iunpun, "ATOM" // TRIM( iotk_index( i ) ), attr )
         CALL iotk_scan_attr( attr, "SPECIES", lab )
         CALL iotk_scan_attr( attr, "INDEX",   ityp(i) )
         CALL iotk_scan_attr( attr, "tau",     tau(:,i) )
         CALL iotk_scan_attr( attr, "if_pos",  if_pos(:,i) )
         !
      END DO
      !
      CALL iotk_scan_end( iunpun, "IONS" )
      !
      RETURN
      !
    END SUBROUTINE read_ions
    !
    !
    !
    SUBROUTINE write_gk( iun, ik, mill, filename )
       !
       USE gvecw,                    ONLY : ngw, ngwt
       USE control_flags,            ONLY : do_wf_cmplx, gamma_only !added:giovanni do_wf_cmplx
       USE reciprocal_vectors,       ONLY : ig_l2g, mill_l
       USE mp,                       ONLY : mp_sum
       USE mp_global,                ONLY : intra_image_comm
       USE io_global,                ONLY : ionode
       !
       IMPLICIT NONE
       !
       INTEGER,            INTENT(IN) :: iun, ik
       INTEGER,            INTENT(IN) :: mill(:,:)
       CHARACTER(LEN=256), INTENT(IN) :: filename
       !
       INTEGER, ALLOCATABLE :: igwk(:)
       INTEGER, ALLOCATABLE :: itmp1(:)
       INTEGER  :: npwx_g, npw_g, ig, ngg
       REAL(DP) :: xk(3)

       xk     = 0.0d0
       npwx_g = ngwt
       npw_g  = ngwt

       ALLOCATE( igwk( npwx_g ) )
       ! 
       igwk = 0
       !
       ALLOCATE( itmp1( npw_g ) )
       !
       itmp1 = 0
       !
       !
       DO ig = 1, ngw
          !
          itmp1( ig_l2g( ig ) ) = ig_l2g( ig )
          !
       END DO
       !
       CALL mp_sum( itmp1, intra_image_comm )
       !
       ngg = 0
       !
       DO ig = 1, npw_g
          !
          IF ( itmp1(ig) == ig ) THEN
             !
             ngg = ngg + 1
             !
             igwk( ngg ) = ig
             !
          END IF
          !
       END DO

       DEALLOCATE( itmp1 )
       !
       IF ( ionode ) THEN
          !
          CALL iotk_open_write( iun, FILE = TRIM( filename ), &
                                ROOT="GK-VECTORS", BINARY = .TRUE. )
          !
          CALL iotk_write_dat( iun, "NUMBER_OF_GK-VECTORS", npw_g )
          CALL iotk_write_dat( iun, "MAX_NUMBER_OF_GK-VECTORS", npwx_g )
          CALL iotk_write_dat( iun, "DO_WF_CMPLX", do_wf_cmplx ) !added:giovanni do_wf_cmplx
          CALL iotk_write_dat( iun, "GAMMA_ONLY", gamma_only.and..not.do_wf_cmplx )!modified:giovanni for post-processing
          !
          CALL iotk_write_attr ( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
          CALL iotk_write_dat( iun, "K-POINT_COORDS", xk(:), ATTR = attr )
          !
          CALL iotk_write_dat( iun, "INDEX", igwk( 1:npw_g ) )
          CALL iotk_write_dat( iun, "GRID", mill( 1:3, igwk( 1:npw_g ) ), COLUMNS = 3 )
          !
          CALL iotk_close_write( iun )
          !
       END IF
       !
       DEALLOCATE( igwk )

       RETURN

    END SUBROUTINE write_gk
    !
    !
    !
!     SUBROUTINE write_translation(trans_matrix, trans_vec, )!added:giovanni
! 
! 
!                CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin,        &
!                                ctot( :, ib : ib + nbnd_tot - 1 ), ngwt, do_wf_cmplx, gamma_only,& !added:giovanni do_wf_cmplx
!                                nbnd_tot, ig_l2g, ngw, filename, scalef )
! 	      ib = iupdwn(iss)
! 	      nb = nupdwn(iss)
!       !
! 	      CALL read_wfc( iunout, ik_eff, nk, kunit, iss_, nspin_, &
!                      c2(:,ib:ib+nb-1), ngwt_, nbnd_, ig_l2g, ngw,  &
!                      filename, scalef )
!       ! 
!     END SUBROUTINE write_translation

END MODULE cp_restart
