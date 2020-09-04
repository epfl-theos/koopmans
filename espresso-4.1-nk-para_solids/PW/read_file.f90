!
! Copyright (C) 2001-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE read_file()
  !----------------------------------------------------------------------------
  !
  ! ... This routine allocates space for all quantities already computed
  ! ... in the pwscf program and reads them from the data file.
  ! ... All quantities that are initialized in subroutine "setup" when
  ! ... starting from scratch should be initialized here when restarting
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, nsp, ityp, tau, if_pos, extfor
  USE basis,                ONLY : natomwfc
  USE cell_base,            ONLY : tpiba2, alat,omega, at, bg, ibrav
  USE force_mod,            ONLY : force
  USE klist,                ONLY : nkstot, nks, xk, wk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, nbndx, et, wg, npwx
  USE symme,                ONLY : irt, nsym, ftau, s, d1, d2, d3
  USE ktetra,               ONLY : tetra, ntetra 
  USE extfield,             ONLY : forcefield, tefield
  USE cellmd,               ONLY : cell_factor, lmovecell
  USE gvect,                ONLY : gg, ecutwfc, ngm, g, nr1, nr2, nr3, nrxx,&
                                   nrx1, nrx2, nrx3, eigts1, eigts2, eigts3, &
                                   nl, gstart
  USE gsmooth,              ONLY : ngms, nls, nrx1s, nr1s, nr2s, nr3s
  USE spin_orb,             ONLY : lspinorb, domag
  USE scf,                  ONLY : rho, rho_core, rhog_core, v
  USE wavefunctions_module, ONLY : psic
  USE vlocal,               ONLY : strf
  USE io_files,             ONLY : tmp_dir, prefix, iunpun, nwordwfc, iunwfc, qexml_version
  USE buffers,              ONLY : open_buffer, close_buffer
  USE uspp_param,           ONLY : upf
  USE noncollin_module,     ONLY : noncolin, npol, nspin_lsda, nspin_mag, nspin_gga
  USE pw_restart,           ONLY : pw_readfile
  USE xml_io_base,          ONLY : pp_check_file
  USE uspp,                 ONLY : okvan, becsum
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_onecenter,        ONLY : paw_potential
  USE paw_init,             ONLY : paw_init_onecenter, allocate_paw_internals
  USE ldaU,                 ONLY : eth
  USE realus,               ONLY : qpointlist,betapointlist,init_realspace_vars,real_space
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER  :: i, is, ik, ibnd, nb, nt, ios, isym, ierr
  REAL(DP) :: rdum(1,1), ehart, etxc, vtxc, etotefield, charge
  REAL(DP) :: sr(3,3,48)
  LOGICAL  :: exst
  !
  !
  ! ... first we get the version of the qexml file
  !     if not already read
  !
  CALL pw_readfile( 'header', ierr )
  CALL errore( 'read_file ', 'unable to determine qexml version', ABS(ierr) )
  !
  ! ... then we check if the file can be used for post-processing
  !
  IF ( .NOT. pp_check_file() ) &
     CALL infomsg( 'read_file', 'file ' // TRIM( tmp_dir ) // TRIM( prefix ) &
               & // '.save not guaranteed to be safe for post-processing' )
  !
  ! ... here we read the variables that dimension the system
  ! ... in parallel execution, only root proc reads the file
  ! ... and then broadcasts the values to all other procs
  !
  ! ... a reset of the internal flags is necessary because some codes call
  ! ... read_file() more than once
  !
  CALL pw_readfile( 'reset', ierr )
  CALL pw_readfile( 'dim',   ierr )
  !
  CALL errore( 'read_file ', 'problem reading file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.save', ierr )
  !
  ! ... allocate space for atomic positions, symmetries, forces, tetrahedra
  !
  IF ( nat < 0 ) &
     CALL errore( 'read_file', 'wrong number of atoms', 1 )
  !
  ! ... allocation
  !
  ALLOCATE( ityp( nat ) )
  !
  ALLOCATE( tau(    3, nat ) )
  ALLOCATE( if_pos( 3, nat ) )
  ALLOCATE( force(  3, nat ) )
  ALLOCATE( extfor(  3, nat ) )
  !
  IF ( tefield ) ALLOCATE( forcefield( 3, nat ) )
  !
  ALLOCATE( irt( 48, nat ) )
  ALLOCATE( tetra( 4, MAX( ntetra, 1 ) ) )
  !
  ! ... here we read all the variables defining the system
  ! ... in parallel execution, only root proc read the file
  ! ... and then broadcast the values to all other procs
  !
  !-------------------------------------------------------------------------------
  ! ... XML punch-file
  !-------------------------------------------------------------------------------
  !
  CALL set_dimensions()
  !
  ! ... check whether LSDA
  !
  IF ( lsda ) THEN
     !
     nspin = 2
     npol  = 1
     !
  ELSE IF ( noncolin ) THEN
     !
     nspin        = 4
     npol         = 2
     current_spin = 1
     !
  ELSE
     !
     nspin        = 1
     npol         = 1
     current_spin = 1
     !
  END IF
  !
  cell_factor = 1.D0
  lmovecell = .FALSE.
  !
  ! ... allocate memory for eigenvalues and weights (read from file)
  !
  nbndx = nbnd
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ) )
  !
  CALL pw_readfile( 'nowave', ierr )
  !
  ! ... distribute across pools k-points and related variables.
  ! ... nks is defined by the following routine as the number 
  ! ... of k-points in the current pool
  !
  CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  !
  CALL poolscatter( nbnd, nkstot, et, nks, et )
  CALL poolscatter( nbnd, nkstot, wg, nks, wg )
  !
  ! ... check on symmetry
  !
  IF (nat > 0) CALL checkallsym( nsym, s, nat, tau, &
                    ityp, at, bg, nr1, nr2, nr3, irt, ftau, alat, omega )
  !
  !  Set the different spin indices
  !
  nspin_mag  = nspin
  nspin_lsda = nspin
  nspin_gga  = nspin
  IF (nspin==4) THEN
     nspin_lsda=1
     IF (domag) THEN
        nspin_gga=2
     ELSE
        nspin_gga=1
        nspin_mag=1
     ENDIF
  ENDIF
  !
  ! ... read pseudopotentials
  !
  CALL pw_readfile( 'pseudo', ierr )
  !
  CALL readpp()
  !
  okvan = ANY ( upf(:)%tvanp )
  okpaw = ANY ( upf(1:nsp)%tpawp )
  !
  IF ( .NOT. lspinorb ) CALL average_pp ( nsp )
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL pre_init()
  CALL allocate_fft()
  CALL ggen()
  !
  ! ... allocate the potential and wavefunctions
  !
  CALL allocate_locpot()
  CALL allocate_nlpot()
  IF (okpaw) THEN
     CALL allocate_paw_internals()
     CALL paw_init_onecenter()
     CALL d_matrix(d1,d2,d3)
  ENDIF
  CALL allocate_wfc()
  !
  ! ... read the charge density
  !
  CALL pw_readfile( 'rho', ierr )
  !
  ! ... re-calculate the local part of the pseudopotential vltot
  ! ... and the core correction charge (if any) - This is done here
  ! ... for compatibility with the previous version of read_file
  !
  CALL init_vloc()
  !
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, &
                   nr1, nr2, nr3, strf, eigts1, eigts2, eigts3 )
  !
  CALL setlocal()
  !
  CALL set_rhoc()
  !
  ! ... bring rho to G-space
  !
  DO is = 1, nspin
     !
     psic(:) = rho%of_r(:,is)
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
     !
     rho%of_g(:,is) = psic(nl(:))
     !
  END DO
  !
  ! ... recalculate the potential
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, v )

  !
  ! ... reads the wavefunctions and writes them in 'distributed' form 
  ! ... to unit iunwfc (for compatibility)
  !
  nwordwfc = nbnd*npwx*npol
  !
  CALL open_buffer ( iunwfc, 'wfc', nwordwfc, nks, exst )
  !
  CALL pw_readfile( 'wave', ierr )
  !
  CALL init_us_1()
  IF (okpaw) then
     becsum = rho%bec
     CALL PAW_potential(rho%bec, ddd_PAW)
  ENDIF 
  if ( real_space ) THEN !initialisation of real space related stuff
    !OBM - correct parellism issues
    !call qpointlist()
    call betapointlist()
    call init_realspace_vars()
    !call betapointlist_v2()
    write(stdout,'(5X,"Real space initialisation completed")')
  endif
  CALL newd()
  CALL close_buffer  ( iunwfc, 'KEEP' )
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_dimensions()
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : pi
      USE cell_base, ONLY : alat, tpiba, tpiba2
      USE gvect,     ONLY : ecutwfc, dual, gcutm
      USE gsmooth,   ONLY : gcutms, doublegrid
      !
      !
      ! ... Set the units in real and reciprocal space
      !
      tpiba  = 2.D0 * pi / alat
      tpiba2 = tpiba**2
      !
      ! ... Compute the cut-off of the G vectors
      !
      gcutm = dual * ecutwfc / tpiba2
      !
      doublegrid = ( dual > 4.D0 )
      !
      IF ( doublegrid ) THEN
         !
         gcutms = 4.D0 * ecutwfc / tpiba2
         !
      ELSE
         !
         gcutms = gcutm
         !
      END IF
      !
    END SUBROUTINE set_dimensions
    !
END SUBROUTINE read_file
