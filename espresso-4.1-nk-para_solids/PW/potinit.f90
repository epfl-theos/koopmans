!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE potinit()
  !----------------------------------------------------------------------------
  !
  ! ... This routine initializes the self consistent potential in the array
  ! ... vr. There are three possible cases:
  !
  ! ... a) the code is restarting from a broken run:
  ! ...    read rho from data stored during the previous run
  ! ... b) the code is performing a non-scf calculation following a scf one:
  ! ...    read rho from the file produced by the scf calculation
  ! ... c) the code starts a new calculation:
  ! ...    calculate rho as a sum of atomic charges
  ! 
  ! ... In all cases the scf potential is recalculated and saved in vr
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : pi
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : alat, omega
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE basis,                ONLY : starting_pot
  USE klist,                ONLY : nelec
  USE lsda_mod,             ONLY : lsda, nspin
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   ngm, gstart, nl, g, gg
  USE gsmooth,              ONLY : doublegrid
  USE control_flags,        ONLY : lscf
  USE scf,                  ONLY : rho, rho_core, rhog_core, &
                                   vltot, v, vrs, kedtau
  USE funct,                ONLY : dft_is_meta
  USE wavefunctions_module, ONLY : psic
  USE ener,                 ONLY : ehart, etxc, vtxc, epaw
  USE ldaU,                 ONLY : niter_with_fixed_ns
  USE ldaU,                 ONLY : lda_plus_u, Hubbard_lmax, eth
  USE noncollin_module,     ONLY : noncolin, report
  USE io_files,             ONLY : tmp_dir, prefix, iunocc, input_drho
  USE spin_orb,             ONLY : domag
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE mp_global,            ONLY : intra_image_comm, inter_pool_comm, intra_pool_comm
  USE io_global,            ONLY : ionode, ionode_id
  USE pw_restart,           ONLY : pw_readfile
  USE io_rho_xml,           ONLY : read_rho
  USE xml_io_base,          ONLY : check_file_exst
  !
  USE uspp,                 ONLY : becsum
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_init,             ONLY : PAW_atomic_becsum
  USE paw_onecenter,        ONLY : PAW_potential
  !
  IMPLICIT NONE
  !
  REAL(DP)              :: charge           ! the starting charge
  REAL(DP)              :: etotefield       !
  REAL(DP)              :: fact
  INTEGER               :: is, ios
  LOGICAL               :: exst 
  CHARACTER(LEN=256)    :: filename
  !
  CALL start_clock('potinit')
  !
  ! check for both .dat and .xml extensions (compatibility reasons) 
  !
  filename =  TRIM( tmp_dir ) // TRIM( prefix ) // '.save/charge-density.dat'
  exst     =  check_file_exst( TRIM(filename) )
  !
  IF ( .NOT. exst ) THEN
      !
      filename =  TRIM( tmp_dir ) // TRIM( prefix ) // '.save/charge-density.xml'
      exst     =  check_file_exst( TRIM(filename) )
      !
  ENDIF
  !
  !
  IF ( starting_pot == 'file' .AND. exst ) THEN
     !
     ! ... Cases a) and b): the charge density is read from file
     ! ... this also reads rho%ns if lda+U and rho%bec if PAW
     !
     CALL pw_readfile( 'rho', ios )
     !
     IF ( ios /= 0 ) THEN
        !
        WRITE( stdout, '(/5X,"Error reading from file :"/5X,A,/)' ) &
            TRIM( filename )
        !
        CALL errore ( 'potinit' , 'reading starting density', ios)
        !
     ELSE IF ( lscf ) THEN
        !
        WRITE( stdout, '(/5X, &
             & "The initial density is read from file :"/5X,A,/)' ) &
            TRIM( filename )
        !
     ELSE
        !
        WRITE( stdout, '(/5X, &
             & "The potential is recalculated from file :"/5X,A,/)' ) &
            TRIM( filename )
        !
     END IF
     !
  ELSE
     !
     ! ... Case c): the potential is built from a superposition 
     ! ... of atomic charges contained in the array rho_at
     !
     IF ( starting_pot == 'file' .AND. .NOT. exst ) &
        WRITE( stdout, '(5X,"Cannot read rho : file not found")' )
     !
     WRITE( UNIT = stdout, &
            FMT = '(/5X,"Initial potential from superposition of free atoms")' )
     !
     !
     CALL atomic_rho( rho%of_r, nspin )
     ! ... in the lda+U case set the initial value of ns
     !
     IF ( lda_plus_u ) CALL init_ns()
     ! ... in the paw case uses atomic becsum
     IF ( okpaw )      CALL PAW_atomic_becsum()
     !
     IF ( input_drho /= ' ' ) THEN
        !
        IF ( nspin > 1 ) CALL errore &
             ( 'potinit', 'spin polarization not allowed in drho', 1 )
        !
        CALL read_rho ( v%of_r, 1, input_drho )
        !
        WRITE( UNIT = stdout, &
               FMT = '(/5X,"a scf correction to at. rho is read from",A)' ) &
            TRIM( input_drho )
        !
        rho%of_r = rho%of_r + v%of_r
        !
     END IF
     !
  END IF
  !
  ! ... check the integral of the starting charge
  !
  IF ( nspin == 2 ) THEN
     !
     charge = SUM ( rho%of_r(:,1:nspin) )*omega / ( nr1*nr2*nr3 )
     !
  ELSE
     !
     charge = SUM ( rho%of_r(:,1) )*omega / ( nr1*nr2*nr3 )
     !
  END IF
  !
  CALL mp_sum(  charge , intra_pool_comm )
  !
  IF ( lscf .AND. ABS( charge - nelec ) / charge > 1.D-7 ) THEN
     !
     WRITE( stdout, &
            '(/,5X,"starting charge ",F10.5,", renormalised to ",F10.5)') &
         charge, nelec
     !
     IF (nat>0) THEN
        rho%of_r = rho%of_r / charge * nelec
     ELSE
        rho%of_r = nelec / omega
     ENDIF
     !
  ELSE IF ( .NOT. lscf .AND. ABS( charge - nelec ) / charge > 1.D-3 ) THEN
     !
     CALL errore( 'potinit', 'starting and expected charges differ', 1 )
     !
  END IF
  !
  ! ... bring starting rho to G-space
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
  if ( dft_is_meta()) then
     ! ... define a starting (TF) guess for rho%kin_r and rho%kin_g
     fact = (3.d0*pi*pi)**(2.0/3.0)
     DO is = 1, nspin
        rho%kin_r(:,is) = fact * abs(rho%of_r(:,is)*nspin)**(5.0/3.0)/nspin
        psic(:) = rho%kin_r(:,is)
        CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
        rho%kin_g(:,is) = psic(nl(:))
     END DO
     !
  end if
  !
  ! ... compute the potential and store it in vr
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, v )
  IF (okpaw) CALL PAW_potential(rho%bec, ddd_PAW, epaw)
  !
  ! ... define the total local potential (external+scf)
  !
  CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, nrxx, nspin, doublegrid )
  !
  ! ... write on output the parameters used in the lda+U calculation
  !
  IF ( lda_plus_u ) THEN
     !
     WRITE( stdout, '(/5X,"Parameters of the lda+U calculation:")')
     WRITE( stdout, '(5X,"Number of iteration with fixed ns =",I3)') &
         niter_with_fixed_ns
     WRITE( stdout, '(5X,"Starting ns and Hubbard U :")')
     !
     CALL write_ns()
     !
  END IF
  !
  IF ( report /= 0 .AND. &
       noncolin .AND. domag .AND. lscf ) CALL report_mag()
  !
  CALL stop_clock('potinit')
  !
  RETURN
  !
END SUBROUTINE potinit
