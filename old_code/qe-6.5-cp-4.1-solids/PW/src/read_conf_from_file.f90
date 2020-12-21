!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE read_conf_from_file( lmovecell, at_old, omega_old, ierr)
  !-----------------------------------------------------------------------
  ! FIXME: half of the variables are passed as arguments, half in modules
  ! FIXME: this routines does two different things
  !
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout, ionode, ionode_id
  USE io_files,        ONLY : restart_dir, xmlfile, &
                              psfile, pseudo_dir, pseudo_dir_cur
  USE ions_base,       ONLY : nat, nsp, ityp, amass, atm, tau
  USE cell_base,       ONLY : alat, ibrav, at, bg, omega
  USE mp,              ONLY : mp_bcast
  USE mp_images,       ONLY : intra_image_comm
  USE qexsd_module,    ONLY : qexsd_readschema
  USE qexsd_copy,      ONLY : qexsd_copy_atomic_species, &
                              qexsd_copy_atomic_structure
  USE qes_types_module,ONLY : output_type
  USE qes_libs_module, ONLY : qes_reset
  USE qes_bcast_module,ONLY : qes_bcast
  !
  IMPLICIT NONE
  !
  LOGICAL,INTENT(in)     :: lmovecell
  REAL(DP),INTENT(inout) :: at_old(3,3), omega_old
  INTEGER, INTENT(out)   :: ierr
  !
  TYPE ( output_type) :: output_obj
  INTEGER :: nat_
  !
  pseudo_dir_cur = restart_dir () 
  WRITE( stdout, '(/5X,"Atomic positions and unit cell read from directory:", &
                &  /,5X,A)') pseudo_dir_cur
  !
  ! ... check if restart file is present, if so read config parameters
  !
  IF (ionode) CALL qexsd_readschema ( xmlfile(), ierr, output_obj )
  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
  IF ( ierr > 0 ) CALL errore ( 'read_conf_from_file', &
       'fatal error reading xml file', ierr ) 
  CALL qes_bcast(output_obj, ionode_id, intra_image_comm)
  !
  IF (ierr == 0 ) THEN
     !
     CALL qexsd_copy_atomic_species ( output_obj%atomic_species, &
          nsp, atm, amass, PSFILE=psfile, PSEUDO_DIR=pseudo_dir )
     IF ( pseudo_dir == ' ' ) pseudo_dir=pseudo_dir_cur
     CALL qexsd_copy_atomic_structure (output_obj%atomic_structure, nsp, &
          atm, nat_, tau, ityp, alat, at(:,1), at(:,2), at(:,3), ibrav )
     CALL qes_reset (output_obj)
     IF ( nat_ /= nat ) CALL errore('read_conf_from_file','bad atom number',1)
     at(:,:) = at(:,:) / alat
     tau(:,1:nat) = tau(:,1:nat)/alat  
     CALL volume (alat,at(:,1),at(:,2),at(:,3),omega)
     CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
     !
  ELSE
     !
     WRITE( stdout, '(5X,"Nothing found: ", &
                       & "using input atomic positions and unit cell",/)' )
     RETURN
     !
  END IF
  !
  WRITE( stdout, * )
  !
  IF ( lmovecell ) THEN
     !
     ! ... input value of at and omega (currently stored in xxx_old variables)
     ! ... must be used to initialize G vectors and other things
     ! ... swap xxx and xxx_old variables and scale the atomic position to the
     ! ... input cell shape in order to check the symmetry.
     !
     CALL cryst_to_cart( nat, tau, bg, - 1 )
     CALL dswap( 9, at, 1, at_old,1  )
     CALL dswap( 1, omega, 1, omega_old, 1 )
     CALL cryst_to_cart( nat, tau, at, + 1 )
     !
     CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE read_conf_from_file
