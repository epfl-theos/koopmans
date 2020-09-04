!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE metadyn_vars
  !----------------------------------------------------------------------------
  !
  ! ... this module contains the variables necessary for the implementation of
  ! ... meta-dynamics and for the calculation of free-energy barriers by means 
  ! ... of the fourier string method
  !
  ! ... code written by Carlo Sbraccia (2005)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER :: ncolvar
  !
  INTEGER  :: fe_nstep, sw_nstep, eq_nstep
  !
  REAL(DP), ALLOCATABLE :: dfe_acc(:)
  REAL(DP), ALLOCATABLE :: fe_grad(:)
  REAL(DP), ALLOCATABLE :: new_target(:)
  REAL(DP), ALLOCATABLE :: to_target(:)
  !
  LOGICAL  :: to_new_target
  !
  INTEGER  :: max_metadyn_iter
  !
  REAL(DP), ALLOCATABLE :: fe_step(:)
  REAL(DP), ALLOCATABLE :: gaussian_pos(:)
  !
  REAL(DP), ALLOCATABLE :: metadyn_history(:,:)
  !
  REAL(DP) :: g_amplitude
  !
  REAL(DP) :: etot_av
  !
  INTEGER :: first_metadyn_iter
  !
  CHARACTER(LEN=80) :: metadyn_fmt
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE init_metadyn_vars()
      !------------------------------------------------------------------------
      !
      USE input_parameters, ONLY : ncolvar_inp, &
                                   g_amplitude_ => g_amplitude, &
                                   fe_step_     => fe_step,  &
                                   fe_nstep_    => fe_nstep, &
                                   sw_nstep_    => sw_nstep, &
                                   eq_nstep_    => eq_nstep
      USE control_flags,    ONLY : lmetadyn, nstep
      !
      IMPLICIT NONE
      !
      !
      ncolvar = ncolvar_inp
      !
      ALLOCATE( fe_step(    ncolvar ) )
      ALLOCATE( dfe_acc(    ncolvar ) )
      ALLOCATE( fe_grad(    ncolvar ) )
      ALLOCATE( new_target( ncolvar ) )
      ALLOCATE( to_target(  ncolvar ) )
      !
      IF ( lmetadyn ) THEN
         !
         ALLOCATE( gaussian_pos( ncolvar ) )
         ALLOCATE( metadyn_history( ncolvar, nstep ) )
         !
      END IF
      !
      fe_nstep    = fe_nstep_
      sw_nstep    = sw_nstep_
      eq_nstep    = eq_nstep_
      g_amplitude = g_amplitude_
      fe_step(:)  = fe_step_(1:ncolvar)
      !
      RETURN
      !
    END SUBROUTINE init_metadyn_vars
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_metadyn_vars()
      !------------------------------------------------------------------------
      !
      IF ( ALLOCATED( fe_step ) )          DEALLOCATE( fe_step )
      IF ( ALLOCATED( dfe_acc ) )          DEALLOCATE( dfe_acc )
      IF ( ALLOCATED( fe_grad ) )          DEALLOCATE( fe_grad )
      IF ( ALLOCATED( new_target ) )       DEALLOCATE( new_target )
      IF ( ALLOCATED( to_target ) )        DEALLOCATE( to_target )
      IF ( ALLOCATED( metadyn_history ) )  DEALLOCATE( metadyn_history )
      IF ( ALLOCATED( gaussian_pos ) )     DEALLOCATE( gaussian_pos )
      !
      RETURN
      !
    END SUBROUTINE deallocate_metadyn_vars
    !
END MODULE metadyn_vars
