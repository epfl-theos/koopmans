!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! written by Carlo Cavazzoni

!-----------------------------------------------------------------------

   SUBROUTINE writefile_cp_real                                         &
     &     ( h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,           &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,nhpdim,ekincm,&
     &       xnhh0,xnhhm,vnhh,velh, fion, tps, mat_z, occ_f, rho )
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!
      USE kinds,            ONLY: DP
      USE ions_base,        ONLY: cdmi, taui
      USE cell_base,        ONLY: s_to_r
      USE cp_restart,       ONLY: cp_writefile
      USE cp_interfaces,    ONLY: set_evtot, set_eitot
      USE electrons_base,   ONLY: nspin, iupdwn, nupdwn
      USE electrons_module, ONLY: ei, nupdwn_emp
      USE io_files,         ONLY: outdir
      USE ensemble_dft,     ONLY: tens, tsmear
      USE mp,               ONLY: mp_bcast
      USE control_flags,    ONLY: tksw, ndw, evc_restart
!
      implicit none
      integer, INTENT(IN) ::  nfi
      REAL(DP), INTENT(IN) :: h(3,3), hold(3,3)
      complex(DP), INTENT(IN) :: c0(:,:), cm(:,:)
      REAL(DP), INTENT(IN) :: tausm(:,:), taus(:,:), fion(:,:)
      REAL(DP), INTENT(IN) :: vels(:,:), velsm(:,:)
      REAL(DP), INTENT(IN) :: acc(:), lambda(:,:,:), lambdam(:,:,:)
      REAL(DP), INTENT(IN) :: xnhe0, xnhem, vnhe, ekincm
      REAL(DP), INTENT(IN) :: xnhp0(:), xnhpm(:), vnhp(:)
      integer,      INTENT(in) :: nhpcl, nhpdim
      REAL(DP), INTENT(IN) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
      REAL(DP), INTENT(in) :: tps
      REAL(DP), INTENT(in) :: rho(:,:)
      REAL(DP), INTENT(in) :: occ_f(:)
      REAL(DP), INTENT(in) :: mat_z(:,:,:)

      REAL(DP) :: ht(3,3), htm(3,3), htvel(3,3), gvel(3,3)
      INTEGER  :: nk = 1
      REAL(DP) :: xk(3,1) = 0.0d0, wk(1) = 2.0d0
      COMPLEX(DP), ALLOCATABLE :: ctot(:,:)
      REAL(DP),    ALLOCATABLE :: eitot(:,:)
      INTEGER  :: nupdwn_tot( 2 ), iupdwn_tot( 2 )

      write(6,*) "writefile_cp_real: this subroutine is not working in this implementation" !added:giovanni:debug
      stop
      if ( ndw < 1 ) then
         !
         ! Do not write restart file if the unit number 
         ! is negative, this is used mainly for benchmarks and tests
         !
         return
         !
      end if

      ht     = TRANSPOSE( h ) 
      htm    = TRANSPOSE( hold ) 
      htvel  = TRANSPOSE( velh ) 
      gvel   = 0.0d0
      
      nupdwn_tot = nupdwn + nupdwn_emp
      iupdwn_tot(1) = iupdwn(1)
      iupdwn_tot(2) = nupdwn_tot(1) + 1  !! NlN check if it's correct through all the routine
      !
      ALLOCATE( eitot( nupdwn_tot(1), nspin ) )
      !
      CALL set_eitot( eitot )
      !
      IF( tksw .or. evc_restart ) THEN
         !
         ALLOCATE( ctot( SIZE( c0, 1 ), nupdwn_tot(1) * nspin ) )
         !
         CALL set_evtot( c0, ctot, lambda, iupdwn_tot, nupdwn_tot )
         !
      END IF
      !
      !  Sincronize lambdas, whose replicas could diverge on
      !  different processors
      !
      IF( tens ) THEN
        !
        CALL cp_writefile( ndw, outdir, .TRUE., nfi, tps, acc, nk, xk, wk,   &
             ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi , taus,        &
             vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim, occ_f , &
             occ_f , lambda, lambdam, xnhe0, xnhem, vnhe, ekincm, ei,&
             rho, c0, cm, ctot, iupdwn, nupdwn, iupdwn, nupdwn, mat_z = mat_z  )
        !
      ELSE IF ( tsmear ) THEN
        ! 
        CALL cp_writefile( ndw, outdir, .TRUE., nfi, tps, acc, nk, xk, wk,      &
             ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi , taus,        &
             vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim, occ_f , &
             occ_f , lambda, lambdam, xnhe0, xnhem, vnhe, ekincm, eitot,&
             rho, c0, cm, ctot, iupdwn, nupdwn, iupdwn_tot, nupdwn_tot, mat_z = mat_z )
        !
      ELSE
        ! 
        CALL cp_writefile( ndw, outdir, .TRUE., nfi, tps, acc, nk, xk, wk,      &
             ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi , taus,        &
             vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim, occ_f , &
             occ_f , lambda, lambdam, xnhe0, xnhem, vnhe, ekincm, eitot,&
             rho, c0, cm, ctot, iupdwn, nupdwn, iupdwn_tot, nupdwn_tot )
        !
      END IF

      DEALLOCATE( eitot )
      !
      IF( tksw .or. evc_restart ) DEALLOCATE( ctot )

      return
      end subroutine writefile_cp_real

   SUBROUTINE writefile_cp_twin                                         &
     &     ( h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,           &
     &       lambda,lambdam,lambda_bare,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,& 
     &       nhpcl,nhpdim,ekincm,&
     &       xnhh0,xnhhm,vnhh,velh, fion, tps, mat_z, occ_f, rho )
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!
      USE kinds,            ONLY: DP
      USE ions_base,        ONLY: cdmi, taui
      USE cell_base,        ONLY: s_to_r
      USE cp_restart,       ONLY: cp_writefile
      USE cp_interfaces,    ONLY: set_evtot, set_eitot
      USE electrons_base,   ONLY: nspin, iupdwn, nupdwn
      USE electrons_module, ONLY: ei, nupdwn_emp
      USE io_files,         ONLY: outdir
      USE ensemble_dft,     ONLY: tens, tsmear
      USE mp,               ONLY: mp_bcast
      USE control_flags,    ONLY: tksw, ndw, evc_restart
      USE electrons_module, ONLY: wfc_spreads
      USE nksic,            ONLY: do_orbdep
      USE twin_types
!
      implicit none
      integer, INTENT(IN) ::  nfi
      REAL(DP), INTENT(IN) :: h(3,3), hold(3,3)
      complex(DP), INTENT(IN) :: c0(:,:), cm(:,:)
      REAL(DP), INTENT(IN) :: tausm(:,:), taus(:,:), fion(:,:)
      REAL(DP), INTENT(IN) :: vels(:,:), velsm(:,:)
      REAL(DP), INTENT(IN) :: acc(:)
      TYPE(twin_matrix), dimension(:), INTENT(IN) :: lambda, lambdam, lambda_bare
      REAL(DP), INTENT(IN) :: xnhe0, xnhem, vnhe, ekincm
      REAL(DP), INTENT(IN) :: xnhp0(:), xnhpm(:), vnhp(:)
      integer,      INTENT(in) :: nhpcl, nhpdim
      REAL(DP), INTENT(IN) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
      REAL(DP), INTENT(in) :: tps
      REAL(DP), INTENT(in) :: rho(:,:)
      REAL(DP), INTENT(in) :: occ_f(:)
      TYPE(twin_matrix), dimension(:), INTENT(IN) :: mat_z

      REAL(DP) :: ht(3,3), htm(3,3), htvel(3,3), gvel(3,3)
      INTEGER  :: nk = 1
      REAL(DP) :: xk(3,1) = 0.0d0, wk(1) = 2.0d0
      COMPLEX(DP), ALLOCATABLE :: ctot(:,:)
      REAL(DP),    ALLOCATABLE :: eitot(:,:)
      INTEGER  :: nupdwn_tot( 2 ), iupdwn_tot( 2 )


      if ( ndw < 1 ) then
         !
         ! Do not write restart file if the unit number 
         ! is negative, this is used mainly for benchmarks and tests
         !
         return
         !
      end if

      ht     = TRANSPOSE( h ) 
      htm    = TRANSPOSE( hold ) 
      htvel  = TRANSPOSE( velh ) 
      gvel   = 0.0d0
      
      nupdwn_tot = nupdwn + nupdwn_emp
      iupdwn_tot(1) = iupdwn(1)
      iupdwn_tot(2) = nupdwn_tot(1) + 1  !! NlN check if it's correct through all the routine
      !
      ALLOCATE( eitot( nupdwn_tot(1), nspin ) )
      !
      CALL set_eitot( eitot )
      !
      IF( tksw .or. evc_restart ) THEN
         !
         ALLOCATE( ctot( SIZE( c0, 1 ), nupdwn_tot(1) * nspin ) )
         !
         CALL set_evtot( c0, ctot, lambda, iupdwn_tot, nupdwn_tot )
         !
      END IF
      !
      !  Sincronize lambdas, whose replicas could diverge on
      !  different processors
      !
      IF(do_orbdep) THEN !Sort wavefunctions with respect to spread !!added:giovanni
         !
         IF(allocated(wfc_spreads)) THEN
            !
            !call spread_sort(c0, ngw, nspin, nbsp, nudx, nupdwn, iupdwn, wfc_spreads)
            !
         ENDIF
         !
      ENDIF

      IF( tens ) THEN
        !
        CALL cp_writefile( ndw, outdir, .TRUE., nfi, tps, acc, nk, xk, wk,   &
             ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi , taus,        &
             vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim, occ_f , &
             occ_f , lambda, lambdam, lambda_bare, xnhe0, xnhem, vnhe, ekincm, ei,&
             rho, c0, cm, ctot, iupdwn, nupdwn, iupdwn, nupdwn, mat_z = mat_z  )
        !
      ELSE IF ( tsmear ) THEN
        ! 
        CALL cp_writefile( ndw, outdir, .TRUE., nfi, tps, acc, nk, xk, wk,      &
             ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi , taus,        &
             vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim, occ_f , &
             occ_f , lambda, lambdam, lambda_bare, xnhe0, xnhem, vnhe, ekincm, eitot,&
             rho, c0, cm, ctot, iupdwn, nupdwn, iupdwn_tot, nupdwn_tot, mat_z = mat_z )
        !
      ELSE
        ! 
        CALL cp_writefile( ndw, outdir, .TRUE., nfi, tps, acc, nk, xk, wk,      &
             ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi , taus,        &
             vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim, occ_f , &
             occ_f , lambda, lambdam, lambda_bare, xnhe0, xnhem, vnhe, ekincm, eitot,&
             rho, c0, cm, ctot, iupdwn, nupdwn, iupdwn_tot, nupdwn_tot )
        !
      END IF

      DEALLOCATE( eitot )
      !
      IF( tksw .or. evc_restart ) DEALLOCATE( ctot )

      return
      end subroutine writefile_cp_twin

!-----------------------------------------------------------------------
      subroutine readfile_cp_real                                        &
     &     ( flag, h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,    &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,nhpdim,ekincm,&
     &       xnhh0,xnhhm,vnhh,velh,&
     &       fion, tps, mat_z, occ_f )
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!
      USE kinds,          ONLY : DP
      USE io_files,       ONLY : outdir
      USE electrons_base, ONLY : nspin, keep_occ
      USE ions_base,      ONLY : cdmi, taui
      USE cp_restart,     ONLY : cp_readfile, cp_read_cell, cp_read_wfc
      USE ensemble_dft,   ONLY : tens, tsmear
      USE autopilot,      ONLY : event_step, event_index, max_event_step
      USE cp_autopilot,   ONLY : employ_rules
      USE control_flags,  ONLY : ndr
!
      implicit none
      INTEGER, INTENT(in) :: flag
      integer ::  nfi
      REAL(DP) :: h(3,3), hold(3,3)
      complex(DP) :: c0(:,:), cm(:,:)
      REAL(DP) :: tausm(:,:),taus(:,:), fion(:,:)
      REAL(DP) :: vels(:,:), velsm(:,:)
      REAL(DP) :: acc(:),lambda(:,:,:), lambdam(:,:,:)
      REAL(DP) :: xnhe0,xnhem,vnhe
      REAL(DP) :: xnhp0(:), xnhpm(:), vnhp(:)
      integer, INTENT(inout) :: nhpcl,nhpdim
      REAL(DP) :: ekincm
      REAL(DP) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
      REAL(DP), INTENT(OUT) :: tps
      REAL(DP), INTENT(INOUT) :: mat_z(:,:,:), occ_f(:)
      !
      REAL(DP) :: ht(3,3), htm(3,3), htvel(3,3), gvel(3,3)
      integer :: nk = 1, ispin
      REAL(DP) :: xk(3,1) = 0.0d0, wk(1) = 2.0d0
      REAL(DP), ALLOCATABLE :: occ_ ( : )
      REAL(DP) :: b1(3) , b2(3), b3(3)
       

      IF( flag == -1 ) THEN
        CALL cp_read_cell( ndr, outdir, .TRUE., ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh )
        h     = TRANSPOSE( ht )
        hold  = TRANSPOSE( htm )
        velh  = TRANSPOSE( htvel )
        RETURN
      ELSE IF ( flag == 0 ) THEN
        DO ispin = 1, nspin
          CALL cp_read_wfc( ndr, outdir, 1, 1, ispin, nspin, c2 = cm(:,:), tag = 'm' )
        END DO
        RETURN
      END IF

      ALLOCATE( occ_ ( SIZE( occ_f ) ) )

      IF( tens ) THEN
         CALL cp_readfile( ndr, outdir, .TRUE., nfi, tps, acc, nk, xk, wk, &
                ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi, taus, &
                vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ_ , &
                occ_ , lambda, lambdam, b1, b2, b3, &
                xnhe0, xnhem, vnhe, ekincm, c0, cm, mat_z = mat_z )
      ELSE IF ( tsmear ) THEN
         CALL cp_readfile( ndr, outdir, .TRUE., nfi, tps, acc, nk, xk, wk, &
                ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi, taus, &
                vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ_ , &
                occ_ , lambda, lambdam, b1, b2, b3, &
                xnhe0, xnhem, vnhe, ekincm, c0, cm, mat_z = mat_z )
      ELSE
         CALL cp_readfile( ndr, outdir, .TRUE., nfi, tps, acc, nk, xk, wk, &
                ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi, taus, &
                vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ_ , &
                occ_ , lambda, lambdam, b1, b2, b3, &
                xnhe0, xnhem, vnhe, ekincm, c0, cm )
      END IF
      ! AutoPilot (Dynamic Rules) Implementation
      event_index = 1

      do while (event_step(event_index) <= nfi)
         ! Assuming that the remaining dynamic parm values are set correctly by reading 
         ! the the restart file.
         ! if this is not true, employ rules as events are updated right here as:
         call employ_rules()
         event_index = event_index + 1
         if(  event_index > max_event_step ) then
            CALL errore( ' readfile ' , ' maximum events exceeded for dynamic rules ' , 1 )
         end if
      enddo
      
      IF( .NOT. keep_occ ) THEN
         occ_f( : ) = occ_ ( : )
      END IF

      DEALLOCATE( occ_ )

      return
      end subroutine readfile_cp_real

!-----------------------------------------------------------------------
      subroutine readfile_cp_twin                                        &
     &     ( flag, h,hold,nfi,c0,cm,taus,tausm,vels,velsm,acc,    &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,nhpdim,ekincm,&
     &       xnhh0,xnhhm,vnhh,velh,&
     &       fion, tps, mat_z, occ_f )
!-----------------------------------------------------------------------
!
! read from file and distribute data calculated in preceding iterations
!
      USE kinds,          ONLY : DP
      USE io_files,       ONLY : outdir
      USE electrons_base, ONLY : nspin, keep_occ
      USE ions_base,      ONLY : cdmi, taui
      USE cp_restart,     ONLY : cp_readfile, cp_read_cell, cp_read_wfc
      USE ensemble_dft,   ONLY : tens, tsmear
      USE autopilot,      ONLY : event_step, event_index, max_event_step
      USE cp_autopilot,   ONLY : employ_rules
      USE control_flags,  ONLY : ndr
      USE twin_types
!
      implicit none
      INTEGER, INTENT(in) :: flag
      integer ::  nfi
      REAL(DP) :: h(3,3), hold(3,3)
      complex(DP) :: c0(:,:), cm(:,:)
      REAL(DP) :: tausm(:,:),taus(:,:), fion(:,:)
      REAL(DP) :: vels(:,:), velsm(:,:)
      REAL(DP) :: acc(:)
      TYPE(twin_matrix), DIMENSION(:) :: lambda, lambdam
      REAL(DP) :: xnhe0,xnhem,vnhe
      REAL(DP) :: xnhp0(:), xnhpm(:), vnhp(:)
      integer, INTENT(inout) :: nhpcl,nhpdim
      REAL(DP) :: ekincm
      REAL(DP) :: xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
      REAL(DP), INTENT(OUT) :: tps
      REAL(DP), INTENT(INOUT) :: occ_f(:)
      TYPE(twin_matrix), DIMENSION(:) :: mat_z
      !
      REAL(DP) :: ht(3,3), htm(3,3), htvel(3,3), gvel(3,3)
      integer :: nk = 1, ispin
      REAL(DP) :: xk(3,1) = 0.0d0, wk(1) = 2.0d0
      REAL(DP), ALLOCATABLE :: occ_ ( : )
      REAL(DP) :: b1(3) , b2(3), b3(3)
        

      IF( flag == -1 ) THEN
        CALL cp_read_cell( ndr, outdir, .TRUE., ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh )
        h     = TRANSPOSE( ht )
        hold  = TRANSPOSE( htm )
        velh  = TRANSPOSE( htvel )
        RETURN
      ELSE IF ( flag == 0 ) THEN
        DO ispin = 1, nspin
          CALL cp_read_wfc( ndr, outdir, 1, 1, ispin, nspin, c2 = cm(:,:), tag = 'm' )
        END DO
        RETURN
      END IF

      ALLOCATE( occ_ ( SIZE( occ_f ) ) )

      IF( tens ) THEN
         CALL cp_readfile( ndr, outdir, .TRUE., nfi, tps, acc, nk, xk, wk, &
                ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi, taus, &
                vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ_ , &
                occ_ , lambda, lambdam, b1, b2, b3, &
                xnhe0, xnhem, vnhe, ekincm, c0, cm, mat_z = mat_z )
      ELSE IF ( tsmear ) THEN
         CALL cp_readfile( ndr, outdir, .TRUE., nfi, tps, acc, nk, xk, wk, &
                ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi, taus, &
                vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ_ , &
                occ_ , lambda, lambdam, b1, b2, b3, &
                xnhe0, xnhem, vnhe, ekincm, c0, cm, mat_z = mat_z )
      ELSE
         CALL cp_readfile( ndr, outdir, .TRUE., nfi, tps, acc, nk, xk, wk, &
                ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh, taui, cdmi, taus, &
                vels, tausm, velsm, fion, vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ_ , &
                occ_ , lambda, lambdam, b1, b2, b3, &
                xnhe0, xnhem, vnhe, ekincm, c0, cm )
      END IF

      ! AutoPilot (Dynamic Rules) Implementation
      event_index = 1

      do while (event_step(event_index) <= nfi)
         ! Assuming that the remaining dynamic parm values are set correctly by reading 
         ! the the restart file.
         ! if this is not true, employ rules as events are updated right here as:
         call employ_rules()
         event_index = event_index + 1
         if(  event_index > max_event_step ) then
            CALL errore( ' readfile ' , ' maximum events exceeded for dynamic rules ' , 1 )
         end if
      enddo
      
      IF( .NOT. keep_occ ) THEN
         occ_f( : ) = occ_ ( : )
      END IF

      DEALLOCATE( occ_ )

      return
      end subroutine readfile_cp_twin

!=----------------------------------------------------------------------------=!


   SUBROUTINE writefile_fpmd &
     ( nfi, trutime, c0, cm, occ, atoms_0, atoms_m, acc, taui, cdmi, ht_m, &
       ht_0, rho, lambda, tlast )
                                                                        
        USE kinds,             ONLY: DP
        USE cell_base,         ONLY: boxdimensions, r_to_s
        USE control_flags,     ONLY: ndw
        USE control_flags,     ONLY: tksw, evc_restart
        USE atoms_type_module, ONLY: atoms_type
        USE electrons_nose,    ONLY: xnhe0, xnhem, vnhe
        USE electrons_base,    ONLY: nspin, iupdwn, nupdwn
        USE cell_nose,         ONLY: xnhh0, xnhhm, vnhh
        USE ions_nose,         ONLY: vnhp, xnhp0, xnhpm, nhpcl, nhpdim
        USE cp_restart,        ONLY: cp_writefile
        USE electrons_module,  ONLY: nupdwn_emp
        USE io_files,          ONLY: outdir
        USE cp_interfaces,     ONLY: set_evtot, set_eitot
        USE kohn_sham_states,  ONLY: print_all_states

        IMPLICIT NONE 
 
        INTEGER,              INTENT(IN)    :: nfi
        COMPLEX(DP),          INTENT(IN)    :: c0(:,:), cm(:,:) 
        REAL(DP),             INTENT(IN)    :: occ(:)
        TYPE (boxdimensions), INTENT(IN)    :: ht_m, ht_0
        TYPE (atoms_type),    INTENT(IN)    :: atoms_0, atoms_m
        REAL(DP),             INTENT(IN)    :: rho(:,:)
        REAL(DP),             INTENT(IN)    :: taui(:,:)
        REAL(DP),             INTENT(IN)    :: acc(:), cdmi(:) 
        REAL(DP),             INTENT(IN)    :: trutime
        REAL(DP),             INTENT(IN)    :: lambda(:,:,:)
        LOGICAL,              INTENT(IN)    :: tlast

        REAL(DP) :: ekincm
        COMPLEX(DP), ALLOCATABLE :: ctot(:,:)
        REAL(DP),    ALLOCATABLE :: eitot(:,:)
        INTEGER  :: nupdwn_tot( 2 ), iupdwn_tot( 2 )

        INTEGER  :: nkpt = 1
        REAL(DP) :: xk(3,1) = 0.0d0, wk(1) = 2.0d0
             
        IF( ndw < 1 ) RETURN
        !
        !   this is used for benchmarking and debug
        !   if ndw < 1 Do not save wave functions and other system
        !   properties on the writefile subroutine

        ekincm = 0.0d0
        !
        nupdwn_tot = nupdwn + nupdwn_emp
        iupdwn_tot(1) = iupdwn(1)
!        iupdwn_tot(2) = nupdwn(1) + 1
        iupdwn_tot(2) = nupdwn_tot(1) + 1  !! NlN check if it's correct through all the routine
        !
        ALLOCATE( eitot( nupdwn_tot(1), nspin ) )
        !
        CALL set_eitot( eitot )

        IF( tksw .or. evc_restart ) THEN
           !
           ALLOCATE( ctot( SIZE( c0, 1 ), nupdwn_tot(1) * nspin ) )
           !
           CALL set_evtot( c0, ctot, lambda, iupdwn_tot, nupdwn_tot )
           !
           IF( tlast ) CALL print_all_states( ctot, iupdwn_tot, nupdwn_tot )
           !
        END IF
        !
        CALL cp_writefile( ndw, outdir, .TRUE., nfi, trutime, acc, nkpt, xk, wk, &
          ht_0%a, ht_m%a, ht_0%hvel, ht_0%gvel, xnhh0, xnhhm, vnhh, taui, cdmi,   &
          atoms_0%taus, atoms_0%vels, atoms_m%taus, atoms_m%vels, atoms_0%for,    &
          vnhp, xnhp0, xnhpm, nhpcl, nhpdim, occ, occ, lambda, lambda,            &
          xnhe0, xnhem, vnhe, ekincm, eitot, rho, c0, cm, ctot, iupdwn, nupdwn,   &
          iupdwn_tot, nupdwn_tot )

        DEALLOCATE( eitot )

        IF( tksw .or. evc_restart ) DEALLOCATE( ctot )

     RETURN 
   END SUBROUTINE writefile_fpmd


!=----------------------------------------------------------------------------=!


!------------------------------------------------------------------------------!
   SUBROUTINE set_eitot_x( eitot )
!------------------------------------------------------------------------------!
      USE kinds,            ONLY: DP
      USE electrons_base,   ONLY: nupdwn, nspin
      USE electrons_module, ONLY: nupdwn_emp, ei, ei_emp, n_emp
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(OUT) :: eitot(:,:)
      !
      INTEGER :: n
      !
      eitot = 0.0d0
      !
      eitot( 1:nupdwn(1), 1 ) = ei( 1:nupdwn(1), 1 )
      IF( nspin == 2 ) eitot( 1:nupdwn(2), 2 ) = ei( 1:nupdwn(2), 2 )
      !  
      IF( n_emp > 0 ) THEN
         ! 
         n = nupdwn(1)
         eitot( n+1 : n+nupdwn_emp(1), 1 ) = ei_emp( 1 : nupdwn_emp(1), 1 )
         IF( nspin == 2 ) THEN
            n = nupdwn(2)
            eitot( n+1 : n+nupdwn_emp(2), 2 ) = ei_emp( 1 : nupdwn_emp(2), 2 )
         END IF
         !
      END IF
      !
      RETURN
   END SUBROUTINE set_eitot_x

!------------------------------------------------------------------------------!
   SUBROUTINE set_evtot_real_x( c0, ctot, lambda, iupdwn_tot, nupdwn_tot )
!------------------------------------------------------------------------------!
      USE kinds,             ONLY: DP
      USE electrons_base,    ONLY: nupdwn, nspin, iupdwn, nudx
      USE electrons_module,  ONLY: nupdwn_emp, n_emp, iupdwn_emp
      USE cp_interfaces,     ONLY: readempty, crot, readempty_twin
      USE cp_main_variables, ONLY: collect_lambda, descla
      USE control_flags,     ONLY: ndw
      USE input_parameters,  ONLY: print_evc0_occ_empty 
      USE wavefunctions_module, ONLY: c0_occ_emp_aux
      !
      IMPLICIT NONE
      !
      COMPLEX(DP), INTENT(IN)  :: c0(:,:)
      COMPLEX(DP), INTENT(OUT) :: ctot(:,:)
      REAL(DP),    INTENT(IN)  :: lambda(:,:,:)
      INTEGER,     INTENT(IN)  :: iupdwn_tot(2), nupdwn_tot(2)
      !
      COMPLEX(DP), ALLOCATABLE :: cemp(:,:)
      REAL(DP),    ALLOCATABLE :: eitmp(:)
      REAL(DP),    ALLOCATABLE :: lambda_repl(:,:)
      LOGICAL                  :: t_emp
      !
      ALLOCATE( eitmp( nudx ) )
      ALLOCATE( lambda_repl( nudx, nudx ) )
      !
      ctot = 0.0d0
      !
      CALL collect_lambda( lambda_repl, lambda(:,:,1), descla(:,1) )
      !
      CALL crot( ctot, c0, SIZE( c0, 1 ), nupdwn(1), iupdwn_tot(1), iupdwn(1), lambda_repl, nudx, eitmp )
      !
      IF( nspin == 2 ) THEN
         CALL collect_lambda( lambda_repl, lambda(:,:,2), descla(:,2) )
         CALL crot( ctot, c0, SIZE( c0, 1 ), nupdwn(2), iupdwn_tot(2), iupdwn(2), lambda_repl, nudx, eitmp )
      END IF
      !
      DEALLOCATE( lambda_repl )
      !
      t_emp = .FALSE.
      !
      IF( n_emp > 0 ) THEN
          !
          ALLOCATE( cemp( SIZE( c0, 1 ), n_emp * nspin ) )
          cemp = 0.0d0
          t_emp = readempty( cemp, n_emp * nspin, ndw )
          !
      END IF
      !
      IF( t_emp ) THEN
         ctot( :, nupdwn( 1 )+1 : nupdwn_tot(1) ) = cemp( :, 1:nupdwn_emp( 1 ) )
         IF( nspin == 2 ) THEN
             ctot( :, iupdwn_tot(2) + nupdwn(2) : iupdwn_tot(2) + nupdwn_tot(1) - 1 ) = &
                cemp( :, iupdwn_emp(2) : iupdwn_emp(2) + nupdwn_emp(2) - 1 )
         END IF
      END IF
      !
      IF( n_emp > 0 ) DEALLOCATE( cemp )
      !
      ! print evc0 of occ and empty in xml_io format
      ! 
      IF (print_evc0_occ_empty .and. (n_emp > 0)) THEN
         !
         ALLOCATE( c0_occ_emp_aux( SIZE( c0, 1 ), nupdwn_tot(1) * nspin ) )
         !
         t_emp = .FALSE.
         !
         ALLOCATE( cemp( SIZE( c0, 1 ), n_emp * nspin ) )
         cemp = 0.0d0
         t_emp = readempty_twin( cemp,  n_emp * nspin, ndw )
         !
         IF (t_emp) THEN
            c0_occ_emp_aux(:, iupdwn_tot(1) : nupdwn(1)) = c0(:, iupdwn(1):nupdwn(1))
            c0_occ_emp_aux(:, nupdwn( 1 )+1 : nupdwn_tot(1)) = cemp( :, 1:nupdwn_emp(1)) 
            IF ( nspin == 2 ) THEN
               c0_occ_emp_aux(:, iupdwn_tot(2) : iupdwn_tot(2) + nupdwn(2) -1 ) = c0(:, iupdwn(2):iupdwn(2) + nupdwn(2)-1)
               c0_occ_emp_aux(:, iupdwn_tot(2) + nupdwn(2) : iupdwn_tot(2) + nupdwn_tot(2) - 1 ) = &
               cemp( :, iupdwn_emp(2) : iupdwn_emp(2) + nupdwn_emp(2) - 1 )
            ENDIF
         ENDIF
         !
         DEALLOCATE( cemp )
      ENDIF
      ! 
      DEALLOCATE( eitmp )
      !
      RETURN

   END SUBROUTINE set_evtot_real_x

!------------------------------------------------------------------------------!
   SUBROUTINE set_evtot_twin_x( c0, ctot, lambda, iupdwn_tot, nupdwn_tot )
!------------------------------------------------------------------------------!
      USE kinds,             ONLY: DP
      USE electrons_base,    ONLY: nupdwn, nspin, iupdwn, nudx
      USE electrons_module,  ONLY: nupdwn_emp, n_emp, iupdwn_emp
      USE cp_interfaces,     ONLY: readempty, crot, readempty_twin
      USE cp_main_variables, ONLY: collect_lambda, descla
      USE control_flags,     ONLY: ndw
      USE input_parameters,  ONLY: print_evc0_occ_empty
      USE wavefunctions_module, ONLY: c0_occ_emp_aux
      USE twin_types
      !
      IMPLICIT NONE
      !
      COMPLEX(DP), INTENT(IN)  :: c0(:,:)
      COMPLEX(DP), INTENT(OUT) :: ctot(:,:)
      TYPE(twin_matrix), dimension(:) :: lambda
      INTEGER,     INTENT(IN)  :: iupdwn_tot(2), nupdwn_tot(2)
      !
      COMPLEX(DP), ALLOCATABLE :: cemp(:,:)
      REAL(DP),    ALLOCATABLE :: eitmp(:)
      REAL(DP),    ALLOCATABLE :: lambda_repl(:,:)
      COMPLEX(DP),    ALLOCATABLE :: lambda_repl_c(:,:)
      LOGICAL                  :: t_emp
      !
      ALLOCATE( eitmp( nudx ) )
      !
      ctot = 0.0d0
      !
      IF(.not.lambda(1)%iscmplx) THEN
          ALLOCATE( lambda_repl( nudx, nudx ) )
          CALL collect_lambda( lambda_repl, lambda(1)%rvec(:,:), descla(:,1) )
          CALL crot( ctot, c0, SIZE( c0, 1 ), nupdwn(1), iupdwn_tot(1), iupdwn(1), lambda_repl, nudx, eitmp )
      ELSE
          ALLOCATE( lambda_repl_c( nudx, nudx ) )
          CALL collect_lambda( lambda_repl_c, lambda(1)%cvec(:,:), descla(:,1) )
          CALL crot( ctot, c0, SIZE( c0, 1 ), nupdwn(1), iupdwn_tot(1), iupdwn(1), lambda_repl_c, nudx, eitmp )
      ENDIF
      !
      IF( nspin == 2 ) THEN
          IF(.not.lambda(1)%iscmplx) THEN
              CALL collect_lambda( lambda_repl, lambda(2)%rvec(:,:), descla(:,2) )
              CALL crot( ctot, c0, SIZE( c0, 1 ), nupdwn(2), iupdwn_tot(2), iupdwn(2), lambda_repl, nudx, eitmp )
          ELSE
              CALL collect_lambda( lambda_repl_c, lambda(2)%cvec(:,:), descla(:,2) )
              CALL crot( ctot, c0, SIZE( c0, 1 ), nupdwn(2), iupdwn_tot(2), iupdwn(2), lambda_repl_c, nudx, eitmp )
          ENDIF
      END IF
      !
      IF(.not.lambda(1)%iscmplx) THEN
          DEALLOCATE( lambda_repl )
      ELSE
          DEALLOCATE( lambda_repl_c )
      ENDIF
      !
      t_emp = .FALSE.
      !
      IF( n_emp > 0 ) THEN
          !
          ALLOCATE( cemp( SIZE( c0, 1 ), n_emp * nspin ) )
          cemp = 0.0d0
          t_emp = readempty( cemp, n_emp * nspin, ndw )
          !
      END IF
      !
      IF( t_emp ) THEN
         ctot( :, nupdwn( 1 )+1 : nupdwn_tot(1) ) = cemp( :, 1:nupdwn_emp( 1 ) )
         IF( nspin == 2 ) THEN
             ctot( :, iupdwn_tot(2) + nupdwn(2) : iupdwn_tot(2) + nupdwn_tot(1) - 1 ) = &
                cemp( :, iupdwn_emp(2) : iupdwn_emp(2) + nupdwn_emp(2) - 1 )
         END IF
      END IF
      !
      IF( n_emp > 0 ) DEALLOCATE( cemp )
      !
      DEALLOCATE( eitmp )
      !
      ! print evc0 of occ and empty in xml_io format
      ! 
      IF (print_evc0_occ_empty .and. (n_emp > 0)) THEN
         !
         ALLOCATE( c0_occ_emp_aux( SIZE( c0, 1 ), nupdwn_tot(1) * nspin ) )
         !
         t_emp = .FALSE.
         !
         ALLOCATE( cemp( SIZE( c0, 1 ), n_emp * nspin ) )
         cemp = 0.0d0
         t_emp = readempty_twin( cemp,  n_emp * nspin, ndw )
         !
         IF (t_emp) THEN
            c0_occ_emp_aux(:, iupdwn_tot(1) : nupdwn(1)) = c0(:, iupdwn(1):nupdwn(1))
            c0_occ_emp_aux(:, nupdwn( 1 )+1 : nupdwn_tot(1)) = cemp( :, 1:nupdwn_emp(1))
            IF ( nspin == 2 ) THEN
               c0_occ_emp_aux(:, iupdwn_tot(2) : iupdwn_tot(2) + nupdwn(2) -1 ) = c0(:, iupdwn(2):iupdwn(2)+ nupdwn(2)-1)
               c0_occ_emp_aux(:, iupdwn_tot(2) + nupdwn(2) : iupdwn_tot(2) + nupdwn_tot(2) - 1 ) = &
                                    cemp( :, iupdwn_emp(2) : iupdwn_emp(2) + nupdwn_emp(2) - 1 )
            ENDIF
         ENDIF
         !
         DEALLOCATE( cemp )
         !  
      ENDIF
      !
      RETURN

   END SUBROUTINE set_evtot_twin_x
