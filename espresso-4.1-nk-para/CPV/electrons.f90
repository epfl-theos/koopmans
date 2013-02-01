!
! Copyright (C) 2002-2009 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
  MODULE electrons_module
!=----------------------------------------------------------------------------=!
#include "f_defs.h"
        USE kinds
        USE dspev_module,       ONLY: pdspev_drv, dspev_drv
        USE electrons_base,     ONLY: nbnd, nbndx, nbsp, nbspx, nspin, nel, nelt, &
                                      nupdwn, iupdwn, telectrons_base_initval, f, &
                                      nudx
        USE ensemble_dft,       ONLY: fmat0_diag
        USE cp_electronic_mass, ONLY: ecutmass => emass_cutoff, emass, emass_precond
        use constants,            only : e2, fpi, hartree_si, electronvolt_si

        IMPLICIT NONE
        SAVE

        PRIVATE 

!  ...  declare module-scope variables

        INTEGER, PARAMETER :: nspinx  = 2
        LOGICAL :: band_first = .TRUE.

        INTEGER :: n_emp               =  0  ! number of empty states
        INTEGER :: nudx_emp            =  0  ! maximum number of empty states per spin
        INTEGER :: nupdwn_emp(nspinx)  =  0  ! number of empty states
        INTEGER :: iupdwn_emp(nspinx)  =  0  ! number of empty states

        INTEGER :: nb_l(nspinx)    =  0  ! local number of states ( for each spin components )
        INTEGER :: n_emp_l(nspinx) =  0
        !
        INTEGER  :: max_emp = 0    !  maximum number of iterations for empty states
        REAL(DP) :: ethr_emp, etot_emp, eodd_emp !  threshold for convergence
        !
        INTEGER, ALLOCATABLE :: ib_owner(:)
        INTEGER, ALLOCATABLE :: ib_local(:)

        REAL(DP), ALLOCATABLE :: ei(:,:)
        REAL(DP), ALLOCATABLE :: ei_emp(:,:)
        REAL(DP), ALLOCATABLE :: wfc_centers(:,:,:) !added:giovanni wfc_centers
        REAL(DP), ALLOCATABLE :: wfc_centers_emp(:,:,:) !added:giovanni wfc_centers_emp
        REAL(DP), ALLOCATABLE :: wfc_spreads(:,:,:) !added:giovanni wfc_spreads
        REAL(DP), ALLOCATABLE :: wfc_spreads_emp(:,:,:) !added:giovanni wfc_spreads_emp
        INTEGER, ALLOCATABLE  :: sort_spreads(:,:) !added:giovanni wfc_spreads_emp
        INTEGER, ALLOCATABLE  :: sort_spreads_emp(:,:) !added:giovanni wfc_spreads_emp

        LOGICAL :: icompute_spread=.false. !added:giovanni 

!  ...  Fourier acceleration

        LOGICAL :: toccrd = .FALSE.  ! read occupation number from standard input

        PUBLIC :: electrons_setup
        PUBLIC :: bmeshset, occn_info
        PUBLIC :: deallocate_electrons
        PUBLIC :: n_emp, ei_emp, n_emp_l, ib_owner, ib_local, nb_l
        PUBLIC :: ei, nupdwn_emp, iupdwn_emp, nudx_emp
        PUBLIC :: print_eigenvalues, print_centers_spreads
        PUBLIC :: max_emp, ethr_emp
        PUBLIC :: empty_print_info, empty_init
        PUBLIC :: sort_spreads, sort_spreads_emp,wfc_centers, wfc_spreads, icompute_spread !added:giovanni 
        PUBLIC :: wfc_centers_emp, wfc_spreads_emp !added:giovanni 
        PUBLIC :: etot_emp, eodd_emp !added:giovanni 


!
!  end of module-scope declarations
!
!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!



   SUBROUTINE occn_info( occ )
     !
     !   This subroutine prints occupation numbers to stdout
     !
     USE io_global, ONLY: stdout, ionode
     !
     REAL(DP) :: occ(:)
     INTEGER  :: i, iss
     !
     IF( ionode ) THEN
       WRITE( stdout, fmt="(3X,'Occupation number from init')" )
       IF( nspin == 1 ) THEN
         WRITE( stdout, fmt = " (3X, 'nbnd = ', I5 ) " ) nbnd
         WRITE( stdout, fmt = " (3X,10F5.2)" ) ( occ( i ), i = 1, nbnd )
       ELSE
         DO iss = 1, nspin
           WRITE( stdout, fmt = " (3X,'spin = ', I3, ' nbnd = ', I5 ) " ) iss, nupdwn( iss )
           WRITE( stdout, fmt = " (3X,10F5.2)" ) ( occ( i+iupdwn(iss)-1 ), i = 1, nupdwn( iss ) )
         END DO
       END IF
     END IF
     !
     RETURN
   END SUBROUTINE occn_info

!  ----------------------------------------------
!  ----------------------------------------------

   SUBROUTINE bmeshset

     !   This subroutine initialize the variables for the 
     !   distribution across processors of the overlap matrixes 
     !   of sizes ( nx, nx )

     USE mp_global, ONLY: me_image, nproc_image

     IMPLICIT NONE

     INTEGER :: i, ierr

     IF( band_first ) THEN
       CALL errore(' bmeshset ',' module not initialized ',0)
     END IF

     DO i = 1, nspin 
       !
       IF( i > nspinx ) CALL errore( ' bmeshset ',' spin too large ', i)
       !
       nb_l( i ) = nupdwn( i ) / nproc_image
       IF( me_image < MOD( nupdwn( i ), nproc_image ) ) nb_l( i ) = nb_l( i ) + 1
       !
       n_emp_l( i ) = nupdwn_emp( i ) / nproc_image
       IF( me_image < MOD( nupdwn_emp( i ), nproc_image ) ) n_emp_l( i ) = n_emp_l( i ) + 1
       !
     END DO

     IF( ALLOCATED( ib_owner ) ) DEALLOCATE( ib_owner )
     ALLOCATE( ib_owner( MAX( n_emp, nbndx ) ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' bmeshset ',' allocating ib_owner ', ierr)
     IF( ALLOCATED( ib_local ) ) DEALLOCATE( ib_local )
     ALLOCATE( ib_local( MAX( n_emp, nbndx ) ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' bmeshset ',' allocating ib_local ', ierr)

     !  here define the association between processors and electronic states
     !  round robin distribution is used

     ib_local =  0
     ib_owner = -1
     DO i = 1, MAX( n_emp, nbndx )
       ib_local( i ) = ( i - 1 ) / nproc_image        !  local index of the i-th band 
       ib_owner( i ) = MOD( ( i - 1 ), nproc_image )  !  owner of th i-th band
       IF( me_image <= ib_owner( i ) ) THEN
         ib_local( i ) = ib_local( i ) + 1
       END IF
     END DO

     RETURN
   END SUBROUTINE bmeshset

!  ----------------------------------------------
!
!
!
!  ----------------------------------------------


   SUBROUTINE electrons_setup( n_emp_ , emass_inp, ecutmass_inp )

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: n_emp_
     REAL(DP),  INTENT(IN) :: emass_inp, ecutmass_inp
     INTEGER :: ierr, i,j
 

     IF( .NOT. telectrons_base_initval ) &
       CALL errore( ' electrons_setup ', ' electrons_base not initialized ', 1 )

     n_emp = n_emp_
     !
     nupdwn_emp(1) = n_emp
     iupdwn_emp(1) = 1
     nudx_emp=n_emp

     IF( nspin == 2 ) THEN
        nupdwn_emp(2) = n_emp
        iupdwn_emp(2) = 1 + n_emp
     END IF

     IF( ALLOCATED( ei ) ) DEALLOCATE( ei )
     ALLOCATE( ei( nudx, nspin ), STAT=ierr)
     IF( ierr/=0 ) CALL errore( ' electrons ',' allocating ei ',ierr)
     ei = 0.0_DP

     IF( ALLOCATED( ei_emp ) ) DEALLOCATE( ei_emp )
     IF( n_emp > 0 ) THEN
       ALLOCATE( ei_emp( n_emp, nspin ), STAT=ierr)
       IF( ierr/=0 ) CALL errore( ' electrons ',' allocating ei_emp ',ierr)
       ei_emp = 0.0_DP
     END IF

!begin_added:giovanni
     IF( ALLOCATED( sort_spreads ) ) DEALLOCATE( sort_spreads )
     IF(nudx > 0) THEN
        ALLOCATE( sort_spreads(nudx,nspin ), STAT=ierr)
        IF( ierr/=0 ) CALL errore( ' electrons ',' allocating sort_spreads ',ierr)
        do j=1,nspin
           do i=1,size(sort_spreads(:,j))
              sort_spreads(i,j)=i
           enddo
        enddo
     ENDIF

     IF( ALLOCATED( sort_spreads_emp ) ) DEALLOCATE( sort_spreads_emp )
     IF(nudx_emp > 0) THEN
        ALLOCATE( sort_spreads_emp(nudx_emp,nspin ), STAT=ierr)
        IF( ierr/=0 ) CALL errore( ' electrons ',' allocating sort_spreads_emp ',ierr)
        do j=1,nspin
           do i=1,size(sort_spreads_emp(:,j))
              sort_spreads_emp(i,j)=i
           enddo
        enddo
     ENDIF
     
     IF( ALLOCATED( wfc_centers ) ) DEALLOCATE( wfc_centers )
     IF(nudx > 0) THEN
        ALLOCATE( wfc_centers(4,nudx,nspin ), STAT=ierr)
        IF( ierr/=0 ) CALL errore( ' electrons ',' allocating wfc_centers ',ierr)
        wfc_centers = 0.0_DP
     ENDIF

     IF( ALLOCATED( wfc_spreads ) ) DEALLOCATE( wfc_spreads )
     IF(nudx > 0) THEN
        ALLOCATE( wfc_spreads( nudx, nspin, 2 ), STAT=ierr)
        IF( ierr/=0 ) CALL errore( ' electrons ',' allocating wfc_spreads ',ierr)
        wfc_spreads = 0.0_DP
     ENDIF
!end_added:giovanni

!begin_added:giovanni
     IF( ALLOCATED( wfc_centers_emp ) ) DEALLOCATE( wfc_centers_emp )
     IF( nudx_emp > 0 ) THEN
        ALLOCATE( wfc_centers_emp(4, nudx_emp,nspin ), STAT=ierr)
        IF( ierr/=0 ) CALL errore( ' electrons ',' allocating wfc_centers_emp ',ierr)
        wfc_centers_emp = 0.0_DP
     ENDIF

     IF( ALLOCATED( wfc_spreads_emp ) ) DEALLOCATE( wfc_spreads_emp )
     IF( nudx_emp > 0 ) THEN
        ALLOCATE( wfc_spreads_emp( nudx_emp, nspin, 2 ), STAT=ierr)
        IF( ierr/=0 ) CALL errore( ' electrons ',' allocating wfc_spreads_emp ',ierr)
        wfc_spreads_emp = 0.0_DP
     ENDIF
!end_added:giovanni

     ecutmass = ecutmass_inp
     emass    = emass_inp
     IF ( ecutmass < 0.0_DP ) &
       CALL errore(' electrons ',' ecutmass out of range ' , 0)

     band_first = .FALSE.

     RETURN
   END SUBROUTINE electrons_setup

!----------------------------------------------------------------------

        SUBROUTINE empty_print_info(iunit)
          !
          USE kinds,            ONLY: DP
          INTEGER, INTENT(IN) :: iunit
          !
          IF ( n_emp > 0 ) THEN
              WRITE (iunit,"(//, 3X, 'Empty states minimization')")
              WRITE (iunit,"(    3X, '--------------------------')")
              WRITE (iunit,"(3X, '   states = ',i8)") n_emp
              WRITE (iunit,"(3X, '  maxiter = ',i8)") max_emp
              WRITE (iunit,"(3X, '     ethr = ',D12.4)") ethr_emp
          ENDIF
          !
          RETURN
        END SUBROUTINE empty_print_info

!----------------------------------------------------------------------

        SUBROUTINE empty_init( max_emp_ , ethr_emp_ )

          USE kinds,            ONLY: DP

          INTEGER, INTENT(IN) :: max_emp_
          REAL(DP), INTENT(IN) :: ethr_emp_

          max_emp   = max_emp_
          ethr_emp  = ethr_emp_
          !
          RETURN
        END SUBROUTINE empty_init


!  ----------------------------------------------


   SUBROUTINE print_eigenvalues( ei_unit, tfile, tstdout, nfi, tps )
      !
      use constants,      only : autoev 
      USE io_global,      ONLY : stdout, ionode
      USE ensemble_dft,   ONLY : tens, tsmear
      !
      INTEGER,  INTENT(IN) :: ei_unit
      LOGICAL,  INTENT(IN) :: tfile, tstdout
      INTEGER,  INTENT(IN) :: nfi
      REAL(DP), INTENT(IN) :: tps
      !
      INTEGER :: i, j, ik
      !
      IF ( tfile ) THEN
          WRITE(ei_unit,30) nfi, tps
      END IF
      !
      ik = 1
      !
      IF(tstdout) THEN
         IF(nspin==1) THEN
            WRITE( stdout,1101) 
            WRITE( stdout, 1444) MAXVAL(ei(1:nupdwn(1),1)*autoev, nupdwn(1))
         ELSE
            WRITE( stdout,1101) 
            WRITE( stdout, 1444) MAX(MAXVAL(ei(1:nupdwn(1),1)*autoev, nupdwn(1)), MAXVAL(ei(1:nupdwn(2),2)*autoev, nupdwn(2)))
         ENDIF
      
         IF(n_emp.gt.0) THEN 
            IF(nspin==1) THEN
               WRITE( stdout,1201) 
               WRITE( stdout, 1444) MINVAL(ei_emp(1:n_emp,1)*autoev, n_emp)
            ELSE
               WRITE( stdout,1201) 
               WRITE( stdout, 1444) MIN(MINVAL(ei_emp(1:n_emp,1)*autoev, n_emp), MINVAL(ei_emp(1:n_emp,2)*autoev, n_emp))
            ENDIF
         ENDIF
      ENDIF

      DO j = 1, nspin
         !
         IF( tstdout ) THEN
            WRITE( stdout,1002) ik, j
            WRITE( stdout,1004) ( ei( i, j ) * autoev, i = 1, nupdwn(j) )

            !
            IF ( tens .OR. tsmear ) THEN
                WRITE( stdout,1082) ik, j
                WRITE( stdout,1084) ( f( i ), i = iupdwn(j), iupdwn(j)+nupdwn(j)-1 )
                !
                WRITE( stdout,1092) ik, j
                WRITE( stdout,1084) ( fmat0_diag( i ), i = iupdwn(j), iupdwn(j)+nupdwn(j)-1 )
            ENDIF
            !
            IF( n_emp .GT. 0 ) THEN
               WRITE( stdout,1005) ik, j
               WRITE( stdout,1004) ( ei_emp( i, j ) * autoev , i = 1, n_emp )
               IF(nupdwn(j)>0) &
               WRITE( stdout,1006) ( ei_emp( 1, j ) - ei( nupdwn(j), j ) ) * autoev
            END IF
         END IF
         !
         IF( tfile ) THEN
            WRITE(ei_unit,1010) ik, j
            WRITE(ei_unit,1020) ( ei( i, j ) * autoev, i = 1, nupdwn(j) )
            IF( n_emp .GT. 0 ) THEN
               WRITE(ei_unit,1011) ik, j
               WRITE(ei_unit,1020) ( ei_emp( i, j ) * autoev , i = 1, n_emp )
               IF(nupdwn(j)>0) &
                WRITE(ei_unit,1021) (( ei_emp( 1, j ) - ei( nupdwn(j), j ) ) * autoev)
            END IF
         END IF
         !
      END DO
      !
  30  FORMAT(2X,'STEP:',I7,1X,F10.2)
 1002 FORMAT(/,3X,'Eigenvalues (eV), kp = ',I3, ' , spin = ',I2,/)
 1101 FORMAT(/,3X,'HOMO Eigenvalue (eV)',/) !added_giovanni
 1022 FORMAT(/,3X,'Centers (Bohr), kp = ',I3, ' , spin = ',I2,/) !added_giovanni
 1222 FORMAT(/,3X,'Spreads (Bohr^2), kp = ',I3, ' , spin = ',I2,/) !added_giovanni
 1082 FORMAT(/,3X,'Occupations,      kp = ',I3, ' , spin = ',I2,/)
 1092 FORMAT(/,3X,'DensityMat diag,  kp = ',I3, ' , spin = ',I2,/)
 1005 FORMAT(/,3X,'Empty States Eigenvalues (eV), kp = ',I3, ' , spin = ',I2,/)
 1201 FORMAT(/,3X,'LUMO Eigenvalue (eV)',/) !added_giovanni
 1004 FORMAT(10F8.2)
 1044 FORMAT(4F8.2)
 1444 FORMAT(1F8.2)
 1084 FORMAT(10F8.4)
 1006 FORMAT(/,3X,'Electronic Gap (eV) = ',F8.2,/)
 1010 FORMAT(3X,'Eigenvalues (eV), kp = ',I3, ' , spin = ',I2)
 1011 FORMAT(3X,'Empty States Eigenvalues (eV), kp = ',I3, ' , spin = ',I2)
 1020 FORMAT(10F8.2)
 1021 FORMAT(3X,'Electronic Gap (eV) = ',F8.2)
 1030 FORMAT(3X,'nfill = ', I4, ', nempt = ', I4, ', kp = ', I3, ', spin = ',I2)
      !
      RETURN
   END SUBROUTINE print_eigenvalues

   SUBROUTINE print_centers_spreads( spread_unit, tfile, tstdout, nfi, tps )
      !
      use constants,      only : autoev 
      USE io_global,      ONLY : stdout, ionode
      USE ensemble_dft,   ONLY : tens, tsmear
      use nksic,          ONLY : complexification_index, pink, pink_emp, do_orbdep, pzalpha=> odd_alpha, pzalpha_emp=>odd_alpha_emp
      !
      INTEGER,  INTENT(IN) :: spread_unit
      LOGICAL,  INTENT(IN) :: tfile, tstdout
      INTEGER,  INTENT(IN) :: nfi
      REAL(DP), INTENT(IN) :: tps
      !
      INTEGER :: i, j, ik
      !
      !
      ik = 1
      !
      DO j = 1, nspin
         !
         IF( tstdout ) THEN
            
            WRITE( stdout,1222) ik, j
            !
            IF(do_orbdep) THEN
               WRITE( stdout,1444) ( wfc_centers(1:4, i, j ),  wfc_spreads( i, j , 1), wfc_spreads( i, j, 2), pink(iupdwn(j)-1+sort_spreads(i,j))*hartree_si/electronvolt_si, pzalpha(iupdwn(j)-1+sort_spreads(i,j)), i = 1, nupdwn(j) )
            ELSE
               WRITE( stdout,1445) ( wfc_centers(1:4, i, j ),  wfc_spreads( i, j , 1), wfc_spreads( i, j, 2), i = 1, nupdwn(j) )
            ENDIF
            !
            IF( n_emp .GT. 0 ) THEN
               WRITE( stdout,12224) ik, j
!               WRITE( stdout,1005) ik, j
               IF(do_orbdep) THEN
                  !
                  WRITE( stdout,1444) ( wfc_centers_emp(1:4, i, j ),  wfc_spreads_emp( i, j , 1), wfc_spreads_emp( i, j, 2), pink_emp(iupdwn_emp(j)-1+sort_spreads_emp(i,j))*hartree_si/electronvolt_si, i = 1, nupdwn_emp(j) )
                  !
               ELSE
                  WRITE( stdout,1445) ( wfc_centers(1:4, i, j ),  wfc_spreads( i, j , 1), wfc_spreads( i, j, 2), i = 1, nupdwn(j) )
               ENDIF
               
            END IF
         END IF
         !
         !
      END DO
      !
      IF(tstdout) THEN
         !
         WRITE(stdout, *) " Complexification index "
         WRITE(stdout, *) complexification_index
         !
      ENDIF
      !
  30  FORMAT(2X,'STEP:',I7,1X,F10.2)
 1022 FORMAT(/,3X,'Centers (Bohr), kp = ',I3, ' , spin = ',I2,/)
 1222 FORMAT(/,3X,'Charge  ---   Centers xyz (Bohr)  ---  Spreads (Bohr^2) - SH(eV), kp = ',I3, ' , spin = ',I2,/)
 12224 FORMAT(/,3X,'Empty Charge  ---   Centers xyz (Bohr)  ---  Spreads (Bohr^2) - SH(eV), kp = ',I3, ' , spin = ',I2,/)
 1005 FORMAT(/,3X,'Empty States Eigenvalues (eV), kp = ',I3, ' , spin = ',I2,/)
 1444 FORMAT(F8.2,'   ---',3F8.2,'   ---',4F8.3)
 1445 FORMAT(F8.2,'   ---',3F8.2,'   ---',2F8.2)
 1121 FORMAT(/3X,'Manifold complexification index = ',2F8.4/)
 1084 FORMAT(10F8.4)
      !
      RETURN
   END SUBROUTINE print_centers_spreads


!  ----------------------------------------------

   SUBROUTINE deallocate_electrons
      INTEGER :: ierr
      IF(ALLOCATED(ei))       THEN
            DEALLOCATE(ei, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating ei ',ierr )
      END IF
      IF(ALLOCATED(ei_emp))   THEN
            DEALLOCATE(ei_emp, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating ei_emp ',ierr )
      END IF
      IF(ALLOCATED(ib_owner)) THEN
            DEALLOCATE(ib_owner, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating ib_owner ',ierr )
      END IF
      IF(ALLOCATED(ib_local)) THEN
            DEALLOCATE(ib_local, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating ib_local ',ierr )
      END IF
      IF(ALLOCATED(sort_spreads))       THEN
            DEALLOCATE(sort_spreads, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating sort_spreads ',ierr )
      END IF
      IF(ALLOCATED(sort_spreads_emp))       THEN
            DEALLOCATE(sort_spreads_emp, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating sort_spreads ',ierr )
      END IF
      IF(ALLOCATED(wfc_centers))       THEN
            DEALLOCATE(wfc_centers, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating wfc_centers ',ierr )
      END IF
      IF(ALLOCATED(wfc_spreads))       THEN
            DEALLOCATE(wfc_spreads, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating wfc_spreads ',ierr )
      END IF
!       write(6,*) "deallocating empty", ubound(wfc_centers_emp), ubound(wfc_spreads_emp)
      IF(ALLOCATED(wfc_centers_emp))       THEN
!       write(6,*) "deallocating wfc_centers_emp"
            DEALLOCATE(wfc_centers_emp, STAT=ierr)
!       write(6,*) "deallocated wfc_centers_emp"
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating wfc_centers_emp ',ierr )
      END IF
      IF(ALLOCATED(wfc_spreads_emp))       THEN
            DEALLOCATE(wfc_spreads_emp, STAT=ierr)
            IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating wfc_spreads_emp ',ierr )
      END IF
!             write(6,*) "deallocated empty"

      RETURN
   END SUBROUTINE deallocate_electrons
        



!=----------------------------------------------------------------------------=!
  END MODULE electrons_module
!=----------------------------------------------------------------------------=!
