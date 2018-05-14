!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! Written and Revised by Carlo Cavazzoni

!=----------------------------------------------------------------------------------=!

   SUBROUTINE runcp_uspp_x &
      ( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, fromscra, restart, tprint_ham )
      !
      !  This subroutine performs a Car-Parrinello or Steepest-Descent step
      !  on the electronic variables, computing forces on electrons
      ! 
      !  on input:
      !  c0  wave functions at time t
      !  cm  wave functions at time t - dt 
      !
      !  on output:
      !  cm  wave functions at time t + dt, not yet othogonalized 
      !
      USE parallel_include
      USE kinds,               ONLY : DP
      USE mp_global,           ONLY : nogrp, ogrp_comm, me_image, nolist
      USE mp,                  ONLY : mp_sum
      USE fft_base,            ONLY : dffts
      use wave_base,           only : wave_steepest, wave_verlet
      use control_flags,       only : lwf, tsde, use_task_groups, program_name, &
                                                       gamma_only, do_wf_cmplx !added:giovanni do_wf_cmplx
      use uspp,                only : deeq, vkb
      use reciprocal_vectors,  only : gstart
      use electrons_base,      only : n=>nbsp, ispin, f, nspin, nupdwn, iupdwn
      USE electrons_base,      ONLY:  nx=>nbspx
      use wannier_subroutines, only : ef_potential
      use efield_module,       only : dforce_efield, tefield, dforce_efield2, tefield2
      use gvecw,               only : ngw, ngwx
      USE cp_main_variables,   ONLY : hamilt, iprint_stdout
      USE cp_interfaces,       ONLY : dforce
      USE task_groups,         ONLY : tg_gather
      USE ldaU,                ONLY : vupsi, lda_plus_u
      use hfmod,               only : do_hf, vxxpsi
      use nksic,               only : do_orbdep, vsic, vsicpsi, deeq_sic, & 
                                      f_cutoff, valpsi
      use ensemble_dft,        only : tens, tsmear
      use twin_types !added:giovanni
      use input_parameters,    only : odd_nkscalfact 
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nfi
      REAL(DP) :: fccc, ccc
      REAL(DP) :: ema0bg(:), dt2bye
      REAL(DP) :: rhos(:,:)
      type(twin_matrix) :: bec!(:,:) !modified:giovanni
      COMPLEX(DP) :: c0(:,:), cm(:,:)
      LOGICAL, OPTIONAL, INTENT(IN) :: fromscra
      LOGICAL, OPTIONAL, INTENT(IN) :: restart
      LOGICAL, OPTIONAL, INTENT(IN) :: tprint_ham
      !
      !
      real(DP) ::  verl1, verl2, verl3
      real(DP),    allocatable :: emadt2(:)
      real(DP),    allocatable :: emaver(:)
      real(DP),    allocatable :: faux(:)
      complex(DP), allocatable :: c2(:), c3(:)
      REAL(DP),    ALLOCATABLE :: tg_rhos(:,:)
      integer :: nsiz, incr, idx, idx_in, ierr
      integer :: iwfc, nwfc, is, tg_rhos_siz, c2_siz
      integer :: isp, j, jj, j0, i, ii, i0
      integer :: iflag
      logical :: ttsde
      LOGICAL :: tprint_ham_ 
      LOGICAL :: lgam
      integer :: iss

     lgam=gamma_only.and..not.do_wf_cmplx

     allocate(faux(nx))

     iflag = 0
     !
     IF( PRESENT( fromscra ) ) THEN
         IF( fromscra ) iflag = 1
     ENDIF
     !
     IF( PRESENT( restart ) ) THEN
         IF( restart ) iflag = 2
     ENDIF
     !
     tprint_ham_ = .FALSE.

     DO iss=1,size(hamilt)
      if(.not. hamilt(iss)%iscmplx) then
	hamilt(iss)%rvec = 0.0d0
      else
	hamilt(iss)%cvec = CMPLX(0.0d0,0.d0)
      endif
     END DO

     IF ( PRESENT( tprint_ham ) ) THEN
         !
         tprint_ham_ = tprint_ham
     ENDIF
     !
     IF( use_task_groups ) THEN
         tg_rhos_siz = nogrp * dffts%nnrx
         c2_siz      = nogrp * ngwx
     ELSE
         tg_rhos_siz = 1
         c2_siz      = ngw 
     ENDIF

     !
     ! ...  set verlet variables 
     !
     verl1 = 2.0d0 * fccc
     verl2 = 1.0d0 - verl1
     verl3 = 1.0d0 * fccc
     ALLOCATE( emadt2( ngw ) )
     ALLOCATE( emaver( ngw ) )

     ccc    = fccc * dt2bye
     emadt2 = dt2bye * ema0bg
     emaver = emadt2 * verl3

     IF( iflag == 0 ) THEN
       ttsde  = tsde
     ELSE IF( iflag == 1 ) THEN
       ttsde = .TRUE.
     ELSE IF( iflag == 2 ) THEN
       ttsde = .FALSE.
     END IF

     IF( lwf ) THEN

        call ef_potential( nfi, rhos, bec%rvec, deeq, vkb, c0, cm, emadt2, emaver, verl1, verl2 ) !warning:giovanni not yet modified

     ELSE

        allocate( c2( c2_siz ), c3( c2_siz ) )
        allocate( tg_rhos( tg_rhos_siz, nspin ) )
        !
        c2      = 0D0
        c3      = 0D0
       
        ! size check
        IF ( do_orbdep .AND. use_task_groups ) &
            CALL errore('runcp_uspp','NKSIC and task_groups incompatible',10)

        IF( use_task_groups ) THEN
           !
           !  The potential in rhos is distributed accros all processors
           !  We need to redistribute it so that it is completely contained in the
           !  processors of an orbital TASK-GROUP
           !
           DO i = 1, nspin
              CALL tg_gather( dffts, rhos(:,i), tg_rhos(:,i) )
           END DO
           
!            IF(lgam) THEN
             incr = 2 * nogrp
!            ELSE
!              incr = nogrp
!            ENDIF

        ELSE

!            IF(lgam) THEN
             incr = 2 
!            ELSE
!              incr = 1
!            ENDIF

        END IF


        !
        ! faux takes into account spin multiplicity.
        !
        faux(:) = f(:) * DBLE( nspin ) / 2.0d0
        !
        DO j = 1, n
            faux(j) = max( f_cutoff, faux(j) )
        ENDDO
        

        DO i = 1, n, incr
           !
           ! ... apply preconditioning on occupation
           ! ... (this preconditioning must be taken into account when
           ! ... calculating eigenvalues in eigs0.f90)
           !

           !
           IF( use_task_groups ) THEN
              !
              !The input coefficients to dforce cover eigenstates i:i+2*NOGRP-1
              !Thus, in dforce the dummy arguments for c0(1,i) and
              !c0(1,i+1) hold coefficients for eigenstates i,i+2*NOGRP-2,2
              !and i+1,i+2*NOGRP...for example if NOGRP is 4 then we would have
              !1-3-5-7 and 2-4-6-8
              !

              if( tefield .OR. tefield2 ) then
                 CALL errore( ' runcp_uspp ', ' electric field with task group not implemented yet ', 1 )
              end if
              if( lda_plus_u ) then
                 CALL errore( ' runcp_uspp ', ' lda_plus_u with task group not implemented yet ', 1 )
              end if

              CALL dforce( i, bec, vkb, c0, c2, c3, tg_rhos, tg_rhos_siz, ispin, faux, n, nspin )

           ELSE
!begin_added:giovanni:debug ------------ FORCES
!               write(6,*) "c0, debug, before dforce", i
!               write(6,*) c0(1,i), c0(2,i), c0(3,i), rhos(1,1)
!               write(6,*) "c2, debug, before dforce"
!               write(6,*) c2(1), c2(2), c2(3), rhos(2,1)
!               write(6,*) "c3, debug, before dforce"
!               write(6,*) c3(1), c3(2), c3(3), rhos(3,1)
!end_added:giovanni:debug ------------ FORCES
              CALL dforce( i, bec, vkb, c0, c2, c3, rhos, SIZE(rhos,1), ispin, faux, n, nspin )
!begin_added:giovanni:debug ------------ FORCES
!               write(6,*) "c0, debug, after dforce"
!               write(6,*) c0(1,i), c0(2,i), c0(3,i), rhos(1,1)
!               write(6,*) "c2, debug, after dforce"
!               write(6,*) c2(1), c2(2), c2(3), rhos(2,1)
!               write(6,*) "c3, debug, after dforce"
!               write(6,*) c3(1), c3(2), c3(3), rhos(3,1)
!end_added:giovanni:debug ------------ FORCES
           END IF


           IF ( lda_plus_u ) THEN
               !
               IF ( tens .OR. tsmear ) THEN
                   !
                   c2(:) = c2(:) - vupsi(:,i) 
                   c3(:) = c3(:) - vupsi(:,i+1) 
                   !
               ELSE
                   !
                   c2(:) = c2(:) - vupsi(:,i) * faux(i)
                   c3(:) = c3(:) - vupsi(:,i+1) * faux(i+1)
                   !
               ENDIF
               !
           ENDIF
           !
           IF ( do_orbdep ) THEN
               !
               IF ( odd_nkscalfact ) THEN
                   !
                   IF ( tens .OR. tsmear ) THEN
                      !
                      c2(:) = c2(:) - valpsi(i, :)
                      c3(:) = c3(:) - valpsi(i+1, :)
                      !
                   ELSE
                      !
                      c2(:) = c2(:) - valpsi(i,:)   * faux(i)
                      c3(:) = c3(:) - valpsi(i+1, :) * faux(i+1)
                      !
                   ENDIF
                   !
               ENDIF
               !
               ! faux takes into account spin multiplicity.
               !
               CALL nksic_eforce( i, n, nx, vsic, deeq_sic, bec, ngw, c0(:,i), c0(:,i+1), vsicpsi, lgam )
               !
               IF ( tens .OR. tsmear ) THEN
                   !
                   c2(:) = c2(:) - vsicpsi(:,1)
                   c3(:) = c3(:) - vsicpsi(:,2)
                   !
               ELSE
                   !
                   c2(:) = c2(:) - vsicpsi(:,1) * faux(i)   
                   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1) 
                   !
               ENDIF

               !
               ! build the matrix elements of 
               ! the  Delta h^SIC hamiltonian
               !
               ! (for the sake of plotting the evolution 
               !  of its imaginary part )
               !
               IF ( tprint_ham_ ) THEN
                   !
                   DO ii = 0, 1
                       ! 
                       IF ( i+ii > n ) CYCLE
                       !
                       isp = ispin( i+ii )
                       i0  = i+ii
                       !
                       IF ( nspin==2 ) i0 = i0 -iupdwn(isp) +1 
                       ! NOTE: iupdwn(isp) is set to zero if nspin=1
                       !
                       DO j0 = 1, nupdwn(isp)
                           !
                           jj = j0
                           IF ( nspin==2 ) jj = jj +iupdwn(isp) -1 
                           !
                           IF(.not.hamilt(isp)%iscmplx) THEN
			      hamilt(isp)%rvec(j0, i0) = 2.0d0 * DOT_PRODUCT( c0(:,jj), vsicpsi(:,ii+1) )
                              !
			      IF ( gstart == 2 ) THEN
				  hamilt(isp)%rvec( j0, i0) =  hamilt(isp)%rvec( j0, i0) -c0(1,jj)*vsicpsi(1,ii+1) 
			      ENDIF
                           ELSE
			      hamilt(isp)%cvec(j0, i0) = DOT_PRODUCT( c0(:,jj), vsicpsi(:,ii+1) ) !warning:giovanni put conjugate??
                           ENDIF
                           !
                       ENDDO
                       !
                   ENDDO
                   !
               ENDIF
               !
           ENDIF
           
           !
           ! HF exchange 
           ! 
           IF ( do_hf ) THEN
               !
               IF ( tens .OR. tsmear ) THEN
                   !
                   c2(:) = c2(:) - vxxpsi(:,i)
                   c3(:) = c3(:) - vxxpsi(:,i+1)
                   !
               ELSE
                   !
                   c2(:) = c2(:) - vxxpsi(:,i)   * faux(i)
                   c3(:) = c3(:) - vxxpsi(:,i+1) * faux(i+1)
                   !
               ENDIF
               !
           ENDIF
           
           !
           ! spin multiplicity and occupation factors 
           ! are taken into account inside the calls
           !
           IF( tefield ) THEN
               CALL dforce_efield ( bec%rvec, i, c0, c2, c3, rhos)
           ENDIF
           !
           IF( tefield2 ) THEN
               CALL dforce_efield2 ( bec%rvec, i, c0, c2, c3, rhos)
           ENDIF

           IF( iflag == 2 ) THEN
              DO idx = 1, incr, 2
                 IF( i + idx - 1 <= n ) THEN
                    cm( :, i+idx-1) = c0(:,i+idx-1)
                    cm( :, i+idx  ) = c0(:,i+idx  )
                 END IF
              ENDDO
           END IF
           ! 
           idx_in = 1
           DO idx = 1, incr, 2
              IF( i + idx - 1 <= n ) THEN
                 IF (tsde) THEN 
                    CALL wave_steepest( cm(:, i+idx-1 ), c0(:, i+idx-1 ), emaver, c2, ngw, idx_in )
                    CALL wave_steepest( cm(:, i+idx   ), c0(:, i+idx   ), emaver, c3, ngw, idx_in )
                 ELSE
                    CALL wave_verlet( cm(:, i+idx-1 ), c0(:, i+idx-1 ), verl1, verl2, emaver, c2, ngw, idx_in )
                    CALL wave_verlet( cm(:, i+idx   ), c0(:, i+idx   ), verl1, verl2, emaver, c3, ngw, idx_in )
                 ENDIF
                 IF ( gstart == 2 ) THEN
                    IF(lgam) THEN
                       cm(1,i+idx-1) = CMPLX(DBLE(cm(1,i+idx-1)),0.0d0)
                       cm(1,i+idx  ) = CMPLX(DBLE(cm(1,i+idx  )),0.0d0)
                    ENDIF
                 END IF
              END IF
              !
              idx_in = idx_in + 1
              !
           END DO
           ! 
        end do

        DEALLOCATE( c2 )
        DEALLOCATE( c3 )
        DEALLOCATE( tg_rhos )

     END IF
     !
     IF ( tprint_ham_ ) THEN
       DO iss=1, size(hamilt)
         IF(.not.hamilt(iss)%iscmplx) THEN
	    CALL mp_sum( hamilt(iss)%rvec )
         ELSE
	    CALL mp_sum( hamilt(iss)%cvec )
         ENDIF
       END DO
     ENDIF
     !
     DEALLOCATE( emadt2 )
     DEALLOCATE( emaver )
!
   END SUBROUTINE runcp_uspp_x
!
!
!=----------------------------------------------------------------------------=!
!
!

!=----------------------------------------------------------------------------=!

    SUBROUTINE runcp_uspp_force_pairing_x  & !warning:giovanni still to be modified
       ( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, intermed, fromscra, &
         restart )
  !
!  same as runcp, except that electrons are paired forcedly
!  i.e. this handles a state dependant Hamiltonian for the paired and unpaired electrons
!  unpaired is assumed to exist, to be unique, and located in highest index band

      USE kinds,               ONLY : DP
      USE wave_base,           ONLY : wave_steepest, wave_verlet
      USE control_flags,       ONLY : lwf, tsde, program_name, use_task_groups
      USE uspp,                ONLY : deeq, vkb
      USE reciprocal_vectors,  ONLY : gstart
      USE wannier_subroutines, ONLY : ef_potential
      USE efield_module,       ONLY : dforce_efield, tefield
      USE electrons_base,      ONLY : ispin, nspin, f, n=>nbsp
      USE cp_interfaces,       ONLY : dforce
  !
      USE gvecw, ONLY: ngw
  !
  !
      USE electrons_base,   ONLY: nx=>nbnd, nupdwn, iupdwn, nbspx, nbsp
      USE mp, ONLY: mp_sum 
      USE mp_global, ONLY: intra_image_comm 
!@@@@
      USE ldaU
!@@@@
  !
      IMPLICIT NONE
      INTEGER, INTENT(in) :: nfi
      REAL(DP) :: fccc, ccc
      REAL(DP) :: ema0bg(:), dt2bye
      REAL(DP) :: rhos(:,:)
      REAL(DP) :: bec(:,:)
      COMPLEX(DP) :: c0(:,:), cm(:,:)
      REAL(DP)    :: intermed
      LOGICAL, OPTIONAL, INTENT(in) :: fromscra
      LOGICAL, OPTIONAL, INTENT(in) :: restart
!
      REAL(DP) ::  verl1, verl2, verl3
      REAL(DP), ALLOCATABLE:: emadt2(:)
      REAL(DP), ALLOCATABLE:: emaver(:)
      COMPLEX(DP), ALLOCATABLE:: c2(:), c3(:)
      INTEGER :: i
      INTEGER :: iflag
      LOGICAL :: ttsde
!
       INTEGER     :: ierr,  nb, np_dw, is_dw, npair, n_unp, n_dwn, n_pair 
       REAL(DP)    :: ei_unp_mem, ei_unp_wfc
       COMPLEX(DP) :: intermed3
       REAL(DP),    ALLOCATABLE :: occ(:)
       COMPLEX(DP), ALLOCATABLE :: c4(:), c5(:)
!
! ... Controlling on sic applicability
!
       IF( lwf ) CALL errore('runcp_uspp_force_pairing', &
                           'Wannier function and sic are not compatibile',1)
       IF( tefield ) CALL errore('runcp_uspp_force_pairing', &
                           'Electric field and sic are not implemented',2)
       IF( nspin == 1 ) CALL errore(' runcp_force_pairing ',' inconsistent nspin ', 1)

       IF( use_task_groups ) CALL errore(' runcp_force_pairing ',' task_groups not implemented ', 1)
!       
       ALLOCATE( emadt2( ngw ) )
       ALLOCATE( emaver( ngw ) )      
!
       iflag = 0
       IF( PRESENT( fromscra ) ) THEN
          IF( fromscra ) iflag = 1
       END IF
       IF( PRESENT( restart ) ) THEN
          IF( restart ) iflag = 2
       END IF
!       
       IF( iflag == 0 ) THEN
          ttsde  = tsde
       ELSE IF( iflag == 1 ) THEN
          ttsde = .TRUE.
       ELSE IF( iflag == 2 ) THEN
          ttsde = .FALSE.
       END IF
!
       ALLOCATE( c2(ngw), c3(ngw), c4(ngw), c5(ngw) )
       !
       ! ...  set verlet variables
       !
       verl1 = 2.0d0 * fccc
       verl2 = 1.0d0 - verl1
       verl3 = 1.0d0 * fccc 
!
       ccc    = fccc * dt2bye
       emadt2 = dt2bye * ema0bg
       emaver = emadt2 * verl3
!
       IF( nupdwn(1) /= (nupdwn(2) + 1) ) &
          CALL errore(' runcp_force_pairing ',' inconsistent number of states ', 1)

       n_unp = nupdwn(1)
       n_dwn = nupdwn(2)
       is_dw = iupdwn(2) 
       np_dw = nbsp 
!
       ALLOCATE( occ( nbspx ) )
!
       occ( 1:np_dw )  = 1.0d0
       occ( nbspx   )  = 0.0d0
!
! c0(dwn_paired) == c0(up_paired)
! cm(dwn_paired) == cm(up_paired)
! the nbspx dwn state has to be empty
!
!
      c0(:, is_dw:np_dw ) = c0(:, 1:n_dwn )
      cm(:, is_dw:np_dw ) = cm(:, 1:n_dwn )
!
      c0(:, nbspx ) = (0.d0, 0.d0)
      cm(:, nbspx ) = (0.d0, 0.d0)
!
     IF( MOD(n_unp, 2) == 0 ) npair = n_unp - 2
     IF( MOD(n_unp, 2) /= 0 ) npair = n_unp - 1

      DO i = 1, npair, 2 
         !
         CALL dforce(i,bec,vkb,c0,c2,c3,rhos(:,1:1),SIZE(rhos,1),ispin,f,n,nspin)
         CALL dforce(i,bec,vkb,c0,c4,c5,rhos(:,2:2),SIZE(rhos,1),ispin,f,n,nspin)
         !
         c2 = occ( i )*(c2 + c4)  
         c3 = occ(i+1)*(c3 + c5) 
         !
         IF( iflag == 2 ) THEN
            cm(:,i)        = c0(:,i)
            cm(:,i+1)      = c0(:,i+1)
         END IF
         !
         IF( ttsde ) THEN
             CALL wave_steepest( cm(:, i  ), c0(:, i  ), emaver, c2 )
             CALL wave_steepest( cm(:, i+1), c0(:, i+1), emaver, c3 )
         ELSE
             CALL wave_verlet( cm(:, i  ), c0(:, i  ), verl1, verl2, emaver, c2 )
             CALL wave_verlet( cm(:, i+1), c0(:, i+1), verl1, verl2, emaver, c3 )
         END IF
         !
         IF ( gstart == 2 ) THEN
                cm(1,  i)    = CMPLX(DBLE(cm(1,  i)),0.d0)
                cm(1, i+1)   = CMPLX(DBLE(cm(1,  i+1)),0.d0)
         END IF
      !
      END DO
      !
      IF( MOD(n_unp, 2) == 0 ) THEN

         npair = n_unp - 1 
!
         CALL dforce(npair,bec,vkb,c0,c2,c3,rhos(:,1:1),SIZE(rhos,1),ispin,f,n,nspin)
         CALL dforce(npair,bec,vkb,c0,c4,c5,rhos(:,2:2),SIZE(rhos,1),ispin,f,n,nspin)
!
         c2 = c2 + c4
         !
         IF( iflag == 2 ) cm( :, npair ) = c0( :, npair )
!
         IF( ttsde ) THEN
           CALL wave_steepest( cm(:, npair  ), c0(:, npair  ), emaver, c2 )
         ELSE
           CALL wave_verlet( cm(:, npair), c0(:, npair), verl1, verl2, emaver, c2 )
         ENDIF
!
         IF ( gstart == 2 ) cm(1, npair) = CMPLX(DBLE(cm(1, npair)),0.d0)

      ENDIF
!
      c0(:, is_dw:np_dw ) = c0(:, 1:n_dwn )
      cm(:, is_dw:np_dw ) = cm(:, 1:n_dwn )
!
      c0(:, nbspx ) = (0.d0, 0.d0)
      cm(:, nbspx ) = (0.d0, 0.d0)
!

!
! The electron unpaired is signed by n_unp and spin up 
! for the unpaired electron the ei_unp is the value of lambda
! "TRUE" ONLY WHEN THE POT is NORM_CONSERVING
!

      CALL dforce( n_unp, bec, vkb, c0, c2, c3, rhos, SIZE(rhos,1), ispin,f,n,nspin )
      !
      intermed  = - 2.d0 * sum(c2 * conjg(c0(:,n_unp)))
      IF ( gstart == 2 ) THEN
        intermed  = intermed + 1.d0 * c2(1) * conjg(c0(1,n_unp))
      END IF
      CALL mp_sum ( intermed, intra_image_comm )
      !           
      IF( iflag == 2 ) cm(:, n_unp) = c0(:, n_unp) 
      !
      IF( ttsde ) THEN
        CALL wave_steepest( cm(:, n_unp), c0(:, n_unp), emaver, c2 )
      ELSE
        CALL wave_verlet( cm(:, n_unp), c0(:, n_unp), verl1, verl2, emaver, c2 )
      ENDIF 
      !
      IF ( gstart == 2 ) cm(1, n_unp) = CMPLX(DBLE(cm(1, n_unp)),0.d0)
      !
      DEALLOCATE( occ )
      DEALLOCATE( emadt2 )
      DEALLOCATE( emaver )
      DEALLOCATE(c2, c4)
      DEALLOCATE(c3, c5)

   END SUBROUTINE runcp_uspp_force_pairing_x

