!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE move_electrons_x( nfi, tfirst, tlast, b1, b2, b3, fion, &
                           enthal, enb, enbi, fccc, ccc, dt2bye, stress, tprint_ham )
  !----------------------------------------------------------------------------
  !
  ! ... this routine updates the electronic degrees of freedom
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : lwf, tfor, tprnfor, thdyn, use_task_groups
  USE cg_module,            ONLY : tcg
  USE cp_main_variables,    ONLY : eigr, bec, irb, eigrb, rhog, rhos, rhor, &
                                   ei1, ei2, ei3, sfac, ema0bg, becdr, &
                                   taub, lambda, lambdam, lambdap, vpot
  USE wavefunctions_module, ONLY : c0, cm, phi => cp
  USE cell_base,            ONLY : omega, ibrav, h, press, a1, a2, a3
  use grid_dimensions,      only : nr1, nr2, nr3, nr1x, nr2x, nr3x, nnrx
  USE uspp,                 ONLY : becsum, vkb, nkb
  USE energies,             ONLY : ekin, enl, entropy, etot
  USE grid_dimensions,      ONLY : nnrx
  USE electrons_base,       ONLY : nbsp, nbspx, nspin, f, nudx
  USE core,                 ONLY : nlcc_any, rhoc
  USE ions_positions,       ONLY : tau0
  USE ions_base,            ONLY : nat
  USE dener,                ONLY : detot, denl, dekin6
  USE efield_module,        ONLY : tefield, ipolp, qmat, gqq, evalue, &
                                   tefield2, ipolp2, qmat2, gqq2, evalue2
  !
  USE wannier_subroutines,  ONLY : get_wannier_center, wf_options, &
                                   write_charge_and_exit, ef_tune
  USE ensemble_dft,         ONLY : compute_entropy2, z0t, c0diag, becdiag, tens, tsmear, fmat0_diag
  USE efield_module,        ONLY : berry_energy, berry_energy2
  USE cp_interfaces,        ONLY : runcp_uspp, runcp_uspp_force_pairing, &
                                   interpolate_lambda
  USE gvecw,                ONLY : ngw
  USE orthogonalize_base,   ONLY : calphi
  USE control_flags,        ONLY : force_pairing
  USE cp_interfaces,        ONLY : rhoofr, compute_stress, invfft
  USE electrons_base,       ONLY : ispin, iupdwn, nupdwn 
  USE mp_global,            ONLY : me_image, intra_image_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  USE efield_mod,           ONLY : do_efield
  USE fft_base,             ONLY : dfftp
  USE io_global,            ONLY : ionode, ionode_id
  USE hfmod,                ONLY : do_hf, vxxpsi, exx
  USE nksic,                ONLY : do_orbdep, vsic, wtot, fsic, fion_sic, deeq_sic, pink
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: nfi
  LOGICAL,  INTENT(IN)    :: tfirst, tlast
  REAL(DP), INTENT(IN)    :: b1(3), b2(3), b3(3)
  REAL(DP)                :: fion(:,:)
  REAL(DP), INTENT(IN)    :: dt2bye
  REAL(DP)                :: fccc, ccc
  REAL(DP)                :: enb, enbi
  REAL(DP)                :: enthal
  REAL(DP)                :: ei_unp
  REAL(DP)                :: stress(3,3)
  LOGICAL, OPTIONAL, INTENT(IN) :: tprint_ham
  !
  INTEGER :: i, j, is, n2
!$$ The following local variables are for the inner-loop, i.e., unitary rotation
  REAL(DP)                :: vsicah(nbspx,nbspx),vsich(nbspx,nbspx)
  COMPLEX(DP)             :: psi1(nnrx), psi2(nnrx)
  REAL(DP)                   :: vsicah2sum, vsicahtmp, vsich2sum, vsichtmp
  INTEGER :: nbnd1,nbnd2,no_iter
  REAL(DP)                   :: Omat(nbspx,nbspx), Omattmp(nbspx,nbspx), Omattot(nbspx,nbspx), Wmat(nbspx,nbspx)
  REAL(DP), save             :: ncalled = 0
!$$ for a test
  REAL(DP)                   :: Omatim(nbspx,nbspx)
!$$
  COMPLEX(DP)                :: Hmat(nbspx,nbspx), Umat(nbspx,nbspx), Cmattmp(nbspx,nbspx)
  REAL(DP)                   :: Heig(nbspx), alpha_sd, dwfnnorm, dtmp
  COMPLEX(DP)                :: exp_iHeig(nbspx)
  COMPLEX(DP)                :: ci
  COMPLEX(DP)                :: wfn_ctmp(ngw,nbspx), wfn_ctmp2(ngw,nbspx)
  REAL(DP)                   :: dtmp1,dtmp2
!$$
  !
  !
  ci = (0.d0,1.d0)
  dwfnnorm = 1.0/(dble(nr1x)*dble(nr2x)*dble(nr3x))
!$$
  electron_dynamic: IF ( tcg ) THEN
     !
     CALL runcg_uspp( nfi, tfirst, tlast, eigr, bec, irb, eigrb, &
                      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac, &
                      fion, ema0bg, becdr, lambdap, lambda, vpot  )
     !
     CALL compute_stress( stress, detot, h, omega )
     !
  ELSE
     !
     IF ( lwf ) &
          CALL get_wannier_center( tfirst, cm, bec, eigr, &
                                   eigrb, taub, irb, ibrav, b1, b2, b3 )
     !
     IF ( .NOT. tsmear ) THEN
         !
         ! std implementation
         !
         CALL rhoofr( nfi, c0, irb, eigrb, bec, &
                         becsum, rhor, rhog, rhos, enl, denl, ekin, dekin6 )
         !
     ELSE
         !
         ! take into account of the proper density matrix
         ! and rotate orbitals back to theie diagonal representation
         ! to compute rho and other physical quantities, like the kinetic energy
         !
         ! rotates the wavefunctions c0 and the overlaps bec
         ! (the occupation matrix f_ij becomes diagonal f_i)
         !
         CALL rotate( z0t, c0, bec, c0diag, becdiag, .false. )
         !
         CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, &
                         becsum, rhor, rhog, rhos, enl, denl, ekin, dekin6 )
         !
     ENDIF

     !
     ! ... put core charge (if present) in rhoc(r)
     !
     IF ( nlcc_any ) CALL set_cc( irb, eigrb, rhoc )
     !
     IF ( lwf ) THEN
        !
        CALL write_charge_and_exit( rhog )
        CALL ef_tune( rhog, tau0 )
        !
     END IF
     !
     vpot = rhor
     !
     CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
                  ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion )

     !
     ! compute auxiliary potentials
     !
     if( do_orbdep ) then
         !
         if ( tens .or. tsmear) then
             fsic = fmat0_diag
         else
             fsic = f
         endif
         !

         call nksic_potential( nbsp, nbspx, c0, fsic, bec, becsum, deeq_sic, &
                    ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )
         !
!$$ We should update etot only once at the end of this do_orbdep routine
!$$         etot = etot + sum(pink(1:nbsp))
         !

!$$ the following is for tests now.
!         if(ionode) write(1034,*) sum(pink(:)), etot-sum(pink(:))

!$$         if( ncalled .ge. 0 ) then
!$$ Currently, the inner loop routine is disabled: will be updated later
!$$ Also, it is far from optimization.
!$$
         if( .false. ) then
!$$ pseudo code
!$$ if(.true.)
!$$   do while (.true.) loop
!$$     (1) calculate the anti-hermitian part Wmat of the vsic which is real.
!$$     (2) if converged, exit
!$$     (3) diagonalize Hmat = iWmat which is hermitian
!$$     (4) within the steepest descent, calculate Omat = exp(coeff * i * Hmat)
!$$     (5) Omattot -> Omattot x Omat
!$$     (6) rotate c0 by Omat
!$$     (7) call nksic_potential
!$$   enddo
!$$ endif
!$$ (8) outside "if (do_orbdep)", rotate cm by the converged Omat


!$$ Initialization of Omattot
           Omattot(:,:)=0.d0
           do nbnd1=1,nbspx
             Omattot(nbnd1,nbnd1)=1.d0
           enddo

           no_iter = 0

!$$
!$$
!$$ alpha_sd = size of the step for steepest descent matrix change
!$$            will be made an external input variable later.

!$$ Delta W = - del_W E_SIC

!$$
           alpha_sd = 1.d0
!$$           alpha_sd = 0.0
!$$
           do nbnd1=1,nbspx
             wfn_ctmp(:,nbnd1)=c0(:,nbnd1)
           enddo

           do while (.true.)
!$$             if(no_iter .ge. 300) then
             if(no_iter .ge. 1) then
!               if(ionode) write(1033,*)
               exit
             endif

             no_iter = no_iter + 1

!$$ (1) This part calculates the anti-hermitian Wmat and see whether a convergence has been achieved
             vsicah2sum  = 0.d0
             vsich2sum   = 0.d0
             vsicah(:,:) = 0.d0
             vsich(:,:)  = 0.d0

             do nbnd1=1,nbspx
               CALL c2psi( psi1, nnrx, wfn_ctmp(:,nbnd1), (0.d0,0.d0), ngw, 1)
               CALL invfft('Dense', psi1, dfftp )

               do nbnd2=1,nbspx
                 CALL c2psi( psi2, nnrx, wfn_ctmp(:,nbnd2), (0.d0,0.d0), ngw, 1)
                 CALL invfft('Dense', psi2, dfftp )

                 vsichtmp = 0.d0
                 vsicahtmp = 0.d0

                 do i=1,nnrx
                   vsichtmp  = vsichtmp  + dble( conjg(psi1(i)) * (vsic(i,nbnd2)+vsic(i,nbnd1)) * psi2(i) ) * dwfnnorm
                   vsicahtmp = vsicahtmp + dble( conjg(psi1(i)) * (vsic(i,nbnd2)-vsic(i,nbnd1)) * psi2(i) ) * dwfnnorm
                 enddo

                 CALL mp_sum(vsichtmp,intra_image_comm)
                 CALL mp_sum(vsicahtmp,intra_image_comm)

                 vsich2sum = vsich2sum + vsichtmp*vsichtmp
                 vsicah2sum = vsicah2sum + vsicahtmp*vsicahtmp

                 vsich(nbnd1,nbnd2) = vsichtmp
                 vsicah(nbnd1,nbnd2) = vsicahtmp
               enddo
             enddo
             vsich2sum  = vsich2sum
             vsicah2sum = vsicah2sum

             if(ionode) write(1023,*) vsich2sum, vsicah2sum

!$$ print out ESIC part & other total energy
!             if(ionode) write(1033,*) sum(pink(:)), etot
!$$

!$$ (2) if converged, exit
             if(vsicah2sum.le.1.d-8) then
!               if(ionode) write(1023,*) '# inner-loop converged.'
!               if(ionode) write(1023,*)
!               if(ionode) write(1033,*)
!$$ For just one time rotation, we should not exit this loop: I know it is confusing.
!$$ will come back here soon.
!$$               exit
             endif

!$$ (3) Now this part diagonalizes Hmat = iWmat
             Hmat(:,:) = ci * vsicah(:,:)
!$$ diagonalize Hmat
             if(ionode) then
               CALL zdiag(nbspx,nbspx,Hmat(1,1),Heig(1),Umat(1,1),1)
             endif

             CALL mp_bcast(Umat, ionode_id, intra_image_comm)
             CALL mp_bcast(Heig, ionode_id, intra_image_comm)

             do nbnd1=1,nbspx
               dtmp = alpha_sd*Heig(nbnd1)
               exp_iHeig(nbnd1) = cos(dtmp) + ci*sin(dtmp)
             enddo

!$$ (4) Cmattmp = exp(i * alpha * Heig) * Umat   ; Omat = Umat^dagger * Cmattmp
             do nbnd1=1,nbspx
               Cmattmp(:,nbnd1) = Umat(:,nbnd1) * exp_iHeig(nbnd1)
             enddo
             Omat = dble ( MATMUL(Cmattmp, conjg(transpose(Umat)) ) )

!$$ The following is to check whether Omat is actually real or not
             Omatim = aimag ( MATMUL(Cmattmp, conjg(transpose(Umat)) ) )

             dtmp1 = 0.d0
             dtmp2 = 0.d0

             do nbnd1=1,nbspx
               do nbnd2=1,nbspx
                 dtmp1 = dtmp1 + Omat(nbnd1,nbnd2)**2
                 dtmp2 = dtmp2 + Omatim(nbnd1,nbnd2)**2
               enddo
               if(ionode) write(1027,*) (Omat(nbnd1,nbnd2),nbnd2=1,nbspx)
             enddo

             if(ionode) write(1027,*)

             if(ionode) write(1025,*) dtmp1,dtmp2


!$$ (5) Update Omattot
             Omattmp = MATMUL(Omattot,Omat)
             Omattot = Omattmp

!$$ (6) Wavefunction c0 rotation using according to Omat
             wfn_ctmp2(:,:) = (0.d0,0.d0)
             do nbnd1=1,nbspx
               do nbnd2=1,nbspx
                 wfn_ctmp2(:,nbnd1)=wfn_ctmp2(:,nbnd1) + wfn_ctmp(:,nbnd2) * Omat(nbnd2,nbnd1)
               enddo
             enddo
             wfn_ctmp(:,1:nbspx) = wfn_ctmp2(:,1:nbspx)

!$$ (7) recalculate vsic according to the new wavefunction
             call nksic_potential( nbsp, nbspx, wfn_ctmp, fsic, bec, becsum, deeq_sic, &
                        ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )
             !
!$$             etot = etot + sum(pink(1:nbsp))
             !
           enddo  !$$ do while (.true.)

!$$ for tests
           if(ionode) write(1031,*) tfirst, tlast, no_iter
!$$


           do nbnd1=1,nbspx
             if(ionode) write(1029,*) (Omattot(nbnd1,nbnd2),nbnd2=1,nbspx)
           enddo

           dtmp1 = 0.0
           do nbnd1=1,nbspx
             do nbnd2=1,nbspx
               dtmp1 = dtmp1 + Omattot(nbnd1,nbnd2)**2
             enddo
           enddo
           if(ionode) write(1029,*) dtmp1

           if(ionode) write(1029,*)



!$$ (8) Wavefunction cm rotation according to Omattot
           wfn_ctmp2(:,:) = (0.d0,0.d0)
           do nbnd1=1,nbspx
             do nbnd2=1,nbspx
               wfn_ctmp2(:,nbnd1)=wfn_ctmp2(:,nbnd1) + cm(:,nbnd2) * Omattot(nbnd2,nbnd1)
!$$
!$$ The following is to revert c0 to the wavefunction before rotation
!$$               wfn_ctmp2(:,nbnd1)=wfn_ctmp2(:,nbnd1) + wfn_ctmp(:,nbnd2) * Omattot(nbnd1,nbnd2)
             enddo
           enddo
           cm(:,1:nbspx) = wfn_ctmp2(:,1:nbspx)
!$$
!$$ The following is to revert c0 to the wavefunction before rotation
!$$           c0(:,1:nbspx) = wfn_ctmp2(:,1:nbspx)
           call nksic_potential( nbsp, nbspx, c0, fsic, bec, becsum, deeq_sic, &
                      ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )
           !
!$$           etot = etot + sum(pink(1:nbsp))
           !

           wfn_ctmp2(:,:) = (0.d0,0.d0)
           do nbnd1=1,nbspx
             do nbnd2=1,nbspx
               wfn_ctmp2(:,nbnd1)=wfn_ctmp2(:,nbnd1) + c0(:,nbnd2) * Omattot(nbnd2,nbnd1)
             enddo
           enddo
           c0(:,1:nbspx) = wfn_ctmp2(:,1:nbspx)
!$$
!$$ The following is for a test
!$$           c0(:,1:nbspx) = wfn_ctmp(:,1:nbspx)
!$$

         endif !$$ if(.true.) for tests now it is if(.false.)

         ncalled = ncalled + 1

         !
!$$ fsic has been defined at the beginning of do_orbdep routine
!$$         if ( tens .or. tsmear) then
!$$             fsic = fmat0_diag
!$$         else
!$$             fsic = f
!$$         endif
         !
         call nksic_potential( nbsp, nbspx, c0, fsic, bec, becsum, deeq_sic, &
                    ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )
         !
         etot = etot + sum(pink(1:nbsp))
         !
     endif !if( do_orbdep )
     !
     if( do_hf ) then
         !
         call hf_potential( nbsp, nbspx, c0, f, ispin, iupdwn, nupdwn, &
                            nbsp, nbspx, c0, f, ispin, iupdwn, nupdwn, &
                            rhor, rhog, vxxpsi, exx)
         !
         etot = etot + sum(exx(1:nbsp))
         !
     endif
     !
     if( do_efield ) then
         !
         call calc_dipole(c0, h)
         !
     endif
     !
     IF ( lwf ) CALL wf_options( tfirst, nfi, cm, becsum, bec, &
                                 eigr, eigrb, taub, irb, ibrav, b1,   &
                                 b2, b3, vpot, rhog, rhos, enl, ekin  )
     !
     CALL compute_stress( stress, detot, h, omega )
     !
     enthal = etot + press * omega
     !
     IF( tefield )  THEN
        !
        CALL berry_energy( enb, enbi, bec, c0, fion )
        !
        etot = etot + enb + enbi
        !
     END IF
     IF( tefield2 )  THEN
        !
        CALL berry_energy2( enb, enbi, bec, c0, fion )
        !
        etot = etot + enb + enbi
        !
     END IF

     !
     !=======================================================================
     !
     !              verlet algorithm
     !
     !     loop which updates electronic degrees of freedom
     !     cm=c(t+dt) is obtained from cm=c(t-dt) and c0=c(t)
     !     the electron mass rises with g**2
     !
     !=======================================================================
     !
     ! This call must be done after the call to nksic_potential
     ! or nksic_inner_loop
     !
     CALL newd( vpot, irb, eigrb, becsum, fion )
     !
     if( do_orbdep ) then
         !
         fion = fion + fion_sic
         !
     endif
     !
     CALL prefor( eigr, vkb )
     !
     IF( force_pairing ) THEN
        !
        CALL runcp_uspp_force_pairing( nfi, fccc, ccc, ema0bg, dt2bye, &
                      rhos, bec, c0, cm, ei_unp )
        !
     ELSE
        !
        CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, &
                         tprint_ham = tprint_ham )
        !
     ENDIF
     !
     !----------------------------------------------------------------------
     !                 contribution to fion due to lambda
     !----------------------------------------------------------------------
     !
     ! ... nlfq needs deeq bec
     !
     IF ( tfor .OR. tprnfor ) CALL nlfq( c0, eigr, bec, becdr, fion )
     !
     IF ( (tfor.or.tprnfor) .AND. tefield ) &
        CALL bforceion( fion, .TRUE. , ipolp, qmat, bec, becdr, gqq, evalue )
     IF ( (tfor.or.tprnfor) .AND. tefield2 ) &
        CALL bforceion( fion, .TRUE. , ipolp2, qmat2, bec, becdr, gqq2, evalue2 )
     !
     IF( force_pairing ) THEN
        lambda( :, :, 2 ) =  lambda(:, :, 1 )
        lambdam( :, :, 2 ) = lambdam(:, :, 1 )
     ENDIF
     ! 
     IF ( tfor .OR. thdyn ) then
        CALL interpolate_lambda( lambdap, lambda, lambdam )
     ELSE
        ! take care of the otherwise uninitialized lambdam
        lambdam = lambda
     END IF
     !
     ! ... calphi calculates phi
     ! ... the electron mass rises with g**2
     !
     CALL calphi( c0, ngw, bec, nkb, vkb, phi, nbsp, ema0bg )
     !
     ! ... begin try and error loop (only one step!)
     !
     ! ... nlfl and nlfh need: lambda (guessed) becdr
     !
     IF ( tfor .OR. tprnfor ) CALL nlfl( bec, becdr, lambda, fion )
     !
  END IF electron_dynamic
  !
  RETURN
  !
END SUBROUTINE move_electrons_x
