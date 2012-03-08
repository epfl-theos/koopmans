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
!$$
!$$ CHP (August 10 / 2011)
!$$
!$$ An optimal unitary rotation among occupied states is calculated by setting
!$$ 'do_innerloop=.true.' in the namelist SYSTEM of the input file.
!$$ Also, either CG or SD with parabolic minimization or steepest-descent method without
!$$ parabolic minimization can be chosen: 'do_innerloop_cg' being .true. or .false.
!$$ (Default: do_innerloop_cg=.false.)
!$$ Convergence for this inner loop minimization and (conventional) outer loop minimization
!$$ can be found in the files convg_inner.dat and convg_outer.dat, respectively.
!$$ esic_conv_thr is the threshold of convergence for the inner loop minimization.
!$$
!$$ When do_innerloop_cg=.true., one could also set innerloop_cg_nsd (the number of
!$$ initial steepest-descent steps) and innerloop_cg_nreset (the number of CG step running
!$$ before resetting the CG direction, i.e., when this variable is set to 1, the calculation
!$$ becomes SD).  Initially, SD is better than CG but when near the energy minimum, CG
!$$ works better.
!$$
!$$ When the OUTERLOOP dynamics is damped dynamics, we should set innerloop_dd_nstep,
!$$ which is the number of outerloop steps between each inner loop minimization. When it is
!$$ set to 1, we do inner loop minimization at every outer loop step.
!$$
!$$ Nota Bene:
!$$ 1. When using methods without energy functional such as nk0, do_innerloop_cg
!$$ should be set to .false.  In general, it is faster to set this one to .false.
!$$ 2. Fractional occupation is not supported yet.  This is the case for the outer loop's
!$$ CG routine as well - even without SIC.
!$$

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
!$$
  USE nksic,                ONLY : do_pz, do_innerloop,do_innerloop_cg, innerloop_dd_nstep
  use ions_base,            only : nsp
  use electrons_base,       only : nel, nelt
!$$
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
  INTEGER :: ninner
  REAL(DP)                   :: Omattot(nbspx,nbspx)
  INTEGER , save             ::  nouter = 0
  REAL(DP)                   :: ene0
!$$
  !
  !

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

         ene0 = sum(pink(1:nbsp))
         !

         nouter = nouter + 1
         !
!$$#ifdef __DEBUG
         if( ionode .and.( nouter == 1) ) then
           !
           open(1032,file='convg_outer.dat',status='unknown')
           write(1032,'("#   ninner    nouter     non-sic energy (Ha)         sic energy (Ha)")')

           if(do_innerloop) then
             open(1031,file='convg_inner.dat',status='unknown')
             write(1031,'("#   ninner    nouter     non-sic energy (Ha)         sic energy (Ha)    RMS force eigenvalue")')

             if(do_innerloop_cg) then
               open(1037,file='cg_convg.dat',status='unknown')!for debug and tuning purposes
             endif

           endif
         endif
!$$#endif

!$$ Inner loop convergence is performed only once at the first iteration:
!$$ therefore it is better to start from LDA wavefunctions which are good
!$$ except for the unitary rotation.

         ninner = 0
!$$
!$$ For Benzene, it has been checked that pz and nk both need only one
!$$ inner loop optimization.
!$$

         if( do_innerloop .and. ( nouter == 1 .or. mod(nouter,innerloop_dd_nstep) == 0 ) ) then
             !
             ! if( do_innerloop .and. ( nouter.eq.1) ) then
             ! if( do_innerloop ) then
             ! if(do_innerloop .and. nouter.eq.1) then
             ! if(.false.) then
             !
             if(.not.do_innerloop_cg) then
                 call nksic_rot_emin(nouter,ninner,etot,Omattot)
             else
                 call nksic_rot_emin_cg(nouter,ninner,etot,Omattot)
             endif
             !
             ene0 = sum(pink(:))
             !
         endif

         !
#ifdef __DEBUG
!$$      ! to see the outer loop energy convergence
         if(ionode) write(1032,'(2I10,2F24.13)') ninner, nouter,etot,sum(pink(:))
         !
#endif
!$$
         !
         etot = etot + ene0
         !
!$$
!         if(nouter.eq.1.and.ionode) then
!           write(1033,*) 'fsic',fsic
!           write(1033,*) 'nupdwn',nupdwn
!           write(1033,*) 'iupdwn',iupdwn
!           write(1033,*) 'nbsp,nbspx',nbsp,nbspx
!         endif
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
