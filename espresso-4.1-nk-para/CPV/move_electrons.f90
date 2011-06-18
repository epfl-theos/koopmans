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
!$$ CHP (June 18 / 2011)
!$$
!$$ When damped dynamics is used, an optimal unitary rotation among occupied states is
!$$ calculated by setting 'do_innerloop=.true.' in the namelist SYSTEM of the input file.
!$$ Convergence for this inner loop minimization and (conventional) outer loop minimization
!$$ can be found in the files convg_inner.dat and convg_outer.dat, respectively.
!$$
!$$ Nota Bene: This inner-loop minimization cannot be used for the case where an energy
!$$ functional cannot be defined like nk0.
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
  USE nksic,                ONLY : do_pz, do_innerloop
  use ions_base, only: nsp
  use electrons_base, only: nel,nelt,nupdwn,iupdwn
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
  REAL(DP)                :: bec1(nkb,nbsp)
  REAL(DP), allocatable   :: vsicah(:,:),overlap(:,:)
  COMPLEX(DP)             :: psi1(nnrx), psi2(nnrx)
  REAL(DP)                   :: vsicahtmp,overlaptmp
  INTEGER :: nbnd1,nbnd2,ninner
  REAL(DP), allocatable      :: Omat1(:,:),Omatim(:,:)
  REAL(DP)                   :: Omat1tot(nbspx,nbspx)
  REAL(DP)                   :: Omattot(nbspx,nbspx)
  REAL(DP)                   :: vsic1(nnrx,nbspx)
  INTEGER , save             ::  nouter = 0
  REAL(DP)                   :: ene0
  COMPLEX(DP), allocatable   :: Hmat(:,:), Umat(:,:), Cmattmp(:,:)
  REAL(DP), allocatable      :: Heig(:)
  REAL(DP)                   :: dwfnnorm, dtmp
  COMPLEX(DP), allocatable   :: exp_iHeig(:)
  COMPLEX(DP)                :: ci
  COMPLEX(DP)                :: wfn_ctmp(ngw,nbspx), wfn_ctmp2(ngw,nbspx), wfn_ctmp1(ngw,nbspx)
  REAL(DP)                   :: passof,passoprodmin
  INTEGER                    :: npassofailmax
  INTEGER, SAVE              :: npassofail=0
  REAL(DP), SAVE             :: pinksumprev=1.d8, passoprod=0.3d0
  REAL(DP)                   :: pink1(nbspx)
  INTEGER                    :: isp
!  COMPLEX(DP)                :: Htest(2,2), Utest(2,2)
!  REAL(DP)                   :: Eigtest(2)
!$$
  !
  !
!$$
  ci = (0.d0,1.d0)
  dwfnnorm = 1.0/(dble(nr1x)*dble(nr2x)*dble(nr3x))
  npassofailmax = 5 ! when to stop dividing passoprod by 2
!$$

!$$ LAPACK TEST
!  Htest(1,1) =  0.0
!  Htest(1,2) = -ci
!  Htest(2,1) = ci
!  Htest(2,2) = 0.0
!  if(ionode) then
!    CALL zdiag(2,2,Htest(1,1),Eigtest(1),Utest(1,1),1)
!  endif
!
!  CALL mp_bcast(Utest, ionode_id, intra_image_comm)
!  CALL mp_bcast(Htest, ionode_id, intra_image_comm)
!
!  if(ionode) then
!    write(555,*) 'Printing out Htest ...'
!    write(555,*) Htest(1,1),Htest(1,2)
!    write(555,*) Htest(2,1),Htest(2,2)
!
!    write(555,*) 'Printing out Utest ...'
!    write(555,*) Utest(1,1),Utest(1,2)
!    write(555,*) Utest(2,1),Utest(2,2)
!
!    write(555,*) 'Printing out Eigtest ...'
!    write(555,*) Eigtest(1),Eigtest(2)
!
!    write(555,*) 'nbsp,nbspx,nspin,nudx', nbsp,nbspx,nspin,nudx
!    write(555,*) 'nel(2),nelt,nupdwn(2),iupdwn(2)',nel,nelt,nupdwn,iupdwn
!  endif
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

         ene0 = sum(pink(1:nbsp))
         !

         nouter = nouter + 1
         if(ionode.and.( nouter.eq.1)) then
           open(1032,file='convg_outer.dat',status='unknown')
           write(1032,'("#   ninner    nouter     non-sic energy (Ha)         sic energy (Ha)")')
         endif

!$$ Inner loop convergence is performed only once at the first iteration:
!$$ therefore it is better to start from LDA wavefunctions which are good
!$$ except for the unitary rotation.

         ninner = 0

         if( do_innerloop .and. ( nouter.eq.1) ) then

           Omattot(:,:)=0.d0
           do nbnd1=1,nbspx
             Omattot(nbnd1,nbnd1)=1.d0
           enddo

           do nbnd1=1,nbspx
             wfn_ctmp(:,nbnd1)=c0(:,nbnd1)
           enddo

           if(ionode) then
             open(1031,file='convg_inner.dat',status='unknown')
             write(1031,'("#   ninner    nouter     non-sic energy (Ha)         sic energy (Ha)")')
           endif

           do while (.true.)

             ninner = ninner + 1

!$$ print out ESIC part & other total energy
             if(ionode) write(1031,'(2I10,2F24.13)') ninner, nouter,etot,ene0

             Omat1tot(:,:)=0.d0
             do nbnd1=1,nbspx
               Omat1tot(nbnd1,nbnd1)=1.d0
             enddo


!$$ This part calculates the anti-hermitian Wmat and see whether a convergence has been achieved

             wfn_ctmp1(:,:) = (0.d0,0.d0)

             do isp=1,nspin

               allocate(vsicah(nupdwn(isp),nupdwn(isp)))
               allocate(overlap(nupdwn(isp),nupdwn(isp)))
               allocate(Hmat(nupdwn(isp),nupdwn(isp)))
               allocate(Umat(nupdwn(isp),nupdwn(isp)))
               allocate(Heig(nupdwn(isp)))
               allocate(exp_iHeig(nupdwn(isp)))
               allocate(Cmattmp(nupdwn(isp),nupdwn(isp)))
               allocate(Omat1(nupdwn(isp),nupdwn(isp)))

               vsicah(:,:) = 0.d0

               do nbnd1=1,nupdwn(isp)
                 CALL c2psi( psi1, nnrx, wfn_ctmp(:,iupdwn(isp)-1 + nbnd1), (0.d0,0.d0), ngw, 1)
                 CALL invfft('Dense', psi1, dfftp )

                 do nbnd2=1,nupdwn(isp)
                   if(nbnd2.lt.nbnd1) then
                     vsicahtmp = -vsicah(nbnd2,nbnd1)
                     overlaptmp = overlap(nbnd2,nbnd1)
                   else
                     CALL c2psi( psi2, nnrx, wfn_ctmp(:,iupdwn(isp)-1 + nbnd2), (0.d0,0.d0), ngw, 1)
                     CALL invfft('Dense', psi2, dfftp )

                     vsicahtmp = 0.d0
                     overlaptmp = 0.d0

                     do i=1,nnrx
!$$ Imposing Pederson condition
                       vsicahtmp = vsicahtmp &
                           + 2.d0 * dble( conjg(psi1(i)) * (vsic(i,nbnd2)-vsic(i,nbnd1)) * psi2(i) ) * dwfnnorm
!                           + 2.d0 * dble( vsic(i,nbnd2)-vsic(i,nbnd1) ) * dwfnnorm
                       overlaptmp = overlaptmp + dble( conjg(psi1(i)) * psi2(i) ) * dwfnnorm
                     enddo

                     CALL mp_sum(vsicahtmp,intra_image_comm)
                     CALL mp_sum(overlaptmp,intra_image_comm)
                   endif ! if(nbnd2.lt.nbnd1)

                   vsicah(nbnd1,nbnd2) = vsicahtmp
                   overlap(nbnd1,nbnd2) = overlaptmp

                 enddo ! nbnd2=1,nupdwn(isp)
               enddo ! nbnd1=1,nupdwn(isp)

!$$ Now this part diagonalizes Hmat = iWmat
               Hmat(:,:) = ci * vsicah(:,:)
!$$ diagonalize Hmat
               if(ionode) then
                 CALL zdiag(nupdwn(isp),nupdwn(isp),Hmat(1,1),Heig(1),Umat(1,1),1)
               endif

               CALL mp_bcast(Umat, ionode_id, intra_image_comm)
               CALL mp_bcast(Heig, ionode_id, intra_image_comm)

!$$ We set the step size in such a way that the phase change
!$$ of the wavevector with largest eigenvalue upon rotation is fixed
               passof = passoprod/max(abs(Heig(1)),abs(Heig(nupdwn(isp))))

               do nbnd1=1,nupdwn(isp)
                 dtmp =  passof * Heig(nbnd1)
                 exp_iHeig(nbnd1) = cos(dtmp) + ci*sin(dtmp)
               enddo

!$$ Cmattmp = exp(i * passof * Heig) * Umat   ; Omat = Umat^dagger * Cmattmp
               do nbnd1=1,nupdwn(isp)
                 Cmattmp(:,nbnd1) = Umat(:,nbnd1) * exp_iHeig(nbnd1)
               enddo

               Omat1 = dble ( MATMUL(Cmattmp, conjg(transpose(Umat)) ) )

!$$ Wavefunction c0 rotation using according to Omat
               do nbnd1=1,nupdwn(isp)
                 do nbnd2=1,nupdwn(isp)
                   wfn_ctmp1(:,iupdwn(isp)-1 + nbnd1)=wfn_ctmp1(:,iupdwn(isp)-1 + nbnd1) &
                       + wfn_ctmp(:,iupdwn(isp)-1 + nbnd2) * Omat1(nbnd2,nbnd1)
                 enddo
               enddo

! Assigning the rotation matrix for a specific spin isp
               Omat1tot(iupdwn(isp):iupdwn(isp)-1+nupdwn(isp),iupdwn(isp):iupdwn(isp)-1+nupdwn(isp)) = Omat1(:,:)

!               if(ionode) write(1041,*) ninner, nouter,isp
!               if(ionode) write(1042,*) ninner, nouter,isp
!               if(ionode) write(1043,*) ninner, nouter,isp
!!               if(ionode) write(1041,*) pink  ! results are okay
!!               if(ionode) write(1041,*) pink1 ! 15th result is NaN
!               do nbnd1=1,nupdwn(isp)
!                 if(ionode) write(1041,'(15E9.1)') (overlap(nbnd1,nbnd2),nbnd2=1,nupdwn(isp))
!                 if(ionode) write(1042,'(15E9.1)') (Omat1(nbnd1,nbnd2),nbnd2=1,nupdwn(isp))
!                 if(ionode) write(1043,'(15E9.1)') (vsicah(nbnd1,nbnd2),nbnd2=1,nupdwn(isp))
!               enddo
!!               do nbnd1=1,nbspx
!!                 if(ionode) write(1041,'(30E9.1)') (Omat1tot(nbnd1,nbnd2),nbnd2=1,nbspx)
!!               enddo

               deallocate(vsicah)
               deallocate(overlap)
               deallocate(Umat)
               deallocate(Hmat)
               deallocate(Heig)
               deallocate(exp_iHeig)
               deallocate(Cmattmp)
               deallocate(Omat1)

             enddo ! do isp=1,nspin

!$$ recalculate bec & vsic according to the new wavefunction
             call calbec(1,nsp,eigr,wfn_ctmp1,bec1)

             vsic1(:,:) = 0.d0
             pink1(:) = 0.d0

             call nksic_potential( nbsp, nbspx, wfn_ctmp1, fsic, bec1, becsum, deeq_sic, &
                        ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic1, pink1 )
!$$ Currently, this nksic_potential call messes up pink1 and vsic1 !!!
!             call nksic_potential( nbsp, nbspx, c0, fsic, bec1, becsum, deeq_sic, &
!                        ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic1, pink1 )


!$$ if converged, exit
             if(sum(pink1(:)).gt.ene0) then
               npassofail = npassofail+1

               if(ionode) then
                 write(1031,'("# procedure  ",I2," / ",I2," is finished.")') npassofail,npassofailmax
                 write(1031,*)
               endif

               if(npassofail.eq.npassofailmax) then
!$$ if we reach at the maximum allowed npassofail number, we exit without further update
                 exit
               endif
               passoprod = passoprod * 0.5d0
               cycle
             endif

!$$ we keep track of all the rotations to rotate cm later
             Omattot = MATMUL(Omattot,Omat1tot)

             pink(:) = pink1(:)
             ene0 = sum(pink1(:))
             vsic(:,:) = vsic1(:,:)
             wfn_ctmp(:,:) = wfn_ctmp1(:,:)
             bec(:,:) = bec1(:,:)
           enddo  !$$ do while (.true.)

!$$ Now update c0 (wfn_ctmp is in fact not necessary but ...)
           c0(:,1:nbspx) = wfn_ctmp(:,1:nbspx)

!$$ Wavefunction cm rotation according to Omattot
           wfn_ctmp2(:,:) = (0.d0,0.d0)
           do nbnd1=1,nbspx
             do nbnd2=1,nbspx
               wfn_ctmp2(:,nbnd1)=wfn_ctmp2(:,nbnd1) + cm(:,nbnd2) * Omattot(nbnd2,nbnd1)
             enddo
           enddo
           cm(:,1:nbspx) = wfn_ctmp2(:,1:nbspx)

           if(ionode) close(1031)

         endif !$$ if( do_innerloop .and. ( nouter.eq.1) )

!$$ to see the outer loop energy convergence
         if(ionode) write(1032,'(2I10,2F24.13)') ninner, nouter,etot,ene0
!$$
         !
         etot = etot + ene0
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
