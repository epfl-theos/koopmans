!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!=======================================================================
   subroutine runcg_uspp( nfi, tfirst, tlast, eigr, bec, irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac, &
                          fion, ema0bg, becdr, lambdap, lambda, lambda_bare, vpot  )
!=======================================================================

      use kinds,                    only : dp
      use control_flags,            only : iprint, thdyn, tpre, iprsta, &
                                           tfor, taurdr, tprnfor, gamma_only, do_wf_cmplx !added:giovanni gamma_only, do_wf_cmplx
      use control_flags,            only : ndr, ndw, nbeg, nomore, tsde, tortho, tnosee, &
                                           tnosep, trane, tranp, tsdp, tcp, tcap, ampre, &
                                           amprp, tnoseh, non_ortho
      use core,                     only : nlcc_any
      !---ensemble-DFT
      use energies,                 only : eht, epseu, exc, etot, eself, enl, ekin,&
                                           atot, entropy, egrand, eodd
      use electrons_base,           only : f, nspin, nel, iupdwn, nupdwn, nudx, nelt, &
                                           nbspx, nbsp, ispin
      use ensemble_dft,             only : tens, tsmear,   ef,  z0t, c0diag,  &
                                           becdiag, fmat0, fmat0_diag, e0,  id_matrix_init
      !---
      use gvecp,                    only : ngm
      use gvecs,                    only : ngs
      use gvecb,                    only : ngb
      use gvecw,                    only : ngw, ngwx
      use reciprocal_vectors,       only : ng0 => gstart
      use cvan,                     only : nvb, ish
      use ions_base,                only : na, nat, pmass, nax, nsp, rcmax
      use grid_dimensions,          only : nnr => nnrx, nr1, nr2, nr3
      use cell_base,                only : ainv, a1, a2, a3
      use cell_base,                only : omega, alat
      use cell_base,                only : h, hold, deth, wmass, tpiba2
      use smooth_grid_dimensions,   only : nnrsx, nr1s, nr2s, nr3s
      use smallbox_grid_dimensions, only : nnrb => nnrbx, nr1b, nr2b, nr3b
      use local_pseudo,             only : vps, rhops
      use io_global,                ONLY : io_global_start, stdout, ionode, ionode_id
      use mp_global,                ONLY : intra_image_comm, np_ortho, me_ortho, ortho_comm, me_image
      use dener
      use cdvan
      use constants,                only : pi, au_gpa
      use io_files,                 only : psfile, pseudo_dir
      USE io_files,                 ONLY : outdir, prefix
      use uspp,                     only : nhsa=> nkb, nhsavb=> nkbus, betae => vkb, rhovan => becsum, deeq,qq
      use uspp_param,               only : nh
      use cg_module,                only : ene_ok,  maxiter,niter_cg_restart, &
                                           conv_thr, passop, enever, itercg
      use ions_positions,           only : tau0
      use wavefunctions_module,     only : c0, cm, phi => cp, cdual, cmdual, cstart
      use efield_module,            only : tefield, evalue, ctable, qmat, detq, ipolp, &
                                           berry_energy, ctabin, gqq, gqqm, df, pberryel, &
                                           tefield2, evalue2, ctable2, qmat2, detq2, ipolp2, &
                                           berry_energy2, ctabin2, gqq2, gqqm2, pberryel2
      use mp,                       only : mp_sum, mp_bcast
      use cp_electronic_mass,       ONLY : emass_cutoff
      use orthogonalize_base,       ONLY : calphi
      use cp_interfaces,            ONLY : rhoofr, dforce, compute_stress, nlfl, set_x_minus1, xminus1
      USE cp_main_variables,        ONLY : nlax, collect_lambda, distribute_lambda, descla, nrlx, nlam
      USE descriptors,              ONLY : la_npc_ , la_npr_ , la_comm_ , la_me_ , la_nrl_ , ldim_cyclic
      USE mp_global,                ONLY : me_image, my_image_id
      !
      use nksic,                    only : do_orbdep, do_innerloop, do_innerloop_cg, innerloop_cg_nsd, &
                                           innerloop_cg_nreset, innerloop_init_n, innerloop_cg_ratio, &
                                           vsicpsi, vsic, wtot, fsic, fion_sic, deeq_sic, f_cutoff, & 
                                           pink, do_wxd, sizwtot, do_bare_eigs
      use hfmod,                    only : do_hf, vxxpsi, exx
      use twin_types !added:giovanni
      use control_flags,            only : non_ortho
      use cp_main_variables,        only : becdual, becmdual, overlap, ioverlap, becstart
      use electrons_module,         only : wfc_spreads, wfc_centers, icompute_spread, manifold_overlap
      use ldau,                     only : lda_plus_u, vupsi
      use printout_base,            only : printout_base_open, printout_base_unit, &
                                       printout_base_close
      use control_flags,            only : iprint_manifold_overlap, iprint_spreads
!
      implicit none
!
      CHARACTER(LEN=80) :: uname
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      integer     :: nfi
      logical     :: tfirst , tlast
      complex(dp) :: eigr(ngw,nat)
      type(twin_matrix)    :: bec !modified:giovanni
      type(twin_tensor)    :: becdr!(nhsa,nspin*nlax,3) !modified:giovanni
      integer     :: irb(3,nat)
      complex(dp) :: eigrb(ngb,nat)
      real(dp)    :: rhor(nnr,nspin)
      real(dp)    :: vpot(nnr,nspin)
      complex(dp) :: rhog(ngm,nspin)
      real(dp)    :: rhos(nnrsx,nspin)
      real(dp)    :: rhoc(nnr)
      complex(dp) :: ei1(-nr1:nr1,nat)
      complex(dp) :: ei2(-nr2:nr2,nat)
      complex(dp) :: ei3(-nr3:nr3,nat)
      complex(dp) :: sfac( ngs, nsp )
      real(dp)    :: fion(3,nat)
      real(dp)    :: ema0bg(ngw)
      type(twin_matrix) :: lambdap(nspin)!(nlam,nlam,nspin) !modified:giovanni
      type(twin_matrix) :: lambda(nspin)!(nlam,nlam,nspin)   !modified:giovanni
      type(twin_matrix) :: lambda_bare(nspin)     !(nlam,nlam,nspin)   !modified:giovanni
!
!
      integer     :: i, j, ig, k, is, iss,ia, iv, jv, il, ii, jj, kk, ip, isp
      integer     :: inl, jnl, niter, istart, nss, nrl, me_rot, np_rot , comm
      real(dp)    :: enb, enbi, x
      real(dp)    :: entmp, sta
      complex(dp) :: gamma_c  !warning_giovanni, is it real anyway?
      complex(dp), allocatable :: c2(:), c3(:), c2_bare(:), c3_bare(:)
      complex(dp), allocatable :: hpsi(:,:), hpsi0(:,:), gi(:,:), hi(:,:), gi_bare(:,:)
!       real(DP),    allocatable :: s_minus1(:,:)    !factors for inverting US S matrix
!       real(DP),    allocatable :: k_minus1(:,:)    !factors for inverting US preconditioning matrix
      type(twin_matrix) :: s_minus1!(:,:)    !factors for inverting US S matrix
      type(twin_matrix) :: k_minus1!(:,:)    !factors for inverting US preconditioning matrix
      real(DP),    allocatable :: lambda_repl(:,:) ! replicated copy of lambda
      real(DP),    allocatable :: lambda_dist(:,:) ! replicated copy of lambda
      complex(DP),    allocatable :: lambda_repl_c(:,:) ! replicated copy of lambda
      complex(DP),    allocatable :: lambda_dist_c(:,:) ! replicated copy of lambda

      !
      real(dp)    :: sca, dumm(1)
      logical     :: newscheme, firstiter
      integer     :: maxiter3
      !
!       real(kind=DP), allocatable :: bec0(:,:), becm(:,:), becdrdiag(:,:,:)
      type(twin_tensor) :: becdrdiag !modified:giovanni
      type(twin_matrix) :: bec0, becm !modified:giovanni
      real(kind=DP), allocatable :: ave_ene(:)!average kinetic energy for preconditioning
      real(kind=DP), allocatable :: fmat_(:,:)!average kinetic energy for preconditioning
      complex(kind=DP), allocatable :: fmat_c_(:,:)!average kinetic energy for preconditioning
      ! 
      logical     :: pre_state!if .true. does preconditioning state by state
      !
      complex(DP)    :: esse_c,essenew_c !factors in c.g.
      logical     :: ltresh!flag for convergence on energy
      real(DP)    :: passo!step to minimum
      real(DP)    :: etotnew, etotold!energies
      real(DP)    :: eoddnew, eoddold!odd energies
      real(DP)    :: spasso!sign of small step
      logical     :: restartcg!if .true. restart again the CG algorithm, performing a SD step
      integer     :: numok!counter on converged iterations
      integer     :: iter3
      real(DP)    :: passof,passov !step to minimum: effective, estimated
      real(DP)    :: ene0,ene1,dene0,enesti !energy terms for linear minimization along hi
      !
      real(DP),    allocatable :: faux(:) ! takes into account spin multiplicity
      real(DP),    allocatable :: hpsinorm(:), hpsinosicnorm(:)
      complex(DP), allocatable :: hpsinosic(:,:)
      complex(DP), allocatable :: hitmp(:,:)
      integer     :: ninner,nbnd1,nbnd2,itercgeff
      complex(DP)    :: Omattot(nbspx,nbspx)
      real(DP)    :: dtmp, temp
      real(dp)    :: tmppasso, ene_save(100), ene_save2(100), ene_lda
      !
      logical :: lgam, switch=.false., ortho_switch=.false.
      complex(DP) :: phase
      integer :: ierr, northo_flavor
      real(DP) :: deltae,sic_coeff1, sic_coeff2 !coefficients which may change according to the flavour of SIC
      integer :: me, iunit_manifold_overlap, iunit_spreads
      !
      ! Initialize some basic variables
      me = me_image+1
      lgam = gamma_only.and..not.do_wf_cmplx
      deltae = 2.d0*conv_thr
      etotnew=0.d0
      etotold=0.d0
      eoddnew=0.d0
      eoddold=0.d0
      ! Initialize printing
      IF( ionode ) THEN
         iunit_spreads = printout_base_unit( "sha" )
         CALL printout_base_open( "sha" )
         WRITE(iunit_spreads, *) " "
         iunit_manifold_overlap = printout_base_unit( "ovp" )
         CALL printout_base_open( "ovp" )
         WRITE(iunit_manifold_overlap, *) " "
      ENDIF
      !
      northo_flavor=1
      !
      IF(northo_flavor==2) THEN
         !
         sic_coeff1=1.d0
         sic_coeff2=1.d0
         !
      ELSE
         !
         sic_coeff1=0.5d0
         sic_coeff2=0.5d0
         !
      ENDIF
      !
      allocate (faux(nbspx))
      !
      allocate (ave_ene(nbsp))
      allocate (c2(ngw),c3(ngw))

      !begin_added:giovanni
      call init_twin(bec0, lgam)
      call allocate_twin(bec0, nhsa, nbsp, lgam)
      call init_twin(becm, lgam)
      call allocate_twin(becm, nhsa, nbsp, lgam)
      call init_twin(becdrdiag, lgam)
      call allocate_twin(becdrdiag, nhsa, nspin*nlax,3, lgam)

!       do iss=1,nspin
! 	call init_twin(lambda(iss), lgam)
! 	call allocate_twin(lambda(iss), nlam, nlam, lgam)
! 	call init_twin(lambdap(iss), lgam)
! 	call allocate_twin(lambdap(iss), nlam, nlam, lgam)
!       enddo
      !end_added:giovanni

      call start_clock('runcg_uspp')
      newscheme=.false.
      firstiter=.true.

      pre_state=.false.!normally is disabled

      maxiter3=12
!$$
      if(do_orbdep) maxiter3=10
!$$
      ninner=0

      !
      ! the following is just a beginning; many things to be done...
      !
      if(do_orbdep) then
          !
          if ( tens .or. tsmear) then
              fsic = fmat0_diag
          else
              fsic = f
          endif
          !
      endif


#ifdef __DEBUG
      if(ionode) then
         uname = TRIM( outdir ) // "/" // trim(prefix) // '.' &
                 // trim(int_to_char( my_image_id )) // '_' // trim(int_to_char( me_image))
         open(37,file=uname,status='unknown')!for debug and tuning purposes
         open(1037,file='cg_convg.dat',status='unknown')!for debug and tuning purposes
      endif
#endif

      if( tfirst .and. ionode ) &
         write(stdout,"(/,a,/)") 'PERFORMING CONJUGATE GRADIENT MINIMIZATION OF EL. STATES'
      
!set tpa preconditioning

      call  emass_precond_tpa( ema0bg, tpiba2, emass_cutoff )
     
      call prefor(eigr,betae) 

      ltresh    = .false.
      itercg    = 1
      etotold   = 1.d8
      eoddold   = 1.d8
      restartcg = .true.
      passof = passop
      ene_ok = .false.

!$$
      itercgeff = 1
!$$

      !orthonormalize c0
      IF(non_ortho) THEN
         call calbec(1,nsp,eigr,c0,bec)
         call compute_duals(c0,cdual,nbspx,1)
         call calbec(1,nsp,eigr,cdual,becdual)
      ELSE
         IF(do_orbdep.and.ortho_switch) THEN
            call lowdin(c0, lgam)
            call calbec(1,nsp,eigr,c0,bec)
         ELSE
            call calbec(1,nsp,eigr,c0,bec)
            call gram(betae,bec,nhsa,c0,ngw,nbsp)
         ENDIF
      ENDIF
      !call calbec(1,nsp,eigr,c0,bec)
      IF(non_ortho) THEN
         CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, nbsp, lgam)
      ELSE
         CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, nbsp, lgam)
      ENDIF
      !
      ! calculates the factors for S and K inversion in US case
      !
      if ( nvb > 0 ) then
          !
          call init_twin(s_minus1,lgam)
          call allocate_twin(s_minus1, nhsavb, nhsavb, lgam)
          call init_twin(k_minus1,lgam)
          call allocate_twin(k_minus1, nhsavb, nhsavb, lgam)
          call  set_x_minus1_twin(betae,s_minus1,dumm,.false.)
          call  set_x_minus1_twin(betae,k_minus1,ema0bg,.true.)
          !
      else
          !
          call init_twin(s_minus1,lgam)
          call allocate_twin(s_minus1, 1, 1, lgam)
          call init_twin(k_minus1,lgam)
          call allocate_twin(k_minus1, 1, 1, lgam)
          !
      endif  
      !
      ! set index on number of converged iterations
      !
      numok = 0

      allocate( hpsi(ngw,nbsp) )
      allocate( hpsi0(ngw,nbsp) )
      allocate( gi(ngw,nbsp), hi(ngw,nbsp) )
      !
      allocate(hitmp(ngw,nbsp))
      hitmp(:,:) = CMPLX(0.d0,0.d0)
      ! allocate(hpsinosic(ngw,n))
      !
      gi(:,:)=CMPLX(0.d0,0.d0)
      hi(:,:)=CMPLX(0.d0,0.d0)


      !=======================================================================
      !                 begin of the main loop
      !=======================================================================
      !
      OUTER_LOOP: &
      do while ( itercg < maxiter .and. (.not. ltresh) )
        !
        call start_clock( "outer_loop" )

!$$$$
!$$$$        if(itercg.ge.10) do_innerloop=.false.
!$$$$

!$$
#ifdef __DEBUG
        if( do_orbdep .and. ionode .and.( itercg == 1) ) then

          open(1032,file='convg_outer.dat',status='unknown')
          write(1032,'("#   ninner    nouter     non-sic energy (Ha)         sic energy (Ha)")')

          if(do_innerloop) then
            open(1031,file='convg_inner.dat',status='unknown')
            write(1031,'("#   ninner    nouter     non-sic energy (Ha)         sic energy (Ha)    RMS force eigenvalue")')
          endif
        endif
#endif
!$$

        ENERGY_CHECK: &
        if(.not. ene_ok ) then

          call calbec(1,nsp,eigr,c0,bec)
          IF(non_ortho) THEN
             call calbec(1,nsp,eigr,cdual,becdual)
          ENDIF
          
          if(.not.tens) then
             !
             if(non_ortho) then
                call rhoofr(nfi,c0(:,:),cdual,irb,eigrb,bec,becdual,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
             else
                call rhoofr(nfi,c0(:,:),irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
             endif
          else

            if(newscheme.or.firstiter) then 
               call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,c0,bec,firstiter,vpot)
               firstiter=.false.
            endif
            !     calculation of the rotated quantities

            call rotate_twin( z0t, c0(:,:), bec, c0diag, becdiag, .false. )
            !     calculation of rho corresponding to the rotated wavefunctions
            call rhoofr(nfi,c0diag,irb,eigrb,becdiag                        &
                     &                    ,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
          endif
           
          !
          ! when cycle is restarted go to diagonal representation
          !
          ! CHP: do we need to do the following even if when we do not use ensemble dft?
          !      I have added this additional constraint.
          !
          if( tens .and. mod(itercg,niter_cg_restart) ==1 .and. itercg >= 2 ) then
              !
              call rotate_twin( z0t, c0(:,:), bec, c0diag, becdiag, .false. )
              c0(:,:)=c0diag(:,:)
              call copy_twin(bec,becdiag) !modified:giovanni
!               bec(:,:)=becdiag(:,:)
              !
              call id_matrix_init( descla, nspin )
              !
          endif
        
          !
          ! calculates the potential
          !
          !     put core charge (if present) in rhoc(r)
          !
          if (nlcc_any) call set_cc(irb,eigrb,rhoc)

          !
          !---ensemble-DFT

          vpot = rhor

!$$
          CALL start_clock( 'vofrho1' )
!$$
          call vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast,             &
                 &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)

          ene_lda = etot
!$$
          CALL stop_clock( 'vofrho1' )
!$$

!$$
          if( do_orbdep ) then
              !
              if ( tens .or. tsmear) then
                  fsic = fmat0_diag
              else
                  fsic = f
              endif
              !
              IF(non_ortho) THEN
                 call nksic_potential_non_ortho( nbsp, nbspx, c0, cdual, fsic, bec, becdual, &
                                    rhovan, deeq_sic, &
                                    ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, wfc_centers, &
                                    wfc_spreads, icompute_spread, .false. )
              ELSE
                 !
!                  write(6,*) "checkbounds", ubound(c0)
!                  write(6,*) "checkbounds", ubound(fsic)
!                  write(6,*) "checkbounds", ubound(rhovan)
!                  write(6,*) "checkbounds", ubound(wfc_centers), nudx
!                  write(6,*) "checkbounds", ubound(wfc_spreads)
!                  write(6,*) "checkbounds", ubound(deeq_sic)
!                  write(6,*) "checkbounds", ubound(ispin)
!                  write(6,*) "checkbounds", ubound(iupdwn)
!                  write(6,*) "checkbounds", ubound(nupdwn)
!                  write(6,*) "checkbounds", ubound(rhor), "rhor"
!                  write(6,*) "checkbounds", ubound(rhog), "rhog"
!                  write(6,*) "checkbounds", ubound(wtot), "wtot"
!                  write(6,*) "checkbounds", ubound(vsic), "vsic"
!                  write(6,*) "checkbounds", ubound(pink), "pink", nudx
                 !
                 call nksic_potential( nbsp, nbspx, c0, fsic, bec, rhovan, deeq_sic, &
                                    ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, wfc_centers, &
                                    wfc_spreads, icompute_spread, .false. )
!                  write(6,*) "checkbounds", ubound(c0)
!                  write(6,*) "checkbounds", ubound(fsic)
!                  write(6,*) "checkbounds", ubound(rhovan)
!                  write(6,*) "checkbounds", ubound(wfc_centers), nudx
!                  write(6,*) "checkbounds", ubound(wfc_spreads)
!                  write(6,*) "checkbounds", ubound(deeq_sic)
!                  write(6,*) "checkbounds", ubound(ispin)
!                  write(6,*) "checkbounds", ubound(iupdwn)
!                  write(6,*) "checkbounds", ubound(nupdwn)
!                  write(6,*) "checkbounds", ubound(rhor), "rhor"
!                  write(6,*) "checkbounds", ubound(rhog), "rhog"
!                  write(6,*) "checkbounds", ubound(wtot), "wtot"
!                  write(6,*) "checkbounds", ubound(vsic), "vsic"
!                  write(6,*) "checkbounds", ubound(pink), "pink", nudx
              ENDIF

              eodd=sum(pink(1:nbsp))
!               write(6,*) eodd, etot, "EODD0", etot+eodd
              etot = etot + eodd
              !
          endif
!$$
          if( do_hf ) then
              !
              call hf_potential( nbsp, nbspx, c0, f, ispin, iupdwn, nupdwn, &
                                 nbsp, nbspx, c0, f, ispin, iupdwn, nupdwn, &
                                 rhor, rhog, vxxpsi, exx)
              !
              etot = etot + sum(exx(1:nbsp))
              !
          endif

          if (.not.tens) then
              etotnew=etot
          else
              etotnew=etot+entropy
          end if

          if(tefield  ) then!just in this case calculates elfield stuff at zeo field-->to be bettered
            
             call berry_energy( enb, enbi, bec%rvec, c0(:,:), fion )
             etot=etot+enb+enbi
          endif
          if(tefield2  ) then!just in this case calculates elfield stuff at zeo field-->to be bettered

             call berry_energy2( enb, enbi, bec%rvec, c0(:,:), fion )
             etot=etot+enb+enbi
          endif

        else

          etot=enever
          if(.not.tens) then 
             etotnew=etot
          else
             etotnew=etot+entropy
          endif
          ene_ok=.false.

        end if ENERGY_CHECK

!$$
        if( do_orbdep ) then

#ifdef __DEBUG
          if( ionode .and. itercg == 1 ) then
             write(1032,'(2I10,2F24.13)') 0,0,etot-eodd,eodd
          endif
#endif

          if(do_innerloop) then
!$$$$          if(do_innerloop.and.itercg.le.20) then
!$$$$
             !
             !call start_clock( "inner_loop" )
             !
             eodd    = sum(pink(1:nbsp))
             etot    = etot - eodd
             etotnew = etotnew - eodd
             ninner  = 0

             if(.not.do_innerloop_cg) then
                 call nksic_rot_emin(itercg,ninner,etot,Omattot, lgam)
             else
                 call nksic_rot_emin_cg(itercg,innerloop_init_n,ninner,etot,Omattot,deltae*innerloop_cg_ratio,lgam)
             endif

!$$ Now rotate hi(:,:) according to Omattot!
!$$ It seems that not rotating hi gives us better convergence.
!$$ So, we do not perform the following routine.
!$$
!            if(ninner.ge.2) then
!              hitmp(:,:) = CMPLX(0.d0,0.d0)
!              do nbnd1=1,n
!                do nbnd2=1,n
!                  hitmp(:,nbnd1)=hitmp(:,nbnd1) + hi(:,nbnd2) * Omattot(nbnd2,nbnd1)
!                enddo
!              enddo
!              hi(:,:) = hitmp(:,:)
!            endif
!$$
             eodd    = sum(pink(1:nbsp))
             etot    = etot + eodd
             etotnew = etotnew + eodd
             eoddnew = eodd

             !call stop_clock( "inner_loop" )
             !
           endif
           !
        endif
!$$

!$$     
#ifdef __DEBUG
        ! for debug and tuning purposes
        if ( ionode ) write(37,*)itercg, itercgeff, etotnew
        if ( ionode ) write(1037,'("iteration =",I4,"  eff iteration =",I4,"   Etot (Ha) =",F22.14)')&
            itercg, itercgeff, etotnew 
#endif
        if ( ionode ) write(stdout,'(5x,"iteration =",I4,"  eff iteration =",I4,"   Etot (Ha) =",F22.14," delta_E=",E22.14)')&
            itercg, itercgeff, etotnew, deltae

        if ( ionode .and. mod(itercg,10) == 0 ) write(stdout,"()" )

        IF(iprint_spreads>0) THEN
           !
           IF ( ionode .and. mod(itercg, iprint_spreads)==0) THEN
              !
              IF(nspin==1) THEN
                 !
                 WRITE( iunit_spreads, '(400f20.14)') wfc_spreads(:,1,2)               
                 !
              ELSE IF(nspin==2) THEN
                 !
                 WRITE( iunit_spreads, '(2(400f20.14)(3x))') wfc_spreads(:,1,2), wfc_spreads(:,2,2)
                 !
              ENDIF
              !
           ENDIF
           !
        ENDIF

        IF(iprint_manifold_overlap>0) THEN
           !
           IF(mod(itercg, iprint_manifold_overlap)==0) THEN
              !
              CALL compute_manifold_overlap( cstart, c0, becstart, bec, ngwx, nbspx, manifold_overlap ) !added:giovanni
              !
              IF ( ionode ) THEN
                 !
                 IF(nspin==1) THEN
                    !
                    WRITE( iunit_manifold_overlap, '(2f20.14)') manifold_overlap(1)
                    !
                 ELSE IF(nspin==2) THEN
                    !
                    WRITE( iunit_manifold_overlap, '(2(2f20.14)(3x))') manifold_overlap(1), manifold_overlap(2)
                    !
                 ENDIF
                 !
              ENDIF
              !
           ENDIF
           !
        ENDIF

!$$


!$$ to see the outer loop energy convergence
        if (do_orbdep) then
            !
            eodd = sum(pink(1:nbsp))
#ifdef __DEBUG
            if(ionode) write(1032,'(2I10,2F24.13)') ninner,itercg,etot-eodd,eodd
#endif
            !
        endif
!$$
        deltae=abs(etotnew-etotold) 
        if(do_orbdep) then
           deltae=deltae+abs(eoddnew-eoddold)
        endif
        !
        if( deltae < conv_thr ) then
           numok=numok+1
        else 
           numok=0
        endif
        !
        if( numok >= 4 ) ltresh=.true.
        !
        if(ltresh.or.itercg==maxiter-1) icompute_spread=.true.
        !
        etotold=etotnew
        eoddold=eoddnew
        ene0=etot
        !
        if( tens .and. newscheme ) ene0=ene0+entropy
        
        

!$$$$ For a test: Calculates wavefunctions very close to c0.
!    if(.false.) then
!      cm(1:ngw,1:n)=c0(1:ngw,1:n)
!      if(ng0.eq.2) then
!        cm(1,:)=0.5d0*(cm(1,:)+CONJG(cm(1,:)))
!      endif
!
!      call lowdin(cm)
!      call calbec(1,nsp,eigr,cm,becm)
!
!      call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
!      vpot = rhor
!
!      call vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast,             &
!                  &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
!
!      if(do_orbdep) then
!        call nksic_potential( n, nbspx, cm, fsic, bec, rhovan, deeq_sic, &
!                 ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )
!        etot = etot + sum(pink(:))
!      endif
!
!      ene0 = etot
!    endif
!$$$$


!!$$$$ For a test: Calculates wavefunction very close to c0.
!    if(.false.) then
!      !
!#ifdef __DEBUG
!      if(ionode) write(1000,*) 'Now entering the routine...'
!      if(ionode) write(1000,*) itercg
!#endif
!      !
!      cm(1:ngw,1:nbsp)=c0(1:ngw,1:nbsp)
!      if(ng0.eq.2) then
!        cm(1,:)=0.5d0*(cm(1,:)+CONJG(cm(1,:)))
!      endif
!
!      call lowdin(cm)
!      call calbec(1,nsp,eigr,cm,becm)
!!      becm=bec
!
!      call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
!      vpot = rhor
!
!      call vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast,             &
!                  &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
!
!      ene_save2(1)=etot
!
!      if(do_orbdep) then
!          !
!          call nksic_potential( nbsp, nbspx, cm, fsic, bec, rhovan, deeq_sic, &
!                                ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )
!          etot = etot + sum(pink(1:nbsp))
!          !
!      endif
!      !
!      if( do_hf ) then
!          !
!          call hf_potential( nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
!                             nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
!                             rhor, rhog, vxxpsi, exx)
!          !
!          etot = etot + sum(exx(1:nbsp))
!          !
!      endif
!
!#ifdef __DEBUG
!      if(ionode) then
!          write(1000,'(3e30.20)')  ene0,etot,etot-ene0
!          write(1000,'(3e30.20)')  esic,sum(pink(:)), sum(pink(:))-esic
!          write(1000,*)
!      endif
!#endif
!      !
!    endif
!!$$$$


        !update d

        call newd(vpot,irb,eigrb,rhovan,fion)


        call prefor(eigr,betae)!ATTENZIONE

!$$
        ! faux takes into account spin multiplicity.
        !
        faux(1:nbspx)=0.d0
        faux(1:nbsp) = max(f_cutoff,f(1:nbsp)) * DBLE( nspin ) / 2.0d0
!$$
        IF(non_ortho) THEN
           !
           hpsi0(:,:) = CMPLX(0.d0,0.d0)
           hitmp(:,:) = CMPLX(0.d0,0.d0)
           !
        ENDIF
        !
        do i=1,nbsp,2
!$$  FIRST CALL TO DFORCE
          CALL start_clock( 'dforce1' )
!$$          call dforce( i, bec, betae, c0,c2,c3,rhos, nnrsx, ispin,f,n,nspin)
          IF(non_ortho) THEN
             call dforce( i, bec, betae, c0,c2,c3,rhos, nnrsx, ispin, faux, nbsp, nspin)
          ELSE
             call dforce( i, bec, betae, c0,c2,c3,rhos, nnrsx, ispin, faux, nbsp, nspin)
          ENDIF
          !
          CALL stop_clock( 'dforce1' )
!$$

          if(tefield .and. (evalue.ne.0.d0)) then
            call dforceb(c0, i, betae, ipolp, bec%rvec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            c2(1:ngw)=c2(1:ngw)+evalue*df(1:ngw)
            call dforceb(c0, i+1, betae, ipolp, bec%rvec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            c3(1:ngw)=c3(1:ngw)+evalue*df(1:ngw)
          endif
          !
          if(tefield2 .and. (evalue2.ne.0.d0)) then
            call dforceb(c0, i, betae, ipolp2, bec%rvec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
            c2(1:ngw)=c2(1:ngw)+evalue2*df(1:ngw)
            call dforceb(c0, i+1, betae, ipolp2, bec%rvec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
            c3(1:ngw)=c3(1:ngw)+evalue2*df(1:ngw)
          endif
          
          IF ( lda_plus_u ) THEN
               !
               IF ( tens .OR. tsmear ) THEN
                   !
                   c2(:) = c2(:) - vupsi(:,i) 
                   if( i+1 <= nbsp ) c3(:) = c3(:) - vupsi(:,i+1) 
                   !
               ELSE
                   !
                   c2(:) = c2(:) - vupsi(:,i) * faux(i)
                   if( i+1 <= nbsp ) c3(:) = c3(:) - vupsi(:,i+1) * faux(i+1)
                   !
               ENDIF
               !
           ENDIF

!$$
!          hpsinosic(1:ngw,  i)=c2(1:ngw)
!          if(i+1 <= n) then
!            hpsinosic(1:ngw,i+1)=c3(1:ngw)
!          endif
!          if (ng0.eq.2) then
!            hpsinosic(1,  i)=CMPLX(DBLE(hpsinosic(1,  i)), 0.d0)
!            if(i+1 <= n) then
!              hpsinosic(1,i+1)=CMPLX(DBLE(hpsinosic(1,i+1)), 0.d0)
!            endif
!          end if
!$$

!$$
          if ( do_orbdep ) then
              !
              ! faux takes into account spin multiplicity.
              !
              IF(non_ortho) THEN
                 !
                 IF(northo_flavor==2) THEN
                    CALL compute_overlap( c0, ngw, nbsp, overlap)
                    !FLAVOUR2_NONORTHO_SIC
                 ENDIF
                 !
                 CALL nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, bec, ngw, c0(:,i), c0(:,i+1), vsicpsi, lgam )
                 !
              ELSE
                 CALL nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, bec, ngw, c0(:,i), c0(:,i+1), vsicpsi, lgam )
              ENDIF

              IF(non_ortho) THEN
                 !
                 IF(northo_flavor==1) THEN
                    !
                    c2(:) = c2(:) - 0.5d0 * vsicpsi(:,1) * faux(i)
                    !
                    if( i+1 <= nbsp )   c3(:) = c3(:) - 0.5d0 * vsicpsi(:,2) * faux(i+1)
                    !
                    CALL nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, becdual, ngw, cdual(:,i), cdual(:,i+1), vsicpsi, lgam )
                    !
                 ELSE IF(northo_flavor==2) THEN
                    !
                    vsicpsi(:,1) = vsicpsi(:,1)/overlap(i-iupdwn(ispin(i))+1,i-iupdwn(ispin(i))+1,ispin(i))!FLAVOUR2_NONORTHO_SIC
                    !
                 ENDIF
                 !
                 hitmp(:,i) = hitmp(:,i) - sic_coeff1 * vsicpsi(:,1) * faux(i)
                 !
                 if( i+1 <= nbsp ) then
                    !
                    IF(northo_flavor==2) THEN
                       !
                       vsicpsi(:,2) = vsicpsi(:,2)/overlap(i-iupdwn(ispin(i+1))+2,i-iupdwn(ispin(i+1))+2,ispin(i+1))!FLAVOUR2_NONORTHO_SIC
                       !
                    ENDIF
                    !
                    hitmp(:,i+1) = hitmp(:,i+1) - sic_coeff1 * vsicpsi(:,2) * faux(i+1)
                    !
                 endif
                 !
              ELSE
                 !
                 c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
                 !
                 if( i+1 <= nbsp )   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
                 !
              ENDIF
              !
          endif
!$$
          if ( do_hf ) then
              !
              c2(:) = c2(:) - vxxpsi(:,i) * faux(i)
              !
              if( i+1 <= nbsp )   c3(:) = c3(:) - vxxpsi(:,i+1) * faux(i+1)
              !
          endif


          hpsi(1:ngw,  i)=c2(1:ngw)
          if(i+1 <= nbsp) then
              hpsi(1:ngw,i+1)=c3(1:ngw)
          endif
          !
          IF(lgam) THEN
	    if (ng0.eq.2) then
		hpsi(1,  i)=CMPLX(DBLE(hpsi(1,  i)), 0.d0)
		if(i+1 <= nbsp) then
		    hpsi(1,i+1)=CMPLX(DBLE(hpsi(1,i+1)), 0.d0)
		endif
	    endif
          ENDIF
          !
        enddo
        
        IF(non_ortho) THEN
           !
           IF(northo_flavor==1) THEN
              !
              call times_overlap(c0, hitmp, gi, nbsp, 1)
              hpsi(:,:) = hpsi(:,:)+gi(:,:)
              call pc3nc_non_ortho(c0,cdual,hpsi, lgam)
              !
           ELSE IF(northo_flavor==2) THEN
              !
              call pc4nc_non_ortho(gi,c0,c0,hitmp,lgam)
              !
              do i=1,nbsp
                 !
                 gi(:,i) = gi(:,i)/overlap(i-iupdwn(ispin(i))+1,i-iupdwn(ispin(i))+1,ispin(i))
                 !
              enddo
              !
              hitmp(:,:) = hitmp(:,:)+gi(:,:)
              call times_overlap(c0, hitmp, hpsi0, nbsp, 1)
!             !
              call pc3nc_non_ortho(c0,cdual,hpsi, lgam)
              !
              hpsi(:,:) = hpsi(:,:)+hpsi0(:,:)
           ENDIF
           !
        ENDIF
!$$
!        if(.not.tens) then
!          do i=1,n
!            hpsinorm(i) = 0.d0
!            hpsinosicnorm(i) = 0.d0
!            do ig=1,ngw
!              hpsinorm(i)=hpsinorm(i)+DBLE(CONJG(hpsi(ig,i))*hpsi(ig,i))
!              hpsinosicnorm(i)=hpsinosicnorm(i)+DBLE(CONJG(hpsinosic(ig,i))*hpsinosic(ig,i))
!            enddo
!          end do
!        endif
!        call mp_sum(hpsinorm(1:n),intra_image_comm)
!        call mp_sum(hpsinosicnorm(1:n),intra_image_comm)
!        if(ionode) write(100,*) 'hpsinorm is ',(hpsinorm(i),i=1,n)
!        if(ionode) write(100,*) 'hpsinosicnorm is ',(hpsinosicnorm(i),i=1,n)
!$$

        if(pre_state) call ave_kin(c0,SIZE(c0,1),nbsp,ave_ene)

!$$        call pcdaga2(c0,phi,hpsi)
!$$     HPSI IS ORTHOGONALIZED TO  c0
        IF(non_ortho) THEN
!            call pcdaga2(c0,cdual,hpsi, lgam)
!            call pc3nc_non_ortho(c0,cdual,hpsi, lgam)
!            call pc2_non_ortho(c0,cdual,bec,becdual,hpsi,becm,lgam)
        ELSE
           if(switch.or.(.not.do_orbdep)) then
             call pcdaga2(c0,phi,hpsi, lgam)
           else
!           call calbec(1,nsp,eigr,hpsi,becm)
            call pc3nc(c0,hpsi,lgam)
!           call pc3us(c0,bec,hpsi,becm,lgam)
!           call pcdaga3(c0,phi,hpsi, lgam)
           endif
           !
        ENDIF
!$$

!begin_added:giovanni debug, check orthonormality
!        temp=0.d0
!        do ig=1,ngw
!        temp=temp+2.d0*DBLE(CONJG(c0(ig,1)+hpsi(ig,1))*(c0(ig,1)+hpsi(ig,1)))
!        enddo
!        if(ng0==2.and.lgam) then
!        temp=temp-DBLE(CONJG((c0(1,1)+hpsi(1,1)))*(c0(1,1)+hpsi(1,1)))
!        endif
!        call mp_sum(temp,intra_image_comm)
!        write(6,*) "debug", temp
!end_added:giovanni

!$$
!        if(ionode) then
!          do i=1,n
!            write(701,*) sum(phi(1:ngw,i)),sum(c0(1:ngw,i))
!          enddo
!          write(701,*) 'nhsa ',nhsa
!          write(701,*)
!        endif
!$$
        !TWO VECTORS INITIALIZED TO HPSI
        IF(non_ortho) THEN

!            call times_overlap(c0, hpsi, hpsi0, nbsp, 1)
!            gi(1:ngw,1:nbsp) = hpsi(1:ngw,1:nbsp)
!            hpsi(1:ngw,1:nbsp)    = hpsi0(1:ngw,1:nbsp)
!            call pc3nc_non_ortho(c0,cdual,hpsi,lgam)          
           call times_overlap(c0, hpsi, hpsi0, nbsp, -1)
           gi(1:ngw,1:nbsp)    = hpsi0(1:ngw,1:nbsp)
           hpsi(1:ngw,1:nbsp)    = hpsi0(1:ngw,1:nbsp)
        ELSE
           hpsi0(1:ngw,1:nbsp) = hpsi(1:ngw,1:nbsp)
           gi(1:ngw,1:nbsp)    = hpsi(1:ngw,1:nbsp)
        ENDIF

	!COMPUTES ULTRASOFT-PRECONDITIONED HPSI, non kinetic-preconditioned, is the subsequent reorthogonalization necessary in the norm conserving case???: giovanni
        call calbec(1,nsp,eigr,hpsi,becm)
        IF(non_ortho) THEN
        !
        ELSE
           call xminus1_twin(hpsi,betae,dumm,becm,s_minus1,.false.)
        ENDIF
!        call sminus1(hpsi,becm,betae)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!look if the following two lines are really needed
        call calbec(1,nsp,eigr,hpsi,becm)
!$$        call pc2(c0,bec,hpsi,becm)
!$$     THIS ORTHOGONALIZED PRECONDITIONED VECTOR HPSI

        IF(non_ortho) THEN
!             call pc2_non_ortho(c0,cdual,bec,becdual,hpsi,becm,lgam)
!             call pc3nc_non_ortho(cdual,c0,hpsi,lgam)
        ELSE
           if(switch.or.(.not.do_orbdep)) then
             call pc2(c0,bec,hpsi,becm, lgam)
           else
             call pc3nc(c0,hpsi,lgam)
!           call pc3us(c0,bec,hpsi,becm, lgam)
           endif
        ENDIF
!$$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       COMPUTES ULTRASOFT+KINETIC-preconditioned GI
!        call kminus1(gi,betae,ema0bg)
        IF(do_orbdep) THEN ! preconditioning with respect to spreads.. the gradient along wavefunctions with the largest localization is a bit shortened
           !
           !do i=1,nbsp
           !   gi(:,i) = gi(:,i)*(1.d0+1.d0/sqrt(wfc_spreads(1+i-iupdwn(ispin(i)),ispin(i),2)))
           !enddo
           !
        ENDIF
        !
        IF(non_ortho) THEN
           call calbec(1,nsp,eigr,gi,becm)
           call xminus1_twin(gi,betae,ema0bg,becm,k_minus1,.true.)
           call times_overlap(c0, gi, cmdual, nbsp, 1)
           gi(:,:)=cmdual(:,:)
        ELSE
           if(.not.pre_state) then
              call xminus1_twin(gi,betae,ema0bg,becm,k_minus1,.true.)
           else
              call xminus1_state(gi,betae,ema0bg,becm,k_minus1,.true.,ave_ene) !warning:giovanni not yet implemented
           endif
        ENDIF
        !
        call calbec(1,nsp,eigr,gi,becm)
!$$        call pc2(c0,bec,gi,becm)
!$$     !ORTHOGONALIZES GI to c0
        IF(non_ortho) THEN
!            call pc2_non_ortho(c0,cdual,bec,becdual,gi,becm,lgam)
!            call pc3nc_non_ortho(c0,cdual,gi,lgam)
        ELSE
           if(switch.or.(.not.do_orbdep)) then
             call pc2(c0,bec,gi,becm, lgam)
           else
             call pc3nc(c0,gi, lgam)
!           call pc3us(c0,bec, gi,becm, lgam)
           endif
        ENDIF
!$$
        if(tens) call calcmt_twin( f, z0t, fmat0, .false. )
        call calbec(1,nsp,eigr,hpsi,bec0) 
!  calculates gamma
        gamma_c=CMPLX(0.d0,0.d0)
        
        if(.not.tens) then
           do i=1,nbsp
              IF(lgam) THEN
		do ig=1,ngw
		  gamma_c=gamma_c+2.d0*DBLE(CONJG(gi(ig,i))*hpsi(ig,i))
		enddo
		if (ng0.eq.2) then
		  gamma_c=gamma_c-DBLE(CONJG(gi(1,i))*hpsi(1,i))
		endif
              ELSE
		do ig=1,ngw
		  gamma_c=gamma_c+CONJG(gi(ig,i))*hpsi(ig,i)
		enddo
              ENDIF
           enddo
           

	  call mp_sum( gamma_c, intra_image_comm )
           
	   if (nvb.gt.0) then
            if(.not.becm%iscmplx) then
		do i=1,nbsp
		  do is=1,nvb
		      do iv=1,nh(is)
			do jv=1,nh(is)
			    do ia=1,na(is)
			      inl=ish(is)+(iv-1)*na(is)+ia
			      jnl=ish(is)+(jv-1)*na(is)+ia
			      gamma_c=gamma_c+ qq(iv,jv,is)*becm%rvec(inl,i)*bec0%rvec(jnl,i)
			    end do
			end do
		      end do
		  end do
		enddo
            else
		do i=1,nbsp
		  do is=1,nvb
		      do iv=1,nh(is)
			do jv=1,nh(is)
			    do ia=1,na(is)
			      inl=ish(is)+(iv-1)*na(is)+ia
			      jnl=ish(is)+(jv-1)*na(is)+ia
			      gamma_c=gamma_c+ qq(iv,jv,is)*CONJG(becm%cvec(inl,i))*(bec0%cvec(jnl,i)) !warning:giovanni CONJG
			    end do
			end do
		      end do
		  end do
		enddo
	    endif
          endif

        else

           do iss=1,nspin
              nss=nupdwn(iss)
              istart=iupdwn(iss)
              me_rot = descla( la_me_ , iss )
              np_rot = descla( la_npc_ , iss ) * descla( la_npr_ , iss )
              if(.not.fmat0(iss)%iscmplx) then
		allocate( fmat_ ( nrlx, nudx ) )
		do ip = 1, np_rot
		  if( me_rot == ( ip - 1 ) ) then
		      fmat_ = fmat0(iss)%rvec(:,:)
		  end if
		  nrl = ldim_cyclic( nss, np_rot, ip - 1 )
		  CALL mp_bcast( fmat_ , ip - 1 , intra_image_comm )
		  do i=1,nss
		      jj = ip
		      do j=1,nrl
			do ig=1,ngw
			    gamma_c=gamma_c+2.d0*DBLE(CONJG(gi(ig,i+istart-1))*hpsi(ig,jj+istart-1))*fmat_(j,i)
			enddo
			if (ng0.eq.2) then
			    gamma_c=gamma_c-DBLE(CONJG(gi(1,i+istart-1))*hpsi(1,jj+istart-1))*fmat_(j,i)
			endif
			jj = jj + np_rot
		      enddo
		  enddo
		enddo
		deallocate( fmat_ )
              else
		allocate( fmat_c_ ( nrlx, nudx ) ) !warning:giovanni, need to put some conjugates somewhere?
		do ip = 1, np_rot
		  if( me_rot == ( ip - 1 ) ) then
		      fmat_c_ = fmat0(iss)%cvec(:,:)
		  end if
		  nrl = ldim_cyclic( nss, np_rot, ip - 1 )
		  CALL mp_bcast( fmat_c_ , ip - 1 , intra_image_comm )
		  do i=1,nss
		      jj = ip
		      do j=1,nrl
			do ig=1,ngw
			    gamma_c=gamma_c+CONJG(gi(ig,i+istart-1))*hpsi(ig,jj+istart-1)*fmat_c_(j,i)
			enddo
			jj = jj + np_rot
		      enddo
		  enddo
		enddo
		deallocate( fmat_c_ )
              endif
           enddo
           if(nvb.gt.0) then
              do iss=1,nspin
                 nss=nupdwn(iss)
                 istart=iupdwn(iss)
                 me_rot = descla( la_me_ , iss )
                 np_rot = descla( la_npc_ , iss ) * descla( la_npr_ , iss )
                 if(.not.fmat0(iss)%iscmplx) then
		  allocate( fmat_ ( nrlx, nudx ) )
		  do ip = 1, np_rot
		      if( me_rot == ( ip - 1 ) ) then
			fmat_ = fmat0(iss)%rvec(:,:)
		      end if
		      nrl = ldim_cyclic( nss, np_rot, ip - 1 )
		      CALL mp_bcast( fmat_ , ip - 1 , intra_image_comm )

		      do i=1,nss
			jj = ip 
			do j=1,nrl
			    do is=1,nvb
			      do iv=1,nh(is)
				  do jv=1,nh(is)
				    do ia=1,na(is)
					inl=ish(is)+(iv-1)*na(is)+ia
					jnl=ish(is)+(jv-1)*na(is)+ia
					gamma_c=gamma_c+ qq(iv,jv,is)*becm%rvec(inl,i+istart-1)*bec0%rvec(jnl,jj+istart-1)*fmat_(j,i)
				    end do
				  end do
			      end do
			    enddo
			    jj = jj + np_rot
			enddo
		      enddo
		  end do
		  deallocate( fmat_ )
                 else
		  allocate( fmat_c_ ( nrlx, nudx ) )
		  do ip = 1, np_rot
		      if( me_rot == ( ip - 1 ) ) then
			fmat_c_(:,:) = fmat0(iss)%cvec(:,:)
		      end if
		      nrl = ldim_cyclic( nss, np_rot, ip - 1 )
		      CALL mp_bcast( fmat_c_ , ip - 1 , intra_image_comm )

		      do i=1,nss
			jj = ip 
			do j=1,nrl
			    do is=1,nvb
			      do iv=1,nh(is)
				  do jv=1,nh(is)
				    do ia=1,na(is)
					inl=ish(is)+(iv-1)*na(is)+ia
					jnl=ish(is)+(jv-1)*na(is)+ia
					gamma_c=gamma_c+ qq(iv,jv,is)*CONJG(becm%cvec(inl,i+istart-1)) &
                                  *(bec0%cvec(jnl,jj+istart-1))*fmat_c_(j,i)
				    end do
				  end do
			      end do
			    enddo
			    jj = jj + np_rot
			enddo
		      enddo
		  end do
		  deallocate( fmat_c_ )
                 endif
              enddo
           endif
	    call mp_sum( gamma_c, intra_image_comm )
        endif
        !case of first iteration

! 	IF(lgam) THEN
	  gamma_c=CMPLX(DBLE(gamma_c),0.d0)
! 	ENDIF

!$$        if(itercg==1.or.(mod(itercg,niter_cg_restart).eq.1).or.restartcg) then
        if( itercg==1 .or. mod(itercg,niter_cg_restart)==0 .or. restartcg) then
!$$

          restartcg=.false.
!$$  We do not have to reset passof every exception of CG!
!warning:giovanni if we do not reset we may have fake convergences!!!!
          passof=passop

          hi(1:ngw,1:nbsp)=gi(1:ngw,1:nbsp)!hi is the search direction

          esse_c=gamma_c

        else

          !find direction hi for general case 
          !calculates gamma for general case, not using Polak Ribiere


          essenew_c=gamma_c
          gamma_c=gamma_c/esse_c
          esse_c=essenew_c

          hi(1:ngw,1:nbsp)=gi(1:ngw,1:nbsp)+(gamma_c)*hi(1:ngw,1:nbsp)

        endif
!note that hi, is saved  on gi, because we need it before projection on conduction states

        !find minimum along direction hi:

        !project hi on conduction sub-space

        call calbec(1,nsp,eigr,hi,bec0)
!$$        call pc2(c0,bec,hi,bec0)
!$$
        IF(non_ortho) THEN
!            call pc2_non_ortho(c0, cdual, bec, becdual, hi, bec0, lgam)
!            call pc3nc_non_ortho(c0,cdual,hi, lgam)
        ELSE
           if(switch.or.(.not.do_orbdep)) then
              call pc2(c0,bec,hi,bec0, lgam)
           else
              call pc3nc(c0,hi,lgam)
!           call pc3us(c0,bec,hi,bec0, lgam)
           endif
        ENDIF
!$$
        !do quadratic minimization
        !             
        !calculate derivative with respect to  lambda along direction hi

        dene0=0.
        if(.not.tens) then
          do i=1,nbsp
            IF(lgam) THEN              
	      do ig=1,ngw
		dene0=dene0-4.d0*DBLE(CONJG(hi(ig,i))*hpsi0(ig,i))
	      enddo
	      if (ng0.eq.2) then
		dene0=dene0+2.d0*DBLE(CONJG(hi(1,i))*hpsi0(1,i))
	      endif
            ELSE
	      do ig=1,ngw
		dene0=dene0-2.d0*DBLE(CONJG(hi(ig,i))*hpsi0(ig,i))
	      enddo
            ENDIF
          end do
!$$ We need the following because n for spin 2 is double that for spin 1!
          dene0 = dene0 *2.d0/nspin
!$$          dene0 = dene0 *4.d0/nspin
!$$
        else
          !in the ensamble case the derivative is Sum_ij (<hi|H|Psi_j>+ <Psi_i|H|hj>)*f_ji
          !     calculation of the kinetic energy x=xmin    
         call calcmt_twin( f, z0t, fmat0, .false. )
         do iss = 1, nspin
            nss    = nupdwn(iss)
            istart = iupdwn(iss)!warning:giovanni this is a bug for a fully polarized system
            me_rot = descla( la_me_ , iss )
            np_rot = descla( la_npc_ , iss ) * descla( la_npr_ , iss )
            if(.not. fmat0(iss)%iscmplx) then
	      allocate( fmat_ ( nrlx, nudx ) )
	      do ip = 1, np_rot
		if( me_rot == ( ip - 1 ) ) then
		    fmat_(:,:) = fmat0(iss)%rvec(:,:)
		end if
		nrl = ldim_cyclic( nss, np_rot, ip - 1 )
		CALL mp_bcast( fmat_ , ip - 1 , intra_image_comm )
		do i=1,nss
		    jj = ip
		    do j=1,nrl
		      do ig=1,ngw
			  dene0=dene0-2.d0*DBLE(CONJG(hi(ig,i+istart-1))*hpsi0(ig,jj+istart-1))*fmat_(j,i)
			  dene0=dene0-2.d0*DBLE(CONJG(hpsi0(ig,i+istart-1))*hi(ig,jj+istart-1))*fmat_(j,i)
		      enddo
		      if (ng0.eq.2) then
			  dene0=dene0+DBLE(CONJG(hi(1,i+istart-1))*hpsi0(1,jj+istart-1))*fmat_(j,i)
			  dene0=dene0+DBLE(CONJG(hpsi0(1,i+istart-1))*hi(1,jj+istart-1))*fmat_(j,i)
		      end if
		      jj = jj + np_rot
		    enddo
		enddo
	      end do
	      deallocate( fmat_ )
            else
	      allocate( fmat_c_ ( nrlx, nudx ) )
	      do ip = 1, np_rot
		if( me_rot == ( ip - 1 ) ) then
		    fmat_c_(:,:) = fmat0(iss)%cvec(:,:)
		end if
		nrl = ldim_cyclic( nss, np_rot, ip - 1 )
		CALL mp_bcast( fmat_ , ip - 1 , intra_image_comm )
		do i=1,nss
		    jj = ip
		    do j=1,nrl
		      do ig=1,ngw
			  dene0=dene0-CONJG(hi(ig,i+istart-1))*hpsi0(ig,jj+istart-1)*fmat_c_(j,i)
			  dene0=dene0-CONJG(hpsi0(ig,i+istart-1))*hi(ig,jj+istart-1)*fmat_c_(j,i)
		      enddo
		      jj = jj + np_rot
		    enddo
		enddo
	      end do
	      deallocate( fmat_c_ )
            endif
         enddo
      endif

      call mp_sum( dene0, intra_image_comm )

        !if the derivative is positive, search along opposite direction
      if(dene0.gt.0.d0) then
         spasso=-1.D0
      else
         spasso=1.d0
      endif

!$$$$ Calculates wavefunction at very close to c0.
!    if(.false.) then
!      tmppasso=0.d-8
!      if(ionode) write(8000,*) itercg
!      do i=1,5
!        cm(1:ngw,1:n)=c0(1:ngw,1:n)+spasso * tmppasso * hi(1:ngw,1:n)
!        if(ng0.eq.2) then
!          cm(1,:)=0.5d0*(cm(1,:)+CONJG(cm(1,:)))
!        endif
!
!        call lowdin(cm)
!        call calbec(1,nsp,eigr,cm,becm)
!
!        call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
!        vpot = rhor
!
!        call vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast,             &
!                    &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
!
!        ene_save2(i)=etot
!
!        if(do_orbdep) then
!          call nksic_potential( nbsp, nbspx, cm, fsic, bec, rhovan, deeq_sic, &
!                   ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )
!          etot = etot + sum(pink(:))
!        endif
!
!        if(ionode) then
!          write(8000,'(2e30.20,3e20.10)')  ene0,etot,dene0,tmppasso,((etot-ene0)+1.d-10)/(tmppasso+1.d-10)/dene0
!        endif
!
!        ene_save(i)=etot
!
!        tmppasso=tmppasso+1.d-8
!      enddo
!
!      if(ionode) then
!        write(8000,'(2e30.20,3e20.10)')  ene_save(1),ene_save(2),dene0,1.d-8,(ene_save(2)-ene_save(1))/(1.d-8)/dene0
!        write(8000,*)
!        write(9000,'(3e30.20)')  ene_lda,ene_save2(1), ene_lda-ene_save2(1)
!        write(9000,*)
!      endif
!
!    endif
!$$$$


! open(file="~/marzari/debug.txt", unit=8000)
!$$$$ Calculates wavefunction at very close to c0.
!    if(.true.) then
!      tmppasso=1.d-4
!      !
! #ifdef __DEBUG
!      if(ionode) write(6,*) "debug", itercg
! #endif
!      do i=1,5
!        cm(1:ngw,1:nbsp)=c0(1:ngw,1:nbsp)+spasso * tmppasso * hi(1:ngw,1:nbsp)
!        if(ng0.eq.2) then
!          cm(1,:)=0.5d0*(cm(1,:)+CONJG(cm(1,:)))
!        endif
! 
!        call lowdin(cm)
!        call calbec(1,nsp,eigr,cm,becm)
! 
!        call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
!        vpot = rhor
! 
!        call vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast,             &
!                    &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
! 
!        ene_save2(i)=etot
! 
!        if(do_orbdep) then
!            !
!            call nksic_potential( nbsp, nbspx, cm, fsic, bec, rhovan, deeq_sic, &
!                                  ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )
!            !
!            etot = etot + sum(pink(1:nbsp))
!            !
!        endif
!        !
!        if( do_hf ) then
!            !
!            call hf_potential( nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
!                               nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
!                               rhor, rhog, vxxpsi, exx)
!            !
!            etot = etot + sum(exx(1:nbsp))
!            !
!        endif
!        !
! #ifdef __DEBUG
!        if(ionode) then
!            write(6,'(2e30.20,3e20.10)')  ene0,etot,dene0,tmppasso,(etot-ene0)/tmppasso/dene0
!        endif
! #endif
!        !if(ionode) then
!        !    write(stdout,'(2e30.20,3e20.10)')  ene0,etot,dene0,tmppasso,(etot-ene0)/tmppasso/dene0
!        !endif
! 
!        ene_save(i)=etot
! 
!        tmppasso=tmppasso*0.1d0
!        !
!      enddo
! 
! #ifdef __DEBUG
!      if(ionode) write(6,*) "debug"
! #endif
!      !
!    endif
! close(8000)
!$$$$



      !
      ! calculates wave-functions on a point on direction hi
      !
      cm(1:ngw,1:nbsp)=c0(1:ngw,1:nbsp)+spasso*passof*hi(1:ngw,1:nbsp)
      !
!$$   ! I do not know why the following 3 lines 
      ! were not in the original code (CHP)
      !
      if(lgam.and.ng0 == 2)  cm(1,:)=0.5d0*(cm(1,:)+CONJG(cm(1,:)))
!$$

      !orthonormalize

      !
      IF(non_ortho) THEN
         call calbec(1,nsp,eigr,cm,becm)
         call compute_duals(cm,cmdual,nbspx,1)
         call calbec(1,nsp,eigr,cmdual,becmdual)
         write(6,*) "checkdual", cdual(1:2,1)
      ELSE
         if(do_orbdep.and.ortho_switch) then
            call lowdin(cm, lgam)
            call calbec(1,nsp,eigr,cm,becm)
         else
            call calbec(1,nsp,eigr,cm,becm)
            call gram(betae,becm,nhsa,cm,ngw,nbsp)
         endif
      ENDIF
        !call calbec(1,nsp,eigr,cm,becm)

        !****calculate energy ene1
        if(.not.tens) then
           if(non_ortho) then
              call rhoofr(nfi,cm(:,:),cmdual, irb,eigrb,becm,becmdual,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
           else
              call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
           endif
        else
          if(newscheme) then 
              call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                        rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,.false., vpot  )  
          endif

          !     calculation of the rotated quantities
          call rotate_twin( z0t, cm(:,:), becm, c0diag, becdiag, .false. )
          !     calculation of rho corresponding to the rotated wavefunctions
          call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
        endif

        !calculate potential
        !
        !     put core charge (if present) in rhoc(r)
        !
        if (nlcc_any) call set_cc(irb,eigrb,rhoc)
        !
        vpot = rhor
        !
!$$
!        if(ionode) write(*,*) 'Now doing vofrho2'
        CALL start_clock( 'vofrho2' )
!$$
        call vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast,             &
                      &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
!$$
        CALL stop_clock( 'vofrho2' )
!$$

!$$
        if( do_orbdep ) then
            !warning:giovanni don't we need becm down here??? otherwise problems with ultrasoft!!
            IF(non_ortho) THEN
               call nksic_potential_non_ortho( nbsp, nbspx, cm, cmdual, fsic, bec, becdual, rhovan, deeq_sic, &
                                  ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx,&
                                  wfc_centers, wfc_spreads, &
                                  icompute_spread, .false.)
            ELSE
               call nksic_potential( nbsp, nbspx, cm, fsic, bec, rhovan, deeq_sic, &
                                  ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, &
                                  wfc_centers, wfc_spreads, &
                                  icompute_spread, .false.)
            ENDIF
            !
            eodd=sum(pink(1:nbsp))
!             write(6,*) eodd, etot, "EODD2", etot+eodd !debug:giovanni
            etot = etot + eodd
            !
        endif
!$$
        if( do_hf ) then
            !
            call hf_potential( nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
                               nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
                               rhor, rhog, vxxpsi, exx)
            !
            etot = etot + sum(exx(1:nbsp))
            !
        endif

        if( tefield  ) then!to be bettered
          call berry_energy( enb, enbi, becm, cm(:,:), fion )
          etot=etot+enb+enbi
        endif
        !
        if( tefield2  ) then!to be bettered
          call berry_energy2( enb, enbi, becm, cm(:,:), fion )
          etot=etot+enb+enbi
        endif

        ene1=etot
        if( tens .and. newscheme) ene1=ene1+entropy
              
            
        !find the minimum

        call minparabola(ene0,spasso*dene0,ene1,passof,passo,enesti)

        if( ionode .and. iprsta > 1 ) write(stdout,"(6f20.12)") ene0,dene0,ene1,passo, DBLE(gamma_c), esse_c

        !set new step

        passov=passof
!         passof=2.d0*passo
!$$ doing the following makes the convergence better...
        passof=passo
!$$$$
              
        !calculates wave-functions at minimum

        cm(1:ngw,1:nbsp) = c0(1:ngw,1:nbsp) +spasso*passo*hi(1:ngw,1:nbsp)
        !
        if(lgam.and. ng0 == 2 )  THEN
          cm(1,:) = 0.5d0*(cm(1,:)+CONJG(cm(1,:)))
        ELSE !warning:giovanni this would fix the phase of the new position.. should
             !        not influence the calculation
        ! do i=1,nbsp
        !  phase=0.d0
        !  IF(ng0 == 2 ) THEN
        !   phase = cm(1,i)/(abs(cm(1,i))+1.d-10)
        !  ENDIF
        !  call mp_sum(phase, intra_image_comm)
        !  cm(:,i) = cm(:,i)*CONJG(phase)
        ! enddo
        ENDIF
      
        IF(non_ortho) THEN
           call calbec(1,nsp,eigr,cm,becm)
           call compute_duals(cm,cmdual,nbspx,1)
           call calbec(1,nsp,eigr,cmdual,becmdual)
        ELSE
           IF(do_orbdep.and.ortho_switch) THEN
              call lowdin(cm, lgam)
              call calbec(1,nsp,eigr,cm,becm)
           ELSE
              call calbec(1,nsp,eigr,cm,becm)
              call gram(betae,becm,nhsa,cm,ngw,nbsp)
           ENDIF
        ENDIF

        !test on energy: check the energy has really diminished

        !call calbec(1,nsp,eigr,cm,becm)
        if(.not.tens) then
          !
          if(non_ortho) then
             call rhoofr(nfi,cm(:,:),cmdual,irb,eigrb,becm,becmdual,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
          else
             call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
          endif
          !
        else
          if(newscheme)  then
              call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,.false., vpot  )
          endif
          !     calculation of the rotated quantities
          call rotate_twin( z0t, cm(:,:), becm, c0diag, becdiag, .false. )
          !     calculation of rho corresponding to the rotated wavefunctions
          call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
        endif

        !calculates the potential
        !
        !     put core charge (if present) in rhoc(r)
        !
        if (nlcc_any) call set_cc(irb,eigrb,rhoc)
        !
        vpot = rhor
!$$
!        if(ionode) write(*,*) 'Now doing vofrho3'
        CALL start_clock( 'vofrho3' )
!$$
        !
        call vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast,             &
                       &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
!$$
        CALL stop_clock( 'vofrho3' )
!$$

        if( tefield )  then!to be bettered
          call berry_energy( enb, enbi, becm, cm(:,:), fion )
          etot=etot+enb+enbi
        endif
        if( tefield2 )  then!to be bettered
          call berry_energy2( enb, enbi, becm, cm(:,:), fion )
          etot=etot+enb+enbi
        endif

!$$
        if(do_orbdep) then
            !warning:giovanni... don't we need becm down here?? otherwise problem with ultrasoft!!
            IF(non_ortho) THEN
               call nksic_potential_non_ortho( nbsp, nbspx, cm, cmdual, fsic, becm, becmdual, rhovan, deeq_sic, &
                                  ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx,&
                                  wfc_centers, wfc_spreads, &
                                  icompute_spread, .false.)
            ELSE
               call nksic_potential( nbsp, nbspx, cm, fsic, bec, rhovan, deeq_sic, &
                                  ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx,&
                                  wfc_centers, wfc_spreads, &
                                  icompute_spread, .false.)
            ENDIF
            eodd = sum(pink(1:nbsp))
            etot = etot + eodd
            !
        endif
!$$ 
        if( do_hf ) then
            !
            call hf_potential( nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
                               nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
                               rhor, rhog, vxxpsi, exx)
            !
            etot = etot + sum(exx(1:nbsp))
            !
        endif

        enever=etot
        if( tens .and. newscheme) enever=enever+entropy
        !
        if( tens .and. newscheme .and. ionode ) then
#ifdef __DEBUG
            write(37,'(a3,4f20.10)')   'CG1',ene0,ene1,enesti,enever
            write(37,'(a3,4f10.7,/)')  'CG2',spasso,passov,passo,(enever-ene0)/passo/dene0
#endif
            !write(stdout,'(a3,4f20.10)')   'CG1',ene0,ene1,enesti,enever
            !write(stdout,'(a3,4f10.7,/)')  'CG2',spasso,passov,passo,(enever-ene0)/passo/dene0
        else
            !
#ifdef __DEBUG
            if(ionode) then
                write(37,'(a3,4f20.10)') 'CG1',ene0+entropy,ene1+entropy,enesti+entropy,enever+entropy
                write(37,'(a3,3f12.7,e20.10,f12.7)')  'CG2',spasso,passov,passo,dene0,(enever-ene0)/passo/dene0
                write(37,"()")
                write(1037,'(a3,4f20.10)') 'CG1',ene0+entropy,ene1+entropy,enesti+entropy,enever+entropy
                write(1037,'(a3,3f12.7,e20.10,f12.7)')  'CG2',spasso,passov,passo,dene0,(enever-ene0)/passo/dene0
                write(1037, "()")
            endif
#endif
            !write(stdout,'(a3,4f20.10)') 'CG1',ene0+entropy,ene1+entropy,enesti+entropy,enever+entropy
            !write(stdout,'(a3,3f12.7,e20.10,f12.7)')  'CG2',spasso,passov,passo,dene0,(enever-ene0)/passo/dene0
            !write(stdout, "()")
            !
        endif

        !
        !check with  what supposed
        !
        if(ionode .and. iprsta > 1 ) then
            write(stdout,"(2x,a,f20.12)") 'cg_sub: estimate :'  , (enesti-enever)/(ene0-enever)
            write(stdout,"(2x,a,3f20.12)") 'cg_sub: minmum   :'  , enever,passo,passov
        endif

        !
        !if the energy has diminished with respect to  ene0 and ene1 , everything ok
        !
        if( ((enever.lt.ene0) .and. (enever.lt.ene1)).or.(tefield.or.tefield2)) then
          c0(:,:)=cm(:,:)
          call copy_twin(bec,becm) !modified:giovanni
          ene_ok=.true.
          if(non_ortho) then
             cdual(:,:)=cmdual(:,:)
             call copy_twin(becdual,becmdual)
          endif
        elseif( (enever.ge.ene1) .and. (enever.lt.ene0)) then
          if(ionode) then
             write(stdout,"(2x,a,i5,f20.12)") 'cg_sub: missed minimum, case 1, iteration',itercg, passof
          endif
          c0(1:ngw,1:nbsp)=c0(1:ngw,1:nbsp)+spasso*passov*hi(1:ngw,1:nbsp)
!$$
          passof=2.d0*passov
!$$
          restartcg=.true.
          !
          IF(non_ortho) THEN
             call calbec(1,nsp,eigr,c0,bec)
             call compute_duals(c0,cdual,nbspx,1)
             call calbec(1,nsp,eigr,cdual,becdual)
          ELSE
             IF(do_orbdep.and.ortho_switch) THEN
                call lowdin(c0, lgam)
                call calbec(1,nsp,eigr,c0,bec)
             ELSE
                call calbec(1,nsp,eigr,c0,bec)
                call gram(betae,bec,nhsa,c0,ngw,nbsp)
             ENDIF
          ENDIF
          !
          ene_ok=.false.
          !if  ene1 << energy <  ene0; go to  ene1
        else if( (enever.ge.ene0).and.(ene0.gt.ene1)) then
          if(ionode) then
             write(stdout,"(2x,a,i5)") 'cg_sub: missed minimum, case 2, iteration',itercg
          endif  
          c0(1:ngw,1:nbsp)=c0(1:ngw,1:nbsp)+spasso*passov*hi(1:ngw,1:nbsp)
!$$
          passof=1.d0*passov
!$$
          restartcg=.true.!ATTENZIONE
          !
          IF(non_ortho) THEN
             call calbec(1,nsp,eigr,c0,bec)
             call compute_duals(c0,cdual,nbspx,1)
             call calbec(1,nsp,eigr,cdual,becdual)
          ELSE
             IF(do_orbdep.and.ortho_switch) THEN
                call lowdin(c0, lgam)
                call calbec(1,nsp,eigr,c0,bec)
             ELSE
                call calbec(1,nsp,eigr,c0,bec)
                call gram(betae,bec,nhsa,c0,ngw,nbsp)
             ENDIF
          ENDIF
          !
          !if ene > ene0,en1 do a steepest descent step
          ene_ok=.false.
        else if((enever.ge.ene0).and.(ene0.le.ene1)) then
          if(ionode) then
             write(stdout,"(2x,a,i5)") 'cg_sub: missed minimum, case 3, iteration',itercg
          endif

          iter3=0
          do while(enever.ge.ene0 .and. iter3.lt.maxiter3)
            iter3=iter3+1

            passov=passov*0.5d0
            cm(1:ngw,1:nbsp)=c0(1:ngw,1:nbsp)+spasso*passov*hi(1:ngw,1:nbsp)
!$$
            passof=1.d0*passov
            itercgeff=itercgeff+1
!$$
            ! chenge the searching direction
            spasso=spasso*(-1.d0)

            IF(non_ortho) THEN 
                 call calbec(1,nsp,eigr,cm,becm)
                 call compute_duals(cm,cmdual,nbspx,1)
                 call calbec(1,nsp,eigr,cmdual,becmdual)
            ELSE
               IF(do_orbdep.and.ortho_switch) THEN
                  call lowdin(cm, lgam)
                  call calbec(1,nsp,eigr,cm,becm)
               ELSE
                  call calbec(1,nsp,eigr,cm,becm)
                  call gram(betae,bec,nhsa,cm,ngw,nbsp)
               ENDIF
            ENDIF

            if(.not.tens) then
              if(non_ortho) then
                 call rhoofr(nfi,cm(:,:),cmdual,irb,eigrb,becm,becmdual,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
              else
                 call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
              endif
            else
              if(newscheme)  then
                  call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,.false., vpot  )
              endif
              !     calculation of the rotated quantities
              call rotate_twin( z0t, cm(:,:), becm, c0diag, becdiag, .false. )
              !     calculation of rho corresponding to the rotated wavefunctions
              call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
            endif
  
            !calculates the potential
            !
            !     put core charge (if present) in rhoc(r)
            !
            if (nlcc_any) call set_cc(irb,eigrb,rhoc)
            !
            vpot = rhor
            !
!$$
            CALL start_clock( 'vofrho4' )
!$$
            call vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
                         ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion)
!$$
            CALL stop_clock( 'vofrho4' )
!$$

            if( tefield)  then !to be bettered
                !
                call berry_energy( enb, enbi, becm, cm(:,:), fion )
                etot=etot+enb+enbi
                !
            endif
            if( tefield2)  then !to be bettered
                !
                call berry_energy2( enb, enbi, becm, cm(:,:), fion )
                etot=etot+enb+enbi
                !
            endif

!$$
            if(do_orbdep) then
                !warning:giovanni don't we need becm down here??? otherwise problems with ultrasoft
                IF(non_ortho) THEN
                   call nksic_potential_non_ortho( nbsp, nbspx, cm,cmdual, fsic, becm,becmdual, rhovan, deeq_sic, &
                                      ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, &
                                      wfc_centers, wfc_spreads, &
                                      icompute_spread, .false.)
                ELSE
                   call nksic_potential( nbsp, nbspx, cm, fsic, bec, rhovan, deeq_sic, &
                                      ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, &
                                      wfc_centers, wfc_spreads, &
                                      icompute_spread, .false.)
                ENDIF
                !
                eodd = sum(pink(1:nbsp))
                etot = etot + eodd
                !
            endif
!$$
            if( do_hf ) then
                !
                call hf_potential( nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
                                   nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
                                   rhor, rhog, vxxpsi, exx)
                !
                etot = etot + sum(exx(1:nbsp))
                !
            endif
            !
            enever=etot
            !
            if( tens .and. newscheme) enever=enever+entropy
            !
          enddo
!$$
          if (ionode) write(stdout,"(2x,a,i5)") 'iter3 = ',iter3
!$$

!$$
          !if(.not.do_orbdep) then
              if(iter3 == maxiter3 .and. enever.gt.ene0) then
                write(stdout,"(2x,a)") 'missed minimum: iter3 = maxiter3'
                write(stdout,*) enever, ene0
!                 if(non_ortho) then
!                    call compute_duals(c0,cdual,nbspx,1)
!                    call calbec(1,nsp,eigr,cdual,becdual)
!                    write(6,*) "checkdual", cdual(1:2,1)
!                 endif
              else if(enever.le.ene0) then
                c0(:,:)=cm(:,:)
                call copy_twin(bec,becm)
                if(non_ortho) then
                   cdual(:,:)=cmdual(:,:)
                   call copy_twin(becdual,becmdual)
                 endif
              endif

          !endif
!$$

          restartcg=.true.
          ene_ok=.false.

!$$
          if(iter3 == maxiter3) then
            passof=passop
          endif
!$$
        end if
        
        if(tens.and.newscheme) enever=enever-entropy
 
        if(.not. ene_ok) call calbec (1,nsp,eigr,c0,bec)

        !calculates phi for pc_daga
        IF(non_ortho) THEN
           CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, nbsp, lgam )
        ELSE
           CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, nbsp, lgam )
        ENDIF
  
        !=======================================================================
        !
        !                 start of the inner loop
        !                 (Uij degrees of freedom)
        !
        !=======================================================================
        !
        if(tens.and. .not.newscheme) then
            !
            call start_clock( "inner_loop" )
            !
            call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                                   rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac, &
                                   c0, bec, firstiter, vpot  )
            ! the following sets up the new energy
            enever=etot
            !
            call stop_clock( "inner_loop" )
            !
        endif
        ! 
        !=======================================================================
        !                 end of the inner loop
        !=======================================================================
        !
!        if ( ( mod( itercg, isave ) == 0 ) ) then
!            !
!            CALL writefile( h, hold ,nfi, c0, c0old, taus, tausm,  &
!                            vels, velsm, acc, lambda, lambdam, xnhe0, xnhem,     &
!                            vnhe, xnhp0, xnhpm, vnhp, nhpcl,nhpdim,ekincm, xnhh0,&
!                            xnhhm, vnhh, velh, fion, tps, z0t, f, rhor )
!            !
!        endif
        !
        !=======================================================================
        !                 end write to file
        !=======================================================================
  
        itercg=itercg+1

!$$
        itercgeff=itercgeff+1
!$$
        !
        call stop_clock( "outer_loop" )

      enddo OUTER_LOOP

#ifdef __DEBUG
        ! for debug and tuning purposes
        if ( ionode ) write(37,*)itercg, itercgeff, etotnew
        if ( ionode ) write(1037,'("iteration =",I4,"  eff iteration =",I4,"   Etot (Ha) =",F22.14)')&
            itercg, itercgeff, etotnew
#endif
      ! 
      !=======================================================================
      !                 end of the main loop
      !=======================================================================

      !
      !calculates atomic forces and lambda
      !

      !
      ! if pressure is need the following is written because of caldbec
      !
      if(tpre) then
         !
         call  calbec(1,nsp,eigr,c0,bec)
         !
         if(.not.tens) then
             call  caldbec( ngw, nhsa, nbsp, 1, nsp, eigr, c0, dbec )
             if(non_ortho) then
                call compute_duals(c0,cdual,nbspx,1)
                call calbec(1,nsp,eigr,cdual,becdual)
                call rhoofr(nfi,c0(:,:),cdual,irb,eigrb,bec,becdual,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
             else
                call rhoofr(nfi,c0(:,:),irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
             endif
         else
             !
             !     calculation of the rotated quantities
             call rotate_twin( z0t, c0(:,:), bec, c0diag, becdiag, .false. )
             !
             !     calculation of rho corresponding to the rotated wavefunctions
             call caldbec( ngw, nhsa, nbsp, 1, nsp, eigr, c0diag, dbec )
             call rhoofr( nfi, c0diag, irb, eigrb, becdiag,      &
                          rhovan, rhor, rhog, rhos, enl, denl, ekin, dekin6)
             !
         endif

         !calculates the potential
         !
         !     put core charge (if present) in rhoc(r)
         !
         if (nlcc_any) call set_cc(irb,eigrb,rhoc)

         !
         !---ensemble-DFT
         !
         vpot = rhor
!$$
         CALL start_clock( 'vofrho5' )
!$$
         call vofrho(nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
                     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion)
!$$
         CALL stop_clock( 'vofrho5' )
!$$

!$$
!$$ Why there are not other terms here???
!$$

!$$
         if(do_orbdep) then
             !
             IF(non_ortho) THEN
                call nksic_potential_non_ortho( nbsp, nbspx, c0, cdual, fsic, bec, becdual, rhovan, &
                                   deeq_sic, &
                                   ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, &
                                   wfc_centers, wfc_spreads, &
                                   icompute_spread, .false.)
             ELSE
                call nksic_potential( nbsp, nbspx, c0, fsic, bec, rhovan, deeq_sic, &
                                   ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, &
                                   wfc_centers, wfc_spreads, &
                                   icompute_spread, .false.)
             ENDIF
             eodd = sum(pink(1:nbsp))
             etot = etot + eodd
             !
         endif
!$$
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
     endif


     if(tens) call calcmt_twin( f, z0t, fmat0, .false. )

     call newd(vpot,irb,eigrb,rhovan,fion)
     if (.not.tens) then
        if (tfor .or. tprnfor) call nlfq(c0,eigr,bec,becdr,fion, lgam)
     else
        if (tfor .or. tprnfor) call nlfq(c0diag,eigr,becdiag,becdrdiag,fion)
     endif
  
     call prefor(eigr,betae)
!$$
     ! faux takes into account spin multiplicity.
     !
     faux(1:nbsp) = max(f_cutoff,f(1:nbsp)) * DBLE( nspin ) / 2.0d0
     !
!$$
     IF(non_ortho) THEN
        hpsi0(:,:)=CMPLX(0.d0,0.d0)
     ENDIF
     !
     IF(do_bare_eigs) THEN
        !
        allocate(c2_bare(ngw), c3_bare(ngw))
        allocate(gi_bare(ngw,nbsp))
        c2_bare=0.d0
        c3_bare=0.d0
        gi_bare=0.d0
        !
     ENDIF
     !
     do i=1,nbsp,2
!$$
         CALL start_clock( 'dforce2' )
!$$          call dforce( i, bec, betae, c0,c2,c3,rhos, nnrsx, ispin,f,n,nspin)
         IF(non_ortho) THEN
            call dforce(i,bec,betae,c0,c2,c3,rhos,nnrsx,ispin,faux,nbsp,nspin)
         ELSE
            call dforce(i,bec,betae,c0,c2,c3,rhos,nnrsx,ispin,faux,nbsp,nspin)
         ENDIF
         !
         CALL start_clock( 'dforce2' )
!$$
         if(tefield.and.(evalue .ne. 0.d0)) then

            call dforceb &
               (c0, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            do ig=1,ngw
              c2(ig)=c2(ig)+evalue*df(ig)
            enddo
            call dforceb &
               (c0, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            do ig=1,ngw
              c3(ig)=c3(ig)+evalue*df(ig)
            enddo
            !
         endif

         if(tefield2.and.(evalue2 .ne. 0.d0)) then

            call dforceb &
               (c0, i, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
            do ig=1,ngw
              c2(ig)=c2(ig)+evalue2*df(ig)
            enddo
            call dforceb &
               (c0, i+1, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
            do ig=1,ngw
              c3(ig)=c3(ig)+evalue2*df(ig)
            enddo

         endif

         IF(do_bare_eigs) THEN
            !
            c2_bare(:) = c2(:)
            c3_bare(:) = c3(:)
            !
         ENDIF
!$$
         if ( do_orbdep ) then
             !
             ! faux takes into account spin multiplicity.
             !
             IF(non_ortho) THEN
                CALL nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, bec, ngw, c0(:,i), c0(:,i+1), vsicpsi, lgam )
             ELSE
                CALL nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, bec, ngw, c0(:,i), c0(:,i+1), vsicpsi, lgam )
             ENDIF
             !
             IF(non_ortho) THEN
                !
                c2(:) = c2(:) - 0.5d0 * vsicpsi(:,1) * faux(i)
                !
                if( i+1 <= nbsp )   c3(:) = c3(:) - 0.5d0 * vsicpsi(:,2) * faux(i+1)
                !
                CALL nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, becdual, ngw, cdual(:,i), cdual(:,i+1), vsicpsi, lgam )
                hpsi0(:,i) = hpsi0(:,i) - 0.5d0 * vsicpsi(:,1) * faux(i)
                !
                IF(i+1<=nbsp) THEN
                   hpsi0(:,i+1) = hpsi0(:,i+1) - 0.5d0 * vsicpsi(:,2)*faux(i+1)
                ENDIF
                !
             ELSE
                !
                c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
                !
                if( i+1 <= nbsp )   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
                !
             ENDIF
             !
         endif
!$$
          IF ( lda_plus_u ) THEN
               !
               IF ( tens .OR. tsmear ) THEN
                   !
                   c2(:) = c2(:) - vupsi(:,i) 
                   if( i+1 <= nbsp ) c3(:) = c3(:) - vupsi(:,i+1) 
                   !
               ELSE
                   !
                   c2(:) = c2(:) - vupsi(:,i) * faux(i)
                   if( i+1 <= nbsp ) c3(:) = c3(:) - vupsi(:,i+1) * faux(i+1)
                   !
               ENDIF
               !
           ENDIF

         if ( do_hf ) then
             !
             c2(:) = c2(:) - vxxpsi(:,i) * faux(i)
             !
             if( i+1 <= nbsp )   c3(:) = c3(:) - vxxpsi(:,i+1) * faux(i+1)
             !
         endif
 
         do ig=1,ngw
            gi(ig,  i)=c2(ig)
            if(i+1 <= nbsp) gi(ig,i+1)=c3(ig)
         enddo
         !
         if (lgam.and.ng0.eq.2) then
            gi(1,  i)=CMPLX(DBLE(gi(1,  i)),0.d0)
            if(i+1 <= nbsp) gi(1,i+1)=CMPLX(DBLE(gi(1,i+1)),0.d0)
         endif
         
         IF(do_bare_eigs) THEN
            !
            do ig=1,ngw
               gi_bare(ig,  i)=c2_bare(ig)
               if(i+1 <= nbsp) gi_bare(ig,i+1)=c3_bare(ig)
            enddo
            !
            if (lgam.and.ng0.eq.2) then
               gi_bare(1,  i)=CMPLX(DBLE(gi_bare(1,  i)),0.d0)
               if(i+1 <= nbsp) gi_bare(1,i+1)=CMPLX(DBLE(gi_bare(1,i+1)),0.d0)
            endif
            !
         ENDIF

     enddo
     
     IF(non_ortho) THEN
        call times_overlap(c0, hpsi0, hpsi, nbsp, 1)
        gi(:,:) = gi(:,:)+hpsi(:,:)
!         call pc3nc_non_ortho(c0,cdual,gi, lgam)
     ENDIF

     IF(.not.lambda(1)%iscmplx) THEN
        allocate(lambda_repl(nudx,nudx))
     ELSE
        allocate(lambda_repl_c(nudx,nudx))
     ENDIF
     !
     IF(non_ortho) THEN
        hitmp(1:ngw,1:nbsp) = cdual(1:ngw,1:nbsp)
     ELSE
        hitmp(1:ngw,1:nbsp) = c0(1:ngw,1:nbsp)
     ENDIF
     !
     do is = 1, nspin
        !
        nss = nupdwn(is)
        istart = iupdwn(is)
           
        IF(.not.lambda(1)%iscmplx) THEN
           lambda_repl = 0.d0
        ELSE
           lambda_repl_c = CMPLX(0.d0,0.d0)
        ENDIF
        !
        IF(non_ortho) THEN
           !
           do i = 1, nss
              do j = 1, nss
                 ii = i + istart - 1
                 jj = j + istart - 1
                 IF(.not.lambda(1)%iscmplx) THEN
                    do ig = 1, ngw
                       lambda_repl( i, j ) = lambda_repl( i, j ) - &
                       2.d0 * DBLE( CONJG( hitmp( ig, ii ) ) * gi( ig, jj) )
                    enddo
                    if( ng0 == 2 ) then
                       lambda_repl( i, j ) = lambda_repl( i, j ) + &
                       DBLE( CONJG( hitmp( 1, ii ) ) * gi( 1, jj ) )
                    endif
                 ELSE
                    do ig = 1, ngw
                       lambda_repl_c( i, j ) = lambda_repl_c( i, j ) - &
                       CONJG( hitmp( ig, ii ) ) * gi( ig, jj)
                    enddo
                 ENDIF
              enddo
           enddo
           !
        ELSE
           !
           do i = 1, nss
              do j = i, nss
                 ii = i + istart - 1
                 jj = j + istart - 1
                 IF(.not.lambda(1)%iscmplx) THEN
                    do ig = 1, ngw
                       lambda_repl( i, j ) = lambda_repl( i, j ) - &
                       2.d0 * DBLE( CONJG( hitmp( ig, ii ) ) * gi( ig, jj) )
                    enddo
                    if( ng0 == 2 ) then
                       lambda_repl( i, j ) = lambda_repl( i, j ) + &
                       DBLE( CONJG( hitmp( 1, ii ) ) * gi( 1, jj ) )
                    endif
                    lambda_repl( j, i ) = lambda_repl( i, j )
                 ELSE
                    do ig = 1, ngw
                       lambda_repl_c( i, j ) = lambda_repl_c( i, j ) - &
                       CONJG( hitmp( ig, ii ) ) * gi( ig, jj)
                    enddo
                    lambda_repl_c( j, i ) = CONJG(lambda_repl_c( i, j ))
                 ENDIF
              enddo
           enddo
           !
         ENDIF
        !
        IF(.not.lambda(1)%iscmplx) THEN
           CALL mp_sum( lambda_repl, intra_image_comm )
           CALL distribute_lambda( lambda_repl, lambda(is)%rvec( :, :), descla( :, is ) )
        ELSE
           CALL mp_sum( lambda_repl_c, intra_image_comm )
           CALL distribute_lambda( lambda_repl_c, lambda(is)%cvec( :, :), descla( :, is ) )
        ENDIF
        !
        !
     end do

     IF(do_bare_eigs) call compute_lambda_bare()

     IF(.not.lambda(1)%iscmplx) THEN
        DEALLOCATE( lambda_repl )
     ELSE
        DEALLOCATE( lambda_repl_c )
     ENDIF
  
     if ( tens ) then
        !
        ! in the ensemble case matrix labda must be multiplied with f
        IF(.not.lambda(1)%iscmplx) THEN
           ALLOCATE( lambda_dist( nlam, nlam ) )
        ELSE
           ALLOCATE( lambda_dist_c( nlam, nlam ) )
        ENDIF

        do iss = 1, nspin
           !
           nss    = nupdwn( iss )
           !
           call set_twin(lambdap(iss), CMPLX(0.d0,0.d0)) !modified:giovanni
           !
           IF(.not.lambdap(iss)%iscmplx) THEN
              CALL cyc2blk_redist( nss, fmat0(iss)%rvec(1,1), nrlx, SIZE(fmat0(iss)%rvec,2), lambda_dist, nlam, nlam, descla(1,iss) )
           ELSE
              CALL cyc2blk_zredist( nss, fmat0(iss)%cvec(1,1), nrlx, SIZE(fmat0(iss)%cvec,2), lambda_dist_c, nlam, nlam, descla(1,iss) )
           ENDIF
           !
           ! Perform lambdap = lambda * fmat0
           !
           IF(.not. lambdap(iss)%iscmplx) then !modified:giovanni
              CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, lambda(iss)%rvec(1,1), nlam, lambda_dist, nlam, &
                                  0.0d0, lambdap(iss)%rvec(1,1), nlam, descla(1,iss) )
           ELSE
              CALL sqr_zmm_cannon( 'N', 'N', nss, (1.0d0,0.d0), lambda(iss)%cvec(1,1), nlam, lambda_dist_c, nlam, &
                                  (0.0d0,0.d0), lambdap(iss)%cvec(1,1), nlam, descla(1,iss) ) !warning:giovanni C or N?
           ENDIF
           !
           !begin_modified:giovanni
           IF(.not.lambdap(iss)%iscmplx) THEN
             lambda_dist(:,:) = lambda(iss)%rvec(:,:)
           ELSE
             lambda_dist_c(:,:) = lambda(iss)%cvec(:,:)
           ENDIF

           call copy_twin(lambda(iss), lambdap(iss))

           IF(.not.lambdap(iss)%iscmplx) THEN
              lambdap(iss)%rvec(:,:) = lambda_dist(:,:)
           ELSE
              lambdap(iss)%cvec(:,:) = lambda_dist_c(:,:)
           ENDIF
           !end_modified:giovanni
           !
        end do
        !
        IF(.not.lambdap(iss)%iscmplx) THEN
           DEALLOCATE( lambda_dist )
        ELSE
          DEALLOCATE( lambda_dist_c )
        ENDIF
        !
        call nlsm2(ngw,nhsa,nbsp,nspin,eigr,c0(:,:),becdr, lgam)
        !
     endif
     !
     call nlfl_twin(bec,becdr,lambda,fion, lgam)
     ! bforceion adds the force term due to electronic berry phase
     ! only in US-case
          
     if( tefield.and.(evalue .ne. 0.d0) ) then
        call bforceion(fion,tfor.or.tprnfor,ipolp, qmat,bec,becdr,gqq,evalue)
     endif
     !
     if( tefield2.and.(evalue2 .ne. 0.d0) ) then
        call bforceion(fion,tfor.or.tprnfor,ipolp2, qmat2,bec,becdr,gqq2,evalue2)
     endif
     deallocate(hpsi0,hpsi,gi,hi)
     deallocate(hitmp, STAT=ierr)
     !        
     call deallocate_twin(s_minus1)
     call deallocate_twin(k_minus1)
#ifdef __DEBUG
     !
     !for debug and tuning purposes
     !
     if(ionode) close(37)
     if(ionode) close(1037)
#endif
     call stop_clock('runcg_uspp')

!         deallocate(bec0,becm,becdrdiag)

     !begin_modified:giovanni
     call deallocate_twin(bec0)
     call deallocate_twin(becm)
     call deallocate_twin(becdrdiag)
     !
!         do i=1,nspin
! 	  call deallocate_twin(lambda(i))
! 	  call deallocate_twin(lambdap(i))
!         enddo
     !
     !end_modified:giovanni

     deallocate(ave_ene)
     deallocate(c2,c3)
     IF(allocated(c2_bare)) deallocate(c2_bare)
     IF(allocated(c3_bare)) deallocate(c3_bare)
     IF(allocated(gi_bare)) deallocate(gi_bare)
     !
     IF(ionode) THEN
        ! 
        CALL printout_base_close( "sha" )
        CALL printout_base_close( "ovp" )
        !
     ENDIF
     !
     return

     contains

     subroutine compute_lambda_bare()
     
        IF(non_ortho) THEN
           hitmp(:,:) = cdual(:,:)
        ELSE
           hitmp(:,:) = c0(:,:)
        ENDIF

        do is = 1, nspin
           !
           nss = nupdwn(is)
           istart = iupdwn(is)

           IF(.not.lambda(1)%iscmplx) THEN
              lambda_repl = 0.d0
           ELSE
              lambda_repl_c = CMPLX(0.d0,0.d0)
           ENDIF
           !
           IF(non_ortho) THEN
              !
              do i = 1, nss
                 do j = 1, nss
                    ii = i + istart - 1
                    jj = j + istart - 1
                    IF(.not.lambda(1)%iscmplx) THEN
                       do ig = 1, ngw
                          lambda_repl( i, j ) = lambda_repl( i, j ) - &
                          2.d0 * DBLE( CONJG( hitmp( ig, ii ) ) * gi_bare( ig, jj) )
                       enddo
                       if( ng0 == 2 ) then
                          lambda_repl( i, j ) = lambda_repl( i, j ) + &
                          DBLE( CONJG( hitmp( 1, ii ) ) * gi_bare( 1, jj ) )
                       endif
                    ELSE
                       do ig = 1, ngw
                          lambda_repl_c( i, j ) = lambda_repl_c( i, j ) - &
                          CONJG( hitmp( ig, ii ) ) * gi_bare( ig, jj)
                       enddo
                    ENDIF
                 enddo
              enddo
              !
           ELSE
              !
              do i = 1, nss
                 do j = i, nss
                    ii = i + istart - 1
                    jj = j + istart - 1
                    IF(.not.lambda(1)%iscmplx) THEN
                       do ig = 1, ngw
                          lambda_repl( i, j ) = lambda_repl( i, j ) - &
                          2.d0 * DBLE( CONJG( hitmp( ig, ii ) ) * gi_bare( ig, jj) )
                       enddo
                       if( ng0 == 2 ) then
                          lambda_repl( i, j ) = lambda_repl( i, j ) + &
                          DBLE( CONJG( hitmp( 1, ii ) ) * gi_bare( 1, jj ) )
                       endif
                       lambda_repl( j, i ) = lambda_repl( i, j )
                    ELSE
                       do ig = 1, ngw
                          lambda_repl_c( i, j ) = lambda_repl_c( i, j ) - &
                          CONJG( hitmp( ig, ii ) ) * gi_bare( ig, jj)
                       enddo
                       lambda_repl_c( j, i ) = CONJG(lambda_repl_c( i, j ))
                    ENDIF
                 enddo
              enddo
              !
           ENDIF
           !
           IF(.not.lambda_bare(1)%iscmplx) THEN
              CALL mp_sum( lambda_repl, intra_image_comm )
              CALL distribute_lambda( lambda_repl, lambda_bare(is)%rvec( :, :), descla( :, is ) )
           ELSE
              CALL mp_sum( lambda_repl_c, intra_image_comm )
              CALL distribute_lambda( lambda_repl_c, lambda_bare(is)%cvec( :, :), descla( :, is ) )
           ENDIF
           !
           !
        end do

     return
     !
     end subroutine compute_lambda_bare

     END SUBROUTINE runcg_uspp
!
!=======================================================================
   subroutine runcg_uspp_emp( c0_emp, cm_emp, bec_emp, f_emp, fsic_emp, n_empx,&
                          n_emps, ispin_emp, iupdwn_emp, nupdwn_emp, phi_emp, lambda_emp, &
                          maxiter_emp, wxd_emp, vsic_emp, sizvsic_emp, pink_emp, nnrx, rhovan_emp, &
                          deeq_sic_emp, nudx_emp, eodd_emp, etot_emp, &
                          filledstates_potential, &
                          nfi, tfirst, tlast, eigr, bec, irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ema0bg, desc_emp)
!=======================================================================

      use kinds,                    only : dp
      use control_flags,            only : iprint, thdyn, tpre, iprsta, &
                                           tfor, taurdr, tprnfor, gamma_only, do_wf_cmplx, tstress !added:giovanni gamma_only, do_wf_cmplx
      use control_flags,            only : ndr, ndw, nbeg, nomore, tsde, tortho, tnosee, &
                                           tnosep, trane, tranp, tsdp, tcp, tcap, ampre, &
                                           amprp, tnoseh, non_ortho
      use core,                     only : nlcc_any
      !---ensemble-DFT
      use energies,                 only : eht, epseu, exc, etot, eself, enl, ekin,&
                                           atot, entropy, egrand, eodd
      use electrons_base,           only : f, nspin, nel, iupdwn, nupdwn, nudx, nelt, &
                                           nbspx, nbsp, ispin
      use ensemble_dft,             only : tens, tsmear,   ef,  z0t, c0diag,  &
                                           becdiag, fmat0, fmat0_diag, e0,  id_matrix_init
      !---
      use gvecp,                    only : ngm
      use gvecs,                    only : ngs
      use gvecb,                    only : ngb
      use gvecw,                    only : ngw
      use reciprocal_vectors,       only : ng0 => gstart
      use cvan,                     only : nvb, ish
      use ions_base,                only : na, nat, pmass, nax, nsp, rcmax
      use grid_dimensions,          only : nnr => nnrx, nr1, nr2, nr3
      use cell_base,                only : ainv, a1, a2, a3
      use cell_base,                only : omega, alat
      use cell_base,                only : h, hold, deth, wmass, tpiba2
      use smooth_grid_dimensions,   only : nnrsx, nr1s, nr2s, nr3s
      use smallbox_grid_dimensions, only : nnrb => nnrbx, nr1b, nr2b, nr3b
      use local_pseudo,             only : vps, rhops
      use io_global,                ONLY : io_global_start, stdout, ionode, ionode_id
      use mp_global,                ONLY : intra_image_comm, np_ortho, me_ortho, ortho_comm, me_image
!       use dener
      use cdvan
      use constants,                only : pi, au_gpa
      use io_files,                 only : psfile, pseudo_dir
      USE io_files,                 ONLY : outdir, prefix
      use uspp,                     only : nhsa=> nkb, nhsavb=> nkbus, betae => vkb, rhovan => becsum, deeq,qq
      use uspp_param,               only : nh
      use cg_module,                only : ene_ok,  maxiter,niter_cg_restart, &
                                           conv_thr, passop, enever, itercg
      use ions_positions,           only : tau0
      use wavefunctions_module,     only : c0,cm, phi => cp, cdual, cmdual
      use efield_module,            only : tefield, evalue, ctable, qmat, detq, ipolp, &
                                           berry_energy, ctabin, gqq, gqqm, df, pberryel, &
                                           tefield2, evalue2, ctable2, qmat2, detq2, ipolp2, &
                                           berry_energy2, ctabin2, gqq2, gqqm2, pberryel2
      use mp,                       only : mp_sum, mp_bcast
      use cp_electronic_mass,       ONLY : emass_cutoff
      use orthogonalize_base,       ONLY : calphi
      use cp_interfaces,            ONLY : rhoofr, dforce, compute_stress, nlfl, set_x_minus1, xminus1
      USE cp_main_variables,        ONLY : nlax, collect_lambda, distribute_lambda, descla, nrlx, nlam
      USE descriptors,              ONLY : la_npc_ , la_npr_ , la_comm_ , la_me_ , la_nrl_ , ldim_cyclic
      USE mp_global,                ONLY : me_image,my_image_id
      !
      use nksic,                    only : do_orbdep, do_innerloop, do_innerloop_cg, innerloop_cg_nsd, &
                                           innerloop_cg_nreset, innerloop_init_n, innerloop_cg_ratio, &
                                           vsicpsi, vsic, wtot, fsic, fion_sic, deeq_sic, f_cutoff, pink, &
                                           do_wxd, sizwtot
      use hfmod,                    only : do_hf, vxxpsi, exx
      use twin_types !added:giovanni
      use control_flags,            only : non_ortho
      use cp_main_variables,        only : becdual, becmdual, overlap, ioverlap
      use electrons_module,         only : wfc_spreads_emp, wfc_centers_emp, icompute_spread, &
                                           wfc_centers, wfc_spreads
      use ldau,                     only : lda_plus_u, vupsi
      use cp_interfaces,            only : gram_empty, nlsm1
      use uspp_param,               only : nhm
      USE descriptors,          ONLY : descla_siz_
!
      implicit none
!
      CHARACTER(LEN=80) :: uname
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      integer     :: nfi
      logical     :: tfirst , tlast
      integer :: sizvsic_emp
      complex(dp) :: eigr(ngw,nat)
      type(twin_matrix)    :: bec !modified:giovanni
      integer     :: irb(3,nat)
      complex(dp) :: eigrb(ngb,nat)
      real(dp)    :: rhor(nnr,nspin)
      complex(dp) :: rhog(ngm,nspin)
      real(dp)    :: rhos(nnrsx,nspin)
      real(dp)    :: rhoc(nnr)
      real(dp)    :: ema0bg(ngw)
      integer :: n_emps, n_empx, iupdwn_emp(nspin), nupdwn_emp(nspin), maxiter_emp, nnrx, &
                 nudx_emp, ispin_emp(n_empx)
      real(dp) :: f_emp(n_empx), fsic_emp(n_empx), wxd_emp(sizvsic_emp,2), vsic_emp(sizvsic_emp, n_empx), &
                  pink_emp(n_empx), rhovan_emp(nhm*(nhm+1)/2, nat, nspin), &
                  deeq_sic_emp(nhm,nhm,nat,n_empx), eodd_emp, etot_emp, & 
                  filledstates_potential(nnrsx,nspin)
      complex(dp) :: c0_emp(ngw, n_empx), cm_emp(ngw, n_empx), phi_emp(ngw, n_empx)
      type(twin_matrix) :: bec_emp, lambda_emp(nspin)
      integer :: desc_emp(descla_siz_, 2)

      !
      ! local variables
      !
      integer     :: i, j, ig, k, is, iss,ia, iv, jv, il, ii, jj, kk, ip
      integer     :: inl, jnl, niter, istart, nss, nrl, me_rot, np_rot , comm
      real(dp)    :: enb, enbi, x
      real(dp)    :: entmp, sta
      complex(dp) :: gamma_c  !warning_giovanni, is it real anyway?
      complex(dp), allocatable :: c2(:), c3(:)
      complex(dp), allocatable :: hpsi(:,:), hpsi0(:,:), gi(:,:), hi(:,:)
!       real(DP),    allocatable :: s_minus1(:,:)    !factors for inverting US S matrix
!       real(DP),    allocatable :: k_minus1(:,:)    !factors for inverting US preconditioning matrix
      type(twin_matrix) :: s_minus1!(:,:)    !factors for inverting US S matrix
      type(twin_matrix) :: k_minus1!(:,:)    !factors for inverting US preconditioning matrix
      real(DP),    allocatable :: lambda_repl(:,:) ! replicated copy of lambda
      real(DP),    allocatable :: lambda_dist(:,:) ! replicated copy of lambda
      complex(DP),    allocatable :: lambda_repl_c(:,:) ! replicated copy of lambda
      complex(DP),    allocatable :: lambda_dist_c(:,:) ! replicated copy of lambda

      !
      real(dp)    :: sca, dumm(1)
      logical     :: newscheme, firstiter
      integer     :: maxiter3
      !
!       real(kind=DP), allocatable :: bec0(:,:), becm(:,:), becdrdiag(:,:,:)
      type(twin_matrix) :: bec0_emp, becm_emp !modified:giovanni
      real(kind=DP), allocatable :: ave_ene(:)!average kinetic energy for preconditioning
      real(kind=DP), allocatable :: fmat_(:,:)!average kinetic energy for preconditioning
      complex(kind=DP), allocatable :: fmat_c_(:,:)!average kinetic energy for preconditioning
      ! 
      logical     :: pre_state!if .true. does preconditioning state by state
      !
      complex(DP)    :: esse_c,essenew_c !factors in c.g.
      logical     :: ltresh!flag for convergence on energy
      real(DP)    :: passo!step to minimum
      real(DP)    :: etotnew, etotold!energies
      real(DP)    :: spasso!sign of small step
      logical     :: restartcg!if .true. restart again the CG algorithm, performing a SD step
      integer     :: numok!counter on converged iterations
      integer     :: iter3
      real(DP)    :: passof,passov !step to minimum: effective, estimated
      real(DP)    :: ene0,ene1,dene0,enesti !energy terms for linear minimization along hi
      !
      real(DP),    allocatable :: faux(:) ! takes into account spin multiplicity
      real(DP),    allocatable :: hpsinorm(:), hpsinosicnorm(:)
      complex(DP), allocatable :: hpsinosic(:,:)
      complex(DP), allocatable :: hitmp(:,:)
      integer     :: ninner,nbnd1,nbnd2,itercgeff
      complex(DP)    :: Omattot(n_empx,n_empx)
      real(DP)    :: dtmp, temp
      real(dp)    :: tmppasso, ene_save(100), ene_save2(100), ene_lda
      !
      logical :: lgam, switch=.false., ortho_switch=.false.
      complex(DP) :: phase
      integer :: ierr, northo_flavor
      real(DP) :: deltae,sic_coeff1, sic_coeff2 !coefficients which may change according to the flavour of SIC
      real(DP) :: ekin_emp, enl_emp, dekin_emp(6), denl_emp(3,3), epot_emp
      real(DP), allocatable :: rhor_emp(:,:), rhos_emp(:,:), rhoc_emp(:)
      complex(DP), allocatable :: rhog_emp(:,:)
      real(DP), allocatable :: faux_emp(:)
      integer :: ndwwf=0, in, in_emp, issw
      !
      lgam = gamma_only.and..not.do_wf_cmplx
      !
      deltae = 2.d0*conv_thr
      !
      allocate (faux(n_empx), faux_emp(n_empx))
      !
      faux_emp=0.d0
      !
      allocate (c2(ngw),c3(ngw))

      call init_twin(bec0_emp, lgam)
      call allocate_twin(bec0_emp, nhsa, n_emps, lgam)
      call init_twin(becm_emp, lgam)
      call allocate_twin(becm_emp, nhsa, n_emps, lgam)
      allocate(rhor_emp(nnr,nspin), rhos_emp(nnrsx,nspin), rhog_emp(ngm,nspin))
      !
      if(nlcc_any) then
         !
         allocate(rhoc_emp(nnr))
         !
      endif
      !
      ! prototype call to gram_empty: use throughout the coding
      !
      ! CALL gram_empty(.false. , eigr, betae, bec_emp, bec, nhsa, &
      !                       c0_emp( :, iupdwn_emp(iss): ), c0( :, iupdwn(iss): ), ngw, nupdwn_emp(iss), nupdwn(iss), iupdwn_emp(iss), iupdwn(iss))

      call start_clock('runcg_uspp')
      
      firstiter=.true.

      pre_state=.false.!normally is disabled

      maxiter3=12
!$$
      if(do_orbdep) maxiter3=10
!$$
      ninner=0

      !
      ! the following is just a beginning; many things to be done...
      !
      if(do_orbdep) then
          !
              fsic = f
          !
      endif

      if( tfirst .and. ionode ) &
         write(stdout,"(/,a,/)") 'PERFORMING CONJUGATE GRADIENT MINIMIZATION OF EL. STATES'
      
!set tpa preconditioning

      call  emass_precond_tpa( ema0bg, tpiba2, emass_cutoff )
     
      call prefor(eigr,betae) 

      ltresh    = .false.
      itercg    = 1
      etotold   = 1.d8
      restartcg = .true.
      passof = passop
      ene_ok = .false.

!$$
      itercgeff = 1
!$$
      CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
      !write(6,*) "checkbec", bec_emp%cvec
      !orthonormalize c0
      do iss=1,nspin
         !
         in     = iupdwn(iss)
         in_emp = iupdwn_emp(iss)
         !
         issw   = iupdwn( iss )
         !
         if(nupdwn(iss)>0.and. nupdwn_emp(iss)>0) &
         CALL gram_empty(.true. , eigr, betae, bec_emp, bec, nhsa, &
                             c0_emp( :, in_emp:), c0( :, issw: ), ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, issw)
         !
      enddo
      !write(6,*) "checkbec2", bec_emp%cvec, nsp, n_emps
      CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
      !write(6,*) "checkbec3", bec_emp%cvec
      !
!          IF(do_orbdep.and.ortho_switch) THEN
!             !
!             call lowdin(c0, lgam) !WARNING, FOR EMPTY STATES IT DOES NOT WORK
!             CALL nlsm1 ( n_emps, 1, nvb, eigr, c0_emp, bec_emp, 1, lgam )
!             !call calbec(1,nsp,eigr,c0,bec)
!          ELSE
!             CALL nlsm1 ( n_emps, 1, nvb, eigr, c0_emp, bec_emp, 1, lgam )
!             !call calbec(1,nsp,eigr,c0,bec)
!             call gram(betae, bec_emp, nhsa, c0_emp, ngw, n_emps)
!          ENDIF
      !call calbec(1,nsp,eigr,c0,bec)
         CALL calphi( c0_emp, ngw, bec_emp, nhsa, betae, phi_emp, n_emps, lgam, ema0bg )
         !CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, nbsp, lgam)
      !
      ! calculates the factors for S and K inversion in US case
      !
      if ( nvb > 0 ) then
          !
          call init_twin(s_minus1,lgam)
          call allocate_twin(s_minus1, nhsavb, nhsavb, lgam)
          call init_twin(k_minus1,lgam)
          call allocate_twin(k_minus1, nhsavb, nhsavb, lgam)
          call  set_x_minus1_twin(betae,s_minus1,dumm,.false.)
          call  set_x_minus1_twin(betae,k_minus1,ema0bg,.true.)
          !
      else
          !
          call init_twin(s_minus1,lgam)
          call allocate_twin(s_minus1, 1, 1, lgam)
          call init_twin(k_minus1,lgam)
          call allocate_twin(k_minus1, 1, 1, lgam)
          !
      endif  
      !
      ! set index on number of converged iterations
      !
      numok = 0

      allocate( hpsi(ngw,n_emps) )
      allocate( hpsi0(ngw,n_emps) )
      allocate( gi(ngw,n_emps), hi(ngw,n_emps) )
      !
      allocate(hitmp(ngw,n_emps))
      hitmp(:,:) = CMPLX(0.d0,0.d0)
      ! allocate(hpsinosic(ngw,n))
      !
      gi(:,:)=CMPLX(0.d0,0.d0)
      hi(:,:)=CMPLX(0.d0,0.d0)


      !=======================================================================
      !                 begin of the main loop
      !=======================================================================
      !
      OUTER_LOOP: &
      do while ( itercg < maxiter .and. (.not. ltresh) )
        !
        call start_clock( "outer_loop" )

!$$$$
!$$$$        if(itercg.ge.10) do_innerloop=.false.
!$$$$

        ENERGY_CHECK: &
        !
        if(.not. ene_ok ) then

          !call calbec(1,nsp,eigr,c0,bec)
          CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
                !write(6,*) "checkbec3b", bec_emp%cvec
          !
!           write(6,*) "checkbounds", ubound(f_emp)
!           write(6,*) "checkbounds", ubound(ispin_emp), ispin_emp
!           write(6,*) "checkbounds", ubound(iupdwn_emp)
!           write(6,*) "checkbounds", ubound(nupdwn_emp)
!           write(6,*) "checkbounds", ubound(c0_emp)
!           write(6,*) "checkbounds", ubound(bec_emp%cvec)
!           write(6,*) "checkbounds", ubound(rhovan_emp)
!           write(6,*) "checkbounds", ubound(rhor_emp)
!           write(6,*) "checkbounds", ubound(rhog_emp)
!           write(6,*) "checkbounds", ubound(rhos_emp)

          !
          call rhoofr_cp_ortho_new &
          ( n_empx, n_emps, nudx_emp, f_emp, ispin_emp, iupdwn_emp, &
          nupdwn_emp, nspin, nfi, c0_emp, irb, eigrb, bec_emp, &
          rhovan_emp, rhor_emp, rhog_emp, rhos_emp, enl_emp, denl_emp, &
          ekin_emp, dekin_emp, tstress, ndwwf)
          
          call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp) 
          
          etot_emp=epot_emp+enl_emp+ekin_emp
          write(6,*) "etot", etot_emp, "epot", epot_emp, "enl_emp", enl_emp, "ekin", ekin_emp
          
          ! potential energy of filled states times density of empty states
          
          !call rhoofr(nfi,c0(:,:),irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
           
          !
          ! when cycle is restarted go to diagonal representation
          !
          ! CHP: do we need to do the following even if when we do not use ensemble dft?
          !      I have added this additional constraint.
          !
        
          !
          ! calculates the potential
          !
          !     put core charge (if present) in rhoc(r)
          !
          
          !if (nlcc_any) call set_cc(irb,eigrb,rhoc) !warning:I comment this for core charge, for the moment

          !
          !---ensemble-DFT

          !vpot = rhor

!$$
          !CALL start_clock( 'vofrho1' )
!$$
          !call vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast,             &
          !       &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)

          !ene_lda = etot
!$$
          CALL stop_clock( 'vofrho1' )
!$$

!$$
          if( do_orbdep ) then
              !
!                                write(6,*) "check1bounds", ubound(c0_emp)
!                  write(6,*) "checkbounds", ubound(fsic_emp)
!                  write(6,*) "checkbounds", ubound(rhovan_emp)
!                  write(6,*) "checkbounds", ubound(wfc_centers_emp), nudx_emp
!                  write(6,*) "checkbounds", ubound(wfc_spreads_emp)
!                  write(6,*) "checkbounds", ubound(deeq_sic_emp)
!                  write(6,*) "checkbounds", ubound(ispin_emp)
!                  write(6,*) "checkbounds", ubound(iupdwn_emp)
!                  write(6,*) "checkbounds", ubound(nupdwn_emp)
!                  write(6,*) "checkbounds", ubound(rhor), "rhor"
!                  write(6,*) "checkbounds", ubound(rhog), "rhog"
!                  write(6,*) "checkbounds", ubound(wtot), "wtot"
!                  write(6,*) "checkbounds", ubound(vsic_emp), "vsic"
!                  write(6,*) "checkbounds", ubound(pink_emp), "pink", nudx_emp
                 call nksic_potential( n_emps, n_empx, c0_emp, faux_emp, bec_emp, rhovan_emp, deeq_sic_emp, &
                 ispin_emp, iupdwn_emp, nupdwn_emp, rhor, rhoc, wtot, sizwtot, vsic_emp, .false., & 
                 pink_emp, nudx_emp, wfc_centers_emp, wfc_spreads_emp, icompute_spread, .true. )

              eodd_emp=sum(pink_emp(1:n_empx))
!               write(6,*) eodd_emp, etot_emp, "EODD0", pink_emp
              etot_emp = etot_emp + eodd_emp
              !
          endif
!$$
!           if( do_hf ) then
!               !
!               call hf_potential( nbsp, nbspx, c0, f, ispin, iupdwn, nupdwn, &
!                                  nbsp, nbspx, c0, f, ispin, iupdwn, nupdwn, &
!                                  rhor, rhog, vxxpsi, exx)
!               !
!               etot_emp = etot_emp + sum(exx(1:nbsp))
!               !
!           endif

              etotnew=etot_emp

        else

          etot_emp=enever
          etotnew=etot_emp
          ene_ok=.false.

        end if ENERGY_CHECK

!$$
        if( do_orbdep ) then

          if(do_innerloop) then
!$$$$          if(do_innerloop.and.itercg.le.20) then
!$$$$
             !
             !call start_clock( "inner_loop" )
             !
             eodd_emp    = sum(pink_emp(1:n_emps))
             etot_emp    = etot_emp - eodd_emp
             etotnew = etotnew - eodd_emp
             ninner  = 0

             if(.not.do_innerloop_cg) then !warning, this does not work with empty states
                 call nksic_rot_emin(itercg,ninner,etot,Omattot, lgam)
             else
                 call nksic_rot_emin_cg_new(c0_emp, cm_emp, vsic_emp, ngw, nnrx, &
                      bec_emp, itercg, innerloop_init_n, ninner, etot, Omattot, &
                      deltae*innerloop_cg_ratio, n_empx, nudx_emp, nspin, iupdwn_emp, &
                      nupdwn_emp, pink_emp, wfc_centers_emp, wfc_spreads_emp, lgam)
             endif

!$$ Now rotate hi(:,:) according to Omattot!
!$$ It seems that not rotating hi gives us better convergence.
!$$ So, we do not perform the following routine.
!$$
!            if(ninner.ge.2) then
!              hitmp(:,:) = CMPLX(0.d0,0.d0)
!              do nbnd1=1,n
!                do nbnd2=1,n
!                  hitmp(:,nbnd1)=hitmp(:,nbnd1) + hi(:,nbnd2) * Omattot(nbnd2,nbnd1)
!                enddo
!              enddo
!              hi(:,:) = hitmp(:,:)
!            endif
!$$
             eodd_emp    = sum(pink_emp(1:n_emps))
!              write(6,*) eodd, etot, "EODD_inn", etot+eodd
             etot_emp    = etot_emp + eodd_emp
             etotnew = etotnew + eodd_emp

             !call stop_clock( "inner_loop" )
             !
           endif
           !
        endif
!$$

!$$     
        if ( ionode ) write(stdout,'(5x,"iteration =",I4,"  eff iteration =",I4,"   Etot (Ha) =",F22.14," delta_E=",E22.14)')&
            itercg, itercgeff, etotnew, deltae

        if ( ionode .and. mod(itercg,10) == 0 ) write(stdout,"()" )
!$$


!$$ to see the outer loop energy convergence
        if (do_orbdep) then
            !
            eodd_emp = sum(pink_emp(1:n_empx))
            !
        endif
!$$
        deltae=abs(etotnew-etotold) 
        !
        if( deltae < conv_thr ) then
           numok=numok+1
        else 
           numok=0
        endif
        !
        if( numok >= 4 ) ltresh=.true.
        !
        if(ltresh.or.itercg==maxiter-1) icompute_spread=.true.
        !
        etotold=etotnew
        ene0=etot_emp
        !

       !update d

!         call newd(vpot,irb,eigrb,rhovan,fion) !warning:giovanni I comment this

        call prefor(eigr,betae)!ATTENZIONE
!$$
        ! faux takes into account spin multiplicity.
        !
        faux(1:n_empx)=0.d0
        faux(1:n_empx) = f_emp(1:n_empx) * DBLE( nspin ) / 2.0d0
!$$
        !
        do i=1, n_emps, 2
!$$  FIRST CALL TO DFORCE
          CALL start_clock( 'dforce1' )
!$$          call dforce( i, bec, betae, c0,c2,c3,rhos, nnrsx, ispin,f,n,nspin)
             call dforce( i, bec_emp, betae, c0_emp, c2, c3, filledstates_potential, nnrsx, ispin_emp, faux, n_emps, nspin)
          !
          CALL stop_clock( 'dforce1' )
!$$
          
          IF ( lda_plus_u ) THEN
               !
               IF ( tens .OR. tsmear ) THEN
                   !
                   c2(:) = c2(:) - vupsi(:,i) 
                   if( i+1 <= nbsp ) c3(:) = c3(:) - vupsi(:,i+1) 
                   !
               ELSE
                   !
                   c2(:) = c2(:) - vupsi(:,i) * faux(i)
                   if( i+1 <= nbsp ) c3(:) = c3(:) - vupsi(:,i+1) * faux(i+1)
                   !
               ENDIF
               !
           ENDIF

!$$
          if ( do_orbdep ) then
              !
              ! faux takes into account spin multiplicity.
              !
                 CALL nksic_eforce( i, n_emps, n_empx, vsic_emp, deeq_sic_emp, bec_emp, ngw, c0_emp(:,i), c0_emp(:,i+1), vsicpsi, lgam )
                 !
                 c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
                 !
                 if( i+1 <= n_emps )   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
                 !
              !
          endif
!$$
!           if ( do_hf ) then
!               !
!               c2(:) = c2(:) - vxxpsi(:,i) * faux(i)
!               !
!               if( i+1 <= n_emps )   c3(:) = c3(:) - vxxpsi(:,i+1) * faux(i+1)
!               !
!           endif

          hpsi(1:ngw,  i)=c2(1:ngw)
          if(i+1 <= n_emps) then
              hpsi(1:ngw,i+1)=c3(1:ngw)
          endif
          !
          IF(lgam) THEN
	    if (ng0.eq.2) then
		hpsi(1,  i)=CMPLX(DBLE(hpsi(1,  i)), 0.d0)
		if(i+1 <= n_emps) then
		    hpsi(1,i+1)=CMPLX(DBLE(hpsi(1,i+1)), 0.d0)
		endif
	    endif
          ENDIF
          !
        enddo
        
!$$

!         if(pre_state) call ave_kin(c0,SIZE(c0,1),nbsp,ave_ene) ! warning:giovanni I commented this

!$$        call pcdaga2(c0,phi,hpsi)
!$$     HPSI IS ORTHOGONALIZED TO the filled manifold and to c0
        CALL nlsm1 ( n_emps, 1, nsp, eigr, hpsi, becm_emp, 1, lgam )

        do iss=1,nspin
           !
           in     = iupdwn(iss)
           in_emp = iupdwn_emp(iss)
           !
           issw   = iupdwn( iss )
           !
           CALL gram_empty(.true. , eigr, betae, becm_emp, bec, nhsa, &
                             hpsi( :, in_emp: ), c0( :, issw: ), &
                             ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, issw)
           !
        enddo
        
           if(switch.or.(.not.do_orbdep)) then
             call pcdaga2_new(c0_emp,phi_emp,hpsi,n_emps, ispin_emp ,lgam)
           else
!           call calbec(1,nsp,eigr,hpsi,becm)
            call pc3nc_new(c0_emp,hpsi,n_empx, ispin_emp,lgam)
!           call pc3us(c0,bec,hpsi,becm,lgam)
!           call pcdaga3(c0,phi,hpsi, lgam)
           endif
!$$
        !TWO VECTORS INITIALIZED TO HPSI
           hpsi0(1:ngw,1:n_emps) = hpsi(1:ngw,1:n_emps)
           gi(1:ngw,1:n_emps)    = hpsi(1:ngw,1:n_emps)

	!COMPUTES ULTRASOFT-PRECONDITIONED HPSI, non kinetic-preconditioned, is the subsequent reorthogonalization necessary in the norm conserving case???: giovanni
!         call calbec(1,nsp,eigr,hpsi,becm_emp) !warning:giovanni substitute with nlsm1
!         CALL nlsm1 ( n_emps, 1, nvb, eigr, hpsi, becm_emp, 1, lgam )
        !
        call xminus1_twin_new(hpsi,n_emps,betae,dumm,becm_emp,s_minus1,.false.)
!        call sminus1(hpsi,becm,betae)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!look if the following two lines are really needed
!         call calbec(1,nsp,eigr,hpsi,becm_emp) !warning:giovanni substitute with nlsm1
        CALL nlsm1 ( n_emps, 1, nsp, eigr, hpsi, becm_emp, 1, lgam )
        !write(6,*) "becm_emp", becm_emp%cvec
!$$        call pc2(c0,bec,hpsi,becm)
!$$     THIS ORTHOGONALIZES PRECONDITIONED VECTOR HPSI again to filled states and to c0
        do iss=1,nspin
           !
           in     = iupdwn(iss)
           in_emp = iupdwn_emp(iss)
           !
           issw   = iupdwn( iss )
           !
           CALL gram_empty(.true. , eigr, betae, becm_emp, bec, nhsa, &
                             hpsi( :, in_emp: ), c0( :, issw: ), &
                             ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, issw)
           !
        enddo

           if(switch.or.(.not.do_orbdep)) then
             call pc2_new(c0_emp,bec_emp,hpsi,becm_emp, n_emps, &
                      nupdwn_emp, iupdwn_emp, ispin_emp, lgam)
           else
             call pc3nc_new(c0_emp,hpsi,n_emps,ispin_emp, lgam)
!           call pc3us(c0,bec,hpsi,becm, lgam)
           endif
!$$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       COMPUTES ULTRASOFT+KINETIC-preconditioned GI
!        call kminus1(gi,betae,ema0bg)
        IF(do_orbdep) THEN ! preconditioning with respect to spreads.. the gradient along wavefunctions with the largest localization is a bit shortened
           !
           !do i=1,nbsp
           !   gi(:,i) = gi(:,i)*(1.d0+1.d0/sqrt(wfc_spreads(1+i-iupdwn(ispin(i)),ispin(i),2)))
           !enddo
           !
        ENDIF
        !
        CALL nlsm1 ( n_emps, 1, nsp, eigr, gi, bec0_emp, 1, lgam )
        
           if(.not.pre_state) then
              call xminus1_twin_new(gi,n_emps, betae,ema0bg,bec0_emp,k_minus1,.true.)
           else
              call xminus1_state(gi,betae,ema0bg,bec0_emp,k_minus1,.true.,ave_ene) !warning:giovanni not yet implemented
           endif
        !
!         call calbec(1,nsp,eigr,gi,becm_emp) !warning:giovanni substitute with nlsm1
        CALL nlsm1 ( n_emps, 1, nsp, eigr, gi, bec0_emp, 1, lgam )
!$$        call pc2(c0,bec,gi,becm)
!$$     !ORTHOGONALIZES GI to filled states and to c0 ... check if bec0_emp is computed correctly

        do iss=1,nspin
           !
           in     = iupdwn(iss)
           in_emp = iupdwn_emp(iss)
           !
           issw   = iupdwn( iss )
           !
           CALL gram_empty(.true. , eigr, betae, bec0_emp, bec, nhsa, &
                             gi( :, in_emp: ), c0( :, issw:), &
                             ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, issw)
           !
        enddo

           if(switch.or.(.not.do_orbdep)) then
             call pc2_new(c0_emp,bec_emp,gi,bec0_emp,n_emps, &
                       nupdwn_emp, iupdwn_emp, ispin_emp, lgam)
           else
             call pc3nc_new(c0_emp,gi, n_empx,ispin_emp, lgam)
!           call pc3us(c0,bec, gi,becm, lgam)
           endif
!$$
!         call calbec(1,nsp,eigr,hpsi,bec0_emp) !warning:giovanni substitute with nlsm1
!         CALL nlsm1 ( n_emps, 1, nvb, eigr, hpsi, becm_emp, 1, lgam )
!  calculates gamma
        gamma_c=CMPLX(0.d0,0.d0)
        
           do i=1,n_emps
              IF(lgam) THEN
		do ig=1,ngw
		  gamma_c=gamma_c+2.d0*DBLE(CONJG(gi(ig,i))*hpsi(ig,i))
		enddo
		if (ng0.eq.2) then
		  gamma_c=gamma_c-DBLE(CONJG(gi(1,i))*hpsi(1,i))
		endif
              ELSE
		do ig=1,ngw
		  gamma_c=gamma_c+CONJG(gi(ig,i))*hpsi(ig,i)
		enddo
              ENDIF
           enddo
           

	  call mp_sum( gamma_c, intra_image_comm )
           
	   if (nvb.gt.0) then
            if(.not.becm_emp%iscmplx) then
		do i=1,n_emps
		  do is=1,nvb
		      do iv=1,nh(is)
			do jv=1,nh(is)
			    do ia=1,na(is)
			      inl=ish(is)+(iv-1)*na(is)+ia
			      jnl=ish(is)+(jv-1)*na(is)+ia
			      gamma_c=gamma_c+ qq(iv,jv,is)*becm_emp%rvec(inl,i)*bec0_emp%rvec(jnl,i)
			    end do
			end do
		      end do
		  end do
		enddo
            else
		do i=1,n_emps
		  do is=1,nvb
		      do iv=1,nh(is)
			do jv=1,nh(is)
			    do ia=1,na(is)
			      inl=ish(is)+(iv-1)*na(is)+ia
			      jnl=ish(is)+(jv-1)*na(is)+ia
			      gamma_c=gamma_c+ qq(iv,jv,is)*(becm_emp%cvec(inl,i))*CONJG(bec0_emp%cvec(jnl,i)) !warning:giovanni CONJG
			    end do
			end do
		      end do
		  end do
		enddo
	    endif
          endif

        !case of first iteration

! 	IF(lgam) THEN
	  gamma_c=CMPLX(DBLE(gamma_c),0.d0)
! 	  write(6,*) "gamma", gamma_c
! 	ENDIF

!$$        if(itercg==1.or.(mod(itercg,niter_cg_restart).eq.1).or.restartcg) then
        if( itercg==1 .or. mod(itercg,niter_cg_restart)==0 .or. restartcg) then
!$$

          restartcg=.false.
!$$  We do not have to reset passof every exception of CG!
!warning:giovanni if we do not reset we may have fake convergences!!!!
          passof=passop

          hi(1:ngw,1:n_emps)=gi(1:ngw,1:n_emps)!hi is the search direction

          esse_c=gamma_c

        else

          !find direction hi for general case 
          !calculates gamma for general case, not using Polak Ribiere

          essenew_c=gamma_c
          gamma_c=gamma_c/esse_c
          esse_c=essenew_c

          hi(1:ngw,1:n_emps)=gi(1:ngw,1:n_emps)+(gamma_c)*hi(1:ngw,1:n_emps)

        endif
!note that hi, is saved  on gi, because we need it before projection on conduction states

        !find minimum along direction hi:

        !project hi on conduction sub-space

!         call calbec(1,nsp,eigr,hi,bec0_emp)
        CALL nlsm1 ( n_emps, 1, nsp, eigr, hi, bec0_emp, 1, lgam )
        !write(6,*) "bec0", bec0_emp%cvec

!$$        call pc2(c0,bec,hi,bec0)
!$$
           if(switch.or.(.not.do_orbdep)) then
              call pc2_new(c0_emp,bec_emp,hi,bec0_emp,n_emps, &
                      nupdwn_emp, iupdwn_emp, ispin_emp, lgam)
           else
              call pc3nc_new(c0_emp,hi,n_emps,ispin_emp, lgam)
!           call pc3us(c0,bec,hi,bec0, lgam)
           endif
!$$
        !do quadratic minimization
        !             
        !calculate derivative with respect to  lambda along direction hi

        dene0=0.
          do i=1,n_emps
            IF(lgam) THEN              
	      do ig=1,ngw
		dene0=dene0-4.d0*DBLE(CONJG(hi(ig,i))*hpsi0(ig,i))
	      enddo
	      if (ng0.eq.2) then
		dene0=dene0+2.d0*DBLE(CONJG(hi(1,i))*hpsi0(1,i))
	      endif
            ELSE
	      do ig=1,ngw
		dene0=dene0-2.d0*DBLE(CONJG(hi(ig,i))*hpsi0(ig,i))
	      enddo
            ENDIF
          end do
!$$ We need the following because n for spin 2 is double that for spin 1!
          dene0 = dene0 *2.d0/nspin
!$$          dene0 = dene0 *4.d0/nspin
!$$
      call mp_sum( dene0, intra_image_comm )

        !if the derivative is positive, search along opposite direction
      if(dene0.gt.0.d0) then
         spasso=-1.D0
      else
         spasso=1.d0
      endif
      !
      ! calculates wave-functions on a point on direction hi
      !
      cm_emp(1:ngw,1:n_emps)=c0_emp(1:ngw,1:n_emps)+spasso*passof*hi(1:ngw,1:n_emps)
      !
!$$   ! I do not know why the following 3 lines 
      ! were not in the original code (CHP)
      !
      if(lgam.and.ng0 == 2)  cm_emp(1,:)=0.5d0*(cm_emp(1,:)+CONJG(cm_emp(1,:)))
!$$

      !orthonormalize

      !
         if(do_orbdep.and.ortho_switch) then
            !
            call lowdin(cm_emp, lgam)
!             call calbec(1,nsp,eigr,cm_emp,becm_emp) !warning:giovanni substitute with nlsm1
            CALL nlsm1 ( n_emps, 1, nsp, eigr, cm_emp, becm_emp, 1, lgam )
            !
         else
            !
!             call calbec(1,nsp,eigr,cm_emp,becm_emp) !warning:giovanni substitute with nlsm1
            CALL nlsm1 ( n_emps, 1, nsp, eigr, cm_emp, becm_emp, 1, lgam )
            !
            call gram(betae,becm_emp,nhsa,cm_emp,ngw,n_emps)
            !
         endif
            !call calbec(1,nsp,eigr,cm,becm)
            !
        !****calculate energy ene1
       call rhoofr_cp_ortho_new &
       ( n_empx, n_emps, nudx_emp, f_emp, ispin_emp, iupdwn_emp, &
         nupdwn_emp, nspin, nfi, cm_emp, irb, eigrb, becm_emp, &
         rhovan_emp, rhor_emp, rhog_emp, rhos_emp, enl_emp, denl_emp, &
         ekin_emp, dekin_emp, tstress, ndwwf)

       call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp)
       
       etot_emp=epot_emp+enl_emp+ekin_emp
                 write(6,*) "etotm", etot_emp, "epot", epot_emp, "enl_emp", enl_emp, "ekin", ekin_emp

        !calculate potential
        !
        !     put core charge (if present) in rhoc(r)
        !
!         if (nlcc_any) call set_cc(irb,eigrb,rhoc)  !warning:giovanni comment this for empty states
        !
!         vpot = rhor
        !
!$$
!        if(ionode) write(*,*) 'Now doing vofrho2'
        CALL start_clock( 'vofrho2' )
!$$
!         call vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast,             &
!                       &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
!$$
        CALL stop_clock( 'vofrho2' )
!$$

!$$
        if( do_orbdep ) then
            !warning:giovanni don't we need becm down here??? otherwise problems with ultrasoft!!
               call nksic_potential( n_emps, n_empx, cm_emp, faux_emp, becm_emp, &
                                  rhovan_emp, deeq_sic, ispin_emp, iupdwn_emp, &
                                  nupdwn_emp, rhor, rhoc, wtot, sizwtot, vsic_emp, .false., &
                                  pink_emp, nudx_emp, wfc_centers_emp, &
                                  wfc_spreads_emp, icompute_spread, .true.)
            !
            eodd_emp=sum(pink_emp(1:n_emps))
!             write(6,*) eodd_emp, etot_emp, "EODD2", pink_emp !debug:giovanni
            etot_emp = etot_emp + eodd_emp
            !
        endif
!$$
!         if( do_hf ) then
!             !
!             call hf_potential( nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
!                                nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
!                                rhor, rhog, vxxpsi, exx)
!             !
!             etot = etot + sum(exx(1:nbsp))
!             !
!         endif

        ene1=etot_emp
              
        !find the minimum

        call minparabola(ene0,spasso*dene0,ene1,passof,passo,enesti)

        if( ionode .and. iprsta > 1 ) write(stdout,"(6f20.12)") ene0,dene0,ene1,passo, DBLE(gamma_c), esse_c

        !set new step

        passov=passof
!         passof=2.d0*passo
!$$ doing the following makes the convergence better...
        passof=passo
!$$$$
              
        !calculates wave-functions at minimum

        cm_emp(1:ngw,1:n_emps) = c0_emp(1:ngw,1:n_emps) +spasso*passo*hi(1:ngw,1:n_emps)
        !
        if(lgam.and. ng0 == 2 )  THEN
          cm_emp(1,:) = 0.5d0*(cm_emp(1,:)+CONJG(cm_emp(1,:)))
        ELSE !warning:giovanni this would fix the phase of the new position.. should
             !        not influence the calculation
        ENDIF
      
           IF(do_orbdep.and.ortho_switch) THEN
              call lowdin(cm_emp, lgam)
!               call calbec(1,nsp,eigr,cm_emp,becm_emp) !warning:giovanni substitute with nlsm1
              CALL nlsm1 ( n_emps, 1, nsp, eigr, cm_emp, becm_emp, 1, lgam )
           ELSE
!               call calbec(1,nsp,eigr,cm_emp,becm_emp)  !warning:giovanni substitute with nlsm1
              CALL nlsm1 ( n_emps, 1, nsp, eigr, cm_emp, becm_emp, 1, lgam )
              call gram(betae,becm_emp,nhsa,cm_emp,ngw,n_emps)
           ENDIF

        !test on energy: check the energy has really diminished

        !call calbec(1,nsp,eigr,cm,becm)
          !
          call rhoofr_cp_ortho_new &
          ( n_empx, n_emps, nudx_emp, f_emp, ispin_emp, iupdwn_emp, &
          nupdwn_emp, nspin, nfi, cm_emp, irb, eigrb, becm_emp, &
          rhovan_emp, rhor_emp, rhog_emp, rhos_emp, enl_emp, denl_emp, &
          ekin_emp, dekin_emp, tstress, ndwwf)
          
          call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp) 
          
          etot_emp=epot_emp+enl_emp+ekin_emp
          write(6,*) "etotm", etot_emp, "epot", epot_emp, "enl_emp", enl_emp, "ekin", ekin_emp
          
!           call rhoofr(nfi,c0(:,:),irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
          !

        !calculates the potential
        !
        !     put core charge (if present) in rhoc(r)
        !
!         if (nlcc_any) call set_cc(irb,eigrb,rhoc) !warning:giovanni commenting this
        !
!         vpot = rhor
!$$
!        if(ionode) write(*,*) 'Now doing vofrho3'
!         CALL start_clock( 'vofrho3' )
!$$
        !
!         call vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast,             &
!                        &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
!$$
!         CALL stop_clock( 'vofrho3' )
!$$

!$$
        if(do_orbdep) then
            !warning:giovanni... don't we need becm down here?? otherwise problem with ultrasoft!!
            call nksic_potential( n_emps, n_empx, cm_emp, faux_emp, becm_emp, rhovan_emp, deeq_sic_emp, &
            ispin_emp, iupdwn_emp, nupdwn_emp, rhor, rhoc, wtot, sizwtot, vsic_emp, .false., & 
            pink_emp, nudx_emp, wfc_centers_emp, &
            wfc_spreads_emp, icompute_spread, .true. )
            
!                call nksic_potential( nbsp, nbspx, cm, fsic, bec, rhovan, deeq_sic, &
! !                                   ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, do_wxd, pink, nudx,&
!                                   wfc_centers, wfc_spreads, &
!                                   icompute_spread)
            eodd_emp = sum(pink_emp(1:n_emps))
            write(6,*) eodd_emp, etot_emp,"EODD3"
            etot_emp = etot_emp + eodd_emp
            !
        endif
!$$ 
!         if( do_hf ) then
!             !
!             call hf_potential( nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
!                                nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
!                                rhor, rhog, vxxpsi, exx)
!             !
!             etot = etot + sum(exx(1:nbsp))
!             !
!         endif
        enever=etot_emp
        !
        !
        !check with  what supposed
        !
        if(ionode .and. iprsta > 1 ) then
            write(stdout,"(2x,a,f20.12)") 'cg_sub: estimate :'  , (enesti-enever)/(ene0-enever)
            write(stdout,"(2x,a,3f20.12)") 'cg_sub: minmum   :'  , enever,passo,passov
        endif

        !
        !if the energy has diminished with respect to  ene0 and ene1 , everything ok
        !
!         write(6,*) "ENERGYCHECK", ene0, enever, ene1
        if( ((enever.lt.ene0) .and. (enever.lt.ene1)).or.(tefield.or.tefield2)) then
          c0_emp(:,:) = cm_emp(:,:)
          call copy_twin(bec_emp, becm_emp) !modified:giovanni
          ene_ok=.true.
        elseif( (enever.ge.ene1) .and. (enever.lt.ene0)) then
          if(ionode) then
             write(stdout,"(2x,a,i5,f20.12)") 'cg_sub: missed minimum, case 1, iteration',itercg, passof
!              write(6,*) "checkenergies",ene0,enever,ene1
          endif
          c0_emp(1:ngw,1:n_emps)=c0_emp(1:ngw,1:n_emps)+spasso*passov*hi(1:ngw,1:n_emps)
!$$
          passof=2.d0*passov
!$$
          restartcg=.true.
          !
             IF(do_orbdep.and.ortho_switch) THEN
                call lowdin(c0_emp, lgam)
!                 call calbec(1,nsp,eigr,c0_emp,bec_emp)
                CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
             ELSE
!                 call calbec(1,nsp,eigr,c0_emp,bec_emp)
                CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
                call gram(betae,bec_emp,nhsa,c0_emp,ngw,n_emps)
             ENDIF
          !
          ene_ok=.false.
          !if  ene1 << energy <  ene0; go to  ene1
        else if( (enever.ge.ene0).and.(ene0.gt.ene1)) then
          if(ionode) then
             write(stdout,"(2x,a,i5)") 'cg_sub: missed minimum, case 2, iteration',itercg
             write(6,*) "checkenergies",ene0,enever,ene1
          endif  
          c0_emp(1:ngw,1:n_emps)=c0_emp(1:ngw,1:n_emps)+spasso*passov*hi(1:ngw,1:n_emps)
!$$
          passof=1.d0*passov
!$$
          restartcg=.true.!ATTENZIONE
          !
             IF(do_orbdep.and.ortho_switch) THEN
                call lowdin(c0, lgam)
!                 call calbec(1,nsp,eigr,c0_emp,bec_emp)
                CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
             ELSE
!                 call calbec(1,nsp,eigr,c0,bec)
                CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
                call gram(betae,bec_emp,nhsa,c0_emp,ngw,n_emps)
             ENDIF
          !
          !if ene > ene0,en1 do a steepest descent step
          ene_ok=.false.
        else if((enever.ge.ene0).and.(ene0.le.ene1)) then
          if(ionode) then
             write(stdout,"(2x,a,i5)") 'cg_sub: missed minimum, case 3, iteration',itercg
             write(6,*) "checkenergies",ene0,enever,ene1
          endif

          iter3=0
          do while(enever.ge.ene0 .and. iter3.lt.maxiter3)
            iter3=iter3+1

            passov=passov*0.5d0
            cm_emp(1:ngw,1:n_emps)=c0_emp(1:ngw,1:n_emps)+spasso*passov*hi(1:ngw,1:n_emps)
!$$
            passof=1.d0*passov
            itercgeff=itercgeff+1
!$$
            ! chenge the searching direction
            spasso=spasso*(-1.d0)

               IF(do_orbdep.and.ortho_switch) THEN
                  call lowdin(cm_emp, lgam)
!                   call calbec(1, nsp, eigr, cm_emp, becm_emp)
                  CALL nlsm1 ( n_emps, 1, nsp, eigr, cm_emp, becm_emp, 1, lgam )
               ELSE
!                   call calbec(1, nsp, eigr, cm_emp, becm_emp)
                  CALL nlsm1 ( n_emps, 1, nsp, eigr, cm_emp, becm_emp, 1, lgam )
                  call gram(betae, becm_emp, nhsa, cm_emp, ngw, n_emps)
               ENDIF
               
               call rhoofr_cp_ortho_new &
          ( n_empx, n_emps, nudx_emp, f_emp, ispin_emp, iupdwn_emp, &
          nupdwn_emp, nspin, nfi, cm_emp, irb, eigrb, becm_emp, &
          rhovan_emp, rhor_emp, rhog_emp, rhos_emp, enl_emp, denl_emp, &
          ekin_emp, dekin_emp, tstress, ndwwf)
!                call rhoofr(nfi,c0(:,:),irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
          call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp) 
          
          etot_emp=epot_emp+enl_emp+ekin_emp
          write(6,*) "etotm", etot_emp, "epot", epot_emp, "enl_emp", enl_emp, "ekin", ekin_emp
            !calculates the potential
            !
            !     put core charge (if present) in rhoc(r)
            !
!             if (nlcc_any) call set_cc(irb,eigrb,rhoc) ! warning:giovanni commented this
            !
!             vpot = rhor
            !
!$$
!             CALL start_clock( 'vofrho4' )
!$$
!             call vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
!                          ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion)
!$$
!             CALL stop_clock( 'vofrho4' )
!$$

!$$
            if(do_orbdep) then
                !warning:giovanni don't we need becm down here??? otherwise problems with ultrasoft
                call nksic_potential( n_emps, n_empx, cm_emp, faux_emp, becm_emp, rhovan_emp, deeq_sic_emp, &
                ispin_emp, iupdwn_emp, nupdwn_emp, rhor, rhoc, wtot, sizwtot, vsic_emp, .false., & 
                pink_emp, nudx_emp, wfc_centers_emp, &
                wfc_spreads_emp, icompute_spread, .true. )
!                    call nksic_potential( nbsp, nbspx, cm, fsic, bec, rhovan, deeq_sic, &
!                                       ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, do_wxd, pink, nudx, &
!                                       wfc_centers, wfc_spreads, &
!                                       icompute_spread)
                !
                eodd_emp = sum(pink_emp(1:n_emps))
!                 write(6,*) eodd_emp, etot_emp, "EODD4", pink_emp
                etot_emp = etot_emp + eodd_emp
                !
            endif
!$$
!             if( do_hf ) then
!                 !
!                 call hf_potential( nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
!                                    nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
!                                    rhor, rhog, vxxpsi, exx)
!                 !
!                 etot = etot + sum(exx(1:nbsp))
!                 !
!             endif
            !
            enever=etot_emp
            !
          enddo
!$$
          if (ionode) write(stdout,"(2x,a,i5)") 'iter3 = ',iter3
!$$

!$$
          !if(.not.do_orbdep) then
              if(iter3 == maxiter3 .and. enever.gt.ene0) then
                write(stdout,"(2x,a)") 'missed minimum: iter3 = maxiter3'
                write(stdout,*) enever, ene0
!                 if(non_ortho) then
!                    call compute_duals(c0,cdual,nbspx,1)
!                    call calbec(1,nsp,eigr,cdual,becdual)
!                    write(6,*) "checkdual", cdual(1:2,1)
!                 endif
              else if(enever.le.ene0) then
                c0_emp(:,:)=cm_emp(:,:)
                call copy_twin(bec_emp,becm_emp)
              endif

          !endif
!$$

          restartcg=.true.
          ene_ok=.false.

!$$
          if(iter3 == maxiter3) then
            passof=passop
          endif
!$$
        end if
        
 
        if(.not. ene_ok) THEN
!         call calbec (1,nsp,eigr,c0,bec)
           CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
        ENDIF

        !calculates phi for pc_daga
        CALL calphi( c0_emp, SIZE(c0_emp,1), bec_emp, nhsa, betae, phi_emp, n_emps, lgam )
  
        !=======================================================================
        !                 end of the inner loop
        !=======================================================================
        !
!        if ( ( mod( itercg, isave ) == 0 ) ) then
!            !
!            CALL writefile( h, hold ,nfi, c0, c0old, taus, tausm,  &
!                            vels, velsm, acc, lambda, lambdam, xnhe0, xnhem,     &
!                            vnhe, xnhp0, xnhpm, vnhp, nhpcl,nhpdim,ekincm, xnhh0,&
!                            xnhhm, vnhh, velh, fion, tps, z0t, f, rhor )
!            !
!        endif
        !
        !=======================================================================
        !                 end write to file
        !=======================================================================
  
        itercg=itercg+1
!$$
        itercgeff=itercgeff+1
!$$
        !
        call stop_clock( "outer_loop" )

      enddo OUTER_LOOP

      ! 
      !=======================================================================
      !                 end of the main loop
      !=======================================================================

      !
      !calculates atomic forces and lambda
      !

      !
      ! if pressure is need the following is written because of caldbec
      !
!       if(tpre) then
!          !
!          call  calbec(1,nsp,eigr,c0,bec)
!          !
!              call  caldbec( ngw, nhsa, nbsp, 1, nsp, eigr, c0, dbec )
!                 call rhoofr(nfi,c0(:,:),irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
! 
!          !calculates the potential
!          !
!          !     put core charge (if present) in rhoc(r)
!          !
!          if (nlcc_any) call set_cc(irb,eigrb,rhoc)
! 
!          !
!          !---ensemble-DFT
!          !
!          vpot = rhor
! !$$
!          CALL start_clock( 'vofrho5' )
! !$$
!          call vofrho(nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
!                      ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion)
! !$$
!          CALL stop_clock( 'vofrho5' )
! !$$
! 
! !$$
! !$$ Why there are not other terms here???
! !$$
! 
! !$$
         if(do_orbdep) then
!              !
                call nksic_potential( n_emps, n_empx, c0_emp, faux_emp, bec_emp, rhovan_emp, deeq_sic_emp, &
                                   ispin_emp, iupdwn_emp, nupdwn_emp, rhor, rhoc, wtot, sizwtot, vsic_emp, .false., &
                                   pink_emp, nudx_emp, wfc_centers_emp, wfc_spreads_emp, &
                                   icompute_spread, .false.)
             eodd_emp = sum(pink_emp(1:n_empx))
             write(6,*) eodd_emp, etot_emp, "EODD5", pink_emp
             etot_emp = etot_emp + eodd_emp
!              !
         endif
! !$$
!          if( do_hf ) then
!              !
!              call hf_potential( nbsp, nbspx, c0, f, ispin, iupdwn, nupdwn, &
!                                 nbsp, nbspx, c0, f, ispin, iupdwn, nupdwn, &
!                                 rhor, rhog, vxxpsi, exx)
!              !
!              etot = etot + sum(exx(1:nbsp))
!              !
!          endif
!          !
!      endif



!      call newd(vpot,irb,eigrb,rhovan,fion)
!         if (tfor .or. tprnfor) call nlfq(c0,eigr,bec,becdr,fion)
  
!      call prefor(eigr,betae)
!$$
     ! faux takes into account spin multiplicity.
     !
     faux(1:n_empx) = f_emp(1:n_empx) * DBLE( nspin ) / 2.0d0
     !
!$$
     !
     do i=1,n_emps,2
!$$
         CALL start_clock( 'dforce2' )
!$$          call dforce( i, bec, betae, c0,c2,c3,rhos, nnrsx, ispin,f,n,nspin)
         call dforce(i,bec_emp,betae,c0_emp,c2,c3,filledstates_potential,nnrsx,ispin_emp,faux,n_emps,nspin)
         !
         CALL start_clock( 'dforce2' )
!$$

!$$
         if ( do_orbdep ) then
             !
             ! faux takes into account spin multiplicity.
             !
             CALL nksic_eforce( i, n_emps, n_empx, vsic_emp, deeq_sic_emp, bec_emp, ngw, c0_emp(:,i), c0_emp(:,i+1), vsicpsi, lgam )
!                 CALL nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, bec, ngw, c0(:,i), c0(:,i+1), vsicpsi, lgam )
             !
                !
                c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
                !
                if( i+1 <= n_emps )   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
                !
             !
         endif
!$$
!          if ( do_hf ) then
!              !
!              c2(:) = c2(:) - vxxpsi(:,i) * faux(i)
!              !
!              if( i+1 <= nbsp )   c3(:) = c3(:) - vxxpsi(:,i+1) * faux(i+1)
!              !
!          endif

         do ig=1,ngw
            gi(ig,  i)=c2(ig)
            if(i+1 <= n_emps) gi(ig,i+1)=c3(ig)
         enddo
         !
         if (lgam.and.ng0.eq.2) then
            gi(1,  i)=CMPLX(DBLE(gi(1,  i)),0.d0)
            if(i+1 <= n_emps) gi(1,i+1)=CMPLX(DBLE(gi(1,i+1)),0.d0)
         endif

     enddo

     CALL nlsm1 ( n_emps, 1, nsp, eigr, gi, becm_emp, 1, lgam )

     do iss=1,nspin
        !
        in     = iupdwn(iss)
        in_emp = iupdwn_emp(iss)
        !
        issw   = iupdwn( iss )
        !
        CALL gram_empty(.true. , eigr, betae, becm_emp, bec, nhsa, &
                          gi( :, in_emp: ), c0( :, issw: ), &
                          ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, issw)
        
     enddo
          
     IF(.not.lambda_emp(1)%iscmplx) THEN
        allocate(lambda_repl(nudx_emp,nudx_emp))
     ELSE
        allocate(lambda_repl_c(nudx_emp,nudx_emp))
     ENDIF
     !
        hitmp(:,:) = c0_emp(:,:)
     !
     do is = 1, nspin
        !
        nss = nupdwn_emp(is)
        istart = iupdwn_emp(is)
        
        IF(.not.lambda_emp(1)%iscmplx) THEN
           lambda_repl = 0.d0
        ELSE
           lambda_repl_c = CMPLX(0.d0,0.d0)
        ENDIF
        !
        !
           !
           do i = 1, nss
              do j = i, nss
                 ii = i + istart - 1
                 jj = j + istart - 1
                 IF(.not.lambda_emp(1)%iscmplx) THEN
                    do ig = 1, ngw
                       lambda_repl( i, j ) = lambda_repl( i, j ) - &
                       2.d0 * DBLE( CONJG( hitmp( ig, ii ) ) * gi( ig, jj) )
                    enddo
                    if( ng0 == 2 ) then
                       lambda_repl( i, j ) = lambda_repl( i, j ) + &
                       DBLE( CONJG( hitmp( 1, ii ) ) * gi( 1, jj ) )
                    endif
                    lambda_repl( j, i ) = lambda_repl( i, j )
                 ELSE
                    do ig = 1, ngw
                       lambda_repl_c( i, j ) = lambda_repl_c( i, j ) - &
                       CONJG( hitmp( ig, ii ) ) * gi( ig, jj)
                    enddo
                    lambda_repl_c( j, i ) = CONJG(lambda_repl_c( i, j ))
                 ENDIF
              enddo
           enddo
           !
           
        !
        IF(.not.lambda_emp(1)%iscmplx) THEN
           CALL mp_sum( lambda_repl, intra_image_comm )
           CALL distribute_lambda( lambda_repl, lambda_emp(is)%rvec( :, :), desc_emp( :, is ) )
        ELSE
           CALL mp_sum( lambda_repl_c, intra_image_comm )
           CALL distribute_lambda( lambda_repl_c, lambda_emp(is)%cvec( :, :), desc_emp( :, is ) )
        ENDIF
        !
        !
     end do

!     write(6,*) lambda_emp(1)%cvec, "ldambdaemp", lambda_emp(2)%cvec
     
     IF(.not.lambda_emp(1)%iscmplx) THEN
        DEALLOCATE( lambda_repl )
     ELSE
        DEALLOCATE( lambda_repl_c )
     ENDIF
        !
! write(6,*) "nlfl_twin"
!      call nlfl_twin(bec,becdr,lambda,fion, lgam)
     ! bforceion adds the force term due to electronic berry phase
     ! only in US-case
          
     !
     deallocate(hpsi0,hpsi,gi,hi)
     deallocate(hitmp, STAT=ierr)
     write(6,*) "deallocated hitmp", ierr
     !        
     call deallocate_twin(s_minus1)
     call deallocate_twin(k_minus1)
     
     call stop_clock('runcg_uspp')

!         deallocate(bec0,becm,becdrdiag)

     !begin_modified:giovanni
     call deallocate_twin(bec0_emp)
     call deallocate_twin(becm_emp)
     !
     deallocate(rhor_emp, rhos_emp, rhog_emp)
     !
     if(allocated(rhoc_emp)) then
        !
        deallocate(rhoc_emp)
        !
     endif
!         do i=1,nspin
! 	  call deallocate_twin(lambda(i))
! 	  call deallocate_twin(lambdap(i))
!         enddo
     !
     !end_modified:giovanni

!      deallocate(ave_ene)
     deallocate(c2,c3, faux, faux_emp)

          return

     contains
     
     subroutine v_times_rho(v, nspin, rhos_emp, epot_emp)
     
        use kinds, only: DP
        use mp,                   only : mp_sum
        use mp_global,            only : intra_image_comm
        use smooth_grid_dimensions,      only : nnrsx, nr1s, nr2s, nr3s
        use cell_base,            only : omega
     
        implicit none
        
        integer, intent(in) :: nspin
        real(DP), intent(in) :: v(nnrsx,nspin), rhos_emp(nnrsx,nspin)
        real(DP), intent(out) :: epot_emp
        
        integer :: i
        real(DP) :: etemp, fact, rhosum(2)
        
        etemp=0.d0
        rhosum=0.d0
        fact=omega/DBLE(nr1s*nr2s*nr3s)
        
        do i=1,nspin
           !
           etemp = etemp+sum( v(1:nnrsx,i) * rhos_emp(1:nnrsx,i) )
           rhosum(i) =  sum( rhos_emp(1:nnrsx,i) )
           !
        enddo
        
        call mp_sum(etemp, intra_image_comm)
        call mp_sum(rhosum, intra_image_comm)
        
        rhosum=rhosum*fact
        epot_emp=etemp*fact
        
!         write(6,*) "integrated charge of empty", rhosum
!         write(6,*) "integrated energy", epot_emp
        
        return
     
     end subroutine v_times_rho
     
     END SUBROUTINE runcg_uspp_emp

