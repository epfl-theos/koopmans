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
      !
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
      use constants,                only : pi, au_gpa, e2
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
                                           pink, do_wxd, sizwtot, do_bare_eigs, innerloop_until, &
                                           valpsi, odd_alpha
      use hfmod,                    only : do_hf, vxxpsi, exx
      use twin_types !added:giovanni
      use control_flags,            only : non_ortho
      use cp_main_variables,        only : becdual, becmdual, overlap, ioverlap, becstart
      use electrons_module,         only : wfc_spreads, wfc_centers, icompute_spread, manifold_overlap
      use ldau,                     only : lda_plus_u, vupsi
      use printout_base,            only : printout_base_open, printout_base_unit, &
                                           printout_base_close
      use control_flags,            only : iprint_manifold_overlap, iprint_spreads
      use input_parameters,         only : fixed_state, fixed_band, odd_nkscalfact, do_outerloop, &
                                           finite_field_introduced, finite_field_for_empty_state
      !
      implicit none
      !
      CHARACTER(LEN=80) :: uname
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      integer, EXTERNAL :: get_clock
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
      integer     :: i, j, ig, k, is, iss,ia, iv, jv, il, ii, jj, kk, ip, isp
      integer     :: inl, jnl, niter, istart, nss, nrl, me_rot, np_rot , comm
      real(dp)    :: enb, enbi, x
      real(dp)    :: entmp, sta
      complex(dp) :: gamma_c  !warning_giovanni, is it real anyway?
      complex(dp), allocatable :: c2(:), c3(:), c2_bare(:), c3_bare(:)
      complex(dp), allocatable :: hpsi(:,:), hpsi0(:,:), gi(:,:), hi(:,:), gi_bare(:,:)
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
      complex(DP) :: Omattot(nbspx,nbspx)
      real(DP)    :: dtmp, temp
      real(DP)    :: etot_tmp1, etot_tmp2,  tmppasso
      !
      logical :: lgam, switch=.false., ortho_switch=.false., okvan, steepest=.false.
      complex(DP) :: phase
      integer :: ierr, northo_flavor
      real(DP) :: deltae,sic_coeff1, sic_coeff2 !coefficients which may change according to the flavour of SIC
      integer :: me, iunit_manifold_overlap, iunit_spreads
      character(len=10) :: tcpu_cg_here
      real(DP) :: charge
      !
      real(dp) :: rPi, uPi, eff_finite_field
      real(dp), allocatable :: rho_init(:,:), dvpot(:)
      complex(dp), allocatable :: dvpotpsi(:,:)
      real(dp) :: exxdiv, mp1
      !
      real(dp), external :: exx_divergence
      !
      !
      call do_allocation_initialization()
      !
      ! Initializing clock for minimization
      !
      call start_clock('runcg_uspp')

      if( tfirst .and. ionode ) &
         write(stdout,"(/,a,/)") 'PERFORMING CONJUGATE GRADIENT MINIMIZATION OF EL. STATES'
      !         
      ! set tpa mass preconditioning
      !
      call emass_precond_tpa( ema0bg, tpiba2, emass_cutoff )
      ! 
      call prefor(eigr,betae) 
      !
      ! orthonormalize c0
      !
      call orthogonalize_wfc_only(c0, bec)
      !
      ! recompute phi (the augmented wave function) from the new c0
      !
      CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, nbsp, lgam)
      !
      ! calculates the factors for S and K inversion in US case -- they are important for preconditioning
      ! see paper by Hasnip and Pickard, Computer Physics Communications 174 (2006) 24–29
      !
      if ( okvan ) then
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
      !
      allocate( hpsi(ngw,nbsp) )
      allocate( hpsi0(ngw,nbsp) )
      allocate( gi(ngw,nbsp), hi(ngw,nbsp) )
      !
      allocate(hitmp(ngw,nbsp))
      hitmp(:,:) = CMPLX(0.d0,0.d0)
      !
      gi(:,:)=CMPLX(0.d0,0.d0)
      hi(:,:)=CMPLX(0.d0,0.d0)
      ene_ok=.false.
      !
      ! set occupation for the fixed state ... linh
      !
      if (fixed_state) then
         !
         do i = 1, nbspx
            if (i == fixed_band) f(i) = f_cutoff
         enddo
         ! 
      endif
      !    
      !=======================================================================
      !                 begin of the main loop
      !=======================================================================
      !
      OUTER_LOOP: &
      do while ( itercg < maxiter .and. (.not. ltresh) )
        !
        call start_clock( "outer_loop" )
        ! 
        ENERGY_CHECK: &
        if(.not. ene_ok ) then
          ! 
          call calbec(1,nsp,eigr,c0,bec)
          !  
          call rhoofr(nfi,c0(:,:),irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
          !
          ! put core charge (if present) in rhoc(r)
          !
          if (nlcc_any) call set_cc(irb,eigrb,rhoc)
          !
          ! ensemble-DFT
          !
          vpot = rhor
          !
          CALL start_clock( 'vofrho1' )
          !
          call vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast, &
                      ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
          !
          CALL stop_clock( 'vofrho1' )
          !
          if( tefield  ) then
             ! 
             call berry_energy( enb, enbi, bec, c0(:,:), fion )
             !
             etot=etot+enb+enbi
             !
          endif
          !
          if( do_orbdep ) then
              !
              if (odd_nkscalfact) then
                 !  
                 valpsi(:,:) = (0.0_DP, 0.0_DP)
                 odd_alpha(:) = 0.0_DP
                 !
                 call odd_alpha_routine(c0, nbsp, nbspx, lgam, .false.)
                 !
              endif
              !
              fsic = f
              !
              call nksic_potential( nbsp, nbspx, c0, fsic, bec, rhovan, deeq_sic, &
                                    ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, wfc_centers, &
                                    wfc_spreads, icompute_spread, .false. )
              ! 
              eodd=sum(pink(1:nbsp))
              etot=etot + eodd
              !
          endif
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
          if (finite_field_introduced) then 
              !
              if (itercg == 1) then
                  !
                  allocate (rho_init(nnr,nspin))
                  allocate (dvpotpsi(ngw,2))
                  allocate (dvpot(nnr))
                  !
                  ! 1) compute dvpot only for the first time

                 !call perturbing_pot_nic(fixed_band, ispin(fixed_band), rhor, dvpot)
                  if (.false.) then ! here for Koopmans
                     !
                     call perturbing_pot(fixed_band, ispin(fixed_band), rhor, dvpot, uPi,lgam, finite_field_for_empty_state)
                     !
                  endif
                  !
                  if (.true.) then ! here for BSE
                     ! 
                     !call perturbing_pot_bse(fixed_band, fixed_band_1, ispin(fixed_band), ispin(fixed_band_1), rhor, dvpot, uPi,lgam)
                     !
                  endif
                  !
                  ! 2) save rhor to rhor_starting
                  rho_init(:,:) = rhor(:,:)
                  !
              endif
              !
#ifdef __TOBE_FIXED
              ! compute int dvpot rhor 
              call compute_effective_energy(dvpot, rhor, eff_finite_field)
              !
              etot = etot + eff_finite_field
#endif
              !
          endif
          ! 
          etotnew=etot
          !
        else
          !
          etot=enever
          !
          etotnew=etot
          ene_ok=.false.
          !
        endif ENERGY_CHECK
        !
        if( do_orbdep ) then
           !
           call do_innerloop_subroutine()
           !
        endif
        !
        call print_out_observables()
        !
        IF ( .not. do_outerloop ) THEN
           EXIT OUTER_LOOP
        ENDIF
        !
        ! here we store the etot in ene0, to keep track of the energy of the initial point
        !
        call check_convergence_cg()        
        !
        call newd(vpot,irb,eigrb,rhovan,fion)
        call prefor(eigr,betae)!ATTENZIONE
        !
        call compute_hpsi()
        !
        if(pre_state) call ave_kin(c0,SIZE(c0,1),nbsp,ave_ene)
        ! 
        ! HPSI IS ORTHOGONALIZED TO c0
        !
        if(switch.or.(.not.do_orbdep)) then
           !
           if (fixed_state) then
              !
              hpsi(:, fixed_band) = cmplx(0.d0, 0.d0)  
              ! 
           endif
           !
           call pcdaga2(c0,phi,hpsi, lgam)
           !
        else
           !
           if(.not.okvan) then
             !
             if (fixed_state) then
                !
                call pc3nc_fixed(c0, hpsi, lgam)
                !
             endif
             ! 
             call pc3nc(c0,hpsi,lgam)
             !
           else
             !
             call pc3us(c0,bec,hpsi,becm,lgam)
             !
           endif
           !
        endif
        !
        ! TWO VECTORS INITIALIZED TO HPSI
        ! 
        hpsi0(1:ngw,1:nbsp) = hpsi(1:ngw,1:nbsp) 
        !
        gi(1:ngw,1:nbsp)    = hpsi(1:ngw,1:nbsp)
        !
	! COMPUTES ULTRASOFT-PRECONDITIONED HPSI,
        ! non kinetic-preconditioned, 
        ! is the subsequent reorthogonalization necessary 
        ! in the norm conserving case???: giovanni
        ! 
        call calbec(1,nsp,eigr,hpsi,becm)
        call xminus1_twin(hpsi,betae,dumm,becm,s_minus1,.false.)
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! look if the following two lines are really needed
        !
        call orthogonalize(c0,hpsi,becm,bec)
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  COMPUTES ULTRASOFT+KINETIC-preconditioned GI
        ! 
        IF (do_orbdep) THEN 
           !
           ! preconditioning with respect to spreads.
           ! the gradient along wavefunctions with the largest 
           ! localization is a bit shortened
           !
           !do i=1,nbsp
           !   gi(:,i) = gi(:,i)*(1.d0+1.d0/sqrt(wfc_spreads(1+i-iupdwn(ispin(i)),ispin(i),2)))
           !enddo
           !
        ENDIF
        !
        if (.not.pre_state) then
           !
           call xminus1_twin(gi,betae,ema0bg,becm,k_minus1,.true.)
           !
        else
           !
           ! warning:giovanni not yet implemented
           !
           call xminus1_state(gi,betae,ema0bg,becm,k_minus1,.true.,ave_ene) 
           ! 
        endif
        !
        call orthogonalize(c0,gi,becm,bec)
        !    
        !  calculates gamma
        !
        gamma_c=CMPLX(0.d0,0.d0)
        !
        DO i=1,nbsp
           !
           IF (lgam) THEN
              !
              do ig=1,ngw
                ! 
                gamma_c=gamma_c+2.d0*DBLE(CONJG(gi(ig,i))*hpsi(ig,i))
                ! 
              enddo
              !
              if (ng0.eq.2) then
                 !
                 gamma_c=gamma_c-DBLE(CONJG(gi(1,i))*hpsi(1,i))
                 !
              endif
              !
           ELSE
              !
              do ig=1,ngw
                 !
                 gamma_c=gamma_c+CONJG(gi(ig,i))*hpsi(ig,i)
                 ! 
              enddo
              !
           ENDIF
           !
        ENDDO
        !        
        call mp_sum( gamma_c, intra_image_comm )
        ! 
        if (nvb.gt.0) then
           !
           if (.not.becm%iscmplx) then
              ! 
              do i=1,nbsp
                 ! 
                 do is=1,nvb
                    !
                    do iv=1,nh(is)
                       !
                       do jv=1,nh(is)
                          !
                          do ia=1,na(is)
                             !
                             inl=ish(is)+(iv-1)*na(is)+ia
                             !
                             jnl=ish(is)+(jv-1)*na(is)+ia
                             ! 
                             gamma_c=gamma_c+ qq(iv,jv,is)*becm%rvec(inl,i)*bec0%rvec(jnl,i)
                             !
                          enddo
                          !
                       enddo
                       !
                    enddo
                    !
                 enddo
                 !
              enddo
              !
           else
              !
              do i=1,nbsp
                 !
                 do is=1,nvb
                    ! 
                    do iv=1,nh(is)
                       !
                       do jv=1,nh(is)
                          !
                          do ia=1,na(is)
                             !
                             inl=ish(is)+(iv-1)*na(is)+ia
                             !    
                             jnl=ish(is)+(jv-1)*na(is)+ia
                             !
                             gamma_c=gamma_c+ qq(iv,jv,is)*CONJG(becm%cvec(inl,i))*(bec0%cvec(jnl,i)) 
                             !
                          enddo
                          !
                       enddo
                       !
                    enddo
                    !
                 enddo
                 !
              enddo
              !
           endif
           !
        endif
        ! 
        ! case of first iteration
        ! 
	gamma_c=CMPLX(DBLE(gamma_c),0.d0)
        !
        if (steepest) then
           !
           ! steepest descent
           !  
           gamma_c=0.d0
           !
        endif
        !
        if (itercg==1 .or. mod(itercg,niter_cg_restart)==0 .or. restartcg) then
           !
           restartcg=.false.
           !
           !  We do not have to reset passof every exception of CG!
           ! warning:giovanni if we do not reset we may have fake convergences!!!!
           !
           passof=passop
           !
           ! hi is the search direction
           !
           hi(1:ngw,1:nbsp)=gi(1:ngw,1:nbsp)
           ! 
           esse_c=gamma_c
           ! 
        else
           !
           ! find direction hi for general case 
           ! calculates gamma for general case, not using Polak Ribiere
           !
           if (.not.steepest) then
              !
              essenew_c=gamma_c
              !
              gamma_c=gamma_c/esse_c
              ! 
              esse_c=essenew_c
              !
           else
              !    
              esse_c=0.d0
              !
              essenew_c=0.d0
              !
           endif
           !   
           hi(1:ngw,1:nbsp)=gi(1:ngw,1:nbsp)+(gamma_c)*hi(1:ngw,1:nbsp)
           !   
        endif
        !
        ! note that hi is saved on gi, because we need it before projection on conduction states
        !     
        ! ... find minimum along direction hi:
        !
        ! project hi on conduction sub-space
        !
        call orthogonalize(c0,hi,bec0,bec)
        !    
        ! do quadratic minimization
        !             
        ! calculate derivative with respect to lambda along direction hi
        !
        dene0=0.d0
        ! 
        do i=1,nbsp
           !  
           IF (lgam) THEN
              !              
              do ig=1,ngw
                 !
                 dene0=dene0-4.d0*DBLE(CONJG(hi(ig,i))*hpsi0(ig,i))
                 !
              enddo
              !
              if (ng0.eq.2) then
                 ! 
                 dene0=dene0+2.d0*DBLE(CONJG(hi(1,i))*hpsi0(1,i))
                 !   
              endif
              !
           ELSE
              !
              do ig=1,ngw
                 !  
                 dene0=dene0-2.d0*DBLE(CONJG(hi(ig,i))*hpsi0(ig,i))
                 !
              enddo
              ! 
           ENDIF
           !
        enddo
        !
        ! We need the following because n for spin 2 is double that for spin 1!
        !
        dene0 = dene0 *2.d0/nspin
        !
        call mp_sum( dene0, intra_image_comm )
        !
        ! if the derivative is positive, search along opposite direction
        ! 
        if (dene0.gt.0.d0) then
           !
           spasso=-1.D0
           ! 
        else
           !
           spasso=1.d0
           !
        endif
        !
        ! Thank Nicola for providing it. 
        !
        if (.false.) then
           !
           etot_tmp1 = 0.0
           etot_tmp2 = 0.0
           !
           tmppasso = 0.0
           do i=1,2
              !
              if (i==1)  tmppasso= 1.d-5
              if (i==2)  tmppasso=-1.d-5
              !tmppasso =tmppasso+0.00001
              !
              cm(:,:) = c0(:,:) + spasso * tmppasso * hi(:,:)
              !
              if (lgam.and.ng0 == 2)  cm(1,:)=0.5d0*(cm(1,:)+CONJG(cm(1,:)))
              !
              ! orthonormalize
              !
              call orthogonalize_wfc_only(cm,becm)
              ! 
              ! **** calculate energy ene1
              !  
              call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
              !
              ! calculate potential
              !
              !  put core charge (if present) in rhoc(r)
              !
              if (nlcc_any) call set_cc(irb,eigrb,rhoc)
              !
              vpot = rhor
              !
              CALL vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast,  &
                          ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
              !
              if ( tefield  ) then
                 !
                 call berry_energy( enb, enbi, becm, cm(:,:), fion )
                 ! 
                 etot=etot+enb+enbi
                 ! 
              endif
              !
              if (do_orbdep) then
                 !
                 if (odd_nkscalfact) then
                    !  
                    valpsi(:,:) = (0.0_DP, 0.0_DP)
                    odd_alpha(:) = 0.0_DP
                    !
                    call odd_alpha_routine(cm, nbsp, nbspx, lgam, .false.)
                    !
                 endif
                 !
                 ! warning:giovanni don't we need becm down here??? otherwise
                 ! problems
                 ! with ultrasoft!!
                 !
                 call nksic_potential( nbsp, nbspx, cm, fsic, becm, rhovan, deeq_sic,&
                                       ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, &
                                       wfc_centers, wfc_spreads, &
                                       icompute_spread, .false.)
                 !
                 eodd=sum(pink(1:nbsp))
                 !
                 etot = etot + eodd
                 !
              endif
              !
              if ( do_hf ) then
                 !
                 call hf_potential( nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
                                    nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
                                    rhor, rhog, vxxpsi, exx)
                 !
                 etot = etot + sum(exx(1:nbsp))
                 !
              endif
              !
#ifdef __TOBE_FIXED
              if (finite_field_introduced) then
                 !
                 ! compute int dvpot rhor 
                 call compute_effective_energy(dvpot, rhor, eff_finite_field)
                 !
                 etot = etot + eff_finite_field
                 !
              endif
#endif
              !
              write(stdout,*) "etot: ", i, "=", etot
              !
              if (i==1) etot_tmp1 = etot
              if (i==2) etot_tmp2 = etot
              !  
          enddo
          !
          write(stdout,*) "here is numerical derivative vs analytic derivative at step", itercg
          write(stdout,*) "(etot_emp_tmp1-etot_emp_tmp2)/tmppasso, dene0, tmppasso,  &
                           ((etot_emp-etot_emp_tmp2)/(2*tmppasso)/dene0)"
          write(stdout,'(2e25.15,4e20.10)') (etot_tmp1-etot_tmp2)/(2.0*tmppasso), dene0, &
                tmppasso, ((etot_tmp1 -etot_tmp2)/(2.0*tmppasso)/dene0)
          !
        endif
        !
        ! calculates wave-functions on a point on direction hi
        !
        cm(1:ngw,1:nbsp)=c0(1:ngw,1:nbsp)+spasso*passof*hi(1:ngw,1:nbsp)
        !
        !
        ! I do not know why the following 3 lines 
        ! were not in the original code (CHP)
        !
        if (lgam.and.ng0 == 2)  cm(1,:)=0.5d0*(cm(1,:)+CONJG(cm(1,:)))
        !  
        ! orthonormalize
        !
        call orthogonalize_wfc_only(cm,becm)
        ! 
        ! **** calculate energy ene1
        !  
        call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
        !
        ! calculate potential
        !
        !  put core charge (if present) in rhoc(r)
        !
        if (nlcc_any) call set_cc(irb,eigrb,rhoc)
        !
        vpot = rhor
        !
        CALL start_clock( 'vofrho2' )
        !
        CALL vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast,             &
                     ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
        !
        CALL stop_clock( 'vofrho2' )
        !
        if ( tefield  ) then
           !
           call berry_energy( enb, enbi, becm, cm(:,:), fion )
           ! 
           etot=etot+enb+enbi
           ! 
        endif
        !
        if (do_orbdep) then
           !
           if (odd_nkscalfact) then
              !  
              valpsi(:,:) = (0.0_DP, 0.0_DP)
              odd_alpha(:) = 0.0_DP
              !
              call odd_alpha_routine(cm, nbsp, nbspx, lgam, .false.)
              !
           endif
           !
           ! warning:giovanni don't we need becm down here??? otherwise problems with ultrasoft!!
           !
           call nksic_potential( nbsp, nbspx, cm, fsic, becm, rhovan, deeq_sic, &
                               ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, &
                               wfc_centers, wfc_spreads, &
                               icompute_spread, .false.)
           !
           eodd=sum(pink(1:nbsp))
           !
           etot = etot + eodd
           !
        endif
        !
        if ( do_hf ) then
           !
           call hf_potential( nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
                              nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
                              rhor, rhog, vxxpsi, exx)
           !
           etot = etot + sum(exx(1:nbsp))
           !
        endif
        !
#ifdef __TOBE_FIXED
        if (finite_field_introduced) then
           !
           ! compute int dvpot rhor 
           call compute_effective_energy(dvpot, rhor, eff_finite_field)
           !
           etot = etot + eff_finite_field
           !
        endif
#endif
        !
        ene1=etot
        !    
        ! find the minimum
        !
        call minparabola(ene0,spasso*dene0,ene1,passof,passo,enesti)
        !
        if( ionode .and. iprsta > 1 ) write(stdout,"(6f20.12)") ene0,dene0,ene1,passo, DBLE(gamma_c), esse_c
        !  
        ! set new step
        !
        passov=passof
        !
        ! doing the following makes the convergence better...
        !
        passof=passo
        !      
        ! calculates wave-functions at minimum
        !
        cm(1:ngw,1:nbsp) = c0(1:ngw,1:nbsp) +spasso*passo*hi(1:ngw,1:nbsp)
        !
        IF (lgam.and. ng0 == 2 ) cm(1,:) = 0.5d0*(cm(1,:)+CONJG(cm(1,:)))
        !
        call orthogonalize_wfc_only(cm, becm)
        !
        call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
        !
        ! calculates the potential
        !
        ! put core charge (if present) in rhoc(r)
        !
        if (nlcc_any) call set_cc(irb,eigrb,rhoc)
        !
        vpot = rhor
        !
        CALL start_clock( 'vofrho3' )
        !
        call vofrho(nfi,vpot,rhog,rhos,rhoc,tfirst,tlast, &
                    ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
        !
        CALL stop_clock( 'vofrho3' )
        !
        if (tefield  ) then
           ! 
           call berry_energy( enb, enbi, becm, cm(:,:), fion )
           !  
           etot=etot+enb+enbi
           !
        endif
        !
        if(do_orbdep) then
            !
            if (odd_nkscalfact) then
               !  
               valpsi(:,:) = (0.0_DP, 0.0_DP)
               odd_alpha(:) = 0.0_DP
               !
               call odd_alpha_routine(cm, nbsp, nbspx, lgam, .false.)
               !
            endif
            !
            ! warning:giovanni... don't we need becm down here?? otherwise problem with ultrasoft!!
            !
            call nksic_potential( nbsp, nbspx, cm, fsic, becm, rhovan, deeq_sic, &
                                  ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx,&
                                  wfc_centers, wfc_spreads, &
                                  icompute_spread, .false.)
            eodd = sum(pink(1:nbsp))
            etot = etot + eodd
            !
        endif
        ! 
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
#ifdef __TOBE_FIXED
        if (finite_field_introduced) then
           !
           ! compute int dvpot rhor 
           call compute_effective_energy(dvpot, rhor, eff_finite_field)
           !
           write(stdout,*) "eff_finite_field", eff_finite_field
           etot = etot + eff_finite_field
           !
        endif
#endif
        !
        enever=etot
        !
        ! check with  what supposed
        !
        write(stdout,*) 'ene0, dene0, ene1, enesti,enever, passo, passov, passof'
        write(stdout,"(8f18.12)") ene0, dene0, ene1, enesti,enever, passo, passov, passof

        if(ionode .and. iprsta > 1 ) then
            write(stdout,"(2x,a,f20.12)") 'cg_sub: estimate :'  , (enesti-enever)/(ene0-enever)
            write(stdout,"(2x,a,3f20.12)") 'cg_sub: minmum   :'  , enever,passo,passov
        endif
        !
        ! if the energy has diminished with respect to ene0 and ene1 , everything ok
        !
        if (((enever.lt.ene0) .and. (enever.lt.ene1)).or.(tefield.or.tefield2)) then
           !
           c0(:,:)=cm(:,:)
           call copy_twin(bec,becm) !modified:giovanni
           ene_ok=.true.
           !
           ! if  ene1 << energy <  ene0; go to  ene1
           !
        elseif((enever.ge.ene1) .and. (enever.lt.ene0)) then
           !
           if (ionode) then
              ! 
              write(stdout,"(2x,a,i5,f20.12)") 'cg_sub: missed minimum, case 1, iteration',itercg, passof
              ! 
           endif
           ! 
           c0(1:ngw,1:nbsp)=c0(1:ngw,1:nbsp)+spasso*passov*hi(1:ngw,1:nbsp)
           ! 
           passof=2.d0*passov
           !
           restartcg=.true.
           !
           call orthogonalize_wfc_only(c0,bec)
           !
           ene_ok=.false.
           !
           ! if  ene1 << ene0 <= energy; go to  ene1
           !
        elseif((enever.ge.ene0).and.(ene0.gt.ene1)) then
           !   
           if (ionode) then
              write(stdout,"(2x,a,i5)") 'cg_sub: missed minimum, case 2, iteration',itercg
           endif
           ! 
           c0(1:ngw,1:nbsp)=c0(1:ngw,1:nbsp)+spasso*passov*hi(1:ngw,1:nbsp)
           !
           passof=1.d0*passov
           !
           restartcg=.true.!ATTENZIONE
           !
           call orthogonalize_wfc_only(c0,bec)
           !
           ! if ene > ene0, en1 do a steepest descent step
           ! 
           ene_ok=.false.
           !
        elseif((enever.ge.ene0).and.(ene0.le.ene1)) then
           !
           if(ionode) then
             write(stdout,"(2x,a,i5)") 'cg_sub: missed minimum, case 3, iteration',itercg
           endif
           !
           iter3=0
           !  
           do while(enever.ge.ene0 .and. iter3.lt.maxiter3)
             ! 
             iter3=iter3+1
             !  
             passov=passov*0.5d0
             !
             cm(1:ngw,1:nbsp)=c0(1:ngw,1:nbsp)+spasso*passov*hi(1:ngw,1:nbsp)
             !
             passof=1.d0*passov
             ! 
             itercgeff=itercgeff+1
             !
             ! change the searching direction
             !
             spasso=spasso*(-1.d0)
             !
             call orthogonalize_wfc_only(cm,becm)
             !
             call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
             !
             ! calculates the potential
             !
             ! put core charge (if present) in rhoc(r)
             !
             if (nlcc_any) call set_cc(irb,eigrb,rhoc)
             !
             vpot = rhor
             !
             CALL start_clock( 'vofrho4' )
             !
             call vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
                         ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion)
             !
             CALL stop_clock( 'vofrho4' )
             !
             if ( tefield  ) then!to be bettered
                !
                call berry_energy( enb, enbi, becm, cm(:,:), fion )
                !
                etot=etot+enb+enbi
                !
             endif
             !
             if (do_orbdep) then
                !
                if (odd_nkscalfact) then
                   !  
                   valpsi(:,:) = (0.0_DP, 0.0_DP)
                   odd_alpha(:) = 0.0_DP
                   !
                   call odd_alpha_routine(cm, nbsp, nbspx, lgam, .false.)
                   !
                endif
                !
                ! warning:giovanni don't we need becm down here??? otherwise problems with ultrasoft
                ! 
                call nksic_potential( nbsp, nbspx, cm, fsic, becm, rhovan, deeq_sic, &
                                      ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, &
                                      wfc_centers, wfc_spreads, &
                                      icompute_spread, .false.)
                !
                eodd = sum(pink(1:nbsp))
                !
                etot = etot + eodd
                !
             endif
             ! 
             if ( do_hf ) then
                !
                call hf_potential( nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
                                   nbsp, nbspx, cm, f, ispin, iupdwn, nupdwn, &
                                   rhor, rhog, vxxpsi, exx)
                !
                etot = etot + sum(exx(1:nbsp))
                !
             endif
             !
#ifdef __TOBE_FIXED
             if (finite_field_introduced) then
                !
                ! compute int dvpot rhor 
                call compute_effective_energy(dvpot, rhor, eff_finite_field)
                !
                etot = etot + eff_finite_field
                !
             endif
#endif
             !
             enever=etot
             !
             write(stdout, *) iter3, spasso*passov, enever
             !
           enddo
           !
           if (ionode) write(stdout,"(2x,a,i5)") 'iter3 = ',iter3
           !
           if (iter3 == maxiter3 .and. enever.gt.ene0) then
              ! 
              write(stdout,"(2x,a)") 'missed minimum: iter3 = maxiter3'
              write(stdout,*) enever, ene0
              !
           elseif (enever.le.ene0) then
              !
              c0(:,:)=cm(:,:)
              ! 
              call copy_twin(bec,becm)
              !
           endif
           !
           restartcg=.true.
           ene_ok=.false.
           !
           if(iter3 == maxiter3) then
             passof=passop
           endif
           !
        endif
        !  
        if(.not. ene_ok) call calbec (1,nsp,eigr,c0,bec)
        !
        ! calculates phi for pc_daga
        !
        CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, nbsp, lgam )
        ! 
        !=======================================================================
        !                 end of the outer loop
        !=======================================================================
        !
        itercg=itercg+1
        !
        itercgeff=itercgeff+1
        !
        call stop_clock( "outer_loop" )
        !
      enddo OUTER_LOOP
      !
#ifdef __DEBUG
        ! for debug and tuning purposes
        if ( ionode ) write(37,*)itercg, itercgeff, etotnew
        if ( ionode ) write(1037,'("iteration =",I4,"  eff iteration =",I4,"   Etot (Ha) =",F22.14)')&
            itercg, itercgeff, etotnew
#endif
      !
      ! Computes the leading term of the Makov-Payne corrective energy
      ! E = q^2 * madelung_const / ( 2 * L )
      !
      IF ( fixed_state ) THEN
        !
        WRITE( stdout, '(//,A)' ) " -----------------------"
        WRITE( stdout, '(A)' ) " MAKOV-PAYNE CORRECTIONS"
        WRITE( stdout, '(A)' ) " -----------------------"
        exxdiv = exx_divergence()
        !
        ! The following IF loop determines the system charge assuming always 
        ! an even number of electrons for the neutral system. When the number
        ! of electrons is odd (nupdwn(1) .ne. nupdwn(2)), we guess to deal with
        ! an N+1 calculation and the charge is calculated consequently 
        IF ( nupdwn(1) == nupdwn(2) .AND. fixed_band .LE. nupdwn(1) ) THEN
          ! Case N-1
          charge = 1 - f_cutoff
          WRITE( stdout, '(A,F10.6)' ) " N-1 CASE --- q =", charge
          !
        ELSE IF ( fixed_band .GT. nupdwn(1) .OR. fixed_band .GT. nupdwn(2) ) THEN
          ! Case N+1
          charge = - f_cutoff
          WRITE( stdout, '(A,F10.6)' ) " N+1 CASE --- q =", charge
          !
        ELSE
          !
          charge = 1.D0
          WRITE( stdout, '(A)' ) " Cannot understand which case --- q set to 1"
          !
        ENDIF
        !
        mp1 = - exxdiv / omega * charge**2 / 2
        mp1 = mp1 / e2       ! Ry to Ha conversion
        WRITE( stdout, '(/,2X,A,ES20.8)' ) " Makov-Payne 1-order energy : ", mp1
        !
      ENDIF
      !
      ! OBSOLETE: old version for calculating correction for charged systems
      !           using Makov-Payne corrections (simple cubic systems only!)
      !IF (fixed_state .and. do_orbdep) THEN
      !   !
      !   write(stdout, *) "NLN: This is 2nd term in MP formular, computed with localized orbital"
      !   write(stdout, *) "NLN: Use for extended system only where tcc does not have the correction"
      !   !
      !   charge = 1.0_dp
      !   !
      !   call makov_payne_correction_2nd_term ( charge, wfc_spreads(fixed_band, 1, 1))
      !   !
      !ENDIF
      ! 
      !=======================================================================
      !                 end of the main loop
      !=======================================================================
      !
      ! calculates atomic forces and lambda matrix elements
      !
      ! if pressure is need the following is written because of caldbec
      !
      if(tpre) then
         !
         call calbec(1,nsp,eigr,c0,bec)
         !
         call caldbec( ngw, nhsa, nbsp, 1, nsp, eigr, c0, dbec )
         call rhoofr(nfi,c0(:,:),irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
         !
         !calculates the potential
         !
         !     put core charge (if present) in rhoc(r)
         !
         if (nlcc_any) call set_cc(irb,eigrb,rhoc)
         !
         !---ensemble-DFT
         !
         vpot = rhor
         !
         CALL start_clock( 'vofrho5' )
         !
         call vofrho(nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
                     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion)
         !
         CALL stop_clock( 'vofrho5' )
         !
         if(do_orbdep) then
             !
             if (odd_nkscalfact) then
                !  
                valpsi(:,:) = (0.0_DP, 0.0_DP)
                odd_alpha(:) = 0.0_DP
                !
                call odd_alpha_routine(cm, nbsp, nbspx, lgam, .false.)
                !
             endif
             !
             call nksic_potential( nbsp, nbspx, c0, fsic, bec, rhovan, deeq_sic, &
                                ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, &
                                wfc_centers, wfc_spreads, &
                                icompute_spread, .false.)
             !
             eodd = sum(pink(1:nbsp))
             etot = etot + eodd
             !
         endif
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
#ifdef __TOBE_FIXED
         if (finite_field_introduced) then
            !
            ! compute int dvpot rhor 
            call compute_effective_energy(dvpot, rhor, eff_finite_field)
            !
            etot = etot + eff_finite_field
            !
         endif
#endif
         !
     endif
     !
#ifdef __TOBE_FIXED
     if (finite_field_introduced) then
        !
        ! 1) compute drho = rho_final - rho_init
        ! 2) compute rPi = int{ dvpot * drho()}
        ! 3) compute alpha_i = 1 + rPi/uPi
        !
        rho_init(:,:) = rhor(:,:) - rho_init(:,:)
        !
        call compute_effective_energy(dvpot, rho_init, rPi)
        !     
        WRITE(stdout,*) "Here is uPi", uPi 
        WRITE(stdout,*) "Here is rPi", rPi 
        WRITE(stdout,*) "Here is 1+rPi/uPi",  1 + rPi/uPi 
        !
        deallocate(rho_init)
        ! 
     endif
#endif
     !
     call newd(vpot,irb,eigrb,rhovan,fion)
     !
     if (tfor .or. tprnfor) call nlfq(c0,eigr,bec,becdr,fion, lgam)
     ! 
     call prefor(eigr,betae)
     ! 
     ! faux takes into account spin multiplicity.
     !
     faux(1:nbsp) = max(f_cutoff,f(1:nbsp)) * DBLE( nspin ) / 2.0d0
     !
     ! This condition to ensure that the orbital(fixed_band) is frozen
     ! 
     IF (fixed_state) THEN
        faux(fixed_band) = f_cutoff
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
         !
         call dforce(i,bec,betae,c0,c2,c3,rhos,nnrsx,ispin,faux,nbsp,nspin)
         !
         CALL start_clock( 'dforce2' )
!$$
         IF(do_bare_eigs) THEN
            !
            c2_bare(:) = c2(:)
            c3_bare(:) = c3(:)
            !
         ENDIF
!$$
        if(tefield.and.(evalue .ne. 0.d0)) then
            !
            call dforceb &
               (c0, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            !
            c2(:)=c2(:)+evalue*df(:)
            !
            call dforceb &
               (c0, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            !
            c3(:)=c3(:)+evalue*df(:)
            !
        endif

         if ( do_orbdep ) then
            !
            IF ( odd_nkscalfact ) THEN
               !
               c2(:) = c2(:) - valpsi(i,:) * faux(i)
               !
               if( i+1 <= nbsp ) c3(:) = c3(:) - valpsi(i+1,:) * faux(i+1)
               ! 
            ENDIF
            !
            ! faux takes into account spin multiplicity.
            !
            CALL nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, bec, ngw, c0(:,i), c0(:,i+1), vsicpsi, lgam )
            !
            !
            c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
            !
            if( i+1 <= nbsp )   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
            !
         endif
!$$
         IF ( lda_plus_u ) THEN
             !
             c2(:) = c2(:) - vupsi(:,i) * faux(i)
             if( i+1 <= nbsp ) c3(:) = c3(:) - vupsi(:,i+1) * faux(i+1)
             !
         ENDIF

         if ( do_hf ) then
             !
             c2(:) = c2(:) - vxxpsi(:,i) * faux(i)
             !
             if( i+1 <= nbsp )   c3(:) = c3(:) - vxxpsi(:,i+1) * faux(i+1)
             !
         endif
!$$
#ifdef __TOBE_FIXED
         if (finite_field_introduced) then
            !
            ! faux takes into account spin multiplicity.
            !
            CALL finite_field_force ( i, nbsp, nbspx, dvpot, ngw, c0(:,i), c0(:,i+1), dvpotpsi, lgam )
            !
            c2(:) = c2(:) + dvpotpsi(:,1) * faux(i)
            !
            if( i+1 <= nbsp )   c3(:) = c3(:) + dvpotpsi(:,2) * faux(i+1)
            !
         endif
#endif
         ! 
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
         !
     enddo
     !
     IF(.not.lambda(1)%iscmplx) THEN
        allocate(lambda_repl(nudx,nudx))
     ELSE
        allocate(lambda_repl_c(nudx,nudx))
     ENDIF
     !
     hitmp(1:ngw,1:nbsp) = c0(1:ngw,1:nbsp)
     !
     do is = 1, nspin
        !
        nss = nupdwn(is)
        istart = iupdwn(is)
        ! 
        IF(.not.lambda(1)%iscmplx) THEN
           lambda_repl = 0.d0
        ELSE
           lambda_repl_c = CMPLX(0.d0,0.d0)
        ENDIF
        !
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
        IF(.not.lambda(1)%iscmplx) THEN
           CALL mp_sum( lambda_repl, intra_image_comm )
           CALL distribute_lambda( lambda_repl, lambda(is)%rvec( :, :), descla( :, is ) )
        ELSE
           CALL mp_sum( lambda_repl_c, intra_image_comm )
           CALL distribute_lambda( lambda_repl_c, lambda(is)%cvec( :, :), descla( :, is ) )
        ENDIF
        !
     end do

     IF(do_bare_eigs) call compute_lambda_bare()

     IF(.not.lambda(1)%iscmplx) THEN
        DEALLOCATE( lambda_repl )
     ELSE
        DEALLOCATE( lambda_repl_c )
     ENDIF
     !
     call nlfl_twin(bec,becdr,lambda,fion, lgam)
     !
     ! bforceion adds the force term due to electronic berry phase
     ! only in US-case          
     ! 
     if( tefield.and.(evalue .ne. 0.d0) ) then
        call bforceion(fion,tfor.or.tprnfor,ipolp, qmat,bec,becdr,gqq,evalue)
     endif
     !
     call do_deallocation()
     !
     return
     !
     contains

     subroutine do_allocation_initialization()
     
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! INITIALIZATION PART (variables, allocation of arrays, minimization parameters)
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !!! Initialize some basic variables
         !
         me = me_image+1
         lgam = gamma_only.and..not.do_wf_cmplx
         okvan=(nvb>0)
         !
         deltae = 2.d0*conv_thr
         etotnew=0.d0
         etotold=0.d0
         eoddnew=0.d0
         eoddold=0.d0
         !
         !!! Initialize printing
         !
         IF( ionode ) THEN
            !
            iunit_spreads = printout_base_unit( "sha" )
            CALL printout_base_open( "sha" )
            WRITE(iunit_spreads, *) " "
            iunit_manifold_overlap = printout_base_unit( "ovp" )
            CALL printout_base_open( "ovp" )
            WRITE(iunit_manifold_overlap, *) " "
            !
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
         if(innerloop_until<0) then
            !
            innerloop_until=2*maxiter
            !
         endif
         !
         allocate (faux(nbspx))
         !
         allocate (ave_ene(nbsp))
         allocate (c2(ngw),c3(ngw))
         !
         ! Allocation of twin-type variables
         !
         call init_twin(bec0, lgam)
         call allocate_twin(bec0, nhsa, nbsp, lgam)

         call init_twin(becm, lgam)
         call allocate_twin(becm, nhsa, nbsp, lgam)

         call init_twin(becdrdiag, lgam)
         call allocate_twin(becdrdiag, nhsa, nspin*nlax,3, lgam)

         !
         ! initializing variables for checking iterations and convergence
         !

         newscheme=.false.
         firstiter=.true.

         pre_state=.false.!normally is disabled

         maxiter3=12
         !$$
         if(do_orbdep) maxiter3=5
         !$$
         ninner=0

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
     
     end subroutine do_allocation_initialization
     
     subroutine do_deallocation()
     
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
   !         call deallocate_twin(lambda(i))
   !         call deallocate_twin(lambdap(i))
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
        if (finite_field_introduced) deallocate(dvpot, dvpotpsi)
        ! 
     end subroutine do_deallocation
     
     subroutine do_innerloop_subroutine()

      if(do_innerloop .and. innerloop_until>=itercgeff) then
!$$$$          if(do_innerloop.and.itercg.le.20) then
!$$$$
         !
         call start_clock( "inner_loop" )
         !
         eodd    = sum(pink(1:nbsp))
         etot    = etot - eodd
         etotnew = etotnew - eodd
         ninner  = 0

         if(.not.do_innerloop_cg) then
           call nksic_rot_emin(itercg,ninner,etot,Omattot, lgam)
         else
           !call nksic_rot_emin_cg(itercg,innerloop_init_n,ninner,etot,Omattot,deltae*innerloop_cg_ratio,lgam)
           call nksic_rot_emin_cg_general(itercg,innerloop_init_n,ninner,etot,deltae*innerloop_cg_ratio,lgam, &
                                     nbsp, nbspx, nudx, iupdwn, nupdwn, ispin, c0, rhovan, bec, rhor, rhoc, &
                                     vsic, pink, deeq_sic, wtot, fsic, sizwtot, do_wxd, wfc_centers, wfc_spreads, .false.)


         endif

         eodd    = sum(pink(1:nbsp))
         etot    = etot + eodd
         etotnew = etotnew + eodd
         eoddnew = eodd

         call stop_clock( "inner_loop" )
         !
      endif
      !
      eodd = sum(pink(1:nbsp))
      
     end subroutine do_innerloop_subroutine

     subroutine print_out_observables()
     
#ifdef __DEBUG
        ! for debug and tuning purposes
        if ( ionode ) write(37,*)itercg, itercgeff, etotnew
        if ( ionode ) write(1037,'("iteration =",I4,"  eff iteration =",I4,"   Etot (Ha) =",F22.14)')&
            itercg, itercgeff, etotnew 
#endif
!         tcpu_cg_here=get_clock('nk_corr')
        call print_clock('CP')
!       #
        if ( ionode ) then
            if (itercg>2) then
               write(stdout,'(5x,"iteration =",I4,"  eff iteration =",I4,"   Etot (Ha) =",F22.14," delta_E=",E22.14)')&
               itercg, itercgeff, etotnew, deltae
            else
               write(stdout,'(5x,"iteration =",I4,"  eff iteration =",I4,"   Etot (Ha) =",F22.14)')&
               itercg, itercgeff, etotnew
            endif
            !
        endif
        !

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
     
     end subroutine print_out_observables

     subroutine check_convergence_cg()
        !
        deltae=abs(etotnew-etotold)
        !
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
     end subroutine check_convergence_cg

     subroutine compute_hpsi()
        ! 
        ! faux takes into account spin multiplicity.
        !
        faux(1:nbspx)=0.d0
        faux(1:nbsp) = max(f_cutoff,f(1:nbsp)) * DBLE( nspin ) / 2.0d0
        !
        !
        ! This condition to ensure that the orbital(fixed_band) is frozen
        ! 
        IF (fixed_state) THEN
           faux(fixed_band) = f_cutoff
        ENDIF
        !
        do i=1,nbsp,2
          ! 
          ! FIRST CALL TO DFORCE
          !
          CALL start_clock( 'dforce1' )
          !  
          CALL dforce( i, bec, betae, c0, c2, c3,rhos, nnrsx, ispin, faux, nbsp, nspin)
          !
          CALL stop_clock( 'dforce1' )
          !
          ! COMPUTE DFORCE FROM BERRY PHASE
          !
          if ( tefield .and. (evalue.ne.0.d0)) then
             !
             call dforceb(c0, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
             !
             c2(:)=c2(:) + evalue*df(:)
             !  
             call dforceb(c0, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
             !
             c3(:)=c3(:) + evalue*df(:)
             !
          endif

          IF ( lda_plus_u ) THEN
               !
               !
               c2(:) = c2(:) - vupsi(:,i) * faux(i)
               if( i+1 <= nbsp ) c3(:) = c3(:) - vupsi(:,i+1) * faux(i+1)
               !               
               !
           ENDIF

!$$
          if ( do_orbdep ) then
              !
              IF ( odd_nkscalfact ) THEN
                 !
                 c2(:) = c2(:) - valpsi(i,:) * faux(i)
                 !
                 if( i+1 <= nbsp ) c3(:) = c3(:) - valpsi(i+1,:) * faux(i+1)
                 ! 
              ENDIF
              !
              ! faux takes into account spin multiplicity.
              !
              CALL nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, bec, ngw, c0(:,i), c0(:,i+1), vsicpsi, lgam )
              !
              !
              c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
              !
              if( i+1 <= nbsp )   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
              !
              
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
!$$

#ifdef __TOBE_FIXED
          if (finite_field_introduced) then
             !
             ! faux takes into account spin multiplicity.
             !
             CALL finite_field_force ( i, nbsp, nbspx, dvpot, ngw, c0(:,i), c0(:,i+1), dvpotpsi, lgam )
             !
             c2(:) = c2(:) - dvpotpsi(:,1) * faux(i)
             !
             if( i+1 <= nbsp )   c3(:) = c3(:) - dvpotpsi(:,2) * faux(i+1)
             !
          endif
#endif

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
     
     end subroutine compute_hpsi

     subroutine orthogonalize(wfc0, wfc, becwfc, bec0)
        ! 
        type(twin_matrix) :: becwfc, bec0
        complex(DP) :: wfc(:,:), wfc0(:,:)
        ! 
        if (switch.or.(.not.do_orbdep)) then
           !
           if (fixed_state) then
              !
              wfc(:, fixed_band) = cmplx(0.d0, 0.d0)
              ! 
           endif
           !
           call calbec(1,nsp,eigr,wfc,becwfc)
           !  
           call pc2(wfc0,bec0,wfc,becwfc, lgam)
           !
        else
           !
           if (fixed_state) then
              !
              call pc3nc_fixed(wfc0, wfc, lgam)
              !
           endif
           !   
           call calbec(1,nsp,eigr,wfc,becwfc)
           ! 
           if (.not.okvan) then
              !
              call pc3nc(wfc0,wfc,lgam)
              !
           else
              ! 
              call pc3us(wfc0,bec0,wfc,becwfc, lgam)
              !
           endif
           !
        endif
        !  
     end subroutine orthogonalize
     !
     subroutine orthogonalize_wfc_only(wfc,becwfc)
        !
        type(twin_matrix) :: becwfc
        complex(DP) :: wfc(:,:)
        complex(DP) :: s(nbspx,nbspx)
        integer :: nbnd1,nbnd2,ndim,i,j,isp
        !
        call calbec(1,nsp,eigr,wfc,becwfc)
        !
        IF (do_orbdep.and.ortho_switch) THEN
           !            
           if(.not. okvan) then
              !
              call lowdin(wfc, lgam)
              !
           else
              !
              call lowdin_uspp(wfc,becwfc,lgam)
              !
           endif
           !
        ELSE
           !
           if (fixed_state) then
             !
             call gram_swap(betae,becwfc,nhsa,wfc,ngw,nbsp, fixed_band)
             ! 
           else
             ! 
             call gram(betae,becwfc,nhsa,wfc,ngw,nbsp)
             ! 
           endif 
            !
        ENDIF
        !
        call calbec(1,nsp,eigr,wfc,becwfc)
        !
        s(:,:)=CMPLX(0.d0,0.d0)
        !
        DO isp=1,nspin
           ndim=nupdwn(isp)
           DO i=1,ndim
              !
              nbnd1=iupdwn(isp)-1+i
              !
              DO j=1,i
                 !
                 nbnd2=iupdwn(isp)-1+j
                 !
                 call dotcsv( s(j,i), nbspx, nbsp, wfc, becwfc, wfc, becwfc, ngw, iupdwn(isp)+j-1, iupdwn(isp)+i-1, lgam)
                 s(i,j)=CONJG(s(j,i))
                 !
              ENDDO
              !
           ENDDO
        ENDDO
        !
     end subroutine orthogonalize_wfc_only
     !                     
     subroutine compute_lambda_bare ()
        ! 
        hitmp(:,:) = c0(:,:)
        !
        DO is = 1, nspin
           !
           nss = nupdwn(is)
           istart = iupdwn(is)
           !
           IF(.not.lambda(1)%iscmplx) THEN
              lambda_repl = 0.d0
           ELSE
              lambda_repl_c = CMPLX(0.d0,0.d0)
           ENDIF
           !
           DO i = 1, nss
              ! 
              DO j = i, nss
                 ii = i + istart - 1
                 jj = j + istart - 1
                 IF (.not.lambda(1)%iscmplx) THEN
                    !
                    DO ig = 1, ngw
                       lambda_repl( i, j ) = lambda_repl( i, j ) - &
                       2.d0 * DBLE( CONJG( hitmp( ig, ii ) ) * gi_bare( ig, jj) )
                    ENDDO
                    !
                    IF( ng0 == 2 ) THEN
                       lambda_repl( i, j ) = lambda_repl( i, j ) + &
                       DBLE( CONJG( hitmp( 1, ii ) ) * gi_bare( 1, jj ) )
                    ENDIF
                    ! 
                    lambda_repl( j, i ) = lambda_repl( i, j )
                    !
                 ELSE
                    !
                    DO ig = 1, ngw
                       lambda_repl_c( i, j ) = lambda_repl_c( i, j ) - &
                       CONJG( hitmp( ig, ii ) ) * gi_bare( ig, jj)
                    ENDDO
                    !
                    lambda_repl_c( j, i ) = CONJG(lambda_repl_c( i, j ))
                    !    
                 ENDIF
                 !
              ENDDO
              ! 
           ENDDO
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
        ENDDO
        ! 
     return
     !
     end subroutine compute_lambda_bare
     !
     subroutine makov_payne_correction_2nd_term (charge, quadrupole)
       !
       real (DP) :: charge, quadrupole
       real (DP) :: corr2, corr1 
       ! 
       ! 1 / 2 Ry -> a.u.
       corr1 = - 2.8373D0 / alat * charge**2 / 2.0D0
       !
       corr2 = ( 2.D0 / 3.D0 * pi )*( charge*quadrupole)/ omega
       !
       write(stdout, *) "Test MP:", charge, omega, quadrupole
       write(stdout, *) "Makov-Payne 1st energy ", corr1
       write(stdout, *) "Makov-Payne 2nd energy ", corr2
       !
       return  
       !
     end subroutine makov_payne_correction_2nd_term
     ! 
END SUBROUTINE runcg_uspp
