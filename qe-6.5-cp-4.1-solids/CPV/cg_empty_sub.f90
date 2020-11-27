
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
#define DEBUG_NUMDER
!
!=======================================================================
subroutine runcg_uspp_emp( c0_emp, cm_emp, bec_emp, f_emp, fsic_emp, n_empx,&
                          n_emps, ispin_emp, iupdwn_emp, nupdwn_emp, phi_emp, lambda_emp, &
                          maxiter_emp, wxd_emp, vsic_emp, sizvsic_emp, pink_emp, nnrx, rhovan_emp, &
                          deeq_sic_emp, nudx_emp, eodd_emp, etot_emp, &
                          filledstates_potential, nfi, tfirst, tlast, eigr, bec, irb, eigrb, &
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
                                           ispin
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
      ! use dener
      use cdvan
      use constants,                only : pi, au_gpa
      use io_files,                 only : psfile, pseudo_dir
      USE io_files,                 ONLY : outdir, prefix
      use uspp,                     only : nhsa=> nkb, nhsavb=> nkbus, betae => vkb, rhovan => becsum, deeq,qq
      use uspp_param,               only : nh
      use cg_module,                only : ene_ok,  maxiter,niter_cg_restart, &
                                           conv_thr, passop, enever, itercg
      use ions_positions,           only : tau0
      use wavefunctions_module,     only : c0 => cp
      use mp,                       only : mp_sum, mp_bcast
      use cp_electronic_mass,       ONLY : emass_cutoff
      use orthogonalize_base,       ONLY : calphi
      use cp_interfaces,            ONLY : rhoofr, dforce, compute_stress, nlfl, set_x_minus1, xminus1
      USE cp_main_variables,        ONLY : nlax, collect_lambda, distribute_lambda, nrlx, nlam
      USE descriptors,              ONLY : la_npc_ , la_npr_ , la_comm_ , la_me_ , la_nrl_ , ldim_cyclic
      USE mp_global,                ONLY : me_image,my_image_id
      !
      use twin_types !added:giovanni
      use control_flags,            only : non_ortho, iprint_spreads
      use printout_base,            only : printout_base_open, printout_base_unit, &
                                           printout_base_close
      use nksic,                    only : odd_alpha, valpsi, nkscalfact, do_orbdep, wtot, vsicpsi, sizwtot, & 
                                           do_innerloop_empty, do_innerloop_cg, innerloop_cg_nsd, &
                                           innerloop_cg_nreset, innerloop_init_n, innerloop_cg_ratio, &
                                           innerloop_until

      use cp_main_variables,        only : becdual, becmdual, overlap, ioverlap
      use electrons_module,         only : wfc_spreads_emp, wfc_centers_emp, icompute_spread
      use cp_interfaces,            only : gram_empty, nlsm1
      use uspp_param,               only : nhm
      use descriptors,              only : descla_siz_
      use input_parameters,         only : odd_nkscalfact_empty, wo_odd_in_empty_run, odd_nkscalfact
      !
      implicit none
      !
      integer     :: nfi
      logical     :: tfirst , tlast
      integer     :: sizvsic_emp
      complex(dp) :: eigr(ngw,nat)
      type(twin_matrix)    :: bec 
      type(twin_matrix )   :: bec_emp
      type(twin_matrix )   :: lambda_emp(nspin)
      integer     :: irb(3,nat)
      complex(dp) :: eigrb(ngb,nat)
      real(dp)    :: rhor(nnr,nspin)
      complex(dp) :: rhog(ngm,nspin)
      real(dp)    :: rhos(nnrsx,nspin)
      real(dp)    :: rhoc(nnr)
      real(dp)    :: ema0bg(ngw)
      integer     :: n_emps, n_empx, iupdwn_emp(nspin), nupdwn_emp(nspin), maxiter_emp, nnrx, &
                     nudx_emp, ispin_emp(n_empx)
      real(dp)    :: f_emp(n_empx), fsic_emp(n_empx), wxd_emp(sizvsic_emp,2), vsic_emp(sizvsic_emp, n_empx), &
                     pink_emp(n_empx), rhovan_emp(nhm*(nhm+1)/2, nat, nspin), &
                     deeq_sic_emp(nhm,nhm,nat,n_empx), eodd_emp, etot_emp, & 
                     filledstates_potential(nnrsx,nspin)
      complex(dp) :: c0_emp(ngw, n_empx), cm_emp(ngw, n_empx), phi_emp(ngw, n_empx)
      integer, intent(in)    :: desc_emp(descla_siz_, 2)
      !
      ! local variables
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
      type(twin_matrix) :: bec0, becm !modified:giovanni
      real(kind=DP), allocatable :: fmat_(:,:)!average kinetic energy for preconditioning
      complex(kind=DP), allocatable :: fmat_c_(:,:)!average kinetic energy for preconditioning
      ! 
      complex(DP) :: esse_c, essenew_c !factors in c.g.
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
      real(DP)    :: dtmp, temp
      real(dp)    :: tmppasso, ene_save(100), ene_save2(100)
      !
      logical     :: lgam, switch=.false., ortho_switch=.false., okvan, steepest=.false.
      complex(DP) :: phase
      integer     :: ierr, northo_flavor
      real(DP)    :: deltae, sic_coeff1, sic_coeff2 !coefficients which may change according to the flavour of SIC
      integer     :: me, iunit_manifold_overlap, iunit_spreads
      character(len=10) :: tcpu_cg_here
      real(DP):: ekin_emp, enl_emp, dekin_emp(6), denl_emp(3,3), epot_emp
      real(DP), allocatable :: rhor_emp(:,:), rhos_emp(:,:), rhoc_emp(:)
      complex(DP), allocatable :: rhog_emp(:,:)
      real(DP), allocatable    :: faux_emp(:)
      integer                  :: in_emp, issw
      COMPLEX(DP), PARAMETER   :: c_zero=CMPLX(0.d0,0.d0)
      !
      ! var for numerical derivatives
      REAl(DP):: etot_emp_tmp1, etot_emp_tmp2, etot_emp_tmp
      !
      call do_allocation_initialization()
      !
      ! Initializing clock for minimization
      !
      call start_clock('runcg_uspp')

      if( tfirst .and. ionode ) &
         write(stdout,"(/,a,/)") 'PERFORMING CONJUGATE GRADIENT MINIMIZATION OF EMPTY STATES'
      !         
      ! set tpa mass preconditioning
      !
      call emass_precond_tpa ( ema0bg, tpiba2, emass_cutoff )
      ! 
      call prefor(eigr, betae) 
      !
      ! orthonormalize c0_empty
      !
      !call orthogonalize_wfc_only(c0_emp, bec_emp)
      call nlsm1( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
      !
      ! recompute phi (the augmented wave function) from the new c0_empty
      !
      CALL calphi( c0_emp, SIZE(c0_emp,1), bec_emp, nhsa, betae, phi_emp, n_emps, lgam)
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
      allocate( hpsi(ngw, n_emps) )
      allocate( hpsi0(ngw, n_emps) )
      allocate( gi(ngw, n_emps), hi(ngw, n_emps) )
      !
      allocate(hitmp(ngw, n_emps))
      hitmp(:,:) = CMPLX(0.d0,0.d0)
      !
      gi(:,:)=CMPLX(0.d0,0.d0)
      hi(:,:)=CMPLX(0.d0,0.d0)
      ene_ok=.false.
      !
      !=======================================================================
      !                 begin of the main loop
      !=======================================================================
      !
      OUTER_LOOP: &
      do while ( ((itercg < maxiter_emp).and.(.not. ltresh)).or. & 
                 ((maxiter_emp==1).and.(itercg==1).and.(do_orbdep.and.(.not.wo_odd_in_empty_run))) )
        !
        call start_clock( "outer_loop" )
        ! 
        ENERGY_CHECK: &
        if(.not. ene_ok ) then
          !
          call nlsm1( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam ) 
          !  
          call rhoofr_cp_ortho_new &
                         ( n_empx, n_emps, nudx_emp, f_emp, ispin_emp, iupdwn_emp, &
                           nupdwn_emp, nspin, nfi, c0_emp, irb, eigrb, bec_emp, &
                           rhovan_emp, rhor_emp, rhog_emp, rhos_emp, enl_emp, denl_emp, &
                           ekin_emp, dekin_emp, tstress, 0)
          !
          etot_emp = enl_emp + ekin_emp                                 
          !
          call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp) 
          !
          etot_emp = etot_emp + epot_emp
          !
          if ( do_orbdep .and. (.not. wo_odd_in_empty_run)  ) then
             !
             if (odd_nkscalfact_empty) then
                !
                valpsi(:,:)  = (0.0_DP, 0.0_DP)
                odd_alpha(:) = 0.0_DP
                !
                CALL odd_alpha_routine(c0_emp, n_emps, n_empx, lgam, .true.)
                !
             else
                !
                ! here, we want to use only one value alpha for all empty
                ! states,
                ! that value alpha is defined from in input file. 
                ! This require to deactive the odd_nkscalfact here so that 
                ! it does not do odd_alpha in nksic_potential.
                !  
                odd_nkscalfact = .false.
                ! 
             endif
             !
             icompute_spread = .true.
             !
             call nksic_potential( n_emps, n_empx, c0_emp, fsic_emp, &
                                   bec_emp, rhovan_emp, deeq_sic_emp, &
                                   ispin_emp, iupdwn_emp, nupdwn_emp, rhor, rhoc, &
                                   wtot, sizwtot, vsic_emp, .false., pink_emp, nudx_emp, &
                                   wfc_centers_emp, wfc_spreads_emp, &
                                   icompute_spread, .true.)
             !
             ! Print spreads infor
             !
             WRITE( stdout, *) "sum spreads:1", sum(wfc_spreads_emp(1:nupdwn_emp(1), 1, 1)), &
                                                sum(wfc_spreads_emp(1:nupdwn_emp(2), 2, 1))
             WRITE( stdout, *) "sum spreads:2", sum(wfc_spreads_emp(1:nupdwn_emp(1), 1, 2)), &
                                                sum(wfc_spreads_emp(1:nupdwn_emp(2), 2, 2))       
             !
             IF ( (maxiter_emp==1).and.(itercg==1) ) THEN 
                write(stdout, *) "Localization of orbitals from PZS localization"
                do i = 1, nupdwn_emp(2)
                   write(stdout, *) i, wfc_spreads_emp(i, 2, 2), pink_emp (nupdwn_emp(1)+i)
                enddo
                !
             ENDIF 
             !
             do i = 1, n_emps
                !  
                ! Here wxd_emp <-> wtot that computed from nksic_potential of occupied states.
                ! wtot is scaled with nkscalfact constant, we thus need to rescaled it here with
                ! odd_alpha
                !
                if (odd_nkscalfact_empty) wxd_emp(:,:) = wxd_emp(:,:)*odd_alpha(i)/nkscalfact
                !
                vsic_emp(:,i) = vsic_emp(:,i) + wxd_emp(:, ispin_emp(i))
                !
             enddo
             !
             eodd_emp=sum(pink_emp(1:n_empx))
             ! 
             etot_emp = etot_emp + eodd_emp
             !
          endif
          ! 
          etotnew=etot_emp
          !
        else
          !
          etot_emp=enever
          !
          etotnew=etot_emp
          ene_ok=.false.
          !
        endif ENERGY_CHECK
        !
!!! NLN  check begin
        IF ( (maxiter_emp==1).and.(itercg==1) ) THEN 
           if( do_orbdep ) then
             !
             call do_innerloop_subroutine()
             ! 
           endif
           EXIT
        ENDIF
        !!! orgin
        !if( do_orbdep ) then
          !
        !  call do_innerloop_subroutine()
          ! 
        !endif
!!!! NLN check end
        ! 
        call print_out_observables()
        !
        ! here we store the etot in ene0, to keep track of the energy of the initial point
        !
        call check_convergence_cg()        
        !
        !call prefor(eigr, betae)!ATTENZIONE
        !
        call compute_hpsi()
        !
        ! HPSI IS ORTHOGONALIZED TO c0
        !
        ! comp. <beta|hpsi>
        !  
        if(switch.or.(.not.do_orbdep).or.(do_orbdep.and.wo_odd_in_empty_run) ) then
           !
           call pcdaga2_new(c0_emp, phi_emp, hpsi, n_emps, ispin_emp, lgam)
           !
        else
           !
           call pc3nc_new(c0_emp, hpsi, n_emps, ispin_emp, lgam)
           !
        endif
        !
        CALL nlsm1 ( n_emps, 1, nsp, eigr, hpsi, becm, 1, lgam )
        !
        do iss=1,nspin
           !
           in_emp = iupdwn_emp( iss )
           issw   = iupdwn( iss )
           !
           CALL gram_empty(.true. , eigr, betae, becm, bec, nhsa, &
                            hpsi( :, in_emp: ), c0( :, issw: ), ngw, &
                            nupdwn_emp(iss), nupdwn(iss), in_emp, issw)
           !
        enddo
        ! 
        call nlsm1( n_emps, 1, nsp, eigr, hpsi, becm, 1, lgam )  
        !
        ! TWO VECTORS INITIALIZED TO HPSI
        ! 
        hpsi0(1:ngw,:) = hpsi(1:ngw,:) 
        !
        gi(1:ngw,:)    = hpsi(1:ngw,:)
        !
	! COMPUTES ULTRASOFT-PRECONDITIONED HPSI,
        ! non kinetic-preconditioned, 
        ! is the subsequent reorthogonalization necessary 
        ! in the norm conserving case???: giovanni
        !
        call xminus1_twin_new(hpsi, n_emps, betae, dumm, becm, s_minus1,.false.)
        !
        call orthogonalize(c0_emp, hpsi, becm, bec_emp)
        !
        call xminus1_twin_new(gi, n_emps, betae, ema0bg, becm, k_minus1, .true.)
        !
        call orthogonalize(c0_emp, gi, becm, bec_emp)
        !    
        !  calculates gamma
        !
        gamma_c=CMPLX(0.d0,0.d0)
        !
        DO i=1, n_emps
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
              do i=1, n_emps
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
              do i=1, n_emps
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
        if (itercg==1 .or. mod(itercg, niter_cg_restart)==0 .or. restartcg) then
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
           hi(1:ngw, : )= gi(1:ngw, :) 
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
           hi(:,:) = gi(:, :) + gamma_c * hi(:, : )
           !   
        endif
        !
        ! note that hi is saved on gi, because we need it before projection on conduction states
        !     
        ! ... find minimum along direction hi:
        !
        ! project hi on conduction sub-space
        !
        call orthogonalize(c0_emp, hi, bec0, bec_emp)
        !    
        ! do quadratic minimization
        !             
        ! calculate derivative with respect to lambda along direction hi
        !
        dene0=0.d0
        ! 
        do i=1, n_emps
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
        dene0 = dene0 * 2.d0/nspin
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
        ! below is the debug part, that computes numerical derivative.
        ! Thank Nicola for providing it. 
        !
        if (.false.) then
           !
           etot_emp_tmp1 = 0.0
           etot_emp_tmp2 = 0.0
           etot_emp_tmp  = 0.0
           !
           tmppasso = 0.0
           do i=1, 30 
           tmppasso = tmppasso + 0.1
           !
           !if (i==1)  tmppasso=1.d-5
           !if (i==2)  tmppasso=-1.d-5
           !
           cm_emp(:,:) = c0_emp(:,:) + spasso * tmppasso * hi(:,:)
           !
           if (lgam.and.ng0 == 2)  cm_emp(1,:) = 0.5d0*(cm_emp(1,:) + CONJG(cm_emp(1,:)))
           !
           ! orthonormalize
           !
           call orthogonalize_wfc_only(cm_emp, becm)
           !
           call rhoofr_cp_ortho_new &
                          ( n_empx, n_emps, nudx_emp, f_emp, ispin_emp, iupdwn_emp, &
                            nupdwn_emp, nspin, nfi, cm_emp, irb, eigrb, becm, &
                            rhovan_emp, rhor_emp, rhog_emp, rhos_emp, enl_emp, denl_emp, &
                            ekin_emp, dekin_emp, tstress, 0)
           !
           etot_emp_tmp =  enl_emp + ekin_emp
           !
           call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp)
           !
           etot_emp_tmp =  etot_emp_tmp + epot_emp
           !
           if ( do_orbdep .and. (.not. wo_odd_in_empty_run)  ) then
              !
              if (odd_nkscalfact_empty) then
                 !
                 valpsi(:,:)  = (0.0_DP, 0.0_DP)
                 odd_alpha(:) = 0.0_DP
                 !
                 CALL odd_alpha_routine(c0_emp, n_emps, n_empx, lgam, .true.)
                 !
              else
                 !
                 ! here, we want to use only one value alpha for all empty
                 ! states,
                 ! that value alpha is defined from in input file. 
                 ! This require to deactive the odd_nkscalfact here so that 
                 ! it does not do odd_alpha in nksic_potential.
                 !  
                 odd_nkscalfact = .false.
                 ! 
              endif
              !
              call nksic_potential( n_emps, n_empx, cm_emp, fsic_emp, &
                                    becm, rhovan_emp, deeq_sic_emp, &
                                    ispin_emp, iupdwn_emp, nupdwn_emp, rhor, rhoc, &
                                    wtot, sizwtot, vsic_emp, .false., pink_emp, nudx_emp, &
                                    wfc_centers_emp, wfc_spreads_emp, &
                                    icompute_spread, .true.)
             !
             eodd_emp=sum(pink_emp(:))
             !  
             etot_emp_tmp = etot_emp_tmp + eodd_emp
             !
          endif
          !
          write(stdout,*) "etot_emp: ", i, tmppasso, " = ", etot_emp_tmp
          !
          !if (i==1) etot_emp_tmp1 = etot_emp_tmp
          !if (i==2) etot_emp_tmp2 = etot_emp_tmp
          !  
          enddo
          !
          !write(stdout,*) "here is numerical derivative vs analytic derivative at step", itercg
          !write(stdout,*) "(etot_emp_tmp1-etot_emp_tmp2)/tmppasso, dene0, tmppasso, ((etot_emp-ene0)/tmppasso)/dene0" 
          !write(stdout,'(2e25.15,4e20.10)') (etot_emp_tmp1 - etot_emp_tmp2)/(2.0*tmppasso), dene0, tmppasso, ((etot_emp_tmp1 - etot_emp_tmp2)/(2.0*tmppasso)/dene0)
          !
        endif  
        !
        ! calculates wave-functions on a point on direction hi
        !
        cm_emp(:,:) = c0_emp(:,:) + spasso * passof * hi(:,:)
        !
        if (lgam.and.ng0 == 2)  cm_emp(1,:) = 0.5d0*(cm_emp(1,:) + CONJG(cm_emp(1,:)))
        ! 
        ! orthonormalize
        !
        call orthogonalize_wfc_only(cm_emp, becm)
        ! 
        call rhoofr_cp_ortho_new &
                         ( n_empx, n_emps, nudx_emp, f_emp, ispin_emp, iupdwn_emp, &
                           nupdwn_emp, nspin, nfi, cm_emp, irb, eigrb, becm, &
                           rhovan_emp, rhor_emp, rhog_emp, rhos_emp, enl_emp, denl_emp, &
                           ekin_emp, dekin_emp, tstress, 0)        
        !
        etot_emp =  enl_emp + ekin_emp 
        !
        call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp)
        !
        etot_emp =  etot_emp + epot_emp
        !
        if ( do_orbdep .and. (.not. wo_odd_in_empty_run)  ) then
           !
           if (odd_nkscalfact_empty) then
              !
              valpsi(:,:)  = (0.0_DP, 0.0_DP)
              odd_alpha(:) = 0.0_DP
              !
              CALL odd_alpha_routine(c0_emp, n_emps, n_empx, lgam, .true.)
              !
           else
              !
              ! here, we want to use only one value alpha for all empty
              ! states,
              ! that value alpha is defined from in input file. 
              ! This require to deactive the odd_nkscalfact here so that 
              ! it does not do odd_alpha in nksic_potential.
              !  
              odd_nkscalfact = .false.
              ! 
           endif
           !
           call nksic_potential( n_emps, n_empx, cm_emp, fsic_emp, &
                                      becm, rhovan_emp, deeq_sic_emp, &
                                      ispin_emp, iupdwn_emp, nupdwn_emp, rhor, rhoc, &
                                      wtot, sizwtot, vsic_emp, .false., pink_emp, nudx_emp, &
                                      wfc_centers_emp, wfc_spreads_emp, &
                                      icompute_spread, .true.)
           !
           eodd_emp=sum(pink_emp(:))
           !  
           etot_emp = etot_emp + eodd_emp
           !
        endif
        !  
        ene1=etot_emp
        !    
        ! find the minimum
        !
        call minparabola(ene0, spasso*dene0, ene1, passof, passo, enesti)
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
        cm_emp(:,:) = c0_emp(:,:) + spasso * passo * hi(:,:)
        !
        IF (lgam.and. ng0 == 2 ) cm_emp(1,:) = 0.5d0*(cm_emp(1,:)+CONJG(cm_emp(1,:)))
        !
        call orthogonalize_wfc_only(cm_emp, becm)
        !
        call rhoofr_cp_ortho_new &
                         ( n_empx, n_emps, nudx_emp, f_emp, ispin_emp, iupdwn_emp, &
                           nupdwn_emp, nspin, nfi, cm_emp, irb, eigrb, becm, &
                           rhovan_emp, rhor_emp, rhog_emp, rhos_emp, enl_emp, denl_emp, &
                           ekin_emp, dekin_emp, tstress, 0)
        !
        etot_emp =  enl_emp + ekin_emp
        !
        call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp)
        !
        etot_emp =  etot_emp + epot_emp
        !
        if ( do_orbdep .and. (.not. wo_odd_in_empty_run)  ) then
           !
           if (odd_nkscalfact_empty) then
              !
              valpsi(:,:)  = (0.0_DP, 0.0_DP)
              odd_alpha(:) = 0.0_DP
              !
              CALL odd_alpha_routine(cm_emp, n_emps, n_empx, lgam, .true.)
              !
           else
              !
              ! here, we want to use only one value alpha for all empty
              ! states,
              ! that value alpha is defined from in input file. 
              ! This require to deactive the odd_nkscalfact here so that 
              ! it does not do odd_alpha in nksic_potential.
              !  
              odd_nkscalfact = .false.
              ! 
           endif
           !
           icompute_spread = .true.
           call nksic_potential( n_emps, n_empx, cm_emp, fsic_emp, &
                                      becm, rhovan_emp, deeq_sic_emp, &
                                      ispin_emp, iupdwn_emp, nupdwn_emp, rhor, rhoc, &
                                      wtot, sizwtot, vsic_emp, .false., pink_emp, nudx_emp, &
                                      wfc_centers_emp, wfc_spreads_emp, &
                                      icompute_spread, .true.)
           !
           ! Print spreads infor
           !
           WRITE( stdout, *) "sum spreads:1", sum(wfc_spreads_emp(1:nupdwn_emp(1), 1, 1)), & 
                                              sum(wfc_spreads_emp(1:nupdwn_emp(2), 2, 1))
           WRITE( stdout, *) "sum spreads:2", sum(wfc_spreads_emp(1:nupdwn_emp(1), 1, 2)), &
                                              sum(wfc_spreads_emp(1:nupdwn_emp(2), 2, 2))
           !
           do i = 1, n_emps 
              !  
              ! Here wxd_emp <-> wtot that computed from nksic_potential of
              ! occupied states.
              ! wtot is scaled with nkscalfact constant, we thus need to
              ! rescaled it here with
              ! odd_alpha
              !
              if (odd_nkscalfact_empty) wxd_emp(:,:) = wxd_emp(:,:)*odd_alpha(i)/nkscalfact
              !
              vsic_emp(:,i) = vsic_emp(:,i) + wxd_emp(:, ispin_emp(i))
              !
           enddo
           !
           eodd_emp=sum(pink_emp(:))
           !  
           etot_emp = etot_emp + eodd_emp
           !
        endif
        !
        enever=etot_emp
        !
        ! check with  what supposed
        !
        write(stdout,*) 'ene0, dene0, ene1, enesti,enever, passo, passov, passof'
        write(stdout,"(7f18.12)") ene0, dene0, ene1, enesti,enever, passo, passov, passof
        !
        ! if the energy has diminished with respect to ene0 and ene1 , everything ok
        !
        if (((enever.lt.ene0) .and. (enever.lt.ene1))) then
           !
           c0_emp(:,:) = cm_emp(:,:)
           !
           call copy_twin(bec_emp, becm) !modified:giovanni
           !
           ene_ok=.true.
           !
           ! if  ene1 << energy < ene0; go to  ene1
           !
        elseif((enever.ge.ene1) .and. (enever.lt.ene0)) then
           !
           if (ionode) then
              ! 
              write(stdout,"(2x,a,i5,f20.12)") 'cg_sub: missed minimum, case 1, iteration',itercg, passof
              ! 
           endif
           ! 
           c0_emp(:,:) = c0_emp(:,:) + spasso*passov*hi(:,:)
           ! 
           passof = 2.d0*passov
           !
           restartcg = .true.
           !
           call orthogonalize_wfc_only(c0_emp, bec_emp)
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
           c0_emp(:,:) = c0_emp(:,:) + spasso * passov * hi(:,:)
           !
           passof = 1.d0*passov
           !
           restartcg = .true.
           !
           call orthogonalize_wfc_only(c0_emp, bec_emp)
           !
           ene_ok=.false.
           !  
           ! if ene > ene0, en1 do a steepest descent step
           !
        elseif((enever.ge.ene0).and.(ene0.le.ene1)) then
           !
           if(ionode) then
             write(stdout,"(2x,a,i5)") 'cg_sub: missed minimum, case 3, iteration, doing steepest descent',itercg
           endif
           !
           iter3=0
           !  
           do while(enever.ge.ene0 .and. iter3.lt.maxiter3)
             ! 
             iter3 = iter3 + 1
             !  
             passov = passov*0.5d0
             !
             cm_emp(:,:) = c0_emp(:,:) + spasso*passov*hi(:,:)
             !
             passof = 1.d0*passov
             ! 
             itercgeff = itercgeff+1
             !
             ! change the searching direction
             !
             spasso=spasso*(-1.d0)
             !
             call orthogonalize_wfc_only(cm_emp, becm)
             !
             call rhoofr_cp_ortho_new &
                         ( n_empx, n_emps, nudx_emp, f_emp, ispin_emp, iupdwn_emp, &
                           nupdwn_emp, nspin, nfi, cm_emp, irb, eigrb, becm, &
                           rhovan_emp, rhor_emp, rhog_emp, rhos_emp, enl_emp, denl_emp, &
                           ekin_emp, dekin_emp, tstress, 0)
             !
             etot_emp = enl_emp + ekin_emp
             !
             call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp) 
             ! 
             etot_emp = etot_emp + epot_emp 
             !
             if ( do_orbdep .and. (.not. wo_odd_in_empty_run)  ) then
                !
                if (odd_nkscalfact_empty) then
                   !
                   valpsi(:,:)  = (0.0_DP, 0.0_DP)
                   odd_alpha(:) = 0.0_DP
                   !
                   call odd_alpha_routine(cm_emp, n_emps, n_empx, lgam, .true.)
                   !
                else
                   !
                   ! here, we want to use only one value alpha for all empty
                   ! states,
                   ! that value alpha is defined from in input file. 
                   ! This require to deactive the odd_nkscalfact here so that 
                   ! it does not do odd_alpha in nksic_potential.
                   !  
                   odd_nkscalfact = .false.
                   ! 
                endif
                !
                icompute_spread = .false.
                call nksic_potential( n_emps, n_empx, cm_emp, fsic_emp, &
                                      becm, rhovan_emp, deeq_sic_emp, &
                                      ispin_emp, iupdwn_emp, nupdwn_emp, rhor, rhoc, &
                                      wtot, sizwtot, vsic_emp, .false., pink_emp, nudx_emp, &
                                      wfc_centers_emp, wfc_spreads_emp, &
                                      icompute_spread, .true.)
                !
                eodd_emp = sum(pink_emp(:))
                !
                etot_emp = etot_emp + eodd_emp
                !
             endif
             ! 
             enever=etot_emp
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
              c0_emp(:,:)=cm_emp(:,:)
              ! 
              call copy_twin(bec_emp, becm)
              !
           endif
           !
           restartcg=.true.
           ene_ok=.false.
           !
           if (iter3 == maxiter3) then
              !
              passof=passop
              !
           endif
           !
        endif
        !  
        if (.not.ene_ok) call nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
        !
        ! calculates phi for pc_daga
        !
        call calphi( c0_emp, SIZE(c0_emp,1), bec_emp, nhsa, betae, phi_emp, n_emps, lgam )
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
      !=======================================================================
      !                 end of the main loop
      !=======================================================================
      !
      ! faux takes into account spin multiplicity.
      !
      faux(:) = f_emp(:) * DBLE( nspin ) / 2.0d0
      !
      do i = 1, n_emps, 2
         !
         call start_clock( 'dforce2' )
         !
         call dforce(i, bec_emp, betae, c0_emp, c2, c3, filledstates_potential, nnrsx, ispin_emp, faux, n_emps, nspin)
         !
         call start_clock( 'dforce2' )
         ! 
         if ( do_orbdep .and. (.not. wo_odd_in_empty_run) ) then
            !
            ! faux takes into account spin multiplicity.
            !
            call nksic_eforce( i, n_emps, n_empx, vsic_emp, deeq_sic_emp, bec_emp, ngw, c0_emp(:,i), c0_emp(:,i+1), vsicpsi, lgam )
            !
            c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
            !
            if( i+1 <= n_emps)   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
            !
         endif
         !
         do ig=1, ngw
            !
            gi(ig, i)=c2(ig)
            !
            if(i+1 <= n_emps) gi(ig,i+1)=c3(ig)
            ! 
         enddo
         !
         if (lgam.and.ng0.eq.2) then
            ! 
            gi(1,  i)=CMPLX(DBLE(gi(1,  i)),0.d0)
            !
            if(i+1 <= n_emps) gi(1,i+1)=CMPLX(DBLE(gi(1,i+1)),0.d0)
            !    
         endif
         !
      enddo
      !
      if(.not.lambda_emp(1)%iscmplx) then
         allocate(lambda_repl(nudx_emp,nudx_emp))
      else
         allocate(lambda_repl_c(nudx_emp,nudx_emp))
      endif
      !
      hitmp(:,:) = c0_emp(:,:)
      !
      do is = 1, nspin
         !
         nss = nupdwn_emp(is)
         istart = iupdwn_emp(is)
         ! 
         if (.not.lambda_emp(1)%iscmplx) then
            lambda_repl = 0.d0
         else
            lambda_repl_c = CMPLX(0.d0,0.d0)
         endif
         !
         do i = 1, nss
            ! 
            do j = i, nss
               !
               ii = i + istart - 1
               jj = j + istart - 1
               !
               if (.not.lambda_emp(1)%iscmplx) then
                  !
                  do ig = 1, ngw
                     !
                     lambda_repl( i, j ) = lambda_repl( i, j ) - &
                     2.d0 * DBLE( CONJG( hitmp( ig, ii ) ) * gi( ig, jj) )
                     !
                  enddo
                  !
                  if ( ng0 == 2 ) then
                     !  
                     lambda_repl( i, j ) = lambda_repl( i, j ) + &
                    DBLE( CONJG( hitmp( 1, ii ) ) * gi( 1, jj ) )
                    !
                 endif
                 !
                 lambda_repl( j, i ) = lambda_repl( i, j )
                 !
              else
                 !
                 do ig = 1, ngw
                    !
                    lambda_repl_c( i, j ) = lambda_repl_c( i, j ) - &
                    CONJG( hitmp( ig, ii ) ) * gi( ig, jj)
                    !
                 enddo
                 !  
                 lambda_repl_c( j, i ) = CONJG(lambda_repl_c( i, j ))
                 !
              endif
              !
           enddo
           !
        enddo
        !
        if (.not.lambda_emp(1)%iscmplx) then
           !  
           call mp_sum( lambda_repl, intra_image_comm )
           call distribute_lambda( lambda_repl, lambda_emp(is)%rvec( :, :),  desc_emp( :, is ) )
           !
        else
           !
           call mp_sum( lambda_repl_c, intra_image_comm )
           call distribute_lambda( lambda_repl_c, lambda_emp(is)%cvec( :, :), desc_emp( :, is ) )
           ! 
        endif
        !
     enddo
     !
     if (.not.lambda_emp(1)%iscmplx) then
        ! 
        deallocate( lambda_repl )
        !
     else
        !  
        deallocate( lambda_repl_c )
        !
     endif
     !
     call do_deallocation()
     !
     return
     !
     contains
     !
     subroutine do_allocation_initialization()
         !  
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
         allocate (faux(n_empx))
         allocate (faux_emp(n_empx))
         !
         faux_emp=0.d0
         !
         allocate (c2(ngw),c3(ngw))
         !
         allocate(rhor_emp(nnr,nspin), rhos_emp(nnrsx,nspin), rhog_emp(ngm,nspin))
         !
         if (nlcc_any) allocate(rhoc_emp(nnr))
         !  
         ! Allocation of twin-type variables
         !
         call init_twin(bec0, lgam)
         call allocate_twin(bec0, nhsa, n_emps, lgam)
         !
         call init_twin(becm, lgam)
         call allocate_twin(becm, nhsa, n_emps, lgam)
         !
         ! initializing variables for checking iterations and convergence
         !
         newscheme=.false.
         firstiter=.true.
         maxiter3=12
         !
         if(do_orbdep) maxiter3=10
         !
         ninner=0
         !
         ltresh    = .false.
         itercg    = 1
         etotold   = 1.d8
         restartcg = .true.
         passof = passop
         ene_ok = .false.
         !
         itercgeff = 1
         !
     end subroutine do_allocation_initialization
     
     subroutine do_deallocation()
        !
        deallocate(hpsi0,hpsi,gi,hi)
        deallocate(hitmp, STAT=ierr)
        !        
        call deallocate_twin(s_minus1)
        call deallocate_twin(k_minus1)
        !
        call stop_clock('runcg_uspp')
        !
        call deallocate_twin(bec0)
        call deallocate_twin(becm)
        !
        deallocate(c2,c3)
        deallocate(faux)
        deallocate(faux_emp)
        deallocate(rhor_emp,rhos_emp,rhog_emp)
        !
        if (nlcc_any) deallocate(rhoc_emp)
        !
        IF(ionode) THEN
           ! 
           CALL printout_base_close( "sha" )
           CALL printout_base_close( "ovp" )
           !
        ENDIF
        ! 
     end subroutine do_deallocation
     ! 
     subroutine do_innerloop_subroutine()
        !
	if (do_innerloop_empty .and. innerloop_until>=itercgeff) then
           !
           call start_clock( "inner_loop" )
           !
           eodd_emp= sum(pink_emp(:))
           etot_emp= etot_emp - eodd_emp
           etotnew = etotnew  - eodd_emp
           ninner  = 0
           !
           if (.not.do_innerloop_cg) then
              ! 
              write(stdout,*)  "WARNING, do_innerloop_cg should be .true."
              ! 
           else
              !
              call nksic_rot_emin_cg_general(itercg,innerloop_init_n,ninner,etot_emp,deltae*innerloop_cg_ratio,lgam, &
                                     n_emps, n_empx, nudx_emp, iupdwn_emp, nupdwn_emp, ispin_emp, c0_emp, rhovan_emp, bec_emp, rhor, rhoc, &
                                     vsic_emp, pink_emp, deeq_sic_emp, wtot, fsic_emp, sizwtot, .false.,  wfc_centers_emp, wfc_spreads_emp, .true.) 
              !
           endif
           !
           eodd_emp= sum(pink_emp(:)) 
           etot_emp= etot_emp + eodd_emp 
           etotnew = etotnew  + eodd_emp
           ! 
           call stop_clock( "inner_loop" )
           !
        endif
        !
     endsubroutine do_innerloop_subroutine
     !   
     subroutine print_out_observables()
        ! 
        call print_clock('CP')
        !
        if ( ionode ) then
           !
           if (itercg>2) then
              write(stdout,'(5x,"iteration =",I4,"  eff iteration =",I4,"   Etot (Ha) =",F22.14," delta_E=",E22.14)')&
              itercg, itercgeff, etotnew, deltae
           else
              write(stdout,'(5x,"iteration =",I4,"  eff iteration =",I4,"   Etot (Ha) =",F22.14)')&
              itercg, itercgeff, etotnew
           endif
           !
           write(stdout,'(5x,  "Ekin (Ha) = ",F22.14 , " Enl (Ha) = ",F22.14, " Eloc (Ha) =" , F22.14)' )&
              ekin_emp, enl_emp, epot_emp 
           !
           if ( do_orbdep .and. (.not. wo_odd_in_empty_run)  ) then 
              write(stdout,'(1x,  "Fake EODD (Ha) = ",F22.14) ') eodd_emp
           endif
           !
        endif
        !
        if ( ionode .and. mod(itercg,10) == 0 ) write(stdout,"()" )
        !
        !if ( ionode .and. mod(itercg, iprint_spreads)==0) then
        if ( .false.) then
           !
           if(nspin==1) then
             !
             write( iunit_spreads, '(400f20.14)') wfc_spreads_emp(:,1,2)               
             !
           elseif(nspin==2) then
             !
             write( iunit_spreads, '(2(400f20.14)(3x))') wfc_spreads_emp(:,1,2), wfc_spreads_emp(:,2,2)
             !
           endif
           !
        endif
        !
     end subroutine print_out_observables
     !
     subroutine check_convergence_cg()
        !
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
     end subroutine check_convergence_cg
     !
     subroutine compute_hpsi()
        ! 
        ! faux takes into account spin multiplicity.
        !
        faux(:)=0.d0
        faux(1:n_emps) = f_emp(1:n_emps) * DBLE( nspin ) / 2.0d0
        !
        do i=1, n_emps ,2
           ! 
           call dforce( i, bec_emp, betae, c0_emp, c2, c3, filledstates_potential, nnrsx, ispin_emp, faux, n_emps, nspin) 
           !
           ! ODD terms
           !
           if ( do_orbdep .and. (.not. wo_odd_in_empty_run) ) then
              !   
              CALL nksic_eforce( i, n_emps, n_empx, vsic_emp, deeq_sic_emp, bec_emp, ngw, &
                                       c0_emp(:,i), c0_emp(:,i+1), vsicpsi, lgam )
              !
              c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
              !
              if( i+1 <= n_emps )   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
              !
           endif
           !
           hpsi(1:ngw, i)=c2(1:ngw)
           !
           if (i+1 <= n_emps ) then
              hpsi(1:ngw, i+1)=c3(1:ngw)
           endif
           !
           if (lgam) then
              !
              if (ng0.eq.2) then
                 ! 
                 hpsi(1, i)=CMPLX(DBLE(hpsi(1, i)), 0.d0)
                 !
                 if (i+1 <= n_emps) then
                    !
                    hpsi(1,i+1)=CMPLX(DBLE(hpsi(1,i+1)), 0.d0)
                    ! 
                 endif
                 !
              endif
              !
           endif
           !
        enddo
        !  
     end subroutine compute_hpsi

     subroutine orthogonalize(wfc0, wfc, becwfc, bec0)
        ! 
        type(twin_matrix) :: becwfc, bec0
        complex(DP) :: wfc(:,:), wfc0(:,:)
        !
        if (switch.or.(.not.do_orbdep).or.(do_orbdep .and.wo_odd_in_empty_run)) then
           !
           call pc2_new(wfc0, bec0, wfc, becwfc, n_emps, &
                        nupdwn_emp, iupdwn_emp, ispin_emp, lgam)
           !
        else
           !
           if (.not.okvan) then
              !
              call pc3nc_new(wfc0, wfc, n_empx, ispin_emp, lgam)
              !
           else
              ! 
              call pc3nc_new(wfc0, wfc, n_empx, ispin_emp, lgam)
              !
           endif
           !
        endif
        ! 
        CALL nlsm1 ( n_emps, 1, nsp, eigr, wfc, becwfc, 1, lgam )
        !
        do iss=1,nspin
           !
           in_emp = iupdwn_emp( iss )
           issw   = iupdwn( iss )
           !
           CALL gram_empty(.true., eigr, betae, becwfc, bec, nhsa, &
                            wfc( :, in_emp: ), c0( :, issw: ), &
                            ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, issw)
           !
        enddo
        !  
     end subroutine orthogonalize
     !
     subroutine orthogonalize_wfc_only(wfc,becwfc)
        !
        type(twin_matrix) :: becwfc
        complex(DP) :: wfc(:,:)
        !
        call nlsm1 ( n_emps, 1, nsp, eigr, wfc, becwfc, 1, lgam )
        !
        do iss=1,nspin
           !
           issw   = iupdwn(iss)
           in_emp = iupdwn_emp(iss)
           !
           CALL gram_empty(.false., eigr, betae, becwfc, bec, nhsa, &
                            wfc( :, in_emp: ), c0( :, issw: ), &
                            ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, issw)           
           !
        enddo
        ! 
        call nlsm1 ( n_emps, 1, nsp, eigr, wfc, becwfc, 1, lgam ) 
        !
     end subroutine orthogonalize_wfc_only
     !
     subroutine v_times_rho(v, nspin, rhos_emp, epot_emp)
        !
        use kinds, only: DP
        use mp,                   only : mp_sum
        use mp_global,            only : intra_image_comm
        use cell_base,            only : omega
        use smooth_grid_dimensions,      only : nnrsx, nr1s, nr2s, nr3s
        !
        implicit none
        ! 
        integer, intent(in)   :: nspin
        real(DP), intent(in)  :: v(nnrsx,nspin), rhos_emp(nnrsx,nspin)
        real(DP), intent(out) :: epot_emp
        !
        ! local vars
        !   
        integer  :: i
        real(DP) :: etemp, fact, rhosum(2)
        ! 
        etemp=0.d0
        rhosum=0.d0
        fact=omega/DBLE(nr1s*nr2s*nr3s)
        ! 
        do i=1,nspin
           !
           etemp = etemp + sum(v(1:nnrsx,i) * rhos_emp(1:nnrsx,i))
           ! 
           rhosum(i) =  sum(rhos_emp(1:nnrsx,i))
           !
        enddo
        !  
        call mp_sum(etemp, intra_image_comm)
        !  
        call mp_sum(rhosum, intra_image_comm)
        !
        epot_emp = etemp*fact
        !
        rhosum = rhosum*fact
        !
        return
        !
     end subroutine v_times_rho
     !                     
END SUBROUTINE runcg_uspp_emp
