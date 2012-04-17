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
                          fion, ema0bg, becdr, lambdap, lambda, vpot  )
!=======================================================================

      use kinds,                    only : dp
      use control_flags,            only : iprint, thdyn, tpre, iprsta, &
                                           tfor, taurdr, tprnfor, gamma_only, do_wf_cmplx !added:giovanni gamma_only, do_wf_cmplx
      use control_flags,            only : ndr, ndw, nbeg, nomore, tsde, tortho, tnosee, &
                                           tnosep, trane, tranp, tsdp, tcp, tcap, ampre, &
                                           amprp, tnoseh
      use core,                     only : nlcc_any
      !---ensemble-DFT
      use energies,                 only : eht, epseu, exc, etot, eself, enl, ekin,&
                                           atot, entropy, egrand
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
      use wavefunctions_module,     only : c0, cm, phi => cp
      use efield_module,            only : tefield, evalue, ctable, qmat, detq, ipolp, &
                                           berry_energy, ctabin, gqq, gqqm, df, pberryel, &
                                           tefield2, evalue2, ctable2, qmat2, detq2, ipolp2, &
                                           berry_energy2, ctabin2, gqq2, gqqm2, pberryel2
      use mp,                       only : mp_sum, mp_bcast
      use cp_electronic_mass,       ONLY : emass_cutoff
      use orthogonalize_base,       ONLY : calphi
      use cp_interfaces,            ONLY : rhoofr, dforce, compute_stress, nlfl
      USE cp_main_variables,        ONLY : nlax, collect_lambda, distribute_lambda, descla, nrlx, nlam
      USE descriptors,              ONLY : la_npc_ , la_npr_ , la_comm_ , la_me_ , la_nrl_ , ldim_cyclic
      USE mp_global,                ONLY : me_image,my_image_id
      !
      use nksic,                    only : do_orbdep, do_innerloop, do_innerloop_cg, innerloop_cg_nsd, &
                                           innerloop_cg_nreset, &
                                           vsicpsi, vsic, wtot, fsic, fion_sic, deeq_sic, f_cutoff, pink
      use hfmod,                    only : do_hf, vxxpsi, exx
      use twin_types !added:giovanni

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
      type(twin_matrix), dimension(nspin)    :: lambdap!(nlam,nlam,nspin) !modified:giovanni
      type(twin_matrix), dimension(nspin)    :: lambda!(nlam,nlam,nspin)   !modified:giovanni
!
!
      integer     :: i, j, ig, k, is, iss,ia, iv, jv, il, ii, jj, kk, ip
      integer     :: inl, jnl, niter, istart, nss, nrl, me_rot, np_rot , comm
      real(dp)    :: enb, enbi, x
      real(dp)    :: gamma, entmp, sta
      complex(dp), allocatable :: c2(:)
      complex(dp), allocatable :: c3(:)
      complex(dp), allocatable :: hpsi(:,:), hpsi0(:,:), gi(:,:), hi(:,:)
      real(DP),    allocatable :: s_minus1(:,:)    !factors for inverting US S matrix
      real(DP),    allocatable :: k_minus1(:,:)    !factors for inverting US preconditioning matrix 
!       real(DP),    allocatable :: lambda_repl(:,:) ! replicated copy of lambda
!       real(DP),    allocatable :: lambda_dist(:,:) ! replicated copy of lambda
      type(twin_matrix) :: lambda_repl!(:,:) ! replicated copy of lambda !modified:giovanni
      type(twin_matrix) :: lambda_dist!(:,:) ! replicated copy of lambda  !modified:giovanni
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
      ! 
      logical     :: pre_state!if .true. does preconditioning state by state
      !
      real(DP)    :: esse,essenew !factors in c.g.
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
      real(DP)    :: esic
      real(DP)    :: Omattot(nbspx,nbspx)
      complex(DP) :: hi_tmp(ngw,nbsp)
      real(DP)    :: dtmp
      real(dp),    allocatable :: vsicah(:,:)
      real(dp)    :: vsicah2sum
      real(dp)    :: tmppasso, ene_save(100), ene_save2(100), ene_lda
      !
      logical :: lgam
      integer :: ispin
      !
      lgam=gamma_only.and..not.do_wf_cmplx
      !
      allocate (faux(nbspx))
      !allocate(hpsinorm(n))
      !allocate(hpsinosicnorm(n))
      !
!       allocate (bec0(nhsa,nbsp),becm(nhsa,nbsp), becdrdiag(nhsa,nspin*nlax,3))
      allocate (ave_ene(nbsp))
      allocate (c2(ngw),c3(ngw))

      !begin_added:giovanni
      call init_twin(bec0, lgam)
      call allocate_twin(bec0, nhsa, nbsp, lgam)
      call init_twin(becm, lgam)
      call allocate_twin(becm, nhsa, nbsp, lgam)
      call init_twin(becdrdiag, lgam)
      call allocate_twin(becdrdiag, nhsa, nspin*nlax,3, lgam)
      do ispin=1,nspin
	call init_twin(lambda, lgam)
	call allocate_twin(lambda, nlam, nlam, lgam)
	call init_twin(lambdap, lgam)
	call allocate_twin(lambdap, nlam, nlam, lgam)
      enddo
      !end_added:giovanni

      call start_clock('runcg_uspp')
      newscheme=.false.
      firstiter=.true.

      pre_state=.false.!normally is disabled

      maxiter3=250
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
      restartcg = .true.
      passof = passop
      ene_ok = .false.

!$$
      itercgeff = 1
!$$

      !orthonormalize c0
      call calbec(1,nsp,eigr,c0,bec)

      call gram(betae,bec,nhsa,c0,ngw,nbsp)

      !call calbec(1,nsp,eigr,c0,bec)

      ! calculates phi for pcdaga

      !call calphiid(c0,bec,betae,phi)
      CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, nbsp, lgam)
      !
      ! calculates the factors for S and K inversion in US case
      !
      if ( nvb > 0 ) then
          !
          allocate( s_minus1(nhsavb,nhsavb))
          allocate( k_minus1(nhsavb,nhsavb))
          call  set_x_minus1(betae,s_minus1,dumm,.false.)
          call  set_x_minus1(betae,k_minus1,ema0bg,.true.)
          !
      else
          !
          allocate( s_minus1(1,1))
          allocate( k_minus1(1,1))
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
      hitmp(:,:) = (0.d0,0.d0)
      ! allocate(hpsinosic(ngw,n))
      !
      gi(:,:)=(0.d0,0.d0)
      hi(:,:)=(0.d0,0.d0)


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
          if(.not.tens) then
             call rhoofr(nfi,c0(:,:),irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
          else

            if(newscheme.or.firstiter) then 
               call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,c0,bec,firstiter,vpot)
               firstiter=.false.
            endif
            !     calculation of the rotated quantities

            call rotate( z0t, c0(:,:), bec, c0diag, becdiag, .false. )
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
              call rotate( z0t, c0(:,:), bec, c0diag, becdiag, .false. )
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
              call nksic_potential( nbsp, nbspx, c0, fsic, bec, rhovan, deeq_sic, &
                                    ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )

              esic=sum(pink(1:nbsp))
              etot = etot + esic
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
            
             call berry_energy( enb, enbi, bec, c0(:,:), fion )
             etot=etot+enb+enbi
          endif
          if(tefield2  ) then!just in this case calculates elfield stuff at zeo field-->to be bettered

             call berry_energy2( enb, enbi, bec, c0(:,:), fion )
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
             write(1032,'(2I10,2F24.13)') 0,0,etot-esic,esic
          endif
#endif

          if(do_innerloop) then
!$$$$          if(do_innerloop.and.itercg.le.20) then
!$$$$
             !
             !call start_clock( "inner_loop" )
             !
             esic    = sum(pink(1:nbsp))
             etot    = etot - esic
             etotnew = etotnew - esic
             ninner  = 0

             if(.not.do_innerloop_cg) then
                 call nksic_rot_emin(itercg,ninner,etot,Omattot)
             else
                 call nksic_rot_emin_cg(itercg,ninner,etot,Omattot)
             endif

!$$ Now rotate hi(:,:) according to Omattot!
!$$ It seems that not rotating hi gives us better convergence.
!$$ So, we do not perform the following routine.
!$$
!            if(ninner.ge.2) then
!              hitmp(:,:) = (0.d0,0.d0)
!              do nbnd1=1,n
!                do nbnd2=1,n
!                  hitmp(:,nbnd1)=hitmp(:,nbnd1) + hi(:,nbnd2) * Omattot(nbnd2,nbnd1)
!                enddo
!              enddo
!              hi(:,:) = hitmp(:,:)
!            endif
!$$

             esic    = sum(pink(1:nbsp))
             etot    = etot + esic
             etotnew = etotnew + esic
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
        if ( ionode ) write(stdout,'(5x,"iteration =",I4,"  eff iteration =",I4,"   Etot (Ha) =",F22.14)')&
            itercg, itercgeff, etotnew
        if ( ionode .and. mod(itercg,10) == 0 ) write(stdout,"()" )
!$$


!$$ to see the outer loop energy convergence
        if (do_orbdep) then
            !
            esic = sum(pink(1:nbsp))
#ifdef __DEBUG
            if(ionode) write(1032,'(2I10,2F24.13)') ninner,itercg,etot-esic,esic
#endif
            !
        endif
!$$

        if( abs(etotnew-etotold) < conv_thr ) then
           numok=numok+1
        else 
           numok=0
        endif
        !
        if( numok >= 4 ) ltresh=.true.
        !
        etotold=etotnew
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

        do i=1,nbsp,2
!$$
          CALL start_clock( 'dforce1' )
!$$          call dforce( i, bec, betae, c0,c2,c3,rhos, nnrsx, ispin,f,n,nspin)
          call dforce( i, bec, betae, c0,c2,c3,rhos, nnrsx, ispin, faux, nbsp, nspin)
          CALL stop_clock( 'dforce1' )
!$$

          if(tefield .and. (evalue.ne.0.d0)) then
            call dforceb(c0, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            c2(1:ngw)=c2(1:ngw)+evalue*df(1:ngw)
            call dforceb(c0, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            c3(1:ngw)=c3(1:ngw)+evalue*df(1:ngw)
          endif
          !
          if(tefield2 .and. (evalue2.ne.0.d0)) then
            call dforceb(c0, i, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
            c2(1:ngw)=c2(1:ngw)+evalue2*df(1:ngw)
            call dforceb(c0, i+1, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
            c3(1:ngw)=c3(1:ngw)+evalue2*df(1:ngw)
          endif

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
              CALL nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, bec, ngw, c0(:,i), c0(:,i+1), vsicpsi )
              !
              c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
              !
              if( i+1 <= nbsp )   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
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
          if (ng0.eq.2) then
              hpsi(1,  i)=CMPLX(DBLE(hpsi(1,  i)), 0.d0)
              if(i+1 <= nbsp) then
                  hpsi(1,i+1)=CMPLX(DBLE(hpsi(1,i+1)), 0.d0)
              endif
          endif
          !
        enddo

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
!$$
        if(.not.do_orbdep) then
          call pcdaga2(c0,phi,hpsi)
        else
          call pc3(c0,hpsi)
        endif
!$$

!$$
!        if(ionode) then
!          do i=1,n
!            write(701,*) sum(phi(1:ngw,i)),sum(c0(1:ngw,i))
!          enddo
!          write(701,*) 'nhsa ',nhsa
!          write(701,*)
!        endif
!$$

        hpsi0(1:ngw,1:nbsp) = hpsi(1:ngw,1:nbsp)
        gi(1:ngw,1:nbsp)    = hpsi(1:ngw,1:nbsp)
        
        call calbec(1,nsp,eigr,hpsi,becm)
        call xminus1(hpsi,betae,dumm,becm,s_minus1,.false.)
!        call sminus1(hpsi,becm,betae)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!look if the following two lines are really needed
        call calbec(1,nsp,eigr,hpsi,becm)
!$$        call pc2(c0,bec,hpsi,becm)
!$$     
        if(.not.do_orbdep) then
          call pc2(c0,bec,hpsi,becm)
        else
          call pc3(c0,hpsi)
        endif
!$$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        call kminus1(gi,betae,ema0bg)
        if(.not.pre_state) then
           call xminus1(gi,betae,ema0bg,becm,k_minus1,.true.)
        else
           call xminus1_state(gi,betae,ema0bg,becm,k_minus1,.true.,ave_ene)
        endif
        call calbec(1,nsp,eigr,gi,becm)
!$$        call pc2(c0,bec,gi,becm)
!$$     
        if(.not.do_orbdep) then
          call pc2(c0,bec,gi,becm)
        else
          call pc3(c0,gi)
        endif
!$$

        
        if(tens) call calcmt( f, z0t, fmat0, .false. )

        call calbec(1,nsp,eigr,hpsi,bec0) 

!  calculates gamma
        gamma=0.d0
        
        if(.not.tens) then
           do i=1,nbsp
              do ig=1,ngw
                 gamma=gamma+2.d0*DBLE(CONJG(gi(ig,i))*hpsi(ig,i))
              enddo
              if (ng0.eq.2) then
                 gamma=gamma-DBLE(CONJG(gi(1,i))*hpsi(1,i))
              endif
           enddo
           
           call mp_sum( gamma, intra_image_comm )
           
           if (nvb.gt.0) then
              do i=1,nbsp
                 do is=1,nvb
                    do iv=1,nh(is)
                       do jv=1,nh(is)
                          do ia=1,na(is)
                             inl=ish(is)+(iv-1)*na(is)+ia
                             jnl=ish(is)+(jv-1)*na(is)+ia
                             gamma=gamma+ qq(iv,jv,is)*becm(inl,i)*bec0(jnl,i)
                          end do
                       end do
                    end do
                 end do
              enddo
           endif

        else

           do iss=1,nspin
              nss=nupdwn(iss)
              istart=iupdwn(iss)
              me_rot = descla( la_me_ , iss )
              np_rot = descla( la_npc_ , iss ) * descla( la_npr_ , iss )
              allocate( fmat_ ( nrlx, nudx ) )
              do ip = 1, np_rot
                 if( me_rot == ( ip - 1 ) ) then
                    fmat_ = fmat0(:,:,iss)
                 end if
                 nrl = ldim_cyclic( nss, np_rot, ip - 1 )
                 CALL mp_bcast( fmat_ , ip - 1 , intra_image_comm )
                 do i=1,nss
                    jj = ip
                    do j=1,nrl
                       do ig=1,ngw
                          gamma=gamma+2.d0*DBLE(CONJG(gi(ig,i+istart-1))*hpsi(ig,jj+istart-1))*fmat_(j,i)
                       enddo
                       if (ng0.eq.2) then
                          gamma=gamma-DBLE(CONJG(gi(1,i+istart-1))*hpsi(1,jj+istart-1))*fmat_(j,i)
                       endif
                       jj = jj + np_rot
                    enddo
                 enddo
              enddo
              deallocate( fmat_ )
           enddo
           if(nvb.gt.0) then
              do iss=1,nspin
                 nss=nupdwn(iss)
                 istart=iupdwn(iss)
                 me_rot = descla( la_me_ , iss )
                 np_rot = descla( la_npc_ , iss ) * descla( la_npr_ , iss )
                 allocate( fmat_ ( nrlx, nudx ) )
                 do ip = 1, np_rot
                    if( me_rot == ( ip - 1 ) ) then
                       fmat_ = fmat0(:,:,iss)
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
                                      gamma=gamma+ qq(iv,jv,is)*becm(inl,i+istart-1)*bec0(jnl,jj+istart-1)*fmat_(j,i)
                                   end do
                                end do
                             end do
                          enddo
                          jj = jj + np_rot
                       enddo
                    enddo
                 end do
                 deallocate( fmat_ )
              enddo
           endif
           call mp_sum( gamma, intra_image_comm )
        endif




        !case of first iteration

!$$        if(itercg==1.or.(mod(itercg,niter_cg_restart).eq.1).or.restartcg) then
        if( itercg==1 .or. mod(itercg,niter_cg_restart)==0 .or. restartcg) then
!$$

          restartcg=.false.
!$$  We do not have to reset passof every exception of CG!
!$$$$          passof=passop

          hi(1:ngw,1:nbsp)=gi(1:ngw,1:nbsp)!hi is the search direction
          esse=gamma

        else

          !find direction hi for general case 
          !calculates gamma for general case, not using Polak Ribiere
          
          essenew=gamma
          gamma=gamma/esse
          esse=essenew

          hi(1:ngw,1:nbsp)=gi(1:ngw,1:nbsp)+gamma*hi(1:ngw,1:nbsp)

        endif
!note that hi, is saved  on gi, because we need it before projection on conduction states

        !find minimum along direction hi:

        !project hi on conduction sub-space

        call calbec(1,nsp,eigr,hi,bec0)
!$$        call pc2(c0,bec,hi,bec0)
!$$     
        if(.not.do_orbdep) then
          call pc2(c0,bec,hi,bec0)
        else
          call pc3(c0,hi)
        endif
!$$
        

        !do quadratic minimization
        !             
        !calculate derivative with respect to  lambda along direction hi

        dene0=0.
        if(.not.tens) then
          do i=1,nbsp              
            do ig=1,ngw
              dene0=dene0-4.d0*DBLE(CONJG(hi(ig,i))*hpsi0(ig,i))
            enddo
            if (ng0.eq.2) then
              dene0=dene0+2.d0*DBLE(CONJG(hi(1,i))*hpsi0(1,i))
            endif
          end do
!$$ We need the following because n for spin 2 is double that for spin 1!
          dene0 = dene0 *2.d0/nspin
!$$          dene0 = dene0 *4.d0/nspin
!$$
        else
          !in the ensamble case the derivative is Sum_ij (<hi|H|Psi_j>+ <Psi_i|H|hj>)*f_ji
          !     calculation of the kinetic energy x=xmin      
         call calcmt( f, z0t, fmat0, .false. )
         do iss = 1, nspin
            nss    = nupdwn(iss)
            istart = iupdwn(iss)
            me_rot = descla( la_me_ , iss )
            np_rot = descla( la_npc_ , iss ) * descla( la_npr_ , iss )
            allocate( fmat_ ( nrlx, nudx ) )
            do ip = 1, np_rot
               if( me_rot == ( ip - 1 ) ) then
                  fmat_ = fmat0(:,:,iss)
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



!!$$$$ Calculates wavefunction at very close to c0.
!    if(.false.) then
!      tmppasso=1.d-4
!      !
!#ifdef __DEBUG
!      if(ionode) write(8000,*) itercg
!#endif
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
!#ifdef __DEBUG
!        if(ionode) then
!            write(8000,'(2e30.20,3e20.10)')  ene0,etot,dene0,tmppasso,(etot-ene0)/tmppasso/dene0
!        endif
!#endif
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
!#ifdef __DEBUG
!      if(ionode) write(8000,*)
!#endif
!      !
!    endif
!!$$$$



      !
      ! calculates wave-functions on a point on direction hi
      !
      cm(1:ngw,1:nbsp)=c0(1:ngw,1:nbsp)+spasso*passof*hi(1:ngw,1:nbsp)

      !
!$$   ! I do not know why the following 3 lines 
      ! were not in the original code (CHP)
      !
      if(ng0 == 2)  cm(1,:)=0.5d0*(cm(1,:)+CONJG(cm(1,:)))
!$$

      !orthonormalize

      call calbec(1,nsp,eigr,cm,becm)
      call gram(betae,becm,nhsa,cm,ngw,nbsp)

        !call calbec(1,nsp,eigr,cm,becm)

        !calculate energy
        if(.not.tens) then
          call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
        else
          if(newscheme) then 
              call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                        rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,.false., vpot  )  
          endif

          !     calculation of the rotated quantities
          call rotate( z0t, cm(:,:), becm, c0diag, becdiag, .false. )
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
            !
            call nksic_potential( nbsp, nbspx, cm, fsic, bec, rhovan, deeq_sic, &
                                  ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )
            !
            etot = etot + sum(pink(1:nbsp))
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

        if( ionode .and. iprsta > 1 ) write(stdout,"(6f20.12)") ene0,dene0,ene1,passo, gamma, esse

        !set new step

        passov=passof
!$$$$        passof=2.d0*passo
!$$ doing the following makes the convergence better...
        passof=passo
!$$$$
              
        !calculates wave-functions at minimum

        cm(1:ngw,1:nbsp) = c0(1:ngw,1:nbsp) +spasso*passo*hi(1:ngw,1:nbsp)
        !
        if( ng0 == 2 )  cm(1,:) = 0.5d0*(cm(1,:)+CONJG(cm(1,:)))

        call calbec(1,nsp,eigr,cm,becm)
        call gram(betae,becm,nhsa,cm,ngw,nbsp)

        !test on energy: check the energy has really diminished

        !call calbec(1,nsp,eigr,cm,becm)
        if(.not.tens) then
          call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
        else
          if(newscheme)  then
              call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,.false., vpot  )
          endif
          !     calculation of the rotated quantities
          call rotate( z0t, cm(:,:), becm, c0diag, becdiag, .false. )
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
            !
            call nksic_potential( nbsp, nbspx, cm, fsic, bec, rhovan, deeq_sic, &
                                  ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )
            etot = etot + sum(pink(1:nbsp))
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
!           bec(:,:)=becm(:,:)
          call copy_twin(bec,becm) !modified:giovanni
          ene_ok=.true.
        elseif( (enever.ge.ene1) .and. (enever.lt.ene0)) then
          if(ionode) then
             write(stdout,"(2x,a,i5)") 'cg_sub: missed minimum, case 1, iteration',itercg
          endif
          c0(1:ngw,1:nbsp)=c0(1:ngw,1:nbsp)+spasso*passov*hi(1:ngw,1:nbsp)
!$$
          passof=1.d0*passov
!$$
          restartcg=.true.
          call calbec(1,nsp,eigr,c0,bec)
          call gram(betae,bec,nhsa,c0,ngw,nbsp)

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
          call calbec(1,nsp,eigr,c0,bec)
          call gram(betae,bec,nhsa,c0,ngw,nbsp)

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
            call calbec(1,nsp,eigr,cm,becm)
            call gram(betae,bec,nhsa,cm,ngw,nbsp)
            call calbec(1,nsp,eigr,cm,becm)

            if(.not.tens) then
              call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
            else
              if(newscheme)  then
                  call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,.false., vpot  )
              endif
              !     calculation of the rotated quantities
              call rotate( z0t, cm(:,:), becm, c0diag, becdiag, .false. )
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
                !
                call nksic_potential( nbsp, nbspx, cm, fsic, bec, rhovan, deeq_sic, &
                                      ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )
                !
                etot = etot + sum(pink(1:nbsp))
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
          if(.not.do_orbdep) then
              if(iter3 == maxiter3) write(stdout,"(2x,a)") 'missed minimun: iter3 = maxiter3'
              c0(:,:)=cm(:,:)
          endif
!$$

          restartcg=.true.
          ene_ok=.false.

!$$
          if(iter3 == maxiter3) passof=passop
!$$
        end if
        
        if(tens.and.newscheme) enever=enever-entropy
 
        if(.not. ene_ok) call calbec (1,nsp,eigr,c0,bec)

        !calculates phi for pc_daga
        CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, nbsp )
  
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
             call rhoofr(nfi,c0(:,:),irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
         else
             !
             !     calculation of the rotated quantities
             call rotate( z0t, c0(:,:), bec, c0diag, becdiag, .false. )
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
             call nksic_potential( nbsp, nbspx, c0, fsic, bec, rhovan, deeq_sic, &
                                   ispin, iupdwn, nupdwn, rhor, rhog, wtot, vsic, pink )
             etot = etot + sum(pink(1:nbsp))
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


     call calcmt( f, z0t, fmat0, .false. )

     call newd(vpot,irb,eigrb,rhovan,fion)
     if (.not.tens) then
        if (tfor .or. tprnfor) call nlfq(c0,eigr,bec,becdr,fion)
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

     do i=1,nbsp,2
!$$
         CALL start_clock( 'dforce2' )
!$$          call dforce( i, bec, betae, c0,c2,c3,rhos, nnrsx, ispin,f,n,nspin)
         call dforce(i,bec,betae,c0,c2,c3,rhos,nnrsx,ispin,faux,nbsp,nspin)
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

!$$
         if ( do_orbdep ) then
             !
             ! faux takes into account spin multiplicity.
             !
             CALL nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, bec, ngw, c0(:,i), c0(:,i+1), vsicpsi )
             !
             c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
             !
             if( i+1 <= nbsp )   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
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

         do ig=1,ngw
            gi(ig,  i)=c2(ig)
            if(i+1 <= nbsp) gi(ig,i+1)=c3(ig)
         enddo
         !
         if (ng0.eq.2) then
            gi(1,  i)=CMPLX(DBLE(gi(1,  i)),0.d0)
            if(i+1 <= nbsp) gi(1,i+1)=CMPLX(DBLE(gi(1,i+1)),0.d0)
         endif

        enddo

        call init_twin(lambda_repl,lgam) !modified:giovanni
        call allocate_twin(lambda_repl, nudx, nudx, lgam) !modified:giovanni
!         ALLOCATE( lambda_repl( nudx, nudx ) )
        !
        do is = 1, nspin
           !
           nss = nupdwn(is)
           istart = iupdwn(is)
           lambda_repl = 0.d0
           !
           !
           do i = 1, nss
              do j = i, nss
                 ii = i + istart - 1
                 jj = j + istart - 1
                 do ig = 1, ngw
                    lambda_repl( i, j ) = lambda_repl( i, j ) - &
                       2.d0 * DBLE( CONJG( c0( ig, ii ) ) * gi( ig, jj) )
                 enddo
                 if( ng0 == 2 ) then
                    lambda_repl( i, j ) = lambda_repl( i, j ) + &
                       DBLE( CONJG( c0( 1, ii ) ) * gi( 1, jj ) )
                 endif
                 lambda_repl( j, i ) = lambda_repl( i, j )
              enddo
           enddo
           !
           CALL mp_sum( lambda_repl, intra_image_comm )
           !
           CALL distribute_lambda( lambda_repl, lambda( :, :, is ), descla( :, is ) )
           !
        end do

        DEALLOCATE( lambda_repl )
  
        if ( tens ) then
           !
           ! in the ensemble case matrix labda must be multiplied with f

!            ALLOCATE( lambda_dist( nlam, nlam ) )
           call init_twin(lambda_dist, lgam) !modified:giovanni
           call allocate_twin(lambda_dist,nlam,nlam,lgam) !modified:giovanni
 
           do iss = 1, nspin
              !
              nss    = nupdwn( iss )
              !
              lambdap(:,:,iss) = 0.0d0
              !
              CALL cyc2blk_redist( nss, fmat0(1,1,iss), nrlx, SIZE(fmat0,2), lambda_dist, nlam, nlam, descla(1,iss) )
              !
              ! Perform lambdap = lambda * fmat0
              !
              CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, lambda(1,1,iss), nlam, lambda_dist, nlam, &
                                  0.0d0, lambdap(1,1,iss), nlam, descla(1,iss) )
              !
              !begin_modified:giovanni
              call copy_twin(lambda_dist, lambda(iss))
              call copy_twin(lambda(iss), lambdap(iss))
              call copy_twin(lambdap(iss), lambda_dist)
              !end_modified:giovanni
              !
           end do
           !
           call deallocate_twin(lambda_dist)
           DEALLOCATE( lambda_dist )
           !
           call nlsm2(ngw,nhsa,nbsp,nspin,eigr,c0(:,:),becdr)
           !
        endif
        !
  
        !
        call nlfl(bec,becdr,lambda,fion, lgam)
          
        ! bforceion adds the force term due to electronic berry phase
        ! only in US-case
          
        if( tefield.and.(evalue .ne. 0.d0) ) then
           call bforceion(fion,tfor.or.tprnfor,ipolp, qmat,bec,becdr,gqq,evalue)

        endif
        if( tefield2.and.(evalue2 .ne. 0.d0) ) then
           call bforceion(fion,tfor.or.tprnfor,ipolp2, qmat2,bec,becdr,gqq2,evalue2)
        endif
        deallocate(hpsi0,hpsi,gi,hi)
        deallocate(hitmp)
        deallocate( s_minus1,k_minus1)

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
        do i=1,nspin
	  call deallocate_twin(lambda(i))
	  call deallocate_twin(lambdap(i))
        enddo
        !
        !end_modified:giovanni

        deallocate(ave_ene)
        deallocate(c2,c3)

        return

     END SUBROUTINE runcg_uspp

