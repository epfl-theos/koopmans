!
! Copyright (C) 2002-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!     n     = total number of electronic states
!     nx    = if n is even, nx=n ; if it is odd, nx=n+1
!             nx is used only to dimension arrays

!     tpiba   = 2*pi/alat
!     tpiba2  = (2*pi/alat)**2
!     ng      = number of G vectors for density and potential
!     ngl     = number of shells of G

!     G-vector quantities for the thick grid - see also doc in ggen 
!     g       = G^2 in increasing order (in units of tpiba2=(2pi/a)^2)
!     gl      = shells of G^2           ( "   "   "    "      "      )
!     gx      = G-vectors               ( "   "   "  tpiba =(2pi/a)  )
!
!     g2_g    = all G^2 in increasing order, replicated on all procs
!     mill_g  = miller index of G vecs (increasing order), replicated on all procs
!     mill_l  = miller index of G vecs local to the processors
!     ig_l2g  = "l2g" means local to global, this array convert a local
!               G-vector index into the global index, in other words
!               the index of the G-v. in the overall array of G-vectors
!     bi?     = base vector used to generate the reciprocal space
!
!     np      = fft index for G>
!     nm      = fft index for G<
!     mill_l  = G components in crystal axis
!


!
!  lqmax:  maximum angular momentum of Q (Vanderbilt augmentation charges)
! 

!  nbeta    number of beta functions (sum over all l)
!  kkbeta   last radial mesh point used to describe functions
!                 which vanish outside core
!  nqf      coefficients in Q smoothing
!  nqlc     angular momenta present in Q smoothing
!  lll      lll(j) is l quantum number of j'th beta function
!  lmaxq      highest angular momentum that is present in Q functions
!  lmaxkb   highest angular momentum that is present in beta functions
!  dion     bare pseudopotential D_{\mu,\nu} parameters
!              (ionic and screening parts subtracted out)
!  betar    the beta function on a r grid (actually, r*beta)
!  qqq      Q_ij matrix
!  qfunc    Q_ij(r) function (for r>rinner)
!  rinner   radius at which to cut off partial core or Q_ij
!
!  qfcoef   coefficients to pseudize qfunc for different total
!              angular momentum (for r<rinner)
!  vloc_at  local potential for each atom


module local_pseudo
  use kinds, only: DP
  implicit none
  save
  !
  !    rhops = ionic pseudocharges (for Ewald term)
  !    vps   = local pseudopotential in G space for each species
  !
  real(DP), allocatable:: rhops(:,:), vps(:,:)
  !
  !    drhops = derivative of rhops respect to G^2
  !    dvps   = derivative of vps respect to G^2
  !
  real(DP),allocatable:: dvps(:,:), drhops(:,:)
  !
  !    vps0  = correction factors needed to align V(0) to the "traditional"
  !            value used by other plane-wave codes - one per species
  !
  real(DP),allocatable:: vps0(:)
  !
contains
  !
  subroutine allocate_local_pseudo( ng, nsp )
      integer, intent(in) :: ng, nsp
      call deallocate_local_pseudo()
      ALLOCATE( rhops( ng, nsp ) )
      ALLOCATE( vps( ng, nsp ) )
      ALLOCATE( drhops( ng, nsp ) )
      ALLOCATE( dvps( ng, nsp ) )
      ALLOCATE( vps0( nsp ) )
  end subroutine
  !
  subroutine deallocate_local_pseudo
      IF( ALLOCATED( vps0 ) ) DEALLOCATE( vps0 )
      IF( ALLOCATED( dvps ) ) DEALLOCATE( dvps )
      IF( ALLOCATED( drhops ) ) DEALLOCATE( drhops )
      IF( ALLOCATED( vps ) ) DEALLOCATE( vps )
      IF( ALLOCATED( rhops ) ) DEALLOCATE( rhops )
  end subroutine
  !
end module local_pseudo

module qgb_mod
  USE kinds, ONLY: DP
  implicit none
  save
  complex(DP), allocatable :: qgb(:,:,:)
contains
  subroutine deallocate_qgb_mod
      IF( ALLOCATED( qgb ) ) DEALLOCATE( qgb )
  end subroutine deallocate_qgb_mod
end module qgb_mod

module qradb_mod
  USE kinds, ONLY: DP
  implicit none
  save
  real(DP), allocatable:: qradb(:,:,:,:)
contains
  subroutine deallocate_qradb_mod
      IF( ALLOCATED( qradb ) ) DEALLOCATE( qradb )
  end subroutine deallocate_qradb_mod
end module qradb_mod


MODULE metagga  !metagga
  USE kinds, ONLY: DP
  implicit none
  !the variables needed for meta-GGA
  REAL(DP), ALLOCATABLE :: &
       kedtaus(:,:), &! KineticEnergyDensity in real space,smooth grid
       kedtaur(:,:), &! real space, density grid
       crosstaus(:,:,:), &!used by stress tensor,in smooth grid
       dkedtaus(:,:,:,:)  !derivative of kedtau wrt h on smooth grid
  COMPLEX(DP) , ALLOCATABLE :: &
       kedtaug(:,:),    & !KineticEnergyDensity in G space
       gradwfc(:,:)    !used by stress tensor
contains
  subroutine deallocate_metagga
      IF( ALLOCATED(crosstaus))DEALLOCATE(crosstaus)
      IF( ALLOCATED(dkedtaus)) DEALLOCATE(dkedtaus)
      IF( ALLOCATED(gradwfc))  DEALLOCATE(gradwfc)
  end subroutine deallocate_metagga
END MODULE metagga  !end metagga

MODULE dener
  USE kinds, ONLY: DP
  IMPLICIT NONE
  SAVE
  REAL(DP) :: dekin(3,3)
  REAL(DP) :: dh(3,3)
  REAL(DP) :: dps(3,3)
  REAL(DP) :: denl(3,3)
  REAL(DP) :: dxc(3,3)
  REAL(DP) :: dsr(3,3)
  REAL(DP) :: detot(3,3)
  REAL(DP) :: dekin6(6)
  REAL(DP) :: dh6(6)
  REAL(DP) :: dps6(6)
  REAL(DP) :: denl6(6)
  REAL(DP) :: dxc6(6)
  REAL(DP) :: dsr6(6)
  REAL(DP) :: detot6(6)
END MODULE dener

module dqgb_mod
  USE kinds, ONLY: DP
  implicit none
  save
  complex(DP),allocatable:: dqgb(:,:,:,:,:)
contains
  subroutine deallocate_dqgb_mod
      IF( ALLOCATED( dqgb ) ) DEALLOCATE( dqgb )
  end subroutine deallocate_dqgb_mod
end module dqgb_mod

MODULE cdvan
  USE kinds, ONLY: DP
  IMPLICIT NONE
  SAVE
  REAL(DP), ALLOCATABLE :: dbeta(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: dbec(:,:,:,:)     
    ! Warning dbec is distributed over row and column processors of the ortho group
  REAL(DP), ALLOCATABLE :: drhovan(:,:,:,:,:)
CONTAINS
  SUBROUTINE deallocate_cdvan
      IF( ALLOCATED( dbeta ) ) DEALLOCATE( dbeta )
      IF( ALLOCATED( dbec ) ) DEALLOCATE( dbec )
      IF( ALLOCATED( drhovan ) ) DEALLOCATE( drhovan )
  END SUBROUTINE deallocate_cdvan
END MODULE cdvan



module cvan

  ! this file contains common subroutines and modules between
  ! CP and FPMD

  !     ionic pseudo-potential variables
  use parameters, only: nsx
  implicit none
  save
  integer nvb, ish(nsx)
  !     nvb    = number of species with Vanderbilt PPs
  !     ish(is)= used for indexing the nonlocal projectors betae
  !              with contiguous indices inl=ish(is)+(iv-1)*na(is)+1
  !              where "is" is the species and iv=1,nh(is)
  !
  !     indlm: indlm(ind,is)=Y_lm for projector ind
  integer, allocatable:: indlm(:,:)
contains

  subroutine allocate_cvan( nind, ns )
    integer, intent(in) :: nind, ns
    allocate( indlm( nind, ns ) )
  end subroutine allocate_cvan

  subroutine deallocate_cvan( )
    if( allocated(indlm) ) deallocate( indlm )
  end subroutine deallocate_cvan

end module cvan


MODULE stress_param

   USE kinds, ONLY : DP

   IMPLICIT NONE
   SAVE

   INTEGER, DIMENSION(6), PARAMETER :: alpha = (/ 1,2,3,2,3,3 /)
   INTEGER, DIMENSION(6), PARAMETER :: beta  = (/ 1,1,1,2,2,3 /)

   REAL(DP),  DIMENSION(3,3), PARAMETER :: delta = reshape &
         ( (/ 1.0_DP, 0.0_DP, 0.0_DP, &
              0.0_DP, 1.0_DP, 0.0_DP, &
              0.0_DP, 0.0_DP, 1.0_DP  &
            /), (/ 3, 3 /) )

   ! ...  dalbe(:) = delta(alpha(:),beta(:))
   !
   REAL(DP),  DIMENSION(6), PARAMETER :: dalbe = &
         (/ 1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP /)

END MODULE



MODULE core
   !
   USE kinds
   ! 
   IMPLICIT NONE
   SAVE
   !     nlcc_any = 0 no core correction on any atom
   !     rhocb  = core charge in G space (box grid)
   !     rhoc   = core charge in real space  (dense grid)
   !     rhocg  = core charge in G space  (dense grid)
   !     drhocg = derivative of core charge in G space (used for stress)
   !
   LOGICAL :: nlcc_any
   REAL(DP), ALLOCATABLE:: rhocb(:,:)
   REAL(DP), ALLOCATABLE:: rhoc(:)
   REAL(DP), ALLOCATABLE:: rhocg(:,:)
   REAL(DP), ALLOCATABLE:: drhocg(:,:)
   !
CONTAINS
   !
   SUBROUTINE allocate_core( nnrx, ngm, ngb, nsp ) 
     INTEGER, INTENT(IN) :: nnrx, ngm, ngb, nsp
     IF ( nlcc_any ) THEN    
        !
        ALLOCATE( rhoc( nnrx ) )
        ALLOCATE( rhocb( ngb, nsp ) )
        ALLOCATE( rhocg( ngm, nsp ) )
        ALLOCATE( drhocg( ngm, nsp ) )
        !
     ELSE
        !
        ! ... dummy allocation required because this array appears in the
        ! ... list of arguments of some routines
        !
        ALLOCATE( rhoc( 1 ) )
        !
     END IF
   END SUBROUTINE allocate_core
   !
   SUBROUTINE deallocate_core()
      IF( ALLOCATED( rhocb  ) ) DEALLOCATE( rhocb )
      IF( ALLOCATED( rhoc   ) ) DEALLOCATE( rhoc  )
      IF( ALLOCATED( rhocg  ) ) DEALLOCATE( rhocg  )
      IF( ALLOCATED( drhocg ) ) DEALLOCATE( drhocg )
   END SUBROUTINE deallocate_core
   !
END MODULE core
!
module ldaU
  use parameters, only: nsx
  USE kinds
  implicit none
  complex(DP), allocatable :: atomwfc(:,:)
  complex(DP), allocatable :: swfcatom(:,:)
  real(DP) :: Hubbard_U(nsx), Hubbard_lambda(nsx,2), ns0(nsx,2),   &
     & Hubbard_alpha(nsx)
  real(DP) :: e_hubbard = 0.d0, e_lambda = 0.d0
  real(DP), allocatable :: ns(:,:,:,:)
  integer :: Hubbard_l(nsx), Hubbard_lmax=0, n_atomic_wfc
  logical lda_plus_u
  COMPLEX(DP), allocatable::  vupsi(:,:) !@@@@
contains
  !
  subroutine deallocate_lda_plus_u()
     !
     IF( ALLOCATED( atomwfc ) ) DEALLOCATE( atomwfc )
     IF( ALLOCATED( swfcatom ) ) DEALLOCATE( swfcatom )
     IF( ALLOCATED( ns ) ) DEALLOCATE( ns )
     IF( ALLOCATED( vupsi ) ) DEALLOCATE( vupsi )
     !
     !
  end subroutine
  !
end module ldaU
!
module nksic
  !
  use kinds
  implicit none
  save
  !
  real(dp) :: fref
  real(dp) :: rhobarfact
  real(dp) :: nkscalfact
  real(dp) :: vanishing_rho_w
  real(dp) :: f_cutoff
  !
  real(dp) :: etxc
  !
  real(dp),    allocatable :: fsic(:)
  real(dp),    allocatable :: vsic(:,:)
  real(dp),    allocatable :: fion_sic(:,:)
  real(dp),    allocatable :: deeq_sic(:,:,:,:)
  real(dp),    allocatable :: pink(:)
  real(dp),    allocatable :: vxcsic(:,:)
  real(dp),    allocatable :: wxdsic(:,:)
  real(dp),    allocatable :: orb_rhor(:,:)
  real(dp),    allocatable :: rhoref(:,:)
  real(dp),    allocatable :: rhobar(:,:)
  real(dp),    allocatable :: grhobar(:,:,:)
  real(dp),    allocatable :: wrefsic(:)
  !
  complex(dp), allocatable :: vsicpsi(:,:)
  !
  integer :: nknmax
  logical :: do_orbdep
  logical :: do_nk
  logical :: do_pz
  logical :: do_nkpz
  logical :: do_nki
  logical :: do_spinsym
  logical :: do_wxd
  logical :: do_wref

contains
  !
  subroutine allocate_nksic( nnrx, ngw, nspin, nx, nat)
      !
      use funct,          only : dft_is_gradient
      USE uspp_param,     only : nhm
      !
      implicit none
      integer, intent(in):: nx, nspin
      integer, intent(in):: nat
      integer, intent(in):: ngw
      integer, intent(in):: nnrx
      !
      allocate( fsic(nx) )
      allocate( vsic(nnrx,nx) )
      allocate( fion_sic(3,nat) )
      allocate( deeq_sic(nhm, nhm, nat, nx) )
      allocate( pink(nx) )
      allocate( vsicpsi(ngw,2) )
      allocate( vxcsic(nnrx,2) )
      allocate( orb_rhor(nnrx,2) )
      allocate( rhoref(nnrx,2) )
      allocate( rhobar(nnrx,2) )
      !
      if ( do_nk .or. do_nkpz ) then
          allocate( wxdsic(nnrx,2) )
          allocate( wrefsic(nnrx) )
      else if ( do_nki ) then
          allocate( wxdsic(nnrx,2) )
      endif
      !
      if ( dft_is_gradient() ) then 
          allocate( grhobar(nnrx,3,2) )
      else
          allocate( grhobar(1,1,1) )
      endif
      !
      fsic     = 0.0d0
      pink     = 0.0d0
      vsic     = 0.0d0
      vsicpsi  = 0.0d0
      !
      !vxcsic   = 0.0d0
      !wxdsic   = 0.0d0
      !wrefsic  = 0.0d0
      !orb_rhor = 0.0d0
      !rhobar   = 0.0d0
      !rhoref   = 0.0d0
      !
  end subroutine allocate_nksic
  !
  real(dp) function nksic_memusage( )
      ! output in MB (according to 4B integers and 8B reals)  
      real(dp) :: cost
      !
      cost = 0.0_dp
      if ( allocated(fsic) )       cost = cost + real( size(fsic) )       *  8.0_dp 
      if ( allocated(vsic) )       cost = cost + real( size(vsic) )       *  8.0_dp 
      if ( allocated(fion_sic) )   cost = cost + real( size(fion_sic) )   *  8.0_dp 
      if ( allocated(deeq_sic) )   cost = cost + real( size(deeq_sic) )   *  8.0_dp 
      if ( allocated(pink) )       cost = cost + real( size(pink) )       *  8.0_dp 
      if ( allocated(vsicpsi) )    cost = cost + real( size(vsicpsi) )    * 16.0_dp 
      if ( allocated(vxcsic) )     cost = cost + real( size(vxcsic) )     *  8.0_dp 
      if ( allocated(wxdsic) )     cost = cost + real( size(wxdsic) )     *  8.0_dp 
      if ( allocated(orb_rhor))    cost = cost + real( size(orb_rhor))    *  8.0_dp
      if ( allocated(rhoref) )     cost = cost + real( size(rhoref) )     *  8.0_dp
      if ( allocated(rhobar) )     cost = cost + real( size(rhobar) )     *  8.0_dp
      if ( allocated(grhobar) )    cost = cost + real( size(grhobar) )    *  8.0_dp
      if ( allocated(wrefsic) )    cost = cost + real( size(wrefsic) )    *  8.0_dp
      !
      nksic_memusage = cost / 1000000.0_dp
      !   
  end function nksic_memusage
  !
  subroutine deallocate_nksic
      !
      if(allocated(vsic))        deallocate(vsic)
      if(allocated(fion_sic))    deallocate(fion_sic)
      if(allocated(deeq_sic))    deallocate(deeq_sic)
      if(allocated(pink))        deallocate(pink)
      if(allocated(wxdsic))      deallocate(wxdsic)
      if(allocated(vxcsic))      deallocate(vxcsic)
      if(allocated(vsicpsi))     deallocate(vsicpsi)
      if(allocated(wrefsic))     deallocate(wrefsic)
      if(allocated(orb_rhor))    deallocate(orb_rhor)
      if(allocated(grhobar))     deallocate(grhobar)
      if(allocated(rhobar))      deallocate(rhobar)
      if(allocated(rhoref))      deallocate(rhoref)
      !
  end subroutine deallocate_nksic
  !
end module nksic


module hfmod
  !
  use kinds
  implicit none
  complex(dp), allocatable :: vxxpsi(:,:)
  real(dp), allocatable :: exx(:)
  real(dp), allocatable :: dvxchf(:,:)
  real(dp) :: detothf = 0.d0
  real(dp) :: hfscalfact = 1.d0
  logical  :: do_hf
contains
  !
  subroutine allocate_hf(ngw,nx)
  implicit none
  integer, intent(in):: nx
  integer, intent(in):: ngw
  allocate(exx(nx))
  allocate(vxxpsi(ngw,nx))
  !
  exx(:) = 0.0
  !
  end subroutine
  !
  subroutine deallocate_hf
  if(allocated(exx)) deallocate(exx)
  if(allocated(vxxpsi)) deallocate(vxxpsi)
  end subroutine
  !
end module hfmod

module eecp_mod
  USE kinds
  implicit none
  real(dp), allocatable :: gcorr(:)
  complex(dp), allocatable :: gcorr_fft(:)
  real(dp), allocatable :: gcorr1d(:)
  complex(dp), allocatable :: gcorr1d_fft(:)
  real(dp), allocatable :: vcorr(:)
  complex(dp), allocatable :: vcorr_fft(:)
  character (len=256) :: which_compensation
  logical :: do_comp
  real(dp) :: ecomp
contains
  !
  subroutine allocate_ee(nnrx,ngm)
  implicit none
  integer, intent(in):: nnrx
  integer, intent(in):: ngm
  allocate(gcorr(nnrx))
  allocate(gcorr_fft(nnrx))
  allocate(gcorr1d(nnrx))
  allocate(gcorr1d_fft(nnrx))
  allocate(vcorr(nnrx))
  allocate(vcorr_fft(ngm))
  end subroutine
  !
  subroutine deallocate_ee
  if(allocated(gcorr)) deallocate(gcorr)
  if(allocated(gcorr_fft)) deallocate(gcorr_fft)
  if(allocated(gcorr1d)) deallocate(gcorr1d)
  if(allocated(gcorr1d_fft)) deallocate(gcorr1d_fft)
  if(allocated(vcorr)) deallocate(vcorr)
  if(allocated(vcorr_fft)) deallocate(vcorr_fft)
  end subroutine
  !
end module eecp_mod

module efield_mod
  USE kinds
  implicit none
  real(dp) :: ampfield(3)
  logical :: do_efield
  real(dp), allocatable :: efieldpot(:)
  complex(dp), allocatable :: efieldpotg(:)
contains
  !
  subroutine allocate_efield(nnrx,ngm)
  implicit none
  integer, intent(in):: nnrx
  integer, intent(in):: ngm
  allocate(efieldpot(nnrx))
  allocate(efieldpotg(ngm))
  end subroutine
  !
  subroutine deallocate_efield
  if(allocated(efieldpot)) deallocate(efieldpot)
  if(allocated(efieldpotg)) deallocate(efieldpotg)
  end subroutine
  !
end module efield_mod
!
!
! Occupation constraint ...to be implemented...
!
module step_constraint
  USE kinds
  implicit none
  integer, parameter :: natx_ = 5000
  real(DP) :: E_con
  real(DP) :: A_con(natx_,2), sigma_con(natx_), alpha_con(natx_)
  logical :: step_con
  ! complex(DP), allocatable:: vpsi_con(:,:)
  complex(DP) :: vpsi_con(1,1)
end module step_constraint

