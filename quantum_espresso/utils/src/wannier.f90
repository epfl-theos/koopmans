!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Riccardo De Gennaro: module containing some env variables of wann2kcp
!
!
!------------------------------------------------------------------------
MODULE wannier
  !
  !
  USE kinds,        ONLY : DP
  !
  !
  INTEGER :: nnb                                  ! #b
  INTEGER, ALLOCATABLE :: kpb(:,:)                ! k+b (ik,ib)
  INTEGER, ALLOCATABLE :: g_kpb(:,:,:)            ! G_k+b (ipol,ik,ib)
  INTEGER, ALLOCATABLE :: ig_(:,:)                ! G_k+b (ipol,ik,ib)
  LOGICAL, ALLOCATABLE :: excluded_band(:)
  INTEGER :: iun_nnkp, nnbx, nexband
  INTEGER :: n_wannier                            ! number of WF
  INTEGER :: n_proj                               ! number of projection
  INTEGER :: ispinw, ikstart, ikstop, iknum
  CHARACTER(LEN=15) :: wan_mode                   ! running mode
  LOGICAL :: scdm_proj
  LOGICAL :: regular_mesh = .true.
  ! input data from nnkp file
  REAL(DP), ALLOCATABLE :: center_w(:,:)          ! center_w(3,n_wannier)
  INTEGER,  ALLOCATABLE :: spin_eig(:)
  REAL(DP), ALLOCATABLE :: spin_qaxis(:,:)
  INTEGER, ALLOCATABLE  :: l_w(:), mr_w(:)        ! l and mr of wannier (n_wannier) as from table 3.1,3.2 of spec.
  INTEGER, ALLOCATABLE  :: r_w(:)                 ! index of radial function (n_wannier) as from table 3.3 of spec.
  REAL(DP), ALLOCATABLE :: xaxis(:,:),zaxis(:,:)  ! xaxis and zaxis(3,n_wannier)
  REAL(DP), ALLOCATABLE :: alpha_w(:)             ! alpha_w(n_wannier) ( called zona in wannier spec)
  !
  CHARACTER(len=256) :: seedname  = 'wannier'     ! prepended to file names in wannier90
  ! For implementation of wannier_lib
  INTEGER               :: mp_grid(3)             ! dimensions of MP k-point grid
  REAL(DP)              :: rlatt(3,3),glatt(3,3)  ! real and recip lattices (Cartesian co-ords, units of Angstrom)
  REAL(DP), ALLOCATABLE :: kpt_latt(:,:)          ! k-points in crystal co-ords. kpt_latt(3,iknum)
  INTEGER               :: num_bands              ! number of bands left after exclusions
  LOGICAL                  :: old_spinor_proj     ! for compatability for nnkp files prior to W90v2.0
  LOGICAL,ALLOCATABLE :: zerophase(:,:)
  ! wannier2kcp additional variables
  LOGICAL :: wannier_plot
  CHARACTER(LEN=255) :: wannier_plot_list
  INTEGER, ALLOCATABLE :: wann_to_plot(:)
  LOGICAL :: gamma_trick                          ! determines whether or not using SC real wfc (wannier2kcp wan_mode)
  LOGICAL :: print_rho                            ! determines whether or not writing the supercell charge density to file
  !
  !
END MODULE wannier
