!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
MODULE cp_main_variables
  !----------------------------------------------------------------------------
  !
  USE kinds,             ONLY : DP
  USE control_flags,     ONLY : program_name
  USE funct,             ONLY : dft_is_meta
  USE metagga,           ONLY : kedtaur, kedtaus, kedtaug
  USE cell_base,         ONLY : boxdimensions
  USE wave_types,        ONLY : wave_descriptor, wave_descriptor_init
  USE energies,          ONLY : dft_energy_type
  USE pres_ai_mod,       ONLY : abivol, abisur, jellium, t_gauss, rho_gaus, &
                                v_vol, posv, f_vol
  USE twin_types
  !
  IMPLICIT NONE
  SAVE
  !
  ! ... structure factors e^{-ig*R}
  !
  ! ...  G = reciprocal lattice vectors
  ! ...  R_I = ionic positions
  !
  COMPLEX(DP), ALLOCATABLE :: eigr(:,:)        ! exp (i G   dot R_I)
  COMPLEX(DP), ALLOCATABLE :: ei1(:,:)         ! exp (i G_x dot x_I)
  COMPLEX(DP), ALLOCATABLE :: ei2(:,:)         ! exp (i G_y dot y_I)
  COMPLEX(DP), ALLOCATABLE :: ei3(:,:)         ! exp (i G_z dot z_I)
  !
  ! ... structure factors (summed over atoms of the same kind)
  !
  ! S( s, G ) = sum_(I in s) exp( i G dot R_(s,I) )
  ! s       = index of the atomic specie
  ! R_(s,I) = position of the I-th atom of the "s" specie
  !
  COMPLEX(DP), ALLOCATABLE:: sfac(:,:)
  !
  ! ... indexes, positions, and structure factors for the box grid
  !
  REAL(DP), ALLOCATABLE :: taub(:,:)
  COMPLEX(DP), ALLOCATABLE :: eigrb(:,:)
  INTEGER,     ALLOCATABLE :: irb(:,:)
  ! 
  ! ... nonlocal projectors:
  ! ...    bec   = scalar product of projectors and wave functions
  ! ...    betae = nonlocal projectors in g space = beta x e^(-ig.R) 
  ! ...    becdr = <betae|g|psi> used in force calculation
  ! ...    rhovan= \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
  ! ...    deeq  = \int V_eff(r) q_lm(r) dr
  !
  type(twin_tensor) :: becdr!(:,:,:)!,bec(:,:) !modified:giovanni
!   REAL(DP), ALLOCATABLE :: bephi(:,:)!, becp(:,:) !removed:giovanni
  TYPE(twin_matrix)     :: bec, becp, bephi, becdual !added:giovanni
  TYPE(twin_matrix), ALLOCATABLE :: lambda(:)
  TYPE(twin_matrix), ALLOCATABLE :: lambdam(:)
  TYPE(twin_matrix), ALLOCATABLE :: lambdap(:)
  !
  ! For non-orthogonal SIC, naive implementation :giovanni
  !
  COMPLEX(DP), ALLOCATABLE :: overlap(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: ioverlap(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: kinetic_mat(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: pseudopot_mat(:,:,:)
  !
  ! ... mass preconditioning
  !
  REAL(DP), ALLOCATABLE :: ema0bg(:)
  !
  ! ... constraints (lambda at t, lambdam at t-dt, lambdap at t+dt)
  !
!   REAL(DP), ALLOCATABLE :: hamilt(:,:,:) !removed:giovanni
  type(twin_matrix), dimension(:), allocatable :: hamilt(:)
  !
  INTEGER,  ALLOCATABLE :: descla(:,:) ! descriptor of the lambda distribution
                                       ! see descriptors_module
  INTEGER :: nlax = 0                  ! leading dimension of the distribute (by block) lambda matrix 
  INTEGER :: nlam = 1                  ! dimension of lambda matrix, can be 1 or nlax depending on la_proc
  INTEGER :: nrlx = 0                  ! leading dimension of the distribute (by row  ) lambda matrix
  LOGICAL :: la_proc = .FALSE.         ! indicate if a proc own a block of lambda
  !
  INTEGER, PARAMETER :: nacx = 10      ! max number of averaged
                                       ! quantities saved to the restart
  REAL(DP) :: acc(nacx)
  REAL(DP) :: acc_this_run(nacx)
  !
  ! cell geometry
  !
  TYPE (boxdimensions) :: htm, ht0, htp  ! cell metrics
  !
  ! charge densities and potentials
  !
  ! rhog  = charge density in g space
  ! rhor  = charge density in r space (dense grid)
  ! rhos  = charge density in r space (smooth grid)
  ! vpot  = potential in r space (dense grid)
  !
  COMPLEX(DP), ALLOCATABLE :: rhog(:,:)
  REAL(DP),    ALLOCATABLE :: rhor(:,:), rhos(:,:)
  REAL(DP),    ALLOCATABLE :: vpot(:,:)
  !
  ! derivative wrt cell
  !
  COMPLEX(DP), ALLOCATABLE :: drhog(:,:,:,:)
  REAL(DP),    ALLOCATABLE :: drhor(:,:,:,:)

  TYPE (wave_descriptor) :: wfill     ! wave function descriptor for filled
  !
  TYPE(dft_energy_type) :: edft
  !
  INTEGER :: nfi             ! counter on the electronic iterations
  INTEGER :: nprint_nfi=-1   ! counter indicating the last time data have been
                             ! printed on file ( prefix.pos, ... )
  INTEGER :: nfi_run=0       ! counter on the electronic iterations,
                             ! for the present run
  INTEGER :: iprint_stdout=1 ! define how often CP writes verbose information to stdout
  !
  INTERFACE distribute_bec
     module procedure distribute_bec_cmplx, distribute_bec_real
  END INTERFACE

  INTERFACE collect_lambda
      module procedure collect_lambda_real, collect_lambda_cmplx
  END INTERFACE

  INTERFACE collect_zmat
      module procedure collect_zmat_real, collect_zmat_cmplx
  END INTERFACE

  INTERFACE distribute_lambda
      module procedure distribute_lambda_real, distribute_lambda_cmplx
  END INTERFACE

  INTERFACE distribute_zmat
      module procedure distribute_zmat_real, distribute_zmat_cmplx
  END INTERFACE

  INTERFACE setval_lambda
      module procedure setval_lambda_real, setval_lambda_cmplx
  END INTERFACE

  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE allocate_mainvar( ngw, ngwt, ngb, ngs, ng, nr1, nr2, nr3, &
                                 nr1x, nr2x, npl, nnr, nnrsx, nat, nax,  &
                                 nsp, nspin, n, nx, n_emp, nupdwn, nhsa, &
                                 gzero, nudx, tpre )
      !------------------------------------------------------------------------
      !
      USE mp_global,   ONLY: np_ortho, me_ortho, intra_image_comm, ortho_comm, &
                             me_image, ortho_comm_id
      USE mp,          ONLY: mp_max, mp_min
      USE descriptors, ONLY: descla_siz_ , descla_init , nlax_ , la_nrlx_ , lambda_node_
      USE control_flags, ONLY: do_wf_cmplx, gamma_only, non_ortho! added:giovanni
      USE twin_types
      !
      INTEGER,           INTENT(IN) :: ngw, ngwt, ngb, ngs, ng, nr1, nr2, nr3, &
                                       nnr, nnrsx, nat, nax, nsp, nspin, &
                                       n, nx, n_emp, nhsa, nr1x, nr2x, npl
      INTEGER,           INTENT(IN) :: nupdwn(:)
      LOGICAL,           INTENT(IN) :: gzero
      INTEGER,           INTENT(IN) :: nudx
      LOGICAL,           INTENT(IN) :: tpre
      !
      INTEGER  :: iss, nhsa_l
      LOGICAL  :: lgam !added:giovanni

      lgam=gamma_only.and..not.do_wf_cmplx !added:giovanni
      !
      ! ... allocation of all arrays not already allocated in init and nlinit
      !
      ALLOCATE( eigr( ngw, nat ) )
      ALLOCATE( sfac( ngs, nsp ) )
      ALLOCATE( ei1( -nr1:nr1, nat ) )
      ALLOCATE( ei2( -nr2:nr2, nat ) )
      ALLOCATE( ei3( -nr3:nr3, nat ) )
      ALLOCATE( eigrb( ngb, nat ) )
      ALLOCATE( irb( 3, nat ) )
      !
      IF ( dft_is_meta() ) THEN
         !
         ! ... METAGGA
         !
         ALLOCATE( kedtaur( nnr,   nspin ) )
         ALLOCATE( kedtaus( nnrsx, nspin ) )
         ALLOCATE( kedtaug( ng,    nspin ) )
         !
      ELSE
         !
         ! ... dummy allocation required because this array appears in the
         ! ... list of arguments of some routines
         !
         ALLOCATE( kedtaur( 1, nspin ) )
         ALLOCATE( kedtaus( 1, nspin ) )
         ALLOCATE( kedtaug( 1, nspin ) )
         !
      END IF
      !
      ALLOCATE( ema0bg( ngw ) )
      !
      ALLOCATE( rhor( nnr, nspin ) )
      ALLOCATE( vpot( nnr, nspin ) )
      ALLOCATE( rhos( nnrsx, nspin ) )
      ALLOCATE( rhog( ng,    nspin ) )
      IF( program_name == 'CP90' .AND. tpre ) THEN
         IF ( tpre ) THEN
            ALLOCATE( drhog( ng,  nspin, 3, 3 ) )
            ALLOCATE( drhor( nnr, nspin, 3, 3 ) )
         END IF
      END IF
      !
      !  Compute local dimensions for lambda matrixes
      !

      ALLOCATE( descla( descla_siz_ , nspin ) )
      !
      nlax = 0
      nrlx = 0
      DO iss = 1, nspin
         CALL descla_init( descla( :, iss ), nupdwn( iss ), nudx, np_ortho, me_ortho, ortho_comm, ortho_comm_id )
         nlax = MAX( nlax, descla( nlax_ , iss ) )
         nrlx = MAX( nrlx, descla( la_nrlx_ , iss ) )
         IF( descla( lambda_node_ , iss ) > 0 ) la_proc = .TRUE.
      END DO
      !
      nlam = 1
      IF( la_proc ) nlam = nlax
      !
      !  ... End with lambda dimensions
      !
      !
      IF( program_name == 'CP90' ) THEN
         !
         if ( abivol.or.abisur ) then
            !
            allocate(rho_gaus(nnr))
            allocate(v_vol(nnr))
            if (jellium.or.t_gauss) allocate(posv(3,nr1*nr2*nr3))
            if (t_gauss) allocate(f_vol(3,nax,nsp))
            !
         end if
         !
      END IF
      !
!!!! begin_modified:giovanni 
      ALLOCATE( lambda(nspin ) )
      ALLOCATE( lambdam(nspin ) )
      ALLOCATE( lambdap(nspin ) )
      DO iss=1,nspin
          lambda(iss)%iscmplx=.not.lgam
          call init_twin(lambda(iss), lgam)
          call allocate_twin(lambda(iss), nlam, nlam, lgam)
          call init_twin(lambdap(iss), lgam)
          call allocate_twin(lambdap(iss), nlam, nlam, lgam)
          call init_twin(lambdam(iss), lgam)
          call allocate_twin(lambdam(iss), nlam, nlam, lgam)
      ENDDO
      
      IF(non_ortho) THEN
         allocate(ioverlap(nudx,nudx,nspin), &
           overlap(nudx,nudx,nspin), kinetic_mat(nudx,nudx,nspin), &
         pseudopot_mat(nudx,nudx,nspin))
      ENDIF
!!!! end_modified:giovanni
      !
      ! becdr, distributed over row processors of the ortho group
      !
      ! AF: this avoids problems when nhsa, i.e. nkb, is zero
      !nhsa_l = MAX( nhsa, 1) 
      nhsa_l = nhsa
      !
      call init_twin(becdr, lgam)
      call allocate_twin(becdr, nhsa_l, nspin*nlax, 3, lgam) !added:giovanni
      !
      call init_twin(bec, lgam)
      call allocate_twin(bec,nhsa_l,n, lgam)!added:giovanni
      call init_twin(becp, lgam)
      call allocate_twin(becp,nhsa_l,n, lgam)!added:giovanni
      call init_twin(bephi, lgam)
      call allocate_twin(bephi, nhsa_l, nspin*nlax, lgam)!added:giovanni
      
      IF(non_ortho) THEN
         call init_twin(becdual, lgam)
         call allocate_twin(becdual,nhsa_l,n, lgam)!added:giovanni
      ENDIF
!       ALLOCATE( bec( nhsa_l,n ) )!removed:giovanni
      !
!       ALLOCATE( bephi( nhsa_l, nspin*nlax ) ) !removed:giovanni
!       ALLOCATE( becp(  nhsa_l, n ) ) !removed:giovanni
      !
      CALL wave_descriptor_init( wfill, ngw, ngwt, nupdwn,  nupdwn, &
            1, 1, nspin, 'gamma', gzero )
      !
      RETURN
      !
    END SUBROUTINE allocate_mainvar
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_mainvar()
      use io_global, only: ionode
      IMPLICIT NONE

      INTEGER :: iss
      !------------------------------------------------------------------------
      !
      IF( ALLOCATED( ei1 ) )     DEALLOCATE( ei1 )
      IF( ALLOCATED( ei2 ) )     DEALLOCATE( ei2 )
      IF( ALLOCATED( ei3 ) )     DEALLOCATE( ei3 )
      IF( ALLOCATED( eigr ) )    DEALLOCATE( eigr )
      IF( ALLOCATED( sfac ) )    DEALLOCATE( sfac )
      IF( ALLOCATED( eigrb ) )   DEALLOCATE( eigrb )
      IF( ALLOCATED( irb ) )     DEALLOCATE( irb )
      IF( ALLOCATED( rhor ) )    DEALLOCATE( rhor )
      IF( ALLOCATED( rhos ) )    DEALLOCATE( rhos )
      IF( ALLOCATED( rhog ) )    DEALLOCATE( rhog )
      IF( ALLOCATED( drhog ) )   DEALLOCATE( drhog )
      IF( ALLOCATED( drhor ) )   DEALLOCATE( drhor )
      IF( ALLOCATED( ema0bg ) )  DEALLOCATE( ema0bg )
      IF( ALLOCATED( kedtaur ) ) DEALLOCATE( kedtaur )
      IF( ALLOCATED( kedtaus ) ) DEALLOCATE( kedtaus )
      IF( ALLOCATED( kedtaug ) ) DEALLOCATE( kedtaug )
      IF( ALLOCATED( vpot ) )    DEALLOCATE( vpot )
      IF( ALLOCATED( taub ) )    DEALLOCATE( taub )
      IF( ALLOCATED( descla ) )  DEALLOCATE( descla )
      !added:giovanni -- deallocation of structured types -- the check is inside deallocate
      return
if(ionode) then
!       write(0,*) "debug1"
endif
      CALL deallocate_twin(bec)
      CALL deallocate_twin(becp)
      CALL deallocate_twin(becdr) 
      CALL deallocate_twin(bephi)
      ! for non orthogonal case
      CALL deallocate_twin(becdual)

if(ionode) then
! write(0,*) "debug2"
endif
      IF(allocated(hamilt)) THEN 
	  DO iss=1, size(hamilt)
	      CALL deallocate_twin(hamilt(iss))
	  END DO
          DEALLOCATE(hamilt)
      ENDIF
if(ionode) then
! write(0,*) "debug3"
endif
      IF(allocated(lambda)) THEN 
	  DO iss=1, size(lambda)
	      CALL deallocate_twin(lambda(iss))
	  END DO
          DEALLOCATE(lambda)
      ENDIF
! if(ionode) then
! write(0,*) "debug4"
! endif
      IF(allocated(lambdam)) THEN 
	  DO iss=1, size(lambdam)
	      CALL deallocate_twin(lambdam(iss))
	  END DO
          DEALLOCATE(lambdam)
      ENDIF
! if(ionode) then
! write(0,*) "debug5"
! endif
      IF(allocated(lambdap)) THEN 
	  DO iss=1, size(lambdap)
	      CALL deallocate_twin(lambdap(iss))
	  END DO
          DEALLOCATE(lambdap)
      ENDIF

! if(ionode) then
! write(0,*) "debug6"
! endif
      ! -- deallocation of structured types -------------------------------------------------------------
      
      IF(allocated(overlap)) THEN
         deallocate(overlap)
      ENDIF

      IF(allocated(ioverlap)) THEN
         deallocate(ioverlap)
      ENDIF

      IF(allocated(kinetic_mat)) THEN
         deallocate(kinetic_mat)
      ENDIF

      IF(allocated(pseudopot_mat)) THEN
         deallocate(pseudopot_mat)
      ENDIF

      RETURN
      !
    END SUBROUTINE deallocate_mainvar
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE distribute_lambda_real( lambda_repl, lambda_dist, desc )
       USE descriptors, ONLY: lambda_node_ , ilar_ , ilac_ , nlac_ , nlar_
       REAL(DP), INTENT(IN)  :: lambda_repl(:,:)
       REAL(DP), INTENT(OUT) :: lambda_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:)
       INTEGER :: i, j, ic, ir
       IF( desc( lambda_node_ ) > 0 ) THEN
          ir = desc( ilar_ )       
          ic = desc( ilac_ )       
          DO j = 1, desc( nlac_ )
             DO i = 1, desc( nlar_ )
                lambda_dist( i, j ) = lambda_repl( i + ir - 1, j + ic - 1 )
             END DO
          END DO
       END IF
       RETURN
    END SUBROUTINE distribute_lambda_real

    SUBROUTINE distribute_lambda_cmplx( lambda_repl, lambda_dist, desc )
       USE descriptors, ONLY: lambda_node_ , ilar_ , ilac_ , nlac_ , nlar_
       COMPLEX(DP), INTENT(IN)  :: lambda_repl(:,:)
       COMPLEX(DP), INTENT(OUT) :: lambda_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:)
       INTEGER :: i, j, ic, ir
       IF( desc( lambda_node_ ) > 0 ) THEN
          ir = desc( ilar_ )       
          ic = desc( ilac_ )       
          DO j = 1, desc( nlac_ )
             DO i = 1, desc( nlar_ )
                lambda_dist( i, j ) = lambda_repl( i + ir - 1, j + ic - 1 )
             END DO
          END DO
       END IF
       RETURN
    END SUBROUTINE distribute_lambda_cmplx
    !
    !------------------------------------------------------------------------
    SUBROUTINE distribute_bec_cmplx( bec_repl, bec_dist, desc, nspin )
       USE descriptors, ONLY: lambda_node_ , ilar_ , nlar_ , la_n_ , nlax_
       COMPLEX(DP), INTENT(IN)  :: bec_repl(:,:)
       COMPLEX(DP) :: bec_dist(:,:) !modified:giovanni
       INTEGER,  INTENT(IN)  :: desc(:,:)
       INTEGER,  INTENT(IN)  :: nspin
       INTEGER :: i, ir, n, nlax
       CHARACTER(len=21) :: subname = "distribute_bec_cmplx"
       !
!        IF(.not.bec_dist%iscmplx) THEN
!            call errore(subname, "incompatible types", 1)
!        ENDIF

       IF( desc( lambda_node_ , 1 ) > 0 ) THEN
          !
          bec_dist = CMPLX(0.0d0,0.d0)
          !
          ir = desc( ilar_ , 1 )
          DO i = 1, desc( nlar_ , 1 )
             bec_dist( :, i ) = bec_repl( :, i + ir - 1 )
          END DO
          !
          IF( nspin == 2 ) THEN
             n     = desc( la_n_ , 1 )  !  number of states with spin 1 ( nupdw(1) )
             nlax  = desc( nlax_ , 1 )   !  array elements reserved for each spin ( bec(:,2*nlax) )
             ir = desc( ilar_ , 2 )
             DO i = 1, desc( nlar_ , 2 )
                bec_dist( :, i + nlax ) = bec_repl( :, i + ir - 1 + n )
             END DO
          END IF
          !
       END IF
       RETURN
    END SUBROUTINE distribute_bec_cmplx
    !
    !-----------------------------------------------------------------------------------------------
    SUBROUTINE distribute_bec_real( bec_repl, bec_dist, desc, nspin )
       USE descriptors, ONLY: lambda_node_ , ilar_ , nlar_ , la_n_ , nlax_
       REAL(DP), INTENT(IN)  :: bec_repl(:,:)
       REAL(DP) :: bec_dist(:,:) !modified:giovanni
       INTEGER,  INTENT(IN)  :: desc(:,:)
       INTEGER,  INTENT(IN)  :: nspin
       INTEGER :: i, ir, n, nlax
       CHARACTER(len=21) :: subname = "distribute_bec_real"
       !
!        IF(bec_dist%iscmplx) THEN
!            call errore(subname, "incompatible types", 1)
!        ENDIF

       IF( desc( lambda_node_ , 1 ) > 0 ) THEN
          !
          bec_dist = 0.0d0
          !
          ir = desc( ilar_ , 1 )
          DO i = 1, desc( nlar_ , 1 )
             bec_dist( :, i ) = bec_repl( :, i + ir - 1 )
          END DO
          !
          IF( nspin == 2 ) THEN
             n     = desc( la_n_ , 1 )  !  number of states with spin 1 ( nupdw(1) )
             nlax  = desc( nlax_ , 1 )   !  array elements reserved for each spin ( bec(:,2*nlax) )
             ir = desc( ilar_ , 2 )
             DO i = 1, desc( nlar_ , 2 )
                bec_dist( :, i + nlax ) = bec_repl( :, i + ir - 1 + n )
             END DO
          END IF
          !
       END IF
       RETURN
    END SUBROUTINE distribute_bec_real
    !
    !------------------------------------------------------------------------
    SUBROUTINE distribute_zmat_real( zmat_repl, zmat_dist, desc )
       USE descriptors, ONLY: lambda_node_ , la_nrl_ , la_me_ , la_npr_ , la_npc_ , la_n_
       REAL(DP), INTENT(IN)  :: zmat_repl(:,:)
       REAL(DP), INTENT(OUT) :: zmat_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:)
       INTEGER :: i, ii, j, me, np
       me = desc( la_me_ )
       np = desc( la_npc_ ) * desc( la_npr_ )
       IF( desc( lambda_node_ ) > 0 ) THEN
          DO j = 1, desc( la_n_ )
             ii = me + 1
             DO i = 1, desc( la_nrl_ )
                zmat_dist( i, j ) = zmat_repl( ii, j )
                ii = ii + np
             END DO
          END DO
       END IF
       RETURN
    END SUBROUTINE distribute_zmat_real

    SUBROUTINE distribute_zmat_cmplx( zmat_repl, zmat_dist, desc )
       USE descriptors, ONLY: lambda_node_ , la_nrl_ , la_me_ , la_npr_ , la_npc_ , la_n_
       COMPLEX(DP), INTENT(IN)  :: zmat_repl(:,:)
       COMPLEX(DP), INTENT(OUT) :: zmat_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:)
       INTEGER :: i, ii, j, me, np
       me = desc( la_me_ )
       np = desc( la_npc_ ) * desc( la_npr_ )
       IF( desc( lambda_node_ ) > 0 ) THEN
          DO j = 1, desc( la_n_ )
             ii = me + 1
             DO i = 1, desc( la_nrl_ )
                zmat_dist( i, j ) = zmat_repl( ii, j )
                ii = ii + np
             END DO
          END DO
       END IF
       RETURN
    END SUBROUTINE distribute_zmat_cmplx
    !
    !------------------------------------------------------------------------
    SUBROUTINE collect_lambda_real( lambda_repl, lambda_dist, desc )
       USE mp_global,   ONLY: intra_image_comm
       USE mp,          ONLY: mp_sum
       USE descriptors, ONLY: lambda_node_ , ilar_ , ilac_ , nlac_ , nlar_
       REAL(DP), INTENT(OUT) :: lambda_repl(:,:)
       REAL(DP), INTENT(IN)  :: lambda_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:)
       INTEGER :: i, j, ic, ir
       lambda_repl = 0.0d0
       IF( desc( lambda_node_ ) > 0 ) THEN
          ir = desc( ilar_ )       
          ic = desc( ilac_ )       
          DO j = 1, desc( nlac_ )
             DO i = 1, desc( nlar_ )
                lambda_repl( i + ir - 1, j + ic - 1 ) = lambda_dist( i, j )
             END DO
          END DO
       END IF
       CALL mp_sum( lambda_repl, intra_image_comm )
       RETURN
    END SUBROUTINE collect_lambda_real

    SUBROUTINE collect_lambda_cmplx( lambda_repl, lambda_dist, desc )
       USE mp_global,   ONLY: intra_image_comm
       USE mp,          ONLY: mp_sum
       USE descriptors, ONLY: lambda_node_ , ilar_ , ilac_ , nlac_ , nlar_
       COMPLEX(DP), INTENT(OUT) :: lambda_repl(:,:)
       COMPLEX(DP), INTENT(IN)  :: lambda_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:)
       INTEGER :: i, j, ic, ir

       lambda_repl = CMPLX(0.0d0,0.d0)
       IF( desc( lambda_node_ ) > 0 ) THEN
          ir = desc( ilar_ )       
          ic = desc( ilac_ )       
          DO j = 1, desc( nlac_ )
             DO i = 1, desc( nlar_ )
                lambda_repl( i + ir - 1, j + ic - 1 ) = lambda_dist( i, j )
             END DO
          END DO
       END IF
       CALL mp_sum( lambda_repl, intra_image_comm )
       RETURN
    END SUBROUTINE collect_lambda_cmplx
    !
    !------------------------------------------------------------------------
    SUBROUTINE collect_bec( bec_repl, bec_dist, desc, nspin )
       USE mp_global,   ONLY: intra_image_comm
       USE mp,          ONLY: mp_sum
       USE descriptors, ONLY: lambda_node_ , ilar_ , nlar_ , la_myc_ , nlax_ , la_n_
       REAL(DP), INTENT(OUT) :: bec_repl(:,:)
       REAL(DP), INTENT(IN)  :: bec_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:,:)
       INTEGER,  INTENT(IN)  :: nspin
       INTEGER :: i, ir, n, nlax, iss
       !
       bec_repl = 0.0d0
       !
       !  bec is distributed across row processor, the first column is enough
       !
       IF( ( desc( lambda_node_ , 1 ) > 0 ) .AND. ( desc( la_myc_ , 1 ) == 1 ) ) THEN
          ir = desc( ilar_ , 1 )
          DO i = 1, desc( nlar_ , 1 )
             bec_repl( :, i + ir - 1 ) = bec_dist( :, i )
          END DO
          IF( nspin == 2 ) THEN
             n  = desc( la_n_ , 1 )   ! number of states with spin==1 ( nupdw(1) )
             nlax = desc( nlax_ , 1 ) ! array elements reserved for each spin ( bec(:,2*nlax) )
             ir = desc( ilar_ , 2 )
             DO i = 1, desc( nlar_ , 2 )
                bec_repl( :, i + ir - 1 + n ) = bec_dist( :, i + nlax )
             END DO
          END IF
       END IF
       !
       CALL mp_sum( bec_repl, intra_image_comm )
       !
       RETURN
    END SUBROUTINE collect_bec
    !
    !------------------------------------------------------------------------
    SUBROUTINE collect_zmat_real( zmat_repl, zmat_dist, desc )
       USE mp_global,   ONLY: intra_image_comm
       USE mp,          ONLY: mp_sum
       USE descriptors, ONLY: lambda_node_ , la_nrl_ , la_me_ , la_npr_ , la_npc_ , la_n_
       REAL(DP), INTENT(OUT) :: zmat_repl(:,:)
       REAL(DP), INTENT(IN)  :: zmat_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:)
       INTEGER :: i, ii, j, me, np, nrl
       zmat_repl = 0.0d0
       me = desc( la_me_ )
       np = desc( la_npc_ ) * desc( la_npr_ )
       nrl = desc( la_nrl_ )
       IF( desc( lambda_node_ ) > 0 ) THEN
          DO j = 1, desc( la_n_ )
             ii = me + 1
             DO i = 1, nrl
                zmat_repl( ii, j ) = zmat_dist( i, j )
                ii = ii + np
             END DO
          END DO
       END IF
       CALL mp_sum( zmat_repl, intra_image_comm )
       RETURN
    END SUBROUTINE collect_zmat_real
    !
    !------------------------------------------------------------------------
    SUBROUTINE collect_zmat_cmplx( zmat_repl, zmat_dist, desc )
       USE mp_global,   ONLY: intra_image_comm
       USE mp,          ONLY: mp_sum
       USE descriptors, ONLY: lambda_node_ , la_nrl_ , la_me_ , la_npr_ , la_npc_ , la_n_
       COMPLEX(DP), INTENT(OUT) :: zmat_repl(:,:)
       COMPLEX(DP), INTENT(IN)  :: zmat_dist(:,:)
       INTEGER,  INTENT(IN)  :: desc(:)
       INTEGER :: i, ii, j, me, np, nrl
       zmat_repl = CMPLX(0.0d0,0.d0)
       me = desc( la_me_ )
       np = desc( la_npc_ ) * desc( la_npr_ )
       nrl = desc( la_nrl_ )
       IF( desc( lambda_node_ ) > 0 ) THEN
          DO j = 1, desc( la_n_ )
             ii = me + 1
             DO i = 1, nrl
                zmat_repl( ii, j ) = zmat_dist( i, j )
                ii = ii + np
             END DO
          END DO
       END IF
       CALL mp_sum( zmat_repl, intra_image_comm )
       RETURN
    END SUBROUTINE collect_zmat_cmplx
    !
    !------------------------------------------------------------------------
    SUBROUTINE setval_lambda_real( lambda_dist, i, j, val, desc )
       USE descriptors, ONLY: lambda_node_ , ilar_ , ilac_ , nlac_ , nlar_
       REAL(DP), INTENT(OUT) :: lambda_dist(:,:)
       INTEGER,  INTENT(IN)  :: i, j
       REAL(DP), INTENT(IN)  :: val
       INTEGER,  INTENT(IN)  :: desc(:)
       IF( desc( lambda_node_ ) > 0 ) THEN
          IF( ( i >= desc( ilar_ ) ) .AND. ( i - desc( ilar_ ) + 1 <= desc( nlar_ ) ) ) THEN
             IF( ( j >= desc( ilac_ ) ) .AND. ( j - desc( ilac_ ) + 1 <= desc( nlac_ ) ) ) THEN
                lambda_dist( i - desc( ilar_ ) + 1, j - desc( ilac_ ) + 1 ) = val
             END IF
          END IF
       END IF
       RETURN
    END SUBROUTINE setval_lambda_real

    SUBROUTINE setval_lambda_cmplx( lambda_dist, i, j, val, desc )
       USE descriptors, ONLY: lambda_node_ , ilar_ , ilac_ , nlac_ , nlar_
       COMPLEX(DP), INTENT(OUT) :: lambda_dist(:,:)
       INTEGER,  INTENT(IN)  :: i, j
       COMPLEX(DP), INTENT(IN)  :: val
       INTEGER,  INTENT(IN)  :: desc(:)
       IF( desc( lambda_node_ ) > 0 ) THEN
          IF( ( i >= desc( ilar_ ) ) .AND. ( i - desc( ilar_ ) + 1 <= desc( nlar_ ) ) ) THEN
             IF( ( j >= desc( ilac_ ) ) .AND. ( j - desc( ilac_ ) + 1 <= desc( nlac_ ) ) ) THEN
                lambda_dist( i - desc( ilar_ ) + 1, j - desc( ilac_ ) + 1 ) = val
             END IF
          END IF
       END IF
       RETURN
    END SUBROUTINE setval_lambda_cmplx
    !
    !
END MODULE cp_main_variables
