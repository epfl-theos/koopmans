!
! Copyright (C) 2002-2008 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!=----------------------------------------------------------------------------=!
     SUBROUTINE interpolate_lambda_real_x( lambdap, lambda, lambdam )
!=----------------------------------------------------------------------------=!
       USE kinds, ONLY: DP
       IMPLICIT NONE
       REAL(DP) :: lambdap(:,:,:), lambda(:,:,:), lambdam(:,:,:) 
       !
       ! interpolate new lambda at (t+dt) from lambda(t) and lambda(t-dt):
       !
       lambdap= 2.d0*lambda - lambdam
       lambdam=lambda 
       lambda =lambdap
       RETURN
     END SUBROUTINE interpolate_lambda_real_x
!=----------------------------------------------------------------------------=!
     SUBROUTINE interpolate_lambda_twin_x( lambdap, lambda, lambdam )
!=----------------------------------------------------------------------------=!
       USE kinds, ONLY: DP
       USE twin_types
   
       IMPLICIT NONE
       TYPE(twin_matrix), dimension(:) :: lambdap, lambda, lambdam
       !
       INTEGER :: i
       !
       ! interpolate new lambda at (t+dt) from lambda(t) and lambda(t-dt):
       !
       DO i=1,size(lambda)
         IF(.not.lambda(i)%iscmplx) THEN
	    lambdap(i)%rvec= 2.d0*lambda(i)%rvec - lambdam(i)%rvec
	    lambdam(i)%rvec=lambda(i)%rvec 
	    lambda(i)%rvec =lambdap(i)%rvec
         ELSE
	    lambdap(i)%cvec= 2.d0*lambda(i)%cvec - lambdam(i)%cvec
	    lambdam(i)%cvec=lambda(i)%cvec 
	    lambda(i)%cvec =lambdap(i)%cvec
         ENDIF
       ENDDO
       RETURN
     END SUBROUTINE interpolate_lambda_twin_x

!=----------------------------------------------------------------------------=!
     SUBROUTINE update_lambda_x( i, lambda, c0, c2, n, noff, tdist )
!=----------------------------------------------------------------------------=!
       USE kinds,              ONLY: DP
       USE electrons_module,   ONLY: ib_owner, ib_local
       USE mp_global,          ONLY: me_image, intra_image_comm
       USE mp,                 ONLY: mp_sum
       USE wave_base,          ONLY: hpsi
       USE reciprocal_vectors, ONLY: gzero
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n, noff
       REAL(DP)            :: lambda(:,:)
       COMPLEX(DP)         :: c0(:,:), c2(:)
       INTEGER, INTENT(IN) :: i
       LOGICAL, INTENT(IN) :: tdist   !  if .true. lambda is distributed
       !
       REAL(DP), ALLOCATABLE :: prod(:)
       INTEGER :: ibl
       !
       ALLOCATE( prod( n ) )
       prod = hpsi( gzero, c0, SIZE( c0, 1 ), c2, n, noff )
       CALL mp_sum( prod, intra_image_comm )
       IF( tdist ) THEN
          IF( me_image == ib_owner( i ) ) THEN
             ibl = ib_local( i )
             lambda( ibl, : ) = prod( : )
          END IF
       ELSE
          lambda( i, : ) = prod( : )
       END IF
       DEALLOCATE( prod )
       RETURN
     END SUBROUTINE update_lambda_x

!=----------------------------------------------------------------------------=!
  subroutine elec_fakekine_x( ekincm, ema0bg, emass, c0, cm, ngw, n, noff, delt )
!=----------------------------------------------------------------------------=!
    !
    !  This subroutine computes the CP(fake) wave functions kinetic energy
    
    USE kinds,              only : DP
    use mp,                 only : mp_sum
    use mp_global,          only : intra_image_comm
    use reciprocal_vectors, only : gstart
    use wave_base,          only : wave_speed2
    use control_flags,      only: gamma_only, do_wf_cmplx!added:giovanni
    !
    IMPLICIT NONE
    !
    integer, intent(in)      :: ngw    !  number of plane wave coeff.
    integer, intent(in)      :: n      !  number of bands
    integer, intent(in)      :: noff   !  offset for band index
    real(DP), intent(out)    :: ekincm
    real(DP), intent(in)     :: ema0bg( ngw ), delt, emass
    complex(DP), intent(in)  :: c0( ngw, n ), cm( ngw, n )
    !
    real(DP), allocatable :: emainv(:)
    real(DP) :: ftmp
    integer  :: i
    logical :: lgam!added:giovanni

    lgam=gamma_only.and..not.do_wf_cmplx!added:giovanni
    ALLOCATE( emainv( ngw ) )
    emainv = 1.0d0 / ema0bg
    ftmp = 1.0d0
    if( gstart == 2 ) ftmp = 0.5d0

    ekincm=0.0d0
    do i = noff, n + noff - 1
      ekincm = ekincm + 2.0d0 * wave_speed2( c0(:,i), cm(:,i), emainv, ftmp, lgam )!added:giovanni lgam
    end do
    ekincm = ekincm * emass / ( delt * delt )

    CALL mp_sum( ekincm, intra_image_comm )
    DEALLOCATE( emainv )

    return
  end subroutine elec_fakekine_x

!=----------------------------------------------------------------------------=!
   SUBROUTINE protate_real_x ( c0, bec, c0rot, becrot, ngwl, nss, noff, lambda, &
                        na, nsp, ish, nh, np_rot, me_rot )
!=----------------------------------------------------------------------------=!

      !  this routine rotates the wave functions using the matrix lambda
      !  it works with a block-like distributed matrix
      !  of the Lagrange multipliers ( lambda ).
      !  no replicated data are used, allowing scalability for large problems.
      !  the layout of lambda is as follows :
      !
      !  (PE 0)                 (PE 1)               ..  (PE NPE-1)
      !  lambda(1      ,1:nx)   lambda(2      ,1:nx) ..  lambda(NPE      ,1:nx)
      !  lambda(1+  NPE,1:nx)   lambda(2+  NPE,1:nx) ..  lambda(NPE+  NPE,1:nx)
      !  lambda(1+2*NPE,1:nx)   lambda(2+2*NPE,1:nx) ..  lambda(NPE+2*NPE,1:nx)
      !
      !  distributes lambda's rows across processors with a blocking factor
      !  of 1, ( row 1 to PE 1, row 2 to PE 2, .. row nproc_image+1 to PE 1 and
      !  so on).
      !  nrl = local number of rows
      !  ----------------------------------------------

      ! ... declare modules

      USE kinds,            ONLY: DP
      USE mp,               ONLY: mp_bcast
      USE mp_global,        ONLY: intra_image_comm
      USE dspev_module,     ONLY: pdspev_drv, dspev_drv

      IMPLICIT NONE

      ! ... declare subroutine arguments

      INTEGER, INTENT(IN) :: ngwl, nss, noff
      INTEGER, INTENT(IN) :: na(:), nsp, ish(:), nh(:)
      INTEGER, INTENT(IN) :: np_rot, me_rot
      COMPLEX(DP), INTENT(IN) :: c0(:,:)
      COMPLEX(DP), INTENT(OUT) :: c0rot(:,:)
      REAL(DP), INTENT(IN) :: lambda(:,:)
      REAL(DP), INTENT(IN) :: bec(:,:)
      REAL(DP), INTENT(OUT) :: becrot(:,:)

      ! ... declare other variables
      INTEGER   :: i, j, ip
      INTEGER   :: jl, nrl_ip, is, ia, jv, jnl
      REAL(DP), ALLOCATABLE :: uu(:,:)

      IF( nss < 1 ) THEN
        RETURN
      END IF

      CALL start_clock('protate')

      DO i = 1, nss
         c0rot( :, i+noff-1 ) = 0.0d0
         becrot(:,i+noff-1 ) = 0.0d0
      END DO



!      becrot = 0.0d0
!      c0rot  = 0.0d0

         DO ip = 1, np_rot

            nrl_ip = nss / np_rot
            IF( (ip-1) < mod( nss, np_rot ) ) THEN
              nrl_ip = nrl_ip + 1
            END IF
 
            ALLOCATE( uu( nrl_ip, nss ) )
            IF( me_rot .EQ. (ip-1) ) THEN
              uu = lambda( 1:nrl_ip, 1:nss )
            END IF
            CALL mp_bcast( uu, (ip-1), intra_image_comm)
 
            j      = ip
            DO jl = 1, nrl_ip
              DO i = 1, nss
                CALL DAXPY(2*ngwl,uu(jl,i),c0(1,j+noff-1),1,c0rot(1,i+noff-1),1)
              END DO

              do is=1,nsp
                 do jv=1,nh(is)
                    do ia=1,na(is)
                       jnl=ish(is)+(jv-1)*na(is)+ia
                       do i = 1, nss
                          becrot(jnl,i+noff-1) = becrot(jnl,i+noff-1)+ uu(jl, i) * bec( jnl, j+noff-1 )
                       end do
                    end do
                 end do
              end do

              j = j + np_rot
            END DO

            DEALLOCATE(uu)
 
         END DO

      CALL stop_clock('protate')

      RETURN
   END SUBROUTINE protate_real_x

!=----------------------------------------------------------------------------=!
   SUBROUTINE protate_cmplx_x ( c0, bec, c0rot, becrot, ngwl, nss, noff, lambda, &
                        na, nsp, ish, nh, np_rot, me_rot  )
!=----------------------------------------------------------------------------=!

      !  this routine rotates the wave functions using the matrix lambda
      !  it works with a block-like distributed matrix
      !  of the Lagrange multipliers ( lambda ).
      !  no replicated data are used, allowing scalability for large problems.
      !  the layout of lambda is as follows :
      !
      !  (PE 0)                 (PE 1)               ..  (PE NPE-1)
      !  lambda(1      ,1:nx)   lambda(2      ,1:nx) ..  lambda(NPE      ,1:nx)
      !  lambda(1+  NPE,1:nx)   lambda(2+  NPE,1:nx) ..  lambda(NPE+  NPE,1:nx)
      !  lambda(1+2*NPE,1:nx)   lambda(2+2*NPE,1:nx) ..  lambda(NPE+2*NPE,1:nx)
      !
      !  distributes lambda's rows across processors with a blocking factor
      !  of 1, ( row 1 to PE 1, row 2 to PE 2, .. row nproc_image+1 to PE 1 and
      !  so on).
      !  nrl = local number of rows
      !  ----------------------------------------------

      ! ... declare modules

      USE kinds,            ONLY: DP
      USE mp,               ONLY: mp_bcast
      USE mp_global,        ONLY: intra_image_comm
      USE dspev_module,     ONLY: pdspev_drv, dspev_drv

      IMPLICIT NONE

      ! ... declare subroutine arguments

      INTEGER, INTENT(IN) :: ngwl, nss, noff
      INTEGER, INTENT(IN) :: na(:), nsp, ish(:), nh(:)
      INTEGER, INTENT(IN) :: np_rot, me_rot
      COMPLEX(DP), INTENT(IN) :: c0(:,:)
      COMPLEX(DP), INTENT(OUT) :: c0rot(:,:)
      COMPLEX(DP), INTENT(IN) :: lambda(:,:)
      COMPLEX(DP), INTENT(IN) :: bec(:,:)
      COMPLEX(DP), INTENT(OUT) :: becrot(:,:)

      ! ... declare other variables
      INTEGER   :: i, j, ip
      INTEGER   :: jl, nrl_ip, is, ia, jv, jnl
      COMPLEX(DP), ALLOCATABLE :: uu(:,:)

      IF( nss < 1 ) THEN
        RETURN
      END IF

      CALL start_clock('protate')

      DO i = 1, nss
         c0rot( :, i+noff-1 ) = 0.0d0
         becrot(:,i+noff-1 ) = 0.0d0
      END DO



!      becrot = 0.0d0
!      c0rot  = 0.0d0

         DO ip = 1, np_rot

            nrl_ip = nss / np_rot
            IF( (ip-1) < mod( nss, np_rot ) ) THEN
              nrl_ip = nrl_ip + 1
            END IF
 
            ALLOCATE( uu( nrl_ip, nss ) )
            IF( me_rot .EQ. (ip-1) ) THEN
              uu = lambda( 1:nrl_ip, 1:nss )
            END IF
            CALL mp_bcast( uu, (ip-1), intra_image_comm)
 
            j      = ip
            DO jl = 1, nrl_ip
              DO i = 1, nss
                CALL ZAXPY(ngwl,CONJG(uu(jl,i)),c0(1,j+noff-1),1,c0rot(1,i+noff-1),1)
              END DO

              do is=1,nsp
                 do jv=1,nh(is)
                    do ia=1,na(is)
                       jnl=ish(is)+(jv-1)*na(is)+ia
                       do i = 1, nss
                          becrot(jnl,i+noff-1) = becrot(jnl,i+noff-1)+ CONJG(uu(jl, i)) * bec( jnl, j+noff-1 )
                       end do
                    end do
                 end do
              end do

              j = j + np_rot
            END DO

            DEALLOCATE(uu)
 
         END DO

      CALL stop_clock('protate')

      RETURN
   END SUBROUTINE protate_cmplx_x

!=----------------------------------------------------------------------------=!
   SUBROUTINE crot_gamma2_real ( c0rot, c0, ngw, n, noffr, noff, lambda, nx, eig )
!=----------------------------------------------------------------------------=!

      !  this routine rotates the wave functions to the Kohn-Sham base
      !  it works with a block-like distributed matrix
      !  of the Lagrange multipliers ( lambda ).
      !
      ! ... declare modules

      USE kinds,            ONLY: DP
      USE mp,               ONLY: mp_bcast
      USE dspev_module,     ONLY: dspev_drv

      IMPLICIT NONE

      ! ... declare subroutine arguments

      INTEGER,     INTENT(IN)    :: ngw, n, nx, noffr, noff
      COMPLEX(DP), INTENT(INOUT) :: c0rot(:,:)
      COMPLEX(DP), INTENT(IN)    :: c0(:,:)
      REAL(DP),    INTENT(IN)    :: lambda(:,:)
      REAL(DP),    INTENT(OUT)   :: eig(:)

      ! ... declare other variables
      !
      REAL(DP), ALLOCATABLE :: vv(:,:), ap(:)
      INTEGER               :: i, j, k

      IF( nx < 1 ) THEN
        RETURN
      END IF

      ALLOCATE( vv( nx, nx ) )

      ! NON distributed lambda

      ALLOCATE( ap( nx * ( nx + 1 ) / 2 ) )

      K = 0
      DO J = 1, n
         DO I = J, n
            K = K + 1
            ap( k ) = lambda( i, j )
         END DO
      END DO

      CALL dspev_drv( 'V', 'L', n, ap, eig, vv, nx )

      DEALLOCATE( ap )

      DO i = 1, n
         c0rot( :, i+noffr-1 ) = 0.0d0
      END DO

      DO j = 1, n
         DO i = 1, n
            CALL DAXPY( 2*ngw, vv(j,i), c0(1,j+noff-1), 1, c0rot(1,i+noffr-1), 1 )
         END DO
      END DO

      DEALLOCATE( vv )

      RETURN
   END SUBROUTINE crot_gamma2_real

!=----------------------------------------------------------------------------=!
   SUBROUTINE crot_gamma2_cmplx ( c0rot, c0, ngw, n, noffr, noff, lambda, nx, eig )
!=----------------------------------------------------------------------------=!

      !  this routine rotates the wave functions to the Kohn-Sham base
      !  it works with a block-like distributed matrix
      !  of the Lagrange multipliers ( lambda ).
      !
      ! ... declare modules

      USE kinds,            ONLY: DP
      USE mp,               ONLY: mp_bcast
      USE zhpev_module,     ONLY: zhpev_drv

      IMPLICIT NONE

      ! ... declare subroutine arguments

      INTEGER,     INTENT(IN)    :: ngw, n, nx, noffr, noff
      COMPLEX(DP), INTENT(INOUT) :: c0rot(:,:)
      COMPLEX(DP), INTENT(IN)    :: c0(:,:)
      COMPLEX(DP),    INTENT(IN)    :: lambda(:,:)
      REAL(DP),    INTENT(OUT)   :: eig(:)

      ! ... declare other variables
      !
      COMPLEX(DP), ALLOCATABLE :: vv(:,:), ap(:)
      INTEGER               :: i, j, k

      IF( nx < 1 ) THEN
        RETURN
      END IF

      ALLOCATE( vv( nx, nx ) )

      ! NON distributed lambda

      ALLOCATE( ap( nx * ( nx + 1 ) / 2 ) )

      K = 0
      DO J = 1, n
         DO I = J, n
            K = K + 1
            ap( k ) = lambda( i, j )
         END DO
      END DO

      CALL zhpev_drv( 'V', 'L', n, ap, eig, vv, nx )

      DEALLOCATE( ap )

      DO i = 1, n
         c0rot( :, i+noffr-1 ) = 0.0d0
      END DO

      DO j = 1, n
         DO i = 1, n
            CALL ZAXPY( ngw, vv(j,i), c0(1,j+noff-1), 1, c0rot(1,i+noffr-1), 1 )
         END DO
      END DO

      DEALLOCATE( vv )

      RETURN
   END SUBROUTINE crot_gamma2_cmplx


!=----------------------------------------------------------------------------=!
   SUBROUTINE proj_gamma( a, b, ngw, n, noff, lambda)
!=----------------------------------------------------------------------------=!

        !  projection A=A-SUM{B}<B|A>B
        !  no replicated data are used, allowing scalability for large problems.
        !  The layout of lambda is as follows :
        !
        !  (PE 0)                 (PE 1)               ..  (PE NPE-1)
        !  lambda(1      ,1:nx)   lambda(2      ,1:nx) ..  lambda(NPE      ,1:nx)
        !  lambda(1+  NPE,1:nx)   lambda(2+  NPE,1:nx) ..  lambda(NPE+  NPE,1:nx)
        !  lambda(1+2*NPE,1:nx)   lambda(2+2*NPE,1:nx) ..  lambda(NPE+2*NPE,1:nx)
        !
        !  distribute lambda's rows across processors with a blocking factor
        !  of 1, ( row 1 to PE 1, row 2 to PE 2, .. row nproc_image+1 to PE 1 and so on).
        !  ----------------------------------------------
         
! ...   declare modules
        USE kinds,              ONLY: DP
        USE mp_global,          ONLY: nproc_image, me_image
        USE wave_base,          ONLY: dotp
        USE reciprocal_vectors, ONLY: gzero

        IMPLICIT NONE

! ...   declare subroutine arguments
        INTEGER,     INTENT( IN )  :: ngw, n, noff
        COMPLEX(DP), INTENT(INOUT) :: a(:,:), b(:,:)
        REAL(DP),    OPTIONAL      :: lambda(:,:)

! ...   declare other variables
        REAL(DP), ALLOCATABLE :: ee(:)
        INTEGER :: i, j
        COMPLEX(DP) :: alp

! ... end of declarations
!  ----------------------------------------------

        IF( n < 1 ) THEN
          RETURN
        END IF

        ALLOCATE( ee( n ) )
        DO i = 1, n
          DO j = 1, n
            ee(j) = -dotp( gzero, ngw, b(:,j+noff-1), a(:,i+noff-1) )
          END DO
          IF( PRESENT(lambda) ) THEN
            IF( MOD( (i-1), nproc_image ) == me_image ) THEN
              DO j = 1, n
                lambda( (i-1) / nproc_image + 1, j ) = ee(j)
              END DO
            END IF
          END IF
          DO j = 1, n
            alp = CMPLX(ee(j),0.0d0)
            CALL ZAXPY( ngw, alp, b(1,j+noff-1), 1, a(1,i+noff-1), 1 )
          END DO
        END DO
        DEALLOCATE(ee)

        RETURN
   END SUBROUTINE proj_gamma

!=----------------------------------------------------------------------------=!
   SUBROUTINE wave_atom_init_x( cm, n, noff )
!=----------------------------------------------------------------------------=!

      !  this routine sets the initial atomic wavefunctions

! ... declare modules
      USE kinds,              ONLY: DP
      USE mp,                 ONLY: mp_sum
      USE mp_wave,            ONLY: splitwf
      USE reciprocal_vectors, ONLY: ngw
      USE random_numbers,     ONLY: randy
      USE control_flags,      ONLY: ampre, tatomicwfc, trane
      USE uspp_param,         ONLY: n_atom_wfc
      USE ions_base,          ONLY: nat, ityp
      USE cp_main_variables,  ONLY: eigr
      USE reciprocal_vectors, ONLY: ngw, gx, ngwx
      USE constants,          ONLY: tpi      

      IMPLICIT NONE

      ! ... declare subroutine arguments 
      INTEGER,     INTENT(IN)  :: n, noff
      COMPLEX(DP), INTENT(INOUT) :: cm(:,:)

      ! ... declare other variables
      INTEGER :: ig, natomwfc
      REAL(DP) :: rr, arg
      COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:) ! atomic wfcs for initialization
      INTEGER :: ibnd, n_starting_wfc, n_starting_atomic_wfc
      CHARACTER(len=15) :: subname="wave_atom_init"

      ! ... Check array dimensions
      IF( SIZE( cm, 1 ) < ngw ) THEN 
        CALL errore(' wave_atom_init ', ' wrong dimensions ', 3)
      END IF

      ! ... Reset them to zero
 
      cm( :, noff : noff + n - 1 ) = 0.0d0

      ! ... initialize the wave functions in such a way that the values
      ! ... of the components are independent on the number of processors

      ampre = 0.01d0
      write(6,*) "computing n_atom_wfc"
      natomwfc = n_atom_wfc(nat, ityp, .false.) ! third value is noncolin, which is not yet implemented
      !
      IF ( tatomicwfc ) THEN
         !
         n_starting_wfc = MAX( natomwfc, n )
         n_starting_atomic_wfc = natomwfc
         !
      ENDIF

      ALLOCATE( wfcatom( ngwx, n_starting_wfc ) )

      IF ( tatomicwfc ) THEN
         !
      write(6,*) "calling atomic_wfc"
         CALL atomic_wfc( eigr, n_starting_atomic_wfc, wfcatom(:,1:n_starting_atomic_wfc) )
      write(6,*) "called atomic_wfc"
         !
         IF ( trane .AND. &
            n_starting_wfc == n_starting_atomic_wfc ) THEN
            !
            ! ... in this case, introduce a small randomization of wavefunctions
            ! ... to prevent possible "loss of states"
            !
            DO ibnd = 1, n_starting_atomic_wfc
               !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig, rr, arg)
               DO ig = 1, ngw
                  !           
                  rr  = randy()
                  arg = tpi * randy()
                  !
                  wfcatom(ig,ibnd) = wfcatom(ig,ibnd) * &
                     ( 1.0_DP + 0.05_DP * CMPLX( rr*COS(arg), rr*SIN(arg)) )
                  !
               END DO
!$OMP END PARALLEL DO
               !
            END DO
            !
         END IF
         !
      END IF 
      !      
      ! ... if not enough atomic wfc are available,
      ! ... fill missing wfcs with random numbers
      !      
      DO ibnd = n_starting_atomic_wfc + 1, n_starting_wfc
         !      
         !      
         wfcatom(:,ibnd) = (0.0_dp, 0.0_dp)
         !      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig, rr, arg)
         DO ig = 1, ngw
            !      
            rr  = randy()
            arg = tpi * randy()
            !
            wfcatom(ig,ibnd) = &
                CMPLX( rr*COS( arg ), rr*SIN( arg ) ) / &
                       ( ( gx(1,ig) )**2 + &
                         ( gx(2,ig) )**2 + &
                         ( gx(3,ig) )**2 + 1.0_DP )
         END DO
!$OMP END PARALLEL DO
         !
      END DO

      IF(n_starting_wfc.lt.n) THEN
         call errore(subname, "warning: wavefunction not fully initialized", 1)
      ENDIF
      cm( :, noff : noff + min(n,n_starting_wfc) - 1 ) = wfcatom(:,1:min(n,n_starting_wfc))
      deallocate(wfcatom)

      RETURN

    END SUBROUTINE wave_atom_init_x

!=----------------------------------------------------------------------------=!
   SUBROUTINE wave_rand_init_x( cm, n, noff )
!=----------------------------------------------------------------------------=!

      !  this routine sets the initial wavefunctions at random

! ... declare modules
      USE kinds,              ONLY: DP
      USE mp,                 ONLY: mp_sum
      USE mp_wave,            ONLY: splitwf
      USE reciprocal_vectors, ONLY: ig_l2g, ngw, gzero, ngwt
      USE random_numbers,     ONLY: randy
      
      IMPLICIT NONE

      ! ... declare subroutine arguments 
      INTEGER,     INTENT(IN)  :: n, noff
      COMPLEX(DP), INTENT(OUT) :: cm(:,:)

      ! ... declare other variables
      INTEGER :: ntest, ig, ib
      REAL(DP) ::  rranf1, rranf2, ampre
      COMPLEX(DP), ALLOCATABLE :: pwt( : )

      ! ... Check array dimensions
      IF( SIZE( cm, 1 ) < ngw ) THEN 
        CALL errore(' wave_rand_init ', ' wrong dimensions ', 3)
      END IF

      ! ... Reset them to zero
 
      cm( :, noff : noff + n - 1 ) = 0.0d0

      ! ... initialize the wave functions in such a way that the values
      ! ... of the components are independent on the number of processors

      ampre = 0.01d0
      ALLOCATE( pwt( ngwt ) )

      ntest = ngwt / 4
      IF( ntest < SIZE( cm, 2 ) ) THEN
         ntest = ngwt
      END IF
      !
      ! ... assign random values to wave functions
      !
      DO ib = noff, noff + n - 1
        pwt( : ) = 0.0d0
        DO ig = 3, ntest
          rranf1 = 0.5d0 - randy()
          rranf2 = randy()
          pwt( ig ) = ampre * CMPLX(rranf1, rranf2)
        END DO
        DO ig = 1, ngw
          cm( ig, ib ) = pwt( ig_l2g( ig ) )
        END DO
      END DO
      IF ( gzero ) THEN
        cm( 1, noff : noff + n - 1 ) = (0.0d0, 0.0d0)
      END IF

      DEALLOCATE( pwt )

      RETURN
    END SUBROUTINE wave_rand_init_x

!=----------------------------------------------------------------------------=!
   SUBROUTINE wave_sine_init_x( cm, n, noff ) !added:giovanni
!=----------------------------------------------------------------------------=!

      !  this routine sets the initial wavefunctions at random

! ... declare modules
      USE kinds,              ONLY: DP
      USE mp,                 ONLY: mp_sum
      USE mp_wave,            ONLY: splitwf
      USE reciprocal_vectors, ONLY: ngw, gx
      USE random_numbers,     ONLY: randy
      
      IMPLICIT NONE

      ! ... declare subroutine arguments 
      INTEGER,     INTENT(IN)  :: n, noff
      COMPLEX(DP), INTENT(OUT) :: cm(:,:)

      ! ... declare other variables
      INTEGER :: ntest, ig, ib
      REAL(DP) ::  rranf1, rranf2, ampre
      COMPLEX(DP), ALLOCATABLE :: pwt( : )

      ! ... Check array dimensions
      IF( SIZE( cm, 1 ) < ngw ) THEN 
        CALL errore(' wave_rand_init ', ' wrong dimensions ', 3)
      END IF

      ! ... Reset them to zero
 
      cm( :, noff : noff + n - 1 ) = 0.0d0

      ! ... initialize the wave functions in such a way that the values
      ! ... of the components are independent of the number of processors

      ampre = 0.01d0
      ALLOCATE( pwt( ngw ) )

      ntest = ngw
      IF( ntest < SIZE( cm, 2 ) ) THEN
         ntest = ngw
      END IF
      !
      ! ... assign real values to wave functions
      !
      DO ib = noff, noff + n - 1
        pwt( : ) = CMPLX(0.0d0,0.d0)
        DO ig = 1, ntest
          rranf1 = dcos(ib*(gx(1, ig )))+dcos(ib*(gx(2, ig )))+dcos(ib*(gx(3, ig )))
          rranf2 = (dsin(ib*(gx(1,ig)))+dsin(ib*(gx(2, ig)))+dsin(ib*(gx(3, ig))))
          pwt( ig ) = ampre * CMPLX(rranf1, rranf2)
        END DO
        DO ig = 1, ngw
          cm( ig, ib ) = pwt(  ig  )
        END DO
      END DO

      DEALLOCATE( pwt )

      RETURN
    END SUBROUTINE wave_sine_init_x
