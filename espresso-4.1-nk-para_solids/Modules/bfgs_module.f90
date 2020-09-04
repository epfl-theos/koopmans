!
! Copyright (C) 2003-2007 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE bfgs_module
   !----------------------------------------------------------------------------
   !
   ! ... Ionic relaxation through the Newton-Raphson optimization scheme
   ! ... based on the Broyden-Fletcher-Goldfarb-Shanno algorithm for the
   ! ... estimate of the inverse Hessian matrix.
   ! ... The ionic relaxation is performed in cartesian coordinates using
   ! ... a "trust radius" line search based on Wolfe conditions.
   !
   ! ... Written by Carlo Sbraccia ( 5/12/2003 )
   ! ... Maintained by Carlo Sbraccia ( 2003-2007 )
   ! ... Modified for variable-cell-shape relaxation by 
   ! ...   Javier Antonio Montoya and Stefano de Gironcoli (Dec 2007)
   !
   ! ... references :
   !
   ! ... 1) Roger Fletcher, Practical Methods of Optimization, John Wiley and
   ! ...    Sons, Chichester, 2nd edn, 1987.
   ! ... 2) Salomon R. Billeter, Alexander J. Turner, Walter Thiel,
   ! ...    Phys. Chem. Chem. Phys. 2, 2177 (2000).
   ! ... 3) Salomon R. Billeter, Alessandro Curioni, Wanda Andreoni,
   ! ...    Comput. Mat. Science 27, 437, (2003).
   ! ... 4) Ren Weiqing, PhD Thesis: Numerical Methods for the Study of Energy
   ! ...    Landscapes and Rare Events.
   !
   !
   USE kinds,     ONLY : DP
   USE io_files,  ONLY : iunbfgs, prefix
   USE constants, ONLY : eps16
   !
   USE basic_algebra_routines
   !
   IMPLICIT NONE
   !
   PRIVATE
   !
   ! ... public methods
   !
   PUBLIC :: bfgs, terminate_bfgs
   !
   ! ... public variables
   !
   PUBLIC :: bfgs_ndim,        &
             trust_radius_max, &
             trust_radius_min, &
             trust_radius_ini, &
             w_1,              &
             w_2
   !
   ! ... global variables
   !
   SAVE
   !
   CHARACTER (len=8) :: fname="energy" ! name of the function to be minimized
   !
   REAL(DP), ALLOCATABLE :: &
      pos(:),                &! positions + cell
      grad(:),              &! gradients + cell_force
      pos_p(:),           &! positions at the previous iteration
      grad_p(:),         &! gradients at the previous iteration
      inv_hess(:,:),   &! inverse hessian matrix ( updated using BFGS formula )
      metric(:,:),      &
      h_block(:,:),         &
      hinv_block(:,:),       &
      step(:),                &! the last bfgs step
      step_old(:),            &! old bfgs steps
      pos_old(:,:),           &! list of m old positions
      grad_old(:,:),          &! list of m old gradients
      pos_best(:)              ! best extrapolated positions
   REAL(DP) :: &
      trust_radius,             &! displacement along the bfgs direction
      trust_radius_old,         &! old displacement along the bfgs direction
      energy_p                   ! energy at the previous iteration
   INTEGER :: &
      scf_iter,              &! number of scf iterations
      bfgs_iter,           &! number of bfgs iterations
      gdiis_iter               ! number of gdiis iterations
   !
   LOGICAL :: &
      tr_min_hit               ! .TRUE. if the trust_radius has already been
                               !  set to the minimum value at the previous step
   !
   LOGICAL :: &
      conv_bfgs                ! .TRUE. when bfgs convergence has been achieved
   !
   ! ... default values for all these variables are set in
   ! ... Modules/read_namelist.f90 (SUBROUTINE ions_defaults)
   !
   INTEGER :: &
      bfgs_ndim                ! dimension of the subspace for GDIIS
                               ! fixed to 1 for standard BFGS algorithm
   REAL(DP)  :: &
      trust_radius_max,       &! maximum allowed displacement
      trust_radius_min,       &! minimum allowed displacement
      trust_radius_ini         ! initial displacement
   REAL(DP)  :: &
      w_1,                    &! parameters for Wolfe conditions
      w_2                      ! parameters for Wolfe conditions
   !
   ! ... Note that trust_radius_max, trust_radius_min, trust_radius_ini,
   ! ... w_1, w_2, bfgs_ndim have a default value, but can also be assigned
   ! ... in the input.
   !
CONTAINS
   !
   !------------------------------------------------------------------------
   SUBROUTINE bfgs( pos_in, h, energy, grad_in, fcell, fixion, scratch, stdout,&
                 energy_thr, grad_thr, cell_thr, energy_error, grad_error,     &
                 cell_error, istep, nstep, step_accepted, stop_bfgs, lmovecell )
      !------------------------------------------------------------------------
      !
      ! ... list of input/output arguments :
      !
      !  pos            : vector containing 3N coordinates of the system ( x )
      !  energy         : energy of the system ( V(x) )
      !  grad           : vector containing 3N components of grad( V(x) )
      !  fixion         : vector used to freeze a deg. of freedom
      !  scratch        : scratch directory
      !  stdout         : unit for standard output
      !  energy_thr     : treshold on energy difference for BFGS convergence
      !  grad_thr       : treshold on grad difference for BFGS convergence
      !                    the largest component of grad( V(x) ) is considered
      !  energy_error   : energy difference | V(x_i) - V(x_i-1) |
      !  grad_error     : the largest component of
      !                         | grad(V(x_i)) - grad(V(x_i-1)) |
      !  cell_error     : the largest component of: omega*(stress-press*I)
      !  nstep          : the maximun nuber of scf-steps
      !  step_accepted  : .TRUE. if a new BFGS step is done
      !  stop_bfgs      : .TRUE. if BFGS convergence has been achieved
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(INOUT) :: pos_in(:)
      REAL(DP),         INTENT(OUT)   :: h(3,3)
      REAL(DP),         INTENT(INOUT) :: energy
      REAL(DP),         INTENT(INOUT) :: grad_in(:)
      REAL(DP),         INTENT(INOUT) :: fcell(3,3)
      INTEGER,          INTENT(IN)    :: fixion(:)
      CHARACTER(LEN=*), INTENT(IN)    :: scratch
      INTEGER,          INTENT(IN)    :: stdout
      REAL(DP),         INTENT(IN)    :: energy_thr, grad_thr, cell_thr
      INTEGER,          INTENT(OUT)   :: istep
      INTEGER,          INTENT(IN)    :: nstep
      REAL(DP),         INTENT(OUT)   :: energy_error, grad_error, cell_error
      LOGICAL,          INTENT(OUT)   :: step_accepted, stop_bfgs
      LOGICAL,          INTENT(IN)    :: lmovecell
      !
      INTEGER  :: n, i, j, k, nat
      LOGICAL  :: lwolfe
      REAL(DP) :: dE0s, den
      ! ... for scaled coordinates
      REAL(DP) :: hinv(3,3),g(3,3),ginv(3,3),garbage
      !
      !
      n = SIZE( pos_in ) + 9
      nat = size (pos_in) / 3
      if (nat*3 /= size (pos_in)) call errore('bfgs',' strange dimension',1)
      !
      ! ... work-space allocation
      !
      ALLOCATE( pos(    n ) )
      ALLOCATE( grad(   n ) )
      !
      ALLOCATE( grad_old( n, bfgs_ndim ) )
      ALLOCATE( pos_old(  n, bfgs_ndim ) )
      !
      ALLOCATE( inv_hess( n, n ) )
      !
      ALLOCATE( pos_p(    n ) )
      ALLOCATE( grad_p(   n ) )
      ALLOCATE( step(     n ) )
      ALLOCATE( step_old( n ) )
      ALLOCATE( pos_best( n ) )
      ! ... scaled coordinates work-space
      ALLOCATE( hinv_block( n-9, n-9 ) )
      ! ... cell related work-space
      ALLOCATE( metric( n , n  ) )
      !
      ! ... the BFGS file read (pos & grad) in scaled coordinates
      !
      call invmat(3, h, hinv, garbage)
      hinv_block = 0.d0
      FORALL ( k=0:nat-1, i=1:3, j=1:3 ) hinv_block(i+3*k,j+3*k) = hinv(i,j)
      !
      ! ... generate metric to work with scaled ionic coordinates
      g = MATMUL(TRANSPOSE(h),h)
      call invmat(3,g,ginv,garbage)
      metric = 0.d0
      FORALL ( k=0:nat-1,   i=1:3, j=1:3 ) metric(i+3*k,j+3*k) = g(i,j)
      FORALL ( k=nat:nat+2, i=1:3, j=1:3 ) metric(i+3*k,j+3*k) = 10.0* ginv(i,j)
      !
      ! ... generate bfgs vectors for the degrees of freedom and their gradients
      pos = 0.0
      pos(1:n-9) = pos_in
      if (lmovecell) FORALL( i=1:3, j=1:3)  pos( n-9 + j+3*(i-1) ) = h(i,j)
      grad = 0.0
      grad(1:n-9) = grad_in
      if (lmovecell) FORALL( i=1:3, j=1:3) grad( n-9 + j+3*(i-1) ) = fcell(i,j)
      !
      ! if the cell moves the quantity to be minimized is the enthalpy
      IF ( lmovecell ) fname="enthalpy"
      !
      CALL read_bfgs_file( pos, grad, fixion, energy, scratch, n, stdout )
      !
      scf_iter = scf_iter + 1
      istep    = scf_iter
      !
      ! ... convergence is checked here
      !
      energy_error = ABS( energy_p - energy )
      grad_error = MAXVAL( ABS( MATMUL( TRANSPOSE(hinv_block), grad(1:n-9)) ) )
      conv_bfgs = energy_error < energy_thr
      conv_bfgs = conv_bfgs .AND. ( grad_error < grad_thr )
      !
      IF( lmovecell) THEN
          cell_error    = MAXVAL( ABS( grad(n-8:n) ) )
          conv_bfgs = conv_bfgs .AND. ( cell_error < cell_thr ) 
      END IF
      !
      stop_bfgs = conv_bfgs .OR. ( scf_iter >= nstep )
      !
      ! ... quick return if possible
      !
      IF ( stop_bfgs ) GOTO 1000
      !
      ! ... some output is written
      !
      WRITE( UNIT = stdout, &
           & FMT = '(/,5X,"number of scf cycles",T30,"= ",I3)' ) scf_iter
      WRITE( UNIT = stdout, &
           & FMT = '(5X,"number of bfgs steps",T30,"= ",I3,/)' ) bfgs_iter
      IF ( scf_iter > 1 ) WRITE( UNIT = stdout, &
           & FMT = '(5X,A," old",T30,"= ",F18.10," Ry")' ) fname,energy_p
      WRITE( UNIT = stdout, &
           & FMT = '(5X,A," new",T30,"= ",F18.10," Ry",/)' ) fname,energy
      !
      ! ... the bfgs algorithm starts here
      !
      IF ( ( energy > energy_p ) .AND. ( scf_iter > 1 ) ) THEN
         !
         ! ... the previous step is rejected, line search goes on
         !
         step_accepted = .FALSE.
         !
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"CASE: ",A,"_new > ",A,"_old",/)' ) fname,fname
         !
         ! ... the new trust radius is obtained by quadratic interpolation
         !
         ! ... E(s) = a*s*s + b*s + c      ( we use E(0), dE(0), E(s') )
         !
         ! ... s_min = - 0.5*( dE(0)*s'*s' ) / ( E(s') - E(0) - dE(0)*s' )
         !
         dE0s = ( grad_p(:) .dot. step_old(:) )
         !
         den = energy - energy_p - dE0s
         !
         IF ( den > eps16 ) THEN
            !
            trust_radius = - 0.5_DP*dE0s*trust_radius_old / den
            !
         ELSE
            !
            ! ... no quadratic interpolation is possible: we use bisection
            !
            trust_radius = 0.5_DP*trust_radius_old
            !
         END IF
         !
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr")' ) &
              trust_radius
         !
         ! ... values from the last succeseful bfgs step are restored
         !
         pos(:)  = pos_p(:)
         energy  = energy_p
         grad(:) = grad_p(:)
         !
         IF ( trust_radius < trust_radius_min ) THEN
            !
            ! ... the history is reset ( the history can be reset at most two
            ! ... consecutive times )
            !
            WRITE( UNIT = stdout, &
                   FMT = '(/,5X,"trust_radius < trust_radius_min")' )
            WRITE( UNIT = stdout, FMT = '(/,5X,"resetting bfgs history",/)' )
            !
            IF ( tr_min_hit ) THEN
               !
               ! ... the history has already been reset at the previous step :
               ! ... something is going wrong
               !
               CALL errore( 'bfgs', &
                            'bfgs history already reset at previous step', 1 )
               !
            END IF
            !
            CALL reset_bfgs( n )
            !
            step(:) = - ( inv_hess(:,:) .times. grad(:) )
            !
            trust_radius = trust_radius_min
            !
            tr_min_hit = .TRUE.
            !
         ELSE
            !
            ! ... old bfgs direction ( normalized ) is recovered
            !
            step(:) = step_old(:) / trust_radius_old
            !
            tr_min_hit = .FALSE.
            !
         END IF
         !
      ELSE
         !
         ! ... a new bfgs step is done
         !
         bfgs_iter = bfgs_iter + 1
         !
         IF ( bfgs_iter == 1 ) THEN
            !
            ! ... first iteration
            !
            IF ( grad_error < 0.01_DP ) &
               trust_radius_ini = MIN( 0.2_DP, trust_radius_ini )
            !
            step_accepted = .FALSE.
            !
         ELSE
            !
            step_accepted = .TRUE.
            !
            WRITE( UNIT = stdout, &
                 & FMT = '(5X,"CASE: ",A,"_new < ",A,"_old",/)' ) fname,fname
            !
            CALL check_wolfe_conditions( lwolfe, energy, grad )
            !
            CALL update_inverse_hessian( pos, grad, n, stdout )
            !
         END IF
         !
         IF ( bfgs_ndim > 1 ) THEN
            !
            ! ... GDIIS extrapolation
            !
            CALL gdiis_step()
            !
         ELSE
            !
            ! ... standard Newton-Raphson step
            !
            step(:) = - ( inv_hess(:,:) .times. grad(:) )
            !
         END IF
         !
         IF ( ( grad(:) .dot. step(:) ) > 0.0_DP ) THEN
            !
            WRITE( UNIT = stdout, &
                   FMT = '(5X,"uphill step: resetting bfgs history",/)' )
            !
            CALL reset_bfgs( n )
            !
            step(:) = - ( inv_hess(:,:) .times. grad(:) )
            !
         END IF
         !
         ! ... the new trust radius is computed
         !
         IF ( bfgs_iter == 1 ) THEN
            !
            trust_radius = trust_radius_ini
            !
            tr_min_hit = .FALSE.
            !
         ELSE
            !
            trust_radius = trust_radius_old
            !
            CALL compute_trust_radius( lwolfe, energy, grad, n, stdout )
            !
         END IF
         !
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr")' ) &
              trust_radius
         !
      END IF
      !
      ! ... step along the bfgs direction
      !
      IF ( scnorm( step(:) ) < eps16 ) &
         CALL errore( 'bfgs', 'NR step-length unreasonably short', 1 )
      !
      step(:) = trust_radius*step(:)/scnorm( step(:) )
      !
      ! ... information required by next iteration is saved here ( this must
      ! ... be done before positions are updated )
      !
      CALL write_bfgs_file( pos, energy, grad, scratch )
      !
      ! ... positions and cell are updated
      !
      pos(:) = pos(:) + step(:)
      !
1000  CONTINUE
      ! ... input ions+cell variables
      IF ( lmovecell ) FORALL( i=1:3, j=1:3) h(i,j) = pos( n-9 + j+3*(i-1) )
      pos_in = pos(1:n-9)
      ! ... update forces
      grad_in = grad(1:n-9)
      !
      ! ... work-space deallocation
      !
      DEALLOCATE( pos )
      DEALLOCATE( grad )
      DEALLOCATE( pos_p )
      DEALLOCATE( grad_p )
      DEALLOCATE( pos_old )
      DEALLOCATE( grad_old )
      DEALLOCATE( inv_hess )
      DEALLOCATE( step )
      DEALLOCATE( step_old )
      DEALLOCATE( pos_best )
      DEALLOCATE( hinv_block )
      DEALLOCATE( metric )
      !
      RETURN
      !
   CONTAINS
      !
      !--------------------------------------------------------------------
      SUBROUTINE gdiis_step()
         !--------------------------------------------------------------------
         USE basic_algebra_routines
         IMPLICIT NONE
         !
         REAL(DP), ALLOCATABLE :: res(:,:), overlap(:,:), work(:)
         INTEGER,  ALLOCATABLE :: iwork(:)
         INTEGER               :: k, k_m, info
         REAL(DP)              :: gamma0
         !
         !
         gdiis_iter = gdiis_iter + 1
         !
         k   = MIN( gdiis_iter, bfgs_ndim )
         k_m = k + 1
         !
         ALLOCATE( res( n, k ) )
         ALLOCATE( overlap( k_m, k_m ) )
         ALLOCATE( work( k_m ), iwork( k_m ) )
         !
         work(:)  = 0.0_DP
         iwork(:) = 0
         !
         ! ... the new direction is added to the workspace
         !
         DO i = bfgs_ndim, 2, -1
            !
            pos_old(:,i)  = pos_old(:,i-1)
            grad_old(:,i) = grad_old(:,i-1)
            !
         END DO
         !
         pos_old(:,1)  = pos(:)
         grad_old(:,1) = grad(:)
         !
         ! ... |res_i> = H^-1 \times |g_i>
         !
         CALL DGEMM( 'N', 'N', n, k, n, 1.0_DP, &
                     inv_hess, n, grad_old, n, 0.0_DP, res, n )
         !
         ! ... overlap_ij = <grad_i|res_j>
         !
         CALL DGEMM( 'T', 'N', k, k, n, 1.0_DP, &
                     res, n, res, n, 0.0_DP, overlap, k_m )
         !
         overlap( :, k_m) = 1.0_DP
         overlap(k_m, : ) = 1.0_DP
         overlap(k_m,k_m) = 0.0_DP
         !
         ! ... overlap is inverted via Bunch-Kaufman diagonal pivoting method
         !
         CALL DSYTRF( 'U', k_m, overlap, k_m, iwork, work, k_m, info )
         CALL DSYTRI( 'U', k_m, overlap, k_m, iwork, work, info )
         CALL errore( 'gdiis_step', 'error in Bunch-Kaufman inversion', info )
         !
         ! ... overlap is symmetrised
         !
         FORALL( i = 1:k_m, j = 1:k_m, j > i ) overlap(j,i) = overlap(i,j)
         !
         pos_best(:) = 0.0_DP
         step(:)     = 0.0_DP
         !
         DO i = 1, k
            !
            gamma0 = overlap(k_m,i)
            !
            pos_best(:) = pos_best(:) + gamma0*pos_old(:,i)
            !
            step(:) = step(:) - gamma0*res(:,i)
            !
         END DO
         !
         ! ... the step must be consistent with the last positions
         !
         step(:) = step(:) + ( pos_best(:) - pos(:) )
         !
         IF ( ( grad(:) .dot. step(:) ) > 0.0_DP ) THEN
            !
            ! ... if the extrapolated direction is uphill use only the
            ! ... last gradient and reset gdiis history
            !
            step(:) = - ( inv_hess(:,:) .times. grad(:) )
            !
            gdiis_iter = 0
            !
         END IF
         !
         DEALLOCATE( res, overlap, work, iwork )
         !
      END SUBROUTINE gdiis_step
      !
   END SUBROUTINE bfgs
   !
   !------------------------------------------------------------------------
   SUBROUTINE reset_bfgs( n )
      !------------------------------------------------------------------------
      ! ... inv_hess in re-initalized to the initial guess 
      ! ... defined as the inverse metric 
      !
      INTEGER, INTENT(IN) :: n
      !
      REAL(DP) :: garbage
      !
      call invmat(n, metric, inv_hess, garbage)
      !
      gdiis_iter = 0
      !
   END SUBROUTINE reset_bfgs
   !
   !------------------------------------------------------------------------
   SUBROUTINE read_bfgs_file( pos, grad, fixion, energy, scratch, n, stdout )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(INOUT) :: pos(:)
      REAL(DP),         INTENT(INOUT) :: grad(:)
      INTEGER,          INTENT(IN)    :: fixion(:)
      CHARACTER(LEN=*), INTENT(IN)    :: scratch
      INTEGER,          INTENT(IN)    :: n
      INTEGER,          INTENT(IN)    :: stdout
      REAL(DP),         INTENT(INOUT) :: energy
      !
      CHARACTER(LEN=256) :: bfgs_file
      LOGICAL            :: file_exists
      REAL(DP) :: garbage
      !
      !
      bfgs_file = TRIM( scratch ) // TRIM( prefix ) // '.bfgs'
      !
      INQUIRE( FILE = TRIM( bfgs_file ) , EXIST = file_exists )
      !
      IF ( file_exists ) THEN
         !
         ! ... bfgs is restarted from file
         !
         OPEN( UNIT = iunbfgs, FILE = TRIM( bfgs_file ), &
               STATUS = 'UNKNOWN', ACTION = 'READ' )
         !
         READ( iunbfgs, * ) pos_p
         READ( iunbfgs, * ) grad_p
         READ( iunbfgs, * ) scf_iter
         READ( iunbfgs, * ) bfgs_iter
         READ( iunbfgs, * ) gdiis_iter
         READ( iunbfgs, * ) energy_p
         READ( iunbfgs, * ) pos_old
         READ( iunbfgs, * ) grad_old
         READ( iunbfgs, * ) inv_hess
         READ( iunbfgs, * ) tr_min_hit
         !
         CLOSE( UNIT = iunbfgs )
         !
         trust_radius_old = scnorm( pos(:) - pos_p(:) )
         !
         step_old = ( pos(:) - pos_p(:) ) / trust_radius_old
         !
      ELSE
         !
         ! ... bfgs initialization
         !
         WRITE( UNIT = stdout, FMT = '(/,5X,"BFGS Geometry Optimization")' )
         !
         ! initialize the inv_hess to the inverse of the metric
         call invmat(n, metric, inv_hess, garbage)
         !
         pos_p      = 0.0_DP
         grad_p     = 0.0_DP
         scf_iter   = 0
         bfgs_iter  = 0
         gdiis_iter = 0
         energy_p   = energy
         step_old   = 0.0_DP
         !
         trust_radius_old = trust_radius_ini
         !
         pos_old(:,:)  = 0.0_DP
         grad_old(:,:) = 0.0_DP
         !
         tr_min_hit = .FALSE.
         !
      END IF
      !
   END SUBROUTINE read_bfgs_file
   !
   !------------------------------------------------------------------------
   SUBROUTINE write_bfgs_file( pos, energy, grad, scratch )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(IN) :: pos(:)
      REAL(DP),         INTENT(IN) :: energy
      REAL(DP),         INTENT(IN) :: grad(:)
      CHARACTER(LEN=*), INTENT(IN) :: scratch
      !
      !
      OPEN( UNIT = iunbfgs, FILE = TRIM( scratch )//TRIM( prefix )//'.bfgs', &
            STATUS = 'UNKNOWN', ACTION = 'WRITE' )
      !
      WRITE( iunbfgs, * ) pos
      WRITE( iunbfgs, * ) grad
      WRITE( iunbfgs, * ) scf_iter
      WRITE( iunbfgs, * ) bfgs_iter
      WRITE( iunbfgs, * ) gdiis_iter
      WRITE( iunbfgs, * ) energy
      WRITE( iunbfgs, * ) pos_old
      WRITE( iunbfgs, * ) grad_old
      WRITE( iunbfgs, * ) inv_hess
      WRITE( iunbfgs, * ) tr_min_hit
      !
      CLOSE( UNIT = iunbfgs )
      !
   END SUBROUTINE write_bfgs_file
   !
   !------------------------------------------------------------------------
   SUBROUTINE update_inverse_hessian( pos, grad, n, stdout )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN)  :: pos(:)
      REAL(DP), INTENT(IN)  :: grad(:)
      INTEGER,  INTENT(IN)  :: n
      INTEGER,  INTENT(IN)  :: stdout
      !
      REAL(DP), ALLOCATABLE :: y(:), s(:)
      REAL(DP), ALLOCATABLE :: Hy(:), yH(:)
      REAL(DP)              :: sdoty
      !
      ALLOCATE( y( n ), s( n ), Hy( n ), yH( n ) )
      !
      s(:) = pos(:)  - pos_p(:)
      y(:) = grad(:) - grad_p(:)
      !
      sdoty = ( s(:) .dot. y(:) )
      !
      IF ( ABS( sdoty ) < eps16 ) THEN
         !
         ! ... the history is reset
         !
         WRITE( stdout, '(/,5X,"WARNING: unexpected ", &
                         &     "behaviour in update_inverse_hessian")' )
         WRITE( stdout, '(  5X,"         resetting bfgs history",/)' )
         !
         CALL reset_bfgs( n )
         !
         RETURN
         !
      END IF
      !
      Hy(:) = ( inv_hess .times. y(:) )
      yH(:) = ( y(:) .times. inv_hess )
      !
      ! ... BFGS update
      !
      inv_hess = inv_hess + 1.0_DP / sdoty * &
                 ( ( 1.0_DP + ( y .dot. Hy ) / sdoty ) * matrix( s, s ) - &
                  ( matrix( s, yH ) +  matrix( Hy, s ) ) )
      !
      DEALLOCATE( y, s, Hy, yH )
      !
      RETURN
      !
   END SUBROUTINE update_inverse_hessian
   !
   !------------------------------------------------------------------------
   SUBROUTINE check_wolfe_conditions( lwolfe, energy, grad )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(DP), INTENT(IN)  :: energy
      REAL(DP), INTENT(IN)  :: grad(:)
      LOGICAL,  INTENT(OUT) :: lwolfe
      !
      !
      lwolfe = ( energy - energy_p ) < w_1 * ( grad_p .dot. step_old )
      !
      lwolfe = lwolfe .AND. &
               ABS( grad .dot. step_old ) > - w_2 * ( grad_p .dot. step_old )
      !
   END SUBROUTINE check_wolfe_conditions
   !
   !------------------------------------------------------------------------
   SUBROUTINE compute_trust_radius( lwolfe, energy, grad, n, stdout )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      LOGICAL,  INTENT(IN)  :: lwolfe
      REAL(DP), INTENT(IN)  :: energy
      REAL(DP), INTENT(IN)  :: grad(:)
      INTEGER,  INTENT(IN)  :: n
      INTEGER,  INTENT(IN)  :: stdout
      !
      REAL(DP) :: a
      LOGICAL  :: ltest
      !
      !
      ltest = ( energy - energy_p ) < w_1 * ( grad_p .dot. step_old )
      ltest = ltest .AND. ( scnorm( step ) > trust_radius_old )
      !
      IF ( ltest ) THEN
         !
         a = 1.5_DP
         !
      ELSE
         !
         a = 1.1_DP
         !
      END IF
      !
      IF ( lwolfe ) THEN
         !
         trust_radius = MIN( trust_radius_max, 2.0_DP*a*trust_radius_old )
         !
      ELSE
         !
         trust_radius = MIN( trust_radius_max, &
                             a*trust_radius_old, scnorm( step ) )
         !
      END IF
      !
      IF ( trust_radius < trust_radius_min ) THEN
         !
         ! ... the history is reset
         !
         IF ( tr_min_hit ) THEN
            !
            ! ... the history has already been reset at the previous step :
            ! ... something is going wrong
            !
            CALL errore( 'bfgs', 'history already reset at previous step', 1 )
            !
         END IF
         !
         WRITE( UNIT = stdout, &
                FMT = '(5X,"small trust_radius: resetting bfgs history",/)' )
         !
         CALL reset_bfgs( n )
         !
         step(:) = - ( inv_hess(:,:) .times. grad(:) )
         !
         trust_radius = trust_radius_min
         !
         tr_min_hit = .TRUE.
         !
      ELSE
         !
         tr_min_hit = .FALSE.
         !
      END IF
      !
   END SUBROUTINE compute_trust_radius
   !
   !----------------------------------------------------------------------- 
   REAL(DP) FUNCTION scnorm( vect )
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: vect(:)
      !
      scnorm = SQRT( DOT_PRODUCT( vect  ,  MATMUL( metric, vect ) ) )
      !
   END FUNCTION scnorm
   !
   !------------------------------------------------------------------------
   SUBROUTINE terminate_bfgs( energy, energy_thr, grad_thr, cell_thr, &
                              lmovecell, stdout, scratch )
      !------------------------------------------------------------------------
      !
      USE io_files, ONLY : prefix, delete_if_present
      !
      IMPLICIT NONE
      REAL(DP),         INTENT(IN) :: energy, energy_thr, grad_thr, cell_thr
      LOGICAL,          INTENT(IN) :: lmovecell
      INTEGER,          INTENT(IN) :: stdout
      CHARACTER(LEN=*), INTENT(IN) :: scratch
      !
      IF ( conv_bfgs ) THEN
         !
         WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"bfgs converged in ",I3," scf cycles and ", &
              &         I3," bfgs steps")' ) scf_iter, bfgs_iter
         IF ( lmovecell ) THEN
            WRITE( UNIT = stdout, &
              & FMT = '(5X,"(criteria: energy < ",E8.2,", force < ",E8.2, &
              &       ", cell < ",E8.2,")")') energy_thr, grad_thr, cell_thr
         ELSE
            WRITE( UNIT = stdout, &
              & FMT = '(5X,"(criteria: energy < ",E8.2,", force < ",E8.2, &
              &                        ")")') energy_thr, grad_thr
         END IF
         WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"End of BFGS Geometry Optimization")' )
         WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"Final ",A," = ",F18.10," Ry")' ) fname, energy
         !
         CALL delete_if_present( TRIM( scratch ) // TRIM( prefix ) // '.bfgs' )
         !
      ELSE
         !
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"The maximum number of steps has been reached.")' )
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"End of BFGS Geometry Optimization")' )
         !
      END IF
      !
   END SUBROUTINE terminate_bfgs
   !
END MODULE bfgs_module
