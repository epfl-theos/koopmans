!
! Copyright (C) 2002-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!====================================================================
   SUBROUTINE inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ei1, ei2, ei3, &
                          sfac, c0, bec, firstiter, vpot, lgam )
!====================================================================
      !
      ! minimizes the total free energy with respect to the
      ! occupation matrix f_ij at fixed Kohn-Sham orbitals
      ! Cf. Marzari, Vanderbilt, Payne PRL 79, 1337 (1997)
      !

      ! declares modules
      USE kinds,          ONLY: dp
      USE core,           ONLY: nlcc_any
      USE energies,       ONLY: etot, enl, &
                                ekin, atot, entropy
      USE electrons_base, ONLY: f, nspin, iupdwn, nupdwn, nudx, &
                                nelt, nx => nbspx, n => nbsp, ispin 

      USE ensemble_dft,   ONLY: ninner, ismear, etemp, &
                                z0t, c0diag, becdiag, &
                                fmat0,  &
                                e0,   &
                                 compute_entropy2, &
                                compute_entropy_der, compute_entropy
      USE gvecp,          ONLY: ngm
      USE gvecs,          ONLY: ngs
      USE gvecb,          ONLY: ngb
      USE gvecw,          ONLY: ngw
      USE reciprocal_vectors, &
                          ONLY: ng0 => gstart
      USE ions_base,      ONLY: nat, nsp
      USE grid_dimensions, &
                          ONLY: nnr => nnrx, nr1, nr2, nr3
      USE smooth_grid_dimensions, &
                          ONLY: nnrsx
      USE io_global,      ONLY: io_global_start, ionode, &
                                ionode_id
      USE mp_global,      ONLY: intra_image_comm
      USE dener
      ! USE derho
      USE cdvan
      USE uspp,           ONLY: nhsa=> nkb, betae => vkb, &
                                rhovan => becsum
      USE cg_module,      ONLY:  ene_ok, &
                                enever
      USE ions_positions, ONLY: tau0
      USE mp,             ONLY: mp_sum,mp_bcast
      use cp_interfaces,  only: rhoofr, dforce
      USE twin_types !added:giovanni

      !
      IMPLICIT NONE

      ! declares local variables and counters
      INTEGER                :: nfi
      LOGICAL                :: tfirst 
      LOGICAL                :: tlast
      COMPLEX(DP)            :: eigr( ngw, nat )
      COMPLEX(DP)            :: c0( ngw, n )
      REAL(DP)               :: bec( nhsa, n )
      LOGICAL                :: firstiter


      INTEGER                :: irb( 3, nat )
      COMPLEX (DP)           :: eigrb( ngb, nat )
      REAL(DP)               :: rhor( nnr, nspin )
      REAL(DP)               :: vpot( nnr, nspin )
      COMPLEX(DP)            :: rhog( ngm, nspin )
      REAL(DP)               :: rhos( nnrsx, nspin )
      REAL(DP)               :: rhoc( nnr )
      COMPLEX(DP)            :: ei1( nr1:nr1, nat )
      COMPLEX(DP)            :: ei2( nr2:nr2, nat )
      COMPLEX(DP)            :: ei3( nr3:nr3, nat )
      COMPLEX(DP)            :: sfac( ngs, nsp )
      LOGICAL                 :: lgam

      INTEGER                :: i
      INTEGER                :: j
      INTEGER                :: ig
      INTEGER                :: k
      INTEGER                :: is
      REAL(DP)               :: entmp
      INTEGER                :: npt
      REAL(DP)               :: deltax
      REAL(DP)               :: deltaxmin
      REAL(DP)               :: xinit
      REAL(kind=DP), ALLOCATABLE :: dval(:), ex(:)!, ztmp(:,:)
      COMPLEX(DP), ALLOCATABLE :: dx_c(:)
      REAL(kind=DP), ALLOCATABLE :: fion2(:,:)
      type(twin_matrix), dimension(:), allocatable :: c0hc0, z1, zx, zxt, zaux
      COMPLEX(kind=DP), ALLOCATABLE :: h0c0(:,:)
      REAL(kind=DP), ALLOCATABLE :: f0(:),f1(:),fx(:), faux(:)
      type(twin_matrix), dimension(:), allocatable :: fmat1, fmatx, dfmat
!       REAL(kind=DP), ALLOCATABLE :: epsi0(:,:,:)
      type(twin_matrix), dimension(:), allocatable :: epsi0
      REAL(kind=DP) :: atot0,atot1,atotmin,etot0,etot1, ef1, enocc
      REAL(kind=DP) :: eqa,eqb,eqc, etot2, entropy2
      COMPLEX(DP) :: dentdx1_c, dedx1_c, dadx1_c
      REAL(kind=DP) :: f2,x,xmin
      INTEGER ::  niter,nss,istart,il

      CALL errore( " inner_loop ", " sub. not updated ", 1 )
   
      CALL start_clock( 'inner_loop' )
      ! initializes variables
      allocate(fion2(3,nat))
      allocate(h0c0(ngw,nx))
      allocate(dval(nx),ex(nx),dx_c(nx),f0(nx),f1(nx),fx(nx),faux(nx))
      allocate(z1(nspin),zx(nspin),zxt(nspin),zaux(nspin))
      allocate(fmat1(nspin),fmatx(nspin),dfmat(nspin))
      allocate(c0hc0(nspin), epsi0(nspin))
!       allocate(epsi0(nudx,nudx,nspin))
! 
      do is=1,nspin
	call init_twin(z1(is),lgam)
        call allocate_twin(z1(is),nudx,nudx,lgam)
	call init_twin(zx(is),lgam)
        call allocate_twin(zx(is),nudx,nudx,lgam)
	call init_twin(zxt(is),lgam)
        call allocate_twin(zxt(is),nudx,nudx,lgam)
	call init_twin(zaux(is),lgam)
        call allocate_twin(zaux(is),nudx,nudx,lgam)
	call init_twin(c0hc0(is),lgam)
        call allocate_twin(c0hc0(is),nudx,nudx,lgam)
	call init_twin(epsi0(is),lgam)
        call allocate_twin(epsi0(is),nudx,nudx,lgam)
	call init_twin(fmat1(is),lgam)
        call allocate_twin(fmat1(is),nudx,nudx,lgam)
	call init_twin(fmatx(is),lgam)
        call allocate_twin(fmatx(is),nudx,nudx,lgam)
	call init_twin(dfmat(is),lgam)
        call allocate_twin(dfmat(is),nudx,nudx,lgam)
      enddo
! 
      fion2( :, : )= 0.D0
      npt=10
      deltaxmin=1.D-8
      ! calculates the initial free energy if necessary
      IF( .not. ene_ok ) THEN
        ! calculates the overlaps bec between the wavefunctions c0
        ! and the beta functions
        CALL calbec( 1, nsp, eigr, c0, bec )
 
        ! rotates the wavefunctions c0 and the overlaps bec
        ! (the occupation matrix f_ij becomes diagonal f_i)      
        CALL rotate_twin( z0t, c0(:,:), bec, c0diag, becdiag, firstiter)
  
        ! calculates the electronic charge density
        CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, rhovan, &
                     rhor, rhog, rhos, enl, denl, ekin, dekin6 )
        IF(nlcc_any) CALL set_cc( irb, eigrb, rhoc )
  
        vpot = rhor

        ! calculates the SCF potential, the total energy
        ! and the ionic forces
        CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, &
                     tlast, ei1, ei2, ei3, irb, eigrb, sfac, &
                     tau0, fion2 )
        
        ! calculates the entropy
        CALL compute_entropy2( entropy, f, n, nspin )

      END IF
      ene_ok=.FALSE.
      atot0=etot+entropy
      etot0=etot

      ! calculates the occupation matrix 
      ! fmat_ij = \sum_k z_ki * f_k * z_kj
      CALL calcmt_twin( f, z0t, fmat0,firstiter )
   
      ! calculateas the energy contribution associated with 
      ! the augmentation charges and the 
      ! corresponding contribution to the ionic force
      CALL newd( vpot, irb, eigrb, rhovan, fion2 )
      CALL prefor( eigr, betae ) ! ATTENZIONE
  
      ! iterates on niter
      INNERLOOP : DO niter= 1, ninner

        ! operates the Hamiltonian on the wavefunction c0
        h0c0( :, : )= 0.D0
          DO i= 1, n, 2                      
            CALL dforce( i, bec, betae, c0, h0c0(:,i), h0c0(:,i+1), rhos, nnrsx, ispin, f, n, nspin )
          END DO
        
        ! calculates the Hamiltonian matrix in the basis {c0}           
        DO is=1,nspin
          call set_twin(c0hc0(is), CMPLX(0.d0,0.d0))
        ENDDO        

        DO is= 1, nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
            DO i= 1, nss
              DO k= 1, nss
                IF(.not.c0hc0(is)%iscmplx) THEN
		  DO ig= 1, ngw
		    c0hc0(is)%rvec(k,i)= c0hc0(is)%rvec(k,i) &
		      - 2.0d0*DBLE( CONJG( c0( ig,k+istart-1 ) ) &
		      * h0c0( ig, i+istart-1 ) )
		  END DO
		  IF( ng0 .eq. 2 ) THEN
		    c0hc0(is)%rvec(k,i)= c0hc0(is)%rvec(k,i) &
		      + DBLE( CONJG( c0( 1, k+istart-1 ) ) &
		      * h0c0( 1, i+istart-1 ) )
		  END IF
                ELSE
                  DO ig= 1, ngw
		    c0hc0(is)%cvec(k,i)= c0hc0(is)%cvec(k,i) &
		      - ( CONJG( c0( ig,k+istart-1 ) ) &
		      * h0c0( ig, i+istart-1 ) )
		  END DO
                ENDIF
              END DO
            END DO
            IF(.not.c0hc0(is)%iscmplx) THEN
	      CALL mp_sum( c0hc0(is)%rvec( 1:nss, 1:nss), intra_image_comm )
            ELSE
	      CALL mp_sum( c0hc0(is)%cvec( 1:nss, 1:nss), intra_image_comm )
            ENDIF
        END DO
 
        DO is= 1, nspin
          nss= nupdwn( is )
          call copy_twin(epsi0(is), c0hc0(is))
!           epsi0( 1:nss, 1:nss, is )= c0hc0( 1:nss, 1:nss, is ) ! ATTENZIONE
        END DO
          
        ! diagonalizes the Hamiltonian matrix
        !e0( : )= 0.D0 This is not set to 0 anymore ecause of blockocc
        DO is= 1, nspin
          istart= iupdwn( is )
          nss= nupdwn( is ) 
            IF( ionode ) THEN
              if(epsi0(is)%iscmplx) then
		CALL ddiag( nss, nss, epsi0(is)%rvec(1,1), dval(1), &
                          z1(is)%rvec(1,1), 1 )
              else
		CALL ddiag( nss, nss, epsi0(is)%cvec(1,1), dval(1), &
                          z1(is)%cvec(1,1), 1 )
              endif
            END IF
            CALL mp_bcast( dval, ionode_id, intra_image_comm )
            CALL twin_mp_bcast( z1(is), ionode_id)
            DO i= 1, nss
              e0( i+istart-1 )= dval( i )
            END DO
        END DO
  

        ! calculates the occupations and the fermi energy at the 
        ! end of the search direction
        CALL efermi( nelt, n, etemp, 1, f1, ef1, e0, enocc, ismear, & 
                     nspin )

        ! fmat1_ij = \sum_k z_ik * f_k * z_jk
        CALL calcm( f1, z1, fmat1, firstiter)
             
        ! calculates of dfmat_ij
        ! ( dfmat defines the search direction in occupation space)
        DO is= 1, nspin
          nss= nupdwn( is )
          if(.not.dfmat(is)%iscmplx) then
            dfmat(is)%rvec(1:nss,1:nss) = - fmat0(is)%rvec(1:nss,1:nss) &
                                       + fmat1(is)%rvec(1:nss,1:nss)
          else
            dfmat(is)%cvec(1:nss,1:nss) = - fmat0(is)%cvec(1:nss,1:nss) &
                                       + fmat1(is)%cvec(1:nss,1:nss)
          endif
        END DO
            
        ! 
        f0( 1:n )= f( 1:n )
       
        ! calculates fmatx= fmat0 + x* dfmat
        ! (here x=1)      
        x=1.D0
        DO is= 1, nspin
          nss= nupdwn( is )
          IF(.not.fmatx(is)%iscmplx) THEN
            fmatx(is)%rvec(1:nss,1:nss) = fmat0(is)%rvec(1:nss,1:nss) &
                                    + x * dfmat(is)%rvec( 1:nss, 1:nss)
          ELSE
            fmatx(is)%cvec(1:nss,1:nss) = fmat0(is)%cvec(1:nss,1:nss) &
                                    + x * dfmat(is)%cvec(1:nss,1:nss)
          ENDIF
        END DO
                      
        ! diagonalizes fmatx 
        fx( : ) = 0.0d0
        DO is=1, nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
            IF( ionode ) THEN
              if(lgam) then
		CALL ddiag( nss, nss, fmatx(is)%rvec(1,1), dval(1), &
                          zaux(is)%rvec(1,1), 1 )
              else
		CALL zdiag( nss, nss, fmatx(is)%cvec(1,1), dval(1), &
                          zaux(is)%cvec(1,1), 1 )
              endif
            END IF
            CALL mp_bcast( dval, ionode_id, intra_image_comm )
            CALL twin_mp_bcast( zaux(is), ionode_id)
            DO i= 1, nss
              faux( i+istart-1 )= dval( i )
            END DO
            DO i= 1, nss
              fx( i+istart-1 )= faux( nss-i+istart )
              IF(.not.zx(is)%iscmplx) then
		DO j=1, nss
		  zx(is)%rvec( i, j)= zaux(is)%rvec( i, nss-j+1 )
		END DO
              ELSE
		DO j=1, nss
		  zx(is)%cvec( i, j )= zaux(is)%cvec( i, nss-j+1)
		END DO
              ENDIF
            END DO
            IF(.not.zxt(is)%iscmplx) then
	      DO i= 1, nss
		DO k= 1, nss
		  zxt(is)%rvec( k, i)= zx(is)%rvec( i, k )
		END DO
	      END DO
            ELSE
	      DO i= 1, nss
		DO k= 1, nss
		  zxt(is)%cvec(k, i)= zx(is)%cvec(i, k)
		END DO
	      END DO
            ENDIF
        END DO

        ! updates f
        f( 1:n )= fx( 1:n )

        ! re-calculates fmatx
        CALL calcmt_twin( f, zxt, fmatx, firstiter)
              
        ! calculates the entropy and its derivatives with respect
        ! to each occupation at x
        CALL compute_entropy2( entropy, fx, n, nspin )
        CALL compute_entropy_der( ex, fx, n, nspin )
  
        ! calculates the free energy at x    
        CALL rotate( zxt, c0(:,:), bec, c0diag, becdiag, firstiter)
        CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, rhovan, &
                     rhor, rhog, rhos, enl, denl, ekin, dekin6 ) 
        IF(nlcc_any) CALL set_cc( irb, eigrb, rhoc )
        vpot = rhor
        CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
                     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion2 )
        CALL newd( vpot, irb, eigrb, rhovan, fion2 )
        CALL prefor( eigr, betae )
        atot1=etot+entropy
        etot1=etot
  
        ! calculates the Hamiltonian matrix
        h0c0( :, : )= 0.D0
          DO i= 1, n, 2
            CALL dforce( i, bec, betae, c0, h0c0(:,i), h0c0(:,i+1), rhos, nnrsx, ispin, f, n, nspin )
          END DO
        
        do is=1,nspin
	  call set_twin(c0hc0(is), CMPLX(0.d0,0.d0))
        enddo
!         c0hc0(:,:,:)=0.d0

        DO is= 1, nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
            DO i= 1, nss
              DO k= 1, nss
                IF(lgam) THEN
		  DO ig= 1, ngw
		    c0hc0(is)%rvec( k, i)= c0hc0(is)%rvec(k, i) &
		      - 2.0d0*DBLE( CONJG( c0( ig,k+istart-1 ) ) &
		      * h0c0( ig, i+istart-1 ) )
		  END DO
		  IF( ng0 .eq. 2 ) THEN
		    c0hc0(is)%rvec( k, i ) = c0hc0(is)%rvec(k, i) &
		      + DBLE(CONJG(c0( 1, k+istart-1)) &
		      * h0c0(1, i+istart-1))
		  END IF
                ELSE
		  DO ig= 1, ngw
		    c0hc0(is)%cvec( k, i)= c0hc0(is)%cvec(k, i) &
		      - ( CONJG( c0( ig,k+istart-1 ) ) &
		      * h0c0( ig, i+istart-1 ) )
		  END DO
                ENDIF
              END DO
            END DO
            CALL twin_mp_sum( c0hc0(is))
      ENDDO


        DO is= 1, nspin
          nss= nupdwn( is )
          call copy_twin(epsi0(is), c0hc0(is))
        END DO
                   
        !     calculates 
        !     (1) the energy derivative 
        !         dE / dx (1) = dE / df_ji * (f1_ji - f0_ji) 
        !     (2) the entropy derivative
        !         d Tr S(f) / df_ij = S'(f)_ji
        !         ( d( -T S )/dx is calculated as 
        !           (ex)_j [(zt)_ji (dfmat)_ik (zt)_jk] 
        !           instead of as
        !           [(zt)_jk (ex)_j (zt)_ji] (dfmat)_ik )
        !     (3) the free energy derivative
        dedx1_c= CMPLX(0.D0,0.d0)
        dentdx1_c= CMPLX(0.D0, 0.d0)
        DO is= 1,nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
            DO i= 1, nss
              dx_c( i+istart-1 )= CMPLX(0.D0,0.D0)
              IF(zxt(is)%iscmplx) THEN
		DO k= 1, nss
		  DO j= 1, nss
		    dx_c( i+istart-1 )= dx_c( i+istart-1 ) &
				    - zxt(is)%rvec(i,k) * fmat0(is)%rvec(k,j) &
				    * zxt(is)%rvec(i,j)
		  END DO
		END DO
              ELSE
		DO k= 1, nss
		  DO j= 1, nss
		    dx_c( i+istart-1 )= dx_c( i+istart-1 ) &
				    - zxt(is)%cvec(i,k) * fmat0(is)%cvec(k,j) &
				    * CONJG(zxt(is)%cvec(i,j))
		  END DO
		END DO
              ENDIF
              dx_c( i+istart-1 )= dx_c( i+istart-1 ) + fx( i+istart-1 )
            END DO
        END DO
        DO is= 1, nspin
          nss= nupdwn( is )
          istart= iupdwn( is )
            DO i= 1, nss
              dentdx1_c= dentdx1_c - etemp * dx_c( i+istart-1 ) &
                       * ex(i+istart-1)
              IF(.not.dfmat(is)%iscmplx) THEN
		DO k= 1, nss
		  dedx1_c= dedx1_c + dfmat(is)%rvec(i, k) * epsi0(is)%rvec(k, i)
		END DO
              ELSE
		DO k= 1, nss
		  dedx1_c= dedx1_c + dfmat(is)%cvec(i, k) * epsi0(is)%cvec(k, i)
		END DO
              ENDIF
            END DO
        END DO
        dadx1_c = dedx1_c + dentdx1_c
                   
        ! performs the line minimization
        ! (a) the energy contribution is approximated 
        !     by a second order polynomial
        ! (b) the entropic term is approximated by a function 
        !     of the form \sum_i s(a_i*x**2+b_i*x+c_i)
        !    (where s(f)=-f*ln(f)-(1-f)*ln(1-f) ).
        !    The coefficients a_i, b_i and c_i are calculated
        !    by first-order perturbation
        eqc= etot0
        eqa= dedx1_c - etot1 + etot0 !note that we expect dedx1_c to be real
        eqb= etot1 - etot0 - eqa
        atotmin= atot0
        xmin= 0.D0
        xinit=0.D0
        deltax= 1.D0 / DBLE( npt )
        DO WHILE ( deltax .gt. deltaxmin )
          SAMPLING : DO il= 0, npt
            x= xinit + deltax * DBLE( il )
            IF( x .gt. 1.D0 ) EXIT SAMPLING
            entropy2= 0.D0
            DO is=1,nspin
              nss= nupdwn( is )
              istart= iupdwn( is )
              DO i= 1, nss
                f2= fx( i+istart-1 ) + ( x-1 ) * dx_c( i+istart-1 ) &
                  + ( - fx( i+istart-1 ) + dx_c( i+istart-1 ) + &
                  f0( i+istart-1 ) ) * ( x-1 )**2
                CALL compute_entropy( entmp, f2, nspin )
                entropy2 = entropy2 + entmp
              END DO
            END DO
            etot2= eqa * x ** 2 + eqb * x + eqc
            IF ( ( etot2 + entropy2 ) .lt. atotmin ) THEN
              xmin= x
              atotmin= etot2+ entropy2
            END IF
          END DO SAMPLING
          xinit= MAX( 0.D0, xmin - deltax )
          deltax= 2.D0 * deltax / DBLE( npt )
        END DO

        IF( ionode ) THEN
          WRITE(37,'(a5,3f15.10)') &
                             'XMIN', xmin, atotmin, atotmin-atot0
          IF ( atotmin-atot0 .gt. 0.D0 ) &
          WRITE(37,*) "INNER LOOP, WARNING : increasing free energy"
        END IF
       
        !
        IF ( xmin .eq. 1.0D0 ) GOTO 300
  
        ! calculates the occupation matrix at xmin
        DO is= 1, nspin
          nss= nupdwn( is )
            IF(.not.fmat0(is)%iscmplx) THEN
	      DO i= 1, nss
		DO j= 1, nss
		  fmatx(is)%rvec( i, j )= fmat0(is)%rvec( i, j) &
				    + xmin * dfmat(is)%rvec( i, j)
		END DO
	      END DO
            ELSE
	      DO i= 1, nss
		DO j= 1, nss
		  fmatx(is)%cvec( i, j )= fmat0(is)%cvec( i, j) &
				    + xmin * dfmat(is)%cvec( i, j )
		END DO
	      END DO
            ENDIF
        END DO
               
        ! diagonalizes the occupation matrix at xmin
        fx( : )= 0.D0
        DO is= 1, nspin
          nss= nupdwn( is ) 
          istart= iupdwn( is )
! 
            IF(lgam) THEN
	      IF(ionode) CALL ddiag( nss, nss, fmatx(is)%rvec(1,1), &
                                dval(1), zaux(is)%rvec(1,1), 1 )  
            ELSE
	      IF(ionode) CALL ddiag( nss, nss, fmatx(is)%cvec(1,1), &
                                dval(1), zaux(is)%cvec(1,1), 1 )  
            ENDIF
! 
            CALL mp_bcast( dval, ionode_id, intra_image_comm )
            CALL twin_mp_bcast( zaux(is), ionode_id)
            DO i= 1, n
              faux( i+istart-1 )= dval( i )
            END DO
            DO i= 1, nss
              fx( i+istart-1 )= faux( nss-i+istart )
              IF(.not.zx(is)%iscmplx) THEN
		DO j= 1, nss
		  zx(is)%rvec( i, j ) = zaux(is)%rvec( i, nss-j+1 )
		END DO
              ELSE
		DO j= 1, nss
		  zx(is)%cvec( i, j) = zaux(is)%cvec( i, nss-j+1 )
		END DO
              ENDIF
            END DO
        END DO
    
        ! updates f
        f( 1:n )= fx( 1:n )
  
 300    CONTINUE
         
        ! updates z0t
        DO is= 1, nspin
          nss= nupdwn( is )
          if(.not.z0t(is)%iscmplx) THEN
	    DO i= 1, nss
	      DO k= 1, nss
		z0t(is)%rvec(k, i)= zx(is)%rvec(k, i)
	      END DO
	    END DO
          ELSE
	    DO i= 1, nss
	      DO k= 1, nss
		z0t(is)%cvec(k, i)= zx(is)%cvec(k, i)
	      END DO
	    END DO
          ENDIF
        END DO
        
        ! calculates the total free energy
        CALL calcmt_twin( f, z0t, fmat0, firstiter)
        CALL rotate( z0t, c0(:,:), bec, c0diag, becdiag, firstiter )
        CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, &
                     rhovan, rhor, rhog, rhos, enl, denl, ekin, dekin6 ) 
        IF(nlcc_any) CALL set_cc(irb,eigrb,rhoc)
        vpot = rhor
        CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
                     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion2 )
        CALL newd( vpot, irb, eigrb, rhovan, fion2 )
        CALL compute_entropy2( entropy, f, n, nspin )
        CALL prefor( eigr, betae )
        ene_ok= .TRUE.
        atotmin = etot + entropy      
        IF (ionode) write(37,'(a3,i2,2f15.10)') 'CI',niter,atot0,atotmin
        atot0=atotmin
        atot=atotmin
        etot0=etot
        enever=etot
        if(xmin==0.d0) exit 
     END DO INNERLOOP
! 
     do is=1,nspin
      call deallocate_twin(z1(is))
      call deallocate_twin(zx(is))
      call deallocate_twin(zxt(is))
      call deallocate_twin(zaux(is))
      call deallocate_twin(c0hc0(is))
      call deallocate_twin(epsi0(is))
      call deallocate_twin(fmat1(is))
      call deallocate_twin(fmatx(is))
      call deallocate_twin(dfmat(is))
     enddo

     deallocate(fion2,c0hc0,h0c0,z1)
     deallocate(zx,zxt,zaux,dval,ex,dx_c)
     deallocate(f0,f1,fx,faux)
     deallocate(dfmat,epsi0)
! 
     CALL stop_clock( 'inner_loop' )
     return
!====================================================================      
   END SUBROUTINE inner_loop
!====================================================================
