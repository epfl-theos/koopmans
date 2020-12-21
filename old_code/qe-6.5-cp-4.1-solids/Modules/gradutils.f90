!
! Copyright (C) 2018 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
! Various routines computing gradient and similar quantities via FFT
!--------------------------------------------------------------------
! FIXME: there is a dependency upon "cell_base" via variable tpiba
!        (2\pi/a) that maybe should be taken out from here?
!--------------------------------------------------------------------
SUBROUTINE external_gradient( a, grada )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing gradients of a real function in real space,
  ! to be called by an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : g
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: grada( 3, dfftp%nnr )

! A in real space, grad(A) in real space
  CALL fft_gradient_r2r( dfftp, a, g, grada )

  RETURN

END SUBROUTINE external_gradient
!----------------------------------------------------------------------------
SUBROUTINE fft_gradient_r2r( dfft, a, g, ga )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a
  ! ... input : dfft     FFT descriptor
  ! ...         a(:)     a real function on the real-space FFT grid
  ! ...         g(3,:)   G-vectors, in 2\pi/a units
  ! ... output: ga(3,:)  \grad a, real, on the real-space FFT grid
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba
  USE fft_interfaces,ONLY : fwfft, invfft 
  USE fft_types, ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  REAL(DP), INTENT(IN)  :: a(dfft%nnr), g(3,dfft%ngm)
  REAL(DP), INTENT(OUT) :: ga(3,dfft%nnr)
  !
  INTEGER  :: ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
  !
  ALLOCATE(  aux( dfft%nnr ) )
  ALLOCATE( gaux( dfft%nnr ) )
  !
  aux = CMPLX( a(:), 0.0_dp, kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Rho', aux, dfft)
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  DO ipol = 1, 3
     !
     gaux(:) = (0.0_dp, 0.0_dp)
     !
     gaux(dfft%nl(:)) = g(ipol,:) * CMPLX( -AIMAG( aux(dfft%nl(:)) ), &
                                             REAL( aux(dfft%nl(:)) ), kind=DP)
     !
     IF ( dfft%lgamma ) THEN
        !
        gaux(dfft%nlm(:)) = CMPLX(  REAL( gaux(dfft%nl(:)) ), &
                                  -AIMAG( gaux(dfft%nl(:)) ), kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Rho', gaux, dfft)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = tpiba * DBLE( gaux(:) )
     !
  END DO
  !
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE fft_gradient_r2r
!
!--------------------------------------------------------------------
SUBROUTINE fft_qgradient (dfft, a, xq, g, ga)
  !--------------------------------------------------------------------
  !
  ! Like fft_gradient_r2r, for complex arrays having a e^{iqr} behavior
  ! ... input : dfft     FFT descriptor
  ! ...         a(:)     a complex function on the real-space FFT grid
  ! ...         xq(3)    q-vector, in 2\pi/a units
  ! ...         g(3,:)   G-vectors, in 2\pi/a units
  ! ... output: ga(3,:)  \grad a, complex, on the real-space FFT grid
  !
  USE kinds,     ONLY: dp
  USE cell_base, ONLY: tpiba
  USE fft_types, ONLY : fft_type_descriptor
  USE fft_interfaces, ONLY: fwfft, invfft
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  !
  COMPLEX(DP), INTENT(IN)  :: a(dfft%nnr)
  REAL(DP), INTENT(IN):: xq(3), g(3,dfft%ngm)
  COMPLEX(DP), INTENT(OUT) :: ga(3,dfft%nnr)

  INTEGER  :: n, ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)

  ALLOCATE (gaux(dfft%nnr))
  ALLOCATE (aux (dfft%nnr))

  ! bring a(r) to G-space, a(G) ...
  aux (:) = a(:)

  CALL fwfft ('Rho', aux, dfft)
  ! multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
  DO ipol = 1, 3
     gaux (:) = (0.0_dp, 0.0_dp)
     DO n = 1, dfft%ngm
        gaux(dfft%nl(n)) = CMPLX( 0.0_dp, xq (ipol) + g(ipol,n), kind=DP ) * &
             aux (dfft%nl(n))
        IF ( dfft%lgamma ) gaux(dfft%nlm(n)) = CONJG( gaux (dfft%nl(n)) )
     END DO
     ! bring back to R-space, (\grad_ipol a)(r) ...

     CALL invfft ('Rho', gaux, dfft)
     ! ...and add the factor 2\pi/a  missing in the definition of q+G
     DO n = 1, dfft%nnr
        ga (ipol,n) = gaux (n) * tpiba
     END DO
  END DO

  DEALLOCATE (aux)
  DEALLOCATE (gaux)

  RETURN

END SUBROUTINE fft_qgradient
!
!----------------------------------------------------------------------------
SUBROUTINE fft_gradient_g2r( dfft, a, g, ga )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a - like fft_gradient with a(G) instead of a(r)
  ! ... input : dfft     FFT descriptor
  ! ...         a(:)     a(G), a complex function in G-space
  ! ...         g(3,:)   G-vectors, in 2\pi/a units
  ! ... output: ga(3,:)  \grad a, real, on the real-space FFT grid
  !
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE fft_interfaces,ONLY : invfft
  USE fft_types, ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  COMPLEX(DP), INTENT(IN)  :: a(dfft%ngm)
  REAL(DP),    INTENT(IN)  :: g(3,dfft%ngm)
  REAL(DP),    INTENT(OUT) :: ga(3,dfft%nnr)
  !
  INTEGER                  :: ipol, n
  COMPLEX(DP), ALLOCATABLE :: gaux(:)
  !
  !
  ALLOCATE( gaux( dfft%nnr ) )
  ga(:,:) = 0.D0
  !
  IF ( dfft%lgamma) THEN
     !
     ! ... Gamma tricks: perform 2 FFT's in a single shot
     ! x and y
     ipol = 1
     gaux(:) = (0.0_dp,0.0_dp)
     !
     ! ... multiply a(G) by iG to get the gradient in real space
     !
     DO n = 1, dfft%ngm
        gaux(dfft%nl (n)) = CMPLX( 0.0_dp, g(ipol,  n), kind=DP )* a(n) - &
                                           g(ipol+1,n) * a(n)
        gaux(dfft%nlm(n)) = CMPLX( 0.0_dp,-g(ipol,  n), kind=DP )*CONJG(a(n)) +&
                                            g(ipol+1,n) * CONJG(a(n))
     ENDDO
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Rho', gaux, dfft)
     !
     ! ... bring back to R-space, (\grad_ipol a)(r)
     ! ... add the factor 2\pi/a  missing in the definition of q+G
     !
     DO n = 1, dfft%nnr
        ga (ipol  , n) =  REAL( gaux(n) ) * tpiba
        ga (ipol+1, n) = AIMAG( gaux(n) ) * tpiba
     ENDDO
     ! z
     ipol = 3
     gaux(:) = (0.0_dp,0.0_dp)
     !
     ! ... multiply a(G) by iG to get the gradient in real space
     !
     gaux(dfft%nl (:)) = g(ipol,:) * CMPLX( -AIMAG(a(:)), REAL(a(:)), kind=DP)
     gaux(dfft%nlm(:)) = CONJG( gaux(dfft%nl(:)) )
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Rho', gaux, dfft)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = tpiba * REAL( gaux(:) )
     !
  ELSE
     !
     DO ipol = 1, 3
        !
        gaux(:) = (0.0_dp,0.0_dp)
        !
        ! ... multiply a(G) by iG to get the gradient in real space
        !
        gaux(dfft%nl(:)) = g(ipol,:) * CMPLX( -AIMAG(a(:)), REAL(a(:)), kind=DP)
        !
        ! ... bring back to R-space, (\grad_ipol a)(r) ...
        !
        CALL invfft ('Rho', gaux, dfft)
        !
        ! ...and add the factor 2\pi/a  missing in the definition of G
        !
        ga(ipol,:) = tpiba * REAL( gaux(:) )
        !
     END DO
     !
  END IF
  !
  DEALLOCATE( gaux )
  !
  RETURN
  !
END SUBROUTINE fft_gradient_g2r

!----------------------------------------------------------------------------
SUBROUTINE fft_graddot( dfft, a, g, da )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates da = \sum_i \grad_i a_i in R-space
  ! ... input : dfft     FFT descriptor
  ! ...         a(3,:)   a real function on the real-space FFT grid
  ! ...         g(3,:)   G-vectors, in 2\pi/a units
  ! ... output: ga(:)    \sum_i \grad_i a_i, real, on the real-space FFT grid
  !
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE fft_interfaces,ONLY : fwfft, invfft
  USE fft_types, ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  REAL(DP), INTENT(IN)     :: a(3,dfft%nnr), g(3,dfft%ngm)
  REAL(DP), INTENT(OUT)    :: da(dfft%nnr)
  !
  INTEGER                  :: n, ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
  COMPLEX(DP) :: fp, fm, aux1, aux2
  !
  ALLOCATE( aux(dfft%nnr) )
  ALLOCATE( gaux(dfft%nnr) )
  !
  gaux(:) = (0.0_dp,0.0_dp)
  !
  IF ( dfft%lgamma ) THEN
     !
     ! Gamma tricks: perform 2 FFT's in a single shot
     ! x and y
     ipol = 1
     aux(:) = CMPLX( a(ipol,:), a(ipol+1,:), kind=DP)
     !
     ! ... bring a(ipol,r) to G-space, a(G) ...
     !
     CALL fwfft ('Rho', aux, dfft)
     !
     ! ... multiply by iG to get the gradient in G-space
     !
     DO n = 1, dfft%ngm
        !
        fp = (aux(dfft%nl(n)) + aux (dfft%nlm(n)))*0.5_dp
        fm = (aux(dfft%nl(n)) - aux (dfft%nlm(n)))*0.5_dp
        aux1 = CMPLX( REAL(fp), AIMAG(fm), kind=DP)
        aux2 = CMPLX(AIMAG(fp), -REAL(fm), kind=DP)
        gaux (dfft%nl(n)) = &
             CMPLX(0.0_dp, g(ipol  ,n),kind=DP) * aux1 + &
             CMPLX(0.0_dp, g(ipol+1,n),kind=DP) * aux2
     ENDDO
     ! z
     ipol = 3
     aux(:) = CMPLX( a(ipol,:), 0.0_dp, kind=DP)
     !
     ! ... bring a(ipol,r) to G-space, a(G) ...
     !
     CALL fwfft ('Rho', aux, dfft)
     !
     ! ... multiply by iG to get the gradient in G-space
     ! ... fill both gaux(G) and gaux(-G) = gaux*(G)
     !
     DO n = 1, dfft%ngm
        gaux(dfft%nl(n)) = gaux(dfft%nl(n)) + g(ipol,n) * &
             CMPLX( -AIMAG( aux(dfft%nl(n)) ), &
                      REAL( aux(dfft%nl(n)) ), kind=DP)
        gaux(dfft%nlm(n)) = CONJG( gaux(dfft%nl(n)) )
     END DO
     !
  ELSE
     !
     DO ipol = 1, 3
        !
        aux = CMPLX( a(ipol,:), 0.0_dp, kind=DP)
        !
        ! ... bring a(ipol,r) to G-space, a(G) ...
        !
        CALL fwfft ('Rho', aux, dfft)
        !
        ! ... multiply by iG to get the gradient in G-space
        !
        DO n = 1, dfft%ngm
           gaux(dfft%nl(n)) = gaux(dfft%nl(n)) + g(ipol,n) * &
                CMPLX( -AIMAG( aux(dfft%nl(n)) ), &
                         REAL( aux(dfft%nl(n)) ), kind=DP)
        END DO
        !
     END DO
     !
  END IF
  !
  ! ... bring back to R-space, (\grad_ipol a)(r) ...
  !
  CALL invfft ('Rho', gaux, dfft)
  !
  ! ... add the factor 2\pi/a  missing in the definition of G and sum
  !
  da(:) = tpiba * REAL( gaux(:) )
  !
  DEALLOCATE( aux, gaux )
  !
  RETURN
  !
END SUBROUTINE fft_graddot

!--------------------------------------------------------------------
SUBROUTINE fft_qgraddot ( dfft, a, xq, g, da)
  !--------------------------------------------------------------------
  !
  ! Like fft_graddot, for complex arrays having a e^{iqr} dependency
  ! ... input : dfft     FFT descriptor
  ! ...         a(3,:)   a complex function on the real-space FFT grid
  ! ...         xq(3)    q-vector, in 2\pi/a units
  ! ...         g(3,:)   G-vectors, in 2\pi/a units
  ! ... output: ga(:)    \sum_i \grad_i a_i, complex, on the real-space FFT grid
  !
  USE kinds,          ONLY : DP
  USE cell_base,      ONLY : tpiba
  USE fft_interfaces, ONLY : fwfft, invfft
  USE fft_types, ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  COMPLEX(DP), INTENT(IN)  :: a(3,dfft%nnr)
  REAL(DP), INTENT(IN)     :: xq(3), g(3,dfft%ngm)
  COMPLEX(DP), INTENT(OUT) :: da(dfft%nnr)
  
  INTEGER :: n, ipol
  COMPLEX(DP), allocatable :: aux (:)

  ALLOCATE (aux (dfft%nnr))
  da(:) = (0.0_dp, 0.0_dp)
  DO ipol = 1, 3
     ! copy a(ipol,r) to a complex array...
     DO n = 1, dfft%nnr
        aux (n) = a (ipol, n)
     END DO
     ! bring a(ipol,r) to G-space, a(G) ...
     CALL fwfft ('Rho', aux, dfft)
     ! multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
     DO n = 1, dfft%ngm
        da (dfft%nl(n)) = da (dfft%nl(n)) + &
             CMPLX(0.0_dp, xq (ipol) + g (ipol, n),kind=DP) * aux(dfft%nl(n))
     END DO
  END DO
  IF ( dfft%lgamma ) THEN
     DO n = 1, dfft%ngm
        da (dfft%nlm(n)) = CONJG( da (dfft%nl(n)) )
     END DO
  END IF
  !  bring back to R-space, (\grad_ipol a)(r) ...
  CALL invfft ('Rho', da, dfft)
  ! ...add the factor 2\pi/a  missing in the definition of q+G
  da (:) = da (:) * tpiba

  DEALLOCATE(aux)

  RETURN
 
END SUBROUTINE fft_qgraddot

!--------------------------------------------------------------------
! Routines computing laplacian via FFT
!--------------------------------------------------------------------

!--------------------------------------------------------------------
SUBROUTINE external_laplacian( a, lapla )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing laplacian in real space, to be called by 
  ! an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : gg
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: lapla( dfftp%nnr )

! A in real space, lapl(A) in real space
  CALL fft_laplacian( dfftp, a, gg, lapla )

  RETURN

END SUBROUTINE external_laplacian

!--------------------------------------------------------------------
SUBROUTINE fft_laplacian( dfft, a, gg, lapla )
!--------------------------------------------------------------------
  !
  ! ... Calculates lapla = laplacian(a)
  ! ... input : dfft     FFT descriptor
  ! ...         a(:)     a real function on the real-space FFT grid
  ! ...         gg(:)    square modules of G-vectors, in (2\pi/a)^2 units
  ! ... output: lapla(:) \nabla^2 a, real, on the real-space FFT grid
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba2
  USE fft_types, ONLY : fft_type_descriptor
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  REAL(DP), INTENT(IN)  :: a(dfft%nnr), gg(dfft%ngm)
  REAL(DP), INTENT(OUT) :: lapla(dfft%nnr)
  !
  INTEGER                  :: ig
  COMPLEX(DP), ALLOCATABLE :: aux(:), laux(:)
  !
  !
  ALLOCATE(  aux( dfft%nnr ) )
  ALLOCATE( laux( dfft%nnr ) )
  !
  aux = CMPLX( a(:), 0.0_dp, kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Rho', aux, dfft)
  !
  ! ... Compute the laplacian
  !
  laux(:) = (0.0_dp, 0.0_dp)
  !
  DO ig = 1, dfft%ngm
     !
     laux(dfft%nl(ig)) = -gg(ig)*aux(dfft%nl(ig))
     !
  END DO
  !
  IF ( dfft%lgamma ) THEN
     !
     laux(dfft%nlm(:)) = CMPLX( REAL(laux(dfft%nl(:)) ), &
                              -AIMAG(laux(dfft%nl(:)) ), kind=DP)
     !
  ENDIF
  !
  ! ... bring back to R-space, (\lapl a)(r) ...
  !
  CALL invfft ('Rho', laux, dfft)
  !
  ! ... add the missing factor (2\pi/a)^2 in G
  !
  lapla = tpiba2 * REAL( laux )
  !
  DEALLOCATE( laux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE fft_laplacian
!
!--------------------------------------------------------------------
! Routines computing hessian via FFT
!--------------------------------------------------------------------
!
!----------------------------------------------------------------------
SUBROUTINE fft_hessian_g2r ( dfft, a, g, ha )
!----------------------------------------------------------------------
  !
  ! ... Calculates ha = hessian(a)
  ! ... input : dfft     FFT descriptor
  ! ...         a(:)     a real function on the real-space FFT grid
  ! ...         g(3,:)   G-vectors, in (2\pi/a)^2 units
  ! ... output: ha(6,:)  hessian(a), real, on the real-space FFT grid
  ! ...                  lower-packed matrix indeces 1-6 correspond to:
  ! ...                  1 = xx, 2 = yx, 3 = yy, 4 = zx, 5 = zy, 6 = zz
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba
  USE fft_types, ONLY : fft_type_descriptor
  USE fft_interfaces,ONLY : fwfft, invfft
  USE fft_helper_subroutines, ONLY: fftx_oned2threed
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  REAL(DP), INTENT(IN)  :: g(3,dfft%ngm)
  COMPLEX(DP), INTENT(IN)  :: a(dfft%ngm)
  REAL(DP), INTENT(OUT) :: ha( 6, dfft%nnr )
  !
  INTEGER                  :: ig, ir
  COMPLEX(DP), ALLOCATABLE :: aux(:), haux(:,:)
  !
  IF ( .NOT. dfft%lgamma ) CALL errore ('fft_hessian_g2r',&
       'only gamma case is implemented',1)
  ALLOCATE ( aux(dfft%nnr))
  ALLOCATE (haux(dfft%ngm,2))
  ! xx, yx
  DO ig=1,dfft%ngm
     haux(ig,1) = -tpiba**2*g(1,ig)**2     *a(ig)
     haux(ig,2) = -tpiba**2*g(1,ig)*g(2,ig)*a(ig)
  END DO
  CALL fftx_oned2threed( dfft, aux, haux(:,1), haux(:,2) )
  CALL invfft('Rho', aux, dfft)
  DO ir=1,dfft%nnr
     ha(1,ir) = DBLE(aux(ir))
     ha(2,ir) =AIMAG(aux(ir))
  END DO
  ! yy, zx
  DO ig=1,dfft%ngm
     haux(ig,1) = -tpiba**2*g(2,ig)**2     *a(ig)
     haux(ig,2) = -tpiba**2*g(1,ig)*g(3,ig)*a(ig)
  END DO
  CALL fftx_oned2threed( dfft, aux, haux(:,1), haux(:,2) )
  CALL invfft('Rho', aux, dfft)
  DO ir=1,dfft%nnr
     ha(3,ir) = DBLE(aux(ir))
     ha(4,ir) =AIMAG(aux(ir))
  END DO
  ! zy, zz
  DO ig=1,dfft%ngm
     haux(ig,1) = -tpiba**2*g(2,ig)*g(3,ig)*a(ig)
     haux(ig,2) = -tpiba**2*g(3,ig)**2     *a(ig)
  END DO
  CALL fftx_oned2threed( dfft, aux, haux(:,1), haux(:,2) )
  CALL invfft('Rho', aux, dfft)
  DO ir=1,dfft%nnr
     ha(5,ir) = DBLE(aux(ir))
     ha(6,ir) =AIMAG(aux(ir))
  END DO
  !
  DEALLOCATE(aux)
  DEALLOCATE(haux)
  
END SUBROUTINE fft_hessian_g2r
!--------------------------------------------------------------------
SUBROUTINE fft_hessian( dfft, a, g, ga, ha )
!--------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a and ha = hessian(a)
  ! ... input : dfft     FFT descriptor
  ! ...         a(:)     a real function on the real-space FFT grid
  ! ...         g(3,:)   G-vectors, in (2\pi/a)^2 units
  ! ... output: ga(3,:)  \grad a, real, on the real-space FFT grid
  ! ...         ha(3,3,:)  hessian(a), real, on the real-space FFT grid
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba
  USE fft_types, ONLY : fft_type_descriptor
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  REAL(DP), INTENT(IN)  :: a(dfft%nnr), g(3,dfft%ngm)
  REAL(DP), INTENT(OUT) :: ga( 3, dfft%nnr )
  REAL(DP), INTENT(OUT) :: ha( 3, 3, dfft%nnr )
  !
  INTEGER                  :: ipol, jpol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:), haux(:)
  !
  !
  ALLOCATE(  aux( dfft%nnr ) )
  ALLOCATE( gaux( dfft%nnr ) )
  ALLOCATE( haux( dfft%nnr ) )
  !
  aux = CMPLX( a(:), 0.0_dp, kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Rho', aux, dfft)
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  DO ipol = 1, 3
     !
     gaux(:) = (0.0_dp,0.0_dp)
     !
     gaux(dfft%nl(:)) = g(ipol,:) * CMPLX( -AIMAG( aux(dfft%nl(:)) ), &
                                             REAL( aux(dfft%nl(:)) ), kind=DP )
     !
     IF ( dfft%lgamma ) THEN
        !
        gaux(dfft%nlm(:)) = CMPLX(  REAL( gaux(dfft%nl(:)) ), &
                                  -AIMAG( gaux(dfft%nl(:)) ), kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Rho', gaux, dfft)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = tpiba * REAL( gaux(:) )
     !
     ! ... compute the second derivatives
     !
     DO jpol = 1, ipol
        !
        haux(:) = (0.0_dp,0.0_dp)
        !
        haux(dfft%nl(:)) = - g(ipol,:) * g(jpol,:) * &
             CMPLX( REAL( aux(dfft%nl(:)) ), &
                   AIMAG( aux(dfft%nl(:)) ), kind=DP)
        !
        IF ( dfft%lgamma ) THEN
           !
           haux(dfft%nlm(:)) = CMPLX(  REAL( haux(dfft%nl(:)) ), &
                                     -AIMAG( haux(dfft%nl(:)) ), kind=DP)
           !
        END IF
        !
        ! ... bring back to R-space, (\grad_ipol a)(r) ...
        !
        CALL invfft ('Rho', haux, dfft)
        !
        ! ...and add the factor 2\pi/a  missing in the definition of G
        !
        ha(ipol, jpol, :) = tpiba * tpiba * REAL( haux(:) )
        !
        ha(jpol, ipol, :) = ha(ipol, jpol, :) 
        !
     END DO
     !
  END DO
  !
  DEALLOCATE( haux )
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE fft_hessian
!--------------------------------------------------------------------
SUBROUTINE external_hessian( a, grada, hessa )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing hessian in real space, to be called by
  ! an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : g
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: grada( 3, dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: hessa( 3, 3, dfftp%nnr )

! A in real space, grad(A) and hess(A) in real space
  CALL fft_hessian( dfftp, a, g, grada, hessa )

  RETURN

END SUBROUTINE external_hessian

!--------------------------------------------------------------------
