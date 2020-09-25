!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Written by Riccardo De Gennaro, EPFL (Sept 2020).
!
!
!-----------------------------------------------------------------------
MODULE scell_wfc
  !---------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  USE fft_support,         ONLY : good_fft_dimension
  USE fft_base,            ONLY : dffts
  !
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: extend_wfc
  !
  CONTAINS
  !
  !---------------------------------------------------------------------
  SUBROUTINE extend_wfc( psic, psicx, dfftx, kvec )
    !-------------------------------------------------------------------
    !
    ! ...  Here we extend in real space the periodic part of the Bloch
    ! ...  functions from the PW unit cell to the whole extension of 
    ! ...  BVK boundary conditions (supercell).
    !
    ! ...  psic  : input wfc defined on the unit cell
    ! ...  psicx : output wfc extended to the supercell
    ! ...  dfftx : fft descriptor of the supercell
    !
    USE io_global,           ONLY : stdout
    USE mp,                  ONLY : mp_sum
    USE mp_bands,            ONLY : intra_bgrp_comm
    USE fft_types,           ONLY : fft_type_descriptor
    USE cell_base,           ONLY : at
    USE constants,           ONLY : eps8, tpi
    !
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: psic(:)     ! input pcell wfc
    COMPLEX(DP), INTENT(OUT) :: psicx(:)   ! output supercell wfc
    TYPE( fft_type_descriptor ) :: dfftx   ! fft descriptor for the supercell
    REAL(DP), INTENT(IN) :: kvec(3)        ! k-vector for the phase factor 
    !
    COMPLEX(DP), ALLOCATABLE :: psicg(:)   ! global pcell wfc
    COMPLEX(DP) :: phase                   ! phase factor e^(ikr)
    INTEGER :: nnrg
    INTEGER :: i, j, k, ir, ir_end, idx, j0, k0
    INTEGER :: ip, jp, kp
    REAL(DP) :: r(3)
    REAL(DP) :: dot_prod
    !
    !
    ! ... broadcast the unit cell wfc to all the procs 
    !
    nnrg = dffts%nnr
    CALL mp_sum( nnrg, intra_bgrp_comm )
    ALLOCATE( psicg(nnrg) )
    CALL bcast_psic( psic, psicg )
    !
    !
#if defined (__MPI)
    j0 = dfftx%my_i0r2p
    k0 = dfftx%my_i0r3p
    IF( dfftx%nr1x == 0 ) dfftx%nr1x = good_fft_dimension( dfftx%nr1 )
    ir_end = MIN( dfftx%nnr, dfftx%nr1x*dfftx%my_nr2p*dfftx%my_nr3p )
#else
    j0 = 0
    k0 = 0
    ir_end = dfftx%nnr
#endif
    !
    !
    DO ir = 1, ir_end
      !
      ! ... three dimensional indexes
      !
      idx = ir - 1
      k = idx / ( dfftx%nr1x * dfftx%my_nr2p )
      idx = idx - ( dfftx%nr1x * dfftx%my_nr2p ) * k
      k = k + k0
      IF ( k .GE. dfftx%nr3 ) CYCLE
      j = idx / dfftx%nr1x
      idx = idx - dfftx%nr1x * j
      j = j + j0
      IF ( j .GE. dfftx%nr2 ) CYCLE
      i = idx
      IF ( i .GE. dfftx%nr1 ) CYCLE
      !
      ! ... ip, jp and kp represent the indexes folded into the 
      ! ... reference unit cell
      !
      ip = MOD( i, dffts%nr1 )
      jp = MOD( j, dffts%nr2 )
      kp = MOD( k, dffts%nr3 )
      !
      psicx(ir) = psicg( ip + jp*dffts%nr1x + kp*dffts%nr1x*dffts%my_nr2p + 1 )
      !
      !
      ! ... calculate the phase factor e^(ikr) and applies it to psicx
      !
      r(1) = DBLE(i) / dffts%nr1
      r(2) = DBLE(j) / dffts%nr2
      r(3) = DBLE(k) / dffts%nr3
      !
      CALL cryst_to_cart( 1, r, at, 1 )
      !
      dot_prod = tpi * SUM( kvec(:) * r(:) )
      phase = CMPLX( COS(dot_prod), SIN(dot_prod), KIND=DP )
      psicx(ir) = phase * psicx(ir)
      !
    ENDDO
    !
    !
  END SUBROUTINE extend_wfc
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE bcast_psic( psic, psicg )
    !-------------------------------------------------------------------
    !
    ! ...  This routine broadcasts the local wavefunction (psic)
    ! ...  to all the procs into a global variable psicg indexed
    ! ...  following the global ordering of the points in the
    ! ...  FFT grid
    !
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: psic(:)
    COMPLEX(DP), INTENT(OUT) :: psicg(:)
    !
    INTEGER :: i, j, k, ir, ir_end, idx, j0, k0, irg
    !
    !
#if defined (__MPI)
    j0 = dffts%my_i0r2p
    k0 = dffts%my_i0r3p
    IF( dffts%nr1x == 0 ) dffts%nr1x = good_fft_dimension( dffts%nr1 )
    ir_end = MIN( dffts%nnr, dffts%nr1x*dffts%my_nr2p*dffts%my_nr3p )
#else
    j0 = 0
    k0 = 0
    ir_end = dffts%nnr
#endif
    !
    !
    DO ir = 1, ir_end
      !
      ! ... three dimensional indexes
      !
      idx = ir - 1
      k = idx / ( dffts%nr1x * dffts%my_nr2p )
      idx = idx - ( dffts%nr1x * dffts%my_nr2p ) * k
      k = k + k0
      IF ( k .GE. dffts%nr3 ) CYCLE
      j = idx / dffts%nr1x
      idx = idx - dffts%nr1x * j
      j = j + j0
      IF ( j .GE. dffts%nr2 ) CYCLE
      i = idx
      IF ( i .GE. dffts%nr1 ) CYCLE
      !
      ! ... defining global index and saving psicg
      !
      irg = i + j*dffts%nr1x + k*dffts%nr1x*dffts%my_nr2p + 1
      psicg(irg) = psic(ir)
      !
    ENDDO    
    !
    !
  END SUBROUTINE bcast_psic
  !
  !
END MODULE scell_wfc
