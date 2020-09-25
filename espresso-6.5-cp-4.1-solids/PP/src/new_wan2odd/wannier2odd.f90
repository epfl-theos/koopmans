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
MODULE wannier2odd
  !---------------------------------------------------------------------
  !
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: wan2odd
  !
  INTEGER, PUBLIC :: iunwann = 8        ! unit for supercell WFs file
  INTEGER, PUBLIC :: nwordwfcx          ! nword for supercell wavefunctions
  INTEGER, PUBLIC :: nrtot              ! total num of R-vectors
  !
  CONTAINS
  !
  !---------------------------------------------------------------------
  SUBROUTINE wan2odd
    !-------------------------------------------------------------------
    !
    ! ...  This routine:
    !
    ! ...  1) reads the KS states u_nk(G) from PW and the rotation 
    ! ...     matrices U(k) from Wannier90
    !
    ! ...  2) Fourier transforms u_nk to real-space and extends them 
    ! ...     to the supercell defined by the k-points sampling
    !
    ! ...  3) applies the matrices U(k) and realizes the Wannier 
    ! ...     functions in G-space
    !
    ! ...  4) Wannier functions are finally written in a CP-readable
    ! ...     file
    !
    USE kinds,               ONLY : DP
    USE io_global,           ONLY : stdout, ionode
    USE io_files,            ONLY : nwordwfc, iunwfc, restart_dir
    USE io_base,             ONLY : write_rhog
    USE mp_pools,            ONLY : my_pool_id
    USE mp_bands,            ONLY : my_bgrp_id, root_bgrp_id, &
                                    root_bgrp, intra_bgrp_comm
    USE mp,                  ONLY : mp_sum
    USE wavefunctions,       ONLY : evc, psic
    USE fft_base,            ONLY : dffts
    USE fft_interfaces,      ONLY : invfft, fwfft
    USE buffers,             ONLY : open_buffer, close_buffer, &
                                    save_buffer, get_buffer
    USE lsda_mod,            ONLY : nspin
    USE klist,               ONLY : xk, ngk, igk_k, nkstot, nelec
    USE gvect,               ONLY : ngm
    USE wvfct,               ONLY : wg, npwx
    USE cell_base,           ONLY : tpiba, omega, at
    USE control_flags,       ONLY : gamma_only
    USE constants,           ONLY : eps6, tpi
    USE scf,                 ONLY : rho
    USE noncollin_module,    ONLY : npol
    USE wannier,             ONLY : iknum, ikstart, mp_grid, num_bands, n_wannier
    USE input_parameters,    ONLY : nkstot
    USE scell_wfc,           ONLY : extend_wfc
    USE read_wannier,        ONLY : read_u_matrices
    USE fft_supercell,       ONLY : dfftcp, setup_scell_fft, bg_cp, at_cp, omega_cp, &
                                    gamma_only_cp, npwxcp, ngmcp, gstart_cp, &
                                    mill_cp, ig_l2g_cp, check_fft, g_cp
    !
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256) :: dirname
    INTEGER :: ik, ikevc, ibnd, iw
    INTEGER :: i, j, k, ir, ipw
    INTEGER :: npw
    INTEGER :: iunwfcx = 24                       ! unit for supercell wfc file
    INTEGER :: io_level = 1
    LOGICAL :: exst
    COMPLEX(DP), ALLOCATABLE :: evcx(:,:)
    COMPLEX(DP), ALLOCATABLE :: wann(:,:), wann_r(:,:)
    COMPLEX(DP), ALLOCATABLE :: psicx(:)
    COMPLEX(DP), ALLOCATABLE :: rhor(:)           ! supercell density in real space
    COMPLEX(DP), ALLOCATABLE :: rhog(:,:)         ! supercell density in G-space
    COMPLEX(DP), ALLOCATABLE :: u_mat(:,:,:)      ! Wannier U-matrix
    COMPLEX(DP), ALLOCATABLE :: u_mat_opt(:,:,:)  ! Wannier optimal U-matrix
    REAL(DP) :: kvec(3), rvec(3), gvec(3)
    REAL(DP) :: nelec_, charge
    REAL(DP) :: dot_prod
    COMPLEX(DP) :: phase
    !
    !
    WRITE(*,*) 'RICCARDO nkstot = ', nkstot
    CALL start_clock( 'wannier2odd' )
    CALL setup_scell_fft
    !
    nwordwfcx = num_bands*npwxcp*npol
    ALLOCATE( evcx(npwxcp*npol,num_bands) )
    ALLOCATE( wann(npwxcp*npol,num_bands) )
    ALLOCATE( wann_r(npwxcp*npol,num_bands) )
    ALLOCATE( psicx(dfftcp%nnr) )
    ALLOCATE( rhor(dfftcp%nnr) )
    ALLOCATE( rhog(ngmcp,nspin) )
    WRITE(stdout,'(2x, "RICCARDO - check on parallelization - dffts%nnr:", i10)') dffts%nnr
    !
    rhor(:) = ( 0.D0, 0.D0 )
    rhog(:,:) = ( 0.D0, 0.D0 )
    !
    !
    ! ... open buffer for direct-access to the extended wavefunctions
    !
    CALL open_buffer( iunwfcx, 'wfcx', nwordwfcx, io_level, exst )
    !
    ! ... loop to read the primitive cell wavefunctions and 
    ! ... extend them to the supercell
    ! 
    DO ik = 1, iknum
      !
      ikevc = ik + ikstart - 1
      CALL davcio( evc, 2*nwordwfc, iunwfc, ikevc, -1 )
      npw = ngk(ik)
      kvec(:) = xk(:,ik)
      !
      DO ibnd = 1, num_bands     !!! SPIN TO CHECK!!!!
        !
        psic(:) = (0.D0, 0.D0)
        psicx(:) = (0.D0, 0.D0)
        psic( dffts%nl(igk_k(1:npw,ik)) ) = evc(1:npw, ibnd)
        IF( gamma_only ) psic( dffts%nlm(igk_k(1:npw,ik)) ) = CONJG(evc(1:npw, ibnd))
        CALL invfft( 'Wave', psic, dffts )
        !
        ! ... here we extend the wfc to the whole supercell
        !
        ! ... NB: the routine extend_wfc applies also the phase factor
        ! ...     e^(ikr) so the output wfc (psicx) is a proper Bloch
        ! ...     function and not just its periodic part
        !
        CALL extend_wfc( psic, psicx, dfftcp, kvec )
        !
        ! ... now we calculate the total density in the supercell
        !
        rhor(:) = rhor(:) + ( DBLE( psicx(:) )**2 + &
                             AIMAG( psicx(:) )**2 ) * wg(ibnd,ik) / omega
        !
        CALL fwfft( 'Wave', psicx, dfftcp )
        IF ( gamma_only_cp ) THEN
          evcx(1:npwxcp,ibnd) = psicx( dfftcp%nlm(1:npwxcp) )  ! THIS DOES NOT WORK!!! TO BE FIXED
        ELSE
          evcx(1:npwxcp,ibnd) = psicx( dfftcp%nl(1:npwxcp) )
        ENDIF
        !
      ENDDO ! ibnd
      !
      ! ... save the extended wavefunctions into the buffer
      !
      CALL save_buffer( evcx, nwordwfcx, iunwfcx, ik )
      !
    ENDDO ! ik
    !
    !
    CALL fwfft( 'Rho', rhor, dfftcp )
    IF ( gamma_only_cp ) THEN
      rhog(1:ngmcp,1) = rhor( dfftcp%nlm(1:ngmcp) )            ! THIS DOES NOT WORK!!! TO BE FIXED
    ELSE
      rhog(1:ngmcp,1) = rhor( dfftcp%nl(1:ngmcp) )
    ENDIF
    !
    ! ... check on the total charge
    !
    charge = 0.D0
    IF ( gstart_cp == 2 ) THEN
      charge = rhog(1,1) * omega_cp
    ENDIF
    CALL mp_sum( charge, intra_bgrp_comm )
    nelec_ = nelec * iknum
    IF ( check_fft ) nelec_ = nelec
    IF ( ABS( charge - nelec_ ) > 1.D-3 * charge ) &
         CALL errore( 'wan2odd', 'wrong total charge', 1 )
    !
    ! ... check on rho(G) when dfftcp is taken equal to dffts
    !
    IF ( check_fft ) THEN
      DO ik = 1, ngmcp
        IF ( ABS( DBLE(rhog(ik,1) - rho%of_g(ik,1)) ) .ge. eps6 .or. &
             ABS( AIMAG(rhog(ik,1) - rho%of_g(ik,1)) ) .ge. eps6 ) THEN
          CALL errore( 'wan2odd', 'rhog and rho%of_g differ', 1 )
        ENDIF
      ENDDO
    ENDIF
    !
    ! ... write G-space density to file
    !
    dirname = restart_dir()
    IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
         CALL write_rhog( TRIM(dirname) // "charge-density-x", &
         root_bgrp, intra_bgrp_comm, &
         bg_cp(:,1)*tpiba, bg_cp(:,2)*tpiba, bg_cp(:,3)*tpiba, &
         gamma_only_cp, mill_cp, ig_l2g_cp, rhog(:,:) )
    !
    !
    ! ... read the Wannier gauge matrices
    !
    ALLOCATE( u_mat(n_wannier,n_wannier,iknum) )
    ALLOCATE( u_mat_opt(n_wannier,n_wannier,iknum) )
    !
    CALL read_u_matrices( u_mat, u_mat_opt )
    !
    ! ... here the Wannier functions are realized
    !
    CALL open_buffer( iunwann, 'wann', nwordwfcx, io_level, exst )
    !
    ! ... first the Wannier functions in the reference cell 
    ! ... (R=0) are obtained
    !
    wann(:,:) = ( 0.D0, 0.D0 )
    !
    DO iw = 1, num_bands
      DO ik = 1, iknum
        !
        evcx(:,:) = ( 0.D0, 0.D0 )
        CALL get_buffer( evcx, nwordwfcx, iunwfcx, ik )
        !
        DO ibnd = 1, num_bands
          !
          wann(:,iw) = wann(:,iw) + u_mat(ibnd,iw,ik) * evcx(:,ibnd)
          !
        ENDDO
        !
      ENDDO
    ENDDO
    !    
    ! ... then the other Wannier functions are obtained by
    ! ... symmetry, w_Rn(g) = e^(igR) * w_0n(g)
    ! ... and saved to file
    !
    ir = 0
    !
    DO i = 1, mp_grid(1)
      DO j = 1, mp_grid(2)
        DO k = 1, mp_grid(3)
          !
          IF ( i==1 .and. j==1 .and. k==1 ) THEN
            wann_r = wann
            GOTO 100
          ENDIF
          !
          rvec(:) = (/ i-1, j-1, k-1 /)
          CALL cryst_to_cart( 1, rvec, at_cp, 1 )
          !
          DO ipw = 1, npwxcp
            !
            gvec(:) = g_cp(:,dfftcp%nl(ipw))
            dot_prod = tpi * SUM( rvec(:) * gvec )
            phase = CMPLX( COS(dot_prod), SIN(dot_prod), KIND=DP )
            !
            wann_r(ipw,:) = phase * wann(ipw,:)
            !
          ENDDO
          !
100       ir = ir + 1
          CALL save_buffer( wann_r, nwordwfcx, iunwann, ir )
          !
        ENDDO
      ENDDO
    ENDDO
    !
    nrtot = ir
    !
    ! ... write the WFs to a CP-Koopmans-readable file
    !
    !CALL write_wannier_cp( 'occupied' )
    !
    CALL close_buffer( iunwfcx, 'delete' )
    CALL close_buffer( iunwann, 'delete' )
    !
    !
    CALL stop_clock( 'wannier2odd' )
    !
    !
  END SUBROUTINE wan2odd
  !
  !
END MODULE wannier2odd
