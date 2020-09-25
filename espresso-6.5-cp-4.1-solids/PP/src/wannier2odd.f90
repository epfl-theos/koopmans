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
  USE kinds,               ONLY : DP
  !
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: wan2odd
  !
  CONTAINS
  !
  !---------------------------------------------------------------------
  SUBROUTINE wan2odd( seedname, ikstart )
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
    USE io_global,           ONLY : stdout, ionode
    USE io_files,            ONLY : nwordwfc, iunwfc, restart_dir
    USE io_base,             ONLY : write_rhog
    USE mp_pools,            ONLY : my_pool_id
    USE mp_bands,            ONLY : my_bgrp_id, root_bgrp_id, &
                                    root_bgrp, intra_bgrp_comm
    USE wavefunctions,       ONLY : evc, psic
    USE fft_base,            ONLY : dffts
    USE fft_interfaces,      ONLY : invfft, fwfft
    USE buffers,             ONLY : open_buffer, close_buffer, &
                                    save_buffer, get_buffer
    USE lsda_mod,            ONLY : nspin
    USE klist,               ONLY : xk, ngk, igk_k
    USE gvect,               ONLY : ngm
    USE wvfct,               ONLY : wg, npwx
    USE cell_base,           ONLY : tpiba, omega, at
    USE control_flags,       ONLY : gamma_only
    USE constants,           ONLY : tpi
    USE noncollin_module,    ONLY : npol
    USE scell_wfc,           ONLY : extend_wfc
    USE read_wannier,        ONLY : read_wannier_chk, num_bands, num_wann, num_kpts, & 
                                    kgrid, u_mat, u_mat_opt
    USE fft_supercell,       ONLY : dfftcp, setup_scell_fft, bg_cp, at_cp, &
                                    gamma_only_cp, npwxcp, ngmcp, &
                                    mill_cp, ig_l2g_cp, check_fft, g_cp
    USE cp_files,            ONLY : write_wannier_cp
    !
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), INTENT(IN) :: seedname
    INTEGER, INTENT(IN) :: ikstart
    !
    CHARACTER(LEN=256) :: dirname
    INTEGER :: ik, ikevc, ibnd, iw
    INTEGER :: i, j, k, ir, ipw
    INTEGER :: npw
    INTEGER :: nwordwfcx
    INTEGER :: iunwfcx = 24                       ! unit for supercell wfc file
    INTEGER :: io_level = 1
    LOGICAL :: exst
    COMPLEX(DP), ALLOCATABLE :: evcx(:,:)
    COMPLEX(DP), ALLOCATABLE :: ewan(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: psicx(:)
    COMPLEX(DP), ALLOCATABLE :: rhor(:)           ! supercell density in real space
    COMPLEX(DP), ALLOCATABLE :: rhog(:,:)         ! supercell density in G-space
!!! RICCARDO debug >>>
    COMPLEX(DP), ALLOCATABLE :: rhor_(:)
    COMPLEX(DP), ALLOCATABLE :: rhog_(:,:)
!!! RICCARDO debug <<<
    REAL(DP) :: kvec(3), rvec(3), gvec(3)
    REAL(DP) :: dot_prod
    COMPLEX(DP) :: phase
    !
    !
    CALL start_clock( 'wannier2odd' )
    !
    !
    ! ... initialize the supercell and read from Wannier90
    !
    CALL read_wannier_chk( seedname )
    CALL setup_scell_fft
    !
    !
    nwordwfcx = num_bands*npwxcp*npol
    ALLOCATE( evcx(npwxcp*npol,num_bands) )
    ALLOCATE( ewan(npwxcp*npol,num_bands,num_kpts) )
    ALLOCATE( psicx(dfftcp%nnr) )
    ALLOCATE( rhor(dfftcp%nnr) )
    ALLOCATE( rhog(ngmcp,nspin) )
    WRITE(stdout,'(2x, "RICCARDO - check on parallelization - dffts%nnr:", i10, /)') dffts%nnr
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
    DO ik = 1, num_kpts
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
    CALL check_rho( rhog )
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
    ! ... here the Wannier functions are realized:
    ! ... first the Wannier functions in the reference cell 
    ! ... (R=0) are obtained
    !
    ewan(:,:,:) = ( 0.D0, 0.D0 )
    !
!    DO iw = 1, num_bands
!      DO ik = 1, num_kpts
!        !
!        evcx(:,:) = ( 0.D0, 0.D0 )
!        CALL get_buffer( evcx, nwordwfcx, iunwfcx, ik )
!        !
!        DO ibnd = 1, num_bands
!          !
!          ewan(:,iw,1) = ewan(:,iw,1) + u_mat(ibnd,iw,ik) * evcx(:,ibnd) * SQRT( wg(ibnd,ik) )&
!                                                         / SQRT( DBLE(num_kpts) )
!          !
!        ENDDO
!        !
!      ENDDO
!    ENDDO
    ! 
!    CALL close_buffer( iunwfcx, 'delete' )
    !
    ! ... then the other Wannier functions are obtained by
    ! ... symmetry, w_Rn(g) = e^(igR) * w_0n(g)
    !
    ir = 0
    !
    DO i = 1, kgrid(1)
      DO j = 1, kgrid(2)
        DO k = 1, kgrid(3)
          !
          ir = ir + 1
!          IF ( i==1 .and. j==1 .and. k==1 ) CYCLE
          !
          rvec(:) = (/ i-1, j-1, k-1 /)
          CALL cryst_to_cart( 1, rvec, at, 1 )
          !
!!! RICCARDO debug >>>
          DO iw = 1, num_bands
            DO ik = 1, num_kpts
              !
              kvec(:) = xk(:,ik)
              dot_prod = tpi * SUM( kvec(:) * rvec(:) )
              phase = CMPLX( COS(dot_prod), -SIN(dot_prod), KIND=DP )
              evcx(:,:) = ( 0.D0, 0.D0 )
              CALL get_buffer( evcx, nwordwfcx, iunwfcx, ik )
              !
              DO ibnd = 1, num_bands
                !
                ewan(:,iw,ir) = ewan(:,iw,ir) + phase * u_mat(ibnd,iw,ik) * &
                           evcx(:,ibnd) * SQRT( wg(ibnd,ik) ) / SQRT( DBLE(num_kpts) )
                !
              ENDDO
              !
            ENDDO
          ENDDO
!!! RICCARDO debug <<<
!          DO ipw = 1, npwxcp
!            !
!            gvec(:) = g_cp(:,ipw) 
!            dot_prod = tpi * SUM( rvec(:) * gvec )
!            phase = CMPLX( COS(dot_prod), SIN(dot_prod), KIND=DP )
!!            phase = ( 1.0, 0.0 )
!            !
!            ewan(ipw,:,ir) = phase * ewan(ipw,:,1)
!            !
!          ENDDO
          !
          !
        ENDDO
      ENDDO
    ENDDO
    !
    CALL close_buffer( iunwfcx, 'delete' )
    !
!!! RICCARDO debug >>>
    ALLOCATE( rhor_(dfftcp%nnr) )
    ALLOCATE( rhog_(ngmcp,nspin) )
    rhor_(:) = ( 0.D0, 0.D0 )
    rhog_(:,:) = ( 0.D0, 0.D0 )
    do ik=1,num_kpts
      do ibnd=1,num_bands
        psicx = (0.d0,0.d0)
        psicx(dfftcp%nl(1:npwxcp)) = ewan(1:npwxcp,ibnd,ik)
        call invfft('Wave',psicx,dfftcp)
        rhor_(:) = rhor_(:) + ( DBLE( psicx(:) )**2 + &
                               AIMAG( psicx(:) )**2 ) / omega
      enddo
    enddo
    CALL fwfft( 'Rho', rhor_, dfftcp )
    IF ( gamma_only_cp ) THEN
      rhog_(1:ngmcp,1) = rhor_( dfftcp%nlm(1:ngmcp) )            ! THIS DOES NOT WORK!!! TO BE FIXED
    ELSE
      rhog_(1:ngmcp,1) = rhor_( dfftcp%nl(1:ngmcp) )
    ENDIF
    !
    CALL check_rho( rhog_, rhog )
    CALL plot_density( 'total', rhog_(:,1) )
!!! RICCARDO debug <<<
    !
    ! ... write the WFs to a CP-Koopmans-readable file
    !
    CALL plot_density( 'single', ewan(:,1,1), 1, 1 )
    CALL write_wannier_cp( 'occupied', ewan, npwxcp, num_bands, num_kpts, ig_l2g_cp )
    !
    !
    CALL stop_clock( 'wannier2odd' )
    !
    !
  END SUBROUTINE wan2odd
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE check_rho( rhog, rhogref )
    !-------------------------------------------------------------------
    !
    ! ...  this routine performs some checks on the total density:
    ! ...  1) the total charge
    ! ...  2) if rhoref is present, it checks if it matches with rho
    !
    USE mp,                  ONLY : mp_sum
    USE mp_bands,            ONLY : intra_bgrp_comm
    USE klist,               ONLY : nelec
    USE scf,                 ONLY : rho
    USE constants,           ONLY : eps6
    USE fft_supercell,       ONLY : gstart_cp, omega_cp, check_fft, ngmcp
    USE read_wannier,        ONLY : num_kpts
    !
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: rhog(:,:) 
    COMPLEX(DP), INTENT(IN), OPTIONAL :: rhogref(:,:)
    !
    REAL(DP) :: nelec_, charge
    INTEGER :: ik
    !
    !
    ! ... check the total charge
    !
    charge = 0.D0
    IF ( gstart_cp == 2 ) THEN
      charge = rhog(1,1) * omega_cp
    ENDIF
    !
    CALL mp_sum( charge, intra_bgrp_comm )
    nelec_ = nelec * num_kpts
    IF ( check_fft ) nelec_ = nelec
    IF ( ABS( charge - nelec_ ) > 1.D-3 * charge ) &
         CALL errore( 'wan2odd', 'wrong total charge', 1 )
    !
    !
    ! ... check rho(G) when dfftcp is taken equal to dffts
    !
    IF ( check_fft ) THEN
      DO ik = 1, ngmcp
        IF ( ABS( DBLE(rhog(ik,1) - rho%of_g(ik,1)) ) .ge. eps6 .or. &
             ABS( AIMAG(rhog(ik,1) - rho%of_g(ik,1)) ) .ge. eps6 ) THEN
          CALL errore( 'wan2odd', 'rhog and rho%of_g differ', ik )
        ENDIF
      ENDDO
    ENDIF
    !
    !
    ! ... when present, rhogref is compared to rhog
    !
    IF ( PRESENT(rhogref) ) THEN
      DO ik = 1, ngmcp
        IF ( ABS( DBLE(rhog(ik,1) - rhogref(ik,1)) ) .ge. eps6 .or. &
             ABS( AIMAG(rhog(ik,1) - rhogref(ik,1)) ) .ge. eps6 ) THEN
          CALL errore( 'wan2odd', 'rhog and rhogref differ', ik )
        ENDIF
      ENDDO 
    ENDIF
    !
    !
  END SUBROUTINE check_rho
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE plot_density( typ, evc, ir, ibnd )
    !-------------------------------------------------------------------
    !
    ! ...  This routine generates a XSF file, in a format readable
    ! ...  by XCrySDen, with the plot of the density of the input
    ! ...  (Wannier) function extended to the supercell
    !
    !      typ  : 'total' for plot of the total charge density
    !             'single' for plot of a single orbital density
    !
    !      evc  : input charge density ('total') or Wannier function ('single').
    !             These are expected in plane waves
    !
    !      ir   : R-vector index, expected only when typ='single',
    !             otherwise it is ignored
    !
    !      ibnd : band index, expected only when typ='single',
    !             otherwise it is ignored
    !
    USE io_global,           ONLY : ionode, stdout
    USE fft_interfaces,      ONLY : invfft
    USE cell_base,           ONLY : alat, omega, at
    USE constants,           ONLY : BOHR_RADIUS_ANGS
    USE ions_base,           ONLY : atm
    USE parameters,          ONLY : ntypx
    USE fft_supercell,       ONLY : dfftcp, at_cp, nat_cp, tau_cp, ityp_cp, &
                                    ngmcp, npwxcp
    !
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: typ         ! type of plot
    COMPLEX(DP), INTENT(IN) :: evc(:)           ! Wannier function to plot (plane waves)
    INTEGER, INTENT(IN), OPTIONAL :: ir         ! R-vector index of the WF
    INTEGER, INTENT(IN), OPTIONAL :: ibnd       ! band index of the WF
    !
    CHARACTER(LEN=30) :: filename
    INTEGER :: fileunit=224
    INTEGER :: i, rr
    REAL(DP) :: alang
    REAL(DP) :: orig(3), dirs(3,3)
    COMPLEX(DP) :: psic(dfftcp%nnr)
    REAL(DP) :: orb_dens(dfftcp%nnr)
    !
    !
    psic(:) = ( 0.D0, 0.D0 )
    !
    IF ( typ == 'total' ) THEN
      !
      IF ( ionode ) WRITE(stdout,'(5x, "Plot of the total charge density")')
      !
      WRITE( filename, 100 )
      psic(dfftcp%nl(1:ngmcp)) = evc(1:ngmcp)
      CALL invfft( 'Rho', psic, dfftcp )
      orb_dens(:) = psic(:)
      !
    ELSE IF ( typ == 'single' ) THEN
      !
      IF ( ( .not. PRESENT(ir) ) .or. ( .not. PRESENT(ibnd) ) ) &
        CALL errore( 'plot_density', 'Missing optional input parameter', 1 )
      !
      IF ( ionode ) THEN
        WRITE(stdout,'(5x, "Plot of the Wannier orbital density: ", 2I6)') ir, ibnd
      ENDIF
      !
      !
      WRITE( filename, 101 ) ir, ibnd
      psic(dfftcp%nl(1:npwxcp)) = evc(1:npwxcp) 
      CALL invfft( 'Wave', psic, dfftcp )
      orb_dens(:) = ( DBLE( psic(:) )**2 + AIMAG( psic(:) )**2 ) / omega
      !
    ELSE
      !
      CALL errore( 'plot_density', 'Invalid value for typ', 1 )
      !
    ENDIF
    !  
    !
    alang = alat * BOHR_RADIUS_ANGS
    !
    orig(:) = (/ 0.0, 0.0, 0.0 /)      ! origin of the datagrid
    dirs(:,1) = at_cp(:,1) * alat      ! 1st spanning vector datagrid
    dirs(:,2) = at_cp(:,2) * alat      ! 2nd spanning vector datagrid
    dirs(:,3) = at_cp(:,3) * alat      ! 3rd spanning vector datagrid
    !
    !
    IF ( ionode ) THEN
      !
      OPEN( UNIT=fileunit, FILE=trim(filename), STATUS='unknown', FORM='formatted' ) 
      !
      WRITE( fileunit, 201 ) at_cp(:,1)*alang, at_cp(:,2)*alang, at_cp(:,3)*alang
      WRITE( fileunit, 202 ) at_cp(:,1)*alang, at_cp(:,2)*alang, at_cp(:,3)*alang
      WRITE( fileunit, 203 ) nat_cp
      WRITE( fileunit, 204 ) ( atm(ityp_cp(i)), tau_cp(:,i)*alat*BOHR_RADIUS_ANGS, i=1,nat_cp )
      WRITE( fileunit, 205 )
      WRITE( fileunit, 206 ) dfftcp%nr1, dfftcp%nr2, dfftcp%nr3, &
                             orig, dirs(:,1), dirs(:,2), dirs(:,3) 
      WRITE( fileunit, 207 ) ( DBLE(psic(rr)), rr=1,dfftcp%nnr )
      WRITE( fileunit, 208 ) 
      !
      CLOSE( fileunit )
      !
    ENDIF
    !
    !
100 FORMAT( 'total_density.xsf' )    
101 FORMAT( 'WF_R', I4.4, '_B', I4.4, '.xsf' )    ! ex: WF_R1_B2.xsf
    !
201 FORMAT( 'CRYSTAL', /,'PRIMVEC', /, 3F12.7, /, 3F12.7, /, 3F12.7 )
202 FORMAT( 'CONVVEC', /, 3F12.7, /, 3F12.7, /, 3F12.7 )
203 FORMAT( 'PRIMCOORD', /, I6, '  1' )
204 FORMAT( A2, 3X, 3F12.7 )
205 FORMAT( /, 'BEGIN_BLOCK_DATAGRID_3D', /, '3D_field', /, 'BEGIN_DATAGRID_3D_UNKNOWN' )
206 FORMAT( 3I6, /, 3F12.6, /, 3F12.7, /, 3F12.7, /, 3F12.7 )
207 FORMAT( 6E13.5 )
208 FORMAT( 'END_DATAGRID_3D', /, 'END_BLOCK_DATAGRID_3D' )
    !
    !
  END SUBROUTINE plot_density
  !
  !
END MODULE wannier2odd
