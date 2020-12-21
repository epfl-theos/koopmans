
! Copyright (C) 2002-2007 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE h_psi( lda, n, m, psi, hpsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the product of the Hamiltonian
  ! ... matrix with m wavefunctions contained in psi
  !
  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  ! ...    hpsi  H*psi
  !
  USE kinds,    ONLY : DP
  USE bp,       ONLY : lelfield,l3dstring,gdir, efield, efield_cry
  USE becmod,   ONLY : becp, rbecp, becp_nc, calbec
  USE lsda_mod, ONLY : current_spin
  USE scf,      ONLY : vrs  
  USE wvfct,    ONLY : g2kin
  USE uspp,     ONLY : vkb, nkb
  USE ldaU,     ONLY : lda_plus_u
  USE gvect,    ONLY : gstart
  USE funct,    ONLY : dft_is_meta
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY: npol, noncolin
  USE realus,  ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                           bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, v_loc_psir, check_fft_orbital_gamma
  USE mp_global,            ONLY : nogrp
  USE control_flags,        ONLY : use_task_groups

#ifdef EXX
  USE exx,      ONLY : vexx
  USE funct,    ONLY : exx_is_active
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: lda, n, m
  COMPLEX(DP), INTENT(IN)  :: psi(lda*npol,m) 
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda*npol,m)   
  !
  INTEGER     :: ipol, ibnd, incr
  !
  CALL start_clock( 'h_psi' )
  !  
  ! ... Here we apply the kinetic energy (k+G)^2 psi
  !
  DO ibnd = 1, m
     IF ( noncolin ) THEN
        hpsi (1:n, ibnd) = g2kin (1:n) * psi (1:n, ibnd)
        hpsi (lda+1:lda+n, ibnd) = g2kin (1:n) * psi (lda+1:lda+n, ibnd)
     ELSE
        hpsi(1:n,ibnd) = g2kin(1:n) * psi(1:n,ibnd)
     END IF
  END DO
  !
  if (dft_is_meta()) call h_psi_meta (lda, n, m, psi, hpsi)
  !
  ! ... Here we add the Hubbard potential times psi
  !
  IF ( lda_plus_u ) CALL vhpsi( lda, n, m, psi, hpsi )
  !
  ! ... the local potential V_Loc psi
  !
  CALL start_clock( 'h_psi:vloc' )
  IF ( gamma_only ) THEN
     ! 
     IF (( use_task_groups ) .AND. ( m >= nogrp )) then 
      incr = 2 * nogrp
     else
      incr = 2
     endif

        IF (  real_space .and. nkb > 0  ) then !fixme: real_space without beta functions does not make sense
         do ibnd = 1 , m , incr
          !call check_fft_orbital_gamma(psi,ibnd,m)
          call fft_orbital_gamma(psi,ibnd,m,.true.) !transform the psi real space, saved in temporary memory
          call calbec_rs_gamma(ibnd,m,rbecp) !rbecp on psi
          call fft_orbital_gamma(hpsi,ibnd,m) ! psi is now replaced by hpsi
          call v_loc_psir(ibnd,m) ! hpsi -> hpsi + psi*vrs  (psi read from temporary memory)
          call add_vuspsir_gamma(ibnd,m) ! hpsi -> hpsi + vusp
          call bfft_orbital_gamma(hpsi,ibnd,m,.true.) !transform back hpsi, clear psi in temporary memory
         enddo
         !
        ELSE
         !not real space
         !CALL vloc_psi( lda, n, m, psi, vrs(1,current_spin), hpsi )
         CALL vloc_psi_gamma ( lda, n, m, psi, vrs(1,current_spin), hpsi ) 
        ENDIF 
     !
  ELSE IF ( noncolin ) THEN 
     !
     CALL vloc_psi_nc ( lda, n, m, psi, vrs, hpsi )
     !
  ELSE  
     !
     CALL vloc_psi_k ( lda, n, m, psi, vrs(1,current_spin), hpsi )
     !
  END IF  
  CALL stop_clock( 'h_psi:vloc' )
  !
  ! ... Here the product with the non local potential V_NL psi
  !
  IF ( nkb > 0 .and. .not. real_space) THEN !since the real space stuff has to be treated differently
     !
     CALL start_clock( 'h_psi:vnl' )
     IF ( gamma_only) THEN
        CALL calbec ( n, vkb, psi, rbecp, m )
        CALL add_vuspsi( lda, n, m, psi, hpsi )
     ELSE IF ( noncolin) THEN
        CALL calbec ( n, vkb, psi, becp_nc, m )
        CALL add_vuspsi_nc (lda, n, m, psi, hpsi )
     ELSE 
        CALL calbec( n, vkb, psi, becp, m )
        CALL add_vuspsi( lda, n, m, psi, hpsi )
     END IF
     CALL stop_clock( 'h_psi:vnl' )
     !
  END IF
#ifdef EXX
  IF ( exx_is_active() ) CALL vexx( lda, n, m, psi, hpsi )
#endif
  !
  ! ... electric enthalpy if required
  !
  IF ( lelfield ) THEN
     !
     IF ( .NOT.l3dstring ) THEN
        CALL h_epsi_her_apply( lda, n, m, psi, hpsi,gdir, efield )
     ELSE
        DO ipol=1,3
           CALL h_epsi_her_apply( lda, n, m, psi, hpsi,ipol,efield_cry(ipol) )
        END DO
     END IF
     !
  END IF
  !
  ! ... Gamma-only trick: set to zero the imaginary part of hpsi at G=0
  !
  IF ( gamma_only .AND. gstart == 2 ) &
      hpsi(1,1:m) = CMPLX( DBLE( hpsi(1,1:m) ), 0.D0 )
  !
  CALL stop_clock( 'h_psi' )
  !
  RETURN
  !
END SUBROUTINE h_psi
