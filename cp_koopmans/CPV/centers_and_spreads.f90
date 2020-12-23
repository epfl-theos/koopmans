!
! Copyright (C) 2002-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE centers_and_spreads
  !---------------------------------------------------------------------
  !
  ! ... This contains all the routines and variables important for the
  ! ... calculation of centers and spreads of the variational orbitals
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE reciprocal_vectors,   ONLY : g2_g, g, gx
  USE constants,            ONLY : BOHR_RADIUS_ANGS
  !
  !
  IMPLICIT NONE
  !
  SAVE
  !
  REAL(DP), ALLOCATABLE :: centers_occ(:,:), centers_emp(:,:)
  REAL(DP), ALLOCATABLE :: spreads_occ(:), spreads_emp(:)
  REAL(DP), ALLOCATABLE :: bstar(:,:)    ! first shell of b-vectors
  REAL(DP), ALLOCATABLE :: wb(:)         ! weight of the b-vectors 
  REAL(DP) :: b2                         ! squared modulus of the first shell
  INTEGER :: nbtot                       ! num of vectors in the first shell
  !
  ! ...  end of module-scope declarations
  !
  !---------------------------------------------------------------------
  !
CONTAINS
  !
  !
  SUBROUTINE get_centers_spreads( wfc, num_states, typ ) 
    !---------------------------------------------------------------------
    !
    ! ...  This routine calculates the centers and spreads of the variational
    ! ...  orbitals by using the Gamma-point algorithm of Wannier90
    ! ...  (see Comput. Phys. Commun. 178, 685 (2008))
    !
    ! ...  typ = 'occ' : occupied states
    ! ...  typ = 'emp' : empty states
    !
    !
    USE electrons_base,         ONLY : nspin
    USE cell_base,              ONLY : tpiba
    USE gvecw,                  ONLY : ngw
    USE cp_interfaces,          ONLY : invfft, fwfft
    USE fft_base,               ONLY : dffts
    USE smooth_grid_dimensions, ONLY : nnrsx
    USE control_flags,          ONLY : do_wf_cmplx
    USE mp,                     ONLY : mp_sum
    USE mp_global,              ONLY : intra_image_comm
    USE constants,              ONLY : tpi
    !
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: wfc(:,:)
    INTEGER, INTENT(IN) :: num_states
    CHARACTER(3), INTENT(IN) :: typ
    !
    REAL(DP), ALLOCATABLE :: centers(:,:)
    REAL(DP), ALLOCATABLE :: spreads(:)
    COMPLEX(DP) :: ZDOTC
    COMPLEX(DP), ALLOCATABLE :: Mnn(:)
    COMPLEX(DP), ALLOCATABLE :: psi(:)
    COMPLEX(DP), ALLOCATABLE :: wfcp(:), wfcm(:)
    COMPLEX(DP), ALLOCATABLE :: ph_aux(:)
    COMPLEX(DP), ALLOCATABLE :: phase(:)
    INTEGER :: norb
    INTEGER :: ib, bindex, ig, n
    !
    !
    norb = num_states / nspin
    ALLOCATE( centers(3,norb) )
    ALLOCATE( spreads(norb) )
    centers(:,:) = 0.D0
    spreads(:) = 0.D0
    !
    !
    IF ( do_wf_cmplx ) THEN
      WRITE( stdout, * ) '# WARNING: centers to check when do_wf_complex=.true.'
    ENDIF
    !
    !
    ALLOCATE( ph_aux(ngw) )
    ALLOCATE( psi(nnrsx), phase(nnrsx) )
    ALLOCATE( wfcp(ngw), wfcm(ngw) )
    ALLOCATE( Mnn(nbtot) )
    !
    !
    IF ( .not. ALLOCATED( bstar ) ) CALL star_of_k()
    !
    WRITE( stdout, 301 )
    !
    !
    DO n = 1, norb 
      !
      Mnn(:) = ( 0.D0, 0.D0 )
      !
      DO ib = 1, nbtot
        !
        ! FFT wfc to real-space -> psi
        !
        psi(:) = ( 0.D0, 0.D0 )
        CALL c2psi( psi, nnrsx, wfc(:,n), wfc(:,n), ngw, 1 ) ! what if wfc is complex ????
        CALL invfft( 'Wave', psi, dffts )
        !
        ! u_b(r) = e^(-ibr) * u(r)
        ! so here we calculate the phase factor e^(ibr)
        !
        bindex = -1
        ig = 1
        DO WHILE ( g(ig) .lt. b2 + 1.E-6 )
          !
          IF ( ( ABS( gx(1,ig) - bstar(1,ib) ) .lt. 1.E-6 ) .and. &
               ( ABS( gx(2,ig) - bstar(2,ib) ) .lt. 1.E-6 ) .and. &
               ( ABS( gx(3,ig) - bstar(3,ib) ) .lt. 1.E-6 ) ) bindex = ig
          !
          ig = ig + 1
          !
        ENDDO
        !
        ph_aux(:) = ( 0.D0, 0.D0 )
        ph_aux(bindex) = ( 1.D0, 0.D0 )
        CALL c2psi( phase, nnrsx, ph_aux, ph_aux, ngw, 0 )
        CALL invfft( 'Wave', phase, dffts )
        !
        psi(1:nnrsx) = psi(1:nnrsx) * phase(1:nnrsx)
        !
        ! FFT back to G-space
        !
        CALL fwfft( 'Wave', psi, dffts )
        CALL psi2c( psi, nnrsx, wfcp(:), wfcm(:), ngw, 2 ) ! what if wfc is complex ????
        IF ( g2_g(1) .lt. 1.E-6 ) wfcm(1) = (0.D0, 0.D0)
        wfcm(:) = CONJG(wfcm(:))
        !
        ! here we calculate Mb(n,n) = < u_n | u_nb >
        !
        Mnn(ib) = ZDOTC( ngw, wfcp, 1, wfc(:,n), 1 ) + ZDOTC( ngw, wfc(:,n), 1, wfcm, 1 )
        CALL mp_sum( Mnn(ib), intra_image_comm )
        !
        !
        centers(:,n) = centers(:,n) - wb(ib) * AIMAG( LOG( Mnn(ib) ) ) * bstar(:,ib)
        !!! RICCARDO debug >>>
        !WRITE(*,'("RICCARDO centers(",i2,",",i2,") = ",3f12.6)') ib,n,centers(:,n)
        !WRITE(*,'("RICCARDO M(",i2,",",i2,") = ",2f12.6,4x,"ImLogM = ",f12.6)') Mnn(ib),AIMAG( LOG( Mnn(ib) ) )
        !WRITE(*,*) aimag(log(Mnn(ib)))*wb(n)
        !!! RICCARDO debug <<<
        !
        spreads(n) = spreads(n) + wb(ib) * ( ( 1 - Mnn(ib) * CONJG( Mnn(ib) ) ) + &
                                         ( AIMAG( LOG( Mnn(ib) ) ) )**2 )
        !
      ENDDO ! ib
      !
      centers(:,n) = centers(:,n) / ( tpiba / BOHR_RADIUS_ANGS )
      spreads(n) = spreads(n) / ( tpiba / BOHR_RADIUS_ANGS ) ** 2
      !
      WRITE( stdout, 401 ) n, centers(:,n), spreads(n)
      !
    ENDDO ! n 
    !
    ! 
    301 FORMAT( //, 3x, 'Centers and Spreads of variational orbitals', /  &
                    
                    3x, 75('-'), / &
                    3x, 'Orbital #', 3x, '|', 3x, '  x [Ang]      y [Ang]      z [Ang]  ', 3x, '|', 3x, '< r^2 > [Ang^2]', / &
                    3x, 75('-') )
    !
    401 FORMAT( 3x, i5, 7x, 3f13.6, 7x, f14.8 ) 
    !
    !
    SELECT CASE ( typ )
      !
      CASE ( 'occ' )
        !
        IF ( .not. ALLOCATED( centers_occ ) ) THEN
          ALLOCATE( centers_occ(3,norb) )
          ALLOCATE( spreads_occ(norb) )
        ENDIF
        !
        centers_occ(:,:) = centers(:,:)
        spreads_occ(:) = spreads(:)
        !
      CASE ( 'emp' )
        !
        IF ( .not. ALLOCATED( centers_emp ) ) THEN
          ALLOCATE( centers_emp(3,norb) )
          ALLOCATE( spreads_emp(norb) )
        ENDIF
        !
        centers_emp(:,:) = centers(:,:)
        spreads_emp(:) = spreads(:)
        !
      CASE DEFAULT
        !
        CALL errore( 'get_centers_spreads', 'typ must be occ or emp', 1 )
        !
    END SELECT
    !
    !
  END SUBROUTINE get_centers_spreads
  !
  !
  SUBROUTINE star_of_k( ) 
    !---------------------------------------------------------------------
    !
    ! ...  This routine finds the stars (shells) of nearest-neighbor
    ! ...  b-vectors. The calculated vectors (in cartesian coordinates)
    ! ...  are saved in the global variable bstar.
    !
    ! ...  NB: for the moment only the first star, i.e. first nearest
    ! ...      neighbors, is calculated which is enough for cubic systems. 
    !
    !
    USE mp,                   ONLY : mp_sum, mp_gather
    USE mp_global,            ONLY : mpime, nproc, intra_image_comm
    USE cell_base,            ONLY : ibrav, at, tpiba
    !
    !
    IMPLICIT NONE
    !
    INTEGER, ALLOCATABLE :: gtot_g(:)
    INTEGER :: gcumul
    INTEGER :: i, ig
    REAL(DP) :: tpiba_ang
    !
    !
    ! first b-shell
    !
    b2 = MINVAL( g2_g, MASK = g2_g .gt. 1.E-6 )    ! skipping G=0
    !
    ! here we determine how many (gtot) G-vectors in the first shell each
    ! process has and we put all these values in a global array (gtot_g)
    !
    i = 1
    nbtot = 0
    ALLOCATE( gtot_g(nproc) )
    gtot_g(:) = 0
    !
    DO WHILE ( g(i) .lt. b2 + 1.E-6 )
      IF ( ABS( g(i) - b2 ) .lt. 1.E-6 ) nbtot = nbtot + 1
      i = i + 1
    ENDDO
    !
    gtot_g(mpime+1) = nbtot
    CALL mp_sum( gtot_g, intra_image_comm )
    CALL mp_sum( nbtot, intra_image_comm )
    ALLOCATE( bstar(3,nbtot) )
    bstar(:,:) = 0.D0
    !
    !
    ! check on nbtot - for the moment only for cubic systems (ibrav = 1,2,3) -
    ! keeping in mind that it has to be equal to half of the coordination
    ! number of the underlying reciprocal lattice
    !
    IF ( ibrav .eq. 1 .and. nbtot .ne. 3 ) CALL errore( 'star_of_k', &
    'nbtot not matching (1/2) the coordination number of the reciprocal lattice', nbtot*2 ) 
    IF ( ibrav .eq. 2 .and. nbtot .ne. 4 ) CALL errore( 'star_of_k', &
    'nbtot not matching (1/2) the coordination number of the reciprocal lattice', nbtot*2 ) 
    IF ( ibrav .eq. 3 .and. nbtot .ne. 6 ) CALL errore( 'star_of_k', &
    'nbtot not matching (1/2) the coordination number of the reciprocal lattice', nbtot*2 ) 
    !
    ! 
    ! now we put in bstar_ all the G-vectors in the first shell
    !
    gcumul = 0
    IF ( mpime .gt. 0 ) gcumul = SUM( gtot_g(:mpime) )
    !  
    i = 1
    ig = 1
    DO WHILE ( g(i) .lt. b2 + 1.E-6 )
      !
      IF ( ABS( g(i) - b2 ) .lt. 1.E-6 ) THEN
        bstar(:,gcumul+ig) = gx(:,i)
        ig = ig + 1
      ENDIF
      !
      i = i + 1
      !
    ENDDO
    !
    CALL mp_sum( bstar, intra_image_comm )
    !
    ALLOCATE( wb(nbtot) )
    wb(:) = 3.D0 / nbtot / b2
    !
    tpiba_ang = tpiba / BOHR_RADIUS_ANGS
    !
    !
    WRITE( stdout, 101 )
    WRITE( stdout, * )
    WRITE( stdout, 201 ) ( i, bstar(:,i)*tpiba_ang, wb(i)/tpiba_ang**2, i = 1, nbtot )
    !
    !
    101 FORMAT( //, 3x, 'First shell of b-vectors (b-vectors [Ang^-1] and weights [Ang^2])', / &
                    3x, '-----------------------------------------------------------------' )
    !
    201 FORMAT( 5x,       'b(', i2, ') =  [', 3f12.6, ' ] ,', 6x, 'wb = ', f8.6 ) 
    !
    !
  END SUBROUTINE star_of_k
  !
  !
END MODULE centers_and_spreads
