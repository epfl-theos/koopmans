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
  ! ... This module contains all the routines and variables important for
  ! ... the calculation of centers and spreads of the variational orbitals
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE reciprocal_vectors,   ONLY : g2_g, g, gx
  USE electrons_base,       ONLY : nspin
  USE constants,            ONLY : BOHR_RADIUS_ANGS
  !
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  REAL(DP), ALLOCATABLE, PUBLIC :: centers_occ(:,:,:), centers_emp(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: spreads_occ(:,:), spreads_emp(:,:)
  !
  REAL(DP), ALLOCATABLE :: bstar(:,:)    ! first shell of b-vectors
  REAL(DP), ALLOCATABLE :: wb(:)         ! weight of the b-vectors 
  REAL(DP) :: b2                         ! squared modulus of the first shell
  INTEGER :: nbtot                       ! num of vectors in the first shell
  INTEGER :: norb
  !
  PUBLIC :: get_centers_spreads, write_centers_and_spreads, &
            read_wannier_centers, read_wannier_spreads
  !
  ! ...  end of module-scope declarations
  !
  !---------------------------------------------------------------------
  !
CONTAINS
  !
  !
  SUBROUTINE get_centers_spreads( wfc, num_states, typ, units, verbose ) 
    !---------------------------------------------------------------------
    !
    ! ...  This routine calculates the centers and spreads of the variational
    ! ...  orbitals by using the Gamma-point algorithm of the Wannier90 code
    ! ...  (see Comput. Phys. Commun. 178, 685 (2008))
    !
    ! ...  typ = 'occ' : occupied states
    ! ...  typ = 'emp' : empty states
    !
    ! ...  verbose = .true. : at each call the routine prints out the current 
    ! ...                     centers and spreads
    !
    ! ...  units : 'bohr', 'ang', 'crystal', 'alat'
    ! ...          useful only when verbose is .true., defines in which units 
    ! ...          centers and spreads must be expressed (the stored global
    ! ...          variables are not affected and they are always expressed
    ! ...          in Bohr)
    !
    !
    USE io_global,              ONLY : ionode
    USE cell_base,              ONLY : tpiba, alat, bg, at
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
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: units
    LOGICAL, OPTIONAL, VALUE :: verbose
    !
    REAL(DP), ALLOCATABLE :: centers(:,:)
    REAL(DP), ALLOCATABLE :: spreads(:)
    REAL(DP) :: centers_(3), spreads_
    COMPLEX(DP) :: ZDOTC
    COMPLEX(DP), ALLOCATABLE :: Mnn(:)
    COMPLEX(DP), ALLOCATABLE :: psi(:)
    COMPLEX(DP), ALLOCATABLE :: wfcp(:), wfcm(:)
    COMPLEX(DP), ALLOCATABLE :: ph_aux(:)
    COMPLEX(DP), ALLOCATABLE :: phase(:)
    INTEGER :: iss
    INTEGER :: ib, bindex, ig, i, n, ns
    CHARACTER(LEN=3) :: typ_c
    CHARACTER(LEN=9) :: units_, units_c
    !
    CHARACTER(1), EXTERNAL :: capital
    !
    !
    IF ( do_wf_cmplx ) THEN
      WRITE( stdout, * ) '# WARNING: centers to check when do_wf_complex=.true.'
    ENDIF
    !
    IF ( .not. PRESENT(verbose) ) verbose = .false.
    !
    IF ( .not. PRESENT(units) ) THEN
      units_ = 'bohr'
    ELSE
      units_ = units
    ENDIF
    IF ( trim(units_) .ne. 'bohr' .and. trim(units_) .ne. 'ang' .and. &
         trim(units_) .ne. 'alat' .and. trim(units_) .ne. 'crystal' ) THEN
      CALL errore( 'get_centers_spreads', 'Invalid value assigned to units', 1 )
    ENDIF
    units_c = trim(units_)
    units_c(1:1) = capital( units_c(1:1) )
    !
    norb = num_states / nspin
    ALLOCATE( centers(3,norb) )
    ALLOCATE( spreads(norb) )
    CALL alloc_centers_spreads( typ )
    !
    IF ( .not. ALLOCATED( bstar ) ) CALL star_of_k()
    !
    ALLOCATE( ph_aux(ngw) )
    ALLOCATE( psi(nnrsx), phase(nnrsx) )
    ALLOCATE( wfcp(ngw), wfcm(ngw) )
    ALLOCATE( Mnn(nbtot) )
    !
    DO i = 1, LEN_TRIM(typ)
      typ_c(i:i) = capital( typ(i:i) )
    ENDDO
    !
    !
    ns = 0
    DO iss = 1, nspin
      !
      centers(:,:) = 0.D0
      spreads(:) = 0.D0
      !
      IF ( ionode .and. verbose ) THEN
        IF ( trim(units_) .eq. 'bohr' .or. trim(units_) .eq. 'ang' ) THEN
          WRITE( stdout, 301 ) typ_c, iss, (trim(units_c), i=1,4)
        ELSE
          WRITE( stdout, 302 ) typ_c, iss, (trim(units_c), i=1,3)
        ENDIF
      ENDIF
      !
      DO n = 1, norb 
        !
        ns = ns + 1
        Mnn(:) = ( 0.D0, 0.D0 )
        !
        DO ib = 1, nbtot
          !
          ! FFT wfc to real-space -> psi
          !
          psi(:) = ( 0.D0, 0.D0 )
          CALL c2psi( psi, nnrsx, wfc(:,ns), wfc(:,ns), ngw, 1 ) ! what if wfc is complex ????
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
          IF ( bindex .ne. -1 ) ph_aux(bindex) = ( 1.D0, 0.D0 )
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
          Mnn(ib) = ZDOTC( ngw, wfcp, 1, wfc(:,ns), 1 ) + ZDOTC( ngw, wfc(:,ns), 1, wfcm, 1 )
          CALL mp_sum( Mnn(ib), intra_image_comm )
          !
          !
          centers(:,n) = centers(:,n) - wb(ib) * AIMAG( LOG( Mnn(ib) ) ) * bstar(:,ib)
          !
          spreads(n) = spreads(n) + wb(ib) * ( ( 1 - Mnn(ib) * CONJG( Mnn(ib) ) ) + &
                                           ( AIMAG( LOG( Mnn(ib) ) ) )**2 )
          !
        ENDDO ! ib
        !
        ! centers and spreads are stored in Bohr units
        centers(:,n) = centers(:,n) / tpiba
        spreads(n) = spreads(n) / tpiba**2
        !
        SELECT CASE( trim(units_) )
          CASE ( 'ang' )
            centers_(:) = centers(:,n) * BOHR_RADIUS_ANGS
            spreads_ = spreads(n) * BOHR_RADIUS_ANGS**2
          CASE ( 'alat' )
            centers_(:) = centers(:,n) / alat
            spreads_ = spreads(n) / alat**2
          CASE ( 'crystal' )
            centers_(:) = centers(:,n) / alat
            CALL cryst_to_cart( 1, centers, bg, -1 )
            spreads_ = spreads(n) / alat**2
          CASE DEFAULT ! units=bohr
            centers_(:) = centers(:,n)
            spreads_ = spreads(n)
        END SELECT
        !
        IF ( ionode .and. verbose ) WRITE( stdout, 401 ) n, centers_, spreads_
        !
      ENDDO ! n
      !
      !
      IF ( typ == 'occ' ) THEN
        centers_occ(:,:,iss) = centers(:,:)
        spreads_occ(:,iss) = spreads(:)
      ELSE IF ( typ == 'emp' ) THEN
        centers_emp(:,:,iss) = centers(:,:)
        spreads_emp(:,iss) = spreads(:)
      ELSE
        CALL errore( 'get_centers_spreads', 'Invalid value assigned to typ', 1 )
      ENDIF
      !
      !
    ENDDO ! iss
    !
    ! 
301 FORMAT( //, 3x, 'Centers and Spreads of ', a3, ' variational orbitals, spin = ', i1, /  &
                3x, 76('-'), / &
                3x, 'Orbital #', 2x, '|', 3x, '  x [', a4, ']     y [', a4, ']     z [', a4, ']  ', 3x, '|', 3x, '< r^2 > [', a4, '^2]', / &
                3x, 76('-') )
302 FORMAT( //, 3x, 'Centers and Spreads of ', a3, ' variational orbitals, spin = ', i1, /  &
                3x, 76('-'), / &
                3x, 'Orbital #', 2x, '|', 3x, '  x [', a4, ']     y [', a4, ']     z [', a4, ']  ', 3x, '|', 3x, '< r^2 > [alat^2]', / &
                3x, 76('-') )
401 FORMAT( 3x, i5, 7x, 3f13.6, 8x, f14.8 ) 
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
!    WRITE( stdout, 101 )
!    WRITE( stdout, * )
!    WRITE( stdout, 201 ) ( i, bstar(:,i)*tpiba_ang, wb(i)/tpiba_ang**2, i = 1, nbtot )
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
  SUBROUTINE write_centers_and_spreads( num_wann, centers, spreads, emp )
    !---------------------------------------------------------------------
    !
    ! ...  This routine prints the centers and spreads of Wannier functions
    ! ...  (occupied or empty) to file. The units are Ang.
    !
    ! ...  emp = .true. when printing empty states
    !
    USE kinds,                ONLY : DP
    USE io_files,             ONLY : prefix
    USE io_global,            ONLY : ionode, stdout
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: num_wann        ! number of Wannier functions
    REAL(DP), INTENT(IN) :: centers(3,num_wann) 
    REAL(DP), INTENT(IN) :: spreads(num_wann) 
    LOGICAL, INTENT(IN) :: emp
    !
    CHARACTER(LEN=268) :: filename
    CHARACTER(LEN=9)  :: cdate, ctime
    CHARACTER(LEN=100) :: header
    INTEGER :: n
    !
    !
    IF ( emp ) THEN
      filename = TRIM(prefix)//'_cp_centers_emp.xyz'
    ELSE
      filename = TRIM(prefix)//'_cp_centers.xyz'
    ENDIF
    !
    CALL date_and_tim( cdate, ctime )
    header = 'Wannier centers and spreads, written by CP on '//cdate//' at '//ctime
    !
    IF ( ionode ) THEN
      !
      OPEN( 300, file=filename, status='unknown' )
      !
      WRITE( 300, '(i6)' ) num_wann
      WRITE( 300, * ) header
      !
      WRITE( 300, 55 ) ( centers(:,n), spreads(n), n = 1, num_wann )
      !
      CLOSE( 300 )
      !
    ENDIF
    !
    !
    55 FORMAT( 'X', 5x, 3f17.8, 2x, f17.8 )
    !
    !
  END SUBROUTINE write_centers_and_spreads
  !
  !
  SUBROUTINE read_wannier_centers( centers, num_wann, emp )
    !---------------------------------------------------------------------
    !
    ! ...  This routine reads the centers of Wannier functions from .xyz
    ! ...  file print out by Wannier90, fold them into the R=0 primitive
    ! ...  cell and gives them in output (in crystal units)
    !
    ! ...  emp = .true. when reading empty states
    !
    USE kinds,                ONLY : DP
    USE cell_base,            ONLY : bg, alat
    USE constants,            ONLY : BOHR_RADIUS_ANGS
    USE io_files,             ONLY : prefix
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: num_wann        ! number of Wannier functions 
    LOGICAL, INTENT(IN) :: emp             
    !
    REAL(DP), INTENT(OUT) :: centers(3,num_wann)
    !
    LOGICAL :: exst
    INTEGER :: n 
    CHARACTER(LEN=268) :: filename
    CHARACTER(LEN=256) :: input_line
    !
    !
    IF ( emp ) THEN
      filename = trim(prefix)//'_emp_centres.xyz'
    ELSE
      filename = trim(prefix)//'_centres.xyz'
    ENDIF
    !
    INQUIRE( file=filename, exist=exst )
    !
    IF ( .not. exst ) CALL errore( 'read_wannier_centers', 'File not found', 1 )
    !
    OPEN( 100, file=filename, form='formatted', status='old' )
    !
    READ( 100, *, end=10, err=20 )    ! skip 1st line
    READ( 100, *, end=10, err=20 )    ! skip 2nd line
         !
    DO n = 1, num_wann
      !
      READ( 100, '(a256)', end=10, err=20 ) input_line
      !
      IF ( input_line(1:1) .ne. 'X' ) CALL errore( 'read_wannier_centers', &
              'X must precede each Wannier center line', 1 )
      !
      READ( input_line(2:), *, end=10, err=20 ) centers(:,n)
      !
    ENDDO
    !
    READ( 100, * ) input_line
    IF ( input_line(1:1) == 'X' ) CALL errore( 'read_wannier_centers', &
            'Missing some center!', 1 )
    !
    CLOSE( 100 )
    !
    !
    centers = centers / ( alat * BOHR_RADIUS_ANGS )
    !
    DO n = 1, num_wann
      !
      CALL cryst_to_cart( 1, centers(:,n), bg, -1 )
      !
    ENDDO
    !
    RETURN
    !
  10  CALL errore ( 'read_wannier_centers', 'end of file while reading', 1 )
  20  CALL errore ( 'read_wannier_centers', 'error while reading', 1 )
    !
    !
  END SUBROUTINE read_wannier_centers
  !
  !
  SUBROUTINE read_wannier_spreads( spreads, num_wann, emp )
    !---------------------------------------------------------------------
    !
    ! ...  This routine reads the spreads of Wannier functions from .wout
    ! ...  file print out by Wannier90, gives them in output (in Ang^2)
    !
    ! ...  emp = .true. when reading empty states
    !
    USE io_files,             ONLY : prefix
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: num_wann        ! number of Wannier functions 
    LOGICAL, INTENT(IN) :: emp
    !
    REAL(DP), INTENT(OUT) :: spreads(num_wann)
    !
    LOGICAL :: exst
    INTEGER :: n 
    CHARACTER(LEN=268) :: filename
    CHARACTER(LEN=256) :: input_line
    !
    !
    IF ( emp ) THEN
      filename = trim(prefix)//'_emp.wout'
    ELSE
      filename = trim(prefix)//'.wout'
    ENDIF
    !
    INQUIRE( file=filename, exist=exst )
    !
    IF ( .not. exst ) CALL errore( 'read_wannier_spreads', 'File not found', 1 )
    !
    OPEN( 200, file=filename, form='formatted', status='old' )
    !
    READ( 200, '(a256)', end=10, err=20 ) input_line
    DO WHILE ( input_line .ne. ' Final State' )
      READ( 200, '(a256)', end=10, err=20 ) input_line
    ENDDO
    !
    DO n = 1, num_wann
      !
      READ( 200, '(a256)', end=10, err=20 ) input_line
      READ( input_line(65:), * ) spreads(n)
      !
    ENDDO
    !
    CLOSE( 200 )
    !
    !
    RETURN
    !
  10  CALL errore ( 'read_wannier_spreads', 'end of file while reading', 1 )
  20  CALL errore ( 'read_wannier_spreads', 'error while reading', 1 )
    !
    !
  END SUBROUTINE read_wannier_spreads
  !
  !
  SUBROUTINE alloc_centers_spreads( typ )
    !---------------------------------------------------------------------
    !
    ! ...  allocates centers_occ, centers_emp, spreads_occ and spreads_emp
    !
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=3), INTENT(IN) :: typ
    !
    !
    SELECT CASE ( typ )
      !
      CASE ( 'occ' )
        !
        IF ( .not. ALLOCATED( centers_occ ) ) THEN
          ALLOCATE( centers_occ(3,norb,nspin) )
          ALLOCATE( spreads_occ(norb,nspin) )
        ENDIF
        !
      CASE ( 'emp' )
        !
        IF ( .not. ALLOCATED( centers_emp ) ) THEN
          ALLOCATE( centers_emp(3,norb,nspin) )
          ALLOCATE( spreads_emp(norb,nspin) )
        ENDIF
        !
      CASE DEFAULT
        !
        CALL errore( 'alloc_centers_spreads', 'typ must be occ or emp', 1 )
        !
    END SELECT
    !
    !
  END SUBROUTINE alloc_centers_spreads
  !
  !
END MODULE centers_and_spreads
