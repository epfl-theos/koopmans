!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE ldaU
  !--------------------------------------------------------------------------
  !! The quantities needed in lda+U calculations.
  !
  USE kinds,        ONLY : DP
  USE parameters,   ONLY : lqmax, ntypx
  USE basis,        ONLY : natomwfc
  USE ions_base,    ONLY : nat, ntyp => nsp, ityp
  !
  SAVE
  !
  INTEGER, PARAMETER :: nspinx=2
  COMPLEX(DP), ALLOCATABLE :: wfcU(:,:)
  !! atomic wfcs with U term
  COMPLEX(DP), ALLOCATABLE :: d_spin_ldau(:,:,:)
  !! the rotations in spin space for all symmetries
  REAL(DP) :: eth
  !! the Hubbard contribution to the energy
  REAL(DP) :: Hubbard_U(ntypx)
  !! the Hubbard U
  REAL(DP) :: Hubbard_J0(ntypx)
  !! the Hubbard J, in simplified LDA+U
  REAL(DP) :: Hubbard_J(3,ntypx)
  !! extra Hubbard parameters for full LDA+U:  
  !! * p: J(1)=J  
  !! * d: J(1)=J, J(2)=B   
  !! * f: J(1)=J, J(2)=E2, J(3)=E3   
  REAL(DP) :: Hubbard_alpha(ntypx)
  !! the Hubbard alpha (used to calculate U)
  REAL(DP) :: Hubbard_beta(ntypx)
  !! the Hubbard beta (used to calculate J0)
  REAL(DP) :: starting_ns(lqmax,nspinx,ntypx)
  !! starting ns
  INTEGER :: nwfcU
  !! total no. of atomic wavefunctions having U term
  INTEGER :: niter_with_fixed_ns
  !! no. of iterations with fixed ns
  INTEGER :: lda_plus_u_kind
  !! 1/0 --> full/simplified(old) LDA+U calculation
  INTEGER :: Hubbard_l(ntypx)
  !! the angular momentum of Hubbard states
  INTEGER :: Hubbard_lmax = 0 
  !! maximum angular momentum of Hubbard states
  LOGICAL :: is_hubbard(ntypx)
  !! .TRUE. if this atom species has U correction
  LOGICAL :: lda_plus_u
  !! .TRUE. if lda+u calculation is performed
  LOGICAL :: conv_ns
  !! .TRUE. if ns are converged
  CHARACTER(len=30) :: U_projection
  !! 'atomic', 'ortho-atomic', 'file'
  INTEGER, ALLOCATABLE :: oatwfc(:)
  !! specifies how input coordinates are given
  INTEGER, ALLOCATABLE :: offsetU(:) 
  !! offset of atomic wfcs used for projections
  REAL(DP), ALLOCATABLE :: q_ae(:,:,:)
  !! coefficients for projecting onto beta functions
  REAL(DP), ALLOCATABLE :: q_ps(:,:,:)
  !! (matrix elements on AE and PS atomic wfcs)
  !
CONTAINS
  !
  SUBROUTINE init_lda_plus_u( psd, noncolin )
    !
    !! NOTE: average_pp must be called before init_lda_plus_u
    !
    IMPLICIT NONE
    !
    CHARACTER (LEN=2), INTENT(IN) :: psd(:)
    LOGICAL, INTENT(IN) :: noncolin
    !
    INTEGER, EXTERNAL :: set_Hubbard_l
    INTEGER :: na, nt
    !
    !
    IF ( .NOT. lda_plus_u ) THEN
       Hubbard_lmax = 0
       is_hubbard(:) = .FALSE.
       RETURN
    ENDIF
    !
    Hubbard_lmax = -1
    ! Set the default of Hubbard_l for the species which have
    ! Hubbard_U=0 (in that case set_Hubbard_l will not be called)
    Hubbard_l(:) = -1
    !
    IF ( lda_plus_u_kind == 0 ) THEN
       !
       DO nt = 1, ntyp
          !
          is_hubbard(nt) = Hubbard_U(nt) /= 0.0_DP     .OR. &
                           Hubbard_alpha(nt) /= 0.0_DP .OR. &
                           Hubbard_J0(nt) /= 0.0_DP    .OR. &
                           Hubbard_beta(nt) /= 0.0_DP
          !
          IF ( is_hubbard(nt) ) THEN
             Hubbard_l(nt) = set_Hubbard_l( psd(nt) )
             Hubbard_lmax = MAX( Hubbard_lmax, Hubbard_l(nt) )
          ENDIF
          !
       ENDDO
       !
    ELSEIF ( lda_plus_u_kind == 1 ) THEN
       !
       IF ( U_projection == 'pseudo' ) CALL errore( 'init_lda_plus_u', &
            & 'full LDA+U not implemented with pseudo projection type', 1 )
       !
       IF (noncolin) THEN
          IF ( .NOT. ALLOCATED (d_spin_ldau) ) ALLOCATE( d_spin_ldau(2,2,48) )
          CALL comp_dspinldau()
       ENDIF
       !
       DO nt = 1, ntyp
          IF (Hubbard_alpha(nt)/=0.d0 ) CALL errore( 'init_lda_plus_u', &
               'full LDA+U does not support Hubbard_alpha calculation', 1 )

          is_hubbard(nt) = Hubbard_U(nt)/= 0.0_dp .OR. &
                         ANY( Hubbard_J(:,nt) /= 0.0_dp )
        
          IF ( is_hubbard(nt) ) THEN
             !
             Hubbard_l(nt) = set_Hubbard_l( psd(nt) )
             Hubbard_lmax = MAX( Hubbard_lmax, Hubbard_l(nt) )
             !
             IF (Hubbard_U(nt) == 0.0_dp) Hubbard_U(nt) = 1.d-14
             !
             IF ( Hubbard_l(nt) == 2 ) THEN
                IF ( Hubbard_J(2,nt) == 0.d0 ) &
                     Hubbard_J(2,nt) = 0.114774114774d0 * Hubbard_J(1,nt)
             ELSEIF ( Hubbard_l(nt) == 3 ) THEN
                IF ( Hubbard_J(2,nt) == 0.d0 ) &
                     Hubbard_J(2,nt) = 0.002268d0 * Hubbard_J(1,nt)
                IF ( Hubbard_J(3,nt)==0.d0 )   &
                     Hubbard_J(3,nt) = 0.0438d0 * Hubbard_J(1,nt)
             ENDIF
          ENDIF
          !
       ENDDO
       !
    ELSE
       !
       CALL errore( 'init_lda_plus_u', 'lda_plus_u_kind should be 0 or 1', 1 )
       !
    ENDIF
    !
    IF ( Hubbard_lmax == -1 ) CALL errore( 'init_lda_plus_u', &
         'lda_plus_u calculation but Hubbard_l not set', 1 )
    !
    IF ( Hubbard_lmax > 3 ) &
         CALL errore( 'init_lda_plus_u', 'Hubbard_l should not be > 3 ', 1 )
    !
    ! compute index of atomic wfcs used as projectors
    IF ( .NOT.ALLOCATED(oatwfc)) ALLOCATE( oatwfc(nat) )
    CALL offset_atom_wfc ( .FALSE., oatwfc, nwfcU )
    !
    ! nwfcU is set to natomwfc by the routine above
    IF ( nwfcU /= natomwfc ) &
         CALL errore( 'offset_atom_wfc', 'wrong number of wavefunctions', 1 )
    !
    ! for each atom, compute index of its projectors (among projectors only)
    IF ( .NOT.ALLOCATED(offsetU)) ALLOCATE( offsetU(nat) )
    CALL offset_atom_wfc( .TRUE., offsetU, nwfcU )
    !
  END SUBROUTINE init_lda_plus_u
  !
  !
  SUBROUTINE deallocate_ldaU( flag )
  !
  !
  LOGICAL, INTENT(IN) :: flag
  !
  IF ( flag ) THEN
     IF ( ALLOCATED( oatwfc ) )     DEALLOCATE( oatwfc )
     IF ( ALLOCATED( offsetU ) )    DEALLOCATE( offsetU )
     IF ( ALLOCATED( q_ae ) )       DEALLOCATE( q_ae )
     IF ( ALLOCATED( q_ps ) )       DEALLOCATE( q_ps )
     IF ( ALLOCATED( d_spin_ldau )) DEALLOCATE( d_spin_ldau )
  ENDIF
  IF ( ALLOCATED( wfcU ) )    DEALLOCATE( wfcU )
  !
  END SUBROUTINE deallocate_ldaU
  !
  !
  SUBROUTINE copy_U_wfc( swfcatom, noncolin )
  !
  !  Copy (orthogonalized) atomic wavefunctions "swfcatom"
  !  having a Hubbard U correction to array "wfcU"
  !
  IMPLICIT NONE
  COMPLEX(KIND=DP), INTENT(IN) :: swfcatom(:,:)
  LOGICAL, INTENT(IN), OPTIONAL :: noncolin
  LOGICAL :: twice
  INTEGER :: na, nt, m1, m2

  IF ( PRESENT(noncolin) ) THEN
     twice = noncolin
  ELSE
     twice = .FALSE.
  ENDIF
  DO na=1,nat
     nt = ityp(na)
     if ( is_hubbard(nt) ) THEN
        m1 = 1
        m2 = 2*hubbard_l(nt)+1
        IF ( twice ) m2 = 2*m2
        wfcU(:,offsetU(na)+m1:offsetU(na)+m2) = swfcatom(:,oatwfc(na)+m1:oatwfc(na)+m2)
     ENDIF
  ENDDO
  !
  END SUBROUTINE copy_U_wfc
  !
  !
END MODULE ldaU
!
