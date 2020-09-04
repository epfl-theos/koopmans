
! Copyright (C) 2004 PWSCF group 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 

!------------------------------
 MODULE grid_module
!------------------------------
  USE kinds,        ONLY : DP
  IMPLICIT NONE
  PRIVATE

  !
  ! general purpose vars 
  !
  REAL(DP), ALLOCATABLE  :: focc(:,:), wgrid(:)
  REAL(DP)               :: alpha 
  !
  !
  PUBLIC :: grid_build, grid_destroy
  PUBLIC :: focc, wgrid, alpha

CONTAINS

!---------------------------------------------
  SUBROUTINE grid_build(nw, wmax, wmin)
  !-------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE wvfct,     ONLY : nbnd, wg
  USE klist,     ONLY : nks, wk, nelec
  USE lsda_mod,  ONLY : nspin
  USE uspp,      ONLY : okvan
  !
  IMPLICIT NONE
  !
  ! input vars
  INTEGER,  INTENT(IN) :: nw
  REAL(DP), INTENT(IN) :: wmax ,wmin
  !
  ! local vars
  INTEGER         :: iw,ik,i,ierr
 
  !
  ! check on the number of bands: we need to include empty bands in order to allow
  ! to write the transitions
  !
  ! IF ( REAL(nbnd, DP)  <= nelec / 2.0_DP ) CALL errore('epsilon', 'ban band number', 1)

  !
  ! spin is not implemented
  !
  IF( nspin > 2 ) CALL errore('grid_build','Non collinear spin  calculation not implemented',1)

  !
  ! USPP are not implemented (dipole matrix elements are not trivial at all)
  !
  !IF ( okvan ) CALL errore('grid_build','USPP are not implemented',1)

  ALLOCATE ( focc( nbnd, nks), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating focc', ABS(ierr))
  !
  ALLOCATE( wgrid( nw ), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating wgrid', ABS(ierr))

  !
  ! check on k point weights, no symmetry operations are allowed 
  !
  DO ik = 2, nks
     !
     IF ( ABS( wk(1) - wk(ik) ) > 1.0d-8 ) &
        CALL errore('grid_build','non unifrom kpt grid', ik )
     !
  ENDDO
  !
  ! occupation numbers, to be normalized differently 
  ! whether we are spin resolved or not
  !
  IF(nspin==1) THEN
    DO ik = 1,nks
    DO i  = 1,nbnd
         focc(i,ik)= wg(i, ik ) * 2.0_DP / wk( ik )
    ENDDO
    ENDDO
  ELSE IF(nspin==2) THEN
    DO ik = 1,nks
    DO i  = 1,nbnd
         focc(i,ik)= wg(i, ik ) * 1.0_DP / wk( ik )
    ENDDO
    ENDDO
  ENDIF
  !
  ! set the energy grid
  !
  alpha = (wmax - wmin) / REAL(nw, DP)
  !
  DO iw = 1, nw 
      wgrid(iw) = wmin + iw * alpha
  ENDDO
  !
END SUBROUTINE grid_build
!
!
!----------------------------------
  SUBROUTINE grid_destroy
  !----------------------------------
  IMPLICIT NONE
  INTEGER :: ierr
  !
  IF ( ALLOCATED( focc) ) THEN
      !
      DEALLOCATE ( focc, wgrid, STAT=ierr)
      CALL errore('grid_destroy','deallocating grid stuff',ABS(ierr))
      !
  ENDIF
  !
END SUBROUTINE grid_destroy

END MODULE grid_module


!------------------------------
PROGRAM epsilon
!------------------------------
  !
  ! Compute the complex macroscopic dielectric function,
  ! at the RPA level, neglecting local field effects.
  ! Eps is computed both on the real or immaginary axis
  !
  ! Authors: Andrea Benassi, Andrea Ferretti, Carlo Cavazzoni
  !
  ! NOTE: Part of the basic implementation is taken from pw2gw.f90;
  !
  !
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : stdout, ionode, ionode_id
  USE mp,          ONLY : mp_bcast
  USE iotk_module
  USE xml_io_base
  USE io_files,    ONLY : nd_nmbr, tmp_dir, prefix, outdir, trimcheck
  USE constants,   ONLY : RYTOEV
  USE ener,        ONLY : ef
  USE klist,       ONLY : lgauss  
  USE ktetra,      ONLY : ltetra
  USE wvfct,       ONLY : nbnd
  USE lsda_mod,    ONLY : nspin
  ! 
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER                 :: nw,nbndmin,nbndmax
  REAL(DP)                :: intersmear,intrasmear,wmax,wmin,shift, & 
                             ekin_eout, ekin_error , photon_ener, & 
                             polar_angle, azimuthal_angle, & 
                             photon_angle, e_fermi
  CHARACTER(10)           :: calculation,smeartype
  LOGICAL                 :: metalcalc, homo_gas, wfc_real, & 
                             modified_pw, othor_pw
  !
  NAMELIST / inputpp / prefix, outdir, calculation
  NAMELIST / energy_grid / smeartype,intersmear,intrasmear,wmax,wmin, & 
                           nbndmin,nbndmax,nw,shift,ekin_eout,ekin_error, &
                           polar_angle, azimuthal_angle, photon_angle, &
                           photon_ener, homo_gas, wfc_real, e_fermi,   &
                           modified_pw, othor_pw                          
  
  !
  ! local variables
  !
  INTEGER :: ios

!--------------------------------------------- 
! program body
!--------------------------------------------- 
!

  CALL init_clocks( .TRUE. )
  CALL start_clock( 'epsilon' )
  !
  !
  CALL start_postproc(nd_nmbr)

  !
  ! Set default values for variables in namelist 
  !
  calculation  = 'eps'
  prefix       = 'pwscf'
  shift        = 0.0d0
  outdir       = './'
  intersmear   = 0.136
  wmin         = 0.0d0
  wmax         = 30.0d0
  nbndmin      = 1
  nbndmax      = 0 
  nw           = 600
  smeartype    = 'gauss'
  intrasmear   = 0.0d0 
  metalcalc    = .FALSE. 
  !
  ekin_eout    = 30.0d0 
  ekin_error   = 1.0d0 
  polar_angle  = 30
  azimuthal_angle  = 30
  photon_angle  = 50
  photon_ener  = 80
  e_fermi      = 4.0
  wfc_real     = .true.
  homo_gas     = .true.
  modified_pw  = .true.
  othor_pw     = .false.
  !
  !
  ! this routine allows the user to redirect the input using -input 
  ! instead of <  
  !
  CALL input_from_file( )

  !
  ! read input file
  !
  IF (ionode) WRITE( stdout, "( 2/, 5x, 'Reading input file...' ) " )
  ios = 0 
  !
  IF ( ionode ) READ (5, inputpp, IOSTAT=ios)
  !
  CALL mp_bcast ( ios, ionode_id ) 
  IF (ios/=0) CALL errore('epsilon', 'reading namelist INPUTPP', ABS(ios))
  !
  IF ( ionode ) THEN 
     !
     READ (5, energy_grid, IOSTAT=ios)
     !
     tmp_dir = trimcheck(outdir)
     !
  ENDIF
  !
  CALL mp_bcast ( ios, ionode_id ) 
  IF (ios/=0) CALL errore('epsilon', 'reading namelist ENERGY_GRID', ABS(ios))
  ! 
  ! ... Broadcast variables 
  !
  IF (ionode) WRITE( stdout, "( 5x, 'Broadcasting variables...' ) " )

  CALL mp_bcast( smeartype, ionode_id ) 
  CALL mp_bcast( calculation, ionode_id )
  CALL mp_bcast( prefix, ionode_id ) 
  CALL mp_bcast( tmp_dir, ionode_id ) 
  CALL mp_bcast( shift, ionode_id ) 
  CALL mp_bcast( outdir, ionode_id ) 
  CALL mp_bcast( intrasmear, ionode_id )
  CALL mp_bcast( intersmear, ionode_id) 
  CALL mp_bcast( wmax, ionode_id ) 
  CALL mp_bcast( wmin, ionode_id ) 
  CALL mp_bcast( nw, ionode_id ) 
  CALL mp_bcast( nbndmin, ionode_id ) 
  CALL mp_bcast( nbndmax, ionode_id ) 
  !
  CALL mp_bcast( ekin_eout, ionode_id ) 
  CALL mp_bcast( ekin_error, ionode_id ) 
  CALL mp_bcast( polar_angle, ionode_id ) 
  CALL mp_bcast( azimuthal_angle, ionode_id ) 
  CALL mp_bcast( photon_angle, ionode_id ) 
  CALL mp_bcast( photon_ener, ionode_id ) 
  CALL mp_bcast( homo_gas, ionode_id ) 
  CALL mp_bcast( wfc_real, ionode_id ) 
  CALL mp_bcast( e_fermi, ionode_id ) 
  CALL mp_bcast( modified_pw, ionode_id ) 
  CALL mp_bcast( othor_pw, ionode_id ) 
  !
  ! read PW simulation parameters from prefix.save/data-file.xml 
  !
  IF (ionode) WRITE( stdout, "( 5x, 'Reading PW restart file...' ) " )
  
  CALL read_file
  CALL openfil_pp

  !
  ! few conversions
  !
  
  IF (ionode) WRITE(stdout,"(2/, 5x, 'Fermi energy [eV] is: ',f8.5)") ef *RYTOEV

  IF (lgauss .OR. ltetra) THEN 
      metalcalc=.true. 
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a metal...' ) " )
  ELSE
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a dielectric...' ) " )
  ENDIF

  IF (nbndmax == 0) nbndmax = nbnd

  !
  ! ... run the specific pp calculation
  !
  IF (ionode) WRITE(stdout,"(/, 5x, 'Performing ',a,' calculation...')") TRIM(calculation)

  CALL start_clock( 'calculation' )
  !
  SELECT CASE ( TRIM(calculation) )
  !
  CASE ( 'eps' )
      !
      CALL eps_calc ( intersmear,intrasmear,nw,wmax,wmin,nbndmin,nbndmax,shift,metalcalc,nspin )
      !
  CASE ( 'jdos' )
      !
      CALL jdos_calc ( smeartype,intersmear,nw,wmax,wmin,nbndmin,nbndmax,shift,nspin )
      !
  CASE ( 'photospec' )
      !
      CALL photoemission_spectr_pw ( smeartype,intersmear,intrasmear,nw,wmax,wmin,nbndmin,nbndmax, &
                                     shift,nspin,metalcalc,polar_angle, azimuthal_angle,&
                                     photon_angle,photon_ener,homo_gas,wfc_real,e_fermi, &
                                     modified_pw, othor_pw)
      !
  CASE ( 'offdiag' )
      !
      CALL offdiag_calc ( intersmear,intrasmear,nw,wmax,wmin,nbndmin,nbndmax,shift,metalcalc,nspin )
      !
  CASE ( 'occ' )
      !
      CALL occ_calc ()
      !
  CASE DEFAULT
      !
      CALL errore('epsilon','invalid CALCULATION = '//TRIM(calculation),1)
      !
  END SELECT
  !
  CALL stop_clock( 'calculation' )
  
  !
  ! few info about timing
  !
  CALL stop_clock( 'epsilon' )
  !
  IF ( ionode ) WRITE( stdout , "(/)" )
  !
  CALL print_clock( 'epsilon' )
  CALL print_clock( 'calculation' )
  CALL print_clock( 'dipole_calc' )
  !
  IF ( ionode ) WRITE( stdout, *  )

  !
  !
  CALL stop_pp ()

END PROGRAM epsilon 


!-----------------------------------------------------------------------------
SUBROUTINE eps_calc ( intersmear,intrasmear, nw, wmax, wmin, nbndmin, nbndmax, shift, &
                      metalcalc , nspin)
  !-----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE cell_base,            ONLY : tpiba2, omega
  USE wvfct,                ONLY : nbnd, et
  USE ener,                 ONLY : efermi => ef
  USE klist,                ONLY : nks, nkstot, degauss
  USE io_global,            ONLY : ionode, stdout
  !
  USE grid_module,          ONLY : alpha, focc, wgrid, grid_build, grid_destroy                
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum
  ! 
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(IN) :: nw,nbndmin,nbndmax,nspin
  REAL(DP),        INTENT(IN) :: wmax, wmin, intersmear,intrasmear, shift
  LOGICAL,         INTENT(IN) :: metalcalc 
  !
  ! local variables
  !
  INTEGER       :: i, ik, iband1, iband2,is 
  INTEGER       :: iw, iwp, ierr
  REAL(DP)      :: etrans, const, w, renorm(3)
  !
  REAL(DP), ALLOCATABLE    :: epsr(:,:), epsi(:,:), epsrc(:,:,:), epsic(:,:,:)
  REAL(DP), ALLOCATABLE    :: ieps(:,:), eels(:,:), iepsc(:,:,:), eelsc(:,:,:)
  REAL(DP), ALLOCATABLE    :: dipole(:,:,:)
  COMPLEX(DP),ALLOCATABLE  :: dipole_aux(:,:,:)
!
!--------------------------
! main routine body
!--------------------------
!
  !
  ! perform some consistency checks, calculate occupation numbers and setup w grid
  !
  CALL grid_build(nw, wmax, wmin)
  !   
  ! allocate main spectral and auxiliary quantities   
  !   
  ALLOCATE( dipole(3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole', ABS(ierr) )
  !
  ALLOCATE( dipole_aux(3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole_aux', ABS(ierr) )

!
! spin unresolved calculation
!
IF (nspin == 1) THEN 
  !
  ALLOCATE( epsr( 3, nw), epsi( 3, nw), eels( 3, nw), ieps(3,nw ), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating eps', ABS(ierr))

  !
  ! initialize response functions
  !
  epsr(:,:)  = 0.0_DP
  epsi(:,:)  = 0.0_DP
  ieps(:,:)  = 0.0_DP

  !
  ! main kpt loop
  !
  kpt_loop: &
  DO ik = 1, nks
     !
     ! For every single k-point: order k+G for
     !                           read and distribute wavefunctions
     !                           compute dipole matrix 3 x nbnd x nbnd parallel over g 
     !                           recover g parallelism getting the total dipole matrix
     ! 
     CALL dipole_calc( ik, dipole_aux, metalcalc , nbndmin, nbndmax)
     !
     dipole(:,:,:)= tpiba2 * REAL( dipole_aux(:,:,:) * CONJG(dipole_aux(:,:,:)), DP ) 
      
     !
     ! Calculation of real and immaginary parts 
     ! of the macroscopic dielettric function from dipole
     ! approximation. 
     ! 'intersmear' is the brodening parameter  
     !
     !Interband
     ! 
     DO iband2 = nbndmin,nbndmax
         !
         IF ( focc(iband2,ik) < 2.0d0) THEN
     DO iband1 = nbndmin,nbndmax 
         !
         IF (iband1==iband2) CYCLE
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
         IF (ABS(focc(iband2,ik)-focc(iband1,ik))< 1e-3) CYCLE
               !
               ! transition energy
               !
               etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV + shift
               !
               ! loop over frequencies
               !
               DO iw = 1, nw 
                   !
                   w = wgrid(iw)
                   !
                   epsi(:,iw) = epsi(:,iw) + dipole(:,iband1,iband2) * intersmear * w* &
                                             RYTOEV**3 * (focc(iband1,ik))/  &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans ) 
                                   
                   epsr(:,iw) = epsr(:,iw) + dipole(:,iband1,iband2) * RYTOEV**3 * &
                                             (focc(iband1,ik)) * &
                                             (etrans**2 - w**2 ) / &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans ) 
               ENDDO

         ENDIF
     ENDDO
         ENDIF
     ENDDO
     !
     !Intraband (only if metalcalc is true)
     !
     IF (metalcalc) THEN
     DO iband1 = nbndmin,nbndmax
         !
         IF ( focc(iband1,ik) < 2.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
               !
               ! loop over frequencies
               !
               DO iw = 1, nw 
                   !
                   w = wgrid(iw)
                   !
                  epsi(:,iw) = epsi(:,iw) +  dipole(:,iband1,iband1) * intrasmear * w* &
                                RYTOEV**2 * (EXP((et(iband1,ik)-efermi)/degauss ))/  &
                    (( w**4 + intrasmear**2 * w**2 )*(1+EXP((et(iband1,ik)-efermi)/ & 
                    degauss))**2*degauss ) 
                                   
                  epsr(:,iw) = epsr(:,iw) - dipole(:,iband1,iband1) * RYTOEV**2 * &
                                            (EXP((et(iband1,ik)-efermi)/degauss )) * w**2 / & 
                    (( w**4 + intrasmear**2 * w**2 )*(1+EXP((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss )
               ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDIF
  ENDDO kpt_loop

  !
  ! recover over kpt parallelization (inter_pool)
  !
  CALL mp_sum( epsr, inter_pool_comm )
  CALL mp_sum( epsi, inter_pool_comm )
  
  !
  ! impose the correct normalization
  !
  const = 64.0d0 * PI / ( omega * REAL(nkstot, DP) )
  epsr(:,:) = 1.0_DP + epsr(:,:) * const  
  epsi(:,:) =          epsi(:,:) * const 

  !
  ! Calculation of eels spectrum 
  !
  DO iw = 1, nw
      !
      eels(:,iw) = epsi(:,iw) / ( epsr(:,iw)**2 + epsi(:,iw)**2 )
      !
  ENDDO
  
  !
  !  calculation of dielectric function on the immaginary frequency axe
  !                   

  DO iw = 1, nw
  DO iwp = 2, nw 
      !
      ieps(:,iw) = ieps(:,iw) + wgrid(iwp) * epsi(:,iwp) / ( wgrid(iwp)**2 + wgrid(iw)**2) 
      !
  ENDDO
  ENDDO               

  ieps(:,:) = 1.0d0 + 2 / PI * ieps(:,:) * alpha 

  !
  ! check  dielectric function  normalizzation via sumrule  
  !
 DO i=1,3  
     renorm(i) = alpha * SUM( epsi(i,:) * wgrid(:) ) 
 ENDDO
  !
  IF ( ionode ) THEN
      !
      WRITE(stdout,"(/,5x, 'The bulk xx plasmon frequency [eV] is: ',f15.9 )")  SQRT(renorm(1) * 2.0d0 / PI)
      WRITE(stdout,"(5x, 'The bulk yy plasmon frequency [eV] is: ',f15.9 )")  SQRT(renorm(2) * 2.0d0 / PI)
      WRITE(stdout,"(5x, 'The bulk zz plasmon frequency [eV] is: ',f15.9 )")  SQRT(renorm(3) * 2.0d0 / PI)
      WRITE(stdout,"(/,5x, 'Writing output on file...' )")
      !
      ! write results on data files
      !

      OPEN (30, FILE='epsr.dat', FORM='FORMATTED' )
      OPEN (40, FILE='epsi.dat', FORM='FORMATTED' )
      OPEN (41, FILE='eels.dat', FORM='FORMATTED' )
      OPEN (42, FILE='ieps.dat', FORM='FORMATTED' )
      !
      WRITE(30, "(2x,'# energy grid [eV]     epsr_x  epsr_y  epsr_z')" )
      WRITE(40, "(2x,'# energy grid [eV]     epsi_x  epsi_y  epsi_z')" )
      WRITE(41, "(2x,'# energy grid [eV]  eels components [arbitrary units]')" )
      WRITE(42, "(2x,'# energy grid [eV]     ieps_x  ieps_y  ieps_z ')" )
      !
      DO iw =1, nw
          !
          WRITE(30,"(4f15.6)") wgrid(iw), epsr(1:3, iw)
          WRITE(40,"(4f15.6)") wgrid(iw), epsi(1:3, iw)
          WRITE(41,"(4f15.6)") wgrid(iw), eels(1:3, iw)
          WRITE(42,"(4f15.6)") wgrid(iw), ieps(1:3, iw)
          !
      ENDDO
      !
      CLOSE(30)
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      !
  ENDIF

  DEALLOCATE ( epsr, epsi, eels, ieps)
!
! collinear spin calculation
!
ELSE IF (nspin == 2 ) THEN
  !
  ALLOCATE( epsrc( 0:1, 3, nw), epsic( 0:1,3, nw), eelsc( 0:1,3, nw), iepsc(0:1,3,nw ), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating eps', ABS(ierr))

  !
  ! initialize response functions
  !
  epsrc(:,:,:)  = 0.0_DP
  epsic(:,:,:)  = 0.0_DP
  iepsc(:,:,:)  = 0.0_DP

  !
  ! main kpt loop
  !

spin_loop: &
DO is=0,1
  kpt_loopspin: &
! if nspin=2 the number of nks must be even (even if the calculation 
! is performed at gamma point only), so nks must be always a multiple of 2
  DO ik = 1 + is * INT(nks/2), INT(nks/2) +  is * INT(nks/2)
     !
     ! For every single k-point: order k+G for
     !                           read and distribute wavefunctions
     !                           compute dipole matrix 3 x nbnd x nbnd parallel over g 
     !                           recover g parallelism getting the total dipole matrix
     ! 
     CALL dipole_calc( ik, dipole_aux, metalcalc , nbndmin, nbndmax)
     !
     dipole(:,:,:)= tpiba2 * REAL( dipole_aux(:,:,:) * CONJG(dipole_aux(:,:,:)), DP ) 
      
     !
     ! Calculation of real and immaginary parts 
     ! of the macroscopic dielettric function from dipole
     ! approximation. 
     ! 'intersmear' is the brodening parameter  
     !
     !Interband
     ! 
     DO iband2 = nbndmin,nbndmax
         !
         IF ( focc(iband2,ik) < 1.0d0) THEN
     DO iband1 = nbndmin,nbndmax 
         !
         IF (iband1==iband2) CYCLE
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
         IF (ABS(focc(iband2,ik)-focc(iband1,ik))< 1e-3) CYCLE
               !
               ! transition energy
               !
               etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV + shift
               !
               ! loop over frequencies
               !
               DO iw = 1, nw 
                   !
                   w = wgrid(iw)
                   !
                   epsic(is,:,iw) = epsic(is,:,iw) + dipole(:,iband1,iband2) * intersmear * w* &
                                             RYTOEV**3 * (focc(iband1,ik))/  &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans ) 
                                   
                   epsrc(is,:,iw) = epsrc(is,:,iw) + dipole(:,iband1,iband2) * RYTOEV**3 * &
                                             (focc(iband1,ik)) * &
                                             (etrans**2 - w**2 ) / &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans ) 
               ENDDO

         ENDIF
     ENDDO
         ENDIF
     ENDDO
     !
     !Intraband (only if metalcalc is true)
     !
     IF (metalcalc) THEN
     DO iband1 = nbndmin,nbndmax
         !
         IF ( focc(iband1,ik) < 1.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
               !
               ! loop over frequencies
               !
               DO iw = 1, nw 
                   !
                   w = wgrid(iw)
                   !
                  epsic(is,:,iw) = epsic(is,:,iw) +  dipole(:,iband1,iband1) * intrasmear * w* &
                                RYTOEV**2 * (EXP((et(iband1,ik)-efermi)/degauss ))/  &
                    (( w**4 + intrasmear**2 * w**2 )*(1+EXP((et(iband1,ik)-efermi)/ & 
                    degauss))**2*degauss ) 
                                   
                  epsrc(is,:,iw) = epsrc(is,:,iw) - dipole(:,iband1,iband1) * RYTOEV**2 * &
                                            (EXP((et(iband1,ik)-efermi)/degauss )) * w**2 / & 
                    (( w**4 + intrasmear**2 * w**2 )*(1+EXP((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss )
               ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDIF
  ENDDO kpt_loopspin
ENDDO spin_loop
  !
  ! recover over kpt parallelization (inter_pool)
  !
  CALL mp_sum( epsr, inter_pool_comm )
  CALL mp_sum( epsi, inter_pool_comm )
  
  !
  ! impose the correct normalization
  !
  const = 128.0d0 * PI / ( omega * REAL(nkstot, DP) )
  epsrc(:,:,:) = 1.0_DP + epsrc(:,:,:) * const  
  epsic(:,:,:) =          epsic(:,:,:) * const 

  !
  ! Calculation of eels spectrum 
  !
  DO iw = 1, nw
      !
      eelsc(:,:,iw) = epsic(:,:,iw) / ( epsrc(:,:,iw)**2 + epsic(:,:,iw)**2 )
      !
  ENDDO
  
  !
  !  calculation of dielectric function on the immaginary frequency axe
  !                   

  DO iw = 1, nw
  DO iwp = 2, nw 
      !
      iepsc(:,:,iw) = iepsc(:,:,iw) + wgrid(iwp) * epsic(:,:,iwp) / ( wgrid(iwp)**2 + wgrid(iw)**2) 
      !
  ENDDO
  ENDDO               

  iepsc(:,:,:) = 1.0d0 + 2.0_DP / PI * iepsc(:,:,:) * alpha 

  IF (ionode) THEN
      WRITE(stdout,"(/,5x, 'Writing output on file...' )")
      !
      ! write results on data files
      !

      OPEN (30, FILE='uepsr.dat', FORM='FORMATTED' )
      OPEN (40, FILE='uepsi.dat', FORM='FORMATTED' )
      OPEN (41, FILE='ueels.dat', FORM='FORMATTED' )
      OPEN (42, FILE='uieps.dat', FORM='FORMATTED' )
      OPEN (43, FILE='depsr.dat', FORM='FORMATTED' )
      OPEN (44, FILE='depsi.dat', FORM='FORMATTED' )
      OPEN (45, FILE='deels.dat', FORM='FORMATTED' )
      OPEN (46, FILE='dieps.dat', FORM='FORMATTED' )
      OPEN (47, FILE='epsr.dat', FORM='FORMATTED' )
      OPEN (48, FILE='epsi.dat', FORM='FORMATTED' )
      OPEN (49, FILE='eels.dat', FORM='FORMATTED' )
      OPEN (50, FILE='ieps.dat', FORM='FORMATTED' )
      !
      WRITE(30, "(2x,'# energy grid [eV]     epsr_x  epsr_y  epsr_z')" )
      WRITE(40, "(2x,'# energy grid [eV]     epsi_x  epsi_y  epsi_z')" )
      WRITE(41, "(2x,'# energy grid [eV]  eels components [arbitrary units]')" )
      WRITE(42, "(2x,'# energy grid [eV]     ieps_x  ieps_y  ieps_z ')" )
      WRITE(43, "(2x,'# energy grid [eV]     epsr_x  epsr_y  epsr_z')" )
      WRITE(44, "(2x,'# energy grid [eV]     epsi_x  epsi_y  epsi_z')" )
      WRITE(45, "(2x,'# energy grid [eV]  eels components [arbitrary units]')" )
      WRITE(46, "(2x,'# energy grid [eV]     ieps_x  ieps_y  ieps_z ')" )
      WRITE(47, "(2x,'# energy grid [eV]     epsr_x  epsr_y  epsr_z')" )
      WRITE(48, "(2x,'# energy grid [eV]     epsi_x  epsi_y  epsi_z')" )
      WRITE(49, "(2x,'# energy grid [eV]  eels components [arbitrary units]')" )
      WRITE(50, "(2x,'# energy grid [eV]     ieps_x  ieps_y  ieps_z ')" )
      !
      DO iw =1, nw
          !
          WRITE(30,"(4f15.6)") wgrid(iw), epsrc(0,1:3, iw)
          WRITE(40,"(4f15.6)") wgrid(iw), epsic(0,1:3, iw)
          WRITE(41,"(4f15.6)") wgrid(iw), eelsc(0,1:3, iw)
          WRITE(42,"(4f15.6)") wgrid(iw), iepsc(0,1:3, iw)
          WRITE(43,"(4f15.6)") wgrid(iw), epsrc(1,1:3, iw)
          WRITE(44,"(4f15.6)") wgrid(iw), epsic(1,1:3, iw)
          WRITE(45,"(4f15.6)") wgrid(iw), eelsc(1,1:3, iw)
          WRITE(46,"(4f15.6)") wgrid(iw), iepsc(1,1:3, iw)
          WRITE(47,"(4f15.6)") wgrid(iw), epsrc(1,1:3, iw)+epsrc(0,1:3, iw)
          WRITE(48,"(4f15.6)") wgrid(iw), epsic(1,1:3, iw)+epsic(0,1:3, iw)
          WRITE(49,"(4f15.6)") wgrid(iw), eelsc(1,1:3, iw)+eelsc(0,1:3, iw)
          WRITE(50,"(4f15.6)") wgrid(iw), iepsc(1,1:3, iw)+iepsc(0,1:3, iw)
          !
      ENDDO
      !
      CLOSE(30)
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      CLOSE(43)
      CLOSE(44)
      CLOSE(45)
      CLOSE(46)
      CLOSE(47)
      CLOSE(48)
      CLOSE(49)
      CLOSE(50)
      !
  ENDIF
  DEALLOCATE ( epsrc, epsic, eelsc, iepsc)
ENDIF
  !
  ! local cleaning
  !
  CALL grid_destroy()
  !
  DEALLOCATE (  dipole, dipole_aux )

END SUBROUTINE eps_calc
 
!----------------------------------------------------------------------------------------
SUBROUTINE jdos_calc ( smeartype,intersmear,nw,wmax,wmin,nbndmin,nbndmax,shift,nspin )
  !--------------------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE wvfct,                ONLY : nbnd, et
  USE klist,                ONLY : nks
  USE io_global,            ONLY : ionode, stdout
  USE grid_module,          ONLY : alpha, focc, wgrid, grid_build, grid_destroy                
  ! 
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(IN) :: nw,nbndmin,nbndmax,nspin
  REAL(DP),        INTENT(IN) :: wmax, wmin, intersmear, shift
  CHARACTER(*),    INTENT(IN) :: smeartype
  !
  ! local variables
  !
  INTEGER       :: ik, is, iband1, iband2 
  INTEGER       :: iw, ierr
  REAL(DP)      :: etrans, w, renorm, count, srcount(0:1), renormzero,renormuno
  !
  REAL(DP), ALLOCATABLE    :: jdos(:),srjdos(:,:)
!
!--------------------------
! main routine body
!--------------------------
!
! No wavefunctions are needed in order to compute jdos, only eigenvalues, 
! they are distributed to each task so
! no mpi calls are necessary in this routine
  !
  ! perform some consistency checks, calculate occupation numbers and setup w grid
  !
  CALL grid_build(nw, wmax, wmin )

!
! spin unresolved calculation
!
IF (nspin == 1) THEN 
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( jdos(nw), STAT=ierr )
      IF (ierr/=0) CALL errore('epsilon','allocating jdos',ABS(ierr))
  !
  ! initialize jdos
  !
  jdos(:)=0.0_DP

  ! Initialising a counter for the number of transition
  count=0.0_DP

  !
  ! main kpt loop
  !

  IF (smeartype=='lorentz') THEN

    kpt_lor: &
    DO ik = 1, nks
       !
       ! Calculation of joint density of states
       ! 'intersmear' is the brodening parameter
       !
       DO iband2 = 1,nbnd
           IF ( focc(iband2,ik) <  2.0d0) THEN
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV  + shift
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 count = count + (focc(iband1,ik)-focc(iband2,ik))
                 !
                 ! loop over frequencies
                 !
                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     jdos(iw) = jdos(iw) + intersmear * (focc(iband1,ik)-focc(iband2,ik)) &
                                  / ( PI * ( (etrans -w )**2 + (intersmear)**2 ) )

                 ENDDO

           ENDIF
       ENDDO
           ENDIF
       ENDDO

    ENDDO kpt_lor

  ELSE IF (smeartype=='gauss') THEN

    kpt_gauss: &
    DO ik = 1, nks

       !
       ! Calculation of joint density of states
       ! 'intersmear' is the brodening parameter
       !
       DO iband2 = 1,nbnd
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband2,ik) <  2.0d0) THEN
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV  + shift
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 ! loop over frequencies
                 !

                 count=count+ (focc(iband1,ik)-focc(iband2,ik))

                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     jdos(iw) = jdos(iw) + (focc(iband1,ik)-focc(iband2,ik)) * &
                                EXP(-(etrans-w)**2/intersmear**2) &
                                  / (intersmear * SQRT(PI))

                 ENDDO

           ENDIF
           ENDIF
       ENDDO
       ENDDO

    ENDDO kpt_gauss

  ELSE

    CALL errore('epsilon', 'invalid SMEARTYPE = '//TRIM(smeartype), 1)

  ENDIF

  !
  ! jdos normalizzation
  !

  jdos(:)=jdos(:)/count

  !
  ! check jdos normalization
  !

  renorm = alpha * SUM( jdos(:) )
  !
  ! write results on data files
  !
  IF (ionode) THEN
     WRITE(stdout,"(/,5x, 'Integration over JDOS gives: ',f15.9,' instead of 1.0d0' )") renorm
     WRITE(stdout,"(/,5x, 'Writing output on file...' )")
                                               
     OPEN (30, FILE='jdos.dat', FORM='FORMATTED' )
     !
     WRITE(30, "(2x,'# energy grid [eV]     JDOS [1/eV] ')" )
     !
     DO iw =1, nw
         !
         WRITE(30,"(4f15.6)") wgrid(iw), jdos(iw)
         !
     ENDDO
     !
     CLOSE(30)
  ENDIF
  !
  ! local cleaning
  !
  DEALLOCATE ( jdos )

!
! collinear spin calculation
!
ELSE IF(nspin==2) THEN
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( srjdos(0:1,nw), STAT=ierr )
      IF (ierr/=0) CALL errore('epsilon','allocating spin resolved jdos',ABS(ierr))
  !
  ! initialize jdos
  !
  srjdos(:,:)=0.0_DP

  ! Initialising a counter for the number of transition
  srcount(:)=0.0_DP

  !
  ! main kpt loop
  !

  IF (smeartype=='lorentz') THEN

  DO is=0,1
    ! if nspin=2 the number of nks must be even (even if the calculation 
    ! is performed at gamma point only), so nks must be always a multiple of 2
    DO ik = 1 + is * INT(nks/2), INT(nks/2) +  is * INT(nks/2)
       !
       ! Calculation of joint density of states
       ! 'intersmear' is the brodening parameter
       !
       DO iband2 = 1,nbnd
           IF ( focc(iband2,ik) <  2.0d0) THEN
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV  + shift
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 ! loop over frequencies
                 !
                 srcount(is)=srcount(is)+ (focc(iband1,ik)-focc(iband2,ik))

                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     srjdos(is,iw) = srjdos(is,iw) + intersmear * (focc(iband1,ik)-focc(iband2,ik)) &
                                  / ( PI * ( (etrans -w )**2 + (intersmear)**2 ) )

                 ENDDO

           ENDIF
       ENDDO
           ENDIF
       ENDDO

    ENDDO 
 ENDDO

  ELSE IF (smeartype=='gauss') THEN

  DO is=0,1
    ! if nspin=2 the number of nks must be even (even if the calculation 
    ! is performed at gamma point only), so nks must be always a multiple of 2
    DO ik = 1 + is * INT(nks/2), INT(nks/2) +  is * INT(nks/2)
       !
       ! Calculation of joint density of states
       ! 'intersmear' is the brodening parameter
       !
       DO iband2 = 1,nbnd
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband2,ik) <  2.0d0) THEN
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV  + shift
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 ! loop over frequencies
                 !

                 srcount(is)=srcount(is)+ (focc(iband1,ik)-focc(iband2,ik))

                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     srjdos(is,iw) = srjdos(is,iw) + (focc(iband1,ik)-focc(iband2,ik)) * &
                                EXP(-(etrans-w)**2/intersmear**2) &
                                  / (intersmear * SQRT(PI))

                 ENDDO

           ENDIF
           ENDIF
       ENDDO
       ENDDO

    ENDDO 
 ENDDO

  ELSE

    CALL errore('epsilon', 'invalid SMEARTYPE = '//TRIM(smeartype), 1)

  ENDIF

  !
  ! jdos normalizzation
  !
  DO is = 0,1
    srjdos(is,:)=srjdos(is,:)/srcount(is)
  ENDDO
  !
  ! check jdos normalization
  !

  renormzero = alpha * SUM( srjdos(0,:) )
  renormuno = alpha * SUM( srjdos(1,:) )
  !
  ! write results on data files
  !
  IF (ionode) THEN
     WRITE(stdout,"(/,5x, 'Integration over spin UP JDOS gives: ',f15.9,' instead of 1.0d0' )") renormzero
     WRITE(stdout,"(/,5x, 'Integration over spin DOWN JDOS gives: ',f15.9,' instead of 1.0d0' )") renormuno
     WRITE(stdout,"(/,5x, 'Writing output on file...' )")
                                               
     OPEN (30, FILE='jdos.dat', FORM='FORMATTED' )
     !
     WRITE(30, "(2x,'# energy grid [eV]     UJDOS [1/eV]      DJDOS[1:eV]')" )
     !
     DO iw =1, nw
         !
         WRITE(30,"(4f15.6)") wgrid(iw), srjdos(0,iw), srjdos(1,iw)
         !
     ENDDO
     !
     CLOSE(30)
  ENDIF

  DEALLOCATE ( srjdos )
ENDIF
  !
  ! local cleaning
  !
  CALL grid_destroy()

END SUBROUTINE jdos_calc

!-----------------------------------------------------------------------------
SUBROUTINE offdiag_calc ( intersmear,intrasmear, nw, wmax, wmin, nbndmin, nbndmax,&
                          shift, metalcalc, nspin )
  !-----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE cell_base,            ONLY : tpiba2, omega
  USE wvfct,                ONLY : nbnd, et
  USE ener,                 ONLY : efermi => ef
  USE klist,                ONLY : nks, nkstot, degauss
  USE grid_module,          ONLY : focc, wgrid, grid_build, grid_destroy
  USE io_global,            ONLY : ionode, stdout
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum

  ! 
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(IN) :: nw,nbndmin,nbndmax,nspin
  REAL(DP),        INTENT(IN) :: wmax, wmin, intersmear,intrasmear, shift
  LOGICAL,         INTENT(IN) :: metalcalc 
  !
  ! local variables
  !
  INTEGER       :: ik, iband1, iband2 
  INTEGER       :: iw, ierr, it1, it2
  REAL(DP)      :: etrans, const, w
  !
  COMPLEX(DP), ALLOCATABLE :: dipole_aux(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: epstot(:,:,:),dipoletot(:,:,:,:)
  !
  !--------------------------
  ! main routine body
  !--------------------------
  !
  ! perform some consistency checks, calculate occupation numbers and setup w grid
  !
  CALL grid_build(nw, wmax, wmin )
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( dipoletot(3,3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipoletot', ABS(ierr) )
  !
  ALLOCATE( dipole_aux(3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole_aux', ABS(ierr) )
  !
  ALLOCATE(epstot( 3,3, nw),STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating epstot', ABS(ierr))

   !
   ! initialize response functions
   !
   epstot  = (0.0_DP,0.0_DP)
   !
   ! main kpt loop
   !
   DO ik = 1, nks
     !
     ! For every single k-point: order k+G for
     !                           read and distribute wavefunctions
     !                           compute dipole matrix 3 x nbnd x nbnd parallel over g 
     !                           recover g parallelism getting the total dipole matrix
     ! 
     CALL dipole_calc( ik, dipole_aux, metalcalc, nbndmin, nbndmax)
     !
     DO it2 = 1, 3
        DO it1 = 1, 3
           dipoletot(it1,it2,:,:) = tpiba2 * dipole_aux(it1,:,:) * CONJG( dipole_aux(it2,:,:) )
        ENDDO
     ENDDO
     !
     ! Calculation of real and immaginary parts
     ! of the macroscopic dielettric function from dipole
     ! approximation.
     ! 'intersmear' is the brodening parameter
     !
     DO iband2 = 1,nbnd
         IF ( focc(iband2,ik) <  2.0d0) THEN
     DO iband1 = 1,nbnd
         !
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
             !
             ! transition energy
             !
             etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV + shift
             !
             IF (ABS(focc(iband2,ik)-focc(iband1,ik))< 1e-4) CYCLE
             !
             ! loop over frequencies
             !
             DO iw = 1, nw
                  !
                  w = wgrid(iw)
                  !
                  epstot(:,:,iw) = epstot(:,:,iw) + dipoletot(:,:,iband1,iband2)*RYTOEV**3/(etrans) *&
                                   focc(iband1,ik)/(etrans**2 - w**2 - (0,1)*intersmear*w)
             ENDDO
             !
         ENDIF
     ENDDO
         ENDIF
     ENDDO
     !
     !Intraband (only if metalcalc is true)
     !
     IF (metalcalc) THEN
     DO iband1 = 1,nbnd
         !
         IF ( focc(iband1,ik) < 2.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
               !
               ! loop over frequencies
               !
               DO iw = 1, nw
                   !
                   w = wgrid(iw)
                   !
                   epstot(:,:,iw) = epstot(:,:,iw) - dipoletot(:,:,iband1,iband1)* &
                                RYTOEV**2 * (EXP((et(iband1,ik)-efermi)/degauss ))/  &
                    (( w**2 + (0,1)*intrasmear*w)*(1+EXP((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss ) 
               ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDIF
  ENDDO

  !
  ! recover over kpt parallelization (inter_pool)
  !
  CALL mp_sum( epstot, inter_pool_comm )

  !
  ! impose the correct normalization
  !
  const = 64.0d0 * PI / ( omega * REAL(nkstot, DP) )
  !
  epstot(:,:,:) = 1.0_DP + epstot(:,:,:) * const
  !
  ! write results on data files
  !
  IF (ionode) THEN
      !
      WRITE(stdout,"(/,5x, 'Writing output on file...' )")
      !
      OPEN (41, FILE='epsxx.dat', FORM='FORMATTED' )
      OPEN (42, FILE='epsxy.dat', FORM='FORMATTED' )
      OPEN (43, FILE='epsxz.dat', FORM='FORMATTED' )
      OPEN (44, FILE='epsyx.dat', FORM='FORMATTED' )
      OPEN (45, FILE='epsyy.dat', FORM='FORMATTED' )
      OPEN (46, FILE='epsyz.dat', FORM='FORMATTED' )
      OPEN (47, FILE='epszx.dat', FORM='FORMATTED' )
      OPEN (48, FILE='epszy.dat', FORM='FORMATTED' )
      OPEN (49, FILE='epszz.dat', FORM='FORMATTED' )
      ! 
      WRITE(41, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(42, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(43, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(44, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(45, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(46, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(47, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(48, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(49, "(2x,'# energy grid [eV]     epsr     epsi')" )
      !
      DO iw =1, nw
         !
         WRITE(41,"(4f15.6)") wgrid(iw), REAL(epstot(1,1, iw)), AIMAG(epstot(1,1, iw))
         WRITE(42,"(4f15.6)") wgrid(iw), REAL(epstot(1,2, iw)), AIMAG(epstot(1,2, iw))
         WRITE(43,"(4f15.6)") wgrid(iw), REAL(epstot(1,3, iw)), AIMAG(epstot(1,3, iw))
         WRITE(44,"(4f15.6)") wgrid(iw), REAL(epstot(2,1, iw)), AIMAG(epstot(2,1, iw))
         WRITE(45,"(4f15.6)") wgrid(iw), REAL(epstot(2,2, iw)), AIMAG(epstot(2,2, iw))
         WRITE(46,"(4f15.6)") wgrid(iw), REAL(epstot(2,3, iw)), AIMAG(epstot(2,3, iw))
         WRITE(47,"(4f15.6)") wgrid(iw), REAL(epstot(3,1, iw)), AIMAG(epstot(3,1, iw))
         WRITE(48,"(4f15.6)") wgrid(iw), REAL(epstot(3,2, iw)), AIMAG(epstot(3,2, iw))
         WRITE(49,"(4f15.6)") wgrid(iw), REAL(epstot(3,3, iw)), AIMAG(epstot(3,3, iw))
         !
      ENDDO
      !
      CLOSE(30)
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      !
  ENDIF

  !
  ! local cleaning
  !
  CALL grid_destroy()
  DEALLOCATE ( dipoletot, dipole_aux, epstot )

END SUBROUTINE offdiag_calc


!--------------------------------------------------------------------
SUBROUTINE dipole_calc( ik, dipole_aux, metalcalc, nbndmin, nbndmax )
  !------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : npw, nbnd, igk, g2kin
  USE wavefunctions_module, ONLY : evc
  USE klist,                ONLY : xk
  USE cell_base,            ONLY : tpiba2
  USE gvect,                ONLY : ngm, g, ecutwfc
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE grid_module,          ONLY : focc
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum

IMPLICIT NONE
  !
  ! global variables
  INTEGER, INTENT(IN)        :: ik,nbndmin,nbndmax
  COMPLEX(DP), INTENT(INOUT) :: dipole_aux(3,nbnd,nbnd)
  LOGICAL, INTENT(IN)        :: metalcalc
  !
  ! local variables
  INTEGER :: iband1,iband2,ig
  COMPLEX(DP)   :: caux 

  !
  ! Routine Body
  ! 
  CALL start_clock( 'dipole_calc' )

  !
  ! setup k+G grids for each kpt
  !
  CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)  
  ! 
  ! read wfc for the given kpt
  !
  CALL davcio (evc, nwordwfc, iunwfc, ik, - 1) 
  !
  ! compute matrix elements
  !
  dipole_aux(:,:,:) = (0.0_DP,0.0_DP)
  !
  DO iband2 = nbndmin,nbndmax
      IF ( focc(iband2,ik) <  2.0d0) THEN
  DO iband1 = nbndmin,nbndmax 
      !
      IF ( iband1==iband2 ) CYCLE 
      IF ( focc(iband1,ik) >= 1e-4 ) THEN
            !
            DO  ig=1,npw
                 !
                 caux= CONJG(evc(ig,iband1))*evc(ig,iband2) 
                 !
                 dipole_aux(:,iband1,iband2) = dipole_aux(:,iband1,iband2) + &
                       ( g(:,igk(ig)) ) * caux
                 !
            ENDDO
      ENDIF
      !
  ENDDO
      ENDIF
  ENDDO
  !
  ! The diagonal terms are taken into account only if the system is treated like a metal, not
  ! in the intraband therm. Because of this we can recalculate the diagonal component of the dipole 
  ! tensor directly as we need it for the intraband therm, without interference with interband one.
  !
  IF (metalcalc) THEN
     !
     DO iband1 = nbndmin,nbndmax
        DO  ig=1,npw
          !
          caux= CONJG(evc(ig,iband1))*evc(ig,iband1) 
          !
          dipole_aux(:,iband1,iband1) = dipole_aux(:,iband1,iband1) + &
                                        ( g(:,igk(ig))+ xk(:,ik) ) * caux
          !
        ENDDO
     ENDDO
     !
  ENDIF
  !
  ! recover over G parallelization (intra_pool)
  !
  CALL mp_sum( dipole_aux, intra_pool_comm ) 
  !
  CALL stop_clock( 'dipole_calc' )
  !
END SUBROUTINE dipole_calc


!-------------------------------------------------
SUBROUTINE occ_calc ()
  !-------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE klist,     ONLY : nkstot, wk, degauss
  USE wvfct,     ONLY : nbnd, wg, et
  USE ener,      ONLY : ef
  USE mp_global, ONLY : me_pool
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE  :: focc(:,:),foccp(:,:)
  CHARACTER(25)          :: filename 
  INTEGER                :: ierr, i, ik
  ! 
  ALLOCATE ( focc( nbnd, nkstot), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating focc', ABS(ierr))
  !
  ALLOCATE ( foccp( nbnd, nkstot), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating foccp', ABS(ierr))

  IF (me_pool==0) THEN
      !
      filename = 'occupations.dat'
!      WRITE(filename,"(I3,'.occupation.dat')")me_pool
      OPEN (unit=50, file=TRIM(filename))
      WRITE(50,*) '#energy (Ry)      occupation factor       derivative'
    
      DO ik = 1,nkstot
      DO i  = 1,nbnd
           focc(i,ik)= wg(i, ik ) * 2.0_DP/wk( ik )
           foccp(i,ik)= 2* EXP((et(i,ik)-ef)/degauss)/((1+EXP((et(i,ik)-ef)/degauss))**2*degauss)
           WRITE(50,*)et(i,ik),focc(i,ik),foccp(i,ik) 
      ENDDO
      ENDDO
    
      CLOSE (50) 
      !
  ENDIF
  !
  DEALLOCATE ( focc, STAT=ierr)
  CALL errore('grid_destroy','deallocating grid stuff',ABS(ierr))
  !
  DEALLOCATE ( foccp, STAT=ierr)
  CALL errore('grid_destroy','deallocating grid stuff',ABS(ierr))

END SUBROUTINE occ_calc

!----------------------------------------------------------------------------------------
SUBROUTINE photoemission_spectr_pw ( smeartype,intersmear,intrasmear,nw,wmax,wmin,nbndmin,nbndmax, &
                                     shift,nspin,metalcalc,polar_angle, azimuthal_angle,&
                                     photon_angle,photon_ener,homo_gas,wfc_real,e_fermi,&
                                     modified_pw,othor_pw)
  !--------------------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE wvfct,                ONLY : npw, nbnd, igk
  USE wavefunctions_module, ONLY : evc
  USE cell_base,            ONLY : tpiba2
  USE wvfct,                ONLY : nbnd, et
  USE gvect,                ONLY : ngm, g, ecutwfc
  USE klist,                ONLY : nks, xk
  USE io_global,            ONLY : ionode, stdout
  USE grid_module,          ONLY : alpha, focc, wgrid, grid_build, grid_destroy               
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
  ! 
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(IN) :: nw,nbndmin,nbndmax,nspin
  REAL(DP),        INTENT(IN) :: wmax, wmin, intersmear, shift, photon_ener,  & 
                                 polar_angle, azimuthal_angle, photon_angle,  &
                                 e_fermi, intrasmear
  CHARACTER(*),    INTENT(IN) :: smeartype
  LOGICAL,         INTENT(IN) :: metalcalc, homo_gas, wfc_real, modified_pw,  &
                                 othor_pw
  !
  ! local variables
  !
  INTEGER       :: ik, is, iband1, iband2, ipol, ig
  INTEGER       :: iw, ierr, ie
  REAL(DP)      :: etrans, w, renorm, count, srcount(0:1), renormzero,renormuno
  REAL(DP)      :: polar_angle_radial, azimuthal_angle_radial, photon_angle_radial
  REAL(DP)      :: ekin_eout, ekin_eout_list (nbndmax)
  REAL(DP)      :: ekin_plus, ekin_minus, module_k, kx, ky, kz, delta_ecut_G, &
                   delta_kx_Gx, delta_ky_Gy, delta_kz_Gz, max_num, constant,  & 
                   max_sigma_tot
  REAL(DP)      :: sigma_tot (nbndmax)
  REAL(DP)      :: ssigma_tot (2, nbndmax)
  REAL(DP)      :: sigma_2   (nbndmax)
  REAL(DP)      :: ssigma_2   (2,nbndmax)
  REAL(DP)      :: beta      (nbndmax)
  REAL(DP)      :: sbeta      (2,nbndmax)
  REAL(DP)      :: eigen_mol (nbndmax)
  REAL(DP)      :: seigen_mol (2,nbndmax)
  !
  REAL(DP), ALLOCATABLE    :: photospec(:), srphotospec(:,:)
  REAL(DP), ALLOCATABLE    :: g2kin (:)
  REAL(DP), ALLOCATABLE    :: dipole_2(:,:)
  REAL(DP), ALLOCATABLE    :: gamma_2_opw(:,:)
  REAL(DP), ALLOCATABLE    :: lambda_2_opw(:,:)
  COMPLEX(DP),ALLOCATABLE  :: dipole_aux(:,:,:)
  COMPLEX(DP),ALLOCATABLE  :: dipole_opw_gamma(:,:,:)
  COMPLEX(DP),ALLOCATABLE  :: scala_opw_lambda(:,:)
  !
  !--------------------------
  ! main routine body
  !--------------------------
  !
  ! No wavefunctions are needed in order to compute jdos, only eigenvalues, 
  ! they are distributed to each task so
  ! no mpi calls are necessary in this routine
  !
  ! perform some consistency checks, calculate occupation numbers and setup w
  ! grid
  !
  CALL grid_build(nw, wmax, wmin )
  !
  ! change angles from degree to radial unit
  !
  polar_angle_radial = (polar_angle/180.0)*PI
  azimuthal_angle_radial = (azimuthal_angle/180.0)*PI
  photon_angle_radial = (photon_angle/180.0)*PI
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( g2kin(npw), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating photoemission_spectr',ABS(ierr))
  !
  ALLOCATE( dipole_2(nbndmax, npw), STAT=ierr ) 
  IF (ierr/=0) CALL errore('epsilon','allocating dipole', ABS(ierr) )
  !
  ALLOCATE( dipole_aux(3, nbndmax, npw), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole_aux', ABS(ierr) )
  !
  IF (othor_pw) THEN
    ALLOCATE( dipole_opw_gamma(3, nbndmax, npw) )
    ALLOCATE( scala_opw_lambda(nbndmax, npw) )
    ALLOCATE( gamma_2_opw(nbndmax, npw) )
    ALLOCATE( lambda_2_opw(nbndmax, npw) )
  ENDIF
  !
  ! spin unresolved calculation
  !
IF (nspin == 1) THEN
  !
  ALLOCATE( photospec(nw), STAT=ierr )
      IF (ierr/=0) CALL errore('epsilon','allocating photoemission_spectr',ABS(ierr))
  !
  ! initialize  photoemission
  !
  constant = 1.0_DP
  !
  photospec(:)=0.0_DP
  !
  ! Initialising a counter for the number of transition
  count=0.0_DP
  !
  ! main kpt loop
  !
  IF (smeartype=='lorentz') THEN
    !
    kpt_lor: &
    DO ik = 1, nks
       !
       ! For every single k-point: order k+G for
       !                           read and distribute wavefunctions
       !                           compute dipole matrix 3 x nbnd x nbnd
       !                           parallel over g 
       !                           recover g parallelism getting the total
       !                           dipole matrix
       !
       CALL dipole_calc_pw( ik, dipole_aux, dipole_opw_gamma, scala_opw_lambda, metalcalc, nbndmin, nbndmax, &
                            photon_angle_radial, azimuthal_angle_radial, &
                            homo_gas, othor_pw)
       !
       dipole_2(:,:) = 0.0_DP 
       DO ipol = 1, 3  
             dipole_2(:,:) = dipole_2(:,:) + tpiba2 * &
                           ( (REAL(dipole_aux(ipol,:,:)))**2 + (AIMAG(dipole_aux(ipol,:,:)))**2 )  
       ENDDO
       !
       IF (othor_pw) THEN
          ! 
          gamma_2_opw(:,:) = 0.0_DP
          lambda_2_opw(:,:) = 0.0_DP
          DO ipol = 1, 3   
             gamma_2_opw(:,:) = gamma_2_opw(:,:) + tpiba2 * & 
                          ( (REAL(dipole_opw_gamma(ipol,:,:)))**2 + (AIMAG(dipole_opw_gamma(ipol,:,:)))**2 )
          ENDDO
          !
          lambda_2_opw(:,:) = lambda_2_opw(:,:)  + tpiba2 * & 
                            scala_opw_lambda(:,:) * conjg(scala_opw_lambda(:,:)) 
          !
       ENDIF
       !
       ! Calculation of photoemission spectra
       ! 'intersmear' is the brodening parameter
       !
       IF (homo_gas) THEN
         sigma_tot(:) = 0.0_DP
         sigma_2(:)   = 0.0_DP
       ENDIF
       !
       DO iband1 = nbndmin, nbndmax 
          !
          IF (homo_gas) THEN   
             !
             eigen_mol(iband1) = (0.0_DP - et(iband1, ik) )*RYTOEV
             !
             IF (modified_pw ) THEN
                ekin_eout = photon_ener
             ELSE   
                ekin_eout = photon_ener - (0.0_DP - et(iband1, ik) )*RYTOEV
             ENDIF
             !  
          ELSE
             ekin_eout = photon_ener -((0.0_DP - et(iband1, ik) )*RYTOEV + e_fermi)
          ENDIF
          ! 
          IF (ekin_eout > photon_ener .or. ekin_eout <=0.0_DP ) CYCLE
          !  
          module_k = sqrt ((ekin_eout/13.6056923)/tpiba2)
          ! 
          kx = module_k * cos(azimuthal_angle_radial) * sin(polar_angle_radial) 
          ky = module_k * sin(azimuthal_angle_radial) * sin(polar_angle_radial)
          kz = module_k * cos(polar_angle_radial)  
          !  
          ! compute the total cross-section
          ! and ansymmetry parameter
          !
          IF (homo_gas) THEN
            !
            DO ig = 1, npw
              !
              g2kin(ig) = tpiba2*((xk(1,ik)+g(1,igk(ig)))**2 + (xk(2,ik)+g(2,igk(ig))) **2 + (xk(3,ik)+g(3,igk(ig)))**2)
              delta_ecut_G =  intrasmear/( PI * (( g2kin(ig)*13.6056923 - ekin_eout )**2 + (intrasmear)**2 )) 
              !
              sigma_tot(iband1) = sigma_tot(iband1) + constant*module_k*dipole_2(iband1, ig)*delta_ecut_G
              !
              IF (othor_pw ) THEN
                sigma_2(iband1) = sigma_2(iband1)   + 1.5*constant*module_k*(gamma_2_opw(iband1, ig) - lambda_2_opw(iband1, ig))*delta_ecut_G
              ENDIF
              !
            ENDDO
            !
          ENDIF
          ! 
          DO ig = 1, npw
             !
             g2kin(ig) = tpiba2*((xk(1,ik)+g(1,igk(ig)))**2 + (xk(2,ik)+g(2,igk(ig))) **2 + (xk(3,ik)+g(3,igk(ig)))**2)
             delta_ecut_G =  intrasmear/( PI * (( g2kin(ig)*13.6056923 - ekin_eout  )**2 + (intrasmear)**2 )) 
             ! 
             delta_kx_Gx  =  intrasmear/( PI * (( ( kx - ( xk(1,ik) + g(1,igk(ig)))) *sqrt(tpiba2))**2 + (intrasmear)**2 ) )
             delta_ky_Gy  =  intrasmear/( PI * (( ( ky - ( xk(2,ik) + g(2,igk(ig)))) *sqrt(tpiba2))**2 + (intrasmear)**2 ) )
             delta_kz_Gz  =  intrasmear/( PI * (( ( kz - ( xk(3,ik) + g(3,igk(ig)))) *sqrt(tpiba2))**2 + (intrasmear)**2 ) )
             !
             ! transition energy 
             !
             etrans = (0.0d0 - et(iband1, ik) ) * RYTOEV  + shift ! g2kin were called in the dipole_pw routine
             !
             IF( etrans < 1.0d-10 ) CYCLE
             !
             ! loop over frequencies
             !
             DO iw = 1, nw
                !
                w = wgrid(iw)
                !
                IF (homo_gas) THEN
                   ! 
                   photospec(iw) = photospec(iw) + module_k * intersmear * dipole_2(iband1, ig) &
                                     / ( PI * ( (etrans - w )**2 + (intersmear)**2 )  ) * delta_ecut_G 
                   !
                   IF (wfc_real) then 
                      !
                      photospec(iw) = photospec(iw) + 2.0D0 * module_k * intersmear * dipole_2(iband1, ig) & 
                                     / ( PI * ( (etrans - w )**2 + (intersmear)**2 )  ) * delta_ecut_G
                      !
                   ENDIF
                   !
                ELSE 
                   !
                   photospec(iw) = photospec(iw) + dipole_2(iband1, ig) * delta_ecut_G &
                                 * delta_kx_Gx * delta_ky_Gy * delta_kz_Gz &
                                 * intersmear / ( PI * ( (etrans - w )**2 + (intersmear)**2 )  ) 
                   !
                ENDIF 
                ! 
             ENDDO
             ! 
          ENDDO 
          !
       ENDDO
       !  
    ENDDO kpt_lor
    !
  ENDIF
  !
  ! recover over G parallelization (intra_pool)
  !
  CALL mp_sum( photospec(:), intra_pool_comm )
  !
  IF(homo_gas) THEN
    !  
    CALL mp_sum( sigma_tot(:), intra_pool_comm )
    CALL mp_sum( sigma_2(:), intra_pool_comm )
    max_sigma_tot = maxval(sigma_tot(:))
    !
    IF (othor_pw) THEN
       DO iband1 = nbndmin, nbndmax
         beta(iband1)=2.0_DP*(1.0_DP-sigma_2(iband1)/sigma_tot(iband1) )
       ENDDO 
    ELSE
       beta(:) = 2.0_DP 
    ENDIF
    !
    IF (ionode) THEN
      WRITE(stdout,"(/,5x, 'Writing the molecule gas properties' )")
      DO iband1 = nbndmin, nbndmax
         WRITE(stdout,"(4f15.6)") eigen_mol(iband1), sigma_tot(iband1)/max_sigma_tot, beta(iband1)
      ENDDO
    ENDIF 
    !
  ENDIF
  !
  ! write results on data files
  !
  max_num = maxval(photospec(:))
  !
  IF (ionode) THEN
     WRITE(stdout,"(/,5x, 'Writing output on file...' )")
                                               
     OPEN (30, FILE='photospectra.dat', FORM='FORMATTED' )
     !
     WRITE(30, "(2x,'# energy grid [eV]     photospectra [ab.] ')" )
     !
     DO iw =1, nw
         !
         WRITE(30,"(4f15.6)") wgrid(iw), photospec(iw), photospec(iw)/max_num
         !
     ENDDO
     !
     CLOSE(30)
  ENDIF
  !
  ! local cleaning
  !
  DEALLOCATE ( g2kin, photospec , dipole_2, dipole_aux)
  !
ELSEIF (nspin==2) THEN
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( srphotospec(0:1,nw), STAT=ierr )
      IF (ierr/=0) CALL errore('epsilon','allocating spin resolved photospectra',ABS(ierr))
  !
  ! initialize 
  !
  srphotospec(:,:)=0.0_DP
  !
  ! Initialising a counter for the number of transition
  srcount(:)=0.0_DP
  !
  constant = 1.0_DP
  !
  IF (homo_gas) THEN
    ssigma_tot(:, :) = 0.0_DP
    ssigma_2(:, :)   = 0.0_DP
  ENDIF
  !
  ! main kpt loop
  !
  IF (smeartype=='lorentz') THEN
  !
  DO is = 0, 1
     !
     ! if nspin=2 the number of nks must be even (even if the calculation 
     ! is performed at gamma point only), so nks must be always a multiple of 2
     ! 
     DO ik = 1 + is * INT(nks/2), INT(nks/2) +  is * INT(nks/2)
        !
	! For every single k-point: order k+G for
	!                           read and distribute wavefunctions
	!                           compute dipole matrix 3 x nbnd x nbnd
	!                           parallel over g 
	!                           recover g parallelism getting the total
	!                           dipole matrix
	!
        CALL dipole_calc_pw( ik, dipole_aux, dipole_opw_gamma, scala_opw_lambda, metalcalc, nbndmin, nbndmax, &
                             photon_angle_radial, azimuthal_angle_radial, &
                             homo_gas, othor_pw)
	!
	dipole_2(:,:) = 0.0_DP
	DO ipol = 1, 3
	    dipole_2(:,:) = dipole_2(:,:) + tpiba2 * ( (REAL(dipole_aux(ipol,:,:)))**2 + (AIMAG(dipole_aux(ipol,:,:)))**2 )
	ENDDO
	!
	IF (othor_pw) THEN
	  ! 
	  gamma_2_opw(:,:) = 0.0_DP
	  lambda_2_opw(:,:) = 0.0_DP
	  DO ipol = 1, 3
	    gamma_2_opw(:,:) = gamma_2_opw(:,:) + tpiba2 * & 
			( (REAL(dipole_opw_gamma(ipol,:,:)))**2 + (AIMAG(dipole_opw_gamma(ipol,:,:)))**2 )
	  ENDDO
	  !
	  lambda_2_opw(:,:) = lambda_2_opw(:,:)  + tpiba2 * & 
			    scala_opw_lambda(:,:) * conjg(scala_opw_lambda(:,:))
	  !
	ENDIF
	!
	! Calculation of photoemission spectra
	! 'intersmear' is the brodening parameter
	!
	DO iband1 = nbndmin, nbndmax 
	  !
	  IF (homo_gas) THEN
	    !
	    seigen_mol(is, iband1) = (0.0_DP - et(iband1, ik) )*RYTOEV
	    ! 
	    IF (modified_pw ) THEN
		ekin_eout = photon_ener
	    ELSE
		ekin_eout = photon_ener - (0.0_DP - et(iband1, ik) )*RYTOEV
	    ENDIF
	    !  
	  ELSE
            ! 
	    ekin_eout = photon_ener -((0.0_DP - et(iband1, ik) )*RYTOEV + e_fermi)
            !
	  ENDIF
	  ! 
	  IF (ekin_eout > photon_ener .or. ekin_eout <=0.0_DP ) CYCLE
	  !  
	  module_k = sqrt ((ekin_eout/13.6056923)/tpiba2)
	  ! 
	  kx = module_k * cos(azimuthal_angle_radial) * sin(polar_angle_radial)
	  ky = module_k * sin(azimuthal_angle_radial) * sin(polar_angle_radial)
	  kz = module_k * cos(polar_angle_radial)
	  !
	  ! compute the total cross-section
	  ! and ansymmetry parameter
	  !
	  IF (homo_gas) THEN
	    !
	    DO ig = 1, npw
	      !
	      g2kin(ig) = tpiba2*((xk(1,ik)+g(1,igk(ig)))**2 + (xk(2,ik)+g(2,igk(ig))) **2 + (xk(3,ik)+g(3,igk(ig)))**2)
	      delta_ecut_G =  intrasmear/( PI * (( g2kin(ig)*13.6056923 - ekin_eout )**2 + (intrasmear)**2 ))
	      !
	      ssigma_tot(is, iband1) = ssigma_tot(is, iband1) + constant*module_k*dipole_2(iband1, ig)*delta_ecut_G
	      !
	      IF (othor_pw ) THEN
		ssigma_2(is, iband1) = ssigma_2(is, iband1)   + 1.5*constant*module_k &
					* (gamma_2_opw(iband1,ig) - lambda_2_opw(iband1,ig))*delta_ecut_G 
	      ENDIF
	      !
	    ENDDO
	    !
	  ENDIF
	  ! 
	  DO ig = 1, npw
	      g2kin(ig) = tpiba2*((xk(1,ik)+g(1,igk(ig)))**2 + (xk(2,ik)+g(2,igk(ig))) **2 + (xk(3,ik)+g(3,igk(ig)))**2)
	      delta_ecut_G =  intrasmear/( PI * (( g2kin(ig)*13.6056923 - ekin_eout  )**2 + (intrasmear)**2 ))
	      ! 
	      delta_kx_Gx  =  intrasmear/( PI * (( ( kx - ( xk(1,ik) + g(1,igk(ig)))) *sqrt(tpiba2))**2 + (intrasmear)**2 ) )
	      delta_ky_Gy  =  intrasmear/( PI * (( ( ky - ( xk(2,ik) + g(2,igk(ig)))) *sqrt(tpiba2))**2 + (intrasmear)**2 ) )
	      delta_kz_Gz  =  intrasmear/( PI * (( ( kz - ( xk(3,ik) + g(3,igk(ig)))) *sqrt(tpiba2))**2 + (intrasmear)**2 ) )
	      !
	      ! transition energy 
	      !
	      etrans = (0.0d0 - et(iband1, ik) ) * RYTOEV  + shift ! g2kin were called in the dipole_pw routine
	      !
	      IF( etrans < 1.0d-10 ) CYCLE
	      !
	      ! loop over frequencies
	      !
	      DO iw = 1, nw
		!
		w = wgrid(iw)
		!
		IF (homo_gas) THEN
		  !
		  IF (wfc_real) THEN  
		    !
		    srphotospec(is, iw) = srphotospec(is, iw) + 2.0D0 * module_k* intersmear * dipole_2(iband1, ig) * delta_ecut_G &
					  / ( PI * ( (etrans - w )**2 + (intersmear)**2 )  ) 
		    !
		  ELSE
		    !
		    srphotospec(is, iw) = srphotospec(is, iw) + module_k * intersmear * dipole_2(iband1, ig) * delta_ecut_G &
					/ ( PI * ( (etrans - w )**2 + (intersmear)**2 )  ) 
		    !
		  ENDIF
		  !
		ELSE
		  !
		  srphotospec(is, iw) = srphotospec(is, iw) + dipole_2(iband1, ig) * delta_ecut_G &
				    * delta_kx_Gx * delta_ky_Gy * delta_kz_Gz * intersmear & 
				    / ( PI * ( (etrans - w )**2 + (intersmear)**2 )  )
		  !
		ENDIF
		!
	      ENDDO
	      ! 
	  ENDDO
	  !  
	ENDDO
	!
      ENDDO ! kpt_lor
      !
  ENDDO 
  !
  ENDIF
  !
  ! recover over G parallelization (intra_pool)
  !
  DO is = 0, 1
    !
    CALL mp_sum( srphotospec(is,:), intra_pool_comm )
    !
    IF(homo_gas) THEN
      !  
      CALL mp_sum( ssigma_tot(is,:), intra_pool_comm )
      CALL mp_sum( ssigma_2(is,:), intra_pool_comm )
      !
      max_sigma_tot = maxval(ssigma_tot(is, :))
      ! 
      IF (othor_pw) THEN
        !
        DO iband1 = nbndmin, nbndmax
           sbeta(is, iband1)= 2.0_DP*(1.0_DP-ssigma_2(is, iband1)/ssigma_tot(is, iband1) )
        ENDDO
        ! 
      ELSE
        ! 
        sbeta(is,:) = 2.0_DP
        !
      ENDIF
      !
      IF (ionode) THEN
        IF (is == 0) WRITE(stdout,"(/,5x, 'Writing the molecule gas properties with spin up' )") 
        IF (is == 1) WRITE(stdout,"(/,5x, 'Writing the molecule gas properties with spin down' )") 
        DO iband1 = nbndmin, nbndmax
          WRITE(stdout,"(4f15.6)") seigen_mol(is, iband1), ssigma_tot(is, iband1)/max_sigma_tot, sbeta(is, iband1)
        ENDDO
      ENDIF
      !
    ENDIF
    !
  ENDDO
  !
  ! write results on data files
  !
  IF (ionode) THEN
     WRITE(stdout,"(/,5x, 'Writing output on file...' )")

     OPEN (30, FILE='spin_resolved_photospectra.dat', FORM='FORMATTED' )
     !
     WRITE(30, "(2x,'# energy grid [eV]   photospectra spin_up [ab.]   photospectra spin_dw [ab.] ')" )
     !
     DO iw =1, nw
         !
         WRITE(30,"(4f15.10)") wgrid(iw), srphotospec(0, iw), srphotospec(1,iw) 
         !
     ENDDO
     !
     CLOSE(30)
     !
  ENDIF
  !
  ! local cleaning
  !
  DEALLOCATE ( g2kin, srphotospec, dipole_2, dipole_aux)
  ! 
ENDIF
  !
  IF (othor_pw) THEN
    DEALLOCATE( dipole_opw_gamma )
    DEALLOCATE( scala_opw_lambda )
    DEALLOCATE( gamma_2_opw )
    DEALLOCATE( lambda_2_opw )
  ENDIF
  !
  CALL grid_destroy()
  !
  return
END SUBROUTINE photoemission_spectr_pw

!--------------------------------------------------------------------
SUBROUTINE dipole_calc_pw ( ik, dipole_aux, dipole_opw_gamma, scala_opw_lambda, metalcalc, &
                            nbndmin, nbndmax, photon_angle, azimuthal_angle, homo_gas, othor_pw) 
  !------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : npw, nbnd, igk, g2kin, et
  USE io_global,            ONLY : ionode, stdout
  USE wavefunctions_module, ONLY : evc
  USE klist,                ONLY : xk
  USE cell_base,            ONLY : tpiba2
  USE gvect,                ONLY : ngm, g, ecutwfc
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE grid_module,          ONLY : focc
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
  USE io_global,            ONLY : ionode, stdout

IMPLICIT NONE
  !
  ! global variables
  INTEGER, INTENT(IN)        :: ik,nbndmin, nbndmax
  REAL(DP),INTENT(IN)        :: photon_angle, azimuthal_angle 
  COMPLEX(DP), INTENT(INOUT) :: dipole_aux(3, nbndmax, npw)
  COMPLEX(DP), INTENT(INOUT) :: dipole_opw_gamma(3, nbndmax, npw)
  COMPLEX(DP), INTENT(INOUT) :: scala_opw_lambda(nbndmax, npw)
  LOGICAL, INTENT(IN)        :: metalcalc, homo_gas, othor_pw
  !
  ! local variables
  INTEGER :: iband1, iband2, ig
  REAL(DP):: sqrtk2
  COMPLEX(DP) :: caux1, caux2
  COMPLEX(DP) :: sumx, sumy, sumz  
  COMPLEX(DP) :: ax(npw), ay(npw), az(npw) 
  COMPLEX(DP) :: ax_add(npw), ay_add(npw), az_add(npw) 
  !
  ! Routine Body
  ! 
  CALL start_clock( 'dipole_calc' )
  !
  ! setup k+G grids for each kpt
  !
  CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)  
  !
  ! read wfc for the given kpt
  !
  CALL davcio (evc, nwordwfc, iunwfc, ik, - 1) 
  !
  ! compute matrix elements
  !
  dipole_aux(:,:,:) = (0.0_DP,0.0_DP)
  !
  IF (othor_pw) THEN
     dipole_opw_gamma(:,:,:) = (0.0_DP,0.0_DP)
     scala_opw_lambda(:,:) = (0.0_DP,0.0_DP)
  ENDIF
  !
  DO iband1 = nbndmin, nbndmax 
     !
     DO ig = 1, npw
        !
        caux1 = evc(ig,iband1)
        !
        IF (homo_gas) then 
           dipole_aux(1,iband1,ig) = dipole_aux(1,iband1,ig) + &
                                     ( g(1,igk(ig)) ) * caux1
           dipole_aux(2,iband1,ig) = dipole_aux(2,iband1,ig) + &
                                     ( g(2,igk(ig)) ) * caux1
           dipole_aux(3,iband1,ig) = dipole_aux(3,iband1,ig) + &
                                     ( g(3,igk(ig)) ) * caux1
        ELSE
           dipole_aux(1,iband1,ig) = dipole_aux(1,iband1,ig) + &
                                     ( g(1,igk(ig)) )* cos(photon_angle) &
                                      * cos(azimuthal_angle) * caux1
           dipole_aux(2,iband1,ig) = dipole_aux(2,iband1,ig) + &
                                     ( g(2,igk(ig)) )* cos(photon_angle) &
                                      * sin(azimuthal_angle) * caux1
           dipole_aux(3,iband1,ig) = dipole_aux(3,iband1,ig) + &
                                     ( g(3,igk(ig)) )* sin(photon_angle) * caux1
        ENDIF
        ! 
     ENDDO
     !
     IF (othor_pw) THEN
        !
        ax(:) = (0.0_DP,0.0_DP) 
        ay(:) = (0.0_DP,0.0_DP) 
        az(:) = (0.0_DP,0.0_DP)
        ax_add(:) = (0.0_DP,0.0_DP)
        ay_add(:) = (0.0_DP,0.0_DP)
        az_add(:) = (0.0_DP,0.0_DP) 
        !
        DO iband2 = nbndmin, nbndmax
          !
          sumx = (0.0_DP,0.0_DP)
          sumy = (0.0_DP,0.0_DP)
          sumz = (0.0_DP,0.0_DP)
          ! 
          DO ig = 1, npw
            caux1 = evc(ig,iband1) 
            caux2 = evc(ig,iband2) 
            if (homo_gas) then
               sumx = sumx +  g(1,igk(ig)) * conjg(caux1)*caux2    
               sumy = sumy +  g(2,igk(ig)) * conjg(caux1)*caux2    
               sumz = sumz +  g(3,igk(ig)) * conjg(caux1)*caux2    
            else
               sumx = sumx +  g(1,igk(ig)) * conjg(caux1)*caux2 * &
                           cos(photon_angle) * cos(azimuthal_angle)
               sumy = sumy +  g(2,igk(ig)) * conjg(caux1)*caux2 * &
                           cos(photon_angle) * sin(azimuthal_angle) 
               sumz = sumz +  g(3,igk(ig)) * conjg(caux1)*caux2 * &
                           sin(photon_angle)
            endif 
          ENDDO
          !
          CALL mp_sum( sumx, intra_pool_comm )
          CALL mp_sum( sumy, intra_pool_comm )
          CALL mp_sum( sumz, intra_pool_comm )
          ! 
          DO ig = 1, npw
            caux2  = evc(ig,iband2)
            ax(ig) = ax(ig) + conjg(caux2)*sumx 
            ay(ig) = ay(ig) + conjg(caux2)*sumy 
            az(ig) = az(ig) + conjg(caux2)*sumz 
          ENDDO
          !
          DO ig = 1, npw
            caux2  = evc(ig,iband2)
            ax_add(ig) = ax_add(ig) + conjg(caux2)*sumx*g(1,igk(ig))
            ay_add(ig) = ay_add(ig) + conjg(caux2)*sumy*g(2,igk(ig))
            az_add(ig) = az_add(ig) + conjg(caux2)*sumz*g(3,igk(ig))
          ENDDO
          !
        ENDDO
        !
        DO ig = 1, npw
          !
          ! gradient component of total sigma
          ! 
          dipole_aux(1,iband1,ig) = dipole_aux(1,iband1,ig) - ax(ig)
          dipole_aux(2,iband1,ig) = dipole_aux(2,iband1,ig) - ay(ig)
          dipole_aux(3,iband1,ig) = dipole_aux(3,iband1,ig) - az(ig)
          !
          ! gradient component of Gamma
          !   
          dipole_opw_gamma(1,iband1,ig) = dipole_opw_gamma(1,iband1,ig) + (ax(ig))
          dipole_opw_gamma(2,iband1,ig) = dipole_opw_gamma(2,iband1,ig) + (ay(ig))
          dipole_opw_gamma(3,iband1,ig) = dipole_opw_gamma(3,iband1,ig) + (az(ig))
          !   
          sqrtk2 = sqrt((xk(1,ik)+g(1,igk(ig)))**2 + (xk(2,ik)+g(2,igk(ig))) **2 + (xk(3,ik)+g(3,igk(ig)))**2)
          !
          IF (sqrtk2 > 1.0E-05) THEN
            !
            ! scala component of Lambda
            !
            scala_opw_lambda(iband1,ig) = scala_opw_lambda(iband1,ig) + & 
                                          (ax_add(ig) + ay_add(ig) + az_add(ig))/sqrtk2      
            !    
          ENDIF
          !
        ENDDO
        !
     ENDIF
     ! 
  ENDDO 
  !
  CALL stop_clock( 'dipole_calc' )
  !
  return
  !
END SUBROUTINE dipole_calc_pw
