!
! Copyright (C) 2003-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!          25/06/2009 (Riccardo Sabatini)
!               reformulation using a unique saw(x) function (included in 
!               cell_base) in all e-field related routines and inclusion of 
!               a macroscopic electronic dipole contribution in the mixing 
!               scheme. 
!
!   the calculation of the dipole is split in the ionic (compute_ion_dip)
!   and electronic (compute_el_dip) contributions.
!
SUBROUTINE compute_ion_dip(emaxpos, eopreg, edir, ion_dipole)
  !
  !
  !---------------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout, ionode
  USE ions_base,  ONLY : nat, ityp, tau, zv
  USE constants, ONLY : fpi
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg, omega, alat, saw
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: emaxpos, eopreg
  INTEGER, INTENT(IN)  :: edir
  REAL(DP), INTENT(OUT) ::  ion_dipole
  !
  REAL(DP) :: bmod
  INTEGER  :: na
  REAL(DP) :: sawarg, tvectb, zvia

  !--------------------------
  !  Fix some values for later calculations
  !--------------------------
  bmod=SQRT(bg(1,edir)**2+bg(2,edir)**2+bg(3,edir)**2)

  !--------------------------
  !  Calculate IONIC dipole
  !--------------------------
  
  !
  ! P_{ion} = \sum^{nat}_{s} z_{v} Saw\left( \vec{t_{s}}\cdot\vec{b_{edir}}} 
  !                            \right) \frac{alat}{bmod} \frac{4\pi}{\Omega}
  !

  ion_dipole=0.d0
  
  DO na = 1, nat
     !
     ! Ion charge
     zvia = zv(ityp(na))
     
     ! Position vector
     tvectb = tau(1,na)*bg(1,edir) + tau(2,na)*bg(2,edir) + tau(3,na)*bg(3,edir)

     ion_dipole = ion_dipole + zvia* saw(emaxpos,eopreg, tvectb ) &
                                                * (alat/bmod) * (fpi/omega)
 
  END DO

  
  RETURN
  
END SUBROUTINE compute_ion_dip
!
SUBROUTINE compute_el_dip(emaxpos, eopreg, edir, charge, e_dipole)
  !
  !
  !---------------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout, ionode
  USE lsda_mod,     ONLY : nspin
  USE constants, ONLY : fpi
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg, omega, alat, saw
  USE gvect,      ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
  USE fft_base,   ONLY : dfftp
  USE mp_global,  ONLY : me_pool, intra_pool_comm
  USE mp,         ONLY : mp_sum
  USE fft_base,  ONLY : grid_gather  
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: emaxpos, eopreg
  REAL(DP), INTENT(IN), DIMENSION(nrxx,nspin) :: charge
  INTEGER, INTENT(IN)  :: edir
  REAL(DP), INTENT(OUT) ::  e_dipole
  !
  REAL(DP), ALLOCATABLE :: rho_all(:), aux(:)
  REAL(DP) :: rhoir,bmod
  INTEGER  :: i, k, j, ip, ir, index, index0, na
  REAL(DP) :: sawarg, tvectb

  !--------------------------
  !  Fix some values for later calculations
  !--------------------------
  bmod=SQRT(bg(1,edir)**2+bg(2,edir)**2+bg(3,edir)**2)

  !
  !--------------------------
  !  Calculate ELECTRONIC dipole
  !--------------------------  
  
  !
  ! Case with edir = 3 (in the formula changes only tha rgument of saw, i for
  ! edir=1 and j for edir = 2)
  !
  ! P_{ele} = \sum_{ijk} \rho_{r_{ijk}} Saw\left( \frac{k}{nr3} \right) 
  !                   \frac{alat}{bmod} \frac{\Omega}{nrxx} \frac{4\pi}{\Omega}
  !

  e_dipole  = 0.D0
  
  !
  ! Procedure for parallel summation
  !

  index0 = 0
  !
#if defined (__PARA)
  !
  DO i = 1, me_pool
     index0 = index0 + nrx1*nrx2*dfftp%npp(i)
  END DO
  !
#endif

  !
  ! Loop in the charge array
  !
  DO ir = 1, nrxx
     !
     ! ... three dimensional indexes
     !
     index = index0 + ir - 1
     k     = index / (nrx1*nrx2)
     index = index - (nrx1*nrx2)*k
     j     = index / nrx1
     index = index - nrx1*j
     i     = index
     
     !
     ! Define the argument for the saw function     
     !
     if (edir.eq.1) sawarg = (i*1.0)/(nrx1*1.0)
     if (edir.eq.2) sawarg = (j*1.0)/(nrx2*1.0)
     if (edir.eq.3) sawarg = (k*1.0)/(nrx3*1.0)
     
     rhoir = charge(ir,1)
     !
     IF ( nspin == 2 ) rhoir = rhoir + charge(ir,2)
          
     e_dipole = e_dipole + rhoir * saw(emaxpos,eopreg, sawarg) &
                                   * (alat/bmod) * (fpi/(nr1*nr2*nr3))

  END DO

  CALL mp_sum(  e_dipole , intra_pool_comm )
  
  RETURN
  
END SUBROUTINE compute_el_dip
