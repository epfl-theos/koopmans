!
! Copyright (C) 2003-2004 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... written by J. Tobik
!
! Changes 30/06/2003 (ADC) : 
!               Calculation of corrections to energy and forces due
!               to the field.
!               Added possibility to subtract the dipole field 
!               for slab or molecule calculation.
!               (See Bengtsson PRB 59, 12 301 (1999) and
!                    Meyer and Vanderbilt, PRB 63, 205426 (2001).)
!
!          25/06/2009 (Riccardo Sabatini)
!               reformulation using a unique saw(x) function (included in 
!               cell_base) in all e-field related routines and inclusion of 
!               a macroscopic electronic dipole contribution in the mixing 
!               scheme. 
!
#include "f_defs.h"

!
!--------------------------------------------------------------------------
SUBROUTINE add_efield(vpoten,etotefield,rho,iflag)
  !--------------------------------------------------------------------------
  !
  !   This routine adds an electric field to the local potential. The
  !   field is made artificially periodic by introducing a saw-tooth
  !   potential. The field is parallel to a reciprocal lattice vector bg, 
  !   according to the index edir.
  !
  !   if dipfield is false the electric field correction is added to the
  !   potential given as input (the bare local potential) only
  !   at the first call to this routine. In the following calls
  !   the routine exit.
  !
  !   if dipfield is true the dipole moment per unit surface is calculated
  !   and used to cancel the electric field due to periodic boundary
  !   conditions. This potential is added to the Hartree and xc potential
  !   in v_of_rho. NB: in this case the electric field contribution to the 
  !   band energy is subtracted by deband.
  !
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8, e2
  USE ions_base,     ONLY : nat, ityp, zv
  USE cell_base,     ONLY : alat, at, omega, bg, saw
  USE extfield,      ONLY : tefield, dipfield, edir, eamp, emaxpos, &
                            eopreg, forcefield
  USE force_mod,     ONLY : lforce
  USE gvect,         ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
  USE io_global,     ONLY : stdout,ionode
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_global,     ONLY : intra_image_comm, me_pool, intra_pool_comm
  USE fft_base,      ONLY : dfftp
  USE mp,            ONLY : mp_bcast, mp_sum
  USE control_flags, ONLY : iverbosity
  
  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP),INTENT(INOUT) :: vpoten(nrxx)    ! the ef is added to this potential
  REAL(DP),INTENT(OUT)   :: etotefield      ! the contribution to etot due to ef
  REAL(DP),INTENT(IN)    :: rho(nrxx,nspin) ! the density whose dipole is computed
  LOGICAL,INTENT(IN)     :: iflag ! set to true to force recalculation of field
  !
  ! local variables
  !
  INTEGER :: index, index0, i, j, k
  INTEGER :: ir, na, ipol
  REAL(DP) :: length, vamp, value, sawarg, e_dipole, ion_dipole
  REAL(DP) :: tot_dipole, bmod, debye

  LOGICAL :: first=.TRUE.
  SAVE first
  
  !---------------------
  !  Execution control
  !---------------------

  IF (.NOT.tefield) RETURN
  ! efiled only needs to be added on the first iteration, if dipfield
  ! is not used. note that for relax calculations it has to be added
  ! again on subsequent relax steps.
  IF ((.NOT.dipfield).AND.(.NOT.first) .AND..NOT. iflag) RETURN
  first=.FALSE.

  IF ((edir.lt.1).or.(edir.gt.3)) THEN
     CALL errore('add_efield',' wrong edir',1)
  ENDIF

  !---------------------
  !  Variable initialization
  !---------------------

  bmod=SQRT(bg(1,edir)**2+bg(2,edir)**2+bg(3,edir)**2)
  debye = 2.54176D0

  tot_dipole=0.d0
  e_dipole=0.d0
  ion_dipole=0.d0
  
  !---------------------
  !  Calculate dipole
  !---------------------
  
  if (dipfield) then
  !
  ! dipole correction is active 
  !
     CALL compute_el_dip(emaxpos, eopreg, edir, rho, e_dipole)
     CALL compute_ion_dip(emaxpos, eopreg, edir, ion_dipole)
    
     tot_dipole  = -e_dipole + ion_dipole

#ifdef __PARA
     CALL mp_bcast(tot_dipole, 0, intra_image_comm)
#endif
  !  
  !  E_{TOT} = -e^{2} \left( eamp - dip \right) dip \frac{\Omega}{4\pi} 
  !
     etotefield=-e2*(eamp-tot_dipole/2.d0)*tot_dipole*omega/fpi 

  !---------------------
  !  Define forcefield
  !  
  !  F_{s} = e^{2} \left( eamp - dip \right) z_{v}\cross\frac{\vec{b_{3}}}{bmod} 
  !---------------------
    
     IF (lforce) THEN
        DO na=1,nat
           DO ipol=1,3
              forcefield(ipol,na)= e2 *(eamp - tot_dipole) &
                               *zv(ityp(na))*bg(ipol,edir)/bmod
           ENDDO
        ENDDO
     ENDIF

  else
  !
  ! dipole correction is not active
  !

     CALL compute_ion_dip(emaxpos, eopreg, edir, ion_dipole)

  !  
  !  E_{TOT} = -e^{2} eamp * iondip \frac{\Omega}{4\pi} 
  !
     etotefield=-e2*eamp*ion_dipole*omega/fpi 

  !---------------------
  !  Define forcefield
  !  
  !  F_{s} = e^{2}  eamp z_{v}\cross\frac{\vec{b_{3}}}{bmod} 
  !---------------------
    
     IF (lforce) THEN
        DO na=1,nat
           DO ipol=1,3
              forcefield(ipol,na)= e2 *eamp &
                               *zv(ityp(na))*bg(ipol,edir)/bmod
           ENDDO
        ENDDO
     ENDIF

  end if

  !
  !  Calcualte potential and print values 
  !   
  
  length=(1.0-eopreg)*(alat*SQRT(at(1,edir)**2+at(2,edir)**2+at(3,edir)**2))
  
  vamp=e2*(eamp-tot_dipole)*length

  IF (ionode) THEN
       !
       ! Output data
       !
       WRITE( stdout,*)
       WRITE( stdout,'(5x,"Adding external electric field":)')

       IF (dipfield) then
          WRITE( stdout,'(/5x,"Computed dipole along edir(",i1,") : ")' ) edir

          !
          !  If verbose prints also the different components
          !
          IF (iverbosity>0) THEN
              WRITE( stdout, '(8X,"Elec. dipole ",1F15.4," au,  ", 1F15.4," Debye")' ) &
                                            e_dipole, (e_dipole*debye)
              WRITE( stdout, '(8X,"Ion. dipole  ",1F15.4," au,  ", 1F15.4," Debye")' ) &
                                          ion_dipole, (ion_dipole*debye)
          ENDIF

          WRITE( stdout, '(8X,"Dipole       ",1F15.4," au,  ", 1F15.4," Debye")' ) &
                                            (tot_dipole* (omega/fpi)),   &
                                            ((tot_dipole* (omega/fpi))*debye)  

          WRITE( stdout, '(8x,"Dipole field     ", f11.4," au")') tot_dipole
          WRITE( stdout,*)

       ENDIF

       IF (abs(eamp)>0.d0) WRITE( stdout,'(8x,"E field amplitude [a.u.]: ", es11.4)') eamp 
        
       WRITE( stdout,'(8x,"Potential amp.   ", f11.4," Ry")') vamp 
       WRITE( stdout,'(8x,"Total length     ", f11.4," bhor")') length
       WRITE( stdout,*)     
  ENDIF


  !
  !------------------------------
  !  Add potential
  !  
  !  V\left(ijk\right) = e^{2} \left( eamp - dip \right) z_{v} 
  !          Saw\left( \frac{k}{nr3} \right) \frac{alat}{bmod} 
  !          
  !---------------------

  ! Index for parallel summation
  !
  index0 = 0
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
     
     if (edir.eq.1) sawarg = (i*1.0)/(nrx1*1.0)
     if (edir.eq.2) sawarg = (j*1.0)/(nrx2*1.0)
     if (edir.eq.3) sawarg = (k*1.0)/(nrx3*1.0)
     
     value = e2*(eamp - tot_dipole)*saw(emaxpos,eopreg,sawarg) * (alat/bmod)

     vpoten(ir) = vpoten(ir) + value

  END DO
  
  
  RETURN

END SUBROUTINE add_efield
