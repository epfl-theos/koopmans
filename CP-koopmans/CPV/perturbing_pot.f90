!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!#define DEBUG
!
!-----------------------------------------------------------------------
      SUBROUTINE perturbing_pot(current_bnd, current_spin, rhor, pot)
!-----------------------------------------------------------------------
! ... Calculate pot(r) = \int dr' f_Hxr(r,r')n_i(r') ...
! take as an input:
! 1) the band and spin index of the chosen orbital
! 2) the total charge density (to compute the xc kernel)
! give as OUTPUT:
! the perturbing potential in pot
! 
!
      USE kinds,                 ONLY : DP
      USE constants,             ONLY : fpi, hartree_si, electronvolt_si
      USE cell_base,             ONLY : tpiba2, omega
      USE gvecw,                 ONLY : ngw
      USE gvecp,                 ONLY : ngm
      USE recvecs_indexes,       ONLY : np, nm
      USE grid_dimensions,       ONLY : nnrx, nr1, nr2, nr3
      USE electrons_base,        ONLY : nspin
      USE control_flags,         ONLY : gamma_only, do_wf_cmplx
      USE reciprocal_vectors,    ONLY : gstart, g
      USE cp_interfaces,         ONLY : fwfft, invfft
      USE fft_base,              ONLY : dfftp
      USE mp,                    ONLY : mp_sum
      USE mp_global,             ONLY : intra_image_comm
      USE cp_interfaces,         ONLY : nksic_get_orbitalrho
      USE io_global,             ONLY : stdout
      USE eecp_mod,              ONLY : do_comp
      USE wavefunctions_module,  ONLY : c0
      USE funct,                 ONLY : dmxc, dmxc_spin
      USE uspp,                  ONLY : okvan
      USE twin_types
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN):: current_spin, current_bnd
      ! band and spin index of the orbital 
      REAL(DP), INTENT(IN):: rhor(nnrx,nspin)
      ! INPUT: total charge density to compute the xc kernel
      REAL(DP) :: orb_rhor(nnrx), dmuxc_r(nnrx, nspin, nspin)
      ! orbital density as an input 
      REAL(DP), INTENT(OUT):: pot(nnrx)
      ! OUTPUT: potential \int f_Hxc(r,r')*orb_rho_r(r')
      COMPLEX(DP), ALLOCATABLE :: orb_rhog(:,:), vhaux(:), vtmp(:), vcorr(:)
      ! ... orbital density in G space, auxiliary potential in real and G space
      COMPLEX(DP) :: c1(ngw), psi1(nnrx)
      ! ... auxiliary variable for the wavefunction ...
      INTEGER ig, ir, is, is1
      ! ... counter  
      REAL(DP) sum, uPi
      !
      !
      !
      ALLOCATE ( orb_rhog(ngm,1) )
      ALLOCATE ( vhaux(nnrx) )
      ALLOCATE ( vtmp(ngm) )
      ALLOCATE ( vcorr(ngm) )
      !
      IF (okvan) call errore('perturbing_pot','USPP not implemented yet',1)
      !
      ! ... xc kernel initialization ...
      !
      dmuxc_r = 0.D0
      DO ir =1 , nnrx
         !
         IF (nspin == 2) THEN
            CALL dmxc_spin (rhor(ir,1), rhor(ir,2), dmuxc_r(ir,1,1), dmuxc_r(ir,1,2), dmuxc_r(ir,2,1), dmuxc_r(ir,2,2) )
            dmuxc_r(ir,:,:)=dmuxc_r(ir,:,:)*0.5D0 ! From Rydberg to Hartree
#ifdef DEBUG
            IF (MOD(ir,200) == 0) WRITE(stdout,'(I8, 6F15.6)') ir, rhor(ir,1), rhor(ir,2), dmuxc_r(ir,1,1), dmuxc_r(ir,1,2), dmuxc_r(ir,2,1), dmuxc_r(ir,2,2)
#endif
         ELSE
            dmuxc_r(ir,1,1) = 0.5D0*dmxc(rhor(ir,1)) ! factor 0.5 from Ry to Ha
         ENDIF
         !
      ENDDO
      !
      pot=0.D0
      orb_rhog = (0.D0,0.D0)
      !
      ! ... Orbital density ...
      ! 
      c1=c0(:, current_bnd)
      CALL c2psi ( psi1, nnrx, c1, c1, ngw, 0 )
      CALL invfft('Dense', psi1, dfftp ) 
      !
      orb_rhor= 0.D0
      DO ir = 1, nnrx
         !
         orb_rhor(ir) = (( abs(psi1(ir)) ))**2/omega
         !
      ENDDO
      !
#ifdef DEBUG
      sum=0.D0
      DO ir = 1, nnrx
         sum =sum + orb_rhor(ir)
      ENDDO
      CALL mp_sum (sum, intra_image_comm)
      WRITE(stdout,'(2x, "orbital charge", 2F18.12)') sum/( nr1*nr2*nr3 )*omega
#endif 
      ! ... Hartree potential ...
      ! 
      ! orbital density in G spce 
      vhaux(:) = (0.D0, 0.D0)
      vhaux=CMPLX(orb_rhor(:),0.D0)
      CALL fwfft('Dense',vhaux,dfftp )
      DO ig = 1,ngm; orb_rhog(ig,1) = vhaux( np(ig) ); ENDDO
      ! 
      ! compute hartree like potential 
      IF ( gstart == 2 ) vtmp(1)=(0.d0,0.d0)
      DO ig=gstart,ngm
          vtmp(ig) = orb_rhog(ig,1) * fpi/( tpiba2*g(ig) )
      ENDDO
      ! 
      ! compute periodic corrections
      IF ( do_comp ) THEN
          !
          call calc_compensation_potential( vcorr, orb_rhog(:,1),.true.)
          vtmp(:) = vtmp(:) + vcorr(:)
          !
      ENDIF
      !
      ! Go back to real space
      vhaux = (0.D0,0.D0)
      DO ig = 1, ngm
         vhaux(np(ig)) = vtmp(ig)
         vhaux(nm(ig)) = CONJG( vtmp(ig) )
      ENDDO
      !
      CALL invfft('Dense', vhaux, dfftp)
      !
      ! update the potential
      pot(1:nnrx) = DBLE(vhaux(1:nnrx))
      !
#ifdef DEBUG
      sum=0.D0 
      DO ir = 1, nnrx
         sum = sum + orb_rhor(ir)*pot(ir)
      ENDDO
      CALL mp_sum (sum, intra_image_comm)
      sum = sum *0.5D0*omega/(nr1*nr2*nr3)
      WRITE(stdout,*) "hartree energy", sum, "Ha" 
      !
      DO is = 1 ,nspin
         DO is1 =1, nspin
            WRITE(6,'(i3,i3,3F18.6,/)') is,is1,dmuxc_r(1:3,is,is1)
         ENDDO
       ENDDO
#endif
      !
      ! ... Add xc contribution
      DO ir = 1, nnrx
          pot(ir) = pot(ir) + dmuxc_r(ir,current_spin,current_spin) * orb_rhor(ir)
      ENDDO
      !
#ifdef DEBUG
      uPi=0.D0
      DO ir = 1, nnrx
         uPi = uPi + ( orb_rhor(ir) *  pot(ir) )
      ENDDO
      CALL mp_sum(  uPi , intra_image_comm )
      uPi = uPi / ( nr1*nr2*nr3 )*omega
      WRITE(stdout,*) "unrelaxed Koopmans uPi=", uPi
#endif
      !
      DEALLOCATE ( orb_rhog )
      DEALLOCATE ( vhaux )
      DEALLOCATE ( vtmp )
      !
!---------------------------------------------------------------
      END SUBROUTINE perturbing_pot
!---------------------------------------------------------------
