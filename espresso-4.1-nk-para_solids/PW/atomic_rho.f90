!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine atomic_rho (rhoa, nspina)
  !-----------------------------------------------------------------------
  ! This routine calculates rhoa as the superposition of atomic charges.
  !
  ! nspina is the number of spin components to be calculated
  !
  ! if nspina = 1 the total atomic charge density is calculated
  ! if nspina = 2 the spin up and spin down atomic charge densities are
  !               calculated assuming an uniform atomic spin-polarization
  !               equal to starting_magnetization(nt)
  ! if nspina = 4 noncollinear case. The total density is calculated
  !               in the first component and the magnetization vector 
  !               in the other three.
  !
  ! NB: nspina may not be equal to nspin because in some cases (as in update)
  ! the total charge only could be needed, even in a LSDA calculation.
  !
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE atom,                 ONLY : rgrid, msh
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : tpiba, omega
  USE gvect,                ONLY : ngm, ngl, nrxx, nr1, nr2, nr3, nrx1, nrx2, &
                                   nrx3, gstart, nl, nlm, gl, igtongl
  USE lsda_mod,             ONLY : starting_magnetization, lsda
  USE vlocal,               ONLY : strf
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE noncollin_module,     ONLY : angle1, angle2
  USE uspp_param,           ONLY : upf
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  implicit none
  !
  integer :: nspina
  ! the number of spin polarizations
  real(DP) :: rhoa (nrxx, nspina)
  ! the output atomic charge
  !
  ! local variables
  !
  real(DP) :: rhoneg, rhoima, gx
  real(DP), allocatable :: rhocgnt (:), aux (:)
  complex(DP), allocatable :: rhocg (:,:)
  integer :: ir, is, ig, igl, nt, ndm
  !
  ! superposition of atomic charges contained in the array rho_at
  ! (read from pseudopotential files)
  !
  ! allocate work space (psic must already be allocated)
  !
  allocate (rhocg(  ngm, nspina))    
  ndm = MAXVAL ( msh(1:ntyp) )
  allocate (aux(ndm))    
  allocate (rhocgnt( ngl))    
  rhoa(:,:) = 0.d0
  rhocg(:,:) = (0.d0,0.d0)

  do nt = 1, ntyp
     !
     ! Here we compute the G=0 term
     !
     if (gstart == 2) then
        do ir = 1, msh (nt)
           aux (ir) = upf(nt)%rho_at (ir)
        enddo
        call simpson (msh (nt), aux, rgrid(nt)%rab, rhocgnt (1) )
     endif
     !
     ! Here we compute the G<>0 term
     !
     do igl = gstart, ngl
        gx = sqrt (gl (igl) ) * tpiba
        do ir = 1, msh (nt)
           if (rgrid(nt)%r(ir) < 1.0d-8) then
              aux(ir) = upf(nt)%rho_at(ir)
           else
              aux(ir) = upf(nt)%rho_at(ir) * &
                        sin(gx*rgrid(nt)%r(ir)) / (rgrid(nt)%r(ir)*gx)
           endif
        enddo
        call simpson (msh (nt), aux, rgrid(nt)%rab, rhocgnt (igl) )
     enddo
     !
     ! we compute the 3D atomic charge in reciprocal space
     !
     if (nspina == 1) then
        do ig = 1, ngm
           rhocg(ig,1) = rhocg(ig,1) + &
                         strf(ig,nt) * rhocgnt(igtongl(ig)) / omega
        enddo
     else if (nspina == 2) then
        do ig = 1, ngm
           rhocg(ig,1) = rhocg(ig,1) + &
                         0.5d0 * ( 1.d0 + starting_magnetization(nt) ) * &
                         strf(ig,nt) * rhocgnt(igtongl(ig)) / omega
           rhocg(ig,2) = rhocg(ig,2) + &
                         0.5d0 * ( 1.d0 - starting_magnetization(nt) ) * &
                         strf(ig,nt) * rhocgnt(igtongl(ig)) / omega
        enddo
     else
!
!    Noncolinear case
!
        do ig = 1,ngm
           rhocg(ig,1) = rhocg(ig,1) + &
                strf(ig,nt)*rhocgnt(igtongl(ig))/omega

           ! Now, the rotated value for the magnetization

           rhocg(ig,2) = rhocg(ig,2) + &
                starting_magnetization(nt)* &
                sin(angle1(nt))*cos(angle2(nt))* &
                strf(ig,nt)*rhocgnt(igtongl(ig))/omega
           rhocg(ig,3) = rhocg(ig,3) + &
                starting_magnetization(nt)* &
                sin(angle1(nt))*sin(angle2(nt))* &
                strf(ig,nt)*rhocgnt(igtongl(ig))/omega
           rhocg(ig,4) = rhocg(ig,4) + &
                starting_magnetization(nt)* &
                cos(angle1(nt))* &
                strf(ig,nt)*rhocgnt(igtongl(ig))/omega
        end do
     endif
  enddo

  deallocate (rhocgnt)
  deallocate (aux)

  do is = 1, nspina
     !
     ! and we return to real space
     !
     psic(:) = (0.d0,0.d0)
     psic (nl (:) ) = rhocg (:, is)
     if (gamma_only) psic ( nlm(:) ) = CONJG( rhocg (:, is) )
     call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
     !
     ! we check that everything is correct
     !
     rhoneg = 0.d0
     rhoima = 0.d0
     do ir = 1, nrxx
        rhoneg = rhoneg + MIN (0.d0,  DBLE (psic (ir)) )
        rhoima = rhoima + abs (AIMAG (psic (ir) ) )
     enddo
     rhoneg = omega * rhoneg / (nr1 * nr2 * nr3)
     rhoima = omega * rhoima / (nr1 * nr2 * nr3)
#ifdef __PARA
     call mp_sum(  rhoneg, intra_pool_comm )
     call mp_sum(  rhoima, intra_pool_comm )
#endif
     IF ( rhoima > 1.0d-4 ) THEN
        WRITE( stdout,'(5x,"Check: imaginary charge or magnetization=",&
          & f12.6," (component ",i1,") set to zero")') rhoima, is
     END IF
     IF ( (is == 1) .OR. lsda ) THEN
        !
        IF ( (rhoneg < -1.0d-4) ) THEN
           IF ( lsda ) THEN 
              WRITE( stdout,'(5x,"Check: negative starting charge=", &
                   &"(component",i1,"):",f12.6)') is, rhoneg
           ELSE
              WRITE( stdout,'(5x,"Check: negative starting charge=", &
          &          f12.6)') rhoneg
           END IF
        END IF
     END IF
     !
     ! set imaginary terms to zero - negative terms are not set to zero
     ! because it is basically useless to do it in real space: negative
     ! charge will re-appear when Fourier-transformed back and forth
     !
     DO ir = 1, nrxx
        rhoa (ir, is) =  DBLE (psic (ir))
     END DO
     !
  enddo

  deallocate (rhocg)
  return
end subroutine atomic_rho

