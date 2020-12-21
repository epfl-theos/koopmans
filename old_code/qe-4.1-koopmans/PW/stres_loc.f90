!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
subroutine stres_loc (sigmaloc)
  !----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE atom,                 ONLY : msh, rgrid
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : omega, tpiba2
  USE gvect,                ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, &
                                   nrx3, nrxx, nl, g, ngl, gl, igtongl
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho
  USE vlocal,               ONLY : strf, vloc
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE uspp_param,           ONLY : upf
  USE noncollin_module,     ONLY : nspin_lsda
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  implicit none
  !
  real(DP) :: sigmaloc (3, 3)
  real(DP) , allocatable :: dvloc(:)
  real(DP) :: evloc, fact
  integer :: ng, nt, l, m, is
  ! counter on g vectors
  ! counter on atomic type
  ! counter on angular momentum
  ! counter on spin components

  allocate(dvloc(ngl))
  sigmaloc(:,:) = 0.d0
  psic(:)=(0.d0,0.d0)
  do is = 1, nspin_lsda
     call DAXPY (nrxx, 1.d0, rho%of_r (1, is), 1, psic, 2)
  enddo

  call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  ! psic contains now the charge density in G space
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  evloc = 0.0d0
  do nt = 1, ntyp
     if (gstart==2) evloc = evloc + &
          psic (nl (1) ) * strf (1, nt) * vloc (igtongl (1), nt)
     do ng = gstart, ngm
        evloc = evloc +  DBLE (CONJG(psic (nl (ng) ) ) * strf (ng, nt) ) &
             * vloc (igtongl (ng), nt) * fact
     enddo
  enddo
  !
  !      WRITE( 6,*) ' evloc ', evloc, evloc*omega   ! DEBUG
  !
  do nt = 1, ntyp
     IF ( .NOT. ASSOCIATED ( upf(nt)%vloc ) ) THEN
        !
        ! special case: pseudopotential is coulomb 1/r potential
        !
        call dvloc_coul (upf(nt)%zp, tpiba2, ngl, gl, omega, dvloc)
        !
     ELSE
        !
        ! normal case: dvloc contains dV_loc(G)/dG
        !
        call dvloc_of_g (rgrid(nt)%mesh, msh (nt), rgrid(nt)%rab, rgrid(nt)%r,&
          upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngl, gl, omega, dvloc)
        !
     END IF
     ! no G=0 contribution
     do ng = 1, ngm
        do l = 1, 3
           do m = 1, l
              sigmaloc(l, m) = sigmaloc(l, m) +  DBLE( CONJG( psic(nl(ng) ) ) &
                    * strf (ng, nt) ) * 2.0d0 * dvloc (igtongl (ng) ) &
                    * tpiba2 * g (l, ng) * g (m, ng) * fact
           enddo
        enddo
     enddo
  enddo
  !
  do l = 1, 3
     sigmaloc (l, l) = sigmaloc (l, l) + evloc
     do m = 1, l - 1
        sigmaloc (m, l) = sigmaloc (l, m)
     enddo
  enddo
#ifdef __PARA
  call mp_sum(  sigmaloc, intra_pool_comm )
#endif
  deallocate(dvloc)
  return
end subroutine stres_loc

