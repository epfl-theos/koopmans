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
SUBROUTINE atomic_wfc_nc_proj (ik, wfcatom)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the superposition of atomic wavefunctions
  ! for k-point "ik" - output in "wfcatom" - noncolinear case only
  ! If lspinorb=.TRUE. it makes linear combinations of eigenstates of 
  ! the atomic total angular momenta j and j_z; otherwise, of eigenstates of 
  ! the orbital angular momenta l, l_z and of s_z (the z-component of the spin).
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi, fpi, pi
  USE cell_base,  ONLY : omega, tpiba
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,      ONLY : natomwfc
  USE gvect,      ONLY : ig1, ig2, ig3, eigts1, eigts2, eigts3, g
  USE klist,      ONLY : xk
  USE wvfct,      ONLY : npwx, npw, nbnd, igk
  USE us,         ONLY : tab_at, dq
  USE uspp_param, ONLY : upf
  USE noncollin_module, ONLY : noncolin, npol, angle1, angle2
  USE spin_orb,   ONLY : lspinorb, rot_ylm, fcoef, lmaxx
  !
  implicit none
  !
  integer, intent(in) :: ik
  complex(DP), intent(out) :: wfcatom (npwx, npol, natomwfc)
  !
  integer :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3, nwfcm
  real(DP), allocatable :: qg(:), ylm (:,:), chiq (:,:,:), gk (:,:)
  complex(DP), allocatable :: sk (:), aux(:)
  complex(DP) :: kphase, lphase
  real(DP) :: arg, px, ux, vx, wx

  call start_clock ('atomic_wfc')

  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = 0
  do nt = 1, ntyp
     lmax_wfc = MAX ( lmax_wfc, MAXVAL (upf(nt)%lchi(1:upf(nt)%nwfc) ) )
  enddo
  !
  nwfcm = MAXVAL ( upf(1:ntyp)%nwfc )
  !
  allocate ( ylm (npw,(lmax_wfc+1)**2), chiq(npw,nwfcm,ntyp), &
             sk(npw), gk(3,npw), qg(npw) )
  !
  do ig = 1, npw
     gk (1,ig) = xk(1, ik) + g(1, igk(ig) )
     gk (2,ig) = xk(2, ik) + g(2, igk(ig) )
     gk (3,ig) = xk(3, ik) + g(3, igk(ig) )
     qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  !  ylm = spherical harmonics
  !
  call ylmr2 ((lmax_wfc+1)**2, npw, gk, qg, ylm)
  !
  ! set now q=|k+G| in atomic units
  !
  do ig = 1, npw
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  do nt = 1, ntyp
     do nb = 1, upf(nt)%nwfc
        if ( upf(nt)%oc (nb) >= 0.d0) then
           do ig = 1, npw
              px = qg (ig) / dq - int (qg (ig) / dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT( qg (ig) / dq ) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq (ig, nb, nt) = &
                     tab_at (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                     tab_at (i1, nb, nt) * px * vx * wx / 2.d0 - &
                     tab_at (i2, nb, nt) * px * ux * wx / 2.d0 + &
                     tab_at (i3, nb, nt) * px * ux * vx / 6.d0
           enddo
        endif
     enddo
  enddo

  deallocate (qg, gk)
  allocate ( aux(npw) )

  do na = 1, nat
     arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
     kphase = CMPLX (cos (arg), - sin (arg) )
     !
     !     sk is the structure factor
     !
     do ig = 1, npw
        iig = igk (ig)
        sk (ig) = kphase * eigts1 (ig1 (iig), na) * eigts2 (ig2 (iig), na) * &
                           eigts3 (ig3 (iig), na)
     enddo
     !
     nt = ityp (na)
     do nb = 1, upf(nt)%nwfc
        if (upf(nt)%oc(nb) >= 0.d0) then
           l = upf(nt)%lchi(nb)
           lphase = (0.d0,1.d0)**l
           !
           !  the factor i^l MUST BE PRESENT in order to produce
           !  wavefunctions for k=0 that are real in real space
           !
           IF ( lspinorb ) THEN
              !
              IF ( upf(nt)%has_so ) THEN
                 !
                 call atomic_wfc_so ( )
                 !
              ELSE
                 !
                 call atomic_wfc_so2 ( )
                 !
              ENDIF
              !
           ELSE
              !
              call atomic_wfc_nc_z ( )
              !
           END IF
           !
        END IF
        !
     END DO
     !
  END DO

  if (n_starting_wfc /= natomwfc) call errore ('atomic_wfc_nc_proj', &
       'internal error: some wfcs were lost ', 1)

  deallocate(aux, sk, chiq, ylm)

  call stop_clock ('atomic_wfc')
  return

CONTAINS

  SUBROUTINE atomic_wfc_so ( )
   !
   ! ... spin-orbit case
   !
   real(DP) :: fact(2), j
   real(DP), external :: spinor
   integer :: ind, ind1, n1, is, sph_ind
   !
   j = upf(nt)%jchi(nb)
   do m = -l-1, l
      fact(1) = spinor(l,j,m,1)
      fact(2) = spinor(l,j,m,2)
      if (abs(fact(1)) > 1.d-8 .or. abs(fact(2)) > 1.d-8) then
         n_starting_wfc = n_starting_wfc + 1
         if (n_starting_wfc > natomwfc) call errore &
              ('atomic_wfc_so', 'internal error: too many wfcs', 1)
         DO is=1,2
            IF (abs(fact(is)) > 1.d-8) THEN
               ind=lmaxx+1+sph_ind(l,j,m,is)
               aux=(0.d0,0.d0)
               DO n1=1,2*l+1
                  ind1=l**2+n1
                  if (abs(rot_ylm(ind,n1)) > 1.d-8) &
                      aux(:)=aux(:)+rot_ylm(ind,n1)*ylm(:,ind1)
               ENDDO
               DO ig=1,npw
                  wfcatom (ig,is,n_starting_wfc) = lphase*fact(is)*&
                        sk(ig)*aux(ig)*chiq (ig, nb, nt)
               END DO
            ELSE
                wfcatom (:,is,n_starting_wfc) = (0.d0,0.d0)
            END IF
         END DO
      END IF
   END DO
   !
   END SUBROUTINE atomic_wfc_so
   ! 
   SUBROUTINE atomic_wfc_so2 ( )
   !
   ! ... spin-orbit case with no spin-orbit PP
   !
   real(DP) :: fact(2), j
   real(DP), external :: spinor
   integer :: ind, ind1, n1, n2, is, sph_ind
   !
   DO n2 = l, l + 1
      j = n2 - 0.5_dp
      IF (j > 0.0_dp)  THEN 
         DO m = -l-1, l
            fact(1) = spinor(l,j,m,1)
            fact(2) = spinor(l,j,m,2)
            IF (abs(fact(1)) > 1.d-8 .or. abs(fact(2)) > 1.d-8) THEN
               n_starting_wfc = n_starting_wfc + 1
               IF (n_starting_wfc > natomwfc) CALL errore &
                  ('atomic_wfc_so2', 'internal error: too many wfcs', 1)
               DO is=1,2
                  IF (abs(fact(is)) > 1.d-8) THEN
                     ind=lmaxx+1+sph_ind(l,j,m,is)
                     aux=(0.0_dp,0.0_dp)
                     DO n1=1,2*l+1
                        ind1=l**2+n1
                        IF (abs(rot_ylm(ind,n1)) > 1.d-8) &
                           aux(:)=aux(:)+rot_ylm(ind,n1)*ylm(:,ind1)
                     ENDDO
                     DO ig=1,npw
                        wfcatom (ig,is,n_starting_wfc) = lphase * &
                           fact(is)*sk(ig)*aux(ig)*chiq(ig,nb,nt)
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO 
      ENDIF
   ENDDO 
   !
   END SUBROUTINE atomic_wfc_so2
   !
   SUBROUTINE atomic_wfc_nc_z ( )
   !
   ! ... noncolinear case, magnetization along z
   !
   DO m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
      if (n_starting_wfc + 2*l + 1 > natomwfc) call errore &
            ('atomic_wfc_nc', 'internal error: too many wfcs', 1)
      DO ig=1,npw
         aux(ig) = sk(ig)*ylm(ig,lm)*chiq(ig,nb,nt)
      END DO
!
      DO ig=1,npw
         wfcatom(ig,1,n_starting_wfc) = aux(ig)
         wfcatom(ig,2,n_starting_wfc) = (0.0_dp, 0.0_dp)
!
         wfcatom(ig,1,n_starting_wfc+2*l+1) = (0.0_dp, 0.0_dp)
         wfcatom(ig,2,n_starting_wfc+2*l+1) = aux(ig)
      END DO
   END DO
   n_starting_wfc = n_starting_wfc + 2*l+1
   !
   END SUBROUTINE atomic_wfc_nc_z

END SUBROUTINE atomic_wfc_nc_proj
