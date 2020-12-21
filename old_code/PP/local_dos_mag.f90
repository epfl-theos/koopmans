!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE local_dos_mag(spin_component, kpoint, kband, raux)
  !----------------------------------------------------------------------------
  !
  ! ... calculates the symmetrized charge density and sum of occupied
  ! ... eigenvalues.
  ! ... this version works also for metals (gaussian spreading technique)  
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : omega,tpiba2
  USE gvect,                ONLY : nrxx, ngm, g, ecutwfc 
  USE gsmooth,              ONLY : nls, nr1s, nr2s, nr3s, &
                                   nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
  USE klist,                ONLY : nks, xk
  USE scf,                  ONLY : rho
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE uspp,                 ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions_module, ONLY : evc, psic_nc
  USE noncollin_module,     ONLY : noncolin, npol
  USE spin_orb,             ONLY : lspinorb, fcoef
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, g2kin
  USE becmod,               ONLY : calbec
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: spin_component, kpoint, kband
  REAL(DP) :: raux(nrxx)

  INTEGER :: ikb, jkb, ijkb0, ih, jh, ijh, na, np
  ! counters on beta functions, atoms, pseudopotentials  
  INTEGER :: ir, is, ig, ibnd, ik
  ! counter on 3D r points
  ! counter on spin polarizations
  ! counter on g vectors
  ! counter on bands
  ! counter on k points  
  !
  REAL(DP) :: w1
  ! weights
  COMPLEX(DP), ALLOCATABLE :: becp_nc(:,:,:)
  ! contains <beta|psi>
  !
  COMPLEX(DP), ALLOCATABLE :: be1(:,:), be2(:,:)
  !
  INTEGER :: ipol, kh, kkb, is1, is2

  becsum(:,:,:) = 0.D0
  rho%of_r(:,:)      = 0.D0
  w1=1.D0/omega

  ALLOCATE( becp_nc( nkb, npol, nbnd ) )
  IF (lspinorb) ALLOCATE(be1(nhm,2), be2(nhm,2))
  !
  ! ... here we sum for each k point the contribution
  ! ... of the wavefunctions to the charge
  !
  DO ik = 1, nks
     IF (ik == kpoint) THEN
        CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)
        IF (nkb > 0) CALL init_us_2 (npw, igk, xk (1, ik), vkb)
        CALL calbec ( npw, vkb, evc, becp_nc)
        !
        !
        DO ibnd = 1, nbnd
           !
           IF (ibnd == kband) then
              psic_nc = (0.D0,0.D0)
              DO ig = 1, npw
                 psic_nc(nls(igk(ig)),1)=evc(ig     ,ibnd)
                 psic_nc(nls(igk(ig)),2)=evc(ig+npwx,ibnd)
              END DO
              DO ipol=1,npol
                 call cft3s (psic_nc(1,ipol), nr1s, nr2s, nr3s, nrx1s, &
                                                          nrx2s, nrx3s, 2)
              END DO
              IF (spin_component==1) THEN
                 DO ir = 1,nrxxs
                    rho%of_r(ir,2) = rho%of_r(ir,2) + 2.D0*w1* &
                         (DBLE(psic_nc(ir,1))* DBLE(psic_nc(ir,2)) + &
                         AIMAG(psic_nc(ir,1))*AIMAG(psic_nc(ir,2)))
                 END DO
              END IF
              IF (spin_component==2) THEN
                 DO ir = 1,nrxxs
                    rho%of_r(ir,3) = rho%of_r(ir,3) + 2.D0*w1* &
                         (DBLE(psic_nc(ir,1))*AIMAG(psic_nc(ir,2)) - &
                          DBLE(psic_nc(ir,2))*AIMAG(psic_nc(ir,1)))
                 END DO
              END IF
              IF (spin_component==3) THEN
                 DO ir = 1,nrxxs
                    rho%of_r(ir,4) = rho%of_r(ir,4) + w1* &
                        (DBLE(psic_nc(ir,1))**2+AIMAG(psic_nc(ir,1))**2 &
                        -DBLE(psic_nc(ir,2))**2-AIMAG(psic_nc(ir,2))**2)
                 END DO
              END IF

              ijkb0 = 0
              DO np = 1, ntyp
                 !
                 IF ( upf(np)%tvanp ) THEN
                    !
                    DO na = 1, nat
                       !
                       IF (ityp(na)==np) THEN
                          !
                          IF (upf(np)%has_so) THEN
                             be1=(0.d0,0.d0)
                             be2=(0.d0,0.d0)
                             DO ih = 1, nh(np)
                                ikb = ijkb0 + ih
                                DO kh = 1, nh(np)
                                   IF ((nhtol(kh,np)==nhtol(ih,np)).AND. &
                                   (nhtoj(kh,np)==nhtoj(ih,np)).AND.     &
                                   (indv(kh,np)==indv(ih,np))) THEN
                                      kkb=ijkb0 + kh
                                      DO is1=1,2
                                         DO is2=1,2
                                            be1(ih,is1)=be1(ih,is1)+  &
                                               fcoef(ih,kh,is1,is2,np)*  &
                                                    becp_nc(kkb,is2,ibnd)
                                            be2(ih,is1)=be2(ih,is1)+ &
                                               fcoef(kh,ih,is2,is1,np)* &
                                               CONJG(becp_nc(kkb,is2,ibnd))
                                         END DO
                                      END DO
                                   END IF
                                END DO
                             END DO
                          END IF
                          ijh = 1
                          !
                          DO ih = 1, nh(np)
                             !
                             ikb = ijkb0 + ih
                             !
                             IF (upf(np)%has_so) THEN
                                IF (spin_component==1) &
                                  becsum(ijh,na,2)=becsum(ijh,na,2)+ &
                                  (be1(ih,2)*be2(ih,1)+ be1(ih,1)*be2(ih,2))
                                IF (spin_component==2) &
                                   becsum(ijh,na,3)=becsum(ijh,na,3)+ &
                                         (0.d0,-1.d0)*      &  
                                  (be1(ih,2)*be2(ih,1)-be1(ih,1)*be2(ih,2))
                                IF (spin_component==3) &
                                   becsum(ijh,na,4)=becsum(ijh,na,4)+ &
                                  (be1(ih,1)*be2(ih,1)-be1(ih,2)*be2(ih,2))
                             ELSE
                                IF (spin_component==1) &
                                   becsum(ijh,na,2)=becsum(ijh,na,2)  &
                                     + (CONJG(becp_nc(ikb,2,ibnd))   &
                                                 *becp_nc(ikb,1,ibnd)   &
                                     +     CONJG(becp_nc(ikb,1,ibnd))   &
                                                 *becp_nc(ikb,2,ibnd) )
                                IF (spin_component==2) &
                                   becsum(ijh,na,3)=becsum(ijh,na,3)+2.d0   &
                                      *AIMAG(CONJG(becp_nc(ikb,1,ibnd))* &
                                                   becp_nc(ikb,2,ibnd) )
                                IF (spin_component==3) &
                                   becsum(ijh,na,4) = becsum(ijh,na,4)    &
                                          + ( CONJG(becp_nc(ikb,1,ibnd)) &
                                                   *becp_nc(ikb,1,ibnd)  &
                                          -      CONJG(becp_nc(ikb,2,ibnd)) &
                                                      *becp_nc(ikb,2,ibnd) )
                             END IF
                             !
                             ijh = ijh + 1
                             !
                             DO jh = ( ih + 1 ), nh(np)
                                !
                                jkb = ijkb0 + jh
                                !
                                IF (upf(np)%has_so) THEN
                                   IF (spin_component==1) &
                                   becsum(ijh,na,2)=becsum(ijh,na,2)+( &
                                     (be1(jh,2)*be2(ih,1)+be1(jh,1)*be2(ih,2))+&
                                     (be1(ih,2)*be2(jh,1)+be1(ih,1)*be2(jh,2)))
                                   IF (spin_component==2) &
                                   becsum(ijh,na,3)=becsum(ijh,na,3)+ &
                                      (0.d0,-1.d0)*((be1(jh,2)*&
                                      be2(ih,1)-be1(jh,1)*be2(ih,2))+ &
                                      (be1(ih,2)*be2(jh,1)-be1(ih,1)*be2(jh,2)))
                                   IF (spin_component==3) &
                                   becsum(ijh,na,4)=becsum(ijh,na,4)+ &
                                      ((be1(jh,1)*be2(ih,1)- &
                                        be1(jh,2)*be2(ih,2))+  &
                                       (be1(ih,1)*be2(jh,1)-  &
                                        be1(ih,2)*be2(jh,2)) )
                                ELSE
                                   IF (spin_component==1) &
                                   becsum(ijh,na,2)=becsum(ijh,na,2)+ 2.d0* &
                                           DBLE(CONJG(becp_nc(ikb,2,ibnd))* &
                                                      becp_nc(jkb,1,ibnd) + &
                                                CONJG(becp_nc(ikb,1,ibnd))* &
                                                      becp_nc(jkb,2,ibnd) )
                                   IF (spin_component==2) &
                                   becsum(ijh,na,3)=becsum(ijh,na,3)+ &
                                                       2.d0* &
                                            AIMAG(CONJG(becp_nc(ikb,1,ibnd))* &
                                                        becp_nc(jkb,2,ibnd) + &
                                                  CONJG(becp_nc(ikb,1,ibnd))* &
                                                        becp_nc(jkb,2,ibnd) )
                                   IF (spin_component==3) &
                                   becsum(ijh,na,4)=becsum(ijh,na,4)+ 2.d0* &
                                            DBLE(CONJG(becp_nc(ikb,1,ibnd))* &
                                                       becp_nc(jkb,1,ibnd) - &
                                                 CONJG(becp_nc(ikb,2,ibnd))* &
                                                       becp_nc(jkb,2,ibnd) )
                                END IF
                                !            
                                ijh = ijh + 1
                                !
                             END DO
                             !
                          END DO
                          !
                          ijkb0 = ijkb0 + nh(np)
                          !
                       END IF
                       !
                    END DO
                    !
                 ELSE
                    !
                    DO na = 1, nat
                       !
                       IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
                       !
                    END DO
                    !
                 END IF
                 !
              END DO
              !
           END IF
           !
        END DO
        !
     END IF
     !
  END DO 
     !
  IF ( doublegrid ) THEN
    is=spin_component+1
    CALL interpolate( rho%of_r(1,is), rho%of_r(1,is), 1 )
  END IF
  !
  ! ... Here we add the Ultrasoft contribution to the charge and magnetization
  !
  IF ( okvan ) CALL addusdens()

  DO ir=1,nrxx
     raux(ir)=rho%of_r(ir,spin_component+1)
  END DO
  !
  IF (lspinorb) DEALLOCATE(be1, be2)
  DEALLOCATE( becp_nc )
  RETURN
  !
END SUBROUTINE local_dos_mag
