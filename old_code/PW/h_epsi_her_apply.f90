!
! Copyright (C) 2005 Paolo Umari
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"
!-----------------------------------------------------------------------
subroutine h_epsi_her_apply(lda, n,nbande, psi, hpsi, pdir, e_field)
  !-----------------------------------------------------------------------
  !
  ! this subroutine applies w_k+w_k* on psi, 
  ! (as in Souza et al.  PRB B 69, 085106 (2004))
  ! the output is put into hpsi
  !
  ! evcel must contain the wavefunctions from previous iteration
  ! spin polarized systems supported only with fixed occupations

  USE kinds,    ONLY : DP
  USE us
  USE wvfct,    ONLY : igk, g2kin, npwx, npw, nbnd, ik => current_k
  USE gsmooth,  ONLY : nls, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE ldaU,     ONLY : lda_plus_u
  USE lsda_mod, ONLY : current_spin, nspin
  USE scf,      ONLY : vrs  
  USE gvect
  USE uspp
  USE uspp_param, ONLY: nh, nhm, nbetam
  USE bp
  USE basis
  USE klist
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
  USE ions_base, ONLY: ityp, tau, nat,ntyp => nsp
  USE constants, ONLY : e2, pi, tpi, fpi
  USE fixed_occ
  USE io_global, ONLY : stdout
  USE becmod,    ONLY : calbec
  USE mp_global, ONLY : intra_pool_comm
  USE mp,        ONLY : mp_sum
  !
  implicit none
  INTEGER, INTENT(in) :: pdir!direction on which the polarization is calculated
  REAL(DP) :: e_field!electric field along pdir    

  !
  INTEGER :: lda !leading dimension
  INTEGER ::  n! total number of wavefunctions 
  INTEGER :: nbande!number of wavefunctions to be calculated 
  
  COMPLEX(DP) :: psi (lda, nbande ), hpsi (lda,nbande)


  COMPLEX(DP), EXTERNAL :: zdotc
  
  COMPLEX(DP), ALLOCATABLE  :: evct(:,:)!temporary array
  COMPLEX(DP) :: ps(nkb,nbnd)

  COMPLEX(DP) :: becp0(nkb,nbnd)

  INTEGER :: nkbtona(nkb)
   INTEGER :: nkbtonh(nkb)
  REAL(DP) :: g2kin_bp(npwx)
  INTEGER :: igk1(npwx)
  COMPLEX(DP) :: sca, sca1, pref
  INTEGER nb,mb, jkb, nhjkb, na, np, nhjkbm,jkb1,i,j
  INTEGER :: jkb_bp,nt,ig, ijkb0,ibnd,jh,ih,ikb
  

  if(e_field==0.d0) return
  
  ALLOCATE( evct(npwx,nbnd))

  if(okvan) then
!  --- Initialize arrays ---
     jkb_bp=0
      DO nt=1,ntyp
         DO na=1,nat
            IF (ityp(na)== nt) THEN
               DO i=1, nh(nt)
                  jkb_bp=jkb_bp+1
                  nkbtona(jkb_bp) = na
                  nkbtonh(jkb_bp) = i
               END DO
            END IF
         END DO
      END DO
   endif



  if(okvan) THEN
     CALL calbec ( npw, vkb, psi, becp0, nbande )
  endif



  do nb=1,nbande
      
!apply w_k     
     do mb=1,nbnd!index on states of evcel
        sca = zdotc(npw,evcel(1,mb),1,psi(1,nb),1)
        call mp_sum( sca, intra_pool_comm )

        

        if(okvan) then 
           pref = (0.d0,0.d0)
           
           DO jkb=1,nkb
              nhjkb = nkbtonh(jkb)
              na = nkbtona(jkb)
              np = ityp(na)
              nhjkbm = nh(np)
              jkb1 = jkb - nhjkb
              DO j = 1,nhjkbm
                 ! bec_evcel is relative to ik
                 pref = pref+CONJG(bec_evcel(jkb,mb))*becp0(jkb1+j,nb) &
                      *qq(nhjkb,j,np)
              ENDDO
           ENDDO
           sca= sca + pref
        endif


        do ig=1,npw

           hpsi(ig,nb) = hpsi(ig,nb) + &
                &     fact_hepsi(ik,pdir)*sca*(evcelm(ig,mb,pdir)-evcelp(ig,mb,pdir))
        enddo
     enddo
!apply w_k*

     if(.not.okvan) then
        do mb=1,nbnd!index on states of evcel        
           sca = zdotc(npw,evcelm(1,mb,pdir),1,psi(1,nb),1)
           sca1 = zdotc(npw,evcelp(1,mb,pdir),1,psi(1,nb),1)
           call mp_sum( sca, intra_pool_comm )
           call mp_sum(  sca1, intra_pool_comm )
           
           do ig=1,npw

              hpsi(ig,nb) = hpsi(ig,nb) + &
                   &     CONJG(fact_hepsi(ik,pdir))*evcel(ig,mb)*(sca-sca1)
           enddo
        enddo
   
     else ! US case

! copy evcel into evct
        do ig=1,npwx*nbnd
           evct(ig,1)=evcel(ig,1)
        enddo
!  calculate S|evct>
 
        ps (:,:) = (0.d0, 0.d0)
        ijkb0 = 0
        do nt = 1, ntyp
           do na = 1, nat
              if (ityp (na) == nt) then
                 do ibnd = 1, nbnd
                    do jh = 1, nh (nt)
                       jkb = ijkb0 + jh
                       do ih = 1, nh (nt)
                          ikb = ijkb0 + ih
                           ps (ikb, ibnd) = ps (ikb, ibnd) + &
                               qq(ih,jh,nt)* bec_evcel(jkb,ibnd)

                       enddo
                    enddo
                 enddo
                 ijkb0 = ijkb0 + nh (nt)
              endif
           enddo
        enddo
        call ZGEMM ('N', 'N', npw, nbnd , nkb, (1.d0, 0.d0) , vkb, &!vkb is relative to the last ik read
             npwx, ps, nkb, (1.d0, 0.d0) , evct, npwx)
        do mb=1,nbnd!index on states of evcel       
           sca = zdotc(npw,evcelm(1,mb,pdir),1,psi(1,nb),1)
           sca1 = zdotc(npw,evcelp(1,mb,pdir),1,psi(1,nb),1)         
           call mp_sum( sca, intra_pool_comm )
           call mp_sum( sca1, intra_pool_comm )

           do ig=1,npw

              hpsi(ig,nb) = hpsi(ig,nb) + &
                   &     CONJG(fact_hepsi(ik,pdir))*evct(ig,mb)*(sca-sca1)
           enddo
        enddo

     endif
  ENDDO

  DEALLOCATE( evct)



  
!  --
!------------------------------------------------------------------------------!
   return
 END SUBROUTINE h_epsi_her_apply
