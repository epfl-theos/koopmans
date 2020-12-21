!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine qqberry2( gqq,gqqm, ipol)

!  this subroutine computes the array gqq and gqqm
!  gqq=int_dr qq(r)exp(iGr)=<Beta_r|exp(iGr)|Beta_r'>
!  gqqm=int_dr qq(r)exp(-iGr)=<Beta_r|exp(-iGr)|Beta_r'>

!   gqq output: as defined above

  use uspp_param,         only: upf, lmaxq, nbetam, nh, nhm, oldvan
  use uspp,               only: indv, lpx, lpl, ap,nhtolm
  use atom,               only: rgrid
  use core
  use gvecw,              only: ngw
  use reciprocal_vectors, only: mill_l
  use constants
  use cvan,               only: nvb
  use ions_base
  use ions_base,          only: nas => nax
  use cell_base,          only: a1, a2, a3
  use reciprocal_vectors, only: ng0 => gstart, gx, g
  use mp,                 only: mp_sum
  use mp_global,          only: intra_image_comm
  use cp_interfaces,      only: fill_qrl
  
  implicit none

  complex(8) gqq(nhm,nhm,nas,nsp)
  complex(8) gqqm(nhm,nhm,nas,nsp)
  real(8) gmes
  integer :: ipol

! local variables

  integer :: ndm, ig, is, iv, jv, i, istart, il,l,ir, igi,ia
  real(8), allocatable:: fint(:),jl(:)
  real(8), allocatable:: qrl(:,:,:), qradb2(:,:,:,:) 
  real(8) c, xg
  complex(8) qgbs,sig
  integer :: ivs, jvs, ivl, jvl, lp, ijv
  real(8), allocatable:: ylm(:,:)

  IF( .NOT. ALLOCATED( rgrid ) ) &
     CALL errore( ' qqberry2 ', ' rgrid not allocated ', 1 )
  IF( .NOT. ALLOCATED( upf ) ) &
     CALL errore( ' qqberry2 ', ' upf not allocated ', 1 )

  ndm = MAXVAL (upf(1:nsp)%kkbeta)
  allocate( fint( ndm), jl(ndm))
  allocate( qradb2(nbetam,nbetam,lmaxq,nsp))
  allocate( ylm(ngw, lmaxq*lmaxq))

  CALL ylmr2( lmaxq*lmaxq, ngw, gx, g, ylm )

  qradb2 = 0.0d0
     
  do is=1,nsp
     do ia=1,nas
        do jv=1,nhm
           do iv=1,nhm
              gqq(iv,jv,ia,is)=(0.d0,0.d0)
              gqqm(iv,jv,ia,is)=(0.d0,0.d0)
           enddo
        enddo
     enddo
  enddo
  
  if(ipol.eq.1) then
     gmes=a1(1)**2+a1(2)**2+a1(3)**2
     gmes=2*pi/SQRT(gmes)
  endif
  if(ipol.eq.2) then
     gmes=a2(1)**2+a2(2)**2+a2(3)**2
     gmes=2*pi/SQRT(gmes)
  endif
  if(ipol.eq.3) then
     gmes=a3(1)**2+a3(2)**2+a3(3)**2
     gmes=2*pi/SQRT(gmes)
  endif    
  ! only for Vanderbilt species 
  do is=1,nvb
     c=fpi                 !/omegab
     !
     ALLOCATE ( qrl( upf(is)%kkbeta, upf(is)%nbeta*(upf(is)%nbeta+1)/2, &
                     upf(is)%nqlc ) )
     !
     call fill_qrl ( is, qrl )
     ! now the radial part
     do l=1,upf(is)%nqlc
        xg= gmes !only orthorombic cells
        !!!call bess(xg,l,upf(is)%kkbeta,rgrid(is)%r,jl)
        call sph_bes ( upf(is)%kkbeta, rgrid(is)%r, xg, l-1, jl )
        do iv= 1,upf(is)%nbeta
           do jv=iv,upf(is)%nbeta
              ijv = (jv-1)*jv/2 + iv
!     
!     note qrl(r)=r^2*q(r)
!
              do ir=1,upf(is)%kkbeta
                 fint(ir)=qrl(ir,ijv,l)*jl(ir)
              end do
              if (oldvan(is)) then
                 call herman_skillman_int ( upf(is)%kkbeta,fint,rgrid(is)%rab,&
                                            qradb2(iv,jv,l,is) )
              else
                 call simpson ( upf(is)%kkbeta,fint,rgrid(is)%rab,&
                                qradb2(iv,jv,l,is) )
              endif
              qradb2(iv,jv,l,is)=  c*qradb2(iv,jv,l,is)
              if ( iv /= jv ) qradb2(jv,iv,l,is)=  qradb2(iv,jv,l,is)
           end do
        end do
     end do
     DEALLOCATE ( qrl )    
  enddo

  igi=-1
  do ig=1,ngw
     if(ipol.eq.1 ) then
        if(mill_l(1,ig).eq.1 .and. mill_l(2,ig).eq.0  .and. mill_l(3,ig).eq. 0) igi=ig
     endif
     if(ipol.eq.2 ) then
        if(mill_l(1,ig).eq.0 .and. mill_l(2,ig).eq.1  .and. mill_l(3,ig).eq. 0) igi=ig
     endif
     if(ipol.eq.3 ) then
        if(mill_l(1,ig).eq.0 .and. mill_l(2,ig).eq.0   .and. mill_l(3,ig).eq. 1) igi=ig
     endif
  enddo
  if( igi.ne.-1) then

!setting array beigr
             
     do is=1,nvb
        do iv= 1,nh(is)
           do jv=iv,nh(is)
              ivs=indv(iv,is)
              jvs=indv(jv,is)
              ivl=nhtolm(iv,is)
              jvl=nhtolm(jv,is)
!
!     lpx = max number of allowed y_lm
!     lp  = composite lm to indentify them
!
              qgbs=(0.d0,0.d0)
              do i=1,lpx(ivl,jvl)
                 lp=lpl(ivl,jvl,i)
!
!     extraction of angular momentum l from lp:  
!
                 if (lp.eq.1) then
                    l=1         
                 else if ((lp.ge.2) .and. (lp.le.4)) then
                    l=2
                 else if ((lp.ge.5) .and. (lp.le.9)) then
                    l=3
                 else if ((lp.ge.10).and.(lp.le.16)) then
                    l=4
                 else if ((lp.ge.17).and.(lp.le.25)) then
                    l=5
                 else if (lp.ge.26) then 
                    call errore(' qvanb ',' lp.ge.26 ',lp)
                 endif
!
!       sig= (-i)^l
!
                 sig=(0.d0,-1.d0)**(l-1)
                  
                 sig=sig*ap(lp,ivl,jvl)
                 qgbs=qgbs+sig*ylm(igi,lp)*qradb2(ivs,jvs,l,is)
                
              end do
              
              do ia=1,na(is)
                     
                 gqqm(iv,jv,ia,is)=qgbs
                 gqqm(jv,iv,ia,is)=qgbs
                 gqq(iv,jv,ia,is)=CONJG(gqqm(iv,jv,ia,is))
                 gqq(jv,iv,ia,is)=CONJG(gqqm(iv,jv,ia,is))
              end do
           end do
        enddo
     enddo
  endif

  call mp_sum(gqq(:,:,:,:),intra_image_comm)
  call mp_sum(gqqm(:,:,:,:),intra_image_comm)

  deallocate( fint)
  deallocate( jl)
  deallocate(qradb2)
  deallocate(ylm)
  
  return
end subroutine qqberry2






! this subroutine updates gqq and gqqm to the 
! (new) atomic position


subroutine qqupdate(eigr, gqqm0, gqq, gqqm, ipol)

!   gqq output: as defined above

  use cvan
  use gvecw, only: ngw
  use ions_base, only : nas => nax, nat, na, nsp
  use reciprocal_vectors, only: mill_l
  use uspp_param, only: nh, nhm
  use mp, only: mp_sum
  use mp_global, only: intra_image_comm

  implicit none

 
  complex(8) eigr(ngw,nat)
  complex(8) gqq(nhm,nhm,nas,nsp)
  complex(8) gqqm(nhm,nhm,nas,nsp)
  complex(8) gqqm0(nhm,nhm,nas,nsp)

  integer ipol
  
  integer igi,ig,is,iv,jv,ia,isa


  do is=1,nsp
     do ia=1,nas
        do jv=1,nhm
           do iv=1,nhm
              gqq(iv,jv,ia,is)=(0.d0,0.d0)
              gqqm(iv,jv,ia,is)=(0.d0,0.d0)
           enddo
        enddo
     enddo
  enddo

  igi=-1
  do ig=1,ngw
     if(ipol.eq.1 ) then
        if(mill_l(1,ig).eq.1 .and. mill_l(2,ig).eq.0  .and. mill_l(3,ig).eq. 0) igi=ig
     endif
     if(ipol.eq.2 ) then
        if(mill_l(1,ig).eq.0 .and. mill_l(2,ig).eq.1  .and. mill_l(3,ig).eq. 0) igi=ig
     endif
     if(ipol.eq.3 ) then
        if(mill_l(1,ig).eq.0 .and. mill_l(2,ig).eq.0  .and. mill_l(3,ig).eq. 1) igi=ig
     endif
  enddo
  if( igi.ne.-1) then

  
     isa = 1
     do is=1,nvb
        do ia=1,na(is)
           do iv= 1,nh(is)
              do jv=iv,nh(is)
                 gqqm(iv,jv,ia,is)= gqqm0(iv,jv,ia,is)*eigr(igi,isa)
                 gqqm(jv,iv,ia,is)= gqqm0(iv,jv,ia,is)*eigr(igi,isa)
                 gqq(iv,jv,ia,is)=CONJG(gqqm(iv,jv,ia,is))
                 gqq(jv,iv,ia,is)=CONJG(gqqm(iv,jv,ia,is))
              enddo
           enddo
           isa = isa + 1
        enddo
     enddo
  endif
  call mp_sum(gqq(:,:,:,:),intra_image_comm)
  call mp_sum(gqqm(:,:,:,:),intra_image_comm)
  return
end subroutine qqupdate
