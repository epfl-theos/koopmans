!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine qmatrixd(c0, bec0,ctable, gqq, qmat, detq, ipol)

  ! this subroutine computes the inverse of the matrix Q
  ! Q_ij=<Psi_i^0|e^iG_ipol.r|Psi_j^0>
  ! and det Q
  ! Matrix Q is symmetric, and we make us of it

  !    c0 input: the unperturbed wavefunctions
  !    bec0 input: the coefficients <Phi_Rj|Psi_i^0>
  !    ctable input: the coorespondence array
  !    gqq input: the intqq(r) exp(iG_ipol*r) array
  !    qmat output: the inverse q matrix
  !    detq output: det Q
  !    ipol: electric field direction

  use kinds, only : DP
  use gvecw, only: ngw
  use cvan,  only: nvb, ish
  use ions_base, only : nax, nsp, na
  use reciprocal_vectors, only: gstart
  use uspp_param, only: nh, nhm
  use electrons_base, only: nx => nbspx, n => nbsp, ispin
  use mp, only: mp_sum, mp_alltoall
  use mp_global, only: intra_image_comm, nproc_image
  USE efield_module, ONLY : ctable_missing_1,ctable_missing_2, n_g_missing_p,&
       &      ctable_missing_rev_1,ctable_missing_rev_2
  use twin_types !added:giovanni
  
  implicit none
  
  type(twin_matrix) :: bec0 !modified:giovanni (nkb,n)
  complex(DP) ::  gqq(nhm,nhm,nax,nsp)
  complex(DP) :: c0(ngw,nx),  qmat(nx,nx), detq
  integer :: ctable(ngw,2)
  integer, intent(in) :: ipol
  ! local variables
  integer ig,ix,jx, iv,jv,is,ia, inl,jnl, ip
  complex(DP) :: sca
  integer :: info
  integer, allocatable :: ipiv(:,:)
  complex(DP), allocatable :: work(:)
  complex(DP), allocatable :: sndbuf(:,:,:),rcvbuf(:,:,:)


  qmat(:,:)=CMPLX(0.d0,0.d0)

  ALLOCATE( ipiv( nx, nx ), work( nx ) )

  do ix=1,n
     do jx=ix,n

! first the local part

        sca=CMPLX(0.d0,0.d0)
        if(ispin(ix) == ispin(jx) ) then
       
           do ig=1,ngw
              if(ctable(ig,1).ne.(ngw+1))then
                 if(ctable(ig,1).ge.0) then
                    sca=sca+CONJG(c0(ctable(ig,1),ix))*c0(ig,jx)
                 endif
              endif
           enddo

          do ig=1,ngw
             if(ctable(ig,1).ne.(ngw+1))then
                if(ctable(ig,1).lt. 0) then
                   sca=sca+c0(-ctable(ig,1),ix)*c0(ig,jx)
                endif
             endif
          enddo

          do ig=gstart,ngw
             if(ctable(ig,2).ne.(ngw+1)) then
                if(ctable(ig,2).lt.0) then
                   sca=sca+c0(-ctable(ig,2),ix)*CONJG(c0(ig,jx))
                endif
             endif
          enddo
          
          do ig=gstart,ngw
             if(ctable(ig,2).ne.(ngw+1)) then
                if(ctable(ig,2).ge.0) then
                   sca=sca+CONJG(c0(ctable(ig,2),ix))*conjg(c0(ig,jx))
                endif
             endif
          enddo

        

#ifdef __PARA

          if(ipol /= 3) then
             !
             allocate(sndbuf(n_g_missing_p(ipol),2,nproc_image))
             sndbuf(:,:,:)=(0.d0,0.d0)
             allocate(rcvbuf(n_g_missing_p(ipol),2,nproc_image))
!copy arrays to snd buf
             do ip=1,nproc_image
                do ig=1,n_g_missing_p(ipol)
                   if(ipol==1) then
                      if(ctable_missing_1(ig,1,ip)/=0) then
                         sndbuf(ig,1,ip)=c0(ctable_missing_1(ig,1,ip),jx)
                      endif
                   else
                      if(ctable_missing_2(ig,1,ip)/=0) then
                         sndbuf(ig,1,ip)=c0(ctable_missing_2(ig,1,ip),jx)
                      endif
                   endif
                enddo
                do ig=1,n_g_missing_p(ipol)
                   if(ipol==1) then
                      if(ctable_missing_1(ig,2,ip)/=0) then
                         sndbuf(ig,2,ip)=conjg(c0(ctable_missing_1(ig,2,ip),jx))
                      endif
                   else
                      if(ctable_missing_2(ig,2,ip)/=0) then
                         sndbuf(ig,2,ip)=conjg(c0(ctable_missing_2(ig,2,ip),jx))
                      endif
                   endif
                enddo
             enddo


             CALL mp_alltoall( sndbuf, rcvbuf, intra_image_comm )

!update sca
             do ip=1,nproc_image
                do ig=1,n_g_missing_p(ipol)
                   if(ipol==1) then
                      if(ctable_missing_rev_1(ig,1,ip) >0) then
                         sca=sca+conjg(c0(ctable_missing_rev_1(ig,1,ip),ix))*rcvbuf(ig,1,ip)
                      else if(ctable_missing_rev_1(ig,1,ip)< 0) then
                         sca=sca+c0(-ctable_missing_rev_1(ig,1,ip),ix)*rcvbuf(ig,1,ip)
                      endif
                   else
                      if(ctable_missing_rev_2(ig,1,ip) >0) then
                         sca=sca+conjg(c0(ctable_missing_rev_2(ig,1,ip),ix))*rcvbuf(ig,1,ip)
                      else if(ctable_missing_rev_2(ig,1,ip)< 0) then
                         sca=sca+c0(-ctable_missing_rev_2(ig,1,ip),ix)*rcvbuf(ig,1,ip)
                      endif
                   endif
                enddo
                do ig=1,n_g_missing_p(ipol)
                   if(ipol==1) then
                      if(ctable_missing_rev_1(ig,2,ip) >0) then
                         sca=sca+conjg(c0(ctable_missing_rev_1(ig,2,ip),ix))*rcvbuf(ig,2,ip)
                      else if(ctable_missing_rev_1(ig,2,ip)< 0) then
                         sca=sca+c0(-ctable_missing_rev_1(ig,2,ip),ix)*rcvbuf(ig,2,ip)
                      endif
                   else
                      if(ctable_missing_rev_2(ig,2,ip) >0) then
                         sca=sca+conjg(c0(ctable_missing_rev_2(ig,2,ip),ix))*rcvbuf(ig,2,ip)
                      else if(ctable_missing_rev_2(ig,2,ip)< 0) then
                         sca=sca+c0(-ctable_missing_rev_2(ig,2,ip),ix)*rcvbuf(ig,2,ip)
                      endif
                   endif
                enddo

             enddo
            
             deallocate(rcvbuf,sndbuf)

          endif

#endif

        
          call mp_sum( sca, intra_image_comm )
       endif
       qmat(ix,jx)=sca
       
       !  now the non local vanderbilt part
       
       sca = CMPLX(0.d0,0.d0)
       if(ispin(ix)==ispin(jx)) then
          if(.not.bec0%iscmplx) then
	    do is=1,nvb!loop on vanderbilt species
	      do ia=1,na(is)!loop on atoms
		  do iv=1,nh(is)!loop on projectors
		    do jv=1,nh(is)
			inl=ish(is)+(iv-1)*na(is)+ia
			jnl=ish(is)+(jv-1)*na(is)+ia                
			sca=sca+gqq(iv,jv,ia,is)*bec0%rvec(inl,ix)*bec0%rvec(jnl,jx)
		    enddo
		  enddo
	      enddo
	    enddo
          else
            do is=1,nvb!loop on vanderbilt species
	      do ia=1,na(is)!loop on atoms
		  do iv=1,nh(is)!loop on projectors
		    do jv=1,nh(is)
			inl=ish(is)+(iv-1)*na(is)+ia
			jnl=ish(is)+(jv-1)*na(is)+ia                
			sca=sca+gqq(iv,jv,ia,is)*CONJG(bec0%cvec(inl,ix))*bec0%cvec(jnl,jx)
		    enddo
		  enddo
	      enddo
	    enddo
          endif

          qmat(ix,jx)=qmat(ix,jx)+sca
       endif
       qmat(jx,ix)=qmat(ix,jx)       
     enddo
  enddo

  !LAPACK
  call zgetrf (n,n,qmat,nx,ipiv,info)
  !  write(6,*) 'info trf', info
  detq=(1.d0,0.d0)
 
  do ix=1,n
     if(ix.ne.ipiv(ix,1)) detq=-detq
  enddo
  do ix=1,n
     detq = detq*qmat(ix,ix)
  enddo

  call zgetri (n,qmat,nx,ipiv,work,nx,info)

! force qmat to be symmetric

  do ix=1,n
     do jx=ix+1,n
        qmat(jx,ix)=0.5d0*(qmat(ix,jx)+qmat(jx,ix))
        qmat(ix,jx)=qmat(jx,ix)
     enddo
  enddo

  deallocate( ipiv, work )

  return
end subroutine qmatrixd

subroutine qmatrixd_old(c0, bec0,ctable, gqq, qmat, detq, ipol)

  ! this subroutine computes the inverse of the matrix Q
  ! Q_ij=<Psi_i^0|e^iG_ipol.r|Psi_j^0>
  ! and det Q
  ! Matrix Q is symmetric, and we make us of it

  !    c0 input: the unperturbed wavefunctions
  !    bec0 input: the coefficients <Phi_Rj|Psi_i^0>
  !    ctable input: the coorespondence array
  !    gqq input: the intqq(r) exp(iG_ipol*r) array
  !    qmat output: the inverse q matrix
  !    detq output: det Q
  !    ipol: electric field direction

  use kinds, only : DP
  use gvecw, only: ngw
  use cvan,  only: nvb, ish
  use ions_base, only : nax, nsp, na
  use reciprocal_vectors, only: gstart
  use uspp_param, only: nh, nhm
  use uspp, only : nkb
  use electrons_base, only: nx => nbspx, n => nbsp, ispin
  use mp, only: mp_sum, mp_alltoall
  use mp_global, only: intra_image_comm, nproc_image
  USE efield_module, ONLY : ctable_missing_1,ctable_missing_2,n_g_missing_p,&
       &      ctable_missing_rev_1,ctable_missing_rev_2
  
  implicit none
  
  real(DP) :: bec0(nkb,n)
  complex(DP) ::  gqq(nhm,nhm,nax,nsp)
  complex(DP) :: c0(ngw,nx),  qmat(nx,nx), detq
  integer :: ctable(ngw,2)
  integer, intent(in) :: ipol
  ! local variables
  integer ig,ix,jx, iv,jv,is,ia, inl,jnl, ip
  complex(DP) :: sca
  integer :: info
  integer, allocatable :: ipiv(:,:)
  complex(DP), allocatable :: work(:)
  complex(DP), allocatable :: sndbuf(:,:,:),rcvbuf(:,:,:)


  qmat(:,:)=CMPLX(0.d0,0.d0)

  ALLOCATE( ipiv( nx, nx ), work( nx ) )

  do ix=1,n
     do jx=ix,n

! first the local part

        sca=CMPLX(0.d0,0.d0)
        if(ispin(ix) == ispin(jx) ) then
       
           do ig=1,ngw
              if(ctable(ig,1).ne.(ngw+1))then
                 if(ctable(ig,1).ge.0) then
                    sca=sca+CONJG(c0(ctable(ig,1),ix))*c0(ig,jx)
                 endif
              endif
           enddo

          do ig=1,ngw
             if(ctable(ig,1).ne.(ngw+1))then
                if(ctable(ig,1).lt. 0) then
                   sca=sca+c0(-ctable(ig,1),ix)*c0(ig,jx)
                endif
             endif
          enddo

          do ig=gstart,ngw
             if(ctable(ig,2).ne.(ngw+1)) then
                if(ctable(ig,2).lt.0) then
                   sca=sca+c0(-ctable(ig,2),ix)*CONJG(c0(ig,jx))
                endif
             endif
          enddo
          
          do ig=gstart,ngw
             if(ctable(ig,2).ne.(ngw+1)) then
                if(ctable(ig,2).ge.0) then
                   sca=sca+CONJG(c0(ctable(ig,2),ix))*conjg(c0(ig,jx))
                endif
             endif
          enddo

        

#ifdef __PARA

          if(ipol /= 3) then
             !
             allocate(sndbuf(n_g_missing_p(ipol),2,nproc_image))
             sndbuf(:,:,:)=(0.d0,0.d0)
             allocate(rcvbuf(n_g_missing_p(ipol),2,nproc_image))
!copy arrays to snd buf
             do ip=1,nproc_image
                do ig=1,n_g_missing_p(ipol)
                   if(ipol==1) then
                      if(ctable_missing_1(ig,1,ip)/=0) then
                         sndbuf(ig,1,ip)=c0(ctable_missing_1(ig,1,ip),jx)
                      endif
                   else
                      if(ctable_missing_2(ig,1,ip)/=0) then
                         sndbuf(ig,1,ip)=c0(ctable_missing_2(ig,1,ip),jx)
                      endif
                   endif
                enddo
                do ig=1,n_g_missing_p(ipol)
                   if(ipol==1) then
                      if(ctable_missing_1(ig,2,ip)/=0) then
                         sndbuf(ig,2,ip)=conjg(c0(ctable_missing_1(ig,2,ip),jx))
                      endif
                   else
                      if(ctable_missing_2(ig,2,ip)/=0) then
                         sndbuf(ig,2,ip)=conjg(c0(ctable_missing_2(ig,2,ip),jx))
                      endif
                   endif
                enddo
             enddo


             CALL mp_alltoall( sndbuf, rcvbuf, intra_image_comm )

!update sca
             do ip=1,nproc_image
                do ig=1,n_g_missing_p(ipol)
                   if(ipol==1) then
                      if(ctable_missing_rev_1(ig,1,ip) >0) then
                         sca=sca+conjg(c0(ctable_missing_rev_1(ig,1,ip),ix))*rcvbuf(ig,1,ip)
                      else if(ctable_missing_rev_1(ig,1,ip)< 0) then
                         sca=sca+c0(-ctable_missing_rev_1(ig,1,ip),ix)*rcvbuf(ig,1,ip)
                      endif
                   else
                      if(ctable_missing_rev_2(ig,1,ip) >0) then
                         sca=sca+conjg(c0(ctable_missing_rev_2(ig,1,ip),ix))*rcvbuf(ig,1,ip)
                      else if(ctable_missing_rev_2(ig,1,ip)< 0) then
                         sca=sca+c0(-ctable_missing_rev_2(ig,1,ip),ix)*rcvbuf(ig,1,ip)
                      endif
                   endif
                enddo
                do ig=1,n_g_missing_p(ipol)
                   if(ipol==1) then
                      if(ctable_missing_rev_1(ig,2,ip) >0) then
                         sca=sca+conjg(c0(ctable_missing_rev_1(ig,2,ip),ix))*rcvbuf(ig,2,ip)
                      else if(ctable_missing_rev_1(ig,2,ip)< 0) then
                         sca=sca+c0(-ctable_missing_rev_1(ig,2,ip),ix)*rcvbuf(ig,2,ip)
                      endif
                   else
                      if(ctable_missing_rev_2(ig,2,ip) >0) then
                         sca=sca+conjg(c0(ctable_missing_rev_2(ig,2,ip),ix))*rcvbuf(ig,2,ip)
                      else if(ctable_missing_rev_2(ig,2,ip)< 0) then
                         sca=sca+c0(-ctable_missing_rev_2(ig,2,ip),ix)*rcvbuf(ig,2,ip)
                      endif
                   endif
                enddo

             enddo
            
             deallocate(rcvbuf,sndbuf)

          endif

#endif

        
          call mp_sum( sca, intra_image_comm )
       endif
       qmat(ix,jx)=sca
       
       !  now the non local vanderbilt part
       
       sca = CMPLX(0.d0,0.d0)
       if(ispin(ix)==ispin(jx)) then
          do is=1,nvb!loop on vanderbilt species
             do ia=1,na(is)!loop on atoms
                do iv=1,nh(is)!loop on projectors
                   do jv=1,nh(is)
                      inl=ish(is)+(iv-1)*na(is)+ia
                      jnl=ish(is)+(jv-1)*na(is)+ia                
                      sca=sca+gqq(iv,jv,ia,is)*bec0(inl,ix)*bec0(jnl,jx)
                   enddo
                enddo
             enddo
          enddo

          qmat(ix,jx)=qmat(ix,jx)+sca
       endif
       qmat(jx,ix)=qmat(ix,jx)       
     enddo
  enddo

  !LAPACK
  call zgetrf (n,n,qmat,nx,ipiv,info)
  !  write(6,*) 'info trf', info
  detq=(1.d0,0.d0)
 
  do ix=1,n
     if(ix.ne.ipiv(ix,1)) detq=-detq
  enddo
  do ix=1,n
     detq = detq*qmat(ix,ix)
  enddo

  call zgetri (n,qmat,nx,ipiv,work,nx,info)

! force qmat to be symmetric

  do ix=1,n
     do jx=ix+1,n
        qmat(jx,ix)=0.5d0*(qmat(ix,jx)+qmat(jx,ix))
        qmat(ix,jx)=qmat(jx,ix)
     enddo
  enddo

  deallocate( ipiv, work )

  return
end subroutine qmatrixd_old
