!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!-----------------------------------------------------------------------
   subroutine calcmt( fdiag, zmat, fmat, firstiter)
!-----------------------------------------------------------------------
!
!  constructs fmat=z0^t.fdiag.z0    zmat = z0^t
!
      USE kinds,             ONLY: DP
      use electrons_base,    ONLY: nudx, nspin, nupdwn, iupdwn, nx => nbspx
      USE cp_main_variables, ONLY: descla, nrlx
      USE descriptors,       ONLY: la_npc_ , la_npr_ , la_comm_ , la_me_ , la_nrl_ , &
                                   lambda_node_ , ldim_cyclic
      USE mp,                ONLY: mp_sum, mp_bcast

      implicit none
      logical firstiter

      real(DP) :: zmat( nrlx, nudx, nspin ), fmat( nrlx, nudx, nspin ), fdiag( nx )
                  !  NOTE: zmat and fmat are distributed by row across processors
                  !        fdiag is replicated

      integer  :: iss, nss, istart, i, j, k, ii
      integer  :: np_rot, me_rot, nrl, comm_rot, ip, nrl_ip

      real(DP), ALLOCATABLE :: mtmp(:,:)
      real(DP) :: f_z0t


      call start_clock('calcmt')

      fmat = 0.0d0

      DO iss = 1, nspin

         nss      = nupdwn( iss )
         istart   = iupdwn( iss )
         np_rot   = descla( la_npr_  , iss ) * descla( la_npc_ , iss )
         me_rot   = descla( la_me_   , iss )
         nrl      = descla( la_nrl_  , iss )
         comm_rot = descla( la_comm_ , iss )

         IF( descla( lambda_node_ , iss ) > 0 ) THEN

            ALLOCATE( mtmp( nrlx, nudx ) )

            DO ip = 1, np_rot

               IF( me_rot == ( ip - 1 ) ) THEN
                  mtmp = zmat(:,:,iss)
               END IF
               nrl_ip = ldim_cyclic( nss, np_rot, ip - 1 )
               CALL mp_bcast( mtmp , ip - 1 , comm_rot )

               DO j = 1, nss
                  ii = ip
                  DO i = 1, nrl_ip
                     f_z0t = fdiag( j + istart - 1 ) * mtmp( i, j )
                     DO k = 1, nrl
                        fmat( k, ii, iss ) = fmat( k, ii, iss )+ zmat( k, j, iss ) * f_z0t 
                     END DO
                     ii = ii + np_rot
                  END DO
               END DO

            END DO

            DEALLOCATE( mtmp )

         END IF

      END DO

      call stop_clock('calcmt')

      RETURN
      END SUBROUTINE calcmt

!-----------------------------------------------------------------------
   subroutine calcmt_twin( fdiag, zmat, fmat, firstiter)
!-----------------------------------------------------------------------
!
!  constructs fmat=z0^t.fdiag.z0    zmat = z0^t
!
      USE kinds,             ONLY: DP
      use electrons_base,    ONLY: nudx, nspin, nupdwn, iupdwn, nx => nbspx
      USE cp_main_variables, ONLY: descla, nrlx
      USE descriptors,       ONLY: la_npc_ , la_npr_ , la_comm_ , la_me_ , la_nrl_ , &
                                   lambda_node_ , ldim_cyclic
      USE mp,                ONLY: mp_sum, mp_bcast
      USE twin_types

      implicit none
      logical firstiter

      real(DP) ::  fdiag( nx )
                  !  NOTE: zmat and fmat are distributed by row across processors
                  !        fdiag is replicated
      type(twin_matrix) :: zmat(:), fmat(:)

      integer  :: iss, nss, istart, i, j, k, ii
      integer  :: np_rot, me_rot, nrl, comm_rot, ip, nrl_ip

      real(DP), ALLOCATABLE :: mtmp(:,:)
      complex(DP), ALLOCATABLE :: mtmp_c(:,:)
      real(DP) :: f_z0t
      complex(DP) :: f_z0t_c


      call start_clock('calcmt')

      do iss=1,nspin
	call set_twin(fmat(iss), CMPLX(0.d0,0.d0))
      enddo

      DO iss = 1, nspin

         nss      = nupdwn( iss )
         istart   = iupdwn( iss )
         np_rot   = descla( la_npr_  , iss ) * descla( la_npc_ , iss )
         me_rot   = descla( la_me_   , iss )
         nrl      = descla( la_nrl_  , iss )
         comm_rot = descla( la_comm_ , iss )

         IF( descla( lambda_node_ , iss ) > 0 ) THEN
            
            IF(.not.fmat(iss)%iscmplx) THEN
	      ALLOCATE( mtmp( nrlx, nudx ) )

	      DO ip = 1, np_rot

		IF( me_rot == ( ip - 1 ) ) THEN
		    mtmp = zmat(iss)%rvec(:,:)
		END IF
		nrl_ip = ldim_cyclic( nss, np_rot, ip - 1 )
		CALL mp_bcast( mtmp , ip - 1 , comm_rot )

		DO j = 1, nss
		    ii = ip
		    DO i = 1, nrl_ip
		      f_z0t = fdiag( j + istart - 1 ) * mtmp( i, j )
		      DO k = 1, nrl
			  fmat(iss)%rvec( k, ii) = fmat(iss)%rvec( k, ii)+ zmat(iss)%rvec(k, j) * f_z0t 
		      END DO
		      ii = ii + np_rot
		    END DO
		END DO

	      END DO

	      DEALLOCATE( mtmp )
            ELSE
	      ALLOCATE( mtmp_c( nrlx, nudx ) )

	      DO ip = 1, np_rot

		IF( me_rot == ( ip - 1 ) ) THEN
		    mtmp_c(:,:) = zmat(iss)%cvec(:,:)
		END IF
		nrl_ip = ldim_cyclic( nss, np_rot, ip - 1 )
		CALL mp_bcast( mtmp_c , ip - 1 , comm_rot )

		DO j = 1, nss
		    ii = ip
		    DO i = 1, nrl_ip
		      f_z0t_c = fdiag( j + istart - 1 ) * mtmp_c( i, j )
		      DO k = 1, nrl
			  fmat(iss)%cvec( k, ii) = fmat(iss)%cvec(k, ii)+ zmat(iss)%cvec( k, j) * f_z0t_c 
		      END DO
		      ii = ii + np_rot
		    END DO
		END DO

	      END DO

	      DEALLOCATE( mtmp_c )
            ENDIF

         END IF

      END DO

      call stop_clock('calcmt')

      RETURN
      END SUBROUTINE calcmt_twin

!-----------------------------------------------------------------------
      subroutine rotate( z0, c0, bec, c0diag, becdiag, firstiter )
!-----------------------------------------------------------------------
      use kinds, only: dp
      use cvan
      use electrons_base, only: nudx, nspin, nupdwn, iupdwn, nx => nbspx, n => nbsp
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb
      use gvecw, only: ngw
      use ions_base, only: nsp, na
      USE cp_main_variables, ONLY: descla, nrlx
      USE descriptors,       ONLY: la_npc_ , la_npr_ , la_comm_ , la_me_ , la_nrl_
      USE cp_interfaces,     ONLY: protate

      implicit none
      integer iss, nss, istart
      integer :: np_rot, me_rot, nrl, comm_rot
      real(kind=DP)    z0( nrlx, nudx, nspin )
      real(kind=DP)    bec( nhsa, n ), becdiag( nhsa, n )
      complex(kind=DP) c0( ngw, nx ), c0diag( ngw, nx )
      logical firstiter
  
      CALL start_clock( 'rotate' )

      DO iss = 1, nspin
         istart   = iupdwn( iss )
         nss      = nupdwn( iss )
         np_rot   = descla( la_npr_  , iss ) * descla( la_npc_ , iss )
         me_rot   = descla( la_me_   , iss )
         nrl      = descla( la_nrl_  , iss )
         comm_rot = descla( la_comm_ , iss )
         CALL protate ( c0, bec, c0diag, becdiag, ngw, nss, istart, z0(:,:,iss), &
                        na, nsp, ish, nh, np_rot, me_rot)
      END DO

      CALL stop_clock( 'rotate' )
      return
      end subroutine rotate

!-----------------------------------------------------------------------
      subroutine rotate_twin( z0, c0, bec, c0diag, becdiag, firstiter )
!-----------------------------------------------------------------------
      use kinds, only: dp
      use cvan
      use electrons_base, only: nspin, nupdwn, iupdwn, nx => nbspx
      use uspp_param, only: nh
      use gvecw, only: ngw
      use ions_base, only: nsp, na
      USE cp_main_variables, ONLY: descla
      USE descriptors,       ONLY: la_npc_ , la_npr_ , la_comm_ , la_me_ , la_nrl_
      USE cp_interfaces,     ONLY: protate
      USE twin_types

      implicit none
      integer iss, nss, istart
      integer :: np_rot, me_rot, nrl, comm_rot
!       real(kind=DP)    z0( nrlx, nudx, nspin )
!       real(kind=DP)    bec( nhsa, n ), becdiag( nhsa, n )
      type(twin_matrix) :: z0(:)
      type(twin_matrix) :: bec, becdiag
      complex(kind=DP) c0( ngw, nx ), c0diag( ngw, nx )
      logical firstiter
  
      CALL start_clock( 'rotate' )

      DO iss = 1, nspin
         istart   = iupdwn( iss )
         nss      = nupdwn( iss )
         np_rot   = descla( la_npr_  , iss ) * descla( la_npc_ , iss )
         me_rot   = descla( la_me_   , iss )
         nrl      = descla( la_nrl_  , iss )
         comm_rot = descla( la_comm_ , iss )
         IF(.not.z0(iss)%iscmplx) THEN
	  CALL protate ( c0, bec%rvec, c0diag, becdiag%rvec, ngw, nss, istart, z0(iss)%rvec, &
                        na, nsp, ish, nh, np_rot, me_rot)
         ELSE
	  CALL protate ( c0, bec%cvec, c0diag, becdiag%cvec, ngw, nss, istart, z0(iss)%cvec, &
                        na, nsp, ish, nh, np_rot, me_rot )
         ENDIF
      END DO

      CALL stop_clock( 'rotate' )
      return
      end subroutine rotate_twin

!-----------------------------------------------------------------------
      subroutine ddiag(nx,n,amat,dval,dvec,iflag)
!-----------------------------------------------------------------------
!
      use dspev_module, only: dspev_drv
      use kinds , only : dp

      implicit none

      integer nx,n,ndim,iflag,k,i,j
      real(dp)   dval(n)
      real(dp) amat(nx,n), dvec(nx,n)
      real(dp), allocatable::  ap(:)

      ndim=(n*(n+1))/2
      allocate(ap(ndim))
      ap(:)=0.d0

      k=0
      do j=1,n
       do i=1,j
        k=k+1
        ap(k)=amat(i,j)
       end do
      end do

      CALL dspev_drv( 'V', 'U', n, ap, dval, dvec, nx )

      deallocate(ap)

      return
    end subroutine ddiag


!$$
!-----------------------------------------------------------------------
      subroutine zdiag(nx,n,amat,dval,cvec,iflag)
!-----------------------------------------------------------------------
!
      use zhpev_module, only: zhpev_drv
      use kinds , only : dp

      implicit none

      integer nx,n,ndim,iflag,k,i,j
      real(dp)   dval(n)
      complex(dp) amat(nx,n), cvec(nx,n)
      complex(dp), allocatable::  ap(:)

      ndim=(n*(n+1))/2
      allocate(ap(ndim))
      ap(:)=CMPLX(0.d0,0.d0)

      k=0
      do j=1,n
       do i=1,j
        k=k+1
        ap(k)=amat(i,j)
       end do
      end do

      CALL zhpev_drv( 'V', 'U', n, ap, dval, cvec, nx )

      deallocate(ap)

      return
    end subroutine zdiag
!$$


!-----------------------------------------------------------------------
      subroutine calcm(fdiag,zmat,fmat,firstiter)
!-----------------------------------------------------------------------
!
!  constructs fmat=zmat.fdiag.zmat^t
!
      use electrons_base, only: nudx, nspin, nupdwn, iupdwn, nx => nbspx
      use kinds, only : dp

      implicit none

      logical firstiter


      integer iss, nss, istart, i, j, k
      real(dp) zmat(nudx,nudx,nspin), fmat(nudx,nudx,nspin),         &
    &   fdiag(nx)

      call errore(" calcm ", " subroutine not updated ", 1)

      call start_clock('calcm')


        do iss=1,nspin
         nss=nupdwn(iss)
         istart=iupdwn(iss)
         do i=1,nss
          do k=1,nss
           fmat(k,i,iss)=0.0d0
           do j=1,nss
            fmat(k,i,iss)=fmat(k,i,iss)+                                  &
    &            zmat(k,j,iss)*fdiag(j+istart-1)*zmat(i,j,iss)
           end do
          end do
         end do
        end do

      call stop_clock('calcm')
      return
      end subroutine calcm

    subroutine minparabola(ene0,dene0,ene1,passop,passo,stima)
!this subroutines finds the minimum of a quadratic real function
      
      use kinds, only : dp

      implicit none
      real(dp) ene0,dene0,ene1,passop,passo,stima
      real(dp) a,b,c!a*x^2+b*x+c
      
      c=ene0
      b=dene0
      a=(ene1-b*passop-c)/(passop**2.d0)
      
      passo = -b/(2.d0*a)
      if( a.lt.0.d0) then
         if(ene1.lt.ene0) then
            passo=passop
!$$            passo=passop*2.d0
         else 
!$$ This case never happens.
!$$ (1) b is always negative and hence -b*passop is positive.
!$$ (2) in order for a to be negative, therefore, ene1-c(=ene0).
!$$     should be very negative.
            passo=0.5d0*passop
         endif
      endif


      stima=a*passo**2.d0+b*passo+c


      return
    end subroutine minparabola

subroutine pc2_non_ortho(a, adual, beca, becadual, b,becb, lgam)      
               
! this function applies the operator Pc
            
!    this subroutine applies the Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i-a_j><a_j|S|b_i>
    
      use kinds, only: dp 
      use ions_base, only: na
      use mp_global, only: intra_image_comm
      use cvan 
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, nupdwn, iupdwn, nspin
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb
      use uspp, only :qq
      use parallel_toolkit, only : rep_matmul_drv
      use twin_types ! added:giovanni
      
                           
      implicit none        
                           
      complex(kind=DP) a(ngw,n), adual(ngw,n), b(ngw,n)
                     
      type(twin_matrix) ::   beca,becb,becadual!(nhsa,n) !modified:giovanni
      logical :: lgam
! local variables
      integer is, iv, jv, ia, inl, jnl, i, j
      real(DP), allocatable :: bectmp(:,:)
      complex(DP), allocatable :: bectmp_c(:,:)
      real(DP), allocatable :: qq_tmp(:,:), qqb_tmp(:,:)
      complex(DP), allocatable :: qqb_tmp_c(:,:), qq_tmp_c(:,:)
      complex(DP), allocatable :: zbectmp(:,:)
      integer :: nl_max
      integer :: nss,iss, istart

      logical :: mat_par=.true.!if true uses parallel routines

      CALL start_clock( 'pc2' )

      do iss= 1, nspin
         nss= nupdwn( iss )
!          write(6,*) "nupdwn", iss, nupdwn(iss), iupdwn(iss)
         if(nss>0) THEN
            
            istart= iupdwn( iss )

            if(lgam) then
               allocate(bectmp(nss,nss))
               bectmp(:,:)=0.d0
            else
               allocate(bectmp_c(nss,nss))
               bectmp_c(:,:)=CMPLX(0.d0,0.d0)
            endif
   ! 
            allocate(zbectmp(nss,nss))
            call zgemm('C','N',nss,nss,ngw,(1.d0,0.d0),a(:,istart),ngw,b(:,istart),ngw,(0.d0,0.d0),zbectmp,nss)

            if(lgam) then
               do j=1,nss
                   do i=1,nss
                     bectmp(i,j)=2.d0*DBLE(zbectmp(i,j))
                     if(ng0.eq.2) bectmp(i,j)=bectmp(i,j)-DBLE(a(1,j))*DBLE(b(1,i))
                   enddo
               enddo
               call mp_sum( bectmp(:,:), intra_image_comm)
            else
               do j=1,nss
                   do i=1,nss
                     bectmp_c(i,j)=zbectmp(i,j)
                   enddo
               enddo
               call mp_sum( bectmp_c(:,:), intra_image_comm)
            endif
            deallocate(zbectmp)
            if(nvb >= 0) then

               nl_max=0
               do is=1,nvb
                  nl_max=nl_max+nh(is)*na(is)
               enddo
               if(lgam) then
                   allocate (qqb_tmp(nl_max,nss))
                   allocate (qq_tmp(nl_max,nl_max))
                   qq_tmp(:,:)=0.d0
                   do is=1,nvb
                     do iv=1,nh(is)
                         do jv=1,nh(is)
                           do ia=1,na(is)
                               inl=ish(is)+(iv-1)*na(is)+ia
                               jnl=ish(is)+(jv-1)*na(is)+ia
                               qq_tmp(inl,jnl)=qq(iv,jv,is)
                           enddo
                         enddo
                     enddo
                   enddo
               else
                   allocate (qqb_tmp_c(nl_max,nss))
                   allocate (qq_tmp_c(nl_max,nl_max))
                   qq_tmp_c(:,:)=CMPLX(0.d0,0.d0)
                   do is=1,nvb
                     do iv=1,nh(is)
                         do jv=1,nh(is)
                           do ia=1,na(is)
                               inl=ish(is)+(iv-1)*na(is)+ia
                               jnl=ish(is)+(jv-1)*na(is)+ia
                               qq_tmp_c(inl,jnl)=CMPLX(qq(iv,jv,is),0.d0)
                           enddo
                         enddo
                     enddo
                   enddo
               endif
               !
               if(lgam) then
                   if( nhsa > 0 .and. .not. mat_par)  then
                     call dgemm('N','N',nl_max,nss,nl_max,1.d0,qq_tmp,nl_max,becb%rvec(:,istart),nhsa,0.d0,qqb_tmp,nl_max)
                     call dgemm('T','N',nss,nss,nl_max,1.d0,beca%rvec(:,istart),nhsa,qqb_tmp,nl_max,1.d0,bectmp,nss)
                   else if ( nhsa > 0 ) then
                     call para_dgemm ('N','N',nl_max,nss,nl_max,1.d0,qq_tmp,nl_max,&
                                       becb%rvec(:,istart),nhsa,0.d0,qqb_tmp,nl_max, intra_image_comm)
                     call para_dgemm ('T','N',nss,nss,nl_max,1.d0,beca%rvec(:,istart),nhsa, &
                                       qqb_tmp,nl_max,1.d0,bectmp,nss, intra_image_comm)
                   endif
                  deallocate(qq_tmp,qqb_tmp)
               else
                   if( nhsa > 0 .and. .not. mat_par)  then
                     call zgemm('N','N',nl_max,nss,nl_max,(1.d0,0.d0),qq_tmp_c, &
                          nl_max,becb%cvec(:,istart),nhsa,(0.d0,0.d0), qqb_tmp_c,nl_max)
                     call zgemm('C','N',nss,nss,nl_max,(1.d0,0.d0), & 
                          beca%cvec(:,istart),nhsa,qqb_tmp_c,nl_max,(1.d0,0.d0),bectmp_c,nss)
                   else if ( nhsa > 0 ) then
                     call para_zgemm ('N','N',nl_max,nss,nl_max,(1.d0,0.d0),qq_tmp_c,nl_max,&
                                       becb%cvec(:,istart),nhsa,(0.d0,0.d0),qqb_tmp_c,nl_max, intra_image_comm)
                     call para_zgemm ('C','N',nss,nss,nl_max,(1.d0,0.d0),beca%cvec(:,istart),nhsa, &
                                       qqb_tmp_c,nl_max,(1.d0,0.d0),bectmp_c,nss, intra_image_comm)
                   endif
                  deallocate(qq_tmp_c,qqb_tmp_c)
               endif   
               !
            endif
            allocate(zbectmp(nss,nss))
            if(lgam) then
               do i=1,nss
                   do j=1,nss
                     zbectmp(i,j)=CMPLX(bectmp(i,j),0.d0)
                   enddo
               enddo
            else
               do i=1,nss
                   do j=1,nss
                     zbectmp(i,j)=bectmp_c(i,j)
                   enddo
               enddo
            endif
            call zgemm('N','N',ngw,nss,nss,(-1.d0,0.d0),adual(:,istart),ngw,zbectmp,nss,(1.d0,0.d0),b(:,istart),ngw)
            deallocate(zbectmp)

            ! this computes the new bec
            if(lgam) then
               if ( nhsa > 0 ) then
                   call dgemm('N','N',nhsa,nss,nss,1.0d0,becadual%rvec(:,istart), &
                              nhsa,bectmp,nss,1.0d0,becb%rvec(:,istart),nhsa)
               endif
               deallocate(bectmp)
            else
               if ( nhsa > 0 ) then
                   call zgemm('N','N',nhsa,nss,nss,(1.0d0,0.d0),becadual%cvec(:,istart), &
                              nhsa,bectmp_c,nss,(1.0d0,0.d0),becb%cvec(:,istart),nhsa)
               endif
               deallocate(bectmp_c)
            endif
            !
         ENDIF
      enddo!on spin
      CALL stop_clock( 'pc2' )
      return
    end subroutine pc2_non_ortho

subroutine pc2(a,beca,b,becb, lgam)      
               
! this function applies the operator Pc
            
!    this subroutine applies the Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i-a_j><a_j|S|b_i>
    
      use kinds, only: dp 
      use ions_base, only: na
      use mp_global, only: intra_image_comm
      use cvan 
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, nupdwn, iupdwn, nspin
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb
      use uspp, only :qq
      use parallel_toolkit, only : rep_matmul_drv
      use twin_types ! added:giovanni
      
                           
      implicit none        
                           
      complex(kind=DP) a(ngw,n), b(ngw,n)
                     
      type(twin_matrix) ::   beca,becb!(nhsa,n) !modified:giovanni
      logical :: lgam
! local variables
      integer is, iv, jv, ia, inl, jnl, i, j
      real(DP), allocatable :: bectmp(:,:)
      complex(DP), allocatable :: bectmp_c(:,:)
      real(DP), allocatable :: qq_tmp(:,:), qqb_tmp(:,:)
      complex(DP), allocatable :: qqb_tmp_c(:,:), qq_tmp_c(:,:)
      complex(DP), allocatable :: zbectmp(:,:)
      integer :: nl_max
      integer :: nss,iss, istart

      logical :: mat_par=.true.!if true uses parallel routines

      CALL start_clock( 'pc2' )

      do iss= 1, nspin
         nss= nupdwn( iss )
!          write(6,*) "nupdwn", iss, nupdwn(iss), iupdwn(iss)
         if(nss>0) THEN
            
            istart= iupdwn( iss )

            if(lgam) then
               allocate(bectmp(nss,nss))
               bectmp(:,:)=0.d0
            else
               allocate(bectmp_c(nss,nss))
               bectmp_c(:,:)=CMPLX(0.d0,0.d0)
            endif
   ! 
            allocate(zbectmp(nss,nss))
            call zgemm('C','N',nss,nss,ngw,(1.d0,0.d0),a(:,istart),ngw,b(:,istart),ngw,(0.d0,0.d0),zbectmp,nss)

            if(lgam) then
               do j=1,nss
                   do i=1,nss
                     bectmp(i,j)=2.d0*DBLE(zbectmp(i,j))
                     if(ng0.eq.2) bectmp(i,j)=bectmp(i,j)-DBLE(a(1,j))*DBLE(b(1,i))
                   enddo
               enddo
               call mp_sum( bectmp(:,:), intra_image_comm)
            else
               do j=1,nss
                   do i=1,nss
                     bectmp_c(i,j)=zbectmp(i,j)
                   enddo
               enddo
               call mp_sum( bectmp_c(:,:), intra_image_comm)
            endif
            deallocate(zbectmp)
            if(nvb >= 0) then

               nl_max=0
               do is=1,nvb
                  nl_max=nl_max+nh(is)*na(is)
               enddo
               if(lgam) then
                   allocate (qqb_tmp(nl_max,nss))
                   allocate (qq_tmp(nl_max,nl_max))
                   qq_tmp(:,:)=0.d0
                   do is=1,nvb
                     do iv=1,nh(is)
                         do jv=1,nh(is)
                           do ia=1,na(is)
                               inl=ish(is)+(iv-1)*na(is)+ia
                               jnl=ish(is)+(jv-1)*na(is)+ia
                               qq_tmp(inl,jnl)=qq(iv,jv,is)
                           enddo
                         enddo
                     enddo
                   enddo
               else
                   allocate (qqb_tmp_c(nl_max,nss))
                   allocate (qq_tmp_c(nl_max,nl_max))
                   qq_tmp_c(:,:)=CMPLX(0.d0,0.d0)
                   do is=1,nvb
                     do iv=1,nh(is)
                         do jv=1,nh(is)
                           do ia=1,na(is)
                               inl=ish(is)+(iv-1)*na(is)+ia
                               jnl=ish(is)+(jv-1)*na(is)+ia
                               qq_tmp_c(inl,jnl)=CMPLX(qq(iv,jv,is),0.d0)
                           enddo
                         enddo
                     enddo
                   enddo
               endif
               !
               if(lgam) then
                   if( nhsa > 0 .and. .not. mat_par)  then
                     call dgemm('N','N',nl_max,nss,nl_max,1.d0,qq_tmp,nl_max,becb%rvec(:,istart),nhsa,0.d0,qqb_tmp,nl_max)
                     call dgemm('T','N',nss,nss,nl_max,1.d0,beca%rvec(:,istart),nhsa,qqb_tmp,nl_max,1.d0,bectmp,nss)
                   else if ( nhsa > 0 ) then
                     call para_dgemm ('N','N',nl_max,nss,nl_max,1.d0,qq_tmp,nl_max,&
                                       becb%rvec(:,istart),nhsa,0.d0,qqb_tmp,nl_max, intra_image_comm)
                     call para_dgemm ('T','N',nss,nss,nl_max,1.d0,beca%rvec(:,istart),nhsa, &
                                       qqb_tmp,nl_max,1.d0,bectmp,nss, intra_image_comm)
                   endif
                  deallocate(qq_tmp,qqb_tmp)
               else
                   if( nhsa > 0 .and. .not. mat_par)  then
                     call zgemm('N','N',nl_max,nss,nl_max,(1.d0,0.d0),qq_tmp_c, &
                                nl_max,becb%cvec(:,istart),nhsa,(0.d0,0.d0), qqb_tmp_c,nl_max)
                     call zgemm('C','N',nss,nss,nl_max,(1.d0,0.d0),beca%cvec(:,istart),& 
                                nhsa,qqb_tmp_c,nl_max,(1.d0,0.d0),bectmp_c,nss)
                   else if ( nhsa > 0 ) then
                     call para_zgemm ('N','N',nl_max,nss,nl_max,(1.d0,0.d0),qq_tmp_c,nl_max,&
                                       becb%cvec(:,istart),nhsa,(0.d0,0.d0),qqb_tmp_c,nl_max, intra_image_comm)
                     call para_zgemm ('C','N',nss,nss,nl_max,(1.d0,0.d0),beca%cvec(:,istart),nhsa, &
                                       qqb_tmp_c,nl_max,(1.d0,0.d0),bectmp_c,nss, intra_image_comm)
                   endif
                  deallocate(qq_tmp_c,qqb_tmp_c)
               endif   
               !
            endif
            allocate(zbectmp(nss,nss))
            if(lgam) then
               do i=1,nss
                   do j=1,nss
                     zbectmp(i,j)=CMPLX(bectmp(i,j),0.d0)
                   enddo
               enddo
            else
               do i=1,nss
                   do j=1,nss
                     zbectmp(i,j)=bectmp_c(i,j)
                   enddo
               enddo
            endif
            call zgemm('N','N',ngw,nss,nss,(-1.d0,0.d0),a(:,istart),ngw,&
                       zbectmp,nss,(1.d0,0.d0),b(:,istart),ngw)
            deallocate(zbectmp)

            ! this computes the new bec
            if(lgam) then
               if ( nhsa > 0 ) then
                   call dgemm('N','N',nhsa,nss,nss,1.0d0,beca%rvec(:,istart),nhsa,bectmp,nss,1.0d0,becb%rvec(:,istart),nhsa)
               endif
               deallocate(bectmp)
            else
               if ( nhsa > 0 ) then
                   call zgemm('N','N',nhsa,nss,nss,(1.0d0,0.d0),beca%cvec(:,istart),&
                               nhsa,bectmp_c,nss,(1.0d0,0.d0),becb%cvec(:,istart),nhsa)
               endif
               deallocate(bectmp_c)
            endif
            !
         ENDIF
      enddo!on spin
      CALL stop_clock( 'pc2' )
      return
    end subroutine pc2

    subroutine pc2_new(a,beca,b,becb,n,nupdwn,iupdwn,ispin, lgam)      
               
! this function applies the operator Pc
            
!    this subroutine applies the Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i-a_j><a_j|S|b_i>
    
      use kinds, only: dp 
      use ions_base, only: na
      use mp_global, only: intra_image_comm
      use cvan 
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only:  nspin
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb
      use uspp, only :qq
      use parallel_toolkit, only : rep_matmul_drv
      use twin_types ! added:giovanni
      
                           
      implicit none        

      integer, intent(in) :: n, nupdwn(nspin), iupdwn(nspin), ispin(n)
      complex(kind=DP) a(ngw,n), b(ngw,n)
                     
      type(twin_matrix) ::   beca,becb!(nhsa,n) !modified:giovanni
      logical :: lgam
! local variables
      integer is, iv, jv, ia, inl, jnl, i, j
      real(DP), allocatable :: bectmp(:,:)
      complex(DP), allocatable :: bectmp_c(:,:)
      real(DP), allocatable :: qq_tmp(:,:), qqb_tmp(:,:)
      complex(DP), allocatable :: qqb_tmp_c(:,:), qq_tmp_c(:,:)
      complex(DP), allocatable :: zbectmp(:,:)
      integer :: nl_max
      integer :: nss,iss, istart

      logical :: mat_par=.true.!if true uses parallel routines

      CALL start_clock( 'pc2' )

      do iss= 1, nspin
         nss= nupdwn( iss )
!          write(6,*) "nupdwn", iss, nupdwn(iss), iupdwn(iss)
         if(nss>0) THEN
            
            istart= iupdwn( iss )

            if(lgam) then
               allocate(bectmp(nss,nss))
               bectmp(:,:)=0.d0
            else
               allocate(bectmp_c(nss,nss))
               bectmp_c(:,:)=CMPLX(0.d0,0.d0)
            endif
   ! 
            allocate(zbectmp(nss,nss))
            call zgemm('C','N',nss,nss,ngw,(1.d0,0.d0),a(:,istart), &
                        ngw,b(:,istart),ngw,(0.d0,0.d0),zbectmp,nss)

            if(lgam) then
               do j=1,nss
                   do i=1,nss
                     bectmp(i,j)=2.d0*DBLE(zbectmp(i,j))
                     if(ng0.eq.2) bectmp(i,j)=bectmp(i,j)-DBLE(a(1,j))*DBLE(b(1,i))
                   enddo
               enddo
               call mp_sum( bectmp(:,:), intra_image_comm)
            else
               do j=1,nss
                   do i=1,nss
                     bectmp_c(i,j)=zbectmp(i,j)
                   enddo
               enddo
               call mp_sum( bectmp_c(:,:), intra_image_comm)
            endif
            deallocate(zbectmp)
            if(nvb >= 0) then

               nl_max=0
               do is=1,nvb
                  nl_max=nl_max+nh(is)*na(is)
               enddo
               if(lgam) then
                   allocate (qqb_tmp(nl_max,nss))
                   allocate (qq_tmp(nl_max,nl_max))
                   qq_tmp(:,:)=0.d0
                   do is=1,nvb
                     do iv=1,nh(is)
                         do jv=1,nh(is)
                           do ia=1,na(is)
                               inl=ish(is)+(iv-1)*na(is)+ia
                               jnl=ish(is)+(jv-1)*na(is)+ia
                               qq_tmp(inl,jnl)=qq(iv,jv,is)
                           enddo
                         enddo
                     enddo
                   enddo
               else
                   allocate (qqb_tmp_c(nl_max,nss))
                   allocate (qq_tmp_c(nl_max,nl_max))
                   qq_tmp_c(:,:)=CMPLX(0.d0,0.d0)
                   do is=1,nvb
                     do iv=1,nh(is)
                         do jv=1,nh(is)
                           do ia=1,na(is)
                               inl=ish(is)+(iv-1)*na(is)+ia
                               jnl=ish(is)+(jv-1)*na(is)+ia
                               qq_tmp_c(inl,jnl)=CMPLX(qq(iv,jv,is),0.d0)
                           enddo
                         enddo
                     enddo
                   enddo
               endif
               !
               if(lgam) then
                   if( nhsa > 0 .and. .not. mat_par)  then
                     call dgemm('N','N',nl_max,nss,nl_max,1.d0,qq_tmp,nl_max,becb%rvec(:,istart),nhsa,0.d0,qqb_tmp,nl_max)
                     call dgemm('T','N',nss,nss,nl_max,1.d0,beca%rvec(:,istart),nhsa,qqb_tmp,nl_max,1.d0,bectmp,nss)
                   else if ( nhsa > 0 ) then
                     call para_dgemm ('N','N',nl_max,nss,nl_max,1.d0,qq_tmp,nl_max,&
                                       becb%rvec(:,istart),nhsa,0.d0,qqb_tmp,nl_max, intra_image_comm)
                     call para_dgemm ('T','N',nss,nss,nl_max,1.d0,beca%rvec(:,istart),nhsa, &
                                       qqb_tmp,nl_max,1.d0,bectmp,nss, intra_image_comm)
                   endif
                  deallocate(qq_tmp,qqb_tmp)
               else
                   if( nhsa > 0 .and. .not. mat_par)  then
                     call zgemm('N','N',nl_max,nss,nl_max,(1.d0,0.d0),qq_tmp_c,&
                                 nl_max,becb%cvec(:,istart),nhsa,(0.d0,0.d0), qqb_tmp_c,nl_max)
                     call zgemm('C','N',nss,nss,nl_max,(1.d0,0.d0),beca%cvec(:,istart),&
                                 nhsa,qqb_tmp_c,nl_max,(1.d0,0.d0),bectmp_c,nss)
                   else if ( nhsa > 0 ) then
                     call para_zgemm ('N','N',nl_max,nss,nl_max,(1.d0,0.d0),qq_tmp_c,nl_max,&
                                       becb%cvec(:,istart),nhsa,(0.d0,0.d0),qqb_tmp_c,nl_max, intra_image_comm)
                     call para_zgemm ('C','N',nss,nss,nl_max,(1.d0,0.d0),beca%cvec(:,istart),nhsa, &
                                       qqb_tmp_c,nl_max,(1.d0,0.d0),bectmp_c,nss, intra_image_comm)
                   endif
                  deallocate(qq_tmp_c,qqb_tmp_c)
               endif   
               !
            endif
            allocate(zbectmp(nss,nss))
            if(lgam) then
               do i=1,nss
                   do j=1,nss
                     zbectmp(i,j)=CMPLX(bectmp(i,j),0.d0)
                   enddo
               enddo
            else
               do i=1,nss
                   do j=1,nss
                     zbectmp(i,j)=bectmp_c(i,j)
                   enddo
               enddo
            endif
            call zgemm('N','N',ngw,nss,nss,(-1.d0,0.d0),a(:,istart), &
                        ngw,zbectmp,nss,(1.d0,0.d0),b(:,istart),ngw)
            deallocate(zbectmp)

            ! this computes the new bec
            if(lgam) then
               if ( nhsa > 0 ) then
                   call dgemm('N','N',nhsa,nss,nss,1.0d0,beca%rvec(:,istart), & 
                               nhsa,bectmp,nss,1.0d0,becb%rvec(:,istart),nhsa)
               endif
               deallocate(bectmp)
            else
               if ( nhsa > 0 ) then
                   call zgemm('N','N',nhsa,nss,nss,(1.0d0,0.d0),beca%cvec(:,istart), & 
                               nhsa,bectmp_c,nss,(1.0d0,0.d0),becb%cvec(:,istart),nhsa)
               endif
               deallocate(bectmp_c)
            endif
            !
         ENDIF
      enddo!on spin
      CALL stop_clock( 'pc2' )
      return
    end subroutine pc2_new
    
    subroutine pcdaga2_non_ortho(a,adual, as ,b, lgam )

! this function applies the operator Pc

!    this subroutine applies the Pc^dagerr operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i - S|a_j><a_j|b_i>

      use kinds
      use mp_global, only: intra_image_comm
      use cvan
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, ispin

      implicit none

      complex(dp) a(ngw,n), adual(ngw,n), b(ngw,n), as(ngw,n)
      logical :: lgam
      ! local variables
      integer i, j,ig
      complex(dp) sca
      complex(DP), allocatable:: scar(:)
      !
      call start_clock('pcdaga2')
      allocate(scar(n))
      do j=1,n
         do i=1,n
            sca=0.0d0
            if(ispin(i) == ispin(j)) then
               IF(lgam) THEN
		  if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
		  do  ig=1,ngw           !loop on g vectors
		      sca=sca+DBLE(CONJG(a(ig,j))*b(ig,i))
		  enddo
		  sca = sca*2.0d0  !2. for real weavefunctions
		  if (ng0.eq.2) sca = sca - DBLE(a(1,j))*DBLE(b(1,i))
               ELSE
		  do  ig=1,ngw           !loop on g vectors
		      sca=sca+CONJG(a(ig,j))*b(ig,i)
		  enddo
               ENDIF
            endif
            scar(i) = sca
         enddo
                   
         call mp_sum( scar, intra_image_comm )

         do i=1,n
            if(ispin(i) == ispin(j)) then
               sca = scar(i)
	       do ig=1,ngw
		  b(ig,i)=b(ig,i)-sca*as(ig,j)
	       enddo
		! this to prevent numerical errors
               IF(lgam) THEN
		  if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
               ENDIF
            endif
         enddo
      enddo
      deallocate(scar)
      call stop_clock('pcdaga2')
      return
      end subroutine pcdaga2_non_ortho


    subroutine pcdaga2(a,as ,b, lgam )

! this function applies the operator Pc

!    this subroutine applies the Pc^dagerr operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i - S|a_j><a_j|b_i>

      use kinds
      use mp_global, only: intra_image_comm
      use cvan
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, ispin

      implicit none

      complex(dp) a(ngw,n), b(ngw,n), as(ngw,n)
      logical :: lgam
      ! local variables
      integer i, j,ig
      complex(dp) sca
      complex(DP), allocatable:: scar(:)
      !
      call start_clock('pcdaga2')
      allocate(scar(n))
      do j=1,n
         do i=1,n
            sca=0.0d0
            if(ispin(i) == ispin(j)) then
               IF(lgam) THEN
		  if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
		  do  ig=1,ngw           !loop on g vectors
		      sca=sca+DBLE(CONJG(a(ig,j))*b(ig,i))
		  enddo
		  sca = sca*2.0d0  !2. for real weavefunctions
		  if (ng0.eq.2) sca = sca - DBLE(a(1,j))*DBLE(b(1,i))
               ELSE
		  do  ig=1,ngw           !loop on g vectors
		      sca=sca+CONJG(a(ig,j))*b(ig,i)
		  enddo
               ENDIF
            endif
            scar(i) = sca
         enddo
                   
         call mp_sum( scar, intra_image_comm )

         do i=1,n
            if(ispin(i) == ispin(j)) then
               sca = scar(i)
	       do ig=1,ngw
		  b(ig,i)=b(ig,i)-sca*as(ig,j)
	       enddo
		! this to prevent numerical errors
               IF(lgam) THEN
		  if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
               ENDIF
            endif
         enddo
      enddo
      deallocate(scar)
      call stop_clock('pcdaga2')
      return
      end subroutine pcdaga2

      subroutine pcdaga2_new(a,as ,b, n, ispin, lgam )

! this function applies the operator Pc

!    this subroutine applies the Pc^dagerr operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i - S|a_j><a_j|b_i>

      use kinds
      use mp_global, only: intra_image_comm
      use cvan
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum

      implicit none

      integer, intent(in) :: n, ispin(n)
      complex(dp) a(ngw,n), b(ngw,n), as(ngw,n)
      logical :: lgam
      ! local variables
      integer  i, j,ig
      complex(dp) sca
      complex(DP), allocatable:: scar(:)
      !
      call start_clock('pcdaga2')
      allocate(scar(n))
      do j=1,n
         do i=1,n
            sca=0.0d0
            if(ispin(i) == ispin(j)) then
               IF(lgam) THEN
                  if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
                  do  ig=1,ngw           !loop on g vectors
                      sca=sca+DBLE(CONJG(a(ig,j))*b(ig,i))
                  enddo
                  sca = sca*2.0d0  !2. for real weavefunctions
                  if (ng0.eq.2) sca = sca - DBLE(a(1,j))*DBLE(b(1,i))
               ELSE
                  do  ig=1,ngw           !loop on g vectors
                      sca=sca+CONJG(a(ig,j))*b(ig,i)
                  enddo
               ENDIF
            endif
            scar(i) = sca
         enddo
                   
         call mp_sum( scar, intra_image_comm )

         do i=1,n
            if(ispin(i) == ispin(j)) then
               sca = scar(i)
               do ig=1,ngw
                  b(ig,i)=b(ig,i)-sca*as(ig,j)
               enddo
                ! this to prevent numerical errors
               IF(lgam) THEN
                  if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
               ENDIF
            endif
         enddo
      enddo
      deallocate(scar)
      call stop_clock('pcdaga2')
      return
      end subroutine pcdaga2_new
      
    subroutine pcdaga3(a,as ,b, lgam )

 !For LOWDIN orthogonalization

!    this subroutine applies the Pc^dagerr operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i - S|a_j>(<a_j|b_i>+<a_i|b_j>)/2

      use kinds
      use mp_global, only: intra_image_comm
      use cvan
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, ispin

      implicit none

      complex(dp) a(ngw,n), b(ngw,n), as(ngw,n)
      logical :: lgam
      ! local variables
      integer i, j,ig
      complex(dp) sca
      complex(DP) :: bold(ngw,n)
      complex(DP), allocatable:: scar(:)
      !
      call start_clock('pcdaga2')
      allocate(scar(n))
      bold(:,:) = b(:,:)
      !
      do j=1,n
         do i=1,n
            sca=0.0d0
            if(ispin(i) == ispin(j)) then
               IF(lgam) THEN
		  if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
		  do  ig=1,ngw           !loop on g vectors
		      sca=sca+DBLE( CONJG(a(ig,j))*bold(ig,i)+a(ig,i)*CONJG(bold(ig,j)))
		  enddo
		  sca = sca*2.0d0  !2. for real weavefunctions
		  if (ng0.eq.2) sca = sca - DBLE(CONJG(a(1,j))*bold(1,i)+a(1,i)*CONJG(bold(1,j)))
               ELSE
		  do  ig=1,ngw           !loop on g vectors
		      sca=sca+CONJG(a(ig,j))*bold(ig,i)+a(ig,i)*CONJG(bold(ig,j))
		  enddo
               ENDIF
            endif
            scar(i) = sca*0.5d0
         enddo
                   
         call mp_sum( scar, intra_image_comm )

         do i=1,n
            if(ispin(i) == ispin(j)) then
               sca = scar(i)
	       do ig=1,ngw
		  b(ig,i)=b(ig,i)-sca*as(ig,j)
	       enddo
		! this to prevent numerical errors
               IF(lgam) THEN
		  if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
               ENDIF
            endif
         enddo
      enddo
      deallocate(scar)
      call stop_clock('pcdaga2')
      return
      end subroutine pcdaga3

!$$
!----------------------------------------------------------------------
     SUBROUTINE lowdin(a, lgam)
!----------------------------------------------------------------------

      use kinds
      use io_global, only: ionode
      use mp_global, only: intra_image_comm
      use gvecw, only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use cvan, only: nvb
      use electrons_base, only: n => nbsp, nspin, nupdwn, iupdwn

      implicit none

      complex(dp) a(ngw,n), aold(ngw,n)
      integer i, j,k,ig, isp,ndim,nbnd1,nbnd2
      real(dp) sqrt_seig(n)
      complex(DP) :: sca
      real(DP), allocatable :: seig(:)
      complex(DP), allocatable :: s(:,:), omat(:,:), sqrt_s(:,:)
      logical :: lgam, okvan
      !
      okvan=nvb>0
      aold(:,:)=a(:,:)

      do isp=1,nspin

        ndim=nupdwn(isp)

        allocate(s(ndim,ndim))
        allocate(omat(ndim,ndim))
        allocate(seig(ndim))
        allocate(sqrt_s(ndim,ndim))

        s(:,:)=CMPLX(0.d0,0.d0)

        do i=1,ndim
           !
           nbnd1=iupdwn(isp)-1+i
           !
           do j=1,i
              !
              nbnd2=iupdwn(isp)-1+j
              !
              sca=CMPLX(0.0d0,0.d0)
              !
              IF(lgam) THEN
                 !
                 if (ng0.eq.2) aold(1,nbnd1) = CMPLX(DBLE(a(1,nbnd1)),0.0d0)
                 !
                 do  ig=1,ngw           !loop on g vectors
                    !
                    sca=sca+DBLE(CONJG(a(ig,nbnd2))*a(ig,nbnd1))
                    !
                 enddo
                 !
                 sca = sca*2.0d0  !2. for real weavefunctions
                 if (ng0.eq.2) sca = sca - DBLE(CONJG(a(1,nbnd2))*a(1,nbnd1))
                 s(i,j) = CMPLX(DBLE(sca),0.d0)
                 s(j,i) = s(i,j)
                 !
              ELSE
                 !
                 do  ig=1,ngw           !loop on g vectors
                    !
                    sca=sca+CONJG(a(ig,nbnd2))*a(ig,nbnd1)
                    s(i,j) = sca
                    s(j,i) = CONJG(sca)
                    !
                 enddo
                 !
              ENDIF
              !
           enddo
           !
        enddo
        
        call mp_sum( s, intra_image_comm )
        call zdiag(ndim,ndim,s,seig,omat,1)

        do i=1,ndim
           !
           if(seig(i).lt.0.d0.and.ionode) write(*,*) 'seig is negative ',seig(:)
           !
        enddo

        sqrt_seig(:)=1.d0/DSQRT(seig(:))

        sqrt_s(:,:)=CMPLX(0.d0,0.d0)

        do i=1,ndim
           !
           do j=1,i
!             if(j.lt.i) then
!               sqrt_s(i,j)=sqrt_s(j,i)
!             else
              sca=0.d0
              do k=1,ndim
                 !
                 sca=sca+sqrt_seig(k) * omat(i,k)*CONJG(omat(j,k))
                 !
              enddo
              sqrt_s(i,j) = sca
              sqrt_s(j,i) = CONJG(sca)
              !
          enddo
          !
        enddo

        do i=1,ndim
           !
           nbnd1=iupdwn(isp)-1+i
           a(:,nbnd1) = CMPLX(0.d0,0.d0)
           !
           do j=1,ndim
              !
              nbnd2=iupdwn(isp)-1+j
              a(:,nbnd1) = a(:,nbnd1) + sqrt_s(i,j) * aold(:,nbnd2)
              !
           enddo
           !
        enddo

        deallocate(s)
        deallocate(omat)
        deallocate(seig)
        deallocate(sqrt_s)

      enddo

     END SUBROUTINE lowdin
!$$

!----------------------------------------------------------------------
     SUBROUTINE lowdin_uspp(a, beca, lgam)
!----------------------------------------------------------------------

      use kinds
      use io_global, only: ionode
      use gvecw, only: ngw
      use mp, only: mp_sum
      use cvan, only: nvb
      use twin_types
      use electrons_base, only: n => nbsp, nbspx, nspin, nupdwn, iupdwn

      implicit none

      complex(dp) a(ngw,n), aold(ngw,n)
      integer i, j,k,isp,ndim,nbnd1,nbnd2
      real(dp) sqrt_seig(n)
      complex(DP) :: sca
      type(twin_matrix) :: beca
      real(DP), allocatable :: seig(:)
      complex(DP), allocatable :: s(:,:), omat(:,:), sqrt_s(:,:)
      logical :: lgam, okvan
      !
      okvan=nvb>0
      aold(:,:)=a(:,:)

      do isp=1,nspin

        ndim=nupdwn(isp)

         if (ndim>0) then
            allocate(s(ndim,ndim))
            allocate(omat(ndim,ndim))
            allocate(seig(ndim))
            allocate(sqrt_s(ndim,ndim))

            s(:,:)=CMPLX(0.d0,0.d0)

            do i=1,ndim
               !
               nbnd1=iupdwn(isp)-1+i
               !
               do j=1,i
                  !
                  nbnd2=iupdwn(isp)-1+j
                  !
                  call dotcsv( s(j,i), nbspx, n, a, beca, a, beca, ngw, iupdwn(isp)+j-1, iupdwn(isp)+i-1, lgam)
                  s(i,j)=CONJG(s(j,i))
                  !
               enddo
               !
            enddo
            do i=1,ndim
               write(111,*) s(i,:)
               write(111,*) i
            enddo
               write(111,*) "END"
            !         call mp_sum( s, intra_image_comm )
            call zdiag(ndim,ndim,s,seig,omat,1)

            do i=1,ndim
               !
               if(seig(i).lt.0.d0.and.ionode) write(*,*) 'seig is negative ',seig(:)
               !
            enddo

            sqrt_seig(:)=1.d0/DSQRT(seig(:))

            sqrt_s(:,:)=CMPLX(0.d0,0.d0)

            do i=1,ndim
               !
               do j=1,i
                  !
                  sca=0.d0
                  do k=1,ndim
                     !
                     sca=sca+sqrt_seig(k) * omat(i,k)*CONJG(omat(j,k))
                     !
                  enddo
                  sqrt_s(i,j) = sca
                  sqrt_s(j,i) = CONJG(sca)
                  !
               enddo
               !
            enddo

            do i=1,ndim
               !
               nbnd1=iupdwn(isp)-1+i
               a(:,nbnd1) = CMPLX(0.d0,0.d0)
               !
               do j=1,ndim
                  !
                  nbnd2=iupdwn(isp)-1+j
                  a(:,nbnd1) = a(:,nbnd1) + sqrt_s(j,i) * aold(:,nbnd2)
                  !
               enddo
               !
            enddo

            deallocate(s)
            deallocate(omat)
            deallocate(seig)
            deallocate(sqrt_s)
            !
         endif

      enddo

     END SUBROUTINE lowdin_uspp
!$$

!$$
    subroutine pc3us(a, beca, b, becb, lgam)

! this function applies the modified Pc operator which is
! equivalent to Lowdin orthonormalization of the revised wavefunctions.
! this subroutine works for pseudopotentials. 
               
! this function applies the operator Pc
            
!    this subroutine applies the Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i-a_j>(<a_j|S|b_i>+<b_j|S|a_i>)/2
    
      use kinds, only: dp 
      use ions_base, only: na
      use mp_global, only: intra_image_comm
      use cvan 
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, nupdwn, iupdwn, nspin
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb
      use uspp, only :qq
      use parallel_toolkit, only : rep_matmul_drv
      use twin_types      
                           
      implicit none        
                           
      complex(kind=DP) a(ngw,n), b(ngw,n)
                     
      type(twin_matrix) ::   beca,becb!(nhsa,n) !modified:giovanni
      logical :: lgam
! local variables
      integer is, iv, jv, ia, inl, jnl, i, j
      real(DP), allocatable :: bectmp(:,:)
      complex(DP), allocatable :: bectmp_c(:,:)
      real(DP), allocatable :: qq_tmp(:,:), qqb_tmp(:,:)
      complex(DP), allocatable :: qqb_tmp_c(:,:), qq_tmp_c(:,:)
      complex(DP), allocatable :: zbectmp(:,:)
      integer :: nl_max
      integer :: nss,iss, istart

      logical :: mat_par=.true.!if true uses parallel routines

      CALL start_clock( 'pc3us' )

      do iss= 1, nspin
         !
         nss= nupdwn( iss )
         istart= iupdwn( iss )

         if(lgam) then
            !
            allocate(bectmp(nss,nss))
            bectmp(:,:)=0.d0
            !
         else
            !
            allocate(bectmp_c(nss,nss))
            bectmp_c(:,:)=CMPLX(0.d0,0.d0)
            !
         endif
         !
         allocate(zbectmp(nss,nss))
         !
         call zgemm('C','N',nss,nss,ngw,(1.d0,0.d0),a(:,istart),ngw,b(:,istart),ngw,(0.d0,0.d0),zbectmp,nss)

         if(lgam) then
            !
            do j=1,nss
               !
               do i=1,nss
                  !
                  bectmp(i,j)=2.d0*DBLE(zbectmp(i,j))
                  if(ng0.eq.2) bectmp(i,j)=bectmp(i,j)-DBLE(CONJG(a(1,j))*(b(1,i)))
                  !
               enddo
               !
            enddo
            !
            call mp_sum( bectmp(:,:), intra_image_comm)
            !
         else
            !
            do j=1,nss
               !
               do i=1,nss
                  !
                  bectmp_c(i,j)=zbectmp(i,j)
                  !
               enddo
               !
            enddo
            !
            call mp_sum( bectmp_c(:,:), intra_image_comm)
            !
         endif
         !
         deallocate(zbectmp)
         !
         if(nvb >= 0) then
            !
            nl_max=0
            !
            do is=1,nvb
               !
               nl_max=nl_max+nh(is)*na(is)
               !
            enddo
            !
            if(lgam) then
               !
               allocate (qqb_tmp(nl_max,nss))
               allocate (qq_tmp(nl_max,nl_max))
               qq_tmp(:,:)=0.d0
               !
               do is=1,nvb
                  !
                  do iv=1,nh(is)
                     !
                     do jv=1,nh(is)
                        !
                        do ia=1,na(is)
                           !
                           inl=ish(is)+(iv-1)*na(is)+ia
                           jnl=ish(is)+(jv-1)*na(is)+ia
                           qq_tmp(inl,jnl)=qq(iv,jv,is)
                           !
                        enddo
                        !
                     enddo
                     !
                  enddo
                  !
               enddo
               !
            else
               !
               allocate (qqb_tmp_c(nl_max,nss))
               allocate (qq_tmp_c(nl_max,nl_max))
               qq_tmp_c(:,:)=CMPLX(0.d0,0.d0)
               !
               do is=1,nvb
                  !
                  do iv=1,nh(is)
                     !
                     do jv=1,nh(is)
                        !
                        do ia=1,na(is)
                           !
                           inl=ish(is)+(iv-1)*na(is)+ia
                           jnl=ish(is)+(jv-1)*na(is)+ia
                           qq_tmp_c(inl,jnl)=CMPLX(qq(iv,jv,is),0.d0)
                           !
                        enddo
                        !
                     enddo
                     !
                  enddo
                  !
               enddo
               !
            endif
            !
            if(lgam) then
               !
               if( nhsa > 0 .and. .not. mat_par)  then
                  !
                  call dgemm('N','N',nl_max,nss,nl_max,1.d0,qq_tmp,nl_max,becb%rvec(:,istart),nhsa,0.d0,qqb_tmp,nl_max)
                  call dgemm('T','N',nss,nss,nl_max,1.d0,beca%rvec(:,istart),nhsa,qqb_tmp,nl_max,1.d0,bectmp,nss)
                  !
               else if ( nhsa > 0 ) then
                  !
                  call para_dgemm ('N','N',nl_max,nss,nl_max,1.d0,qq_tmp,nl_max,&
                  becb%rvec(:,istart),nhsa,0.d0,qqb_tmp,nl_max, intra_image_comm)
                  call para_dgemm ('T','N',nss,nss,nl_max,1.d0,beca%rvec(:,istart),nhsa, &
                  qqb_tmp,nl_max,1.d0,bectmp,nss, intra_image_comm)
                  !
               endif
               !
               deallocate(qq_tmp,qqb_tmp)
               !
            else
               !
               if( nhsa > 0 .and. .not. mat_par)  then
                  !
                  call zgemm('N','N',nl_max,nss,nl_max,(1.d0,0.d0),qq_tmp_c, & 
                              nl_max,becb%cvec(:,istart),nhsa,(0.d0,0.d0), qqb_tmp_c,nl_max)
                  call zgemm('C','N',nss,nss,nl_max,(1.d0,0.d0),beca%cvec(:,istart), &
                              nhsa,qqb_tmp_c,nl_max,(1.d0,0.d0),bectmp_c,nss)
                  !
               else if ( nhsa > 0 ) then
                  !
                  call para_zgemm ('N','N',nl_max,nss,nl_max,(1.d0,0.d0),qq_tmp_c,nl_max,&
                  becb%cvec(:,istart),nhsa,(0.d0,0.d0),qqb_tmp_c,nl_max, intra_image_comm)
                  call para_zgemm ('C','N',nss,nss,nl_max,(1.d0,0.d0),beca%cvec(:,istart),nhsa, &
                  qqb_tmp_c,nl_max,(1.d0,0.d0),bectmp_c,nss, intra_image_comm)
                  !
               endif
               !
               deallocate(qq_tmp_c,qqb_tmp_c)
               !
            endif   
            !
         endif

         allocate(zbectmp(nss,nss))

         if(lgam) then
            !
            do i=1,nss
               !
               do j=1,nss
                  !
                  zbectmp(i,j)=0.5d0*(CMPLX(bectmp(i,j),0.d0) + CMPLX(bectmp(j,i),0.d0))
                  !
               enddo
               !
            enddo
            !
            bectmp(:,:) = DBLE(zbectmp(:,:))
            !
         else
            !
            do i=1,nss
               !
               do j=1,nss
                  !
                  zbectmp(i,j)=0.5d0*(bectmp_c(i,j)+CONJG(bectmp_c(j,i)))
                  !
               enddo
               !
            enddo
            !
            bectmp_c(:,:) = CONJG(zbectmp(:,:))
            !
         endif

         call zgemm('N','N',ngw,nss,nss,(-1.d0,0.d0),a(:,istart),ngw,zbectmp,nss,(1.d0,0.d0),b(:,istart),ngw)
         deallocate(zbectmp)

         ! this computes the new bec
         if(lgam) then
            !
            if ( nhsa > 0 ) then
               !
               call dgemm('N','N',nhsa,nss,nss,1.0d0,beca%rvec(:,istart),nhsa,bectmp,nss,1.0d0,becb%rvec(:,istart),nhsa)
               !
            endif
            !
            deallocate(bectmp)
            !
         else
            !
            if ( nhsa > 0 ) then
               !
               call zgemm('N','N',nhsa,nss,nss,(1.0d0,0.d0),beca%cvec(:,istart), & 
                           nhsa,bectmp_c,nss,(1.0d0,0.d0),becb%cvec(:,istart),nhsa)
               !
            endif
            !
            deallocate(bectmp_c)
            !
         endif
         !
      enddo!on spin
      !
      CALL stop_clock( 'pc3us' )
      !
      return
      !
      end subroutine pc3us

    subroutine pc3nc(a,b, lgam)

! this function applies the modified Pc operator which is
! equivalent to Lowdin orthonormalization of the revised wavefunctions.
! currently implemented only for norm-conserving pseudopotentials. 

!    this subroutine applies the modified Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i - |a_j>(<a_j|b_i>+<b_j|a_i>)/2

      use kinds
      use mp_global, only: intra_image_comm
      use gvecw, only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, ispin

      implicit none

      complex(dp) a(ngw,n), b(ngw,n)
      logical :: lgam
      ! local variables
      complex(DP) :: bold(ngw,n)
      integer i, j,ig
!       real(dp) sca
      complex(DP) :: sca_c
      complex(DP), allocatable:: scar_c(:)
      !
      call start_clock('pc3')

      allocate(scar_c(n))

      bold(:,:)=b(:,:)

      do j=1,n
         !
         do i=1,n
            !
            sca_c=CMPLX(0.0d0,0.d0)
            !
            if(ispin(i) == ispin(j)) then
               !
               if(lgam) then
                  !
                  if (ng0.eq.2) bold(1,i) = CMPLX(DBLE(bold(1,i)),0.0d0)
                  !
               endif
               !
               do  ig=1,ngw           !loop on g vectors
                  !
                  sca_c=sca_c+CONJG(a(ig,j))*bold(ig,i) !uncomment this for lowdin ortho
                  sca_c=sca_c+(a(ig,i))*CONJG(bold(ig,j)) !remove the 2.d0 for lowdin ortho
                  !
               enddo
               !sca = sca*2.0d0  !2. for real weavefunctions
               !$$ not necessary: sca = sca*2.0d0  !2. for real weavefunctions
               if(lgam) then
                  !
                  if (ng0.eq.2) then
                     !
                     sca_c = CMPLX(DBLE(sca_c),0.d0) - &
                             CMPLX(0.5d0*DBLE(CONJG(a(1,j))*(bold(1,i))+(a(1,i))*CONJG(bold(1,j))),0.d0) !use this one for lowdin ortho
                   !sca_c = CMPLX(DBLE(sca_c),0.d0) - CMPLX(DBLE((a(1,i))*CONJG(bold(1,j))),0.d0) !comment this one for lowdin ortho
                  else
                     !
                     sca_c = CMPLX(DBLE(sca_c), 0.d0)
                     !
                  endif
                  !
               else
                  !
                  sca_c=0.5d0*sca_c
                  !
               endif
               !
	            scar_c(i) = sca_c
               !
            endif
            !
         enddo

         call mp_sum( scar_c, intra_image_comm )

         do i=1,n
            !
            if(ispin(i) == ispin(j)) then
               !
               sca_c = scar_c(i)
               !
               do ig=1,ngw
                  !
                  b(ig,i)=b(ig,i)-sca_c*a(ig,j)
                  !
               enddo
               ! this to prevent numerical errors
               if(lgam) then
                  ! 
                  if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
                  !
               endif
               !
            endif
            !
         enddo
         !
      enddo
      !
      deallocate(scar_c)
      call stop_clock('pc3')
      return
      end subroutine pc3nc

    subroutine pc3nc_new(a,b,n,ispin, lgam)

! this function applies the modified Pc operator which is
! equivalent to Lowdin orthonormalization of the revised wavefunctions.
! currently implemented only for norm-conserving pseudopotentials. 

!    this subroutine applies the modified Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i - |a_j>(<a_j|b_i>+<b_j|a_i>)/2

      use kinds
      use mp_global, only: intra_image_comm
      use gvecw, only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum

      implicit none

      integer, intent(in) :: n, ispin(n)
      complex(dp) :: a(ngw,n), b(ngw,n)
      logical :: lgam
      ! local variables
      complex(DP) :: bold(ngw,n)
      integer i, j,ig
!       real(dp) sca
      complex(DP) :: sca_c
      complex(DP), allocatable:: scar_c(:)
      !
      call start_clock('pc3')

      allocate(scar_c(n))

      bold(:,:)=b(:,:)

      do j=1,n
         !
         do i=1,n
            !
            sca_c=CMPLX(0.0d0,0.d0)
            !
            if(ispin(i) == ispin(j)) then
                !
                if(lgam) then
                   !
                   if (ng0.eq.2) bold(1,i) = CMPLX(DBLE(bold(1,i)),0.0d0)
                   !
                endif
                !
                do  ig=1,ngw
                    !
                    sca_c=sca_c+CONJG(a(ig,j))*bold(ig,i) !uncomment this for lowdin ortho
                    sca_c=sca_c+(a(ig,i))*CONJG(bold(ig,j)) !remove the 2.d0 for lowdin ortho
                    !
                enddo
                !sca = sca*2.0d0  !2. for real weavefunctions
                !$$ not necessary: sca = sca*2.0d0  !2. for real weavefunctions
                if(lgam) then
                   !
                   if (ng0.eq.2) then
                      !
                      sca_c = CMPLX(DBLE(sca_c),0.d0) - &
                              CMPLX(0.5d0*DBLE(CONJG(a(1,j))*(bold(1,i))+(a(1,i))*CONJG(bold(1,j))),0.d0) !use this one for lowdin ortho
                   !sca_c = CMPLX(DBLE(sca_c),0.d0) - CMPLX(DBLE((a(1,i))*CONJG(bold(1,j))),0.d0) !comment this one for lowdin ortho
                   else
                      !
                      sca_c = CMPLX(DBLE(sca_c), 0.d0)
                      !
                   endif
                   !
                else
                   !
                   sca_c=0.5d0*sca_c
                   !
                endif
                !
                scar_c(i) = sca_c
                !
            endif
            !
         enddo
         !
         call mp_sum( scar_c, intra_image_comm )
         !
         do i=1,n
            !
            if(ispin(i) == ispin(j)) then
               !
               sca_c = scar_c(i)
               !
               do ig=1,ngw
                  !
                  b(ig,i)=b(ig,i)-sca_c*a(ig,j)
                  !
               enddo
               ! this to prevent numerical errors
               if(lgam) then
                  ! 
                  if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
                  !
               endif
               !
            endif
            !
         enddo
         !
      enddo
      !
      deallocate(scar_c)
      call stop_clock('pc3')
      return
      end subroutine pc3nc_new

    subroutine pc3nc_both(a, b, n_emp, c0, n, ispin_emp, ispin, lgam)

! this function applies the modified Pc operator which is
! equivalent to Lowdin orthonormalization of the revised wavefunctions.
! currently implemented only for norm-conserving pseudopotentials. 

!    this subroutine applies the modified Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i - |a_j>(<a_j|b_i>+<b_j|a_i>)/2

      use kinds
      use mp_global, only: intra_image_comm
      use gvecw, only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum

      implicit none

      integer, intent(in) :: n_emp, ispin_emp(n_emp), n, ispin(n)
      complex(dp) :: a(ngw,n_emp), b(ngw,n_emp), c0(ngw,n)
      logical :: lgam
      ! local variables
      complex(DP) :: bold(ngw,n_emp)
      integer i, j, ig, ispin_tot(n+n_emp)
      complex(DP) :: sca_c
      complex(DP), allocatable:: scar_c(:)
      !
      call start_clock('pc3')

      allocate(scar_c(n_emp))
      !
      ispin_tot(1:n_emp) = ispin_emp(1:n_emp)
      ispin_tot(n_emp+1:n_emp+n) = ispin(n)
      !
      bold(:,:)=b(:,:)
      !
      do j=1,n_emp+n
         !
         do i=1,n_emp
            !
            sca_c=CMPLX(0.0d0,0.d0)
            !
            if(ispin_tot(i) == ispin_tot(j)) then
               !
               if(lgam) then
                  !
                  if (ng0.eq.2) bold(1,i) = CMPLX(DBLE(bold(1,i)),0.0d0)
                  !
               endif
               !
               IF(j<=n_emp) THEN
                  !
                  do  ig=1,ngw
                     !
                     sca_c=sca_c+CONJG(a(ig,j))*bold(ig,i) !uncomment this for lowdin ortho
                     sca_c=sca_c+(a(ig,i))*CONJG(bold(ig,j)) !remove the 2.d0 for lowdin ortho
                     !
                  enddo
                  !
               ELSE
                  !
                  do  ig=1,ngw
                     !
                     sca_c=sca_c+CONJG(c0(ig,j-n_emp))*bold(ig,i) !uncomment this for lowdin ortho
!                      sca_c=sca_c+(a(ig,i))*CONJG(bold(ig,j)) !remove the 2.d0 for lowdin ortho
                     !
                  enddo
                  !
               ENDIF
                !sca = sca*2.0d0  !2. for real weavefunctions
                !$$ not necessary: sca = sca*2.0d0  !2. for real weavefunctions
               IF(lgam) then
                  !
                  IF (ng0.eq.2) then
                     !
                     IF(j<=n_emp) THEN
                        !
                        sca_c = CMPLX(DBLE(sca_c),0.d0) - &
                                CMPLX(0.5d0*DBLE(CONJG(a(1,j))*(bold(1,i))+(a(1,i))*CONJG(bold(1,j))),0.d0) !use this one for lowdin ortho
                        !
                     ELSE
                        !
                        sca_c = CMPLX(DBLE(sca_c),0.d0) - & 
                                CMPLX(0.5d0*DBLE(CONJG(c0(1,j-n_emp))*(bold(1,i))),0.d0) !use this one for lowdin ortho
                        !
                     ENDIF
                     !
                   !sca_c = CMPLX(DBLE(sca_c),0.d0) - CMPLX(DBLE((a(1,i))*CONJG(bold(1,j))),0.d0) !comment this one for lowdin ortho
                  ELSE
                     !
                     sca_c = CMPLX(DBLE(sca_c), 0.d0)
                     !
                  ENDIF
                 !
               ELSE
                  !
                  sca_c=0.5d0*sca_c
                  !
               ENDIF
               !
               scar_c(i) = sca_c
               !
            ENDIF
            !
         ENDDO

         call mp_sum( scar_c, intra_image_comm )

         do i=1,n_emp
            !
            if(ispin_tot(i) == ispin_tot(j)) then
            
               sca_c = scar_c(i)
               !
               IF(j<=n_emp) THEN
                  !
                  do ig=1,ngw
                     !
                     b(ig,i)=b(ig,i)-sca_c*a(ig,j)
                     !
                  enddo
                  !
               ELSE
                  !
                  do ig=1,ngw
                     !
                     b(ig,i)=b(ig,i)-sca_c*c0(ig,j-n_emp)
                     !
                  enddo
                  !
               ENDIF
               ! this to prevent numerical errors
               IF(lgam) THEN
                  !
                  if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
                  !
               ENDIF
               !
            ENDIF
            !
         ENDDO
         !
      ENDDO
      !
      deallocate(scar_c)
      !
      call stop_clock('pc3')
      !
      return
      !
    end subroutine pc3nc_both
      
    subroutine pc4nc(a,b, lgam)

! this function applies the modified Pc operator which is
! equivalent to Lowdin orthonormalization of the revised wavefunctions.
! currently implemented only for norm-conserving pseudopotentials. 

!    this subroutine applies the modified Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i - |a_j>(<b_j|a_i>)

      use kinds
      use mp_global, only: intra_image_comm
      use gvecw, only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, ispin

      implicit none

      complex(dp) a(ngw,n), b(ngw,n)
      logical :: lgam
      ! local variables
      complex(DP) :: bold(ngw,n)
      integer i, j,ig
!       real(dp) sca
      complex(DP) :: sca_c
      complex(DP), allocatable:: scar_c(:)
      !
      call start_clock('pc3')

      allocate(scar_c(n))

      bold(:,:)=b(:,:)

      do j=1,n
         do i=1,n
	    sca_c=CMPLX(0.0d0,0.d0)
	    if(ispin(i) == ispin(j)) then
		if(lgam) then
                 if (ng0.eq.2) bold(1,i) = CMPLX(DBLE(bold(1,i)),0.0d0)
               endif
		do  ig=1,ngw           !loop on g vectors
		    !sca_c=sca_c+CONJG(a(ig,j))*bold(ig,i) !uncomment this for lowdin ortho
		    sca_c=sca_c+2.d0*(a(ig,i))*CONJG(bold(ig,j)) !remove the 2.d0 for lowdin ortho
		enddo
		!sca = sca*2.0d0  !2. for real weavefunctions
		!$$ not necessary: sca = sca*2.0d0  !2. for real weavefunctions
               if(lgam) then
		if (ng0.eq.2) then
                   !sca_c = CMPLX(DBLE(sca_c),0.d0) - CMPLX(0.5d0*DBLE(CONJG(a(1,j))*(bold(1,i))+(a(1,i))*CONJG(bold(1,j))),0.d0) !use this one for lowdin ortho
                   sca_c = CMPLX(DBLE(sca_c),0.d0) - CMPLX(DBLE((a(1,i))*CONJG(bold(1,j))),0.d0) !comment this one for lowdin ortho
                 else
                   sca_c = CMPLX(DBLE(sca_c), 0.d0)
                 endif
               else
                 sca_c=0.5d0*sca_c
               endif
	      scar_c(i) = sca_c
            endif
         enddo

         call mp_sum( scar_c, intra_image_comm )

         do i=1,n
            if(ispin(i) == ispin(j)) then
               sca_c = scar_c(i)
               do ig=1,ngw
                  b(ig,i)=b(ig,i)-sca_c*a(ig,j)
               enddo
               ! this to prevent numerical errors
               if(lgam) then 
                if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
               endif
            endif
         enddo
      enddo
      deallocate(scar_c)
      call stop_clock('pc3')
      return
      end subroutine pc4nc


      subroutine pc3nc_non_ortho(a,adual, b, lgam)

! this function applies the modified Pc operator which is
! equivalent to Lowdin orthonormalization of the revised wavefunctions.
! currently implemented only for norm-conserving pseudopotentials. 

!    this subroutine applies the modified Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i - |adual_j>(<a_j|b_i>+<b_j|a_i>)/2

      use kinds
      use mp_global, only: intra_image_comm
      use gvecw, only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, ispin

      implicit none

      complex(dp) a(ngw,n), adual(ngw,n), b(ngw,n)
      logical :: lgam
      ! local variables
      complex(DP) :: bold(ngw,n)
      integer i, j,ig
!       real(dp) sca
      complex(DP) :: sca_c
      complex(DP), allocatable:: scar_c(:)
      !
      call start_clock('pc3')

      allocate(scar_c(n))

      bold(:,:)=b(:,:)

      do j=1,n
         do i=1,n
            sca_c=CMPLX(0.0d0,0.d0)
            if(ispin(i) == ispin(j)) then
                if(lgam) then
                 if (ng0.eq.2) bold(1,i) = CMPLX(DBLE(bold(1,i)),0.0d0)
               endif
                do  ig=1,ngw           !loop on g vectors
                    sca_c=sca_c+CONJG(a(ig,j))*bold(ig,i)
                    sca_c=sca_c+(a(ig,i))*CONJG(bold(ig,j))
                enddo
                !sca = sca*2.0d0  !2. for real weavefunctions
                !$$ not necessary: sca = sca*2.0d0  !2. for real weavefunctions
               if(lgam) then
                if (ng0.eq.2) then
                   sca_c = CMPLX(DBLE(sca_c),0.d0) - & 
                           CMPLX(0.5d0*DBLE(CONJG(a(1,j))*(bold(1,i))+(a(1,i))*CONJG(bold(1,j))),0.d0)
                 else
                   sca_c = CMPLX(DBLE(sca_c), 0.d0)
                 endif
               else
                 sca_c=0.5d0*sca_c
               endif
              scar_c(i) = sca_c
            endif
         enddo

         call mp_sum( scar_c, intra_image_comm )

         do i=1,n
            if(ispin(i) == ispin(j)) then
               sca_c = scar_c(i)
               do ig=1,ngw
                  b(ig,i)=b(ig,i)-sca_c*adual(ig,j)
               enddo
               ! this to prevent numerical errors
               if(lgam) then 
                if (ng0.eq.2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
               endif
            endif
         enddo
      enddo
      deallocate(scar_c)
      call stop_clock('pc3')
      return
      end subroutine pc3nc_non_ortho
      
      subroutine pc4nc_non_ortho(c, a, adual, b, lgam)

! this function applies the modified Pc operator which is
! equivalent to Lowdin orthonormalization of the revised wavefunctions.
! currently implemented only for norm-conserving pseudopotentials. 

!    this subroutine applies the modified Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:c_i = - |adual_i>(<a_i|b_i>)

      use kinds
      use mp_global, only: intra_image_comm
      use gvecw, only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, ispin

      implicit none

      complex(dp) a(ngw,n), adual(ngw,n), b(ngw,n), c(ngw,n)
      logical :: lgam
      ! local variables
      complex(DP) :: bold(ngw,n)
      integer i, j,ig
!       real(dp) sca
      complex(DP) :: sca_c
      complex(DP), allocatable:: scar_c(:)
      !
      call start_clock('pc4')

      allocate(scar_c(n))

      bold(:,:)=b(:,:)

!       do j=1,n
         do i=1,n
            j=i
            sca_c=CMPLX(0.0d0,0.d0)
            if(ispin(i) == ispin(j)) then
                if(lgam) then
                 if (ng0.eq.2) bold(1,i) = CMPLX(DBLE(bold(1,i)),0.0d0)
               endif
                do  ig=1,ngw           !loop on g vectors
                    sca_c=sca_c+CONJG(a(ig,j))*bold(ig,i)
                    sca_c=sca_c+(a(ig,i))*CONJG(bold(ig,j))
                enddo
                !sca = sca*2.0d0  !2. for real weavefunctions
                !$$ not necessary: sca = sca*2.0d0  !2. for real weavefunctions
               if(lgam) then
                if (ng0.eq.2) then
                   sca_c = CMPLX(DBLE(sca_c),0.d0) - & 
                           CMPLX(0.5d0*DBLE(CONJG(a(1,j))*(bold(1,i))+(a(1,i))*CONJG(bold(1,j))),0.d0)
                 else
                   sca_c = CMPLX(DBLE(sca_c), 0.d0)
                 endif
               else
                 sca_c=0.5d0*sca_c
               endif
              scar_c(i) = sca_c
            endif
         enddo

         call mp_sum( scar_c, intra_image_comm )

         do i=1,n
            j=i
            if(ispin(i) == ispin(j)) then
               sca_c = scar_c(i)
               do ig=1,ngw
                  c(ig,i)=-sca_c*adual(ig,j)
               enddo
               ! this to prevent numerical errors
               if(lgam) then 
                if (ng0.eq.2) c(1,i) = CMPLX(DBLE(c(1,i)),0.0d0)
               endif
            endif
         enddo
!       enddo
      deallocate(scar_c)
      call stop_clock('pc4')
      return
      end subroutine pc4nc_non_ortho
      
     subroutine set_x_minus1_real(betae,m_minus1,ema0bg,use_ema)
     !
     ! this function calculates the factors for the inverse of the US K  matrix
     ! it takes care of the preconditioning
     ! see paper by Hasnip and Pickard, Computer Physics Communications 174 (2006) 24–29
     !

      use kinds, only: dp
      use ions_base, only: na
      use io_global, only: stdout
      use mp_global, only: intra_image_comm
      use cvan
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum, mp_bcast
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb,qq,nhsavb=>nkbus
      use io_global, ONLY: ionode, ionode_id

      implicit none

      complex(DP) :: betae(ngw,nhsa)
      real(DP)    :: m_minus1(nhsavb,nhsavb)
      real(DP)    :: ema0bg(ngw)
      logical     :: use_ema


! local variables
      real(DP),allocatable :: q_matrix(:,:), c_matrix(:,:)
      integer is, iv, jv, ia, inl, jnl, i, ig, js, ja
      real(DP) sca
      integer info, lwork
      integer, allocatable :: ipiv(:)
      real(dp),allocatable :: work(:)

      call start_clock('set_x_minus1')
      allocate(ipiv(nhsavb))
      allocate(work(nhsavb))

      lwork=nhsavb

      allocate(q_matrix(nhsavb,nhsavb),c_matrix(nhsavb,nhsavb))
!construct q matrix
      q_matrix(:,:) = 0.d0

      do is=1,nvb
         do iv=1,nh(is)
            do jv=1,nh(is)
               do ia=1,na(is)
                    inl=ish(is)+(iv-1)*na(is)+ia
                    jnl=ish(is)+(jv-1)*na(is)+ia
                    q_matrix(inl,jnl)= qq(iv,jv,is)
               enddo
            enddo
         enddo
      enddo

!construct C matrix
! m_minus1 used to be C matrix
      m_minus1(:,:) = 0.d0
      do is=1,nvb
         do ia=1,na(is)
            do iv=1,nh(is)
               do js=1,nvb
                  do ja=1,na(js)
                     do jv=1,nh(js)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(js)+(jv-1)*na(js)+ja
                        sca=0.d0
                        if (use_ema) then
                              ! k_minus case
                           do  ig=1,ngw           !loop on g vectors
                              sca=sca+ema0bg(ig)*DBLE(CONJG(betae(ig,inl))*betae(ig,jnl))
                           enddo
                           sca = sca*2.0d0  !2. for real weavefunctions
                           if (ng0.eq.2) sca = sca - ema0bg(1)*DBLE(CONJG(betae(1,inl))*betae(1,jnl))
                           else
                              ! s_minus case
                           do  ig=1,ngw           !loop on g vectors
                              sca=sca+DBLE(CONJG(betae(ig,inl))*betae(ig,jnl))
                           enddo
                           sca = sca*2.0d0  !2. for real weavefunctions
                           if (ng0.eq.2) sca = sca - DBLE(CONJG(betae(1,inl))*betae(1,jnl))
                        endif
                        m_minus1(inl,jnl)=sca
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      call mp_sum( m_minus1, intra_image_comm )

!calculate -(1+QC)**(-1) * Q
      CALL DGEMM('N','N',nhsavb,nhsavb,nhsavb,1.0d0,q_matrix,nhsavb,m_minus1,nhsavb,0.0d0,c_matrix,nhsavb)

      do i=1,nhsavb
         c_matrix(i,i)=c_matrix(i,i)+1.d0
      enddo

      if(ionode) then
        call dgetrf(nhsavb,nhsavb,c_matrix,nhsavb,ipiv,info)
        if(info .ne. 0) write(stdout,*) 'set_k_minus1 Problem with dgetrf :', info
        call dgetri(nhsavb,c_matrix,nhsavb,ipiv,work,lwork,info)
        if(info .ne. 0) write(stdout,*) 'set_k_minus1 Problem with dgetri :', info
      endif
      call mp_bcast( c_matrix, ionode_id, intra_image_comm )


      CALL DGEMM('N','N',nhsavb,nhsavb,nhsavb,-1.0d0,c_matrix,nhsavb,q_matrix,nhsavb,0.0d0,m_minus1,nhsavb)

      deallocate(q_matrix,c_matrix)
      deallocate(ipiv,work)
      call stop_clock('set_x_minus1')
      return
      
    end subroutine set_x_minus1_real

     subroutine set_x_minus1_twin(betae,m_minus1,ema0bg,use_ema)

     !
     ! this function calculates the factors for the inverse of the US K  matrix
     ! it takes care of the preconditioning
     ! see paper by Hasnip and Pickard, Computer Physics Communications 174 (2006) 24–29
     !
     ! this subroutine stores in m_mins1 the matrix R of the above paper
     !

      use kinds, only: dp
      use ions_base, only: na
      use io_global, only: stdout
      use mp_global, only: intra_image_comm
      use cvan
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum, mp_bcast
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb,qq,nhsavb=>nkbus
      use io_global, ONLY: ionode, ionode_id
      use twin_types

      implicit none

      complex(DP) :: betae(ngw,nhsa)
!       real(DP)    :: m_minus1(nhsavb,nhsavb)
      type(twin_matrix) :: m_minus1
      real(DP)    :: ema0bg(ngw)
      logical     :: use_ema


! local variables
      real(DP),allocatable :: q_matrix(:,:) ,c_matrix(:,:)
      complex(DP), allocatable :: q_matrix_c(:,:), c_matrix_c(:,:)
      integer is, iv, jv, ia, inl, jnl, i, ig, js, ja
      complex(DP) :: sca
      integer info, lwork
      integer, allocatable :: ipiv(:)
      real(dp),allocatable :: work(:)
      complex(dp),allocatable :: work_c(:)
      complex(dp), parameter :: c_zero=CMPLX(0.d0,0.d0), c_one=CMPLX(1.d0,0.d0)
      complex(dp), parameter :: c_mone=CMPLX(-1.d0,0.d0)

      call start_clock('set_x_minus1')
      allocate(ipiv(nhsavb))


      lwork=nhsavb

      IF(.not.m_minus1%iscmplx) THEN
         allocate(q_matrix(nhsavb,nhsavb),c_matrix(nhsavb,nhsavb))
         allocate(work(nhsavb))
      ELSE
	 allocate(q_matrix_c(nhsavb,nhsavb),c_matrix_c(nhsavb,nhsavb))
	 allocate(work_c(nhsavb))
      ENDIF
!construct q matrix
      IF(.not.m_minus1%iscmplx) THEN
	q_matrix(:,:) = 0.d0

	do is=1,nvb
	  do iv=1,nh(is)
	      do jv=1,nh(is)
		do ia=1,na(is)
		      inl=ish(is)+(iv-1)*na(is)+ia
		      jnl=ish(is)+(jv-1)*na(is)+ia
		      q_matrix(inl,jnl)= qq(iv,jv,is)
		enddo
	      enddo
	  enddo
	enddo
      ELSE
	q_matrix_c(:,:) = CMPLX(0.d0,0.d0)

	do is=1,nvb
	  do iv=1,nh(is)
	      do jv=1,nh(is)
		do ia=1,na(is)
		      inl=ish(is)+(iv-1)*na(is)+ia
		      jnl=ish(is)+(jv-1)*na(is)+ia
		      q_matrix_c(inl,jnl)= CMPLX(qq(iv,jv,is),0.d0)
		enddo
	      enddo
	  enddo
	enddo
      ENDIF
!construct b matrix
! m_minus1 used to be b matrix
      call set_twin(m_minus1, CMPLX(0.d0,0.d0))
!       m_minus1(:,:) = 0.d0
      do is=1,nvb
         do ia=1,na(is)
            do iv=1,nh(is)
               do js=1,nvb
                  do ja=1,na(js)
                     do jv=1,nh(js)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(js)+(jv-1)*na(js)+ja
                        sca=CMPLX(0.d0,0.d0)
                        if (use_ema) then
                           !k_minus case
                           IF(.not.m_minus1%iscmplx) THEN
                              !
                              do  ig=1,ngw           !loop on g vectors
                                 sca=sca+ema0bg(ig)*DBLE(CONJG(betae(ig,inl))*betae(ig,jnl))
                              enddo
                              !
                              sca = sca*2.0d0  !2. for real weavefunctions
                              !
                              if (ng0.eq.2) sca = sca - ema0bg(1)*DBLE(CONJG(betae(1,inl))*betae(1,jnl))
                              !
                           ELSE
                              !
                              do  ig=1,ngw           !loop on g vectors
                                 !
                                 sca=sca+ema0bg(ig)*CONJG(betae(ig,inl))*(betae(ig,jnl))
                                 !
                              enddo
                              !
                           ENDIF
                        else
                           ! s_minus case
			  IF(.not.m_minus1%iscmplx) THEN
			    do  ig=1,ngw           !loop on g vectors
			      sca=sca+DBLE(CONJG(betae(ig,inl))*betae(ig,jnl))
			    enddo
			    sca = sca*2.0d0  !2. for real weavefunctions
			    if (ng0.eq.2) sca = sca - DBLE(CONJG(betae(1,inl))*betae(1,jnl))
                          ELSE
			    do  ig=1,ngw           !loop on g vectors
			      sca=sca+CONJG(betae(ig,inl))*(betae(ig,jnl))
			    enddo
                          ENDIF
                        endif
                        !
                        call set_twin(m_minus1,inl,jnl,sca)
!                         write(6,*) "sca", sca, ema0bg(1), use_ema, betae(3,1)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      call twin_mp_sum( m_minus1)
!calculate -(1+QB)**(-1) * Q
      IF(.not.m_minus1%iscmplx) THEN
	CALL DGEMM('N','N',nhsavb,nhsavb,nhsavb,1.0d0,q_matrix,nhsavb,m_minus1%rvec,nhsavb,0.0d0,c_matrix,nhsavb)
	do i=1,nhsavb
         c_matrix(i,i)=c_matrix(i,i)+1.d0
	enddo
      ELSE
	CALL ZGEMM('N','N',nhsavb,nhsavb,nhsavb,c_one,q_matrix_c,nhsavb,m_minus1%cvec,nhsavb,c_zero,c_matrix_c,nhsavb) !warning:giovanni conjugate?
	do i=1,nhsavb
         c_matrix_c(i,i)=c_matrix_c(i,i)+CMPLX(1.d0,0.d0)
	enddo
      ENDIF
      if(ionode) then
	IF(.not.m_minus1%iscmplx) THEN
	  call dgetrf(nhsavb,nhsavb,c_matrix,nhsavb,ipiv,info)
	  if(info .ne. 0) write(stdout,*) 'set_k_minus1 Problem with dgetrf :', info
	  call dgetri(nhsavb,c_matrix,nhsavb,ipiv,work,lwork,info)
	  if(info .ne. 0) write(stdout,*) 'set_k_minus1 Problem with dgetri :', info
        ELSE
	  call zgetrf(nhsavb,nhsavb,c_matrix_c,nhsavb,ipiv,info)
	  if(info .ne. 0) write(stdout,*) 'set_k_minus1 Problem with dgetrf :', info
	  call zgetri(nhsavb,c_matrix_c,nhsavb,ipiv,work_c,lwork,info)
	  if(info .ne. 0) write(stdout,*) 'set_k_minus1 Problem with dgetri :', info
        ENDIF
      endif
      IF(.not.m_minus1%iscmplx) THEN
	call mp_bcast( c_matrix, ionode_id, intra_image_comm )
	CALL DGEMM('N','N',nhsavb,nhsavb,nhsavb,-1.0d0,c_matrix,nhsavb,q_matrix,nhsavb,0.0d0,m_minus1%rvec,nhsavb) 
      ELSE
	call mp_bcast( c_matrix_c, ionode_id, intra_image_comm )
	CALL ZGEMM('N','N',nhsavb,nhsavb,nhsavb,c_mone,c_matrix_c,nhsavb,q_matrix_c,nhsavb,c_zero,m_minus1%cvec,nhsavb) !warning:giovanni put a conjugate?
      ENDIF
      IF(.not.m_minus1%iscmplx) THEN
	deallocate(q_matrix,c_matrix, work)
      ELSE
	deallocate(q_matrix_c,c_matrix_c, work_c)
      ENDIF
      !
      deallocate(ipiv)
      !
!       call set_twin(m_minus1,CMPLX(0.d0,0.d0))
      call stop_clock('set_x_minus1')
      return
    end subroutine set_x_minus1_twin
!
      subroutine xminus1_real(c0,betae,ema0bg,beck,m_minus1,do_k)
! if (do_k) then
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = K^{-1}|c0>
! else
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = s^{-1}|c0> 
! endif
      use kinds, only: dp
      use ions_base, only: na
      use mp_global, only: intra_image_comm
      use cvan
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb, nhsavb=>nkbus
      use electrons_base, only: n => nbsp
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use mp, only: mp_sum
      use reciprocal_vectors, only: ng0 => gstart
!
      implicit none
      complex(dp) c0(ngw,n), betae(ngw,nhsa)
      real(dp)    beck(nhsa,n), ema0bg(ngw)
      real(DP)    :: m_minus1(nhsavb,nhsavb)
      logical :: do_k
! local variables
      complex(dp), allocatable :: phi(:,:)
      real(dp) , allocatable   :: qtemp(:,:)
      integer is, iv, ia, inl, i, j ,ig
      real(dp) becktmp

      
      logical :: mat_par=.true.!if true uses parallel routines      

      call start_clock('xminus1')
!$$
!      if(ionode) write(700,*) 'nvb is',nvb
!$$
      if (nvb.gt.0) then
!calculates beck
         if (do_k) then
            beck(:,:) = 0.d0

            do is=1,nvb
               do iv=1,nh(is)
                  do ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     do i=1,n
                        becktmp = 0.0d0
                        do ig=1,ngw
                           becktmp=becktmp+ema0bg(ig)*DBLE(CONJG(betae(ig,inl))*c0(ig,i))
                        enddo
                        becktmp = becktmp*2.0d0
                        if (ng0.eq.2) becktmp = becktmp-ema0bg(1)*DBLE(CONJG(betae(1,inl))*c0(1,i)) 
                        beck(inl,i) = beck(inl,i) + becktmp
                     enddo
                  enddo
               enddo
            enddo
            call mp_sum( beck, intra_image_comm )
         endif
!
!
      allocate(phi(ngw,n))
      allocate(qtemp(nhsavb,n))
      phi(1:ngw,1:n) = 0.0d0
      qtemp(:,:) = 0.0d0
      if(.not.mat_par) then
         call dgemm( 'N', 'N', nhsavb, n, nhsavb, 1.0d0, m_minus1,nhsavb ,    &
                    beck, nhsa, 0.0d0, qtemp,nhsavb )
      else
         call para_dgemm( 'N', 'N', nhsavb, n, nhsavb, 1.0d0, m_minus1,nhsavb ,    &
                    beck, nhsa, 0.0d0, qtemp,nhsavb,intra_image_comm )
      endif

!NB  nhsavb is the total number of US projectors
!    it works because the first pseudos are the vanderbilt's ones

         CALL DGEMM( 'N', 'N', 2*ngw, n, nhsavb, 1.0d0, betae, 2*ngw,    &
                    qtemp, nhsavb, 0.0d0, phi, 2*ngw )
         if (do_k) then
            do j=1,n
               do ig=1,ngw
                  c0(ig,j)=(phi(ig,j)+c0(ig,j))*ema0bg(ig)
               end do
            end do
         else
            do j=1,n
               do i=1,ngw
                  c0(i,j)=(phi(i,j)+c0(i,j))
               end do
            end do
         endif
      deallocate(qtemp,phi)

      else
         if (do_k) then
            do j=1,n
               do ig=1,ngw
                  c0(ig,j)=c0(ig,j)*ema0bg(ig)
               end do
            end do
         endif
      endif
      call stop_clock('xminus1')
      return
     end subroutine xminus1_real

!
      subroutine xminus1_twin(c0,betae,ema0bg,beck,m_minus1,do_k)
! if (do_k) then
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = K^{-1}|c0>
! else
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = s^{-1}|c0> 
! endif
      use kinds, only: dp
      use ions_base, only: na
      use mp_global, only: intra_image_comm
      use cvan
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb, nhsavb=>nkbus
      use electrons_base, only: n => nbsp
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use mp, only: mp_sum
      use reciprocal_vectors, only: ng0 => gstart
      use twin_types
!
      implicit none
      complex(dp) c0(ngw,n), betae(ngw,nhsa)
      real(dp) ::    ema0bg(ngw)
      type(twin_matrix) :: beck
!       complex(DP)    :: m_minus1(nhsavb,nhsavb)
      type(twin_matrix)    :: m_minus1 !(nhsavb,nhsavb)
      logical :: do_k
! local variables
      complex(dp), allocatable :: phi(:,:)
      real(dp) , allocatable   :: qtemp(:,:)
      complex(dp) , allocatable   :: qtemp_c(:,:)
      integer is, iv, ia, inl, i, j,ig
      real(dp) becktmp
      complex(dp) becktmp_c      
      logical :: mat_par=.true.!if true uses parallel routines
      complex(DP), parameter :: c_one=CMPLX(1.d0,0.d0)
      complex(DP), parameter :: c_zero=CMPLX(0.d0,0.d0)

      call start_clock('xminus1')
!$$
!      if(ionode) write(700,*) 'nvb is',nvb
!$$
      if (nvb.gt.0) then
!calculates beck
         if (do_k) then
            call set_twin(beck,CMPLX(0.d0,0.d0))
            do is=1,nvb
               do iv=1,nh(is)
                  do ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     do i=1,n
                        IF(.not.beck%iscmplx) THEN
			  becktmp = 0.0d0
			  do ig=1,ngw
			    becktmp=becktmp+ema0bg(ig)*DBLE(CONJG(betae(ig,inl))*c0(ig,i))
			  enddo
			  becktmp = becktmp*2.0d0
			  if (ng0.eq.2) becktmp = becktmp-ema0bg(1)*DBLE(CONJG(betae(1,inl))*c0(1,i)) 
			  beck%rvec(inl,i) = beck%rvec(inl,i) + becktmp
                        ELSE
			  becktmp_c = CMPLX(0.0d0, 0.d0)
			  do ig=1,ngw
			    becktmp_c=becktmp_c+ema0bg(ig)*(CONJG(betae(ig,inl))*c0(ig,i))
			  enddo
			  beck%cvec(inl,i) = beck%cvec(inl,i) + becktmp_c
                        ENDIF
                     enddo
                  enddo
               enddo
            enddo
	    call twin_mp_sum( beck )
         endif
!
!
	  allocate(phi(ngw,n))
	  phi(1:ngw,1:n) = 0.0d0
	  IF(.not.m_minus1%iscmplx) THEN
	    allocate(qtemp(nhsavb,n))
	    qtemp(:,:) = 0.0d0
	  ELSE
	    allocate(qtemp_c(nhsavb,n))
	    qtemp_c(:,:) = CMPLX(0.0d0, 0.d0)
	  ENDIF
	  if(.not.mat_par) then
	    IF(.not.m_minus1%iscmplx) THEN
	      call dgemm( 'N', 'N', nhsavb, n, nhsavb, 1.0d0, m_minus1%rvec,nhsavb ,    &
			  beck%rvec, nhsa, 0.0d0, qtemp,nhsavb )
	    ELSE
	      call zgemm( 'N', 'N', nhsavb, n, nhsavb, c_one, m_minus1%cvec,nhsavb ,    &
			  beck%cvec, nhsa, c_zero, qtemp_c, nhsavb )
	    ENDIF
	  else
	    IF(.not.m_minus1%iscmplx) THEN
	      call para_dgemm( 'N', 'N', nhsavb, n, nhsavb, 1.0d0, m_minus1%rvec,nhsavb ,    &
			  beck%rvec, nhsa, 0.0d0, qtemp,nhsavb,intra_image_comm )
	    ELSE
	      call para_zgemm( 'N', 'N', nhsavb, n, nhsavb, (1.0d0,0.d0), m_minus1%cvec,nhsavb ,    &
			  beck%cvec, nhsa, (0.0d0,0.d0), qtemp_c,nhsavb,intra_image_comm )
	    ENDIF
	  endif
!NB  nhsavb is the total number of US projectors
!    it works because the first pseudos are the vanderbilt's ones
	  IF(.not.m_minus1%iscmplx) THEN
	    CALL DGEMM( 'N', 'N', 2*ngw, n, nhsavb, 1.0d0, betae, 2*ngw,    &
			qtemp, nhsavb, 0.0d0, phi, 2*ngw )
	  ELSE
	    CALL ZGEMM( 'N', 'N', ngw, n, nhsavb, (1.0d0,0.d0), betae, ngw,    &
			qtemp_c, nhsavb, (0.0d0,0.d0), phi, ngw ) !warning:giovanni is it like this??
	  ENDIF
          if (do_k) then
            do j=1,n
               do ig=1,ngw
                  c0(ig,j)=(phi(ig,j)+c0(ig,j))*ema0bg(ig)
               end do
            end do
          else
            do j=1,n
               do i=1,ngw
                  c0(i,j)=(phi(i,j)+c0(i,j))
               end do
            end do
          endif
	  deallocate(phi)
	  IF(.not.m_minus1%iscmplx) THEN
	    deallocate(qtemp)
	  ELSE
	    deallocate(qtemp_c)
	  ENDIF
      else
         if (do_k) then
            do j=1,n
               do ig=1,ngw
                  c0(ig,j)=c0(ig,j)*ema0bg(ig)
               end do
            end do
         endif
      endif
      call stop_clock('xminus1')
      return
     end subroutine xminus1_twin

     subroutine xminus1_twin_new(c0,n,betae,ema0bg,beck,m_minus1,do_k)
     ! 
     ! if (do_k) then
     !-----------------------------------------------------------------------
     !     input: c0 , bec=<c0|beta>, betae=|beta>
     !     computes the matrix phi (with the old positions)
     !     where  |phi> = K^{-1}|c0>
     ! else
     !-----------------------------------------------------------------------
     !     input: c0 , bec=<c0|beta>, betae=|beta>
     !     computes the matrix phi (with the old positions)
     !       where  |phi> = s^{-1}|c0> 
     ! endif
      !
      use kinds, only: dp
      use ions_base, only: na
      use mp_global, only: intra_image_comm
      use cvan
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb, nhsavb=>nkbus
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use mp, only: mp_sum
      use reciprocal_vectors, only: ng0 => gstart
      use twin_types
      !
      implicit none
      !
      integer, intent(in) :: n
      real(dp)            :: ema0bg(ngw)
      complex(dp)         :: c0(ngw,n), betae(ngw,nhsa)
      type(twin_matrix)   :: beck
      type(twin_matrix)   :: m_minus1
      logical :: do_k
      !
      ! local variables
      !
      complex(dp), allocatable   :: phi(:,:)
      real(dp) ,   allocatable   :: qtemp(:,:)
      complex(dp), allocatable   :: qtemp_c(:,:)
      integer is, iv, ia, inl, i, j, ig
      real(dp) becktmp
      complex(dp) becktmp_c      
      logical :: mat_par=.true.!if true uses parallel routines
      complex(DP), parameter :: c_one=CMPLX(1.d0,0.d0)
      complex(DP), parameter :: c_zero=CMPLX(0.d0,0.d0)
      ! 
      call start_clock('xminus1')
      !
      if (nvb.gt.0) then
         if (do_k) then
            call set_twin(beck,CMPLX(0.d0,0.d0))
            do is=1,nvb
               do iv=1,nh(is)
                  do ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     do i=1,n
                        IF(.not.beck%iscmplx) THEN
                          becktmp = 0.0d0
                          do ig=1,ngw
                            becktmp=becktmp+ema0bg(ig)*DBLE(CONJG(betae(ig,inl))*c0(ig,i))
                          enddo
                          becktmp = becktmp*2.0d0
                          if (ng0.eq.2) becktmp = becktmp-ema0bg(1)*DBLE(CONJG(betae(1,inl))*c0(1,i)) 
                          beck%rvec(inl,i) = beck%rvec(inl,i) + becktmp
                        ELSE
                          becktmp_c = CMPLX(0.0d0, 0.d0)
                          do ig=1,ngw
                            becktmp_c=becktmp_c+ema0bg(ig)*(CONJG(betae(ig,inl))*c0(ig,i))
                          enddo
                          beck%cvec(inl,i) = beck%cvec(inl,i) + becktmp_c
                        ENDIF
                     enddo
                  enddo
               enddo
            enddo
            call twin_mp_sum( beck )
         endif
         !
         !
         allocate(phi(ngw,n))
          phi(1:ngw,1:n) = 0.0d0
          IF(.not.m_minus1%iscmplx) THEN
            allocate(qtemp(nhsavb,n))
            qtemp(:,:) = 0.0d0
          ELSE
            allocate(qtemp_c(nhsavb,n))
            qtemp_c(:,:) = CMPLX(0.0d0, 0.d0)
          ENDIF
          if(.not.mat_par) then
            IF(.not.m_minus1%iscmplx) THEN
              call dgemm( 'N', 'N', nhsavb, n, nhsavb, 1.0d0, m_minus1%rvec,nhsavb ,    &
                          beck%rvec, nhsa, 0.0d0, qtemp,nhsavb )
            ELSE
              call zgemm( 'N', 'N', nhsavb, n, nhsavb, c_one, m_minus1%cvec,nhsavb ,    &
                          beck%cvec, nhsa, c_zero, qtemp_c, nhsavb )
            ENDIF
          else
            IF(.not.m_minus1%iscmplx) THEN
              call para_dgemm( 'N', 'N', nhsavb, n, nhsavb, 1.0d0, m_minus1%rvec,nhsavb ,    &
                          beck%rvec, nhsa, 0.0d0, qtemp,nhsavb,intra_image_comm )
            ELSE
              call para_zgemm( 'N', 'N', nhsavb, n, nhsavb, (1.0d0,0.d0), m_minus1%cvec,nhsavb ,    &
                          beck%cvec, nhsa, (0.0d0,0.d0), qtemp_c,nhsavb,intra_image_comm )
            ENDIF
          endif
!NB  nhsavb is the total number of US projectors
!    it works because the first pseudos are the vanderbilt's ones
          IF(.not.m_minus1%iscmplx) THEN
            CALL DGEMM( 'N', 'N', 2*ngw, n, nhsavb, 1.0d0, betae, 2*ngw,    &
                        qtemp, nhsavb, 0.0d0, phi, 2*ngw )
          ELSE
            CALL ZGEMM( 'N', 'N', ngw, n, nhsavb, (1.0d0,0.d0), betae, ngw,    &
                        qtemp_c, nhsavb, (0.0d0,0.d0), phi, ngw ) !warning:giovanni is it like this??
          ENDIF
          if (do_k) then
            do j=1,n
               do ig=1,ngw
                  c0(ig,j)=(phi(ig,j)+c0(ig,j))*ema0bg(ig)
               end do
            end do
          else
            do j=1,n
               do i=1,ngw
                  c0(i,j)=(phi(i,j)+c0(i,j))
               end do
            end do
          endif
          deallocate(phi)
          IF(.not.m_minus1%iscmplx) THEN
            deallocate(qtemp)
          ELSE
            deallocate(qtemp_c)
          ENDIF
      else
         if (do_k) then
            do j=1,n
               do ig=1,ngw
                  c0(ig,j)=c0(ig,j)*ema0bg(ig)
               end do
            end do
         endif
      endif
      call stop_clock('xminus1')
      return
     end subroutine xminus1_twin_new     
     
     SUBROUTINE emass_precond_tpa( ema0bg, tpiba2, emaec )
       !
       ! kinetic energy preconditioning is computed here:
       ! (1+T')^{-1}
       !
       use kinds, ONLY : dp
       use gvecw, ONLY : ggp, ngw
       
       IMPLICIT NONE
       
       REAL(DP), INTENT(OUT) :: ema0bg(ngw)
       REAL(DP), INTENT(IN) ::  tpiba2, emaec
       INTEGER :: i

       real(DP) :: x

       call start_clock('emass_p_tpa')
       do i = 1, ngw
          !
          x=0.5d0*tpiba2*ggp(i)/emaec
          ema0bg(i) = 1.d0/(1.d0+(16.d0*x**4)/(27.d0+18.d0*x+12.d0*x**2+8.d0*x**3))
          !
       end do
       call stop_clock('emass_p_tpa')
       
     RETURN
     
     END SUBROUTINE emass_precond_tpa

      subroutine ave_kin( c, ngwx, n, ene_ave )
!this subroutine calculates the average kinetic energy of
!each state , to be used for preconditioning


      USE kinds,              ONLY: DP
      USE constants,          ONLY: pi, fpi
      USE gvecw,              ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart
      USE gvecw,              ONLY: ggp
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_image_comm
      USE cell_base,          ONLY: tpiba2
                                                                                                                             
      IMPLICIT NONE
                                                                                                                            
                                                                                                                             
      ! input
                                                                                                                             
      INTEGER,     INTENT(IN) :: ngwx, n
      COMPLEX(kind=DP), INTENT(IN) :: c( ngwx, n )
      REAL(kind=DP), INTENT(out) :: ene_ave(n)!average kinetic energy to be calculated
      !
      ! local
                                                                                                                             
      INTEGER  :: ig, i

      !
      DO i=1,n
         ene_ave(i)=0.d0
         DO ig=gstart,ngw
            ene_ave(i)=ene_ave(i)+DBLE(CONJG(c(ig,i))*c(ig,i))*ggp(ig)
         END DO
      END DO


      CALL mp_sum( ene_ave(1:n), intra_image_comm )
      ene_ave(:)=ene_ave(:)*tpiba2
                                                                                                                             
      RETURN
    END subroutine ave_kin



      subroutine xminus1_state(c0,betae,ema0bg,beck,m_minus1,do_k,ave_kin)
! if (do_k) then
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = K^{-1}|c0>
! else
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = s^{-1}|c0>
! endif
!adapted for state by state
      use kinds, only: dp
      use ions_base, only: na
      use mp_global, only: intra_image_comm
      use cvan
      use uspp_param, only: nh
      use uspp, only :nhsa=>nkb, nhsavb=>nkbus
      use electrons_base, only: n => nbsp
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use mp, only: mp_sum
      use reciprocal_vectors, only: ng0 => gstart
      USE gvecw,              ONLY: ggp
      USE cell_base,          ONLY: tpiba2


!
      implicit none
      complex(dp) c0(ngw,n), betae(ngw,nhsa)
      real(dp)    beck(nhsa,n), ema0bg(ngw)
      real(DP)    :: m_minus1(nhsavb,nhsavb)
      logical :: do_k
      real(kind=DP) :: ave_kin(n)!average kinetic energy per state
! local variables
      complex(dp), allocatable :: phi(:,:)
      real(dp) , allocatable   :: qtemp(:,:)
      integer is, iv, ia, inl, i, j, ig
      real(dp) becktmp
      real(kind=DP) :: prec_fact, x

 
      call start_clock('xminus1')
      if (nvb.gt.0) then
!calculates beck
         if (do_k) then
            beck(:,:) = 0.d0
 
            do is=1,nvb
               do iv=1,nh(is)
                  do ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     do i=1,n
                        becktmp = 0.0d0
                        do ig=1,ngw
                           becktmp=becktmp+ema0bg(ig)*DBLE(CONJG(betae(ig,inl))*c0(ig,i))
                        enddo
                        becktmp = becktmp*2.0d0
                        if (ng0.eq.2) becktmp = becktmp-ema0bg(1)*DBLE(CONJG(betae(1,inl))*c0(1,i))
                        beck(inl,i) = beck(inl,i) + becktmp
                     enddo
                  enddo
               enddo
            enddo
            call mp_sum( beck, intra_image_comm )
         endif
!
!
      allocate(phi(ngw,n))
      allocate(qtemp(nhsavb,n))
      phi(1:ngw,1:n) = 0.0d0
      qtemp(:,:) = 0.0d0
      call dgemm( 'N', 'N', nhsavb, n, nhsavb, 1.0d0, m_minus1,nhsavb ,    &
                    beck, nhsa, 0.0d0, qtemp,nhsavb )
 
 
 
!NB nhsavb is the total number of US projectors, it works because the first pseudos are the vanderbilt's ones

      CALL DGEMM( 'N', 'N', 2*ngw, n, nhsavb, 1.0d0, betae, 2*ngw,    &
           qtemp, nhsavb, 0.0d0, phi, 2*ngw )
      if (do_k) then
         do j=1,n
            do ig=1,ngw
               x=tpiba2*ggp(i)/ave_kin(j)
               prec_fact = 1.d0/(1.d0+(16.d0*x**4)/(27.d0+18.d0*x+12.d0*x**2+8.d0*x**3))
                c0(ig,j)=c0(ig,j)*prec_fact
               !c0(ig,j)=(phi(ig,j)+c0(ig,j))*ema0bg(ig)
            end do
         end do
      else
         do j=1,n
            do i=1,ngw
               c0(i,j)=(phi(i,j)+c0(i,j))
            end do
         end do
      endif
      deallocate(qtemp,phi)
      
   else
      if (do_k) then
         do j=1,n
            do ig=1,ngw
               x=tpiba2*ggp(ig)/ave_kin(j)
               prec_fact = 1.d0/(1.d0+(16.d0*x**4)/(27.d0+18.d0*x+12.d0*x**2+8.d0*x**3))
               c0(ig,j)=c0(ig,j)*prec_fact
            end do
         end do
      endif
   endif
   call stop_clock('xminus1')
   return
 end subroutine xminus1_state
!
! ... some simple routines for parallel linear algebra (the matrices are
! ... always replicated on all the cpus)
!
! ... written by carlo sbraccia ( 2006 )
!
!----------------------------------------------------------------------------
SUBROUTINE para_dgemm( transa, transb, m, n, k, &
                       alpha, a, lda, b, ldb, beta, c, ldc, comm )
  !----------------------------------------------------------------------------
  !
  ! ... trivial parallelization (splitting matrix B by columns) of DGEMM 
  !
  USE kinds, ONLY : DP
  USE parallel_toolkit
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=1), INTENT(IN)    :: transa, transb
  INTEGER,          INTENT(IN)    :: m, n, k
  REAL(DP),         INTENT(IN)    :: alpha, beta
  INTEGER,          INTENT(IN)    :: lda, ldb, ldc
  REAL(DP),         INTENT(INOUT) :: a(lda,*), b(ldb,*), c(ldc,*)
  INTEGER,          INTENT(IN)    :: comm
  !
  ! ... quick return if possible
  !
  IF ( m == 0 .OR. n == 0 .OR. &
       ( ( alpha == 0.0_DP .OR. k == 0 ) .AND. beta == 1.0_DP ) ) RETURN
  !
!write(*,*) 'DEBUG: para_dgemm'
  !
  CALL rep_matmul_drv( transa, transb, m, n, k, &
                       alpha, a, lda, b, ldb, beta, c, ldc, comm )
  RETURN
  !
END SUBROUTINE para_dgemm

SUBROUTINE para_zgemm( transa, transb, m, n, k, &
                       alpha, a, lda, b, ldb, beta, c, ldc, comm )
  !----------------------------------------------------------------------------
  !
  ! ... trivial parallelization (splitting matrix B by columns) of DGEMM 
  !
  USE kinds, ONLY : DP
  USE parallel_toolkit
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=1), INTENT(IN)    :: transa, transb
  INTEGER,          INTENT(IN)    :: m, n, k
  COMPLEX(DP),         INTENT(IN)    :: alpha, beta
  INTEGER,          INTENT(IN)    :: lda, ldb, ldc
  COMPLEX(DP),         INTENT(INOUT) :: a(lda,*), b(ldb,*), c(ldc,*)
  INTEGER,          INTENT(IN)    :: comm
  !
  ! ... quick return if possible
  !
  IF ( m == 0 .OR. n == 0 .OR. &
       ( ( alpha == 0.0_DP .OR. k == 0 ) .AND. beta == 1.0_DP ) ) RETURN
  !
!write(*,*) 'DEBUG: para_dgemm'
  !
  CALL zrep_matmul_drv( transa, transb, m, n, k, &
                       alpha, a, lda, b, ldb, beta, c, ldc, comm )
  RETURN
  !
END SUBROUTINE para_zgemm
