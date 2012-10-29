!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
   subroutine nlsm1_real ( n, nspmn, nspmx, eigr, c, becp )
!-----------------------------------------------------------------------

      !     computes: the array becp
      !     becp(ia,n,iv,is)=
      !         = sum_g [(-i)**l beta(g,iv,is) e^(-ig.r_ia)]^* c(g,n)
      !         = delta_l0 beta(g=0,iv,is) c(g=0,n)
      !          +sum_g> beta(g,iv,is) 2 re[(i)**l e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of c*(g)=c(-g)  (g> see routine ggen)
      !     input : beta(ig,l,is), eigr, c
      !     output: becp as parameter
      !
      USE kinds,      ONLY : DP
      USE mp,         ONLY : mp_sum
      USE mp_global,  ONLY : nproc_image, intra_image_comm
      USE ions_base,  only : na, nat
      USE gvecw,      only : ngw
      USE uspp,       only : nkb, nhtol, beta
      USE cvan,       only : ish
      USE uspp_param, only : nh
      !
      USE reciprocal_vectors, ONLY : gstart
!
      implicit none

      integer,   intent(in)  :: n, nspmn, nspmx
      complex(DP), intent(in)  :: eigr( ngw, nat ), c( ngw, n )
      real(DP), intent(out) :: becp( nkb, n )
      !
      integer   :: isa, ig, is, iv, ia, l, ixr, ixi, inl, i, nhx
      real(DP)  :: signre, signim, arg
      real(DP), allocatable :: becps( :, : )
      complex(DP), allocatable :: wrk2( :, : )
      complex(DP), parameter :: c_one=CMPLX(1.d0,0.d0), c_zero=CMPLX(0.d0,0.d0)
      complex(DP), parameter :: ci=CMPLX(0.d0,1.d0)
      complex(DP) :: cl, arg_c
      !
      call start_clock( 'nlsm1' )

      isa = 0
      do is = 1, nspmn - 1
        isa = isa + na(is)
      end do

      do is = nspmn, nspmx
         !
         IF( nh( is ) < 1 ) THEN
            isa = isa + na(is)
            CYCLE
         END IF
         !
         allocate( wrk2(  ngw, na( is ) ) )
         !
         IF( nproc_image > 1 ) THEN
            nhx = nh( is ) * na( is )
            IF( MOD( nhx, 2 ) /= 0 ) nhx = nhx + 1
            ALLOCATE( becps( nhx, n ) )
            becps = 0.0d0
         END IF
         !
         do iv = 1, nh( is )
            !
!$omp parallel default(shared), private(l,ixr,ixi,signre,signim,ig,arg,ia)
            l = nhtol( iv, is )
            cl = (-ci)**l
!             !
!             if (l == 0) then
!                ixr = 1
!                ixi = 2
!                signre =  1.0d0
!                signim =  1.0d0
!             else if (l == 1) then
!                ixr = 2
!                ixi = 1
!                signre =  1.0d0
!                signim = -1.0d0
!             else if (l == 2) then
!                ixr = 1
!                ixi = 2
!                signre = -1.0d0
!                signim = -1.0d0
!             else if (l == 3) then
!                ixr = 2
!                ixi = 1
!                signre = -1.0d0
!                signim =  1.0d0
!             endif
!
!$omp do
            do ia=1,na(is)
               !
               !  q = 0   component (with weight 1.0)
               !
               if (gstart == 2) then
                  wrk2(1, ia ) = cl*beta(1,iv,is)*eigr(1,ia+isa)
!                   wrk2( 2, 1, ia ) = signim*beta(1,iv,is)*eigr(ixi,1,ia+isa)
               end if
               !
               !   q > 0   components (with weight 2.0)
               !
               do ig = gstart, ngw
                  arg_c = CMPLX(2.0d0 * beta(ig,iv,is), 0.d0)*cl
                  wrk2( ig, ia ) = arg_c*eigr(ig,ia+isa)
!                   wrk2( 2, ig, ia ) = signim*arg*eigr(ixi,ig,ia+isa)
               end do
               !
            end do
!$omp end do
            
!$omp end parallel
            !
            IF( nproc_image > 1 ) THEN
               inl=(iv-1)*na(is)+1
               CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becps( inl, 1 ), nhx )
            ELSE
               inl=ish(is)+(iv-1)*na(is)+1
               CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becp( inl, 1 ), nkb )
            END IF

         end do

         deallocate( wrk2 )


         IF( nproc_image > 1 ) THEN
            !
            inl = ish(is) + 1
            !
            CALL mp_sum( becps, intra_image_comm )

            do i = 1, n
               do iv = inl , ( inl + na(is) * nh(is) - 1 )
                  becp( iv, i ) = becps( iv - inl + 1, i )
               end do
            end do

            DEALLOCATE( becps )

         END IF

         isa = isa + na(is)

      end do

      call stop_clock( 'nlsm1' )

      return
   end subroutine nlsm1_real
!--------------------------------------------------------l---------------
!

!-----------------------------------------------------------------------
   subroutine nlsm1_twin(n, nspmn, nspmx, eigr, c, becp, lbound_bec, lgam2)!added:giovanni lgam
!-----------------------------------------------------------------------

      !     computes: the array becp
      !     becp(ia,n,iv,is)=
      !         = sum_g [(-i)**l beta(g,iv,is) e^(-ig.r_ia)]^* c(g,n)
      !         = delta_l0 beta(g=0,iv,is) c(g=0,n)
      !          +sum_g> beta(g,iv,is) 2 re[(i)**l e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of c*(g)=c(-g)  (g> see routine ggen)
      !     input : beta(ig,l,is), eigr, c
      !     output: becp as parameter
      !
      USE kinds,      ONLY : DP
      USE mp,         ONLY : mp_sum
      USE mp_global,  ONLY : nproc_image, intra_image_comm
      USE ions_base,  only : na, nat
      USE gvecw,      only : ngw
      USE uspp,       only : nkb, nhtol, beta
      USE cvan,       only : ish
      USE uspp_param, only : nh
      USE twin_types !added:giovanni
      !
      USE reciprocal_vectors, ONLY : gstart
!
      implicit none

      integer,   intent(in)  :: n, nspmn, nspmx, lbound_bec
      complex(DP), intent(in) :: eigr(ngw,nat),c(ngw, n )!modified:giovanni
!       real(DP), intent(out) :: becp
      type(twin_matrix) :: becp!( nkb, n ) !modified:giovanni
      logical :: lgam2!added:giovanni
      !
      integer   :: isa, ig, is, iv, ia, l, ixr, ixi, inl, i, nhx
      real(DP)  :: signre, signim, arg
      real(DP), allocatable :: becps( :, : )
      complex(DP), allocatable :: becps_c( :, : )
      complex(DP), allocatable :: wrk2_c( :, : )
      complex(DP), parameter :: c_one=CMPLX(1.d0,0.d0), c_zero=CMPLX(0.d0,0.d0)
      complex(DP), parameter :: ci=CMPLX(0.d0,1.d0)
      complex(DP) :: cl, arg_c
      logical :: lgam!added:giovanni
      integer :: i1,i2,i3
      !
      lgam=lgam2
      call start_clock( 'nlsm1' )
      ! isa 
      isa = 0
      do is = 1, nspmn - 1
        isa = isa + na(is)
      end do

      IF((.not.lgam).and.nproc_image==1) THEN
	ALLOCATE( becps_c( nkb, n ))
	becps_c = CMPLX(0.0d0,0.d0)
      ENDIF

      do is = nspmn, nspmx
         !
         IF( nh( is ) < 1 ) THEN
            isa = isa + na(is)
            CYCLE
         END IF
         !
         allocate( wrk2_c( ngw, na( is ) ) )
         wrk2_c=CMPLX(0.d0,0.d0)
         !
         IF( nproc_image > 1 ) THEN
            nhx = nh( is ) * na( is )
            IF( MOD( nhx, 2 ) /= 0 ) nhx = nhx + 1 
            IF(lgam) THEN
	      ALLOCATE( becps( nhx, n ) )
	      becps = 0.0d0
	    ELSE
	      ALLOCATE( becps_c( nhx, n ) )
	      becps_c = CMPLX(0.0d0,0.d0)
	    ENDIF
         END IF


         !
         IF(lgam) THEN !added:giovanni
	  do iv = 1, nh( is )
	      !
  !$omp parallel default(shared), private(l,ixr,ixi,signre,signim,ig,arg,ia)
	      l = nhtol( iv, is )
!               write(6,*) "check_l_giovanni", l, iv, is
              cl =(-ci)**l
! 	      write(6,'(2((F20.13)(3x)))') cl
	      !
  !
  !$omp do
	      do ia=1,na(is)
		!
		!  q = 0   component (with weight 1.0)
		!
		if (gstart == 2) then
		    wrk2_c( 1, ia ) = CMPLX(beta(1,iv,is),0.d0)*cl*eigr(1,ia+isa)
		end if
		!
		!   q > 0   components (with weight 2.0)
		!
		do ig = gstart, ngw
		    arg_c = CMPLX(2.0d0*beta(ig,iv,is), 0.d0)*cl
		    wrk2_c( ig, ia ) = arg_c*eigr(ig,ia+isa)
		end do
		!
	      end do
  !$omp end do
	      
  !$omp end parallel
	      !
	      IF( nproc_image > 1 ) THEN
		inl=(iv-1)*na(is)+1
		CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2_c, 2*ngw, c, 2*ngw, 0.0d0, becps( inl, 1 ), nhx )
	      ELSE
                inl=ish(is)+(iv-1)*na(is)+1
                CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2_c, 2*ngw, c, 2*ngw, 0.0d0, becp%rvec( inl, lbound_bec ), nkb )		
	      END IF

	  end do
         ELSE
!begin_added:giovanni
	  do iv = 1, nh( is )
	      !
  !$omp parallel default(shared), private(l,ixr,ixi,signre,signim,ig,arg,ia)
	      l = nhtol( iv, is )
              cl=(-ci)**l
  !
  !$omp do
	      do ia=1,na(is)
                !
		do ig = 1, ngw
		    arg_c = cl*CMPLX(beta(ig,iv,is),0.d0)
		    wrk2_c( ig, ia ) = (arg_c*eigr(ig,ia+isa))
		end do
		!
	      end do
  !$omp end do

  !$omp end parallel
	      !
	      IF( nproc_image > 1 ) THEN
		inl=(iv-1)*na(is)+1
		CALL ZGEMM( 'C', 'N', na(is), n, ngw, c_one, wrk2_c, ngw, c, ngw, c_zero, becps_c( inl, 1 ), nhx )
	      ELSE
		inl=ish(is)+(iv-1)*na(is)+1
		CALL ZGEMM( 'C', 'N', na(is), n, ngw, c_one, wrk2_c, ngw, c, ngw,c_zero, becps_c( inl, lbound_bec ), nkb)
	      END IF

	  end do
!end_added:giovanni
         ENDIF

         deallocate( wrk2_c )

	 IF(nproc_image>1) THEN 

	      inl = ish(is) + 1

	      IF(lgam) THEN
		  !
		  CALL mp_sum( becps, intra_image_comm )

		  do i = 1, n
		    do iv = inl , ( inl + na(is) * nh(is) - 1 )
			becp%rvec( iv, i +lbound_bec -1) = becps( iv - inl + 1, i )
		    end do
		  end do

		  DEALLOCATE( becps )
		!
	      ELSE IF(.not.lgam) THEN

		CALL mp_sum( becps_c, intra_image_comm )

		    do i = 1, n
		      do iv = inl , ( inl + na(is) * nh(is) - 1 )
			  becp%cvec( iv, i +lbound_bec -1) = (becps_c( iv - inl + 1, i ))
		      end do
		    end do

		DEALLOCATE( becps_c )

	    ENDIF

	 END IF

         isa = isa + na(is)

      end do
!begin_added:giovanni
      IF(nproc_image==1.and.(.not.lgam)) THEN
          becp%cvec(1:nkb,lbound_bec:lbound_bec + n-1)=(becps_c(1:nkb,1:n))
          deallocate(becps_c)         
      endif
!end_added:giovanni
      call stop_clock( 'nlsm1' )
      
      return
   end subroutine nlsm1_twin
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
   subroutine nlsm1_dist_real ( n, nspmn, nspmx, eigr, c, becp, nlax, nspin, desc )
!-----------------------------------------------------------------------
      !  
      ! This version is for becp distributed over procs  
      !  
      !     computes: the array becp
      !     becp(ia,n,iv,is)=
      !         = sum_g [(-i)**l beta(g,iv,is) e^(-ig.r_ia)]^* c(g,n)
      !         = delta_l0 beta(g=0,iv,is) c(g=0,n)
      !          +sum_g> beta(g,iv,is) 2 re[(i)**l e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of c*(g)=c(-g)  (g> see routine ggen)
      !     input : beta(ig,l,is), eigr, c
      !     output: becp as parameter
      !
      USE kinds,      ONLY : DP
      USE mp,         ONLY : mp_sum
      USE mp_global,  ONLY : nproc_image, intra_image_comm
      USE ions_base,  only : na, nat
      USE gvecw,      only : ngw
      USE uspp,       only : nkb, nhtol, beta
      USE cvan,       only : ish
      USE uspp_param, only : nh
      !
      USE reciprocal_vectors, ONLY : gstart
      USE descriptors,        ONLY : descla_siz_ , lambda_node_ , nlar_ , ilar_ , la_n_
!
      implicit none

      integer,  intent(in)  :: n, nspmn, nspmx, nlax, nspin
      integer,  intent(in)  :: desc( descla_siz_ , nspin )
      real(DP), intent(in)  :: eigr( 2, ngw, nat ), c( 2, ngw, n )
      real(DP), intent(out) :: becp( nkb, nlax*nspin )
      !
      integer   :: isa, ig, is, iv, ia, l, ixr, ixi, inl, i, nhx
      integer   :: nr, ir, nup
      real(DP)  :: signre, signim, arg
      real(DP), allocatable :: becps( :, : )
      real(DP), allocatable :: wrk2( :, :, : )
      !
      call start_clock( 'nlsm1' )

      isa = 0
      do is = 1, nspmn - 1
        isa = isa + na(is)
      end do

      do is = nspmn, nspmx
         !
         IF( nh( is ) < 1 ) THEN
            isa = isa + na(is)
            CYCLE
         END IF
         !
         allocate( wrk2( 2, ngw, na( is ) ) )
         !
         IF( nproc_image > 1 ) THEN
            nhx = nh( is ) * na( is )
            IF( MOD( nhx, 2 ) /= 0 ) nhx = nhx + 1
            ALLOCATE( becps( nhx, n ) )
            becps = 0.0d0
         END IF
         !
         do iv = 1, nh( is )
            !
!$omp parallel default(shared), private(l,ixr,ixi,signre,signim,ig,arg,ia)
            l = nhtol( iv, is )
            !
            if (l == 0) then
               ixr = 1
               ixi = 2
               signre =  1.0d0
               signim =  1.0d0
            else if (l == 1) then
               ixr = 2
               ixi = 1
               signre =  1.0d0
               signim = -1.0d0
            else if (l == 2) then
               ixr = 1
               ixi = 2
               signre = -1.0d0
               signim = -1.0d0
            else if (l == 3) then
               ixr = 2
               ixi = 1
               signre = -1.0d0
               signim =  1.0d0
            endif
!
!$omp do
            do ia=1,na(is)
               !
               !  q = 0   component (with weight 1.0)
               !
               if (gstart == 2) then
                  wrk2( 1, 1, ia ) = signre*beta(1,iv,is)*eigr(ixr,1,ia+isa)
                  wrk2( 2, 1, ia ) = signim*beta(1,iv,is)*eigr(ixi,1,ia+isa)
               end if
               !
               !   q > 0   components (with weight 2.0)
               !
               do ig = gstart, ngw
                  arg = 2.0d0 * beta(ig,iv,is)
                  wrk2( 1, ig, ia ) = signre*arg*eigr(ixr,ig,ia+isa)
                  wrk2( 2, ig, ia ) = signim*arg*eigr(ixi,ig,ia+isa)
               end do
               !
            end do
!$omp end do
            
!$omp end parallel
            
            !
            IF( nproc_image > 1 ) THEN
               inl=(iv-1)*na(is)+1
               CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becps( inl, 1 ), nhx )
            ELSE
               inl=ish(is)+(iv-1)*na(is)+1
               CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becp( inl, 1 ), nkb )
            END IF

         end do

         deallocate( wrk2 )


         IF( nproc_image > 1 ) THEN
            !
            inl = ish(is) + 1
            !
            CALL mp_sum( becps, intra_image_comm )

            IF( desc( lambda_node_ , 1 ) > 0 ) THEN
               ir = desc( ilar_ , 1 )
               nr = desc( nlar_ , 1 )
               do i = 1, nr
                  do iv = inl , ( inl + na(is) * nh(is) - 1 )
                     becp( iv, i ) = becps( iv - inl + 1, i + ir - 1 )
                  end do
               end do
            END IF
            !
            IF( nspin == 2 ) THEN
               IF( desc( lambda_node_ , 2 ) > 0 ) THEN
                  nup = desc( la_n_ , 1 )
                  ir = desc( ilar_ , 2 )
                  nr = desc( nlar_ , 2 )
                  do i = 1, nr
                     do iv = inl , ( inl + na(is) * nh(is) - 1 )
                        becp( iv, i + nlax ) = becps( iv - inl + 1, i + ir - 1 + nup )
                     end do
                  end do
               END IF
            END IF

            DEALLOCATE( becps )

         END IF

         isa = isa + na(is)

      end do

      call stop_clock( 'nlsm1' )

      return
   end subroutine nlsm1_dist_real
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
   subroutine nlsm1_dist_twin ( n, nspmn, nspmx, eigr, c, becp, nlax, nspin, desc, lgam2 )
!-----------------------------------------------------------------------
      !  
      ! This version is for becp distributed over procs  
      !  
      !     computes: the array becp
      !     becp(ia,n,iv,is)=
      !         = sum_g [(-i)**l beta(g,iv,is) e^(-ig.r_ia)]^* c(g,n)
      !         = delta_l0 beta(g=0,iv,is) c(g=0,n)
      !          +sum_g> beta(g,iv,is) 2 re[(i)**l e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of c*(g)=c(-g)  (g> see routine ggen)
      !     input : beta(ig,l,is), eigr, c
      !     output: becp as parameter
      !
      USE kinds,      ONLY : DP
      USE mp,         ONLY : mp_sum
      USE mp_global,  ONLY : nproc_image, intra_image_comm
      USE ions_base,  only : na, nat
      USE gvecw,      only : ngw
      USE uspp,       only : nkb, nhtol, beta
      USE cvan,       only : ish
      USE uspp_param, only : nh
      !
      USE reciprocal_vectors, ONLY : gstart
      USE descriptors,        ONLY : descla_siz_ , lambda_node_ , nlar_ , ilar_ , la_n_
      USE twin_types
!
      implicit none

      integer,  intent(in)  :: n, nspmn, nspmx, nlax, nspin
      integer,  intent(in)  :: desc( descla_siz_ , nspin )
      complex(DP), intent(in)  :: eigr( ngw, nat ), c( ngw, n )
      type(twin_matrix), intent(out) :: becp !( nkb, nlax*nspin )
      logical, intent(IN) :: lgam2
      !
      integer   :: isa, ig, is, iv, ia, l, ixr, ixi, inl, i, nhx
      integer   :: nr, ir, nup
      real(DP)  :: signre, signim, arg
      real(DP), allocatable :: becps( :, : )
      complex(DP), allocatable :: becps_c( :, : )
      complex(DP), allocatable :: wrk2_c( :, : )
      complex(DP), parameter :: c_one=CMPLX(1.d0,0.d0), c_zero=CMPLX(0.d0,0.d0)
      complex(DP), parameter :: ci=CMPLX(0.d0,1.d0)
      complex(DP) :: cl, arg_c
      logical :: lgam !added:giovanni
      integer :: i1,i2,i3
      !
      lgam=lgam2
      call start_clock( 'nlsm1' )
      ! isa 
      isa = 0
      do is = 1, nspmn - 1
        isa = isa + na(is)
      end do

      IF(nproc_image==1.and.(.not.lgam)) THEN
	ALLOCATE( becps_c( nkb, n ))
	becps_c = CMPLX(0.0d0,0.d0)
      ENDIF

      do is = nspmn, nspmx
         !
         IF( nh( is ) < 1 ) THEN
            isa = isa + na(is)
            CYCLE
         END IF
         !
         allocate( wrk2_c( ngw, na( is ) ) )
         wrk2_c=CMPLX(0.d0,0.d0)
         !
         IF( nproc_image > 1 ) THEN
            nhx = nh( is ) * na( is )
            IF( MOD( nhx, 2 ) /= 0 ) nhx = nhx + 1 
            IF(lgam) THEN
	      ALLOCATE( becps( nhx, n ) )
	      becps = 0.0d0
	    ELSE
	      ALLOCATE( becps_c( nhx, n ) )
	      becps_c = CMPLX(0.0d0,0.d0)
	    ENDIF
         END IF
         !
         IF(lgam) THEN !added:giovanni
	  do iv = 1, nh( is )
	      !
  !$omp parallel default(shared), private(l,ixr,ixi,signre,signim,ig,arg,ia)
	      l = nhtol( iv, is )
              cl =(-ci)**l
	      !
  !
  !$omp do
	      do ia=1,na(is)
		!
		!  q = 0   component (with weight 1.0)
		!
		if (gstart == 2) then
		    wrk2_c( 1, ia ) =CMPLX(beta(1,iv,is),0.d0)*cl*eigr(1,ia+isa)
		end if
		!
		!   q > 0   components (with weight 2.0)
		!
		do ig = gstart, ngw
		    arg_c = CMPLX(2.0d0*beta(ig,iv,is), 0.d0)*cl
		    wrk2_c( ig, ia ) = arg_c*eigr(ig,ia+isa)
		end do
		!
	      end do
  !$omp end do
	      
  !$omp end parallel
	      !
	      IF( nproc_image > 1 ) THEN
		inl=(iv-1)*na(is)+1
		CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2_c, 2*ngw, c, 2*ngw, 0.0d0, becps( inl, 1 ), nhx )
	      ELSE
                inl=ish(is)+(iv-1)*na(is)+1
                CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2_c, 2*ngw, c, 2*ngw, 0.0d0, becp%rvec( inl, 1 ), nkb )		
	      END IF
	  end do

         ELSE
!begin_added:giovanni
	  do iv = 1, nh( is )
	      !
  !$omp parallel default(shared), private(l,ixr,ixi,signre,signim,ig,arg,ia)
	      l = nhtol( iv, is )
              cl=(-ci)**l
  !
  !$omp do
	      do ia=1,na(is)
                !
		do ig = 1, ngw
		    arg_c = cl*CMPLX(beta(ig,iv,is),0.d0)
		    wrk2_c( ig, ia ) = arg_c*eigr(ig,ia+isa)
		end do
		!
	      end do
  !$omp end do

  !$omp end parallel
	      !
	      IF( nproc_image > 1 ) THEN
		inl=(iv-1)*na(is)+1
		CALL ZGEMM( 'C', 'N', na(is), n, ngw, c_one, wrk2_c, ngw, c, ngw, c_zero, becps_c( inl, 1 ), nhx )
	      ELSE
		inl=ish(is)+(iv-1)*na(is)+1
		CALL ZGEMM( 'C', 'N', na(is), n, ngw, c_one, wrk2_c, ngw, c, ngw,c_zero, becps_c( inl, 1 ), nkb)
	      END IF
	  end do
!end_added:giovanni
         ENDIF

         deallocate( wrk2_c )

         IF( nproc_image > 1 ) THEN
            !
            inl = ish(is) + 1
            !
            IF(lgam) THEN

	      CALL mp_sum( becps, intra_image_comm )

	      IF( desc( lambda_node_ , 1 ) > 0 ) THEN
		ir = desc( ilar_ , 1 )
		nr = desc( nlar_ , 1 )
		do i = 1, nr
		    do iv = inl , ( inl + na(is) * nh(is) - 1 )
		      becp%rvec( iv, i ) = becps( iv - inl + 1, i + ir - 1 )
		    end do
		end do
	      END IF

	      IF( nspin == 2 ) THEN
		IF( desc( lambda_node_ , 2 ) > 0 ) THEN
		    nup = desc( la_n_ , 1 )
		    ir = desc( ilar_ , 2 )
		    nr = desc( nlar_ , 2 )
		    do i = 1, nr
		      do iv = inl , ( inl + na(is) * nh(is) - 1 )
			  becp%rvec( iv, i + nlax ) = becps( iv - inl + 1, i + ir - 1+  nup )
		      end do
		    end do
		END IF
	      END IF

              DEALLOCATE( becps )

            ELSE  !if not lgam

	      CALL mp_sum( becps_c, intra_image_comm )
	      IF( desc( lambda_node_ , 1 ) > 0 ) THEN
		ir = desc( ilar_ , 1 )
		nr = desc( nlar_ , 1 )
		do i = 1, nr
		    do iv = inl , ( inl + na(is) * nh(is) - 1 )
		      becp%cvec( iv, i ) = (becps_c( iv - inl + 1, i + ir - 1 ))
		    end do
		end do
	      END IF

	      IF( nspin == 2 ) THEN
		IF( desc( lambda_node_ , 2 ) > 0 ) THEN
		    nup = desc( la_n_ , 1 )
		    ir = desc( ilar_ , 2 )
		    nr = desc( nlar_ , 2 )
		    do i = 1, nr
		      do iv = inl , ( inl + na(is) * nh(is) - 1 )
			  becp%cvec( iv, i + nlax ) = (becps_c( iv - inl + 1, i +ir - 1+  nup ))
		      end do
		    end do
		END IF
	      END IF            

              DEALLOCATE( becps_c )
            ENDIF
         END IF

         isa = isa + na(is)

      end do
!       write(0,*) "inloop_giovanni", nproc_image==1.and.(.not.lgam).and.(.not.becp%iscmplx)
!begin_added:giovanni
      IF(nproc_image==1.and.(.not.lgam)) THEN
	  becp%cvec=(becps_c)
	  DEALLOCATE(becps_c)
      ENDIF

      call stop_clock( 'nlsm1' )

      return
   end subroutine nlsm1_dist_twin
!-----------------------------------------------------------------------


!-------------------------------------------------------------------------
   subroutine nlsm2( ngw, nkb, n, nspin, eigr, c, becdr, lgam2 )
!-----------------------------------------------------------------------

      !     computes: the array becdr
      !     becdr(ia,n,iv,is,k)
      !      =2.0 sum_g> g_k beta(g,iv,is) re[ (i)**(l+1) e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of  c*(g)=c(-g)  (g> see routine ggen)
      !     input : eigr, c
      !     output: becdr
      !
 
      USE kinds,      ONLY : DP
      use ions_base,  only : nsp, na, nat
      use uspp,       only : nhtol, beta  !, nkb
      use cvan,       only : ish
      use uspp_param, only : nh
      use cell_base,  only : tpiba
      use mp,         only : mp_sum
      use mp_global,  only : nproc_image, intra_image_comm
      use cp_main_variables,  only : nlax, descla, distribute_bec
      use reciprocal_vectors, only : gx, gstart
      use twin_types !added:giovanni
!
      implicit none
    
      integer,  intent(in)  :: ngw, nkb, n, nspin
      complex(DP), intent(in)  :: eigr(ngw,nat), c(ngw,n)
      type(twin_tensor), intent(out) :: becdr!(nkb,nspin*nlax,3) !modified:giovanni
      logical :: lgam2
      !
      real(DP), allocatable :: gk(:)
      complex(DP), allocatable :: wrk2_c(:,:)
      real(DP), allocatable :: becdr_repl(:,:)
      complex(DP), allocatable :: becdr_repl_c(:,:)
      !
      integer   :: ig, is, iv, ia, k, l, ixr, ixi, inl, isa, i
      real(DP) :: signre, signim, arg
      logical :: lgam ! added:giovanni
      complex(DP), parameter :: c_one=CMPLX(1.d0,0.d0), c_zero=CMPLX(0.d0,0.d0)
      complex(DP), parameter :: ci=CMPLX(0.d0,1.d0) !added:giovanni
      complex(DP) :: cl, arg_c !added:giovanni
      real(DP) :: fact
!
      lgam=lgam2

      call start_clock( 'nlsm2' )

      allocate( gk( ngw ) )

      IF(lgam) THEN
         becdr%rvec = 0.d0
         allocate( becdr_repl( nkb, n ) )
      ELSE
         becdr%cvec = c_zero
         allocate( becdr_repl_c( nkb, n ) )
      ENDIF
!
      do k = 1, 3

         IF(lgam) then
	    becdr_repl = 0.d0
         ELSE
	    becdr_repl_c = CMPLX(0.d0,0.d0)
         ENDIF

         do ig=1,ngw
            gk(ig)=gx(k,ig)*tpiba
         end do
!
         isa = 0

         do is=1,nsp

            allocate( wrk2_c( ngw, na( is ) ) )
            wrk2_c=CMPLX(0.d0,0.d0)
            
            IF(lgam) THEN 
                fact=2.d0
            ELSE
                fact=1.d0
            ENDIF

	    do iv=1,nh(is)
	      !
	      !     order of states:  s_1  p_x1  p_z1  p_y1  s_2  p_x2  p_z2  p_y2
	      !
!$omp parallel default(shared), private(l,ixr,ixi,signre,signim,ig,arg,ia)
	      l=nhtol(iv,is)
	      cl=(-ci)**(l+1)
!    
!$omp do
	      do ia=1,na(is)
		  !    q = 0   component (with weight 1.0)
		  if (gstart == 2) then
		    wrk2_c(1,ia) = cl*CMPLX(gk(1)*beta(1,iv,is),0.d0)*eigr(1,ia+isa)
		  end if
		  !    q > 0   components (with weight 2.0)
		  do ig=gstart,ngw
		    arg_c = CMPLX(fact*gk(ig)*beta(ig,iv,is),0.d0)
		    wrk2_c(ig,ia) = cl*arg_c*eigr(ig,ia+isa)
		  end do
	      end do
!$omp end do
!$omp end parallel 
	      inl=ish(is)+(iv-1)*na(is)+1
              IF(lgam) THEN
	         CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2_c, 2*ngw, c, 2*ngw, 0.0d0, becdr_repl( inl, 1 ), nkb )
              ELSE
                 CALL ZGEMM( 'C', 'N', na(is), n, ngw, c_one, wrk2_c, ngw, c, ngw, c_zero, becdr_repl_c( inl, 1 ), nkb )
              ENDIF
	    end do

            deallocate( wrk2_c )

            isa = isa + na(is)

         end do

         IF( nproc_image > 1 ) THEN
            IF(lgam) THEN
               CALL mp_sum( becdr_repl(:,:), intra_image_comm )
            ELSE
               CALL mp_sum( becdr_repl_c(:,:), intra_image_comm )
            ENDIF
         END IF

         IF(lgam) THEN
	    CALL distribute_bec( becdr_repl, becdr%rvec(:,:,k), descla, nspin )
         ELSE
            CALL distribute_bec((becdr_repl_c(:,:)), becdr%cvec(:,:,k), descla, nspin )
         ENDIF

      end do

      deallocate( gk )
      IF(lgam) THEN
	deallocate( becdr_repl )
      ELSE
	deallocate( becdr_repl_c )
      ENDIF

      call stop_clock( 'nlsm2' )
!
      return
   end subroutine nlsm2
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------
   subroutine nlsm2_repl( ngw, nkb, n, eigr, c, becdr, lgam )
!-----------------------------------------------------------------------

      !     computes: the array becdr
      !     becdr(ia,n,iv,is,k)
      !      =2.0 sum_g> g_k beta(g,iv,is) re[ (i)**(l+1) e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of  c*(g)=c(-g)  (g> see routine ggen)
      !     input : eigr, c
      !     output: becdr
      !
 
      USE kinds,      ONLY : DP
      use ions_base,  only : nsp, na, nat
      use uspp,       only : nhtol, beta  !, nkb
      use cvan,       only : ish
      use uspp_param, only : nh
      use cell_base,  only : tpiba
      use mp,         only : mp_sum
      use mp_global,  only : nproc_image, intra_image_comm
      use reciprocal_vectors, only : gx, gstart
      use twin_types !added:giovanni
!
      implicit none
    
      integer,  intent(in)  :: ngw, nkb, n
      complex(DP), intent(in)  :: eigr(ngw,nat), c(ngw,n)
      type(twin_tensor)     :: becdr
!      real(DP), intent(out) :: becdr(nkb,n,3)
      logical :: lgam
      !
      real(DP), allocatable :: gk(:)
      complex(DP), allocatable :: wrk2_c(:,:)
      !
      integer   :: ig, is, iv, ia, k, l, ixr, ixi, inl, isa, i
      real(DP) :: fact
      complex(DP), parameter :: c_one=CMPLX(1.d0,0.d0), c_zero=CMPLX(0.d0,0.d0)
      complex(DP), parameter :: ci=CMPLX(0.d0,1.d0) !added:giovanni
      complex(DP) :: cl, arg_c !added:giovanni
      !
      call start_clock( 'nlsm2' )

      allocate( gk( ngw ) )

      IF(lgam) THEN
         becdr%rvec = 0.d0
      ELSE
         becdr%cvec = c_zero
      ENDIF
!
      do k = 1, 3

         do ig=1,ngw
            gk(ig)=gx(k,ig)*tpiba
         end do
!
         isa = 0

         do is=1,nsp

            allocate( wrk2_c( ngw, na( is ) ) )
            !
            IF(lgam) THEN 
                fact=2.d0
            ELSE
                fact=1.d0
            ENDIF
            !
            do iv=1,nh(is)
               !
               !     order of states:  s_1  p_x1  p_z1  p_y1  s_2  p_x2  p_z2  p_y2
               !
!$omp parallel default(shared), private(l,ixr,ixi,signre,signim,ig,arg,ia)
               l=nhtol(iv,is)
               cl=(-ci)**(l+1)
!    
!$omp do
               do ia=1,na(is)
                  !    q = 0   component (with weight 1.0)
                  if (gstart == 2) then
                     wrk2_c(1,ia) = cl*CMPLX(gk(1)*beta(1,iv,is),0.d0)*eigr(1,ia+isa)
                  end if
                  !    q > 0   components (with weight 2.0)
                  do ig=gstart,ngw
                     arg_c = CMPLX(fact*gk(ig)*beta(ig,iv,is),0.d0)
                     wrk2_c(ig,ia) = cl*arg_c*eigr(ig,ia+isa)
                  end do
               end do
!$omp end do
!$omp end parallel 
               inl=ish(is)+(iv-1)*na(is)+1
               IF(lgam) THEN
                  CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2_c, 2*ngw, c, 2*ngw, 0.0d0, becdr%rvec( inl, 1, k ), nkb )
               ELSE
                  CALL ZGEMM( 'C', 'N', na(is), n, ngw, c_one, wrk2_c, ngw, c, ngw, c_zero, becdr%cvec( inl, 1, k ), nkb )
               ENDIF
            end do

            deallocate( wrk2_c )

            isa = isa + na(is)

         end do

         IF( nproc_image > 1 ) THEN
            IF(lgam) THEN
               CALL mp_sum( becdr%rvec(:,:,k), intra_image_comm )
            ELSE
               CALL mp_sum( becdr%cvec(:,:,k), intra_image_comm )
            ENDIF
         END IF
      end do

      deallocate( gk )

      call stop_clock( 'nlsm2' )
!
      return
   end subroutine nlsm2_repl
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
   real(8) function ennl_non_ortho( rhovan, bec, becdual )!added:giovanni lgam
!-----------------------------------------------------------------------
      !
      ! calculation of nonlocal potential energy term and array rhovan
      !
      use kinds,          only : DP
      use cvan,           only : ish
      use uspp_param,     only : nhm, nh
      use uspp,           only : nkb, dvan
      use electrons_base, only : n => nbsp, nspin, ispin, f
      use ions_base,      only : nsp, nat, na
      use twin_types
      use control_flags,  only : gamma_only, do_wf_cmplx
      !
      implicit none
      !
      ! input
      !
      type(twin_matrix) :: bec, becdual!( nkb, n )!modified:giovanni
      real(DP) :: rhovan( nhm*(nhm+1)/2, nat, nspin )
      !
      ! local
      !
      real(DP) :: sumt, sums(2), ennl_t
      complex(DP) :: sumt_c, sums_c(2), ennl_tc
      integer  :: is, iv, jv, ijv, inl, jnl, isa, isat, ism, ia, iss, i
      logical :: lgam!added:giovanni lgam
      !
      lgam=gamma_only.and..not.do_wf_cmplx
      !
      ennl_t = 0.d0 
      ennl_tc = CMPLX(0.d0,0.d0) 
      !
      !  xlf does not like name of function used for OpenMP reduction
      !
!$omp parallel default(shared), &
!$omp private(is,iv,jv,ijv,isa,isat,ism,ia,inl,jnl,sums,i,iss,sumt), reduction(+:ennl_t)
      if(.not.bec%iscmplx) then
        do is = 1, nsp
          do iv = 1, nh(is)
              do jv = iv, nh(is)
                ijv = (jv-1)*jv/2 + iv
                isa = 0
                do ism = 1, is - 1
                    isa = isa + na(ism)
                end do
!$omp do
                do ia = 1, na(is)
                    inl = ish(is)+(iv-1)*na(is)+ia
                    jnl = ish(is)+(jv-1)*na(is)+ia
                    isat = isa+ia
                    sums = 0.d0
                    do i = 1, n
                      iss = ispin(i)
                      sums(iss) = sums(iss) + 0.5d0*f(i) * (becdual%rvec(inl,i) * bec%rvec(jnl,i) +becdual%rvec(jnl,i) * bec%rvec(inl,i))
                    end do
                    sumt = 0.d0
                    do iss = 1, nspin
                      rhovan( ijv, isat, iss ) = sums( iss )
                      sumt = sumt + sums( iss )
                    end do
                    if( iv .ne. jv ) sumt = 2.d0 * sumt
                    ennl_t = ennl_t + sumt * dvan( jv, iv, is)
                end do
!$omp end do
              end do
          end do
        end do
      else
        do is = 1, nsp
          do iv = 1, nh(is)
              do jv = iv, nh(is)
                ijv = (jv-1)*jv/2 + iv
                isa = 0
                do ism = 1, is - 1
                    isa = isa + na(ism)
                end do
!$omp do
                do ia = 1, na(is)
                    inl = ish(is)+(iv-1)*na(is)+ia
                    jnl = ish(is)+(jv-1)*na(is)+ia
                    isat = isa+ia
                    sums_c = CMPLX(0.d0,0.d0)
                    do i = 1, n
                      iss = ispin(i)
                      sums_c(iss) = sums_c(iss) + CMPLX(0.5d0*f(i),0.d0)  &
      &                * (bec%cvec(inl,i) * CONJG(becdual%cvec(jnl,i))+bec%cvec(jnl,i) * CONJG(becdual%cvec(inl,i)))
                    end do
                    sumt_c = CMPLX(0.d0,0.d0)
                    do iss = 1, nspin
                      rhovan( ijv, isat, iss ) = DBLE(sums_c( iss ))
                      sumt_c = sumt_c + sums_c( iss )
                    end do
                    if( iv .ne. jv ) sumt_c = CMPLX(2.d0,0.d0) * sumt_c
                    ennl_tc = ennl_tc + sumt_c * CMPLX(dvan( jv, iv, is),0.d0)
                end do
!$omp end do
              end do
          end do
        end do
      endif
!$omp end parallel
      !
      if(.not.bec%iscmplx) then
        ennl_non_ortho = ennl_t
      else
        ennl_non_ortho = DBLE(ennl_tc)
      endif
      !
      return
   end function ennl_non_ortho
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
   real(8) function ennl( rhovan, bec )!added:giovanni lgam
!-----------------------------------------------------------------------
      !
      ! calculation of nonlocal potential energy term and array rhovan
      !
      use kinds,          only : DP
      use cvan,           only : ish
      use uspp_param,     only : nhm, nh
      use uspp,           only : nkb, dvan
      use electrons_base, only : n => nbsp, nspin, ispin, f
      use ions_base,      only : nsp, nat, na
      use twin_types
      use control_flags,  only : gamma_only, do_wf_cmplx
      !
      implicit none
      !
      ! input
      !
      type(twin_matrix) :: bec!( nkb, n )!modified:giovanni
      real(DP) :: rhovan( nhm*(nhm+1)/2, nat, nspin )
      !
      ! local
      !
      real(DP) :: sumt, sums(2), ennl_t
      complex(DP) :: sumt_c, sums_c(2), ennl_tc
      integer  :: is, iv, jv, ijv, inl, jnl, isa, isat, ism, ia, iss, i
      logical :: lgam!added:giovanni lgam
      !
      lgam=gamma_only.and..not.do_wf_cmplx
      !
      ennl_t = 0.d0 
      ennl_tc = CMPLX(0.d0,0.d0) 
      !
      !  xlf does not like name of function used for OpenMP reduction
      !
!$omp parallel default(shared), &
!$omp private(is,iv,jv,ijv,isa,isat,ism,ia,inl,jnl,sums,i,iss,sumt), reduction(+:ennl_t)
      if(.not.bec%iscmplx) then
        do is = 1, nsp
           do iv = 1, nh(is)
              do jv = iv, nh(is)
                 ijv = (jv-1)*jv/2 + iv
                 isa = 0
                 do ism = 1, is - 1
                    isa = isa + na(ism)
                 end do
!$omp do
                 do ia = 1, na(is)
                    inl = ish(is)+(iv-1)*na(is)+ia
                    jnl = ish(is)+(jv-1)*na(is)+ia
                    isat = isa+ia
                    sums = 0.d0
                    do i = 1, n
                      iss = ispin(i)
                      sums(iss) = sums(iss) + f(i) * bec%rvec(inl,i) * bec%rvec(jnl,i)
                    end do
                    sumt = 0.d0
                    do iss = 1, nspin
                       rhovan( ijv, isat, iss ) = sums( iss )
                       sumt = sumt + sums( iss )
                    end do
                    if( iv .ne. jv ) sumt = 2.d0 * sumt
                    ennl_t = ennl_t + sumt * dvan( jv, iv, is)
                 end do
!$omp end do
	      end do
	  end do
	end do
      else
        do is = 1, nsp
	  do iv = 1, nh(is)
	      do jv = iv, nh(is)
		ijv = (jv-1)*jv/2 + iv
		isa = 0
		do ism = 1, is - 1
		    isa = isa + na(ism)
		end do
!$omp do
		do ia = 1, na(is)
		    inl = ish(is)+(iv-1)*na(is)+ia
		    jnl = ish(is)+(jv-1)*na(is)+ia
		    isat = isa+ia
		    sums_c = CMPLX(0.d0,0.d0)
	            do i = 1, n
		      iss = ispin(i)
		      sums_c(iss) = sums_c(iss) + CMPLX(f(i),0.d0)  &
                      * ((bec%cvec(inl,i)) * CONJG(bec%cvec(jnl,i)))
		    end do
		    sumt_c = CMPLX(0.d0,0.d0)
		    do iss = 1, nspin
		      rhovan( ijv, isat, iss ) = DBLE(sums_c( iss ))
		      sumt_c = sumt_c + sums_c( iss )
		    end do
		    if( iv .ne. jv ) sumt_c = CMPLX(2.d0,0.d0) * sumt_c
		    ennl_tc = ennl_tc + sumt_c * CMPLX(dvan( jv, iv, is),0.d0)
		end do
!$omp end do
	      end do
	  end do
	end do
      endif
!$omp end parallel
      !
      if(.not.bec%iscmplx) then
        ennl = ennl_t
      else
        ennl = DBLE(ennl_tc)
      endif
      !
      return
   end function ennl
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
   real(8) function ennl_new( n, nspin, ispin, f, rhovan, bec )!added:giovanni lgam
!-----------------------------------------------------------------------
      !
      ! calculation of nonlocal potential energy term and array rhovan
      !
      use kinds,          only : DP
      use cvan,           only : ish
      use uspp_param,     only : nhm, nh
      use uspp,           only : nkb, dvan
      use ions_base,      only : nsp, nat, na
      use twin_types
      use control_flags,  only : gamma_only, do_wf_cmplx
      !
      implicit none
      !
      ! input
      !
      integer, intent(in) :: n, nspin, ispin(n), f(n)
      type(twin_matrix) :: bec!( nkb, n )!modified:giovanni
      real(DP) :: rhovan( nhm*(nhm+1)/2, nat, nspin )
      !
      ! local
      !
      real(DP) :: sumt, sums(2), ennl_t
      complex(DP) :: sumt_c, sums_c(2), ennl_tc
      integer  :: is, iv, jv, ijv, inl, jnl, isa, isat, ism, ia, iss, i
      logical :: lgam!added:giovanni lgam
      !
      lgam=gamma_only.and..not.do_wf_cmplx
      !
      ennl_t = 0.d0 
      ennl_tc = CMPLX(0.d0,0.d0) 
      !
      !  xlf does not like name of function used for OpenMP reduction
      !
!$omp parallel default(shared), &
!$omp private(is,iv,jv,ijv,isa,isat,ism,ia,inl,jnl,sums,i,iss,sumt), reduction(+:ennl_t)
      if(.not.bec%iscmplx) then
        do is = 1, nsp
           do iv = 1, nh(is)
              do jv = iv, nh(is)
                 ijv = (jv-1)*jv/2 + iv
                 isa = 0
                 do ism = 1, is - 1
                    isa = isa + na(ism)
                 end do
!$omp do
                 do ia = 1, na(is)
                    inl = ish(is)+(iv-1)*na(is)+ia
                    jnl = ish(is)+(jv-1)*na(is)+ia
                    isat = isa+ia
                    sums = 0.d0
                    do i = 1, n
                      iss = ispin(i)
                      sums(iss) = sums(iss) + f(i) * bec%rvec(inl,i) * bec%rvec(jnl,i)
                    end do
                    sumt = 0.d0
                    do iss = 1, nspin
                       rhovan( ijv, isat, iss ) = sums( iss )
                       sumt = sumt + sums( iss )
                    end do
                    if( iv .ne. jv ) sumt = 2.d0 * sumt
                    ennl_t = ennl_t + sumt * dvan( jv, iv, is)
                 end do
!$omp end do
              end do
          end do
        end do
      else
        do is = 1, nsp
          do iv = 1, nh(is)
              do jv = iv, nh(is)
                ijv = (jv-1)*jv/2 + iv
                isa = 0
                do ism = 1, is - 1
                    isa = isa + na(ism)
                end do
!$omp do
                do ia = 1, na(is)
                    inl = ish(is)+(iv-1)*na(is)+ia
                    jnl = ish(is)+(jv-1)*na(is)+ia
                    isat = isa+ia
                    sums_c = CMPLX(0.d0,0.d0)
                    do i = 1, n
                      iss = ispin(i)
                      sums_c(iss) = sums_c(iss) + CMPLX(f(i),0.d0)  &
                      * ((bec%cvec(inl,i)) * CONJG(bec%cvec(jnl,i)))
                    end do
                    sumt_c = CMPLX(0.d0,0.d0)
                    do iss = 1, nspin
                      rhovan( ijv, isat, iss ) = DBLE(sums_c( iss ))
                      sumt_c = sumt_c + sums_c( iss )
                    end do
                    if( iv .ne. jv ) sumt_c = CMPLX(2.d0,0.d0) * sumt_c
                    ennl_tc = ennl_tc + sumt_c * CMPLX(dvan( jv, iv, is),0.d0)
                end do
!$omp end do
              end do
          end do
        end do
      endif
!$omp end parallel
      !
      if(.not.bec%iscmplx) then
        ennl_new = ennl_t
      else
        ennl_new = DBLE(ennl_tc)
      endif
      !
      return
   end function ennl_new
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
   subroutine calrhovan_real( rhovan, bec, iwf )
!-----------------------------------------------------------------------
      !
      ! calculation of rhovan relative to state iwf
      !
      use kinds,          only : DP
      use cvan,           only : ish
      use uspp_param,     only : nhm, nh
      use uspp,           only : nkb, dvan
      use electrons_base, only : n => nbsp, nspin, ispin, f
      use ions_base,      only : nsp, nat, na
      !
      implicit none
      !
      ! input
      !
      real(DP) :: bec( nkb, n )
      real(DP) :: rhovan( nhm*(nhm+1)/2, nat, nspin )
      integer, intent(in) :: iwf
      !
      ! local
      !
      integer   :: is, iv, jv, ijv, inl, jnl, isa, ism, ia, iss
      !
      do is = 1, nsp
         do iv = 1, nh(is)
            do jv = iv, nh(is)
               ijv = (jv-1)*jv/2 + iv
               isa = 0
               do ism = 1, is - 1
                  isa = isa + na(ism)
               end do
               do ia = 1, na(is)
                  inl = ish(is)+(iv-1)*na(is)+ia
                  jnl = ish(is)+(jv-1)*na(is)+ia
                  isa = isa+1
                  iss = ispin(iwf)
                  rhovan( ijv, isa, iss ) = f(iwf) * bec(inl,iwf) * bec(jnl,iwf)
               end do
            end do
         end do
      end do
      !
      return
   end subroutine calrhovan_real
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine calrhovan_twin( rhovan, bec, iwf )
!-----------------------------------------------------------------------
      !
      ! calculation of rhovan relative to state iwf
      !
      use kinds,          only : DP
      use cvan,           only : ish
      use uspp_param,     only : nhm, nh
      use uspp,           only : nkb, dvan
      use electrons_base, only : n => nbsp, nspin, ispin, f
      use ions_base,      only : nsp, nat, na
      use twin_types
      use cp_main_variables, only : becdual
      use control_flags,     only : non_ortho
      use wavefunctions_module, only : cdual
      !
      implicit none
      !
      ! input
      !
      type(twin_matrix) :: bec !( nkb, n )
      real(DP) :: rhovan( nhm*(nhm+1)/2, nat, nspin )
      integer, intent(in) :: iwf
      !
      ! local
      !
      integer   :: is, iv, jv, ijv, inl, jnl, isa, ism, ia, iss
      !
      if(.not.bec%iscmplx) then
	  do is = 1, nsp
	    do iv = 1, nh(is)
		do jv = iv, nh(is)
		  ijv = (jv-1)*jv/2 + iv
		  isa = 0
		  do ism = 1, is - 1
		      isa = isa + na(ism)
		  end do
		  do ia = 1, na(is)
		      inl = ish(is)+(iv-1)*na(is)+ia
		      jnl = ish(is)+(jv-1)*na(is)+ia
		      isa = isa+1
		      iss = ispin(iwf)
		      IF(non_ortho) THEN
		         rhovan( ijv, isa, iss ) = f(iwf) * becdual%rvec(inl,iwf) * bec%rvec(jnl,iwf)
		      ELSE
		         rhovan( ijv, isa, iss ) = f(iwf) * bec%rvec(inl,iwf) * bec%rvec(jnl,iwf)
		      ENDIF
		  end do
		end do
	    end do
	  end do
       else
	  do is = 1, nsp
	    do iv = 1, nh(is)
		do jv = iv, nh(is)
		  ijv = (jv-1)*jv/2 + iv
		  isa = 0
		  do ism = 1, is - 1
		      isa = isa + na(ism)
		  end do
		  do ia = 1, na(is)
		      inl = ish(is)+(iv-1)*na(is)+ia
		      jnl = ish(is)+(jv-1)*na(is)+ia
		      isa = isa+1
		      iss = ispin(iwf)
		      IF(non_ortho) THEN
                          rhovan( ijv, isa, iss ) = f(iwf) * DBLE(CONJG(becdual%cvec(inl,iwf)) * &                                                        (bec%cvec(jnl,iwf)))
		      ELSE
                          rhovan( ijv, isa, iss ) = f(iwf) * DBLE(CONJG(bec%cvec(inl,iwf)) * &
                                                              (bec%cvec(jnl,iwf)))
                       ENDIF
		  end do
		end do
	    end do
	  end do
       endif
      !
      return
   end subroutine calrhovan_twin
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine calbec ( nspmn, nspmx, eigr, c, bec)
!-----------------------------------------------------------------------
      
      !     this routine calculates array bec
      !
      !        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
      !                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
      !
      !     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
      !
      
      USE kinds,          ONLY : DP
      use ions_base,      only : na, nat
      use io_global,      only : stdout
      use cvan,           only : ish
      use electrons_base, only : n => nbsp
      use gvecw,          only : ngw
      use control_flags,  only : iprint, iprsta, gamma_only, &
  &                              do_wf_cmplx
      use uspp_param,     only : nh
      use uspp,           only : nkb
      use twin_types
!
      implicit none
      !
      integer,     intent(in)  :: nspmn, nspmx
      type(twin_matrix)        :: bec!( nkb, n )
      complex(DP), intent(in)  :: c( ngw, n ), eigr( ngw,nat )

      ! local variables

      integer :: is, ia, i , iv
      logical :: lgam
!
      lgam=gamma_only.and..not.do_wf_cmplx

      call start_clock( 'calbec' )
      call nlsm1_twin( n, nspmn, nspmx, eigr, c, bec, 1, lgam )
!
      if ( iprsta > 2 ) then
         WRITE( stdout,*)
         do is=1,nspmx
            if(nspmx.gt.1) then
               WRITE( stdout,'(33x,a,i4)') ' calbec: bec (is)',is
               if(.not.bec%iscmplx) then
                 WRITE( stdout,'(8f9.4)')                                       &
     &              ((bec%rvec(ish(is)+(iv-1)*na(is)+1,i),iv=1,nh(is)),i=1,n)
               else
                 WRITE( stdout,'(8(2((f9.4),(f9.4))))') &
     &              ((bec%cvec(ish(is)+(iv-1)*na(is)+1,i),iv=1,nh(is)),i=1,n)
               endif
            else
               do ia=1,na(is)
                  WRITE( stdout,'(33x,a,i4)') ' calbec: bec (ia)',ia
                  if(lgam) then
                    WRITE( stdout,'(8f9.4)')                                    &
     &             ((bec%rvec(ish(is)+(iv-1)*na(is)+ia,i),iv=1,nh(is)),i=1,n)
                  else
                    WRITE( stdout,'(8(2((f9.4),(f9.4))))')                                    &
     &             ((bec%cvec(ish(is)+(iv-1)*na(is)+ia,i),iv=1,nh(is)),i=1,n)
                  endif
               end do
            end if
         end do
      endif
      call stop_clock( 'calbec' )
!
      return
   end subroutine calbec
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE caldbec( ngw, nkb, n, nspmn, nspmx, eigr, c, dbec ) !warning:giovanni this still does not work with complex wavefunctions
  !-----------------------------------------------------------------------
  !
  !     this routine calculates array dbec, derivative of bec:
  !
  !        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
  !                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
  !
  !     with respect to cell parameters h
  !
  !     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
  !
  USE kinds,      ONLY : DP
  use mp,         only : mp_sum
  use mp_global,  only : nproc_image, intra_image_comm
  use ions_base,  only : na, nat
  use cvan,       only : ish
  use cdvan,      only : dbeta
  use uspp,       only : nhtol
  use uspp_param, only : nh, nhm
  use reciprocal_vectors, only : gstart
  USE cp_main_variables,  ONLY : descla, la_proc, nlax, nlam
  USE descriptors,        ONLY : nlar_ , nlac_ , ilar_ , ilac_ , nlax_ , la_myr_ , la_myc_
  use electrons_base,     only : nspin, iupdwn, nupdwn
  !
  implicit none
  !
  integer,      intent(in)  :: ngw, nkb, n
  integer,      intent(in)  :: nspmn, nspmx
  complex(DP), intent(in)  :: c(ngw,n)
  real(DP),    intent(in)  :: eigr(2,ngw,nat)
  real(DP),    intent(out) :: dbec( nkb, 2*nlam, 3, 3 )
  !
  real(DP), allocatable :: wrk2(:,:,:), dwrk(:,:)
  !
  integer   :: ig, is, iv, ia, l, ixr, ixi, inl, i, j, ii, isa, nanh, iw, iss, nr, ir, istart, nss
  real(DP) :: signre, signim, arg
  !
  !
  !
  do j=1,3
     do i=1,3

        isa = 0
        do is = 1, nspmn - 1
          isa = isa + na(is)
        end do

        do is=nspmn,nspmx
           allocate( wrk2( 2, ngw, na(is) ) )
           nanh = na(is)*nh(is)
           allocate( dwrk( nanh, n ) )
           do iv=1,nh(is)
              l=nhtol(iv,is)
              if (l == 0) then
                 ixr = 1
                 ixi = 2
                 signre =  1.0d0
                 signim =  1.0d0
              else if (l == 1) then
                 ixr = 2
                 ixi = 1
                 signre =  1.0d0
                 signim = -1.0d0
              else if (l == 2) then
                 ixr = 1
                 ixi = 2
                 signre = -1.0d0
                 signim = -1.0d0
              else if (l == 3) then
                 ixr = 2
                 ixi = 1
                 signre = -1.0d0
                 signim =  1.0d0
              else
                 CALL errore(' caldbec  ', ' l not implemented ', ABS( l ) )
              endif
              !
              do ia=1,na(is)
                 if (gstart == 2) then
                    !     q = 0   component (with weight 1.0)
                    wrk2(1,1,ia)= signre*dbeta(1,iv,is,i,j)*eigr(ixr,1,ia+isa)
                    wrk2(2,1,ia)= signim*dbeta(1,iv,is,i,j)*eigr(ixi,1,ia+isa)
                 end if
                 !     q > 0   components (with weight 2.0)
                 do ig = gstart, ngw
                    arg = 2.0d0*dbeta(ig,iv,is,i,j)
                    wrk2(1,ig,ia) = signre*arg*eigr(ixr,ig,ia+isa)
                    wrk2(2,ig,ia) = signim*arg*eigr(ixi,ig,ia+isa)
                 end do
              end do
              inl=(iv-1)*na(is)+1
              CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, dwrk(inl,1), nanh )
           end do
           deallocate( wrk2 )
           if( nproc_image > 1 ) then
              call mp_sum( dwrk, intra_image_comm )
           end if
           inl=ish(is)+1
           do iss=1,nspin
              IF( la_proc ) THEN
                 nr = descla( nlar_ , iss )
                 ir = descla( ilar_ , iss )
                 istart = iupdwn( iss )
                 nss    = nupdwn( iss )
                 do ii = 1, nr
                    do iw = 1, nanh
                       dbec( iw + inl - 1, ii + (iss-1)*nspin, i, j ) = dwrk( iw, ii + ir - 1 + istart - 1 )
                    end do
                 end do
              END IF
           end do
           deallocate( dwrk )
           isa = isa + na(is)
        end do
     end do
  end do

  !
  return
end subroutine caldbec
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine dennl( bec, dbec, drhovan, denl )
  !-----------------------------------------------------------------------
  !
  !  compute the contribution of the non local part of the
  !  pseudopotentials to the derivative of E with respect to h
  !
  USE kinds,      ONLY : DP
  use cvan,       only : ish
  use uspp_param, only : nh, nhm
  use uspp,       only : nkb, dvan, deeq
  use ions_base,  only : nsp, na, nat
  use cell_base,  only : h
  use io_global,  only : stdout
  use mp,         only : mp_sum
  use mp_global,  only : intra_image_comm
  USE cp_main_variables,  ONLY : descla, la_proc, nlax, nlam
  USE descriptors,        ONLY : nlar_ , nlac_ , ilar_ , ilac_ , nlax_ , la_myr_ , la_myc_
  use electrons_base,     only : n => nbsp, ispin, f, nspin, iupdwn, nupdwn
  use reciprocal_vectors, only : gstart

  implicit none

  real(DP), intent(in)  :: dbec( nkb, 2*nlam, 3, 3 )
  real(DP), intent(in)  :: bec( nkb, n )
  real(DP), intent(out) :: drhovan( nhm*(nhm+1)/2, nat, nspin, 3, 3 )
  real(DP), intent(out) :: denl( 3, 3 )

  real(DP) :: dsum(3,3),dsums(2,3,3), detmp(3,3)
  integer   :: is, iv, jv, ijv, inl, jnl, isa, ism, ia, iss, i,j,k
  integer   :: istart, nss, ii, ir, nr
  !
  denl=0.d0
  drhovan=0.0d0

  IF( la_proc ) THEN


  do is=1,nsp
     do iv=1,nh(is)
        do jv=iv,nh(is)
           ijv = (jv-1)*jv/2 + iv
           isa=0
           do ism=1,is-1
              isa=isa+na(ism)
           end do
           do ia=1,na(is)
              inl=ish(is)+(iv-1)*na(is)+ia
              jnl=ish(is)+(jv-1)*na(is)+ia
              isa=isa+1
              dsums=0.d0
              do iss=1,nspin
                 IF( descla( la_myr_ , iss ) == descla( la_myc_ , iss ) ) THEN
                 nr = descla( nlar_ , iss )
                 ir = descla( ilar_ , iss )
                 istart = iupdwn( iss )
                 nss    = nupdwn( iss )
                 do i=1,nr
                    ii = i+istart-1+ir-1
                    do k=1,3
                       do j=1,3
                          dsums(iss,k,j)=dsums(iss,k,j)+f(ii)*       &
 &                          (dbec(inl,i+(iss-1)*nlam,k,j)*bec(jnl,ii)          &
 &                          + bec(inl,ii)*dbec(jnl,i+(iss-1)*nlam,k,j))
                       enddo
                    enddo
                 end do
                 END IF
              end do
              !
              do iss=1,nspin
                 IF( descla( la_myr_ , iss ) == descla( la_myc_ , iss ) ) THEN
                 dsum=0.d0
                 do k=1,3
                    do j=1,3
                       drhovan(ijv,isa,iss,j,k)=dsums(iss,j,k)
                       dsum(j,k)=dsum(j,k)+dsums(iss,j,k)
                    enddo
                 enddo
                 if(iv.ne.jv) dsum=2.d0*dsum
                 denl = denl + dsum * dvan(jv,iv,is)
                 END IF
              end do
           end do
        end do
     end do
  end do

  END IF

  CALL mp_sum( denl,    intra_image_comm )
  do k=1,3
     do j=1,3
        CALL mp_sum( drhovan(:,:,:,j,k), intra_image_comm )
     end do
  end do

!  WRITE(6,*) 'DEBUG enl (CP) = '
!  detmp = denl
!  detmp = MATMUL( detmp(:,:), TRANSPOSE( h ) )
!  WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)

  !
  return
end subroutine dennl
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine nlfq( c, eigr, bec, becdr, fion, lgam2)
  !-----------------------------------------------------------------------
  !
  !     contribution to fion due to nonlocal part
  !
  USE kinds,          ONLY : DP
  use uspp,           only : nkb, dvan, deeq
  use uspp_param,     only : nhm, nh
  use cvan,           only : ish, nvb
  use ions_base,      only : nax, nat, nsp, na
  use electrons_base, only : n => nbsp, ispin, f, nspin, iupdwn, nupdwn
  use gvecw,          only : ngw
  use constants,      only : pi, fpi
  use mp_global,      only : me_image, intra_image_comm, nproc_image
  use mp,             only : mp_sum
  USE cp_main_variables, ONLY: nlax, descla, la_proc
  USE descriptors,       ONLY: nlar_ , nlac_ , ilar_ , ilac_ , lambda_node_ , &
                               la_myr_ , la_myc_
  USE twin_types !added:giovanni
  !
  implicit none
  !
  type(twin_matrix) :: bec!(nkb,n) !modified:giovanni 
  complex(DP) :: c(ngw, n) !modified:giovanni
  type(twin_tensor) :: becdr!( nkb, nspin*nlax, 3 ) !modified:giovanni
  complex(DP), intent(in)  :: eigr( ngw, nat )
  real(DP),    intent(out) :: fion( 3, nat )
  logical, intent(IN) :: lgam2
  !
  integer   :: k, is, ia, isa, iss, inl, iv, jv, i, ir, nr, nss, istart, ioff
  real(DP) :: temp
  !
  real(DP), allocatable :: tmpbec(:,:), tmpdr(:,:)
  complex(DP), allocatable :: tmpbec_c(:,:), tmpdr_c(:,:) 
  real(DP), allocatable :: fion_loc(:,:)
  logical :: lgam !added:giovanni:debug
  character(len=5) :: subname = "nlfq"
#ifdef __OPENMP 
  INTEGER :: mytid, ntids, omp_get_thread_num, omp_get_num_threads
#endif  
  !
  call start_clock( 'nlfq' )
  !
  lgam=lgam2.or..not.becdr%iscmplx !added:giovanni
  !
  IF(becdr%iscmplx.neqv.bec%iscmplx) THEN !added:giovanni:debug
     call errore(subname, "incompatible twin types", 1)
     stop
  ENDIF
  !
  !     nlsm2 fills becdr
  !
  call nlsm2( ngw, nkb, n, nspin, eigr, c, becdr, lgam2 )
  !
  allocate ( fion_loc( 3, nat ) )
  !
  fion_loc = 0.0d0
  !
  DO k = 1, 3

!$omp parallel default(shared), &
!$omp private(tmpbec,tmpdr,isa,is,ia,iss,nss,istart,ir,nr,ioff,iv,jv,inl,temp,i,mytid,ntids)
#ifdef __OPENMP
     mytid = omp_get_thread_num()  ! take the thread ID
     ntids = omp_get_num_threads() ! take the number of threads
#endif

     IF(lgam) THEN
        allocate ( tmpbec( nhm, nlax ), tmpdr( nhm, nlax ) )
     ELSE
        allocate ( tmpbec_c( nhm, nlax ), tmpdr_c( nhm, nlax ) )
     ENDIF

     isa = 0
     !
     DO is=1,nsp
        DO ia=1,na(is)

           isa=isa+1

#ifdef __OPENMP
           ! distribute atoms round robin to threads
           !
           IF( MOD( isa, ntids ) /= mytid ) CYCLE
#endif  
           DO iss = 1, nspin

              nss = nupdwn( iss )
              istart = iupdwn( iss )

              IF( la_proc .AND. &
                  ( descla( la_myr_ , iss ) == descla( la_myc_ , iss ) ) ) THEN

                 ! only processors on the diagonal of the square proc grid enter here.
                 ! This is to distribute the load among different multi-core nodes,
                 ! and maximize the memory bandwith per core.
                 IF(lgam) THEN
                    tmpbec = 0.d0
                    tmpdr  = 0.d0
                 ELSE
                    tmpbec_c = CMPLX(0.d0, 0.d0)
                    tmpdr_c  = CMPLX(0.d0, 0.d0)
                 ENDIF

                 ir = descla( ilar_ , iss )
                 nr = descla( nlar_ , iss )

                 ioff = istart-1+ir-1
                 
                 IF(lgam) THEN
		    do iv=1,nh(is)
			do jv=1,nh(is)
			  inl=ish(is)+(jv-1)*na(is)+ia
			  temp=dvan(iv,jv,is)+deeq(jv,iv,isa,iss)
			  do i=1,nr
			      tmpbec(iv,i)=tmpbec(iv,i)+temp*bec%rvec(inl,i+ioff)
			  end do
			end do
		    end do
                 ELSE
		    do iv=1,nh(is)
			do jv=1,nh(is)
			  inl=ish(is)+(jv-1)*na(is)+ia
			  temp=dvan(iv,jv,is)+deeq(jv,iv,isa,iss)
			  do i=1,nr
			      tmpbec_c(iv,i)=tmpbec_c(iv,i)+temp*bec%cvec(inl,i+ioff)
			  end do
			end do
		    end do
                 ENDIF

                 IF(lgam) THEN
		    do iv=1,nh(is)
			inl=ish(is)+(iv-1)*na(is)+ia
			do i=1,nr
			  tmpdr(iv,i)=f(i+ioff)*becdr%rvec( inl, i+(iss-1)*nlax, k )
			end do
		    end do
                 ELSE
		    do iv=1,nh(is)
			inl=ish(is)+(iv-1)*na(is)+ia
			do i=1,nr
			  tmpdr_c(iv,i)=f(i+ioff)*becdr%cvec( inl, i+(iss-1)*nlax, k )
			end do
		    end do
                 ENDIF

                 IF(lgam) THEN
		    do i=1,nr
			do iv=1,nh(is)
			  tmpdr(iv,i)=tmpdr(iv,i)*tmpbec(iv,i)
			end do
		    end do
                 ELSE
		    do i=1,nr
			do iv=1,nh(is)
			  tmpdr_c(iv,i)=tmpdr_c(iv,i)*CONJG(tmpbec_c(iv,i))
			end do
		    end do
                 ENDIF

                 IF(lgam) THEN
                    fion_loc(k,isa) = fion_loc(k,isa)-2.d0*SUM(tmpdr)
                 ELSE
                    fion_loc(k,isa) = fion_loc(k,isa)-2.d0*DBLE(SUM(tmpdr_c))
                 ENDIF

              END IF
           END DO
        END DO
     END DO
     
     IF(lgam) THEN
        deallocate (tmpbec, tmpdr)
     ELSE
        deallocate (tmpbec_c, tmpdr_c)
     ENDIF
!$omp end parallel
  END DO
  !
  CALL mp_sum( fion_loc, intra_image_comm )
  !
  fion = fion + fion_loc
  !
  !     end of x/y/z loop
  !
  deallocate ( fion_loc )
  !
  call stop_clock( 'nlfq' )
  !
  return
end subroutine nlfq

