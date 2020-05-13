
subroutine pc3nc_fixed(phi, grad, lgam)
      !
      use kinds
      use io_global, only: stdout
      use mp_global, only: intra_image_comm
      use gvecw, only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      use mp, only: mp_sum
      use electrons_base, only: n => nbsp, ispin
      use input_parameters,only : fixed_band
      !
      implicit none
      !
      complex(dp), intent(in):: phi(ngw,n)
      complex(dp), intent(inout):: grad(ngw, n)
      logical :: lgam
      !
      ! local variables
      !
      integer :: nb,ig, iter
      real(dp):: norm_dgrad, norm_thr
      complex(dp) :: sca_c
      complex(dp) :: d1(ngw), c1(ngw), grad0(ngw), dgrad(ngw)
      !
      call start_clock('pc3nc_fixed')
      !
      norm_thr = 1.0E-12
      grad0(:) = grad(:, fixed_band)  
      ! 
      ! compute c1 = \sum_i |phi_i><grad_i|phi_io > 
      !
      c1(:) = cmplx(0.0d0, 0.0d0)
      do nb=1, n
         ! 
         sca_c=cmplx(0.0d0, 0.0d0)
         ! 
         if (ispin(nb) == ispin(fixed_band)) then
            !  
            if (lgam) then
               !
               if (ng0.eq.2) grad(1,nb) = cmplx(dble(grad(1,nb)),0.0d0)
               !
               do ig=1,ngw  !loop on g vectors
                  sca_c = sca_c + conjg(grad(ig,nb))*phi(ig,fixed_band)
               enddo
               !  
               sca_c = sca_c*2.0d0  !2. for real weavefunctions
               !
               if (ng0.eq.2) then
                   sca_c = cmplx(dble(sca_c), 0.d0) &
                         - cmplx(dble(conjg(grad(1,nb))*phi(1,fixed_band)),0.d0)
               else 
                   sca_c = cmplx(dble(sca_c), 0.d0)
               endif 
               !  
            else
               !
               do ig=1,ngw           !loop on g vectors
                  ! 
                  sca_c = sca_c + conjg(grad(ig,nb)) * phi(ig,fixed_band)
                  !  
               enddo
               !
            endif
            !
         endif
         !
         call mp_sum( sca_c, intra_image_comm )
         !
         c1(:) = c1(:) + sca_c*phi(:, nb) 
         ! 
         if (lgam) then
            ! 
            if (ng0.eq.2) c1(1) = cmplx(dble(c1(1)),0.0d0)
            !
         endif
         !
      enddo
      !
      ! main iteration
      !
      do iter = 1, 100
         !
         ! compute d1 = [\sum_i |phi_i><\phi_i|grad0> ]
         !
         d1(:)=cmplx(0.0d0, 0.0d0)
         do nb=1, n
            ! 
            sca_c=cmplx(0.0d0, 0.0d0)
            ! 
            if (ispin(nb) == ispin(fixed_band)) then
               !  
               if (lgam) then
                  !
                  if (ng0.eq.2) grad0(1) = cmplx(dble(grad0(1)),0.0d0)
                  !
                  do ig=1,ngw  !loop on g vectors
                     sca_c = sca_c + conjg(phi(ig,nb))*grad0(ig)
                  enddo
                  !  
                  sca_c = sca_c*2.0d0  !2. for real weavefunctions
                  !
                  if (ng0.eq.2) then
                      sca_c = cmplx(dble(sca_c), 0.d0) &
                            - cmplx(dble(conjg(phi(1,nb))* grad0(1)),0.d0)
                  else
                      sca_c = cmplx(dble(sca_c), 0.d0)
                  endif
                  !  
               else
                  !
                  do ig=1,ngw           !loop on g vectors
                     ! 
                     sca_c = sca_c + conjg(phi(ig,nb))*grad0(ig)
                     !  
                  enddo
                  !
               endif
               !
            endif
            !
            call mp_sum( sca_c, intra_image_comm )
            !
            d1(:) = d1(:) + sca_c*phi(:, nb)
            ! 
            if (lgam) then
               ! 
               if (ng0.eq.2) d1(1) = cmplx(dble(d1(1)),0.0d0)
               !
            endif
            !
         enddo
         !
         ! compute the gradient
         !
         dgrad(:) = grad0(:) - (1.d0/2.d0)*(d1(:) + c1(:))
         !
         ! compute norm |dgrad| 
         !
         norm_dgrad = 0.0d0
         do ig=1,ngw 
           norm_dgrad = norm_dgrad + conjg(dgrad(ig))*dgrad(ig)
         enddo
         !
         if (lgam) then
            norm_dgrad = 2.0d0*norm_dgrad 
            if (ng0.eq.2) norm_dgrad = norm_dgrad - conjg(dgrad(1))*dgrad(1)
         endif
         !
         call mp_sum (norm_dgrad, intra_image_comm)
         !
         write(stdout,*) 'norm dgrad', iter, norm_dgrad
         !
         if (norm_dgrad <= norm_thr) then
            !
            goto 100
            ! 
         else
            !
            grad0(:) = grad0(:) - dgrad(:)
            ! 
         endif
         !  
      enddo
      !
      100 continue
      !
      grad(:, fixed_band) = grad0(:) 
      !  
      return
      ! 
      end subroutine pc3nc_fixed
