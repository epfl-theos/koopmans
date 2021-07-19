!
! Copyright (C) 2002-2007 Quantum-Espresso groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!----------------------------------------------------------------------
      subroutine vol_clu(rho_real,rho_g,s_fac,flag)
!----------------------------------------------------------------------
! it computes the volume of the cluster (cluster calculations) starting
! from the measure of the region of space occupied by the electronic density
! above a given threshold

      use cell_base
      use electrons_base
      use ions_base
      use ions_positions, only: tau0
      use constants, only: pi
      use parameters
      use reciprocal_vectors
      use gvecs
      use gvecp, only: ngm
      use recvecs_indexes
      use cp_main_variables, only: drhor
      use control_flags, only: tpre
      use local_pseudo
      USE cp_interfaces,    ONLY: fwfft, invfft
      use grid_dimensions, only: nr1, nr2, nr3,                         &
     &                   nr1x, nr2x, nr3x, nnr => nnrx
      use pres_ai_mod, only: rho_thr, n_cntr, cntr, step_rad, fill_vac, &
     &                       delta_eps, delta_sigma, axis,              &
     &                       abisur, dthr, Surf_t, rho_gaus, v_vol,     &
     &                       posv, xc0, weight, volclu, stress_vol,     &
     &                       surfclu, n_ele, jellium, R_j, h_j, e_j,    &
     &                       nelect, P_ext
      use fft_base
#ifdef __PARA
      use mp_global, only: nproc, mpime
      use io_global, only: ionode
      USE mp,                 ONLY: mp_bcast, mp_sum
      USE mp_global,          ONLY: intra_image_comm
#endif

      implicit none

#ifdef __PARA
      include 'mpif.h'
#endif

      real(kind=8) dx, dxx, xcc(4500)
      real(kind=8) weight0, wpiu, wmeno, maxr, minr
      real(kind=8) tau00(3), dist
      real(kind=8) rho_real(nnr,nspin), rhoc
      real(kind=8) alfa(nsx), alfa0, sigma, hgt 
      real(kind=8) pos_cry(3), pos_car(3), pos_aux(3)
      real(kind=8) pos_cry0(3), dpvdh(3,3)
      real(kind=8) v_d(3)
      real(kind=8) mtot, rad0, cm(3)
      real(kind=8) modr, lap
      real(kind=8) prod, aux1
      real(kind=8) gxl, xyr, xzr, yzr
      real(kind=8), allocatable:: vec(:,:,:), aiuto(:,:,:)
      real(kind=8), allocatable:: drho(:,:), d2rho(:,:)
      real(kind=8), allocatable:: dxdyrho(:), dxdzrho(:)
      real(kind=8), allocatable:: dydzrho(:)
      real(kind=8), allocatable:: tauv(:,:,:)

      complex(kind=8) s_fac(ngs,nsp), ci
      complex(kind=8) sum_sf, aux, auxx, fact, rho_g(ngm,nspin) 
      complex(kind=8), allocatable :: psi(:), rhofill(:), rhotmp(:,:)

      integer ir, ir1, ir2, ir3, is, iss, ia, flag, ierr
      integer i, j, k, l, ig, cnt, nmin, nmax, n_at

#ifdef __PARA
      real(kind=8) maxr_p(nproc), minr_p(nproc), maxr_pp, minr_pp
      integer shift(nproc), incr(nproc),  ppp(nproc) 
      integer displs(nproc), ip, me
#endif
      if (abisur) allocate(drho(3,nnr))
      if (abisur) allocate(d2rho(3,nnr))
      if (abisur) allocate(dxdyrho(nnr))
      if (abisur) allocate(dxdzrho(nnr))
      if (abisur) allocate(dydzrho(nnr))
      allocate(psi(nnr))

      call start_clock( 'vol_clu' )

      ci = (0.d0,1.d0)

#ifdef __PARA
      me = mpime + 1
      do ip=1,nproc
         ppp(ip) =  dfftp%nnp  * ( dfftp%npp(ip) )
         if (ip.eq.1) then
            shift(ip)=0
         else
            shift(ip)=shift(ip-1) + ppp(ip-1)
         end if
      end do
#endif

      sigma = rho_thr/3.d0 !3.d0
      hgt = 0.0050d0 !5000.d0*rho_thr
! We smear the step function defining the volume and approximate its derivative
! with a gaussian. Here we sample the integral of this gaussian. It has to 
! be done once for ever
      dx = 5.d0*sigma/60.d0
      if (flag.eq.1) then
         dxx = dx/40.d0
         weight(1) = 0.d0
         xcc(1) = rho_thr - 5.d0*sigma
         xc0(1) = xcc(1)
         cnt = 1
         do i = 2,121
            weight(i) = weight(i-1)
            do j = 1,40
               cnt = cnt + 1
               xcc(cnt) = xcc(cnt-1) + dxx
               if (j.eq.40) then
                  xc0(i) = xcc(cnt)
               end if
               aux1 = xcc(cnt)-dxx/2.d0-rho_thr
               weight(i) = weight(i) + 1.d0/(sigma*dsqrt(pi*2.d0)) *    &
     &                     dxx * dexp(-1.d0*aux1**2/(2.d0*sigma**2))
            end do
         end do
! This doesn't work yet.....
         if (jellium) then
            do ir3 = 1,nr3
               do ir2 = 1,nr2
                  do ir1 = 1,nr1
                     ir = ir1 + (ir2-1)*nr1 + (ir3-1)*nr2*nr1
                     dist = 0.d0
                     do i = 1,3
                        posv(i,ir) = (DBLE(ir1)-1.0d0)*a1(i)/DBLE(nr1) +&
     &                              (DBLE(ir2)-1.0d0)*a2(i)/DBLE(nr2) +&
     &                              (DBLE(ir3)-1.0d0)*a3(i)/DBLE(nr3)
                     end do
                  end do
               end do
            end do
         end if
      end if

      n_at = MAXVAL ( na(1:nsp) )
      allocate ( tauv(3,n_at,nsp) )
      n_at = 0
      do is = 1,nsp
         alfa(is) = step_rad(is)/2.d0
         do ia = 1,na(is)
            n_at = n_at + 1
            do k = 1,3
               tauv(k,ia,is) = tau0(k,n_at)
            end do
         end do
      end do

      stress_vol = 0.d0
      dpvdh = 0.d0
      
! Now we compute the volume and other quantities

      volclu = 0.d0
      n_ele = 0.d0
      surfclu = 0.d0

! Let's add rhops to fill possible holes in the valence charge density on top
! of the ions

      allocate(rhotmp(ngm,nspin))
      rhotmp = CMPLX(0.d0,0.d0)

      if (nspin.eq.1) then
         do ig = 1,ngm
            rhotmp(ig,1)=rho_g(ig,1)
         end do
      else
         do ig = 1,ngm
            do iss = 1,2
               rhotmp(ig,iss) = rho_g(ig,iss) 
            end do
         end do
      end if

! To fill the vacuum inside hollow structures

      if (fill_vac) then
         allocate(rhofill(ngm))
         rhofill = 0.d0
         do k = 1,3
            cm(k) = 0.d0
            mtot = 0.d0
            do is = 1,nsp
               do ia = 1,na(is)
                  cm(k) = cm(k) + tauv(k,ia,is)*pmass(is)
               end do
               mtot = mtot + pmass(is)
            end do
            cm(k) = cm(k)/mtot
         end do
      end if

      if (fill_vac) then
         do i = 1,n_cntr
            do is = 1,nsp
               if (cntr(is)) then
                  rad0 = step_rad(is) + DBLE(i)*delta_sigma
                  alfa0 = rad0/2.d0
                  do ia = 1,na(is)
                     do k = 1,3
                        if (k.ne.axis) then
                           tau00(k) = (tauv(k,ia,is)-cm(k))*              &
     &                                  (1.d0-delta_eps*DBLE(i))+cm(k)
                        else
                           tau00(k) = tauv(k,ia,is)
                        end if
                     end do
                     do ig = 1,ngm
                        prod = 0.d0
                        do k = 1,3
                           prod = prod + gx(k,ig)*tau00(k)
                        end do
                        prod = prod*tpiba
                        fact = CMPLX(dcos(prod),-1.d0*dsin(prod))
                        aux = alfa0*hgt*dexp(-0.50d0*alfa0**2*g(ig)*tpiba2)
                        rhofill(ig) = rhofill(ig) + aux*fact
                     end do 
                  end do
               end if
            end do
         end do   
         if (nspin.eq.1) then
            do ig=1,ngm
               rhotmp(ig,1) = rhotmp(ig,1) + rhofill(ig)
            end do
         else
            do ig = 1,ngm
               do iss = 1,2
                  rhotmp(ig,iss) = rhotmp(ig,iss) + 0.5d0*rhofill(ig)
               end do
            end do
         end if
      end if

      if (fill_vac) then
         deallocate(rhofill)
      end if

      if (abisur)                                                       &
     &   call gradrho(nspin,rhotmp,drho,d2rho,dxdyrho,dxdzrho,dydzrho)

      psi = CMPLX(0.d0,0.d0)
      if (nspin.eq.1) then
         do ig = 1,ngm
            psi(np(ig)) = rhotmp(ig,1)
            psi(nm(ig)) = conjg(rhotmp(ig,1))
         end do
         call invfft('Dense',psi, dfftp )
         do ir = 1,nnr
            rho_gaus(ir) = real(psi(ir))
         end do
      else            
        do ig = 1,ngm
            psi(np(ig)) = rhotmp(ig,1) + ci*rhotmp(ig,2)
            psi(nm(ig)) = conjg(rhotmp(ig,1)) + ci*conjg(rhotmp(ig,2))
         end do
         call invfft('Dense',psi, dfftp )
         do ir = 1,nnr
            rho_gaus(ir) = real(psi(ir))+aimag(psi(ir))
         end do
      end if
      deallocate(psi)
      deallocate(rhotmp)

      e_j = 0.d0

      do ir = 1,nnr
   
         v_vol(ir) = 0.d0

         if (jellium) then
#ifdef __PARA
            do j = 1,3
               pos_aux(j) = posv(j,ir+shift(me))
            end do
#else
            do j = 1,3
               pos_aux(j) = posv(j,ir)
            end do
#endif
            dist = 0.d0
            do j = 1,3
               dist = dist + (pos_aux(j) - 0.5d0*(a1(j)+a2(j)+a3(j)))**2
            end do
            dist = dsqrt(dist)
            if (dist.ge.R_j) then
               v_vol(ir) = - nelect/dist
               v_vol(ir) = 0.d0
            else
! The last term in the internal potential is for its continuity
               v_vol(ir) = + 0.5d0*nelect*dist**2/R_j**3                &
                           - 1.5d0*nelect/R_j
               v_vol(ir) = - h_j
            end if
            if (nspin.eq.1) then
               e_j = e_j + v_vol(ir) * rho_real(ir,1) * omega /         &
     &                                DBLE(nr1*nr2*nr3)
            else
               e_j = e_j + v_vol(ir) *                                  &
                     ( rho_real(ir,1) + rho_real(ir,2) ) * omega /      &
     &                                DBLE(nr1*nr2*nr3)
            end if
         end if

         rhoc = rho_gaus(ir)
! Volume and surface
         if (rhoc.gt.rho_thr+5.d0*sigma) then
            weight0 = 1.d0
            wpiu = 1.d0
            i = int((rhoc-rho_thr-dthr+5.d0*sigma)/dx) + 1
            if (i.gt.120) then
               wmeno = 1.d0
            else
               wmeno = weight(i) + (weight(i+1)-weight(i)) *            &
     &                 (rhoc-rho_thr-dthr-DBLE(i-1)*dx+5.d0*sigma)/dx
            end if
            go to 79
         end if
! Volume and surface
         k = int((rhoc-rho_thr+5.d0*sigma)/dx) + 1
         weight0 = weight(k) + (weight(k+1)-weight(k)) *                &
                   (rhoc-rho_thr+5.d0*sigma-DBLE(k-1)*dx)/dx
         if (abisur) then
            if (rhoc-rho_thr+dthr.gt.5.d0*sigma) then
               wpiu = weight0
               i = int((rhoc-rho_thr-dthr+5.d0*sigma)/dx) + 1
               wmeno = weight(i)+(weight(i+1)-weight(i))*               &
     &            (rhoc-rho_thr-dthr+5.d0*sigma-DBLE(i-1)*dx)/dx
            else if (rho_thr+dthr-rhoc.gt.5.d0*sigma) then
               wmeno = 0.d0
               i = int((rhoc-rho_thr+dthr+5.d0*sigma)/dx) + 1
               wpiu = weight0
            else
               i = int((rhoc-rho_thr+dthr+5.d0*sigma)/dx) + 1
               wpiu = weight0
               i = int((rhoc-rho_thr-dthr+5.d0*sigma)/dx) + 1
               wmeno = weight(i)+(weight(i+1)-weight(i))*               &
     &            (rhoc-rho_thr-dthr+5.d0*sigma-DBLE(i-1)*dx)/dx
            end if
         end if
  79     continue
         if (nspin.eq.1) then
            n_ele = n_ele + weight0 * rho_real(ir,1)
         else
            n_ele = n_ele + weight0 * (rho_real(ir,1) + rho_real(ir,2))
         end if
         volclu = volclu + weight0
         v_vol(ir) = v_vol(ir) + P_ext /(sigma*dsqrt(pi*2.d0)) *     &
     &               dexp(-1.d0*(rhoc-rho_thr)**2/(2.d0*sigma**2))
         if (tpre) then
            do k = 1,3
               do j = 1,3
                  do is = 1,nspin
                     dpvdh(k,j) = dpvdh(k,j) +                       &
     &                    v_vol(ir)*drhor(ir,is,k,j)*omega/          &
     &                    DBLE(nr1*nr2*nr3)
                  end do
               end do
            end do
         end if

         if (abisur) then
            modr = 0.d0
            lap = 0.d0
            gxl = 0.d0
            do j = 1,3
               modr = modr + drho(j,ir)**2
               lap = lap + d2rho(j,ir)
               gxl = gxl + drho(j,ir)**2*d2rho(j,ir)
            end do
            xyr = 2.d0*dxdyrho(ir)*drho(1,ir)*drho(2,ir)
            xzr = 2.d0*dxdzrho(ir)*drho(1,ir)*drho(3,ir)
            yzr = 2.d0*dydzrho(ir)*drho(2,ir)*drho(3,ir)
            modr = dsqrt(modr)
            surfclu = surfclu + (wpiu-wmeno)*modr
            v_vol(ir) = v_vol(ir) -1.d0*Surf_t/dthr * (wpiu-wmeno) *    &
     &                  (lap/modr - (gxl + xyr + xzr + yzr)/modr**3)
         end if

      end do

#ifdef __PARA
      call mp_sum(volclu,intra_image_comm)
      call mp_sum(n_ele,intra_image_comm)
      if (jellium) call mp_sum(e_j,intra_image_comm)
      call mp_sum(surfclu,intra_image_comm)
      call mp_sum(dpvdh,intra_image_comm)
#endif
      volclu = volclu * omega / DBLE(nr1*nr2*nr3)
      n_ele = n_ele * omega / DBLE(nr1*nr2*nr3)
      surfclu = surfclu * omega / DBLE(nr1*nr2*nr3) / dthr
      do i = 1,3
         do j = 1,3
            stress_vol(i,j) =  dpvdh(i,1)*h(j,1) + dpvdh(i,2)*h(j,2) +  &
     &                         dpvdh(i,3)*h(j,3)
         end do
      end do
 
      deallocate( tauv )
      if ( abisur ) deallocate( drho )
      if ( abisur ) deallocate( d2rho )
      if ( abisur ) deallocate( dxdyrho )
      if ( abisur ) deallocate( dxdzrho )
      if ( abisur ) deallocate( dydzrho )

      call stop_clock( 'vol_clu' )

      return
      end

