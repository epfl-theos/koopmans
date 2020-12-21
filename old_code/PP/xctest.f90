!
! Copyright (C) 2005 Quantum-espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
program xctest
  USE mp, ONLY: mp_start, mp_end
  use kinds, only: DP
  use funct, only: set_dft_from_indices
  implicit none
  integer :: nnr = 1000
  integer :: nspin = 2
  real(DP), allocatable :: rhor( :, : )
  real(DP), allocatable :: grhor( :, :, : )
  integer iexch,icorr,igcx,igcc
  
  CALL mp_start()

  iexch=1
  icorr=3
  igcx=1
  igcc=3
  call set_dft_from_indices(iexch,icorr,igcx,igcc)

  open(unit=17,form='unformatted',status='old')
  read(17) nnr, nspin
  allocate(rhor( nnr, nspin ))
  allocate(grhor( nnr, 3, nspin ))
  read(17) rhor
  read(17) grhor
  close(17)

  !CALL test_gcxc( nnr, nspin, rhor, grhor )
  CALL test_xc( nnr, nspin, rhor, grhor )
  
  CALL mp_end()
end program xctest

subroutine test_gcxc( nnr, nspin, rhor, grhor )
  use kinds, only: DP
!  use funct, only: gcxc
  implicit none
  integer, intent(in) :: nnr, nspin
  real(DP) :: rhor( nnr, nspin )
  real(DP) :: grhor( nnr, 3, nspin )
  !
  real(DP), parameter :: epsr = 1.0d-10, epsg = 1.0d-10
  real(DP), parameter :: e2   = 1.0d0
  real(DP) :: grho2( nspin )
  real(DP) :: arho, segno
  real(DP) :: sx_w, sc_w, v1x_w, v2x_w, v1c_w, v2c_w
  real(DP) :: sx, sc, v1x, v2x, v1c, v2c
  real(DP) :: sx_m, sc_m, v1x_m, v2x_m, v1c_m, v2c_m
  real(DP) :: sx_d, sc_d, v1x_d, v2x_d, v1c_d, v2c_d
  integer :: k, is, ipol

    do k = 1, nnr
       !
       !
       do is = 1, nspin
          grho2 (is) = grhor(k, 1, is)**2 + grhor(k, 2, is)**2 + grhor(k, 3, is)**2
       enddo
       !
       !
       if (nspin == 1) then
          !
          !    This is the spin-unpolarised case
          !
          arho = abs (rhor (k, 1) )
          segno = sign (1.d0, rhor (k, 1) )
          if (arho > epsr .and. grho2 (1) > epsg) then

             ! call gcxc (arho, grho2(1), sx, sc, v1x, v2x, v1c, v2c)

             call becke88 (arho, grho2(1), sx, v1x, v2x)
             call wrap_b88 (arho, grho2(1), sx_w, v1x_w, v2x_w)  ! DEBUG
             call glyp (arho, grho2(1), sc, v1c, v2c)
             call wrap_glyp (arho, grho2(1), sc_w, v1c_w, v2c_w)  ! DEBUG

             sx_d = (sx_w - sx) / (abs(sx) + abs(sx_w))
             sc_d = (sc_w - sc) / (abs(sc) + abs(sc_w))
             v1x_d = (v1x_w - v1x) / (abs(v1x) + abs(v1x_w))
             v1c_d = (v1c_w - v1c) / (abs(v1c) + abs(v1c_w))
             v2x_d = (v2x_w - v2x) / (abs(v2x) + abs(v2x_w))
             v2c_d = (v2c_w - v2c) / (abs(v2c) + abs(v2c_w))
             

             write(18,*) arho,grho2(1), sx_d, sc_d
             write(19,*) arho,grho2(1), v1x_d, v1c_d
             write(20,*) arho,grho2(1), v2x_w, v2x, v2x_d
             write(21,*) arho,grho2(1), v2c_w, v2c, v2c_d

             !
             ! first term of the gradient correction : D(rho*Exc)/D(rho)

             ! v (k, 1) = v (k, 1) + e2 * (v1x + v1c)

             ! HERE h contains D(rho*Exc)/D(|grad rho|) / |grad rho|
             !
             ! h (k, 1, 1) = e2 * (v2x + v2c)
             ! etxc = etxc + e2 * (sx + sc) * segno

          else
             ! h (k, 1, 1) = 0.d0
             sx = 0.0d0
             sc = 0.0d0
          endif
          !
       endif
       !
    end do

    return
end subroutine test_gcxc

!
!
!

subroutine test_xc( nnr, nspin, rhor, grhor )
  use kinds, only: DP
  use funct, only: get_iexch, get_icorr, get_igcx, get_igcc

  implicit none
  integer, intent(in) :: nnr, nspin
  real(DP) :: rhor( nnr, nspin )
  real(DP) :: grhor( nnr, 3, nspin )
  !
  real(DP) :: rhon( nnr, nspin )
  real(DP) :: grhon( nnr, 3, nspin )
  real(DP) :: exc, excn, rhod, grhod
  integer :: ir, is, ipol
  integer iexch,icorr,igcx,igcc


  iexch = get_iexch()
  icorr = get_icorr()
  igcx  = get_igcx()
  igcc  = get_igcc()

  rhon  = rhor
  grhon = grhor
  !
  ! original CP xc selection
  !
      if (iexch==1.and.icorr==1.and.igcx==0.and.igcc==0) then
         ! LDA (Perdew-Zunger)
         call expxc(nnr,nspin,rhor,exc)
      else if (iexch==1.and.icorr==4.and.igcx==2.and.igcc==2) then
         ! PW91
         call ggapwold(nnr,nspin,grhor,rhor,exc)
      else if (iexch==1.and.icorr==3.and.igcx==1.and.igcc==3) then
         ! BLYP
         call ggablyp4(nnr,nspin,grhor,rhor,exc)
      else if (iexch==1.and.icorr==4.and.igcx==3.and.igcc==4) then
         ! PBE
         call ggapbe(nnr,nspin,grhor,rhor,exc)
      else
         call errore('exc-cor','no such exch-corr',1)
      end if
  !
  ! Wrapper to PW xc selection
  !
  call exch_corr_cp(nnr,nspin,grhon,rhon,excn)
  !
  write(6,*) 'EXC = ', exc, excn
  do is = 1, nspin
    do ir = 1, nnr
      rhod = abs( rhor( ir, is ) - rhon( ir, is ) ) / ( abs( rhor( ir, is ) ) + abs( rhon( ir, is ) ) )
      WRITE(18,100) ir,is,rhod
    end do
  end do
  do is = 1, nspin
    do ir = 1, nnr
      do ipol = 1, 3
      grhod = abs( grhor( ir, ipol, is ) - grhon( ir, ipol, is ) ) / &
             ( abs( grhor( ir, ipol, is ) ) + abs( grhon( ir, ipol, is ) ) )
      WRITE(19,100) ir,is,grhod
      end do
    end do
  end do
100 FORMAT( I5, I2, 1X, E15.8, 1X, E15.8 )
end subroutine test_xc
