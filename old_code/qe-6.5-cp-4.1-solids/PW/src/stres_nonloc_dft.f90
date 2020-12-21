!
! Copyright (C) 2010- Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
subroutine stres_nonloc_dft( rho, rho_core, nspin, sigma_nonloc_dft )

  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  use funct,            ONLY : get_igcc, get_inlc 
  USE mp,               ONLY : mp_sum
  USE fft_base,         ONLY : dfftp
  USE vdW_DF,           ONLY : vdW_DF_stress
  USE rVV10,            ONLY : stress_rVV10
  USE io_global,        ONLY : stdout
  !
  IMPLICIT NONE
  !
  integer, intent(in) ::nspin
  real(DP), intent(in)     :: rho (dfftp%nnr), rho_core (dfftp%nnr)
  real(DP), intent(inout)  :: sigma_nonloc_dft (3, 3)

  integer :: l, m, inlc


  sigma_nonloc_dft(:,:) = 0.d0
  inlc = get_inlc()

  if ( inlc==1 .or. inlc==2 ) then
     CALL vdW_DF_stress(rho, rho_core, nspin, sigma_nonloc_dft)
  elseif ( inlc == 3 ) then
     CALL stress_rVV10(rho, rho_core, nspin, sigma_nonloc_dft)
  end if

  return

end subroutine stres_nonloc_dft

