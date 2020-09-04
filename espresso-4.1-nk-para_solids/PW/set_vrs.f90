!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
subroutine set_vrs (vrs, vltot, vr, kedtau, kedtaur,nrxx, nspin, doublegrid)
  !--------------------------------------------------------------------
  ! set the total local potential vrs on the smooth mesh to be used in 
  ! h_psi, adding the (spin dependent) scf (H+xc) part and the sum of 
  ! all the local pseudopotential contributions.
  !
  USE kinds
  USE funct, only : dft_is_meta
  USE gsmooth, only : nrxxs
  implicit none

  integer :: nspin, nrxx
  ! input: number of spin components: 1 if lda, 2 if lsd, 4 if noncolinear
  ! input: the fft grid dimension
  real(DP) :: vrs (nrxx, nspin), vltot (nrxx), vr (nrxx, nspin), &
              kedtau(nrxxs,nspin), kedtaur(nrxx,nspin)
  ! output: total local potential on the smooth grid
  !         vrs=vltot+vr
  ! input: the total local pseudopotential
  ! input: the scf(H+xc) part of the local potential
  logical :: doublegrid
  ! input: true if a doublegrid is used

  integer:: is

  do is = 1, nspin
     !
     ! define the total local potential (external + scf) for each spin ...
     !
     if (is > 1 .and. nspin == 4) then
        !
        ! noncolinear case: only the first component contains vltot
        !
        vrs (:, is) = vr (:, is)
     else
        vrs (:, is) = vltot (:) + vr (:, is)
     end if
     !
     ! ... and interpolate it on the smooth mesh if necessary
     !
     if (doublegrid) call interpolate (vrs (1, is), vrs (1, is), - 1)
     if (dft_is_meta()) call interpolate(kedtaur(1,is),kedtau(1,is),-1)
  enddo
  return

end subroutine set_vrs
