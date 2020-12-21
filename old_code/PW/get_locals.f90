!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
      subroutine get_locals(rholoc,magloc, rho)
!---------------------------------------------------------------------------
!
! Here local integrations are carried out around atoms.
! The points and weights for these integrations are determined in the
! subroutine make_pointlists, the result may be printed in the
! subroutine report_mag. If constraints are present, the results of this
! calculation are used in v_of_rho for determining the penalty functional.
!
      USE kinds,      ONLY : DP
      USE ions_base,  ONLY : nat
      USE cell_base,  ONLY : omega
      USE gvect,      ONLY : nr1, nr2, nr3, nrxx
      USE lsda_mod,   ONLY : nspin
      USE mp_global,  ONLY : intra_pool_comm
      USE mp,         ONLY : mp_sum
      use noncollin_module

      implicit none
!
! I/O variables
!
      real(DP) ::   &
          rholoc(nat),   &     ! integrated charge arount the atoms
          magloc(nspin-1,nat)  ! integrated magnetic moment around the atom
      real(DP) :: rho (nrxx, nspin)
!
! local variables
!
      integer i,ipol
      real(DP) :: fact
      real(DP), allocatable :: auxrholoc(:,:)
 
      allocate (auxrholoc(0:nat,nspin))
      auxrholoc(:,:) = 0.d0
      do i=1,nrxx
         auxrholoc(pointlist(i),1:nspin) = auxrholoc(pointlist(i),1:nspin) + &
                                           rho(i,1:nspin) * factlist(i)
      end do
      call mp_sum( auxrholoc( 0:nat, 1:nspin ), intra_pool_comm )
!
      fact =  omega/(nr1*nr2*nr3)
      if (nspin.eq.2) then
         rholoc(1:nat)   = (auxrholoc(1:nat,1)+auxrholoc(1:nat,2)) * fact
         magloc(1,1:nat) = (auxrholoc(1:nat,1)-auxrholoc(1:nat,2)) * fact
      else
         rholoc(1:nat) = auxrholoc(1:nat,1) * fact
         if (noncolin) then
            do ipol=1,3
               magloc(ipol,1:nat) = auxrholoc(1:nat,ipol+1) * fact
            end do
         end if
      endif
!
      deallocate (auxrholoc)

      end subroutine get_locals
