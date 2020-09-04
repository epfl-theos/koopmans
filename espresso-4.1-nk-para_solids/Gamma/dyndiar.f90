!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dyndiar (dyn,nat3,nmodes,u,nat,ityp,amass,w2,dynout)
  !-----------------------------------------------------------------------
  !
  !   diagonalizes the dynamical matrix "dyn", returns energies in "w2"
  !   and mode displacements in "dynout". dyn is unchanged on output.
  !
#include "f_defs.h"
  USE kinds, only : DP
  USE io_global,  ONLY : stdout
  implicit none
  integer :: nmodes, nat3, nat,ityp(nat), iudyn
  real(DP):: dyn(nat3,nmodes), u(nat3,nmodes), amass(*)
  real(DP):: dynout(nat3,nmodes), w2(nat3)
  !
  integer:: nu_i, nu_j, mu, na, nb, nt, i, j
  real(DP), allocatable :: m(:,:), z(:,:)
  real(DP) :: rydthz, rydcm1, w1, unorm, sum, dif
  !
  allocate  ( m  ( nat3, nat3))    
  allocate  ( z  ( nat3, nat3))    
  !
  call DCOPY(nat3*nmodes,dyn,1,dynout,1)
  !
  !  Impose symmetry to the matrix
  !
  dif=0.d0
  do nu_i=1,nmodes
     do nu_j=1,nu_i-1
        dif = dif + abs(dynout(nu_i,nu_j)-dynout(nu_j,nu_i))
        dynout(nu_j,nu_i) = 0.5d0*(dynout(nu_i,nu_j)+dynout(nu_j,nu_i))
        dynout(nu_i,nu_j) = dynout(nu_j,nu_i)
     end do
  end do
  WRITE( stdout,9000) dif
  !
  !  Impose Acoustic Sum Rule
  !
  dif=0.d0
  do i=1,3
     do j=1,3
        do na=1,nat
           sum=0.d0
           do nb=1,nat
              if (na.ne.nb) sum=sum+dynout((na-1)*3+i,(nb-1)*3+j)
           end do
           dif = dif + abs(dynout((na-1)*3+i,(na-1)*3+j) + sum)
           dynout((na-1)*3+i,(na-1)*3+j) = -sum
        end do
     end do
  end do
  WRITE( stdout,9005) dif
  !
  !  fill the mass matrix
  !
  do nu_i = 1,nmodes
     do nu_j = 1,nmodes
        m(nu_i,nu_j) = 0.0d0
        do mu = 1,3*nat
           na = (mu-1)/3+1
           nt = ityp(na)
           m(nu_i,nu_j) = m(nu_i,nu_j) + amass(nt)*u(mu,nu_i)*u(mu,nu_j)
        end do
     end do
  end do
  !
  !  solve the generalized eigenvalue problem w2*(M*z) = (Cz)
  !  Note that z are eigendisplacements in the base of input
  !  modes u and that they are normalized as <z|M|z>=I
  !
  call rdiaghg (nat3, nmodes, dynout, m, nat3, w2, z)
  !
  !  conversion factors ryd=>thz e ryd=>1/cm
  !
  rydthz = 13.6058d0*241.796d0
  rydcm1 = 13.6058d0*8065.5d0
  !
  !  write frequencies
  !
  WRITE( stdout,'(5x,"diagonalizing the dynamical matrix ..."//)')
  WRITE( stdout,'(1x,74("*"))')
  !
  dynout (:,:) = 0.0d0
  do nu_i = 1,nmodes
     w1 = sqrt(abs(w2(nu_i)))
     if (w2(nu_i).lt.0.0) w1 = -w1
     WRITE( stdout,9010) nu_i, w1*rydthz, w1*rydcm1
     !  bring eigendisplacements in cartesian axis
     do mu = 1,3*nat
        do i = 1,nmodes
           dynout(mu,nu_i) = dynout(mu,nu_i) + z(i,nu_i)*u(mu,i)
        end do
     end do
  end do
  WRITE( stdout,'(1x,74("*"))')
  !
  deallocate(z)
  deallocate(m)
  return
  !
9000 format ('  Symmetry violation  sum_ij |D_ij-D_ji| :',f15.6)
9005 format ('  ASR violation  sum_i |D_ij| :',f15.6)
9010 format(5x,'omega(',i3,') =',f10.6,' [THz] =',f11.6,' [cm-1]')
  !
end subroutine dyndiar
