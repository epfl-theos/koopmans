!
! Copyright (C) 2001-2008 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  !----------------------------------------------------------------------------
  !
  ! ... This routine divides the k points across nodes, sets the variable
  ! ... nks equal to the local (on this processors) number of k-points
  ! ... (nkstot on input is the total number of k-points)
  ! ... The distributed has "granularity kunit", that is, kunit consecutive 
  ! ... points stay on the same processor. Usually kunit=1; kunit=2 is used 
  ! ... in phonon calculations, when one has interspersed k_i and k_i+q and
  ! ... it is needed that they stay on the same processor
  !
  USE io_global, only : stdout
  USE kinds,     ONLY : DP
  USE mp_global, ONLY : my_pool_id, npool, kunit
#if defined (EXX)
  USE exx,       ONLY : index_xk, index_xkq, nkqs, nqs
  USE funct,     ONLY : dft_is_hybrid
#endif
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: lsda
    ! logical for local spin density approx.
  INTEGER, INTENT(IN)  :: nkstot
    ! total number of k-points
  INTEGER, INTENT(INOUT) :: isk(nkstot)
    ! spin index of each kpoint (when lsda=.t.)
  INTEGER, INTENT(OUT)  :: nks
    ! number of k-points per pool
  REAL (DP), INTENT(INOUT) :: xk(3,nkstot), wk(nkstot)
    ! k-points
    ! k-point weights
#if defined (EXX)
  INTEGER :: ikk, iq
#endif
  !
#if defined (__PARA)
  !
  INTEGER :: ik, nbase, rest
  !
  !
  IF ( MOD( nkstot, kunit ) /= 0 ) &
     CALL errore( 'd_&_i', ' nkstot/kunit is not an integer', nkstot )
  !
  nks    = kunit * ( nkstot / kunit / npool )
  !
  IF ( nks == 0 ) CALL errore( 'd_&_i', ' some nodes have no k-points', 1 )
  !
  rest = ( nkstot - nks * npool ) / kunit
  !
  IF ( ( my_pool_id + 1 ) <= rest ) nks = nks + kunit
  !
  ! ... calculates nbase = the position in the list of the first point that
  ! ...                    belong to this npool - 1
  !
  nbase = nks * my_pool_id
  !
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit
  !
  ! ... displaces these points in the first positions of the list
  !
  IF ( nbase > 0 ) THEN
     !
     xk(:,1:nks) =  xk(:,nbase+1:nbase+nks)
     !
     wk(1:nks) = wk(nbase+1:nbase+nks)
     !
     IF ( lsda ) isk(1:nks) = isk(nbase+1:nbase+nks)
     !
#if defined (EXX)
     IF ( dft_is_hybrid() ) THEN
        index_xk(1:nkqs) = index_xk(1:nkqs) - nbase
        index_xkq(1:nks,1:nqs) = index_xkq(nbase+1:nbase+nks,1:nqs)
        !
        ! consistency check
        !
        do ik=1,nks
           do iq =1,nqs
              ikk = index_xk(index_xkq(ik,iq))
              if ( ikk < 1 .or. ikk > nks ) then
                 write (stdout,*) ik, iq, index_xkq(ik,iq), ikk
                 call errore ('d_&_i',' error in EXX indexing',1)
              end if
           end do
        end do
        
     END IF
#endif
  !
  END IF
  !
#else
  !
  nks = nkstot
  !
#endif
  !
  RETURN
  !
END SUBROUTINE divide_et_impera
