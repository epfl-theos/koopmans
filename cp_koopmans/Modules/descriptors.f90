!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

   MODULE descriptors
      !
      IMPLICIT NONE
      SAVE

      INTEGER  ldim_block, ldim_cyclic, ldim_block_cyclic, ldim_block_sca
      INTEGER  lind_block, lind_cyclic, lind_block_cyclic, lind_block_sca
      INTEGER  gind_block, gind_cyclic, gind_block_cyclic, gind_block_sca
      EXTERNAL ldim_block, ldim_cyclic, ldim_block_cyclic, ldim_block_sca
      EXTERNAL lind_block, lind_cyclic, lind_block_cyclic, lind_block_sca
      EXTERNAL gind_block, gind_cyclic, gind_block_cyclic, gind_block_sca

      !  Descriptor for Cannon's algorithm
      !
      !  Parameters to define and manage the Descriptor
      !  of square matricxes block distributed on a square grid of processors
      !  to be used with Cannon's algorithm for matrix multiplication
      !
      INTEGER, PARAMETER :: descla_siz_ = 16
      INTEGER, PARAMETER :: ilar_        = 1
      INTEGER, PARAMETER :: nlar_        = 2
      INTEGER, PARAMETER :: ilac_        = 3
      INTEGER, PARAMETER :: nlac_        = 4
      INTEGER, PARAMETER :: nlax_        = 5
      INTEGER, PARAMETER :: lambda_node_ = 6
      INTEGER, PARAMETER :: la_n_        = 7
      INTEGER, PARAMETER :: la_nx_       = 8
      INTEGER, PARAMETER :: la_npr_      = 9
      INTEGER, PARAMETER :: la_npc_      = 10
      INTEGER, PARAMETER :: la_myr_      = 11
      INTEGER, PARAMETER :: la_myc_      = 12
      INTEGER, PARAMETER :: la_comm_     = 13
      INTEGER, PARAMETER :: la_me_       = 14
      INTEGER, PARAMETER :: la_nrl_      = 15
      INTEGER, PARAMETER :: la_nrlx_     = 16
      !
      ! desc( ilar_ )  globla index of the first row in the local block of lambda
      ! desc( nlar_ )  number of row in the local block of lambda ( the "2" accounts for spin)
      ! desc( ilac_ )  global index of the first column in the local block of lambda
      ! desc( nlac_ )  number of column in the local block of lambda
      ! desc( nlax_ )  leading dimension of the distribute lambda matrix
      ! desc( lambda_node_ )  if > 0 the proc holds a block of the lambda matrix
      ! desc( la_n_ )     global dimension of the matrix
      ! desc( la_nx_ )    global leading dimension
      ! desc( la_npr_ )   number of row processors 
      ! desc( la_npc_ )   number of column processors 
      ! desc( la_myr_ )   processor row index
      ! desc( la_myc_ )   processor column index
      ! desc( la_comm_ )  communicator
      ! desc( la_me_ ) processor index ( from 0 to desc( la_npr_ ) * desc( la_npc_ ) - 1 )
      ! desc( la_nrl_ ) number of local row, when the matrix is cyclically distributed across proc
      ! desc( la_nrlx_ ) leading dimension, when the matrix is distributed by row

       
   CONTAINS

   !------------------------------------------------------------------------
   !
   SUBROUTINE descla_local_dims( i2g, nl, n, nx, np, me )
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: i2g  !  global index of the first local element
      INTEGER, INTENT(OUT) :: nl   !  local number of elements
      INTEGER, INTENT(IN)  :: n    !  number of actual element in the global array
      INTEGER, INTENT(IN)  :: nx   !  dimension of the global array (nx>=n) to be distributed
      INTEGER, INTENT(IN)  :: np   !  number of processors
      INTEGER, INTENT(IN)  :: me   !  taskid for which i2g and nl are computed
      !
      !  note that we can distribute a global array larger than the
      !  number of actual elements. This could be required for performance
      !  reasons, and to have an equal partition of matrix having different size
      !  like matrixes of spin-up and spin-down 
      !
#if __SCALAPACK
      nl  = ldim_block_sca( nx, np, me )
      i2g = gind_block_sca( 1, nx, np, me )
#else
      nl  = ldim_block( nx, np, me )
      i2g = gind_block( 1, nx, np, me )
#endif
      ! This is to try to keep a matrix N * N into the same
      ! distribution of a matrix NX * NX, useful to have 
      ! the matrix of spin-up distributed in the same way
      ! of the matrix of spin-down
      !
      IF( i2g + nl - 1 > n ) nl = n - i2g + 1
      IF( nl < 0 ) nl = 0
      RETURN
      !
   END SUBROUTINE descla_local_dims
   !
   !
   SUBROUTINE descla_init( desc, n, nx, np, me, comm, includeme )
      !
      IMPLICIT NONE  
      INTEGER, INTENT(OUT) :: desc(:)
      INTEGER, INTENT(IN)  :: n   !  the size of this matrix
      INTEGER, INTENT(IN)  :: nx  !  the max among different matrixes sharing 
                                  !  this descriptor or the same data distribution
      INTEGER, INTENT(IN)  :: np(2), me(2), comm
      INTEGER, INTENT(IN)  :: includeme
      INTEGER  :: ir, nr, ic, nc, lnode, nlax, nrl, nrlx
      INTEGER  :: ip, npp
      
      IF( np(1) /= np(2) ) &
         CALL errore( ' descla_init ', ' only square grid of proc are allowed ', 2 )
      IF( n < 0 ) &
         CALL errore( ' descla_init ', ' dummy argument n less than 1 ', 3 )
      IF( nx < n ) &
         CALL errore( ' descla_init ', ' dummy argument nx less than n ', 4 )
      IF( np(1) < 1 ) &
         CALL errore( ' descla_init ', ' dummy argument np less than 1 ', 5 )

      ! find the block maximum dimensions

#if __SCALAPACK
      nlax = ldim_block_sca( nx, np(1), 0 )
#else
      nlax = ldim_block( nx, np(1), 0 )
      DO ip = 1, np(1) - 1
         nlax = MAX( nlax, ldim_block( nx, np(1), ip ) )
      END DO
#endif
      !
      ! find local dimensions, if appropriate
      !
      IF( includeme == 1 ) THEN
         !
         CALL descla_local_dims( ir, nr, n, nx, np(1), me(1) )
         CALL descla_local_dims( ic, nc, n, nx, np(2), me(2) )
         !
         lnode = 1
         !
      ELSE
         !
         nr = 0
         nc = 0
         !  
         ir = 0
         ic = 0
         !
         lnode = -1
         !
      END IF

      desc( ilar_ ) = ir
      desc( nlar_ ) = nr
      desc( ilac_ ) = ic
      desc( nlac_ ) = nc
      desc( nlax_ ) = nlax
      desc( lambda_node_ ) = lnode
      desc( la_n_  ) = n
      desc( la_nx_ ) = nx
      desc( la_npr_ ) = np(1)
      desc( la_npc_ ) = np(2)
      desc( la_myr_ ) = me(1)
      desc( la_myc_ ) = me(2)
      desc( la_comm_ ) = comm
      desc( la_me_ )  = desc( la_myc_ ) + desc( la_myr_ ) * desc( la_npr_ )
     
      npp = np(1) * np(2)

      !  Compute local dimension of the cyclically distributed matrix
      !
      IF( includeme == 1 ) THEN
         nrl  = ldim_cyclic( n, npp, desc( la_me_ ) )
      ELSE
         nrl = 0
      END IF
      nrlx = n / npp + 1

      desc( la_nrl_  ) = nrl
      desc( la_nrlx_ ) = nrlx

      IF( nr < 0 .OR. nc < 0 ) &
         CALL errore( ' descla_init ', ' wrong valune for computed nr and nc ', 1 )
      IF( nlax < 1 ) &
         CALL errore( ' descla_init ', ' wrong value for computed nlax ', 2 )
      IF( nlax < nr ) &
         CALL errore( ' descla_init ', ' nlax < nr ', ( nr - nlax ) )
      IF( nlax < nc ) &
         CALL errore( ' descla_init ', ' nlax < nc ', ( nc - nlax ) )
      IF( nrlx < nrl ) &
         CALL errore( ' descla_init ', ' nrlx < nrl ', ( nrl - nrlx ) )
      IF( nrl < 0 ) &
         CALL errore( ' descla_init ', ' nrl < 0 ', ABS( nrl ) )

      ! WRITE(*,*) 'me1,me2,nr,nc,ir,ic= ', me(1), me(2), nr, nc, ir, ic

      RETURN
   END SUBROUTINE descla_init


   END MODULE descriptors
