!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE wavefunctions_module
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP

     IMPLICIT NONE
     SAVE

     !
     COMPLEX(DP), ALLOCATABLE, TARGET :: &
       evc(:,:)     ! wavefunctions in the PW basis set
                    ! noncolinear case: first index
                    ! is a combined PW + spin index
     !
     COMPLEX(DP) , ALLOCATABLE, TARGET :: &
       psic(:), &      ! additional memory for FFT
       psic_nc(:,:)    ! as above for the noncolinear case
     !
     ! electronic wave functions, FPMD code
     !
     COMPLEX(DP), ALLOCATABLE :: c0(:,:)  ! wave functions at time t
     COMPLEX(DP), ALLOCATABLE :: cm(:,:)  ! wave functions at time t-delta t
     COMPLEX(DP), ALLOCATABLE :: cp(:,:)  ! wave functions at time t+delta t
     COMPLEX(DP), ALLOCATABLE :: cstart(:,:)  ! wave functions at start
     COMPLEX(DP), ALLOCATABLE :: c0_fixed(:,:) ! wave functions at start fixed
     COMPLEX(DP), ALLOCATABLE :: c0fixed_emp(:,:) ! empty wave functions at start fixed
     COMPLEX(DP), ALLOCATABLE :: c0_occ_emp_aux(:,:) ! empty wave functions to saved
     COMPLEX(DP), ALLOCATABLE :: c0fixed_aux(:,:) ! empty wave functions to saved
     COMPLEX(DP), ALLOCATABLE :: ctot_aux(:,:)

     ! below dual wavefunctions, allocated only in the non orthogonal case
     COMPLEX(DP), ALLOCATABLE :: cdual(:,:)  ! dual wave functions at time t
     COMPLEX(DP), ALLOCATABLE :: cmdual(:,:)  ! dual wave functions at time t


   CONTAINS

     SUBROUTINE deallocate_wavefunctions
       IF( ALLOCATED( c0 ) ) DEALLOCATE( c0 )
       IF( ALLOCATED( cm ) ) DEALLOCATE( cm )
       IF( ALLOCATED( cp ) ) DEALLOCATE( cp )
       IF( ALLOCATED( psic_nc ) ) DEALLOCATE( psic_nc )
       IF( ALLOCATED( psic ) ) DEALLOCATE( psic )
       IF( ALLOCATED( evc ) ) DEALLOCATE( evc )
       IF( ALLOCATED( cdual ) ) DEALLOCATE( cdual )
       IF( ALLOCATED( cmdual ) ) DEALLOCATE( cmdual )
       IF( ALLOCATED( cstart ) ) DEALLOCATE( cstart )
       IF( ALLOCATED( c0_fixed ) ) DEALLOCATE( c0_fixed )
       IF( ALLOCATED( c0_occ_emp_aux ) ) DEALLOCATE( c0_occ_emp_aux )
       IF( ALLOCATED( c0fixed_aux ) ) DEALLOCATE( c0fixed_aux )
     END SUBROUTINE deallocate_wavefunctions

!=----------------------------------------------------------------------------=!
   END MODULE wavefunctions_module
!=----------------------------------------------------------------------------=!
