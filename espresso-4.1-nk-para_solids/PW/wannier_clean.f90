! Copyright (C) 2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)

#include "f_defs.h"

!----------------------------------------------------------------------
subroutine wannier_clean()
  !----------------------------------------------------------------------
  !    
  ! ... This routine deallocates all dynamically allocated arrays for wannier calc and closes all open files
  !
  USE wannier_new, only: wan_in, wan_pot, wannier_energy, wannier_occ, pp, coef
  USE io_files
  USE buffers
  USE ldaU,       ONLY : swfcatom
  
  IMPLICIT NONE
  LOGICAL :: opnd
  
  if(allocated(wan_in)) deallocate(wan_in)
  if(allocated(wan_pot)) deallocate(wan_pot)
  if(allocated(wannier_energy)) deallocate(wannier_energy)
  if(allocated(wannier_occ)) deallocate(wannier_occ)
  if(allocated(pp)) deallocate(pp)
  if(allocated(coef)) deallocate(coef)

  CALL close_buffer( iunwpp, 'keep' )
  CALL close_buffer( iunwf, 'keep' )
  
  INQUIRE( UNIT = iunat, OPENED = opnd )  
  IF ( opnd ) CALL close_buffer( iunat, 'delete' )
  INQUIRE( UNIT = iunsat, OPENED = opnd )  
  IF ( opnd ) CALL close_buffer( iunsat, 'delete' )
  INQUIRE( UNIT = iunigk, OPENED = opnd )  
  IF ( opnd ) CALL close_buffer( iunigk, 'delete' )
  
  IF(ALLOCATED(swfcatom)) DEALLOCATE(swfcatom)

  return
  !
end subroutine wannier_clean

