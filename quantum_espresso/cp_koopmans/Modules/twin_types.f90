!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! #include "f_defs.h"
!
!----------------------------------------------------------------------------
MODULE twin_types
  !----------------------------------------------------------------------------
  !
  USE kinds,             ONLY : DP
  USE control_flags,     ONLY : gamma_only, do_wf_cmplx
  !

  IMPLICIT NONE

  INTERFACE init_twin
    MODULE PROCEDURE init_twin_matrix, init_twin_tensor
  END INTERFACE
  
  INTERFACE set_twin
    MODULE PROCEDURE set_twin_matrix, set_twin_tensor, set_index_twin_matrix, set_index_twin_tensor
  END INTERFACE

  INTERFACE copy_twin
    MODULE PROCEDURE copy_twin_matrix, copy_twin_tensor
  END INTERFACE

  INTERFACE zaxpy_twin
    MODULE PROCEDURE zaxpy_twin_matrix, zaxpy_twin_tensor
  END INTERFACE

  INTERFACE allocate_twin
    MODULE PROCEDURE allocate_twin_matrix, allocate_twin_tensor
  END INTERFACE

  INTERFACE deallocate_twin
    MODULE PROCEDURE deallocate_twin_matrix, deallocate_twin_tensor
  END INTERFACE

  INTERFACE twin_mp_sum
    MODULE PROCEDURE tmatrix_mp_sum, ttensor_mp_sum
  END INTERFACE

  INTERFACE twin_mp_bcast
    MODULE PROCEDURE tmatrix_mp_bcast, ttensor_mp_bcast
  END INTERFACE

  TYPE :: twin_matrix

    REAL(DP), DIMENSION(:,:), POINTER :: rvec
    COMPLEX(DP), DIMENSION(:,:), POINTER :: cvec
    
    INTEGER :: xdim
    INTEGER :: ydim
    
    LOGICAL :: isalloc
    LOGICAL :: iscmplx

  END TYPE twin_matrix

  TYPE :: twin_tensor

    REAL(DP), DIMENSION(:,:,:), POINTER :: rvec
    COMPLEX(DP), DIMENSION(:,:,:), POINTER :: cvec
    
    INTEGER :: xdim
    INTEGER :: ydim
    INTEGER :: zdim
    
    LOGICAL :: isalloc
    LOGICAL :: iscmplx

  END TYPE twin_tensor
!! rvec-> rdata
!! ierr ovunque
  CONTAINS

  SUBROUTINE init_twin_matrix(tmatrix, lgam)

   type(twin_matrix)   :: tmatrix
   logical, intent(in) :: lgam

   tmatrix%isalloc=.false.
   tmatrix%iscmplx = .not.lgam
   
   return
  END SUBROUTINE init_twin_matrix

  SUBROUTINE set_twin_matrix(tmatrix, value)

   type(twin_matrix)   :: tmatrix
   COMPLEX(DP), INTENT(IN) :: value

   IF(tmatrix%iscmplx) THEN
     tmatrix%cvec=value
   ELSE
     tmatrix%rvec = DBLE(value)
   ENDIF
   
   return
  END SUBROUTINE set_twin_matrix

  SUBROUTINE set_index_twin_matrix(tmatrix,i,j, value)

   type(twin_matrix)   :: tmatrix
   COMPLEX(DP), INTENT(IN) :: value
   INTEGER, INTENT(IN) :: i,j

   IF(tmatrix%iscmplx) THEN
     tmatrix%cvec(i,j)=value
   ELSE
     tmatrix%rvec(i,j) = DBLE(value)
   ENDIF
   
   return
  END SUBROUTINE set_index_twin_matrix

  SUBROUTINE init_twin_tensor(ttensor, lgam)

   type(twin_tensor)   :: ttensor
   logical, intent(in) :: lgam

   ttensor%isalloc=.false.
   ttensor%iscmplx = .not.lgam
   
   return
  END SUBROUTINE init_twin_tensor

  SUBROUTINE set_twin_tensor(ttensor, value)

   type(twin_tensor)   :: ttensor
   COMPLEX(DP), INTENT(IN) :: value

   IF(ttensor%iscmplx) THEN
     ttensor%cvec=value
   ELSE
     ttensor%rvec = DBLE(value)
   ENDIF
   
   return
  END SUBROUTINE set_twin_tensor

  SUBROUTINE set_index_twin_tensor(ttensor,i,j,k, value)

   type(twin_tensor)   :: ttensor
   COMPLEX(DP), INTENT(IN) :: value
   INTEGER, INTENT(IN) :: i,j,k

   IF(ttensor%iscmplx) THEN
     ttensor%cvec(i,j,k)=value
   ELSE
     ttensor%rvec(i,j,k) = DBLE(value)
   ENDIF
   
   return
  END SUBROUTINE set_index_twin_tensor

  SUBROUTINE copy_twin_matrix(tmatrix1,tmatrix2)

   type(twin_matrix)   :: tmatrix1,tmatrix2

   character(len=17)   :: subname="copy_twin_tensor"

   IF(tmatrix1%xdim.ne.tmatrix2%xdim .OR. tmatrix1%ydim.ne.tmatrix2%ydim) THEN
      call errore(subname,"copying twin matrices with incompatible size", 1)
   ENDIF
!    COMPLEX(DP), INTENT(IN) :: value

   IF(tmatrix1%iscmplx) THEN
     IF(tmatrix2%iscmplx) THEN
      tmatrix1%cvec(:,:)=tmatrix2%cvec(:,:)
     ELSE
      call errore(subname,"copying real tensor into complex tensor", 1)
      tmatrix1%cvec(:,:)=tmatrix2%rvec(:,:)
     ENDIF
   ELSE
     IF(.not.tmatrix2%iscmplx) THEN
      tmatrix1%rvec(:,:)=tmatrix2%rvec(:,:)
     ELSE
      call errore(subname,"copying complex tensor into real tensor", 1)
      tmatrix1%rvec(:,:)=tmatrix2%cvec(:,:)
     ENDIF
   ENDIF
   
   return
  END SUBROUTINE copy_twin_matrix

  SUBROUTINE zaxpy_twin_matrix(tmatrix1,tmatrix2, coef)

   type(twin_matrix)   :: tmatrix1,tmatrix2
   COMPLEX(DP)    :: coef
   character(len=17)   :: subname="copy_twin_tensor"
!    COMPLEX(DP), INTENT(IN) :: value

   IF(tmatrix1%xdim.ne.tmatrix2%xdim .OR. tmatrix1%ydim.ne.tmatrix2%ydim) THEN
      call errore(subname,"copying twin matrices with incompatible size", 1)
   ENDIF

   IF(tmatrix1%iscmplx) THEN
     IF(tmatrix2%iscmplx) THEN
      tmatrix1%cvec(:,:)=tmatrix2%cvec(:,:)+coef*tmatrix1%cvec(:,:)
     ELSE
      call errore(subname,"copying real tensor into complex tensor", 1)
      tmatrix1%cvec(:,:)=tmatrix2%rvec(:,:)+DBLE(coef*tmatrix1%cvec(:,:))
     ENDIF
   ELSE
     IF(.not.tmatrix2%iscmplx) THEN
      tmatrix1%rvec(:,:)=tmatrix2%rvec(:,:)+DBLE(coef)*tmatrix1%rvec(:,:)
     ELSE
      call errore(subname,"copying complex tensor into real tensor", 1)
      tmatrix1%rvec(:,:)=tmatrix2%cvec(:,:)+coef*tmatrix1%rvec(:,:)
     ENDIF
   ENDIF
   
   return
  END SUBROUTINE zaxpy_twin_matrix

  SUBROUTINE copy_twin_tensor(ttensor1,ttensor2)

   type(twin_tensor)   :: ttensor1, ttensor2
   character(len=17)   :: subname="copy_twin_tensor"
!    COMPLEX(DP), INTENT(IN) :: value

   IF(ttensor1%xdim.ne.ttensor2%xdim .OR. ttensor1%ydim.ne.ttensor2%ydim &
     .OR. ttensor1%zdim.ne.ttensor2%zdim) THEN
      call errore(subname,"copying twin tensors with incompatible size", 1)
   ENDIF

   IF(ttensor1%iscmplx) THEN
     IF(ttensor2%iscmplx) THEN
      ttensor1%cvec(:,:,:)=ttensor2%cvec(:,:,:)
     ELSE
      call errore(subname,"copying real tensor into complex tensor", 1)
      ttensor1%cvec(:,:,:)=ttensor2%rvec(:,:,:)
     ENDIF
   ELSE
     IF(.not.ttensor2%iscmplx) THEN
      ttensor1%rvec(:,:,:)=ttensor2%rvec(:,:,:)
     ELSE
      call errore(subname,"copying complex tensor into real tensor", 1)
      ttensor1%rvec(:,:,:)=ttensor2%cvec(:,:,:)
     ENDIF
   ENDIF
   
   return
  END SUBROUTINE copy_twin_tensor

  SUBROUTINE zaxpy_twin_tensor(ttensor1,ttensor2, coef)

   type(twin_tensor)   :: ttensor1,ttensor2
   COMPLEX(DP)    :: coef
   character(len=17)   :: subname="copy_twin_tensor"
!    COMPLEX(DP), INTENT(IN) :: value

   IF(ttensor1%xdim.ne.ttensor2%xdim .OR. ttensor1%ydim.ne.ttensor2%ydim &
     .OR. ttensor1%zdim.ne.ttensor2%zdim) THEN
      call errore(subname,"copying twin tensors with incompatible size", 1)
   ENDIF

   IF(ttensor1%iscmplx) THEN
     IF(ttensor2%iscmplx) THEN
      ttensor1%cvec(:,:,:)=ttensor2%cvec(:,:,:)+coef*ttensor1%cvec(:,:,:)
     ELSE
      call errore(subname,"copying real tensor into complex tensor", 1)
      ttensor1%cvec(:,:,:)=ttensor2%rvec(:,:,:)+DBLE(coef*ttensor1%cvec(:,:,:))
     ENDIF
   ELSE
     IF(.not.ttensor2%iscmplx) THEN
      ttensor1%rvec(:,:,:)=ttensor2%rvec(:,:,:)+DBLE(coef)*ttensor1%rvec(:,:,:)
     ELSE
      call errore(subname,"copying complex tensor into real tensor", 1)
      ttensor1%rvec(:,:,:)=ttensor2%cvec(:,:,:)+coef*ttensor1%rvec(:,:,:)
     ENDIF
   ENDIF
   
   return
  END SUBROUTINE zaxpy_twin_tensor

  SUBROUTINE allocate_twin_matrix(tmatrix,xlen,ylen,doreal)
   
   type(twin_matrix)   :: tmatrix
   INTEGER, INTENT(IN) :: xlen,ylen
   LOGICAL, INTENT(IN) :: doreal

   character(len=24)   :: subname="allocate_twin_matrix"
   INTEGER             :: ierr

   IF(tmatrix%isalloc) THEN
     call deallocate_twin(tmatrix)
   ENDIF

   IF(.not.doreal) THEN
   !write(6,*) "TWIN:allocating complex matrix", xlen, ylen
!      nullify(tmatrix%cvec)
     ALLOCATE(tmatrix%cvec(max(1,xlen),max(1,ylen)), STAT=ierr)
     IF(ierr/=0) call errore(subname,"allocating twin_matrix cvec", abs(ierr))
     tmatrix%iscmplx=.true.
     tmatrix%cvec=CMPLX(0.d0,0.d0)
   ELSE
   !write(6,*) "TWIN:allocating real matrix", xlen, ylen
!      nullify(tmatrix%rvec)
     allocate(tmatrix%rvec(max(1,xlen),max(1,ylen)), STAT=ierr)
     IF(ierr/=0) call errore(subname,"allocating twin_matrix rvec", abs(ierr))
     tmatrix%iscmplx=.false.
     tmatrix%rvec=0.d0
   ENDIF

   tmatrix%xdim=max(1,xlen)
   tmatrix%ydim=max(1,ylen)
   tmatrix%isalloc=.true.
   return

  END SUBROUTINE allocate_twin_matrix

  SUBROUTINE deallocate_twin_matrix(tmatrix)

    type(twin_matrix) :: tmatrix

    CHARACTER(len=26) :: subname="deallocate_twin_matrix"
    INTEGER           :: ierr

    IF(.not.tmatrix%iscmplx) THEN
       IF(associated(tmatrix%rvec)) THEN
	  deallocate(tmatrix%rvec, STAT=ierr)
         nullify(tmatrix%rvec)
         IF(ierr/=0) call errore(subname,"deallocating twin_matrix rvec", abs(ierr))       
       ENDIF
    ENDIF

    IF(tmatrix%iscmplx) THEN
       !write(6,*) "deallocating cvec", tmatrix%xdim, tmatrix%ydim, tmatrix%iscmplx, tmatrix%isalloc, associated(tmatrix%rvec)
       IF(associated(tmatrix%cvec)) THEN
	  deallocate(tmatrix%cvec, STAT=ierr)
          nullify(tmatrix%cvec)
          IF(ierr/=0) call errore(subname,"deallocating twin_matrix cvec", abs(ierr))
       ENDIF
    ENDIF

    tmatrix%xdim=0
    tmatrix%ydim=0
    tmatrix%iscmplx=.false.
    tmatrix%isalloc=.false.
    return
  END SUBROUTINE deallocate_twin_matrix

  SUBROUTINE allocate_twin_tensor(ttensor,xlen,ylen,zlen,doreal)
   
   type(twin_tensor) :: ttensor
   INTEGER, INTENT(IN) :: xlen,ylen,zlen
   LOGICAL :: doreal

   character(len=24)   :: subname="allocate_twin_tensor"
   INTEGER             :: ierr   

   IF(ttensor%isalloc) THEN
     call deallocate_twin(ttensor)
   ENDIF
   nullify(ttensor%rvec)
   nullify(ttensor%cvec)

   IF(.not.doreal) THEN
     ALLOCATE(ttensor%cvec(max(1,xlen),max(1,ylen),max(1,zlen)), STAT=ierr)
     IF(ierr/=0) call errore(subname,"allocating twin_tensor cvec", abs(ierr))
     ttensor%iscmplx=.true.
     ttensor%cvec=CMPLX(0.d0,0.d0)
   ELSE
     ALLOCATE(ttensor%rvec(max(1,xlen),max(1,ylen),max(1,zlen)), STAT=ierr)
     IF(ierr/=0) call errore(subname,"allocating twin_tensor rvec", abs(ierr))
     ttensor%iscmplx=.false.
     ttensor%rvec=0.d0
   ENDIF

   ttensor%xdim=max(1,xlen)
   ttensor%ydim=max(1,ylen)
   ttensor%ydim=max(1,zlen)
   ttensor%isalloc=.true.
   return

  END SUBROUTINE allocate_twin_tensor

  SUBROUTINE deallocate_twin_tensor(ttensor)

    type(twin_tensor) :: ttensor
    character(len=26) :: subname="deallocate_twin_tensor"

    INTEGER :: ierr

    IF(.not.ttensor%iscmplx) THEN
       IF(associated(ttensor%rvec)) THEN
	  DEALLOCATE(ttensor%rvec, STAT=ierr)
	  nullify(ttensor%rvec)
	  IF(ierr/=0) call errore(subname,"deallocating twin_tensor rvec", abs(ierr))
       ENDIF
    ENDIF
    
    IF(ttensor%iscmplx) THEN
       IF(associated(ttensor%cvec)) THEN
	  DEALLOCATE(ttensor%cvec, STAT=ierr)
	  nullify(ttensor%cvec)
	  IF(ierr/=0) call errore(subname,"deallocating twin_tensor cvec", abs(ierr))
       ENDIF
    ENDIF

    ttensor%xdim=0
    ttensor%ydim=0
    ttensor%zdim=0

    ttensor%iscmplx=.false.
    ttensor%isalloc=.false.
    nullify(ttensor%rvec)
    nullify(ttensor%cvec)

    return
  END SUBROUTINE deallocate_twin_tensor

  SUBROUTINE tmatrix_mp_sum(tmatrix)

    use mp, only: mp_sum, mp_bcast
    use mp_global,                ONLY : intra_image_comm
    
    IMPLICIT NONE

    type(twin_matrix) :: tmatrix    

    IF(.not.tmatrix%iscmplx) THEN
      call mp_sum(tmatrix%rvec, intra_image_comm)
    ELSE
      call mp_sum(tmatrix%cvec, intra_image_comm)
    ENDIF
   
  END SUBROUTINE tmatrix_mp_sum

  SUBROUTINE tmatrix_mp_bcast(tmatrix, node)

    use mp, only: mp_sum, mp_bcast
    use mp_global,                ONLY : intra_image_comm
    
    IMPLICIT NONE

    type(twin_matrix) :: tmatrix    
    integer :: node

    IF(.not.tmatrix%iscmplx) THEN
      call mp_bcast(tmatrix%rvec, node, intra_image_comm)
    ELSE
      call mp_bcast(tmatrix%cvec, node, intra_image_comm)
    ENDIF
   
  END SUBROUTINE tmatrix_mp_bcast

  SUBROUTINE ttensor_mp_bcast(ttensor, node)

    use mp, only: mp_sum, mp_bcast
    use mp_global,                ONLY : intra_image_comm
    
    IMPLICIT NONE

    type(twin_tensor) :: ttensor
    integer :: node

    IF(.not.ttensor%iscmplx) THEN
      call mp_bcast(ttensor%rvec, node, intra_image_comm)
    ELSE
      call mp_bcast(ttensor%cvec, node, intra_image_comm)
    ENDIF
   
  END SUBROUTINE ttensor_mp_bcast

  SUBROUTINE ttensor_mp_sum(ttensor)

    use mp, only: mp_sum, mp_bcast
    use mp_global,                ONLY : intra_image_comm
    
    IMPLICIT NONE

    type(twin_tensor) :: ttensor

    IF(.not.ttensor%iscmplx) THEN
      call mp_sum(ttensor%rvec, intra_image_comm)
    ELSE
      call mp_sum(ttensor%cvec, intra_image_comm)
    ENDIF
   
  END SUBROUTINE ttensor_mp_sum

  complex(DP) FUNCTION scalar_twin(vec1,vec2,sizevec,gstart,lgam)

   complex(DP), dimension(:) :: vec1,vec2
   integer :: gstart
   logical :: lgam
   integer :: sizevec
   character(len=12) :: subname="scalar_twin"

   complex(DP) :: aid

   if((size(vec1).ne.sizevec) .or. (size(vec2).ne.sizevec)) then
     call errore(subname,"inconsistent vector size", 1)
   endif

   aid=CMPLX(0.d0,0.d0)
   if(lgam) then
     aid=2.d0*DBLE(DOT_PRODUCT(vec1(1:sizevec),vec2(1:sizevec)))
     if(gstart == 2) then
       aid=aid - DBLE(CONJG(vec1(1))*vec2(1))
     endif
   else
     aid=DOT_PRODUCT(CONJG(vec1(1:sizevec)),vec2(1:sizevec))
   endif
   
   scalar_twin=aid   

  END FUNCTION scalar_twin

END MODULE twin_types
