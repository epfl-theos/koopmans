MODULE read_wannier_files
   !
   USE kinds, ONLY: DP
   !
   IMPLICIT NONE
   SAVE
   !
   complex (dp), allocatable :: u_matrix(:,:,:)
   complex (dp), allocatable :: u_matrix_opt (:,:,: )
   real(kind=dp), allocatable :: wannier_centres(:,:)
   real(kind=dp), allocatable :: wannier_spreads(:)
   integer, allocatable :: ndimwin (:)
   logical, allocatable :: lwindow (:,:) 
   integer, allocatable :: exclude_bands(:) 
   logical :: have_disentangled 
   integer :: num_bands, num_wann, num_exclude_bands
   !  
   CONTAINS
   !    
   subroutine conv_read_chkpt (seedname)
    !
    !=======================================!
    ! Read formatted checkpoint file        !
    !=======================================!
    use kinds,                only : dp
    USE mp_global,            ONLY: intra_image_comm
    USE mp,                   ONLY: mp_bcast, mp_sum
    USE io_global,            ONLY: ionode, ionode_id, stdout
    !
    implicit none
    !

    character(len=50), intent (in) :: seedname
    ! 
    integer :: chk_unit,i,j,k,l,nkp,ierr
    character(len=33) :: header
    character(len=20) :: checkpoint
    !
    real(kind=dp), allocatable :: kpt_latt(:,:) !kpoints in lattice vecs
    !
    ! u_matrix_opt gives the num_wann dimension optimal subspace from the
    ! original bloch states
    !
    ! u_matrix gives the unitary rotations from the optimal subspace to the
    ! optimally smooth states. 
    ! m_matrix we store here, because it is needed for restart of wannierise
    !
    complex(kind=dp), allocatable :: m_matrix(:,:,:,:)
    !
    real(kind=dp)  :: real_lattice(3,3)
    real(kind=dp)  :: recip_lattice(3,3)
    integer        :: mp_grid(3) 
    ! 
    integer :: num_kpts, nntot
    real(dp) :: omega_invariant
    !
    write(stdout,'(1x,3a)') 'Reading information from file ',trim(seedname),'.chk :'
    !
    chk_unit=100 ! io_file_unit()
    ! 
    if(ionode) then
      !
      open(unit=chk_unit,file=trim(seedname)//'.chk',status='old',form='unformatted')
      !
      ! Read comment line
      !
      read(chk_unit) header
      !
      ! Consistency checks
      !
      read(chk_unit) num_bands                           ! Number of bands
      read(chk_unit) num_exclude_bands                   ! Number of excluded bands
      !
    endif
    ! 
    call mp_bcast( header, ionode_id, intra_image_comm ) 
    call mp_bcast( num_bands, ionode_id, intra_image_comm ) 
    call mp_bcast( num_exclude_bands, ionode_id, intra_image_comm ) 
    !
    write(stdout,'(1x,a)') trim(header)
    ! 
    write(stdout,'(a,i0)') "Number of bands: ", num_bands
    !
    if (num_exclude_bands < 0) then
       ! 
       call errore('conv_read_chkpt', 'Invalid value for num_exclude_bands', num_exclude_bands)
       !
    endif
    !
    allocate(exclude_bands(num_exclude_bands),stat=ierr)
    ! 
    exclude_bands = .false.
    if(ionode) then
      ! 
      read(chk_unit) (exclude_bands(i),i=1,num_exclude_bands) ! Excluded bands
      ! 
    endif
    !
    call mp_bcast( exclude_bands, ionode_id, intra_image_comm ) 
    ! 
    write(stdout,'(a)',advance='no') "Excluded bands: "
    !
    if (num_exclude_bands == 0) then
       !
       write(stdout,'(a)') "none."
       !
    else
       !
       do i=1,num_exclude_bands-1
          !
          write(stdout,'(I0,a)',advance='no') exclude_bands(i), ','
          !
       end do
       !
       write(stdout,'(I0,a)') exclude_bands(num_exclude_bands), '.'
       !
    end if
    !
    if(ionode) then
      !
      read(chk_unit) ((real_lattice(i,j),i=1,3),j=1,3)  ! Real lattice
      read(chk_unit) ((recip_lattice(i,j),i=1,3),j=1,3) ! Reciprocal lattice
      read(chk_unit) num_kpts                           ! K-points
      read(chk_unit) (mp_grid(i),i=1,3)                 ! M-P grid
      !
    endif
    !
    call mp_bcast( real_lattice, ionode_id, intra_image_comm ) 
    call mp_bcast( recip_lattice, ionode_id, intra_image_comm ) 
    call mp_bcast( num_kpts, ionode_id, intra_image_comm ) 
    call mp_bcast( mp_grid, ionode_id, intra_image_comm ) 
    !
    write(stdout,'(a)') "Real lattice: read."
    !
    write(stdout,'(a)') "Reciprocal lattice: read."
    !
    write(stdout,'(a,I0)') "Num kpts:", num_kpts
    !
    write(stdout,'(a)') "mp_grid: read."
    !
    if (.not.allocated(kpt_latt)) allocate(kpt_latt(3,num_kpts),stat=ierr)
    !
    if(ionode) then
      !
      read(chk_unit) ((kpt_latt(i,nkp),i=1,3),nkp=1,num_kpts)
      read(chk_unit) nntot                ! nntot
      read(chk_unit) num_wann                ! num_wann
      read(chk_unit) checkpoint             ! checkpoint
      read(chk_unit) have_disentangled      ! whether a disentanglement has been performed
      !
    endif  
    !
    call mp_bcast( kpt_latt, ionode_id, intra_image_comm )  
    call mp_bcast( nntot, ionode_id, intra_image_comm )  
    call mp_bcast( num_wann, ionode_id, intra_image_comm )  
    call mp_bcast( checkpoint, ionode_id, intra_image_comm )  
    call mp_bcast( have_disentangled, ionode_id, intra_image_comm )  
    !
    write(stdout,'(a)') "kpt_latt: read."
    !
    write(stdout,'(a,I0)') "nntot:", nntot
    ! 
    write(stdout,'(a,I0)') "num_wann:", num_wann
    !
    checkpoint=adjustl(trim(checkpoint))
    !
    write(stdout,'(a,I0)') "checkpoint: " // trim(checkpoint)
    !
    if (have_disentangled) then
       !
       write(stdout,'(a)') "have_disentangled: TRUE"
       !
       if (ionode) then
          read(chk_unit) omega_invariant     ! omega invariant
       endif
       ! 
       call mp_bcast( omega_invariant, ionode_id, intra_image_comm )  
       !
       write(stdout,'(a)') "omega_invariant: read."
       ! 
       ! lwindow
       if (.not.allocated(lwindow)) then
          !
          allocate(lwindow(num_bands,num_kpts),stat=ierr)
          !
       endif
       !
       lwindow = .false.
       if (ionode) then
          read(chk_unit) ((lwindow(i,nkp),i=1,num_bands),nkp=1,num_kpts)
       endif
       !
       call mp_bcast( lwindow, ionode_id, intra_image_comm )  
       !
       write(stdout,'(a)') "lwindow: read."
       !
       ! ndimwin
       ! 
       if (.not.allocated(ndimwin)) then
          !
          allocate(ndimwin(num_kpts),stat=ierr)
          !
       endif
       !
       if (ionode) then
          read(chk_unit) (ndimwin(nkp),nkp=1,num_kpts)
       endif
       !
       call mp_bcast( ndimwin, ionode_id, intra_image_comm )
       !
       write(stdout,'(a)') "ndimwin: read."
       !
       ! U_matrix_opt
       !
       if (.not.allocated(u_matrix_opt)) then
          ! 
          allocate(u_matrix_opt(num_bands,num_wann,num_kpts),stat=ierr)
          !
       endif
       ! 
       if (ionode) then
          read(chk_unit) (((u_matrix_opt(i,j,nkp),i=1,num_bands),j=1,num_wann),nkp=1,num_kpts)
       endif
       !
       call mp_bcast( u_matrix_opt, ionode_id, intra_image_comm )
       !
       write(stdout,'(a)') "U_matrix_opt: read."
       !
    else
       ! 
       write(stdout,'(a)') "have_disentangled: FALSE"
       !
    endif
    !
    ! U_matrix
    !
    if (.not.allocated(u_matrix)) then
       !
       allocate(u_matrix(num_wann,num_wann,num_kpts),stat=ierr)
       !
    endif
    !
    if (ionode) then
       read(chk_unit) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)
       ! NLN
       !DO i=1,num_wann
       !   DO j=1,num_wann
       !      WRITE(*,*) '========== U matrix =========='
       !      WRITE(*,*) i,j,1
       !      WRITE(*,*) u_matrix(i,j,1)
       !   ENDDO
       !ENDDO
    endif
    !
    call mp_bcast( u_matrix, ionode_id, intra_image_comm )
    !
    write(stdout,'(a)') "U_matrix: read."
    !
    ! M_matrix
    !
    if (.not.allocated(m_matrix)) then
       ! 
       allocate(m_matrix(num_wann,num_wann,nntot,num_kpts),stat=ierr)
       ! 
    endif
    !
    if (ionode) then
       read(chk_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts)
    endif
    !
    call mp_bcast( m_matrix, ionode_id, intra_image_comm ) 
    !
    write(stdout,'(a)') "M_matrix: read."
    !
    ! wannier_centres
    ! 
    if (.not.allocated(wannier_centres)) then
       ! 
       allocate(wannier_centres(3,num_wann),stat=ierr)
       !
    end if
    !
    if (ionode) then
       read(chk_unit) ((wannier_centres(i,j),i=1,3),j=1,num_wann)
    endif
    !
    call mp_bcast( wannier_centres, ionode_id, intra_image_comm )
    !
    write(stdout,'(a)') "wannier_centres: read."    
    !
    ! wannier spreads
    !
    if (.not.allocated(wannier_spreads)) then
       !
       allocate(wannier_spreads(num_wann),stat=ierr)
       !
    end if
    !
    if (ionode) then
       read(chk_unit) (wannier_spreads(i),i=1, num_wann)
    endif
    call mp_bcast( wannier_spreads, ionode_id, intra_image_comm )
    !
    write(stdout,'(a)') "wannier_spreads: read."
    write(stdout,*) wannier_spreads(1)
    !
100 continue
    ! 
    close(chk_unit)
    !
    write(stdout,'(a/)') ' ... done'
    !
    if (allocated(kpt_latt)) deallocate (kpt_latt)
    if (allocated(m_matrix)) deallocate(m_matrix)
    !if (allocated(wannier_centres)) deallocate(wannier_centres)
    !if (allocated(wannier_spreads)) deallocate(wannier_spreads)
    !
    return
    !
   end subroutine conv_read_chkpt

END MODULE read_wannier_files
