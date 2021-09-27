!=----------------------------------------------------------------------------=!
SUBROUTINE wave_init_wannier_pwscf ( c0, nbndx )
!=----------------------------------------------------------------------------=!
      !
      !  this routine read the occupied  wannier wave functions from pw2wannier PP code
      !  and the wfc will be used as starting wfc for ODD functional. This will use
      !  readempty routine to do this job
      !  The same for empty states, but it will be read in the empty.f90 
      ! 
      USE kinds,              ONLY: DP
      USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
      USE io_global,          ONLY: stdout, ionode, ionode_id
      USE mp,                 ONLY: mp_bcast, mp_sum
      USE mp_wave,            ONLY: splitwf
      USE io_files,           ONLY: outdir
      USE io_files,           ONLY: emptyunitc0
      USE reciprocal_vectors, ONLY: ig_l2g
      USE gvecw,              ONLY: ngw
      USE xml_io_base,        ONLY: restart_dir, wfc_filename
      USE control_flags,      ONLY: ndr
      !
      IMPLICIT none
      ! 
      INTEGER,     INTENT(IN)  :: nbndx
      COMPLEX(DP), INTENT(OUT) :: c0(ngw,nbndx)
      ! 
      LOGICAL :: exst
      INTEGER :: ig, ibnd
      INTEGER :: ngw_rd, nbnd_rd, ngw_l, ngw_g
      ! 
      CHARACTER(LEN=256) :: fileocc, dirname
      !
      COMPLEX(DP), ALLOCATABLE :: ctmp(:)
      !
      ! ... Subroutine Body
      !
      ngw_g    = ngw
      ngw_l    = ngw
      !
      CALL mp_sum( ngw_g, intra_image_comm )
      !
      ALLOCATE( ctmp(ngw_g) )
      !
      dirname   = restart_dir( outdir, ndr )
      !
      fileocc = TRIM( wfc_filename( dirname, 'evc_occupied', 1 ) )
      !
      IF ( ionode ) THEN
         !
         INQUIRE( FILE = TRIM(fileocc), EXIST = exst )
         !
         IF ( exst ) THEN
            !
            OPEN( UNIT=emptyunitc0, FILE=TRIM(fileocc), STATUS='OLD', FORM='UNFORMATTED' )
            !
            READ(emptyunitc0) ngw_rd, nbnd_rd
            !
            IF ( ngw_g .ne. ngw_rd ) THEN
               !
               exst = .false.
               WRITE( stdout,10)  TRIM(fileocc) 
               WRITE( stdout,20)  ngw_g, ngw_rd
               !
            ENDIF
            !
         ENDIF
         !
      ENDIF
      !
 10   FORMAT('*** OCCUPIED STATES : wavefunctions dimensions changed  ', A )
 20   FORMAT('*** NGW_G = ', I8, ' NE_READ = ', I4)
      !  
      CALL mp_bcast(exst,   ionode_id, intra_image_comm)
      CALL mp_bcast(nbnd_rd,  ionode_id, intra_image_comm)
      CALL mp_bcast(ngw_rd, ionode_id, intra_image_comm)
      !
      c0(:,:) = (0.0_DP, 0.0_DP)
      !
      IF ( exst ) THEN
         ! 
         DO ibnd = 1, MIN( nbndx, nbnd_rd )
            !
            IF ( ionode ) THEN
               ! 
               READ(emptyunitc0) ( ctmp(ig), ig = 1, MIN( SIZE(ctmp), ngw_rd ) )
               !
            ENDIF
            ! 
            IF ( ibnd <= nbnd_rd ) THEN
               !
               CALL splitwf(c0(:,ibnd), ctmp, ngw_l, ig_l2g, me_image, &
                            nproc_image, ionode_id, intra_image_comm)
               !
            ENDIF
            !
         ENDDO
         !
      ELSE
         !
         IF (.not. exst ) CALL errore( 'wave_init_wannier_pwscf', 'Something wrong with reading evc_occupied file', 1 )
         !
      ENDIF
      ! 
      IF ( ionode .AND. exst ) THEN
         !
         CLOSE(emptyunitc0) 
         !
      ENDIF
      ! 
      DEALLOCATE(ctmp)
      !
      RETURN
      !
END SUBROUTINE wave_init_wannier_pwscf

!=----------------------------------------------------------------------------=!
SUBROUTINE wave_init_wannier_cp ( c0, ngw,  nbndx, occupied_only)
!=----------------------------------------------------------------------------=!
      !
      !  this routine rotates the wave functions to the Kohn-Sham base
      !  it works with a block-like distributed matrix
      !  of the Lagrange multipliers ( lambda ).
      !
      ! ... declare modules
      ! 
      USE kinds,            ONLY: DP
      USE mp,               ONLY: mp_bcast
      USE dspev_module,     ONLY: dspev_drv
      USE cp_interfaces,    ONLY: set_evtot, readempty
      USE electrons_base,   ONLY: nspin, iupdwn, nupdwn
      USE electrons_module, ONLY: iupdwn_emp, nupdwn_emp
      USE wavefunctions_module, ONLY: ctot_aux
      !  
      IMPLICIT NONE
      !
      ! ... declare subroutine arguments
      !
      INTEGER, INTENT(IN) :: ngw, nbndx
      COMPLEX(DP), INTENT(INOUT) :: c0(ngw, nbndx)
      LOGICAL, INTENT(IN) :: occupied_only
      !
      ! ... declare other variables
      !
      COMPLEX(DP), ALLOCATABLE :: vv(:,:)
      COMPLEX(DP), ALLOCATABLE :: ctot(:,:), c0_tmp(:,:)
      INTEGER               :: i, j, k, nx, n_occs, npw
      !  
      IF (occupied_only) THEN
         !
         ! if this routine is called in occupied calculation,
         ! c0 is the variational orbitalsi, and ctot is tranformed 
         ! from c0. 
         !
         ALLOCATE( ctot( SIZE( c0, 1 ), nupdwn(1)*nspin ) )
         ctot = (0.0d0, 0.0d0)
         ctot(:,1:nupdwn(1)*nspin ) = ctot_aux(:,1:nupdwn(1)*nspin ) 
         !
         if (allocated(ctot_aux)) deallocate (ctot_aux)
         !write(stdout,*) "in restart kipz lambda up", lambda(1)%cvec(1:3,1)
         !write(stdout,*) "in restart kipz lambda dwn", lambda(2)%cvec(1:3,1)
         !CALL set_evtot( c0, ctot, lambda, iupdwn, nupdwn )
         !
         nx = nupdwn(1)
         !
         n_occs = nupdwn(1)
         IF( nspin == 2 ) n_occs = n_occs + nupdwn(2)
         !
      ELSE
         !
         ! if this routine is called in the empty calculation
         ! ctot (canonical orbitals) are read from file,
         !
         ALLOCATE( ctot( SIZE( c0, 1 ),  nupdwn_emp(1) * nspin ) )
         ctot = (0.0d0, 0.0d0)
         ctot(:,:) = c0(:,:) 
         !
         nx = nupdwn_emp(1)
         !
         n_occs = nupdwn_emp(1)
         IF( nspin == 2 ) n_occs = n_occs + nupdwn_emp(2)
         !
      END IF         
      !
      !write(stdout,*) "ctot in restart kipz", ctot(1:3, 1)
      !write(stdout,*) "c0 in restart kipz", c0(1:3, 1)
      ALLOCATE( vv( nx, nx ) )
      ! 
      ! read the optimal unitary matrix
      ! 
      CALL conv_read_chkpt(vv, nx)
      !
      ! convert canonical to minimizing orbitals
      ! 
      ALLOCATE( c0_tmp( SIZE( c0, 1 ), nx ) )  
      c0_tmp(:,:) = (0.0d0, 0.0d0)
      !
      npw=SIZE( c0, 1 )
      do k=1, nx
           do i=1, nx
              c0_tmp(:,i) = c0_tmp(:,i) + vv(k,i)*ctot(:,k)
           enddo !ibnd
      enddo !wannier
      !
      ! Symmetry the up and down spin
      !
      c0=(0.0d0, 0.0d0) 
      !
      IF (occupied_only) THEN
         !
         DO i = 1, MIN(nupdwn(1), n_occs)
            !
            c0(:,i) = c0_tmp(:, i)
            !  
         ENDDO
         !
         IF (nspin==2) THEN
            !
            DO i = 1, MIN(nupdwn(1),nupdwn(2))
               !
               j=i+iupdwn(2)-1
               c0(:,j) = c0_tmp(:,i)
               !
            ENDDO
            !
         ENDIF
         !
      ELSE
         !
         DO i = 1, MIN(nupdwn_emp(1), n_occs)
            !
            c0(:,i) = c0_tmp(:, i)
            !  
         ENDDO
         !  
         IF (nspin==2) THEN
            !
            DO i = 1, MIN(nupdwn_emp(1),nupdwn_emp(2))
               !
               j=i+iupdwn_emp(2)-1
               c0(:,j) = c0_tmp(:,i)
               !
            ENDDO
            !
         ENDIF 
         !
      ENDIF

      DEALLOCATE( vv )
      DEALLOCATE( ctot , c0_tmp)
      !
      RETURN
      ! 
END SUBROUTINE wave_init_wannier_cp
!
!
subroutine conv_read_chkpt(u_matrix_pass, nx)
    !
    !=======================================!
    ! Read formatted checkpoint file        !
    !=======================================!
    use kinds,                only : dp
    USE mp_global,            ONLY: intra_image_comm
    USE mp,                   ONLY: mp_bcast, mp_sum
    USE io_global,            ONLY: ionode, ionode_id, stdout
    USE input_parameters,     ONLY: which_file_wannier
    !
    implicit none
    !
    integer, intent (in) :: nx
    complex(dp), intent(inout) :: u_matrix_pass(nx, nx)
    ! 
    integer :: chk_unit,i,j,k,l,nkp,ierr
    character(len=50) :: seedname
    !
    character(len=33) :: header
    character(len=20) :: checkpoint
    !
    integer, allocatable :: exclude_bands(:)
    real(kind=dp), allocatable :: kpt_latt(:,:) !kpoints in lattice vecs
    integer, allocatable :: ndimwin(:)
    logical, allocatable :: lwindow(:,:)
    !
    ! u_matrix_opt gives the num_wann dimension optimal subspace from the
    ! original bloch states
    !
    complex(kind=dp), allocatable :: u_matrix_opt(:,:,:)
    !
    ! u_matrix gives the unitary rotations from the optimal subspace to the
    ! optimally smooth states. 
    ! m_matrix we store here, because it is needed for restart of wannierise
    !
    complex(kind=dp), allocatable :: u_matrix(:,:,:)
    complex(kind=dp), allocatable :: m_matrix(:,:,:,:)
    ! Wannier centres and spreads
    !
    real(kind=dp), allocatable :: wannier_centres(:,:)
    real(kind=dp), allocatable :: wannier_spreads(:)
    !
    real(kind=dp)  :: real_lattice(3,3)
    real(kind=dp)  :: recip_lattice(3,3)
    integer        :: mp_grid(3) 
    logical   :: have_disentangled
    !
    ! 
    integer :: num_bands, num_exclude_bands
    integer :: num_kpts, nntot, num_wann
    real(dp)    :: omega_invariant
    !
    seedname=which_file_wannier
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
    close(chk_unit)
    !
    write(stdout,'(a/)') ' ... done'
    !
    ! RETURN U_MATRIX
    !
    if (nx/=num_wann) call errore('conv_read_chkpt', 'different number nx and num_wann', nx)
    u_matrix_pass (:,:) = u_matrix(:,:, 1)
    !
    if (allocated(exclude_bands)) deallocate (exclude_bands)
    if (allocated(kpt_latt)) deallocate (kpt_latt)
    if (allocated(lwindow))  deallocate (lwindow)
    if (allocated(ndimwin)) deallocate (ndimwin)
    if (allocated(u_matrix_opt)) deallocate(u_matrix_opt)
    if (allocated(u_matrix)) deallocate(u_matrix)
    if (allocated(m_matrix)) deallocate(m_matrix)
    if (allocated(wannier_centres)) deallocate(wannier_centres)
    if (allocated(wannier_spreads)) deallocate(wannier_spreads)
    !
    return
    !
end subroutine conv_read_chkpt
