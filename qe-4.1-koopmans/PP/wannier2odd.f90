subroutine wan2odd ( seedname, nbnd_to_save, iknum, ikstart, split_orbital_twofiles, spin_component, force_spin_symmetry)
  !
  USE io_global,       only : stdout
  USE kinds,           only : DP
  use io_files,        only : iunwfc, iunatsicwfc, nwordwfc, nwordwann
  USE cell_base,       only : omega, tpiba2
  use gvect,           only : g, ngm, ecutwfc, nrxx
  use gsmooth,         only: nls, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
  use wavefunctions_module, only : evc, psic
  use wvfct,           only : nbnd, npwx, npw, igk, g2kin
  use klist,           only : nkstot, xk, wk, nelec
  use lsda_mod,        only : nspin, isk
  use control_flags,    ONLY : gamma_only
  USE mp,        ONLY : mp_bcast, mp_sum
  use mp_global, ONLY : intra_pool_comm
  use read_wannier_files
  !
  !
  character(len=50), intent (in) :: seedname
  integer, intent (in) :: nbnd_to_save, iknum, ikstart
  logical, intent (in) :: split_orbital_twofiles
  !
  ! local variables
  ! 
  integer :: occ_nbnd, emp_nbnd, i, j, nn, ik, ibnd, iw, ikevc, counter, ibnd1
  complex(DP) :: overlap, ZDOTC
  real(DP):: DDOT
  complex (DP), allocatable:: evc_disen(:,:)
  complex (DP), allocatable:: evc_tmp(:,:)
  complex (DP), allocatable:: evc_disen_tmp(:,:)
  complex (DP), allocatable:: c0_tmp(:,:)
  complex (DP), allocatable:: c0_tmp_to_save(:,:)
  complex (DP), allocatable:: c0_tmp_to_save_occ(:,:)
  complex (DP), allocatable:: c0_tmp_to_save_emp(:,:)
  real (DP), allocatable:: wannier_spreads_tmp(:)
  logical, allocatable:: inc_band(:),  exclude_band(:)
  !
  character(len=4) :: spin_component
  CHARACTER(LEN=256) :: filename
  logical :: force_spin_symmetry
  ! 
  ! main program
  !
  call conv_read_chkpt(seedname)
  !
  if (have_disentangled) allocate (evc_disen (npw, maxval(ndimwin)) )
  !
  allocate (evc_tmp (npw, num_bands))
  evc_tmp(:,:) = (0.0, 0.0)
  !
  if (have_disentangled) then 
     allocate (inc_band ( maxval(ndimwin)) )
  else
     allocate (inc_band (num_bands))
  endif
  !
  do ik=1, iknum
     !
     inc_band(:) = .false.
     num_inc=num_wann
     if (have_disentangled) then
        !
        inc_band(:) = lwindow(:,ik)
        num_inc = ndimwin(ik)
        !
     endif
     !
     ikevc = ik + ikstart - 1
     call davcio (evc, nwordwfc, iunwfc, ikevc, -1)
     call gk_sort (xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
     !
     allocate ( exclude_band(nbnd) )
     exclude_band(:) = .false.
     do ibnd = 1, num_exclude_bands
        exclude_band(exclude_bands(ibnd)) = .true.  
     enddo 
     ! 
     counter = 0
     do ibnd = 1, nbnd
        !
        if (.not.exclude_band(ibnd)) then
           counter = counter + 1
           evc_tmp(:,counter) = evc(:, ibnd)
        endif
        ! 
     enddo 
     !
     if ( counter.ne.num_bands ) &
        call errore('wan2odd',' counter.ne.num_bands',1)
     !
     if (have_disentangled) then
        !
        counter=1
        ! 
        do loop_b=1, num_bands
           ! 
           if(counter > num_inc) exit
           !
           if (inc_band(loop_b)) then 
              !
              evc_disen (:, counter) = evc_tmp (:, loop_b)  
              !
              counter=counter + 1
              ! 
           endif 
           !
        enddo
        !
        if ( (counter-1).ne. num_inc ) &
           call errore('wan2odd',' counter.ne.num_inc',1)
        !
        allocate (evc_disen_tmp(npw, num_wann))
        evc_disen_tmp(:,:) = (0.0, 0.0)
        !  
        do loop_w=1, num_wann
           !    
           do loop_b=1, num_inc
              !
              evc_disen_tmp (:,loop_w) = evc_disen_tmp (:,loop_w) + &
                      u_matrix_opt(loop_b,loop_w, ik) * evc_disen (:,loop_b)   
              !
           enddo
           !   
        enddo
        ! 
     endif
     !
     ! Here is the core of this routine
     !   
     allocate (c0_tmp (npw, num_wann))
     c0_tmp(:,:) = (0.0, 0.0)
     ! 
     do ibnd = 1, num_wann
        !
        do jbnd = 1, num_wann
           !
           if (have_disentangled) then
              c0_tmp(:,jbnd) = c0_tmp(:,jbnd) + u_matrix(ibnd,jbnd, ik) * evc_disen_tmp(:,ibnd)
           else
              c0_tmp(:,jbnd) = c0_tmp(:,jbnd) + u_matrix(ibnd,jbnd, ik) * evc_tmp(:,ibnd)
           endif
           !
        enddo
        !
     enddo
     !
     if (have_disentangled) deallocate(evc_disen_tmp)
     ! 
  enddo ! k-points
  !
  write(stdout, *) "Linh: print overlap between wannier and evc"
  !
  do ibnd = 1, num_wann
     ibnd1 = 0
     do jbnd = 1, nbnd
        if (exclude_band(jbnd)) cycle
        if (gamma_only) then
            overlap = 2.0_dp*DDOT(2*npw, evc(1,jbnd),1, c0_tmp(1,ibnd), 1)
            if (gstart==2) overlap = overlap - real(conjg(evc(1,jbnd))*c0_tmp(1,ibnd))
        else
            overlap = ZDOTC(npw, evc(1,jbnd),1, c0_tmp(1,ibnd),1)
        end if
        call mp_sum(overlap, intra_pool_comm)
        ibnd1=ibnd1+1
        if (ibnd1==1) then
            write(stdout, *) 'Linh, ', ibnd, jbnd, real(overlap), aimag(overlap)
        endif
     end do
  end do
  !
  if (split_orbital_twofiles) then
     !
     if (nbnd_to_save .le. int(nelec/2)) then
           call errore('wan2odd','canot split two file nbnd_to_save .le. int(nelec/2)', 1)
     endif
     !
     occ_nbnd = int(nelec/2)    !! FIXME: WHAT IF nele ODD ?????
     emp_nbnd = nbnd_to_save - int(nelec/2) 
     !
     allocate (wannier_spreads_tmp(num_wann))
     
     wannier_spreads_tmp(:) = wannier_spreads(:)
     !
     CALL  Sort(npw, c0_tmp, wannier_spreads_tmp, num_wann )
     !
     deallocate (wannier_spreads_tmp)
     !
     if (nspin == 2) then
        IF (force_spin_symmetry) THEN 
           allocate (c0_tmp_to_save_occ (npw, occ_nbnd*2) )
        ELSE
           allocate (c0_tmp_to_save_occ (npw, occ_nbnd) )
        ENDIF
     else
        allocate (c0_tmp_to_save_occ (npw, occ_nbnd) )
     endif
     !
     if (nspin == 2) then
        !
        IF (force_spin_symmetry) THEN
           allocate (c0_tmp_to_save_emp (npw, emp_nbnd*2) )
        ELSE
           allocate (c0_tmp_to_save_emp (npw, emp_nbnd) )
        ENDIF
        !
     else
        allocate (c0_tmp_to_save_emp (npw, emp_nbnd) )
     endif  
     !
     c0_tmp_to_save_occ(:,:) = (0.0, 0.0)
     c0_tmp_to_save_emp(:,:) = (0.0, 0.0)
     !
     do ibnd=1, occ_nbnd 
        !
        c0_tmp_to_save_occ (:, ibnd) = c0_tmp (:, ibnd)
        !
     enddo
     !
     do ibnd=1, emp_nbnd 
        !
        c0_tmp_to_save_emp (:, ibnd) = c0_tmp (:, occ_nbnd+ibnd)
        !
     enddo
     !
     if (nspin == 2) then
        ! 
        IF (force_spin_symmetry) THEN 
          counter = 0
          !
          do ibnd=occ_nbnd+1, occ_nbnd*2
            !
             counter = counter + 1
             c0_tmp_to_save_occ (:, ibnd) = c0_tmp_to_save_occ(:, counter)
             ! 
          enddo
          !
          counter = 0
          !
          do ibnd=emp_nbnd+1, emp_nbnd*2
             !
             counter = counter + 1
             c0_tmp_to_save_emp (:, ibnd) = c0_tmp_to_save_emp(:, counter)
             ! 
          enddo
          ! 
          call writeempty_state(c0_tmp_to_save_occ, occ_nbnd*2, 'evc_occupied.dat')
          call writeempty_state(c0_tmp_to_save_emp, emp_nbnd*2, 'evc0_empty.dat')
          !
        ELSE
          ! 
          filename='evc_occupied_'//TRIM(spin_component)//'.dat'
          call writeempty_state(c0_tmp_to_save_occ, occ_nbnd,filename)
          filename='evc0_empty_'//TRIM(spin_component)//'.dat'
          call writeempty_state(c0_tmp_to_save_emp, emp_nbnd,filename)

          !
        ENDIF
        !
     else   !!! Nspin = 1
        !
        call writeempty_state(c0_tmp_to_save_occ, occ_nbnd, 'evc_occupied.dat')
        call writeempty_state(c0_tmp_to_save_emp, emp_nbnd, 'evc0_empty.dat')
        !
     endif
     !
     deallocate ( c0_tmp, c0_tmp_to_save_occ, c0_tmp_to_save_emp, evc_tmp )
     !
  else
     !
     if (nspin == 2) then 
        !
        if (force_spin_symmetry) THEN 
            allocate (c0_tmp_to_save (npw, nbnd_to_save*2) )
         else
            allocate (c0_tmp_to_save (npw, nbnd_to_save) )
         endif
         !
     else
        allocate (c0_tmp_to_save (npw, nbnd_to_save) )
     endif
     !
     c0_tmp_to_save(:,:) = (0.0, 0.0)
     !
     do ibnd=1, nbnd_to_save
        !
        c0_tmp_to_save (:, ibnd) = c0_tmp (:, ibnd)
        !
     enddo
     !
     if (nspin == 2) then 
       !
       if (force_spin_symmetry) then 
          !
          counter = 0
          !
          do ibnd=nbnd_to_save+1, nbnd_to_save*2
             !
             counter = counter + 1
             c0_tmp_to_save (:, ibnd) = c0_tmp (:, counter)
             ! 
          enddo
          !
          call writeempty_state(c0_tmp_to_save, nbnd_to_save*2, 'evc_empty.dat')
          !
        else 
          !
          filename='evc_empty_'//TRIM(spin_component)//'.dat'
          call writeempty_state(c0_tmp_to_save, nbnd_to_save, filename)
          !
        endif
        !
     else
        !  
        call writeempty_state(c0_tmp_to_save, nbnd_to_save, 'evc_empty.dat')
        !
     endif 
     !
     deallocate ( c0_tmp, c0_tmp_to_save, evc_tmp )
     !
  endif  
  !
  if (allocated(evc_disen)) deallocate (evc_disen)
  if (allocated(exclude_bands)) deallocate (exclude_bands)
  if (allocated(exclude_band)) deallocate (exclude_band)
  if (allocated(lwindow))  deallocate (lwindow)
  if (allocated(ndimwin)) deallocate (ndimwin)
  if (allocated(u_matrix_opt)) deallocate(u_matrix_opt)
  if (allocated(u_matrix)) deallocate(u_matrix)
  if (allocated(wannier_centres)) deallocate(wannier_centres)
  if (allocated(wannier_spreads)) deallocate(wannier_spreads)
  !
end subroutine wan2odd
!
!-----------------------------------------------------------------------
SUBROUTINE writeempty_state( c_emp, nempty, fileempty)
!-----------------------------------------------------------------------
        !
        ! ...   This subroutine writes empty states to unit emptyunitc0
        ! 
        USE kinds,              ONLY: DP
        USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
        USE mp_wave,            ONLY: mergewf
        USE mp,                 ONLY: mp_sum
        USE io_global,          ONLY: ionode, ionode_id, stdout
        USE reciprocal_vectors, ONLY: ig_l2g
        USE wvfct,              ONLY: npw
        USE xml_io_base,        ONLY: wfc_filename
        !
        IMPLICIT NONE
        ! 
        COMPLEX(DP), INTENT(IN) :: c_emp(npw, nempty)
        INTEGER,     INTENT(IN) :: nempty
        CHARACTER(LEN=256), INTENT(IN) :: fileempty
        !
        INTEGER :: ig, i, ngw_g, iss, ngw_l, emptyunit
        LOGICAL :: exst, ierr
        COMPLEX(DP), ALLOCATABLE :: ctmp(:)
        CHARACTER(LEN=256) ::  dirname
        !
        ! ... Subroutine Body
        !
        emptyunit = 100
        !
        ngw_g    = npw
        ngw_l    = npw
        !
        CALL mp_sum( ngw_g, intra_image_comm )
        !
        ALLOCATE( ctmp( ngw_g ) )
        !
        IF ( ionode ) THEN
           ! 
           ! 
           OPEN( UNIT = emptyunit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'UNFORMATTED' )
           !
           REWIND( emptyunit )
           !
           WRITE ( emptyunit )  ngw_g, nempty
           !
        ENDIF
        !
        DO i = 1, nempty
           !
           ctmp = (0.0d0, 0.0d0)
           !
           CALL mergewf ( c_emp(:,i), ctmp(:), ngw_l, ig_l2g, me_image, &
                         nproc_image, ionode_id, intra_image_comm )
           !
           IF ( ionode ) THEN
              !
              WRITE ( emptyunit ) ( ctmp(ig), ig=1, ngw_g )
              ! 
           ENDIF
           !
        ENDDO
        ! 
        IF ( ionode ) THEN
           !
           CLOSE ( emptyunit )
           !
        ENDIF
        !
        DEALLOCATE(ctmp)
        !
        RETURN
        ! 
END SUBROUTINE writeempty_state

   SUBROUTINE  Sort(ngw, cp, x, Size)
      USE kinds,           only : DP
      USE io_global,       only : stdout
      IMPLICIT  NONE
      INTEGER, INTENT(IN) :: ngw
      INTEGER, INTENT(IN)                   :: Size
      REAL(DP),    INTENT(INOUT) :: x(Size)
      COMPLEX(DP), INTENT(INOUT) :: cp(ngw, Size)
      !
      INTEGER                               :: i
      INTEGER                               :: Location
      INTEGER                               :: FindMinimum

      DO i = 1, Size-1                  ! except for the last
         Location = FindMinimum(x, i, Size)     ! find min from this to last
         CALL Swap(x(i), x(Location))  ! swap this and the minimum
         CALL dswap( 2*ngw*1, cp(:,i), 1, cp(:,Location), 1 )
      END DO
   END SUBROUTINE  Sort

INTEGER FUNCTION  FindMinimum(x, Start, End)
      USE kinds,           only : DP
      IMPLICIT  NONE
      INTEGER, INTENT(IN)                :: Start, End
      REAL(DP), INTENT(IN)               :: x (End)
      REAL(DP)                           :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)               ! assume the first is the min
      Location = Start                  ! record its position
      DO i = Start+1, End               ! start with next elements
         IF (x(i) < Minimum) THEN       !   if x(i) less than the min?
            Minimum  = x(i)             !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location            ! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      USE kinds,           only : DP
      IMPLICIT  NONE
      REAL(DP), INTENT(INOUT) :: a, b
      REAl(DP)                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

