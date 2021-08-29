!--------------------------------------------------------------
subroutine odd_alpha_routine( nbndx, is_empty)
!--------------------------------------------------------------
      !
      ! alpha_v = (\sum (I,l) {alpha0_I <c_v|chi(I,l)><chi(I,l)|c_v>}) / (\sum (I,l) {<c_v|chi(I,l)><chi(I,l)|c_v>}),
      ! where chi(I,l) is orthornormal pseudo atomic wfc
      ! input: evc_wfc
      ! out: odd_alpha, valpsi >>> pass to sharing variable
      ! full valpsi  will be complete in nksiclib routine 
      !
      !
      use kinds,              ONLY: DP        
      USE mp_global,          ONLY: intra_image_comm
      USE mp,                 ONLY: mp_bcast, mp_sum
      USE io_global,          ONLY: ionode, ionode_id
      use orthogonalize_base, ONLY: calphi
      USE nksic,              ONLY: odd_alpha, valpsi, &
                                    call_index_emp, call_index, &
                                    alpha0_ref_emp, alpha0_ref 
      USE twin_types
      !
      implicit none
      !
#ifdef __PARA
      include 'mpif.h'
#endif
      !
      integer,     intent(in)  :: nbndx
      logical :: is_empty
      !
      ! local variables
      !
      integer ::  n_ref_nkscal0  
      integer,     allocatable :: index_ref_nkscal0(:)
      real(DP),    allocatable :: ref_nkscal0(:)
      real(DP),    allocatable :: spead_ref_nkscal0(:)
      ! 
      integer :: i, iv
      !
      ! this part is computed only one at starting calculation
      !
      !
      if (is_empty) then
         ! 
         call_index_emp = call_index_emp + 1  
         !
      else
         !
         call_index = call_index + 1 
         !
      endif
      !
      if ( ((.not. is_empty).and.(call_index .eq. 1)) .or. (is_empty .and. (call_index_emp .eq. 1)) ) then
         !
         if (is_empty .and. (call_index_emp .eq. 1) ) then
            !
            ! alpha0_ref_emp first is allocated here, then will be deallocated 
            ! in the end of cp run in deallocate_nksic module.f90
            !
            allocate (alpha0_ref_emp(nbndx))
            alpha0_ref_emp = 0.d0
            !
         endif
         !   
         if ((.not. is_empty) .and. (call_index .eq. 1)) then
            !  
            ! alpha0_ref_emp first is allocated here, then will be deallocated 
            ! in the end of cp run in deallocate_nksic module.f90
            !   
            allocate (alpha0_ref(nbndx)) 
            alpha0_ref = 0.d0
            !
         endif
         !  
         ! read from file ref. alpha
         !
         if (ionode) then
            !
            if (is_empty .and. (call_index_emp .eq. 1) ) then
               !
               open (unit = 99, file = 'file_alpharef_empty.txt', form = 'formatted', status = 'old' )
               !
               read(99, *), n_ref_nkscal0
               !
            endif 
            ! 
            if ((.not. is_empty) .and. (call_index .eq. 1)) then
               !
               open (unit = 99, file = 'file_alpharef.txt', form = 'formatted', status = 'old' )
               !
               read(99, *), n_ref_nkscal0
               !
            endif
            !
         endif
         ! 
         call mp_bcast( n_ref_nkscal0, ionode_id, intra_image_comm )
         ! 
         allocate(index_ref_nkscal0(n_ref_nkscal0))
         allocate(ref_nkscal0(n_ref_nkscal0))
         allocate(spead_ref_nkscal0(n_ref_nkscal0))
         !
         if (ionode) then
            !  
            do i = 1, n_ref_nkscal0
               ! 
               read (99, * ) index_ref_nkscal0(i), ref_nkscal0(i), spead_ref_nkscal0(i) 
               !
            enddo
            !
            close (99)
            ! 
         endif 
         !
         call mp_bcast( index_ref_nkscal0, ionode_id, intra_image_comm )
         call mp_bcast( ref_nkscal0, ionode_id, intra_image_comm )
         call mp_bcast( spead_ref_nkscal0, ionode_id, intra_image_comm )
         !
         ! first, we assign the refalpha to all fixed orbitals,
         !   
         do iv = 1, n_ref_nkscal0
            !
            if (is_empty) then
               !
               alpha0_ref_emp(iv) = ref_nkscal0(iv)
               ! 
            else
               !
               alpha0_ref(iv) = ref_nkscal0(iv)
               !
            endif
            !  
         enddo
         !
         deallocate(index_ref_nkscal0)
         deallocate(ref_nkscal0)
         deallocate(spead_ref_nkscal0)
         ! 
      endif 
      !
      ! if the calculation does not update alpha wrt minimization wfc
      ! we return from here
      !
      if (is_empty) then
         !
         odd_alpha(:) = alpha0_ref_emp(:)
         !
      else
         !
         odd_alpha(:) = alpha0_ref(:)
         !
      endif
      !   
      valpsi(:,:)  = (0.0_DP, 0.0_DP)
      !
      return
      !  
end subroutine odd_alpha_routine
