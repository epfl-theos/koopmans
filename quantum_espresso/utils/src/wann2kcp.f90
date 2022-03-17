!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------------------
PROGRAM wann2kcp
  !------------------------------------------------------------------------
  !
  !
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE mp_global,        ONLY : mp_startup
  USE mp_pools,         ONLY : npool
  USE mp_bands,         ONLY : nbgrp
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  USE cell_base,        ONLY : at, bg
  USE lsda_mod,         ONLY : nspin, isk
  USE klist,            ONLY : nkstot, xk
  USE io_files,         ONLY : prefix, tmp_dir
  USE noncollin_module, ONLY : noncolin
  USE control_flags,    ONLY : gamma_only
  USE environment,      ONLY : environment_start, environment_end
  USE wannier2kcp,      ONLY : wan2odd
  USE plot_wan2odd,     ONLY : plot_wann
  USE wannier
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios
  CHARACTER(LEN=4) :: spin_component
  CHARACTER(LEN=256) :: outdir
  !
  NAMELIST / inputpp / outdir, prefix, seedname, wan_mode, &
    spin_component, gamma_trick, print_rho, wannier_plot, wannier_plot_list
  !
  ! initialise environment
  !
#if defined( __MPI )
  CALL mp_startup( )
#endif
  !
  CALL environment_start( 'WANN2KCP' )
  !
  CALL start_clock( 'init_wann2kcp' )
  !
  ! Read input on i/o node and broadcast to the rest
  !
  ios = 0
  IF ( ionode ) THEN
    !
    ! Check to see if we are reading from a file
    !
    CALL input_from_file( )
    !
    !   set default values for variables in namelist
    !
    CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
    IF ( trim(outdir) == ' ' ) outdir = './'
    prefix = ' '
    seedname = 'wannier'
    wan_mode = 'wannier2kcp'
    spin_component = 'none'
    gamma_trick = .false.
    print_rho = .false.
    wannier_plot = .false.
    wannier_plot_list = 'all'
    !
    ! Reading the namelist inputpp
    !
    READ( 5, inputpp, iostat=ios )
    !
    ! Check of namelist variables
    !
    tmp_dir = trimcheck( outdir )
    !
  ENDIF
  !
  CALL mp_bcast( ios, ionode_id, world_comm )
  IF ( ios /= 0 ) CALL errore( 'wann2kcp', 'reading inputpp namelist', abs( ios ) )
  !
  ! Broadcast input variable to all nodes
  !
  CALL mp_bcast( outdir, ionode_id, world_comm )
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( seedname, ionode_id, world_comm )
  CALL mp_bcast( wan_mode, ionode_id, world_comm )
  CALL mp_bcast( spin_component, ionode_id, world_comm )
  CALL mp_bcast( gamma_trick, ionode_id, world_comm )
  CALL mp_bcast( print_rho, ionode_id, world_comm )
  CALL mp_bcast( wannier_plot, ionode_id, world_comm )
  CALL mp_bcast( wannier_plot_list, ionode_id, world_comm )
  !
  ! Check: kpoint distribution with pools not implemented
  !
  IF ( npool > 1 ) CALL errore( 'wann2kcp', 'pools not implemented', npool )
  !
  ! Check: bands distribution not implemented
  !
  IF ( nbgrp > 1 ) CALL errore( 'wann2kcp', 'bands (-nb) not implemented', nbgrp )
  !
  ! Now allocate space for pwscf variables, read and check them.
  !
  WRITE( stdout, * )
  WRITE( stdout, *) ' Reading nscf_save data'
  CALL read_file( )
  WRITE( stdout, * )
  !
  SELECT CASE ( trim( spin_component ) )
  CASE ( 'up' )
    WRITE( stdout, * ) ' Spin CASE ( up )'
    ispinw = 1
    ikstart = 1
    ikstop = nkstot / 2
    iknum = nkstot / 2
  CASE ( 'down' )
    WRITE( stdout, * ) 'Spin CASE ( down )'
    ispinw = 2
    ikstart = nkstot / 2 + 1
    ikstop = nkstot
    iknum = nkstot / 2
  CASE DEFAULT
    WRITE( stdout, * ) ' Spin CASE ( default = unpolarized )'
    ispinw = 0
    ikstart = 1
    ikstop = nkstot
    iknum = nkstot
  END SELECT
  !
  CALL stop_clock( 'init_wann2kcp' )
  !
  WRITE( stdout, * )
  WRITE( stdout, * ) ' Mode is: ', wan_mode
  WRITE( stdout, * )
  !
  IF ( wan_mode == 'wannier2kcp' ) THEN
    !
    CALL read_nnkp( )
    CALL get_wannier_to_plot( )
    CALL openfil_pp( )
    !
    CALL wan2odd( ks_only = .false. )
    !
    IF ( wannier_plot ) CALL plot_wann( wann_to_plot, iknum, n_wannier )
    !
    IF ( ionode ) WRITE( stdout, * )
    CALL print_clock( 'init_wann2kcp' )
    CALL print_clock( 'wannier2kcp' )
    IF ( wannier_plot ) CALL print_clock( 'plot_wann' )
    CALL environment_end( 'WANN2KCP' )
    IF ( ionode ) WRITE( stdout, * )
    !
    CALL stop_pp( )
    !
  ELSE IF ( wan_mode == 'ks2kcp' ) THEN
    !
    CALL openfil_pp( )
    CALL mp_grid_ks2kcp( )
    !
    CALL wan2odd( ks_only = .true. )
    !
    IF ( ionode ) WRITE( stdout, *  )
    CALL print_clock( 'init_wannk2kcp' )
    CALL print_clock( 'ks2kcp' )
    CALL environment_end( 'WANN2KCP' )
    IF ( ionode ) WRITE( stdout, * )
    !
    CALL stop_pp( )
    !
  ENDIF
  !
  !
END PROGRAM wann2kcp
!
!
!-------------------------------------------------------------------------
SUBROUTINE find_mp_grid( )
  !-----------------------------------------------------------------------
  !
  USE constants,     ONLY : eps8
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE wannier

  IMPLICIT NONE

  ! <<<local variables>>>
  INTEGER  :: ik,ntemp,ii
  real(DP) :: min_k,temp(3,iknum),mpg1

  min_k=minval(kpt_latt(1,:))
  ii=0
  DO ik=1,iknum
     IF (abs(kpt_latt(1,ik) - min_k) < eps8) THEN
        ii=ii+1
        temp(:,ii)=kpt_latt(:,ik)
     ENDIF
  ENDDO
  ntemp=ii

  min_k=minval(temp(2,1:ntemp))
  ii=0
  DO ik=1,ntemp
     IF (abs(temp(2,ik) - min_k) < eps8) THEN
        ii=ii+1
     ENDIF
  ENDDO
  mp_grid(3)=ii

  min_k=minval(temp(3,1:ntemp))
  ii=0
  DO ik=1,ntemp
     IF (abs(temp(3,ik) - min_k) < eps8) THEN
        ii=ii+1
     ENDIF
  ENDDO
  mp_grid(2)=ii

  IF ( (mp_grid(2)==0) .or. (mp_grid(3)==0) ) &
       CALL errore('find_mp_grid',' one or more mp_grid dimensions is zero', 1)

  mpg1=iknum/(mp_grid(2)*mp_grid(3))

  mp_grid(1) = nint(mpg1)

  WRITE(stdout,*)
  WRITE(stdout,'(3(a,i3))') '  MP grid is ',mp_grid(1),' x',mp_grid(2),' x',mp_grid(3)

  IF (real(mp_grid(1),kind=DP)/=mpg1) &
       CALL errore('find_mp_grid',' determining mp_grid failed', 1)

  RETURN
END SUBROUTINE find_mp_grid
!-----------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------
SUBROUTINE read_nnkp( )
  !-----------------------------------------------------------------------
  !
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE kinds,            ONLY : DP
  USE constants,        ONLY : eps6, tpi, bohr => BOHR_RADIUS_ANGS
  USE cell_base,        ONLY : at, bg, alat
  USE gvect,            ONLY : g, gg
  USE klist,            ONLY : nkstot, xk
  USE mp,               ONLY : mp_bcast, mp_sum
  USE mp_pools,         ONLY : intra_pool_comm
  USE mp_world,         ONLY : world_comm
  USE wvfct,            ONLY : npwx, nbnd
  USE noncollin_module, ONLY : noncolin
  USE wannier
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  real(DP) :: g_(3), gg_
  INTEGER :: ik, ib, ig, ipol, iw, idum, indexb
  INTEGER numk, i, j
  INTEGER, ALLOCATABLE :: ig_check(:,:)
  real(DP) :: xx(3), xnorm, znorm, coseno
  LOGICAL :: have_nnkp,found
  INTEGER :: tmp_auto ! vv: Needed for the selection of projections with SCDM

  IF (ionode) THEN  ! Read nnkp file on ionode only

     INQUIRE(file=trim(seedname)//".nnkp",exist=have_nnkp)
     IF(.not. have_nnkp) THEN
        CALL errore( 'pw2wannier90', 'Could not find the file '&
           &//trim(seedname)//'.nnkp', 1 )
     ENDIF

     iun_nnkp = find_free_unit()
     OPEN (unit=iun_nnkp, file=trim(seedname)//".nnkp",form='formatted', status="old")

  ENDIF

  nnbx=0

  !   check the information from *.nnkp with the nscf_save data
  WRITE(stdout,*) ' Checking info from wannier.nnkp file'
  WRITE(stdout,*)

  IF (ionode) THEN   ! read from ionode only

     CALL scan_file_to('real_lattice',found)
     if(.not.found) then
        CALL errore( 'pw2wannier90', 'Could not find real_lattice block in '&
           &//trim(seedname)//'.nnkp', 1 )
     endif
     DO j=1,3
        READ(iun_nnkp,*) (rlatt(i,j),i=1,3)
        DO i = 1,3
           rlatt(i,j) = rlatt(i,j)/(alat*bohr)
        ENDDO
     ENDDO
     DO j=1,3
        DO i=1,3
           IF(abs(rlatt(i,j)-at(i,j))>eps6) THEN
              WRITE(stdout,*)  ' Something wrong! '
              WRITE(stdout,*)  ' rlatt(i,j) =',rlatt(i,j),  ' at(i,j)=',at(i,j)
              CALL errore( 'pw2wannier90', 'Direct lattice mismatch', 3*j+i )
           ENDIF
        ENDDO
     ENDDO
     WRITE(stdout,*) ' - Real lattice is ok'

     CALL scan_file_to('recip_lattice',found)
     if(.not.found) then
        CALL errore( 'pw2wannier90', 'Could not find recip_lattice block in '&
           &//trim(seedname)//'.nnkp', 1 )
     endif
     DO j=1,3
        READ(iun_nnkp,*) (glatt(i,j),i=1,3)
        DO i = 1,3
           glatt(i,j) = (alat*bohr)*glatt(i,j)/tpi
        ENDDO
     ENDDO
     DO j=1,3
        DO i=1,3
           IF(abs(glatt(i,j)-bg(i,j))>eps6) THEN
              WRITE(stdout,*)  ' Something wrong! '
              WRITE(stdout,*)  ' glatt(i,j)=',glatt(i,j), ' bg(i,j)=',bg(i,j)
              CALL errore( 'pw2wannier90', 'Reciprocal lattice mismatch', 3*j+i )
           ENDIF
        ENDDO
     ENDDO
     WRITE(stdout,*) ' - Reciprocal lattice is ok'

     CALL scan_file_to('kpoints',found)
     if(.not.found) then
        CALL errore( 'pw2wannier90', 'Could not find kpoints block in '&
           &//trim(seedname)//'.nnkp', 1 )
     endif
     READ(iun_nnkp,*) numk
     IF(numk/=iknum) THEN
        WRITE(stdout,*)  ' Something wrong! '
        WRITE(stdout,*)  ' numk=',numk, ' iknum=',iknum
        CALL errore( 'pw2wannier90', 'Wrong number of k-points', numk)
     ENDIF
     IF(regular_mesh) THEN
        DO i=1,numk
           READ(iun_nnkp,*) xx(1), xx(2), xx(3)
           CALL cryst_to_cart( 1, xx, bg, 1 )
           IF(abs(xx(1)-xk(1,i))>eps6.or. &
                abs(xx(2)-xk(2,i))>eps6.or. &
                abs(xx(3)-xk(3,i))>eps6) THEN
              WRITE(stdout,*)  ' Something wrong! '
              WRITE(stdout,*) ' k-point ',i,' is wrong'
              WRITE(stdout,*) xx(1), xx(2), xx(3)
              WRITE(stdout,*) xk(1,i), xk(2,i), xk(3,i)
              CALL errore( 'pw2wannier90', 'problems with k-points', i )
           ENDIF
        ENDDO
     ENDIF ! regular mesh check
     WRITE(stdout,*) ' - K-points are ok'

  ENDIF ! ionode

  ! Broadcast
  CALL mp_bcast(rlatt,ionode_id, world_comm)
  CALL mp_bcast(glatt,ionode_id, world_comm)

  IF (ionode) THEN   ! read from ionode only
     if(noncolin) then
        old_spinor_proj=.false.
        CALL scan_file_to('spinor_projections',found)
        if(.not.found) then
           !try old style projections
           CALL scan_file_to('projections',found)
           if(found) then
              old_spinor_proj=.true.
           else
              CALL errore( 'pw2wannier90', 'Could not find projections block in '&
                 &//trim(seedname)//'.nnkp', 1 )
           endif
        end if
     else
        old_spinor_proj=.false.
        CALL scan_file_to('projections',found)
        if(.not.found) then
           CALL errore( 'pw2wannier90', 'Could not find projections block in '&
              &//trim(seedname)//'.nnkp', 1 )
        endif
     endif
     READ(iun_nnkp,*) n_proj
  ENDIF

  ! Broadcast
  CALL mp_bcast(n_proj,ionode_id, world_comm)
  CALL mp_bcast(old_spinor_proj,ionode_id, world_comm)

  IF(old_spinor_proj)THEN
  WRITE(stdout,'(//," ****** begin WARNING ****** ",/)')
  WRITE(stdout,'(" The pw.x calculation was done with non-collinear spin ")')
  WRITE(stdout,'(" but spinor = T was not specified in the wannier90 .win file!")')
  WRITE(stdout,'(" Please set spinor = T and rerun wannier90.x -pp  ")')
!  WRITE(stdout,'(/," If you are trying to reuse an old nnkp file, you can remove  ")')
!  WRITE(stdout,'(" this check from pw2wannir90.f90 line 870, and recompile. ")')
  WRITE(stdout,'(/," ******  end WARNING  ****** ",//)')
!  CALL errore("pw2wannier90","Spinorbit without spinor=T",1)
  ENDIF

  ! It is not clear if the next instruction is required or not, it probably depend
  ! on the version of wannier90 that was used to generate the nnkp file:
  IF(old_spinor_proj) THEN
     n_wannier=n_proj*2
  ELSE
     n_wannier=n_proj
  ENDIF

  ALLOCATE( center_w(3,n_proj), alpha_w(n_proj), &
       l_w(n_proj), mr_w(n_proj), r_w(n_proj), &
       zaxis(3,n_proj), xaxis(3,n_proj) )
  if(noncolin.and..not.old_spinor_proj) then
     ALLOCATE( spin_eig(n_proj),spin_qaxis(3,n_proj) )
  endif

  IF (ionode) THEN   ! read from ionode only
     DO iw=1,n_proj
        READ(iun_nnkp,*) (center_w(i,iw), i=1,3), l_w(iw), mr_w(iw), r_w(iw)
        READ(iun_nnkp,*) (zaxis(i,iw),i=1,3),(xaxis(i,iw),i=1,3),alpha_w(iw)
        xnorm = sqrt(xaxis(1,iw)*xaxis(1,iw) + xaxis(2,iw)*xaxis(2,iw) + &
             xaxis(3,iw)*xaxis(3,iw))
        IF (xnorm < eps6) CALL errore ('read_nnkp',' |xaxis| < eps ',1)
        znorm = sqrt(zaxis(1,iw)*zaxis(1,iw) + zaxis(2,iw)*zaxis(2,iw) + &
             zaxis(3,iw)*zaxis(3,iw))
        IF (znorm < eps6) CALL errore ('read_nnkp',' |zaxis| < eps ',1)
        coseno = (xaxis(1,iw)*zaxis(1,iw) + xaxis(2,iw)*zaxis(2,iw) + &
             xaxis(3,iw)*zaxis(3,iw))/xnorm/znorm
        IF (abs(coseno) > eps6) &
             CALL errore('read_nnkp',' xaxis and zaxis are not orthogonal !',1)
        IF (alpha_w(iw) < eps6) &
             CALL errore('read_nnkp',' zona value must be positive', 1)
        ! convert wannier center in cartesian coordinates (in unit of alat)
        CALL cryst_to_cart( 1, center_w(:,iw), at, 1 )
        if(noncolin.and..not.old_spinor_proj) then
           READ(iun_nnkp,*) spin_eig(iw),(spin_qaxis(i,iw),i=1,3)
           xnorm = sqrt(spin_qaxis(1,iw)*spin_qaxis(1,iw) + spin_qaxis(2,iw)*spin_qaxis(2,iw) + &
             spin_qaxis(3,iw)*spin_qaxis(3,iw))
           IF (xnorm < eps6) CALL errore ('read_nnkp',' |xaxis| < eps ',1)
           spin_qaxis(:,iw)=spin_qaxis(:,iw)/xnorm
        endif
     ENDDO
  ENDIF

  ! automatic projections
  IF (ionode) THEN
     CALL scan_file_to('auto_projections',found)
     IF (found) THEN
        READ (iun_nnkp, *) n_wannier
        READ (iun_nnkp, *) tmp_auto

        IF (scdm_proj) THEN
           IF (n_proj > 0) THEN
              WRITE(stdout,'(//, " ****** begin Error message ******",/)')
              WRITE(stdout,'(/," Found a projection block, an auto_projections block",/)')
              WRITE(stdout,'(/," and scdm_proj = T in the input file. These three options are inconsistent.",/)')
              WRITE(stdout,'(/," Please refer to the Wannier90 User guide for correct use of these flags.",/)')
              WRITE(stdout,'(/, " ****** end Error message ******",//)')
              CALL errore( 'pw2wannier90', 'Inconsistent options for projections.', 1 )
           ELSE
              IF (tmp_auto /= 0) CALL errore( 'pw2wannier90', 'Second entry in auto_projections block is not 0. ' // &
              'See Wannier90 User Guide in the auto_projections section for clarifications.', 1 )
           ENDIF
        ELSE
           ! Fire an error whether or not a projections block is found
           CALL errore( 'pw2wannier90', 'scdm_proj = F but found an auto_projections block in '&
                &//trim(seedname)//'.nnkp', 1 )
        ENDIF
     ELSE
        IF (scdm_proj) THEN
           ! Fire an error whether or not a projections block is found
           CALL errore( 'pw2wannier90', 'scdm_proj = T but cannot find an auto_projections block in '&
                &//trim(seedname)//'.nnkp', 1 )
        ENDIF
     ENDIF
  ENDIF

  ! Broadcast
  CALL mp_bcast(n_wannier,ionode_id, world_comm)
  CALL mp_bcast(center_w,ionode_id, world_comm)
  CALL mp_bcast(l_w,ionode_id, world_comm)
  CALL mp_bcast(mr_w,ionode_id, world_comm)
  CALL mp_bcast(r_w,ionode_id, world_comm)
  CALL mp_bcast(zaxis,ionode_id, world_comm)
  CALL mp_bcast(xaxis,ionode_id, world_comm)
  CALL mp_bcast(alpha_w,ionode_id, world_comm)
  if(noncolin.and..not.old_spinor_proj) then
     CALL mp_bcast(spin_eig,ionode_id, world_comm)
     CALL mp_bcast(spin_qaxis,ionode_id, world_comm)
  end if

  WRITE(stdout,'("  - Number of wannier functions is ok (",i3,")")') n_wannier

  IF (.not. scdm_proj) WRITE(stdout,*) ' - All guiding functions are given '
  !
  WRITE(stdout,*)
  WRITE(stdout,*) 'Projections:'
  DO iw=1,n_proj
     WRITE(stdout,'(3f12.6,3i3,f12.6)') &
          center_w(1:3,iw),l_w(iw),mr_w(iw),r_w(iw),alpha_w(iw)
  ENDDO

  IF (ionode) THEN   ! read from ionode only
     CALL scan_file_to('nnkpts',found)
     if(.not.found) then
        CALL errore( 'pw2wannier90', 'Could not find nnkpts block in '&
           &//trim(seedname)//'.nnkp', 1 )
     endif
     READ (iun_nnkp,*) nnb
  ENDIF

  ! Broadcast
  CALL mp_bcast(nnb,ionode_id, world_comm)
  !
  nnbx = max (nnbx, nnb )
  !
  ALLOCATE ( kpb(iknum,nnbx), g_kpb(3,iknum,nnbx),&
             ig_(iknum,nnbx), ig_check(iknum,nnbx) )
  ALLOCATE( zerophase(iknum,nnbx) )
  zerophase = .false.

  !  read data about neighbours
  WRITE(stdout,*)
  WRITE(stdout,*) ' Reading data about k-point neighbours '
  WRITE(stdout,*)

  IF (ionode) THEN
     DO ik=1, iknum
        DO ib = 1, nnb
           READ(iun_nnkp,*) idum, kpb(ik,ib), (g_kpb(ipol,ik,ib), ipol =1,3)
        ENDDO
     ENDDO
  ENDIF

  ! Broadcast
  CALL mp_bcast(kpb,ionode_id, world_comm)
  CALL mp_bcast(g_kpb,ionode_id, world_comm)

  DO ik=1, iknum
     DO ib = 1, nnb
        IF ( (g_kpb(1,ik,ib).eq.0) .and.  &
             (g_kpb(2,ik,ib).eq.0) .and.  &
             (g_kpb(3,ik,ib).eq.0) ) zerophase(ik,ib) = .true.
        g_(:) = REAL( g_kpb(:,ik,ib) )
        CALL cryst_to_cart (1, g_, bg, 1)
        gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
        ig_(ik,ib) = 0
        ig = 1
        DO WHILE  (gg(ig) <= gg_ + eps6)
           IF ( (abs(g(1,ig)-g_(1)) < eps6) .and.  &
                (abs(g(2,ig)-g_(2)) < eps6) .and.  &
                (abs(g(3,ig)-g_(3)) < eps6)  ) ig_(ik,ib) = ig
           ig= ig +1
        ENDDO
     ENDDO
  ENDDO
  ig_check(:,:) = ig_(:,:)
  CALL mp_sum( ig_check, intra_pool_comm )
  DO ik=1, iknum
     DO ib = 1, nnb
        IF (ig_check(ik,ib) ==0) &
          CALL errore('read_nnkp', &
                      ' g_kpb vector is not in the list of Gs', 100*ik+ib )
     ENDDO
  ENDDO
  DEALLOCATE (ig_check)

  WRITE(stdout,*) ' All neighbours are found '
  WRITE(stdout,*)

  ALLOCATE( excluded_band(nbnd) )

  IF (ionode) THEN     ! read from ionode only
     CALL scan_file_to('exclude_bands',found)
     if(.not.found) then
        CALL errore( 'pw2wannier90', 'Could not find exclude_bands block in '&
           &//trim(seedname)//'.nnkp', 1 )
     endif
     READ (iun_nnkp,*) nexband
     excluded_band(1:nbnd)=.false.
     DO i=1,nexband
        READ(iun_nnkp,*) indexb
        IF (indexb<1 .or. indexb>nbnd) &
             CALL errore('read_nnkp',' wrong excluded band index ', 1)
        excluded_band(indexb)=.true.
     ENDDO
  ENDIF
  num_bands=nbnd-nexband

  ! Broadcast
  CALL mp_bcast(nexband,ionode_id, world_comm)
  CALL mp_bcast(excluded_band,ionode_id, world_comm)
  CALL mp_bcast(num_bands,ionode_id, world_comm)

  IF (ionode) CLOSE (iun_nnkp)   ! ionode only

  RETURN
END SUBROUTINE read_nnkp
!
!-----------------------------------------------------------------------
SUBROUTINE scan_file_to( keyword, found )
   !-----------------------------------------------------------------------
   !
   USE wannier, ONLY :iun_nnkp
   USE io_global,  ONLY : stdout
   IMPLICIT NONE
   CHARACTER(len=*), intent(in) :: keyword
   logical, intent(out) :: found
   CHARACTER(len=80) :: line1, line2
!
! by uncommenting the following line the file scan restarts every time
! from the beginning thus making the reading independent on the order
! of data-blocks
!   rewind (iun_nnkp)
!
10 CONTINUE
   READ(iun_nnkp,*,end=20) line1, line2
   IF(line1/='begin')  GOTO 10
   IF(line2/=keyword) GOTO 10
   found=.true.
   RETURN
20 found=.false.
   rewind (iun_nnkp)

END SUBROUTINE scan_file_to
!
!
!-------------------------------------------------------------------------
SUBROUTINE get_wannier_to_plot( )
  !-----------------------------------------------------------------------
  !
  ! ... gets the list of Wannier functions to plot
  !
  USE kinds,           ONLY : DP
  USE wannier,         ONLY : n_wannier, iknum, wannier_plot_list, wann_to_plot
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=1), EXTERNAL :: capital
  !
  CHARACTER(LEN=10), PARAMETER :: c_digit = '0123456789'
  CHARACTER(LEN=2), PARAMETER :: c_range = '-:'
  CHARACTER(LEN=3), PARAMETER :: c_sep = ' ,;'
  CHARACTER(LEN=5), PARAMETER :: c_punc = ' ,:-:'
  CHARACTER(LEN=5), PARAMETER :: c_punc_nospace = ',:-:'
  CHARACTER(LEN=15), PARAMETER :: allchar = '0123456789 ,:-:'
  !
  CHARACTER(LEN=255) :: c_wlist, c_iwann
  INTEGER :: i, check, iwann, range_size, counter
  INTEGER :: i_punc, i_digit
  INTEGER :: nwannx
  INTEGER, ALLOCATABLE :: aux(:)
  !
  !
  nwannx = n_wannier * iknum       ! num WFs in the supercell
  ALLOCATE( aux(nwannx) )
  !
  c_wlist = wannier_plot_list
  !
  IF ( len_trim(c_wlist) == 0 ) &
    CALL errore( 'get_wannier_to_plot', 'wannier_plot_list is blank', 1 )
  !
  DO i = 1, len_trim(c_wlist)
    c_wlist(i:i) = capital(c_wlist(i:i))
  ENDDO
  !
  IF ( trim(c_wlist) == 'ALL' ) THEN
    ALLOCATE( wann_to_plot(nwannx) )
    DO i = 1, nwannx
      wann_to_plot(i) = i
    ENDDO
    !
    RETURN
  ELSE
    check = VERIFY( c_wlist, allchar )
    IF ( check .ne. 0 ) CALL errore( 'get_wannier_to_plot', &
                'Unrecognised character in wannier_plot_list', check ) 
  ENDIF
  !
  c_wlist = ADJUSTL( c_wlist )
  !
  IF ( SCAN( c_wlist, c_punc ) == 0 ) THEN
    READ( c_wlist, *, ERR=101, END=101 ) iwann
    IF ( iwann > nwannx ) &
      CALL errore( 'get_wannier_to_plot', 'wannier_plot_list out of range', iwann )
    !
    ALLOCATE( wann_to_plot(1) )
    wann_to_plot(1) = iwann
    !
    RETURN
  ENDIF
  !
  !
  counter = 0
  DO
    i_punc = SCAN( c_wlist, c_punc )
    !
    IF ( i_punc == 1 ) & 
      CALL errore( 'get_wannier_to_plot', 'Error parsing keyword wannier_plot_list', 2 )
    !
    counter = counter + 1
    c_iwann = c_wlist(1:i_punc-1)
    READ( c_iwann, *, ERR=101, END=101 ) iwann
    aux(counter) = iwann
    IF ( iwann > nwannx ) &
      CALL errore( 'get_wannier_to_plot', 'wannier_plot_list out of range', iwann )
    !
    c_wlist = ADJUSTL( c_wlist(i_punc:) )
    !
    IF ( SCAN( c_wlist, c_range ) == 1 ) THEN
      i_digit = SCAN( c_wlist, c_digit )
      IF ( SCAN( ADJUSTL( c_wlist(2:i_digit) ), c_punc_nospace ) /= 0 ) &
        CALL errore( 'get_wannier_to_plot', 'Error parsing keyword wannier_plot_list', 3 )
      c_wlist = ADJUSTL( c_wlist(i_digit:) )
      i_punc = SCAN( c_wlist, c_punc )
      !
      c_iwann = c_wlist(1:i_punc-1)
      READ( c_iwann, *, ERR=101, END=101 ) iwann
      IF ( iwann > nwannx ) &
        CALL errore( 'get_wannier_to_plot', 'wannier_plot_list out of range', iwann )
      !
      range_size = iwann - aux(counter)
      IF ( range_size <= 0 ) CALL errore( 'get_wannier_to_plot', &
                          'Error parsing keyword wannier_plot_list: incorrect range', 1)
      !
      DO i = 1, range_size
        counter = counter + 1
        aux(counter) = aux(counter-1) + 1
      ENDDO
      !
      c_wlist = ADJUSTL( c_wlist(i_punc:) )
    ENDIF
    !
    IF ( SCAN( c_wlist, c_sep ) == 1 ) c_wlist = ADJUSTL( c_wlist(2:) )
    IF ( SCAN( c_wlist, c_range ) == 1 ) &
      CALL errore( 'get_wannier_to_plot', 'Error parsing keyword wannier_plot_list', 4 )
    !
    IF ( INDEX( c_wlist, ' ' ) == 1 ) EXIT
  ENDDO
  !
  ALLOCATE( wann_to_plot(counter) )
  wann_to_plot(:) = aux(1:counter)
  !
  RETURN
  !
101 CALL errore( 'get_wannier_to_plot', 'Error parsing keyword wannier_plot_list', 1 ) 
  !
  ! 
END SUBROUTINE get_wannier_to_plot
!
!
!----------------------------------------------------------------------------
SUBROUTINE mp_grid_ks2kcp( )
   !---------------------------------------------------------------------------------
   !
   ! ...  This routine generate mp_grid for the ks2kcp mode.
   ! ...  It is necessary to momentarily change the definition of
   ! ...  iknum in order to properly define mp_grid.
   !
   USE wannier,           ONLY : kpt_latt, iknum
   USE lsda_mod,          ONLY : nspin
   USE cell_base,         ONLY : at
   USE klist,             ONLY : xk
   !
   !
   IMPLICIT NONE
   !
   LOGICAL :: ks_only
   !
   !
   iknum = iknum / nspin         ! momentarily change the value of iknum (needed by find_mp_grid)
   ALLOCATE( kpt_latt(3,iknum) )
   kpt_latt(:,1:iknum) = xk(:,1:iknum)
   CALL cryst_to_cart( iknum, kpt_latt, at, -1 )
   CALL find_mp_grid( )
   iknum = iknum * nspin         ! restore the initial value of iknum
   !
   !
 END SUBROUTINE mp_grid_ks2kcp
