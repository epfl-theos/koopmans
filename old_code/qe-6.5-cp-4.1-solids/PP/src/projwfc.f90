!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM do_projwfc
  !-----------------------------------------------------------------------
  !
  ! projects wavefunctions onto orthogonalized atomic wavefunctions,
  ! calculates Lowdin charges, spilling parameter, projected DOS
  ! or computes the LDOS in a volume given in input as function of energy
  !
  ! See files INPUT_PROJWFC.* in Doc/ directory for usage
  ! IMPORTANT: since v.5 namelist name is &projwfc and no longer &inputpp
  !
  USE parameters, ONLY : npk
  USE constants,  ONLY : rytoev
  USE kinds,      ONLY : DP
  USE klist,      ONLY : nks, nkstot, xk, degauss, ngauss, lgauss, ltetra
  USE io_files,   ONLY : prefix, tmp_dir
  USE noncollin_module, ONLY : noncolin
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE environment,ONLY : environment_start, environment_end
  USE mp,         ONLY : mp_bcast
  USE mp_global,  ONLY : mp_startup
  USE mp_world,   ONLY : world_comm
  USE mp_images,  ONLY : intra_image_comm
  USE mp_pools,   ONLY : intra_pool_comm
  USE mp_bands,   ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp_diag,    ONLY : mp_start_diag
  USE command_line_options, ONLY : ndiag_
  USE spin_orb,   ONLY : lforcet
  USE wvfct,      ONLY : et, nbnd
  USE paw_variables, ONLY : okpaw
  ! following modules needed for generation of tetrahedra
  USE ktetra,     ONLY : tetra, tetra_type, opt_tetra_init
  USE symm_base,  ONLY : nsym, s, time_reversal, t_rev
  USE cell_base,  ONLY : at, bg
  USE start_k,    ONLY : k1, k2, k3, nk1, nk2, nk3
  USE lsda_mod,   ONLY : lsda
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER (len=256) :: filpdos, filproj, outdir
  REAL (DP), allocatable :: xk_collect(:,:)
  REAL (DP) :: Emin, Emax, DeltaE, degauss1, ef_0
  INTEGER :: nks2, ngauss1, ios
  LOGICAL :: lwrite_overlaps, lbinary_data, needwf = .TRUE.
  LOGICAL :: lsym, kresolveddos, tdosinboxes, plotboxes, pawproj
  INTEGER, PARAMETER :: N_MAX_BOXES = 999
  INTEGER :: n_proj_boxes, irmin(3,N_MAX_BOXES), irmax(3,N_MAX_BOXES)
  LOGICAL :: lgww  !if .true. use GW QP energies from file bands.dat
  !
  NAMELIST / projwfc / outdir, prefix, ngauss, degauss, lsym, &
             Emin, Emax, DeltaE, filpdos, filproj, lgww, &
             kresolveddos, tdosinboxes, n_proj_boxes, irmin, irmax, plotboxes, &
             lwrite_overlaps, lbinary_data, pawproj, lforcet, ef_0
  !
  ! initialise environment
  !
  CALL mp_startup ( )
  CALL mp_start_diag ( ndiag_, world_comm, intra_bgrp_comm, &
          do_distr_diag_inside_bgrp_ = .true. )
  CALL set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, &
       inter_bgrp_comm )
  !
  CALL environment_start ( 'PROJWFC' )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filproj= ' '
  filpdos= ' '
  Emin   =-1000000.d0
  Emax   =+1000000.d0
  DeltaE = 0.01d0
  ngauss = 0
  lsym   = .true.
  degauss= 0.d0
  lgww   = .false.
  pawproj= .false.
  lwrite_overlaps   = .false.
  lbinary_data = .false.
  kresolveddos = .false.
  tdosinboxes = .false.
  plotboxes   = .false.
  n_proj_boxes= 1
  irmin(:,:)  = 1
  irmax(:,:)  = 0
  ef_0 = 0.d0
  lforcet = .false.
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, projwfc, iostat = ios)
     !
     tmp_dir = trimcheck (outdir)
     ! save the value of degauss and ngauss: they are read from file
     degauss1=degauss
     ngauss1 = ngauss
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id, intra_image_comm )
  IF (ios /= 0) CALL errore ('do_projwfc', 'reading projwfc namelist', abs (ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir,   ionode_id, intra_image_comm )
  CALL mp_bcast( prefix,    ionode_id, intra_image_comm )
  CALL mp_bcast( filproj,   ionode_id, intra_image_comm )
  CALL mp_bcast( ngauss1,   ionode_id, intra_image_comm )
  CALL mp_bcast( degauss1,  ionode_id, intra_image_comm )
  CALL mp_bcast( DeltaE,    ionode_id, intra_image_comm )
  CALL mp_bcast( lsym,      ionode_id, intra_image_comm )
  CALL mp_bcast( Emin,      ionode_id, intra_image_comm )
  CALL mp_bcast( Emax,      ionode_id, intra_image_comm )
  CALL mp_bcast( lwrite_overlaps, ionode_id, intra_image_comm )
  CALL mp_bcast( lbinary_data,    ionode_id, intra_image_comm )
  CALL mp_bcast( lgww,      ionode_id, intra_image_comm )
  CALL mp_bcast( pawproj,   ionode_id, intra_image_comm )
  CALL mp_bcast( tdosinboxes,     ionode_id, intra_image_comm )
  CALL mp_bcast( n_proj_boxes,    ionode_id, intra_image_comm )
  CALL mp_bcast( irmin,     ionode_id, intra_image_comm )
  CALL mp_bcast( irmax,     ionode_id, intra_image_comm )
  CALL mp_bcast( ef_0, ionode_id, intra_image_comm )
  CALL mp_bcast( lforcet, ionode_id, intra_image_comm )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file_new ( needwf )
  !
  IF(lgww) CALL get_et_from_gww ( nbnd, et )
  !
  ! Input checks 
  !
  IF (pawproj) THEN
    IF ( .NOT. okpaw ) CALL errore ('projwfc','option pawproj only for PAW',1)
    IF ( noncolin )  CALL errore ('projwfc','option pawproj and noncolinear spin not implemented',2)
    IF ( lforcet ) CALL errore ('projwfc','incompatible options',1)
    IF ( tdosinboxes ) CALL errore ('projwfc','incompatible options',2)
  END IF
  IF ( lforcet .AND. tdosinboxes ) CALL errore ('projwfc','incompatible options',3)
  IF ( lforcet .AND. lsym ) CALL errore ('projwfc','incompatible options',4)
  !
  !   Tetrahedron method
  !
  IF ( ltetra .AND. degauss1==0.d0 ) THEN
     !
     ! info on tetrahedra is no longer saved to file and must be rebuilt
     !
     ! in the lsda case, only the first half of the k points
     ! are needed in the input of "tetrahedra"
     !
     IF ( lsda ) THEN
        nks2 = nkstot / 2
     ELSE
        nks2 = nkstot
     END IF
     IF(tetra_type < 2) THEN
        ! use linear tetrahedron for both tetra_type=0 and tetra_type=1
        IF ( tetra_type == 0 ) tetra_type = 1
        WRITE( stdout,'(/5x,"Linear tetrahedron method (read from file) ")')
     ELSE
        WRITE( stdout,'(/5x,"Optimized tetrahedron method (read from file) ")')
     END IF
     !
     ! not sure this is needed
     !
     ALLOCATE(xk_collect(3,nkstot))
     CALL poolcollect(3, nks, xk, nkstot, xk_collect)
     !
     CALL opt_tetra_init(nsym, s, time_reversal, t_rev, at, bg, npk, k1,k2,k3, &
          &              nk1, nk2, nk3, nks2, xk_collect, 1)
     !
     DEALLOCATE(xk_collect)
     lgauss = .FALSE.
     !
  ELSE IF (degauss1/=0.d0) THEN
     degauss=degauss1
     ngauss =ngauss1
     WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.true.
     ltetra=.false.
  ELSE IF (lgauss) THEN
     WRITE( stdout,'(/5x,"Gaussian broadening (read from file): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     ltetra=.false.
  ELSE
     degauss=DeltaE/rytoev
     ngauss =0
     WRITE( stdout,'(/5x,"Gaussian broadening (default values): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.true.
     ltetra=.false.
  ENDIF
  !
  IF ( filpdos == ' ') filpdos = prefix
  !
  ! Now compute projections
  !
  IF ( tdosinboxes ) THEN
     CALL projwave_boxes (filpdos, filproj, n_proj_boxes, irmin, irmax, plotboxes)
  ELSE IF ( pawproj ) THEN
     CALL projwave_paw (filproj)
  ELSE
     CALL projwave(filproj, lsym, lwrite_overlaps, lbinary_data )
     IF ( lforcet ) CALL force_theorem ( ef_0, filproj )
  ENDIF
  !
  IF ( ionode .AND. .NOT. lforcet ) THEN
     IF ( tdosinboxes ) THEN
        CALL partialdos_boxes (Emin, Emax, DeltaE, kresolveddos, filpdos, n_proj_boxes)
     ELSE IF ( lsym .OR. kresolveddos ) THEN
        IF (noncolin) THEN
           CALL partialdos_nc (Emin, Emax, DeltaE, kresolveddos, filpdos)
        ELSE
           CALL partialdos (Emin, Emax, DeltaE, kresolveddos, filpdos)
        ENDIF
     ENDIF
  ENDIF
  !
  CALL laxlib_free_ortho_group()
  CALL environment_end ( 'PROJWFC' )
  !
  CALL stop_pp
  !
END PROGRAM do_projwfc

SUBROUTINE get_et_from_gww ( nbnd, et )
  !
  USE kinds, ONLY : dp
  USE constants, ONLY: rytoev
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: nbnd
  REAL(dp), INTENT(OUT):: et(nbnd,1)
  !
  INTEGER :: iun, idum, i, ios
  REAL(DP) :: rdum1, rdum2, rdum3
  LOGICAL :: lex
  INTEGER, EXTERNAL :: find_free_unit
  !
  INQUIRE ( file='bands.dat', EXIST=lex )
  WRITE(stdout,*) 'lex=', lex
  FLUSH(stdout)
  !
  IF(lex) THEN
     WRITE(stdout,*) 'Read the file bands.dat => GWA Eigenvalues used.'
     FLUSH(stdout)
     iun = find_free_unit()
     OPEN(unit=iun, file='bands.dat', status='unknown', form='formatted', &
          IOSTAT=ios)
     READ(iun,*) idum
     DO i=1, nbnd
        READ(iun,*) idum,rdum1,rdum2,et(i,1),rdum3
     ENDDO
     et(:,1)=et(:,1)/rytoev !! in bands.dat file, the QP energies are in eV
  ELSE
     WRITE(stdout,*) 'The file bands.dat does not exist.'
     WRITE(stdout,*) 'Eigenergies are not modified'
     FLUSH(stdout)
  ENDIF
END SUBROUTINE get_et_from_gww
!
SUBROUTINE write_lowdin ( filproj, nat, lmax_wfc, nspin, charges, charges_lm )
  !
  USE kinds, ONLY : dp
  USE io_global, ONLY : stdout, ionode
  USE klist, ONLY: nelec
  !
  IMPLICIT NONE
  !
  CHARACTER (len=*), INTENT(in) :: filproj
  INTEGER, INTENT(IN) :: nat, lmax_wfc, nspin
  REAL(DP), INTENT(in) :: charges (nat, 0:lmax_wfc, nspin )
  REAL(DP), INTENT(in), OPTIONAL :: charges_lm (nat, 0:lmax_wfc, 1:2*lmax_wfc+1, nspin )
  INTEGER, EXTERNAL :: find_free_unit
  !
  CHARACTER(len=256) :: filename
  INTEGER :: is, l, m, na, unit
  REAL(DP) :: totcharge(2), psum
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/)
  CHARACTER (len=7)  :: lm_label(1:7,1:3)=reshape( (/ &
    'z      ','x      ','y      ','       ','       ','       ','       ', &
    'z2     ','xz     ','yz     ','x2-y2  ','xy     ','       ','       ', &
    'z3     ','xz2    ','yz2    ','zx2-zy2','xyz    ','x3-3xy2','3yx2-y3' /), (/7,3/) )
  !
  IF ( TRIM(filproj) == ' ') THEN
     filename='lowdin.txt'
  ELSE
     filename = trim(filproj)//'.lowdin'
  END IF
  !
  IF ( ionode ) THEN
     unit = find_free_unit()
     OPEN( unit=unit, file=trim(filename), status='unknown', form='formatted')
     WRITE( stdout, '(/"Lowdin Charges: "/)')
     WRITE( unit, '(/"Lowdin Charges: "/)')
     !
     DO na = 1, nat
        DO is = 1, nspin
          totcharge(is) = SUM(charges(na,0:lmax_wfc,is))
        ENDDO
        IF ( nspin == 1) THEN
           DO l = 0, lmax_wfc
              WRITE(stdout, 2000,advance='no') na, totcharge(1), l_label(l), charges(na,l,1)
              WRITE(unit, 2000,advance='no') na, totcharge(1), l_label(l), charges(na,l,1)
              IF (l /= 0 .AND. present(charges_lm)) THEN
                 DO m = 1, 2*l+1
                    WRITE( stdout,'(A1,A,"=",F8.4,", ")',advance='no') &
                       l_label(l), trim(lm_label(m,l)), charges_lm(na,l,m,1)
                    WRITE( unit,'(A1,A,"=",F8.4,", ")',advance='no') &
                       l_label(l), trim(lm_label(m,l)), charges_lm(na,l,m,1)
                 ENDDO
              ENDIF
              WRITE(stdout,*)
              WRITE(unit,*)
           ENDDO
        ELSEIF ( nspin == 2) THEN
           WRITE( stdout, 2000) na, totcharge(1) + totcharge(2), &
                ( l_label(l), charges(na,l,1) + charges(na,l,2), l=0,lmax_wfc)
          WRITE( unit, 2000) na, totcharge(1) + totcharge(2), &
                ( l_label(l), charges(na,l,1) + charges(na,l,2), l=0,lmax_wfc)
           DO l = 0, lmax_wfc
              WRITE(stdout,2001,advance='no') totcharge(1), l_label(l), charges(na,l,1)
              WRITE(unit,2001,advance='no') totcharge(1), l_label(l), charges(na,l,1)
              IF (l /= 0 .AND. present(charges_lm)) THEN
                 DO m = 1, 2*l+1
                    WRITE( stdout,'(A1,A,"=",F8.4,", ")',advance='no') &
                       l_label(l), trim(lm_label(m,l)), charges_lm(na,l,m,1)
                    WRITE( unit,'(A1,A,"=",F8.4,", ")',advance='no') &
                       l_label(l), trim(lm_label(m,l)), charges_lm(na,l,m,1)
                 ENDDO
              ENDIF
              WRITE(stdout,*)
              WRITE(unit,*)
           ENDDO
           DO l = 0, lmax_wfc
              WRITE(stdout,2002,advance='no') totcharge(2), l_label(l), charges(na,l,2)
              WRITE(unit,2002,advance='no') totcharge(2), l_label(l), charges(na,l,2)
              IF (l /= 0 .AND. present(charges_lm)) THEN
                 DO m = 1, 2*l+1
                    WRITE( stdout,'(A1,A,"=",F8.4,", ")',advance='no') &
                       l_label(l), trim(lm_label(m,l)), charges_lm(na,l,m,2)
                    WRITE( unit,'(A1,A,"=",F8.4,", ")',advance='no') &
                       l_label(l), trim(lm_label(m,l)), charges_lm(na,l,m,2)
                 ENDDO
              ENDIF
              WRITE(stdout,*)
              WRITE(unit,*)
           ENDDO
           WRITE( stdout, 2003) totcharge(1) - totcharge(2), &
                ( l_label(l), charges(na,l,1) - charges(na,l,2), l=0,lmax_wfc)
           WRITE( unit, 2003) totcharge(1) - totcharge(2), &
                ( l_label(l), charges(na,l,1) - charges(na,l,2), l=0,lmax_wfc)
        ENDIF
     ENDDO
     !
     psum = SUM(charges(:,:,:)) / nelec
     WRITE( stdout, '(5x,"Spilling Parameter: ",f8.4)') 1.0d0 - psum
     WRITE( unit, '(5x,"Spilling Parameter: ",f8.4)') 1.0d0 - psum
     !
     ! Sanchez-Portal et al., Sol. State Commun.  95, 685 (1995).
     ! The spilling parameter measures the ability of the basis provided by
     ! the pseudo-atomic wfc to represent the PW eigenstates,
     ! by measuring how much of the subspace of the Hamiltonian
     ! eigenstates falls outside the subspace spanned by the atomic basis
     !
     CLOSE(unit)
  END IF

2000 FORMAT (5x,"Atom # ",i3,": total charge = ",f8.4,4(", ",a1," =",f8.4))
2001 FORMAT (15x,"  spin up      = ",f8.4,4(", ",a1," =",f8.4))
2002 FORMAT (15x,"  spin down    = ",f8.4,4(", ",a1," =",f8.4))
2003 FORMAT (15x,"  polarization = ",f8.4,4(", ",a1," =",f8.4))

END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE sym_proj_g (rproj0, proj_out)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE basis,      ONLY : natomwfc
  USE wvfct,      ONLY : nbnd
  USE symm_base,  ONLY : nsym, irt, t_rev, d1, d2, d3
  USE projections,ONLY : nlmchi
  !
  IMPLICIT NONE
  REAL(DP),    INTENT(IN) ::rproj0 (natomwfc, nbnd)
  REAL   (DP), INTENT(OUT):: proj_out(natomwfc, nbnd)
  !
  INTEGER :: na, nb, n, l, m, m1, isym, nwfc, nwfc1, ibnd
  REAL   (DP), ALLOCATABLE :: rwork1(:)
  !
  ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1
  !
  CALL d_matrix (d1, d2, d3)
  proj_out(:,:) = 0.d0
  !
  ALLOCATE(rwork1 (nbnd) )
  !
  DO nwfc = 1, natomwfc
     !
     !  atomic wavefunction nwfc is on atom na
     !
     na= nlmchi(nwfc)%na
     n = nlmchi(nwfc)%n
     l = nlmchi(nwfc)%l
     m = nlmchi(nwfc)%m
     !
     DO isym = 1, nsym
        nb = irt (isym, na)
        DO nwfc1 =1, natomwfc
           IF (nlmchi(nwfc1)%na == nb             .and. &
                nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. &
                nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. &
                nlmchi(nwfc1)%m == 1 ) GOTO 10
        ENDDO
        CALL errore('sym_proj_g','cannot symmetrize',1)
10      nwfc1=nwfc1-1
        !
        !  nwfc1 is the first rotated atomic wfc corresponding to nwfc
        !
        IF (l == 0) THEN
           rwork1(:) = rproj0 (nwfc1 + 1,:)
        ELSEIF (l == 1) THEN
           rwork1(:) = 0.d0
           DO m1 = 1, 3
              rwork1(:)=rwork1(:)+d1(m1,m,isym)*rproj0(nwfc1+m1,:)
           ENDDO
        ELSEIF (l == 2) THEN
           rwork1(:) = 0.d0
           DO m1 = 1, 5
              rwork1(:)=rwork1(:)+d2(m1,m,isym)*rproj0(nwfc1+m1,:)
           ENDDO
        ELSEIF (l == 3) THEN
           rwork1(:) = 0.d0
           DO m1 = 1, 7
              rwork1(:)=rwork1(:)+d3(m1,m,isym)*rproj0(nwfc1+m1,:)
           ENDDO
        ENDIF
        DO ibnd = 1, nbnd
           proj_out (nwfc, ibnd ) = proj_out (nwfc, ibnd ) + &
                rwork1(ibnd) * rwork1(ibnd) / nsym
        ENDDO
     ENDDO
  ENDDO
  !
  DEALLOCATE(rwork1 )
  !  
END SUBROUTINE sym_proj_g
!
!-----------------------------------------------------------------------
SUBROUTINE sym_proj_k (proj0, proj_out)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE basis,      ONLY : natomwfc
  USE wvfct,      ONLY : nbnd
  USE symm_base,  ONLY : nsym, irt, t_rev, d1, d2, d3
  USE projections,ONLY : nlmchi
  !
  IMPLICIT NONE
  COMPLEX(DP), INTENT(IN) :: proj0 (natomwfc, nbnd)
  REAL   (DP), INTENT(OUT):: proj_out(natomwfc, nbnd)
  !
  INTEGER :: na, nb, n, l, m, m1, isym, nwfc, nwfc1, ibnd
  COMPLEX(DP), ALLOCATABLE ::  work1(:)
  !
  ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1
  !
  CALL d_matrix (d1, d2, d3)
  proj_out(:,:) = 0.d0
  !
  ALLOCATE(work1 (nbnd) )
  !
  DO nwfc = 1, natomwfc
     !
     !  atomic wavefunction nwfc is on atom na
     !
     na= nlmchi(nwfc)%na
     n = nlmchi(nwfc)%n
     l = nlmchi(nwfc)%l
     m = nlmchi(nwfc)%m
     !
     DO isym = 1, nsym
        nb = irt (isym, na)
        DO nwfc1 =1, natomwfc
           IF (nlmchi(nwfc1)%na == nb             .and. &
                nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. &
                nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. &
                nlmchi(nwfc1)%m == 1 ) GOTO 10
        ENDDO
        CALL errore('sym_proj_k','cannot symmetrize',1)
10      nwfc1=nwfc1-1
        !
        !  nwfc1 is the first rotated atomic wfc corresponding to nwfc
        !
        IF (l == 0) THEN
           work1(:) = proj0 (nwfc1 + 1,:)
        ELSEIF (l == 1) THEN
           work1(:) = 0.d0
           DO m1 = 1, 3
              work1(:)=work1(:)+d1(m1,m,isym)*proj0(nwfc1+m1,:)
           ENDDO
        ELSEIF (l == 2) THEN
           work1(:) = 0.d0
           DO m1 = 1, 5
              work1(:)=work1(:)+d2(m1,m,isym)*proj0(nwfc1+m1,:)
           ENDDO
        ELSEIF (l == 3) THEN
           work1(:) = 0.d0
           DO m1 = 1, 7
              work1(:)=work1(:)+d3(m1,m,isym)*proj0(nwfc1+m1,:)
           ENDDO
        ENDIF
        DO ibnd = 1, nbnd
           proj_out (nwfc, ibnd) = proj_out (nwfc, ibnd) + &
                work1(ibnd) * conjg (work1(ibnd)) / nsym
        ENDDO
     ENDDO
  ENDDO
  !
  DEALLOCATE(work1 )
  !  
END SUBROUTINE sym_proj_k
!
!-----------------------------------------------------------------------
SUBROUTINE sym_proj_so (domag, proj0, proj_out  )
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE basis,      ONLY : natomwfc
  USE wvfct,      ONLY : nbnd
  USE symm_base,  ONLY : nsym, irt, t_rev
  USE projections,ONLY : nlmchi
  !
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: domag
  COMPLEX(DP), INTENT(IN) :: proj0 (natomwfc, nbnd)
  REAL   (DP), INTENT(OUT):: proj_out(natomwfc, nbnd)
  !
  INTEGER :: na, nb, n, l, m, m1, ind, ind0, isym, nwfc, nwfc1, ibnd
  REAL(DP) :: jj
  COMPLEX(DP), ALLOCATABLE ::  work1(:)
  COMPLEX(DP) :: d12(2, 2, 48), d32(4, 4, 48), d52(6, 6, 48), d72(8, 8, 48)
  !
  ! initialize D_Sj for j=1/2, j=3/2, j=5/2 and j=7/2
  !
  CALL d_matrix_so (d12, d32, d52, d72)
  !
  proj_out(:,:) = 0.d0
  !
  ALLOCATE(work1 (nbnd) )
  !
  DO nwfc = 1, natomwfc
     !
     !  atomic wavefunction nwfc is on atom na
     !
     na= nlmchi(nwfc)%na
     n = nlmchi(nwfc)%n
     l = nlmchi(nwfc)%l
     m = nlmchi(nwfc)%m
     ind0 = nlmchi(nwfc)%ind
     jj = nlmchi(nwfc)%jj
     !
     DO isym = 1, nsym
        !-- check for the time reversal
        IF (t_rev(isym) == 1) THEN
           ind = 2*jj + 2 - ind0
        ELSE
           ind = ind0
        ENDIF
        !--
        nb = irt (isym, na)
        DO nwfc1 =1, natomwfc
           IF (nlmchi(nwfc1)%na == nb             .and. &
                nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. &
                nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. &
                nlmchi(nwfc1)%jj == nlmchi(nwfc)%jj .and. &
                nlmchi(nwfc1)%ind == 1 ) GOTO 10
        ENDDO
        CALL errore('sym_proj_so','cannot symmetrize',1)
10      nwfc1=nwfc1-1
        !
        !  nwfc1 is the first rotated atomic wfc corresponding to nwfc
        !
        IF (abs(jj-0.5d0)<1.d-8) THEN
           work1(:) = 0.d0
           DO m1 = 1, 2
              work1(:)=work1(:)+d12(m1,ind,isym)*proj0(nwfc1+m1,:)
           ENDDO
        ELSEIF (abs(jj-1.5d0)<1.d-8) THEN
           work1(:) = 0.d0
           DO m1 = 1, 4
              work1(:)=work1(:)+d32(m1,ind,isym)*proj0(nwfc1 + m1,:)
           ENDDO
        ELSEIF (abs(jj-2.5d0)<1.d-8) THEN
           work1(:) = 0.d0
           DO m1 = 1, 6
              work1(:)=work1(:)+d52(m1,ind,isym)*proj0(nwfc1+m1,:)
           ENDDO
        ELSEIF (abs(jj-3.5d0)<1.d-8) THEN
           work1(:) = 0.d0
           DO m1 = 1, 8
              work1(:)=work1(:)+d72(m1,ind,isym)*proj0(nwfc1+m1,:)
           ENDDO
        ENDIF
        DO ibnd = 1, nbnd
           proj_out (nwfc, ibnd) = proj_out (nwfc, ibnd) + &
                work1(ibnd) * conjg (work1(ibnd)) / nsym
        ENDDO
        ! on symmetries
        !--  in a nonmagnetic case - another loop with the time reversal
        IF ( .not.domag .and. ind==ind0 ) THEN
           ind = 2*jj + 2 - ind0
           nwfc1 = nwfc1 + 1
           GOTO 10
        ENDIF
     ENDDO
     !--  in a nonmagnetic case - rescale
     IF (.not.domag) THEN
        DO ibnd = 1, nbnd
           proj_out(nwfc,ibnd) = 0.5_dp*proj_out(nwfc,ibnd)
        ENDDO
     ENDIF
     ! on atomic wavefunctions
  END DO
  !
  DEALLOCATE (work1)
  !
END SUBROUTINE sym_proj_so
!-----------------------------------------------------------------------
SUBROUTINE sym_proj_nc ( proj0, proj_out  )
  !
  USE kinds,      ONLY : DP
  USE basis,      ONLY : natomwfc
  USE wvfct,      ONLY : nbnd
  USE symm_base,  ONLY : nsym, irt, t_rev
  USE projections,ONLY : nlmchi
  !
  IMPLICIT NONE
  COMPLEX(DP), INTENT(IN) :: proj0 (natomwfc, nbnd)
  REAL   (DP), INTENT(OUT):: proj_out(natomwfc, nbnd)
  !
  INTEGER :: na, nb, n, l, m, m1, ind, ind0, jj, isym, nwfc, nwfc1, ibnd
  COMPLEX(DP), ALLOCATABLE ::  work1(:)
  COMPLEX(DP) :: d012(2, 2, 48), d112(6, 6, 48), d212(10, 10, 48), &
       d312(14, 14, 48)
  !
  ! initialize D_Sl for l=0, l=1, l=2 and l=3
  !
  CALL d_matrix_nc (d012, d112, d212, d312)
  !
  proj_out(:,:) = 0.d0
  !
  ALLOCATE(work1 (nbnd) )
  !
  DO nwfc = 1, natomwfc
     na= nlmchi(nwfc)%na
     n = nlmchi(nwfc)%n
     l = nlmchi(nwfc)%l
     m = nlmchi(nwfc)%m
     ind0 = nlmchi(nwfc)%ind
     !
     DO isym = 1, nsym
        !-- check for the time reversal
        IF (t_rev(isym) == 1) THEN
           ind = 2*m - ind0 + 2*l + 1
        ELSE
           ind = ind0
        ENDIF
        !--
        nb = irt (isym, na)
        DO nwfc1 =1, natomwfc
           IF (nlmchi(nwfc1)%na == nb             .and. &
                nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. &
                nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. &
                nlmchi(nwfc1)%m == 1 .and. &
                nlmchi(nwfc1)%ind == 1) GOTO 15
        ENDDO
        CALL errore('sym_proj_nc','cannot symmetrize',1)
15      nwfc1=nwfc1-1
        IF (l == 0) THEN
           work1(:) = 0.d0
           DO m1 = 1, 2
              work1(:) = work1(:) + d012 (m1, ind, isym) * &
                   proj0 (nwfc1 + m1,:)
           ENDDO
        ELSEIF (l == 1) THEN
           work1(:) = 0.d0
           DO m1 = 1, 6
              work1(:) = work1(:) + d112 (m1, ind, isym) * &
                   proj0 (nwfc1 + m1,:)
           ENDDO
        ELSEIF (l == 2) THEN
           work1(:) = 0.d0
           DO m1 = 1, 10
              work1(:) = work1(:) + d212 (m1, ind, isym) * &
                   proj0 (nwfc1 + m1,:)
           ENDDO
        ELSEIF (l == 3) THEN
           work1(:) = 0.d0
           DO m1 = 1, 14
              work1(:) = work1(:) + d312 (m1, ind, isym) * &
                   proj0 (nwfc1 + m1,:)
           ENDDO
        ENDIF
        DO ibnd = 1, nbnd
           proj_out (nwfc, ibnd) = proj_out (nwfc, ibnd) + &
                work1(ibnd) * conjg (work1(ibnd)) / nsym
        ENDDO
        ! on symmetries
     ENDDO
     ! on atomic wavefunctions
  END DO
  !
  DEALLOCATE (work1)
  !
END SUBROUTINE sym_proj_nc
!-----------------------------------------------------------------------
SUBROUTINE write_proj ( lmax_wfc, filproj, proj )
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE constants,  ONLY : rytoev, eps4
  USE basis,      ONLY : natomwfc
  USE lsda_mod,   ONLY : nspin, isk, current_spin
  USE noncollin_module, ONLY : noncolin
  USE spin_orb,   ONLY : lspinorb
  USE klist,      ONLY : nkstot, xk
  USE ions_base,  ONLY : nat, ityp, atm
  USE wvfct,      ONLY : et, wg, nbnd
  USE projections,ONLY : nlmchi
  !
  IMPLICIT NONE
  INTEGER, INTENT(in) :: lmax_wfc
  CHARACTER (len=*), INTENT(in) :: filproj
  REAL(DP), INTENT(IN) :: proj(natomwfc,nbnd,nkstot)
  !
  INTEGER :: nspin0, nwfc, ibnd, i, j, ik, na, l, m
  INTEGER, ALLOCATABLE :: idx(:)
  REAL(DP) :: psum
  REAL(DP), ALLOCATABLE :: proj1 (:)
  REAL(DP), ALLOCATABLE :: charges(:,:,:), charges_lm(:,:,:,:)
  REAL (DP), EXTERNAL :: compute_mj
  CHARACTER (len=1) :: plus
  !
  INTERFACE 
     SUBROUTINE write_lowdin ( filproj, nat, lmax_wfc, nspin, charges, charges_lm )
       IMPORT  :: DP
       CHARACTER (len=*), INTENT(in) :: filproj
       INTEGER, INTENT(IN) :: nat, lmax_wfc, nspin
       REAL(DP), INTENT(in) :: charges (nat, 0:lmax_wfc, nspin )
       REAL(DP), INTENT(in), OPTIONAL :: charges_lm (nat, 0:lmax_wfc, 1:2*lmax_wfc+1, nspin )
     END SUBROUTINE write_lowdin
  END INTERFACE
  !
  !
  WRITE( stdout,'(/5x,"Atomic states used for projection")')
  WRITE( stdout,'( 5x,"(read from pseudopotential files):"/)')
  DO nwfc = 1, natomwfc
     WRITE(stdout,1000, ADVANCE="no" ) &
          nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
          nlmchi(nwfc)%n, nlmchi(nwfc)%l
     IF ( lspinorb ) THEN
        WRITE(stdout,1001) nlmchi(nwfc)%jj, &
             compute_mj(nlmchi(nwfc)%jj,nlmchi(nwfc)%l,nlmchi(nwfc)%m)
     ELSE IF ( noncolin ) THEN
        WRITE(stdout,1002) nlmchi(nwfc)%m, &
             0.5d0-int(nlmchi(nwfc)%ind/(2*nlmchi(nwfc)%l+2))
     ELSE
        WRITE(stdout,1003) nlmchi(nwfc)%m
     END IF
  ENDDO
1000 FORMAT (5x,"state #",i4,": atom ",i3," (",a3,"), wfc ",i2, &
          " (l=",i1)
1001 FORMAT (" j=",f3.1," m_j=",f4.1,")")
1002 FORMAT (" m=",i2," s_z=",f4.1,")")
1003 FORMAT (" m=",i2,")")
  !
  ALLOCATE(idx(natomwfc), proj1 (natomwfc) )
  DO ik = 1, nkstot
     WRITE( stdout, '(/" k = ",3f14.10)') (xk (i, ik) , i = 1, 3)
     DO ibnd = 1, nbnd
        WRITE( stdout, '("==== e(",i4,") = ",f11.5," eV ==== ")') &
             ibnd, et (ibnd, ik) * rytoev
        !
        ! sort projections by magnitude, in decreasing order
        !
        DO nwfc = 1, natomwfc
           idx (nwfc) = 0
           proj1 (nwfc) = - proj (nwfc, ibnd, ik)
        ENDDO
        !
        ! projections differing by less than 1.d-4 are considered equal
        !
        CALL hpsort_eps (natomwfc, proj1, idx, eps4)
        !
        !  only projections that are larger than 0.001 are written
        !
        DO nwfc = 1, natomwfc
           proj1 (nwfc) = - proj1(nwfc)
           IF ( abs (proj1(nwfc)) < 0.001d0 ) GOTO 20
        ENDDO
        nwfc = natomwfc + 1
20      nwfc = nwfc -1
        !
        ! fancy (?!?) formatting
        !
        plus='+'
        DO i = 1, nwfc
           IF ( i == 1 ) THEN
              WRITE( stdout,'(5X,"psi = ",f5.3,"*[#",i4,"]")',advance='no') &
                 proj1 (i), idx(i)
           ELSE
              WRITE( stdout,'(A,f5.3,"*[#",i4,"]")',advance='no') &
                 plus, proj1 (i), idx(i)
              IF ( MOD(i,5) == 0 ) WRITE( stdout,'(/,10X)', advance='no' )
           END IF
        ENDDO
        WRITE( stdout, * )
        psum = SUM ( proj(1:natomwfc, ibnd, ik) )
        WRITE( stdout, '(4x,"|psi|^2 = ",f5.3)') psum
        !
     ENDDO
  ENDDO
  DEALLOCATE (idx, proj1)
  !
  ! estimate partial charges (Loewdin) on each atom
  !
  IF ( lspinorb ) THEN
     nspin0 = 1
  ELSE IF ( noncolin ) THEN
     nspin0 = 2
  ELSE
     nspin0 = nspin
  END IF
  ALLOCATE ( charges (nat, 0:lmax_wfc, nspin0 ) )
  charges = 0.0d0
  IF ( nspin /= 4 ) THEN
     ALLOCATE ( charges_lm (nat, 0:lmax_wfc, 1:2*lmax_wfc+1, nspin ) )
     charges_lm = 0.d0
  END IF
  DO ik = 1, nkstot
     IF ( nspin == 2 ) THEN
        current_spin = isk ( ik )
     ELSE
        current_spin = 1
     ENDIF
     DO ibnd = 1, nbnd
        DO nwfc = 1, natomwfc
           na= nlmchi(nwfc)%na
           l = nlmchi(nwfc)%l
	   IF ( noncolin .AND. .NOT. lspinorb ) THEN
	      IF (nlmchi(nwfc)%ind<=(2*l+1)) THEN
	         current_spin = 1
              ELSE
	         current_spin = 2
              ENDIF
	   END IF
           charges(na,l,current_spin) = charges(na,l,current_spin) + &
                wg (ibnd,ik) * proj (nwfc, ibnd, ik)
           IF ( nspin /= 4 ) THEN
              m = nlmchi(nwfc)%m
              charges_lm(na,l,m,current_spin) = charges_lm(na,l,m,current_spin) + &
                wg (ibnd,ik) * proj (nwfc, ibnd, ik)
           END IF
        ENDDO
     ENDDO
  ENDDO
  !
  IF ( nspin /= 4 ) THEN
     CALL write_lowdin ( filproj, nat, lmax_wfc, nspin, charges, charges_lm )
     DEALLOCATE (charges_lm)
  ELSE
     CALL write_lowdin ( filproj, nat, lmax_wfc, nspin0, charges )
  END IF
  DEALLOCATE (charges)
  !
END SUBROUTINE write_proj
!
SUBROUTINE force_theorem ( ef_0, filproj )
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : rytoev
  USE io_global,  ONLY : ionode
  USE ions_base,  ONLY : nat, ityp, atm
  USE klist,      ONLY : nks, nkstot, xk
  USE basis,      ONLY : natomwfc
  USE wvfct,      ONLY : wg, et, nbnd
  USE mp,         ONLY : mp_sum
  USE mp_pools,   ONLY : inter_pool_comm, intra_pool_comm
  USE projections,ONLY : proj, nlmchi
  !
  !---- Force Theorem -- (AlexS)
  !
  IMPLICIT NONE
  CHARACTER (len=*), INTENT(in) :: filproj
  REAL(DP), INTENT(IN) :: ef_0
  !
  INTEGER :: ik, i, nwfc, na, l
  REAL(DP) :: eband_proj_tot, eband_tot, psum
  REAL(DP), ALLOCATABLE :: eband_proj(:)
  CHARACTER(len=256) :: filename
  !
  CALL weights()
  !   write(6,*) 'ef_0 = ', ef_0
  !   write(6,*) wg
  eband_tot = 0.d0
  ALLOCATE (eband_proj(natomwfc))
  eband_proj = 0.d0
  !
  !    loop on k points
  !
  DO ik = 1, nks
     !
     DO i = 1, nbnd
        psum = wg(i,ik) * (et(i,ik)-ef_0/rytoev)
        eband_tot = eband_tot + psum
        DO nwfc = 1, natomwfc
           eband_proj(nwfc) = eband_proj(nwfc) + psum*proj(nwfc,i,ik)
        ENDDO
     ENDDO
     !
  ENDDO
  !   
  CALL mp_sum( eband_tot,  inter_pool_comm )
  CALL mp_sum( eband_proj, inter_pool_comm )
  !
  IF ( ionode ) THEN
     !
     !-- Output for the Force Theorem (AlexS)
     !
     filename = trim(filproj)
     OPEN (4,file=filename,form='formatted', status='unknown')
     
     eband_proj_tot = 0.d0
     DO na = 1, nat
        
        psum  = 0.d0
        WRITE(4,*) 'Atom   ', na, atm(ityp(na))
        nwfc = 1
        DO WHILE (nwfc.LE.natomwfc)
           IF (nlmchi(nwfc)%na.eq.na) THEN
              l = nlmchi(nwfc)%l
              IF (l.eq.0)  THEN 
                 write(4,*) '... s_up, s_down'
              ELSEIF (l.eq.1) THEN 
                 write(4,*) '... {p_up}, {p_down}'
              ELSEIF (l.eq.2) THEN 
                 write(4,*) '... {d_up}, {d_down}'
              ELSEIF (l.eq.3) THEN 
                 write(4,*) '... {f_up}, {f_down}'
              ELSE
                 call errore('force_theorem','Force Theorem not implemented for l > 2',1)
              ENDIF
              DO i = 1, 2*l + 1
                 WRITE(4,'("eband_l_m_atom ",3i5,2e23.10)') l, i, na, &
                         eband_proj(nwfc-1+i)*rytoev, &
                         eband_proj(nwfc+i+2*l)*rytoev
                 psum  = psum+eband_proj(nwfc-1+i) +  &
                      eband_proj(nwfc+i+2*l)
              ENDDO
              nwfc = nwfc + 2*(2*l+1)
           ELSE
              nwfc = nwfc + 1
           ENDIF
        ENDDO
        eband_proj_tot = eband_proj_tot + psum
        WRITE(4,'("eband_atom (eV) = ",i5,e30.10)') na, psum*rytoev
        
        WRITE(4,*)
        
     ENDDO
     eband_tot = eband_tot*rytoev
     eband_proj_tot = eband_proj_tot*rytoev
     WRITE( 4,'(''eband_tot, eband_proj_tot (eV) = '',2e30.10)') eband_tot, eband_proj_tot
     
     CLOSE(4)
        
     DEALLOCATE (eband_proj)
     !
  ENDIF
  !
END SUBROUTINE FORCE_THEOREM
!--
!-----------------------------------------------------------------------
SUBROUTINE projwave_paw( filproj)
!    8/12/2014 N. A. W. Holzwarth -- attempt to calculate
!      charge within augmentation sphere for pdos
  !-----------------------------------------------------------------------
  !
  USE atom,       ONLY : rgrid, msh
  USE io_global, ONLY : stdout, ionode
  USE ions_base, ONLY : nat, ntyp => nsp, ityp
  USE constants, ONLY: rytoev
  USE klist, ONLY: xk, nks, nkstot, nelec, igk_k, ngk
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE wvfct, ONLY: npwx, nbnd
  USE uspp, ONLY: nkb, vkb
  USE uspp_param, ONLY : upf
  USE becmod,   ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files,  ONLY : restart_dir
  USE pw_restart_new,ONLY : read_collected_wfc
  USE wavefunctions, ONLY: evc
  !
  USE projections
  !
  IMPLICIT NONE
  !
  CHARACTER (len=*) :: filproj
  LOGICAL           :: lwrite_ovp, lbinary
  !
  INTEGER :: npw, ik, ibnd, i, j, k, na, nb, nt, isym, n,  m, l, nwfc,&
       nwfc1, is, ndm, mr,nbp
  REAL(DP), ALLOCATABLE :: e (:), aux(:), pcharge(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:)
  COMPLEX(DP), ALLOCATABLE :: overlap(:,:), work(:,:),work1(:), proj0(:,:)
  ! Some workspace for k-point calculation ...
  REAL   (DP), ALLOCATABLE ::roverlap(:,:), rwork1(:),rproj0(:,:)
  !
  !
  WRITE( stdout, '(/5x,"Calling projwave_paw .... ")')
  !
  !  NAWH 08/12/2014 -- need nkb functions for this case; must reflect
  !     vkb structure
  
  nwfc=0; mr=0
  do nt=1,ntyp
     nwfc=MAX(nwfc,upf(nt)%nbeta)
     mr=MAX(mr,upf(nt)%kkbeta)
  enddo

  ALLOCATE (pcharge(nwfc,nwfc,ntyp), aux(mr))
  pcharge=0.d0

  do nt=1,ntyp
     do i=1,upf(nt)%nbeta
        l=upf(nt)%lll(i)
        do j=1,upf(nt)%nbeta
           if (upf(nt)%lll(j)==l) then
              aux=0
              k=upf(nt)%kkbeta
              aux(1:k)=upf(nt)%aewfc(1:k,i)*upf(nt)%aewfc(1:k,j)
              call simpson(k,aux,rgrid(nt)%rab,pcharge(i,j,nt))
              write(6,*) "pcharge: ", i,j,l,pcharge(i,j,nt) 
           endif
        enddo
     enddo
  enddo     
     
  DEALLOCATE(aux)

  CALL fill_nlmbeta ( nkb, nwfc ) 

  ALLOCATE( proj (nkb, nbnd, nkstot),proj0(nkb,nbnd) )
  proj      = 0.d0
  !
  CALL allocate_bec_type (nkb, nbnd, becp )
  !
  !    loop on k points
  !
  DO ik = 1, nks
     CALL read_collected_wfc ( restart_dir() , ik, evc )
     npw = ngk(ik)
     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)

     proj0=0; 
     CALL calbec ( npw, vkb, evc, proj0)

     do nwfc=1,nkb
        na=nlmchi(nwfc)%na
        nt=ityp(na)
        nb=nlmchi(nwfc)%n
        l=nlmchi(nwfc)%l
        m=nlmchi(nwfc)%m
        do nwfc1=1,nkb
           if (nlmchi(nwfc1)%na==na.AND.nlmchi(nwfc1)%l==l&
&             .AND.nlmchi(nwfc1)%m==m) THEN
              nbp=nlmchi(nwfc1)%n     
              proj(nwfc,:,ik)=proj(nwfc,:,ik) + &
&                proj0(nwfc,:)*CONJG(proj0(nwfc1,:))*pcharge(nb,nbp,nt)
           endif
        enddo
    enddo    
     
 ENDDO
 DEALLOCATE(proj0,pcharge)

 CALL deallocate_bec_type (becp)
 !
 RETURN
 !
END SUBROUTINE projwave_paw
!
!-----------------------------------------------------------------------
FUNCTION compute_mj(j,l,m)
   !-----------------------------------------------------------------------
   USE kinds, ONLY: DP
   IMPLICIT NONE
   !
   REAL(DP) :: compute_mj, j
   INTEGER  :: l, m

   IF (abs(j-l-0.5d0)<1.d-4) THEN
       compute_mj=m+0.5d0
   ELSEIF (abs(j-l+0.5d0)<1.d-4) THEN
      compute_mj=m-0.5d0
   ELSE
      CALL errore('compute_mj','l and j not compatible',1)
   ENDIF

   RETURN
END FUNCTION compute_mj
!
!-----------------------------------------------------------------------
SUBROUTINE  write_proj_iotk (filename, lbinary, projs, lwrite_ovp, ovps )
  !-----------------------------------------------------------------------
  !
  USE kinds
  USE io_files,         ONLY : iun => iunsat, prefix, tmp_dir, postfix
  USE basis,            ONLY : natomwfc
  USE klist,            ONLY : wk, xk, nkstot, nelec
  USE noncollin_module, ONLY : noncolin
  USE lsda_mod,         ONLY : nspin, isk
  USE ener,             ONLY : ef
  USE wvfct,            ONLY : et, nbnd
  USE iotk_module
  IMPLICIT NONE

  CHARACTER(*),  INTENT(IN) :: filename
  LOGICAL,       INTENT(IN) :: lbinary
  COMPLEX(DP),   INTENT(IN) :: projs(natomwfc,nbnd,nkstot)
  LOGICAL,       INTENT(IN) :: lwrite_ovp
  COMPLEX(DP),   INTENT(IN) :: ovps(natomwfc,natomwfc,nkstot)
  !
  CHARACTER(256)          :: tmp
  CHARACTER(iotk_attlenx) :: attr
  INTEGER :: ik, ik_eff, isp, ia, ierr, num_k_points

!
! subroutine body
!

  tmp = trim( tmp_dir ) // trim( prefix ) // postfix //trim(filename)
  !
  IF ( lbinary ) THEN
      tmp = TRIM(tmp) // ".dat"
  ELSE
      tmp = TRIM(tmp) // ".xml"
  ENDIF
  !
  CALL iotk_open_write(iun, FILE=trim(tmp), ROOT="ATOMIC_PROJECTIONS", &
                       BINARY=lbinary, IERR=ierr )
  IF ( ierr /= 0 ) RETURN
  !
  !
  num_k_points = nkstot
  IF ( nspin == 2 ) num_k_points = nkstot / 2
  !
  CALL iotk_write_begin(iun, "HEADER")
  !
  CALL iotk_write_dat(iun, "NUMBER_OF_BANDS", nbnd)
  CALL iotk_write_dat(iun, "NUMBER_OF_K-POINTS", num_k_points )
  CALL iotk_write_dat(iun, "NUMBER_OF_SPIN_COMPONENTS", nspin)
  CALL iotk_write_dat(iun, "NON-COLINEAR_CALCULATION",noncolin)
  CALL iotk_write_dat(iun, "NUMBER_OF_ATOMIC_WFC", natomwfc)
  CALL iotk_write_dat(iun, "NUMBER_OF_ELECTRONS", nelec )
  CALL iotk_write_attr(attr, "UNITS", " 2 pi / a", FIRST=.true.  )
  CALL iotk_write_empty (iun,  "UNITS_FOR_K-POINTS", ATTR=attr)
  CALL iotk_write_attr(attr, "UNITS", "Rydberg", FIRST=.true.  )
  CALL iotk_write_empty (iun,  "UNITS_FOR_ENERGY", ATTR=attr)
  CALL iotk_write_dat(iun, "FERMI_ENERGY", ef )
  !
  CALL iotk_write_end(iun, "HEADER")
  !
  !
  CALL iotk_write_dat(iun, "K-POINTS", xk(:,1:num_k_points) , COLUMNS=3 )
  CALL iotk_write_dat(iun, "WEIGHT_OF_K-POINTS", wk(1:num_k_points), COLUMNS=8 )
  !
  CALL iotk_write_begin(iun, "EIGENVALUES")
  !
  DO ik=1,num_k_points
     !
     CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
     IF ( nspin == 2 ) THEN
        !
        ik_eff = ik + num_k_points
        !
        CALL iotk_write_dat( iun, "EIG.1", et(:,ik) )
        CALL iotk_write_dat( iun, "EIG.2", et(:,ik_eff) )
        !
     ELSE
        !
        CALL iotk_write_dat( iun, "EIG", et(:,ik) )
        !
     ENDIF
     !
     CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
  ENDDO
  !
  CALL iotk_write_end(iun, "EIGENVALUES")

  !
  ! main loop atomic wfc
  !
  CALL iotk_write_begin(iun, "PROJECTIONS")
  !
  DO ik=1,num_k_points
     !
     CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
     IF ( nspin == 2 ) THEN
        !
        CALL iotk_write_begin ( iun, "SPIN.1" )
           !
           DO ia = 1, natomwfc
               CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik)  )
           ENDDO
           !
        CALL iotk_write_end ( iun, "SPIN.1" )
        !
        ik_eff = ik + num_k_points
        !
        CALL iotk_write_begin ( iun, "SPIN.2" )
           !
           DO ia = 1, natomwfc
               CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik_eff)  )
           ENDDO
           !
        CALL iotk_write_end ( iun, "SPIN.2" )
        !
     ELSE
        !
        DO ia = 1,natomwfc
            CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik)  )
        ENDDO
        !
     ENDIF
     !
     CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
  ENDDO
  !
  CALL iotk_write_end(iun, "PROJECTIONS")

  !
  ! overlaps
  !
  IF ( lwrite_ovp ) THEN
      !
      CALL iotk_write_begin(iun, "OVERLAPS")
      !
      DO ik=1,num_k_points
          !
          CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
          !
          DO isp = 1, nspin
              !
              ik_eff = ik + num_k_points * ( isp -1 )
              !
              CALL iotk_write_dat(iun, "OVERLAP"//trim(iotk_index(isp)), ovps(:,:,ik_eff)  )
              !
              !
          ENDDO
          !
          CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
          !
      ENDDO
      !
      CALL iotk_write_end(iun, "OVERLAPS")
      !
  ENDIF
  !
  ! closing the file
  !
  CALL iotk_close_write(iun)

END SUBROUTINE write_proj_iotk
!
!-----------------------------------------------------------------------
SUBROUTINE write_proj_file ( filproj, proj )
  !-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE lsda_mod,  ONLY : nspin
  USE spin_orb,  ONLY : lspinorb
  USE noncollin_module, ONLY : noncolin
  USE fft_base,  ONLY : dfftp
  USE klist,     ONLY : xk, nkstot
  USE run_info,  ONLY : title
  USE cell_base, ONLY : at, ibrav, celldm
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, atm
  USE wvfct,     ONLY : nbnd
  USE basis,     ONLY : natomwfc
  USE gvect,     ONLY : gcutm 
  USE gvecs,     ONLY : dual
  USE gvecw,     ONLY : ecutwfc
  USE projections, ONLY : nlmchi
  !
  IMPLICIT NONE
  CHARACTER (len=*), INTENT(in) :: filproj
  REAL(DP), INTENT(IN) :: proj(natomwfc,nbnd,nkstot)
  !
  CHARACTER(256) :: filename
  REAL (DP), EXTERNAL :: compute_mj
  INTEGER :: is, ik, nwfc, ibnd, nksinit, nkslast, iunproj=33
  !
  IF ( TRIM(filproj) == ' ' ) RETURN
  !
  DO is=1,nspin
     IF (nspin==2) THEN
        IF (is==1) filename=trim(filproj)//'.projwfc_up'
        IF (is==2) filename=trim(filproj)//'.projwfc_down'
        nksinit=(nkstot/2)*(is-1)+1
        nkslast=(nkstot/2)*is
     ELSE
        filename=trim(filproj)//'.projwfc_up'
        nksinit=1
        nkslast=nkstot
     ENDIF
     CALL write_io_header(filename, iunproj, title, dfftp%nr1x, dfftp%nr2x, &
          dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, &
          celldm, at, gcutm, dual, ecutwfc, nkstot/nspin, nbnd, natomwfc)
     DO nwfc = 1, natomwfc
        IF (lspinorb) THEN
           WRITE(iunproj,'(2i5,1x,a4,1x,a2,1x,2i5,f5.1,1x,f5.1)') &
                nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                nlmchi(nwfc)%els, nlmchi(nwfc)%n, nlmchi(nwfc)%l, &
                nlmchi(nwfc)%jj, &
                compute_mj(nlmchi(nwfc)%jj,nlmchi(nwfc)%l, nlmchi(nwfc)%m)
        ELSE IF (noncolin) THEN
           WRITE(iunproj,'(2i5,1x,a4,1x,a2,1x,3i5,1x,f4.1)') &
                nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                nlmchi(nwfc)%els, nlmchi(nwfc)%n, nlmchi(nwfc)%l, &
                nlmchi(nwfc)%m, &
                0.5d0-int(nlmchi(nwfc)%ind/(2*nlmchi(nwfc)%l+2))
        ELSE
           WRITE(iunproj,'(2i5,1x,a4,1x,a2,1x,3i5)') &
                nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                nlmchi(nwfc)%els, nlmchi(nwfc)%n, nlmchi(nwfc)%l,  &
                nlmchi(nwfc)%m
        END IF
        DO ik=nksinit,nkslast
           DO ibnd=1,nbnd
              WRITE(iunproj,'(2i8,f20.10)') ik,ibnd, abs(proj(nwfc,ibnd,ik))
           ENDDO
        ENDDO
     ENDDO
     CLOSE(iunproj)
  ENDDO
  !
END SUBROUTINE write_proj_file
!
!  projwave with distributed matrixes
!
!-----------------------------------------------------------------------
SUBROUTINE projwave( filproj, lsym, lwrite_ovp, lbinary )
  !-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE ions_base, ONLY : nat
  USE basis,     ONLY : natomwfc, swfcatom
  USE klist,     ONLY : xk, nks, nkstot, nelec, ngk, igk_k
  USE lsda_mod,  ONLY : nspin
  USE noncollin_module, ONLY: noncolin, npol
  USE spin_orb,  ONLY : lspinorb, domag, lforcet
  USE wvfct,     ONLY : npwx, nbnd, et
  USE uspp,      ONLY : nkb, vkb
  USE becmod,    ONLY : bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files,  ONLY : prefix, restart_dir, tmp_dir
  USE control_flags, ONLY : gamma_only
  USE pw_restart_new,ONLY : read_collected_wfc
  USE wavefunctions, ONLY : evc
  !
  USE projections, ONLY: nlmchi, fill_nlmchi, proj, proj_aux, ovps_aux
  !
  USE io_files,  ONLY: nd_nmbr
  USE mp,        ONLY: mp_bcast
  USE mp_pools,  ONLY: root_pool, intra_pool_comm
  USE mp_diag,   ONLY: ortho_comm, np_ortho, me_ortho, ortho_comm_id, &
                       leg_ortho, ortho_cntx, nproc_ortho
  USE parallel_toolkit, ONLY : zsqmred, zsqmher, zsqmdst, zsqmcll, dsqmsym
  USE zhpev_module,     ONLY : pzhpev_drv, zhpev_drv
  USE descriptors,      ONLY : la_descriptor, descla_init
  !
  IMPLICIT NONE
  !
  CHARACTER (len=*), INTENT(IN) :: filproj
  LOGICAL, INTENT(IN)    :: lsym
  LOGICAL, INTENT(IN)    :: lbinary
  LOGICAL, INTENT(INOUT) :: lwrite_ovp
  !
  INTEGER :: npw, npw_, ik, ibnd, i, j, k, na, nb, nt, isym, n,  m, l, nwfc,&
       lmax_wfc, is
  REAL(DP),    ALLOCATABLE :: e (:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:), proj0(:,:)
  COMPLEX(DP), ALLOCATABLE :: e_work_d(:,:)
  ! Some workspace for gamma-point calculation ...
  REAL   (DP), ALLOCATABLE :: rproj0(:,:)
  COMPLEX(DP), ALLOCATABLE :: overlap_d(:,:), work_d(:,:), diag(:,:), vv(:,:)
  REAL   (DP), ALLOCATABLE ::roverlap_d(:,:)
  !
  INTEGER  :: nksinit, nkslast
  LOGICAL :: freeswfcatom
  !
  INTEGER :: iunaux
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(len=256) :: auxname
  !
  TYPE(la_descriptor) :: desc
  TYPE(la_descriptor), ALLOCATABLE :: desc_ip( :, : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
  ! matrix distribution descriptors
  INTEGER :: nx, nrl, nrlx
  ! maximum local block dimension
  LOGICAL :: la_para
  ! flag for parallel linear algebra
  LOGICAL :: la_proc
  ! distinguishes active procs in parallel linear algebra
  !
  !
  IF ( natomwfc <= 0 ) CALL errore &
        ('projwave', 'Cannot project on zero atomic wavefunctions!', 1)
  WRITE( stdout, '(/5x,"Calling projwave .... ")')
  la_para = ( nproc_ortho > 1 )
  IF ( la_para ) WRITE( stdout, &
       '(5x,"linear algebra parallelized on ",i3," procs")') nproc_ortho
  IF ( gamma_only ) &
       WRITE( stdout, '(5x,"gamma-point specific algorithms are used")')
  !
  ! fill structure nlmchi
  !
  CALL fill_nlmchi ( natomwfc, nwfc, lmax_wfc )
  !
  ALLOCATE( proj (natomwfc, nbnd, nkstot) )
  !
  IF( la_para ) lwrite_ovp = .FALSE. ! not implemented
  !
  IF (.not. ALLOCATED(swfcatom)) THEN
     ALLOCATE(swfcatom (npwx*npol , natomwfc ) )
     freeswfcatom = .true.
  ELSE
     freeswfcatom = .false.
  ENDIF
  ALLOCATE(wfcatom (npwx*npol, natomwfc) )
  ALLOCATE(e (natomwfc) )
  !
  ! Open file as temporary storage
  !
  iunaux = find_free_unit()
  auxname = TRIM( restart_dir() ) // 'AUX' // TRIM(nd_nmbr)
  OPEN( unit=iunaux, file=trim(auxname), status='unknown', form='unformatted')
  !
  ALLOCATE( desc_ip( np_ortho(1), np_ortho(2) ) )
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ) )
  !
  CALL desc_init( natomwfc, desc, desc_ip )
  la_proc = ( desc%active_node > 0 )
  nx = desc%nrcx
  !
  IF( ionode ) THEN
     WRITE( stdout, * )
     WRITE( stdout, * ) ' Problem Sizes '
     WRITE( stdout, * ) ' natomwfc = ', natomwfc
     IF ( la_para ) WRITE( stdout, * ) ' nx       = ', nx
     WRITE( stdout, * ) ' nbnd     = ', nbnd
     WRITE( stdout, * ) ' nkstot   = ', nkstot
     WRITE( stdout, * ) ' npwx     = ', npwx
     WRITE( stdout, * ) ' nkb      = ', nkb
     WRITE( stdout, * )
  ENDIF
  !
  !    loop on k points
  !
  DO ik = 1, nks
     !
     npw = ngk(ik)
     CALL read_collected_wfc ( restart_dir() , ik, evc )
     !
     wfcatom(:,:) = (0.0_dp, 0.0_dp)
     IF (lforcet) THEN
        !    AlexS - To project on real harmonics, not on spinors.  
        CALL atomic_wfc_nc_updown(ik, wfcatom)
     ELSE IF ( noncolin ) THEN
        CALL atomic_wfc_nc_proj (ik, wfcatom)
     ELSE
        CALL atomic_wfc (ik, wfcatom)
     ENDIF
     !
     CALL allocate_bec_type (nkb, natomwfc, becp )
     !
     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
     CALL calbec ( npw, vkb, wfcatom, becp)
     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)
     !
     CALL deallocate_bec_type (becp)
     !
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i>
     ! calculate overlap matrix O_ij = <phi_i|\hat S|\phi_j>
     !
     IF( la_proc ) THEN
        ALLOCATE(overlap_d (nx, nx) )
     ELSE
        ALLOCATE(overlap_d (1, 1) )
     ENDIF
     overlap_d = (0.d0,0.d0)
     npw_=npw
     IF ( noncolin ) npw_=npol*npwx
     IF ( gamma_only ) THEN
        !
        ! in the Gamma-only case the overlap matrix (real) is copied 
        ! to a complex one as for the general case - easy but wasteful
        !
        IF( la_proc ) THEN
           ALLOCATE(roverlap_d (nx, nx) )
        ELSE
           ALLOCATE(roverlap_d (1, 1) )
        ENDIF
        roverlap_d = 0.d0
        CALL calbec_ddistmat( npw, wfcatom, swfcatom, natomwfc, nx, roverlap_d )
        overlap_d(:,:)=cmplx(roverlap_d(:,:),0.0_dp, kind=dp)
     ELSE 
        CALL calbec_zdistmat( npw_, wfcatom, swfcatom, natomwfc, nx, overlap_d )
     ENDIF
     !
     ! save overlap matrix if required
     !
     IF ( lwrite_ovp ) WRITE( iunaux ) overlap_d
     !
     ! diagonalize the overlap matrix
     !
     IF ( la_proc ) THEN
        !
        ALLOCATE(work_d (nx, nx) )

        nrl  = desc%nrl
        nrlx = desc%nrlx

        ALLOCATE( diag( nrlx, natomwfc ) )
        ALLOCATE( vv( nrlx, natomwfc ) )
        !
        !  re-distribute the overlap matrix for parallel diagonalization
        !
        CALL blk2cyc_zredist( natomwfc, diag, nrlx, natomwfc, overlap_d, nx, nx, desc )
        !
        ! parallel diagonalization
        !
        CALL pzhpev_drv( 'V', diag, nrlx, e, vv, nrlx, nrl, natomwfc, &
           desc%npc * desc%npr, desc%mype, desc%comm )
        !
        !  bring distributed eigenvectors back to original distribution
        !
        CALL cyc2blk_zredist( natomwfc, vv, nrlx, natomwfc, work_d, nx, nx, desc )
        !
        DEALLOCATE( vv )
        DEALLOCATE( diag )
        !
     ELSE
        ALLOCATE(work_d (1, 1) )
     ENDIF

     CALL mp_bcast( e, root_pool, intra_pool_comm )

     ! calculate O^{-1/2} (actually, its transpose)

     DO i = 1, natomwfc
        e (i) = 1.d0 / dsqrt (e (i) )
     ENDDO

     IF ( la_proc ) THEN
        ALLOCATE(e_work_d (nx, nx) )
        DO j = 1, desc%nc
           DO i = 1, desc%nr
              e_work_d( i, j ) = e( j + desc%ic - 1 ) * work_d( i, j )
           ENDDO
        ENDDO
        CALL sqr_zmm_cannon( 'N', 'C', natomwfc, (1.0_dp,0.0_dp), e_work_d, &
            nx, work_d, nx, (0.0_dp, 0.0_dp), overlap_d, nx, desc )
        CALL zsqmher( natomwfc, overlap_d, nx, desc )
        DEALLOCATE( e_work_d )
     ENDIF
     !
     DEALLOCATE( work_d )
     !
     ! calculate wfcatom = O^{-1/2} \hat S | phi>
     !
     IF ( gamma_only ) THEN
        roverlap_d(:,:)=REAL(overlap_d(:,:),DP)
        CALL wf_times_roverlap( nx, npw, swfcatom, roverlap_d, wfcatom )
        DEALLOCATE( roverlap_d )
     ELSE
        CALL wf_times_overlap( nx, npw_, swfcatom, overlap_d, wfcatom )
     ENDIF
     DEALLOCATE( overlap_d )
     !
     ! make the projection <psi_i| O^{-1/2} \hat S | phi_j>,
     ! symmetrize the projections if required
     !
     IF ( gamma_only ) THEN
        !
        ALLOCATE( rproj0(natomwfc,nbnd) )
        CALL calbec ( npw, wfcatom, evc, rproj0)
        WRITE( iunaux ) rproj0
        IF (lsym) THEN
           CALL sym_proj_g (rproj0, proj(:,:,ik))
        ELSE
           proj(:,:,ik)=abs(rproj0(:,:))**2
        ENDIF
        DEALLOCATE (rproj0)
        !
     ELSE
        !
        ALLOCATE( proj0(natomwfc,nbnd) )
        CALL calbec ( npw_, wfcatom, evc, proj0)
        WRITE( iunaux ) proj0
        IF (lsym) THEN
           IF ( lspinorb ) THEN 
              CALL sym_proj_so ( domag, proj0, proj(:,:,ik) )
           ELSE IF (noncolin) THEN
              CALL sym_proj_nc ( proj0, proj(:,:,ik) )
           ELSE
              CALL sym_proj_k (proj0, proj(:,:,ik))
           END IF
        ELSE
           proj(:,:,ik)=abs(proj0(:,:))**2
        END IF
        DEALLOCATE (proj0)
        !
     ENDIF
     ! on k-points
  ENDDO
  !
  DEALLOCATE (e)
  DEALLOCATE (wfcatom)
  IF (freeswfcatom) DEALLOCATE (swfcatom)
  !
  ! closing the file will cause a lot of I/O
  !!! CLOSE( unit=iunaux )
  !
  !   vectors et and proj are distributed across the pools
  !   collect data for all k-points to the first pool
  !   (I think it is not actually needed for et)
  !
  CALL poolrecover (et,       nbnd, nkstot, nks)
  CALL poolrecover (proj,     nbnd * natomwfc, nkstot, nks)
  !
  ! write to standard output and to file filproj (if required)
  !
  IF ( ionode) THEN
     !
     CALL write_proj( lmax_wfc, filproj, proj )
     CALL write_proj_file ( filproj, proj )
     !
  END IF
  !
  !  Recover proj_aux and (if required) overlap matrices for all k-points
  !
  !  See above: no reason to close and re-open the file
  !!! OPEN( unit=iunaux, file=trim(auxname), status='old', form='unformatted')
  REWIND (unit=iunaux)
  !
  IF ( lwrite_ovp ) THEN
     ALLOCATE( ovps_aux(natomwfc, natomwfc, nkstot) )
  ELSE
     ALLOCATE( ovps_aux(1, 1, 1) )
  ENDIF
  ALLOCATE( proj_aux (natomwfc, nbnd, nkstot) )
  proj_aux = (0.d0, 0.d0)
  !
  DO ik = 1, nks
     !
     IF ( lwrite_ovp ) READ( iunaux)  ovps_aux(:,:,ik)
     IF( gamma_only ) THEN
        ALLOCATE( rproj0( natomwfc, nbnd ) )
        READ( iunaux ) rproj0(:,:)
        proj_aux(:,:,ik) = cmplx( rproj0(:,:), 0.00_dp, kind=dp )
        DEALLOCATE ( rproj0 )
     ELSE
        READ( iunaux ) proj_aux(:,:,ik)
     ENDIF
     !
  ENDDO
  !
  CLOSE( unit=iunaux, status='delete' )
  !
  CALL poolrecover (proj_aux, 2 * nbnd * natomwfc, nkstot, nks)
  IF ( lwrite_ovp ) &
      CALL poolrecover (ovps_aux, 2 * natomwfc * natomwfc, nkstot, nks)
  !
  IF ( ionode ) THEN
     !
     ! write projections to file using iotk
     !
     CALL write_proj_iotk( "atomic_proj", lbinary, proj_aux, lwrite_ovp, &
          ovps_aux )
     !
  ENDIF
  !
  DEALLOCATE( proj_aux, ovps_aux )
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE desc_init( nsiz, desc, desc_ip )
     !
     INTEGER, INTENT(in)  :: nsiz
     TYPE(la_descriptor), INTENT(out) :: desc
     TYPE(la_descriptor), INTENT(out) :: desc_ip(:,:)
     INTEGER :: i, j, rank
     INTEGER :: coor_ip( 2 )
     !
     CALL descla_init( desc, nsiz, nsiz, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )
     !
     DO j = 0, desc%npc - 1
        DO i = 0, desc%npr - 1
           coor_ip( 1 ) = i
           coor_ip( 2 ) = j
           CALL descla_init( desc_ip(i+1,j+1), desc%n, desc%nx, np_ortho, coor_ip, ortho_comm, ortho_cntx, 1 )
           CALL GRID2D_RANK( 'R', desc%npr, desc%npc, i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        ENDDO
     ENDDO
     !
     RETURN
  END SUBROUTINE desc_init
  !

  SUBROUTINE calbec_zdistmat( npw, v, w, n, nx, dm )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm
     !
     USE mp, ONLY : mp_root_sum
     !
     IMPLICIT NONE
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root, ldv, ldw
     INTEGER, INTENT(in) :: npw ! local number of plane wave
     INTEGER, INTENT(in) :: n   ! global dimension of matrix dm
     INTEGER, INTENT(in) :: nx  ! local leading dimension of matrix dm
                                ! WARNING: nx is the same on all proc, SIZE( dm, 1 ) NO!
     COMPLEX(DP), INTENT(out) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = (0.0_dp, 0.0_dp)
     !
     ldv = size( v, 1 )
     ldw = size( w, 1 )
     !
     DO ipc = 1, desc%npc !  loop on column procs
        !
        nc = desc_ip( 1, ipc )%nc
        ic = desc_ip( 1, ipc )%ic
        !
        DO ipr = 1, ipc ! desc( la_npr_ ) ! ipc ! use symmetry for the loop on row procs
           !
           nr = desc_ip( ipr, ipc )%nr
           ir = desc_ip( ipr, ipc )%ir
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL ZGEMM( 'C', 'N', nr, nc, npw, (1.0_dp,0.0_dp) , &
                       v(1,ir), ldv, w(1,ic), ldw, (0.0_dp,0.0_dp), work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, intra_pool_comm )

        ENDDO
        !
     ENDDO
     !
     CALL zsqmher( n, dm, nx, desc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE calbec_zdistmat
  !

  SUBROUTINE calbec_ddistmat( npw, v, w, n, nx, dm )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm
     !
     USE mp, ONLY : mp_root_sum
     USE gvect, ONLY : gstart
     !
     IMPLICIT NONE
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root, ldv, ldw, npw2, npwx2
     INTEGER, INTENT(in) :: npw ! local number of plane wave
     INTEGER, INTENT(in) :: n   ! global dimension of matrix dm
     INTEGER, INTENT(in) :: nx  ! local leading dimension of matrix dm
                                ! WARNING: nx is the same on all proc, SIZE( dm, 1 ) NO!
     REAL(DP), INTENT(out) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     REAL(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     npw2  = 2*npw
     npwx2 = 2*npwx
     !
     work = (0.0_dp, 0.0_dp)
     !
     ldv = size( v, 1 )
     ldw = size( w, 1 )
     !
     DO ipc = 1, desc%npc !  loop on column procs
        !
        nc = desc_ip( 1, ipc )%nc
        ic = desc_ip( 1, ipc )%ic
        !
        DO ipr = 1, ipc ! desc( la_npr_ ) ! ipc ! use symmetry for the loop on row procs
           !
           nr = desc_ip( ipr, ipc )%nr
           ir = desc_ip( ipr, ipc )%ir
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           ! use blas subs. on the matrix block

           CALL DGEMM( 'T', 'N', nr, nc, npw2, 2.D0 , &
                       v(1,ir), npwx2, w(1,ic), npwx2, 0.D0, work, nx )

           IF ( gstart == 2 ) &
              CALL DGER( nr, nc, -1.D0, v(1,ir), npwx2, w(1,ic), npwx2, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, intra_pool_comm )

        ENDDO
        !
     ENDDO
     !
     CALL dsqmsym( n, dm, nx, desc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE calbec_ddistmat
  !
  !
  !
  SUBROUTINE wf_times_overlap( nx, npw, swfc, ovr, wfc )
    !
     INTEGER, INTENT(in) :: nx, npw
     COMPLEX(DP) :: swfc( :, : ), ovr( :, : ), wfc( :, : )
     !
     INTEGER :: npwx
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     npwx = SIZE(swfc,1)
     DO ipc = 1, desc%npc
        !
        nc = desc_ip( 1, ipc )%nc
        ic = desc_ip( 1, ipc )%ic
        !
        beta = (0.0_dp, 0.0_dp)

        DO ipr = 1, desc%npr
           !
           nr = desc_ip( ipr, ipc )%nr
           ir = desc_ip( ipr, ipc )%ir
           !
           root = rank_ip( ipr, ipc )

           IF( ipr-1 == desc%myr .and. ipc-1 == desc%myc .and. la_proc ) THEN
              !
              !  this proc sends his block
              !
              CALL mp_bcast( ovr, root, intra_pool_comm )
              CALL ZGEMM( 'N', 'N', npw, nc, nr, (1.0_dp,0.0_dp), &
                          swfc(1,ir), npwx, ovr, nx, beta, wfc(1,ic), npwx )
           ELSE
              !
              !  all other procs receive
              !
              CALL mp_bcast( vtmp, root, intra_pool_comm )
              CALL ZGEMM( 'N', 'N', npw, nc, nr, (1.0_dp,0.0_dp), &
                       swfc(1,ir), npwx, vtmp, nx, beta, wfc(1,ic), npwx )
           ENDIF
           !
           beta = (1.0_dp,0.0_dp)

        ENDDO
        !
     ENDDO
     !
     DEALLOCATE( vtmp )

     RETURN

  END SUBROUTINE wf_times_overlap

  !
  SUBROUTINE wf_times_roverlap( nx, npw, swfc, ovr, wfc )

     USE gvect, ONLY : gstart

     INTEGER, INTENT(in) :: nx, npw
     COMPLEX(DP) :: swfc( :, : ), wfc( :, : )
     REAL(DP)    :: ovr( :, : )
     !
     INTEGER :: ipc, ipr, npw2, npwx2
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     REAL(DP) :: beta

     npw2  = 2*npw
     npwx2 = 2*SIZE(swfc,1)

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc%npc
        !
        nc = desc_ip( 1, ipc )%nc
        ic = desc_ip( 1, ipc )%ic
        !
        beta = 0.0d0

        DO ipr = 1, desc%npr
           !
           nr = desc_ip( ipr, ipc )%nr
           ir = desc_ip( ipr, ipc )%ir
           !
           root = rank_ip( ipr, ipc )

           IF( ipr-1 == desc%myr .and. ipc-1 == desc%myc .and. la_proc ) THEN
              !
              !  this proc sends his block
              !
              CALL mp_bcast( ovr, root, intra_pool_comm )
              CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          swfc(1,ir), npwx2, ovr, nx, beta, wfc(1,ic), npwx2 )
              !
           ELSE
              !
              !  all other procs receive
              !
              CALL mp_bcast( vtmp, root, intra_pool_comm )
              CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          swfc(1,ir), npwx2, vtmp, nx, beta, wfc(1,ic), npwx2 )
              !
           ENDIF
           !
           beta = 1.0d0

        ENDDO
        !
     ENDDO
     !
     DEALLOCATE( vtmp )

     RETURN

  END SUBROUTINE wf_times_roverlap
  !
END SUBROUTINE projwave
!
