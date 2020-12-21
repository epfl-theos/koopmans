!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!----------------------------------------------------------------------- 
PROGRAM projwfc
  !----------------------------------------------------------------------- 
  ! 
  ! projects wavefunctions onto orthogonalized atomic wavefunctions, 
  ! calculates Lowdin charges, spilling parameter, projected DOS 
  ! 
  !
  ! Input (namelist &inputpp ... / ):                         Default value
  !
  !    prefix        prefix of input file produced by pw.x    'pwscf'
  !                    (wavefunctions are needed)
  !    outdir        directory containing the input file       ./
  !    ngauss        type of gaussian broadening (optional)    0
  !            =  0  Simple Gaussian (default)
  !            =  1  Methfessel-Paxton of order 1
  !            = -1  Marzari-Vanderbilt "cold smearing"
  !            =-99  Fermi-Dirac function
  !    degauss       gaussian broadening, Ry (not eV!)          0.0
  !    Emin, Emax    min, max energy (eV) for DOS plot          band extrema
  !    DeltaE        energy grid step (eV)                      none
  !    lsym          if true the projections are symmetrized    .true.
  !    filproj       file containing the projections            none
  !    filpdos       prefix for output files containing PDOS(E) prefix
  !
  ! Output:
  !
  !   Projections are written to standard output,
  !   and also to file filproj if given as input.
  !   The total DOS and the sum of projected DOS are written to file 
  !   "filpdos".pdos_tot.
  !   The format for the collinear, spin-unpolarized case and the
  !   non-collinear, spin-orbit case is
  !        E DOS(E) PDOS(E)
  !   The format for the collinear, spin-polarized case is
  !        E DOSup(E) DOSdw(E)  PDOSup(E) PDOSdw(E) 
  !   The format for the non-collinear, non spin-orbit case is
  !        E DOS(E) PDOSup(E) PDOSdw(E) 
  !   
  !   In the collinear case and the non-collinear, non spin-orbit case
  !   projected DOS are written to file "filpdos".pdos_atm#N(X)_wfc#M(l),
  !   where N = atom number , X = atom symbol, M = wfc number, l=s,p,d,f
  !   (one file per atomic wavefunction found in the pseudopotential file) 
  !   - The format for the collinear, spin-unpolarized case is
  !        E LDOS(E) PDOS_1(E) ... PDOS_2l+1(E)
  !     where LDOS = \sum m=1,2l+1 PDOS_m(E)
  !     and PDOS_m(E) = projected DOS on atomic wfc with component m
  !   - The format for the collinear, spin-polarized case and the 
  !     non-collinear, non spin-orbit case is as above with
  !     two components for both  LDOS(E) and PDOS_m(E)
  !    
  !   In the non-collinear, spin-orbit case (i.e. if there is at least one
  !   fully relativistic pseudopotential) wavefunctions are projected
  !   onto eigenstates of the total angular-momentum. 
  !   Projected DOS are written to file "filpdos".pdos_atm#N(X)_wfc#M(l_j),
  !   where N = atom number , X = atom symbol, M = wfc number, l=s,p,d,f
  !   and j is the value of the total angular momentum.
  !   In this case the format is
  !      E LDOS(E) PDOS_1(E) ... PDOS_2j+1(E)
  !
  !   All DOS(E) are in states/eV plotted vs E in eV
  !
  !  Order of m-components for each l in the output:
  !
  !  1, cos(phi), sin(phi), cos(2*phi), sin(2*phi), .., cos(l*phi), sin(l*phi)
  !
  !  where phi is the polar angle:x=r cos(theta)cos(phi), y=r cos(theta)sin(phi)
  !  This is determined in file flib/ylmr2.f90 that calculates spherical harm.
  !      L=1 :
  !  1 pz     (m=0)
  !  2 px     (real combination of m=+/-1 with cosine)
  !  3 py     (real combination of m=+/-1 with sine)
  !      L=2 :
  !  1 dz2    (m=0)
  !  2 dzx    (real combination of m=+/-1 with cosine)
  !  3 dzy    (real combination of m=+/-1 with sine)
  !  4 dx2-y2 (real combination of m=+/-2 with cosine)
  !  5 dxy    (real combination of m=+/-1 with sine)
  ! 
  ! Important notice:
  !
  !    The tetrahedron method is presently not implemented.
  !    Gaussian broadening is used in all cases:
  !    - if degauss is set to some value in namelist &inputpp, that value
  !      (and the optional value for ngauss) is used
  !    - if degauss is NOT set to any value in namelist &inputpp, the 
  !      value of degauss and of ngauss are read from the input data
  !      file (they will be the same used in the pw.x calculations)
  !    - if degauss is NOT set to any value in namelist &inputpp, AND
  !      there is no value of degauss and of ngauss in the input data
  !      file, degauss=DeltaE (in Ry) and ngauss=0 will be used
  ! Obsolete variables, ignored:
  !   io_choice
  !   smoothing
  ! 
#include "f_defs.h" 
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE constants,  ONLY : rytoev 
  USE kinds,      ONLY : DP 
  USE klist,      ONLY : degauss, ngauss, lgauss
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir, trimcheck
  USE noncollin_module, ONLY : noncolin
  USE mp,               ONLY : mp_bcast      
  USE control_flags,    ONLY : ortho_para
  !
  IMPLICIT NONE 
  CHARACTER (len=256) :: filpdos, filproj, io_choice, outdir
  REAL (DP)      :: Emin, Emax, DeltaE, degauss1, smoothing
  INTEGER :: ngauss1, ios
  LOGICAL :: lsym
  !
  NAMELIST / inputpp / outdir, prefix, ngauss, degauss, lsym, &
             Emin, Emax, DeltaE, io_choice, smoothing, filpdos, filproj
  ! 
  CALL start_postproc (nd_nmbr) 
  ! 
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf' 
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  filproj= ' '
  filpdos= ' '
  Emin   =-1000000.d0 
  Emax   =+1000000.d0 
  DeltaE = 0.01d0 
  ngauss = 0
  lsym   = .true.
  degauss= 0.d0
  !
  ios = 0
  !
  IF ( ionode )  THEN  
     !
     CALL input_from_file ( )
     !
     READ (5, inputpp, iostat = ios) 
     !
     tmp_dir = trimcheck (outdir) 
     ! save the value of degauss and ngauss: they are read from file
     degauss1=degauss
     ngauss1 = ngauss
     ! 
  END IF
  !
  CALL mp_bcast (ios, ionode_id )
  !
  IF (ios /= 0) CALL errore ('projwfc', 'reading inputpp namelist', ABS (ios) ) 
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( tmp_dir, ionode_id ) 
  CALL mp_bcast( prefix,  ionode_id ) 
  CALL mp_bcast( filproj,  ionode_id ) 
  CALL mp_bcast( ngauss1, ionode_id ) 
  CALL mp_bcast( degauss1,ionode_id ) 
  CALL mp_bcast( DeltaE,  ionode_id ) 
  CALL mp_bcast( lsym,  ionode_id ) 
  CALL mp_bcast( Emin, ionode_id ) 
  CALL mp_bcast( Emax, ionode_id ) 
  ! 
  !   Now allocate space for pwscf variables, read and check them. 
  ! 
  CALL read_file ( )
  !
  CALL openfil_pp ( )
  ! 
  !   decide Gaussian broadening
  !
  IF (degauss1.NE.0.d0) THEN
     degauss=degauss1
     ngauss =ngauss1
     WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.TRUE.
  ELSE IF (lgauss) THEN
     WRITE( stdout,'(/5x,"Gaussian broadening (read from file): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
  ELSE
     degauss=DeltaE/rytoev
     ngauss =0
     WRITE( stdout,'(/5x,"Gaussian broadening (default values): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.TRUE.
  END IF
  !
  IF (noncolin) THEN 
     CALL projwave_nc(filproj, lsym )
  ELSE
     IF( ortho_para > 1 ) THEN
        CALL pprojwave (filproj, lsym)
     ELSE
        CALL projwave (filproj, lsym)
     END IF
  ENDIF
  !
  IF ( ionode .and. lsym ) THEN 
     !
     IF ( filpdos == ' ') filpdos = prefix
     !
     IF (noncolin) THEN 
        CALL partialdos_nc (Emin, Emax, DeltaE, filpdos)
     ELSE
        CALL partialdos (Emin, Emax, DeltaE, filpdos)
     ENDIF
     !
  END IF 
  !
  CALL stop_pp 
  ! 
END PROGRAM projwfc
!
MODULE projections
  USE kinds, ONLY : DP 

  TYPE wfc_label 
     INTEGER na, n, l, m 
  END TYPE wfc_label 
  TYPE(wfc_label), ALLOCATABLE :: nlmchi(:) 

  REAL (DP),    ALLOCATABLE :: proj (:,:,:)
  COMPLEX (DP), ALLOCATABLE :: proj_aux (:,:,:)

END MODULE projections
!
MODULE projections_nc
  USE kinds, ONLY : DP

  TYPE wfc_label_nc 
     INTEGER na, n, l, m, ind 
     REAL (DP) jj
  END TYPE wfc_label_nc 
  TYPE(wfc_label_nc), ALLOCATABLE :: nlmchi(:) 

  REAL (DP),    ALLOCATABLE :: proj (:,:,:)
  COMPLEX (DP), ALLOCATABLE :: proj_aux (:,:,:)

END MODULE projections_nc
!
!----------------------------------------------------------------------- 
SUBROUTINE projwave( filproj, lsym )
  !----------------------------------------------------------------------- 
  ! 
  USE io_global, ONLY : stdout, ionode
  USE printout_base, ONLY: title
  USE ions_base, ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE basis,     ONLY : natomwfc
  USE cell_base 
  USE constants, ONLY: rytoev, eps4 
  USE gvect 
  USE klist, ONLY: xk, nks, nkstot, nelec
  USE ldaU 
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE symme, ONLY: nsym, irt, d1, d2, d3
  USE wvfct 
  USE control_flags, ONLY: gamma_only 
  USE uspp, ONLY: nkb, vkb
  USE uspp_param, ONLY: upf
  USE becmod,   ONLY: becp, rbecp, calbec
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc 
  USE spin_orb, ONLY: lspinorb
  USE wavefunctions_module, ONLY: evc 
  !
  USE projections
  !
  IMPLICIT NONE 
  !
  CHARACTER (len=*) :: filproj
  INTEGER :: ik, ibnd, i, j, k, na, nb, nt, isym, n,  m, m1, l, nwfc,& 
       nwfc1, lmax_wfc, is, ios, iunproj
  REAL(DP), ALLOCATABLE :: e (:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:)
  COMPLEX(DP), ALLOCATABLE :: overlap(:,:), work(:,:),work1(:), proj0(:,:)
  ! Some workspace for k-point calculation ... 
  REAL   (DP), ALLOCATABLE ::roverlap(:,:), rwork1(:),rproj0(:,:)
  ! ... or for gamma-point. 
  REAL(DP), ALLOCATABLE :: charges(:,:,:), proj1 (:)
  REAL(DP) :: psum, totcharge(2)
  INTEGER  :: nksinit, nkslast
  CHARACTER(LEN=256) :: filename
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/) 
  INTEGER, ALLOCATABLE :: idx(:) 
  LOGICAL :: lsym
  ! 
  ! 
  ! 
  WRITE( stdout, '(/5x,"Calling projwave .... ")') 
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used")') 
  END IF
  ! 
  ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1 
  ! 
  CALL d_matrix (d1, d2, d3)   
  ! 
  ! fill structure nlmchi 
  ! 
  ALLOCATE (nlmchi(natomwfc)) 
  nwfc=0 
  lmax_wfc = 0 
  DO na = 1, nat 
     nt = ityp (na) 
     DO n = 1, upf(nt)%nwfc
        IF (upf(nt)%oc (n) >= 0.d0) THEN 
           l = upf(nt)%lchi (n) 
           lmax_wfc = MAX (lmax_wfc, l ) 
           DO m = 1, 2 * l + 1 
              nwfc=nwfc+1 
              nlmchi(nwfc)%na = na 
              nlmchi(nwfc)%n  =  n 
              nlmchi(nwfc)%l  =  l 
              nlmchi(nwfc)%m  =  m 
           ENDDO 
        ENDIF 
     ENDDO 
  ENDDO 
  ! 
  IF (lmax_wfc > 3) CALL errore ('projwave', 'l > 3 not yet implemented', 1) 
  IF (nwfc /= natomwfc) CALL errore ('projwave', 'wrong # of atomic wfcs?', 1) 
  ! 
  ALLOCATE( proj (natomwfc, nbnd, nkstot) ) 
  ALLOCATE( proj_aux (natomwfc, nbnd, nkstot) ) 
  proj      = 0.d0 
  proj_aux  = (0.d0, 0.d0) 
  !
  IF (.NOT. lda_plus_u) ALLOCATE(swfcatom (npwx , natomwfc ) ) 
  ALLOCATE(wfcatom (npwx, natomwfc) ) 
  ALLOCATE(overlap (natomwfc, natomwfc) ) 
  overlap= (0.d0,0.d0) 
  !
  IF ( gamma_only ) THEN 
     ALLOCATE(roverlap (natomwfc, natomwfc) ) 
     roverlap= 0.d0 
     ALLOCATE (rbecp (nkb,natomwfc)) 
  ELSE
     ALLOCATE ( becp (nkb,natomwfc)) 
  END IF 
  ALLOCATE(e (natomwfc) ) 
  ! 
  !    loop on k points 
  ! 
  CALL init_us_1 
  CALL init_at_1 
  ! 
  DO ik = 1, nks 
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin) 
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1) 
 
     CALL atomic_wfc (ik, wfcatom) 
 
     CALL init_us_2 (npw, igk, xk (1, ik), vkb) 
 
     IF ( gamma_only ) THEN 
        CALL calbec ( npw, vkb, wfcatom, rbecp) 
     ELSE 
        CALL calbec ( npw, vkb, wfcatom, becp) 
     END IF 
 
     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom) 
     ! 
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i> 
     ! calculate overlap matrix O_ij = <phi_i|\hat S|\phi_j> 
     ! 
     IF ( gamma_only ) THEN 
        CALL calbec ( npw, wfcatom, swfcatom, roverlap ) 
        overlap(:,:)=CMPLX(roverlap(:,:),0.d0)
        ! TEMP: diagonalization routine for real matrix should be used instead 
     ELSE 
        CALL calbec ( npw, wfcatom, swfcatom, overlap ) 
     END IF 
 
     ! 
     ! calculate O^{-1/2} 
     ! 
     ALLOCATE(work (natomwfc, natomwfc) ) 
     CALL cdiagh (natomwfc, overlap, natomwfc, e, work) 
     DO i = 1, natomwfc 
        e (i) = 1.d0 / dsqrt (e (i) ) 
     ENDDO 
     DO i = 1, natomwfc 
        DO j = i, natomwfc 
           overlap (i, j) = (0.d0, 0.d0) 
           DO k = 1, natomwfc 
              overlap (i, j) = overlap (i, j) + e (k) * work (j, k) * CONJG (work (i, k) ) 
           ENDDO 
           IF (j /= i) overlap (j, i) = CONJG (overlap (i, j)) 
        ENDDO 
     ENDDO 
     DEALLOCATE (work) 
     ! 
     ! calculate wfcatom = O^{-1/2} \hat S | phi> 
     ! 
     IF ( gamma_only ) THEN 
        roverlap(:,:)=REAL(overlap(:,:),DP)
        ! TEMP: diagonalization routine for real matrix should be used instead 
        CALL DGEMM ('n', 't', 2*npw, natomwfc, natomwfc, 1.d0 , & 
             swfcatom, 2*npwx,  roverlap, natomwfc, 0.d0, wfcatom, 2*npwx) 
     ELSE 
        CALL ZGEMM ('n', 't', npw, natomwfc, natomwfc, (1.d0, 0.d0) , & 
             swfcatom, npwx,  overlap, natomwfc, (0.d0, 0.d0), wfcatom, npwx) 
     END IF 
 
     ! 
     ! make the projection <psi_i| O^{-1/2} \hat S | phi_j> 
     ! 
     IF ( gamma_only ) THEN 
        !
        ALLOCATE(rproj0(natomwfc,nbnd), rwork1 (nbnd) ) 
        CALL calbec ( npw, wfcatom, evc, rproj0) 
        !
        proj_aux(:,:,ik) = CMPLX( rproj0(:,:), 0.d0 )
        !
     ELSE 
        !
        ALLOCATE(proj0(natomwfc,nbnd), work1 (nbnd) ) 
        CALL calbec ( npw, wfcatom, evc, proj0) 
        !
        proj_aux(:,:,ik) = proj0(:,:)
        !
     END IF 
     ! 
     ! symmetrize the projections 
     ! 
     IF (lsym) THEN
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
                 IF (nlmchi(nwfc1)%na == nb             .AND. & 
                      nlmchi(nwfc1)%n == nlmchi(nwfc)%n .AND. & 
                      nlmchi(nwfc1)%l == nlmchi(nwfc)%l .AND. & 
                      nlmchi(nwfc1)%m == 1 ) go to 10 
              END DO 
              CALL errore('projwave','cannot symmetrize',1) 
10            nwfc1=nwfc1-1 
              !  
              !  nwfc1 is the first rotated atomic wfc corresponding to nwfc 
              ! 
              IF ( gamma_only ) THEN 
                 IF (l == 0) THEN 
                    rwork1(:) = rproj0 (nwfc1 + 1,:) 
                 ELSE IF (l == 1) THEN  
                    rwork1(:) = 0.d0   
                    DO m1 = 1, 3   
                       rwork1(:)=rwork1(:)+d1(m1,m,isym)*rproj0(nwfc1+m1,:) 
                    ENDDO 
                 ELSE IF (l == 2) THEN  
                    rwork1(:) = 0.d0   
                    DO m1 = 1, 5   
                       rwork1(:)=rwork1(:)+d2(m1,m,isym)*rproj0(nwfc1+m1,:) 
                    ENDDO 
                 ELSE IF (l == 3) THEN  
                    rwork1(:) = 0.d0   
                    DO m1 = 1, 7   
                       rwork1(:)=rwork1(:)+d3(m1,m,isym)*rproj0(nwfc1+m1,:) 
                    ENDDO 
                 ENDIF 
                 DO ibnd = 1, nbnd 
                    proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + & 
                         rwork1(ibnd) * rwork1(ibnd) / nsym 
                 ENDDO 
              ELSE 
                 IF (l == 0) THEN 
                    work1(:) = proj0 (nwfc1 + 1,:) 
                 ELSE IF (l == 1) THEN  
                    work1(:) = 0.d0   
                    DO m1 = 1, 3   
                       work1(:)=work1(:)+d1(m1,m,isym)*proj0(nwfc1+m1,:) 
                    ENDDO 
                 ELSE IF (l == 2) THEN  
                    work1(:) = 0.d0   
                    DO m1 = 1, 5   
                       work1(:)=work1(:)+d2(m1,m,isym)*proj0(nwfc1+m1,:) 
                    ENDDO 
                 ELSE IF (l == 3) THEN  
                    work1(:) = 0.d0   
                    DO m1 = 1, 7   
                       work1(:)=work1(:)+d3(m1,m,isym)*proj0(nwfc1+m1,:) 
                    ENDDO 
                 ENDIF 
                 DO ibnd = 1, nbnd 
                    proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + & 
                         work1(ibnd) * CONJG (work1(ibnd)) / nsym 
                 ENDDO 
              END IF 
           END DO 
        END DO 
     ELSE
        IF ( gamma_only ) THEN 
           DO nwfc=1,natomwfc
              DO ibnd=1,nbnd
                 proj(nwfc,ibnd,ik)=ABS(rproj0(nwfc,ibnd))**2
              END DO
           END DO
        ELSE
           DO nwfc=1,natomwfc
              DO ibnd=1,nbnd
                 proj(nwfc,ibnd,ik)=ABS(proj0(nwfc,ibnd))**2
              END DO
           END DO
        END IF
     END IF
     IF ( gamma_only ) THEN 
        DEALLOCATE (rwork1) 
        DEALLOCATE (rproj0) 
     ELSE 
        DEALLOCATE (work1) 
        DEALLOCATE (proj0) 
     END IF 
     ! on k-points 
  ENDDO 
  !
  DEALLOCATE (e)
  IF ( gamma_only ) THEN 
     DEALLOCATE (roverlap) 
     DEALLOCATE (rbecp)
  ELSE
     DEALLOCATE ( becp)
  END IF 
  DEALLOCATE (overlap) 
  DEALLOCATE (wfcatom) 
  IF (.NOT. lda_plus_u) DEALLOCATE (swfcatom) 
  ! 
  !   vectors et and proj are distributed across the pools 
  !   collect data for all k-points to the first pool
  ! 
  CALL poolrecover (et,       nbnd, nkstot, nks) 
  CALL poolrecover (proj,     nbnd * natomwfc, nkstot, nks) 
  CALL poolrecover (proj_aux, 2 * nbnd * natomwfc, nkstot, nks)
  ! 
  IF ( ionode ) THEN 
     ! 
     ! write on the file filproj
     ! 
     IF (filproj.NE.' ') THEN
        DO is=1,nspin
           IF (nspin==2) THEN
              if (is==1) filename=TRIM(filproj)//'.up'
              if (is==2) filename=TRIM(filproj)//'.down'
              nksinit=(nkstot/2)*(is-1)+1
              nkslast=(nkstot/2)*is
           ELSE
              filename=TRIM(filproj)
              nksinit=1
              nkslast=nkstot
           END IF
           iunproj=33
           CALL write_io_header(filename, iunproj, title, nrx1, nrx2, nrx3, &
                nr1, nr2, nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, &
                ecutwfc, nkstot/nspin,nbnd,natomwfc) 
           DO nwfc = 1, natomwfc 
              WRITE(iunproj,'(2i5,a3,3i5)') & 
                  nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), & 
                  nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m 
              DO ik=nksinit,nkslast
                 DO ibnd=1,nbnd
                   WRITE(iunproj,'(2i8,f20.10)') ik,ibnd, &
                                                 abs(proj(nwfc,ibnd,ik)) 
                 END DO
              END DO
           END DO
           CLOSE(iunproj)
        END DO
     END IF

     !
     ! write projections to file using iotk
     !
     CALL write_proj( "atomic_proj.xml", proj_aux )

     ! 
     ! write on the standard output file 
     ! 
     WRITE( stdout,'(/5x,"Atomic states used for projection")') 
     WRITE( stdout,'( 5x,"(read from pseudopotential files):"/)') 
     DO nwfc = 1, natomwfc 
        WRITE(stdout,1000) & 
             nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), & 
             nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m 
     END DO
1000 FORMAT (5x,"state #",i4,": atom ",i3," (",a3,"), wfc ",i2, &
                " (l=",i1," m=",i2,")") 
     ! 
     ALLOCATE(idx(natomwfc), proj1 (natomwfc) ) 
     DO ik = 1, nkstot 
        WRITE( stdout, '(/" k = ",3f14.10)') (xk (i, ik) , i = 1, 3) 
        DO ibnd = 1, nbnd 
           WRITE( stdout, '(5x,"e = ",f11.5," eV")') et (ibnd, ik) * rytoev
           ! 
           ! sort projections by magnitude, in decreasing order 
           ! 
           DO nwfc = 1, natomwfc 
              idx (nwfc) = 0 
              proj1 (nwfc) = - proj (nwfc, ibnd, ik) 
           END DO
           !
           ! projections differing by less than 1.d-4 are considered equal
           !
           CALL hpsort_eps (natomwfc, proj1, idx, eps4) 
           ! 
           !  only projections that are larger than 0.001 are written 
           ! 
           DO nwfc = 1, natomwfc 
              proj1 (nwfc) = - proj1(nwfc) 
              IF ( ABS (proj1(nwfc)) < 0.001d0 ) go to 20 
           END DO
           nwfc = natomwfc + 1 
20         nwfc = nwfc -1 
           ! 
           ! fancy (?!?) formatting 
           ! 
           WRITE( stdout, '(5x,"psi = ",5(f5.3,"*[#",i4,"]+"))') & 
                (proj1 (i), idx(i), i = 1, MIN(5,nwfc)) 
           DO j = 1, (nwfc-1)/5 
              WRITE( stdout, '(10x,"+",5(f5.3,"*[#",i4,"]+"))') & 
                   (proj1 (i), idx(i), i = 5*j+1, MIN(5*(j+1),nwfc)) 
           END DO
           psum = 0.d0 
           DO nwfc = 1, natomwfc 
              psum = psum + proj (nwfc, ibnd, ik) 
           END DO
           WRITE( stdout, '(4x,"|psi|^2 = ",f5.3)') psum 
           ! 
        ENDDO
     ENDDO
     DEALLOCATE (idx, proj1) 
     ! 
     ! estimate partial charges (Loewdin) on each atom 
     ! 
     ALLOCATE ( charges (nat, 0:lmax_wfc, nspin ) ) 
     charges = 0.0d0 
     DO ik = 1, nkstot 
        IF ( nspin == 1 ) THEN
           current_spin = 1
        ELSE IF ( nspin == 2 ) THEN
           current_spin = isk ( ik )
        ELSE
           CALL errore ('projwfc_nc',' called in the wrong case ',1)
        END IF
        DO ibnd = 1, nbnd 
           DO nwfc = 1, natomwfc 
              na= nlmchi(nwfc)%na 
              l = nlmchi(nwfc)%l 
              charges(na,l,current_spin) = charges(na,l,current_spin) + &
                   wg (ibnd,ik) * proj (nwfc, ibnd, ik) 
           ENDDO 
        END DO 
     END DO 
     ! 
     WRITE( stdout, '(/"Lowdin Charges: "/)') 
     ! 
     DO na = 1, nat 
        DO is = 1, nspin
           totcharge(is) = SUM(charges(na,0:lmax_wfc,is))
        END DO
        IF ( nspin == 1) THEN
           WRITE( stdout, 2000) na, totcharge(1), &
                ( l_label(l), charges(na,l,1), l= 0,lmax_wfc)
        ELSE IF ( nspin == 2) THEN 
           WRITE( stdout, 2000) na, totcharge(1) + totcharge(2), &
                ( l_label(l), charges(na,l,1) + charges(na,l,2), l=0,lmax_wfc)
           WRITE( stdout, 2001) totcharge(1), &
                ( l_label(l), charges(na,l,1), l= 0,lmax_wfc) 
           WRITE( stdout, 2002) totcharge(2), &
                ( l_label(l), charges(na,l,2), l= 0,lmax_wfc) 
           WRITE( stdout, 2003) totcharge(1) - totcharge(2), &
                ( l_label(l), charges(na,l,1) - charges(na,l,2), l=0,lmax_wfc)
        END IF
     END DO 
2000 FORMAT (5x,"Atom # ",i3,": total charge = ",f8.4,4(", ",a1," =",f8.4))
2001 FORMAT (15x,"  spin up      = ",f8.4,4(", ",a1," =",f8.4))
2002 FORMAT (15x,"  spin down    = ",f8.4,4(", ",a1," =",f8.4))
2003 FORMAT (15x,"  polarization = ",f8.4,4(", ",a1," =",f8.4))
     !
     psum = SUM(charges(:,:,:)) / nelec 
     WRITE( stdout, '(5x,"Spilling Parameter: ",f8.4)') 1.0d0 - psum 
     ! 
     ! Sanchez-Portal et al., Sol. State Commun.  95, 685 (1995). 
     ! The spilling parameter measures the ability of the basis provided by 
     ! the pseudo-atomic wfc to represent the PW eigenstates, 
     ! by measuring how much of the subspace of the Hamiltonian 
     ! eigenstates falls outside the subspace spanned by the atomic basis 
     ! 
     DEALLOCATE (charges) 
     !
  END IF
  ! 
  RETURN
  !
END SUBROUTINE projwave
!
!----------------------------------------------------------------------- 
SUBROUTINE projwave_nc(filproj, lsym )
  !----------------------------------------------------------------------- 
  ! 
  USE io_global,  ONLY : stdout, ionode
  USE ions_base, ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE basis,     ONLY : natomwfc
  USE printout_base, ONLY: title
  USE cell_base 
  USE constants, ONLY: rytoev, eps4 
  USE gvect 
  USE klist, ONLY: xk, nks, nkstot, nelec
  USE ldaU 
  USE lsda_mod, ONLY: nspin
  USE noncollin_module, ONLY: noncolin, npol, angle1, angle2
  USE symme, ONLY: nsym, irt, t_rev 
  USE wvfct 
  USE control_flags, ONLY: gamma_only 
  USE uspp, ONLY: nkb, vkb
  USE uspp_param, ONLY: upf
  USE becmod,   ONLY: becp_nc, calbec
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc 
  USE wavefunctions_module, ONLY: evc 
  USE mp_global, ONLY : intra_pool_comm
  USE mp,        ONLY : mp_sum
  !
  USE spin_orb,   ONLY: lspinorb, domag
  USE projections_nc
  !
  IMPLICIT NONE 
  !
  CHARACTER(LEN=*) :: filproj
  LOGICAL :: lsym
  INTEGER :: ik, ibnd, i, j, k, na, nb, nt, isym, ind, n, m, m1, n1, &
             n2, l, nwfc, nwfc1, lmax_wfc, is, nspin0, iunproj,    &
             ind0
  REAL(DP) :: jj
  REAL(DP), ALLOCATABLE :: e (:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:)
  COMPLEX(DP), ALLOCATABLE :: overlap(:,:), work(:,:),work1(:), proj0(:,:)
  ! Some workspace for k-point calculation ... 
  REAL(DP), ALLOCATABLE :: charges(:,:,:), proj1 (:)
  REAL(DP) :: psum, totcharge(2), fact(2), spinor, compute_mj
  INTEGER, ALLOCATABLE :: idx(:) 
  !
  COMPLEX(DP) :: d12(2, 2, 48), d32(4, 4, 48), d52(6, 6, 48), &
                      d72(8, 8, 48)
  COMPLEX(DP) :: d012(2, 2, 48), d112(6, 6, 48), d212(10, 10, 48), &
                      d312(14, 14, 48)
  ! 
  ! 
  ! 
  IF (.NOT.noncolin) CALL errore('projwave_nc','called in the wrong case',1) 
  IF (gamma_only) CALL errore('projwave_nc','gamma_only not yet implemented',1)
  WRITE( stdout, '(/5x,"Calling projwave_nc .... ")') 
  ! 
  ! fill structure nlmchi 
  ! 
  ALLOCATE (nlmchi(natomwfc)) 
  nwfc=0 
  lmax_wfc = 0 
  DO na = 1, nat 
     nt = ityp (na) 
     n2 = 0
     DO n = 1, upf(nt)%nwfc
        IF (upf(nt)%oc (n) >= 0.d0) THEN 
           l = upf(nt)%lchi (n) 
           lmax_wfc = MAX (lmax_wfc, l )
           IF (lspinorb) THEN 
              IF (upf(nt)%has_so) THEN     
                jj = upf(nt)%jchi (n)
                ind = 0
                DO m = -l-1, l
                   fact(1) = spinor(l,jj,m,1)
                   fact(2) = spinor(l,jj,m,2)
                   IF (abs(fact(1)) > 1.d-8 .or. abs(fact(2)) > 1.d-8) THEN 
                      nwfc = nwfc + 1 
                      ind = ind + 1
                      nlmchi(nwfc)%na = na 
                      nlmchi(nwfc)%n  =  n 
                      nlmchi(nwfc)%l  =  l 
                      nlmchi(nwfc)%m  =  m 
                      nlmchi(nwfc)%ind  =  ind  
                      nlmchi(nwfc)%jj  =  jj
                   ENDIF
                ENDDO 
              ELSE
                DO n1 = l, l+1
                  jj= DBLE(n1) - 0.5d0
                  ind = 0
                  IF (jj.gt.0.d0)  THEN 
                    n2 = n2 + 1 
                    DO m = -l-1, l
                      fact(1) = spinor(l,jj,m,1)
                      fact(2) = spinor(l,jj,m,2)
                      IF (abs(fact(1)) > 1.d-8 .or. abs(fact(2)) > 1.d-8) THEN 
                        nwfc = nwfc + 1 
                        ind = ind + 1
                        nlmchi(nwfc)%na = na 
                        nlmchi(nwfc)%n  =  n2 
                        nlmchi(nwfc)%l  =  l 
                        nlmchi(nwfc)%m  =  m 
                        nlmchi(nwfc)%ind  =  ind  
                        nlmchi(nwfc)%jj  =  jj
                      ENDIF
                    ENDDO  
                  ENDIF
                ENDDO
              ENDIF
           ELSE 
              DO m = 1, 2 * l + 1 
                 nwfc=nwfc+1 
                 nlmchi(nwfc)%na = na 
                 nlmchi(nwfc)%n  =  n 
                 nlmchi(nwfc)%l  =  l 
                 nlmchi(nwfc)%m  =  m 
                 nlmchi(nwfc)%ind  =  m 
                 nlmchi(nwfc)%jj  =  0.d0
                 nlmchi(nwfc+2*l+1)%na = na 
                 nlmchi(nwfc+2*l+1)%n  =  n 
                 nlmchi(nwfc+2*l+1)%l  =  l 
                 nlmchi(nwfc+2*l+1)%m  =  m 
                 nlmchi(nwfc+2*l+1)%ind  =  m+2*l+1 
                 nlmchi(nwfc+2*l+1)%jj  =  0.d0
              ENDDO 
              nwfc=nwfc+2*l+1
           ENDIF 
        ENDIF
     ENDDO 
  ENDDO  
  !
  IF (lmax_wfc > 3) CALL errore ('projwave_nc', 'l > 3 not yet implemented', 1)
  IF (nwfc /= natomwfc) CALL errore ('projwave_nc','wrong # of atomic wfcs?',1)
  !
  ALLOCATE(wfcatom (npwx*npol,natomwfc) )
  IF (.NOT. lda_plus_u) ALLOCATE(swfcatom (npwx*npol, natomwfc ) )  
  ALLOCATE (becp_nc ( nkb, npol, natomwfc)) 
  ALLOCATE(e (natomwfc) )  
  ALLOCATE(work (natomwfc, natomwfc) )
  !
  ALLOCATE(overlap (natomwfc, natomwfc) )
  ALLOCATE(proj0(natomwfc,nbnd), work1 (nbnd) )
  ALLOCATE(proj (natomwfc, nbnd, nkstot) )
  ALLOCATE(proj_aux (natomwfc, nbnd, nkstot) )
  overlap  = (0.d0,0.d0)
  proj0    = (0.d0,0.d0)
  proj     = 0.d0  
  proj_aux = (0.d0,0.d0)  
  ! 
  !    loop on k points 
  ! 
  CALL init_us_1 
  CALL init_at_1 
  !
  IF (lspinorb) THEN
     ! 
     ! initialize D_Sj for j=1/2, j=3/2, j=5/2 and j=7/2 
     ! 
     CALL d_matrix_so (d12, d32, d52, d72)
     !
  ELSE
     ! 
     ! initialize D_Sl for l=0, l=1, l=2 and l=3 
     ! 
     CALL d_matrix_nc (d012, d112, d212, d312)
     !
  ENDIF   
  ! 
  DO ik = 1, nks  
     wfcatom = (0.d0,0.d0)
     swfcatom= (0.d0,0.d0)
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin) 
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)
     ! 
     CALL atomic_wfc_nc_proj (ik, wfcatom)
     ! 
     CALL init_us_2 (npw, igk, xk (1, ik), vkb) 
 
     CALL calbec ( npw, vkb, wfcatom, becp_nc )
 
     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom) 
     ! 
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i> 
     ! calculate overlap matrix O_ij = <phi_i|\hat S|\phi_j> 
     ! 
     CALL ZGEMM ('C', 'N', natomwfc, natomwfc, npwx*npol, (1.d0, 0.d0), wfcatom, & 
       npwx*npol, swfcatom, npwx*npol, (0.d0, 0.d0), overlap, natomwfc)
     CALL mp_sum ( overlap, intra_pool_comm )
     ! 
     ! calculate O^{-1/2} 
     ! 
     CALL cdiagh (natomwfc, overlap, natomwfc, e, work) 
     DO i = 1, natomwfc 
        e (i) = 1.d0 / dsqrt (e (i) ) 
     ENDDO 
     DO i = 1, natomwfc 
        DO j = i, natomwfc 
          overlap (i, j) = (0.d0, 0.d0) 
          DO k = 1, natomwfc 
            overlap(i, j) = overlap(i, j) + e(k) * work(j, k) * CONJG(work (i, k) ) 
          ENDDO 
          IF (j /= i) overlap (j, i) = CONJG (overlap (i, j)) 
        ENDDO 
     ENDDO
     ! 
     ! calculate wfcatom = O^{-1/2} \hat S | phi> 
     ! 
     CALL ZGEMM ('n', 't', npwx*npol, natomwfc, natomwfc, (1.d0, 0.d0) , & 
     swfcatom, npwx*npol,  overlap, natomwfc, (0.d0, 0.d0), wfcatom, npwx*npol) 
     ! 
     ! make the projection <psi_i| O^{-1/2} \hat S | phi_j> 
     !
     CALL ZGEMM ('C','N',natomwfc, nbnd, npwx*npol, (1.d0, 0.d0), wfcatom, & 
                 npwx*npol, evc, npwx*npol, (0.d0, 0.d0), proj0, natomwfc)
     CALL mp_sum ( proj0( :, 1:nbnd ), intra_pool_comm )
     !
     proj_aux(:,:,ik) = proj0(:,:)

     !
     IF (lsym) THEN
        DO nwfc = 1, natomwfc 
           ! 
           !  atomic wavefunction nwfc is on atom na 
           !
           IF (lspinorb) THEN 
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
                    IF (nlmchi(nwfc1)%na == nb             .AND. & 
                         nlmchi(nwfc1)%n == nlmchi(nwfc)%n .AND. & 
                         nlmchi(nwfc1)%l == nlmchi(nwfc)%l .AND. &
                         nlmchi(nwfc1)%jj == nlmchi(nwfc)%jj .AND. & 
                         nlmchi(nwfc1)%ind == 1 ) go to 10 
                 END DO 
                 CALL errore('projwave_nc','cannot symmetrize',1) 
10               nwfc1=nwfc1-1 
                 ! 
                 !  nwfc1 is the first rotated atomic wfc corresponding to nwfc 
                 ! 
                 IF (abs(jj-0.5d0).lt.1.d-8) THEN 
                    work1(:) = 0.d0  
                    DO m1 = 1, 2   
                       work1(:)=work1(:)+d12(m1,ind,isym)*proj0(nwfc1+m1,:) 
                    ENDDO 
                 ELSE IF (abs(jj-1.5d0).lt.1.d-8) THEN  
                    work1(:) = 0.d0   
                    DO m1 = 1, 4   
                       work1(:)=work1(:)+d32(m1,ind,isym)*proj0(nwfc1 + m1,:) 
                    ENDDO 
                 ELSE IF (abs(jj-2.5d0).lt.1.d-8) THEN  
                    work1(:) = 0.d0   
                    DO m1 = 1, 6   
                       work1(:)=work1(:)+d52(m1,ind,isym)*proj0(nwfc1+m1,:) 
                    ENDDO 
                 ELSE IF (abs(jj-3.5d0).lt.1.d-8) THEN  
                    work1(:) = 0.d0   
                    DO m1 = 1, 8   
                       work1(:)=work1(:)+d72(m1,ind,isym)*proj0(nwfc1+m1,:) 
                    ENDDO 
                 ENDIF 
                 DO ibnd = 1, nbnd 
                    proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + & 
                         work1(ibnd) * CONJG (work1(ibnd)) / nsym 
                 ENDDO 
                 ! on symmetries
!--  in a nonmagnetic case - another loop with the time reversal
                 IF (.not.domag.and.ind.eq.ind0) THEN
                   ind = 2*jj + 2 - ind0
                   nwfc1 = nwfc1 + 1
                   GOTO 10
                 ENDIF
!--
              ENDDO 
!--  in a nonmagnetic case - rescale
              IF (.not.domag) THEN
                DO ibnd = 1, nbnd
                  proj(nwfc,ibnd,ik) = 0.5d0*proj(nwfc,ibnd,ik)
                ENDDO
              ENDIF
!--
           ELSE
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
                    IF (nlmchi(nwfc1)%na == nb             .AND. & 
                        nlmchi(nwfc1)%n == nlmchi(nwfc)%n .AND. & 
                        nlmchi(nwfc1)%l == nlmchi(nwfc)%l .AND. & 
                        nlmchi(nwfc1)%m == 1 .AND. & 
                        nlmchi(nwfc1)%ind == 1) go to 15 
                 END DO 
                 CALL errore('projwave_nc','cannot symmetrize',1) 
15               nwfc1=nwfc1-1 
                 IF (l == 0) THEN  
                    work1(:) = 0.d0  
                    DO m1 = 1, 2   
                       work1(:) = work1(:) + d012 (m1, ind, isym) * & 
                                  proj0 (nwfc1 + m1,:)
                    ENDDO 
                 ELSE IF (l == 1) THEN  
                    work1(:) = 0.d0   
                    DO m1 = 1, 6        
                       work1(:) = work1(:) + d112 (m1, ind, isym) * &
                                  proj0 (nwfc1 + m1,:)
                    ENDDO 
                 ELSE IF (l == 2) THEN  
                    work1(:) = 0.d0   
                    DO m1 = 1, 10     
                       work1(:) = work1(:) + d212 (m1, ind, isym) * &
                                  proj0 (nwfc1 + m1,:) 
                    ENDDO 
                 ELSE IF (l == 3) THEN  
                    work1(:) = 0.d0   
                    DO m1 = 1, 14     
                       work1(:) = work1(:) + d312 (m1, ind, isym) * &
                                  proj0 (nwfc1 + m1,:)
                    END DO 
                 END IF 
                 DO ibnd = 1, nbnd 
                    proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + & 
                         work1(ibnd) * CONJG (work1(ibnd)) / nsym 
                 ENDDO
                 ! on symmetries
              ENDDO
           ENDIF
           ! on atomic wavefunctions  
        ENDDO
     ELSE
        DO nwfc=1,natomwfc
           DO ibnd=1,nbnd
              proj(nwfc,ibnd,ik)=ABS(proj0(nwfc,ibnd))**2
           END DO
        END DO
     END IF
     ! on k-points 
  ENDDO 
  !
  DEALLOCATE (work) 
  DEALLOCATE (work1) 
  DEALLOCATE (proj0) 
  DEALLOCATE (e)
  DEALLOCATE (becp_nc)
  DEALLOCATE (overlap) 
  DEALLOCATE (wfcatom) 
  IF (.NOT. lda_plus_u) DEALLOCATE (swfcatom) 
  ! 
  !   vectors et and proj are distributed across the pools 
  !   collect data for all k-points to the first pool
  ! 
  CALL poolrecover (et,       nbnd, nkstot, nks) 
  CALL poolrecover (proj,     nbnd * natomwfc, nkstot, nks) 
  CALL poolrecover (proj_aux, 2 * nbnd * natomwfc, nkstot, nks) 
  ! 
  IF ( ionode ) THEN 
     ! 
     ! write on the file filproj
     ! 
     IF (filproj.ne.' ') THEN
        iunproj=33
        CALL write_io_header(filproj, iunproj, title, nrx1, nrx2, nrx3, &
           nr1, nr2, nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
           nkstot,nbnd,natomwfc) 
        DO nwfc = 1, natomwfc 
           IF (lspinorb) THEN
              WRITE(iunproj,1000) &
                   nwfc, nlmchi(nwfc)%na,atm(ityp(nlmchi(nwfc)%na)), &
                   nlmchi(nwfc)%n,nlmchi(nwfc)%jj,nlmchi(nwfc)%l, &
                   compute_mj(nlmchi(nwfc)%jj,nlmchi(nwfc)%l,nlmchi(nwfc)%m)
           ELSE
              WRITE(iunproj,1500) &
                   nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                   nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m, &
                   0.5d0-INT(nlmchi(nwfc)%ind/(2*nlmchi(nwfc)%l+2))
           END IF
           DO ik=1,nkstot
              DO ibnd=1,nbnd
                 WRITE(iunproj,'(2i8,f20.10)') ik,ibnd,abs(proj(nwfc,ibnd,ik)) 
              END DO
           END DO
        END DO
        CLOSE(iunproj)
     END IF

     !
     ! write projections to file using iotk
     !
     CALL write_proj( "atomic_proj.xml", proj_aux )

     ! 
     ! write on the standard output file 
     ! 
     WRITE( stdout,'(/5x,"Atomic states used for projection")') 
     WRITE( stdout,'( 5x,"(read from pseudopotential files):"/)')              
     IF (lspinorb) THEN
        DO nwfc = 1, natomwfc 
           WRITE(stdout,1000) & 
             nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), & 
             nlmchi(nwfc)%n, nlmchi(nwfc)%jj, nlmchi(nwfc)%l,   &
             compute_mj(nlmchi(nwfc)%jj,nlmchi(nwfc)%l,nlmchi(nwfc)%m)
        END DO
1000    FORMAT (5x,"state #",i3,": atom ",i3," (",a3,"), wfc ",i2, &
                   " (j=",f3.1," l=",i1," m_j=",f4.1,")")
     ELSE 
        DO nwfc = 1, natomwfc 
           WRITE(stdout,1500) & 
             nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), & 
             nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m, &
             0.5d0-INT(nlmchi(nwfc)%ind/(2*nlmchi(nwfc)%l+2))
        END DO
1500    FORMAT (5x,"state #",i3,": atom ",i3," (",a3,"), wfc ",i2, &
                   " (l=",i1," m=",i2," s_z=",f4.1,")") 
     ENDIF
     ! 
     ALLOCATE(idx (natomwfc), proj1 (natomwfc) ) 
     DO ik = 1, nkstot 
        WRITE( stdout, '(/" k = ",3f14.10)') (xk (i, ik) , i = 1, 3) 
        DO ibnd = 1, nbnd 
           WRITE( stdout, '(5x,"e = ",f11.5," eV")') et (ibnd, ik) * rytoev 
           ! 
           ! sort projections by magnitude, in decreasing order 
           ! 
           DO nwfc = 1, natomwfc 
              idx (nwfc) = 0 
              proj1 (nwfc) = - proj (nwfc, ibnd, ik) 
           END DO
           CALL hpsort_eps (natomwfc, proj1, idx, eps4) 
           ! 
           !  only projections that are larger than 0.001 are written 
           ! 
           DO nwfc = 1, natomwfc 
              proj1 (nwfc) = - proj1(nwfc) 
              IF ( ABS (proj1(nwfc)) < 0.001d0 ) go to 20 
           END DO
           nwfc = natomwfc + 1 
20         nwfc = nwfc -1 
           ! 
           ! fancy (?!?) formatting 
           ! 
           WRITE( stdout, '(5x,"psi = ",5(f5.3,"*[#",i3,"]+"))') & 
                (proj1 (i), idx(i), i = 1, MIN(5,nwfc)) 
           DO j = 1, (nwfc-1)/5 
              WRITE( stdout, '(10x,"+",5(f5.3,"*[#",i3,"]+"))') & 
                   (proj1 (i), idx(i), i = 5*j+1, MIN(5*(j+1),nwfc)) 
           END DO
           psum = 0.d0 
           DO nwfc = 1, natomwfc 
              psum = psum + proj (nwfc, ibnd, ik) 
           END DO
           WRITE( stdout, '(4x,"|psi|^2 = ",f5.3)') psum 
           ! 
        ENDDO
     ENDDO
     DEALLOCATE (idx, proj1) 
     ! 
     ! estimate partial charges (Loewdin) on each atom 
     ! 
     IF (lspinorb) THEN
        nspin0 = 1
     ELSE
        nspin0 = 2
     END IF
     ALLOCATE ( charges (nat, 0:lmax_wfc, nspin0 ) ) 
     charges = 0.0d0 
     IF (lspinorb) THEN
        DO ik = 1, nkstot 
           DO ibnd = 1, nbnd 
              DO nwfc = 1, natomwfc 
                 na= nlmchi(nwfc)%na 
                 l = nlmchi(nwfc)%l 
                 charges(na,l,1) = charges(na,l,1) + &
                      wg (ibnd,ik) * proj (nwfc, ibnd, ik)
              ENDDO 
           END DO 
        END DO 
     ELSE 
        DO ik = 1, nkstot 
           DO ibnd = 1, nbnd 
              DO nwfc = 1, natomwfc 
                 na= nlmchi(nwfc)%na 
                 l = nlmchi(nwfc)%l
                 IF ( nlmchi(nwfc)%ind.le.(2*l+1)) THEN
                    charges(na,l,1) = charges(na,l,1) + &
                         wg (ibnd,ik) * proj (nwfc, ibnd, ik) 
                 ELSE
                    charges(na,l,2) = charges(na,l,2) + &
                         wg (ibnd,ik) * proj (nwfc, ibnd, ik) 
                 ENDIF
              END DO 
           END DO 
        END DO
     END IF 
     ! 
     WRITE( stdout, '(/"Lowdin Charges: "/)') 
     ! 
     DO na = 1, nat 
        DO is = 1, nspin0
           totcharge(is) = SUM(charges(na,0:lmax_wfc,is))
        END DO
        IF ( nspin0 == 1) THEN
           WRITE( stdout, 2000) na, totcharge(1), &
                ( charges(na,l,1), l= 0,lmax_wfc)
        ELSE IF ( nspin0 == 2) THEN 
           WRITE( stdout, 2000) na, totcharge(1) + totcharge(2), &
                ( charges(na,l,1) + charges(na,l,2), l=0,lmax_wfc)
           WRITE( stdout, 2001) totcharge(1), &
                ( charges(na,l,1), l= 0,lmax_wfc) 
           WRITE( stdout, 2002) totcharge(2), &
                ( charges(na,l,2), l= 0,lmax_wfc) 
           WRITE( stdout, 2003) totcharge(1) - totcharge(2), &
                 ( charges(na,l,1) - charges(na,l,2), l=0,lmax_wfc)
        END IF
     END DO 
2000 FORMAT (5x,"Atom # ",i3,": total charge = ",f8.4 ,&
          & ", s, p, d, f = ",4f8.4) 
2001 FORMAT (15x,"  spin up      = ",f8.4 , &
          & ", s, p, d, f = ",4f8.4) 
2002 FORMAT (15x,"  spin down    = ",f8.4 , &
          & ", s, p, d, f = ",4f8.4) 
2003 FORMAT (15x,"  polarization = ",f8.4 , &
          & ", s, p, d, f = ",4f8.4) 
     !
     psum = SUM(charges(:,:,:)) / nelec 
     WRITE( stdout, '(5x,"Spilling Parameter: ",f8.4)') 1.0d0 - psum 
     ! 
     ! Sanchez-Portal et al., Sol. State Commun.  95, 685 (1995). 
     ! The spilling parameter measures the ability of the basis provided by 
     ! the pseudo-atomic wfc to represent the PW eigenstates, 
     ! by measuring how much of the subspace of the Hamiltonian 
     ! eigenstates falls outside the subspace spanned by the atomic basis 
     ! 
     DEALLOCATE (charges) 
     !
  END IF
  ! 
  RETURN
  !
END SUBROUTINE projwave_nc
!
!----------------------------------------------------------------------- 
SUBROUTINE  partialdos (Emin, Emax, DeltaE, filpdos)
  !----------------------------------------------------------------------- 
  ! 
  USE io_global,  ONLY : stdout 
  USE basis, ONLY : natomwfc
  USE ions_base, ONLY : ityp, atm
  USE klist, ONLY: wk, nkstot, degauss, ngauss, lgauss
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE wvfct, ONLY: et, nbnd
  USE constants, ONLY: rytoev 
  !
  USE projections
  !
  IMPLICIT NONE 
  CHARACTER (len=256) :: filpdos
  REAL(DP) :: Emin, Emax, DeltaE
  !
  CHARACTER (len=33) :: filextension 
  CHARACTER (len=256):: fileout
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/) 
  ! 
  INTEGER :: ik, ibnd,  m, &
       c_tab, nwfc, ne, ie_mid, ie_delta, ie, is 
  REAL(DP) :: etev, delta, Elw, Eup 
  REAL(DP), ALLOCATABLE :: dostot(:,:), pdos(:,:,:), pdostot(:,:), &
       ldos(:,:)
  REAL(DP), EXTERNAL :: w0gauss 
  ! 
  ! 
  ! find band extrema 
  ! 
  Elw = et (1, 1) 
  Eup = et (nbnd, 1) 
  DO ik = 2, nkstot 
     Elw = MIN (Elw, et (1, ik) ) 
     Eup = MAX (Eup, et (nbnd, ik) ) 
  ENDDO
  IF (degauss.NE.0.d0) THEN
     Eup = Eup + 3d0 * degauss
     Elw = Elw - 3d0 * degauss
  ENDIF
  Emin = MAX (Emin/rytoev, Elw) 
  Emax = MIN (Emax/rytoev, Eup) 
  DeltaE = DeltaE/rytoev
  ne = NINT ( (Emax - Emin) / DeltaE+0.500001d0)
  !
  ALLOCATE (pdos(0:ne,natomwfc,nspin))
  ALLOCATE (dostot(0:ne,nspin), pdostot(0:ne,nspin), ldos(0:ne,nspin) )
  pdos(:,:,:) = 0.d0 
  dostot(:,:) = 0.d0 
  pdostot(:,:)= 0.d0 
  current_spin = 1
  ie_delta = 5 * degauss / DeltaE + 1 

  DO ik = 1,nkstot 
     IF ( nspin == 2 ) current_spin = isk ( ik ) 
     DO ibnd = 1, nbnd 
        etev = et(ibnd,ik)
        ie_mid = NINT( (etev-Emin)/DeltaE ) 
        DO ie = MAX(ie_mid-ie_delta, 0), MIN(ie_mid+ie_delta, ne) 
           delta = w0gauss((Emin+DeltaE*ie-etev)/degauss,ngauss) &
                 / degauss / rytoev
           !
           ! pdos(:,nwfc,ns) = DOS (states/eV) for spin "ns" 
           !                   projected over atomic wfc "nwfc"
           !
           DO nwfc = 1, natomwfc 
              pdos(ie,nwfc,current_spin) = pdos(ie,nwfc,current_spin) + & 
                   wk(ik) * delta * proj (nwfc, ibnd, ik)  
           END DO
           !
           ! dostot(:,ns) = total DOS (states/eV) for spin "ns" 
           !
           dostot(ie,current_spin) = dostot(ie,current_spin) + & 
                wk(ik) * delta 
        END DO
     END DO
  END DO
  !
  ! pdostot(:,ns) = sum of all projected DOS
  !
  DO is=1,nspin 
     DO ie=0,ne 
        pdostot(ie,is) = SUM(pdos(ie,:,is)) 
     END DO
  END DO
  
  DO nwfc = 1, natomwfc 
     IF (nlmchi(nwfc)%m == 1) THEN 
        filextension='.pdos_atm#' 
        !             12345678901 
        c_tab = 11 
        IF (nlmchi(nwfc)%na < 10) THEN 
           WRITE (filextension( c_tab : c_tab ),'(i1)') nlmchi(nwfc)%na 
           c_tab = c_tab + 1 
        ELSE IF (nlmchi(nwfc)%na < 100) THEN 
           WRITE (filextension( c_tab : c_tab+1 ),'(i2)') nlmchi(nwfc)%na 
           c_tab = c_tab + 2 
        ELSE IF (nlmchi(nwfc)%na < 1000) THEN 
           WRITE (filextension( c_tab : c_tab+2 ),'(i3)') nlmchi(nwfc)%na 
           c_tab = c_tab + 3 
        ELSE IF (nlmchi(nwfc)%na < 10000) THEN 
           WRITE (filextension( c_tab : c_tab+3 ),'(i4)') nlmchi(nwfc)%na 
           c_tab = c_tab + 4 
        ELSE 
           CALL errore('partialdos',& 
                'file extension not supporting so many atoms', nwfc) 
        ENDIF
        WRITE (filextension(c_tab:c_tab+4),'(a1,a)') & 
             '(',TRIM(atm(ityp(nlmchi(nwfc)%na))) 
        c_tab = c_tab + LEN_TRIM(atm(ityp(nlmchi(nwfc)%na))) + 1 
        IF (nlmchi(nwfc)%n >= 10) & 
             CALL errore('partialdos',& 
             'file extension not supporting so many atomic wfc', nwfc) 
        IF (nlmchi(nwfc)%l > 3) & 
             CALL errore('partialdos',& 
             'file extension not supporting so many l', nwfc) 
        WRITE (filextension(c_tab:),'(")_wfc#",i1,"(",a1,")")')  & 
             nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l) 
        fileout = TRIM(filpdos)//TRIM(filextension)
        OPEN (4,file=fileout,form='formatted', & 
             status='unknown') 

        IF (nspin == 1) THEN 
           WRITE (4,'("# E (eV)   ldos(E)  ",$)') 
        ELSE 
           WRITE (4,'("# E (eV)  ldosup(E)  ldosdw(E)",$)') 
        END IF
        DO m=1,2 * nlmchi(nwfc)%l + 1 
           IF (nspin == 1) THEN 
              WRITE(4,'(" pdos(E)   ",$)')  
           ELSE 
              WRITE(4,'(" pdosup(E) ",$)')  
              WRITE(4,'(" pdosdw(E) ",$)')  
           END IF
        END DO
        WRITE(4,*) 
        !
        ! ldos = PDOS summed over m (m=-l:+l)
        !
        ldos  (:,:) = 0.d0
        DO ie= 0, ne 
           DO is=1, nspin
              DO m=1,2 * nlmchi(nwfc)%l + 1 
                 ldos  (ie, is) = ldos  (ie, is) + pdos(ie,nwfc+m-1,is)
              END DO
           END DO
        END DO
        DO ie= 0, ne 
           etev = Emin + ie * DeltaE 
           WRITE (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev,  & 
                (ldos(ie,is), is=1,nspin), & 
                ((pdos(ie,nwfc+m-1,is), is=1,nspin), & 
                m=1,2*nlmchi(nwfc)%l+1)
        END DO
        CLOSE (4) 
     END IF
  END DO
  fileout = TRIM(filpdos)//".pdos_tot"
  OPEN (4,file=fileout,form='formatted', status='unknown') 
  IF (nspin == 1) THEN 
     WRITE (4,'("# E (eV)  dos(E)    pdos(E)")')  
  ELSE 
     WRITE (4,'("# E (eV)  dosup(E)   dosdw(E)  pdosup(E)  pdosdw(E)")')  
  END IF
  DO ie= 0, ne 
     etev = Emin + ie * DeltaE 
     WRITE (4,'(f7.3,4e11.3)') etev*rytoev, (dostot(ie,is), is=1,nspin), & 
          (pdostot(ie,is), is=1,nspin) 
  END DO
  CLOSE (4) 
  DEALLOCATE (ldos, dostot, pdostot)
  DEALLOCATE (pdos)
  !
  DEALLOCATE (nlmchi) 
  DEALLOCATE (proj) 
  DEALLOCATE (proj_aux) 
  !
  RETURN 
END SUBROUTINE partialdos
!
!----------------------------------------------------------------------- 
SUBROUTINE  partialdos_nc (Emin, Emax, DeltaE, filpdos)
  !----------------------------------------------------------------------- 
  ! 
  USE io_global,  ONLY : stdout 
  USE basis, ONLY : natomwfc
  USE ions_base, ONLY : ityp, atm
  USE klist, ONLY: wk, nkstot, degauss, ngauss, lgauss
  USE lsda_mod, ONLY: nspin
  USE wvfct, ONLY: et, nbnd
  USE constants, ONLY: rytoev 
  !
  USE spin_orb,   ONLY: lspinorb
  USE projections_nc
  !
  IMPLICIT NONE 
  CHARACTER (len=256) :: filpdos
  REAL(DP) :: Emin, Emax, DeltaE
  !
  CHARACTER (len=33) :: filextension 
  CHARACTER (len=256):: fileout
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/) 
  ! 
  INTEGER :: ik, ibnd, ind, m, &
       c_tab, nwfc, ne, ie_mid, ie_delta, ie, is, nspin0 
  REAL(DP) :: etev, delta, Elw, Eup, fact(2), spinor
  REAL(DP), ALLOCATABLE :: dostot(:), pdos(:,:,:), pdostot(:,:), &
       ldos(:,:)
  REAL(DP), EXTERNAL :: w0gauss 
  ! 
  ! 
  ! find band extrema 
  ! 
  Elw = et (1, 1) 
  Eup = et (nbnd, 1) 
  DO ik = 2, nkstot 
     Elw = MIN (Elw, et (1, ik) ) 
     Eup = MAX (Eup, et (nbnd, ik) ) 
  ENDDO
  IF (degauss.NE.0.d0) THEN
     Eup = Eup + 3d0 * degauss
     Elw = Elw - 3d0 * degauss
  ENDIF
  Emin = MAX (Emin/rytoev, Elw) 
  Emax = MIN (Emax/rytoev, Eup) 
  DeltaE = DeltaE/rytoev
  ne = NINT ( (Emax - Emin) / DeltaE+0.500001d0)
  !
  IF (lspinorb) THEN
     nspin0 = 1
  ELSE
     nspin0 = 2
  ENDIF
  ALLOCATE (pdos(0:ne,natomwfc,nspin0))
  ALLOCATE (dostot(0:ne), pdostot(0:ne,nspin0), ldos(0:ne,nspin0) )
  pdos(:,:,:) = 0.d0 
  dostot(:) = 0.d0 
  pdostot(:,:)= 0.d0 
  ie_delta = 5 * degauss / DeltaE + 1 

  DO ik = 1,nkstot 
     DO ibnd = 1, nbnd 
        etev = et(ibnd,ik)
        ie_mid = NINT( (etev-Emin)/DeltaE ) 
        DO ie = MAX(ie_mid-ie_delta, 0), MIN(ie_mid+ie_delta, ne) 
           delta = w0gauss((Emin+DeltaE*ie-etev)/degauss,ngauss) &
                 / degauss / rytoev
           !
           ! pdos(:,nwfc,ns) = DOS (states/eV) for spin "ns" 
           !                   projected over atomic wfc "nwfc"
           !
           !
           ! dostot(:) = total DOS (states/eV) 
           !
           IF (lspinorb) THEN
              DO nwfc = 1, natomwfc 
                 pdos(ie,nwfc,1) = pdos(ie,nwfc,1) + & 
                     wk(ik) * delta * proj (nwfc, ibnd, ik)  
              END DO
              dostot(ie) = dostot(ie) + wk(ik) * delta 
           ELSE
              DO nwfc = 1, natomwfc 
                 IF ( nlmchi(nwfc)%ind.le.(2* nlmchi(nwfc)%l+1)) THEN
                    pdos(ie,nwfc,1) = pdos(ie,nwfc,1) + & 
                        wk(ik) * delta * proj (nwfc, ibnd, ik)
                    pdos(ie,nwfc,2) = 0.d0
                 ELSE
                    pdos(ie,nwfc,1) = 0.d0
                    pdos(ie,nwfc,2) = pdos(ie,nwfc,2) + & 
                        wk(ik) * delta * proj (nwfc, ibnd, ik)
                 END IF
              END DO 
              dostot(ie) = dostot(ie) + wk(ik) * delta
           END IF
        END DO
     END DO
  END DO
  !
  ! pdostot(:,ns) = sum of all projected DOS
  !
  DO is=1,nspin0 
     DO ie=0,ne 
        pdostot(ie,is) = SUM(pdos(ie,:,is)) 
     END DO
  END DO
  
  DO nwfc = 1, natomwfc 
     IF (nlmchi(nwfc)%ind == 1) THEN 
        filextension='.pdos_atm#' 
        !             12345678901 
        c_tab = 11 
        IF (nlmchi(nwfc)%na < 10) THEN 
           WRITE (filextension( c_tab : c_tab ),'(i1)') nlmchi(nwfc)%na 
           c_tab = c_tab + 1 
        ELSE IF (nlmchi(nwfc)%na < 100) THEN 
           WRITE (filextension( c_tab : c_tab+1 ),'(i2)') nlmchi(nwfc)%na 
           c_tab = c_tab + 2 
        ELSE IF (nlmchi(nwfc)%na < 1000) THEN 
           WRITE (filextension( c_tab : c_tab+2 ),'(i3)') nlmchi(nwfc)%na 
           c_tab = c_tab + 3 
        ELSE IF (nlmchi(nwfc)%na < 10000) THEN 
           WRITE (filextension( c_tab : c_tab+3 ),'(i4)') nlmchi(nwfc)%na 
           c_tab = c_tab + 4 
        ELSE 
           CALL errore('partialdos_nc',& 
                'file extension not supporting so many atoms', nwfc) 
        ENDIF
        WRITE (filextension(c_tab:c_tab+4),'(a1,a)') & 
             '(',TRIM(atm(ityp(nlmchi(nwfc)%na))) 
        c_tab = c_tab + LEN_TRIM(atm(ityp(nlmchi(nwfc)%na))) + 1 
        IF (nlmchi(nwfc)%n >= 10) & 
             CALL errore('partialdos_nc',& 
             'file extension not supporting so many atomic wfc', nwfc) 
        IF (nlmchi(nwfc)%l > 3) & 
             CALL errore('partialdos_nc',& 
             'file extension not supporting so many l', nwfc) 
        IF (lspinorb) THEN
           WRITE (filextension(c_tab:),'(")_wfc#",i1,"(",a1,"_j",f3.1,")")') & 
             nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l),nlmchi(nwfc)%jj
        ELSE
           WRITE (filextension(c_tab:),'(")_wfc#",i1,"(",a1,")")')  & 
             nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l)
        ENDIF   
        fileout = TRIM(filpdos)//TRIM(filextension)
        OPEN (4,file=fileout,form='formatted', & 
             status='unknown') 

        IF (nspin0 == 1) THEN 
           WRITE (4,'("# E(eV)   ldos(E)   ",$)') 
        ELSE 
           WRITE (4,'("# E(eV)  ldosup(E)  ldosdw(E)",$)') 
        END IF 
        IF (lspinorb) THEN
           ind = 0
           DO m = -nlmchi(nwfc)%l-1, nlmchi(nwfc)%l
              fact(1) = spinor(nlmchi(nwfc)%l,nlmchi(nwfc)%jj,m,1)
              fact(2) = spinor(nlmchi(nwfc)%l,nlmchi(nwfc)%jj,m,2)
              IF (abs(fact(1)).gt.1.d-8.or.abs(fact(2)).gt.1.d-8) THEN          
                 ind = ind + 1
                 WRITE(4,'("pdos(E)_",i1,"   ",$)') ind 
              END IF 
           END DO
        ELSE
           DO ind=1,2 * nlmchi(nwfc)%l + 1 
              WRITE(4,'(" pdosup(E) ",$)')  
              WRITE(4,'(" pdosdw(E) ",$)') 
           END DO
        ENDIF
        WRITE(4,*) 
        !
        ! ldos = PDOS summed over m (m=-l:+l)
        !
        ldos  (:,:) = 0.d0
        IF (lspinorb) THEN
           DO ie= 0, ne 
              IF (abs(nlmchi(nwfc)%jj-nlmchi(nwfc)%l-0.5d0).lt.1.d-8) THEN        
                 DO ind = 1, 2 * nlmchi(nwfc)%l + 2  
                    ldos  (ie, 1) = ldos  (ie, 1) + pdos(ie,nwfc+ind-1,1)
                 END DO
              ELSEIF (abs(nlmchi(nwfc)%jj-nlmchi(nwfc)%l+0.5d0).lt.1.d-8) THEN       
                 DO ind = 1, 2 * nlmchi(nwfc)%l   
                    ldos  (ie, 1) = ldos  (ie, 1) + pdos(ie,nwfc+ind-1,1)
                 END DO
              END IF
           END DO
           DO ie= 0, ne 
              etev = Emin + ie * DeltaE 
              IF (abs(nlmchi(nwfc)%jj-nlmchi(nwfc)%l-0.5d0).lt.1.d-8) THEN 
                 WRITE (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev, ldos(ie,1), & 
                      (pdos(ie,nwfc+ind-1,1), ind=1,2*nlmchi(nwfc)%l+2)
              ELSEIF (abs(nlmchi(nwfc)%jj-nlmchi(nwfc)%l+0.5d0).lt.1.d-8) THEN
                 WRITE (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev, ldos(ie,1), & 
                      (pdos(ie,nwfc+ind-1,1), ind=1,2*nlmchi(nwfc)%l)
              END IF
           END DO
        ELSE
           DO ie= 0, ne 
              DO is=1, nspin0
                 DO ind=1,4 * nlmchi(nwfc)%l + 2 
                    ldos  (ie, is) = ldos  (ie, is) + pdos(ie,nwfc+ind-1,is)
                 END DO
              END DO
           END DO
           DO ie= 0, ne 
              etev = Emin + ie * DeltaE 
              WRITE (4,'(f7.3,2e11.3,14e11.3)') etev*rytoev,  & 
                   (ldos(ie,is), is=1,nspin0), & 
                   ((pdos(ie,nwfc+ind-1+(is-1)*(2*nlmchi(nwfc)%l+1),is), is=1,nspin0), & 
                   ind=1,2*nlmchi(nwfc)%l+1)
           END DO 
        END IF
        CLOSE (4) 
     END IF
  END DO
  fileout = TRIM(filpdos)//".pdos_tot"
  OPEN (4,file=fileout,form='formatted', status='unknown') 
  IF (nspin0 == 1) THEN 
     WRITE (4,'("# E (eV)  dos(E)    pdos(E)")')  
  ELSE 
     WRITE (4,'("# E (eV)  dos(E)   pdosup(E)  pdosdw(E)")')  
  END IF
  DO ie= 0, ne 
     etev = Emin + ie * DeltaE 
     WRITE (4,'(f7.3,4e11.3)') etev*rytoev, dostot(ie), & 
          (pdostot(ie,is), is=1,nspin0) 
  END DO
  CLOSE (4) 
  DEALLOCATE (ldos, dostot, pdostot)
  DEALLOCATE (pdos)
  !
  DEALLOCATE (nlmchi) 
  DEALLOCATE (proj) 
  DEALLOCATE (proj_aux) 
  !
  RETURN 
END SUBROUTINE partialdos_nc
!
!-----------------------------------------------------------------------
SUBROUTINE write_io_header(filplot, iunplot, title, nrx1, nrx2, nrx3, &
           nr1, nr2, nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
           nkstot,nbnd,natomwfc) 
   !-----------------------------------------------------------------------

   USE kinds, ONLY: DP
   USE ions_base, ONLY : zv, atm, tau, ityp
   USE noncollin_module, ONLY: noncolin
   USE spin_orb, ONLY: lspinorb 
   
   IMPLICIT NONE
   CHARACTER (LEN=*) :: filplot
   CHARACTER (LEN=*) :: title
   INTEGER :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nat, ntyp, ibrav
   REAL(DP) :: celldm (6), gcutm, dual, ecutwfc, at(3,3)
   INTEGER :: iunplot, ios, na, nt, i
   INTEGER :: nkstot,nbnd,natomwfc
   !
   if (filplot == ' ') call errore ('write_io_h', 'filename missing', 1)
   
   OPEN (UNIT = iunplot, FILE = filplot, FORM = 'formatted', &
         STATUS = 'unknown', ERR = 101, IOSTAT = ios)
   101     CALL errore ('write_io_h', 'opening file '//TRIM(filplot), abs (ios) )
   WRITE (iunplot, '(a)') title
   WRITE (iunplot, '(8i8)') nrx1, nrx2, nrx3, nr1, nr2, nr3, nat, ntyp
   WRITE (iunplot, '(i6,6f12.8)') ibrav, celldm
   IF  (ibrav == 0) THEN
       WRITE ( iunplot, * ) at(:,1)
       WRITE ( iunplot, * ) at(:,2)
       WRITE ( iunplot, * ) at(:,3)
   END IF
   WRITE (iunplot, '(3f20.10,i6)') gcutm, dual, ecutwfc, 9
   WRITE (iunplot, '(i4,3x,a2,3x,f5.2)') &
         (nt, atm (nt), zv (nt), nt=1, ntyp)
   WRITE (iunplot, '(i4,3x,3f15.9,3x,i2)') (na, &
         (tau (i, na), i = 1, 3), ityp (na), na = 1, nat)
   WRITE (iunplot, '(3i8)') natomwfc, nkstot, nbnd
   WRITE (iunplot, '(2l5)') noncolin, lspinorb
   
   RETURN
END SUBROUTINE write_io_header
!
!-----------------------------------------------------------------------
FUNCTION compute_mj(j,l,m)
   !-----------------------------------------------------------------------
   USE kinds, ONLY: DP
   IMPLICIT NONE
   !
   REAL(DP) :: compute_mj, j
   INTEGER  :: l, m

   IF (ABS(j-l-0.5d0).lt.1.d-4) THEN
       compute_mj=m+0.5d0
   ELSE IF (ABS(j-l+0.5d0).lt.1.d-4) THEN
      compute_mj=m-0.5d0
   ELSE
      call errore('compute_mj','l and j not compatible',1)
   END IF

   RETURN
END FUNCTION compute_mj
!
!-----------------------------------------------------------------------
SUBROUTINE  write_proj (filename, projs)
  !-----------------------------------------------------------------------
  !
  USE kinds
  USE io_files,         ONLY : iun => iunsat, prefix, tmp_dir
  USE basis,            ONLY : natomwfc
  USE cell_base
  USE klist,            ONLY : wk, xk, nkstot, nelec
  USE noncollin_module, ONLY : noncolin
  USE lsda_mod,         ONLY : nspin, isk
  USE ener,             ONLY : ef
  USE wvfct,            ONLY : et, nbnd
  USE iotk_module
  IMPLICIT NONE

  CHARACTER(*),  INTENT(IN) :: filename
  COMPLEX(DP),   INTENT(IN) :: projs(natomwfc,nbnd,nkstot)
  !
  CHARACTER(256)          :: tmp
  CHARACTER(iotk_attlenx) :: attr
  INTEGER :: ik, ik_eff, ia, ierr, num_k_points

!
! subroutine body
!

  tmp = TRIM( tmp_dir ) // TRIM( prefix ) // '.save/' //TRIM(filename)
  !
  CALL iotk_open_write(iun, FILE=TRIM(tmp), ROOT="ATOMIC_PROJECTIONS", IERR=ierr )
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
  CALL iotk_write_attr(attr, "UNITS", " 2 pi / a", FIRST=.TRUE.  )
  CALL iotk_write_empty (iun,  "UNITS_FOR_K-POINTS", ATTR=attr)
  CALL iotk_write_attr(attr, "UNITS", "Rydberg", FIRST=.TRUE.  )
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
     CALL iotk_write_begin( iun, "K-POINT"//TRIM(iotk_index(ik)) )
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
     CALL iotk_write_end( iun, "K-POINT"//TRIM(iotk_index(ik)) )
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
     CALL iotk_write_begin( iun, "K-POINT"//TRIM(iotk_index(ik)) )
     !
     IF ( nspin == 2 ) THEN
        !
        CALL iotk_write_begin ( iun, "SPIN.1" )
           !
           DO ia = 1, natomwfc
               CALL iotk_write_dat(iun, "ATMWFC"//TRIM(iotk_index(ia)), projs(ia,:,ik)  )
           ENDDO
           !
        CALL iotk_write_end ( iun, "SPIN.1" )
        !
        ik_eff = ik + num_k_points
        !
        CALL iotk_write_begin ( iun, "SPIN.2" )
           !
           DO ia = 1, natomwfc
               CALL iotk_write_dat(iun, "ATMWFC"//TRIM(iotk_index(ia)), projs(ia,:,ik_eff)  )
           ENDDO
           !
        CALL iotk_write_end ( iun, "SPIN.2" )
        !
     ELSE
        !
        DO ia = 1,natomwfc
            CALL iotk_write_dat(iun, "ATMWFC"//TRIM(iotk_index(ia)), projs(ia,:,ik)  )
        ENDDO
        !
     ENDIF
     !
     CALL iotk_write_end( iun, "K-POINT"//TRIM(iotk_index(ik)) )
     !
  ENDDO
  !
  CALL iotk_write_end(iun, "PROJECTIONS")

  !
  ! closing the file
  !
  CALL iotk_close_write(iun)

END SUBROUTINE write_proj


!
!  projwave with distributed matrixes
!
!----------------------------------------------------------------------- 
SUBROUTINE pprojwave( filproj, lsym )
  !----------------------------------------------------------------------- 
  ! 
  USE io_global, ONLY : stdout, ionode
  USE printout_base, ONLY: title
  USE ions_base, ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE basis,     ONLY : natomwfc
  USE cell_base 
  USE constants, ONLY: rytoev, eps4 
  USE gvect 
  USE klist, ONLY: xk, nks, nkstot, nelec
  USE ldaU 
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE symme, ONLY: nsym, irt, d1, d2, d3
  USE wvfct 
  USE control_flags, ONLY: gamma_only 
  USE uspp, ONLY: nkb, vkb
  USE uspp_param, ONLY: upf
  USE becmod,   ONLY: becp, rbecp, calbec
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, find_free_unit
  USE spin_orb, ONLY: lspinorb
  USE mp,       ONLY: mp_bcast
  USE mp_global,        ONLY : npool, nproc_pool, me_pool, root_pool, &
                               intra_pool_comm, init_ortho_group, me_image, &
                               ortho_comm, np_ortho, me_ortho, ortho_comm_id, &
                               leg_ortho, mpime
  USE wavefunctions_module, ONLY: evc 
  USE parallel_toolkit, ONLY : zsqmred, zsqmher, zsqmdst, zsqmcll, dsqmsym
  USE zhpev_module,     ONLY : pzhpev_drv, zhpev_drv
  USE descriptors
  USE projections
  !
  IMPLICIT NONE 
  !
  COMPLEX(DP), PARAMETER :: zero = ( 0.0d0, 0.0d0 )
  COMPLEX(DP), PARAMETER :: one  = ( 1.0d0, 0.0d0 )

  CHARACTER (len=*) :: filproj
  INTEGER :: ik, ibnd, i, j, na, nb, nt, isym, n,  m, m1, l, nwfc,& 
       nwfc1, lmax_wfc, is, iunproj, iunaux
  REAL(DP), ALLOCATABLE :: e (:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:)
  COMPLEX(DP), ALLOCATABLE :: work1(:), proj0(:,:)
  COMPLEX(DP), ALLOCATABLE :: overlap_d(:,:), work_d(:,:), diag(:,:), vv(:,:)
  COMPLEX(DP), ALLOCATABLE :: e_work_d(:,:)
  ! Some workspace for k-point calculation ... 
  REAL   (DP), ALLOCATABLE ::rwork1(:),rproj0(:,:)
  REAL   (DP), ALLOCATABLE ::roverlap_d(:,:)
  ! ... or for gamma-point. 
  REAL(DP), ALLOCATABLE :: charges(:,:,:), proj1 (:)
  REAL(DP) :: psum, totcharge(2)
  INTEGER  :: nksinit, nkslast
  CHARACTER(LEN=256) :: filename
  CHARACTER(LEN=256) :: auxname
  CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/) 
  INTEGER, ALLOCATABLE :: idx(:) 
  LOGICAL :: lsym
  INTEGER :: desc( descla_siz_ )
  INTEGER, ALLOCATABLE :: desc_ip( :, :, : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
    ! matrix distribution descriptors
  INTEGER :: nx, nrl, nrlx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  INTEGER, ALLOCATABLE :: notcnv_ip( : )
  INTEGER, ALLOCATABLE :: ic_notcnv( : )
  ! 
  ! 
  WRITE( stdout, '(/5x,"Calling pprojwave .... ")') 
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used")') 
  END IF
  !
  ! Open file as temporary storage
  !
  iunaux = find_free_unit()
  write( auxname, fmt='(I6.1)' ) mpime
  auxname = 'AUX'//TRIM( ADJUSTL( auxname ) )
  OPEN( unit=iunaux, file=TRIM(auxname), status='unknown', form='unformatted')
  ! 
  !
  ALLOCATE( ic_notcnv( np_ortho(2) ) )
  ALLOCATE( notcnv_ip( np_ortho(2) ) )
  ALLOCATE( desc_ip( descla_siz_ , np_ortho(1), np_ortho(2) ) )
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ) )
  !
  CALL desc_init( natomwfc, desc, desc_ip )
  ! 
  ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1 
  ! 
  CALL d_matrix (d1, d2, d3)   
  ! 
  ! fill structure nlmchi 
  ! 
  ALLOCATE (nlmchi(natomwfc)) 
  nwfc=0 
  lmax_wfc = 0 
  DO na = 1, nat 
     nt = ityp (na) 
     DO n = 1, upf(nt)%nwfc
        IF (upf(nt)%oc (n) >= 0.d0) THEN 
           l = upf(nt)%lchi (n) 
           lmax_wfc = MAX (lmax_wfc, l ) 
           DO m = 1, 2 * l + 1 
              nwfc=nwfc+1 
              nlmchi(nwfc)%na = na 
              nlmchi(nwfc)%n  =  n 
              nlmchi(nwfc)%l  =  l 
              nlmchi(nwfc)%m  =  m 
           ENDDO 
        ENDIF 
     ENDDO 
  ENDDO 
  ! 
  IF (lmax_wfc > 3) CALL errore ('projwave', 'l > 3 not yet implemented', 1) 
  IF (nwfc /= natomwfc) CALL errore ('projwave', 'wrong # of atomic wfcs?', 1) 
  ! 
  ! 
  IF( ionode ) THEN
     WRITE( stdout, * ) 
     WRITE( stdout, * ) ' Problem Sizes '
     WRITE( stdout, * ) ' natomwfc = ', natomwfc
     WRITE( stdout, * ) ' nbnd     = ', nbnd
     WRITE( stdout, * ) ' nkstot   = ', nkstot
     WRITE( stdout, * ) ' npwx     = ', npwx
     WRITE( stdout, * ) ' nkb      = ', nkb 
     WRITE( stdout, * ) 
  END IF
  ! 
  ALLOCATE( proj (natomwfc, nbnd, nkstot) ) 
  proj      = 0.d0 
  !
  IF (.NOT. lda_plus_u) ALLOCATE(swfcatom (npwx , natomwfc ) ) 
  ALLOCATE(wfcatom (npwx, natomwfc) ) 
  !
  ALLOCATE(e (natomwfc) ) 
  ! 
  !    loop on k points 
  ! 
  CALL init_us_1 
  CALL init_at_1 
  ! 
  DO ik = 1, nks 
     !
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin) 
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1) 
 
     CALL atomic_wfc (ik, wfcatom) 
 
     CALL init_us_2 (npw, igk, xk (1, ik), vkb) 

     IF ( gamma_only ) THEN 
        ALLOCATE (rbecp (nkb,natomwfc)) 
        CALL calbec ( npw, vkb, wfcatom, rbecp) 
     ELSE 
        ALLOCATE ( becp (nkb,natomwfc)) 
        CALL calbec ( npw, vkb, wfcatom, becp) 
     END IF 
 
     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom) 

     IF ( gamma_only ) THEN 
        DEALLOCATE (rbecp)
     ELSE
        DEALLOCATE ( becp)
     END IF 
     ! 
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i> 
     ! calculate overlap matrix O_ij = <phi_i|\hat S|\phi_j> 
     ! 
     IF( la_proc ) THEN
        ALLOCATE(overlap_d (nx, nx) ) 
     ELSE
        ALLOCATE(overlap_d (1, 1) ) 
     END IF
     overlap_d = (0.d0,0.d0) 
     IF ( gamma_only ) THEN 
        IF( la_proc ) THEN
           ALLOCATE(roverlap_d (nx, nx) ) 
        ELSE
           ALLOCATE(roverlap_d (1, 1) ) 
        END IF
        roverlap_d = 0.d0 
        CALL calbec_ddistmat( npw, wfcatom, swfcatom, natomwfc, nx, roverlap_d )
        overlap_d(:,:)=CMPLX(roverlap_d(:,:),0.d0)
        ! TEMP: diagonalization routine for real matrix should be used instead 
     ELSE 
        CALL calbec_zdistmat( npw, wfcatom, swfcatom, natomwfc, nx, overlap_d )
     END IF 
     ! 
     ! calculate O^{-1/2} 
     ! 
     IF ( desc( lambda_node_ ) > 0 ) THEN
        !
        !  Compute local dimension of the cyclically distributed matrix
        !
        ALLOCATE(work_d (nx, nx) ) 

        nrl  = desc( la_nrl_ )
        nrlx = desc( la_nrlx_ )

        ALLOCATE( diag( nrlx, natomwfc ) )
        ALLOCATE( vv( nrlx, natomwfc ) )
        !
        CALL blk2cyc_zredist( natomwfc, diag, nrlx, natomwfc, overlap_d, nx, nx, desc )
        !
        CALL pzhpev_drv( 'V', diag, nrlx, e, vv, nrlx, nrl, natomwfc, &
           desc( la_npc_ ) * desc( la_npr_ ), desc( la_me_ ), desc( la_comm_ ) )
        !
        CALL cyc2blk_zredist( natomwfc, vv, nrlx, natomwfc, work_d, nx, nx, desc )
        !
        DEALLOCATE( vv )
        DEALLOCATE( diag )
        !
     ELSE
        ALLOCATE(work_d (1, 1) ) 
     END IF

     CALL mp_bcast( e, root_pool, intra_pool_comm )

     DO i = 1, natomwfc 
        e (i) = 1.d0 / dsqrt (e (i) ) 
     ENDDO 

     IF ( desc( lambda_node_ ) > 0 ) THEN
        ALLOCATE(e_work_d (nx, nx) ) 
        DO j = 1, desc( nlac_ )
           DO i = 1, desc( nlar_ )
              e_work_d( i, j ) = e( j + desc( ilac_ ) - 1 ) * work_d( i, j )
           END DO
        END DO
        CALL sqr_zmm_cannon( 'N', 'C', natomwfc, ONE, e_work_d, nx, work_d, nx, ZERO, overlap_d, nx, desc )
        CALL zsqmher( natomwfc, overlap_d, nx, desc )   
        DEALLOCATE( e_work_d )
     END IF
     !
     DEALLOCATE( work_d )
     ! 
     ! calculate wfcatom = O^{-1/2} \hat S | phi> 
     ! 
     IF ( gamma_only ) THEN 
        ! TEMP: diagonalization routine for real matrix should be used instead 
        roverlap_d(:,:)=REAL(overlap_d(:,:),DP)
        CALL wf_times_roverlap( swfcatom, roverlap_d, wfcatom )
        DEALLOCATE( roverlap_d )
     ELSE 
        CALL wf_times_overlap( swfcatom, overlap_d, wfcatom )
        DEALLOCATE( overlap_d )
     END IF 
 
     ! 
     ! make the projection <psi_i| O^{-1/2} \hat S | phi_j> 
     ! 
     IF ( gamma_only ) THEN 
        !
        ALLOCATE( rproj0(natomwfc,nbnd), rwork1 (nbnd) ) 
        CALL calbec ( npw, wfcatom, evc, rproj0) 
        !
        WRITE( iunaux ) rproj0
        !
     ELSE 
        !
        ALLOCATE( proj0(natomwfc,nbnd), work1 (nbnd) ) 
        CALL calbec ( npw, wfcatom, evc, proj0) 
        !
        WRITE( iunaux ) proj0
        !
     END IF 
     ! 
     ! symmetrize the projections 
     ! 
     IF (lsym) THEN
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
              !
              nb = irt (isym, na) 
              DO nwfc1 =1, natomwfc 
                 IF (nlmchi(nwfc1)%na == nb             .AND. & 
                      nlmchi(nwfc1)%n == nlmchi(nwfc)%n .AND. & 
                      nlmchi(nwfc1)%l == nlmchi(nwfc)%l .AND. & 
                      nlmchi(nwfc1)%m == 1 ) go to 10 
              END DO 
              CALL errore('projwave','cannot symmetrize',1) 
10            nwfc1=nwfc1-1 
              !  
              !  nwfc1 is the first rotated atomic wfc corresponding to nwfc 
              ! 
              IF ( gamma_only ) THEN 
                 IF (l == 0) THEN 
                    rwork1(:) = rproj0 (nwfc1 + 1,:) 
                 ELSE IF (l == 1) THEN  
                    rwork1(:) = 0.d0   
                    DO m1 = 1, 3   
                       rwork1(:)=rwork1(:)+d1(m1,m,isym)*rproj0(nwfc1+m1,:) 
                    ENDDO 
                 ELSE IF (l == 2) THEN  
                    rwork1(:) = 0.d0   
                    DO m1 = 1, 5   
                       rwork1(:)=rwork1(:)+d2(m1,m,isym)*rproj0(nwfc1+m1,:) 
                    ENDDO 
                 ELSE IF (l == 3) THEN  
                    rwork1(:) = 0.d0   
                    DO m1 = 1, 7   
                       rwork1(:)=rwork1(:)+d3(m1,m,isym)*rproj0(nwfc1+m1,:) 
                    ENDDO 
                 ENDIF 
                 DO ibnd = 1, nbnd 
                    proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + & 
                         rwork1(ibnd) * rwork1(ibnd) / nsym 
                 ENDDO 
              ELSE 
                 IF (l == 0) THEN 
                    work1(:) = proj0 (nwfc1 + 1,:) 
                 ELSE IF (l == 1) THEN  
                    work1(:) = 0.d0   
                    DO m1 = 1, 3   
                       work1(:)=work1(:)+d1(m1,m,isym)*proj0(nwfc1+m1,:) 
                    ENDDO 
                 ELSE IF (l == 2) THEN  
                    work1(:) = 0.d0   
                    DO m1 = 1, 5   
                       work1(:)=work1(:)+d2(m1,m,isym)*proj0(nwfc1+m1,:) 
                    ENDDO 
                 ELSE IF (l == 3) THEN  
                    work1(:) = 0.d0   
                    DO m1 = 1, 7   
                       work1(:)=work1(:)+d3(m1,m,isym)*proj0(nwfc1+m1,:) 
                    ENDDO 
                 ENDIF 
                 DO ibnd = 1, nbnd 
                    proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + & 
                         work1(ibnd) * CONJG (work1(ibnd)) / nsym 
                 ENDDO 
              END IF 
           END DO 
        END DO 
        !
     ELSE
        !
        IF ( gamma_only ) THEN 
           DO nwfc=1,natomwfc
              DO ibnd=1,nbnd
                 proj(nwfc,ibnd,ik)=ABS(rproj0(nwfc,ibnd))**2
              END DO
           END DO
        ELSE
           DO nwfc=1,natomwfc
              DO ibnd=1,nbnd
                 proj(nwfc,ibnd,ik)=ABS(proj0(nwfc,ibnd))**2
              END DO
           END DO
        END IF
        !
     END IF
     !
     IF ( gamma_only ) THEN 
        DEALLOCATE (rwork1) 
        DEALLOCATE (rproj0) 
     ELSE 
        DEALLOCATE (work1) 
        DEALLOCATE (proj0) 
     END IF 
     !
  ENDDO 
  !
  !
  DEALLOCATE (e)
  !
  DEALLOCATE (wfcatom) 
  IF (.NOT. lda_plus_u) DEALLOCATE (swfcatom) 
  ! 
  CLOSE( unit=iunaux )
  ! 
  ! 
  !   vectors et and proj are distributed across the pools 
  !   collect data for all k-points to the first pool
  ! 
  CALL poolrecover (et,       nbnd, nkstot, nks) 
  CALL poolrecover (proj,     nbnd * natomwfc, nkstot, nks) 

  !
  !  Recover proj_aux
  !
  OPEN( unit=iunaux, file=TRIM(auxname), status='old', form='unformatted')

  ALLOCATE( proj_aux (natomwfc, nbnd, nkstot) ) 
  proj_aux  = (0.d0, 0.d0) 
  DO ik = 1, nks 
     !
     IF( gamma_only ) THEN
        ALLOCATE( rproj0( natomwfc, nbnd ) ) 
        READ( iunaux ) rproj0(:,:)
        proj_aux(:,:,ik) = CMPLX( rproj0(:,:), 0.0d0 )
        DEALLOCATE ( rproj0 ) 
     ELSE
        READ( iunaux ) proj_aux(:,:,ik)
     END IF
     !
  END DO
  !
  CALL poolrecover (proj_aux, 2 * nbnd * natomwfc, nkstot, nks)
  !
  CLOSE( unit=iunaux )
  ! 
  IF ( ionode ) THEN 
     ! 
     ! write on the file filproj
     ! 
     IF (filproj.NE.' ') THEN
        DO is=1,nspin
           IF (nspin==2) THEN
              if (is==1) filename=TRIM(filproj)//'.up'
              if (is==2) filename=TRIM(filproj)//'.down'
              nksinit=(nkstot/2)*(is-1)+1
              nkslast=(nkstot/2)*is
           ELSE
              filename=TRIM(filproj)
              nksinit=1
              nkslast=nkstot
           END IF
           iunproj=33
           CALL write_io_header(filename, iunproj, title, nrx1, nrx2, nrx3, &
                nr1, nr2, nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, &
                ecutwfc, nkstot/nspin,nbnd,natomwfc) 
           DO nwfc = 1, natomwfc 
              WRITE(iunproj,'(2i5,a3,3i5)') & 
                  nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), & 
                  nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m 
              DO ik=nksinit,nkslast
                 DO ibnd=1,nbnd
                   WRITE(iunproj,'(2i8,f20.10)') ik,ibnd, &
                                                 abs(proj(nwfc,ibnd,ik)) 
                 END DO
              END DO
           END DO
           CLOSE(iunproj)
        END DO
     END IF

     !
     ! write projections to file using iotk
     !
     CALL write_proj( "atomic_proj.xml", proj_aux )

     ! 
     ! write on the standard output file 
     ! 
     WRITE( stdout,'(/5x,"Atomic states used for projection")') 
     WRITE( stdout,'( 5x,"(read from pseudopotential files):"/)')              
     DO nwfc = 1, natomwfc 
        WRITE(stdout,1000) & 
             nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), & 
             nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m 
     END DO
1000 FORMAT (5x,"state #",i4,": atom ",i3," (",a3,"), wfc ",i2, &
                " (l=",i1," m=",i2,")") 
     ! 
     ALLOCATE(idx(natomwfc), proj1 (natomwfc) ) 
     !
     DO ik = 1, nkstot 
        WRITE( stdout, '(/" k = ",3f14.10)') (xk (i, ik) , i = 1, 3) 
        DO ibnd = 1, nbnd 
           WRITE( stdout, '(5x,"e = ",f11.5," eV")') et (ibnd, ik) * rytoev
           ! 
           ! sort projections by magnitude, in decreasing order 
           ! 
           DO nwfc = 1, natomwfc 
              idx (nwfc) = 0 
              proj1 (nwfc) = - proj (nwfc, ibnd, ik) 
           END DO
           !
           ! projections differing by less than 1.d-4 are considered equal
           !
           CALL hpsort_eps (natomwfc, proj1, idx, eps4) 
           ! 
           !  only projections that are larger than 0.001 are written 
           ! 
           DO nwfc = 1, natomwfc 
              proj1 (nwfc) = - proj1(nwfc) 
              IF ( ABS (proj1(nwfc)) < 0.001d0 ) go to 20 
           END DO
           nwfc = natomwfc + 1 
20         nwfc = nwfc -1 
           ! 
           ! fancy (?!?) formatting 
           ! 
           WRITE( stdout, '(5x,"psi = ",5(f5.3,"*[#",i4,"]+"))') & 
                (proj1 (i), idx(i), i = 1, MIN(5,nwfc)) 
           DO j = 1, (nwfc-1)/5 
              WRITE( stdout, '(10x,"+",5(f5.3,"*[#",i4,"]+"))') & 
                   (proj1 (i), idx(i), i = 5*j+1, MIN(5*(j+1),nwfc)) 
           END DO
           psum = 0.d0 
           DO nwfc = 1, natomwfc 
              psum = psum + proj (nwfc, ibnd, ik) 
           END DO
           WRITE( stdout, '(4x,"|psi|^2 = ",f5.3)') psum 
           ! 
        ENDDO
     ENDDO
     !
     DEALLOCATE (idx, proj1) 
     ! 
     ! estimate partial charges (Loewdin) on each atom 
     ! 
     ALLOCATE ( charges (nat, 0:lmax_wfc, nspin ) ) 
     charges = 0.0d0 
     DO ik = 1, nkstot 
        IF ( nspin == 1 ) THEN
           current_spin = 1
        ELSE IF ( nspin == 2 ) THEN
           current_spin = isk ( ik )
        ELSE
           CALL errore ('projwfc_nc',' called in the wrong case ',1)
        END IF
        DO ibnd = 1, nbnd 
           DO nwfc = 1, natomwfc 
              na= nlmchi(nwfc)%na 
              l = nlmchi(nwfc)%l 
              charges(na,l,current_spin) = charges(na,l,current_spin) + &
                   wg (ibnd,ik) * proj (nwfc, ibnd, ik) 
           ENDDO 
        END DO 
     END DO 
     ! 
     WRITE( stdout, '(/"Lowdin Charges: "/)') 
     ! 
     DO na = 1, nat 
        DO is = 1, nspin
           totcharge(is) = SUM(charges(na,0:lmax_wfc,is))
        END DO
        IF ( nspin == 1) THEN
           WRITE( stdout, 2000) na, totcharge(1), &
                ( l_label(l), charges(na,l,1), l= 0,lmax_wfc)
        ELSE IF ( nspin == 2) THEN 
           WRITE( stdout, 2000) na, totcharge(1) + totcharge(2), &
                ( l_label(l), charges(na,l,1) + charges(na,l,2), l=0,lmax_wfc)
           WRITE( stdout, 2001) totcharge(1), &
                ( l_label(l), charges(na,l,1), l= 0,lmax_wfc) 
           WRITE( stdout, 2002) totcharge(2), &
                ( l_label(l), charges(na,l,2), l= 0,lmax_wfc) 
           WRITE( stdout, 2003) totcharge(1) - totcharge(2), &
                ( l_label(l), charges(na,l,1) - charges(na,l,2), l=0,lmax_wfc)
        END IF
     END DO 
2000 FORMAT (5x,"Atom # ",i3,": total charge = ",f8.4,4(", ",a1," =",f8.4))
2001 FORMAT (15x,"  spin up      = ",f8.4,4(", ",a1," =",f8.4))
2002 FORMAT (15x,"  spin down    = ",f8.4,4(", ",a1," =",f8.4))
2003 FORMAT (15x,"  polarization = ",f8.4,4(", ",a1," =",f8.4))
     !
     psum = SUM(charges(:,:,:)) / nelec 
     WRITE( stdout, '(5x,"Spilling Parameter: ",f8.4)') 1.0d0 - psum 
     ! 
     ! Sanchez-Portal et al., Sol. State Commun.  95, 685 (1995). 
     ! The spilling parameter measures the ability of the basis provided by 
     ! the pseudo-atomic wfc to represent the PW eigenstates, 
     ! by measuring how much of the subspace of the Hamiltonian 
     ! eigenstates falls outside the subspace spanned by the atomic basis 
     ! 
     DEALLOCATE (charges) 
     !
  END IF
  ! 
  RETURN
  !
CONTAINS
  !
  SUBROUTINE desc_init( nsiz, desc, desc_ip )
     !
     INTEGER, INTENT(IN)  :: nsiz
     INTEGER, INTENT(OUT) :: desc(:)
     INTEGER, INTENT(OUT) :: desc_ip(:,:,:)
     INTEGER :: i, j, rank
     INTEGER :: coor_ip( 2 )
     !
     CALL descla_init( desc, nsiz, nsiz, np_ortho, me_ortho, ortho_comm, ortho_comm_id )
     !
     nx = desc( nlax_ )
     !
     DO j = 0, desc( la_npc_ ) - 1
        DO i = 0, desc( la_npr_ ) - 1
           coor_ip( 1 ) = i
           coor_ip( 2 ) = j
           CALL descla_init( desc_ip(:,i+1,j+1), desc( la_n_ ), desc( la_nx_ ), &
                             np_ortho, coor_ip, ortho_comm, 1 )
           CALL GRID2D_RANK( 'R', desc( la_npr_ ), desc( la_npc_ ), i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        END DO
     END DO
     !
     la_proc = .FALSE.
     IF( desc( lambda_node_ ) > 0 ) la_proc = .TRUE.
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
     INTEGER, INTENT(IN) :: npw ! local number of plane wave
     INTEGER, INTENT(IN) :: n   ! global dimension of matrix dm
     INTEGER, INTENT(IN) :: nx  ! local leading dimension of matrix dm
                                ! WARNING: nx is the same on all proc, SIZE( dm, 1 ) NO! 
     COMPLEX(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = zero
     !
     ldv = SIZE( v, 1 )
     ldw = SIZE( w, 1 )
     !
     DO ipc = 1, desc( la_npc_ ) !  loop on column procs
        !
        nc = desc_ip( nlac_ , 1, ipc )
        ic = desc_ip( ilac_ , 1, ipc )
        !
        DO ipr = 1, ipc ! desc( la_npr_ ) ! ipc ! use symmetry for the loop on row procs
           !
           nr = desc_ip( nlar_ , ipr, ipc )
           ir = desc_ip( ilar_ , ipr, ipc )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL ZGEMM( 'C', 'N', nr, nc, npw, ONE , &
                       v(1,ir), ldv, w(1,ic), ldw, ZERO, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, intra_pool_comm )

        END DO
        !
     END DO
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
     INTEGER, INTENT(IN) :: npw ! local number of plane wave
     INTEGER, INTENT(IN) :: n   ! global dimension of matrix dm
     INTEGER, INTENT(IN) :: nx  ! local leading dimension of matrix dm
                                ! WARNING: nx is the same on all proc, SIZE( dm, 1 ) NO! 
     REAL(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     REAL(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     npw2  = 2*npw
     npwx2 = 2*npwx
     !
     work = zero
     !
     ldv = SIZE( v, 1 )
     ldw = SIZE( w, 1 )
     !
     DO ipc = 1, desc( la_npc_ ) !  loop on column procs
        !
        nc = desc_ip( nlac_ , 1, ipc )
        ic = desc_ip( ilac_ , 1, ipc )
        !
        DO ipr = 1, ipc ! desc( la_npr_ ) ! ipc ! use symmetry for the loop on row procs
           !
           nr = desc_ip( nlar_ , ipr, ipc )
           ir = desc_ip( ilar_ , ipr, ipc )
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

        END DO
        !
     END DO
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
  SUBROUTINE wf_times_overlap( swfc, ovr, wfc )

     COMPLEX(DP) :: swfc( :, : ), ovr( :, : ), wfc( :, : )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc( la_npc_ )
        !
        nc = desc_ip( nlac_ , 1, ipc )
        ic = desc_ip( ilac_ , 1, ipc )
        !
        beta = ZERO

        DO ipr = 1, desc( la_npr_ )
           !
           nr = desc_ip( nlar_ , ipr, ipc )
           ir = desc_ip( ilar_ , ipr, ipc )
           !
           root = rank_ip( ipr, ipc )

           IF( ipr-1 == desc( la_myr_ ) .AND. ipc-1 == desc( la_myc_ ) .AND. la_proc ) THEN
              !
              !  this proc sends his block
              !
              CALL mp_bcast( ovr, root, intra_pool_comm )
              CALL ZGEMM( 'N', 'N', npw, nc, nr, ONE, &
                          swfc(1,ir), npwx, ovr, nx, beta, wfc(1,ic), npwx )
           ELSE
              !
              !  all other procs receive
              !
              CALL mp_bcast( vtmp, root, intra_pool_comm )
              CALL ZGEMM( 'N', 'N', npw, nc, nr, ONE, &
                       swfc(1,ir), npwx, vtmp, nx, beta, wfc(1,ic), npwx )
           END IF
           !
           beta = ONE

        END DO
        !
     END DO
     !
     DEALLOCATE( vtmp )

     RETURN

  END SUBROUTINE wf_times_overlap

  !
  SUBROUTINE wf_times_roverlap( swfc, ovr, wfc )

     USE gvect, ONLY : gstart

     COMPLEX(DP) :: swfc( :, : ), wfc( :, : )
     REAL(DP)    :: ovr( :, : )
     !
     INTEGER :: ipc, ipr, npw2, npwx2
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     REAL(DP) :: beta

     npw2  = 2*npw
     npwx2 = 2*npwx

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc( la_npc_ )
        !
        nc = desc_ip( nlac_ , 1, ipc )
        ic = desc_ip( ilac_ , 1, ipc )
        !
        beta = 0.0d0

        DO ipr = 1, desc( la_npr_ )
           !
           nr = desc_ip( nlar_ , ipr, ipc )
           ir = desc_ip( ilar_ , ipr, ipc )
           !
           root = rank_ip( ipr, ipc )

           IF( ipr-1 == desc( la_myr_ ) .AND. ipc-1 == desc( la_myc_ ) .AND. la_proc ) THEN
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
           END IF
           !
           beta = 1.0d0

        END DO
        !
     END DO
     !
     DEALLOCATE( vtmp )

     RETURN

  END SUBROUTINE wf_times_roverlap
  !
END SUBROUTINE pprojwave
