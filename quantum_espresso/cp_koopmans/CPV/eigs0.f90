!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
   subroutine eigs0( ei, tprint, nspin, nupdwn, iupdwn, lf, f, nx, lambda, nudx, desc )
!-----------------------------------------------------------------------
!     computes eigenvalues (wr) of the real symmetric matrix lambda
!     Note that lambda as calculated is multiplied by occupation numbers
!     so empty states yield zero. Eigenvalues are printed out in eV
!
      use kinds,             only : DP
      use io_global,         only : stdout
      use constants,         only : autoev
      use dspev_module,      only : dspev_drv, pdspev_drv
      USE sic_module,        only : self_interaction
      USE cp_main_variables, only : nlam, la_proc
      USE descriptors,       ONLY : nlar_ , nlac_ , ilar_ , ilac_ , lambda_node_ , la_me_ , la_n_ , &
                                    descla_siz_ , la_npr_ , la_npc_ , la_nrl_ , la_nrlx_ , la_comm_ , &
                                    nlax_ , la_myc_ , la_myr_
      USE mp,                only : mp_sum, mp_bcast
      USE mp_global,         only : intra_image_comm, root_image
      use nksic,             only : f_cutoff

      implicit none
! input
      logical, intent(in) :: tprint, lf
      integer, intent(in) :: nspin, nx, nudx, nupdwn(nspin), iupdwn(nspin)
      integer, intent(in) :: desc( descla_siz_ , 2 )
      real(DP), intent(in) :: lambda( nlam, nlam, nspin ), f( nx )
      real(DP), intent(out) :: ei( nudx, nspin )
! local variables
      real(DP), allocatable :: ap(:), wr(:)
      real(DP) zr(1)
      integer :: iss, j, i, k, n, nspin_eig, npaired
      INTEGER :: ir, nr, np
      logical :: tsic
      CHARACTER(LEN=80) :: msg
!
      tsic = ( ABS( self_interaction) /= 0 )

      IF( tsic ) THEN
         nspin_eig = 1
         npaired   = nupdwn(2)
      ELSE
         nspin_eig = nspin
         npaired   = 0
      END IF


      do iss = 1, nspin_eig

         IF( nudx < nupdwn(iss) ) THEN 
            WRITE( msg, 100 ) nudx, SIZE( ei, 1 ), nupdwn(iss)
100         FORMAT( ' wrong dimension array ei = ', 3I10 )
            CALL errore( ' eigs0 ', msg, 1 )
         END IF

         IF( tsic ) THEN
            n = npaired
         ELSE
            n = nupdwn(iss)
         END IF

         allocate( wr( n ) )

         IF( la_proc ) THEN

            np = desc( la_npc_ , iss ) * desc( la_npr_ , iss )

            IF( np > 1 ) THEN

               !  matrix is distributed

               CALL qe_pdsyevd( .false., n, desc(1,iss), lambda(1,1,iss), SIZE(lambda,1), wr )

            ELSE

               !  matrix is not distributed

               allocate( ap( n * ( n + 1 ) / 2 ) )

               k = 0
               do i = 1, n
                  do j = i, n
                     k = k + 1
                     ap( k ) = lambda( j, i, iss )
                  end do
               end do

               CALL dspev_drv( 'N', 'L', n, ap, wr, zr, 1 )

               deallocate( ap )

            END IF

         END IF

         call mp_bcast( wr, root_image, intra_image_comm )

         if( lf ) then
            do i = 1, n
              wr(i)=wr(i)/max(f(iupdwn(iss)-1+i),f_cutoff)
            end do
         end if
         !
         !     store eigenvalues
         !
         ei( 1:n, iss ) = wr( 1:n )

         IF( tsic ) THEN
            !
            !  store unpaired state
            !
            ei( 1:n,       1 ) = ei( 1:n, 1 ) / 2.0d0
            ei( nupdwn(1), 1 ) = 0.0d0
            if( la_proc ) then
               IF( desc( la_myc_ , iss ) == desc( la_myr_ , iss ) ) THEN
                  ir = desc( ilar_ , iss )
                  nr = desc( nlar_ , iss )
                  IF( nupdwn(1) >= ir .AND. nupdwn(1) < ir + nr ) then
                     ei( nupdwn(1), 1 ) = lambda( nupdwn(1)-ir+1, nupdwn(1)-ir+1, 1 )
                  end if
               END IF
            endif
            call mp_sum( ei( nupdwn(1), 1 ), intra_image_comm )
            !
         END IF

         ! WRITE( stdout,*)  '---- DEBUG ----' ! debug
         ! WRITE( stdout,14) ( wr( i ) * autoev / 2.0d0, i = 1, nupdwn(iss) ) ! debug

         deallocate( wr )

      end do
      !
      !
      do iss = 1, nspin

         IF( tsic .AND. iss == 2 ) THEN
            ei( 1:npaired, 2 ) = ei( 1:npaired, 1 )
         END IF

         IF( tprint ) THEN
            !
            !     print out eigenvalues
            !
            WRITE( stdout,12) 0.d0, 0.d0, 0.d0
            WRITE( stdout,14) ( ei( i, iss ) * autoev, i = 1, nupdwn(iss) )

         ENDIF

      end do

      IF( tprint ) WRITE( stdout,*)

   12 format(//' eigenvalues at k-point: ',3f6.3)
   14 format(10f8.2)
!
      return
   end subroutine eigs0

!-----------------------------------------------------------------------
   subroutine eigs0_twin( ei, tprint, nspin, nupdwn, iupdwn, lf, f, nx, lambda, nudx, desc )
!-----------------------------------------------------------------------
!     computes eigenvalues (wr) of the real symmetric matrix lambda
!     Note that lambda as calculated is multiplied by occupation numbers
!     so empty states yield zero. Eigenvalues are printed out in eV
!
      use kinds,             only : DP
      use io_global,         only : stdout
      use constants,         only : autoev
      use dspev_module,      only : dspev_drv, pdspev_drv
      use zhpev_module,      only : zhpev_drv
      USE sic_module,        only : self_interaction
      USE cp_main_variables, only : la_proc
      USE descriptors,       ONLY : nlar_ , nlac_ , ilar_ , ilac_ , lambda_node_ , la_me_ , la_n_ , &
                                    descla_siz_ , la_npr_ , la_npc_ , la_nrl_ , la_nrlx_ , la_comm_ , &
                                    nlax_ , la_myc_ , la_myr_
      USE mp,                only : mp_sum, mp_bcast
      USE mp_global,         only : intra_image_comm, root_image
      use nksic,             only : f_cutoff
      USE twin_types

      implicit none
! input
      logical, intent(in) :: tprint, lf
      integer, intent(in) :: nspin, nx, nudx, nupdwn(nspin), iupdwn(nspin)
      integer, intent(in) :: desc( descla_siz_ , 2 )
      real(DP), intent(in) :: f( nx )
      type(twin_matrix), dimension(nspin), intent(in) :: lambda
      real(DP), intent(out) :: ei( nudx, nspin )
! local variables
      real(DP), allocatable :: ap(:), wr(:)
      complex(DP), allocatable :: ap_c(:)
      real(DP) zr(1)
      complex(DP) zr_c(1)
      integer :: iss, j, i, k, n, nspin_eig, npaired
      INTEGER :: ir, nr, np
      logical :: tsic
      CHARACTER(LEN=80) :: msg
      !
      tsic = ( ABS( self_interaction) /= 0 )
      ! 
      IF ( tsic ) THEN
         nspin_eig = 1
         npaired   = nupdwn(2)
      ELSE
         nspin_eig = nspin
         npaired   = 0
      ENDIF
      !
      do iss = 1, nspin_eig
         !
         IF ( nudx < nupdwn(iss) ) THEN 
            ! 
            WRITE( msg, 100 ) nudx, SIZE( ei, 1 ), nupdwn(iss)
100         FORMAT( ' wrong dimension array ei = ', 3I10 )
            CALL errore( ' eigs0 ', msg, 1 )
            !
         ENDIF
         !
         IF( tsic ) THEN
            n = npaired
         ELSE
            n = nupdwn(iss)
         END IF
         !
         IF (n.gt.0) THEN
            ! 
            allocate( wr( n ) )
            !  
            IF ( la_proc ) THEN
               !  
               np = desc( la_npc_ , iss ) * desc( la_npr_ , iss )
               !
               IF ( np > 1 ) THEN
                  ! 
                  ! matrix is distributed
                  !  
                  IF (.not.lambda(iss)%iscmplx) THEN
                     !
                     write(6,*) "sizlambda", size(lambda(iss)%rvec,1), lambda(iss)%xdim
                     CALL qe_pdsyevd( .false., n, desc(1,iss), lambda(iss)%rvec(1,1), SIZE(lambda(iss)%rvec,1), wr )
                     !
                  ELSE
                     !
                     write(6,*) "sizlambda", size(lambda(iss)%cvec,1), lambda(iss)%xdim
                     CALL qe_pzheevd( .false., n, desc(1,iss), lambda(iss)%cvec(1,1), SIZE(lambda(iss)%cvec,1), wr )
                     !
                  ENDIF
                  !   
               ELSE
                  !
                  !!  matrix is not distributed
                  !    
                  IF (.not.lambda(1)%iscmplx) THEN
                     ! 
		     allocate( ap( n * ( n + 1 ) / 2 ) )
                     !
		     k = 0
		     do i = 1, n
		        do j = i, n
		           k = k + 1
			   ap( k ) = lambda(iss)%rvec( j, i)
		        enddo
		     enddo
                     !   
		     CALL dspev_drv( 'N', 'L', n, ap, wr, zr, 1 )
                     !  
		     deallocate( ap )
                     !
                  ELSE
                     !
		     allocate( ap_c( n * ( n + 1 ) / 2 ) )
                     !
		     k = 0
		     do i = 1, n
		        do j = i, n
			   k = k + 1
			   ap_c( k ) = lambda(iss)%cvec( j, i)
		        enddo
		     enddo
                     !
		     CALL zhpev_drv( 'N', 'L', n, ap_c, wr, zr_c, 1 )
                     !   
		     deallocate( ap_c )
                     !  
                  ENDIF
                  ! 
               ENDIF
               ! 
            ENDIF
            !
            call mp_bcast( wr, root_image, intra_image_comm )
            !  
            if ( lf ) then
               do i = 1, n
                  wr(i)=wr(i)/max(f(iupdwn(iss)-1+i),f_cutoff)
               end do
            endif
            !
            !  store eigenvalues
            !
            ei( 1:n, iss ) = wr( 1:n )
            ! 
            IF ( tsic ) THEN
               !
               !  store unpaired state
               !
               ei( 1:n,       1 ) = ei( 1:n, 1 ) / 2.0d0
               ei( nupdwn(1), 1 ) = 0.0d0
               ! 
               IF ( la_proc ) THEN
                  ! 
                  IF ( desc( la_myc_ , iss ) == desc( la_myr_ , iss ) ) THEN
                     !
                     ir = desc( ilar_ , iss )
                     nr = desc( nlar_ , iss )
                     !
                     IF ( nupdwn(1) >= ir .AND. nupdwn(1) < ir + nr ) then
                        !
                        IF (.not.lambda(1)%iscmplx) THEN
                           !
                           ei( nupdwn(1), 1 ) = lambda(1)%rvec( nupdwn(1)-ir+1, nupdwn(1)-ir+1)
                           !
                        ELSE
                           ! 
                           ei( nupdwn(1), 1 ) = DBLE(lambda(1)%cvec( nupdwn(1)-ir+1, nupdwn(1)-ir+1))
                           !  
                        ENDIF
                        ! 
                     ENDIF
                     !
                  ENDIF
                  !
               ENDIF
               ! 
               call mp_sum( ei( nupdwn(1), 1 ), intra_image_comm )
               !
            ENDIF
            !  
            ! WRITE( stdout,*)  '---- DEBUG ----' ! debug
            ! WRITE( stdout,14) ( wr( i ) * autoev / 2.0d0, i = 1, nupdwn(iss) ) ! debug
            ! 
            deallocate( wr )
            ! 
         ELSE
            ! 
            ei( 1:n, iss ) = 0.d0
            !
         ENDIF
         !
      ENDDO
      !
      DO iss = 1, nspin
         ! 
         IF ( tsic .AND. iss == 2 ) THEN
            !
            ei( 1:npaired, 2 ) = ei( 1:npaired, 1 )
            !
         ENDIF
         !
         IF ( tprint ) THEN
            !
            !     print out eigenvalues
            !
            WRITE( stdout,12) 0.d0, 0.d0, 0.d0
            WRITE( stdout,14) ( ei( i, iss ) * autoev, i = 1, nupdwn(iss) )
            !
         ENDIF
         !
      ENDDO
      ! 
      IF( tprint ) WRITE( stdout,*)
      !  
   12 format(//' eigenvalues at k-point: ',3f6.3)
   14 format(10f8.2)
      ! 
      return
      !
   end subroutine eigs0_twin

!-----------------------------------------------------------------------
   subroutine eigs0_twin_non_ortho( ei, tprint, nspin, nupdwn, iupdwn, lf, f, nx, lambda, nudx, desc )
!-----------------------------------------------------------------------
!     computes eigenvalues (wr) of the real symmetric matrix lambda
!     Note that lambda as calculated is multiplied by occupation numbers
!     so empty states yield zero. Eigenvalues are printed out in eV
!
      use kinds,             only : DP
      use io_global,         only : stdout
      use constants,         only : autoev
      use dspev_module,      only : dspev_drv, pdspev_drv, dgeev_drv
      use zhpev_module,      only : zhpev_drv, zgeev_drv
      USE sic_module,        only : self_interaction
      USE cp_main_variables, only : la_proc
      USE descriptors,       ONLY : nlar_ , nlac_ , ilar_ , ilac_ , lambda_node_ , la_me_ , la_n_ , &
                                    descla_siz_ , la_npr_ , la_npc_ , la_nrl_ , la_nrlx_ , la_comm_ , &
                                    nlax_ , la_myc_ , la_myr_
      USE mp,                only : mp_sum, mp_bcast
      USE mp_global,         only : intra_image_comm, root_image
      use nksic,             only : f_cutoff
      USE twin_types

      implicit none
! input
      logical, intent(in) :: tprint, lf
      integer, intent(in) :: nspin, nx, nudx, nupdwn(nspin), iupdwn(nspin)
      integer, intent(in) :: desc( descla_siz_ , 2 )
      real(DP), intent(in) :: f( nx )
      type(twin_matrix), dimension(nspin), intent(in) :: lambda
      real(DP), intent(out) :: ei( nudx, nspin )
! local variables
      real(DP), allocatable :: ap(:,:), wr(:), wi(:)
      complex(DP), allocatable :: ap_c(:,:)
      real(DP) zr(1,1)
      complex(DP) zr_c(1,1)
      integer :: iss, j, i, n, nspin_eig, npaired
      INTEGER :: ir, nr, np
      logical :: tsic
      CHARACTER(LEN=80) :: msg
      !
      tsic = ( ABS( self_interaction) /= 0 )
      ! 
      IF ( tsic ) THEN
         nspin_eig = 1
         npaired   = nupdwn(2)
      ELSE
         nspin_eig = nspin
         npaired   = 0
      ENDIF
      !
      do iss = 1, nspin_eig
         !
         IF( nudx < nupdwn(iss) ) THEN 
            WRITE( msg, 100 ) nudx, SIZE( ei, 1 ), nupdwn(iss)
100         FORMAT( ' wrong dimension array ei = ', 3I10 )
            CALL errore( ' eigs0 ', msg, 1 )
         END IF
         !
         IF( tsic ) THEN
            n = npaired
         ELSE
            n = nupdwn(iss)
         END IF
         !
         IF (n.gt.0) THEN
            !  
            allocate( wr( n ), wi( n ) )
            !
            IF ( la_proc ) THEN
               ! 
               np = desc( la_npc_ , iss ) * desc( la_npr_ , iss )
               !
               IF ( np > 1 ) THEN
                  !
                  !  matrix is distributed
                  !
                  IF (.not.lambda(iss)%iscmplx) THEN
                     !
                     write(6,*) "sizlambda", size(lambda(iss)%rvec,1), lambda(iss)%xdim
                     CALL qe_pdsyevd( .false., n, desc(1,iss), lambda(iss)%rvec(1,1), SIZE(lambda(iss)%rvec,1), wr )
                     ! 
                  ELSE
                     ! 
                     write(6,*) "sizlambda", size(lambda(iss)%cvec,1), lambda(iss)%xdim
                     CALL qe_pzheevd( .false., n, desc(1,iss), lambda(iss)%cvec(1,1), SIZE(lambda(iss)%cvec,1), wr )
                     !
                  ENDIF
                  !
               ELSE
                  !
                  !!  matrix is not distributed
                  !
                 IF (.not.lambda(1)%iscmplx) THEN
                    !  
                    allocate( ap( n, n ) )
                    !      
                    do i = 1, n
                       do j = 1, n
                         ap( j, i ) = lambda(iss)%rvec( j, i)
                       end do
                    enddo
                    !
                    CALL dgeev_drv( 'N', 'N', n, ap, n, wr, wi, zr, 1, zr, 1)
                    ! 
                    deallocate( ap )
                    !
                 ELSE
                    !
                    allocate( ap_c( n , n ))
                    do i = 1, n
                      do j = 1, n
                        ap_c( j, i ) = lambda(iss)%cvec( j, i)
                      end do
                    enddo
                    ! 
                    CALL zgeev_drv( 'N', 'N', n, ap_c, n, wr, wi, zr_c, 1, zr_c, 1)
                    ! 
                    deallocate( ap_c )
                    !
                 ENDIF
                 !
               ENDIF
               !
            ENDIF
            !
            call mp_bcast( wr, root_image, intra_image_comm )
            !
            !
            if ( lf ) then
               do i = 1, n
                  wr(i)=wr(i)/max(f(iupdwn(iss)-1+i),f_cutoff)
               enddo
            endif
            !
            ! store eigenvalues
            !
            ei( 1:n, iss ) = wr( 1:n )
            !
            IF ( tsic ) THEN
               !
               ! store unpaired state
               !
               ei( 1:n, 1 ) = ei( 1:n, 1 ) / 2.0d0
               ei( nupdwn(1), 1 ) = 0.0d0
               if ( la_proc ) then
                  IF ( desc( la_myc_ , iss ) == desc( la_myr_ , iss ) ) THEN
                     ir = desc( ilar_ , iss )
                     nr = desc( nlar_ , iss )
                     IF ( nupdwn(1) >= ir .AND. nupdwn(1) < ir + nr ) then
                        IF (.not.lambda(1)%iscmplx) THEN
                           ei( nupdwn(1), 1 ) = lambda(1)%rvec( nupdwn(1)-ir+1, nupdwn(1)-ir+1)
                        ELSE
                           ei( nupdwn(1), 1 ) = DBLE(lambda(1)%cvec( nupdwn(1)-ir+1, nupdwn(1)-ir+1))
                        ENDIF
                    endif
                  ENDIF
               endif
               call mp_sum( ei( nupdwn(1), 1 ), intra_image_comm )
               !
            ENDIF
            !
            deallocate( wr )
            deallocate( wi )
            ! 
         ELSE
            !
            ei( 1:n, iss ) = 0.d0
            !
         ENDIF
         !  
      enddo
      !
      do iss = 1, nspin
         !
         IF( tsic .AND. iss == 2 ) THEN
            ei( 1:npaired, 2 ) = ei( 1:npaired, 1 )
         END IF
         !
         IF( tprint ) THEN
            !
            !     print out eigenvalues
            !
            WRITE( stdout,12) 0.d0, 0.d0, 0.d0
            WRITE( stdout,14) ( ei( i, iss ) * autoev, i = 1, nupdwn(iss) )
            !
         ENDIF
         !
      enddo
      !
      IF( tprint ) WRITE( stdout,*)
   12 format(//' eigenvalues at k-point: ',3f6.3)
   14 format(10f8.2)
      !
      return
      !  
   end subroutine eigs0_twin_non_ortho

!-----------------------------------------------------------------------
   SUBROUTINE rpackgam_x( gam, f, aux )
!-----------------------------------------------------------------------
      USE kinds,            ONLY: DP
      USE mp_global, ONLY: me_image, nproc_image, intra_image_comm
      USE mp, ONLY: mp_sum
      IMPLICIT NONE
          REAL(DP), INTENT(INOUT)  :: gam(:,:)
          REAL(DP), INTENT(OUT), OPTIONAL :: aux(:)
          REAL(DP), INTENT(IN)  :: f(:)
          INTEGER n, nrl, i, j, k, jl
          nrl = SIZE(gam, 1)
          n   = SIZE(gam, 2)
          IF( PRESENT( aux ) ) THEN
            aux = 0.0d0
            IF( me_image < n ) THEN
              DO i = 1, n
                j = me_image + 1
                DO jl = 1, nrl
                  IF( j >= i ) THEN
                    !   maps (j,i) index to low-tri packed (k) index
                    k = (i-1)*n + j - i*(i-1)/2  
                    aux(k) = gam(jl,i) / f(j)
                  END IF
                  j = j + nproc_image
                END DO
              END DO
            END IF
            CALL mp_sum(aux, intra_image_comm)
          ELSE
            IF( me_image < n ) THEN
              DO i = 1, n
                j = me_image + 1
                DO jl = 1, nrl
                  gam(jl,i) = gam(jl,i) / f(j)
                  j = j + nproc_image
                END DO
              END DO
            END IF
          END IF
      RETURN
   END SUBROUTINE rpackgam_x



!-----------------------------------------------------------------------
   SUBROUTINE fermi_energy_x(eig, occ, wke, ef, qtot, temp, sume)
!-----------------------------------------------------------------------

!  this routine computes Fermi energy and weights of occupied states
!  using an improved Gaussian-smearing method
!  refs: C.L.Fu and K.M.Ho, Phys.Rev. B28, 5480 (1983)
!        M.Methfessel and A.T.Paxton Phys.Rev. B40 (15 aug. 89).
!
!  taken from APW code by J. Soler and A. Williams (jk+ss)
!  added computation of occupation numbers without k-point weight

      USE kinds,          ONLY: DP
      USE io_global,      ONLY: stdout
      USE electrons_base, ONLY: nspin, iupdwn

      IMPLICIT NONE

! ... declare subroutine arguments
      REAL(DP) :: occ(:)
      REAL(DP) ef, qtot, temp, sume
      REAL(DP) eig(:,:), wke(:,:)
      REAL(DP), PARAMETER  :: tol = 1.d-10
      INTEGER,   PARAMETER  :: nitmax = 100
      INTEGER ne, nk

! ... declare functions
      REAL(DP) stepf

! ... declare other variables
      REAL(DP) sumq,emin,emax,fac,t,drange
      INTEGER ik,ispin,ie,iter

!  end of declarations
!  ----------------------------------------------

      nk    = 1
      ne    = SIZE( occ, 1)
      sumq=0.d0
      sume=0.d0
      emin=eig(1,1)
      emax=eig(1,1)
      fac=2.d0
      IF (nspin.EQ.2) fac=1.d0

      DO ik=1,nk
        DO ispin=1,nspin
          DO ie=1,ne
            wke(ie,ispin) = fac
            occ(ie+iupdwn(ispin)-1) = fac
            sumq=sumq+wke(ie,ispin)
            sume=sume+wke(ie,ispin)*eig(ie,ispin)
            emin=MIN(emin,eig(ie,ispin))
            emax=MAX(emax,eig(ie,ispin))
          END DO
        END DO
      END DO
      ef=emax
      IF (abs(sumq-qtot).LT.tol) RETURN
      IF (sumq.LT.qtot) THEN
        WRITE( stdout,*) 'FERMIE: NOT ENOUGH STATES'
        WRITE( stdout,*) 'FERMIE: QTOT,SUMQ=',qtot,sumq
        STOP
      END IF
      t = MAX(temp,1.d-6)
      drange = t * SQRT( - LOG( tol*.01d0) )
      emin = emin - drange
      emax = emax + drange
      DO iter = 1, nitmax
        ef   = 0.5d0 * (emin+emax)
        sumq = 0.d0
        sume = 0.d0
        DO ik = 1, nk
          DO ispin = 1, nspin
            DO ie = 1, ne
              wke(ie,ispin) = fac / 2.d0 * stepf((eig(ie,ispin)-ef)/t)
              occ(ie+iupdwn(ispin)-1) = fac / 2.d0 * stepf((eig(ie,ispin)-ef)/t)
              sumq = sumq + wke(ie,ispin)
              sume = sume + wke(ie,ispin) * eig(ie,ispin)
            END DO
          END DO
        END DO
        IF (ABS(sumq-qtot).LT.tol) RETURN
        IF (sumq.LE.qtot) emin=ef
        IF (sumq.GE.qtot) emax=ef
      END DO

      WRITE( stdout,*) 'FERMIE: ITERATION HAS NOT CONVERGED.'
      WRITE( stdout,*) 'FERMIE: QTOT,SUMQ=',qtot,sumq
      STOP

   END SUBROUTINE fermi_energy_x

!
!
!

!-----------------------------------------------------------------------
   SUBROUTINE cp_eigs_real_x( nfi, lambdap, lambda )
!-----------------------------------------------------------------------

      use kinds,             only: DP
      use ensemble_dft,      only: tens, tsmear
      use electrons_base,    only: nx => nbspx, f, nspin
      use electrons_base,    only: iupdwn, nupdwn, nudx
      use electrons_module,  only: ei
      use cp_main_variables, only: descla

      IMPLICIT NONE

      INTEGER  :: nfi
      REAL(DP) :: lambda( :, :, : ), lambdap( :, :, : )
      !
      REAL(DP), ALLOCATABLE :: faux(:)

      ALLOCATE( faux(nx) )
      faux(:) = f(:) * DBLE(nspin) / 2.0d0

      if ( tens ) then 
          !
          call eigs0( ei, .false. , nspin, nupdwn, iupdwn, .false. , faux, nx, lambdap, nudx, descla )
          !
      else if ( tsmear ) then
          !
          call eigs0( ei, .false. , nspin, nupdwn, iupdwn, .false. , faux, nx, lambda, nudx, descla )
          !
      else
          !
          call eigs0( ei, .false. , nspin, nupdwn, iupdwn, .true. , faux, nx, lambda, nudx, descla )
          !
      endif
      !
      DEALLOCATE( faux )
      !
      RETURN
   END SUBROUTINE cp_eigs_real_x

!-----------------------------------------------------------------------
   SUBROUTINE cp_eigs_twin_x( nfi, lambdap, lambda )
!-----------------------------------------------------------------------

      use kinds,             only: DP
      use ensemble_dft,      only: tens, tsmear
      use electrons_base,    only: nx => nbspx, f, nspin
      use electrons_base,    only: iupdwn, nupdwn, nudx
      use electrons_module,  only: ei
      use cp_main_variables, only: descla
      use twin_types

      IMPLICIT NONE

      INTEGER  :: nfi
      type(twin_matrix), dimension(nspin) :: lambda, lambdap
      !
      REAL(DP), ALLOCATABLE :: faux(:)

      ALLOCATE( faux(nx) )
      faux(:) = f(:) * DBLE(nspin) / 2.0d0

      if ( tens ) then 
          !
          call eigs0_twin( ei, .false. , nspin, nupdwn, iupdwn, .false. , faux, nx, lambdap, nudx, descla )
          !
      else if ( tsmear ) then
          !
          call eigs0_twin( ei, .false. , nspin, nupdwn, iupdwn, .false. , faux, nx, lambda, nudx, descla )
          !
      else
          !
          call eigs0_twin( ei, .false. , nspin, nupdwn, iupdwn, .true. , faux, nx, lambda, nudx, descla )
          !
      endif
      !
      DEALLOCATE( faux )
      !
      RETURN
   END SUBROUTINE cp_eigs_twin_x

!-----------------------------------------------------------------------
   SUBROUTINE cp_eigs_twin_non_ortho_x( nfi, lambdap, lambda )
!-----------------------------------------------------------------------

      use kinds,             only: DP
      use ensemble_dft,      only: tens, tsmear
      use electrons_base,    only: nx => nbspx, f, nspin
      use electrons_base,    only: iupdwn, nupdwn, nudx
      use electrons_module,  only: ei
      use cp_main_variables, only: descla
      use twin_types

      IMPLICIT NONE

      INTEGER  :: nfi
      type(twin_matrix), dimension(nspin) :: lambda, lambdap
      !
      REAL(DP), ALLOCATABLE :: faux(:)

      ALLOCATE( faux(nx) )
      faux(:) = f(:) * DBLE(nspin) / 2.0d0

      if ( tens ) then 
          !
          call eigs0_twin_non_ortho( ei, .false. , nspin, nupdwn, iupdwn, .false. , faux, nx, lambdap, nudx, descla )
          !
      else if ( tsmear ) then
          !
          call eigs0_twin_non_ortho( ei, .false. , nspin, nupdwn, iupdwn, .false. , faux, nx, lambda, nudx, descla )
          !
      else
          !
          call eigs0_twin_non_ortho( ei, .false. , nspin, nupdwn, iupdwn, .true. , faux, nx, lambda, nudx, descla )
          !
      endif
      !
      DEALLOCATE( faux )
      !
      RETURN
   END SUBROUTINE cp_eigs_twin_non_ortho_x
