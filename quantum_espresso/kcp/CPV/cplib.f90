!
! Copyright (C) 2002-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE atomic_wfc(eigr, n_atomic_wfc, wfc)
!-----------------------------------------------------------------------
!
! Compute atomic wavefunctions in G-space
!
   USE kinds, ONLY: DP
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: g, gx
   USE ions_base, ONLY: nsp, na, nat
   USE cell_base, ONLY: tpiba
   USE atom, ONLY: rgrid
   USE uspp_param, ONLY: upf
!
   IMPLICIT NONE
   INTEGER, INTENT(in) :: n_atomic_wfc
   COMPLEX(DP), INTENT(in) :: eigr(ngw, nat)
   COMPLEX(DP), INTENT(inout):: wfc(ngw, n_atomic_wfc)
!
   INTEGER :: natwfc, ndm, is, ia, ir, nb, l, m, lm, i, lmax_wfc, isa, ig
   REAL(DP), ALLOCATABLE ::  ylm(:, :), q(:), jl(:), vchi(:), chiq(:)
!
! calculate max angular momentum required in wavefunctions
!
   IF (.NOT. ALLOCATED(rgrid)) &
      CALL errore(' atomic_wfc ', ' rgrid not allocated ', 1)

   lmax_wfc = -1
   DO is = 1, nsp
      lmax_wfc = MAX(lmax_wfc, MAXVAL(upf(is)%lchi(1:upf(is)%nwfc)))
   END DO
   !
   ALLOCATE (ylm(ngw, (lmax_wfc + 1)**2))
   !
   CALL ylmr2((lmax_wfc + 1)**2, ngw, gx, g, ylm)
   ndm = MAXVAL(rgrid(1:nsp)%mesh)
   !
   ALLOCATE (jl(ndm), vchi(ndm))
   ALLOCATE (q(ngw), chiq(ngw))
!
   DO i = 1, ngw
      q(i) = SQRT(g(i))*tpiba
   END DO
!
   !DEALLOCATE(ylm)
   !DEALLOCATE(vchi, jl)
   !write(6,*) ngw,ndm,"deallocating", ubound(q), ubound(chiq), ubound(jl), ubound(ylm), ubound(vchi)
   !DEALLOCATE(chiq)
   !write(6,*) ngw,ndm,"deallocating", ubound(q), ubound(chiq), ubound(jl), ubound(ylm), ubound(vchi)
   !DEALLOCATE(q)
   !stop

   natwfc = 0
   isa = 0
   DO is = 1, nsp
      !
      !   radial fourier transform of the chi functions
      !   NOTA BENE: chi is r times the radial part of the atomic wavefunction
      !
      DO nb = 1, upf(is)%nwfc
         l = upf(is)%lchi(nb)
         DO i = 1, ngw
            CALL sph_bes(rgrid(is)%mesh, rgrid(is)%r, q(i), l, jl)
            DO ir = 1, rgrid(is)%mesh
               vchi(ir) = upf(is)%chi(ir, nb)*rgrid(is)%r(ir)*jl(ir)
            END DO
            CALL simpson_cp90(rgrid(is)%mesh, vchi, rgrid(is)%rab, chiq(i))
         END DO
         !
         !   multiply by angular part and structure factor
         !   NOTA BENE: the factor i^l MUST be present!!!
         !
         DO m = 1, 2*l + 1
            lm = l**2 + m
            DO ia = 1 + isa, na(is) + isa
               natwfc = natwfc + 1
               do ig = 1, ngw
                  wfc(ig, natwfc) = CMPLX(0.d0, 1.d0)**l*eigr(ig, ia)*ylm(ig, lm)*chiq(ig)
               end do
            END DO
         END DO
      END DO
      isa = isa + na(is)
   END DO
!
   IF (natwfc .NE. n_atomic_wfc)                                       &
  &     CALL errore('atomic_wfc', 'unexpected error', natwfc)
!
   DEALLOCATE (ylm)
   DEALLOCATE (q, chiq)
   DEALLOCATE (vchi, jl)
!
   RETURN
END SUBROUTINE atomic_wfc
!

!-----------------------------------------------------------------------
FUNCTION n_atom_wfc_x()
!----------------------------------------------------------------------------
   !
   ! ... Find max number of bands needed
   !
   USE ions_base, ONLY: na, nsp
   USE kinds, ONLY: DP
   USE uspp_param, ONLY: upf
   !
   IMPLICIT NONE
   !
   INTEGER  :: n_atom_wfc_x
   INTEGER  :: is, n
   !
   n_atom_wfc_x = 0
   !
   DO is = 1, nsp
      !
      DO n = 1, upf(is)%nwfc
         !
         IF (upf(is)%oc(n) >= 0.D0) THEN
            !
            n_atom_wfc_x = n_atom_wfc_x + na(is)*(2*upf(is)%lchi(n) + 1)
            !
         END IF
         !
      END DO
      !
   END DO
   !
   RETURN
END FUNCTION

!

!-----------------------------------------------------------------------
FUNCTION cscnorm(bec, nkbx, cp, ngwx, i, n, lgam)
!-----------------------------------------------------------------------
!     requires in input the updated bec(i)
!
   USE ions_base, ONLY: na
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: gstart
   USE cvan, ONLY: ish, nvb
   USE uspp_param, ONLY: nh
   USE uspp, ONLY: qq
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE kinds, ONLY: DP
   USE twin_types !added:giovanni
!
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: i, n
   INTEGER, INTENT(IN) :: ngwx, nkbx
   type(twin_matrix)   :: bec!modified:giovanni
   COMPLEX(DP) :: cp(ngwx, n)
   LOGICAL :: lgam!added:giovanni
   !
   REAL(DP) :: cscnorm
   !
   INTEGER ig, is, iv, jv, ia, inl, jnl
   REAL(DP) rsum
   COMPLEX(DP) csum
   REAL(DP), ALLOCATABLE:: temp(:)
   REAL(DP) :: icoeff !added:giovanni

!!!begin_added:giovanni
   IF (lgam) THEN
      icoeff = 2.d0
   ELSE
      icoeff = 1.d0
   END IF
!!!end_added:giovanni
!
   ALLOCATE (temp(ngw))

   DO ig = 1, ngw
      temp(ig) = DBLE(CONJG(cp(ig, i))*cp(ig, i))
   END DO

   IF (lgam) THEN!added:giovanni
      rsum = 2.d0*SUM(temp)
      IF (gstart == 2) rsum = rsum - temp(1)
   ELSE!added:giovanni
      rsum = SUM(temp) !added:giovanni
   END IF

   CALL mp_sum(rsum, intra_image_comm)

   DEALLOCATE (temp)
!
   IF (.not. bec%iscmplx) THEN
      DO is = 1, nvb
         DO iv = 1, nh(is)
            DO jv = 1, nh(is)
               IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                  DO ia = 1, na(is)
                     inl = ish(is) + (iv - 1)*na(is) + ia
                     jnl = ish(is) + (jv - 1)*na(is) + ia
                     rsum = rsum +                                        &
     &                    qq(iv, jv, is)*bec%rvec(inl, i)*bec%rvec(jnl, i)
                  END DO
               END IF
            END DO
         END DO
      END DO
   ELSE
!begin_added:giovanni
      csum = CMPLX(rsum, 0.d0)
      DO is = 1, nvb
         DO iv = 1, nh(is)
            DO jv = 1, nh(is)
               IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                  DO ia = 1, na(is)
                     inl = ish(is) + (iv - 1)*na(is) + ia
                     jnl = ish(is) + (jv - 1)*na(is) + ia
                     csum = csum +                                        &
       &                  CMPLX(qq(iv, jv, is), 0.d0)*CONJG(bec%cvec(inl, i))*bec%cvec(jnl, i)
                  END DO
               END IF
            END DO
         END DO
      END DO
      rsum = DBLE(csum)
!end_added:giovanni
   END IF
!
   cscnorm = SQRT(rsum)
!
   RETURN
END FUNCTION cscnorm
!
!
!-----------------------------------------------------------------------
SUBROUTINE denlcc(nnr, nspin, vxcr, sfac, drhocg, dcc)
!-----------------------------------------------------------------------
!
! derivative of non linear core correction exchange energy wrt cell
! parameters h
! Output in dcc
!
   USE kinds, ONLY: DP
   USE ions_base, ONLY: nsp
   USE reciprocal_vectors, ONLY: gstart, gx, ngs, g, ngm
   USE recvecs_indexes, ONLY: np
   USE cell_base, ONLY: omega, ainv, tpiba2
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE uspp_param, ONLY: upf
   USE cp_interfaces, ONLY: fwfft
   USE fft_base, ONLY: dfftp

   IMPLICIT NONE

   ! input

   INTEGER, INTENT(IN)   :: nnr, nspin
   REAL(DP)              :: vxcr(nnr, nspin)
   COMPLEX(DP)           :: sfac(ngs, nsp)
   REAL(DP)              :: drhocg(ngm, nsp)

   ! output

   REAL(DP), INTENT(OUT) ::  dcc(3, 3)

   ! local

   INTEGER     :: i, j, ig, is
   COMPLEX(DP) :: srhoc
   REAL(DP)    :: vxcc
   !
   COMPLEX(DP), ALLOCATABLE :: vxc(:)
!
   dcc = 0.0d0
   !
   ALLOCATE (vxc(nnr))
   !
   vxc(:) = vxcr(:, 1)
   !
   IF (nspin > 1) vxc(:) = vxc(:) + vxcr(:, 2)
   !
   CALL fwfft('Dense', vxc, dfftp)
   !
   DO i = 1, 3
      DO j = 1, 3
         DO ig = gstart, ngs
            srhoc = 0.0d0
            DO is = 1, nsp
               IF (upf(is)%nlcc) srhoc = srhoc + sfac(ig, is)*drhocg(ig, is)
            END DO
            vxcc = DBLE(CONJG(vxc(np(ig)))*srhoc)/SQRT(g(ig)*tpiba2)
            dcc(i, j) = dcc(i, j) + vxcc* &
  &                      2.d0*tpiba2*gx(i, ig)*                  &
  &                    (gx(1, ig)*ainv(j, 1) +                         &
  &                     gx(2, ig)*ainv(j, 2) +                         &
  &                     gx(3, ig)*ainv(j, 3))
         END DO
      END DO
   END DO

   DEALLOCATE (vxc)

   dcc = dcc*omega

   CALL mp_sum(dcc(1:3, 1:3), intra_image_comm)

   RETURN
END SUBROUTINE denlcc

!-----------------------------------------------------------------------
SUBROUTINE dotcsc(eigr, cp, ngw, n, lgam)!added:giovanni lgam
!-----------------------------------------------------------------------
!
   USE kinds, ONLY: DP
   USE ions_base, ONLY: na, nat
   USE io_global, ONLY: stdout
   USE reciprocal_vectors, ONLY: gstart
   USE cvan, ONLY: ish, nvb
   USE uspp, ONLY: nkb, qq
   USE uspp_param, ONLY: nh
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE twin_types !added:giovanni
   USE cp_interfaces, ONLY: nlsm1 !added:giovanni
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: ngw, n
   COMPLEX(DP) ::  eigr(ngw, nat), cp(ngw, n)
   LOGICAL, INTENT(IN) :: lgam
! local variables
   REAL(DP) rsum
   COMPLEX(DP) :: csum
   COMPLEX(DP) :: csc(n) ! automatic array
   COMPLEX(DP) :: temp(ngw) ! automatic array

   type(twin_matrix) ::  becp !modified:giovanni
   INTEGER i, kmax, nnn, k, ig, is, ia, iv, jv, inl, jnl
   COMPLEX(DP) :: icoeff

!!!begin_added:giovanni
   IF (lgam) THEN
      icoeff = CMPLX(2.d0, 0.d0)
!         becp%iscmplx=.false.!optional:giovanni
   ELSE
      icoeff = CMPLX(1.d0, 0.d0)
!         becp%iscmplx=.true.!optional:giovanni
   END IF
!!!end_added:giovanni
!
   CALL allocate_twin(becp, nkb, n, lgam)
!
!     < beta | phi > is real. only the i lowest:
!
   nnn = MIN(12, n)

   DO i = nnn, 1, -1
      kmax = i
      CALL nlsm1(i, 1, nvb, eigr, cp, becp, 1, lgam)
!
      DO k = 1, kmax
         DO ig = 1, ngw
            temp(ig) = CONJG(cp(ig, k))*cp(ig, i)
         END DO
         csc(k) = icoeff*SUM(temp)!modified:giovanni
         IF (gstart == 2) csc(k) = csc(k) - (icoeff - CMPLX(1.d0, 0.d0))*(temp(1)) !modified:giovanni
      END DO

      CALL mp_sum(csc(1:kmax), intra_image_comm)

      IF (lgam) THEN
         DO k = 1, kmax
            rsum = 0.d0
            DO is = 1, nvb
               DO iv = 1, nh(is)
                  DO jv = 1, nh(is)
                     DO ia = 1, na(is)
                        inl = ish(is) + (iv - 1)*na(is) + ia
                        jnl = ish(is) + (jv - 1)*na(is) + ia
                        rsum = rsum +                                    &
    &                   qq(iv, jv, is)*becp%rvec(inl, i)*becp%rvec(jnl, k)
                     END DO
                  END DO
               END DO
            END DO
            csc(k) = csc(k) + CMPLX(rsum, 0.d0)
         END DO
!
         WRITE (stdout, '("dotcsc =",12f18.15)') (DBLE(csc(k)), k=1, i)
!
      ELSE
         DO k = 1, kmax
            csum = 0.d0
            DO is = 1, nvb
               DO iv = 1, nh(is)
                  DO jv = 1, nh(is)
                     DO ia = 1, na(is)
                        inl = ish(is) + (iv - 1)*na(is) + ia
                        jnl = ish(is) + (jv - 1)*na(is) + ia
                        csum = csum +                                    &
    &                   qq(iv, jv, is)*CONJG(becp%cvec(inl, i))*becp%cvec(jnl, k)
                     END DO
                  END DO
               END DO
            END DO
            csc(k) = csc(k) + csum
         END DO
!
         WRITE (stdout, '("dotcsc =",12((f18.15)2x(f18.15)))') (csc(k), k=1, i)
!
      END IF

   END DO
   WRITE (stdout, *)
!
   call deallocate_twin(becp)
!
   RETURN
END SUBROUTINE dotcsc

!-----------------------------------------------------------------------
SUBROUTINE dotcsv(csv, nbspx, nbsp, c, bec, v, bev, ngw, ibnd1, ibnd2, lgam)
!-----------------------------------------------------------------------
!
   USE kinds, ONLY: DP
   USE ions_base, ONLY: na
   USE reciprocal_vectors, ONLY: gstart
   USE cvan, ONLY: ish, nvb
   USE uspp, ONLY: qq
   USE uspp_param, ONLY: nh
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE twin_types
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: ngw, nbsp, nbspx, ibnd1, ibnd2
   COMPLEX(DP) ::  c(ngw, nbspx), v(ngw, nbspx)
   COMPLEX(DP), INTENT(OUT) ::  csv
   type(twin_matrix) :: bec, bev
   LOGICAL :: lgam

   ! local variables
   COMPLEX(DP) temp(ngw) ! automatic array

   INTEGER ig, is, ia, iv, jv, inl, jnl
   REAL(DP) :: icoeff

   IF (lgam) THEN
      icoeff = 2.d0
   ELSE
      icoeff = 1.d0
   END IF

!      write(6,*) ubound(c), lbound(c), ubound(v), lbound(v)
!
!     < beta | c > is real. only the i lowest:
!
   DO ig = 1, ngw
      !
      temp(ig) = CONJG(c(ig, ibnd1))*v(ig, ibnd2)
      !
   END DO
   !
   csv = icoeff*(SUM(temp))
   IF (gstart == 2 .and. lgam) csv = csv - DBLE(temp(1))
   IF (lgam) csv = DBLE(csv)

   CALL mp_sum(csv, intra_image_comm)

   IF (lgam) THEN
      DO is = 1, nvb
         DO iv = 1, nh(is)
            DO jv = 1, nh(is)
               DO ia = 1, na(is)
                  inl = ish(is) + (iv - 1)*na(is) + ia
                  jnl = ish(is) + (jv - 1)*na(is) + ia
                  csv = csv + qq(iv, jv, is)*bec%rvec(inl, ibnd1)*bev%rvec(jnl, ibnd2)
               END DO
            END DO
         END DO
      END DO
   ELSE
      DO is = 1, nvb
         DO iv = 1, nh(is)
            DO jv = 1, nh(is)
               DO ia = 1, na(is)
                  inl = ish(is) + (iv - 1)*na(is) + ia
                  jnl = ish(is) + (jv - 1)*na(is) + ia
                  csv = csv + qq(iv, jv, is)*(bec%cvec(inl, ibnd1))*CONJG(bev%cvec(jnl, ibnd2))
               END DO
            END DO
         END DO
      END DO
   END IF
   !
   RETURN
END SUBROUTINE dotcsv
!

!-----------------------------------------------------------------------
SUBROUTINE compute_manifold_overlap(cstart, c, becstart, bec, ngwx, n, sss) !added:giovanni
!-----------------------------------------------------------------------
!
!  computes overlap between the manifold represented by cstart and that of c0
!  useful to check how far we are from the starting guess, during functional oprimization
!  Might even provide a convergence criterion.
!
   USE kinds, ONLY: DP
   USE constants, ONLY: pi, fpi
   USE gvecw, ONLY: ngw
   USE mp, ONLY: mp_sum
   USE control_flags, ONLY: gamma_only, do_wf_cmplx
   USE electrons_base, ONLY: iupdwn, nupdwn, nspin, nudx, nbspx, nbsp
   USE twin_types

   INTEGER, INTENT(IN) :: ngwx, n
   COMPLEX(DP), INTENT(INOUT) :: cstart(ngwx, nbspx), c(ngwx, nbspx)
   COMPLEX(DP), INTENT(OUT) :: sss(nspin)
   TYPE(twin_matrix) :: becstart, bec

   COMPLEX(DP) :: ss(nudx, nudx, nspin)
   LOGICAL :: lgam
   INTEGER :: i, j, isp

   lgam = gamma_only .and. .not. do_wf_cmplx
   !write(6,*) ubound(c), lbound(c), ubound(cstart), lbound(cstart)
   ss = CMPLX(0.d0, 0.d0)
   !
   DO isp = 1, nspin
      !
      DO i = 1, nupdwn(isp)
         !
         DO j = 1, nupdwn(isp)
            !
            call dotcsv(ss(j, i, isp), nbspx, nbsp, cstart, becstart, c, bec, ngw, iupdwn(isp) + j - 1, iupdwn(isp) + i - 1, lgam)
            !
         END DO
         !
      END DO
      !
   END DO
   !
   DO isp = 1, nspin
      !
      sss(isp) = 0.
      !
      DO j = 1, nupdwn(isp)
         !
         DO i = 1, nupdwn(isp)
            !
            sss(isp) = sss(isp) + abs(ss(i, j, isp))**2
!               write(*,*) ss(i,j,isp)-CONJG(ss(j,i,isp))
            !
         END DO
         !
      END DO
      !
      ! renormalize by number of electrons of that spin
      sss(isp) = sss(isp)/nupdwn(isp)
      !
   END DO
   !
END SUBROUTINE compute_manifold_overlap

!-----------------------------------------------------------------------
SUBROUTINE compute_overlap(c, ngwx, n, ss) !added:giovanni
!-----------------------------------------------------------------------
!
!  computes overlap matrix for the non-orthogonal case
!
!
   USE kinds, ONLY: DP
   USE constants, ONLY: pi, fpi
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: gstart
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE control_flags, ONLY: gamma_only, do_wf_cmplx
   USE electrons_base, ONLY: iupdwn, nupdwn, nspin, nudx

   INTEGER, INTENT(IN) :: ngwx, n
   COMPLEX(DP), INTENT(INOUT) :: c(ngwx, n)
   COMPLEX(DP), INTENT(OUT) :: ss(nudx, nudx, nspin)

   LOGICAL :: lgam
   INTEGER :: ig, i, j, isp
   REAL(DP) :: icoeff

   lgam = gamma_only .and. .not. do_wf_cmplx

   IF (lgam) THEN
      icoeff = 2.d0
   ELSE
      icoeff = 1.d0
   END IF

   ss = CMPLX(0.d0, 0.d0)
   !
   DO isp = 1, nspin
      !
      DO i = 1, nupdwn(isp)
         !
!             DO ig=1,ngw
!                !
!                ss(i,i,isp)=ss(i,i,isp)+icoeff*DBLE(CONJG(c(ig,iupdwn(isp)+i-1))* &
!                         c(ig,iupdwn(isp)+i-1))
!                !
!             END DO
!             !
!             IF(gstart==2) THEN
!                !
!                ss(i,i,isp)=ss(i,i,isp)-(icoeff-1.d0)*DBLE(CONJG(c(1,iupdwn(isp)+i-1))* &
!                            c(1,iupdwn(isp)+i-1))
!                !
!             ENDIF
         !
         DO j = 1, nupdwn(isp)
            !
            DO ig = 1, ngw
               !
               ss(j, i, isp) = ss(j, i, isp) + icoeff*CONJG(c(ig, iupdwn(isp) + j - 1))* &
                               c(ig, iupdwn(isp) + i - 1)
               !
            END DO
            !
            IF (gstart == 2) THEN
               ss(j, i, isp) = ss(j, i, isp) - (icoeff - 1.d0)*CONJG(c(1, iupdwn(isp) + j - 1))* &
                               c(1, iupdwn(isp) + i - 1)
            END IF
            !
            !ss(i,j,isp) = CONJG(ss(j,i,isp))
            !
         END DO
         !
      END DO
      !
      !
      !DO i=1,nupdwn(isp)
      !   c(:,iupdwn(isp)+i-1) = c(:,iupdwn(isp)+i-1)/sqrt(abs(ss(i,i,isp)))
      !   DO j=1,i-1
      !      ss(j,i,isp) = ss(j,i,isp)/sqrt(abs(ss(i,i,isp)))
      !      ss(i,j,isp) = CONJG(ss(j,i,isp))
      !   ENDDO
      !   ss(i,i,isp)=CMPLX(1.d0,0.d0)
      !ENDDO
   END DO
   !
   CALL mp_sum(ss(1:nudx, 1:nudx, 1:nspin), intra_image_comm)
   if (lgam) THEN
      ss(1:nudx, 1:nudx, 1:nspin) = DBLE(ss(1:nudx, 1:nudx, 1:nspin))
   END IF
   !
END SUBROUTINE compute_overlap

!-----------------------------------------------------------------------
SUBROUTINE invert_overlap(ss, iss, iflag) !added:giovanni
!-----------------------------------------------------------------------
!
!  contracs part of hamiltonian matrix with inverse overlap matrix
!
!
   USE kinds, ONLY: DP
   USE electrons_base, ONLY: nupdwn, nspin, nudx
   USE io_global, ONLY: ionode

   COMPLEX(DP), INTENT(IN) :: ss(nudx, nudx, nspin)
   COMPLEX(DP), INTENT(OUT) :: iss(nudx, nudx, nspin)
   INTEGER, INTENT(IN) :: iflag
   INTEGER, dimension(:), allocatable :: ipiv
   COMPLEX, DIMENSION(:), allocatable :: work

   INTEGER :: isp, info, lwork
   !
   lwork = nudx*nudx + 10
   allocate (ipiv(nudx))
   allocate (work(lwork))
   iss = CMPLX(0.d0, 0.d0)
   iss(1:nudx, 1:nudx, 1:nspin) = ss(1:nudx, 1:nudx, 1:nspin)
   do isp = 1, nspin
      call zgetrf(nupdwn(isp), nupdwn(isp), iss(1, 1, isp), nudx, ipiv, info)
      call zgetri(nupdwn(isp), iss(1, 1, isp), nudx, ipiv, work, lwork, info)
      if (ionode) THEN
         write (6, *) "inverted?", matmul(ss(1:nupdwn(isp), 1:nupdwn(isp), isp), &
                                          iss(1:nupdwn(isp), 1:nupdwn(isp), isp))
         write (6, *) "overlap", ss(1:nupdwn(isp), 1:nupdwn(isp), isp)
         write (6, *) "inverse_overlap", iss(1:nupdwn(isp), 1:nupdwn(isp), isp)
      END IF
   end do
   !
   deallocate (ipiv)
   deallocate (work)
   !
END SUBROUTINE invert_overlap

!-----------------------------------------------------------------------
SUBROUTINE compute_duals(c, cd, n, iflag) !added:giovanni
!-----------------------------------------------------------------------
!
!  contracs part of hamiltonian matrix with inverse overlap matrix
!
!
   USE kinds, ONLY: DP
   USE constants, ONLY: pi, fpi
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: gstart
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE control_flags, ONLY: gamma_only, do_wf_cmplx
   USE electrons_base, ONLY: iupdwn, nupdwn, nspin, nudx

   INTEGER, INTENT(IN) :: iflag, n
   COMPLEX(DP) :: c(ngw, n)
   COMPLEX(DP), INTENT(OUT) :: cd(ngw, n)
   COMPLEX(DP) :: temp1, temp2
   COMPLEX(DP), allocatable :: ss(:, :, :), iss(:, :, :)
   LOGICAL :: lgam

   INTEGER :: i, j, isp, cdindex

   lgam = gamma_only .and. .not. do_wf_cmplx
   !
   allocate (ss(nudx, nudx, nspin), iss(nudx, nudx, nspin))
   ss = 0.d0
   iss = 0.d0
   !
   call compute_overlap(c, ngw, n, ss)
   call invert_overlap(ss, iss, iflag)
   !
   cd(1:ngw, 1:n) = CMPLX(0.d0, 0.d0)
   !
   do isp = 1, nspin
      do i = 1, nupdwn(isp)
         cdindex = iupdwn(isp) + i - 1
         do j = 1, nupdwn(isp)
            cd(1:ngw, cdindex) = cd(1:ngw, cdindex) + CONJG(iss(i, j, isp))*c(1:ngw, iupdwn(isp) + j - 1)
         end do
      end do
   end do
   !
   temp1 = 0.d0
   temp2 = 0.d0
   !
   do i = 1, ngw
      temp1 = temp1 + 2.d0*CONJG(cd(i, 1))*c(i, 1)
      temp2 = temp2 + 2.d0*CONJG(cd(i, 2))*c(i, 1)
   end do
   if (gstart == 2) THEN
      temp1 = temp1 - CONJG(cd(1, 1))*c(1, 1)
      temp2 = temp2 - CONJG(cd(1, 2))*c(1, 1)
   END IF
   call mp_sum(temp1, intra_image_comm)
   call mp_sum(temp2, intra_image_comm)
!       write(6,*) "checkdiag", temp1, temp2
   deallocate (ss, iss)
END SUBROUTINE compute_duals

!-----------------------------------------------------------------------
SUBROUTINE times_overlap(c, cin, cd, n, iflag) !added:giovanni
!-----------------------------------------------------------------------
!
!  contracs part of hamiltonian matrix with inverse overlap matrix
!
!
   USE kinds, ONLY: DP
   USE constants, ONLY: pi, fpi
   USE gvecw, ONLY: ngw
   USE mp, ONLY: mp_sum
   USE control_flags, ONLY: gamma_only, do_wf_cmplx
   USE electrons_base, ONLY: iupdwn, nupdwn, nspin, nudx

   INTEGER, INTENT(IN) :: iflag, n
   COMPLEX(DP) :: c(ngw, n)
   COMPLEX(DP), INTENT(IN) :: cin(ngw, n)
   COMPLEX(DP), intent(OUT) :: cd(ngw, n)
!       COMPLEX(DP) :: ss(nudx,nudx,nspin)
   COMPLEX(DP), allocatable :: ss(:, :, :), iss(:, :, :)
   LOGICAL :: lgam

   INTEGER :: i, j, isp, cdindex

   lgam = gamma_only .and. .not. do_wf_cmplx
   !
   allocate (ss(nudx, nudx, nspin), iss(nudx, nudx, nspin))
   ss = 0.d0
   iss = 0.d0
   !
   call compute_overlap(c, ngw, n, ss)
   IF (iflag .lt. 0) THEN
      call invert_overlap(ss, iss, iflag)
      ss = CONJG(iss)
   END IF
   !
   cd(1:ngw, 1:n) = CMPLX(0.d0, 0.d0)
   !
   do isp = 1, nspin
      do i = 1, nupdwn(isp)
         cdindex = iupdwn(isp) + i - 1
         do j = 1, nupdwn(isp)
            cd(1:ngw, cdindex) = cd(1:ngw, cdindex) + (ss(i, j, isp))*cin(1:ngw, iupdwn(isp) + j - 1)
         end do
      end do
   end do
   !
!       write(6,*) "checkdiag", temp1, temp2
   deallocate (ss, iss)
END SUBROUTINE times_overlap

!-----------------------------------------------------------------------
FUNCTION enkin_non_ortho(c, cdual, ngwx, f, n)
!-----------------------------------------------------------------------
   !
   ! calculation of kinetic energy term
   !
   USE kinds, ONLY: DP
   USE constants, ONLY: pi, fpi
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: gstart
   USE gvecw, ONLY: ggp
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE cell_base, ONLY: tpiba2
   USE control_flags, ONLY: gamma_only, do_wf_cmplx

   IMPLICIT NONE

   REAL(DP)                :: enkin_non_ortho

   ! input

   INTEGER, INTENT(IN) :: ngwx, n
   COMPLEX(DP), INTENT(IN) :: c(ngwx, n), cdual(ngwx, n)
   REAL(DP), INTENT(IN) :: f(n)
   !
   ! local

   INTEGER  :: ig, i
   REAL(DP) :: sk(n)  ! automatic array
   LOGICAL :: lgam !added:giovanni
   REAL(DP) :: icoeff

   lgam = gamma_only .and. .not. do_wf_cmplx

   IF (lgam) THEN
      icoeff = 1.d0
   ELSE
      icoeff = 0.5d0
   END IF
   !
   ! matrix kk should be a global object, to allocate at the beginning
   !
   !
#ifdef __DEBUG_NONORTHO
   !
   ! compute the kinetic matrix kk, to contract with the
   ! inverse overlap matrix
   !
   kk = CMPLX(0.d0, 0.d0)
   !
   DO isp = 1, nspin
      !
      DO i = 1, nupdwn(isp)
         DO j = 1, i - 1
            DO ig = gstart, ngw
               kk(j, i, isp) = kk(j, i, isp) + CONJG(c(ig, iupdwn(isp) + j - 1))* &
                               c(ig, iupdwn(isp) + i - 1)*ggp(ig)
            END DO
            kk(i, j, isp) = CONJG(kk(j, i, isp))
         END DO
         !
         DO ig = gstart, ngw
            kk(i, i, isp) = kk(i, i, isp) + DBLE(CONJG(c(ig, iupdwn(isp) + i - 1))* &
                                                 c(ig, iupdwn(isp) + i - 1))*ggp(ig)
         END DO
         !
      END DO
      !
   END DO
   !
   CALL mp_sum(kk(1:nudx, 1:nudx, 1:nspin), intra_image_comm)
   !
   kk = kk*icoeff
   !
#endif
   !
   DO i = 1, n
      sk(i) = 0.0d0
      DO ig = gstart, ngw
         sk(i) = sk(i) + DBLE(CONJG(cdual(ig, i))*c(ig, i))*ggp(ig)
      END DO
   END DO
   !
   CALL mp_sum(sk(1:n), intra_image_comm)
   !
   enkin_non_ortho = 0.0d0
   DO i = 1, n
      enkin_non_ortho = enkin_non_ortho + f(i)*sk(i)
   END DO
   !

   ! ... reciprocal-space vectors are in units of alat/(2 pi) so a
   ! ... multiplicative factor (2 pi/alat)**2 is required

   enkin_non_ortho = enkin_non_ortho*tpiba2*icoeff
!
   RETURN
END FUNCTION enkin_non_ortho

!-----------------------------------------------------------------------
FUNCTION enkin(c, ngwx, f, n)
!-----------------------------------------------------------------------
   !
   ! calculation of kinetic energy term
   !
   USE kinds, ONLY: DP
   USE constants, ONLY: pi, fpi
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: gstart
   USE gvecw, ONLY: ggp
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE cell_base, ONLY: tpiba2
   USE control_flags, ONLY: gamma_only, do_wf_cmplx

   IMPLICIT NONE

   REAL(DP)                :: enkin

   ! input

   INTEGER, INTENT(IN) :: ngwx, n
   COMPLEX(DP), INTENT(IN) :: c(ngwx, n)
   REAL(DP), INTENT(IN) :: f(n)
   !
   ! local

   INTEGER  :: ig, i
   REAL(DP) :: sk(n)  ! automatic array
   LOGICAL :: lgam !added:giovanni
   REAL(DP) :: icoeff

   lgam = gamma_only .and. .not. do_wf_cmplx

   IF (lgam) THEN
      icoeff = 1.d0
   ELSE
      icoeff = 0.5d0
   END IF
   !
   ! matrix kk should be a global object, to allocate at the beginning
   !
   !
   !
   DO i = 1, n
      sk(i) = 0.0d0
      DO ig = gstart, ngw
         sk(i) = sk(i) + DBLE(CONJG(c(ig, i))*c(ig, i))*ggp(ig)
      END DO
   END DO
   !
   CALL mp_sum(sk(1:n), intra_image_comm)
   !
   enkin = 0.0d0
   DO i = 1, n
      enkin = enkin + f(i)*sk(i)
   END DO
   !
   ! ... reciprocal-space vectors are in units of alat/(2 pi) so a
   ! ... multiplicative factor (2 pi/alat)**2 is required

   enkin = enkin*tpiba2*icoeff
!
   RETURN
END FUNCTION enkin
!
!-----------------------------------------------------------------------
FUNCTION enkin_new(c, ngwx, f, n, nspin, nudx, iupdwn, nupdwn)
!-----------------------------------------------------------------------
   !
   ! calculation of kinetic energy term
   !
   USE kinds, ONLY: DP
   USE constants, ONLY: pi, fpi
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: gstart
   USE gvecw, ONLY: ggp
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE cell_base, ONLY: tpiba2
   USE control_flags, ONLY: gamma_only, do_wf_cmplx
!       USE electrons_base,     ONLY: iupdwn, nupdwn, nspin, nudx

   IMPLICIT NONE

   REAL(DP)                :: enkin_new

   ! input

   INTEGER, INTENT(IN) :: ngwx, n, nudx, nspin
   INTEGER, INTENT(IN) :: iupdwn(nspin), nupdwn(nspin)
   COMPLEX(DP), INTENT(IN) :: c(ngwx, n)
   REAL(DP), INTENT(IN) :: f(n)
   !
   ! local

   INTEGER  :: ig, i
   REAL(DP) :: sk(n)  ! automatic array
   LOGICAL :: lgam !added:giovanni
   REAL(DP) :: icoeff

   lgam = gamma_only .and. .not. do_wf_cmplx

   IF (lgam) THEN
      icoeff = 1.d0
   ELSE
      icoeff = 0.5d0
   END IF
   !
   ! matrix kk should be a global object, to allocate at the beginning
   !
   !
   !
   DO i = 1, n
      sk(i) = 0.0d0
      DO ig = gstart, ngw
         sk(i) = sk(i) + DBLE(CONJG(c(ig, i))*c(ig, i))*ggp(ig)
      END DO
   END DO
   !
   CALL mp_sum(sk(1:n), intra_image_comm)
   !
   enkin_new = 0.0d0
   DO i = 1, n
      enkin_new = enkin_new + f(i)*sk(i)
   END DO
   !
   ! ... reciprocal-space vectors are in units of alat/(2 pi) so a
   ! ... multiplicative factor (2 pi/alat)**2 is required

   enkin_new = enkin_new*tpiba2*icoeff
!
   RETURN
END FUNCTION enkin_new
!
!-----------------------------------------------------------------------
SUBROUTINE enkin_dens(c, ngwx, f)
!-----------------------------------------------------------------------
   !
   ! calculation of kinetic energy term
   !
   USE kinds, ONLY: DP
   USE constants, ONLY: pi, fpi
!       USE gvecw,              ONLY: ngw
   USE reciprocal_vectors, ONLY: gx
!       USE gvecw,              ONLY: ggp
   USE gvecp, only: ng => ngm
   USE mp, ONLY: mp_sum
   USE control_flags, ONLY: gamma_only, do_wf_cmplx
!       USE electrons_base,     ONLY: iupdwn, nupdwn, nspin, nudx

   IMPLICIT NONE

   ! input

   INTEGER, INTENT(IN) :: ngwx
   COMPLEX(DP), INTENT(INOUT) :: c(ngwx)
   REAL(DP), INTENT(IN) :: f
   !
   ! local

   INTEGER  :: ig
!       COMPLEX(DP), allocatable :: sk(:)  ! automatic array
   LOGICAL :: lgam !added:giovanni
   REAL(DP) :: icoeff

   lgam = gamma_only .and. .not. do_wf_cmplx

   IF (lgam) THEN
      icoeff = 1.d0
   ELSE
      icoeff = 0.5d0
   END IF
   !
   ! matrix kk should be a global object, to allocate at the beginning
   !
   !
!       allocate(sk(ngwx))
   !
   DO ig = 1, ng
      !
      c(ig) = DBLE(CONJG(c(ig))*c(ig))*(gx(1, ig)**2 + gx(2, ig)**2 + gx(3, ig)**2)
      !
   END DO
   !
   RETURN
END SUBROUTINE enkin_dens
!
!-----------------------------------------------------------------------
SUBROUTINE gausin(eigr, cm)
!-----------------------------------------------------------------------
!
! initialize wavefunctions with gaussians - edit to fit your system
!
   USE kinds, ONLY: DP
   USE ions_base, ONLY: na, nat
   USE electrons_base, ONLY: n => nbsp
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: gx, g
!
   IMPLICIT NONE
!
   COMPLEX(DP) eigr(ngw, nat), cm(ngw, n)
   REAL(DP) sigma, auxf
   INTEGER nband, is, ia, ig, isa
!
   sigma = 12.0d0
   nband = 0
!!!      do is=1,nsp
   isa = 0
   is = 1
   DO ia = 1, na(is)
! s-like gaussians
      nband = nband + 1
      DO ig = 1, ngw
         auxf = EXP(-g(ig)/sigma**2)
         cm(ig, nband) = auxf*eigr(ig, ia + isa)
      END DO
! px-like gaussians
      nband = nband + 1
      DO ig = 1, ngw
         auxf = EXP(-g(ig)/sigma**2)
         cm(ig, nband) = auxf*eigr(ig, ia + isa)*gx(1, ig)
      END DO
! py-like gaussians
      nband = nband + 1
      DO ig = 1, ngw
         auxf = EXP(-g(ig)/sigma**2)
         cm(ig, nband) = auxf*eigr(ig, ia + isa)*gx(2, ig)
      END DO
! pz-like gaussians
      nband = nband + 1
      DO ig = 1, ngw
         auxf = EXP(-g(ig)/sigma**2)
         cm(ig, nband) = auxf*eigr(ig, ia + isa)*gx(3, ig)
      END DO
   END DO
   isa = isa + na(is)
   is = 2
   DO ia = 1, na(is)
! s-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)
!            end do
! px-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(1,ig)
!            end do
! py-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(2,ig)
!            end do
! pz-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(3,ig)
!            end do
! dxy-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(1,ig)*gx(2,ig)
!            end do
! dxz-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(1,ig)*gx(3,ig)
!            end do
! dxy-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*gx(2,ig)*gx(3,ig)
!            end do
! dx2-y2-like gaussians
!            nband=nband+1
!            do ig=1,ngw
!               auxf=exp(-g(ig)/sigma**2)
!               cm(ig,nband)=auxf*eigr(ig,ia+isa)*                        &
!     &              (gx(1,ig)**2-gx(2,ig)**2)
!            end do
   END DO
!!!      end do
   RETURN
END SUBROUTINE gausin
!
!-------------------------------------------------------------------------
SUBROUTINE scalar_us(bec, nkbx, betae, cp, ngwx, i, csc, n, lgam) !added:giovanni SUBROUTINE scalar_us
!-----------------------------------------------------------------------
!     requires in input the updated bec(k) for k<i
!     on output: bec(i) is recalculated
!
   USE ions_base, ONLY: na
   USE cvan, ONLY: nvb, ish
   USE uspp, ONLY: nkb, qq
   USE uspp_param, ONLY: nh
   USE electrons_base, ONLY: ispin
   USE gvecw, ONLY: ngw
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE kinds, ONLY: DP
   USE reciprocal_vectors, ONLY: gstart
   USE twin_types !added:giovanni
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: i, nkbx, ngwx, n
   COMPLEX(DP) :: betae(ngwx, nkb)
   COMPLEX(DP)    :: cp(ngwx, n)
   type(twin_matrix) :: bec!( nkbx, n )!modified:giovanni
   COMPLEX(DP)    :: csc(n) !modified:giovanni
   LOGICAL :: lgam !added:giovanni
   INTEGER :: k, kmax, ig, is, iv, jv, ia, inl, jnl
   REAL(DP)    :: rsum
   COMPLEX(DP) :: csum
   REAL(DP), ALLOCATABLE :: temp(:)
   COMPLEX(DP), ALLOCATABLE :: temp_c(:) !added:giovanni

!!!begin_added:giovanni
   IF (lgam) THEN
      ALLOCATE (temp(ngw))
      temp = 0.d0
   ELSE
      ALLOCATE (temp_c(ngw))
      temp_c = CMPLX(0.d0, 0.d0)
   END IF
!!!end_added:giovanni
   !
   !     calculate csc(k)=<cp(i)|cp(k)>,  k<i --->  <cp(k)|cp(i)>
   !
   kmax = i - 1
   !
!$omp parallel default(shared), private( temp, k, ig )

!$omp do
   IF (lgam) THEN
      DO k = 1, kmax
         csc(k) = CMPLX(0.0d0, 0.d0)
         IF (ispin(i) .EQ. ispin(k)) THEN
            DO ig = 1, ngw
               temp(ig) = DBLE(CONJG(cp(ig, i))*cp(ig, k))
            END DO
            csc(k) = CMPLX(2.d0*SUM(temp, ngw), 0.d0)
            IF (gstart == 2) csc(k) = csc(k) - CMPLX(temp(1), 0.d0)
         END IF
      END DO
   ELSE
      !begin_added:giovanni
      DO k = 1, kmax
         csc(k) = CMPLX(0.0d0, 0.d0)
         IF (ispin(i) .EQ. ispin(k)) THEN
            DO ig = 1, ngw
               temp_c(ig) = CONJG(cp(ig, k))*cp(ig, i)
            END DO
            csc(k) = SUM(temp_c, ngw)
         END IF
      END DO
      !end_added:giovanni
   END IF
!$omp end do
   IF (lgam) THEN
      DEALLOCATE (temp)
   ELSE
      DEALLOCATE (temp_c)
   END IF
!$omp end parallel

   CALL mp_sum(csc(1:kmax), intra_image_comm)
   IF (lgam) THEN
      ALLOCATE (temp(ngw))
   ELSE
      ALLOCATE (temp_c(ngw))
   END IF
!
!     calculate csc(k)=<cp(i)|S|cp(k)>,  k<i -----> csc(k) = <cp(k)|S|cp(i)>
!
   IF (.not. bec%iscmplx) THEN
      DO k = 1, kmax
         IF (ispin(i) .EQ. ispin(k)) THEN
            rsum = 0.d0
            DO is = 1, nvb
               DO iv = 1, nh(is)
                  DO jv = 1, nh(is)
                     IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                        DO ia = 1, na(is)
                           inl = ish(is) + (iv - 1)*na(is) + ia
                           jnl = ish(is) + (jv - 1)*na(is) + ia
                           rsum = rsum + qq(iv, jv, is)*bec%rvec(inl, i)*bec%rvec(jnl, k)
                        END DO
                     END IF
                  END DO
               END DO
            END DO
            csc(k) = csc(k) + CMPLX(rsum, 0.d0)
         END IF
      END DO
   ELSE
      DO k = 1, kmax
         IF (ispin(i) .EQ. ispin(k)) THEN
            csum = CMPLX(0.d0, 0.d0)
            DO is = 1, nvb
               DO iv = 1, nh(is)
                  DO jv = 1, nh(is)
                     IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
! write(6,*) "updating via qq" ! added
                        DO ia = 1, na(is)
                           inl = ish(is) + (iv - 1)*na(is) + ia
                           jnl = ish(is) + (jv - 1)*na(is) + ia
                           csum = csum + CMPLX(qq(iv, jv, is), 0.d0)*(bec%cvec(inl, i))*CONJG(bec%cvec(jnl, k))
                        END DO
                     END IF
                  END DO
               END DO
            END DO
            csc(k) = csc(k) + csum
         END IF
      END DO
   END IF

   IF (lgam) THEN
      DEALLOCATE (temp)
   ELSE
      DEALLOCATE (temp_c) !added:giovanni
   END IF
!
!       write(6,*) "bec", bec%rvec
   RETURN
END SUBROUTINE scalar_us

!
!-------------------------------------------------------------------------
SUBROUTINE scalar_character(cp, n, ngwx, csc, nunit, lgam) !added:giovanni SUBROUTINE scalar_character
!-----------------------------------------------------------------------
!     requires in input the updated bec(k) for k<i
!     on output: bec(i) is recalculated
!
   USE reciprocal_vectors, ONLY: gx
   USE cell_base, ONLY: a3
   USE electrons_base, ONLY: ispin
   USE gvecw, ONLY: ngw
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE kinds, ONLY: DP
   USE reciprocal_vectors, ONLY: gstart
   USE twin_types !added:giovanni
   USE constants
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: ngwx, n, nunit
!       COMPLEX(DP) :: betae( ngwx, nkb )
   COMPLEX(DP)    :: cp(ngwx, n), phase
   REAL(DP) :: rvector(3)
!       type(twin_matrix) :: bec!( nkbx, n )!modified:giovanni
   COMPLEX(DP)    :: csc(n, n) !modified:giovanni
   LOGICAL :: lgam !added:giovanni
   INTEGER :: k, i, kmax, ig
   REAL(DP), ALLOCATABLE :: temp(:)
   COMPLEX(DP), ALLOCATABLE :: temp_c(:) !added:giovanni

!!!begin_added:giovanni
   IF (lgam) THEN
      ALLOCATE (temp(ngw))
      temp = 0.d0
   ELSE
      ALLOCATE (temp_c(ngw))
      temp_c = CMPLX(0.d0, 0.d0)
   END IF
!!!end_added:giovanni
   !
   !     calculate csc(k)=-ilog(<cp(k)|T| cp(k)>)
   kmax = n
   rvector(:) = 2.d0*pi*a3(:)/20.d0*nunit
   !
!$omp parallel default(shared), private( temp, k, ig )

!$omp do
   IF (lgam) THEN
      DO i = 1, kmax
         DO k = 1, kmax
            csc(k, i) = CMPLX(0.0d0, 0.d0)
            IF (ispin(i) .EQ. ispin(k)) THEN
               DO ig = 1, ngw
                  phase = exp(CMPLX(0.d0, gx(1, ig)*rvector(1) + gx(2, ig)*rvector(2) + gx(3, ig)*rvector(3)))
                  temp(ig) = DBLE(CONJG(cp(ig, k))*phase*cp(ig, i))
               END DO
               csc(k, i) = CMPLX(2.d0*SUM(temp, ngw), 0.d0)
               IF (gstart == 2) csc(k, i) = csc(k, i) - CMPLX(temp(1), 0.d0)
            END IF
         END DO
      END DO
   ELSE
      !begin_added:giovanni
      DO i = 1, kmax
         DO k = 1, kmax
            csc(k, i) = CMPLX(0.0d0, 0.d0)
            IF (ispin(i) .EQ. ispin(k)) THEN
               DO ig = 1, ngw
                  phase = exp(CMPLX(0.d0, gx(1, ig)*rvector(1) + gx(2, ig)*rvector(2) + gx(3, ig)*rvector(3)))
                  temp_c(ig) = CONJG(cp(ig, k))*phase*cp(ig, i)
               END DO
               csc(k, i) = SUM(temp_c, ngw)
            END IF
         END DO
      END DO
      !end_added:giovanni
   END IF
!$omp end do
   IF (lgam) THEN
      DEALLOCATE (temp)
   ELSE
      DEALLOCATE (temp_c)
   END IF
!$omp end parallel

   CALL mp_sum(csc(1:kmax, 1:kmax), intra_image_comm)
!       write(6,*) "bec", bec%rvec
   RETURN
END SUBROUTINE scalar_character

!-------------------------------------------------------------------------
SUBROUTINE gracsc(bec, nkbx, betae, cp, ngwx, i, csc, n, lgam) !added:giovanni lgam
!-----------------------------------------------------------------------
!     requires in input the updated bec(k) for k<i
!     on output: bec(i) is recalculated
!
   USE ions_base, ONLY: na
   USE cvan, ONLY: nvb, ish
   USE uspp, ONLY: nkb, nhsavb => nkbus, qq
   USE uspp_param, ONLY: nh
   USE electrons_base, ONLY: ispin
   USE gvecw, ONLY: ngw
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE kinds, ONLY: DP
   USE reciprocal_vectors, ONLY: gstart
   USE twin_types !added:giovanni
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: i, nkbx, ngwx, n
   COMPLEX(DP) :: betae(ngwx, nkb)
   COMPLEX(DP)    :: cp(ngwx, n)
   type(twin_matrix) :: bec!( nkbx, n )!modified:giovanni
   COMPLEX(DP)    :: csc(n) !modified:giovanni
   LOGICAL :: lgam !added:giovanni
   INTEGER :: k, kmax, ig, is, iv, jv, ia, inl, jnl
   REAL(DP)    :: rsum
   COMPLEX(DP) :: csum
   REAL(DP), ALLOCATABLE :: temp(:)
   COMPLEX(DP), ALLOCATABLE :: temp_c(:) !added:giovanni

!!!begin_added:giovanni
   IF (lgam) THEN
      ALLOCATE (temp(ngw))
      temp = 0.d0
   ELSE
      ALLOCATE (temp_c(ngw))
      temp_c = CMPLX(0.d0, 0.d0)
   END IF
!!!end_added:giovanni
   !
   !     calculate csc(k)=<cp(i)|cp(k)>,  k<i --->  <cp(k)|cp(i)>
   !
   kmax = i - 1
   !
!$omp parallel default(shared), private( temp, k, ig )

!$omp do
   IF (lgam) THEN
      DO k = 1, kmax
         csc(k) = CMPLX(0.0d0, 0.d0)
         IF (ispin(i) .EQ. ispin(k)) THEN
            DO ig = 1, ngw
               temp(ig) = DBLE(CONJG(cp(ig, i))*cp(ig, k))
            END DO
            csc(k) = CMPLX(2.d0*SUM(temp, ngw), 0.d0)
            IF (gstart == 2) csc(k) = csc(k) - CMPLX(temp(1), 0.d0)
         END IF
      END DO
   ELSE
      !begin_added:giovanni
      DO k = 1, kmax
         csc(k) = CMPLX(0.0d0, 0.d0)
         IF (ispin(i) .EQ. ispin(k)) THEN
            DO ig = 1, ngw
               temp_c(ig) = CONJG(cp(ig, k))*cp(ig, i)
            END DO
            csc(k) = SUM(temp_c, ngw)
         END IF
      END DO
      !end_added:giovanni
   END IF
!$omp end do
   IF (lgam) THEN
      DEALLOCATE (temp)
   ELSE
      DEALLOCATE (temp_c)
   END IF
!$omp end parallel

   CALL mp_sum(csc(1:kmax), intra_image_comm)
   IF (lgam) THEN
      ALLOCATE (temp(ngw))
   ELSE
      ALLOCATE (temp_c(ngw))
   END IF
   !
   !     calculate bec(i)=<cp(i)|beta> NO:giovanni--> <beta|cp(i)>
   !
   IF (lgam) THEN
      DO inl = 1, nhsavb
         DO ig = 1, ngw
            temp(ig) = DBLE(cp(ig, i)*CONJG(betae(ig, inl)))
         END DO
         bec%rvec(inl, i) = 2.d0*SUM(temp)!modified:giovanni
         IF (gstart == 2) bec%rvec(inl, i) = bec%rvec(inl, i) - temp(1)!modified:giovanni
      END DO
      CALL mp_sum(bec%rvec(1:nhsavb, i), intra_image_comm)
   ELSE
!begin_added:giovanni
      DO inl = 1, nhsavb
         DO ig = 1, ngw
            temp_c(ig) = cp(ig, i)*CONJG(betae(ig, inl))
         END DO
         IF (bec%iscmplx) then
            bec%cvec(inl, i) = SUM(temp_c)
         ELSE !added:giovanni:debug
            bec%rvec(inl, i) = DBLE(SUM(temp_c))
         END IF
      END DO
      IF (bec%iscmplx) then
         CALL mp_sum(bec%cvec(1:nhsavb, i), intra_image_comm)
      ELSE
         CALL mp_sum(bec%rvec(1:nhsavb, i), intra_image_comm)
      END IF
!end_added:giovanni
   END IF

!
!     calculate csc(k)=<cp(i)|S|cp(k)>,  k<i -----> csc(k) = <cp(k)|S|cp(i)>
!
   IF (.not. bec%iscmplx) THEN
      DO k = 1, kmax
         IF (ispin(i) .EQ. ispin(k)) THEN
            rsum = 0.d0
            DO is = 1, nvb
               DO iv = 1, nh(is)
                  DO jv = 1, nh(is)
                     IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                        DO ia = 1, na(is)
                           inl = ish(is) + (iv - 1)*na(is) + ia
                           jnl = ish(is) + (jv - 1)*na(is) + ia
                           rsum = rsum + qq(iv, jv, is)*bec%rvec(inl, i)*bec%rvec(jnl, k)
                        END DO
                     END IF
                  END DO
               END DO
            END DO
            csc(k) = csc(k) + CMPLX(rsum, 0.d0)
         END IF
      END DO
   ELSE
      DO k = 1, kmax
         IF (ispin(i) .EQ. ispin(k)) THEN
            csum = CMPLX(0.d0, 0.d0)
            DO is = 1, nvb
               DO iv = 1, nh(is)
                  DO jv = 1, nh(is)
                     IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
! write(6,*) "updating via qq" ! added
                        DO ia = 1, na(is)
                           inl = ish(is) + (iv - 1)*na(is) + ia
                           jnl = ish(is) + (jv - 1)*na(is) + ia
                           csum = csum + CMPLX(qq(iv, jv, is), 0.d0)*(bec%cvec(inl, i))*CONJG(bec%cvec(jnl, k))
                        END DO
                     END IF
                  END DO
               END DO
            END DO
            csc(k) = csc(k) + csum
         END IF
      END DO
   END IF
!
!     orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
!
!     corresponing bec:  bec(i)=<cp(i)|beta>-(csc(k))*<cp(k)|beta>
!
   IF (.not. bec%iscmplx) THEN
      DO k = 1, kmax
         DO inl = 1, nkbx
            bec%rvec(inl, i) = bec%rvec(inl, i) - DBLE(csc(k))*bec%rvec(inl, k)
         END DO
      END DO
   ELSE
!begin_added:giovanni
      DO k = 1, kmax
         DO inl = 1, nkbx
            bec%cvec(inl, i) = bec%cvec(inl, i) - (csc(k))*bec%cvec(inl, k)
         END DO
      END DO
!         write(6,*) "output complex bec", bec%cvec
!end_added:giovanni
   END IF

   IF (lgam) THEN
      DEALLOCATE (temp)
   ELSE
      DEALLOCATE (temp_c) !added:giovanni
   END IF
!
!       write(6,*) "bec", bec%rvec
   RETURN
END SUBROUTINE gracsc

!-------------------------------------------------------------------------
SUBROUTINE gracsc2(bec, nkbx, betae, cp, ngwx, i, k, csc, n, lgam) !added:giovanni lgam
!-----------------------------------------------------------------------
!     requires in input the updated bec(k) for k<i
!     on output: bec(i) is recalculated
!
   USE ions_base, ONLY: na
   USE cvan, ONLY: nvb, ish
   USE uspp, ONLY: nkb, nhsavb => nkbus, qq
   USE uspp_param, ONLY: nh
   USE electrons_base, ONLY: ispin
   USE gvecw, ONLY: ngw
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE kinds, ONLY: DP
   USE reciprocal_vectors, ONLY: gstart
   USE twin_types !added:giovanni
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: i, nkbx, ngwx, n, k
   COMPLEX(DP) :: betae(ngwx, nkb)
   COMPLEX(DP)    :: cp(ngwx, n)
   type(twin_matrix) :: bec!( nkbx, n )!modified:giovanni
   COMPLEX(DP)    :: csc(n) !modified:giovanni
   LOGICAL :: lgam !added:giovanni
   INTEGER :: ig, is, iv, jv, ia, inl, jnl
   REAL(DP)    :: rsum
   COMPLEX(DP) :: csum
   REAL(DP), ALLOCATABLE :: temp(:)
   COMPLEX(DP), ALLOCATABLE :: temp_c(:) !added:giovanni

!!!begin_added:giovanni
   IF (lgam) THEN
      ALLOCATE (temp(ngw))
      temp = 0.d0
   ELSE
      ALLOCATE (temp_c(ngw))
      temp_c = CMPLX(0.d0, 0.d0)
   END IF
!!!end_added:giovanni
   !
   !     calculate csc(k)=<cp(i)|cp(k)>,  k<i
   !
!$omp parallel default(shared), private( temp, k, ig )

!$omp do

   IF (lgam) THEN
      csc(k) = CMPLX(0.0d0, 0.d0)
      IF (ispin(i) .EQ. ispin(k)) THEN
         DO ig = 1, ngw
            temp(ig) = DBLE(CONJG(cp(ig, k))*cp(ig, i))
         END DO
         csc(k) = CMPLX(2.d0*SUM(temp, ngw), 0.d0)
         IF (gstart == 2) csc(k) = csc(k) - CMPLX(temp(1), 0.d0)
      END IF
   ELSE
      !begin_added:giovanni
      csc(k) = CMPLX(0.0d0, 0.d0)
      IF (ispin(i) .EQ. ispin(k)) THEN
         DO ig = 1, ngw
            temp_c(ig) = CONJG(cp(ig, k))*cp(ig, i)
! cmplx(cp(1,ig,k) * cp(1,ig,i) + cp(2,ig,k) * cp(2,ig,i),  &
!                                       cp(1,ig,k) * cp(2,ig,i) - cp(2,ig,k) * cp(1,ig,i), DP)
         END DO
         csc(k) = SUM(temp_c, ngw)
      END IF
      !end_added:giovanni
   END IF
!$omp end do
   IF (lgam) THEN
      DEALLOCATE (temp)
   ELSE
      DEALLOCATE (temp_c)
   END IF
!$omp end parallel

   CALL mp_sum(csc(k), intra_image_comm)
   IF (lgam) THEN
      ALLOCATE (temp(ngw))
   ELSE
      ALLOCATE (temp_c(ngw))
   END IF
   !
   !     calculate bec(i)=<cp(i)|beta>
   !
   IF (lgam) THEN
      DO inl = 1, nhsavb
         DO ig = 1, ngw
            temp(ig) = DBLE(CONJG(cp(ig, i))*betae(ig, inl))
         END DO
         bec%rvec(inl, i) = SUM(temp)!modified:giovanni
         IF (gstart == 2) bec%rvec(inl, i) = bec%rvec(inl, i) - temp(1)!modified:giovanni
      END DO
      CALL mp_sum(bec%rvec(1:nhsavb, i), intra_image_comm)
   ELSE
!begin_added:giovanni
      DO inl = 1, nhsavb
         DO ig = 1, ngw
            temp_c(ig) = CONJG(cp(ig, i))*betae(ig, inl)
! cmplx(cp(1,ig,i)* DBLE(betae(ig,inl)) +      &
!       &               cp(2,ig,i)*AIMAG(betae(ig,inl)),                &
!       &               cp(1,ig,i)*AIMAG(betae(ig,inl))-                &
!       &               cp(2,ig,i)*DBLE(betae(ig,inl)), DP)
         END DO
         IF (bec%iscmplx) then
            bec%cvec(inl, i) = SUM(temp_c)
         ELSE !added:giovanni:debug
            bec%rvec(inl, i) = DBLE(SUM(temp_c))
         END IF
      END DO
      IF (bec%iscmplx) then
         CALL mp_sum(bec%cvec(1:nhsavb, i), intra_image_comm)
      ELSE
         CALL mp_sum(bec%rvec(1:nhsavb, i), intra_image_comm)
      END IF
!end_added:giovanni
   END IF

!
!     calculate csc(k)=<cp(i)|S|cp(k)>,  k<i
!
   IF (.not. bec%iscmplx) THEN
      IF (ispin(i) .EQ. ispin(k)) THEN
         rsum = 0.d0
         DO is = 1, nvb
            DO iv = 1, nh(is)
               DO jv = 1, nh(is)
                  IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                     DO ia = 1, na(is)
                        inl = ish(is) + (iv - 1)*na(is) + ia
                        jnl = ish(is) + (jv - 1)*na(is) + ia
                        rsum = rsum + qq(iv, jv, is)*bec%rvec(inl, i)*bec%rvec(jnl, k)
                     END DO
                  END IF
               END DO
            END DO
         END DO
         csc(k) = csc(k) + CMPLX(rsum, 0.d0)
      END IF
   ELSE
      IF (ispin(i) .EQ. ispin(k)) THEN
         csum = CMPLX(0.d0, 0.d0)
         DO is = 1, nvb
            DO iv = 1, nh(is)
               DO jv = 1, nh(is)
                  IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                     DO ia = 1, na(is)
                        inl = ish(is) + (iv - 1)*na(is) + ia
                        jnl = ish(is) + (jv - 1)*na(is) + ia
                        csum = csum + CMPLX(qq(iv, jv, is), 0.d0)*CONJG(bec%cvec(inl, i))*(bec%cvec(jnl, k))
                     END DO
                  END IF
               END DO
            END DO
         END DO
         csc(k) = csc(k) + csum
      END IF
   END IF
!
!     orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
!
!     corresponing bec:  bec(i)=<cp(i)|beta>-csc(k)<cp(k)|beta>
!
   IF (.not. bec%iscmplx) THEN
      DO inl = 1, nkbx
         bec%rvec(inl, i) = bec%rvec(inl, i) - DBLE(csc(k))*bec%rvec(inl, k)
      END DO
   ELSE
!begin_added:giovanni
      DO inl = 1, nkbx
         bec%cvec(inl, i) = bec%cvec(inl, i) - CONJG(csc(k))*bec%cvec(inl, k)
      END DO
!end_added:giovanni
   END IF

   IF (lgam) THEN
      DEALLOCATE (temp)
   ELSE
      DEALLOCATE (temp_c) !added:giovanni
   END IF
!
   RETURN
END SUBROUTINE gracsc2

!-------------------------------------------------------------------------
SUBROUTINE smooth_csv_real(c, v, ngwx, csv, n)
!-----------------------------------------------------------------------

   USE gvecw, ONLY: ngw
   USE kinds, ONLY: DP
   USE reciprocal_vectors, ONLY: gstart
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: ngwx, n
   REAL(DP)    :: c(2, ngwx)
   REAL(DP)    :: v(2, ngwx, n)
   REAL(DP)    :: csv(n)
   INTEGER     :: k, ig
   REAL(DP), ALLOCATABLE :: temp(:)

   !
   !     calculate csv(k)=<c|v(k)>
   !
   ALLOCATE (temp(ngw))

   DO k = 1, n
      DO ig = 1, ngw
         temp(ig) = v(1, ig, k)*c(1, ig) + v(2, ig, k)*c(2, ig)
      END DO
      csv(k) = 2.0d0*SUM(temp)
      IF (gstart == 2) csv(k) = csv(k) - temp(1)
   END DO

   DEALLOCATE (temp)
!
   RETURN
END SUBROUTINE smooth_csv_real

!-------------------------------------------------------------------------
SUBROUTINE smooth_csv_twin(c, v, ngwx, csv, n)
!-----------------------------------------------------------------------

   USE gvecw, ONLY: ngw
   USE kinds, ONLY: DP
   USE reciprocal_vectors, ONLY: gstart
   USE twin_types
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: ngwx, n
   COMPLEX(DP)    :: c(ngwx)
   COMPLEX(DP)    :: v(ngwx, n)
   type(twin_matrix)    :: csv !( n )
!       LOGICAL :: lgam
   INTEGER     :: k, ig
   REAL(DP), ALLOCATABLE :: temp(:)
   COMPLEX(DP), ALLOCATABLE :: temp_c(:)
   !
   !     calculate csv(k)=<c|v(k)>
   !
   IF (.not. csv%iscmplx) THEN
      ALLOCATE (temp(ngw))
   ELSE
      ALLOCATE (temp_c(ngw))
   END IF

   IF (.not. csv%iscmplx) THEN
      DO k = 1, n
         DO ig = 1, ngw
            temp(ig) = DBLE(CONJG(v(ig, k))*c(ig))
         END DO
         csv%rvec(k, 1) = 2.0d0*SUM(temp)
         IF (gstart == 2) csv%rvec(k, 1) = csv%rvec(k, 1) - temp(1)
      END DO
   ELSE
      DO k = 1, n
         DO ig = 1, ngw
            temp_c(ig) = CONJG(v(ig, k))*c(ig)
         END DO
         csv%cvec(k, 1) = SUM(temp_c)
      END DO
   END IF

   IF (.not. csv%iscmplx) THEN
      DEALLOCATE (temp)
   ELSE
      DEALLOCATE (temp_c)
   END IF
!
   RETURN
END SUBROUTINE smooth_csv_twin

!-------------------------------------------------------------------------
SUBROUTINE grabec_real(becc, nkbx, betae, c, ngwx)
!-----------------------------------------------------------------------
   !
   !     on output: bec(i) is recalculated
   !
   USE uspp, ONLY: nkb, nhsavb => nkbus
   USE gvecw, ONLY: ngw
   USE kinds, ONLY: DP
   USE reciprocal_vectors, ONLY: gstart
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: nkbx, ngwx
   COMPLEX(DP) :: betae(ngwx, nkb)
   REAL(DP)    :: becc(nkbx), c(2, ngwx)
   INTEGER     :: ig, inl
   REAL(DP), ALLOCATABLE :: temp(:)
   !
   ALLOCATE (temp(ngw))
   !
   !     calculate becc=<c|beta>
   !
   DO inl = 1, nhsavb
      DO ig = 1, ngw
         temp(ig) = c(1, ig)*DBLE(betae(ig, inl)) +             &
  &               c(2, ig)*AIMAG(betae(ig, inl))
      END DO
      becc(inl) = 2.d0*SUM(temp)
      IF (gstart == 2) becc(inl) = becc(inl) - temp(1)
   END DO

   DEALLOCATE (temp)

   RETURN
END SUBROUTINE grabec_real

!-------------------------------------------------------------------------
SUBROUTINE grabec_twin(becc, nkbx, betae, c, ngwx, l2_bec)
!-----------------------------------------------------------------------
   !
   !     on output: bec(i) is recalculated
   !
   USE uspp, ONLY: nkb, nhsavb => nkbus
   USE gvecw, ONLY: ngw
   USE kinds, ONLY: DP
   USE reciprocal_vectors, ONLY: gstart
   USE twin_types
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: nkbx, ngwx, l2_bec
   COMPLEX(DP) :: betae(ngwx, nkb)
   !
   COMPLEX(DP)    ::  c(1, ngwx)
   type(twin_matrix) :: becc !( nkbx ),
   INTEGER     :: ig, inl
   REAL(DP), ALLOCATABLE :: temp(:)
   COMPLEX(DP), ALLOCATABLE :: temp_c(:)
   !
   !
   !     calculate becc=<c|beta>
   !
   IF (.not. becc%iscmplx) THEN
      ALLOCATE (temp(ngw))
      DO inl = 1, nhsavb
         DO ig = 1, ngw
            temp(ig) = DBLE(CONJG(c(1, ig))*betae(ig, inl))
         END DO
         becc%rvec(inl, l2_bec) = 2.d0*SUM(temp)
         IF (gstart == 2) becc%rvec(inl, l2_bec) = becc%rvec(inl, l2_bec) - temp(1)
      END DO
      DEALLOCATE (temp)
   ELSE
      ALLOCATE (temp_c(ngw))
      DO inl = 1, nhsavb
         DO ig = 1, ngw
            temp_c(ig) = CONJG(c(1, ig))*betae(ig, inl)
         END DO
         becc%cvec(inl, l2_bec) = SUM(temp_c)
      END DO
      DEALLOCATE (temp_c)
   END IF

   RETURN
END SUBROUTINE grabec_twin

!-------------------------------------------------------------------------
SUBROUTINE bec_csv_real(becc, becv, nkbx, csv, n)
!-----------------------------------------------------------------------
!     requires in input the updated becc and becv(k)
!     on output: csv is updated
!
   USE ions_base, ONLY: na
   USE cvan, ONLY: nvb, ish
   USE uspp, ONLY: qq
   USE uspp_param, ONLY: nh
   USE kinds, ONLY: DP
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: nkbx, n
   REAL(DP)    :: becc(nkbx)
   REAL(DP)    :: becv(nkbx, n)
   REAL(DP)    :: csv(n)
   INTEGER     :: k, is, iv, jv, ia, inl, jnl
   REAL(DP)    :: rsum

!     calculate csv(k) = csv(k) + <c| SUM_nm |beta(n)><beta(m)|v(k)>,  k<i
!
   DO k = 1, n
      rsum = 0.d0
      DO is = 1, nvb
         DO iv = 1, nh(is)
            DO jv = 1, nh(is)
               IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                  DO ia = 1, na(is)
                     inl = ish(is) + (iv - 1)*na(is) + ia
                     jnl = ish(is) + (jv - 1)*na(is) + ia
                     rsum = rsum + qq(iv, jv, is)*becc(inl)*becv(jnl, k)
                  END DO
               END IF
            END DO
         END DO
      END DO
      csv(k) = csv(k) + rsum
   END DO
!
   RETURN
END SUBROUTINE bec_csv_real
!-------------------------------------------------------------------------
SUBROUTINE bec_csv_twin(becc, becv, nkbx, csv, n, lcc)
!-----------------------------------------------------------------------
!     requires in input the updated becc and becv(k)
!     on output: csv is updated
!
   USE ions_base, ONLY: na
   USE cvan, ONLY: nvb, ish
   USE uspp, ONLY: qq
   USE uspp_param, ONLY: nh
   USE kinds, ONLY: DP
   USE twin_types
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: nkbx, n, lcc
   type(twin_matrix)    :: becc!( nkbx )
   type(twin_matrix)    :: becv!( nkbx, n )
   type(twin_matrix)    :: csv!( n )
   INTEGER     :: k, is, iv, jv, ia, inl, jnl
   REAL(DP)    :: rsum
   COMPLEX(DP) :: csum

!     calculate csv(k) = csv(k) + <c| SUM_nm |beta(n)><beta(m)|v(k)>,  k<i
!
   IF (.not. csv%iscmplx) THEN
      DO k = 1, n
         rsum = 0.d0
         DO is = 1, nvb
            DO iv = 1, nh(is)
               DO jv = 1, nh(is)
                  IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                     DO ia = 1, na(is)
                        inl = ish(is) + (iv - 1)*na(is) + ia
                        jnl = ish(is) + (jv - 1)*na(is) + ia
                        rsum = rsum + qq(iv, jv, is)*becc%rvec(inl, lcc)*becv%rvec(jnl, k)
                     END DO
                  END IF
               END DO
            END DO
         END DO
         csv%rvec(k, 1) = csv%rvec(k, 1) + rsum
      END DO
   ELSE
      DO k = 1, n
         csum = CMPLX(0.d0, 0.d0)
         DO is = 1, nvb
            DO iv = 1, nh(is)
               DO jv = 1, nh(is)
                  IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                     DO ia = 1, na(is)
                        inl = ish(is) + (iv - 1)*na(is) + ia
                        jnl = ish(is) + (jv - 1)*na(is) + ia
                        csum = csum + qq(iv, jv, is)*(becc%cvec(inl, lcc))*CONJG(becv%cvec(jnl, k))
                     END DO
                  END IF
               END DO
            END DO
         END DO
         csv%cvec(k, 1) = csv%cvec(k, 1) + csum
      END DO
   END IF
!
   RETURN
END SUBROUTINE bec_csv_twin

!$$
!----------------------------------------------------------------------
SUBROUTINE calc_wfnoverlap(a, b)
!----------------------------------------------------------------------

   use kinds
   use io_global, only: ionode
   use mp_global, only: intra_image_comm
   use gvecw, only: ngw
   use reciprocal_vectors, only: ng0 => gstart
   use mp, only: mp_sum
   use electrons_base, only: n => nbsp, nspin, nupdwn, iupdwn

   implicit none

   complex(dp) a(ngw, n), b(ngw, n)
   integer i, j, ig, isp, ndim, nbnd1, nbnd2
   real(dp) sca
   real(DP), allocatable :: s(:, :)
   !

   s(:, :) = 0.d0

   do isp = 1, nspin

      ndim = nupdwn(isp)

      allocate (s(ndim, ndim))

      s(:, :) = 0.d0

      do i = 1, ndim
         nbnd1 = iupdwn(isp) - 1 + i
         if (ng0 .eq. 2) a(1, nbnd1) = CMPLX(DBLE(a(1, nbnd1)), 0.0d0)
         if (ng0 .eq. 2) b(1, nbnd1) = CMPLX(DBLE(b(1, nbnd1)), 0.0d0)
      end do

      do i = 1, ndim
         nbnd1 = iupdwn(isp) - 1 + i
         do j = 1, ndim
            nbnd2 = iupdwn(isp) - 1 + j
            sca = 0.0d0
            do ig = 1, ngw           !loop on g vectors
               sca = sca + DBLE(CONJG(a(ig, nbnd1))*b(ig, nbnd2))
            end do
            sca = sca*2.0d0  !2. for real weavefunctions
            if (ng0 .eq. 2) sca = sca - DBLE(a(1, nbnd1))*DBLE(b(1, nbnd2))
            s(i, j) = sca
         end do
      end do
      call mp_sum(s, intra_image_comm)

      do i = 1, ndim
         do j = 1, ndim
            s(i, j) = s(i, j)*s(i, j)
         end do
      end do

      if (ionode) then
         write (1235, *) 'spin ', isp, ', sum ', sum(s(:, :))/ndim
         do i = 1, ndim
            write (1235, '(40f5.2)') (s(i, j), j=1, ndim)
         end do
         write (1235, *)
      end if

      deallocate (s)

   end do

END SUBROUTINE calc_wfnoverlap
!$$

!-------------------------------------------------------------------------
SUBROUTINE gram(betae, bec, nkbx, cp, ngwx, n)
!-----------------------------------------------------------------------
!     gram-schmidt orthogonalization of the set of wavefunctions cp
!
   USE uspp, ONLY: nkb
   USE gvecw, ONLY: ngw
   USE kinds, ONLY: DP
   USE control_flags, ONLY: gamma_only, do_wf_cmplx !added:giovanni
   USE twin_types !added:giovanni
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: nkbx, ngwx, n
   type(twin_matrix)   :: bec!( nkbx, n )!modified:giovanni
   COMPLEX(DP)   :: cp(ngwx, n), betae(ngwx, nkb)
!
   REAL(DP) :: anorm, cscnorm
   COMPLEX(DP), ALLOCATABLE :: csc(:) !modified:giovanni
   INTEGER :: i, k
   LOGICAL :: lgam !added:giovanni
   EXTERNAL cscnorm
   !
   lgam = gamma_only .and. .not. do_wf_cmplx !added:giovanni
   !
   CALL start_clock('gram')
!       write(6,*) bec%rvec !added:giovanni:debug
!       stop
   ALLOCATE (csc(n))
   !
   csc = CMPLX(0.d0, 0.d0)
   !
   DO i = 1, n
      !
      CALL gracsc(bec, nkbx, betae, cp, ngwx, i, csc, n, lgam)!added:giovanni lgam
      !
      ! calculate orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
      !
      DO k = 1, i - 1
         CALL ZAXPY(ngw, -csc(k), cp(1, k), 1, cp(1, i), 1)!modified:giovanni
      END DO
      anorm = cscnorm(bec, nkbx, cp, ngwx, i, n, lgam)
      CALL ZSCAL(ngw, CMPLX(1.0d0/anorm, 0.d0), cp(1, i), 1)
      !
      ! these are the final bec's
      !
      IF (nkbx > 0) THEN
         IF (.not. bec%iscmplx) THEN
            CALL DSCAL(nkbx, 1.0d0/anorm, bec%rvec(1:nkbx, i), 1)!modified:giovanni
         ELSE
            CALL ZSCAL(nkbx, CMPLX(1.0d0/anorm, 0.d0), bec%cvec(1:nkbx, i), 1)!added:giovanni
         END IF
      END IF
      !
   END DO
   !
!       write(6,*) "csc_giovanni_debug", csc !added:giovanni:debug
   DEALLOCATE (csc)

   CALL stop_clock('gram')
!
   RETURN
END SUBROUTINE gram
!
!-------------------------------------------------------------------------
SUBROUTINE gram2(betae, bec, nkbx, cp, ngwx, n)
!-----------------------------------------------------------------------
!     gram-schmidt orthogonalization of the set of wavefunctions cp
!
   USE uspp, ONLY: nkb
   USE gvecw, ONLY: ngw
   USE kinds, ONLY: DP
   USE control_flags, ONLY: gamma_only, do_wf_cmplx !added:giovanni
   USE twin_types !added:giovanni
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: nkbx, ngwx, n
   type(twin_matrix)   :: bec!( nkbx, n )!modified:giovanni
   COMPLEX(DP)   :: cp(ngwx, n), betae(ngwx, nkb)
!
   REAL(DP) :: anorm, cscnorm
   COMPLEX(DP), ALLOCATABLE :: csc(:) !modified:giovanni
   INTEGER :: i, k
   LOGICAL :: lgam !added:giovanni
   EXTERNAL cscnorm

   lgam = gamma_only .and. .not. do_wf_cmplx !added:giovanni

!
   CALL start_clock('gram')

   ALLOCATE (csc(n))
   csc = CMPLX(0.d0, 0.d0)
!
   DO i = 1, n
      DO k = 1, i - 1
         !
         CALL gracsc2(bec, nkbx, betae, cp, ngwx, i, k, csc, n, lgam)!added:giovanni lgam
         anorm = cscnorm(bec, nkbx, cp, ngwx, k, n, lgam)
         !
         ! calculate orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
         !
         CALL ZAXPY(ngw, -csc(k)*CMPLX(1/anorm, 0.d0), cp(1, k), 1, cp(1, i), 1)!modified:giovanni
      END DO
      anorm = cscnorm(bec, nkbx, cp, ngwx, i, n, lgam)
      CALL ZSCAL(ngw, CMPLX(1.0d0/anorm, 0.d0), cp(1, i), 1)
      !

      !
      !         these are the final bec's
      !
      IF (nkbx > 0) THEN
         IF (.not. bec%iscmplx) THEN
            CALL DSCAL(nkbx, 1.0d0/anorm, bec%rvec(1:nkbx, i), 1)!modified:giovanni
         ELSE
            CALL ZSCAL(nkbx, CMPLX(1.0d0/anorm, 0.d0), bec%cvec(1:nkbx, i), 1)!added:giovanni
         END IF
      END IF

   END DO
!
   Do i = 1, n
   DO k = 1, i
      write (6, *) dot_product(conjg(cp(:, i)), cp(:, k))
   END DO
   END DO
   DEALLOCATE (csc)

   CALL stop_clock('gram')
!
   RETURN
END SUBROUTINE gram2

!-----------------------------------------------------------------------
SUBROUTINE initbox(tau0, taub, irb, ainv, a1, a2, a3)
!-----------------------------------------------------------------------
!
!     sets the indexes irb and positions taub for the small boxes
!     around atoms
!
   USE kinds, ONLY: DP
   USE ions_base, ONLY: nsp, na, nat
   USE grid_dimensions, ONLY: nr1, nr2, nr3
   USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx
   USE control_flags, ONLY: iprsta
   USE io_global, ONLY: stdout
   USE mp_global, ONLY: nproc_image, me_image
   USE fft_base, ONLY: dfftb, dfftp, fft_dlay_descriptor
   USE fft_types, ONLY: fft_box_set

   IMPLICIT NONE
! input
   REAL(DP), INTENT(in)  :: tau0(3, nat)
! output
   INTEGER, INTENT(out) :: irb(3, nat)
   REAL(DP), INTENT(out) :: taub(3, nat)
! input
   REAL(DP), INTENT(in)  :: ainv(3, 3)
   REAL(DP), INTENT(in)  :: a1(3)
   REAL(DP), INTENT(in)  :: a2(3)
   REAL(DP), INTENT(in)  :: a3(3)
! local
   REAL(DP) :: x(3), xmod
   INTEGER  :: nr(3), nrb(3), xint, is, ia, i, isa
!
   IF (nr1b < 1) CALL errore &
      ('initbox', 'incorrect value for box grid dimensions', 1)
   IF (nr2b < 1) CALL errore &
      ('initbox', 'incorrect value for box grid dimensions', 2)
   IF (nr3b < 1) CALL errore &
      ('initbox', 'incorrect value for box grid dimensions', 3)

   nr(1) = nr1
   nr(2) = nr2
   nr(3) = nr3
   nrb(1) = nr1b
   nrb(2) = nr2b
   nrb(3) = nr3b
!
   isa = 0
   DO is = 1, nsp
      DO ia = 1, na(is)
         isa = isa + 1
!
         DO i = 1, 3
!
! bring atomic positions to crystal axis
!
            x(i) = ainv(i, 1)*tau0(1, isa) +                         &
  &                ainv(i, 2)*tau0(2, isa) +                         &
  &                ainv(i, 3)*tau0(3, isa)
!
! bring x in the range between 0 and 1
!
            x(i) = MOD(x(i), 1.d0)
            IF (x(i) .LT. 0.d0) x(i) = x(i) + 1.d0
!
! case of nrb(i) even
!
            IF (MOD(nrb(i), 2) .EQ. 0) THEN
!
! find irb = index of the grid point at the corner of the small box
!           (the indices of the small box run from irb to irb+nrb-1)
!
               xint = INT(x(i)*nr(i))
               irb(i, isa) = xint + 1 - nrb(i)/2 + 1
               IF (irb(i, isa) .LT. 1) irb(i, isa) = irb(i, isa) + nr(i)
!
! x(i) are the atomic positions in crystal coordinates, where the
! "crystal lattice" is the small box lattice and the origin is at
! the corner of the small box. Used to calculate phases exp(iG*taub)
!
               xmod = x(i)*nr(i) - xint
               x(i) = (xmod + nrb(i)/2 - 1)/nr(i)
            ELSE
!
! case of nrb(i) odd - see above for comments
!
               xint = NINT(x(i)*nr(i))
               irb(i, isa) = xint + 1 - (nrb(i) - 1)/2
               IF (irb(i, isa) .LT. 1) irb(i, isa) = irb(i, isa) + nr(i)
               xmod = x(i)*nr(i) - xint
               x(i) = (xmod + (nrb(i) - 1)/2)/nr(i)
            END IF
         END DO
!
! bring back taub in cartesian coordinates
!
         DO i = 1, 3
            taub(i, isa) = x(1)*a1(i) + x(2)*a2(i) + x(3)*a3(i)
         END DO
      END DO
   END DO

   CALL fft_box_set(dfftb, nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx, &
                    nat, irb, me_image + 1, nproc_image, dfftp%npp, dfftp%ipp)

   IF (iprsta > 2) THEN
      isa = 1
      DO is = 1, nsp
         WRITE (stdout, '( /, 2x, "species= ", i2 )') is
         DO ia = 1, na(is)
            WRITE (stdout, 2000) ia, (irb(i, isa), i=1, 3)
2000        FORMAT(2x, 'atom= ', i3, ' irb1= ', i3, ' irb2= ', i3, ' irb3= ', i3)
            isa = isa + 1
         END DO
      END DO
   END IF

#ifdef __PARA
   !
   ! for processor that do not call fft on the box
   ! artificially start the clock
   !
   CALL start_clock('fftb')
   CALL stop_clock('fftb')
   !
#endif
!
   RETURN
END SUBROUTINE initbox
!
!-------------------------------------------------------------------------
SUBROUTINE newd(vr, irb, eigrb, rhovan, fion)
!-----------------------------------------------------------------------
!
!     this routine calculates array deeq:
!         deeq_i,lm = \int V_eff(r) q_i,lm(r) dr
!     and the corresponding term in forces
!         fion_i = \int V_eff(r) \sum_lm rho_lm (dq_i,lm(r)/dR_i) dr
!     where
!         rho_lm = \sum_j f_j <psi_j|beta_l><beta_m|psi_j>
!
   USE kinds, ONLY: dp
   USE uspp_param, ONLY: nh, nhm
   USE uspp, ONLY: deeq
   USE cvan, ONLY: nvb
   USE ions_base, ONLY: nat, na
   USE constants, ONLY: pi, fpi
   USE grid_dimensions, ONLY: nnr => nnrx
   USE gvecb, ONLY: ngb, npb, nmb, gxb
   USE small_box, ONLY: omegab, tpibab
   USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, &
                                       nnrb => nnrbx
   USE qgb_mod, ONLY: qgb
   USE electrons_base, ONLY: nspin
   USE control_flags, ONLY: thdyn, tfor, tprnfor
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE cp_interfaces, ONLY: invfft
   USE fft_base, ONLY: dfftb
   USE control_flags, ONLY: gamma_only, do_wf_cmplx !added:giovanni
!
   IMPLICIT NONE
! input
   INTEGER irb(3, nat)
   REAL(DP) rhovan(nhm*(nhm + 1)/2, nat, nspin)
   COMPLEX(DP) eigrb(ngb, nat)
   REAL(DP) vr(nnr, nspin)
! output
   REAL(DP) fion(3, nat)
! local
   INTEGER isup, isdw, iss, iv, ijv, jv, ik, nfft, isa, ia, is, ig
   REAL(DP) fvan(3, nat, nvb), fac, fac1, fac2, boxdotgrid
   COMPLEX(DP) ci, facg1, facg2
   COMPLEX(DP), ALLOCATABLE :: qv(:)
   LOGICAL :: lgam !added:giovanni
   INTEGER :: istep !added:giovanni
   EXTERNAL boxdotgrid
!
   CALL start_clock('newd')
   lgam = gamma_only .and. .not. do_wf_cmplx
   ci = (0.d0, 1.d0)
   fac = omegab/DBLE(nr1b*nr2b*nr3b)
   deeq(:, :, :, :) = 0.d0
   fvan(:, :, :) = 0.d0

   ALLOCATE (qv(nnrb))
   IF (lgam) THEN
      istep = 2
   ELSE
      istep = 2
   END IF
!
! calculation of deeq_i,lm = \int V_eff(r) q_i,lm(r) dr
!
   isa = 1
   DO is = 1, nvb
#ifdef __PARA
      nfft = 1
#else
      nfft = istep
#endif
      DO ia = 1, na(is), nfft
         nfft = 1
#ifdef __PARA
         IF (dfftb%np3(isa) <= 0) go to 15
#endif
         IF (ia .EQ. na(is)) nfft = 1
!
! two ffts at the same time, on two atoms (if possible: nfft=2)
!
         DO iv = 1, nh(is)
            DO jv = iv, nh(is)
               ijv = (jv - 1)*jv/2 + iv
               qv(:) = (0.d0, 0.d0)
               IF (nfft .EQ. 2) THEN
                  DO ig = 1, ngb
                     qv(npb(ig)) = eigrb(ig, isa)*qgb(ig, ijv, is)   &
  &                          + ci*eigrb(ig, isa + 1)*qgb(ig, ijv, is)
                     qv(nmb(ig)) = CONJG(                             &
  &                               eigrb(ig, isa)*qgb(ig, ijv, is))  &
  &                          + ci*CONJG(                             &
  &                               eigrb(ig, isa + 1)*qgb(ig, ijv, is))
                  END DO
               ELSE
!                      IF(lgam) THEN
                  DO ig = 1, ngb
                     qv(npb(ig)) = eigrb(ig, isa)*qgb(ig, ijv, is)
                     qv(nmb(ig)) = CONJG(                            &
 &                                eigrb(ig, isa)*qgb(ig, ijv, is))
                  END DO
!                     ELSE
!                      DO ig=1,ngb
!                         qv(npb(ig)) = eigrb(ig,isa)*qgb(ig,ijv,is)
!                           qv(nmb(ig)) = CONJG(                            &
!       &                                eigrb(ig,isa)*qgb(ig,ijv,is))
!                      END DO
!                     ENDIF
               END IF
!
               CALL invfft('Box', qv, dfftb, isa)
!
               DO iss = 1, nspin
                  deeq(iv, jv, isa, iss) = fac*                        &
  &                    boxdotgrid(irb(1, isa), 1, qv, vr(1, iss))
                  IF (iv .NE. jv)                                      &
  &                    deeq(jv, iv, isa, iss) = deeq(iv, jv, isa, iss)
!
                  IF (nfft .EQ. 2) THEN
                     deeq(iv, jv, isa + 1, iss) = fac*                    &
  &                       boxdotgrid(irb(1, isa + 1), 2, qv, vr(1, iss))
                     IF (iv .NE. jv)                                   &
  &                       deeq(jv, iv, isa + 1, iss) = deeq(iv, jv, isa + 1, iss)
                  END IF
               END DO
            END DO
         END DO
15       isa = isa + nfft
      END DO
   END DO

   CALL mp_sum(deeq, intra_image_comm)

   IF (.NOT. (tfor .OR. thdyn .OR. tprnfor)) go to 10
!
! calculation of fion_i = \int V_eff(r) \sum_lm rho_lm (dq_i,lm(r)/dR_i) dr
!
   isa = 1
   IF (nspin .EQ. 1) THEN
!     =================================================================
!     case nspin=1: two ffts at the same time, on two atoms (if possible)
!     -----------------------------------------------------------------
      iss = 1
      isa = 1
      DO is = 1, nvb
#ifdef __PARA
         nfft = 1
#else
         nfft = istep
#endif
         DO ia = 1, na(is), istep
            nfft = 1
#ifdef __PARA
            IF (dfftb%np3(isa) <= 0) go to 20
#endif
            IF (ia .EQ. na(is)) nfft = 1
            DO ik = 1, 3
               qv(:) = (0.d0, 0.d0)
               DO iv = 1, nh(is)
                  DO jv = iv, nh(is)
                     ijv = (jv - 1)*jv/2 + iv
                     IF (iv .NE. jv) THEN
                        fac1 = 2.d0*fac*tpibab*rhovan(ijv, isa, iss)
                        IF (nfft .EQ. 2) fac2 = 2.d0*fac*tpibab*         &
  &                                           rhovan(ijv, isa + 1, iss)
                     ELSE
                        fac1 = fac*tpibab*rhovan(ijv, isa, iss)
                        IF (nfft .EQ. 2) fac2 = fac*tpibab*        &
  &                                           rhovan(ijv, isa + 1, iss)
                     END IF
                     IF (nfft .EQ. 2) THEN
                        DO ig = 1, ngb
                           facg1 = CMPLX(0.d0, -gxb(ik, ig))*         &
  &                                   qgb(ig, ijv, is)*fac1
                           facg2 = CMPLX(0.d0, -gxb(ik, ig))*         &
  &                                   qgb(ig, ijv, is)*fac2
                           qv(npb(ig)) = qv(npb(ig))                 &
  &                                    + eigrb(ig, isa)*facg1  &
  &                                    + ci*eigrb(ig, isa + 1)*facg2
                           qv(nmb(ig)) = qv(nmb(ig))                 &
  &                                + CONJG(eigrb(ig, isa)*facg1)&
  &                                + ci*CONJG(eigrb(ig, isa + 1)*facg2)
                        END DO
                     ELSE
!                            IF(lgam) THEN
                        DO ig = 1, ngb
                           facg1 = CMPLX(0.d0, -gxb(ik, ig))*         &
 &                                   qgb(ig, ijv, is)*fac1
                           qv(npb(ig)) = qv(npb(ig))                 &
 &                                    + eigrb(ig, isa)*facg1
                           qv(nmb(ig)) = qv(nmb(ig))                 &
 &                               + CONJG(eigrb(ig, isa)*facg1)
                        END DO
!                            ELSE
!                                 facg1 = CMPLX(0.d0,-gxb(ik,ig)) *         &
!       &                                   qgb(ig,ijv,is)*fac1
!                                 qv(npb(ig)) = qv(npb(ig))                 &
!       &                                    +    eigrb(ig,isa)*facg1
!                                 qv(nmb(ig)) = qv(nmb(ig))                 &
!       &                               +  CONJG( eigrb(ig,isa)*facg1)
!                            ENDIF
                     END IF
                  END DO
               END DO
!
               CALL invfft('Box', qv, dfftb, isa)
!
               fvan(ik, ia, is) =                                      &
  &                    boxdotgrid(irb(1, isa), 1, qv, vr(1, iss))
!
               IF (nfft .EQ. 2) fvan(ik, ia + 1, is) =                     &
  &                    boxdotgrid(irb(1, isa + 1), 2, qv, vr(1, iss))
            END DO
20          isa = isa + nfft
         END DO
      END DO
   ELSE
!     =================================================================
!     case nspin=2: up and down spin fft's combined into a single fft
!     -----------------------------------------------------------------
      isup = 1
      isdw = 2
      isa = 1
      DO is = 1, nvb
         DO ia = 1, na(is)
#ifdef __PARA
            IF (dfftb%np3(isa) <= 0) go to 25
#endif
            DO ik = 1, 3
               qv(:) = (0.d0, 0.d0)
!
               DO iv = 1, nh(is)
                  DO jv = iv, nh(is)
                     ijv = (jv - 1)*jv/2 + iv
                     IF (iv .NE. jv) THEN
                        fac1 = 2.d0*fac*tpibab*rhovan(ijv, isa, isup)
                        fac2 = 2.d0*fac*tpibab*rhovan(ijv, isa, isdw)
                     ELSE
                        fac1 = fac*tpibab*rhovan(ijv, isa, isup)
                        fac2 = fac*tpibab*rhovan(ijv, isa, isdw)
                     END IF
!                         IF(lgam) THEN
                     DO ig = 1, ngb
                        facg1 = fac1*CMPLX(0.d0, -gxb(ik, ig))*     &
  &                                qgb(ig, ijv, is)*eigrb(ig, isa)
                        facg2 = fac2*CMPLX(0.d0, -gxb(ik, ig))*     &
  &                                qgb(ig, ijv, is)*eigrb(ig, isa)
                        qv(npb(ig)) = qv(npb(ig))                    &
  &                                    + facg1 + ci*facg2
                        qv(nmb(ig)) = qv(nmb(ig))                    &
  &                                    + CONJG(facg1) + ci*CONJG(facg2)
                     END DO
!                         ELSE
!                           DO ig=1,ngb
!                             facg1 = fac1 * CMPLX(0.d0,-gxb(ik,ig)) *     &
!       &                                qgb(ig,ijv,is) * eigrb(ig,isa)
!                             facg2 = fac2 * CMPLX(0.d0,-gxb(ik,ig)) *     &
!       &                                qgb(ig,ijv,is) * eigrb(ig,isa)
!                             qv(npb(ig)) = qv(npb(ig))                    &
!       &                                    + facg1 + ci*facg2
!                            qv(nmb(ig)) = qv(nmb(ig))                    &
!      &                                    +CONJG(facg1)+ci*CONJG(facg2)
!                           END DO
!                         ENDIF
                  END DO
               END DO
!
               CALL invfft('Box', qv, dfftb, isa)
!
               fvan(ik, ia, is) =                                      &
  &                    boxdotgrid(irb(1, isa), isup, qv, vr(1, isup)) + &
  &                    boxdotgrid(irb(1, isa), isdw, qv, vr(1, isdw))
            END DO
25          isa = isa + 1
         END DO
      END DO
   END IF

   CALL mp_sum(fvan, intra_image_comm)

   isa = 0
   DO is = 1, nvb
      DO ia = 1, na(is)
         isa = isa + 1
         fion(:, isa) = fion(:, isa) - fvan(:, ia, is)
      END DO
   END DO

10 CONTINUE
   DEALLOCATE (qv)
!
   CALL stop_clock('newd')
!
   RETURN
END SUBROUTINE newd

!-------------------------------------------------------------------------
SUBROUTINE nlfl_real(bec, becdr, lambda, fion)
!-----------------------------------------------------------------------
!     contribution to fion due to the orthonormality constraint
!
!

   USE kinds, ONLY: DP
   USE ions_base, ONLY: na, nat
   USE uspp, ONLY: nhsa => nkb, qq
   USE uspp_param, ONLY: nhm, nh
   USE cvan, ONLY: ish, nvb
   USE electrons_base, ONLY: nbsp, nspin, iupdwn, nupdwn
   USE constants, ONLY: pi, fpi
   USE cp_main_variables, ONLY: nlam, nlax, descla, la_proc
   USE descriptors, ONLY: nlar_, nlac_, ilar_, ilac_
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
!
   IMPLICIT NONE
   REAL(DP) bec(nhsa, nbsp), becdr(nhsa, nspin*nlax, 3), lambda(nlam, nlam, nspin)
   REAL(DP) fion(3, nat)
!
   INTEGER :: k, is, ia, iv, jv, i, j, inl, isa, iss, nss, istart, ir, ic, nr, nc
   REAL(DP), ALLOCATABLE :: temp(:, :), tmpbec(:, :), tmpdr(:, :)
   REAL(DP), ALLOCATABLE :: fion_tmp(:, :)
   !
   CALL start_clock('nlfl')
   !
   ALLOCATE (fion_tmp(3, nat))
   !
   fion_tmp = 0.0d0
   !

   ALLOCATE (temp(nlax, nlax), tmpbec(nhm, nlax), tmpdr(nlax, nhm))

   DO k = 1, 3
      isa = 0
      DO is = 1, nvb
         DO ia = 1, na(is)
            isa = isa + 1
            !
            DO iss = 1, nspin
               !
               nss = nupdwn(iss)
               istart = iupdwn(iss)
               !
               tmpbec = 0.d0
               tmpdr = 0.d0
!
               IF (la_proc) THEN
                  ! tmpbec distributed by columns
                  ic = descla(ilac_, iss)
                  nc = descla(nlac_, iss)
                  DO iv = 1, nh(is)
                     DO jv = 1, nh(is)
                        inl = ish(is) + (jv - 1)*na(is) + ia
                        IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                           DO i = 1, nc
                              tmpbec(iv, i) = tmpbec(iv, i) + qq(iv, jv, is)*bec(inl, i + istart - 1 + ic - 1)
                           END DO
                        END IF
                     END DO
                  END DO
                  ! tmpdr distributed by rows
                  ir = descla(ilar_, iss)
                  nr = descla(nlar_, iss)
                  DO iv = 1, nh(is)
                     inl = ish(is) + (iv - 1)*na(is) + ia
                     DO i = 1, nr
                        tmpdr(i, iv) = becdr(inl, i + (iss - 1)*nlax, k)
                     END DO
                  END DO
               END IF
!
               IF (nh(is) .GT. 0) THEN
                  !
                  IF (la_proc) THEN
                     ir = descla(ilar_, iss)
                     ic = descla(ilac_, iss)
                     nr = descla(nlar_, iss)
                     nc = descla(nlac_, iss)
                     CALL DGEMM('N', 'N', nr, nc, nh(is), 1.0d0, tmpdr, nlax, tmpbec, nhm, 0.0d0, temp, nlax)
                     DO j = 1, nc
                        DO i = 1, nr
                           fion_tmp(k, isa) = fion_tmp(k, isa) + 2D0*temp(i, j)*lambda(i, j, iss)
                        END DO
                     END DO
                  END IF
!
               END IF

            END DO
!
         END DO
      END DO
   END DO
   !
   DEALLOCATE (temp, tmpbec, tmpdr)
   !
   CALL mp_sum(fion_tmp, intra_image_comm)
   !
   fion = fion + fion_tmp
   !
   DEALLOCATE (fion_tmp)
   !
   CALL stop_clock('nlfl')
   !
   RETURN

END SUBROUTINE nlfl_real

!-------------------------------------------------------------------------
SUBROUTINE nlfl_cmplx(bec, becdr, lambda, fion)
!-----------------------------------------------------------------------
!     contribution to fion due to the orthonormality constraint
!
!

   USE kinds, ONLY: DP
   USE ions_base, ONLY: na, nat
   USE uspp, ONLY: nhsa => nkb, qq
   USE uspp_param, ONLY: nhm, nh
   USE cvan, ONLY: ish, nvb
   USE electrons_base, ONLY: nbsp, nspin, iupdwn, nupdwn
   USE constants, ONLY: pi, fpi
   USE cp_main_variables, ONLY: nlam, nlax, descla, la_proc
   USE descriptors, ONLY: nlar_, nlac_, ilar_, ilac_
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
!
   IMPLICIT NONE
   COMPLEX(DP) bec(nhsa, nbsp), becdr(nhsa, nspin*nlax, 3), lambda(nlam, nlam, nspin)
   REAL(DP) fion(3, nat)
!
   INTEGER :: k, is, ia, iv, jv, i, j, inl, isa, iss, nss, istart, ir, ic, nr, nc
   COMPLEX(DP), ALLOCATABLE :: temp(:, :), tmpbec(:, :), tmpdr(:, :)
   REAL(DP), ALLOCATABLE :: fion_tmp(:, :)
   !
   CALL start_clock('nlfl')
   !
   ALLOCATE (fion_tmp(3, nat))
   !
   fion_tmp = 0.0d0
   !

   ALLOCATE (temp(nlax, nlax), tmpbec(nhm, nlax), tmpdr(nlax, nhm))

   DO k = 1, 3
      isa = 0
      DO is = 1, nvb
         DO ia = 1, na(is)
            isa = isa + 1
            !
            DO iss = 1, nspin
               !
               nss = nupdwn(iss)
               istart = iupdwn(iss)
               !
               tmpbec = 0.d0
               tmpdr = 0.d0
!
               IF (la_proc) THEN
                  ! tmpbec distributed by columns
                  ic = descla(ilac_, iss)
                  nc = descla(nlac_, iss)
                  DO iv = 1, nh(is)
                     DO jv = 1, nh(is)
                        inl = ish(is) + (jv - 1)*na(is) + ia
                        IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                           DO i = 1, nc
                              tmpbec(iv, i) = tmpbec(iv, i) + qq(iv, jv, is)*bec(inl, i + istart - 1 + ic - 1)
                           END DO
                        END IF
                     END DO
                  END DO
                  ! tmpdr distributed by rows
                  ir = descla(ilar_, iss)
                  nr = descla(nlar_, iss)
                  DO iv = 1, nh(is)
                     inl = ish(is) + (iv - 1)*na(is) + ia
                     DO i = 1, nr
                        tmpdr(i, iv) = becdr(inl, i + (iss - 1)*nlax, k)
                     END DO
                  END DO
               END IF
!
               IF (nh(is) .GT. 0) THEN
                  !
                  IF (la_proc) THEN
                     ir = descla(ilar_, iss)
                     ic = descla(ilac_, iss)
                     nr = descla(nlar_, iss)
                     nc = descla(nlac_, iss)
                     CALL ZGEMM('N', 'N', nr, nc, nh(is), (1.0d0, 0.d0), tmpdr, nlax, tmpbec, nhm, (0.0d0, 0.d0), temp, nlax)
                     DO j = 1, nc
                        DO i = 1, nr
                           fion_tmp(k, isa) = fion_tmp(k, isa) + DBLE(2D0*temp(i, j)*lambda(i, j, iss))
                        END DO
                     END DO
                  END IF
!
               END IF

            END DO
!
         END DO
      END DO
   END DO
   !
   DEALLOCATE (temp, tmpbec, tmpdr)
   !
   CALL mp_sum(fion_tmp, intra_image_comm)
   !
   fion = fion + fion_tmp
   !
   DEALLOCATE (fion_tmp)
   !
   CALL stop_clock('nlfl')
   !
   RETURN

END SUBROUTINE nlfl_cmplx

!-------------------------------------------------------------------------
SUBROUTINE nlfl_twin(bec, becdr, lambda, fion, lgam)
!-----------------------------------------------------------------------
!     contribution to fion due to the orthonormality constraint
!
!

   USE kinds, ONLY: DP
   USE ions_base, ONLY: na, nat
   USE uspp, ONLY: qq
   USE uspp_param, ONLY: nhm, nh
   USE cvan, ONLY: ish, nvb
   USE electrons_base, ONLY: nspin, iupdwn, nupdwn
   USE constants, ONLY: pi, fpi
   USE cp_main_variables, ONLY: nlax, descla, la_proc
   USE descriptors, ONLY: nlar_, nlac_, ilar_, ilac_
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE twin_types !added:giovanni
!
   IMPLICIT NONE

   type(twin_matrix) :: lambda(nspin)
   REAL(DP) :: fion(3, nat)
   TYPE(twin_matrix) :: bec
   TYPE(twin_tensor) :: becdr
   LOGICAL :: lgam
!
   INTEGER :: k, is, ia, iv, jv, i, j, inl, isa, iss, nss, istart, ir, ic, nr, nc
   REAL(DP), ALLOCATABLE :: temp(:, :), tmpbec(:, :), tmpdr(:, :)
   COMPLEX(DP), ALLOCATABLE :: temp_c(:, :), tmpbec_c(:, :), tmpdr_c(:, :)
   REAL(DP), ALLOCATABLE :: fion_tmp(:, :)
   COMPLEX(DP), PARAMETER :: c_one = CMPLX(1.d0, 0.d0), c_zero = CMPLX(0.d0, 0.d0)
   !
   CALL start_clock('nlfl')
   !
   ALLOCATE (fion_tmp(3, nat))
   !
   fion_tmp = 0.0d0
   !
   IF (lgam) THEN
      ALLOCATE (temp(nlax, nlax), tmpbec(nhm, nlax), tmpdr(nlax, nhm))
   ELSE
      ALLOCATE (temp_c(nlax, nlax), tmpbec_c(nhm, nlax), tmpdr_c(nlax, nhm))
   END IF

   IF (lgam) THEN
      DO k = 1, 3
         isa = 0
         DO is = 1, nvb
            DO ia = 1, na(is)
               isa = isa + 1
               !
               DO iss = 1, nspin
                  !
                  nss = nupdwn(iss)
                  IF (nss > 0) THEN
                     istart = iupdwn(iss)
                     !
                     tmpbec = 0.d0
                     tmpdr = 0.d0
                     !
                     IF (la_proc) THEN
                        ! tmpbec distributed by columns
                        ic = descla(ilac_, iss)
                        nc = descla(nlac_, iss)
                        DO iv = 1, nh(is)
                           DO jv = 1, nh(is)
                              inl = ish(is) + (jv - 1)*na(is) + ia
                              IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                                 DO i = 1, nc
                                    tmpbec(iv, i) = tmpbec(iv, i) + qq(iv, jv, is)*bec%rvec(inl, i + istart - 1 + ic - 1)
                                 END DO
                              END IF
                           END DO
                        END DO
                        ! tmpdr distributed by rows
                        ir = descla(ilar_, iss)
                        nr = descla(nlar_, iss)
                        DO iv = 1, nh(is)
                           inl = ish(is) + (iv - 1)*na(is) + ia
                           DO i = 1, nr
                              tmpdr(i, iv) = becdr%rvec(inl, i + (iss - 1)*nlax, k)
                           END DO
                        END DO
                     END IF
                     !
                     IF (nh(is) .GT. 0) THEN
                        !
                        IF (la_proc) THEN
                           ir = descla(ilar_, iss)
                           ic = descla(ilac_, iss)
                           nr = descla(nlar_, iss)
                           nc = descla(nlac_, iss)
                           CALL DGEMM('N', 'N', nr, nc, nh(is), 1.0d0, tmpdr, nlax, tmpbec, nhm, 0.0d0, temp, nlax)
                           DO j = 1, nc
                              DO i = 1, nr
                                 fion_tmp(k, isa) = fion_tmp(k, isa) + 2D0*temp(i, j)*lambda(iss)%rvec(i, j)
                              END DO
                           END DO
                        END IF
                        !
                     END IF
                  END IF
               END DO
               !
            END DO
         END DO
      END DO
   ELSE
      DO k = 1, 3
         isa = 0
         DO is = 1, nvb
            DO ia = 1, na(is)
               isa = isa + 1
               !
               DO iss = 1, nspin
                  !
                  nss = nupdwn(iss)
                  istart = iupdwn(iss)
                  !
                  tmpbec_c = CMPLX(0.d0, 0.d0)
                  tmpdr_c = CMPLX(0.d0, 0.d0)
                  !
                  IF (la_proc) THEN
                     ! tmpbec distributed by columns
                     ic = descla(ilac_, iss)
                     nc = descla(nlac_, iss)
                     DO iv = 1, nh(is)
                        DO jv = 1, nh(is)
                           inl = ish(is) + (jv - 1)*na(is) + ia
                           IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                              DO i = 1, nc
                                 tmpbec_c(iv, i) = tmpbec_c(iv, i) + qq(iv, jv, is)*bec%cvec(inl, i + istart - 1 + ic - 1)
                              END DO
                           END IF
                        END DO
                     END DO
                     ! tmpdr distributed by rows
                     ir = descla(ilar_, iss)
                     nr = descla(nlar_, iss)
                     DO iv = 1, nh(is)
                        inl = ish(is) + (iv - 1)*na(is) + ia
                        DO i = 1, nr
                           tmpdr_c(i, iv) = becdr%cvec(inl, i + (iss - 1)*nlax, k)
                        END DO
                     END DO
                  END IF
                  !
                  IF (nh(is) .GT. 0) THEN
                     !
                     IF (la_proc) THEN
                        ir = descla(ilar_, iss)
                        ic = descla(ilac_, iss)
                        nr = descla(nlar_, iss)
                        nc = descla(nlac_, iss)
                        CALL ZGEMM('N', 'N', nr, nc, nh(is), c_one, tmpdr_c, nlax, tmpbec_c, nhm, c_zero, temp_c, nlax) !warning:giovanni:check C
                        DO j = 1, nc
                           DO i = 1, nr
                              fion_tmp(k, isa) = fion_tmp(k, isa) + 2D0*DBLE(CONJG(temp_c(i, j))*lambda(iss)%cvec(i, j))
                           END DO
                        END DO
                     END IF
                     !
                  END IF

               END DO
               !
            END DO
         END DO
      END DO
   END IF
   !
   IF (lgam) THEN
      DEALLOCATE (temp, tmpbec, tmpdr)
   ELSE
      DEALLOCATE (temp_c, tmpbec_c, tmpdr_c)
   END IF
   !
   CALL mp_sum(fion_tmp, intra_image_comm)
   !
   fion = fion + fion_tmp
   !
   DEALLOCATE (fion_tmp)
   !
   CALL stop_clock('nlfl')
   !
   RETURN

END SUBROUTINE nlfl_twin
!
!
!-----------------------------------------------------------------------
SUBROUTINE pbc(rin, a1, a2, a3, ainv, rout)
!-----------------------------------------------------------------------
!
!     brings atoms inside the unit cell
!
   USE kinds, ONLY: DP

   IMPLICIT NONE
! input
   REAL(DP) rin(3), a1(3), a2(3), a3(3), ainv(3, 3)
! output
   REAL(DP) rout(3)
! local
   REAL(DP) x, y, z
!
! bring atomic positions to crystal axis
!
   x = ainv(1, 1)*rin(1) + ainv(1, 2)*rin(2) + ainv(1, 3)*rin(3)
   y = ainv(2, 1)*rin(1) + ainv(2, 2)*rin(2) + ainv(2, 3)*rin(3)
   z = ainv(3, 1)*rin(1) + ainv(3, 2)*rin(2) + ainv(3, 3)*rin(3)
!
! bring x,y,z in the range between -0.5 and 0.5
!
   x = x - NINT(x)
   y = y - NINT(y)
   z = z - NINT(z)
!
! bring atomic positions back in cartesian axis
!
   rout(1) = x*a1(1) + y*a2(1) + z*a3(1)
   rout(2) = x*a1(2) + y*a2(2) + z*a3(2)
   rout(3) = x*a1(3) + y*a2(3) + z*a3(3)
!
   RETURN
END SUBROUTINE pbc

!
!-------------------------------------------------------------------------
SUBROUTINE prefor(eigr, betae)
!-----------------------------------------------------------------------
!
!     input :        eigr =  e^-ig.r_i
!     output:        betae_i,i(g) = (-i)**l beta_i,i(g) e^-ig.r_i
!
   USE kinds, ONLY: DP
   USE ions_base, ONLY: nsp, na, nat
   USE gvecw, ONLY: ngw
   USE cvan, ONLY: ish
   USE uspp, ONLY: nkb, beta, nhtol
   USE uspp_param, ONLY: nh
!
   IMPLICIT NONE
   COMPLEX(DP) :: eigr(ngw, nat)
   COMPLEX(DP) :: betae(ngw, nkb)
!
   INTEGER     :: is, iv, ia, inl, ig, isa
   COMPLEX(DP) :: ci
!
   CALL start_clock('prefor')
   isa = 0
   DO is = 1, nsp
      DO iv = 1, nh(is)
         ci = (0.0d0, -1.0d0)**nhtol(iv, is)
         DO ia = 1, na(is)
            inl = ish(is) + (iv - 1)*na(is) + ia
            DO ig = 1, ngw
               betae(ig, inl) = ci*beta(ig, iv, is)*eigr(ig, ia + isa)
            END DO
         END DO
      END DO
      isa = isa + na(is)
   END DO
   CALL stop_clock('prefor')
!
   RETURN
END SUBROUTINE prefor
!
!-----------------------------------------------------------------------
SUBROUTINE projwfc(c, nx, eigr, betae, n, ei)
!-----------------------------------------------------------------------
   !
   ! Projection on atomic wavefunctions
   !
   USE kinds, ONLY: DP
   USE constants, ONLY: autoev
   USE io_global, ONLY: stdout
   USE mp_global, ONLY: intra_image_comm
   USE mp, ONLY: mp_sum
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: gstart
   USE ions_base, ONLY: nsp, na, nat
   USE uspp, ONLY: nhsa => nkb
   USE uspp_param, ONLY: upf
   USE twin_types !added:giovanni
   USE cp_interfaces, ONLY: nlsm1, s_wfc !added:giovanni
!
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nx, n
   COMPLEX(DP), INTENT(IN) :: c(ngw, nx), eigr(ngw, nat), betae(ngw, nhsa)
   REAL(DP), INTENT(IN) :: ei(nx)
!
   COMPLEX(DP), ALLOCATABLE :: wfc(:, :), swfc(:, :)
   REAL(DP), ALLOCATABLE :: becwfc(:, :)
   REAL(DP), ALLOCATABLE :: overlap(:, :), e(:), z(:, :)
   REAL(DP), ALLOCATABLE :: proj(:, :), temp(:)
   REAL(DP)                 :: somma

   INTEGER :: n_atomic_wfc
   INTEGER :: is, ia, nb, l, m, k, i

   INTERFACE nlsm1_local
      subroutine nlsm1_real(n, nspmn, nspmx, eigr, c, becp) !addded:giovanni
         USE kinds, ONLY: DP

         implicit none

         integer, intent(in)  :: n, nspmn, nspmx
         complex(DP), intent(in)  :: eigr(:, :), c(:, :)
         REAL(DP) :: becp(:, :)
      end subroutine
   END INTERFACE
   !
   ! calculate number of atomic states
   !
   n_atomic_wfc = 0
   !
   DO is = 1, nsp
      DO nb = 1, upf(is)%nwfc
         l = upf(is)%lchi(nb)
         n_atomic_wfc = n_atomic_wfc + (2*l + 1)*na(is)
      END DO
   END DO
   IF (n_atomic_wfc .EQ. 0) RETURN
   !
   ALLOCATE (wfc(ngw, n_atomic_wfc))
   !
   ! calculate wfc = atomic states
   !
   CALL atomic_wfc(eigr, n_atomic_wfc, wfc)
   !
   ! calculate bec = <beta|wfc>
   !
   ALLOCATE (becwfc(nhsa, n_atomic_wfc))
   !
   CALL nlsm1_local(n_atomic_wfc, 1, nsp, eigr, wfc, becwfc)
   !
   ! calculate swfc = S|wfc>
   !
   ALLOCATE (swfc(ngw, n_atomic_wfc))
   !
   CALL s_wfc(n_atomic_wfc, becwfc, betae, wfc, swfc)
   !
   ! calculate overlap(i,j) = <wfc_i|S|wfc_j>
   !
   ALLOCATE (overlap(n_atomic_wfc, n_atomic_wfc))

   CALL DGEMM &
      ('T', 'N', n_atomic_wfc, n_atomic_wfc, 2*ngw, 1.0d0, wfc, 2*ngw, &
       swfc, 2*ngw, 0.0d0, overlap, n_atomic_wfc)

   CALL mp_sum(overlap, intra_image_comm)

   overlap = overlap*2.d0
   IF (gstart == 2) THEN
      DO l = 1, n_atomic_wfc
         DO m = 1, n_atomic_wfc
            overlap(m, l) = overlap(m, l) - DBLE(wfc(1, m))*DBLE(swfc(1, l))
         END DO
      END DO
   END IF
   !
   ! calculate (overlap)^(-1/2)(i,j). An orthonormal set of vectors |wfc_i>
   ! is obtained by introducing |wfc_j>=(overlap)^(-1/2)(i,j)*S|wfc_i>
   !
   ALLOCATE (z(n_atomic_wfc, n_atomic_wfc))
   ALLOCATE (e(n_atomic_wfc))
   !
   CALL rdiag(n_atomic_wfc, overlap, n_atomic_wfc, e, z)
   !
   overlap = 0.d0
   !
   DO l = 1, n_atomic_wfc
      DO m = 1, n_atomic_wfc
         DO k = 1, n_atomic_wfc
            overlap(l, m) = overlap(l, m) + z(m, k)*z(l, k)/SQRT(e(k))
         END DO
      END DO
   END DO
   !
   DEALLOCATE (e)
   DEALLOCATE (z)
   !
   ! calculate |wfc_j>=(overlap)^(-1/2)(i,j)*S|wfc_i>   (note the S matrix!)
   !
   wfc = 0.d0
   DO m = 1, n_atomic_wfc
      DO l = 1, n_atomic_wfc
         wfc(:, m) = wfc(:, m) + overlap(l, m)*swfc(:, l)
      END DO
   END DO
   DEALLOCATE (overlap)
   DEALLOCATE (swfc)
   DEALLOCATE (becwfc)
   !
   ! calculate proj = <c|S|wfc>
   !
   ALLOCATE (proj(n, n_atomic_wfc))

   ALLOCATE (temp(ngw))

   DO m = 1, n
      DO l = 1, n_atomic_wfc
         temp(:) = DBLE(CONJG(c(:, m))*wfc(:, l))
         proj(m, l) = 2.d0*SUM(temp)
         IF (gstart == 2) proj(m, l) = proj(m, l) - temp(1)
      END DO
   END DO

   DEALLOCATE (temp)

   CALL mp_sum(proj, intra_image_comm)

   i = 0
   WRITE (stdout, 90)
   WRITE (stdout, 100)
   DO is = 1, nsp
      DO nb = 1, upf(is)%nwfc
         l = upf(is)%lchi(nb)
         DO m = -l, l
            DO ia = 1, na(is)
               i = i + 1
            END DO
            WRITE (stdout, 110) i - na(is) + 1, i, na(is), is, nb, l, m
         END DO
      END DO
   END DO

   WRITE (stdout, *)
   DO m = 1, n
      somma = 0.d0
      DO l = 1, n_atomic_wfc
         somma = somma + proj(m, l)**2
      END DO
      WRITE (stdout, 120) m, somma, ei(m)*autoev
      WRITE (stdout, 130) (ABS(proj(m, l)), l=1, n_atomic_wfc)
   END DO

90 FORMAT(3X, 'Projection on atomic states')
100 FORMAT(3X, 'atomic state    atom   specie  wfc  l  m')
110 FORMAT(3X, I4, ' - ', I4, 4X, I4, 6X, I3, I5, I4, I3)
120 FORMAT(3X, 'state # ', i4, '    sum c^2 = ', f7.4, ' eV = ', F7.2)
130 FORMAT(3X, 10f7.4)

!
   DEALLOCATE (proj)
   DEALLOCATE (wfc)
   RETURN
END SUBROUTINE projwfc

!
!-----------------------------------------------------------------------
SUBROUTINE rdiag(n, h, ldh, e, v)
!-----------------------------------------------------------------------
!
!   calculates all the eigenvalues and eigenvectors of a complex
!   hermitean matrix H . On output, the matrix H is destroyed
!
   USE kinds, ONLY: DP
   USE dspev_module, ONLY: dspev_drv
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(in)   :: n, ldh
   REAL(DP), INTENT(inout):: h(ldh, n)
   REAL(DP), INTENT(out)  :: e(n)
   REAL(DP), INTENT(out)  :: v(ldh, n)
!
   INTEGER :: i, j, k
   REAL(DP), ALLOCATABLE :: ap(:)
!
   ALLOCATE (ap(n*(n + 1)/2))

   K = 0
   DO J = 1, n
      DO I = J, n
         K = K + 1
         ap(k) = h(i, j)
      END DO
   END DO

   CALL dspev_drv('V', 'L', n, ap, e, v, ldh)

   DEALLOCATE (ap)
!
   RETURN
END SUBROUTINE rdiag

!
!
!-------------------------------------------------------------------------
SUBROUTINE s_wfc_real(n_atomic_wfc1, becwfc, betae, wfc, swfc) !@@@@ Changed n_atomic_wfc to n_atomic_wfc1
!-----------------------------------------------------------------------
!
!     input: wfc, becwfc=<wfc|beta>, betae=|beta>
!     output: swfc=S|wfc>
!
   USE kinds, ONLY: DP
   USE ions_base, ONLY: na
   USE cvan, ONLY: nvb, ish
   USE uspp, ONLY: nhsa => nkb, nhsavb => nkbus, qq
   USE uspp_param, ONLY: nh
   USE gvecw, ONLY: ngw
   USE constants, ONLY: pi, fpi
   IMPLICIT NONE
! input
   INTEGER, INTENT(in)         :: n_atomic_wfc1
   COMPLEX(DP), INTENT(in) :: betae(ngw, nhsa),                   &
  &                               wfc(ngw, n_atomic_wfc1)
   REAL(DP), INTENT(in)    :: becwfc(nhsa, n_atomic_wfc1)
! output
   COMPLEX(DP), INTENT(out):: swfc(ngw, n_atomic_wfc1)
! local
   INTEGER is, iv, jv, ia, inl, jnl, i
   REAL(DP) qtemp(nhsavb, n_atomic_wfc1)
!
   swfc = wfc
!
   IF (nvb .GT. 0) THEN
      qtemp = 0.d0
      DO is = 1, nvb
         DO iv = 1, nh(is)
            DO jv = 1, nh(is)
               IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                  DO ia = 1, na(is)
                     inl = ish(is) + (iv - 1)*na(is) + ia
                     jnl = ish(is) + (jv - 1)*na(is) + ia
                     DO i = 1, n_atomic_wfc1
                        qtemp(inl, i) = qtemp(inl, i) +                &
  &                                    qq(iv, jv, is)*becwfc(jnl, i)
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END DO
!
      CALL DGEMM &
         ('N', 'N', 2*ngw, n_atomic_wfc1, nhsavb, 1.0d0, betae, 2*ngw, &
          qtemp, nhsavb, 1.0d0, swfc, 2*ngw)
!
   END IF
!
!      swfc=swfc+wfc
!
   RETURN
END SUBROUTINE s_wfc_real

!-------------------------------------------------------------------------
SUBROUTINE s_wfc_twin(n_atomic_wfc1, becwfc, betae, wfc, swfc, lgam) !@@@@ Changed n_atomic_wfc to n_atomic_wfc1
!-----------------------------------------------------------------------
!
!     input: wfc, becwfc=<wfc|beta>, betae=|beta>
!     output: swfc=S|wfc>
!
   USE kinds, ONLY: DP
   USE ions_base, ONLY: na
   USE cvan, ONLY: nvb, ish
   USE uspp, ONLY: nhsa => nkb, nhsavb => nkbus, qq
   USE uspp_param, ONLY: nh
   USE gvecw, ONLY: ngw
   USE constants, ONLY: pi, fpi
   USE twin_types
   !
   IMPLICIT NONE
! input
   INTEGER, INTENT(in)         :: n_atomic_wfc1
   COMPLEX(DP), INTENT(in) :: betae(ngw, nhsa),                   &
  &                               wfc(ngw, n_atomic_wfc1)
   TYPE(twin_matrix) :: becwfc
   LOGICAL :: lgam
! output
   COMPLEX(DP), INTENT(out):: swfc(ngw, n_atomic_wfc1)
! local
   INTEGER is, iv, jv, ia, inl, jnl, i
   REAL(DP), ALLOCATABLE :: qtemp(:, :)
   COMPLEX(DP), ALLOCATABLE :: qtemp_c(:, :)
!
   swfc = wfc
!
   IF (nvb .GT. 0) THEN
      IF (lgam) THEN
         allocate (qtemp(nhsavb, n_atomic_wfc1))
         qtemp = 0.d0
         DO is = 1, nvb
            DO iv = 1, nh(is)
               DO jv = 1, nh(is)
                  IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                     DO ia = 1, na(is)
                        inl = ish(is) + (iv - 1)*na(is) + ia
                        jnl = ish(is) + (jv - 1)*na(is) + ia
                        DO i = 1, n_atomic_wfc1
                           qtemp(inl, i) = qtemp(inl, i) +                &
     &                                 qq(iv, jv, is)*becwfc%rvec(jnl, i)
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
         !
         CALL DGEMM &
            ('N', 'N', 2*ngw, n_atomic_wfc1, nhsavb, 1.0d0, betae, 2*ngw, &
             qtemp, nhsavb, 1.0d0, swfc, 2*ngw)
         deallocate (qtemp)
      ELSE
         allocate (qtemp_c(nhsavb, n_atomic_wfc1))
         qtemp_c = CMPLX(0.d0, 0.d0)
         DO is = 1, nvb
            DO iv = 1, nh(is)
               DO jv = 1, nh(is)
                  IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                     DO ia = 1, na(is)
                        inl = ish(is) + (iv - 1)*na(is) + ia
                        jnl = ish(is) + (jv - 1)*na(is) + ia
                        DO i = 1, n_atomic_wfc1
                           qtemp_c(inl, i) = qtemp_c(inl, i) +                &
     &                                    qq(iv, jv, is)*becwfc%cvec(jnl, i)
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
         !
         CALL ZGEMM &
            ('N', 'N', ngw, n_atomic_wfc1, nhsavb, (1.0d0, 0.d0), betae, ngw, &
             qtemp_c, nhsavb, (1.0d0, 0.d0), swfc, ngw)
         deallocate (qtemp_c)
      END IF
   END IF
!
!      swfc=swfc+wfc
!
   RETURN
END SUBROUTINE s_wfc_twin

!-----------------------------------------------------------------------
SUBROUTINE spinsq(c, bec, rhor)
!-----------------------------------------------------------------------
!
!     estimate of <S^2>=s(s+1) in two different ways.
!     1) using as many-body wavefunction a single Slater determinant
!        constructed with Kohn-Sham orbitals:
!
!        <S^2> = (Nup-Ndw)/2 * (Nup-Ndw)/2+1) + Ndw -
!                \sum_up\sum_dw < psi_up | psi_dw >
!
!        where Nup, Ndw = number of up and down states, the sum is over
!        occupied states. Not suitable for fractionary occupancy.
!        In the ultrasoft scheme (c is the smooth part of \psi):
!
!        < psi_up | psi_dw > = \sum_G c*_up(G) c_dw(G) +
!                              \int Q_ij <c_up|beta_i><beta_j|c_dw>
!
!        This is the usual formula, unsuitable for fractionary occupancy.
!     2) using the "LSD model" of Wang, Becke, Smith, JCP 102, 3477 (1995):
!
!        <S^2> = (Nup-Ndw)/2 * (Nup-Ndw)/2+1) + Ndw -
!                \int max(rhoup(r),rhodw(r)) dr
!
!     Requires on input: c=psi, bec=<c|beta>, rhoup(r), rhodw(r)
!     Assumes real psi, with only half G vectors.
!
   USE electrons_base, ONLY: nx => nbspx, n => nbsp, iupdwn, nupdwn, f, nel, nspin
   USE io_global, ONLY: stdout
   USE mp_global, ONLY: intra_image_comm
   USE mp, ONLY: mp_sum
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: gstart
   USE grid_dimensions, ONLY: nr1, nr2, nr3, &
                              nnr => nnrx
   USE cell_base, ONLY: omega
   USE cvan, ONLY: nvb, ish
   USE uspp, ONLY: nhsa => nkb, qq
   USE uspp_param, ONLY: nh
   USE ions_base, ONLY: na
!
   IMPLICIT NONE
! input
   REAL(8) bec(nhsa, n), rhor(nnr, nspin)
   COMPLEX(8) c(ngw, nx)
! local variables
   INTEGER nup, ndw, ir, i, j, jj, ig, ia, is, iv, jv, inl, jnl
   REAL(8) spin0, spin1, spin2, fup, fdw
   REAL(8), ALLOCATABLE:: overlap(:, :), temp(:)
   LOGICAL frac
!
!
   IF (nspin .EQ. 1) RETURN
!
! find spin-up and spin-down states
!
   fup = 0.0d0
   DO i = iupdwn(1), nupdwn(1)
      fup = fup + f(i)
   END DO
   nup = NINT(fup)
   ndw = nel(1) + nel(2) - nup
!
! paranoid checks
!
   frac = ABS(fup - nup) .GT. 1.0d-6
   fup = 0.0d0
   DO i = 1, nup
      fup = fup + f(i)
   END DO
   frac = frac .OR. ABS(fup - nup) .GT. 1.0d-6
   fdw = 0.0d0
   DO j = iupdwn(2), iupdwn(2) - 1 + ndw
      fdw = fdw + f(j)
   END DO
   frac = frac .OR. ABS(fdw - ndw) .GT. 1.0d-6
!
   spin0 = ABS(fup - fdw)/2.d0*(ABS(fup - fdw)/2.d0 + 1.d0) + fdw
!
!     Becke's formula for spin polarization
!
   spin1 = 0.0d0
   DO ir = 1, nnr
      spin1 = spin1 - MIN(rhor(ir, 1), rhor(ir, 2))
   END DO
   CALL mp_sum(spin1, intra_image_comm)
   spin1 = spin0 + omega/(nr1*nr2*nr3)*spin1
   IF (frac) THEN
      WRITE (stdout, '(/'' Spin contamination: s(s+1)='',f5.2,'' (Becke) '',&
  &                             f5.2,'' (expected)'')')              &
  &          spin1, ABS(fup - fdw)/2.d0*(ABS(fup - fdw)/2.d0 + 1.d0)
      RETURN
   END IF
!
!     Slater formula, smooth contribution to  < psi_up | psi_dw >
!
   ALLOCATE (overlap(nup, ndw))
   ALLOCATE (temp(ngw))
   DO j = 1, ndw
      jj = j + iupdwn(2) - 1
      DO i = 1, nup
         overlap(i, j) = 0.d0
         DO ig = 1, ngw
            temp(ig) = DBLE(CONJG(c(ig, i))*c(ig, jj))
         END DO
         overlap(i, j) = 2.d0*SUM(temp)
         IF (gstart == 2) overlap(i, j) = overlap(i, j) - temp(1)
      END DO
   END DO
   DEALLOCATE (temp)
   CALL mp_sum(overlap, intra_image_comm)
   DO j = 1, ndw
      jj = j + iupdwn(2) - 1
      DO i = 1, nup
!
!     vanderbilt contribution to  < psi_up | psi_dw >
!
         DO is = 1, nvb
            DO iv = 1, nh(is)
               DO jv = 1, nh(is)
                  IF (ABS(qq(iv, jv, is)) .GT. 1.e-5) THEN
                     DO ia = 1, na(is)
                        inl = ish(is) + (iv - 1)*na(is) + ia
                        jnl = ish(is) + (jv - 1)*na(is) + ia
                        overlap(i, j) = overlap(i, j) +                &
  &                          qq(iv, jv, is)*bec(inl, i)*bec(jnl, jj)
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END DO
   END DO
!
   spin2 = spin0
   DO j = 1, ndw
      DO i = 1, nup
         spin2 = spin2 - overlap(i, j)**2
      END DO
   END DO
!
   DEALLOCATE (overlap)
!
   WRITE (stdout, '(/" Spin contamination: s(s+1)=",f5.2," (Slater) ",  &
  &          f5.2," (Becke) ",f5.2," (expected)")')              &
  &     spin2, spin1, ABS(fup - fdw)/2.d0*(ABS(fup - fdw)/2.d0 + 1.d0)
!
   RETURN
END SUBROUTINE spinsq

!
!-----------------------------------------------------------------------
SUBROUTINE vofrho(nfi, rhor, rhog, rhos, rhoc, tfirst, tlast,           &
&     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion)
!-----------------------------------------------------------------------
!     computes: the one-particle potential v in real space,
!               the total energy etot,
!               the forces fion acting on the ions,
!               the derivative of total energy to cell parameters h
!     rhor input : electronic charge on dense real space grid
!                  (plus core charge if present)
!     rhog input : electronic charge in g space (up to density cutoff)
!     rhos input : electronic charge on smooth real space grid
!     rhor output: total potential on dense real space grid
!     rhos output: total potential on smooth real space grid
!

   USE kinds, ONLY: dp
   USE control_flags, ONLY: iprint, iprsta, tpre, tfor, &
                            tprnfor, iesr, textfor, gamma_only, do_wf_cmplx !addded:giovanni
   USE io_global, ONLY: stdout
   USE ions_base, ONLY: nsp, na, nat, compute_eextfor
   USE gvecs
   USE gvecp, ONLY: ng => ngm
   USE cell_base, ONLY: omega, r_to_s
   USE cell_base, ONLY: a1, a2, a3, tpiba2, h, ainv
   USE reciprocal_vectors, ONLY: gstart, g, gx
   USE recvecs_indexes, ONLY: np, nm
   USE grid_dimensions, ONLY: nr1, nr2, nr3, &
                              nnr => nnrx
   USE smooth_grid_dimensions, &
      ONLY: nnrsx
   USE electrons_base, ONLY: nspin
   USE constants, ONLY: pi, fpi, au_gpa
   USE energies, ONLY: etot, eself, enl, ekin, epseu, esr, eht, exc, eextfor
   USE local_pseudo, ONLY: vps, rhops
   USE core, ONLY: nlcc_any
   USE gvecb
   USE dener, ONLY: detot, dekin, dps, dh, dsr, dxc, denl, &
                    detot6, dekin6, dps6, dh6, dsr6
   USE cp_main_variables, ONLY: drhog, ht0
   USE mp, ONLY: mp_sum
   USE mp_global, ONLY: intra_image_comm
   USE funct, ONLY: dft_is_meta, dft_is_hybrid
   USE pres_ai_mod, ONLY: abivol, abisur, v_vol, P_ext, volclu, &
                          Surf_t, surfclu
   USE cp_interfaces, ONLY: fwfft, invfft, self_vofhar
   USE sic_module, ONLY: self_interaction
   USE energies, ONLY: self_exc, self_ehte
   USE cp_interfaces, ONLY: compute_gagb, stress_hartree, &
                            add_drhoph, stress_local, force_loc
   USE fft_base, ONLY: dfftp, dffts
   USE ldaU, ONLY: e_hubbard
   USE hfmod, ONLY: do_hf, hfscalfact
   use eecp_mod, only: do_comp, vcorr, &
                       vcorr_fft, ecomp
   USE efield_mod, ONLY: do_efield, efieldpotg
   USE io_global, ONLY: stdout
   IMPLICIT NONE
!
   LOGICAL :: tlast, tfirst
   INTEGER :: nfi
   REAL(DP) rhor(nnr, nspin), rhos(nnrsx, nspin), fion(3, nat)
   REAL(DP) rhoc(nnr), tau0(3, nat)
   COMPLEX(DP) ei1(-nr1:nr1, nat), ei2(-nr2:nr2, nat),     &
  &                ei3(-nr3:nr3, nat), eigrb(ngb, nat),        &
  &                rhog(ng, nspin), sfac(ngs, nsp)
   !
   INTEGER irb(3, nat)
   !
   INTEGER iss, isup, isdw, ig, ir, i, j, k, ij, is
   REAL(DP) vave, ebac, wz, eh
   COMPLEX(DP) fp, fm, ci, zpseu, zh
   COMPLEX(DP), ALLOCATABLE :: rhotmp(:), vtemp(:), aux(:)
   ! COMPLEX(DP), ALLOCATABLE :: drhotmp(:,:,:)
   COMPLEX(DP), ALLOCATABLE :: drhot(:, :)
   COMPLEX(DP), ALLOCATABLE :: v(:), vs(:)
   REAL(DP), ALLOCATABLE    :: gagb(:, :), rhotot(:), v3d(:), v0d(:)
   !
   REAL(DP), ALLOCATABLE :: fion1(:, :)
   REAL(DP), ALLOCATABLE :: stmp(:, :)
   !
   COMPLEX(DP), ALLOCATABLE :: self_vloc(:)
   REAL(DP)                 :: self_ehtet
   LOGICAL                  :: ttsic
   REAL(DP)                 :: detmp(3, 3)
   REAL(DP)                 :: ht(3, 3)
   COMPLEX(DP)              :: screen_coul(1)
!
   INTEGER, DIMENSION(6), PARAMETER :: alpha = (/1, 2, 3, 2, 3, 3/)
   INTEGER, DIMENSION(6), PARAMETER :: beta = (/1, 1, 1, 2, 2, 3/)

   ! ...  dalbe(:) = delta( alpha(:), beta(:) )
   REAL(DP), DIMENSION(6), PARAMETER :: dalbe = &
                                        (/1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP/)
   LOGICAL :: lgam !added:giovanni

   CALL start_clock('vofrho')
   lgam = gamma_only .and. .not. do_wf_cmplx !added:giovanni

   ci = (0.0d0, 1.0d0)
   !
   !     wz = factor for g.neq.0 because of c*(g)=c(-g)
   !
   wz = 2.0d0
   !
   ht = TRANSPOSE(h)
   !
   ALLOCATE (v(nnr))
   ALLOCATE (vs(nnrsx))
   ALLOCATE (vtemp(ng))
   ALLOCATE (rhotmp(ng))
   if (do_comp) then
      allocate (aux(nnr))
      allocate (rhotot(nnr))
      allocate (v3d(nnr))
      allocate (v0d(nnr))
   end if
   !
   IF (tpre) THEN
      ALLOCATE (drhot(ng, 6))
      ALLOCATE (gagb(6, ng))
      CALL compute_gagb(gagb, gx, ng, tpiba2)
   END IF
   !
!
!     ab-initio pressure and surface tension contributions to the potential
!
   if (abivol .or. abisur) call vol_clu(rhor, rhog, nfi)
   !
   ttsic = (ABS(self_interaction) /= 0)
   !
   IF (ttsic) ALLOCATE (self_vloc(ng))
   !
   !     first routine in which fion is calculated: annihilation
   !
   fion = 0.d0
   !
   !     forces on ions, ionic term in real space
   !
   IF (tprnfor .OR. tfor .OR. tfirst .OR. tpre) THEN
      !
      ALLOCATE (stmp(3, nat))
      !
      CALL r_to_s(tau0, stmp, na, nsp, ainv)
      !
      CALL vofesr(iesr, esr, dsr6, fion, stmp, tpre, h)
      !
      call mp_sum(fion, intra_image_comm)
      !
      DEALLOCATE (stmp)
      !
   END IF
!
   rhotmp(1:ng) = rhog(1:ng, 1)
   !

   IF (tpre) THEN
      DO ij = 1, 6
         i = alpha(ij)
         j = beta(ij)
         drhot(:, ij) = 0.0d0
         DO k = 1, 3
            drhot(:, ij) = drhot(:, ij) + drhog(:, 1, i, k)*ht(k, j)
         END DO
      END DO
   END IF
   !
   IF (nspin == 2) THEN
      rhotmp(1:ng) = rhotmp(1:ng) + rhog(1:ng, 2)
      IF (tpre) THEN
         DO ij = 1, 6
            i = alpha(ij)
            j = beta(ij)
            DO k = 1, 3
               drhot(:, ij) = drhot(:, ij) + drhog(:, 2, i, k)*ht(k, j)
            END DO
         END DO
      END IF
   END IF
   !
   !     calculation local potential energy
   !
   zpseu = 0.0d0
   !
!$omp parallel default(shared), private(ig,is)
!$omp do
   DO ig = 1, SIZE(vtemp)
      vtemp(ig) = CMPLX(0.d0, 0.d0)
   END DO
   DO is = 1, nsp
!$omp do
      DO ig = 1, ngs
         vtemp(ig) = vtemp(ig) + CONJG(rhotmp(ig))*sfac(ig, is)*vps(ig, is)
      END DO
   END DO
!$omp do reduction(+:zpseu)
   DO ig = 1, ngs
      zpseu = zpseu + vtemp(ig)
   END DO
!$omp end parallel

   epseu = wz*DBLE(zpseu)
   !
   IF (lgam) THEN
      IF (gstart == 2) epseu = epseu - DBLE(vtemp(1))
   END IF
   !
   CALL mp_sum(epseu, intra_image_comm)

   IF (lgam) THEN
      epseu = epseu*omega
   ELSE
      epseu = 0.5d0*epseu*omega
   END IF
!
   IF (tpre) THEN
      !
      CALL stress_local(dps6, gagb, sfac, rhotmp, drhot, omega)
      !
   END IF
   !
   !
   !     calculation hartree energy
   !
   !
   self_ehtet = 0.d0
   !
   IF (ttsic) self_vloc = 0.d0

   zh = 0.0d0

!$omp parallel default(shared), private(ig,is)

   DO is = 1, nsp
!$omp do
      DO ig = 1, ngs
         rhotmp(ig) = rhotmp(ig) + sfac(ig, is)*rhops(ig, is) !JUST-FOR-NOW
         !rhotmp(ig)=rhotmp(ig)
      END DO
   END DO
   !
!$omp do
   DO ig = gstart, ng
      vtemp(ig) = CONJG(rhotmp(ig))*rhotmp(ig)/g(ig)
   END DO

!$omp do reduction(+:zh)
   DO ig = gstart, ng
      zh = zh + vtemp(ig)
   END DO

!$omp end parallel
   IF (lgam) THEN
      eh = DBLE(zh)*wz*0.5d0*fpi/tpiba2
   ELSE
      eh = DBLE(zh)*wz*0.25d0*fpi/tpiba2
   END IF
!
   CALL mp_sum(eh, intra_image_comm)
   !
   IF (ttsic) THEN
      !
      CALL self_vofhar(.false., self_ehte, self_vloc, rhog, omega, h)
      !
      eh = eh - self_ehte/omega
      !
   END IF
   !
   IF (tpre) THEN
      !
      CALL add_drhoph(drhot, sfac, gagb)
      !
      CALL stress_hartree(dh6, eh*omega, rhotmp, drhot, gagb, omega)
      !
   END IF
   !
   IF (tpre) THEN
      DEALLOCATE (drhot)
   END IF
   !
   !     forces on ions, ionic term in reciprocal space
   !
   ALLOCATE (fion1(3, nat))
   !
   fion1 = 0.d0
   !
   IF (tprnfor .OR. tfor .OR. tpre) THEN
      vtemp(1:ng) = rhog(1:ng, 1)
      IF (nspin == 2) THEN
         vtemp(1:ng) = vtemp(1:ng) + rhog(1:ng, 2)
      END IF
      CALL force_loc(.false., vtemp, fion1, rhops, vps, ei1, ei2, ei3, sfac, omega, screen_coul, lgam)
   END IF
   !
   !     calculation hartree + local pseudo potential
   !
   IF (gstart == 2) vtemp(1) = CMPLX(0.d0, 0.d0)

!$omp parallel default(shared), private(ig,is)
!$omp do
   DO ig = gstart, ng
      vtemp(ig) = rhotmp(ig)*fpi/(tpiba2*g(ig))
   END DO
!
   DO is = 1, nsp
!$omp do
      DO ig = 1, ngs
         vtemp(ig) = vtemp(ig) + sfac(ig, is)*vps(ig, is)  !JUST-FOR-NOW I do not consider the pseudopotential part
      END DO
   END DO
   !
   if (do_comp) then
      !
      call calc_compensation_potential(vcorr_fft, rhotmp, .false.)
      !
      call calc_tcc_energy(ecomp, vcorr_fft, rhotmp, lgam)
      !
      aux = 0.0_dp

!         IF(lgam) THEN !!!### uncomment for k points
      do ig = 1, ng
         aux(np(ig)) = vcorr_fft(ig)
         aux(nm(ig)) = conjg(vcorr_fft(ig))
      end do
!         ELSE !!!### uncomment for k points
!           do ig=1,ng !!!### uncomment for k points
!             aux(np(ig))=vcorr_fft(ig) !!!### uncomment for k points
! !             aux(nm(ig))=conjg(vcorr_fft(ig))
!           end do !!!### uncomment for k points
!         ENDIF
      call invfft('Dense', aux, dfftp)
      vcorr = dble(aux)
      call writetofile(vcorr, nnr, 'vcorrz.dat', dfftp, 'az')
      call writetofile(vcorr, nnr, 'vcorrx.dat', dfftp, 'ax')
      !
      if (tprnfor .or. tfor .or. tfirst .or. tpre) then
         allocate (stmp(3, nat))
         call r_to_s(tau0, stmp, na, nsp, ainv)
         call calc_fcorr(fion, vcorr, stmp, nat, na, nsp, ht0, dfftp)
         deallocate (stmp)
      end if
      !
      aux = 0.0_dp
!         IF(lgam) THEN !!!### uncomment for k points
      do ig = 1, ng
         aux(np(ig)) = vcorr_fft(ig) + vtemp(ig)
         aux(nm(ig)) = conjg(vcorr_fft(ig) + vtemp(ig))
      end do
!         ELSE !!!### uncomment for k points
!           do ig=1,ng !!!### uncomment for k points
!             aux(np(ig))=vcorr_fft(ig)+vtemp(ig) !!!### uncomment for k points
!           end do !!!### uncomment for k points
!         ENDIF !!!### uncomment for k points
      !
      call invfft('Dense', aux, dfftp)
      v0d = dble(aux)
      !
      call writetofile(v0d, nnr, 'v0dz.dat', dfftp, 'az')
      call writetofile(v0d, nnr, 'v0dx.dat', dfftp, 'ax')
      !
      aux = 0.0_dp
!         IF(lgam) THEN !!!### uncomment for k points
      do ig = 1, ng
         aux(np(ig)) = vtemp(ig)
         aux(nm(ig)) = conjg(vtemp(ig))
      end do
!         ELSE !!!### uncomment for k points
!           do ig=1,ng !!!### uncomment for k points
!             aux(np(ig))=vtemp(ig) !!!### uncomment for k points
!           end do !!!### uncomment for k points
!         ENDIF !!!### uncomment for k points
      call invfft('Dense', aux, dfftp)
      v3d = dble(aux)
      call writetofile(v3d, nnr, 'v3dz.dat', dfftp, 'az')
      call writetofile(v3d, nnr, 'v3dx.dat', dfftp, 'ax')
      !
      aux = 0.0_dp

!         IF(lgam) THEN  !!!### uncomment for k points
      do ig = 1, ng
         aux(np(ig)) = rhotmp(ig)
         aux(nm(ig)) = conjg(rhotmp(ig))
      end do
!         ELSE !!!### uncomment for k points
!           do ig=1,ng !!!### uncomment for k points
!             aux(np(ig))=rhotmp(ig) !!!### uncomment for k points
!           end do !!!### uncomment for k points
!         ENDIF !!!### uncomment for k points
      call invfft('Dense', aux, dfftp)
      rhotot = dble(aux)
      call writetofile(rhotot, nnr, 'rhototz.dat', dfftp, 'az')
      call writetofile(rhotot, nnr, 'rhototx.dat', dfftp, 'ax')
      !
      vtemp = vtemp + vcorr_fft
      eh = eh + ecomp/omega
      !
   end if
   if (do_efield) vtemp = vtemp + efieldpotg
!$omp end parallel
!
!     vtemp = v_loc(g) + v_h(g)
!
!     ===================================================================
!      calculation exchange and correlation energy and potential
!     -------------------------------------------------------------------
   IF (nlcc_any) CALL add_cc(rhoc, rhog, rhor)
!
   CALL exch_corr_h(nspin, rhog, rhor, rhoc, sfac, exc, dxc, self_exc)
   call writetofile(rhor, nnr, 'vxc.dat', dfftp, 'az')
   !
   ! correction for traditional use of HF without
   ! reference to the EXX implementation
   !
   IF (do_hf .AND. .NOT. dft_is_hybrid()) THEN
      !
      rhor = rhor*(1.0d0 - hfscalfact)
      exc = exc*(1.0d0 - hfscalfact)
      dxc = dxc*(1.0d0 - hfscalfact)
      self_exc = self_exc*(1.0d0 - hfscalfact)
      !
   END IF
!
!     rhor contains the xc potential in r-space
!
!     ===================================================================
!     fourier transform of xc potential to g-space (dense grid)
!     -------------------------------------------------------------------
!
   IF (nspin == 1) THEN
      iss = 1
      if (abivol .or. abisur) then
!$omp parallel do
         do ir = 1, nnr
            v(ir) = CMPLX(rhor(ir, iss) + v_vol(ir), 0.d0)
         end do
      else
!$omp parallel do
         do ir = 1, nnr
            v(ir) = CMPLX(rhor(ir, iss), 0.d0)
         end do
      end if
      !
      !     v_xc(r) --> v_xc(g)
      !
      CALL fwfft('Dense', v, dfftp)
!
!$omp parallel do
      DO ig = 1, ng
         rhog(ig, iss) = vtemp(ig) + v(np(ig))
      END DO
      !
      !     v_tot(g) = (v_tot(g) - v_xc(g)) +v_xc(g)
      !     rhog contains the total potential in g-space
      !
   ELSE
      isup = 1
      isdw = 2
      if (abivol .or. abisur) then
!$omp parallel do
         do ir = 1, nnr
            v(ir) = CMPLX(rhor(ir, isup) + v_vol(ir), rhor(ir, isdw) + v_vol(ir))
         end do
      else
!$omp parallel do
         do ir = 1, nnr
            v(ir) = CMPLX(rhor(ir, isup), rhor(ir, isdw))
         end do
      end if
      CALL fwfft('Dense', v, dfftp)
!$omp parallel do private(fp,fm)
      DO ig = 1, ng
         fp = v(np(ig)) + v(nm(ig))
         fm = v(np(ig)) - v(nm(ig))
         IF (ttsic) THEN
            rhog(ig, isup) = vtemp(ig) - self_vloc(ig) + 0.5d0*CMPLX(DBLE(fp), AIMAG(fm))
            rhog(ig, isdw) = vtemp(ig) + self_vloc(ig) + 0.5d0*CMPLX(AIMAG(fp), -DBLE(fm))
         ELSE
            rhog(ig, isup) = vtemp(ig) + 0.5d0*CMPLX(DBLE(fp), AIMAG(fm))
            rhog(ig, isdw) = vtemp(ig) + 0.5d0*CMPLX(AIMAG(fp), -DBLE(fm))
         END IF
      END DO
   END IF

!
!     rhog contains now the total (local+Hartree+xc) potential in g-space
!
   IF (tprnfor .OR. tfor) THEN

      IF (nlcc_any) CALL force_cc(irb, eigrb, rhor, fion1)

      CALL mp_sum(fion1, intra_image_comm)
      !
      !    add g-space ionic and core correction contributions to fion
      !
      fion = fion + fion1

   END IF

   DEALLOCATE (fion1)
!
   IF (ttsic) DEALLOCATE (self_vloc)
!
!     ===================================================================
!     fourier transform of total potential to r-space (dense grid)
!     -------------------------------------------------------------------
   v(:) = (0.d0, 0.d0)
   IF (nspin .EQ. 1) THEN
      iss = 1
!$omp parallel do
!         IF(lgam) THEN !!!### uncomment for k points
      DO ig = 1, ng
         v(np(ig)) = rhog(ig, iss)
         v(nm(ig)) = CONJG(rhog(ig, iss))
      END DO
!         ELSE !!!### uncomment for k points
!             DO ig=1,ng !!!### uncomment for k points
!                 v(np(ig))=rhog(ig,iss) !!!### uncomment for k points
!             END DO !!!### uncomment for k points
!         ENDIF !!!### uncomment for k points
!
!     v(g) --> v(r)
!
      CALL invfft('Dense', v, dfftp)
!
!$omp parallel do
      DO ir = 1, nnr
         rhor(ir, iss) = DBLE(v(ir))
      END DO
!
!     calculation of average potential
!
      vave = SUM(rhor(:, iss))/DBLE(nr1*nr2*nr3)
   ELSE
      isup = 1
      isdw = 2
!$omp parallel do
!          IF(lgam) THEN !!!### uncomment for k points
      DO ig = 1, ng
         v(np(ig)) = rhog(ig, isup) + ci*rhog(ig, isdw)
         v(nm(ig)) = CONJG(rhog(ig, isup)) + ci*CONJG(rhog(ig, isdw))
      END DO
!          ELSE !!!### uncomment for k points
!             DO ig=1,ng !!!### uncomment for k points
!                 v(np(ig))=rhog(ig,isup)+ci*rhog(ig,isdw) !!!### uncomment for k points
!                 v(nm(ig))=CONJG(rhog(ig,isup)) +ci*CONJG(rhog(ig,isdw))
!             END DO !!!### uncomment for k points
!          ENDIF !!!### uncomment for k points
!
      CALL invfft('Dense', v, dfftp)
!$omp parallel do
      DO ir = 1, nnr
         rhor(ir, isup) = DBLE(v(ir))
         rhor(ir, isdw) = AIMAG(v(ir))
      END DO
      !
      !     calculation of average potential
      !
      vave = (SUM(rhor(:, isup)) + SUM(rhor(:, isdw)))/2.0d0/DBLE(nr1*nr2*nr3)
   END IF

   CALL mp_sum(vave, intra_image_comm)

   !
   !     fourier transform of total potential to r-space (smooth grid)
   !
   vs(:) = CMPLX(0.d0, 0.d0)
   !
   IF (nspin .EQ. 1) THEN
      !
      iss = 1
!$omp parallel do
!        IF(lgam) THEN !!!### uncomment for k points
      DO ig = 1, ngs
         vs(nms(ig)) = CONJG(rhog(ig, iss))
         vs(nps(ig)) = rhog(ig, iss)
      END DO
!        ELSE !!!### uncomment for k points
!           DO ig=1,ngs !!!### uncomment for k points
!               write(6,*) "debug", nps(ig), nnrsx, ig, ngs  !added:giovanni:debug
!               vs(nps(ig))=rhog(ig,iss) !!!### uncomment for k points
!           END DO !!!### uncomment for k points
!        ENDIF !!!### uncomment for k points
      !
      CALL invfft('Smooth', vs, dffts)
      !
!$omp parallel do
      DO ir = 1, nnrsx
         rhos(ir, iss) = DBLE(vs(ir))
      END DO
      !
   ELSE
      !
      isup = 1
      isdw = 2
!$omp parallel do
!          IF(lgam) THEN !!!### uncomment for k points
      DO ig = 1, ngs
         vs(nps(ig)) = rhog(ig, isup) + ci*rhog(ig, isdw)
         vs(nms(ig)) = CONJG(rhog(ig, isup)) + ci*CONJG(rhog(ig, isdw))
      END DO
!          ELSE !!!### uncomment for k points
!             DO ig=1,ngs !!!### uncomment for k points
!                 vs(nps(ig))=rhog(ig,isup)+ci*rhog(ig,isdw) !!!### uncomment for k points
!                 vs(nms(ig))=CONJG(rhog(ig,isup)) +ci*CONJG(rhog(ig,isdw))
!             END DO !!!### uncomment for k points
!          ENDIF !!!### uncomment for k points
      !
      CALL invfft('Smooth', vs, dffts)
      !
!$omp parallel do
      DO ir = 1, nnrsx
         rhos(ir, isup) = DBLE(vs(ir))
         rhos(ir, isdw) = AIMAG(vs(ir))
      END DO
      !
   END IF
   !

   IF (dft_is_meta()) CALL vofrho_meta(v, vs)  !METAGGA

   ebac = 0.0d0
   !
   eht = eh*omega + esr - eself
   !
   eextfor = 0.0_DP
   IF (textfor) eextfor = compute_eextfor(tau0)
   !
   !     etot is the total energy ; ekin, enl were calculated in rhoofr
   !
   etot = ekin + eht + epseu + enl + exc + ebac + e_hubbard + eextfor
   !
   !     extra contributions
   !
   if (abivol) etot = etot + P_ext*volclu
   if (abisur) etot = etot + Surf_t*surfclu
   !
   IF (tpre) THEN
      !
      detot6 = dekin6 + dh6 + dps6 + dsr6
      !
      call mp_sum(detot6, intra_image_comm)
      !
      DO k = 1, 6
         detmp(alpha(k), beta(k)) = detot6(k)
         detmp(beta(k), alpha(k)) = detmp(alpha(k), beta(k))
      END DO
      !
      detot = MATMUL(detmp(:, :), TRANSPOSE(ainv(:, :)))
      !
      detot = detot + denl + dxc
      !
   END IF
   !
   !
   CALL stop_clock('vofrho')
   !
   !
   IF (tpre) THEN
      !
      DEALLOCATE (gagb)
      !
      IF ((iprsta >= 2) .AND. (MOD(nfi - 1, iprint) == 0)) THEN
         !
         WRITE (stdout, *)
         WRITE (stdout, *) "From vofrho:"
         WRITE (stdout, *) "cell parameters h"
         WRITE (stdout, 5555) (a1(i), a2(i), a3(i), i=1, 3)
         !
         WRITE (stdout, *)
         WRITE (stdout, *) "derivative of e(tot)"
         WRITE (stdout, 5555) ((detot(i, j), j=1, 3), i=1, 3)
         WRITE (stdout, *) "kbar"
         detmp = -1.0d0*MATMUL(detot, TRANSPOSE(h))/omega*au_gpa*10.0d0
         WRITE (stdout, 5555) ((detmp(i, j), j=1, 3), i=1, 3)
         !
         WRITE (stdout, *)
         WRITE (stdout, *) "derivative of e(kin)"
         WRITE (stdout, 5555) ((dekin(i, j), j=1, 3), i=1, 3)
         WRITE (stdout, *) "kbar"
         detmp = -1.0d0*MATMUL(dekin, TRANSPOSE(h))/omega*au_gpa*10.0d0
         WRITE (stdout, 5555) ((detmp(i, j), j=1, 3), i=1, 3)
         !
         WRITE (stdout, *) "derivative of e(h)"
         WRITE (stdout, 5555) ((dh(i, j), j=1, 3), i=1, 3)
         WRITE (stdout, *) "kbar"
         detmp = -1.0d0*MATMUL(dh, TRANSPOSE(h))/omega*au_gpa*10.0d0
         WRITE (stdout, 5555) ((detmp(i, j), j=1, 3), i=1, 3)
         !
         WRITE (stdout, *) "derivative of e(sr)"
         WRITE (stdout, 5555) ((dsr(i, j), j=1, 3), i=1, 3)
         WRITE (stdout, *) "kbar"
         detmp = -1.0d0*MATMUL(dsr, TRANSPOSE(h))/omega*au_gpa*10.0d0
         WRITE (stdout, 5555) ((detmp(i, j), j=1, 3), i=1, 3)
         !
         WRITE (stdout, *) "derivative of e(ps)"
         WRITE (stdout, 5555) ((dps(i, j), j=1, 3), i=1, 3)
         WRITE (stdout, *) "kbar"
         detmp = -1.0d0*MATMUL(dps, TRANSPOSE(h))/omega*au_gpa*10.0d0
         WRITE (stdout, 5555) ((detmp(i, j), j=1, 3), i=1, 3)
         !
         WRITE (stdout, *) "derivative of e(nl)"
         WRITE (stdout, 5555) ((denl(i, j), j=1, 3), i=1, 3)
         WRITE (stdout, *) "kbar"
         detmp = -1.0d0*MATMUL(denl, TRANSPOSE(h))/omega*au_gpa*10.0d0
         WRITE (stdout, 5555) ((detmp(i, j), j=1, 3), i=1, 3)
         !
         WRITE (stdout, *) "derivative of e(xc)"
         WRITE (stdout, 5555) ((dxc(i, j), j=1, 3), i=1, 3)
         WRITE (stdout, *) "kbar"
         detmp = -1.0d0*MATMUL(dxc, TRANSPOSE(h))/omega*au_gpa*10.0d0
         WRITE (stdout, 5555) ((detmp(i, j), j=1, 3), i=1, 3)
      END IF
   END IF

   DEALLOCATE (rhotmp)
   DEALLOCATE (vtemp)
   DEALLOCATE (v)
   DEALLOCATE (vs)
   if (do_comp) then
      deallocate (aux)
      deallocate (rhotot)
      deallocate (v0d)
      deallocate (v3d)
   end if

   RETURN

5555 FORMAT(1x, f12.5, 1x, f12.5, 1x, f12.5/                                &
                                                       &       1x, f12.5, 1x, f12.5, 1x, f12.5/                                &
                                                                &       1x, f12.5, 1x, f12.5, 1x, f12.5//)
!

END SUBROUTINE vofrho

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-----------------------------------------------------------------------
subroutine ldaU_init
!-----------------------------------------------------------------------
!
   USE constants, ONLY: autoev
   use ldaU, ONLY: n_atomic_wfc, atomwfc, lda_plus_u, Hubbard_U
   use ldaU, ONLY: Hubbard_lmax, Hubbard_l, ns, ns_c, vupsi
   use input_parameters, ONLY: atom_label, lda_plus_u_ => lda_plus_u
   use input_parameters, ONLY: Hubbard_U_ => Hubbard_U
   use ions_base, only: na, nsp, nat
   use gvecw, only: ngw
   use electrons_base, only: nspin, nx => nbspx
   USE uspp_param, ONLY: upf
   USE control_flags, ONLY: gamma_only, do_wf_cmplx
   !
   implicit none
   integer is, nb, l
   integer, external :: set_Hubbard_l
   logical :: lgam

   lgam = gamma_only .and. .not. do_wf_cmplx
! allocate vupsi
   lda_plus_u = lda_plus_u_

   allocate (vupsi(ngw, nx))

   vupsi = (0.0d0, 0.0d0)
   ! allocate(vpsi_con(ngw,nx)) ! step_constraint
   n_atomic_wfc = 0

   do is = 1, nsp
      !
      Hubbard_U(is) = Hubbard_U_(is)/autoev
      !
      do nb = 1, upf(is)%nwfc
         l = upf(is)%lchi(nb)
         n_atomic_wfc = n_atomic_wfc + (2*l + 1)*na(is)
      end do
      !
   end do
!
   allocate (atomwfc(ngw, n_atomic_wfc))

   if (lda_plus_u) then
      Hubbard_lmax = -1
      do is = 1, nsp
         if (Hubbard_U(is) .ne. 0.d0) then
!                Hubbard_l(is)=2
            Hubbard_l(is) = set_Hubbard_l(atom_label(is))
            Hubbard_lmax = max(Hubbard_lmax, Hubbard_l(is))
            write (6, *) ' HUBBARD L FOR TYPE ', atom_label(is), ' IS ',&
  &                       Hubbard_l(is)
         end if
      end do
      write (6, *) ' MAXIMUM HUBBARD L IS ', Hubbard_lmax
      if (Hubbard_lmax .eq. -1) call errore                            &
  &        ('setup', 'lda_plus_u calculation but Hubbard_l not set', 1)
   end if
   l = 2*Hubbard_lmax + 1
   !
   IF (lgam) THEN
      allocate (ns(nat, nspin, l, l))
   ELSE
      allocate (ns_c(nat, nspin, l, l))
   END IF
   !
   return
end subroutine ldaU_init
!
!-----------------------------------------------------------------------
subroutine nksic_init
!-----------------------------------------------------------------------
!
   !
   ! this routine is called anyway, even if do_nk=F
   !
   use nksic, ONLY: do_orbdep, do_nk, do_nkipz, do_nkpz, do_pz, &
                    do_nki, do_bare_eigs, do_pz_renorm, kfact, &
                    do_wref, do_wxd, fref, rhobarfact, &
                    vanishing_rho_w, &
                    nknmax, do_spinsym, f_cutoff, &
                    nkscalfact, nksic_memusage, allocate_nksic, odd_alpha
!$$
   use nksic, only: do_innerloop, do_innerloop_empty, do_innerloop_cg, &
                    innerloop_dd_nstep, &
                    innerloop_cg_nsd, innerloop_cg_nreset, innerloop_nmax, &
                    innerloop_cg_ratio, innerloop_init_n, innerloop_until, &
                    innerloop_atleast, l_comp_cmplxfctn_index
!$$
   use input_parameters, ONLY: do_nk_ => do_nk, &
                               do_pz_ => do_pz, &
                               do_nki_ => do_nki, &
                               do_nkpz_ => do_nkpz, &
                               do_nkipz_ => do_nkipz, &
                               do_hf_ => do_hf, &
                               which_orbdep_ => which_orbdep, &
                               fref_ => fref, &
                               rhobarfact_ => rhobarfact, &
                               vanishing_rho_w_ => vanishing_rho_w, &
                               do_wref_ => do_wref, &
                               do_wxd_ => do_wxd, &
                               do_spinsym_ => do_spinsym, &
                               nkscalfact_ => nkscalfact, &
                               nknmax_ => nknmax, &
                               f_cutoff_ => f_cutoff, &
                               do_orbdep_ => do_orbdep, &
                               l_comp_cmplxfctn_index_ => l_comp_cmplxfctn_index
!$$
   use input_parameters, only: do_innerloop_ => do_innerloop, &
                               do_innerloop_empty_ => do_innerloop_empty, &
                               do_innerloop_cg_ => do_innerloop_cg, &
                               innerloop_dd_nstep_ => innerloop_dd_nstep, &
                               innerloop_cg_nsd_ => innerloop_cg_nsd, &
                               innerloop_cg_nreset_ => innerloop_cg_nreset, &
                               innerloop_nmax_ => innerloop_nmax, &
                               innerloop_init_n_ => innerloop_init_n, &
                               innerloop_atleast_ => innerloop_atleast, &
                               innerloop_cg_ratio_ => innerloop_cg_ratio, &
                               innerloop_until_ => innerloop_until, &
                               do_pz_renorm_ => do_pz_renorm, &
                               do_bare_eigs_ => do_bare_eigs, &
                               kfact_ => kfact
!$$
   USE io_global, ONLY: meta_ionode, stdout
   use electrons_base, ONLY: nspin, nbspx
   use gvecw, ONLY: ngw
   use grid_dimensions, ONLY: nnrx
   USE ions_base, ONLY: nat
   !
   implicit none
   !
   logical       :: found, do_hybrid = .FALSE.
   integer       :: i
   character(10) :: subname = 'nksic_init'
   character(1), external :: lowercase

   !
   ! overwriten by which_orbdep, if not empty
   !
   do_nk = do_nk_
   do_pz = do_pz_
   do_nki = do_nki_
   do_nkpz = do_nkpz_
   do_nkipz = do_nkipz_
   !
   do_wxd = do_wxd_
   do_wref = do_wref_
   do_pz_renorm = do_pz_renorm_
   do_bare_eigs = do_bare_eigs_
   kfact = kfact_
   !
   fref = fref_
!$$
   do_innerloop = do_innerloop_
   do_innerloop_empty = do_innerloop_empty_
   l_comp_cmplxfctn_index = l_comp_cmplxfctn_index_
   do_innerloop_cg = do_innerloop_cg_
   innerloop_dd_nstep = innerloop_dd_nstep_
   innerloop_cg_nsd = innerloop_cg_nsd_
   innerloop_cg_nreset = innerloop_cg_nreset_
   innerloop_nmax = innerloop_nmax_
   innerloop_init_n = innerloop_init_n_
   innerloop_atleast = innerloop_atleast_
   innerloop_cg_ratio = innerloop_cg_ratio_
   innerloop_until = innerloop_until_
!$$
   !
   ! use the collective var which_orbdep
   !
   DO i = 1, LEN_TRIM(which_orbdep_)
      which_orbdep_(i:i) = lowercase(which_orbdep_(i:i))
   END DO
   !
   SELECT CASE (TRIM(which_orbdep_))
   CASE ("", "none")
      ! do nothing
   CASE ("hf", "b3lyp", "pbe0")
      do_hybrid = .TRUE.
   CASE ("nk", "non-koopmans")
      do_nk = .TRUE.
      do_wref = .TRUE.
      do_wxd = .TRUE.
   CASE ("nk0")
      do_nk = .TRUE.
      do_wref = .FALSE.
      do_wxd = .FALSE.
   CASE ("nki")
      do_nki = .TRUE.
      do_wxd = .TRUE.
      fref = 1.0
   CASE ("pz", "perdew-zunger")
      do_pz = .TRUE.
   CASE ("nkpz", "pznk")
      do_nkpz = .TRUE.
   CASE ("nkipz", "pznki")
      do_nkipz = .TRUE.
      do_wxd = .TRUE.
      fref = 1.0
   CASE DEFAULT
      call errore(subname, "invalid which_orbdep = "//TRIM(which_orbdep_), 10)
   END SELECT
   !
   IF (.NOT. do_hybrid .AND. do_hf_) do_hybrid = .TRUE.
   !
   do_orbdep = do_orbdep_ .and. .not. do_hybrid

   found = .FALSE.
   !
   IF (do_orbdep) THEN
      !
      if (do_nk .or. do_pz .or. do_nki .or. do_nkpz .or. do_nkipz .or. do_hybrid) found = .true.
      !
      if (.not. found) CALL errore(subname, 'no compatible orbital-dependent scheme specified', 1)
      !
   END IF
   !
   !
   ! check only one orbital dependent scheme is used
   !
   found = .FALSE.
   !
   if (do_nk .and. (do_pz .or. do_nki .or. do_nkpz .or. do_nkipz)) found = .TRUE.
   if (do_nki .and. (do_pz .or. do_nk .or. do_nkpz .or. do_nkipz)) found = .TRUE.
   if (do_pz .and. (do_nk .or. do_nki .or. do_nkpz .or. do_nkipz)) found = .TRUE.
   if (do_nkpz .and. (do_nk .or. do_nki .or. do_pz .or. do_nkipz)) found = .TRUE.
   if (do_nkipz .and. (do_nk .or. do_nki .or. do_pz .or. do_nkpz)) found = .TRUE.
   !
   if (found) CALL errore(subname, 'more than one orb-dependent schme used', 1)
   !
   do_spinsym = do_spinsym_
   vanishing_rho_w = vanishing_rho_w_
   rhobarfact = rhobarfact_
   nkscalfact = nkscalfact_
   nknmax = nknmax_
   f_cutoff = f_cutoff_
   !
   if (do_nki .and. fref /= 1.0) CALL errore(subname, 'nki and fref /= 1.0 ', 1)
   !
   if ((do_nk .or. do_nkpz) .and. meta_ionode) then
      write (stdout, 2000) fref
      write (stdout, 2004) rhobarfact, nkscalfact
   else if (do_pz .and. meta_ionode) then
      write (stdout, 2001) do_pz
   else if ((do_nki .or. do_nkipz) .and. meta_ionode) then
      write (stdout, 2002) do_nki
      write (stdout, 2004) rhobarfact, nkscalfact
   end if
   !
   ! read referece alpha from file, if any | linh
   ! wherein, the information of n_evc0_fixed, ref_alpha0,
   ! broadening of orbitals will be readed.
   !
   ! call readfile_refalpha()
   !
   if (do_orbdep .and. .not. do_hybrid) call allocate_nksic(nnrx, ngw, nspin, nbspx, nat)
   !
   if (do_orbdep) odd_alpha(:) = 1.d0
   !
   if ((do_nk .or. do_nkpz) .and. meta_ionode) then
      write (stdout, 2010) do_wxd, do_wref, do_nkpz
   end if
   !
   if (do_orbdep .and. meta_ionode .and. .not. do_hybrid) then
      !
      write (stdout, 2005) vanishing_rho_w
      if (nknmax > 0) write (stdout, 2030) nknmax
      !
      write (stdout, "(3x, 'NK memusage = ', f10.3, ' MB', /)") &
         nksic_memusage()
   end if
   !
   if (do_orbdep .and. innerloop_until < 1) then
      !
      innerloop_until = -1
      !
   end if
   !
   if (meta_ionode) then
      !
      write (stdout, "()")
      write (stdout, 2006) f_cutoff, do_spinsym
      !
   end if
   !
2000 format(3X, 'NK sic with reference occupation = ', f7.4,/)
2001 format(3X, 'PZ sic = ', l4,/)
2002 format(3X, 'NK sic with integral ref = ', l4,/)
2004 format(3X, 'NK background density factor = ', f7.4, /, &
          3X, 'NK scaling factor = ', f7.4)
2005 format(3X, 'rhothr = ', e8.1)
2006 format(3X, 'f_cutoff = ', f7.4, /, &
          3X, 'do_spinsym   = ', l4)
2010 format(3X, 'NK cross-derivatives = ', l4, /, &
          3X, 'NK reference derivatives = ', l4, /, &
          3X, 'NK on top of PZ = ', l4)
2030 format(3X, 'NK applied up to orbital', i7)

end subroutine nksic_init

!-----------------------------------------------------------------------
subroutine hf_init
!-----------------------------------------------------------------------
!     subroutine introduced by Giovanni Borghi and Andrea Ferretti,
!     following Li, Y. and Dabo, I.
!     Electronic levels and electrical response of
!     periodic molecular structures from plane-wave
!     orbital-dependent calculations.
!     Physical Review B 84, 155127 (2011)
!
   use hfmod, only: do_hf, hfscalfact, allocate_hf
   use nksic, only: f_cutoff
   use input_parameters, only: do_hf_ => do_hf, &
                               which_orbdep_ => which_orbdep, &
                               hfscalfact_ => hfscalfact, &
                               f_cutoff_ => f_cutoff
   use io_global, only: meta_ionode, stdout
   use electrons_base, only: nspin, nbspx
   use gvecw, only: ngw
   use grid_dimensions, only: nnrx
   use funct, only: start_exx, set_dft_from_name, dft_is_hybrid
   !
   implicit none
   !
   logical         :: ishybrid = .false.
   character(30)   :: dft_name

   do_hf = do_hf_
   hfscalfact = hfscalfact_
   f_cutoff = f_cutoff_
   !
   SELECT CASE (TRIM(which_orbdep_))
   CASE ("hf")
      !
      do_hf = .TRUE.
      ishybrid = .TRUE.
      !
   CASE ("b3lyp")
      !
      do_hf = .TRUE.
      hfscalfact = 0.20
      !
   CASE ("pbe0")
      !
      do_hf = .TRUE.
      hfscalfact = 0.25
      !
   END SELECT
   !
   IF (do_hf) ishybrid = .TRUE.
   !
   IF (ishybrid) THEN
      !
      dft_name = TRIM(which_orbdep_)
      IF (LEN_TRIM(dft_name) == 0 .AND. do_hf) dft_name = "hf"
      !
      CALL set_dft_from_name(dft_name)
      !
      IF (meta_ionode) &
         WRITE (stdout, fmt="(/,3X,'Warning XC functionals forced to be: ',A)") dft_name
      !
      CALL start_exx()
      !
   END IF

   !
   !
   if (do_hf .and. meta_ionode) then
      !
      write (stdout, "( 3X,'HF scaling factor = ',f7.4 )") hfscalfact
      write (stdout, "( 3X,'         f_cutoff = ',f7.4 )") f_cutoff
      write (stdout, "( 3X,'    dft_is_hybrid = ',l5 )") dft_is_hybrid()
      !
   end if
   !
   ! allocations
   !
   if (do_hf) call allocate_hf(ngw, nnrx, nspin, nbspx)
   !
end subroutine hf_init

!-----------------------------------------------------------------------
subroutine ee_init
!-----------------------------------------------------------------------
!
   use eecp_mod, only: do_comp, which_compensation, allocate_ee, &
                       tcc_odd
   use input_parameters, only: do_ee, &
                               which_compensation_ => which_compensation, &
                               tcc_odd_ => tcc_odd
   use io_global, only: meta_ionode, stdout
   use grid_dimensions, only: nnrx
   use reciprocal_vectors, only: ngm
   !
   implicit none
   !
   do_comp = do_ee
   which_compensation = which_compensation_
   tcc_odd = tcc_odd_
   if (do_comp .and. meta_ionode) write (stdout, 2010) which_compensation
   !
2010 format(3X, 'EE with periodic-image correction method = ', a20)
   !
   if (do_comp) call allocate_ee(nnrx, ngm)
   !
end subroutine ee_init

!-----------------------------------------------------------------------
subroutine efield_init
!-----------------------------------------------------------------------
!
   use efield_mod, only: do_efield, ampfield, allocate_efield
   use input_parameters, only: do_efield_ => do_efield, &
                               ampfield_ => ampfield
   use io_global, only: meta_ionode, stdout
   use grid_dimensions, only: nnrx
   use reciprocal_vectors, only: ngm
   !
   implicit none
   !
   do_efield = do_efield_
   ampfield = ampfield_
   if (do_efield .and. meta_ionode) write (stdout, 2015) ampfield
2015 format(3X, 'EFIELD with field = ', 3f7.4)
   !
   if (do_efield) call allocate_efield(nnrx, ngm)
   !
end subroutine efield_init
!
!-----------------------------------------------------------------------
integer function set_Hubbard_l(psd) result(hubbard_l)
   !-----------------------------------------------------------------------
   !
   implicit none
   character*3 :: psd
   !
   ! TRANSITION METALS
   !
   if (psd .eq. 'V' .or. psd .eq. 'Cr' .or. psd .eq. 'Mn' .or. psd .eq. 'Fe' .or. &
       psd .eq. 'Co' .or. psd .eq. 'Ni' .or. psd .eq. 'Cu' .or. psd .eq. 'Fe1' .or. &
       psd .eq. 'Fe2') then
      hubbard_l = 2
      !
      ! RARE EARTHS
      !
   elseif (psd .eq. 'Ce') then
      hubbard_l = 3
      !
      ! OTHER ELEMENTS
      !
   elseif (psd .eq. 'H') then
      hubbard_l = 0
   elseif (psd .eq. 'O') then
      hubbard_l = 1
   else
      hubbard_l = -1
      call errore('set_Hubbard_l', 'pseudopotential not yet inserted', 1)
   end if
   return
end function set_Hubbard_l
!
!-----------------------------------------------------------------------
subroutine new_ns_real(c, eigr, betae, hpsi, hpsi_con, forceh)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the on site occupation numbers of the Hubbard ions.
   ! It also calculates the contribution of the Hubbard Hamiltonian to the
   ! electronic potential and to the forces acting on ions.
   !
   use control_flags, ONLY: tfor, tprnfor
   use kinds, ONLY: DP
   use ions_base, only: na, nat, nsp
   use gvecw, only: ngw
   USE uspp, ONLY: nhsa => nkb
   USE uspp_param, ONLY: upf
   use electrons_base, only: nspin, n => nbsp, nx => nbspx, ispin, f
   USE ldaU, ONLY: Hubbard_U, Hubbard_l
   USE ldaU, ONLY: n_atomic_wfc, ns, e_hubbard
   USE cp_interfaces, ONLY: nlsm1, projwfc_hub, s_wfc !added:giovanni
!
   implicit none
#ifdef __PARA
   include 'mpif.h'
#endif
   integer, parameter :: ldmx = 7
   complex(DP), intent(in) :: c(ngw, nx), eigr(ngw, nat), betae(ngw, nhsa)
   complex(DP), intent(out) :: hpsi(ngw, nx), hpsi_con(1, 1)
   real(DP) forceh(3, nat)
   complex(DP), allocatable:: wfc(:, :), swfc(:, :), dphi(:, :, :), spsi(:, :)
   real(DP), allocatable   :: becwfc(:, :), bp(:, :), dbp(:, :, :), wdb(:, :, :)
   real(DP), allocatable   :: dns(:, :, :, :)
   real(DP), allocatable   :: e(:), z(:, :), proj(:, :), temp(:)
   real(DP), allocatable   :: ftemp1(:), ftemp2(:)
   real(DP)                :: lambda(ldmx), somma, ntot, nsum,   &
                              & nsuma, x_value, g_value, step_value
   real(DP) :: f1(ldmx, ldmx), vet(ldmx, ldmx)
   integer is, ia, iat, nb, isp, l, m, m1, m2, k, i, counter, err, ig
   integer iv, jv, inl, jnl, alpha, alpha_a, alpha_s, ipol
   integer, allocatable ::  offset(:, :)
   complex(DP) :: tempsi
   allocate (wfc(ngw, n_atomic_wfc))
   allocate (ftemp1(ldmx))
   allocate (ftemp2(ldmx))
!
! calculate wfc = atomic states
!
!!!      call ewfc(eigr,n_atomic_wfc,wfc)
!
! calculate bec = <beta|wfc>
!
   allocate (becwfc(nhsa, n_atomic_wfc))
!!!      call nlsm1 (n_atomic_wfc,1,nsp,eigr,wfc,becwfc)
!
   allocate (swfc(ngw, n_atomic_wfc))
!!!      call s_wfc(n_atomic_wfc,becwfc,betae,wfc,swfc)
!
! calculate proj = <c|S|wfc>
!
   allocate (proj(n, n_atomic_wfc))
   CALL projwfc_hub(c, nx, eigr, betae, n, n_atomic_wfc,            &
  & wfc, becwfc, swfc, proj) !@@
!
   allocate (offset(nsp, nat))
   counter = 0
   do is = 1, nsp
      do ia = 1, na(is)
         do i = 1, upf(is)%nwfc
            l = upf(is)%lchi(i)
            if (l .eq. Hubbard_l(is)) offset(is, ia) = counter
            counter = counter + 2*l + 1
         end do
      end do
   end do
   if (counter .ne. n_atomic_wfc)                                      &
  &                 call errore('new_ns', 'nstart<>counter', 1)
   ns(:, :, :, :) = 0.d0
   iat = 0
   do is = 1, nsp
      do ia = 1, na(is)
         iat = iat + 1
         if (Hubbard_U(is) .ne. 0.d0) then
            k = offset(is, ia)
            do m1 = 1, 2*Hubbard_l(is) + 1
               do m2 = m1, 2*Hubbard_l(is) + 1
                  do i = 1, n
!                      write(6,*) i,ispin(i),f(i)
                     ns(iat, ispin(i), m1, m2) = ns(iat, ispin(i), m1, m2) + &
    &                               f(i)*proj(i, k + m2)*proj(i, k + m1)
                  end do
!                     ns(iat,:,m2,m1) = ns(iat,:,m1,m2)
                  ns(iat, 1, m2, m1) = ns(iat, 1, m1, m2)
                  ns(iat, 2, m2, m1) = ns(iat, 2, m1, m2)
               end do
            end do
         end if
      end do
   end do
   if (nspin .eq. 1) ns = 0.5d0*ns
! Contributions to total energy
   e_hubbard = 0.d0
   iat = 0
   do is = 1, nsp
      do ia = 1, na(is)
         iat = iat + 1
         if (Hubbard_U(is) .ne. 0.d0) then
            k = offset(is, ia)
            do isp = 1, nspin
               do m1 = 1, 2*Hubbard_l(is) + 1
                  e_hubbard = e_hubbard + 0.5d0*Hubbard_U(is)*    &
  &                           ns(iat, isp, m1, m1)
                  do m2 = 1, 2*Hubbard_l(is) + 1
                     e_hubbard = e_hubbard - 0.5d0*Hubbard_U(is)* &
  &                              ns(iat, isp, m1, m2)*ns(iat, isp, m2, m1)
                  end do
               end do
            end do
         end if
      end do
   end do
   if (nspin .eq. 1) e_hubbard = 2.d0*e_hubbard
!       if (nspin.eq.1) e_lambda = 2.d0*e_lambda
!
!      Calculate the potential and forces on wavefunctions due to U
!
   hpsi(:, :) = CMPLX(0.d0, 0.d0)
   iat = 0
   do is = 1, nsp
      do ia = 1, na(is)
         iat = iat + 1
         if (Hubbard_U(is) .ne. 0.d0) then
            do i = 1, n
               do m1 = 1, 2*Hubbard_l(is) + 1
                  tempsi = proj(i, offset(is, ia) + m1)
                  do m2 = 1, 2*Hubbard_l(is) + 1
                     tempsi = tempsi - 2.d0*ns(iat, ispin(i), m1, m2)*&
  &                                proj(i, offset(is, ia) + m2)
                  end do
                  tempsi = tempsi*Hubbard_U(is)/2.d0*f(i)
                  call ZAXPY(ngw, tempsi, swfc(1, offset(is, ia) + m1), 1, &
  &                           hpsi(1, i), 1)
               end do
            end do
         end if
      end do
   end do
!
!      Calculate the potential and energy due to constraint
!
   hpsi_con(:, :) = 0.d0
!
! Calculate the contribution to forces on ions due to U and constraint
!
   forceh = 0.d0
   if ((tfor) .or. (tprnfor)) then
      allocate (bp(nhsa, n), dbp(nhsa, n, 3), wdb(nhsa, n_atomic_wfc, 3))
      allocate (dns(nat, nspin, ldmx, ldmx))
      allocate (spsi(ngw, n))
!
      call nlsm1(n, 1, nsp, eigr, c, bp)
      call s_wfc(n, bp, betae, c, spsi)
      call nlsm2_repl(ngw, nhsa, n, eigr, c, dbp)
      call nlsm2_repl(ngw, nhsa, n_atomic_wfc, eigr, wfc, wdb)
!
      alpha = 0
      do alpha_s = 1, nsp
         do alpha_a = 1, na(alpha_s)
            alpha = alpha + 1
            do ipol = 1, 3
               call dndtau_real(alpha_a, alpha_s, becwfc, spsi, bp, dbp, wdb,      &
     &                    offset, c, wfc, eigr, betae, proj, ipol, dns)
               iat = 0
               do is = 1, nsp
                  do ia = 1, na(is)
                     iat = iat + 1
                     if (Hubbard_U(is) .ne. 0.d0) then
                        do isp = 1, nspin
                           do m2 = 1, 2*Hubbard_l(is) + 1
                              forceh(ipol, alpha) = forceh(ipol, alpha) -            &
     &                        Hubbard_U(is)*0.5d0*dns(iat, isp, m2, m2)
                              do m1 = 1, 2*Hubbard_l(is) + 1
                                 forceh(ipol, alpha) = forceh(ipol, alpha) +         &
     &                           Hubbard_U(is)*ns(iat, isp, m2, m1)*       &
     &                           dns(iat, isp, m1, m2)
                              end do
                           end do
                        end do
                     end if
! Occupation constraint add here
                  end do
               end do
            end do
         end do
      end do
      if (nspin .eq. 1) then
         forceh = 2.d0*forceh
      end if
!
      deallocate (wfc, becwfc, spsi, proj, offset, swfc, dns, bp, dbp, wdb)
   end if
   return
end subroutine new_ns_real
!
!
!-----------------------------------------------------------------------
subroutine new_ns_twin(c, eigr, betae, hpsi, hpsi_con, forceh, lgam)
!-----------------------------------------------------------------------
!
! This routine computes the on site occupation numbers of the Hubbard ions.
! It also calculates the contribution of the Hubbard Hamiltonian to the
! electronic potential and to the forces acting on ions.
!
   use control_flags, ONLY: tfor, tprnfor
   use kinds, ONLY: DP
   use ions_base, only: na, nat, nsp
   use gvecw, only: ngw
   use reciprocal_vectors, only: ng0 => gstart
   USE uspp, ONLY: nhsa => nkb
   USE uspp_param, ONLY: upf
   use electrons_base, only: nspin, n => nbsp, nx => nbspx, ispin, f
   USE ldaU, ONLY: lda_plus_u, Hubbard_U, Hubbard_l
   USE ldaU, ONLY: n_atomic_wfc, ns, ns_c, e_hubbard
   USE cp_interfaces, ONLY: nlsm1, projwfc_hub, s_wfc !added:giovanni
   USE twin_types
!
   implicit none
#ifdef __PARA
   include 'mpif.h'
#endif
   integer, parameter :: ldmx = 7
   complex(DP), intent(in) :: c(ngw, nx), eigr(ngw, nat),      &
  &                               betae(ngw, nhsa)
   complex(DP), intent(out) :: hpsi(ngw, nx), hpsi_con(1, 1)
   real(DP) forceh(3, nat)
   logical :: lgam
!
   complex(DP), allocatable:: wfc(:, :), swfc(:, :), dphi(:, :, :),   &
  &                               spsi(:, :)
   type(twin_matrix) :: bp
   type(twin_tensor) :: dbp, wdb
   type(twin_matrix)       :: becwfc
   real(DP), allocatable   :: dns(:, :, :, :)
   complex(DP), allocatable   :: dns_c(:, :, :, :)
   real(DP), allocatable   :: e(:), z(:, :),                      &
  &                                temp(:)
   type(twin_matrix)       :: proj
   real(DP), allocatable   :: ftemp1(:), ftemp2(:)
   real(DP)                :: lambda(ldmx), somma, ntot, nsum,   &
  &                           nsuma, x_value, g_value, step_value
   real(DP) :: f1(ldmx, ldmx), vet(ldmx, ldmx)
   integer is, ia, iat, nb, isp, l, m, m1, m2, k, i, counter, err, ig
   integer iv, jv, inl, jnl, alpha, alpha_a, alpha_s, ipol
   integer, allocatable ::  offset(:, :)
   complex(DP) :: tempsi
!
!
   allocate (wfc(ngw, n_atomic_wfc))
   allocate (ftemp1(ldmx))
   allocate (ftemp2(ldmx))
!
! calculate wfc = atomic states
!
!!!      call ewfc(eigr,n_atomic_wfc,wfc)
!
! calculate bec = <beta|wfc>
!
   call init_twin(becwfc, lgam)
   call allocate_twin(becwfc, nhsa, n_atomic_wfc, lgam)
!!!      call nlsm1 (n_atomic_wfc,1,nsp,eigr,wfc,becwfc)
!
   allocate (swfc(ngw, n_atomic_wfc))
!!!      call s_wfc(n_atomic_wfc,becwfc,betae,wfc,swfc)
!
! calculate proj = <c|S|wfc>
!
   call init_twin(proj, lgam)
   call allocate_twin(proj, n, n_atomic_wfc, lgam)

   CALL projwfc_hub(c, nx, eigr, betae, n, n_atomic_wfc,            &
  & wfc, becwfc, swfc, proj, lgam) !@@
!
   allocate (offset(nsp, nat))
   counter = 0
   do is = 1, nsp
      do ia = 1, na(is)
         do i = 1, upf(is)%nwfc
            l = upf(is)%lchi(i)
            if (l .eq. Hubbard_l(is)) offset(is, ia) = counter
            counter = counter + 2*l + 1
         end do
      end do
   end do
   if (counter .ne. n_atomic_wfc)                                      &
  &                 call errore('new_ns', 'nstart<>counter', 1)
   ns(:, :, :, :) = 0.d0
   iat = 0
   do is = 1, nsp
      do ia = 1, na(is)
         iat = iat + 1
         if (Hubbard_U(is) .ne. 0.d0) then
            k = offset(is, ia)
            do m1 = 1, 2*Hubbard_l(is) + 1
               do m2 = m1, 2*Hubbard_l(is) + 1
                  IF (lgam) THEN
                     do i = 1, n
                        !                      write(6,*) i,ispin(i),f(i)
                        ns(iat, ispin(i), m1, m2) = ns(iat, ispin(i), m1, m2) + &
       &                               f(i)*proj%rvec(i, k + m2)*proj%rvec(i, k + m1)
                     end do
!                     ns(iat,:,m2,m1) = ns(iat,:,m1,m2)
                     ns(iat, 1, m2, m1) = ns(iat, 1, m1, m2)
                     ns(iat, 2, m2, m1) = ns(iat, 2, m1, m2)
                  ELSE
                     do i = 1, n
                        !                      write(6,*) i,ispin(i),f(i)
                        ns_c(iat, ispin(i), m1, m2) = ns_c(iat, ispin(i), m1, m2) + &
       &                               f(i)*proj%cvec(i, k + m2)* &
                                       CONJG(proj%cvec(i, k + m1))
                     end do
!                     ns(iat,:,m2,m1) = ns(iat,:,m1,m2)
                     ns_c(iat, 1, m2, m1) = CONJG(ns_c(iat, 1, m1, m2))
                     ns_c(iat, 2, m2, m1) = CONJG(ns_c(iat, 2, m1, m2))
                  END IF
               end do
            end do
         end if
      end do
   end do
   if (nspin .eq. 1) THEN
      IF (lgam) THEN
         ns = 0.5d0*ns
      ELSE
         ns_c = 0.5d0*ns_c
      END IF
   END IF
! Contributions to total energy
   e_hubbard = 0.d0
   iat = 0
   do is = 1, nsp
      do ia = 1, na(is)
         iat = iat + 1
         if (Hubbard_U(is) .ne. 0.d0) then
            k = offset(is, ia)
            do isp = 1, nspin
               IF (lgam) THEN
                  do m1 = 1, 2*Hubbard_l(is) + 1
                     e_hubbard = e_hubbard + 0.5d0*Hubbard_U(is)*    &
     &                           ns(iat, isp, m1, m1)
                     do m2 = 1, 2*Hubbard_l(is) + 1
                        e_hubbard = e_hubbard - 0.5d0*Hubbard_U(is)* &
     &                              ns(iat, isp, m1, m2)*ns(iat, isp, m2, m1)
                     end do
                  end do
               ELSE
                  do m1 = 1, 2*Hubbard_l(is) + 1
                     e_hubbard = e_hubbard + 0.5d0*Hubbard_U(is)*    &
     &                           DBLE(ns_c(iat, isp, m1, m1))
                     do m2 = 1, 2*Hubbard_l(is) + 1
                        e_hubbard = e_hubbard - 0.5d0*Hubbard_U(is)* &
     &                             DBLE(ns_c(iat, isp, m1, m2)*ns_c(iat, isp, m2, m1))
                     end do
                  end do
               END IF
            end do
         end if
      end do
   end do
   if (nspin .eq. 1) e_hubbard = 2.d0*e_hubbard
!       if (nspin.eq.1) e_lambda = 2.d0*e_lambda
!
!      Calculate the potential and forces on wavefunctions due to U
!
   hpsi(:, :) = CMPLX(0.d0, 0.d0)
   iat = 0
   do is = 1, nsp
      do ia = 1, na(is)
         iat = iat + 1
         if (Hubbard_U(is) .ne. 0.d0) then
            IF (lgam) THEN
               do i = 1, n
                  do m1 = 1, 2*Hubbard_l(is) + 1
                     tempsi = proj%rvec(i, offset(is, ia) + m1)
                     do m2 = 1, 2*Hubbard_l(is) + 1
                        tempsi = tempsi - 2.d0*ns(iat, ispin(i), m1, m2)*&
     &                                proj%rvec(i, offset(is, ia) + m2)
                     end do
                     tempsi = tempsi*Hubbard_U(is)/2.d0*f(i)
                     call ZAXPY(ngw, tempsi, swfc(1, offset(is, ia) + m1), 1, &
     &                           hpsi(1, i), 1)
                  end do
               end do
            ELSE
               do i = 1, n
                  do m1 = 1, 2*Hubbard_l(is) + 1
                     tempsi = proj%cvec(i, offset(is, ia) + m1)
                     do m2 = 1, 2*Hubbard_l(is) + 1
                        tempsi = tempsi - 2.d0*ns_c(iat, ispin(i), m1, m2)*&
     &                                proj%cvec(i, offset(is, ia) + m2)
                     end do
                     tempsi = tempsi*Hubbard_U(is)/2.d0*f(i)
                     call ZAXPY(ngw, tempsi, swfc(1, offset(is, ia) + m1), 1, &
     &                           hpsi(1, i), 1)
                  end do
               end do
            END IF
         end if
      end do
   end do
!
!      Calculate the potential and energy due to constraint
!
   hpsi_con(:, :) = 0.d0
!
! Calculate the contribution to forces on ions due to U and constraint
!
   forceh = 0.d0
   if ((tfor) .or. (tprnfor)) then
      call init_twin(bp, lgam)
      call allocate_twin(bp, nhsa, n, lgam)
      call init_twin(dbp, lgam)
      call allocate_twin(dbp, nhsa, n, 3, lgam)
      call init_twin(wdb, lgam)
      call allocate_twin(wdb, nhsa, n_atomic_wfc, 3, lgam)

      IF (lgam) THEN
         allocate (dns(nat, nspin, ldmx, ldmx))
      ELSE
         allocate (dns_c(nat, nspin, ldmx, ldmx))
      END IF
      !
      allocate (spsi(ngw, n))
!
      call nlsm1(n, 1, nsp, eigr, c, bp, 1, lgam)
      call s_wfc(n, bp, betae, c, spsi, lgam)
      call nlsm2_repl(ngw, nhsa, n, eigr, c, dbp)
      call nlsm2_repl(ngw, nhsa, n_atomic_wfc, eigr, wfc, wdb)
!
      alpha = 0
      do alpha_s = 1, nsp
         do alpha_a = 1, na(alpha_s)
            alpha = alpha + 1
            do ipol = 1, 3
               IF (lgam) THEN
                  call dndtau_real(alpha_a, alpha_s, becwfc%rvec, spsi, bp%rvec, dbp%rvec, wdb%rvec,      &
     &                    offset, c, wfc, eigr, betae, proj%rvec, ipol, dns)
               ELSE
                  call dndtau_cmplx(alpha_a, alpha_s, becwfc%cvec, spsi, bp%cvec, dbp%cvec, wdb%cvec,      &
     &                    offset, c, wfc, eigr, betae, proj%cvec, ipol, dns_c)
               END IF
               iat = 0
               do is = 1, nsp
                  do ia = 1, na(is)
                     iat = iat + 1
                     if (Hubbard_U(is) .ne. 0.d0) then
                        do isp = 1, nspin
                           IF (lgam) THEN
                              do m2 = 1, 2*Hubbard_l(is) + 1
                                 forceh(ipol, alpha) = forceh(ipol, alpha) -            &
        &                        Hubbard_U(is)*0.5d0*dns(iat, isp, m2, m2)
                                 do m1 = 1, 2*Hubbard_l(is) + 1
                                    forceh(ipol, alpha) = forceh(ipol, alpha) +         &
        &                           Hubbard_U(is)*ns(iat, isp, m2, m1)*       &
        &                           dns(iat, isp, m1, m2)
                                 end do
                              end do
                           ELSE
                              do m2 = 1, 2*Hubbard_l(is) + 1
                                 forceh(ipol, alpha) = forceh(ipol, alpha) -            &
        &                        Hubbard_U(is)*0.5d0*DBLE(dns_c(iat, isp, m2, m2))
                                 do m1 = 1, 2*Hubbard_l(is) + 1
                                    forceh(ipol, alpha) = forceh(ipol, alpha) +         &
        &                           Hubbard_U(is)*DBLE(ns_c(iat, isp, m2, m1)*       &
        &                           dns_c(iat, isp, m1, m2))
                                 end do
                              end do
                           END IF
                        end do
                     end if
! Occupation constraint add here
                  end do
               end do
            end do
         end do
      end do
      if (nspin .eq. 1) then
         forceh = 2.d0*forceh
      end if
!
      deallocate (wfc, spsi, offset, swfc)

      call deallocate_twin(becwfc)
      call deallocate_twin(proj)
      call deallocate_twin(bp)
      call deallocate_twin(dbp)
      call deallocate_twin(wdb)

      IF (allocated(dns)) deallocate (dns)
      IF (allocated(dns_c)) deallocate (dns_c)
      !
   end if
   return
end subroutine new_ns_twin
!

!-----------------------------------------------------------------------
subroutine write_ns
!-----------------------------------------------------------------------
!
! This routine computes the occupation numbers on atomic orbitals.
! It also write the occupation number in the output file.
!
   USE kinds, only: DP
   USE constants, ONLY: autoev
   use electrons_base, only: nspin
   use electrons_base, only: n => nbsp
   use ions_base, only: na, nat, nsp
   use gvecw, only: ngw
   USE ldaU, ONLY: lda_plus_u, Hubbard_U, Hubbard_l
   USE ldaU, ONLY: n_atomic_wfc, ns, e_hubbard
   USE ldaU, ONLY: Hubbard_lmax
   use dspev_module, only: dspev_drv

   implicit none

   integer :: is, isp, ia, m1, m2, ldim, iat, err, k
! cpunter on atoms type
! counter on spin component
! counter on atoms
! counter on wavefn
! counters on d components
   integer, parameter :: ldmx = 7
   real(DP), allocatable   :: ftemp1(:), ftemp2(:)
   real(DP) :: f1(ldmx*ldmx), vet(ldmx, ldmx)
   real(DP) :: lambda(ldmx), nsum, nsuma
   write (*, *) 'enter write_ns'

   if (2*Hubbard_lmax + 1 .gt. ldmx) &
      call errore('write_ns', 'ldmx is too small', 1)

!  if (step_con) then
!     do isp=1,nspin
!        write (6,'(6(a,i2,a,i2,a,f8.4,6x))') &
!        ('A_con(',is,',',isp,') =', A_con(is,isp),is=1,nsp)
!     enddo
!     write (6,'(6(a,i2,a,f8.4,6x))') &
!           ('sigma_con(',is,') =', sigma_con(is), is=1,nsp)
!     write (6,'(6(a,i2,a,f8.4,6x))') &
!        ('alpha_con(',is,') =', alpha_con(is), is=1,nsp)
!  endif
   write (6, '(6(a,i2,a,f8.4,6x))') &
      ('U(', is, ') =', Hubbard_U(is)*autoev, is=1, nsp)
!  write (6,'(6(a,i2,a,f8.4,6x))') &
!        ('alpha(',is,') =', Hubbard_alpha(is) * autoev, is=1,nsp)
   nsum = 0.d0
   allocate (ftemp1(ldmx))
   allocate (ftemp2(ldmx))
   iat = 0
   write (6, *) 'nsp', nsp
   do is = 1, nsp
      do ia = 1, na(is)
         nsuma = 0.d0
         iat = iat + 1
!        if (iat.eq.1) then
         if (Hubbard_U(is) .ne. 0.d0) then
            do isp = 1, nspin
               do m1 = 1, 2*Hubbard_l(is) + 1
                  nsuma = nsuma + ns(iat, isp, m1, m1)
               end do
            end do
            if (nspin .eq. 1) nsuma = 2.d0*nsuma
            write (6, '(a,x,i2,2x,a,f11.7)') 'atom', iat,              &
  &                                      ' Tr[ns(na)]= ', nsuma
            nsum = nsum + nsuma
!
            do isp = 1, nspin

               k = 0
               do m1 = 1, 2*Hubbard_l(is) + 1
                  do m2 = m1, 2*Hubbard_l(is) + 1
                     k = k + 1
                     f1(k) = ns(iat, isp, m2, m1)
                  end do
               end do

               CALL dspev_drv('V', 'L', 2*Hubbard_l(is) + 1, f1, lambda, vet, ldmx)

               write (6, '(a,x,i2,2x,a,x,i2)') 'atom', iat, 'spin', isp
               write (6, '(a,7f10.7)') 'eigenvalues: ', (lambda(m1), m1=1,&
  &                                2*Hubbard_l(is) + 1)
               write (6, *) 'eigenvectors'
               do m2 = 1, 2*Hubbard_l(is) + 1
                  write (6, '(i2,2x,7(f10.7,x))') m2, (real(vet(m1, m2)),&
  &                            m1=1, 2*Hubbard_l(is) + 1)
               end do
               write (6, *) 'occupations'
               do m1 = 1, 2*Hubbard_l(is) + 1
                  write (6, '(7(f6.3,x))') (ns(iat, isp, m1, m2), m2=1,    &
  &                     2*Hubbard_l(is) + 1)
               end do
            end do
         end if
!        end if
      end do
   end do
   deallocate (ftemp1, ftemp2)
   return
end subroutine write_ns
!-----------------------------------------------------------------------
subroutine genatwfc(n_atomic_wfc, atwfc)
!-----------------------------------------------------------------------
!
! Compute atomic wavefunctions in G-space, in the same order as used in new_ns
!
   use ions_base, only: na, nsp
   use gvecw, only: ngw
   use reciprocal_vectors, only: g, gx, ng0 => gstart
   use cell_base, only: omega, tpiba
   use constants, only: fpi
   USE atom, ONLY: rgrid
   USE uspp_param, ONLY: upf
   USE kinds, ONLY: DP
!
   implicit none
   integer, intent(in) :: n_atomic_wfc
   complex(DP), intent(out):: atwfc(ngw, n_atomic_wfc)
!
   integer natwfc, is, ia, ir, nb, l, m, lm, i, lmax_wfc, ig
   real(DP), allocatable::  ylm(:, :), q(:), jl(:), vchi(:),        &
  &     chiq(:), gxn(:, :)
!
   IF (.NOT. ALLOCATED(rgrid)) &
      CALL errore(' genatwfc ', ' rgrid not allocated ', 1)
!
   allocate (q(ngw))
   allocate (gxn(3, ngw))
   allocate (chiq(ngw))
!
   do ig = 1, ngw
      q(ig) = sqrt(g(ig))*tpiba
   end do
   if (ng0 .eq. 2) gxn(1, :) = 0.0d0
   do ig = ng0, ngw
      gxn(:, ig) = gx(:, ig)/sqrt(g(ig)) !ik<=>ig
   end do
!
   natwfc = 0
!@@@@@
!
! calculate max angular momentum required in wavefunctions
!
   lmax_wfc = -1
   DO is = 1, nsp
      lmax_wfc = MAX(lmax_wfc, MAXVAL(upf(is)%lchi(1:upf(is)%nwfc)))
   END DO
   !
   ALLOCATE (ylm(ngw, (lmax_wfc + 1)**2))
   !
   CALL ylmr2((lmax_wfc + 1)**2, ngw, gx, g, ylm)
!@@@@@

   do is = 1, nsp
      ALLOCATE (jl(rgrid(is)%mesh), vchi(rgrid(is)%mesh))
      do ia = 1, na(is)
!
!   radial fourier transform of the chi functions
!   NOTA BENE: chi is r times the radial part of the atomic wavefunction
!              bess requires l+1, not l, on input
!
         do nb = 1, upf(is)%nwfc
            l = upf(is)%lchi(nb)
            do i = 1, ngw
               call sph_bes(rgrid(is)%mesh, rgrid(is)%r, q(i), l, jl)
               do ir = 1, rgrid(is)%mesh
                  vchi(ir) = upf(is)%chi(ir, nb)*rgrid(is)%r(ir)*jl(ir)
               end do
               call simpson_cp90(rgrid(is)%mesh, vchi, rgrid(is)%rab, chiq(i))
            end do
!
!   multiply by angular part and structure factor
!   NOTA BENE: the factor i^l MUST be present!!!
!
            do m = 1, 2*l + 1
               lm = l**2 + m
!                  call ylmr2b(lm,ngw,ngw,gxn,ylm)
               natwfc = natwfc + 1
               atwfc(:, natwfc) = CMPLX(0.d0, 1.d0)**l*ylm(:, lm)*chiq(:)
            end do
         end do
      end do
      DEALLOCATE (vchi, jl)
   end do
!
   do i = 1, natwfc
      call DSCAL(2*ngw, fpi/sqrt(omega), atwfc(1, i), 1)
   end do
!
   if (natwfc .ne. n_atomic_wfc)                                       &
  &     call errore('atomic_wfc', 'unexpected error', natwfc)
!
   deallocate (ylm)
   deallocate (chiq)
   deallocate (gxn)
   deallocate (q)
!
   return
end subroutine genatwfc
!
!-------------------------------------------------------------------------
subroutine dndtau_real(alpha_a, alpha_s, becwfc, spsi, bp, dbp, wdb,         &
&                  offset, c, wfc,                                   &
&                  eigr, betae,                                     &
&                  proj, ipol, dns)
!-----------------------------------------------------------------------
!
! This routine computes the derivative of the ns with respect to the ionic
! displacement tau(alpha,ipol) used to obtain the Hubbard contribution to the
! atomic forces.
!
   use ions_base, only: na, nat, nsp
   use gvecw, only: ngw
   use electrons_base, only: nspin, n => nbsp, nx => nbspx, ispin, f
   USE uspp, ONLY: nhsa => nkb
   USE ldaU, ONLY: Hubbard_U, Hubbard_l
   USE ldaU, ONLY: n_atomic_wfc, ns
   USE kinds, ONLY: DP
!
   implicit none
   integer, parameter :: ldmx = 7
   integer ibnd, is, i, ia, counter, m1, m2, l, iat, alpha, ldim
! input
   integer, intent(in) :: offset(nsp, nat)
   integer, intent(in) :: alpha_a, alpha_s, ipol
   real(DP), intent(in) :: wfc(2, ngw, n_atomic_wfc), c(2, ngw, nx),  &
  &                            eigr(2, ngw, nat), betae(2, ngw, nhsa),    &
  &                            becwfc(nhsa, n_atomic_wfc),            &
  &                            bp(nhsa, n), dbp(nhsa, n, 3), wdb(nhsa, n_atomic_wfc, 3)
   real(DP), intent(in) :: proj(n, n_atomic_wfc)
   complex(DP), intent(in) :: spsi(ngw, n)
! output
   real(DP), intent(out) :: dns(nat, nspin, ldmx, ldmx)
!
!     dns !derivative of ns(:,:,:,:) w.r.t. tau
!
   real(DP), allocatable :: dproj(:, :)
!
!     dproj(n,n_atomic_wfc) ! derivative of proj(:,:) w.r.t. tau
!
   allocate (dproj(n, n_atomic_wfc))
!
   dns(:, :, :, :) = 0.d0
!
   call dprojdtau_real(c, wfc, becwfc, spsi, bp, dbp, wdb, eigr, alpha_a,     &
&                   alpha_s, ipol, offset(alpha_s, alpha_a), dproj)
!
! compute the derivative of occupation numbers (the quantities dn(m1,m2))
! of the atomic orbitals. They are real quantities as well as n(m1,m2)
!
   iat = 0
   do is = 1, nsp
      do ia = 1, na(is)
         iat = iat + 1
         if (Hubbard_U(is) .ne. 0.d0) then
            ldim = 2*Hubbard_l(is) + 1
            do m1 = 1, ldim
               do m2 = m1, ldim
                  do ibnd = 1, n
                     dns(iat, ispin(ibnd), m1, m2) =                    &
  &                  dns(iat, ispin(ibnd), m1, m2) +                    &
  &                   f(ibnd)*REAL(proj(ibnd, offset(is, ia) + m1)*   &
  &                   (dproj(ibnd, offset(is, ia) + m2)) +              &
  &                         dproj(ibnd, offset(is, ia) + m1)*          &
  &                         (proj(ibnd, offset(is, ia) + m2)))
                  end do
                  dns(iat, :, m2, m1) = dns(iat, :, m1, m2)
               end do
            end do
         end if
      end do
   end do
!
   deallocate (dproj)
   return
end subroutine dndtau_real
!
!
!-------------------------------------------------------------------------
subroutine dndtau_cmplx(alpha_a, alpha_s, becwfc, spsi, bp, dbp, wdb,         &
&                  offset, c, wfc,                                   &
&                  eigr, betae,                                     &
&                  proj, ipol, dns_c)
!-----------------------------------------------------------------------
!
! This routine computes the derivative of the ns with respect to the ionic
! displacement tau(alpha,ipol) used to obtain the Hubbard contribution to the
! atomic forces.
!
   use ions_base, only: na, nat, nsp
   use gvecw, only: ngw
   use electrons_base, only: nspin, n => nbsp, nx => nbspx, ispin, f
   USE uspp, ONLY: nhsa => nkb
   USE ldaU, ONLY: Hubbard_U, Hubbard_l
   USE ldaU, ONLY: n_atomic_wfc, ns_c
   USE kinds, ONLY: DP
!
   implicit none
   integer, parameter :: ldmx = 7
   integer ibnd, is, i, ia, counter, m1, m2, l, iat, alpha, ldim
! input
   integer, intent(in) :: offset(nsp, nat)
   integer, intent(in) :: alpha_a, alpha_s, ipol
   complex(DP), intent(in) :: wfc(ngw, n_atomic_wfc), c(ngw, nx),  &
  &                            eigr(ngw, nat), betae(ngw, nhsa),    &
  &                            becwfc(nhsa, n_atomic_wfc),            &
  &                            bp(nhsa, n), dbp(nhsa, n, 3), wdb(nhsa, n_atomic_wfc, 3)
   complex(DP), intent(in) :: proj(n, n_atomic_wfc)
   complex(DP), intent(in) :: spsi(ngw, n)
! output
   complex(DP), intent(out) :: dns_c(nat, nspin, ldmx, ldmx)
!
!     dns !derivative of ns(:,:,:,:) w.r.t. tau
!
   complex(DP), allocatable :: dproj(:, :)
!
!     dproj(n,n_atomic_wfc) ! derivative of proj(:,:) w.r.t. tau
!
   allocate (dproj(n, n_atomic_wfc))
!
   dns_c(:, :, :, :) = 0.d0
!
   call dprojdtau_cmplx(c, wfc, becwfc, spsi, bp, dbp, wdb, eigr, alpha_a,     &
&                   alpha_s, ipol, offset(alpha_s, alpha_a), dproj)
!
! compute the derivative of occupation numbers (the quantities dn(m1,m2))
! of the atomic orbitals. They are real quantities as well as n(m1,m2)
!
   iat = 0
   do is = 1, nsp
      do ia = 1, na(is)
         iat = iat + 1
         if (Hubbard_U(is) .ne. 0.d0) then
            ldim = 2*Hubbard_l(is) + 1
            do m1 = 1, ldim
               do m2 = m1, ldim
                  do ibnd = 1, n
                     dns_c(iat, ispin(ibnd), m1, m2) =                    &
  &                  dns_c(iat, ispin(ibnd), m1, m2) +                    &
  &                   f(ibnd)*(proj(ibnd, offset(is, ia) + m1)*   &
  &                   (dproj(ibnd, offset(is, ia) + m2)) +              &
  &                         dproj(ibnd, offset(is, ia) + m1)*          &
  &                         (proj(ibnd, offset(is, ia) + m2)))
                  end do
                  dns_c(iat, :, m2, m1) = CONJG(dns_c(iat, :, m1, m2))
               end do
            end do
         end if
      end do
   end do
!
   deallocate (dproj)
   return
end subroutine dndtau_cmplx
!
!
!-----------------------------------------------------------------------
subroutine dprojdtau_real(c, wfc, becwfc, spsi, bp, dbp, wdb, eigr, alpha_a,    &
&                     alpha_s, ipol, offset, dproj)
!-----------------------------------------------------------------------
!
! This routine computes the first derivative of the projection
! <\fi^{at}_{I,m1}|S|\psi_{k,v,s}> with respect to the atomic displacement
! u(alpha,ipol) (we remember that ns_{I,s,m1,m2} = \sum_{k,v}
! f_{kv} <\fi^{at}_{I,m1}|S|\psi_{k,v,s}><\psi_{k,v,s}|S|\fi^{at}_{I,m2}>)
!
   use ions_base, only: na, nat
   use gvecw, only: ngw
   use reciprocal_vectors, only: g, gx, ng0 => gstart
   use electrons_base, only: n => nbsp, nx => nbspx
!      use gvec
!      use constants
   USE uspp, ONLY: nhsa => nkb, qq
   use cvan, ONLY: ish
   USE ldaU, ONLY: Hubbard_U, Hubbard_l
   USE ldaU, ONLY: n_atomic_wfc
   use cell_base, ONLY: tpiba
   USE uspp_param, only: nh !@@@@
   use mp_global, only: intra_image_comm
   use mp, only: mp_sum
   USE kinds, ONLY: DP
!
   implicit none
   integer, parameter :: ldmx = 7
   integer alpha_a, alpha_s, ipol, offset
! input: the displaced atom
! input: the component of displacement
! input: the offset of the wfcs of the atom "alpha_a,alpha_s"
   complex(DP), intent(in) :: spsi(ngw, n),                     &
 &                  c(ngw, nx), eigr(ngw, nat)
! input: the atomic wfc
! input: S|evc>
   real(DP), intent(in) ::becwfc(nhsa, n_atomic_wfc),            &
 &                            wfc(2, ngw, n_atomic_wfc),              &
 &            bp(nhsa, n), dbp(nhsa, n, 3), wdb(nhsa, n_atomic_wfc, 3)
   real(DP), intent(out) :: dproj(n, n_atomic_wfc)
! output: the derivative of the projection
!
   integer i, ig, m1, ibnd, iwf, ia, is, iv, jv, ldim, alpha, l, m, k, inl
!
   real(DP) a1, a2
   real(kind=8), allocatable :: gk(:)
!
   complex(DP), allocatable :: dwfc(:, :)
   real(DP), allocatable :: betapsi(:, :),                       &
  &                              dbetapsi(:, :),                      &
  &                              wfcbeta(:, :), wfcdbeta(:, :), temp(:)
!      dwfc(ngw,ldmx),             ! the derivative of the atomic d wfc
!      betapsi(nh,n),              ! <beta|evc>
!      dbetapsi(nh,n),             ! <dbeta|evc>
!      wfcbeta(n_atomic_wfc,nh),   ! <wfc|beta>
!      wfcdbeta(n_atomic_wfc,nh),  ! <wfc|dbeta>
   ldim = 2*Hubbard_l(alpha_s) + 1
   allocate (dwfc(ngw, ldmx), betapsi(nh(alpha_s), n))
   allocate (dbetapsi(nh(alpha_s), n),                               &
  &           wfcbeta(n_atomic_wfc, nh(alpha_s)))
   allocate (wfcdbeta(n_atomic_wfc, nh(alpha_s)))
   dproj(:, :) = 0.d0
!
! At first the derivative of the atomic wfc is computed
!
!
   allocate (gk(ngw))
   allocate (temp(ngw))
!
   if (Hubbard_U(alpha_s) .ne. 0.d0) then
!
      do ig = 1, ngw
         gk(ig) = gx(ipol, ig)*tpiba
!
         do m1 = 1, ldim
            dwfc(ig, m1) = cmplx(gk(ig)*wfc(2, ig, offset + m1),      &
&                  -1*gk(ig)*wfc(1, ig, offset + m1))
         end do
      end do
!
      do ibnd = 1, n
         do m1 = 1, ldim
            temp(:) = real(conjg(dwfc(:, m1))*spsi(:, ibnd))
            dproj(ibnd, offset + m1) = 2.d0*SUM(temp)
            if (ng0 .eq. 2) dproj(ibnd, offset + m1) = dproj(ibnd, offset + m1) - temp(1)
         end do
      end do
      call mp_sum(dproj, intra_image_comm)
   end if
   do iv = 1, nh(alpha_s)
      inl = ish(alpha_s) + (iv - 1)*na(alpha_s) + alpha_a
      do i = 1, n
         betapsi(iv, i) = bp(inl, i)
         dbetapsi(iv, i) = dbp(inl, i, ipol)
      end do
      do m = 1, n_atomic_wfc
!                 do m1=1,2**Hubbard_l(is) + 1
         wfcbeta(m, iv) = becwfc(inl, m)
         wfcdbeta(m, iv) = wdb(inl, m, ipol)
      end do
   end do
   do ibnd = 1, n
      do iv = 1, nh(alpha_s)
         do jv = 1, nh(alpha_s)
            do m = 1, n_atomic_wfc
!                       do m1=1,2**Hubbard_l(is) + 1
               dproj(ibnd, m) =                                       &
  &                        dproj(ibnd, m) + qq(iv, jv, alpha_s)*       &
  &                         (wfcdbeta(m, iv)*betapsi(jv, ibnd) +      &
  &                           wfcbeta(m, iv)*dbetapsi(jv, ibnd))
            end do
         end do
      end do
   end do
   deallocate (temp, gk)
   deallocate (betapsi)
   deallocate (dwfc)
   deallocate (dbetapsi)
   deallocate (wfcbeta)
   deallocate (wfcdbeta)
   return
end subroutine dprojdtau_real
!
!-----------------------------------------------------------------------
subroutine dprojdtau_cmplx(c, wfc, becwfc, spsi, bp, dbp, wdb, eigr, alpha_a,    &
&                     alpha_s, ipol, offset, dproj)
!-----------------------------------------------------------------------
!
! This routine computes the first derivative of the projection
! <\fi^{at}_{I,m1}|S|\psi_{k,v,s}> with respect to the atomic displacement
! u(alpha,ipol) (we remember that ns_{I,s,m1,m2} = \sum_{k,v}
! f_{kv} <\fi^{at}_{I,m1}|S|\psi_{k,v,s}><\psi_{k,v,s}|S|\fi^{at}_{I,m2}>)
!
   use ions_base, only: na, nat
   use gvecw, only: ngw
   use reciprocal_vectors, only: g, gx, ng0 => gstart
   use electrons_base, only: n => nbsp, nx => nbspx
!      use gvec
!      use constants
   USE uspp, ONLY: nhsa => nkb, qq
   use cvan, ONLY: ish
   USE ldaU, ONLY: Hubbard_U, Hubbard_l
   USE ldaU, ONLY: n_atomic_wfc
   use cell_base, ONLY: tpiba
   USE uspp_param, only: nh !@@@@
   use mp_global, only: intra_image_comm
   use mp, only: mp_sum
   USE kinds, ONLY: DP
!
   implicit none
   integer, parameter :: ldmx = 7
   integer alpha_a, alpha_s, ipol, offset
! input: the displaced atom
! input: the component of displacement
! input: the offset of the wfcs of the atom "alpha_a,alpha_s"
   complex(DP), intent(in) :: spsi(ngw, n),                     &
 &                  c(ngw, nx), eigr(ngw, nat)
! input: the atomic wfc
! input: S|evc>
   complex(DP), intent(in) ::becwfc(nhsa, n_atomic_wfc),            &
 &                            wfc(ngw, n_atomic_wfc),              &
 &            bp(nhsa, n), dbp(nhsa, n, 3), wdb(nhsa, n_atomic_wfc, 3)
   complex(DP), intent(out) :: dproj(n, n_atomic_wfc)
! output: the derivative of the projection
!
   integer i, ig, m1, ibnd, iwf, ia, is, iv, jv, ldim, alpha, l, m, k, inl
!
   real(DP) a1, a2
   real(kind=8), allocatable :: gk(:)
!
   complex(DP), allocatable :: dwfc(:, :)
   complex(DP), allocatable :: betapsi(:, :),                       &
  &                              dbetapsi(:, :),                      &
  &                              wfcbeta(:, :), wfcdbeta(:, :), temp(:)
!      dwfc(ngw,ldmx),             ! the derivative of the atomic d wfc
!      betapsi(nh,n),              ! <beta|evc>
!      dbetapsi(nh,n),             ! <dbeta|evc>
!      wfcbeta(n_atomic_wfc,nh),   ! <wfc|beta>
!      wfcdbeta(n_atomic_wfc,nh),  ! <wfc|dbeta>
   ldim = 2*Hubbard_l(alpha_s) + 1
   allocate (dwfc(ngw, ldmx), betapsi(nh(alpha_s), n))
   allocate (dbetapsi(nh(alpha_s), n),                               &
  &           wfcbeta(n_atomic_wfc, nh(alpha_s)))
   allocate (wfcdbeta(n_atomic_wfc, nh(alpha_s)))
   dproj(:, :) = 0.d0
!
! At first the derivative of the atomic wfc is computed
!
!
   allocate (gk(ngw))
   allocate (temp(ngw))
!
   if (Hubbard_U(alpha_s) .ne. 0.d0) then
!
      do ig = 1, ngw
         gk(ig) = gx(ipol, ig)*tpiba
!
         do m1 = 1, ldim
            dwfc(ig, m1) = (0.d0, -1.d0)*gk(ig)*wfc(ig, offset + m1)
         end do
      end do
!
      do ibnd = 1, n
         do m1 = 1, ldim
            temp(:) = conjg(dwfc(:, m1))*spsi(:, ibnd)
            dproj(ibnd, offset + m1) = SUM(temp)
         end do
      end do
      call mp_sum(dproj, intra_image_comm)
   end if
   do iv = 1, nh(alpha_s)
      inl = ish(alpha_s) + (iv - 1)*na(alpha_s) + alpha_a
      do i = 1, n
         betapsi(iv, i) = bp(inl, i)
         dbetapsi(iv, i) = dbp(inl, i, ipol)
      end do
      do m = 1, n_atomic_wfc
!                 do m1=1,2**Hubbard_l(is) + 1
         wfcbeta(m, iv) = becwfc(inl, m)
         wfcdbeta(m, iv) = wdb(inl, m, ipol)
      end do
   end do
   do ibnd = 1, n
      do iv = 1, nh(alpha_s)
         do jv = 1, nh(alpha_s)
            do m = 1, n_atomic_wfc
!                       do m1=1,2**Hubbard_l(is) + 1
               dproj(ibnd, m) =                                       &
  &                        dproj(ibnd, m) + qq(iv, jv, alpha_s)*       &
  &                         (wfcdbeta(m, iv)*betapsi(jv, ibnd) +      &
  &                           wfcbeta(m, iv)*dbetapsi(jv, ibnd))
            end do
         end do
      end do
   end do
   deallocate (temp, gk)
   deallocate (betapsi)
   deallocate (dwfc)
   deallocate (dbetapsi)
   deallocate (wfcbeta)
   deallocate (wfcdbeta)
   return
end subroutine dprojdtau_cmplx
!
!
!-----------------------------------------------------------------------
subroutine stepfn(A, sigma, x_value, g_value, step_value)
!-----------------------------------------------------------------------
!     This subroutine calculates the value of the gaussian and step
!     functions with a given x_value. A and sigma are given in the
!     input file. ... to be used in occupation_constraint...
!
   USE constants, ONLY: pi
   implicit none
   real(kind=8) A, sigma, x_value, g_value, step_value
   real(kind=8) x
   integer i
   step_value = 0.0d0
   g_value = 0.0d0
!
   do i = 1, 100000
      x = x_value + (i - 100000)/100000.0d0*(x_value + 5.d0*sigma)
!
! Integrate from 5 sigma before the x_value
!
      g_value = A*dexp(-x*x/(2*sigma*sigma))/(sigma*dsqrt(2*pi))
!         write(6,*) 'step', step_value,'g',g_value
!         if (g_value.le.0.0) g_value=0.0
      if ((x_value + 5*sigma) .ge. 0.0d0) then
         step_value = step_value + g_value/100000.0d0*(x_value + 5.d0*sigma)
      end if
   end do
   return
end subroutine stepfn
!
!-----------------------------------------------------------------------
SUBROUTINE projwfc_hub_real(c, nx, eigr, betae, n, n_atomic_wfc,  &
& wfc, becwfc, swfc, proj)
!-----------------------------------------------------------------------
   !
   ! Projection on atomic wavefunctions
   ! Atomic wavefunctions are not orthogonized
   !
   USE kinds, ONLY: DP
   USE constants, ONLY: autoev
   USE io_global, ONLY: stdout
   USE mp_global, ONLY: intra_image_comm
   USE mp, ONLY: mp_sum
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: gstart
   USE ions_base, ONLY: nsp, na, nat
   USE uspp, ONLY: nhsa => nkb
   USE cp_interfaces, ONLY: nlsm1, projwfc_hub, s_wfc !added:giovanni
!
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nx, n, n_atomic_wfc
   COMPLEX(DP), INTENT(IN) :: c(ngw, nx), eigr(ngw, nat), betae(ngw, nhsa)
!
   COMPLEX(DP), INTENT(OUT):: wfc(ngw, n_atomic_wfc),    &
  & swfc(ngw, n_atomic_wfc)
   real(DP), intent(out):: becwfc(nhsa, n_atomic_wfc) !DEBUG
   REAL(DP), ALLOCATABLE :: overlap(:, :), e(:), z(:, :)
   REAL(DP), ALLOCATABLE :: temp(:)
   REAL(DP)                 :: somma, proj(n, n_atomic_wfc)
   INTEGER :: is, ia, nb, l, m, k, i
   !
   ! calculate number of atomic states
   !
   !
   IF (n_atomic_wfc .EQ. 0) RETURN
   !
   !
   ! calculate wfc = atomic states
   !
   CALL atomic_wfc_northo(eigr, n_atomic_wfc, wfc)
   !
   ! calculate bec = <beta|wfc>
   !
   CALL nlsm1(n_atomic_wfc, 1, nsp, eigr, wfc, becwfc)
   !
   ! calculate swfc = S|wfc>
   !
   CALL s_wfc(n_atomic_wfc, becwfc, betae, wfc, swfc)
   !
   ! calculate proj = <c|S|wfc>
   !
   ALLOCATE (temp(ngw))
   DO m = 1, n
      DO l = 1, n_atomic_wfc
         temp(:) = DBLE(CONJG(c(:, m))*swfc(:, l)) !@@@@
         proj(m, l) = 2.d0*SUM(temp)
         IF (gstart == 2) proj(m, l) = proj(m, l) - temp(1)
      END DO
   END DO
   DEALLOCATE (temp)
   CALL mp_sum(proj, intra_image_comm)
!
   RETURN
END SUBROUTINE projwfc_hub_real
!
!-----------------------------------------------------------------------
SUBROUTINE projwfc_hub_twin(c, nx, eigr, betae, n, n_atomic_wfc,  &
& wfc, becwfc, swfc, proj, lgam)
!-----------------------------------------------------------------------
   !
   ! Projection on atomic wavefunctions
   ! Atomic wavefunctions are not orthogonized
   !
   USE kinds, ONLY: DP
   USE constants, ONLY: autoev
   USE io_global, ONLY: stdout
   USE mp_global, ONLY: intra_image_comm
   USE mp, ONLY: mp_sum
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: gstart
   USE ions_base, ONLY: nsp, na, nat
   USE uspp, ONLY: nhsa => nkb
   USE cp_interfaces, ONLY: nlsm1, s_wfc !added:giovanni
   USE twin_types !added:giovanni
!
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nx, n, n_atomic_wfc
   COMPLEX(DP), INTENT(IN) :: c(ngw, nx), eigr(ngw, nat), betae(ngw, nhsa)
   LOGICAL :: lgam
!
   COMPLEX(DP), INTENT(OUT):: wfc(ngw, n_atomic_wfc), &
                              swfc(ngw, n_atomic_wfc)
!       real(DP), intent(out):: becwfc(nhsa,n_atomic_wfc) !DEBUG
   type(twin_matrix) :: becwfc!(nhsa,n_atomic_wfc) !DEBUG
   REAL(DP), ALLOCATABLE :: overlap(:, :), e(:), z(:, :)
   REAL(DP), ALLOCATABLE :: temp(:)
   COMPLEX(DP), ALLOCATABLE :: temp_c(:)
   REAL(DP)                 :: somma
   COMPLEX(DP)              :: somma_c
   TYPE(twin_matrix)        :: proj
   INTEGER :: is, ia, nb, l, m, k, i
   !
   ! calculate number of atomic states
   !
   !
   IF (n_atomic_wfc .EQ. 0) RETURN
   !
   !
   ! calculate wfc = atomic states
   !
   CALL atomic_wfc_northo(eigr, n_atomic_wfc, wfc)
   !
   ! calculate bec = <beta|wfc>
   !
   CALL nlsm1(n_atomic_wfc, 1, nsp, eigr, wfc, becwfc, 1, lgam)
   !
   ! calculate swfc = S|wfc>
   !
   CALL s_wfc(n_atomic_wfc, becwfc, betae, wfc, swfc, lgam)
   !
   ! calculate proj = <c|S|wfc>
   !
   IF (lgam) THEN
      ALLOCATE (temp(ngw))
   ELSE
      ALLOCATE (temp_c(ngw))
   END IF
   !
   DO m = 1, n
      !
      IF (lgam) THEN
         DO l = 1, n_atomic_wfc
            temp(:) = DBLE(CONJG(c(:, m))*swfc(:, l)) !@@@@
            proj%rvec(m, l) = 2.d0*DBLE(SUM(temp))
            IF (gstart == 2) proj%rvec(m, l) = proj%rvec(m, l) - temp(1)
         END DO
      ELSE
         !
         DO l = 1, n_atomic_wfc
            temp_c(:) = CONJG(c(:, m))*swfc(:, l) !@@@@
            proj%cvec(m, l) = SUM(temp_c)
         END DO
         !
      END IF
   END DO
   !
   IF (lgam) THEN
      DEALLOCATE (temp)
      CALL mp_sum(proj%rvec, intra_image_comm)
   ELSE
      DEALLOCATE (temp_c)
      CALL mp_sum(proj%cvec, intra_image_comm)
   END IF
!
   RETURN
END SUBROUTINE projwfc_hub_twin
!
!-----------------------------------------------------------------------
SUBROUTINE atomic_wfc_northo(eigr, n_atomic_wfc, wfc)
!-----------------------------------------------------------------------
!
! Compute atomic wavefunctions in G-space
! Atomic wavefunctions not orthogonalized
!
   USE kinds, ONLY: DP
   USE gvecw, ONLY: ngw
   USE reciprocal_vectors, ONLY: gstart, g, gx
   USE ions_base, ONLY: nsp, na, nat
   USE cell_base, ONLY: tpiba, omega !@@@@
   USE atom, ONLY: rgrid
   USE uspp_param, ONLY: upf
!@@@@@
   USE constants, ONLY: fpi
!@@@@@
!
   IMPLICIT NONE
   INTEGER, INTENT(in) :: n_atomic_wfc
   COMPLEX(DP), INTENT(in) :: eigr(ngw, nat)
   COMPLEX(DP), INTENT(out):: wfc(ngw, n_atomic_wfc)
!
   INTEGER :: natwfc, ndm, is, ia, ir, nb, l, m, lm, i, lmax_wfc, isa
   REAL(DP), ALLOCATABLE ::  ylm(:, :), q(:), jl(:), vchi(:), chiq(:)

   IF (.NOT. ALLOCATED(rgrid)) &
      CALL errore(' atomic_wfc_northo ', ' rgrid not allocated ', 1)
!
! calculate max angular momentum required in wavefunctions
!
   lmax_wfc = -1
   DO is = 1, nsp
      lmax_wfc = MAX(lmax_wfc, MAXVAL(upf(is)%lchi(1:upf(is)%nwfc)))
   END DO
   !
   ALLOCATE (ylm(ngw, (lmax_wfc + 1)**2))
   !
   CALL ylmr2((lmax_wfc + 1)**2, ngw, gx, g, ylm)
   ndm = MAXVAL(rgrid(1:nsp)%mesh)
   !
   ALLOCATE (jl(ndm), vchi(ndm))
   ALLOCATE (q(ngw), chiq(ngw))
!
   DO i = 1, ngw
      q(i) = SQRT(g(i))*tpiba
   END DO
!
   natwfc = 0
   isa = 0
   DO is = 1, nsp
      !
      !   radial fourier transform of the chi functions
      !   NOTA BENE: chi is r times the radial part of the atomic wavefunction
      !
      DO ia = 1 + isa, na(is) + isa
         DO nb = 1, upf(is)%nwfc
            l = upf(is)%lchi(nb)
            DO i = 1, ngw
               CALL sph_bes(rgrid(is)%mesh, rgrid(is)%r, q(i), l, jl)
               DO ir = 1, rgrid(is)%mesh
                  vchi(ir) = upf(is)%chi(ir, nb)*rgrid(is)%r(ir)*jl(ir)
               END DO
               CALL simpson_cp90(rgrid(is)%mesh, vchi, rgrid(is)%rab, chiq(i))
            END DO
            !
            !   multiply by angular part and structure factor
            !   NOTA BENE: the factor i^l MUST be present!!!
            !
            DO m = 1, 2*l + 1
               lm = l**2 + m
               !DO ia = 1 + isa, na(is) + isa
               natwfc = natwfc + 1
               wfc(:, natwfc) = CMPLX(0.d0, 1.d0)**l*eigr(:, ia)*ylm(:, lm)*chiq(:)
               !ENDDO
            END DO
         END DO
      END DO
      isa = isa + na(is)
   END DO
!
   IF (natwfc .NE. n_atomic_wfc)                                       &
  &     CALL errore('atomic_wfc', 'unexpected error', natwfc)
!
!@@@@@
   do i = 1, n_atomic_wfc
      call DSCAL(2*ngw, fpi/sqrt(omega), wfc(1, i), 1)
   end do
!@@@@@
   DEALLOCATE (q, chiq, vchi, jl, ylm)
!
   RETURN
END SUBROUTINE atomic_wfc_northo

!-----------------------------------------------------------------------
SUBROUTINE compute_lambda(c0, gi, lambda, nspin, nbnd, ngw, nudx, desc_emp, nupdwn, iupdwn)
   !-----------------------------------------------------------------------
   !
   ! Compute matrix of lagangian multipliers (i.e. the Hamiltonian on the
   ! variational orbitals)
   !
   USE kinds, ONLY: DP
   USE twin_types
   !USE electrons_module,         ONLY : nupdwn_emp, iupdwn_emp
   USE reciprocal_vectors, ONLY: ng0 => gstart
   USE descriptors, ONLY: descla_siz_
   USE mp_global, ONLY: intra_image_comm
   USE mp, only: mp_sum
   USE cp_main_variables, ONLY: distribute_lambda
   !
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nspin, ngw, nbnd, nudx
   INTEGER, INTENT(IN) :: nupdwn(nspin), iupdwn(nspin)
   INTEGER, INTENT(IN) :: desc_emp(descla_siz_, 2)
   COMPLEX(DP), INTENT(IN) :: c0(ngw, nbnd)
   COMPLEX(DP), INTENT(IN) :: gi(ngw, nbnd)
   TYPE(twin_matrix), INTENT(INOUT)   :: lambda(nspin)
   INTEGER :: nss, is, i, j, ii, jj, istart, ig
   REAL(DP), ALLOCATABLE :: lambda_repl(:, :) ! replicated copy of lambda
   COMPLEX(DP), ALLOCATABLE :: lambda_repl_c(:, :) ! replicated copy of lambda
   !
   if (.not. lambda(1)%iscmplx) then
      allocate (lambda_repl(nudx, nudx))
   else
      allocate (lambda_repl_c(nudx, nudx))
   end if
   !
   do is = 1, nspin
      !
      nss = nupdwn(is)
      istart = iupdwn(is)
      !
      if (.not. lambda(1)%iscmplx) then
         lambda_repl = 0.d0
      else
         lambda_repl_c = CMPLX(0.d0, 0.d0)
      end if
      !
      do i = 1, nss
         !
         do j = i, nss
            !
            ii = i + istart - 1
            jj = j + istart - 1
            !
            if (.not. lambda(1)%iscmplx) then
               !
               do ig = 1, ngw
                  !
                  lambda_repl(i, j) = lambda_repl(i, j) - &
                                      2.d0*DBLE(CONJG(c0(ig, ii))*gi(ig, jj))
                  !
               end do
               !
               if (ng0 == 2) then
                  !
                  lambda_repl(i, j) = lambda_repl(i, j) + &
                                      DBLE(CONJG(c0(1, ii))*gi(1, jj))
                  !
               end if
               !
               lambda_repl(j, i) = lambda_repl(i, j)
               !
            else
               !
               do ig = 1, ngw
                  !
                  lambda_repl_c(i, j) = lambda_repl_c(i, j) - &
                                        CONJG(c0(ig, ii))*gi(ig, jj)
                  !
               end do
               !
               lambda_repl_c(j, i) = CONJG(lambda_repl_c(i, j))
               !
            end if
            !
         end do
         !
      end do
      !
      if (.not. lambda(1)%iscmplx) then
         !
         call mp_sum(lambda_repl, intra_image_comm)
         call distribute_lambda(lambda_repl, lambda(is)%rvec(:, :), desc_emp(:, is))
         !
      else
         !
         call mp_sum(lambda_repl_c, intra_image_comm)
         call distribute_lambda(lambda_repl_c, lambda(is)%cvec(:, :), desc_emp(:, is))
         !
      end if
      !
   end do
   !
   if (.not. lambda(1)%iscmplx) then
      deallocate (lambda_repl)
   else
      deallocate (lambda_repl_c)
   end if
   !
   RETURN
   !
END SUBROUTINE

