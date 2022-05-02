!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!====================================================================
   SUBROUTINE inner_loop_smear( c0, bec, rhos, psihpsi )
!====================================================================
      !
      ! minimizes the total free energy
      ! using cold smearing,
      !
      !

      ! declares modules
      USE kinds,          ONLY: dp
      USE control_flags,  ONLY: iprint, thdyn, tpre, iprsta, &
                                tfor, taurdr, &
                                tprnfor, ndr, ndw, nbeg, nomore, &
                                tsde, tortho, tnosee, tnosep, trane, &
                                tranp, tsdp, tcp, tcap, ampre, &
                                amprp, tnoseh
      USE core,           ONLY: nlcc_any
      USE energies,       ONLY: eht, epseu, exc, etot, eself, enl, &
                                ekin, atot, entropy, egrand
      USE electrons_base, ONLY: f, nspin, nel, iupdwn, nupdwn, nudx, &
                                nelt, nx => nbspx, n => nbsp, ispin 

      USE gvecp,          ONLY: ngm
      USE gvecs,          ONLY: ngs
      USE gvecb,          ONLY: ngb
      USE gvecw,          ONLY: ngw
      USE reciprocal_vectors, &
                          ONLY: gstart
      USE cvan,           ONLY: nvb, ish
      USE ions_base,      ONLY: na, nat, pmass, nax, nsp, rcmax
      USE grid_dimensions, &
                          ONLY: nnr => nnrx, nr1, nr2, nr3
      USE cell_base,      ONLY: ainv, a1, a2, a3
      USE cell_base,      ONLY: omega, alat
      USE cell_base,      ONLY: h, hold, deth, wmass, tpiba2
      USE smooth_grid_dimensions, &
                          ONLY: nnrsx, nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, &
                          ONLY: nnrb => nnrbx, nr1b, nr2b, nr3b
      USE local_pseudo,   ONLY: vps, rhops
      USE io_global,      ONLY: io_global_start, stdout, ionode, &
                                ionode_id
      USE mp_global,      ONLY: intra_image_comm, leg_ortho
      USE dener
      !USE derho
      USE cdvan
      USE io_files,       ONLY: psfile, pseudo_dir, outdir
      USE uspp,           ONLY: nhsa=> nkb, betae => vkb, &
                                rhovan => becsum, deeq
      USE uspp_param,     ONLY: nh
      USE ions_positions, ONLY: tau0
      USE mp,             ONLY: mp_sum,mp_bcast, mp_root_sum

      USE cp_interfaces,  ONLY: rhoofr, dforce
      USE cp_main_variables, ONLY: distribute_lambda, descla, nlax, collect_lambda
      USE descriptors,       ONLY: lambda_node_ , la_npc_ , la_npr_ , descla_siz_ , &
                                   descla_init , la_comm_ , ilar_ , ilac_ , nlar_ , &
                                   nlac_ , la_myr_ , la_myc_ , la_nx_ , la_n_ , la_me_ , la_nrl_


      !
      IMPLICIT NONE

!input variables
      COMPLEX(kind=DP), INTENT(IN)  :: c0( ngw, n )
      REAL(kind=DP),    INTENT(IN)  :: bec( nhsa, n )
      REAL(kind=DP),    INTENT(IN)  :: rhos( nnrsx, nspin )
      REAL(kind=DP),    INTENT(OUT) :: psihpsi( nlax, nlax, nspin )
  

!local variables
      REAL(kind=DP),    ALLOCATABLE :: c0hc0(:,:,:)
      REAL(kind=DP),    ALLOCATABLE :: mtmp(:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: h0c0(:,:)
      INTEGER :: i,k, is, nss, istart, ig, iss

      INTEGER :: np(2), coor_ip(2), ipr, ipc, nr, nc, ir, ic, ii, jj, root, j
      INTEGER :: desc_ip( descla_siz_ )
      INTEGER :: np_rot, me_rot, comm_rot, nrl


      allocate(c0hc0(nlax,nlax,nspin))
      allocate(h0c0(ngw,nx))


      ! operates the Hamiltonian on the wavefunction c0
      h0c0( :, : )= 0.D0
      DO i= 1, n, 2                      
         CALL dforce( i, bec, betae, c0, h0c0(:,i), h0c0(:,i+1), rhos, nnrsx, ispin, f, n, nspin )
      END DO

 
       
      ! calculates the Hamiltonian matrix in the basis {c0}           
      c0hc0(:,:,:)=0.d0
      !
      DO is= 1, nspin

         nss= nupdwn( is )
         istart= iupdwn( is )

         np(1) = descla( la_npr_ , is )
         np(2) = descla( la_npc_ , is )

         DO ipc = 1, np(2)
            DO ipr = 1, np(1)

               coor_ip(1) = ipr - 1
               coor_ip(2) = ipc - 1
               CALL descla_init( desc_ip, descla( la_n_ , is ), &
                                 descla( la_nx_ , is ), np, coor_ip, descla( la_comm_ , is ), 1 )

               nr = desc_ip( nlar_ )
               nc = desc_ip( nlac_ )
               ir = desc_ip( ilar_ )
               ic = desc_ip( ilac_ )

               CALL GRID2D_RANK( 'R', desc_ip( la_npr_ ), desc_ip( la_npc_ ), &
                                 desc_ip( la_myr_ ), desc_ip( la_myc_ ), root )
               !
               root = root * leg_ortho

               ALLOCATE( mtmp( nr, nc ) )
               mtmp = 0.0d0

               CALL DGEMM( 'T', 'N', nr, nc, 2*ngw, - 2.0d0, c0( 1, istart + ir - 1 ), 2*ngw, &
                           h0c0( 1, istart + ic - 1 ), 2*ngw, 0.0d0, mtmp, nr )

               IF (gstart == 2) THEN
                  DO jj = 1, nc
                     DO ii = 1, nr
                        i = ii + ir - 1
                        j = jj + ic - 1
                        mtmp(ii,jj) = mtmp(ii,jj) + DBLE( c0( 1, i + istart - 1 ) ) * DBLE( h0c0( 1, j + istart - 1 ) )
                     END DO
                  END DO
               END IF

               CALL mp_root_sum( mtmp, c0hc0(1:nr,1:nc,is), root, intra_image_comm )

               DEALLOCATE( mtmp )

            END DO
         END DO
      END DO

      
      psihpsi(:,:,:) = c0hc0(:,:,:)


      DEALLOCATE(c0hc0)
      DEALLOCATE(h0c0)

      return
!====================================================================      
    END SUBROUTINE inner_loop_smear
!====================================================================



