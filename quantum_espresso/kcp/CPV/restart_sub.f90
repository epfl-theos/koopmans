!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
SUBROUTINE from_restart( )
   !
   USE kinds,                 ONLY : DP
   USE control_flags,         ONLY : program_name, tbeg, taurdr, tfor, tsdp, tv0rd, &
                                     iprsta, tsde, tzeroe, tzerop, nbeg, tranp, amprp, &
                                     tzeroc, force_pairing, ampre, trane, tpre, dt_old, &
                                     iprint_manifold_overlap !added:giovanni iprint_manifold_overlap
   USE wavefunctions_module,  ONLY : c0, cm, cstart
   USE electrons_module,      ONLY : occn_info
   USE electrons_base,        ONLY : iupdwn, nupdwn, f, nbsp, nbspx
   USE io_global,             ONLY : ionode, stdout
   USE cell_base,             ONLY : ainv, h, hold, deth, r_to_s, s_to_r, boxdimensions, &
                                     velh, a1, a2, a3
   USE ions_base,             ONLY : na, nsp, iforce, vel_srt, nat, randpos
   USE time_step,             ONLY : tps, delt
   USE ions_positions,        ONLY : taus, tau0, tausm, taum, vels, fion, set_velocities
   USE ions_nose,             ONLY : xnhp0, xnhpm
   USE grid_dimensions,       ONLY : nr1, nr2, nr3
   USE reciprocal_vectors,    ONLY : mill_l
   USE printout_base,         ONLY : printout_pos
   USE gvecs,                 ONLY : ngs
   USE gvecw,                 ONLY : ngw
   USE cp_interfaces,         ONLY : phfacs, strucf
   USE energies,              ONLY : eself, dft_energy_type
   USE wave_base,             ONLY : rande_base
   USE efield_module,         ONLY : efield_berry_setup,  tefield, &
                                     efield_berry_setup2, tefield2
   USE small_box,             ONLY : ainvb
   USE uspp,                  ONLY : okvan, vkb, nkb
   USE core,                  ONLY : nlcc_any
   USE cp_main_variables,     ONLY : ht0, htm, lambdap, lambda, lambdam, ei1, ei2, ei3, eigr, &
                                     sfac, bec, taub, irb, eigrb, edft, becstart
   USE cdvan,                 ONLY : dbec
   USE time_step,             ONLY : delt
   USE atoms_type_module,     ONLY : atoms_type
   USE twin_types
   !
   IMPLICIT NONE
   !
   ! ... We are restarting from file recompute ainv
   !
   CALL invmat( 3, h, ainv, deth )
   !
   ! ... Reset total time counter if the run is not strictly 'restart'
   !
   IF ( nbeg < 1 ) tps = 0.D0
   !
   IF ( taurdr ) THEN
      !
      ! ... Input positions read from input file and stored in tau0
      ! ... in readfile, only scaled positions are read
      !
      CALL r_to_s( tau0, taus, na, nsp, ainv )
      !
   END IF
   !
   IF ( ANY( tranp(1:nsp) ) ) THEN
      !
      ! ... Input positions are randomized
      !
      CALL randpos( taus, na, nsp, tranp, amprp, ainv, iforce )
      !
   END IF
   !
   IF ( tzerop .AND. tfor ) THEN
      !
      CALL r_to_s( vel_srt, vels, na, nsp, ainv )
      !
      CALL set_velocities( tausm, taus, vels, iforce, nat, delt )
      !
   END IF
   !
   CALL s_to_r( taus,  tau0, na, nsp, h )
   !
   CALL s_to_r( tausm, taum, na, nsp, h )
   !
   IF ( tzeroc ) THEN
      !
      hold = h
      velh = 0.D0
      !
      htm      = ht0
      ht0%hvel = 0.D0
      !
   END IF
   !
   fion = 0.D0
   !
   IF( force_pairing ) THEN
      cm(:,iupdwn(2):nbsp) = cm(:,1:nupdwn(2))
      c0(:,iupdwn(2):nbsp) = c0(:,1:nupdwn(2))
      IF(.not.lambda(1)%iscmplx) THEN
	  lambdap(2)%rvec( :, :) =  lambdap(1)%rvec( :, :)
	  lambda(2)%rvec( :, :) =  lambda(1)%rvec( :, :)
	  lambdam(2)%rvec( :, :) =  lambdam(1)%rvec( :, :)
      ELSE
	  lambdap(2)%cvec( :, :) =  lambdap(1)%cvec( :, :)
	  lambda(2)%cvec( :, :) =  lambda(1)%cvec( :, :)
	  lambdam(2)%cvec( :, :) =  lambda(1)%cvec( :, :)
      ENDIF
   END IF 
   !
   IF ( tzeroe ) THEN
      !
      lambdam = lambda
      !
      cm = c0
      !
      WRITE( stdout, '(" Electronic velocities set to zero")' )
      !
   END IF
   !
   IF( program_name == "FPMD" ) THEN
      !
      CALL occn_info( f )
      !
      ! ... diagnostics
      !
      IF ( ionode ) THEN
         !
         WRITE( stdout, 100 )
         !
         IF ( .NOT. tbeg ) THEN
            WRITE( stdout, 110 )
         ELSE
            WRITE( stdout, 120 )
         END IF
         !
         IF ( .NOT. taurdr ) THEN
            WRITE( stdout, 130 )
         ELSE
            WRITE( stdout, 140 )
         END IF
         !
         IF ( tfor .AND. ( .NOT. tsdp ) ) THEN
            !
            IF ( .NOT. tv0rd ) THEN
               !
               IF ( .NOT. tzerop ) THEN
                  WRITE( stdout, 150 )
               ELSE
                  WRITE( stdout, 155 )
               END IF
               !
            ELSE
               !
               WRITE( stdout, 160 )
               !
            END IF
          !
         END IF
         !
         IF ( iprsta > 1 ) &
            CALL printout_pos( stdout, taus, nat, head = 'Scaled positions from restart module' )
         !
         IF ( .NOT. tsde ) THEN
            !
            IF ( .NOT. tzeroe ) THEN
               WRITE( stdout, 170 )
            ELSE
               WRITE( stdout, 180 )
            END IF
            !
         END IF
         !
         WRITE( stdout, * )
         !
      END IF
      !
   END IF
   !
   ! ... computes form factors and initializes nl-pseudop. according
   ! ... to starting cell (from ndr or again standard input)
   !
   IF ( okvan .or. nlcc_any ) THEN
      CALL initbox( tau0, taub, irb, ainv, a1, a2, a3 )
      CALL phbox( taub, eigrb, ainvb )
   END IF
   !
   CALL phfacs( ei1, ei2, ei3, eigr, mill_l, taus, nr1, nr2, nr3, nat )
   !
   CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
   !
   CALL prefor( eigr, vkb )
   !
   CALL formf( .TRUE. , eself )
   !
   IF ( trane ) THEN
      !
      WRITE( stdout, 515 ) ampre
      !
515   FORMAT(   3X,'Initial random displacement of el. coordinates',/ &
                3X,'Amplitude = ',F10.6 )
      !
      CALL rande_base( c0, ampre )

      CALL gram( vkb, bec, nkb, c0, ngw, nbsp )
      !
      IF( force_pairing ) c0(:,iupdwn(2):nbsp) = c0(:,1:nupdwn(2))
      !
      cm = c0
      !
   END IF
   !
   CALL calbec( 1, nsp, eigr, c0, bec )
   !
   IF ( tpre     ) CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, c0, dbec )
   !
   IF ( tefield  ) CALL efield_berry_setup( eigr, tau0 )
   IF ( tefield2 ) CALL efield_berry_setup2( eigr, tau0 )
   !
   edft%eself = eself
   !
   IF( tzerop .or. tzeroe .or. tzeroc ) THEN
      IF( .not. ( tzerop .and. tzeroe .and. tzeroc ) ) THEN
         IF( ionode ) THEN
            WRITE( stdout, * ) 'WARNING setting to ZERO ions, electrons and cell velocities without '
            WRITE( stdout, * ) 'setting to ZERO all velocities could generate meaningles trajectories '
         END IF
      END IF
   END IF
   !
   !
   ! dt_old should be -1.0 here if untouched ...
   !
   if ( dt_old > 0.0d0 ) then
      tausm = taus - (taus-tausm)*delt/dt_old
      xnhpm = xnhp0 - (xnhp0-xnhpm)*delt/dt_old
      WRITE( stdout, '(" tausm & xnhpm were rescaled ")' )
   endif
   !
   IF(iprint_manifold_overlap>0) THEN !added:giovanni
      !
      cstart(1:ngw,1:nbspx) = c0(1:ngw,1:nbspx)
      call copy_twin(becstart, bec) 
      !
   ENDIF
   !
   RETURN
   !
100 FORMAT( /,3X,'MD PARAMETERS READ FROM RESTART FILE',/ &
             ,3X,'------------------------------------' )
110 FORMAT(   3X,'Cell variables From RESTART file' )
120 FORMAT(   3X,'Cell variables From INPUT file' )
130 FORMAT(   3X,'Ions positions From RESTART file' )
140 FORMAT(   3X,'Ions positions From INPUT file' )
150 FORMAT(   3X,'Ions Velocities From RESTART file' )
155 FORMAT(   3X,'Ions Velocities set to ZERO' )
160 FORMAT(   3X,'Ions Velocities From STANDARD INPUT' )
170 FORMAT(   3X,'Electronic Velocities From RESTART file' )
180 FORMAT(   3X,'Electronic Velocities set to ZERO' )
   !
END SUBROUTINE from_restart
