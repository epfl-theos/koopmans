!
! Copyright (C) 2009 D. Forrer and M. Pavone
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
!
MODULE london_module
  !
  ! Module for Dispersion Correction
  ! [ V. Barone et al. J. Comp. Chem., 30, 934 (2009) ]
  ! [ S. Grimme, J. Comp. Chem., 27, 1787 (2006) ].
  !
  USE kinds ,           ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  !
  REAL ( DP ) , ALLOCATABLE , PRIVATE :: C6_i  ( : ) ,     &
                                         R_vdw ( : ) ,     &
                                         C6_ij ( : , : ) , &
                                         R_sum ( : , : ) , &
                                         r     ( : , : ) , &
                                         dist2 ( : )
  !
  ! C6_i  ( ntyp )        : atomic C6 coefficient of each atom type
  ! R_vdw ( ntyp )        : Van der Waals Radii of each atom type
  ! C6_ij ( ntyp , ntyp ) : C6 coefficients of each atom type pair: sqrt ( C6i * C6j )
  ! R_sum ( ntyp , ntyp ) : sum of VdW radii
  ! r     ( 3 , mxr )     : ordered distance vectors
  ! dist2 ( mxr )         : ordered distances
  !
  REAL ( DP ) , PUBLIC :: scal6 , lon_rcut
  !
  ! scal6    : global scaling factor
  ! lon_rcut : public cut-off radius
  !
  INTEGER , PRIVATE :: mxr
  !
  ! max number of r ( see rgen)
  !
  REAL ( DP ) , PRIVATE :: r_cut , beta = 20.0_DP
  !
  ! beta  : damping function parameter 
  ! r_cut : cut-off radius in alat units
  !
  CONTAINS
   !
   !---------------------------------------------------------------------------
   ! Initialize parameters
   !---------------------------------------------------------------------------
   !
   SUBROUTINE init_london ( )
      !
      ! extract parameters from database and compute C6_ij and R_sum(i,j)
      !
      USE ions_base ,          ONLY : ntyp => nsp, &
                                      atom_label => atm
      !
      USE cell_base ,          ONLY : alat, omega
      !
      USE io_global,           ONLY : ionode, ionode_id, stdout
      !
#if defined __PARA
      USE mp,                  ONLY : mp_bcast
#endif
      !
      IMPLICIT NONE
      !
      INTEGER, PARAMETER :: maxZ = 54
      REAL (DP) :: vdw_coeffs(2,maxZ)
      ! vdw C6 and radii for the first 54 atoms
      DATA vdw_coeffs / &
         4.857,    1.892,&
         2.775,    1.912,&
        55.853,    1.559,&
        55.853,    2.661,&
       108.584,    2.806,&
        60.710,    2.744,&
        42.670,    2.640,&
        24.284,    2.536,&
        26.018,    2.432,&
        21.855,    2.349,&
       198.087,    2.162,&
       198.087,    2.578,&
       374.319,    3.097,&
       320.200,    3.243,&
       271.980,    3.222,&
       193.230,    3.180,&
       175.885,    3.097,&
       159.927,    3.014,&
       374.666,    2.806,&
       374.666,    2.785,&
       374.666,    2.952,&
       374.666,    2.952,&
       374.666,    2.952,&
       374.666,    2.952,&
       374.666,    2.952,&
       374.666,    2.952,&
       374.666,    2.952,&
       374.666,    2.952,&
       374.666,    2.952,&
       374.666,    2.952,&
       589.405,    3.118,&
       593.221,    3.264,&
       567.896,    3.326,&
       438.498,    3.347,&
       432.600,    3.305,&
       416.642,    3.264,&
       855.833,    3.076,&
       855.833,    3.035,&
       855.833,    3.097,&
       855.833,    3.097,&
       855.833,    3.097,&
       855.833,    3.097,&
       855.833,    3.097,&
       855.833,    3.097,&
       855.833,    3.097,&
       855.833,    3.097,&
       855.833,    3.097,&
       855.833,    3.097,&
      1294.678,    3.160,&
      1342.899,    3.409,&
      1333.532,    3.555,&
      1101.101,    3.575,&
      1092.775,    3.575,&
      1040.391,    3.555/
      !
      INTEGER :: ilab , ata , atb , i
      ! local : counter of atom type
      !         ata , atb : counters of C6_ij matrix
      !         counter
      INTEGER, EXTERNAL :: atomic_number
      !!
      REAL ( DP ) :: R_0, C_0, e_cut , sls
      ! local : buffers
      !
      ! here we allocate parameters
      !
      ALLOCATE ( C6_ij ( ntyp , ntyp ) , &
                 R_sum ( ntyp , ntyp )   )
      !
      IF ( ionode ) THEN
         !
         ! and some buffers on ionode
         !
         ALLOCATE ( C6_i  ( ntyp ) , &
                    R_vdw ( ntyp )   )
         !
         ! here we initialize parameters to unphysical values
         !
         C6_i  ( : )     = -1.d0
         R_vdw ( : )     = -1.d0
         C6_ij ( : , : ) = -1.d0
         R_sum ( : , : ) = -1.d0
         !
         DO ilab = 1 , ntyp
           !
           i = atomic_number ( atom_label ( ilab ) )
           IF ( i > 0 .AND. i < 55 ) THEN
              C6_i  ( ilab )  = vdw_coeffs(1,i)
              R_vdw ( ilab )  = vdw_coeffs(2,i)
           ELSE
             CALL errore ( ' init_london ' ,&
                           'atom ' // atom_label(ilab) //' not found ' , ilab )
           END IF
           !
         END DO
         !
         ! are there all the parameters we need?
         !
         DO ilab = 1 , ntyp
           !
           IF ( ( C6_i  ( ilab ) < 0.d0 ) .or. &
                ( R_vdw ( ilab ) < 0.d0 ) ) THEN
             !
             CALL errore ( ' init_london ' ,&
                           ' one or more parameters not found ' , 4 )
             !
           END IF
           !
         END DO
         !
         ! ...here we store C6_ij parameters of each pair of atom types
         !    into a square matrix C6_ij = sqrt ( C6_i * C6_j )
         !
         DO atb = 1 , ntyp
           !
           DO ata = 1 , ntyp
             !
             C6_ij ( ata , atb ) = sqrt ( C6_i ( ata ) * C6_i ( atb ) )
             !
             R_sum ( ata , atb ) = R_vdw ( ata ) + R_vdw ( atb )
             !
           END DO
           !
         END DO
         !
         WRITE ( stdout ,'( /, 5X, "-------------------------------------" , &
                          & /, 5X, "Parameters for Dispersion Correction:" , &
                          & /, 5X, "-------------------------------------" , &
                          & /, 5X, "  atom      VdW radius       C_6     " , / )' )
         DO ata = 1 , ntyp
            !
            WRITE (stdout , '( 8X, A3 , 6X , F7.3 , 8X , F7.3 )' ) &
                               atom_label ( ata ) , R_vdw ( ata ) , C6_i ( ata )
            !
         END DO
         !
         ! ... atomic parameters are deallocated
         !
         DEALLOCATE ( C6_i , R_vdw )
         !
         ! ... cutoff radius in alat units
         !
         r_cut = lon_rcut / alat
         !
         ! ... define a gross maximum bound of the mxr size
         !
         mxr = INT ( ( 2 * ( lon_rcut + alat ) ) ** 3 / omega )
         !
      END IF
      !
#if defined __PARA
      ! broadcast data to all processors
      !
      CALL mp_bcast ( C6_ij,  ionode_id )
      CALL mp_bcast ( R_sum,  ionode_id )
      CALL mp_bcast ( r_cut,  ionode_id )
      CALL mp_bcast ( mxr  ,  ionode_id )
      !
#endif
      !
      ALLOCATE ( r ( 3 , mxr ) , dist2 ( mxr ) )
      !
      RETURN
      !
   END SUBROUTINE init_london
   !
   !---------------------------------------------------------------------------
   ! Compute dispersion energy
   !---------------------------------------------------------------------------
   !
   FUNCTION energy_london ( alat , nat , ityp , at , bg , tau )
    !
    ! here we compute the dispersion contribution to the total energy
    !
    ! E = - ( C_6^ij / R_ij ** 6 ) * f_damp ( R_ij ) * scal6
    !
    ! where f_damp is the damping function:
    !
    ! f_damp ( R_ij ) = [ 1 + exp ( -beta ( R_ij / (R_i^0+R_j^0) - 1 )) ] ** (-1)
    !
    ! and scal6 is a global scaling factor
    !
#if defined __PARA
    USE mp_global,    ONLY : mpime , nproc , intra_pool_comm
    USE mp,           ONLY : mp_sum
#endif
    !
    IMPLICIT NONE
    !
    !INTEGER , PARAMETER :: mxr = 500000
    ! local:    max number of r ( see rgen )
    !
    INTEGER :: ata , atb , nrm , nr
    ! locals : 
    ! ata , atb : atom counters
    ! nrm :       actual number of vectors computed by rgen
    ! nr :        counter
    !
    INTEGER :: nprocs , first , last , resto , divid
    ! locals :    parallelization stuff
    !
    INTEGER , INTENT ( IN ) :: nat , ityp ( nat ) 
    ! input:
    ! nat :   number of atoms
    ! itype : type of each atom
    !
    REAL ( DP ) :: dist , f_damp , energy_london , dtau ( 3 ) , dist6
    ! locals:
    ! dist          : distance R_ij between the current pair of atoms
    ! f_damp        : damping function
    ! energy_london : the dispersion energy
    ! dtau          : output of rgen ( not used )
    ! dist6         : distance**6
    !
    REAL ( DP ) , INTENT ( IN ) :: alat , tau (3, nat) , &
                                   at ( 3 , 3 ) , bg ( 3 , 3 )
    ! input :
    ! alat : the cell parameter
    ! tau  : atomic positions in alat units
    ! at   : direct lattice vectors
    ! bg   : reciprocal lattice vectors
    !
    energy_london = 0.d0
    !
#if defined __PARA
      !
      ! parallelization
      !
      nprocs = nproc
      !
      resto = MOD ( nat , nproc ) 
      divid = nat / nproc
      !
      IF ( mpime + 1 <= resto ) THEN
         !
         first = ( divid  + 1 ) * mpime + 1
         last  = ( divid  + 1 ) * ( mpime + 1 )
         !
      ELSE
         !
         first = ( ( divid + 1 ) * resto ) + ( divid ) * ( mpime - resto ) + 1
         last  = ( divid  + 1 ) * resto + ( divid ) * ( mpime - resto + 1 )
         !
      END IF
      !
#else
      !
      first = 1
      last  = nat
#endif
      !
      ! ... the dispersion energy
      !
      DO ata = first , last
        !
        DO atb = 1 , nat
          !
          dtau ( : ) = tau ( : , ata ) - tau ( : , atb )
          !
          CALL rgen ( dtau, r_cut, mxr, at, bg, r, dist2, nrm )
          !
          DO nr = 1 , nrm
            !
            dist  = alat * sqrt ( dist2 ( nr ) )
            dist6 = dist ** 6
            !
            f_damp = 1.d0 / ( 1.d0 + &
            !
            exp ( -beta * ( dist / ( R_sum ( ityp ( atb ) , ityp ( ata ) ) ) - 1 )))
            !
            energy_london = energy_london - &
                  ( C6_ij ( ityp ( atb ) , ityp ( ata ) ) / dist6 ) * &
                  f_damp
            !
          END DO
          !
        END DO
        !
      END DO
      !
      energy_london = scal6 * 0.5d0 * energy_london
      !
      !
#if defined (__PARA)
999 CALL mp_sum ( energy_london , intra_pool_comm )
#endif
    !
    RETURN
    !
   END FUNCTION energy_london
   !
   !---------------------------------------------------------------------------
   ! Compute dispersion forces acting on atoms
   !---------------------------------------------------------------------------
   !
   FUNCTION force_london ( alat , nat , ityp , at , bg , tau )
    !
    !
#if defined __PARA
    USE mp_global,    ONLY : mpime , nproc , intra_pool_comm
    USE mp,           ONLY : mp_sum
#endif
    !
    IMPLICIT NONE
    !
    !INTEGER , PARAMETER :: mxr = 500000
    ! local:    max number of r ( see rgen )
    !
    INTEGER :: ata , atb , nrm , nr , ipol
    ! locals :
    ! ata , atb : atom counters
    ! nrm       : actual number of vectors computed by rgen
    ! nr        : counter on neighbours shells
    ! ipol      : counter on coords
    !
    INTEGER :: nprocs , first , last , resto, divid
    ! locals :
    ! nprocs : number of processors in use
    ! first  : lower bound on processor
    ! last   : upper
    !
    INTEGER , INTENT ( IN ) :: nat , ityp ( nat ) 
    ! input:    
    ! nat  : number of atoms
    ! ityp : type of each atom
    !
    REAL ( DP ) :: dist , f_damp , dtau ( 3 ) , force_london ( 3 , nat ) , &
                   dist6 , dist7 , exparg , expval , par , fac , add
    ! locals :
    ! dist         : distance R_ij between the current pair of atoms
    ! f_damp       :  damping function
    ! dtau         :  \vec R_ij
    ! force_london : dispersion forces
    ! dist6        :  dist**6
    ! dist7        :  dist**7
    ! ...  and some buffers
    !
    REAL ( DP ) , INTENT ( IN ) :: alat , tau (3, nat) , &
                                   at ( 3 , 3 ) , bg ( 3 , 3 )
    ! input:
    ! alat : the cell parameter
    ! tau  : atomic positions in alat units
    ! at   : direct lattice vectors
    ! bg   : reciprocal lattice vectors
    !
    !
    force_london ( : , : ) = 0.d0
    !
#if defined __PARA
      !
      ! parallelization
      !
      nprocs = nproc
      !
      resto = MOD ( nat , nproc )
      divid = nat / nproc
      !
      IF ( mpime + 1 <= resto ) THEN
         !
         first = ( divid  + 1 ) * mpime + 1
         last  = ( divid  + 1 ) * ( mpime + 1 )
         !
      ELSE
         !
         first = ( ( divid + 1 ) * resto ) + ( divid ) * ( mpime - resto ) + 1
         last  = ( divid  + 1 ) * resto + ( divid ) * ( mpime - resto + 1 )
         !
      END IF
      !
#else
      !
      first = 1
      last  = nat
#endif
      !
      ! ... the dispersion forces
      !
      DO ata = first , last
        !
        DO atb = 1 , nat
         !
         IF ( ata /= atb ) THEN
           !
           dtau ( : ) = tau ( : , ata ) - tau ( : , atb )
           !
           ! generate neighbours shells
           !
           CALL rgen ( dtau, r_cut, mxr, at, bg, r, dist2, nrm )
           !
           ! compute forces
           !
           par = beta / ( R_sum ( ityp ( atb ) , ityp ( ata ) ) )
           !
           DO nr = 1 , nrm
            !
            dist  = alat * sqrt ( dist2 ( nr ) )
            dist6 = dist ** 6
            dist7 = dist6 * dist
            !
            exparg = - beta * ( dist / ( R_sum ( ityp ( atb ) , ityp ( ata ) ) ) - 1 )
            !
            expval = exp ( exparg )
            !
            fac = C6_ij ( ityp ( atb ) , ityp ( ata ) ) / dist6
            !
            add = 6.d0 / dist
            !
            DO ipol = 1 , 3
              !
              force_london ( ipol , ata ) = force_london ( ipol , ata ) + &
                           ( scal6 / ( 1 + expval ) * fac * &
                           ( - par * expval / ( 1.d0 + expval ) + add ) * &
                           r ( ipol , nr ) * alat / dist )
              !
            END DO
            !
           END DO
           !
         END IF
         !
        END DO
        !
      END DO
      !
#if defined (__PARA)
999 CALL mp_sum ( force_london , intra_pool_comm )
#endif
    !
    RETURN
    !
   END FUNCTION force_london
   !
   !
   !---------------------------------------------------------------------------
   ! Compute dispersion contribution to the stress tensor
   !---------------------------------------------------------------------------
   !
   FUNCTION stres_london ( alat , nat , ityp , at , bg , tau , omega )
    !
    !
#if defined __PARA
    USE mp_global,    ONLY : mpime , nproc , intra_pool_comm
    USE mp,           ONLY : mp_sum
#endif
    !
    IMPLICIT NONE
    !
    !INTEGER , PARAMETER :: mxr = 500000
    ! local:    max number of r ( see rgen )
    !
    INTEGER :: ata , atb , nrm , nr , ipol , lpol , spol
    ! locals :
    ! ata , atb : atom counters
    ! nrm       : actual number of vectors computed by rgen
    ! nr        : counter on neighbours shells
    ! xpol      : coords counters ipol lpol spol
    !
    INTEGER ::  nprocs , first , last , resto, divid
    ! locals : parallelization
    !
    INTEGER , INTENT ( IN ) :: nat , ityp ( nat ) 
    ! input:
    ! nat  : number of atoms
    ! ityp : type of each atom
    !
    REAL ( DP ) :: dist , f_damp , dtau ( 3 ) , stres_london ( 3 , 3 ) , &
                   dist6 , dist7 , exparg , expval , par , fac , add
    ! locals:
    ! dist : distance R_ij of current pair of atoms
    ! f_damp       : damping function
    ! dtau         : \vec R_ij
    ! stres_london : dispersion contribution to stress tensor
    ! dist6        : dist**6
    ! dist7        : dist**7
    !       and some buffers
    !
    REAL ( DP ) , INTENT ( IN ) :: alat , tau (3, nat) , omega , &
                                   at ( 3 , 3 ) , bg ( 3 , 3 )
    ! input :
    ! alat  : the cell parameter
    ! tau   : atomic positions in alat units
    ! omega : unit cell volume
    ! at    : direct lattice vectors
    ! bg    : reciprocal lattice vectors
    !
    !
    !
    stres_london ( : , : ) = 0.d0
    !
    first=0
    last=0
    !
#if defined __PARA
      !
      ! parallelization
      !
      nprocs = nproc
      !
      resto = MOD ( nat , nproc )
      divid = nat / nproc
      !
      IF ( mpime + 1 <= resto ) THEN
         !
         first = ( divid  + 1 ) * mpime + 1
         last  = ( divid  + 1 ) * ( mpime + 1 )
         !
      ELSE
         !
         first = ( ( divid + 1 ) * resto ) + ( divid ) * ( mpime - resto ) + 1
         last  = ( divid  + 1 ) * resto + ( divid ) * ( mpime - resto + 1 )
         !
      END IF
      !
#else
      !
      first = 1
      last  = nat
#endif
      !
      ! ... the dispersion stress tensor
      !
      DO ata = first , last
        !
        DO atb = 1 , nat
           !
           dtau ( : ) = tau ( : , ata ) - tau ( : , atb )
           !
           ! generate neighbours shells
           !
           CALL rgen ( dtau, r_cut, mxr, at, bg, r, dist2, nrm )
           !
           ! compute stress
           !
           par = beta / ( R_sum ( ityp ( atb ) , ityp ( ata ) ) )
           !
           DO nr = 1 , nrm
            !
            dist  = alat * sqrt ( dist2 ( nr ) )
            dist6 = dist ** 6
            dist7 = dist6 * dist
            !
            exparg = - beta * ( dist / ( R_sum ( ityp ( atb ) , ityp ( ata ) ) ) - 1 )
            !
            expval = exp ( exparg )
            !
            fac = C6_ij ( ityp ( atb ) , ityp ( ata ) ) / dist6
            !
            add = 6.d0 / dist
            !
            DO ipol = 1 , 3
              !
              DO lpol = 1 , ipol
                !
                stres_london ( lpol , ipol ) = stres_london ( lpol , ipol ) + &
                           ( scal6 / ( 1 + expval ) * fac * &
                           ( - par * expval / ( 1.d0 + expval ) + add ) * &
                           r ( ipol , nr ) * alat / dist ) * &
                           r ( lpol , nr ) * alat
                !
              END DO
              !
            END DO
            !
           END DO
           !
        END DO
        !
      END DO
      !
      DO ipol = 1 , 3
         !
         DO lpol = ipol + 1 , 3
            !
            stres_london ( lpol , ipol ) = stres_london ( ipol , lpol )
            !
         END DO
         !
      END DO
      !
      stres_london ( : , : ) = - stres_london ( : , : ) / ( 2.d0 * omega )
      !
#if defined (__PARA)
999 CALL mp_sum ( stres_london , intra_pool_comm )
#endif
    !
    RETURN
    !
   END FUNCTION stres_london
   !
   !---------------------------------------------------------------------------
   ! clean memory
   !---------------------------------------------------------------------------
   !
   SUBROUTINE dealloca_london
   !
   ! 
   !
   IMPLICIT NONE
   !
   IF ( ALLOCATED ( C6_i  ) ) DEALLOCATE ( C6_i  )
   IF ( ALLOCATED ( R_vdw ) ) DEALLOCATE ( R_vdw )
   IF ( ALLOCATED ( C6_ij ) ) DEALLOCATE ( C6_ij )
   IF ( ALLOCATED ( R_sum ) ) DEALLOCATE ( R_sum )
   IF ( ALLOCATED ( r     ) ) DEALLOCATE ( r     )
   IF ( ALLOCATED ( dist2 ) ) DEALLOCATE ( dist2 )
   !
   RETURN
   !
   END SUBROUTINE dealloca_london
   !
END MODULE london_module
