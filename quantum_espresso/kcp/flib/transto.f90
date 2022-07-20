!
! Copyright (C) 2001 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "opt_param.h"

!     OPTIMIZED DRIVER FOR MATRIX TRASPOSITION
!   
!     written by Carlo Cavazzoni 
!

      SUBROUTINE mytranspose(x, ldx, y, ldy, n, m)
!
!     x  input matrix (n by m) to be trasposed
!     y  output matrix (m by n), the transpose of x
!

        IMPLICIT NONE


        INTEGER :: ldx, ldy, n, m, what
        REAL(8) :: x(ldx, m), y(ldy, n)
        INTEGER :: i, j, k, d, nb, mb, ib, jb, ioff, joff
        INTEGER :: iind, jind
        INTEGER,  PARAMETER :: bsiz = __BSIZ_VALUE
        REAL(8) :: buf(bsiz, bsiz), bswp

        if( n>ldx ) then
          write(6,fmt='("trasponi: inconsistent ldx and n: ",2I6)') ldx, n
        end if
        if( m>ldy ) then
          write(6,fmt='("trasponi: inconsistent ldy and m: ",2I6)') ldy, m
        end if

        nb = n / bsiz 
        mb = m / bsiz 

        IF( nb < 2 .AND. mb < 2 ) THEN
          what = 1
        ELSE
          what = 2
        END IF

        select case (what)

          case (1)

            do i=1,n
              do j=1,m
                y(j,i) = x(i,j)
              enddo
            enddo

          case (2)

            do ib = 1, nb
              ioff = (ib-1) * bsiz
              do jb = 1, mb
                joff = (jb-1) * bsiz
                do j = 1, bsiz
                  do i = 1, bsiz
                    buf(i,j) = x(i+ioff, j+joff)
                  enddo
                enddo
                do j = 1, bsiz
                  do i = 1, j-1
                    bswp = buf(i,j)
                    buf(i,j) = buf(j,i) 
                    buf(j,i) = bswp
                  enddo
                enddo
                do i=1,bsiz
                  do j=1,bsiz
                    y(j+joff, i+ioff) = buf(j,i)
                  enddo
                enddo
              enddo
            enddo

            IF( MIN(1, MOD(n, bsiz)) > 0 ) THEN
              ioff = nb * bsiz
              do jb = 1, mb
                joff = (jb-1) * bsiz
                do j = 1, bsiz
                  do i = 1, MIN(bsiz, n-ioff)
                    buf(i,j) =  x(i+ioff, j+joff)
                  enddo
                enddo
                do i = 1, MIN(bsiz, n-ioff)
                  do j = 1, bsiz
                    y(j+joff,i+ioff) = buf(i,j)
                  enddo
                enddo
              enddo
            END IF

            IF( MIN(1, MOD(m, bsiz)) > 0 ) THEN
              joff = mb * bsiz
              do ib = 1, nb
                ioff = (ib-1) * bsiz
                do j = 1, MIN(bsiz, m-joff)
                  do i = 1, bsiz
                    buf(i,j) =  x(i+ioff, j+joff)
                  enddo
                enddo
                do i = 1, bsiz
                  do j = 1, MIN(bsiz, m-joff)
                    y(j+joff,i+ioff) = buf(i,j)
                  enddo
                enddo
              enddo
            END IF

            IF( MIN(1,MOD(n,bsiz))>0 .AND. MIN(1,MOD(m,bsiz))>0 ) THEN
              joff = mb * bsiz
              ioff = nb * bsiz
              do j = 1, MIN(bsiz, m-joff)
                do i = 1, MIN(bsiz, n-ioff)
                  buf(i,j) =  x(i+ioff, j+joff)
                enddo
              enddo
              do i = 1, MIN(bsiz, n-ioff)
                do j = 1, MIN(bsiz, m-joff)
                  y(j+joff,i+ioff) = buf(i,j)
                enddo
              enddo
            END IF

#if defined __ESSL
          case (3)

            CALL DGETMO (x, ldx, n, m, y, ldy)
#endif

          case default

            write(6,fmt='("trasponi: undefined method")')

        end select

        RETURN
      END SUBROUTINE  mytranspose


      SUBROUTINE mytransposez(x, ldx, y, ldy, n, m)
!
!     x  input matrix (n by m) to be trasposed
!     y  output matrix (m by n), the transpose of x
!

        IMPLICIT NONE


        INTEGER :: ldx, ldy, n, m, what
        COMPLEX(8) :: x(ldx, m), y(ldy, n)
        INTEGER :: i, j, k, d, nb, mb, ib, jb, ioff, joff
        INTEGER :: iind, jind
        INTEGER,  PARAMETER :: bsiz = __BSIZ_VALUE / 2
        COMPLEX(8) :: buf(bsiz, bsiz), bswp
!         write(6,*) "inside mytransposez" !added:giovanni:debug
        if( n>ldx ) then
          write(6,fmt='("trasponi: inconsistent ldx and n")')
        end if
        if( m>ldy ) then
          write(6,fmt='("trasponi: inconsistent ldy and m")')
        end if

        nb = n / bsiz 
        mb = m / bsiz 

        IF( nb < 2 .AND. mb < 2 ) THEN
          what = 1
        ELSE
          what = 2
        END IF
        write(6,*) "middle mytransposez", what !added:giovanni:debug
        select case (what)

          case (1)

            do i=1,n
              do j=1,m
                y(j,i) = x(i,j)
              enddo
            enddo

          case (2)

            do ib = 1, nb
              ioff = (ib-1) * bsiz
              do jb = 1, mb
                joff = (jb-1) * bsiz
                do j = 1, bsiz
                  do i = 1, bsiz
                    buf(i,j) = x(i+ioff, j+joff)
                  enddo
                enddo
                do j = 1, bsiz
                  do i = 1, j-1
                    bswp = buf(i,j)
                    buf(i,j) = buf(j,i) 
                    buf(j,i) = bswp
                  enddo
                enddo
                do i=1,bsiz
                  do j=1,bsiz
                    y(j+joff, i+ioff) = buf(j,i)
                  enddo
                enddo
              enddo
            enddo

            IF( MIN(1, MOD(n, bsiz)) > 0 ) THEN
              ioff = nb * bsiz
              do jb = 1, mb
                joff = (jb-1) * bsiz
                do j = 1, bsiz
                  do i = 1, MIN(bsiz, n-ioff)
                    buf(i,j) =  x(i+ioff, j+joff)
                  enddo
                enddo
                do i = 1, MIN(bsiz, n-ioff)
                  do j = 1, bsiz
                    y(j+joff,i+ioff) = buf(i,j)
                  enddo
                enddo
              enddo
            END IF

            IF( MIN(1, MOD(m, bsiz)) > 0 ) THEN
              joff = mb * bsiz
              do ib = 1, nb
                ioff = (ib-1) * bsiz
                do j = 1, MIN(bsiz, m-joff)
                  do i = 1, bsiz
                    buf(i,j) =  x(i+ioff, j+joff)
                  enddo
                enddo
                do i = 1, bsiz
                  do j = 1, MIN(bsiz, m-joff)
                    y(j+joff,i+ioff) = buf(i,j)
                  enddo
                enddo
              enddo
            END IF

            IF( MIN(1,MOD(n,bsiz))>0 .AND. MIN(1,MOD(m,bsiz))>0 ) THEN
              joff = mb * bsiz
              ioff = nb * bsiz
              do j = 1, MIN(bsiz, m-joff)
                do i = 1, MIN(bsiz, n-ioff)
                  buf(i,j) =  x(i+ioff, j+joff)
                enddo
              enddo
              do i = 1, MIN(bsiz, n-ioff)
                do j = 1, MIN(bsiz, m-joff)
                  y(j+joff,i+ioff) = buf(i,j)
                enddo
              enddo
            END IF

#if defined __ESSL
          case (3)

            CALL ZGETMO (x, ldx, n, m, y, ldy)
#endif

          case default

            write(6,fmt='("trasponi: undefined method")')

        end select
        RETURN
      END SUBROUTINE mytransposez

      SUBROUTINE mytransposezc(x, ldx, y, ldy, n, m)
!
!     x  input matrix (n by m) to be trasposed
!     y  output matrix (m by n), the transpose of x
!

        IMPLICIT NONE


        INTEGER :: ldx, ldy, n, m, what
        COMPLEX(8) :: x(ldx, m), y(ldy, n)
        INTEGER :: i, j, k, d, nb, mb, ib, jb, ioff, joff
        INTEGER :: iind, jind
        INTEGER,  PARAMETER :: bsiz = __BSIZ_VALUE / 2
        COMPLEX(8) :: buf(bsiz, bsiz), bswp
!         write(6,*) "inside mytransposez" !added:giovanni:debug
        if( n>ldx ) then
          write(6,fmt='("trasponi: inconsistent ldx and n")')
        end if
        if( m>ldy ) then
          write(6,fmt='("trasponi: inconsistent ldy and m")')
        end if

        nb = n / bsiz 
        mb = m / bsiz 

        IF( nb < 2 .AND. mb < 2 ) THEN
          what = 1
        ELSE
          what = 2
        END IF
!         write(6,*) "middle mytransposez", what !added:giovanni:debug
        select case (what)

          case (1)

            do i=1,n
              do j=1,m
                y(j,i) = CONJG(x(i,j))
              enddo
            enddo

          case (2)

            do ib = 1, nb
              ioff = (ib-1) * bsiz
              do jb = 1, mb
                joff = (jb-1) * bsiz
                do j = 1, bsiz
                  do i = 1, bsiz
                    buf(i,j) = x(i+ioff, j+joff)
                  enddo
                enddo
                do j = 1, bsiz
                  do i = 1, j-1
                    bswp = buf(i,j)
                    buf(i,j) = buf(j,i) 
                    buf(j,i) = bswp
                  enddo
                enddo
                do i=1,bsiz
                  do j=1,bsiz
                    y(j+joff, i+ioff) = CONJG(buf(j,i))
                  enddo
                enddo
              enddo
            enddo

            IF( MIN(1, MOD(n, bsiz)) > 0 ) THEN
              ioff = nb * bsiz
              do jb = 1, mb
                joff = (jb-1) * bsiz
                do j = 1, bsiz
                  do i = 1, MIN(bsiz, n-ioff)
                    buf(i,j) =  x(i+ioff, j+joff)
                  enddo
                enddo
                do i = 1, MIN(bsiz, n-ioff)
                  do j = 1, bsiz
                    y(j+joff,i+ioff) = CONJG(buf(i,j))
                  enddo
                enddo
              enddo
            END IF

            IF( MIN(1, MOD(m, bsiz)) > 0 ) THEN
              joff = mb * bsiz
              do ib = 1, nb
                ioff = (ib-1) * bsiz
                do j = 1, MIN(bsiz, m-joff)
                  do i = 1, bsiz
                    buf(i,j) =  x(i+ioff, j+joff)
                  enddo
                enddo
                do i = 1, bsiz
                  do j = 1, MIN(bsiz, m-joff)
                    y(j+joff,i+ioff) = CONJG(buf(i,j))
                  enddo
                enddo
              enddo
            END IF

            IF( MIN(1,MOD(n,bsiz))>0 .AND. MIN(1,MOD(m,bsiz))>0 ) THEN
              joff = mb * bsiz
              ioff = nb * bsiz
              do j = 1, MIN(bsiz, m-joff)
                do i = 1, MIN(bsiz, n-ioff)
                  buf(i,j) =  x(i+ioff, j+joff)
                enddo
              enddo
              do i = 1, MIN(bsiz, n-ioff)
                do j = 1, MIN(bsiz, m-joff)
                  y(j+joff,i+ioff) = CONJG(buf(i,j))
                enddo
              enddo
            END IF

#if defined __ESSL
          case (3)

            CALL ZGETMO (CONJG(x), ldx, n, m, y, ldy)
#endif

          case default

            write(6,fmt='("trasponi: undefined method")')

        end select
!         write(6,*) "exiting mytransposez" !added:giovanni:debug
        RETURN
      END SUBROUTINE mytransposezc

