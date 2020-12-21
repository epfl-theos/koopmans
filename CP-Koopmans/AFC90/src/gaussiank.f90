!
! LibAFCC - Library for auxiliary-function countercharge correction 
! Copyright (c) 2010-2011 I. Dabo
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. See GPL/gpl-3.0.txt. 
! If not, see <http://www.gnu.org/licenses/>.
!
function gaussiank(sigma,x)
  !
  implicit none
  !
  real(8), intent(in) :: sigma
  real(8), intent(in), dimension(:) :: x
  real(8), dimension(size(x)) :: gaussiank
  real(8), parameter :: sigmamax=4.d0
  real(8), parameter :: sigmamin=1.d-6
  real(8), parameter :: tmax=20.d0
  real(8), parameter :: tol=1.d-10
  real(8), parameter :: xthr=1.d-10
  real(8), dimension(2) :: yd,y
  real(8) :: logslope
  real(8) :: temp1,temp2,x1,t1,t0,tx,t,dt
  real(8) :: error
  real(8), allocatable, dimension(:,:) :: z
  integer :: i,n,nt
  !
  interface 
    !
    function besselk(x)
      real(8), intent(in) :: x
      real(8) :: besselk
    end function
    !
    function cylharm(dt,sigma,x)
      real(8), intent(in) :: dt,sigma,x
      real(8) :: cylharm
    end function
    !
    function cylharm0(sigma)
      real(8), intent(in) :: sigma
      real(8) :: cylharm0
    end function
    !
    function cylharmasympt(x)
      real(8), intent(in) :: x
      real(8) :: cylharmasympt
    end function
    !
    function cylharmasymptdot(x)
      real(8), intent(in) :: x
      real(8) :: cylharmasymptdot
    end function
    !
    function cylharmasymptdotdot(x)
      real(8), intent(in) :: x
      real(8) :: cylharmasymptdotdot
    end function
    !
    function cylharmseries(a0,b0,sigma,x)
      real(8), intent(in) :: a0,b0,sigma,x
      real(8) :: cylharmseries
    end function
    !
    function cylharmslope(tol,dt,sigma)
      real(8), intent(in) :: tol,dt,sigma
      real(8) :: cylharmslope
    end function
    !
    function pinterp(x,side,bound)
      real(8), intent(in) :: x
      integer :: side,bound
      real(8) :: pinterp
    end function 
    !
    function qinterp(x,side,bound)
      real(8), intent(in) :: x
      integer :: side,bound
      real(8) :: qinterp
    end function 
    !
    function rungekutta(sigma,t,dt,x)
      real(8), intent(in), dimension(2) :: x
      real(8), intent(in) :: sigma,t,dt
      real(8), dimension(2) :: rungekutta
    end function
    !
    function steprk(sigma)
      real(8), intent(in) :: sigma
      real(8) :: steprk
    end function
    !
    function ydot(sigma,t,x)
      real(8), intent(in) :: sigma,t
      real(8), intent(in), dimension(2) :: x
      real(8), dimension(2) :: ydot
    end function
    !
  end interface
  !
  if (sigma.le.sigmamin) then
    do i=1,size(x)
      gaussiank(i)=besselk(x(i))
    enddo
  elseif (sigma.ge.sigmamax) then
    gaussiank(:)=0.d0
  else
    !
    dt=steprk(sigma)
#ifdef __AFC90_DEBUG
    write(6,*) dt, "time"
    write(0,*) "entering logslope"
#endif
    logslope=cylharmslope(1.d-10,dt,sigma)
!     stop
    !
    x1=max(minval(x),xthr)
    t1=log(x1)
    t0=t1
    nt=1
    y=(/cylharmasympt(t0),cylharmasymptdot(t0)/)
    yd=ydot(sigma,t0,y)
    error=abs(yd(2)-cylharmasymptdotdot(t0))
    searcht0: do while (error>tol)
      if (t0+dt>tmax) then
        print *, 'warning: maximum t0 reached in subroutine gaussiank'
        exit searcht0
      endif
      t0=t0+dt
      nt=nt+1
      y=(/cylharmasympt(t0),cylharmasymptdot(t0)/)
      yd=ydot(sigma,t0,y)
      error=abs(yd(2)-cylharmasymptdotdot(t0))
    enddo searcht0
    !
    allocate(z(2,nt))
    n=nt
    t=t0
    z(1:2,n)=(/cylharmasympt(t),cylharmasymptdot(t)/)
    rk: do while (t>t1)
      z(1:2,n-1)=rungekutta(sigma,t,-dt,z(1:2,n))
      t=t-dt
      n=n-1
    enddo rk
    !
    do i=1,size(x)
      if (x(i)>xthr) then
        tx=log(max(x(i),xthr))
        if (tx>t0) then
          temp1=cylharmasympt(tx)+logslope*besselk(x(i))
        else
          n=int((tx-t1)/dt)+1
          t=t1+(n-1)*dt
          temp1=pinterp((tx-t)/dt,0,0)*z(1,n) &
               +pinterp((tx-t)/dt,1,0)*z(1,n+1) &
               +dt*qinterp((tx-t)/dt,1,0)*z(2,n) &
               +dt*qinterp((tx-t)/dt,1,0)*z(2,n+1) &
               +logslope*besselk(x(i))
          temp2=cylharmseries(cylharm0(sigma),0.d0,sigma,x(i))
          if (abs(temp1-temp2)>tol*abs(temp1).or.isnan(temp2)) then
            gaussiank(i)=temp1
          else
            gaussiank(i)=temp2
          endif
        endif
      else
        gaussiank(i)=cylharm0(sigma)
      endif
    enddo
    deallocate(z)
    !
  endif
  !
end function gaussiank
