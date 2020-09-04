!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine all_electron(ild,ic)
  !---------------------------------------------------------------
  !
  !  this routine is a driver to an all-electron calculation
  !  with the parameters given in input
  !
  !
  use kinds, only : DP
  use io_global, only : stdout, ionode_id, ionode
  use radial_grids, only: ndmx
  use ld1inc, only: isic, grid, zeta, rho, enne, vpot, vxt, enl, &
                     deld, encl, etot, ecxc, evxt, ehrt, epseu, ekin, &
                     vnl, vh, lsd, nspin, nlcc, vdw, nn, ll, oc, nwf, &
                     zed, zval, vxc, exc, excgga, v0, verbosity, &
                     relpert, evel, edar, eso, vsic, vsicnew, vhn1, egc, &
                     wsic, wsictot, w2sic,&
                     do_unrel, nunrel, isw, psi, latt, rytoev_fact, &
                     do_affinity, do_unocc
  implicit none

  integer, intent(in) :: ic
  integer :: ioc,i,n,is
  integer, parameter :: noc = 100
  real(DP) :: frec(0:noc), etotrec(0:noc), enlrec(0:noc)
  real(DP) :: vnew(ndmx,2), rhoc1(ndmx), ze2, ocmin, funrel, delsd(nwf), desic(nwf), ekinnl(nwf)
  real(DP) :: f(nwf), enlunrel, delsdunrel,detotunrel, enlmin, ennemin, enlmax
  real(DP),allocatable :: f1(:), f2(:), f3(:), f4(:), f5(:)
  logical :: ild    ! if true compute log der
  real(DP),parameter :: vanishing_oc=1.0e-6_dp
  !
  !    compute an initial estimate of the potential
  !
  call starting_potential(ndmx,grid%mesh,zval,zed,nwf,oc,nn,ll,grid%r,&
       enl,v0,vxt,vpot,enne,nspin)
  !
  ! allocate variables for SIC, if needed
  !
  if (isic /= 0) then
     allocate(vsic(ndmx,nwf), vsicnew(ndmx), vhn1(ndmx), egc(ndmx))
     allocate(wsic(ndmx,2,nwf), w2sic(ndmx,nwf), wsictot(ndmx,2))
     vsic=0.0_dp
     wsic=0.0_dp
     w2sic=0.0_dp
     wsictot=0.0_dp
  endif
  !
  !     solve the eigenvalue self-consistent equation
  !
  call scf(ic)
  !
  !   compute relativistic corrections to the eigenvalues
  !
  if ( relpert ) call compute_relpert(evel,edar,eso)
  !
  !  compute total energy
  !
  call elsd(zed,grid,rho,vxt,vh,vxc,exc,excgga,nwf,nspin,enl,oc,    &
            etot,ekin,encl,ehrt,ecxc,evxt)
  !
  !   add sic correction if needed
  !
  if(isic /= 0) call esic  
  !
  !   print results
  !
  call write_results
  enlmax = -1e+10_dp
  enlmin = 0.0_dp
  do n=1,nwf
    if(enl(n)>enlmax .and. oc(n)>vanishing_oc) enlmax=enl(n)
  end do
  if( do_unocc ) then
    do n=nwf/2+1,nwf
      if(enl(n)<enlmin .and. oc(n)<vanishing_oc .and. oc(n-nwf/2)<2*ll(n)+1-vanishing_oc) enlmin=enl(n)
    enddo
  else
    do n=1,nwf
      if(enl(n)<enlmin .and. oc(n)<vanishing_oc) enlmin=enl(n)
    end do
  endif
  write(stdout,1202) - enlmin*rytoev_fact
  write(stdout,1201) - enlmax*rytoev_fact
   
1201  format(/5x,'- highest eigenvalue =  ',f10.3, ' eV')
1202  format(/5x,'- lowest eigenvalue  =  ',f10.3, ' eV')
  !
  if (do_unrel) then
    !
    call calc_delsd(zed,grid,vxt,vh,vxc,nwf,nspin,delsd)
    ekinnl(1:nwf)=enl(1:nwf)-delsd(1:nwf)
    if(isic /= 0) then
       call calc_desic(desic)
       ekinnl(1:nwf)=ekinnl(1:nwf)-desic(1:nwf)
    endif
    !
    ocmin=oc(nunrel)
    ennemin=sum(oc(1:nwf))
    if( .not. do_affinity ) then 
      ocmin=ocmin-1.0_dp
      ennemin=ennemin-1.0_dp
    endif
    !
    do ioc = 0, noc
      !
      ! calculate unrelaxed occupation
      !
      funrel = dble(ioc)/dble(noc)
      oc(nunrel) = ocmin+funrel
      enne = ennemin+funrel
      !
      write(stdout,1200) nunrel, oc(nunrel), enne
1200  format(/5x,'occupation( ', i4 ,' ) = ',f6.3,'   number of electrons = ',f10.3)
      !
      call nscf
      !
      ! calculate total energy
      !
      call elsd_unrel(zed,grid,rho,vxt,vh,vxc,exc,excgga,nwf,nspin,enl,oc,    &
             etot,ekin,encl,ehrt,ecxc,evxt,ekinnl)
      !
      ! calculate the energy of orbital
      !
      call calc_delsdunrel(zed,grid,vxt,vh,vxc,nwf,nspin,delsdunrel,nunrel)
      enlunrel = ekinnl(nunrel) + delsdunrel
      !
      ! calculate self-interaction contributions
      !
      if(isic /= 0) call esic_unrel(enlunrel,nunrel,funrel,ekinnl)
      !
      ! write results
      !
      call write_results
      frec(ioc)=funrel
      etotrec(ioc)=etot * rytoev_fact
      enlrec(ioc)=enlunrel * rytoev_fact
      !
    end do
    !
    detotunrel=etotrec(0)-etotrec(noc)
    !
    open(150, file = 'energy.dat')
    write(150,2020) '#','f','Etot(f)','Etot(f)-Etot(0)','eps(f)','Pi(f)'
2020  format(a1,4x,5a25)
    do ioc=0,noc
      write(150,2000) ennemin+frec(ioc),etotrec(ioc),etotrec(ioc)-etotrec(0),enlrec(ioc), &
                      enlrec(ioc)+detotunrel                
2000  format(5x,5f25.6)
    enddo
    close(150)
  end if
  !
  !  compute logarithmic derivative
  !
  if (deld > 0.0_DP .and. ild) call lderiv
  !
  ! compute C6 coefficient if required
  !
  if (vdw) then
     call c6_tfvw ( grid%mesh, zed, grid, rho(1,1) )
     call c6_dft  ( grid%mesh, zed, grid )
  end if
  !
  if (isic /= 0) then
     deallocate(vsic, vsicnew, vhn1, egc)
     deallocate(wsic, w2sic, wsictot)
  endif
  !
  return
  !
end subroutine all_electron
