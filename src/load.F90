subroutine load
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface
!ZHIXIN use w_initial_in
  use control_parameters
  implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

  integer i,m
  integer,dimension(:),allocatable :: mget,mput
  real(wp) :: c0=2.515517,c1=0.802853,c2=0.010328,d1=1.432788,d2=0.189269,d3=0.001308
  real(wp) w_initial,energy,momentum,r,rmi,pi2_inv,delr,tsqrt,cmratio,true_micell
  real(wp) my_ion,my_eon,total_ion,total_eon

  integer j,mp0,mp1,mt0,mt1,ierror,k,mtheta_i,mtheta_load
  real(wp) b,totalv,tdum,zdum,dum,rdum,delp,delz,dell,deltat_i,g,ri,pdum,q,volt,ave0
  real(wp),dimension(:,:,:),allocatable :: volume

! integer(int8) :: mp1_i8,misum_i8
  integer :: mp1_i8,misum_i8

! the following definition is for volume integral using NAG
  integer minpts,ifail
  integer,parameter :: maxpts=1000*3,lenwrk=(3+2)*(1+maxpts/33)
  real(wp) low_lim(3),up_lim(3),wrkstr(lenwrk),acc,error_8,finval
  real(wp),external :: f_jaco,FUN1,dsda,D01AHF,temperature,temperature_e,&
                       u_parallel
  real(wp),dimension(:,:),allocatable :: local_vol
  real(wp) vol_on_node,dum_nag

  real(wp), dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj

  real(wp), external :: object_u0

  real(wp)  fip, dela
  integer   ip

  if(irun==0) ddT_da=0.0_wp
  if(irun==0) sumfluxpsi=0.0_wp      ! sumfluxpsi(a)=int[I(a)dt]
! restart from previous runs
  if(irun /= 0)then
     call restart_read
     return
  endif

  rmi=1.0_wp/real(mi,wp)
  pi2_inv=0.5_wp/pi
  delr=1.0_wp/deltar
! initial weight amplitude
!  w_initial=1.0e-3_wp   ! for ITG
!  w_initial=1.0e-4_wp    ! for TEM
!  w_initial=4.0e-4_wp    ! for TEM
!  w_initial=4.0e-3_wp    ! for TEM

!!!  w_initial=1.0e-6_wp    ! for TEM
!!!  w_initial=1.0e-4_wp    ! for ITG

!  w_initial=1.0e-12_wp    ! for linear growth rate and frequency
!  w_initial=1.0e-18_wp    ! for linear growth rate and frequency
!!!  if(nonlinear<0.5_wp) w_initial=1.0e-20_wp
!ZHIXIN: use w_initial_in to set w_initial
  w_initial=w_initial_in*1.0e0_wp
  if(nonlinear<0.5_wp) w_initial=w_initial_in*1.0e-16_wp
      
! uniform load particle in configuration space 
  mtheta_load=256
!!  allocate(volume(mpsi,mtheta_load))
!!  volume=0.0

  allocate(local_vol(mtheta_load,mpsi))
  local_vol=0.0_wp

!   NOTE:  !!s are replaced by volume multiple integral using d01fcf of NAG

  dell=plenth(mpsi)/real(mtheta_load,wp)
  delz=2.0_wp*pi/real(mzetamax,wp)

  do i=1,mpsi
     mtheta_i=max(2,min(mtheta_load,2*int(plenth(i)/dell/2.0_wp+0.5_wp))) !even # poloidal grid
     deltat_i=2.0_wp*pi/real(mtheta_i,wp)
     do j=1,mtheta_i
        low_lim(1)=psi_grid(i-1)
        up_lim(1)=psi_grid(i)
        low_lim(2)=deltat_i*real(j-1,wp)
        up_lim(2)=deltat_i*real(j,wp)
        low_lim(3)=zetamin
        up_lim(3)=zetamax
        error_8=1.0e-5_wp
        minpts=0
        ifail=1
        call d01fcf(3,low_lim,up_lim,minpts,maxpts,FUN1, &
               error_8,acc,lenwrk,wrkstr,finval,ifail)
        local_vol(j,i)=finval
     enddo
  enddo
  
!!  totalv=sum(volume)
 
  vol_on_node=sum(local_vol)
 
  mp0=0
  do i=1,mpsi
     mtheta_i=max(2,min(mtheta_load,2*int(plenth(i)/dell/2.0+0.5))) !even # poloidal grid
     deltat_i=2.0*pi/real(mtheta_i)

     mt0=mp0
     do j=1,mtheta_i
!!        volt=sum(volume(:,i,j))
!!        dum=real(mi)*volt/totalv

        dum_nag=real(mi,wp)*local_vol(j,i)/vol_on_node  !wx
        dum=dum_nag
 
        mt1=mt0+floor(dum)                     ! integer number of particle
        dum=dum-floor(dum)

        rdum=sprng()

        if(rdum<dum)mt1=mt1+1                  ! fractional number of particle

! uniformly distributed inside poloidal cell

        do k=mt0+1,mt1
           zion(2,k)=sprng()
        end do

        zion(2,mt0+1:mt1)=deltat_i*real(j-1)+deltat_i*zion(2,mt0+1:mt1)

! uniformly distributed in zeta

        do k=mt0+1,mt1
           zion(3,k)=sprng()
        end do

        delz=zetamax-zetamin
        zion(3,mt0+1:mt1)=zetamin+delz*zion(3,mt0+1:mt1)
        mt0=mt1
     enddo

! uniformly distributed inside radial cell
     mp1=mt1
     delp=psi_grid(i)-psi_grid(i-1)

     do k=mp0+1,mp1
         zion(1,k)=sprng()
     end do

     zion(1,mp0+1:mp1)=psi_grid(i-1)+delp*zion(1,mp0+1:mp1)

     mp0=mp1
  enddo

  mp1_i8=mp1   !Copy value of mp1 into 8-byte integer variable  mp1_i8
! Sum up the number of particles from all the processors using an 8-byte
! integer variable since the total may be larger than 2.1 billion.
  call MPI_ALLREDUCE(mp1_i8,misum_i8,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierror)

!     markeri is the inverse of ion marker number at each grid point in term of initial loading
  !  markeri=markeri*(real(misum_i8,wp)/real(ntoroidal,wp))
  markeri=markeri*real(misum_i8,wp)
  markeri=1.0_wp/markeri   ! to avoid divide operation

!     pmarki is the inverse of total particle # (sum over nodes) at an anulus i
  pmarki=pmarki*real(misum_i8,wp)
  pmarki=1.0_wp/pmarki     ! to avoid divide operation

  rtden=rtden*misum/real(misum_i8,wp)
  mi=max(1,mp1)   ! really loaded particle number on my node
  misum=real(misum_i8,wp)  ! really loaded total particle number

  sn0=sn0*misum
  sn0_g=sn0_g*misum
  sn_box=sn_box*misum

  true_micell=misum/real(mzetamax*(mgrid-mpsi),wp)

  cloud_i=cloud_i/misum

  open(988,file='initial_loading.dat',status='unknown')
  if(mype==0)write(988,*)'Average Particle Number Per Cell=', true_micell
  write(988,*)'mype=', mype, 'I have ', mi, ' particles'
  close(988)

  if(mype==0) print*, 'Average Particle Number Per Cell=', true_micell
!  print*, 'mype=', mype, 'I have ', mi, ' particles'

!    ion tracer (ntracer=mi) resides on PE=0 initially.
  ntracer=0
  if(mype==0) ntracer=mi

!!  deallocate(volume)
  deallocate(local_vol)

  do k=1,mi
      zion(4,k)=sprng()
  end do

  do k=1,mi
      zion(5,k)=sprng()
  end do

  do k=1,mi
      zion(6,k)=sprng()
  end do

! Maxwellian distribution in V_para, use zion0(4,:) as temporary storage
!$omp parallel do private(m)
  do m=1,mi
     zion(4,m)=zion(4,m)-0.5_wp
     zion0(4,m)=sign(1.0_wp,zion(4,m))
     zion(4,m)=sqrt(max(1.0e-20_wp,log(1.0_wp/max(1.0e-20_wp,zion(4,m)**2))))
     zion(4,m)=zion(4,m)-(c0+c1*zion(4,m)+c2*zion(4,m)**2)/       &
                  (1.0_wp+d1*zion(4,m)+d2*zion(4,m)**2+d3*zion(4,m)**3)
  end do

  do m=1,mi                
     if(zion(4,m)>umax) zion(4,m)=sprng()
  enddo

!$omp parallel do private(m)
  do m=1,mi

     zion(4,m)=zion0(4,m)*min(umax,zion(4,m))

! initial random weight
     zion(5,m)=2.0_wp*w_initial*(zion(5,m)-0.5_wp)*(1.0_wp+cos(zion(2,m)))

! Maxwellian distribution in V_perp
     zion(6,m)=max(1.0e-20_wp,min(umax*umax,-log(max(1.0e-20_wp,zion(6,m))))) ! =[v_perp/sqrt(2T/m)]^2 ???
  enddo

! transform zion(4,:) to rho_para, zion(6,:) to sqrt(mu) 
  cmratio=qion/aion
!  call esigetb(zion(1,:),zion(2,:),mi)   ! get B
!  call esigetb(zion,mi)   ! get B-->BB
  call esigetbf(zion,mi,7)   ! get B-->BB, I-->F  12/07
  call n_T_u_group                       ! get T, n and omega_t
!  output backu from call n_T_u_group is omega_t
  backu(1:mi)=backu(1:mi)*F(1:mi)/BB(1:mi) ! =u_//=I*omega_t/B

!$omp parallel do private(m)
  do m=1,mi
!     zion(4,m)=cmratio*(zion(4,m)*sqrt(backT(m))+backu(m))/BB(m)
     zion(4,m)=zion(4,m)*sqrt(backT(m))/BB(m)
     if(ifelectron/=0) zeon(4,m)=zion(4,m)  ! neglect electron rotation
     zion(4,m)=zion(4,m)+backu(m)/BB(m)
     zion(6,m)=sqrt(zion(6,m)/BB(m))*sqrt(backT(m))
  end do

! tracer particle initially resides on PE=0
  if(mype == 0)then
     zion(1,ntracer)=0.5_wp*(a0+a1)
     zion(2,ntracer)=0.0_wp
     zion(3,ntracer)=0.5_wp*(zetamin+zetamax)
     pdum=temperature(zion(1,ntracer),0.0_wp)
     backu(ntracer)=u_parallel(zion(1,ntracer),zion(2,ntracer))
     zion(4,ntracer)=(sqrt(pdum)+backu(ntracer))/B_0      ! using v_par=v_th for tracer ion
     zion(5,ntracer)=0.0_wp
     zion(6,ntracer)=sqrt(pdum)/sqrt(2.0_wp*B_0)   ! using v_per=v_th for tracer ion
     zion0(1:5,ntracer)=zion(1:5,ntracer)

     print*,'ntracer=',ntracer,zion(1,ntracer),zion(2,ntracer),zion(3,ntracer),&
                       zion(4,ntracer),zion(5,ntracer),zion(6,ntracer)

     aaa(1)=zion(1,ntracer)
     the(1)=zion(2,ntracer)
     call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
     etracer_i0=0.5_wp*(zion(4,ntracer)*Bj(1))**2+zion(6,ntracer)**2*Bj(1)      ! kinetic energy
  endif

  do m=1,mi
     error_8=1.0e-5_wp
     IFAIL=1
     zion(7,m)=D01AHF(0.0_wp,zion(1,m),error_8,minpts,acc,dsda,-1,IFAIL)
     zion(7,m)=min(s1,max(s0,zion(7,m)))
  enddo

  if(iload==0)then
! for arbitrary loading according to a distribution function g_0(x,v), then zion0(6,*)=f_0/g_0
!         where f_0 is equilibrium distribution with delta_f=f-f_0.
!         Here g_0 contains density <n(r,theta)> rather than n(r,theta), so
!         f_0/g_0=exp(u^2/2T)/<exp(u^2/2T)>.
  dela=(a1-a0)/real(jr,wp)
  dela=1.0_wp/dela
  do m=1,mi
!         the following call takes time
!      call general_flux_average(zion(1,m),object_u0,ave0) ! ave0=<exp(u^2/2T)>
!         replaced by linear interpolation in terms of a pre-generated table flux_A. 3/08
     r=zion(1,m)-a0
     r=r*dela
     ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
     fip=r-real(ip,wp)
     ave0=fip*flux_A(ip+1)+(1.0_wp-fip)*flux_A(ip)       ! ave0=<exp(u^2/2T)> at particles
     zion0(6,m)=exp(backu(m)**2/(2.0*backT(m)))/ave0     ! =1 if u_//=0 (no rotation)
  end do

!!! the following print may hang forever when using large number of processors Jan, 2011
! print*, 'f_0/g_0=', zion0(6,mi/2), mype

! realistic nonuniform loading
!     zion0(6,1:mi)=1.0

  else
!    what else artificial loading, for example uniform loading -- to be added when needed
  endif

  sw=0      ! no use of x-y
     !   sw(i)=1: use x-y for aa<near_axis; if near_axis<0.0=>sw=0: no use of x-y
  
  rdtem0=0.0_wp  ! initial Ti fluctuation

  do m=1,mi
     enion(m)=0.5_wp*(zion(4,m)*BB(m))**2+zion(6,m)**2*BB(m)
  enddo
  !enion(:)=0.5_wp*(zion(4,:)*BB(:))**2+zion(6,:)**2*BB(:)


! set SUM(w)=0 <==> int(delta-n(t=0)*d^3x)=0.0

  dum=0.0_wp
!$omp parallel private(m) reduction(+:dum)
!$omp do
  do m=1,mi
     dum=dum+zion(5,m)
  enddo
!$omp end parallel
  call MPI_ALLREDUCE(dum,pdum,1,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
  pdum=pdum/misum
!$omp parallel do private(m)
  do m=1,mi
     zion(5,m)=zion(5,m)-pdum
  enddo
if(mype==0) print*, 'average weight at t=0 is ', pdum

  
! load electron on top of ion, including both trapped and untrapped electrons
   me=mi
!   call MPI_ALLREDUCE(me,mesum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)
   mesum=misum

  my_ion=0.0_wp
!$omp parallel do private(m) reduction(+:my_ion)
  do m=1,mi
     my_ion=my_ion+zion(5,m)
  end do

  call MPI_Reduce(my_ion,total_ion,1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)

  if(mype==0) print*, 'INITIAL ION WEIGHT=', total_ion/misum



  if(ifelectron/=0) then
!$omp parallel do private(m)
     do m=1,me
        zeon(1,m)=zion(1,m)
        zeon(2,m)=zion(2,m)
        zeon(3,m)=zion(3,m)
        zeon(5,m)=zion(5,m)
        zeon(7,m)=zion(7,m)
!        zeon0(6,m)=zion0(6,m)
        zeon0(6,m)=1.0_wp
     enddo
 
     call n_T_u_group_e               ! get Te  
!$omp parallel do private(m)
     do m=1,me
        zeon(4,m)=zeon(4,m)*msm0*c0cs*sqrt(backTe(m)/backT(m))  !checked; neglect electron rotation
        zeon(6,m)=zion(6,m)*sqrt(msm0*(backTe(m)/backT(m)))     ! =sqrt(mu_e)
     enddo

! the tracer electron is initially on PE=0 with index ntracer_e=me
     ntracer_e=0
     if(mype == 0)then
        ntracer_e=me
        zeon(1,ntracer_e)=zion(1,ntracer)   ! at the same location as ion tracer
        zeon(2,ntracer_e)=zion(2,ntracer)
        zeon(3,ntracer_e)=zion(3,ntracer)
        zeon(7,ntracer_e)=zion(7,ntracer)
        pdum=temperature_e(zeon(1,ntracer_e),0.0_wp)
        zeon(4,ntracer_e)=sqrt(pdum)/B_0*msm0*c0cs  ! using v_par=v_th for tracer electron
        zeon(4,ntracer_e)=zeon(4,ntracer_e)*0.4_wp
!        zeon(5,ntracer_e)=0.0_wp
        zeon(5,ntracer_e)=zeon(5,ntracer)
        zeon(6,ntracer_e)=sqrt(msm0*pdum/(2.0_wp*B_0)) ! using v_per=v_th for tracer electron
        if(ifdphi==2) zeon(8,ntracer_e)=0.0_wp

!!!!!!!!!TEST electron tracer using me/mi=1!!!!!!!!!!!!!!!!!!
!        zeon(:,ntracer_e)=zion(:,ntracer)

        zeon0(1:5,ntracer_e)=zeon(1:5,ntracer_e)
        zeon0(7,ntracer_e)=zeon(7,ntracer_e)
!        zeon0(6,ntracer_e)=zion0(6,ntracer)
        zeon0(6,ntracer_e)=1.0_wp
        if(ifdphi==2) zeon0(8,ntracer_e)=0.0_wp

        print*,'ntracer_e=',ntracer_e,zeon(1,ntracer_e),zeon(2,ntracer_e),zeon(3,ntracer_e),&
                       zeon(4,ntracer_e),zeon(5,ntracer_e),zeon(6,ntracer_e),zeon(7,ntracer_e)
   
        aaa(1)=zeon(1,ntracer_e)
        the(1)=zeon(2,ntracer_e)
        call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
        etracer_e0=0.5*(zeon(4,ntracer_e)*Bj(1))**2*m0ms*csc0**2+zeon(6,ntracer_e)**2*Bj(1)  ! kinetic energy

!        print*, 'TRACER Elec1', zeon(1,ntracer_e)
!        print*, 'TRACER Elec2', zeon(2,ntracer_e)
!        print*, 'TRACER Elec3', zeon(3,ntracer_e)
!        print*, 'TRACER Elec4', zeon(4,ntracer_e)
!        print*, 'TRACER Elec5', zeon(5,ntracer_e)
!        print*, 'TRACER Elec6', zeon(6,ntracer_e)*abs(c0cs)*sqrt(msm0)
     endif

!     zeon(5,:)=zeon(5,:)/10.0
!     zeon(5,:)=zeon(5,:)*1.0d2
!      zeon(5,:)=0.0
 
     markere=markeri*misum/mesum
     pmarke=pmarki*misum/mesum

     cloud_e=cloud_e/mesum

     my_eon=0.0_wp
!$omp parallel do private(m) reduction(+:my_eon)
     do m=1,me
        my_eon=my_eon+zeon(5,m)
     end do

     call MPI_Reduce(my_eon,total_eon,1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)

     if(mype==0) print*, 'INITIAL ELECTRON WEIGHT=', total_eon/mesum

  endif

end subroutine load
