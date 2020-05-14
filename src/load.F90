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
   real(wp) test_u

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
   if(nonlinear<0.5_wp) w_initial=w_initial_in*1.0e-10_wp

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
   acc_newer_H_zion=.true.

   mp1_i8=mp1   !Copy value of mp1 into 8-byte integer variable  mp1_i8
   ! Sum up the number of particles from all the processors using an 8-byte
   ! integer variable since the total may be larger than 2.1 billion.
   call MPI_ALLREDUCE(mp1_i8,misum_i8,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierror)

   iNrNm = iNrNm/real(misum_i8,wp)
   iNmNr = 1.0/iNrNm

   !     markeri is the inverse of ion marker number at each grid point in term of initial loading
   !  markeri=markeri*(real(misum_i8,wp)/real(ntoroidal,wp))
   markeri=markeri*real(misum_i8,wp)
   markeri=1.0_wp/markeri   ! to avoid divide operation

   !     pmarki is the inverse of total particle # (sum over nodes) at an anulus i
   pmarki=pmarki*real(misum_i8,wp)
   pmarki=1.0_wp/pmarki     ! to avoid divide operation

   rtden=rtden*misum/real(misum_i8,wp)
   mi=max(1,mp1)   ! really loaded particle number on my node
   ini_mi = mi
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
   acc_newer_H_zion=.true.

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
      ntracer_energy_ion0=0.5_wp*(zion(4,ntracer)*Bj(1))**2+zion(6,ntracer)**2*Bj(1)      ! kinetic energy
   endif

   do m=1,mi
      error_8=1.0e-5_wp
      IFAIL=1
      zion(7,m)=D01AHF(0.0_wp,zion(1,m),error_8,minpts,acc,dsda,-1,IFAIL)
      zion(7,m)=min(s1,max(s0,zion(7,m)))
   enddo
   acc_newer_H_zion=.true.


   ! Give unique label numbers for every particle
   do m=1,mi
      ! ilabels(m) = mype*10**(floor(log10(real(mi))) + 2) + m
      ilabels(m) = mype + 1e5*m
   enddo
   !$acc update device(ilabels(:mi))

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
      acc_newer_H_zion=.true.

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
   acc_newer_H_zion=.true.

   do m=1,mi
      wi_init(m)  = zion(5,m)
      wi_vpara(m) = 0
      wi_vperp(m) = 0
      wi_accel(m) = 0
      wi_fdamp(m) = 0
   enddo

   if(mype==0) print*, 'average weight at t=0 is ', pdum


   ! load electron on top of ion, including both trapped and untrapped electrons
   me=mi
   ini_me = me
   !   call MPI_ALLREDUCE(me,mesum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)
   mesum=misum

   eNmNr = (mesum/misum)*iNmNr ! Assuming real particle number (particle density) is the same for ions and electrons
   eNrNm = 1.0/eNmNr

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
      acc_newer_H_zeon=.true.

      do m=1,me
         we_init(m)  = zeon(5,m)
         we_vpara(m) = 0
         we_vperp(m) = 0
         we_accel(m) = 0
         we_adiab(m) = 0
         we_fdamp(m) = 0
      enddo



      call n_T_u_group_e               ! get Te
      !$omp parallel do private(m)
      do m=1,me
         zeon(4,m)=zeon(4,m)*msm0*c0cs*sqrt(backTe(m)/backT(m))  !checked; neglect electron rotation
         zeon(6,m)=zion(6,m)*sqrt(msm0*(backTe(m)/backT(m)))     ! =sqrt(mu_e)
      enddo
      acc_newer_H_zeon=.true.

      ! Give unique label numbers for every particle
      do m=1,me
         ! elabels(m) = mype*10**(floor(log10(real(me))) + 2) + m
         elabels(m) = mype +  1e5*m
      enddo
      !$acc update device(elabels(:me))


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
         zeon(4,ntracer_e)=zeon(4,ntracer_e)
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
         ntracer_energy_electron0=0.5*(zeon(4,ntracer_e)*Bj(1))**2*m0ms*csc0**2+zeon(6,ntracer_e)**2*Bj(1)  ! kinetic energy

         !        print*, 'TRACER Elec1', zeon(1,ntracer_e)
         !        print*, 'TRACER Elec2', zeon(2,ntracer_e)
         !        print*, 'TRACER Elec3', zeon(3,ntracer_e)
         !        print*, 'TRACER Elec4', zeon(4,ntracer_e)
         !        print*, 'TRACER Elec5', zeon(5,ntracer_e)
         !        print*, 'TRACER Elec6', zeon(6,ntracer_e)*abs(c0cs)*sqrt(msm0)
      endif
      acc_newer_H_zeon=.true.

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





   if(mype==0) then
      ilabel_of_ntracer    = ilabels(ntracer)
      elabel_of_ntracer_e  = elabels(ntracer_e)
   endif
   call MPI_BCAST(ilabel_of_ntracer,1,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
   call MPI_BCAST(elabel_of_ntracer_e,1,mpi_Rsize,0,MPI_COMM_WORLD,ierror)




   write(stdout,*) "In Loading.F90  I'm=",mype, "ini_mi=", ini_mi, ", ini_me=",ini_me
   call flush(stdout)

   if(mype==0)then
      write(stdout,*) "In Loading.F90, real loaded misum=", misum
   endif

   if(mi/= ini_mi) then
      write(stdout,*) "mi=", mi, ", ini_mi=",ini_mi, " , are different!!! Error!"
      call flush(stdout)
      call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
   endif

   ! (Min-Gu Yoo; 06/25/2018)
   ! For multiple tracers of ions to test Poincare Plot under electromagnetic perturbation



   if (mype==0) then
      ! Assign Electron tracer labels
      do i=1, num_of_etracers
         m = (real(i)/(num_of_etracers+1.0)) * ini_me
         etracers_labels(i) = elabels(m)

         ! !MGYOO (04/02/2019) : to test particles behavior depending on the initial energy
         ! ! test_u from 1 to 10
         ! test_u = 1 + (9/(num_of_etracers-1.0))*(i-1)
         ! zeon(1,m)=0.62
         ! zeon(2,m)=0.0_wp
         ! zeon(3,m)=1e-5
         ! pdum=temperature_e(zeon(1,m),0.0_wp)
         ! zeon(4,m)=test_u*sqrt(pdum)/B_0*msm0*c0cs      ! using v_par=v_th for tracer ion
         ! zeon(5,m)=0.1*w_initial;
         ! zeon(6,m)=test_u*sqrt(msm0*pdum)/sqrt(B_0)   ! using v_per=v_th for tracer ion
         ! zeon(6,m)=0.25*zeon(6,m) ! Reduce magnetic moment to make passing electrons

         ! zeon(7,m)=D01AHF(0.0_wp,zeon(1,m),error_8,minpts,acc,dsda,-1,IFAIL)
         ! zeon(7,m)=min(s1,max(s0,zeon(7,m)))
         ! zeon0(1:5,m)=zeon(1:5,m)
         ! zeon0(7,m)=zeon(7,m)

      enddo

      ! Assign Ion tracer labels
      do i=1, num_of_itracers
         m = (real(i)/(num_of_itracers+1.0)) * ini_mi

         itracers_labels(i) = ilabels(m)

         ! !MGYOO (04/02/2019) : to test particles behavior depending on the initial energy
         ! ! test_u from 1 to 10
         ! test_u = 1 + (9/(num_of_itracers-1.0))*(i-1)
         ! zion(1,m)=0.62
         ! zion(2,m)=0.0_wp
         ! zion(3,m)=1e-5
         ! pdum=temperature(zion(1,m),0.0_wp)
         ! backu(m)=u_parallel(zion(1,m),zion(2,m))
         ! zion(4,m)=test_u*(sqrt(pdum)+backu(m))/B_0      ! using v_par=v_th for tracer ion
         ! zion(5,m)=0.1*w_initial;
         ! zion(6,m)=test_u*sqrt(pdum)/sqrt(B_0)   ! using v_per=v_th for tracer ion
         ! zion(6,m)=0.25*zion(6,m) ! Reduce magnetic moment to make passing electrons

         ! zion(7,m)=D01AHF(0.0_wp,zion(1,m),error_8,minpts,acc,dsda,-1,IFAIL)
         ! zion(7,m)=min(s1,max(s0,zion(7,m)))
         ! zion0(1:5,m)=zion(1:5,m)
         ! zion0(7,m)=zion(7,m)


         ! zion(1,m)= i*((a1-a0)/(num_of_itracers+1.0)) + a0
         ! zion(2,m)=0.0_wp
         ! zion(3,m)=0.5_wp*(zetamin+zetamax)
      end do

      ! call esigetbf(zion,mi,7)   ! get B-->BB, I-->F  12/07
      ! call n_T_u_group   ! get T, n and omega_t
      ! backu(1:mi)=backu(1:mi)*F(1:mi)/BB(1:mi) ! =u_//=I*omega_t/B

      ! do i=1, num_of_itracers
      !   m = (real(i)/(num_of_itracers+1.0)) * ini_mi
      !   pdum=temperature(0.5,0.0_wp) ! same energy just for a test
      !   !pdum=temperature(zion(1,m),0.0_wp)
      !   !zion(4, m) = (sqrt(pdum)+backu(m))/BB(m)      ! using v_par=v_th for tracer ion

      !   ! zion(4, m) = (sqrt(pdum))/BB(m)      ! using v_par=v_th for tracer ion
      !   ! zion(6, m) = 0.1*sqrt(pdum)/sqrt(2.0_wp*BB(m))   ! using v_per=v_th for tracer ion

      !   ! Fixed mu and E test!!
      !   zion(6, m) = sqrt(0.1*(3.0/2.0)*pdum /B_0)  !! 0.1 of fixed total energy E = (3/2)*Ti
      !   zion(4, m) = sign(1.0, sprng()-0.5) * sqrt( 2.0/(BB(m)**2)  * ((3.0/2.0)*pdum  - zion(6,m)**2*BB(m) ) )

      !   zion0(:,m)=zion(:,m)
      ! enddo
   endif


   do m=1,mi
      wi_init(m)  = zion(5,m)
      UnNorm_theta_i(m) =zion(2,m)
      UnNorm_zeta_i(m)  =zion(3,m)
   enddo
   do m=1,me
      we_init(m)  = zeon(5,m)
      UnNorm_theta_e(m) =zeon(2,m)
      UnNorm_zeta_e(m)  =zeon(3,m)
   enddo

   call MPI_BCAST(etracers_labels,num_of_etracers,MPI_Rsize,0,MPI_COMM_WORLD,ierror)
   call MPI_BCAST(itracers_labels,num_of_itracers,MPI_Rsize,0,MPI_COMM_WORLD,ierror)

   ! write(stdout,*) "I'm=", mype,"etracers_labels(1:10)=", etracers_labels(1:10), "itracers_labels(1:10)=", itracers_labels(1:10)
   ! call flush(stdout)


end subroutine load



subroutine load_uniform_MGY
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use field_array
   use diagnosis_array
   use esi_interface
   use control_parameters
   implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

   integer :: minpts, ifail=1, ierror
   real(wp) :: acc, error_8=1e-5
   real(wp), external:: dsda, D01AHF

   real(wp) dela, r, fip, ave0
   real(wp) inv_Nm,Nm ! inverse of total number of Marker particles for each speceis (actually, they are the same)
   real(wp) each_ptl_d6z ! d6z(=da*dgq*dgf*dvpara*dvperp) is the micro phase space volume which each particle occupies
   real(wp) Te_norm(1:mi) !  normalized-Te=Te/E_0, the same unit as normalized Ti
   integer m,ip,i,j,ij
   real(wp) jacobian,density
   real(wp) dpsi,zdum,tdum,pdum
   real(wp) integral_f0
   real(wp) v_abs, v_theta, vpara, vperp ! v_abs is amplitude of total velocity vector, and v_theta is the angle between vpara and vperp

   integer(kind=8) :: sobol_seed=1, sobol_dim_num=5 ! 5D random number (a,gq,gf,v_abs,v_theta)
   real(kind=8) rand_5d(5,mi) ! can be sobol quasi random numbers, or can be sprng() numbers
   real(wp) inv_ds, ff, P_s, P_q  !inv_ds is uniform by definition of s_grids
   integer:: sid ! index for s_grid

   real(wp) vmax,vc,c1,c2,c3,c4,f2val,c5,c6,ratio,gi_val(mi),ge_val(me), Rand !For Flag_Load_and_WeightDef==2 case.

   if(irun==0) ddT_da=0.0_wp
   if(irun==0) sumfluxpsi=0.0_wp      ! sumfluxpsi(a)=int[I(a)dt]
   if(irun /= 0)then
      call restart_read
      return
   endif

   inv_ds = real(mpsi)/(s1-s0)


   me=mi

   Nm = (real(mi)*numberpe)
   inv_Nm=1.0/Nm

   if(mype==0) then
      write(stdout,*) "In load_uniform_MGY, inv_Nm = ",inv_Nm
      call flush(stdout)
   endif

   ! NOTE : markere and markeri definitions are changed
   ! They are now inverse of the (node volume*n0)
   do i=0,mpsi
      pdum=psi_grid(i)
      if(i==0) then
         dpsi=(psi_grid(1)-psi_grid(0))*0.5_wp         ! half cell
      else if(i==mpsi) then
         dpsi=(psi_grid(mpsi)-psi_grid(mpsi-1))*0.5_wp ! half cell
      else
         dpsi=(psi_grid(i+1)-psi_grid(i-1))*0.5_wp     ! full cell
      endif

      do j=0,mtheta(i)
         ij=igrid(i)+j
         zdum=zetamin+deltaz
         tdum=real(j,wp)*deltat(i)+zdum*qtinv(i)

         markeri(ij)=Jacobian(pdum,tdum)*density(pdum,tdum)*dpsi*deltat(i)*deltaz
      enddo
   enddo
   markeri = 1.0/(markeri)
   markere = markeri



   ! =================================================
   ! Ion Part
   ! =================================================

   ! Prepare sobol random numbers covering a whole device within a half-circle velocity space
   ! (a=(a0-a1), gq&gf=(0-2pi), v_abs=(0-umax)*sqrt(T/m), v_theta=(0-pi))
   ! Particles will be distributed into appropriate MPI ranks at the end of the initialization.
   sobol_seed = 1 + mype*mi
   do m=1,mi
      call i8_sobol(sobol_dim_num, sobol_seed, rand_5d(:,m))
   enddo

   ! Follwoing part is to prepare random numbers using sprng
   ! do m=1,mi
   !    rand_5d(1,m)=sprng()
   !    rand_5d(2,m)=sprng()
   !    rand_5d(3,m)=sprng()
   !    rand_5d(4,m)=sprng()
   !    rand_5d(5,m)=sprng()
   ! enddo

   do m=1,mi
      zion(1,m) = a0+ rand_5d(1,m)*(a1-a0)
      zion(7,m)=D01AHF(0.0_wp,zion(1,m),error_8,minpts,acc,dsda,-1,ifail)
      zion(7,m)=min(s1,max(s0,zion(7,m)))

      P_s = zion(7,m)
      sid = (P_s-s0)*inv_ds ! s index of original full grid
      ff  = (P_s-s_grid(sid))*inv_ds
      P_q = ff*q_grid(sid+1) + (1.0-ff)*q_grid(sid)

      ! Pick random points (theta,phi) with Helically aligned along B field
      zion(3,m) = 2.0*pi*rand_5d(2,m)
      zion(2,m) = mod(zion(3,m)/P_q + 2.0*pi*rand_5d(3,m), 2.0*pi)
   enddo
   acc_newer_H_zion=.true.

   call esigetbf(zion,mi,size(zion,1))
   call n_T_u_group  ! get T, n and omega_t, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
   backu(1:mi)=backu(1:mi)*F(1:mi)/BB(1:mi) ! =u_//=I*omega_t/B

   do m=1,mi
      ! ! Do Sampling for the following ranges
      ! ! -umax*sqrt(T/m) < u_para < umax*sqrt(T/m)
      ! ! -umax*sqrt(m*T)/(eB) < rho_para < umax*sqrt(m*T)/(eB)
      ! ! 0 < u_perp < umax*sqrt(2T/m)
      ! ! 0 < sqrt(mu) < sqrt(T/B)*u
      ! zion(4,m) = (2.0*rand_5d(4,m)-1.0)*umax *sqrt(backT(m))/BB(m)  ! u*sqrt(T)/BB + backu/BB
      ! if(ifelectron/=0) zeon(4,m)=zion(4,m)  ! This was for the quiet start..
      ! zion(4,m) = zion(4,m) + backu(m)/BB(m) ! add backu term
      ! zion(6,m) = sqrt(backT(m)/(BB(m)))*rand_5d(5,m)*umax ! sqrt(mu) = sqrt(m*uperp^2/2B)=sqrT(T/B)*u; mu=T/B*u^2
      ! ilabels(m) = mype + 1e5*m

      if(Flag_Load_and_WeightDef==1) then
         ! Find sampling points uniformly distributed within a half circle in the poloar coordinates (vpara,vperp)
         v_abs = umax*sqrt(backT(m))*sqrt(rand_5d(4,m)) ! v_r = amax*sqrt(Rand)
         v_theta = pi*rand_5d(5,m) ! For half-circle, theta should be 0<theta<2pi
      elseif(Flag_Load_and_WeightDef==2) then
         ! First, get v_abs using rand_5d(4,m)
         vmax=umax*sqrt(backT(m))
         vc=2.0*sqrt(backT(m))
         c2=0.01
         f2val=0.2
         c1 = -((-2+c2*vmax+f2val*vmax+ c2*vc-f2val*vc)/( vmax*vc))
         c3 = -((2-2*f2val*vmax-c2*vc+ f2val*vc)/( vmax*(vmax-vc)))
         c4 = -((-2+f2val*vmax+c2*vc)/(vmax-vc))
         Rand = rand_5d(4,m)
         if (Rand < c2*vc+c1*vc**2/2) then
            v_abs = (-c2+sqrt(c2**2+2*c1*Rand))/c1
         else
            Rand = Rand -(c2*vc+c1*vc**2/2)
            v_abs = (-c4 + sqrt(2*c3*Rand+(c4+c3*vc)**2))/c3
         endif

         ! Second, get v_theta using rand_5d(5,m)
         ratio=0.2 ! ratio of (max/min)
         c5=1.0/((1.0+ratio)*0.5*pi)
         c6=(1.0-ratio)/((1.0+ratio)*(0.5*pi)**2)
         Rand = rand_5d(5,m)
         v_theta = 0.5*pi + sign(1.0_wp,Rand-0.5)*((c5-sqrt(c5**2-2*c6*abs(Rand-0.5)))/c6);

         if(v_abs<vc) then
            gi_val(m) = Nm*(c1*v_abs+c2)*(c5-c6*abs(v_theta-0.5*pi))/(Jacobian(zion(1,m),zion(2,m))*v_abs)/((a1-a0)*4.0*pi*pi)
         else
            gi_val(m) = Nm*(c3*v_abs+c4)*(c5-c6*abs(v_theta-0.5*pi))/(Jacobian(zion(1,m),zion(2,m))*v_abs)/((a1-a0)*4.0*pi*pi)
         endif
      endif

      vpara=v_abs*cos(v_theta)
      vperp=v_abs*sin(v_theta)

      zion(4,m) = vpara/BB(m)
      if(ifelectron/=0) zeon(4,m)=zion(4,m)  ! This is for the quiet start..
      zion(4,m) = zion(4,m) + backu(m)/BB(m) ! add backu term
      zion(6,m) = vperp*sqrt(1.0/(2.0*BB(m)))
      ilabels(m) = mype + 1e5*m
   enddo
   acc_newer_H_zion=.true.

   !$acc update device(ilabels(:mi))


   call n_T_u_group  ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))

   ! Set initial weights of ions
   do m=1,mi
      ! zion(5,m) is "df"
      ! zion0(6,m) is "f_tot"


      if(Flag_Load_and_WeightDef==1) then
         !Following line is for the half-circle sampling in the polar coordinates (vpara, vperp)
         ! inv_Nm * (a1-a0) * (2pi) * (2pi) * (0.5*pi*v_abs_max^2)
         each_ptl_d6z = inv_Nm*(a1-a0)*4.0*pi*pi*0.5*pi*umax*umax*backT(m)
         integral_f0 = Jacobian(zion(1,m),zion(2,m))*each_ptl_d6z*back_f0_i(m)
      elseif(Flag_Load_and_WeightDef==2)then
         integral_f0 = back_f0_i(m)/gi_val(m)
      endif

      zion(5,m)  = (2.0*sprng()-1.0)*(1.0+cos(zion(2,m)))*w_initial_in*integral_f0
      zion0(6,m) = integral_f0 + zion(5,m)
   enddo
   acc_newer_H_zion=.true.

   wi_init  = zion(5,:)
   wi_vpara = 0
   wi_vperp = 0
   wi_accel = 0
   wi_fdamp = 0


   ! =================================================
   ! Electron Part
   ! =================================================

   ! Prepare sobol random numbers covering a whole device within a half-circle velocity space
   ! (a=(a0-a1), gq&gf=(0-2pi), v_abs=(0-umax)*sqrt(T/m), v_theta=(0-pi))
   ! Particles will be distributed into appropriate MPI ranks at the end of the initialization.
   sobol_seed = 1 + mype*me
   do m=1,me
      call i8_sobol(sobol_dim_num, sobol_seed, rand_5d(:,m))
   enddo

   ! Follwoing part is to prepare random numbers using sprng
   ! do m=1,me
   !    rand_5d(1,m)=sprng()
   !    rand_5d(2,m)=sprng()
   !    rand_5d(3,m)=sprng()
   !    rand_5d(4,m)=sprng()
   !    rand_5d(5,m)=sprng()
   ! enddo

   if(ifelectron/=0) then
      if(Flag_Quiet_Start==1) then
         do m=1,me
            zeon(1,m)=zion(1,m)
            zeon(2,m)=zion(2,m)
            zeon(3,m)=zion(3,m)
            zeon(7,m)=zion(7,m)
            elabels(m) = mype +  1e5*m
         enddo
      else
         do m=1,me
            zeon(1,m) = a0+ rand_5d(1,m)*(a1-a0)
            zeon(7,m)=D01AHF(0.0_wp,zeon(1,m),error_8,minpts,acc,dsda,-1,ifail)
            zeon(7,m)=min(s1,max(s0,zeon(7,m)))

            P_s = zeon(7,m)
            sid = (P_s-s0)*inv_ds ! s index of original full grid
            ff  = (P_s-s_grid(sid))*inv_ds
            P_q = ff*q_grid(sid+1) + (1.0-ff)*q_grid(sid)

            ! Pick random points (theta,phi) with Helically aligned along B field
            zeon(3,m) = 2.0*pi*rand_5d(2,m)
            zeon(2,m) = mod(zeon(3,m)/P_q + 2.0*pi*rand_5d(3,m), 2.0*pi)

            elabels(m) = mype +  1e5*m
         enddo
      endif
      acc_newer_H_zeon=.true.

      call n_T_u_group_e   ! get Te, note that back_f0_e cannot be calculated here because it requires zion(4) and zion(6)

      do m=1,me
         Te_norm(m) = backTe(m)*msm0
      enddo

      if(Flag_Quiet_Start==1) then
         do m=1,me
            zeon(4,m)=zion(4,m)*sqrt(msm0)*c0cs*sqrt(Te_norm(m)/backT(m))  !checked; neglect electron rotation
            zeon(6,m)=zion(6,m)*sqrt(Te_norm(m)/backT(m))     ! =sqrt(mu_e)
         enddo
      else
         do m=1,me
            ! ! Do Sampling for the following ranges
            ! ! -umax*sqrt(T/m) < u_para < umax*sqrt(T/m)
            ! ! -umax*sqrt(m*T)/(eB) < rho_para < umax*sqrt(m*T)/(eB)
            ! ! 0 < u_perp < umax*sqrt(2T/m)
            ! ! 0 < sqrt(mu) < sqrt(T/B)*u
            ! zeon(4,m) = (2.0*rand_5d(4,m)-1.0)*umax *sqrt(msm0*Te_norm(m))/(csc0*BB(m)) !neglected electron rotation
            ! zeon(6,m) = sqrt(Te_norm(m)/BB(m))*rand_5d(5,m)*umax

            if(Flag_Load_and_WeightDef==1) then
               ! Find sampling points uniformly distributed within a half circle in the poloar coordinates (vpara,vperp)
               v_abs = umax*sqrt(Te_norm(m)*m0ms)*sqrt(rand_5d(4,m)) ! v_r = amax*sqrt(Rand)
               v_theta = pi*rand_5d(5,m) ! For half-circle, theta should be 0<theta<2pi
            elseif(Flag_Load_and_WeightDef==2) then
               ! First, get v_abs using rand_5d(4,m)
               vmax=umax*sqrt(Te_norm(m)*m0ms)
               vc=2.0*sqrt(Te_norm(m)*m0ms)
               c2=0.01
               f2val=0.2
               c1 = -((-2+c2*vmax+f2val*vmax+ c2*vc-f2val*vc)/( vmax*vc))
               c3 = -((2-2*f2val*vmax-c2*vc+ f2val*vc)/( vmax*(vmax-vc)))
               c4 = -((-2+f2val*vmax+c2*vc)/(vmax-vc))
               Rand = rand_5d(4,m)
               if (Rand < c2*vc+c1*vc**2/2) then
                  v_abs = (-c2+sqrt(c2**2+2*c1*Rand))/c1
               else
                  Rand = Rand -(c2*vc+c1*vc**2/2)
                  v_abs = (-c4 + sqrt(2*c3*Rand+(c4+c3*vc)**2))/c3
               endif
               ! Second, get v_theta using rand_5d(5,m)
               ratio=0.2 ! ratio of (max/min)
               c5=1.0/((1.0+ratio)*0.5*pi)
               c6=(1.0-ratio)/((1.0+ratio)*(0.5*pi)**2)
               Rand = rand_5d(5,m)
               v_theta = 0.5*pi + sign(1.0_wp,Rand-0.5)*((c5-sqrt(c5**2-2*c6*abs(Rand-0.5)))/c6);

               if(v_abs<vc) then
                  ge_val(m) = Nm*(c1*v_abs+c2)*(c5-c6*abs(v_theta-0.5*pi))/(Jacobian(zeon(1,m),zeon(2,m))*v_abs)/((a1-a0)*4.0*pi*pi)
               else
                  ge_val(m) = Nm*(c3*v_abs+c4)*(c5-c6*abs(v_theta-0.5*pi))/(Jacobian(zeon(1,m),zeon(2,m))*v_abs)/((a1-a0)*4.0*pi*pi)
               endif
            endif

            vpara=v_abs*cos(v_theta)
            vperp=v_abs*sin(v_theta)

            zeon(4,m) = msm0*vpara/(csc0*BB(m))
            zeon(6,m) = vperp*sqrt(msm0/(2.0*BB(m)))
         enddo
      endif
      acc_newer_H_zeon=.true.

      call n_T_u_group_e   ! This is to get back_f0_e inside. (which requires zeon(4) and zeon(6))

      if(Flag_Quiet_Start==1) then
         zeon(5,1:me) = zion(5,1:me)
         zeon0(6,1:me) = zion0(6,1:me)
      else
         ! Set initial weights of electrons
         do m=1,me
            ! zeon(5,m) is "df"
            ! zeon0(6,m) is "f_tot"

            if(Flag_Load_and_WeightDef==1) then
               !Following line is for the half-circle sampling in the polar coordinates (vpara, vperp)
               ! inv_Nm * (a1-a0) * (2pi) * (2pi) * (0.5*pi*v_abs_max^2)
               each_ptl_d6z = inv_Nm*(a1-a0)*4.0*pi*pi*0.5*pi*umax*umax*Te_norm(m)*m0ms
               integral_f0 = Jacobian(zeon(1,m),zeon(2,m))*each_ptl_d6z*back_f0_e(m)
            elseif(Flag_Load_and_WeightDef==2)then
               integral_f0 = back_f0_e(m)/ge_val(m)
            endif
            zeon(5,m)  = (2.0*sprng()-1.0)*(1.0+cos(zeon(2,m)))*w_initial_in*integral_f0
            zeon0(6,m) = integral_f0 + zeon(5,m)
         enddo
      endif
      acc_newer_H_zeon=.true.

      we_init(:)  = zeon(5,:)
      we_vpara = 0
      we_vperp = 0
      we_accel = 0
      we_adiab = 0
      we_fdamp = 0
   endif

   !$acc update device(elabels(:me))



   ! (Min-Gu Yoo; 06/25/2018)
   ! For multiple tracers of ions to test Poincare Plot under electromagnetic perturbation
   if (mype==0) then
      ! Assign Electron tracer labels
      do i=1, num_of_etracers
         m = (real(i)/(num_of_etracers+1.0)) * me
         etracers_labels(i) = elabels(m)
      enddo

      ! Assign Ion tracer labels
      do i=1, num_of_itracers
         m = (real(i)/(num_of_itracers+1.0)) * mi
         itracers_labels(i) = ilabels(m)
      end do

      ! For diagnostic purpose
      if(mype==0) then
         m = 0.5*ini_mi
         write(stdout,*) m,"-th ion's each_ptl_d6z=", inv_Nm*((a1-a0)*4.0*pi*pi*2.0*umax*sqrt(2.0)*umax*backT(m)), ", back_f0_i(m)=",back_f0_i(m)
         call flush(stdout)

         m = 0.5*ini_me
         write(stdout,*) m,"-th electron's each_ptl_d6z=", inv_Nm*((a1-a0)*4.0*pi*pi*2.0*umax*sqrt(2.0)*umax*Te_norm(m)*m0ms),  ", back_f0_e(m)=",back_f0_e(m)
         call flush(stdout)
      endif
   endif


   do m=1,mi
      wi_init(m)  = zion(5,m)
      UnNorm_theta_i(m) =zion(2,m)
      UnNorm_zeta_i(m)  =zion(3,m)
   enddo
   do m=1,me
      we_init(m)  = zeon(5,m)
      UnNorm_theta_e(m) =zeon(2,m)
      UnNorm_zeta_e(m)  =zeon(3,m)
   enddo

   call MPI_BCAST(etracers_labels,num_of_etracers,MPI_Rsize,0,MPI_COMM_WORLD,ierror)
   call MPI_BCAST(itracers_labels,num_of_itracers,MPI_Rsize,0,MPI_COMM_WORLD,ierror)


   ! Distribute ions and electrons to appropriate MPI ranks
   call shifti
   call shifte

   ini_mi=mi
   ini_me=me

   write(stdout,*) "In load_uniform_MGY, I'm=",mype, "mi=", mi, ", me=",me
   call flush(stdout)

   ntracer=0
   ntracer_e=0
   if(mype==0) then
      ntracer=mi
      ntracer_e=me
      ilabel_of_ntracer    = ilabels(ntracer)
      elabel_of_ntracer_e  = elabels(ntracer_e)
   endif
   call MPI_BCAST(ilabel_of_ntracer,1,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
   call MPI_BCAST(elabel_of_ntracer_e,1,mpi_Rsize,0,MPI_COMM_WORLD,ierror)


end subroutine load_uniform_MGY



module Ptl_locators_partial

   implicit none

   contains

      subroutine electron_locator_partial(m_start,m_end)
         use precision
         use global_parameters
         use particle_array
         use field_array
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0,m1


         integer m,im,ii,j01,j00
         real(wp) rdum,tdum,delr,delt(0:mpsi),delz,wz1,wp1,     &
            pi2_inv,psitmp,thetatmp,zetatmp,psimax,psimin

         m0=1
         m1=me
         if(present(m_start)) m0=m_start
         if(present(m_end))   m1=m_end


         delr=1.0_wp/deltar
         delt=2.0_wp*pi/deltat
         delz=1.0_wp/deltaz
         pi2_inv=0.5_wp/pi
         psimin=s_grid(0)
         psimax=s_grid(mpsi)


         if(acc_newer_H_zeon) then
            !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_H_zeon=.false.
            acc_newer_D_zeon=.false.
         endif

         !$acc parallel loop present(zeon,jteon0,jteon1,wzeon,wpeon,wteon0,wteon1)
         do m=m0,m1
            !     psitmp=zeon(1,m)      ! =a
            psitmp=zeon(7,m)       ! =s
            thetatmp=zeon(2,m)
            zetatmp=zeon(3,m)

            wzeon(m)=(zetatmp-zetamin)*delz

            rdum=psitmp
            rdum=delr*max(0.0_wp,min(s1-s0,rdum-s0))
            ii=max(0,min(mpsi-1,floor(rdum)))
            wp1=rdum-real(ii,wp)
            wpeon(m)=wp1

            ! inner flux surface
            im=ii
            tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
            tdum=(tdum-floor(tdum))*delt(im)
            j00=max(0,min(mtheta(im)-1,floor(tdum)))
            jteon0(m)=igrid(im)+j00  ! pay attention to j00=0 case
            wteon0(m)=tdum-real(j00,wp)

            ! outer flux surface
            im=ii+1
            tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
            tdum=(tdum-floor(tdum))*delt(im)
            j01=max(0,min(mtheta(im)-1,floor(tdum)))
            jteon1(m)=igrid(im)+j01  ! pay attention to j01=0 case
            wteon1(m)=tdum-real(j01,wp)
         enddo

      end subroutine electron_locator_partial

      subroutine ion_locator_at_larmor_points_partial(m_start,m_end)

         use precision
         use global_parameters
         use particle_array
         use field_array
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0,m1

         integer m,i,im,j,k,ip,jt,ii,j11,j10,j01,j00,larmor,ij,ipjt

         real(wp) rdum,tdum,delr,delt(0:mpsi),delz,r,wp1,wp0,wt11,wt10,wt01,wt00,&
            tflr,pi2_inv,psitmp,thetatmp,zetatmp,rhoi,psimax,psimin,pdum,thedum


         integer inner_ip,outter_ip, inner_jt, outter_jt
         real(wp) gyro_phase, intpl_pgyro(larmor_num), intpl_tgyro(larmor_num)


         m0=1
         m1=mi
         if(present(m_start)) m0=m_start
         if(present(m_end))   m1=m_end


         delr=1.0_wp/deltar
         delt=2.0_wp*pi/deltat
         delz=1.0_wp/deltaz
         pi2_inv=0.5_wp/pi
         !  psimin=psi_grid(0)
         !  psimax=psi_grid(mpsi)
         psimin=s_grid(0)
         psimax=s_grid(mpsi)


         if(acc_newer_H_zion) then
            !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
            acc_newer_H_zion=.false.
            acc_newer_D_zion=.false.
         endif
         !$acc parallel loop present(zion,jtion0,jtion1,wzion,wpion,wtion0,wtion1,zion_s_on_GyroOrbit,zion_gq_on_GyroOrbit) private(intpl_pgyro,intpl_tgyro)
         do m=m0,m1
            !psitmp=zion(1,m)   ! =a
            psitmp=zion(7,m)   ! =s
            thetatmp=zion(2,m)
            zetatmp=zion(3,m)
            rhoi=zion(6,m)

            r=psitmp
            ip=max(0,min(mpsi,floor((r-s0)*delr+0.5)))  ! grid size is uniform in s but not a
            !   jt=max(0,min(mtheta(ip),floor(thetatmp*pi2_inv*delt(ip)+0.5)))
            thedum=thetatmp-zetatmp*qtinv(ip)    ! corresponding theta on zeta=0 plane
            thedum=thedum*pi2_inv+100.0_wp
            thedum=2.0_wp*pi*(thedum-floor(thedum))       ! bring back to range [0,2pi]
            jt=max(0,min(mtheta(ip),floor(thedum*pi2_inv*delt(ip)+0.5)))
            ipjt=igrid(ip)+jt     ! note special case of jt=0
            wzion(m)=(zetatmp-zetamin)*delz


            !==========================================
            ! MGYOO (10/24/2018) : to store information of GC position.
            ! GC's lower-left cell-index(jtion) and linear-interpolation-weights (wpion, wtion0, wtion1)
            inner_ip = max(0,min(mpsi-1,floor((r-s0)*delr)))
            thedum=zion(2,m)-zetatmp*qtinv(inner_ip)              ! corresponding theta on zeta=0 plane
            thedum=modulo(thedum,2.0_wp*pi)/deltat(inner_ip)  ! bring back to range [0,2pi] and divde by delta_theta
            inner_jt = floor(thedum)
            wtion0(0,m) = thedum - floor(thedum)
            jtion0(0,m) = igrid(inner_ip) + inner_jt
            wpion(0,m) =  (zion(7,m) - s_grid(inner_ip))/(s_grid(inner_ip+1)-s_grid(inner_ip))

            outter_ip = inner_ip+1
            thedum=zion(2,m)-zetatmp*qtinv(outter_ip)     ! corresponding theta on zeta=0 plane
            thedum=modulo(thedum,2.0_wp*pi)/deltat(outter_ip)  ! bring back to range [0,2pi] and divde by delta_theta
            outter_jt = floor(thedum)
            wtion1(0,m) = thedum - floor(thedum)
            jtion1(0,m) = igrid(outter_ip) + outter_jt
            !==========================================


            !==========================================
            ! MGYOO (10/25/2018) : to get linear-interpolation of pgyro and tgyro
            wp0 = 1.0_wp - wpion(0,m)
            wp1 = wpion(0,m)

            wt00 = 1.0_wp - wtion0(0,m)
            wt10 = wtion0(0,m)

            wt01 = 1.0_wp - wtion1(0,m)
            wt11 = wtion1(0,m)

            !$acc loop seq
            do larmor=1, larmor_num
               intpl_pgyro(larmor) = wp0*(wt00*pgyro(larmor,jtion0(0,m)) + wt10*pgyro(larmor,jtion0(0,m)+1)) + wp1*(wt01*pgyro(larmor,jtion1(0,m)) + wt11*pgyro(larmor,jtion1(0,m)+1))
               intpl_tgyro(larmor) = wp0*(wt00*tgyro(larmor,jtion0(0,m)) + wt10*tgyro(larmor,jtion0(0,m)+1)) + wp1*(wt01*tgyro(larmor,jtion1(0,m)) + wt11*tgyro(larmor,jtion1(0,m)+1))
            enddo
            !==========================================

            !$acc loop seq
            do larmor=1, larmor_num
               !  MODIFIED for FINITE GYRORADIUS CASE
               pdum=max(psimin,min(psimax,psitmp+rhoi*intpl_pgyro(larmor)))
               zion_s_on_GyroOrbit(larmor,m)=pdum
               rdum=delr*max(0.0_wp,min(s1-s0,pdum-s0))
               ii=max(0,min(mpsi-1,floor(rdum)))
               !        wp1=rdum-real(ii)
               !        wpion(larmor,m)=wp1
               wp1=(pdum-s_grid(ii))/(s_grid(ii+1)-s_grid(ii))   !wx  finite rho_i case
               wpion(larmor,m)=wp1                               !wx  finite rho_i case

               ! particle position in theta
               tflr=thetatmp+rhoi*intpl_tgyro(larmor)
               zion_gq_on_GyroOrbit(larmor,m)=tflr

               ! inner flux surface
               im=ii
               tdum=pi2_inv*(tflr-zetatmp*qtinv(im))+10.0_wp
               tdum=(tdum-floor(tdum))*delt(im)
               j00=max(0,min(mtheta(im)-1,floor(tdum)))
               jtion0(larmor,m)=igrid(im)+j00  ! pay attention to j00=0 case
               wtion0(larmor,m)=tdum-real(j00)

               ! outer flux surface
               im=ii+1
               tdum=pi2_inv*(tflr-zetatmp*qtinv(im))+10.0_wp
               tdum=(tdum-floor(tdum))*delt(im)
               j01=max(0,min(mtheta(im)-1,floor(tdum)))
               jtion1(larmor,m)=igrid(im)+j01  ! pay attention to j01=0 case
               wtion1(larmor,m)=tdum-real(j01)
            enddo
         enddo

      end subroutine ion_locator_at_larmor_points_partial


end module Ptl_locators_partial

module Ptl_Src_Loss

   use Ptl_locators_partial
   use precision
   implicit none

   ! Following numbers are local things estimated and controlled by each rank.
   ! They are not global things for all ranks (MPI_COMM_WORLD)
   integer :: e_lost_num_BDY0=0, e_lost_num_BDY1=0 ! number of loss particles across boundary per each call of treat_electrons_out_of_boundary
   integer :: e_cum_lost_num_BDY0=0, e_cum_lost_num_BDY1=0 ! cumulative number of lost particles across boundary during a period of time.

   integer :: i_lost_num_BDY0=0, i_lost_num_BDY1=0 ! number of loss particles across boundary at each istep and per each rank
   integer :: i_cum_lost_num_BDY0=0, i_cum_lost_num_BDY1=0 ! Sum of number of loss particles across boundary over all ranks and isteps.

   type FootPrint
      integer :: dim_gq=50, dim_gf=50 ! Default dimension of gq and gf
      real(wp) inv_dgq, inv_dgf
      real(wp) a, s
      integer :: num_eSD, num_iSD
      real(wp),allocatable, dimension(:) :: gq_1D, gf_1D
      real(wp),allocatable, dimension(:,:) :: a_2D, gq_2D, gf_2D, R_2D, Z_2D, BB_2D
      real(wp),allocatable, dimension(:,:,:) :: eScaData, iScaData ! ScalarData(density, temperautre) of electrons and ions
   end type

   type(FootPrint) FP_BDY1, cum_FP_BDY1

   contains


      subroutine Initialize_FootPrint(FP,target_a)
         use global_parameters
         implicit none

         type(FootPrint), intent(inout) :: FP
         real(wp),intent(in) :: target_a

         integer Nq,Nf ! Nq=FP%dim_gq, Nf=FP%dim_gf
         integer qq,ff

         FP%dim_gq = 200
         FP%dim_gf = 100

         Nq=FP%dim_gq
         Nf=FP%dim_gf

         FP%inv_dgq = real(Nq-1)/(2.0*pi)
         FP%inv_dgf = real(Nf-1)/(2.0*pi)


         FP%num_eSD = 6 !Marker, Particle Flux, Heat Flux
         FP%num_iSD = 6

         allocate(FP%gq_1D(Nq), FP%gf_1D(Nf))
         allocate(FP%a_2D(Nq,Nf), FP%gq_2D(Nq,Nf), FP%gf_2D(Nq,Nf))
         allocate(FP%R_2D(Nq,Nf), FP%Z_2D(Nq,Nf), FP%BB_2D(Nq,Nf))
         allocate(FP%eScaData(Nq,Nf,FP%num_eSD), FP%iScaData(Nq,Nf,FP%num_iSD))

         FP%eScaData=0.0
         FP%iScaData=0.0

         FP%a_2D = target_a

         do qq=1,Nq
            FP%gq_1D(qq) = (2.0*pi)/(Nq-1)*(qq-1)
            FP%gq_2D(qq,:) = FP%gq_1D(qq)
         enddo

         do ff=1,Nf
            FP%gf_1D(ff) = (2.0*pi)/(Nf-1)*(ff-1)
            FP%gf_2D(:,ff) = FP%gf_1D(ff)
         enddo

         call esigetrzb(FP%R_2D, FP%Z_2D, FP%BB_2D, FP%a_2D, FP%gq_2D, Nq*Nf)

      end subroutine

      subroutine Measure_and_Cumulate_Electron_FootPrint_at_BDY1(FP, energy)
         use global_parameters
         use particle_array
         use field_array,only:q_grid
         use mpi_defs
         implicit none

         type(FootPrint), intent(inout) :: FP
         real(wp), intent(in) :: energy(me)

         real(wp) ratio, wgq0,wgq1,wgf0,wgf1
         real(wp) P_gq, P_gf, P_safety_q
         real(wp) target_gq, target_gf
         integer gqid, gfid, qq, ff

         real(wp) zeon_a_arr(me)
         integer lost_mid_arr
         real(wp), allocatable,dimension(:,:) :: out_zeon, out_zeon0
         real(wp), allocatable,dimension(:) :: delta_f, f_tot, out_energy, out_backTe, out_back_f0_e
         integer m,m_out,vid
         real(wp) time0,time1

         allocate(out_zeon(size(zeon,1),me),out_zeon0(size(zeon0,1),me))
         allocate(out_energy(me),out_backTe(me),out_back_f0_e(me))


#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

         !$acc update self(zeon(:,:me),zeon0(:,:me),backTe(:me),back_f0_e(:me))

         m_out=0
         do m=1,me
            if(zeon(1,m)>a1) then
               m_out=m_out+1
               out_zeon(:,m_out)  = zeon(:,m)
               out_zeon0(:,m_out) = zeon0(:,m)
               out_energy(m_out) = energy(m)
               out_backTe(m_out) = backTe(m)
               out_back_f0_e(m_out) = back_f0_e(m)
            endif
         enddo

         allocate(delta_f(m_out),f_tot(m_out))


         if(Flag_Load_and_WeightDef==0)then
            delta_f(:) = out_zeon(5,1:m_out)*eNrNm
            f_tot(:)   = eNrNm
         else
            delta_f(:) = out_zeon(5,1:m_out)
            f_tot(:)   = out_zeon0(6,1:m_out)
         endif


         do m=1,m_out

            P_gq = modulo(out_zeon(2,m),2.0*pi)
            P_gf = modulo(out_zeon(3,m),2.0*pi)

            P_safety_q = q_grid(mpsi)

            gfid  = (P_gf)*FP%inv_dgf +1 ! last +1 is because of Fortran-style indexing
            wgf1  = (P_gf-FP%gf_1D(gfid))*FP%inv_dgf
            wgf0  = 1.0-wgf1


            ! (At the plane of gfid): the nearst left plane in toroidal direction
            target_gf = FP%gf_1D(gfid)
            target_gq = modulo(P_gq - (P_gf-target_gf)/P_safety_q, 2.0*pi)
            gqid  = (target_gq)*FP%inv_dgq +1
            wgq1  = (target_gq-FP%gq_1D(gqid))*FP%inv_dgq
            wgq0  = 1.0-wgq1

            ratio = wgf0*wgq0; ff=gfid; qq=gqid
            vid=1; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backTe(m)*msm0)*delta_f(m) !Temperature definition 1
            vid=3; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot


            ratio = wgf0*wgq1; ff=gfid; qq=gqid+1
            vid=1; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backTe(m)*msm0)*delta_f(m) !Temperature definition 1
            vid=3; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot


            ! (At the plane of gfid+1): the nearst right plane in toroidal direction
            target_gf = FP%gf_1D(gfid+1)
            target_gq = modulo(P_gq - (P_gf-target_gf)/P_safety_q, 2.0*pi)
            gqid  = (target_gq)*FP%inv_dgq +1
            wgq1  = (target_gq-FP%gq_1D(gqid))*FP%inv_dgq
            wgq0  = 1.0-wgq1

            ratio = wgf1*wgq0; ff=gfid+1; qq=gqid
            vid=1; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backTe(m)*msm0)*delta_f(m) !Temperature definition 1
            vid=3; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot

            ratio = wgf1*wgq1; ff=gfid+1; qq=gqid+1
            vid=1; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backTe(m)*msm0)*delta_f(m) !Temperature definition 1
            vid=3; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot
         enddo



         deallocate(delta_f,f_tot)
         deallocate(out_zeon,out_zeon0,out_energy,out_backTe,out_back_f0_e)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Measure_Electron_FootPrint_BDY1 takes", time1-time0, " secs"
         call flush(stdout)
#endif

      end subroutine Measure_and_Cumulate_Electron_FootPrint_at_BDY1

      subroutine Measure_and_Cumulate_Ion_FootPrint_at_BDY1(FP, energy)
         use global_parameters
         use particle_array
         use field_array,only:q_grid
         use mpi_defs
         implicit none

         type(FootPrint), intent(inout) :: FP
         real(wp), intent(in) :: energy(mi)

         real(wp) ratio, wgq0,wgq1,wgf0,wgf1
         real(wp) P_gq, P_gf, P_safety_q
         real(wp) target_gq, target_gf
         integer gqid, gfid, qq, ff

         real(wp) zion_a_arr(mi)
         integer lost_mid_arr
         real(wp), allocatable,dimension(:,:) :: out_zion, out_zion0
         real(wp), allocatable,dimension(:) :: delta_f, f_tot, out_energy, out_backT, out_back_f0_i
         integer m,m_out,vid
         real(wp) time0,time1

         allocate(out_zion(size(zion,1),mi),out_zion0(size(zion0,1),mi))
         allocate(out_energy(mi),out_backT(mi),out_back_f0_i(mi))


#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

         !$acc update self(zion(:,:mi),zion0(:,:mi),backT(:mi),back_f0_i(:mi))

         m_out=0
         do m=1,mi
            if(zion(1,m)>a1) then
               m_out=m_out+1
               out_zion(:,m_out)  = zion(:,m)
               out_zion0(:,m_out) = zion0(:,m)
               out_energy(m_out) = energy(m)
               out_backT(m_out) = backT(m)
               out_back_f0_i(m_out) = back_f0_i(m)
            endif
         enddo

         allocate(delta_f(m_out),f_tot(m_out))


         if(Flag_Load_and_WeightDef==0)then
            delta_f(:) = out_zion(5,1:m_out)*iNrNm
            f_tot(:)   = iNrNm
         else
            delta_f(:) = out_zion(5,1:m_out)
            f_tot(:)   = out_zion0(6,1:m_out)
         endif


         do m=1,m_out

            P_gq = modulo(out_zion(2,m),2.0*pi)
            P_gf = modulo(out_zion(3,m),2.0*pi)

            P_safety_q = q_grid(mpsi)

            gfid  = (P_gf)*FP%inv_dgf +1 ! last +1 is because of Fortran-style indexing
            wgf1  = (P_gf-FP%gf_1D(gfid))*FP%inv_dgf
            wgf0  = 1.0-wgf1


            ! (At the plane of gfid): the nearst left plane in toroidal direction
            target_gf = FP%gf_1D(gfid)
            target_gq = modulo(P_gq - (P_gf-target_gf)/P_safety_q, 2.0*pi)
            gqid  = (target_gq)*FP%inv_dgq +1
            wgq1  = (target_gq-FP%gq_1D(gqid))*FP%inv_dgq
            wgq0  = 1.0-wgq1

            ratio = wgf0*wgq0; ff=gfid; qq=gqid
            vid=1; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backT(m))*delta_f(m) !Temperature definition 1
            vid=3; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot


            ratio = wgf0*wgq1; ff=gfid; qq=gqid+1
            vid=1; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backT(m))*delta_f(m) !Temperature definition 1
            vid=3; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot


            ! (At the plane of gfid+1): the nearst right plane in toroidal direction
            target_gf = FP%gf_1D(gfid+1)
            target_gq = modulo(P_gq - (P_gf-target_gf)/P_safety_q, 2.0*pi)
            gqid  = (target_gq)*FP%inv_dgq +1
            wgq1  = (target_gq-FP%gq_1D(gqid))*FP%inv_dgq
            wgq0  = 1.0-wgq1

            ratio = wgf1*wgq0; ff=gfid+1; qq=gqid
            vid=1; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backT(m))*delta_f(m) !Temperature definition 1
            vid=3; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot

            ratio = wgf1*wgq1; ff=gfid+1; qq=gqid+1
            vid=1; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backT(m))*delta_f(m) !Temperature definition 1
            vid=3; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot
         enddo



         deallocate(delta_f,f_tot)
         deallocate(out_zion,out_zion0,out_energy,out_backT,out_back_f0_i)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Measure_Ion_FootPrint_BDY1 takes", time1-time0, " secs"
         call flush(stdout)
#endif

      end subroutine Measure_and_Cumulate_Ion_FootPrint_at_BDY1

      subroutine Clear_FootPrints_of_Electrons_and_Ions(FP)
         implicit none
         type(FootPrint),intent(inout) :: FP

         FP%eScaData=0
         FP%iScaData=0
      end subroutine

      subroutine MPI_Reduce_FootPrints_of_Electrons_and_Ions(FP,reduce_root)
         use particle_decomp
         use mpi_defs
         implicit none

         type(FootPrint),intent(inout) :: FP
         integer, intent(in) :: reduce_root
         real(wp) tmp_eScaData(FP%dim_gq,FP%dim_gf,FP%num_eSD), tmp_iScaData(FP%dim_gq,FP%dim_gf,FP%num_iSD)
         integer length,ierr

         ! ============================
         ! Electrons
         ! ============================
         tmp_eScaData = FP%eScaData
         length=FP%dim_gq*FP%dim_gf*FP%num_eSD
         call MPI_Reduce(tmp_eScaData,FP%eScaData,length,mpi_Rsize,MPI_SUM,reduce_root,MPI_COMM_WORLD,ierr)

         ! Poloidal (gq) periodic condition
         FP%eScaData(1,:,:) = FP%eScaData(1,:,:) + FP%eScaData(FP%dim_gq,:,:)
         FP%eScaData(FP%dim_gq,:,:) = FP%eScaData(1,:,:)
         ! Toroidal (gf) periodic condition
         FP%eScaData(:,1,:) = FP%eScaData(:,1,:) + FP%eScaData(:,FP%dim_gf,:)
         FP%eScaData(:,FP%dim_gf,:) = FP%eScaData(:,1,:)


         ! ============================
         ! Ions
         ! ============================
         tmp_iScaData = FP%iScaData
         length=FP%dim_gq*FP%dim_gf*FP%num_iSD
         call MPI_Reduce(tmp_iScaData,FP%iScaData,length,mpi_Rsize,MPI_SUM,reduce_root,MPI_COMM_WORLD,ierr)

         ! Poloidal (gq) periodic condition
         FP%iScaData(1,:,:) = FP%iScaData(1,:,:) + FP%iScaData(FP%dim_gq,:,:)
         FP%iScaData(FP%dim_gq,:,:) = FP%iScaData(1,:,:)
         ! Toroidal (gf) periodic condition
         FP%iScaData(:,1,:) = FP%iScaData(:,1,:) + FP%iScaData(:,FP%dim_gf,:)
         FP%iScaData(:,FP%dim_gf,:) = FP%iScaData(:,1,:)

      end subroutine MPI_Reduce_FootPrints_of_Electrons_and_Ions

      subroutine Cumulate_FP_to_cumFP(FP, cum_FP)

         implicit none
         type(FootPrint),intent(inout) :: FP, cum_FP

         cum_FP%eScaData = cum_FP%eScaData + FP%eScaData
         cum_FP%iScaData = cum_FP%iScaData + FP%iScaData

      end subroutine Cumulate_FP_to_cumFP


      subroutine Cumulate_e_lost_num
         implicit none
         e_cum_lost_num_BDY0 = e_cum_lost_num_BDY0 + e_lost_num_BDY0
         e_cum_lost_num_BDY1 = e_cum_lost_num_BDY1 + e_lost_num_BDY1
      end subroutine Cumulate_e_lost_num

      subroutine Cumulate_i_lost_num
         implicit none
         i_cum_lost_num_BDY0 = i_cum_lost_num_BDY0 + i_lost_num_BDY0
         i_cum_lost_num_BDY1 = i_cum_lost_num_BDY1 + i_lost_num_BDY1
      end subroutine Cumulate_i_lost_num

      subroutine Clear_e_cum_lost_num
         implicit none
         e_cum_lost_num_BDY0 = 0
         e_cum_lost_num_BDY1 = 0
      end subroutine Clear_e_cum_lost_num

      subroutine Clear_i_cum_lost_num
         implicit none
         i_cum_lost_num_BDY0 = 0
         i_cum_lost_num_BDY1 = 0
      end subroutine Clear_i_cum_lost_num


      subroutine Inject_electrons_from_Boundaries(delta_t_seconds)
         use global_parameters
         use mpi_defs
         implicit none

         real(wp), intent(in) :: delta_t_seconds

         integer ierror
         integer global_e_cum_lost_num_BDY0
         integer global_e_cum_lost_num_BDY1
         integer num_of_new_ptls_per_rank
         real(wp) time0, time1

#ifdef _MGY_TIMER_
         time0 = MPI_WTIME()
#endif
         call MPI_ALLREDUCE(e_cum_lost_num_BDY0, global_e_cum_lost_num_BDY0, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)
         call MPI_ALLREDUCE(e_cum_lost_num_BDY1, global_e_cum_lost_num_BDY1, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)

         if(mype==0) then
            write(stdout,'(2(A,I8),A,I8,A)') "istep=",istep, "; Inject electrons, global_e_cum_lost_num (BDY0+BDY1)=(",global_e_cum_lost_num_BDY0," + ", global_e_cum_lost_num_BDY1,")"
            call flush(stdout)
         endif

         ! Create Src from inner boundary (BDY1)
         num_of_new_ptls_per_rank = max(0, global_e_cum_lost_num_BDY1/numberpe)
         call load_electrons_In_or_Out_Flux(num_of_new_ptls_per_rank, delta_t_seconds,0.0,a1)

         ! Create Src from inner boundary (BDY0)
         num_of_new_ptls_per_rank = max(0, global_e_cum_lost_num_BDY0/numberpe)
         call load_electrons_In_or_Out_Flux(num_of_new_ptls_per_rank, delta_t_seconds,1.0,a0)

         call Clear_e_cum_lost_num

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "I'm=",mype, ", Inject_electrons_from_Boundaries takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Inject_electrons_from_Boundaries


      subroutine Inject_ions_from_Boundaries(delta_t_seconds)
         use global_parameters
         use mpi_defs

         implicit none

         real(wp), intent(in) :: delta_t_seconds
         integer ierror
         integer global_i_cum_lost_num_BDY0
         integer global_i_cum_lost_num_BDY1
         integer num_of_new_ptls_per_rank
         real(wp) time0, time1

#ifdef _MGY_TIMER_
         time0 = MPI_WTIME()
#endif
         call MPI_ALLREDUCE(i_cum_lost_num_BDY0, global_i_cum_lost_num_BDY0, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)
         call MPI_ALLREDUCE(i_cum_lost_num_BDY1, global_i_cum_lost_num_BDY1, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)

         if(mype==0) then
            write(stdout,'(2(A,I8),A,I8,A)') "istep=",istep, "; Inject ions, global_i_cum_lost_num (BDY0+BDY1)=(",global_i_cum_lost_num_BDY0," + ", global_i_cum_lost_num_BDY1,")"
            call flush(stdout)
         endif

         ! Create Src from inner boundary (BDY1)
         num_of_new_ptls_per_rank = max(0, global_i_cum_lost_num_BDY1/numberpe)
         call load_ions_In_or_Out_Flux(num_of_new_ptls_per_rank, delta_t_seconds,0.0,a1)

         ! Create Src from inner boundary (BDY0)
         num_of_new_ptls_per_rank = max(0, global_i_cum_lost_num_BDY0/numberpe)
         call load_ions_In_or_Out_Flux(num_of_new_ptls_per_rank, delta_t_seconds,1.0,a0)

         call Clear_i_cum_lost_num

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "I'm=",mype, ", Inject_ions_from_Boundaries takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Inject_ions_from_Boundaries



      subroutine load_electrons_In_or_Out_Flux(num_of_new_ptls, delta_t_seconds, FluxFactor_over_f0, target_a)

         use global_parameters
         use particle_array
         use precision
         use field_array
         use esi_interface
         use mpi_defs
         use function_gg_module
         implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

         integer,intent(in) :: num_of_new_ptls
         real(wp),intent(in) :: delta_t_seconds
         real(wp),intent(in) :: FluxFactor_over_f0
         real(wp),intent(in) :: target_a

         integer :: Nm,m0,m,tmp_num
         real(wp) inv_Nm
         real(wp) integral_f0
         real(wp) v_abs, v_theta, vpara, vperp, mean_vabs ! v_abs is amplitude of total velocity vector, and v_theta is the angle between vpara and vperp
         real(wp) vmax,vc,c1,c2,c3,c4,f2val,c5,c6,ratio, Rand !For Flag_Load_and_WeightDef==2 case.


         integer :: minpts, ifail=1, ierror
         real(wp) :: acc, error_8=1e-5
         real(wp), external:: dsda, D01AHF,jacobian

         integer sid
         real(wp) P_s,ff,P_q,inv_ds
         real(wp), allocatable ::Te_norm(:), ge_val(:)

         real(wp), allocatable :: wpi(:,:),ga(:),wga(:,:),gat(:),dgr(:),da(:),dgq(:),dgf(:),Er_particle(:)
         !$acc declare create(dgr,da,dgq,dgf,wpi,Er_particle)
         real(wp) fip,r,dela
         real(wp) e1,e2,e3,wz0,wz1,wp0,wt00,tw10,wp1,wt01,wt11,wt10
         real(wp) intp_phi, intp_FLA_Te0, intp_phi00
         real(wp), allocatable, dimension(:) :: AdiabTerm ! q*phi/T , q is negative charge here

         integer ip,ij
         integer gcmotion_eon, gcmotion_eon_em, gcmotion_eon_em_MGY, gcmotion_eon_em_MGY_vec_form



         if(num_of_new_ptls<=0) then
            return
         endif

         m0 = me+1 ! starting index of newly loaded particles

         if(target_a>a0+1e-3 .and. target_a<a1-1e-3) then
            me = me+ num_of_new_ptls
         else
            ! From these temporay particlse, only the num_of_new_ptls particles in the simulation domain will be picked
            me = me+3*num_of_new_ptls
         endif

         allocate(ge_val(m0:me))


         do m=m0,me
            zeon(1,m)=target_a
            zeon(7,m)=D01AHF(0.0_wp,zeon(1,m),error_8,minpts,acc,dsda,-1,ifail)
            zeon(7,m)=min(s1,max(s0,zeon(7,m)))

            ! Pick random points (theta,phi) with Helically aligned along B field
            ! zeon(3,m) = 2.0*pi*sprng()
            zeon(3,m) = zetamin + sprng()*(zetamax-zetamin)
            zeon(2,m) = 2.0*pi*sprng()

            elabels(m) = -1.0
         enddo
         !$acc update device(zeon(:,m0:me),elabels(m0:me))

         ! call n_T_u_group_e(m0)   ! get Te, note that back_f0_e cannot be calculated here because it requires zion(4) and zion(6)
         call n_T_u_group_e_partial(m0,me)   ! get Te, note that back_f0_e cannot be calculated here because it requires zion(4) and zion(6)
         !$acc update self(BB(m0:me))

         allocate(Te_norm(m0:me))
         do m=m0,me
            Te_norm(m) = backTe(m)*msm0
         enddo


         ! ! Do Sampling for the following ranges
         ! ! -umax*sqrt(T/m) < u_para < umax*sqrt(T/m)
         ! ! -umax*sqrt(m*T)/(eB) < rho_para < umax*sqrt(m*T)/(eB)
         ! ! 0 < u_perp < umax*sqrt(2T/m)
         ! ! 0 < sqrt(mu) < sqrt(T/B)*u
         ! zeon(4,m) = (2.0*rand_5d(4,m)-1.0)*umax *sqrt(msm0*Te_norm(m))/(csc0*BB(m)) !neglected electron rotation
         ! zeon(6,m) = sqrt(Te_norm(m)/BB(m))*rand_5d(5,m)*umax

         if(Flag_Load_and_WeightDef==0) then

         elseif(Flag_Load_and_WeightDef==1) then
            do m=m0,me
               v_abs = umax*sqrt(Te_norm(m)*m0ms)*sqrt(sprng()) ! v_r = amax*sqrt(Rand)
               v_theta = pi*sprng() ! For half-circle, theta should be 0<theta<2pi

               !Following line is for the half-circle sampling in the polar coordinates (vpara, vperp)
               ! Marker Distribution Funcion for full range: ge = Nm / (Jacobian* (2pi) * (2pi) * (0.5*pi*v_abs_max^2))
               ! NOTE : Since Nm is not determined yet, it will be mulitpled later
               ge_val(m) = 1.0/(Jacobian(zeon(1,m),zeon(2,m))*4.0*pi*pi*0.5*pi*umax*umax*Te_norm(m)*m0ms)

               ! NOTE : Direction of vpara is determined by the sign of wall_norm_vec*B_cov_a
               vpara=v_abs*cos(v_theta)
               vperp=v_abs*sin(v_theta)

               zeon(4,m) = msm0*vpara/(csc0*BB(m))
               zeon(6,m) = vperp*sqrt(msm0/(2.0*BB(m)))
            enddo
         elseif(Flag_Load_and_WeightDef==2) then
            do m=m0,me
               ! First, get v_abs using sprng()
               vmax=umax*sqrt(Te_norm(m)*m0ms)
               vc=2.0*sqrt(Te_norm(m)*m0ms)
               c2=0.01
               f2val=0.2
               c1 = -((-2+c2*vmax+f2val*vmax+ c2*vc-f2val*vc)/( vmax*vc))
               c3 = -((2-2*f2val*vmax-c2*vc+ f2val*vc)/( vmax*(vmax-vc)))
               c4 = -((-2+f2val*vmax+c2*vc)/(vmax-vc))
               Rand = sprng()
               if (Rand < c2*vc+c1*vc**2/2) then
                  v_abs = (-c2+sqrt(c2**2+2*c1*Rand))/c1
               else
                  Rand = Rand -(c2*vc+c1*vc**2/2)
                  v_abs = (-c4 + sqrt(2*c3*Rand+(c4+c3*vc)**2))/c3
               endif
               ! Second, get v_theta using sprng()
               ratio=0.2 ! ratio of (max/min)
               c5=1.0/((1.0+ratio)*0.5*pi)
               c6=(1.0-ratio)/((1.0+ratio)*(0.5*pi)**2)
               Rand = sprng()
               v_theta = 0.5*pi + sign(1.0_wp,Rand-0.5)*((c5-sqrt(c5**2-2*c6*abs(Rand-0.5)))/c6);

               if(v_abs<vc) then
                  ge_val(m) = 1.0*(c1*v_abs+c2)*(c5-c6*abs(v_theta-0.5*pi))/(Jacobian(zeon(1,m),zeon(2,m))*v_abs)/(4.0*pi*pi)
               else
                  ge_val(m) = 1.0*(c3*v_abs+c4)*(c5-c6*abs(v_theta-0.5*pi))/(Jacobian(zeon(1,m),zeon(2,m))*v_abs)/(4.0*pi*pi)
               endif

               ! NOTE : Direction of vpara is determined by the sign of wall_norm_vec*B_cov_a
               vpara=v_abs*cos(v_theta)
               vperp=v_abs*sin(v_theta)

               zeon(4,m) = msm0*vpara/(csc0*BB(m))
               zeon(6,m) = vperp*sqrt(msm0/(2.0*BB(m)))
            enddo
         endif
         !$acc update device(zeon(:,m0:me))

         ! call n_T_u_group_e(m0)   ! This is to get back_f0_e inside. (which requires zeon(4) and zeon(6))
         call n_T_u_group_e_partial(m0,me)   ! This is to get back_f0_e inside. (which requires zeon(4) and zeon(6))

         call electron_locator_partial(m0,me) ! This is needed to calculate evector (gradPhi) and its interpolation (wpi)





         allocate(wpi(3,me),ga(me),wga(3,me),gat(me),dgr(me),da(me),dgq(me),dgf(me),Er_particle(me))
         allocate(AdiabTerm(me))



         ! Calculate evector (wpi) first
         !$acc parallel loop present(wzeon,jteon0,wpeon,wteon0,jteon1,wteon1,wpi)
         do m=1,me
            wz1=wzeon(m)
            wz0=1.0_wp-wz1

            ij=jteon0(m)
            wp0=1.0_wp-wpeon(m)
            wt00=1.0_wp-wteon0(m)
            e1=wp0*wt00*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
            e2=wp0*wt00*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
            e3=wp0*wt00*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
            intp_phi = wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))
            intp_FLA_Te0 = wp0*wt00*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            ij=ij+1
            wt10=1.0_wp-wt00
            e1=e1+wp0*wt10*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
            e2=e2+wp0*wt10*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
            e3=e3+wp0*wt10*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
            intp_phi=intp_phi+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))
            intp_FLA_Te0=intp_FLA_Te0+wp0*wt10*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            ij=jteon1(m)
            wp1=1.0_wp-wp0
            wt01=1.0_wp-wteon1(m)
            e1=e1+wp1*wt01*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
            e2=e2+wp1*wt01*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
            e3=e3+wp1*wt01*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
            intp_phi=intp_phi+wp1*wt01*(wz0*phi(0,ij)+wz1*phi(1,ij))
            intp_FLA_Te0=intp_FLA_Te0+wp1*wt01*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            ij=ij+1
            wt11=1.0_wp-wt01
            e1=e1+wp1*wt11*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
            e2=e2+wp1*wt11*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
            e3=e3+wp1*wt11*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
            intp_phi=intp_phi+wp1*wt11*(wz0*phi(0,ij)+wz1*phi(1,ij)) ! --> phi at electrons
            intp_FLA_Te0=intp_FLA_Te0+wp1*wt11*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            wpi(1,m)=e1
            wpi(2,m)=e2
            wpi(3,m)=e3

            r=(zeon(7,m)-s0)/deltar       ! in radial coordinator s
            ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
            fip=r-real(ip,wp)
            intp_phi00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons
            if ( Flag_Def_dne_adiab==0 ) then
               AdiabTerm(m)=csc0*(intp_phi-real(1.0-ZF_SEPARATION,wp)*intp_phi00)/(backTe(m)*msm0) ! Adiabatic term (q*dphi_NZ/Te)
            elseif(Flag_Def_dne_adiab==1) then
               ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) for the New Poisson3D solver
               AdiabTerm(m)=csc0*(intp_phi+real(ZF_SEPARATION,wp)*intp_phi00)/(backTe(m)*msm0) ! Adiabatic term (q*dphi_tot/Te)
            elseif(Flag_Def_dne_adiab==2) then
               ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) with FLA_Te0
               AdiabTerm(m)=csc0*(intp_phi+real(ZF_SEPARATION,wp)*intp_phi00)/(intp_FLA_Te0) ! Adiabatic term (q*dphi_tot/Te)
            else
               AdiabTerm(m)=-1e100
               write(stdout,*) "Error!! Flag_Def_dne_adiab=",Flag_Def_dne_adiab, " is not defined in pushe_new.F90!!"
               STOP
            endif
         enddo

         ! Get Er_particle
         Er401_bp=Er401_bp0
         Er_grid=Er_grid0
         dela=1.0_wp/Ydr
         !$acc parallel loop present(zeon,zeon0,Er_particle)
         do m=1,me
            r=zeon(1,m)-a0
            r=r*dela
            ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
            fip=r-real(ip,wp)
            Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
         enddo

         if(electromagnetic == 1)then
            call Get_MP_alpha_amplitude_arr_at_zeon(istep*tstep, ga(1:me), wga(1:3,1:me), gat(1:me))

            !m=gcmotion_eon_em(dgr,da,dgq,dgf,zeon,wpi,ga,wga,gat,Er_particle,nonlinear,me,zeondim)
            if(Flag_GC_Solver==0) then
               m=gcmotion_eon_em_MGY(dgr,da,dgq,dgf,zeon,wpi,ga,wga,gat,Er_particle,nonlinear,me,zeondim)
            elseif(Flag_GC_Solver==1) then
               m=gcmotion_eon_em_MGY_vec_form(dgr,da,dgq,dgf,zeon,wpi,ga,wga,gat,Er_particle,nonlinear,me,zeondim)
            else
               write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
               call flush(stdout)
               stop
            endif
         else
            m=gcmotion_eon(dgr,da,dgq,dgf,zeon,wpi,Er_particle,nonlinear,me,zeondim)
         endif

         !$acc update self(da(m0:me))

         if(mype==0) then
            write(stdout,'(3(A,ES14.6))') "uveloci=", uveloci,"; utime=", utime, "; multiplied=",utime*uveloci
            write(stdout,'(3(A,ES14.6))') "vpara(m0)/uveloci=", m0ms*csc0*BB(m0)*zeon(4,m0),"; delta_t_seconds=", delta_t_seconds, "; multiplied=",m0ms*csc0*BB(m0)*zeon(4,m0)*uveloci*delta_t_seconds
            write(stdout,'(3(A,ES14.6))') "da(m0)/uveloci=", da(m0),"; delta_t_seconds=", delta_t_seconds, "; multiplied=",da(m0)*uveloci*delta_t_seconds
            call flush(stdout)
         endif

         ! Set initial weights of electrons
         ! zeon(5,m) is "df" or "dh"
         ! zeon0(6,m) is "f_tot"

         do m=m0,me
            vpara = m0ms*csc0*BB(m)*zeon(4,m)  ! vpara = (eB/m)*rho_para
            mean_vabs = sqrt(8.0*Te_norm(m)/(pi*msm0))
            ! integral_f0 = 0.25*mean_vabs*delta_t_seconds*vpara*back_f0_e(m)/ge_val(m)
            ! integral_f0 = vpara*back_f0_e(m)/ge_val(m)
            ! integral_f0 = vpara*delta_t_seconds*back_f0_e(m)/ge_val(m)
            ! integral_f0 = abs(vpara)*uveloci*1e-2*delta_t_seconds*back_f0_e(m)/ge_val(m) ! uveloci*1e-2 = unit velocity in (m/s)
            integral_f0 = abs(da(m))*uveloci*1e-2*delta_t_seconds*back_f0_e(m)/ge_val(m) ! uveloci*1e-2 = unit velocity in (m/s)
            zeon0(6,m) = FluxFactor_over_f0*integral_f0
            zeon(5,m)  = zeon0(6,m)-integral_f0 +int(ifelectron/2)*AdiabTerm(m)*integral_f0 ! For dh, dh=f_tot-f0+q*phi/Te0*f0

            if(da(m)<0) then
               zeon(1,m)=target_a-1e-5
               zeon(7,m)=D01AHF(0.0_wp,zeon(1,m),error_8,minpts,acc,dsda,-1,ifail)
               zeon(7,m)=min(s1,max(s0,zeon(7,m)))
            else
               zeon(1,m)=target_a+1e-5
               zeon(7,m)=D01AHF(0.0_wp,zeon(1,m),error_8,minpts,acc,dsda,-1,ifail)
               zeon(7,m)=min(s1,max(s0,zeon(7,m)))
            endif
         enddo


         if(target_a<a0+1e-3 .or. target_a>a1-1e-3) then ! If target_a is near the boudnaries, check the particles
            tmp_num=0
            do m=m0,me
               if(zeon(1,m)>a0 .and. zeon(1,m)<a1) then
                  tmp_num = tmp_num+1
               endif
               if(tmp_num>=num_of_new_ptls) then
                  me=m
                  exit
               endif
            enddo
         endif

         tmp_num=me-m0+1
         ! Nm is total number of new particles regarding whole cpu ranks (MPI_COMM_WORLD)
         call MPI_ALLREDUCE(tmp_num,Nm,1,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,ierror)
         inv_Nm=1.0/Nm
         !Calculate the correct weight values by considering the number of total marker particles
         zeon0(6,m0:me) = inv_Nm*zeon0(6,m0:me)
         zeon(5,m0:me)  = inv_Nm*zeon(5,m0:me)

         we_init(m0:me)  = zeon(5,m0:me)
         we_vpara(m0:me) = 0
         we_vperp(m0:me) = 0
         we_accel(m0:me) = 0
         we_adiab(m0:me) = 0
         we_fdamp(m0:me) = 0

         zeon0(1:5,m0:me)     = zeon(1:5,m0:me)
         zeon0(7,m0:me)       = zeon(7,m0:me)
         UnNorm_theta_e(m0:me)= zeon(2,m0:me)
         UnNorm_zeta_e(m0:me) = zeon(3,m0:me)


         !$acc update device(we_init(m0:me),we_vpara(m0:me),we_vperp(m0:me),we_accel(m0:me),we_adiab(m0:me),we_fdamp(m0:me))
         !$acc update device(zeon(:,m0:me),zeon0(:,m0:me),elabels(m0:me))
         !$acc update device(UnNorm_theta_e(m0:me),UnNorm_zeta_e(m0:me))

         ! Distribute electrons to appropriate MPI ranks
         call treat_electrons_out_of_boundary(m0,me)
         ! call electron_locator_partial(m0,me)
         ! call n_T_u_group_e_partial(m0,me)


         write(stdout,'(6(A,I8))') "istep=",istep, "; mype=",mype,"; m0=",m0,"; me=",me,"; me-m0=",me-m0,"; memax=",memax
         call flush(stdout)
         call MPI_Barrier(MPI_COMM_WORLD,ierror)

         deallocate(Te_norm)
         deallocate(ge_val)
         deallocate(wpi,ga,wga,gat,dgr,da,dgq,dgf,Er_particle)
         deallocate(AdiabTerm)

      end subroutine load_electrons_In_or_Out_Flux

      subroutine load_ions_In_or_Out_Flux(num_of_new_ptls, delta_t_seconds, FluxFactor_over_f0, target_a)

         use global_parameters
         use particle_array
         use precision
         use field_array
         use esi_interface
         use mpi_defs
         use function_gg_module
         implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

         integer,intent(in) :: num_of_new_ptls
         real(wp),intent(in) :: delta_t_seconds
         real(wp),intent(in) :: FluxFactor_over_f0
         real(wp),intent(in) :: target_a

         integer :: Nm,m0,m,tmp_num
         real(wp) inv_Nm
         real(wp) integral_f0
         real(wp) v_abs, v_theta, vpara, vperp, mean_vabs ! v_abs is amplitude of total velocity vector, and v_theta is the angle between vpara and vperp
         real(wp) vmax,vc,c1,c2,c3,c4,f2val,c5,c6,ratio, Rand !For Flag_Load_and_WeightDef==2 case.


         integer :: minpts, ifail=1, ierror
         real(wp) :: acc, error_8=1e-5
         real(wp), external:: dsda, D01AHF,jacobian

         integer sid
         real(wp) P_s,ff,P_q,inv_ds
         real(wp), allocatable :: gi_val(:)

         real(wp), allocatable :: wpi(:,:),ga(:),wga(:,:),gat(:),dgr(:),da(:),dgq(:),dgf(:),Er_particle(:)
         !$acc declare create(dgr,da,dgq,dgf,wpi,Er_particle)
         real(wp) fip,r,dela
         real(wp) e1,e2,e3,wz0,wz1,wp0,wt00,tw10,wp1,wt01,wt11,wt10
         integer ip,ij,larmor
         integer gcmotion, gcmotion_em, gcmotion_em_MGY, gcmotion_em_MGY_vec_form


         if(num_of_new_ptls<=0) then
            return
         endif

         m0 = mi+1 ! starting index of newly loaded particles

         if(target_a>a0+1e-3 .and. target_a<a1-1e-3) then
            mi = mi+ num_of_new_ptls
         else
            ! From these temporay particlse, only the num_of_new_ptls particles in the simulation domain will be picked
            mi = mi+3*num_of_new_ptls
         endif

         allocate(gi_val(m0:mi))


         do m=m0,mi
            zion(1,m)=target_a
            zion(7,m)=D01AHF(0.0_wp,zion(1,m),error_8,minpts,acc,dsda,-1,ifail)
            zion(7,m)=min(s1,max(s0,zion(7,m)))

            ! Pick random points (theta,phi) with Helically aligned along B field
            ! zion(3,m) = 2.0*pi*sprng()
            zion(3,m) = zetamin + sprng()*(zetamax-zetamin)
            zion(2,m) = 2.0*pi*sprng()

            ilabels(m) = -1.0
         enddo
         !$acc update device(zion(:,m0:mi),ilabels(m0:mi))

         ! call n_T_u_group(m0)   ! get Te, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
         call n_T_u_group_partial(m0,mi)   ! get Te, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
         !$acc update self(BB(m0:mi),backu(m0:mi),F(m0:mi))
         backu(m0:mi)=backu(m0:mi)*F(m0:mi)/BB(m0:mi) ! =u_//=I*omega_t/B
         !$acc update device(backu(m0:mi))



         ! ! Do Sampling for the following ranges
         ! ! -umax*sqrt(T/m) < u_para < umax*sqrt(T/m)
         ! ! -umax*sqrt(m*T)/(eB) < rho_para < umax*sqrt(m*T)/(eB)
         ! ! 0 < u_perp < umax*sqrt(2T/m)
         ! ! 0 < sqrt(mu) < sqrt(T/B)*u
         ! zion(4,m) = (2.0*rand_5d(4,m)-1.0)*umax *sqrt(msm0*Te_norm(m))/(csc0*BB(m)) !neglected electron rotation
         ! zion(6,m) = sqrt(backT(m)/BB(m))*rand_5d(5,m)*umax

         if(Flag_Load_and_WeightDef==0) then

         elseif(Flag_Load_and_WeightDef==1) then
            do m=m0,mi
               v_abs = umax*sqrt(backT(m))*sqrt(sprng()) ! v_r = amax*sqrt(Rand)
               v_theta = pi*sprng() ! For half-circle, theta should be 0<theta<2pi

               !Following line is for the half-circle sampling in the polar coordinates (vpara, vperp)
               ! Marker Distribution Funcion for full range: ge = Nm / (Jacobian* (2pi) * (2pi) * (0.5*pi*v_abs_max^2))
               ! NOTE : Since Nm is not determined yet, it will be mulitpled later
               gi_val(m) = 1.0/(Jacobian(zion(1,m),zion(2,m))*4.0*pi*pi*0.5*pi*umax*umax*backT(m))

               ! NOTE : Direction of vpara is determined by the sign of wall_norm_vec*B_cov_a
               vpara=v_abs*cos(v_theta)
               vperp=v_abs*sin(v_theta)

               zion(4,m) = vpara/(BB(m))
               zion(4,m) = zion(4,m) + backu(m)/BB(m)
               zion(6,m) = vperp*sqrt(1.0/(2.0*BB(m)))
            enddo
         elseif(Flag_Load_and_WeightDef==2) then
            do m=m0,mi
               ! First, get v_abs using sprng()
               vmax=umax*sqrt(backT(m))
               vc=2.0*sqrt(backT(m))
               c2=0.01
               f2val=0.2
               c1 = -((-2+c2*vmax+f2val*vmax+ c2*vc-f2val*vc)/( vmax*vc))
               c3 = -((2-2*f2val*vmax-c2*vc+ f2val*vc)/( vmax*(vmax-vc)))
               c4 = -((-2+f2val*vmax+c2*vc)/(vmax-vc))
               Rand = sprng()
               if (Rand < c2*vc+c1*vc**2/2) then
                  v_abs = (-c2+sqrt(c2**2+2*c1*Rand))/c1
               else
                  Rand = Rand -(c2*vc+c1*vc**2/2)
                  v_abs = (-c4 + sqrt(2*c3*Rand+(c4+c3*vc)**2))/c3
               endif
               ! Second, get v_theta using sprng()
               ratio=0.2 ! ratio of (max/min)
               c5=1.0/((1.0+ratio)*0.5*pi)
               c6=(1.0-ratio)/((1.0+ratio)*(0.5*pi)**2)
               Rand = sprng()
               v_theta = 0.5*pi + sign(1.0_wp,Rand-0.5)*((c5-sqrt(c5**2-2*c6*abs(Rand-0.5)))/c6);

               if(v_abs<vc) then
                  gi_val(m) = 1.0*(c1*v_abs+c2)*(c5-c6*abs(v_theta-0.5*pi))/(Jacobian(zion(1,m),zion(2,m))*v_abs)/(4.0*pi*pi)
               else
                  gi_val(m) = 1.0*(c3*v_abs+c4)*(c5-c6*abs(v_theta-0.5*pi))/(Jacobian(zion(1,m),zion(2,m))*v_abs)/(4.0*pi*pi)
               endif

               ! NOTE : Direction of vpara is determined by the sign of wall_norm_vec*B_cov_a
               vpara=v_abs*cos(v_theta)
               vperp=v_abs*sin(v_theta)

               zion(4,m) = vpara/(BB(m))
               zion(4,m) = zion(4,m) + backu(m)/BB(m)
               zion(6,m) = vperp*sqrt(1.0/(2.0*BB(m)))
            enddo
         endif
         !$acc update device(zion(:,m0:mi))

         ! call n_T_u_group(m0)   ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))
         call n_T_u_group_partial(m0,mi)   ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))

         call ion_locator_at_larmor_points_partial(m0,mi) ! This is needed to calculate evector (gradPhi) and its interpolation (wpi)



         allocate(wpi(3,mi),ga(mi),wga(3,mi),gat(mi),dgr(mi),da(mi),dgq(mi),dgf(mi),Er_particle(mi))



         ! Calculate evector (wpi) first
         !$acc parallel loop present(zion,zion0)
         do m=1,mi
            e1=0.0_wp
            e2=0.0_wp
            e3=0.0_wp
            wz1=wzion(m)
            wz0=1.0_wp-wz1

            !$acc loop seq
            do larmor=1, larmor_num
               ij=jtion0(larmor,m)
               wp0=1.0_wp-wpion(larmor,m)
               wt00=1.0_wp-wtion0(larmor,m)
               e1=e1+wp0*wt00*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
               e2=e2+wp0*wt00*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
               e3=e3+wp0*wt00*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

               ij=ij+1
               wt10=1.0_wp-wt00
               e1=e1+wp0*wt10*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
               e2=e2+wp0*wt10*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
               e3=e3+wp0*wt10*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

               ij=jtion1(larmor,m)
               wp1=1.0_wp-wp0
               wt01=1.0_wp-wtion1(larmor,m)
               e1=e1+wp1*wt01*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
               e2=e2+wp1*wt01*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
               e3=e3+wp1*wt01*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

               ij=ij+1
               wt11=1.0_wp-wt01
               e1=e1+wp1*wt11*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
               e2=e2+wp1*wt11*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
               e3=e3+wp1*wt11*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
            enddo

            wpi(1,m)=e1*inv_larmor_num
            wpi(2,m)=e2*inv_larmor_num
            wpi(3,m)=e3*inv_larmor_num
         enddo

         ! Get Er_particle
         Er401_bp=Er401_bp0
         Er_grid=Er_grid0
         dela=1.0_wp/Ydr
         !$acc parallel loop present(zion,zion0,Er_particle)
         do m=1,mi
            r=zion(1,m)-a0
            r=r*dela
            ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
            fip=r-real(ip,wp)
            Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
         enddo

         if(electromagnetic == 1)then
            call Get_MP_alpha_amplitude_arr_at_zion(istep*tstep, ga(1:mi), wga(1:3,1:mi), gat(1:mi))

            !m=gcmotion_em(dgr,da,dgq,dgf,zion,wpi,ga,wga,gat,Er_particle,nonlinear,mi,ziondim)
            if(Flag_GC_Solver==0) then
               m=gcmotion_em_MGY(dgr,da,dgq,dgf,zion,wpi,ga,wga,gat,Er_particle,nonlinear,mi,ziondim)
            elseif(Flag_GC_Solver==1) then
               m=gcmotion_em_MGY_vec_form(dgr,da,dgq,dgf,zion,wpi,ga,wga,gat,Er_particle,nonlinear,mi,ziondim)
            else
               write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
               call flush(stdout)
               stop
            endif
         else
            m=gcmotion(dgr,da,dgq,dgf,zion,wpi,Er_particle,nonlinear,mi,ziondim)
         endif

         !$acc update self(da(m0:mi))

         if(mype==0) then
            write(stdout,'(3(A,ES14.6))') "uveloci=", uveloci,"; utime=", utime, "; multiplied=",utime*uveloci
            write(stdout,'(3(A,ES14.6))') "vpara(m0)/uveloci=", BB(m0)*zion(4,m0),"; delta_t_seconds=", delta_t_seconds, "; multiplied=",m0ms*csc0*BB(m0)*zion(4,m0)*uveloci*delta_t_seconds
            write(stdout,'(3(A,ES14.6))') "da(m0)/uveloci=", da(m0),"; delta_t_seconds=", delta_t_seconds, "; multiplied=",da(m0)*uveloci*delta_t_seconds
            call flush(stdout)
         endif

         ! Set initial weights of electrons
         ! zion(5,m) is "df" or "dh"
         ! zion0(6,m) is "f_tot"

         do m=m0,mi
            vpara = BB(m)*zion(4,m)  ! vpara = (eB/m)*rho_para
            integral_f0 = abs(da(m))*uveloci*1e-2*delta_t_seconds*back_f0_i(m)/gi_val(m) ! uveloci*1e-2 = unit velocity in (m/s)
            zion0(6,m) = FluxFactor_over_f0*integral_f0
            zion(5,m)  = zion0(6,m)-integral_f0

            if(da(m)<0) then
               zion(1,m)=target_a-1e-5
               zion(7,m)=D01AHF(0.0_wp,zion(1,m),error_8,minpts,acc,dsda,-1,ifail)
               zion(7,m)=min(s1,max(s0,zion(7,m)))
            else
               zion(1,m)=target_a+1e-5
               zion(7,m)=D01AHF(0.0_wp,zion(1,m),error_8,minpts,acc,dsda,-1,ifail)
               zion(7,m)=min(s1,max(s0,zion(7,m)))
            endif
         enddo


         if(target_a<a0+1e-3 .or. target_a>a1-1e-3) then
            tmp_num=0
            do m=m0,mi
               if(zion(1,m)>a0 .and. zion(1,m)<a1) then
                  tmp_num = tmp_num+1
               endif
               if(tmp_num>=num_of_new_ptls) then
                  mi=m
                  exit
               endif
            enddo
         endif

         tmp_num=mi-m0+1
         ! Nm is total number of new particles regarding whole cpu ranks (MPI_COMM_WORLD)
         call MPI_ALLREDUCE(tmp_num,Nm,1,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,ierror)
         inv_Nm=1.0/Nm
         !Calculate the correct weight values by considering the number of total marker particles
         zion0(6,m0:mi) = inv_Nm*zion0(6,m0:mi)
         zion(5,m0:mi)  = inv_Nm*zion(5,m0:mi)


         wi_init(m0:mi)  = zion(5,m0:mi)
         wi_vpara(m0:mi) = 0
         wi_vperp(m0:mi) = 0
         wi_accel(m0:mi) = 0
         wi_fdamp(m0:mi) = 0

         zion0(1:5,m0:mi)     = zion(1:5,m0:mi)
         zion0(7,m0:mi)       = zion(7,m0:mi)
         UnNorm_theta_i(m0:mi)= zion(2,m0:mi)
         UnNorm_zeta_i(m0:mi) = zion(3,m0:mi)


         !$acc update device(wi_init(m0:mi),wi_vpara(m0:mi),wi_vperp(m0:mi),wi_accel(m0:mi),wi_fdamp(m0:mi))
         !$acc update device(zion(:,m0:mi),zion0(:,m0:mi),ilabels(m0:mi))
         !$acc update device(UnNorm_theta_i(m0:mi),UnNorm_zeta_i(m0:mi))

         ! Distribute electrons to appropriate MPI ranks
         call treat_ions_out_of_boundary(m0,mi)
         ! call ion_locator_at_larmor_points_partial(m0,mi)
         ! call n_T_u_group_partial(m0,mi)


         write(stdout,'(6(A,I8))') "istep=",istep, "; mype=",mype,"; m0=",m0,"; mi=",mi,"; mi-m0=",mi-m0,"; mimax=",mimax
         call flush(stdout)
         call MPI_Barrier(MPI_COMM_WORLD,ierror)

         deallocate(gi_val)
         deallocate(wpi,ga,wga,gat,dgr,da,dgq,dgf,Er_particle)

      end subroutine load_ions_In_or_Out_Flux


      subroutine treat_electrons_out_of_boundary(m_start, m_end)

         use global_parameters
         use particle_array
         use field_array
         use function_gg_module
         use Ptl_locators_partial
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0,m1
         integer m,i,j
         real(wp) psimax,psimin
         integer lost_mid_arr(memax) ! index of loss particles
         integer lost_num ! number of loss particles across boundary
         integer, allocatable :: fill_mid_arr(:)
         integer fill_num
         integer mtop, lasth
         integer dim_zeon, dim_zeon0, dim_RK_cum_zeon
         real(wp) zeon_a_arr(me)

         m0=1
         m1=me
         if(present(m_start)) m0=m_start
         if(present(m_end))   m1=m_end


         if(acc_newer_H_zeon) then
            !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_H_zeon=.false.
            acc_newer_D_zeon=.false.
         endif

         dim_zeon = size(zeon,1)
         dim_zeon0 = size(zeon0,1)
         dim_RK_cum_zeon = size(RK_cum_zeon,1)

         psimax=psi_grid(mpsi)
         psimin=psi_grid(0)

         if(Flag_Ptl_BC_type==0) then
            !$acc parallel loop present(zeon,zeon0)
            do m=m0,m1
               if(zeon(1,m) > psimax)then
                  zeon(1,m)=psimax
                  zeon(2,m)=2.0_wp*pi-zeon0(2,m)
                  zeon(3,m)=zeon0(3,m)
                  zeon(4,m)=zeon0(4,m)
                  zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
                  !   zeon(5,m)=-zeon0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
                  !   zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
                  zeon(7,m)=zeon0(7,m)
               elseif(zeon(1,m) < psimin)then
                  zeon(1,m)=psimin
                  zeon(2,m)=2.0_wp*pi-zeon0(2,m)
                  zeon(3,m)=zeon0(3,m)
                  zeon(4,m)=zeon0(4,m)
                  zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
                  !   zeon(5,m)=-zeon0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
                  !   zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
                  zeon(7,m)=zeon0(7,m)
               endif
            enddo
            acc_newer_D_zeon=.true.

            return

         elseif(Flag_Ptl_BC_type==1) then ! delete out-of_boundary particles from the list

            !$acc parallel loop
            do m=m0,m1
               zeon_a_arr(m)=zeon(1,m)
            enddo

            lost_num=0
            e_lost_num_BDY0=0
            e_lost_num_BDY1=0
            ! Find particle indices exiting the boundary
            do m=m0,m1
               if(zeon_a_arr(m) > psimax)then
                  lost_num = lost_num +1
                  lost_mid_arr(lost_num) = m
                  e_lost_num_BDY1 = e_lost_num_BDY1 +1
               elseif(zeon_a_arr(m) < psimin)then
                  lost_num = lost_num +1
                  lost_mid_arr(lost_num) = m
                  e_lost_num_BDY0 = e_lost_num_BDY0 +1
               endif
            enddo


            write(stdout,'(5(A,I8))') "istep=",istep,"; mype=",mype,"; e_lost_num=",lost_num, "=",e_lost_num_BDY0, "+",e_lost_num_BDY1
            ! write(stdout,*) "istep=",istep,"; mype=",mype,"electron lost_mid_arr=",lost_mid_arr(1:lost_num)
            call flush(stdout)


            ! NOTE : mtop always should be mi even if m1!=mi.
            ! fill holes
            mtop = me
            me = me - lost_num ! remained me after loss
            lasth = lost_num

            allocate(fill_mid_arr(lost_num))
            fill_num=0
            do i=1, lost_num
               m = lost_mid_arr(i)
               if(m>me) exit
               do while(mtop == lost_mid_arr(lasth))
                  mtop=mtop-1
                  lasth=lasth-1
               enddo
               fill_mid_arr(i) = mtop
               mtop = mtop -1
               fill_num=i
            enddo

            ! write(stdout,*) "mype=",mype,"fill_num=",fill_num, "fill_mid_arr=",fill_mid_arr(1:fill_num)
            ! call flush(stdout)

            !$acc parallel loop
            do i=1, fill_num
               m = lost_mid_arr(i)
               mtop= fill_mid_arr(i)
               !$acc loop seq
               do j=1,dim_zeon
                  zeon(j, m)   = zeon(j, mtop)
               enddo
               !$acc loop seq
               do j=1,dim_zeon0
                  zeon0(j, m)  = zeon0(j, mtop)
               enddo
               elabels(m)      = elabels(mtop)
               we_init(m)      = we_init(mtop)
               we_vpara(m)     = we_vpara(mtop)
               we_vperp(m)     = we_vperp(mtop)
               we_accel(m)     = we_accel(mtop)
               we_adiab(m)     = we_adiab(mtop)
               we_fdamp(m)     = we_fdamp(mtop)

               UnNorm_theta_e(m) =UnNorm_theta_e(mtop)
               UnNorm_zeta_e(m) =UnNorm_zeta_e(mtop)
               if(Flag_RK4==1) then
                  !$acc loop seq
                  do j=1,dim_RK_cum_zeon
                     RK_cum_zeon(j,m)= RK_cum_zeon(j,mtop)
                  enddo
               endif
            enddo
            deallocate(fill_mid_arr)
            acc_newer_D_zeon=.true.

            call electron_locator_partial(m_start=m0, m_end=min(me,lost_mid_arr(lasth)))
            call n_T_u_group_e_partial(m_start=m0, m_end=min(me,lost_mid_arr(lasth)))
         endif

      end subroutine treat_electrons_out_of_boundary


      subroutine treat_ions_out_of_boundary(m_start, m_end)

         use global_parameters
         use particle_array
         use field_array
         use function_gg_module
         use Ptl_locators_partial
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0, m1

         integer m,i,j
         real(wp) psimax,psimin
         integer lost_mid_arr(mimax) ! index of loss particles
         integer lost_num ! number of loss particles across boundary
         integer, allocatable :: fill_mid_arr(:)
         integer fill_num
         integer mtop, lasth
         integer dim_zion, dim_zion0, dim_RK_cum_zion
         real(wp) zion_a_arr(mi)

         m0=1
         m1=mi
         if(present(m_start)) m0=m_start
         if(present(m_end))   m1=m_end

         if(acc_newer_H_zion) then
            !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
            acc_newer_H_zion=.false.
            acc_newer_D_zion=.false.
         endif

         dim_zion = size(zion,1)
         dim_zion0 = size(zion0,1)
         dim_RK_cum_zion = size(RK_cum_zion,1)

         psimax=psi_grid(mpsi)
         psimin=psi_grid(0)

         if(Flag_Ptl_BC_type==0) then
            !$acc parallel loop present(zion,zion0)
            do m=m0,m1
               if(zion(1,m) > psimax)then
                  zion(1,m)=psimax
                  zion(2,m)=2.0_wp*pi-zion0(2,m)
                  zion(3,m)=zion0(3,m)
                  zion(4,m)=zion0(4,m)
                  zion(5,m)=zion0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
                  !   zion(5,m)=-zion0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
                  !   zion(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
                  zion(7,m)=zion0(7,m)
               elseif(zion(1,m) < psimin)then
                  zion(1,m)=psimin
                  zion(2,m)=2.0_wp*pi-zion0(2,m)
                  zion(3,m)=zion0(3,m)
                  zion(4,m)=zion0(4,m)
                  zion(5,m)=zion0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
                  !   zion(5,m)=-zion0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
                  !   zion(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
                  zion(7,m)=zion0(7,m)
               endif
            enddo
            acc_newer_D_zion=.true.

            return

         elseif(Flag_Ptl_BC_type==1) then ! delete out-of_boundary particles from the list

            !$acc parallel loop
            do m=m0,m1
               zion_a_arr(m)=zion(1,m)
            enddo

            lost_num=0
            i_lost_num_BDY0=0
            i_lost_num_BDY1=0
            ! Find particle indices exiting the boundary
            do m=m0,m1
               if(zion_a_arr(m) > psimax)then
                  lost_num = lost_num +1
                  lost_mid_arr(lost_num) = m
                  i_lost_num_BDY1 = i_lost_num_BDY1 +1
               elseif(zion_a_arr(m) < psimin)then
                  lost_num = lost_num +1
                  lost_mid_arr(lost_num) = m
                  i_lost_num_BDY0 = i_lost_num_BDY0 +1
               endif
            enddo


            write(stdout,'(5(A,I8))') "istep=",istep,"; mype=",mype,"; i_lost_num=",lost_num, "=",i_lost_num_BDY0, "+",i_lost_num_BDY1
            ! write(stdout,*) "istep=",istep,"; mype=",mype,"electron lost_mid_arr=",lost_mid_arr(1:lost_num)
            call flush(stdout)


            ! NOTE : mtop always should be mi even if m1!=mi.
            ! fill holes
            mtop = mi
            mi = mi - lost_num ! remained mi after loss
            lasth = lost_num

            allocate(fill_mid_arr(lost_num))
            fill_num=0
            do i=1, lost_num
               m = lost_mid_arr(i)
               if(m>mi) exit
               do while(mtop == lost_mid_arr(lasth))
                  mtop=mtop-1
                  lasth=lasth-1
               enddo
               fill_mid_arr(i) = mtop
               mtop = mtop -1
               fill_num=i
            enddo

            ! write(stdout,*) "mype=",mype,"fill_num=",fill_num, "fill_mid_arr=",fill_mid_arr(1:fill_num)
            ! call flush(stdout)

            !$acc parallel loop
            do i=1, fill_num
               m = lost_mid_arr(i)
               mtop= fill_mid_arr(i)
               !$acc loop seq
               do j=1,dim_zion
                  zion(j, m)   = zion(j, mtop)
               enddo
               !$acc loop seq
               do j=1,dim_zion0
                  zion0(j, m)  = zion0(j, mtop)
               enddo
               ilabels(m)      = ilabels(mtop)
               wi_init(m)      = wi_init(mtop)
               wi_vpara(m)     = wi_vpara(mtop)
               wi_vperp(m)     = wi_vperp(mtop)
               wi_accel(m)     = wi_accel(mtop)
               wi_fdamp(m)     = wi_fdamp(mtop)

               UnNorm_theta_i(m) =UnNorm_theta_i(mtop)
               UnNorm_zeta_i(m) =UnNorm_zeta_i(mtop)
               if(Flag_RK4==1) then
                  !$acc loop seq
                  do j=1,dim_RK_cum_zion
                     RK_cum_zion(j,m)= RK_cum_zion(j,mtop)
                  enddo
               endif
            enddo
            deallocate(fill_mid_arr)
            acc_newer_D_zion=.true.

            call ion_locator_at_larmor_points_partial(m_start=m0, m_end=min(mi,lost_mid_arr(lasth)))
            call n_T_u_group_partial(m_start=m0, m_end=min(mi,lost_mid_arr(lasth)))
         endif

      end subroutine treat_ions_out_of_boundary




end module Ptl_Src_Loss
