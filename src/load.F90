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
   use function_gg_module
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
   call n_T_u_group                       ! get T, n and omega_t
   !$acc update self(backu(:mi),backT(:mi))
   call esigetbf(zion,mi,ziondim)   ! get B-->BB, I-->F  12/07
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
      ! for arbitrary loading according to a distribution function g_0(x,v), then zion0(8,*)=f_0/g_0
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
         zion0(8,m)=exp(backu(m)**2/(2.0*backT(m)))/ave0     ! =1 if u_//=0 (no rotation)
      end do
      acc_newer_H_zion=.true.

      !!! the following print may hang forever when using large number of processors Jan, 2011
      ! print*, 'f_0/g_0=', zion0(8,mi/2), mype

      ! realistic nonuniform loading
      !     zion0(8,1:mi)=1.0

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
         !        zeon0(8,m)=zion0(8,m)
         zeon0(8,m)=1.0_wp
      enddo
      acc_newer_H_zeon=.true.

      do m=1,me
         we_init(m)  = zeon(5,m)
         we_vpara(m) = 0
         we_vperp(m) = 0
         we_accel(m) = 0
         we_adiab(m) = 0
      enddo



      call n_T_u_group_e               ! get Te
      !$acc update self(backTe(:me))
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
         if(ifdphi==2) zeon(9,ntracer_e)=0.0_wp

         !!!!!!!!!TEST electron tracer using me/mi=1!!!!!!!!!!!!!!!!!!
         !        zeon(:,ntracer_e)=zion(:,ntracer)

         zeon0(1:5,ntracer_e)=zeon(1:5,ntracer_e)
         zeon0(7,ntracer_e)=zeon(7,ntracer_e)
         !        zeon0(8,ntracer_e)=zion0(8,ntracer)
         zeon0(8,ntracer_e)=1.0_wp
         if(ifdphi==2) zeon0(9,ntracer_e)=0.0_wp

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

   if(acc_newer_H_zeon) then
      !$acc update device(we_init(:),we_vpara(:),we_vperp(:),we_accel(:),we_adiab(:))
      !$acc update device(zeon(:,:),zeon0(:,:),elabels(:))
      !$acc update device(UnNorm_theta_e(:),UnNorm_zeta_e(:))
      acc_newer_H_zeon=.false.
      acc_newer_D_zeon=.false.
   endif

   if(acc_newer_H_zion) then
      !$acc update device(wi_init(:),wi_vpara(:),wi_vperp(:),wi_accel(:))
      !$acc update device(zion(:,:),zion0(:,:),ilabels(:))
      !$acc update device(UnNorm_theta_i(:),UnNorm_zeta_i(:))
      acc_newer_H_zion=.false.
      acc_newer_D_zion=.false.
   endif

   call Print_Reference_ntracer_Information

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
   use function_gg_module
   use Structured_Grid
   use Ptl_Src_Loss
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
   real(wp) Te_norm(memax) !  normalized-Te=Te/E_0, the same unit as normalized Ti
   integer m,ip,i,j,ij,mstart,mend
   integer gid, Ngroup
   real(wp), allocatable :: vangle_list(:),uabs_list(:)
   real(wp) jacobian,density
   real(wp) dpsi,zdum,tdum,pdum
   real(wp) integral_f0
   real(wp) v_abs, v_theta, vpara, vperp ! v_abs is amplitude of total velocity vector, and v_theta is the angle between vpara and vperp

   integer(kind=8) :: sobol_seed=1, sobol_dim_num=5 ! 5D random number (a,gq,gf,v_abs,v_theta)
   real(kind=8) rand_5d(5,mi) ! can be sobol quasi random numbers, or can be sprng() numbers
   real(wp) inv_ds, ff, P_s, P_q  !inv_ds is uniform by definition of s_grids
   integer:: sid ! index for s_grid
   integer:: Flag_tracers_at_specific_position=0

   real(wp) vmax,vc,c1,c2,c3,c4,f2val,c5,c6,ratio,gi_val(mi),ge_val(me), Rand !For Flag_Load_and_WeightDef==2 case.

   real(wp) intp_SG_full_phi
   real(wp), allocatable, dimension(:) :: AdiabTerm ! e*phi/T
   real(wp) :: target_w_tot

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

   ! ! Prepare sobol random numbers covering a whole device within a half-circle velocity space
   ! ! (a=(a0-a1), gq&gf=(0-2pi), v_abs=(0-umax)*sqrt(T/m), v_theta=(0-pi))
   ! ! Particles will be distributed into appropriate MPI ranks at the end of the initialization.
   ! sobol_seed = 1 + mype*mi
   ! do m=1,mi
   !    call i8_sobol(sobol_dim_num, sobol_seed, rand_5d(:,m))
   ! enddo

   ! Follwoing part is to prepare random numbers using sprng
   do m=1,mi
      rand_5d(1,m)=sprng()
      rand_5d(2,m)=sprng()
      rand_5d(3,m)=sprng()
      rand_5d(4,m)=sprng()
      rand_5d(5,m)=sprng()
   enddo

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

   call n_T_u_group  ! get T, n and omega_t, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
   !$acc update self(backu(:mi),backT(:mi))
   call esigetbf(zion,mi,size(zion,1))
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
         vpara=v_abs*cos(v_theta)
         vperp=v_abs*sin(v_theta)
         ! Get Marker Distribution function (g_val) at particle position considering Jacobian in the phase space
         gi_val(m) = Nm/((a1-a0)*4.0*pi*pi*0.5*pi*umax*umax*backT(m))/(Jacobian(zion(1,m),zion(2,m))*2.0*pi*vperp)

      elseif(Flag_Load_and_WeightDef==2) then
         ! First, get v_abs using rand_5d(4,m)
         vmax=umax*sqrt(backT(m))
         vc=sqrt(2.0)*sqrt(backT(m))
         c2=0.01
         f2val=0.05
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
         ratio=0.15 ! ratio of (max/min)
         c5=1.0/((1.0+ratio)*0.5*pi)
         c6=(1.0-ratio)/((1.0+ratio)*(0.5*pi)**2)
         Rand = rand_5d(5,m)
         v_theta = 0.5*pi + sign(1.0_wp,Rand-0.5)*((c5-sqrt(c5**2-2*c6*abs(Rand-0.5)))/c6);
         vpara=v_abs*cos(v_theta)
         vperp=v_abs*sin(v_theta)

         ! Get Marker Distribution function (g_val) at particle position considering Jacobian in the phase space
         if(v_abs<vc) then
            gi_val(m) = Nm*(c1*v_abs+c2)*(c5-c6*abs(v_theta-0.5*pi))/(v_abs)/((a1-a0)*4.0*pi*pi)/(Jacobian(zion(1,m),zion(2,m))*2.0*pi*vperp)
         else
            gi_val(m) = Nm*(c3*v_abs+c4)*(c5-c6*abs(v_theta-0.5*pi))/(v_abs)/((a1-a0)*4.0*pi*pi)/(Jacobian(zion(1,m),zion(2,m))*2.0*pi*vperp)
         endif
      endif

      zion(4,m) = vpara/BB(m)
      if(ifelectron/=0) zeon(4,m)=zion(4,m)  ! This is for the quiet start..
      zion(4,m) = zion(4,m) + backu(m)/BB(m) ! add backu term
      zion(6,m) = vperp*sqrt(1.0/(2.0*BB(m)))
      ilabels(m) = mype + 1e5*m
   enddo
   acc_newer_H_zion=.true.

   !$acc update device(ilabels(:mi))


   call n_T_u_group(flag_backf0=.true.)  ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))
   !$acc update self(back_f0_i(:mi))

   ! Set initial weights of ions
   do m=1,mi
      ! zion(5,m) is "df"
      ! zion0(8,m) is "f_tot"

      ! zion0(6,m) represents the phase volume (J(x,v)d5z) that preserves in the Hamiltonian system
      zion0(6,m) = 1.0/gi_val(m)
      integral_f0 = zion0(6,m)*back_f0_i(m)
      zion(5,m)  = (2.0*sprng()-1.0)*(1.0+cos(zion(2,m)))*w_initial_in*integral_f0
      zion0(8,m) = integral_f0 + zion(5,m)
   enddo
   acc_newer_H_zion=.true.

   ! zion(5,:) =0.0

   wi_init  = zion(5,:)
   wi_vpara = 0
   wi_vperp = 0
   wi_accel = 0

   if(mype==0) then
      open(unit=1212,file='test_zion_load.dat',action='write')
      do m=1,mi,100
         write(1212,'(I10,4(ES20.12,1x))') m, zion0(6,m), zion0(8,m), zion0(6,m)*back_f0_i(m), zion(5,m)
      enddo
      call flush(1212)
      close(1212)
   endif

   ! =================================================
   ! Electron Part
   ! =================================================

   ! ! Prepare sobol random numbers covering a whole device within a half-circle velocity space
   ! ! (a=(a0-a1), gq&gf=(0-2pi), v_abs=(0-umax)*sqrt(T/m), v_theta=(0-pi))
   ! ! Particles will be distributed into appropriate MPI ranks at the end of the initialization.
   ! sobol_seed = 1 + mype*me
   ! do m=1,me
   !    call i8_sobol(sobol_dim_num, sobol_seed, rand_5d(:,m))
   ! enddo

   ! Follwoing part is to prepare random numbers using sprng
   do m=1,me
      rand_5d(1,m)=sprng()
      rand_5d(2,m)=sprng()
      rand_5d(3,m)=sprng()
      rand_5d(4,m)=sprng()
      rand_5d(5,m)=sprng()
   enddo

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
      !$acc update self(backTe(:me))
      call esigetb(zeon,me,zeondim)

      do m=1,me
         Te_norm(m) = backTe(m)*msm0
      enddo

      if(Flag_Quiet_Start==1) then
         do m=1,me
            zeon(4,m)=zion(4,m)*sqrt(msm0)*c0cs*sqrt(Te_norm(m)/backT(m))  !checked; neglect electron rotation
            zeon(6,m)=zion(6,m)*sqrt(Te_norm(m)/backT(m))     ! =sqrt(mu_e)

            ! Phase space volume
            ! The power 1.5 is because of (v_abs^2) in area of velocity domain and (v_perp) of velocity Jacobian
            zeon0(6,m)=zion0(6,m)*(Te_norm(m)*m0ms/backT(m))**1.5
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
               vpara=v_abs*cos(v_theta)
               vperp=v_abs*sin(v_theta)
               ! Get Marker Distribution function (g_val) at particle position considering Jacobian in the phase space
               ge_val(m) = Nm/((a1-a0)*4.0*pi*pi*0.5*pi*umax*umax*Te_norm(m)*m0ms)/(Jacobian(zeon(1,m),zeon(2,m))*2.0*pi*vperp)

            elseif(Flag_Load_and_WeightDef==2) then
               ! First, get v_abs using rand_5d(4,m)
               vmax=umax*sqrt(Te_norm(m)*m0ms)
               vc=sqrt(2.0)*sqrt(Te_norm(m)*m0ms)
               c2=0.01
               f2val=0.05
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
               ratio=0.15 ! ratio of (max/min)
               c5=1.0/((1.0+ratio)*0.5*pi)
               c6=(1.0-ratio)/((1.0+ratio)*(0.5*pi)**2)
               Rand = rand_5d(5,m)
               v_theta = 0.5*pi + sign(1.0_wp,Rand-0.5)*((c5-sqrt(c5**2-2*c6*abs(Rand-0.5)))/c6);
               vpara=v_abs*cos(v_theta)
               vperp=v_abs*sin(v_theta)

               ! Get Marker Distribution function (g_val) at particle position considering Jacobian in the phase space
               if(v_abs<vc) then
                  ge_val(m) = Nm*(c1*v_abs+c2)*(c5-c6*abs(v_theta-0.5*pi))/(v_abs)/((a1-a0)*4.0*pi*pi)/(Jacobian(zeon(1,m),zeon(2,m))*2.0*pi*vperp)
               else
                  ge_val(m) = Nm*(c3*v_abs+c4)*(c5-c6*abs(v_theta-0.5*pi))/(v_abs)/((a1-a0)*4.0*pi*pi)/(Jacobian(zeon(1,m),zeon(2,m))*2.0*pi*vperp)
               endif
            else
               write(stdout,*) "Error!! Unknown Flag; Flag_Load_and_WeightDef=",Flag_Load_and_WeightDef
               call flush(stdout)
               stop
            endif

            zeon(4,m) = msm0*vpara/(csc0*BB(m))
            zeon(6,m) = vperp*sqrt(msm0/(2.0*BB(m)))
         enddo
      endif
      acc_newer_H_zeon=.true.

      call n_T_u_group_e(flag_backf0=.true.)   ! This is to get back_f0_e inside. (which requires zeon(4) and zeon(6))
      !$acc update self(back_f0_e(:me))

      ! Calculate Adiab term from initial potential phi
      ! Prepare SG's phi_3D in case
      call Update_SG_phi_3D_by_gathering_SG_phi2D

      allocate(AdiabTerm(me))
      !$acc parallel loop present(backTe)
      do m=1,me
         call Get_Intp_global_phi3D_on_SG(zeon(7,m), zeon(2,m),zeon(3,m), intp_SG_full_phi)

         if(Flag_Def_dne_adiab==1) then
            ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) for the New Poisson3D solver
            AdiabTerm(m)=abs(csc0)*intp_SG_full_phi/(backTe(m)*msm0) ! Adiabatic term (e*dphi_tot/Te)
         else
            AdiabTerm(m)=-1e100
            write(stdout,*) "Error!! Flag_Def_dne_adiab=",Flag_Def_dne_adiab, " is not defined in load_uniform_MGY!!"
            STOP
         endif
      enddo



      if(Flag_Quiet_Start==1) then
         do m=1,me
            integral_f0 = zeon0(6,m)*back_f0_e(m)

            target_w_tot = zion0(8,m) ! same to ion's total weight

            ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
            zeon(9,m) = integral_f0*AdiabTerm(m) ! f_AD
            zeon(8,m) = target_w_tot -int(ifelectron/2)*zeon(9,m)! full-f (ifelectron<2) or full-h (ifelectron==2)
            zeon(5,m)  = zeon(8,m)-integral_f0 ! df or dh
            ! zeon0(9,m) = integral_f0 ! f_0
            zeon0(9,m) = target_w_tot ! f_tot
            zeon0(8,m) = zeon(8,m)
         enddo
      else
         ! Set initial weights of electrons
         do m=1,me
            ! zeon(5,m) is "df"
            ! zeon(8,m) is "full-f or full-h"

            ! zeon0(6,m) represents the phase volume (J(x,v)d5z) that preserves in the Hamiltonian system
            zeon0(6,m) = 1.0/ge_val(m)
            integral_f0 = zeon0(6,m)*back_f0_e(m)

            ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
            zeon(9,m) = integral_f0*AdiabTerm(m) ! f_AD

            target_w_tot = integral_f0 + (2.0*sprng()-1.0)*(1.0+cos(zeon(2,m)))*w_initial_in*integral_f0

            zeon0(9,m) = target_w_tot ! f_tot
            zeon(8,m) = target_w_tot -int(ifelectron/2)*zeon(9,m)! full-f (ifelectron<2) or full-h (ifelectron==2)
            zeon(5,m)  = zeon(8,m)-integral_f0 ! df or dh

            zeon0(8,m) = zeon(8,m)
         enddo
      endif
      acc_newer_H_zeon=.true.

      ! zeon(5,:)=0.0

      if(mype==0) then

         open(unit=1213,file='test_zeon_load.dat',action='write')
         do m=1,me,100
            write(1213,'(I10,4(ES20.12,1x))') m, zeon0(6,m), zeon(8,m), zeon0(6,m)*back_f0_e(m), zeon(5,m)
         enddo
         call flush(1213)
         close(1213)
      endif

      we_init(:)  = zeon(5,:)
      we_vpara = 0
      we_vperp = 0
      we_accel = 0
      we_adiab = 0
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
         m = 0.5*mi
         write(stdout,*) m,"-th ion's each_ptl_d6z=", inv_Nm*((a1-a0)*4.0*pi*pi*2.0*umax*sqrt(2.0)*umax*backT(m)), ", back_f0_i(m)=",back_f0_i(m)
         call flush(stdout)

         m = 0.5*me
         write(stdout,*) m,"-th electron's each_ptl_d6z=", inv_Nm*((a1-a0)*4.0*pi*pi*2.0*umax*sqrt(2.0)*umax*Te_norm(m)*m0ms),  ", back_f0_e(m)=",back_f0_e(m)
         call flush(stdout)
      endif


      Flag_tracers_at_specific_position=1
      if(Flag_tracers_at_specific_position==1) then

         ! Set test etracers
         ! Set positions of trapped particles

         Ngroup=6
         allocate(vangle_list(Ngroup),uabs_list(Ngroup))
         ! vangle_list=[pi/6,pi/3,80*pi/180,pi/6,pi/3,80*pi/180] ! 30, 60, 80, degrees of angle
         ! vangle_list=[30*pi/180,50*pi/180,60*pi/180,30*pi/180,50*pi/180,60*pi/180] ! 30, 50, 70, degrees of angle
         ! uabs_list=[1.0,1.0,1.0,2.5,2.5,2.5]
         ! do gid=1,Ngroup
         !    mstart=(gid-1)*int(num_of_etracers/ngroup)+1
         !    mend  =(gid)*int(num_of_etracers/ngroup)
         !    do i=mstart,mend
         !       m = (real(i)/(num_of_etracers+1.0)) * me
         !       ! zeon(1,m) = 0.64 + 0.02*(i-mstart)/real(mend-mstart)
         !       zeon(1,m) = 0.65 + 0.02*(i-mstart)/real(mend-mstart)
         !       ! zeon(1,m) = 0.79 + 0.02*(i-mstart)/real(mend-mstart)
         !       zeon(2,m) = 0.0
         !       zeon(3,m) = 0.0
         !       zeon(7,m)=min(s1,max(s0,D01AHF(0.0_wp,zeon(1,m),error_8,minpts,acc,dsda,-1,ifail)))
         !    enddo
         ! enddo
         ! acc_newer_H_zeon=.true.
         ! vangle_list=[30*pi/180,56*pi/180,63*pi/180,(180-30)*pi/180,(180-56)*pi/180,(180-63)*pi/180] ! 30, 56, 63, degrees of angle
         ! uabs_list=[1.0,1.0,1.0,1.0,1.0,1.0]

         ! vangle_list=[30*pi/180,56*pi/180,63*pi/180,(180-30)*pi/180,(180-56)*pi/180,(180-63)*pi/180] ! 30, 56, 63, degrees of angle
         uabs_list=[2.0,4.0,6.0,8.0,10.0,12.0]
         do gid=1,Ngroup
            mstart=(gid-1)*int(num_of_etracers/ngroup)+1
            mend  =(gid)*int(num_of_etracers/ngroup)
            do i=mstart,mend
               m = (real(i)/(num_of_etracers+1.0)) * me
               ! zeon(1,m) = 0.64 + 0.02*(i-mstart)/real(mend-mstart)
               zeon(1,m) =0.876
               ! zeon(1,m) = 0.79 + 0.02*(i-mstart)/real(mend-mstart)
               zeon(2,m) = 0.314
               zeon(3,m) = 0.0
               zeon(7,m)=min(s1,max(s0,D01AHF(0.0_wp,zeon(1,m),error_8,minpts,acc,dsda,-1,ifail)))
            enddo
         enddo
         acc_newer_H_zeon=.true.

         call n_T_u_group_e   ! get Te, note that back_f0_e cannot be calculated here because it requires zion(4) and zion(6)
         !$acc update self(backTe(:me))
         call esigetb(zeon,me,zeondim)

         do gid=1,Ngroup
            mstart=(gid-1)*int(num_of_etracers/Ngroup)+1
            mend  =(gid)*int(num_of_etracers/Ngroup)
            do i=mstart,mend
               m = (real(i)/(num_of_etracers+1.0)) * me
               Te_norm(m)= backTe(m)*msm0
               vpara = uabs_list(gid)*sqrt(2*Te_norm(m)*m0ms)*cos(pi*(i-mstart)/real(mend-mstart)) ! Most Probable Veclotiy vp = sqrt(2T/m)
               vperp = uabs_list(gid)*sqrt(2*Te_norm(m)*m0ms)*sin(pi*(i-mstart)/real(mend-mstart)) ! Most Probable Veclotiy vp = sqrt(2T/m)
               zeon(4,m) = msm0*vpara/(csc0*BB(m))
               zeon(6,m) = vperp*sqrt(msm0/(2.0*BB(m)))

               zeon0(6,m) = 0.0
               zeon(8,m)  = 0.0
               zeon(5,m)  = 0.0
               zeon0(8,m) = 0.0
               zeon0(9,m) = 0.0

            enddo
         enddo
         acc_newer_H_zeon=.true.

         call n_T_u_group_e(flag_backf0=.true.)   ! This is to get back_f0_e inside. (which requires zeon(4) and zeon(6))
         !$acc update self(back_f0_e(:me))




         ! Set itracers
         do gid=1,Ngroup
            mstart=(gid-1)*int(num_of_itracers/Ngroup)+1
            mend  =(gid)*int(num_of_itracers/Ngroup)
            do i=mstart,mend
               m = (real(i)/(num_of_itracers+1.0)) * mi
               zion(1,m) = 0.64 + 0.02*(i-mstart)/real(mend-mstart)
               zion(2,m) = 0.0
               zion(3,m) = 0.0
               zion(7,m)=min(s1,max(s0,D01AHF(0.0_wp,zion(1,m),error_8,minpts,acc,dsda,-1,ifail)))
            enddo
         enddo
         acc_newer_H_zion=.true.


         call n_T_u_group   ! get Ti, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
         !$acc update self(backT(:mi))
         call esigetb(zion,mi,ziondim)



         do gid=1,Ngroup
            mstart=(gid-1)*int(num_of_itracers/Ngroup)+1
            mend  =(gid)*int(num_of_itracers/Ngroup)
            do i=mstart,mend
               m = (real(i)/(num_of_itracers+1.0)) * mi
               vpara = uabs_list(gid)*sqrt(2*backT(m))*cos(pi*(i-mstart)/real(mend-mstart)) ! Most Probable Veclotiy vp = sqrt(2T/m)
               vperp = uabs_list(gid)*sqrt(2*backT(m))*sin(pi*(i-mstart)/real(mend-mstart)) ! Most Probable Veclotiy vp = sqrt(2T/m)
               zion(4,m) = vpara/(BB(m))
               zion(6,m) = vperp*sqrt(1.0/(2.0*BB(m)))

               zion0(6,m) = 0.0
               zion(8,m)  = 0.0
               zion(5,m)  = 0.0
               zion0(8,m) = 0.0 ! Full-f Ni=0, otherwise this can make error in the case of Flag_w_Evol=3
            enddo
         enddo
         acc_newer_H_zion=.true.

         call n_T_u_group(flag_backf0=.true.)   ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))
         !$acc update self(back_f0_i(:mi))

         deallocate(vangle_list,uabs_list)

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

   ! To check whether some particles outside the boundary
   ! If Ptl_BC_type=2, ion boundary checking considers the gyro-radius to determien whether it touches the wall
   call treat_ions_out_of_boundary
   !$acc update self(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi),ilabels(:mi))


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


   if(acc_newer_H_zeon) then
      !$acc update device(we_init(:),we_vpara(:),we_vperp(:),we_accel(:),we_adiab(:))
      !$acc update device(zeon(:,:),zeon0(:,:),elabels(:))
      !$acc update device(UnNorm_theta_e(:),UnNorm_zeta_e(:))
      acc_newer_H_zeon=.false.
      acc_newer_D_zeon=.false.
   endif

   if(acc_newer_H_zion) then
      !$acc update device(wi_init(:),wi_vpara(:),wi_vperp(:),wi_accel(:))
      !$acc update device(zion(:,:),zion0(:,:),ilabels(:))
      !$acc update device(UnNorm_theta_i(:),UnNorm_zeta_i(:))
      acc_newer_H_zion=.false.
      acc_newer_D_zion=.false.
   endif


   ! Distribute ions and electrons to appropriate MPI ranks
   ! NOTE: (08/23/2020)
   ! If using OpenACC, these shift routines must be at the last part
   ! because particles got scattered to different PEs in the DEVICE side and not synced with the HOST side.
   call shifti
   call shifte

   ini_mi=mi
   ini_me=me


   call MPI_ALLREDUCE(ini_mi,ini_mi_npartd,1,MPI_INTEGER,MPI_SUM,partd_comm,ierror)
   call MPI_ALLREDUCE(ini_me,ini_me_npartd,1,MPI_INTEGER,MPI_SUM,partd_comm,ierror)


   deallocate(AdiabTerm)

   ! write(stdout,'((A,I4,A,I8,A,I8))') "In load_uniform_MGY, I'm=",mype, ", ini_mi=", ini_mi, ", ini_me=",ini_me
   write(stdout,'(5(A,I0))') "In load_uniform_MGY, I'm=",mype, ", ini_mi=", ini_mi, ", ini_mi_npartd=",ini_mi_npartd, ", ini_me=",ini_me,", ini_me_npartd=",ini_me_npartd
   call flush(stdout)

   call Print_Reference_ntracer_Information

end subroutine load_uniform_MGY


subroutine Adjust_electron_weight_for_initially_zero_potential
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use field_array
   use diagnosis_array
   use esi_interface
   use control_parameters
   use function_gg_module
   use Structured_Grid
   use Ptl_Src_Loss
   use solver_data_class
   use MGY_PETSc_Interface_Data
   implicit none


   integer m,ierror
   real(wp) :: AdiabTerm ! e*phi/T
   real(wp) intp_SG_full_phi, integral_f0
   real(wp) new_f_AD, new_f0

   if(Flag_Poisson3D_on_SG/=1) then
      write(stdout,*) "Error!! Flag_Poisson3D_on_SG=",Flag_Poisson3D_on_SG, " is not defined in Adjust_electron_weight_for_initially_zero_potential!!"
      STOP
   endif


   istep=1

   call n_T_u_group
   call n_T_u_group_e

   call Measure_Ion_Density_Current_on_SG(1)
   call Measure_Electron_Density_Current_on_SG(1)

   if(ifelectron==2 .and. Flag_Def_dne_adiab==2) then
      call Solve_Nonlinear_Poisson3D_on_SG
   else
      call Solve_Poisson3D_on_SG
   endif
   call Cal_gradPhi_on_SG

   call Update_SG_phi_3D_by_gathering_SG_phi2D


   !$acc parallel loop present(backTe,back_f0_e)
   do m=1,me
      call Get_Intp_global_phi3D_on_SG(zeon(7,m), zeon(2,m),zeon(3,m), intp_SG_full_phi)

      if(Flag_Def_dne_adiab==1) then
         ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) for the New Poisson3D solver
         ! considering a conversion from V unit to GTS normalized unit
         AdiabTerm=abs(csc0)*intp_SG_full_phi/(backTe(m)*msm0) ! Adiabatic term (e*dphi_tot/Te)
      else
         AdiabTerm=-1e100
         write(stdout,*) "Error!! Flag_Def_dne_adiab=",Flag_Def_dne_adiab, " is not defined in pushe_new.F90!!"
         STOP
      endif

      new_f0 = zeon0(6,m)*back_f0_e(m)
      new_f_AD = zeon0(6,m)*back_f0_e(m)*AdiabTerm

      ! f_tot = f0 + df (when ifelectron<2)
      ! f_tot = f0 + dh + f_AD (when ifelectron==2)
      zeon0(9,m) = new_f0 +zeon(5,m) + int(ifelectron/2)*new_f_AD

      zeon(8,m) = zeon0(9,m) - (int(ifelectron/2)*new_f_AD) ! full-f or full-h 
   enddo
   acc_newer_D_zeon=.true.

   call Measure_Electron_Density_Current_on_SG(1)
   if(ifelectron==2 .and. Flag_Def_dne_adiab==2) then
      call Solve_Nonlinear_Poisson3D_on_SG
   else
      call Solve_Poisson3D_on_SG
   endif
   call Cal_gradPhi_on_SG
   call Update_SG_phi_3D_by_gathering_SG_phi2D

   !$acc parallel loop present(backTe,back_f0_e)
   do m=1,me
      call Get_Intp_global_phi3D_on_SG(zeon(7,m), zeon(2,m),zeon(3,m), intp_SG_full_phi)

      if(Flag_Def_dne_adiab==1) then
         ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) for the New Poisson3D solver
         AdiabTerm=abs(csc0)*intp_SG_full_phi/(backTe(m)*msm0) ! Adiabatic term (e*dphi_tot/Te)
      else
         AdiabTerm=-1e100
         write(stdout,*) "Error!! Flag_Def_dne_adiab=",Flag_Def_dne_adiab, " is not defined in pushe_new.F90!!"
         STOP
      endif

      integral_f0 = zeon0(6,m)*back_f0_e(m)

      zeon(9,m) = integral_f0*AdiabTerm ! f_AD

      zeon0(9,m) =zeon(8,m) + int(ifelectron/2)*zeon(9,m) ! f_tot
      zeon(5,m) = zeon(8,m)-integral_f0 ! df or dh

      ! if(mype==0 .and. m==1) then
      !    write(stdout,*) 'After adjustment2, zeon0(9,m)=',zeon0(9,m), '; zeon(8,m)=', zeon(8,m),'; zeon(9,m)=', zeon(9,m)
      ! endif
      ! zeon(5,m) = 0.0
   enddo
   acc_newer_D_zeon=.true.


   if(acc_newer_D_zeon) then
      !$acc update self(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
      acc_newer_D_zeon=.false.
      acc_newer_H_zeon=.false.
   endif

   istep=0

end subroutine Adjust_electron_weight_for_initially_zero_potential
