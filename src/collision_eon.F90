subroutine collision_eon
   use mpi_defs
   use global_parameters
   use particle_array
   use field_array
   use esi_interface
   use function_gg_module

   implicit none

   integer i

   !    NOTE: output from call n_T_u_group is omega_t and d(omega_t)/da;  12/07
   !          B(=BB) and I(=F) at particles are needed for calculating u_//=I*omega_t/B;
   !          BB, I can be obtained after call gcmotion, or call esigetbf, esiget2dfunctions.

   call n_T_u_group_e  ! prepare n, T, dlnn/da, dlnT/da, etc. which are also used for PUSHE/PUSHEON

   if(ecollision.gt.0) then   ! include collision
      i=mod(istep,ecollision) 
      if(i==0) then           ! calculate collisions every ecollision-step
         if(irk==1) then      ! calculate collisions at irk=1

            ! Update needed particle information from DEVICE to HOST
            !$acc update self(backne(:me))

            call esigetb(zeon,me,zeondim)   ! get B-->BB

            !           call esigetbf(zion,mi)   ! get B-->BB, I-->F  12/07
            !           backu(1:mi)=backu(1:mi)*F(1:mi)/BB(1:mi) ! =u_//=I*omega_t/B

            call Lorentz_coll

            call CTPandCELL_eon

            !  make sure output backu=omega_t after exit this subroutine
            !           backu(1:mi)=backu(1:mi)*BB(1:mi)/F(1:mi) ! =omega_t

         end if
      end if
   end if   

end subroutine collision_eon



SUBROUTINE Lorentz_coll

   use precision
   use global_parameters
   use particle_array
   use esi_interface

   implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

   integer,dimension(:) :: i1(me),i2(me)
   real*8 sff, pitch, panew,unitf,f_rho,rxi,vperp2,Binv
   integer i,j,k,j1,j2,iss

   real(wp) x(me), x2(me), vz(me)   ! working arrays

   !      C_{pitch-angle}:  Cells are not 
   !              necessary, but taking local collision frequency into 
   !              account through n(r) and T(r)           

   !     dlne=16.0d0
   !     Tee=1.0d0/4.0d0/pi*amasse
   !    &    *(vth0)**3.0d0/(backn0*charge*charge
   !    &    *charge*charge*dlne)*amasse
   !     dt=dtau*Tee

   f_rho=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)

   !$omp parallel do private(i,vperp2)
   do i=1,me
      vz(i)=zeon(4,i)*BB(i)*f_rho   ! normalized v_//
      vperp2=zeon(6,i)*zeon(6,i)*BB(i)*2.0_wp*m0ms    ! normalized vper^2
      x(i)=sqrt(vperp2+vz(i)**2)      ! x is particle velocity in lab frame
      x2(i)=dt_coll_e*backne(i)/(x(i)**3)   ! \nu_e*dt_e
   enddo

   j1=0
   j2=0

   do i=1,me
      if(x2(i).lt.1.0d0) then
         j1=j1+1
         i1(j1)=i
      else
         j2=j2+1
         i2(j2)=i
      endif
   end do

   ! The average value of "j2" at the beginning of the simulation is about "10",
   ! with a maximum around "16", so it is small compared to the total number of
   ! particles "me". We still have to check how this number evolves during the
   ! simulation though.

   do i=1,j1
      zeon(4,i)=sprng()    ! borrow zeon(4,i) to store random number
   end do

   !$omp parallel do private(i,j,pitch,unitf,panew)
   do j=1,j1
      i=i1(j)  !There is a unique "i" for each "j" so no conflicts between threads
      pitch=vz(i)/x(i)
      unitf=dsign(1.0d0,(zeon(4,j)-0.5d0))
      panew=max(-1.0d0,min(1.0d0,    &
         pitch*(1.0d0-x2(i))+unitf*sqrt((1.0d0-pitch**2)*x2(i)) ))

      vz(i)=panew*x(i)
   end do   

   sff=0.4d0
   do j=1,j2
      i=i2(j)
      pitch=vz(i)/x(i)
      iss=x2(i)/sff+1
      iss=min(iss,100000)  ! iss could be very huge for small v particles
      x2(i)=x2(i)/dfloat(iss)

      do k=1,iss
         rxi=sprng()
         unitf=dsign(1.0d0,(rxi-0.5d0))
         panew=max(-1.0d0,min(1.0d0,    &
            pitch*(1.0d0-x2(i))+unitf*sqrt((1.0d0-pitch**2)*x2(i)) ))
         pitch=panew
      end do
      vz(i)=panew*x(i)
   end do

   !  if(me==63) then
   !     do i=1,me
   !       if((x(i)**2-vz(i)**2).lt.0.0d0) then 
   !          print*, '@@@@@@@@@@i=',i, x(i)**2, vz(i)**2, (x(i)**2-vz(i)**2)
   !       end if
   !     end do
   !  end if

   !  update the following quantities
   !$omp parallel do private(i,vperp2,Binv)
   do i=1,me
      vperp2=max(0.0d0,(x(i)**2-vz(i)**2))    ! vper^2
      Binv=1.0_wp/BB(i)
      zeon(4,i)=vz(i)*Binv*for_rho
      zeon(6,i)=sqrt(vperp2*Binv*0.5_wp*msm0)   ! sqrt(mu) of electrons
   enddo

end subroutine Lorentz_coll



SUBROUTINE CTPandCELL_eon

   use global_parameters
   use particle_array
   use field_array
   use esi_interface

   implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

   integer m
   real*8 da_size, dang, ang, rn1, rn2, d1, f_rho, Binv
   real*8 fp, fr, fsr, fsp, fpr, HH, GG, FF, vxy2dum, f0dt, rv0, h1, h2, h3
   integer i, kc, nc, j, k, jx, iii(me)

   real(wp) dgr(me), da(me), dgq(me), dgf(me), mu(me), Er_particle(me)   ! working arrays

   !   real(wp), EXTERNAL :: S15AEF   ! The replacement of nag S15AEF DOES NOT do correct at least on Franklin!!!
   real*8 derf

   real*8, dimension(:), allocatable :: Afsp,Af0dt,Afp,Afsr,Afr,Afpr

   !
   !      C_{TP}: The drag and diffusion operator. Cells are not 
   !              necessary, but taking local collision frequency into 
   !              account through n(r) and T(r)           
   !

   d1=2.0_wp/sqrt(pi)

   f_rho=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)

   !$omp parallel do private(i)
   do i=1,me
      dgr(i)=zeon(4,i)*BB(i)*f_rho                 ! normalized v_//
      da(i)=zeon(6,i)*zeon(6,i)*BB(i)*2.0_wp*m0ms  ! normalized vper^2
      dgq(i)=dgr(i)**2+da(i)                       ! v^2      
      dgf(i)=dgq(i)/backTe(i)*0.5_wp               ! v^2/v_th^2 with v_th^2=2T/m
      mu(i)=sqrt(dgf(i))                           ! v/v_th
      gha(i)=d1*mu(i)*exp(-dgf(i))                 ! 2.0d0/sqrt(pi)*x*dexp(-x^2) with x=v/v_th
      Er_particle(i)=derf(mu(i))-gha(i)
   enddo

   do i=1,me
      zeon(4,i)=sprng()
   end do

   do i=1,me
      zeon(6,i)=sprng()
   end do

   !$omp parallel do private(i,h1,h2,h3,FF,GG,HH,fsp,fsr,fp,fr,fpr,rv0,f0dt)
   do i=1,me
      h1=1.0_wp/dgf(i)
      h2=da(i)
      h3=dgr(i)**2
      FF=2.0_wp*Er_particle(i)
      GG=(1.0_wp-0.5_wp*h1)*Er_particle(i)+gha(i)
      HH=Er_particle(i)*h1
      fsp=dgr(i)*FF
      fsr=2.0_wp*h2*FF-h2*HH-(2.0_wp*h3+h2)*GG
      fp=h3*HH+h2*GG
      fr=4.0_wp*h2*(h2*HH+h3*GG)
      fpr=2.0_wp*h2*dgr(i)*(HH-GG)
      rv0=1.0_wp/(sqrt(dgq(i))*dgq(i))
      f0dt=dt_coll_e*backne(i)*rv0    

      ghq(i)=-fsp*f0dt+(zeon(4,i)-0.5_wp)*sqrt(12.0_wp*fp*f0dt)    ! dvz

      TT(i)=h2-fsr*f0dt+(zeon(6,i)-0.5_wp)*          &
         sqrt(12.0_wp*(fr-fpr*fpr/fp)*f0dt)+     &
         (zeon(4,i)-0.5_wp)*fpr/fp*sqrt(12.0_wp*fp*f0dt) ! new v_perp^2
   end do

   !  particle sorting in 2D grids

   !       dang=topi/dfloat(np_grid)
   !       da_size=1.0d0/dfloat(n_grid)
   !       mm=0
   ! 
   !       do i=1,me
   !         ang=dmod(theta(i),topi)+topi
   !         ang=dmod(ang,topi)
   !         ang=ang/dang
   ! 
   !         nc=idint(ang)+1
   !         jpol(i)=min(np_grid,max(1,nc))
   ! 
   !         kc=idint(aa(i)/da_size)+1
   !         jpsi(i)=min(n_grid,max(1,kc))
   !       end do

   ! Test for electrons that need special treatment
   m=0
   do i=1,me
      vxy2dum=TT(i)
      if(vxy2dum.ge.0.0_wp) then
         dgr(i)=dgr(i)+ghq(i)
         da(i)=TT(i)
      else
         m=m+1
         iii(m)=i
      endif
   end do

   ! The average value of "m" at the beginning of the simulation is about "3",
   ! with a maximum around "6", so it is very small compared to the total number of
   ! particles "me". We still have to check how this number evolves during the
   ! simulation though.

   !
   !       For particles with vxy<0, re-calculate using 
   !       small step size---dt_coll/nss
   !
   if(m.gt.0) then
      allocate(Afsp(m),Af0dt(m),Afp(m),Afsr(m),Afr(m),Afpr(m))
      do 5 k=1,nss
         do i=1,m
            zeon(4,i)=sprng()
         end do
         do i=1,m
            zeon(6,i)=sprng()
         end do

         do j=1,m
            i=iii(j)
            dgq(i)=dgr(i)**2+da(i)
            dgf(i)=dgq(i)/backTe(i)*0.5_wp
            mu(i)=sqrt(dgf(i))
            gha(i)=d1*mu(i)*exp(-dgf(i))
            Er_particle(i)=derf(mu(i))-gha(i)
         end do

         do 4 j=1,m
            i=iii(j)
            h1=1.0_wp/dgf(i)
            h2=da(i)
            h3=dgr(i)**2
            FF=2.0_wp*Er_particle(i)
            GG=(1.0_wp-0.5_wp*h1)*Er_particle(i)+gha(i)
            HH=Er_particle(i)*h1
            Afsp(j)=dgr(i)*FF

            Afsr(j)=2.0_wp*h2*FF-h2*HH-(2.0_wp*h3+h2)*GG
            Afp(j)=h3*HH+h2*GG
            Afr(j)=4.0_wp*h2*(h2*HH+h3*GG)
            Afpr(j)=2.0_wp*h2*dgr(i)*(HH-GG)
            rv0=1.0_wp/(sqrt(dgq(i))*dgq(i))
            Af0dt(j)=dt_coll_e/real(nss,wp)*backne(i)*rv0
            ! smaller step size is dt_coll/dfloat(nss)
            ghq(i)=-Afsp(j)*Af0dt(j)+(zeon(4,j)-0.5_wp)*sqrt(12._wp*Afp(j)*Af0dt(j))  ! dvz
            TT(i)=h2-Afsr(j)*Af0dt(j)+(zeon(6,j)-0.5_wp)*                     &
               sqrt(12.0_wp*(Afr(j)-Afpr(j)*Afpr(j)/Afp(j))*Af0dt(j))+    &
               (zeon(4,j)-0.5_wp)*Afpr(j)/Afp(j)*sqrt(12.0_wp*Afp(j)*Af0dt(j))
            4       continue

            do j=1,m
               i=iii(j)
               do while (TT(i).lt.0.0_wp)
                  rn1=sprng()
                  rn2=sprng()
                  ghq(i)=-Afsp(j)*Af0dt(j)+(rn1-0.5_wp)*sqrt(12.0_wp*Afp(j)*Af0dt(j))
                  TT(i)=da(i)-Afsr(j)*Af0dt(j)+(rn2-0.5_wp)*                  &
                     sqrt(12.0_wp*(Afr(j)-Afpr(j)*Afpr(j)/Afp(j))*Af0dt(j))+   &
                     (rn1-0.5_wp)*Afpr(j)/Afp(j)*sqrt(12.0_wp*Afp(j)*Af0dt(j))
               enddo
               dgr(i)=dgr(i)+ghq(i)
               da(i)=TT(i)
            enddo
            5    continue

            deallocate(Afsp,Af0dt,Afp,Afsr,Afr,Afpr)

         end if

         !$omp parallel do private(i,Binv)
         do i=1,me
            Binv=1.0/BB(i)
            !    rho=vz/B
            zeon(4,i)=dgr(i)*Binv*for_rho               ! rho=v_para/B of electrons
            zeon(6,i)=sqrt(0.5_wp*da(i)*Binv*msm0)      ! sqrt(mu) of electrons
         enddo

      end subroutine CTPandCELL_eon
