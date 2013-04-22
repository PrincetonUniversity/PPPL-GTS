subroutine collision_eon
  use precision
  use mpi_defs
  use global_parameters
  use particle_array
  use field_array
  use esi_interface

  implicit none

  integer i
  real(wp) :: t0,t1,t2,t3,t4
  real(wp),SAVE :: t_n_Tu=0.0_wp,t_esigetb=0.0_wp,t_lorentz=0.0_wp,t_ctp=0.0_wp

!    NOTE: output from call n_T_u_group is omega_t and d(omega_t)/da;  12/07
!          B(=BB) and I(=F) at particles are needed for calculating u_//=I*omega_t/B;
!          BB, I can be obtained after call gcmotion, or call esigetbf, esiget2dfunctions.

  t0=MPI_WTIME()
  call n_T_u_group_e  ! prepare n, T, dlnn/da, dlnT/da, etc. which are also used for PUSHE/PUSHEON
  t1=MPI_WTIME()

  t_n_Tu=t_n_Tu+(t1-t0)
 
  if(ecollision.gt.0) then   ! include collision
     i=mod(istep,ecollision) 
     if(i==0) then           ! calculate collisions every ecollision-step
        if(irk==1) then      ! calculate collisions at irk=1
           call esigetb(zeon,me,zeondim)   ! get B-->BB
          t2=MPI_WTIME()

!           call esigetbf(zion,mi)   ! get B-->BB, I-->F  12/07
!           backu(1:mi)=backu(1:mi)*F(1:mi)/BB(1:mi) ! =u_//=I*omega_t/B

           call Lorentz_coll

          t3=MPI_WTIME()
           call CTPandCELL_eon
          t4=MPI_WTIME()

          t_esigetb=t_esigetb+(t2-t1)
          t_lorentz=t_lorentz+(t3-t2)
          t_ctp=t_ctp+(t4-t3)

!  make sure output backu=omega_t after exit this subroutine
!           backu(1:mi)=backu(1:mi)*BB(1:mi)/F(1:mi) ! =omega_t

        end if
     end if
  end if   

  if(istep==mstep .and. mype==0)then
    write(0,'("in collision_eon: t(n_T_u_eon)=",f12.6,/,&
             &"                  t(esigetb)=",f12.6,/,&
             &"                  t(Lorentz)=",f12.6,/,&
             &"                  t(CTPandCell)=",f12.6,/)')t_n_Tu,t_esigetb,t_lorentz,t_ctp
  endif

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

     real(wp) sff, pitch, panew,unitf,f_rho_inv,rxi,ran_num
     integer i,j,k,j1,j2,iss
 
     real(wp) x, x2, vz

!      C_{pitch-angle}:  Cells are not 
!              necessary, but taking local collision frequency into 
!              account through n(r) and T(r)           
 
!     dlne=16.0d0
!     Tee=1.0d0/4.0d0/pi*amasse
!    &    *(vth0)**3.0d0/(backn0*charge*charge
!    &    *charge*charge*dlne)*amasse
!     dt=dtau*Tee

      f_rho_inv=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)
      sff=0.4d0

!$omp parallel do private(i,vz,x2,x,pitch,rxi,unitf,panew,iss,k)
      do i=1,me
         vz=zeon(4,i)*BB(i)*f_rho_inv   ! normalized v_//
         x2=zeon(6,i)*zeon(6,i)*BB(i)*2.0_wp*m0ms    ! normalized vper^2

         x=sqrt(x2+vz**2)      ! x is particle velocity in lab frame
         x2=dt_coll_e*backne(i)/(x**3)   ! \nu_e*dt_e

         if(x2 < 1.0d0) then
            pitch=vz/x
            rxi=sprng()
            unitf=dsign(1.0d0,(rxi-0.5d0))
            panew=max(-1.0d0,min(1.0d0,    &
                  pitch*(1.0d0-x2)+unitf*sqrt((1.0d0-pitch**2)*x2) ))
            vz=panew*x
         else  ! x2 >= 1
            pitch=vz/x
            iss=int(x2/sff)+1
            iss=min(iss,100000)  ! iss could be very huge for small v particles
            x2=x2/real(iss,wp)
 
            do k=1,iss
               rxi=sprng()
               unitf=dsign(1.0d0,(rxi-0.5d0))
               panew=max(-1.0d0,min(1.0d0,    &
                     pitch*(1.0d0-x2)+unitf*sqrt((1.0d0-pitch**2)*x2) ))
               pitch=panew
            end do
            vz=panew*x
         endif

       ! update zeon(4,i) and zeon(6,i) for new vz and v_perp values
         x2=max(0.0d0,(x**2-vz**2))    ! vper^2
         zeon(4,i)=vz/BB(i)*for_rho
         zeon(6,i)=sqrt(x2/BB(i)*0.5_wp*msm0)   ! sqrt(mu) of electrons
      end do

!  if(me==63) then
!     do i=1,me
!       if((x(i)**2-vz(i)**2).lt.0.0d0) then 
!          print*, '@@@@@@@@@@i=',i, x(i)**2, vz(i)**2, (x(i)**2-vz(i)**2)
!       end if
!     end do
!  end if


end subroutine Lorentz_coll



SUBROUTINE CTPandCELL_eon

  use precision
  use global_parameters
  use particle_array
  use field_array
  use esi_interface
! The esi_interface module contains the following list of 1D vector variables:
!  real(wp), dimension(mimax) :: F,Fa,gFa,gFaa,gYa,gYaa,TT,Ta,PP,Pa, &
!         rr,ra,rq,z,za,zq,BB,Ba,Bq,gh,gha,ghq,sw


  implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

   integer m
   real*8 da_size, dang, ang, rn1, rn2, d1, f_rho_inv
   real*8 fp, fr, fsr, fsp, fpr, HH, GG, FF, vxy2dum, f0dt, rv0, h1, h2, h3
   integer i, kc, nc, j, k, jx

   real(wp) vpara, vperp2, new_vperp2, dvz, dgq, vsq, mu, Er_particle, fMv, Binv

!   real(wp), EXTERNAL :: S15AEF   ! The replacement of nag S15AEF DOES NOT do correct at least on Franklin!!!
   !real*8 derf

!
!      C_{TP}: The drag and diffusion operator. Cells are not 
!              necessary, but taking local collision frequency into 
!              account through n(r) and T(r)           
!

      d1=2.0_wp/sqrt(pi)

      f_rho_inv=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)

!$omp parallel do private(vpara,vperp2,dgq,vsq,mu,fMv,Er_particle,rn1,rn2, &
!$omp& h1,h2,h3,FF,GG,HH,fsp,fp,fr,rv0,f0dt,dvz,new_vperp2)
      do i=1,me
         vpara=zeon(4,i)*BB(i)*f_rho_inv                     ! normalized v_//
         vperp2=zeon(6,i)*zeon(6,i)*BB(i)*2.0_wp*m0ms  ! normalized vper^2
         dgq=vpara**2+vperp2                             ! v^2      
         vsq=dgq/backTe(i)*0.5_wp                  ! v^2/v_th^2 with v_th^2=2T/m
         mu=sqrt(vsq)                              ! v/v_th
         fMv=d1*mu*exp(-vsq)                       ! 2.0d0/sqrt(pi)*x*dexp(-x^2) with x=v/v_th
         Er_particle=derf(mu)-fMv

         rn1=sprng()
         rn2=sprng()

         h1=1.0_wp/vsq
         h2=vperp2
         h3=vpara**2
         FF=2.0_wp*Er_particle
         GG=(1.0_wp-0.5_wp*h1)*Er_particle+fMv
         HH=Er_particle*h1
         fsp=vpara*FF
         fsr=2.0_wp*h2*FF-h2*HH-(2.0_wp*h3+h2)*GG
         fp=h3*HH+h2*GG
         fr=4.0_wp*h2*(h2*HH+h3*GG)
         fpr=2.0_wp*h2*vpara*(HH-GG)
         rv0=1.0_wp/(sqrt(dgq)*dgq)
         f0dt=dt_coll_e*backne(i)*rv0    

         dvz=-fsp*f0dt+(rn1-0.5_wp)*sqrt(12.0_wp*fp*f0dt)    ! dvz

         new_vperp2=h2-fsr*f0dt+(rn2-0.5_wp)*          &
                sqrt(12.0_wp*(fr-fpr*fpr/fp)*f0dt)+     &
                (rn1-0.5_wp)*fpr/fp*sqrt(12.0_wp*fp*f0dt) ! new v_perp^2

         if(new_vperp2 >= 0.0_wp) then
            vpara=vpara+dvz
            vperp2=new_vperp2
         else     ! if vxy < 0, re-calculate using small step size---dt_coll/nss
            do k=1,nss
               rn1=sprng()
               rn2=sprng()
               dgq=vpara**2+vperp2
               vsq=dgq/backTe(i)*0.5_wp
               mu=sqrt(vsq)
               fMv=d1*mu*exp(-vsq)
               Er_particle=derf(mu)-fMv

               h1=1.0_wp/vsq
               h2=vperp2
               h3=vpara**2
               FF=2.0_wp*Er_particle
               GG=(1.0_wp-0.5_wp*h1)*Er_particle+fMv
               HH=Er_particle*h1
               fsp=vpara*FF
               fsr=2.0_wp*h2*FF-h2*HH-(2.0_wp*h3+h2)*GG
               fp=h3*HH+h2*GG
               fr=4.0_wp*h2*(h2*HH+h3*GG)
               fpr=2.0_wp*h2*vpara*(HH-GG)
               rv0=1.0_wp/(sqrt(dgq)*dgq)
               f0dt=dt_coll_e/real(nss,wp)*backne(i)*rv0
                    ! smaller step size is dt_coll/dfloat(nss)
               dvz=-fsp*f0dt+(rn1-0.5_wp)*sqrt(12._wp*fp*f0dt)  ! dvz
               new_vperp2=h2-fsr*f0dt+(rn2-0.5_wp)*                     &
                  sqrt(12.0_wp*(fr-fpr*fpr/fp)*f0dt)+    &
                  (rn1-0.5_wp)*fpr/fp*sqrt(12.0_wp*fp*f0dt)
               do while (new_vperp2 < 0.0_wp)  ! Dangerous loop. Should impose limit...
                  rn1=sprng()
                  rn2=sprng()
                  dvz=-fsp*f0dt+(rn1-0.5_wp)*sqrt(12._wp*fp*f0dt)  ! dvz
                  new_vperp2=h2-fsr*f0dt+(rn2-0.5_wp)*                     &
                     sqrt(12.0_wp*(fr-fpr*fpr/fp)*f0dt)+    &
                     (rn1-0.5_wp)*fpr/fp*sqrt(12.0_wp*fp*f0dt)
               enddo
               vpara=vpara+dvz
               vperp2=new_vperp2
            enddo
         endif
       ! Update particle velocities
       ! rho=vz/B
         Binv=1.0_wp/BB(i)
         zeon(4,i)=vpara*Binv*for_rho               ! rho=v_para/B of electrons
         zeon(6,i)=sqrt(0.5_wp*vperp2*Binv*msm0)      ! sqrt(mu) of electrons

      end do  ! End of loop over electrons

end subroutine CTPandCELL_eon
