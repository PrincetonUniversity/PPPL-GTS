SUBROUTINE CTPandCELL

   use precision
   use global_parameters
   use particle_array
   use field_array
   use esi_interface

   implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

   integer m
   real(wp) da_size, dang, ang, rn1, rn2, d1, Binv
   real(wp) fp, fr, fsr, fsp, fpr, HH, GG, FF, vxy2dum, f0dt, rv0, h1, h2, h3
   integer i, kc, nc, ifail, j, k, jx, iii(mi)

   real(wp) dgr(mi), da(mi), dgq(mi), dgf(mi), mu(mi), Er_particle(mi)   ! working arrays

   !   real(wp), EXTERNAL :: S15AEF   ! The replacement of nag S15AEF DOES NOT do correct at least on Franklin!!!
   real*8 derf

   real(wp), dimension(:), allocatable :: Afsp,Af0dt,Afp,Afsr,Afr,Afpr

   !
   !      C_{TP}: The drag and diffusion operator. Cells are not
   !              necessary, but taking local collision frequency into
   !              account through n(r) and T(r)
   !

   d1=2.0_wp/sqrt(pi)

   !   change to moving frame
   !      vz=vz-backu
   !$omp parallel do private(i)
   do i=1,mi
      dgr(i)=zion(4,i)*BB(i)          ! v_//, BB(m)=B
      dgr(i)=dgr(i)-backu(i)      ! v_// in moving frame
      da(i)=zion(6,i)**2*BB(i)*2.0_wp   ! v_per^2
      dgq(i)=dgr(i)**2+da(i) ! v^2 in moving frame

      !         x2=v2/(2.0d0*backT)            ! v^2/vth^2 in moving frame
      dgf(i)=dgq(i)/backT(i)*0.5_wp  ! v^2/vth^2 in moving frame with vth^2=2T/m

      !         x=sqrt(x2)
      mu(i)=sqrt(dgf(i))

      !         dfi=2.0d0/sqrt(pi)*x*dexp(-x2)
      gha(i)=d1*mu(i)*exp(-dgf(i))
      Er_particle(i)=derf(mu(i))-gha(i)
   enddo

   do i=1,mi
      zion(4,i)=sprng()
   end do

   do i=1,mi
      zion(6,i)=sprng()
   end do

   !$omp parallel do private(i,h1,h2,h3,FF,GG,HH,fsp,fsr,fp,fr,fpr,rv0,f0dt)
   do i=1,mi
      !        h1=1.0d0/x2(i)
      h1=1.0_wp/dgf(i)
      !        h2=vperp2(i)
      h2=da(i)
      !        h3=vz(i)*vz(i)
      h3=dgr(i)**2

      !        FF=2.0d0*fi(i)
      FF=2.0_wp*Er_particle(i)

      !        GG=(1.0d0-0.5d0*h1)*fi(i)+dfi(i)
      GG=(1.0_wp-0.5_wp*h1)*Er_particle(i)+gha(i)

      !        HH=fi(i)*h1
      HH=Er_particle(i)*h1

      !        fsp=vz(i)*FF
      fsp=dgr(i)*FF

      fsr=2.0_wp*h2*FF-h2*HH-(2.0_wp*h3+h2)*GG

      fp=h3*HH+h2*GG

      fr=4.0_wp*h2*(h2*HH+h3*GG)

      !        fpr=2.0d0*h2*vz(i)*(HH-GG)
      fpr=2.0_wp*h2*dgr(i)*(HH-GG)

      !        rv0=vthc/sqrt(v2(i))
      !        rv0=vthc/sqrt(dgq(i))
      !        f0dt=dt_coll*backn(i)/backnc*rv0*rv0*rv0    ! dt_coll=dt/tau_0 where
      !                           ! tau_0 is collision time in terms of backnc and vthc
      rv0=1.0_wp/(sqrt(dgq(i))*dgq(i))
      f0dt=dt_coll*backn(i)*rv0

      ghq(i)=-fsp*f0dt+(zion(4,i)-0.5_wp)*sqrt(12.0_wp*fp*f0dt)    ! dvz

      TT(i)=h2-fsr*f0dt+(zion(6,i)-0.5_wp)*         &
         sqrt(12.0_wp*(fr-fpr*fpr/fp)*f0dt)+     &
         (zion(4,i)-0.5_wp)*fpr/fp*sqrt(12.0_wp*fp*f0dt) ! new v_perp^2

      !  identity: (fr-fpr*fpr/fp)=4.0*h2*HH*GG*(h2+h3)**2/(h3*HH+h2*GG)

   end do

   !  particle sorting in 2D grids

   !      dang=topi/dfloat(np_grid)
   !      da_size=1.0d0/dfloat(n_grid)
   !      mm=0

   !      do i=1,ni
   !        ang=dmod(theta(i),topi)+topi
   !        ang=dmod(ang,topi)
   !        ang=ang/dang

   !        nc=idint(ang)+1
   !        jpol(i)=min(np_grid,max(1,nc))

   !        kc=idint(aa(i)/da_size)+1
   !        jpsi(i)=min(n_grid,max(1,kc))
   !      end do

   ! !MIC$ CASE
   m=0
   !!!!!!     !$omp parallel do private(i,vxy2dum)
   do i=1,mi
      vxy2dum=TT(i)
      if(vxy2dum.ge.0.0_wp) then
         !          vz(i)=vz(i)+dvz(i)
         !          vperp2(i)=vxy2(i)
         dgr(i)=dgr(i)+ghq(i)
         da(i)=TT(i)
      else
         !!!!!!! !$omp atomic
         m=m+1
         iii(m)=i
      endif
   end do

   !  -----------

   !!MIC$ CASE
   !!DIR$ NEXTSCALAR
   !      do i=1,ni
   !         kc=jpsi(i)
   !!         if(kc.lt.(n_grid+1)) then
   !         nc=jpol(i)
   !         mm(nc,kc)=mm(nc,kc)+1
   !         ic(nc,kc,mm(nc,kc))=i
   !!         endif
   !      end do
   !!MIC$ ENDCASE
   !!MIC$ ENDPARALLEL

   !
   !       For particles with vxy<0, re-calculate using
   !       small step size---dt_coll/nss
   !
   if(m.gt.0) then
      allocate(Afsp(m),Af0dt(m),Afp(m),Afsr(m),Afr(m),Afpr(m))
      do 5 k=1,nss

         do i=1,m
            zion(4,i)=sprng()
         end do
         do i=1,m
            zion(6,i)=sprng()
         end do

         do j=1,m
            i=iii(j)
            !        v2(i)=vperp2(i)+vz(i)*vz(i)
            !        x2(i)=v2(i)/(2.0d0*backT(i))
            !        x(i)=sqrt(x2(i))
            !        dfi(i)=2.0d0/dsqrt(pi)*x(i)*dexp(-x2(i))
            dgq(i)=dgr(i)**2+da(i)
            dgf(i)=dgq(i)/backT(i)*0.5_wp
            mu(i)=sqrt(dgf(i))
            gha(i)=d1*mu(i)*exp(-dgf(i))
         end do

         do j=1,m
            i=iii(j)
            ifail=0
            !        fi(i)=S15AEF(x(i),ifail)-dfi(i)
            !        Er_particle(i)=S15AEF(mu(i),ifail)-gha(i)
            ! NOTE!!!!  the replacement of nag S15AEF DOES NOT do correct at least on Franklin!!! Sep. 09
            Er_particle(i)=derf(mu(i))-gha(i)
         end do

         do 4 j=1,m
            i=iii(j)
            !        h1=1.0d0/x2(i)
            h1=1.0_wp/dgf(i)
            !        h2=vperp2(i)
            h2=da(i)
            !        h3=vz(i)*vz(i)
            h3=dgr(i)**2

            !        FF=2.0d0*fi(i)
            FF=2.0_wp*Er_particle(i)
            !        GG=(1.0-0.5*h1)*fi(i)+dfi(i)
            GG=(1.0_wp-0.5_wp*h1)*Er_particle(i)+gha(i)
            !        HH=fi(i)*h1
            HH=Er_particle(i)*h1
            !        Afsp(j)=vz(i)*FF
            Afsp(j)=dgr(i)*FF
            Afsr(j)=2.0_wp*h2*FF-h2*HH-(2.0_wp*h3+h2)*GG
            Afp(j)=h3*HH+h2*GG
            Afr(j)=4.0_wp*h2*(h2*HH+h3*GG)
            !        Afpr(j)=2.0d0*h2*vz(i)*(HH-GG)
            Afpr(j)=2.0_wp*h2*dgr(i)*(HH-GG)
            !        rv0=vthc/sqrt(v2(i))
            !        rv0=vthc/sqrt(dgq(i))
            !        Af0dt(j)=dt_coll/float(nss)*backn(i)/backnc*rv0*rv0*rv0
            !                   ! smaller step size is dt_coll/dfloat(nss)
            rv0=1.0_wp/(sqrt(dgq(i))*dgq(i))
            Af0dt(j)=dt_coll/real(nss,wp)*backn(i)*rv0
            ! smaller step size is dt_coll/dfloat(nss)

            ghq(i)=-Afsp(j)*Af0dt(j)+(zion(4,j)-0.5_wp)*sqrt(12._wp*Afp(j)*Af0dt(j))
            TT(i)=h2-Afsr(j)*Af0dt(j)+(zion(6,j)-0.5_wp)*                     &
               sqrt(12.0_wp*(Afr(j)-Afpr(j)*Afpr(j)/Afp(j))*Af0dt(j))+    &
               (zion(4,j)-0.5_wp)*Afpr(j)/Afp(j)*sqrt(12.0_wp*Afp(j)*Af0dt(j))
            4     continue

            do j=1,m
               i=iii(j)
               11      if(TT(i).ge.0.0_wp) goto 10
               rn1=sprng()
               rn2=sprng()
               ghq(i)=-Afsp(j)*Af0dt(j)+(rn1-0.5_wp)*sqrt(12.0_wp*Afp(j)*Af0dt(j))
               TT(i)=da(i)-Afsr(j)*Af0dt(j)+(rn2-0.5_wp)*                  &
                  sqrt(12.0_wp*(Afr(j)-Afpr(j)*Afpr(j)/Afp(j))*Af0dt(j))+   &
                  (rn1-0.5_wp)*Afpr(j)/Afp(j)*sqrt(12.0_wp*Afp(j)*Af0dt(j))
               goto 11

               10      dgr(i)=dgr(i)+ghq(i)
               da(i)=TT(i)
            end do

            5   continue

            deallocate(Afsp,Af0dt,Afp,Afsr,Afr,Afpr)

         end if

         !$omp parallel do private(i,Binv)
         do i=1,mi
            Binv=1.0/BB(i)
            !     rho=vz/B
            ! moving back to Lab frame after C_TP
            zion(4,i)=(dgr(i)+backu(i))*Binv             ! rho=v_para/B
            zion(6,i)=sqrt(0.5_wp*da(i)*Binv)   ! sqrt(V_perp^2/2B)
         enddo

      end subroutine CTPandCELL
