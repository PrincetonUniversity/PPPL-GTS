subroutine pushion
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use field_array
   use diagnosis_array
   use esi_interface

   implicit none

   !  #define SIMPLE_SPRNG
   !  #define USE_MPI
   !  #include "sprng_f.h"

   !  integer,parameter :: nrbin=5       !wx
   integer i,im,ip,j,jt,kz,m,ii,j11,j10,j01,j00,jtheta,ierror,larmor,ij,ipjt
   real(wp) rdum,dtime,r,grada,vpara,energy,kappa,vdr,wdot,wdrift,&
      vdrenergy,delr,delt(0:mpsi),delz,pi2_inv,adum(0:mpsi),&
      xdum,ydum,xdot,ydot,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,&
      epara,psimax,psimin,pdum,wdrive,ainv,psitmp,thetatmp,&
      zetatmp,rhoi,e1,e2,e3,wpi(3,mi),gaE(mi),wgaE(3,mi),gatE(mi),rdel, &
      dgr(mi),da(mi),dgq(mi),dgf(mi),mu, &
      Er_particle(mi),dela,fip,weqflow,T_inverse,&
      B_inverse,energy1,vdq,wKH,wgradB,omega_t,u_vpara,  &
      v2, avw,   &
      ! following working arrays are defined to use OMP reduction:
   adum_00(0:mpsi),adum_tur(0:mpsi),adum_tot(0:mpsi),adum_11(0:mpsi)
   real(wp) vdgf, vd_para_mag, vd_para_a_cntv, vd_para_gq_cntv, vd_perp_a_cntv, vd_perp_gq_cntv
   real(wp) wdot_vpara, wdot_vperp, wdot_accel, now_w
   real(wp) tot_B_vec_cov(3), tot_B_vec_cntv(3) ! covariant and contravariant components of total B field (B0 + dB) incluidng the magnetic perturbation
   real(wp) grad_lnf0_a_cov, grad_lnf0_gq_cov ! covariant components of grad(ln(f0)) in a and gq axis.
   !(gf component is 0 due to the toroidal axisymmetric assumption )

   !  integer ip_index(mimax)   !wx


   integer gcmotion, gcmotion_em, gcmotion_em_MGY, gcmotion_em_MGY_vec_form

   real(wp) tt0,tt1,cput0,cput1,rtc

   real(wp), external :: D01AHF,dsda,SMOOTHSTEP
   integer minpts,ifail
   real(wp) acc,error_8
   real(wp) E_a,E_t,E_f,timef
   real(wp)  dminv,dmaxv
   real(wp)  wfix     ! =0 in the boundary buffer, and =1 otherwise
   real(wp)  rdm1, rdm2

   real(wp) ddT_da_particle, w_anti_relax, muB, instep,outstep,fdamp,fdamp0, rdum1

   real(wp) wfield(3,mgrid,0:1)   ! weight-field W_0+v_//*W_1+v^2*W_2
   real(wp) FullVelocity_da_dgq_dgf(3,mi) ! ds= a_dot*(Tc/T)**(dsda_power), dgq=gq_dot, dgf=dgf_dot

   delr=1.0_wp/deltar
   delt=1.0_wp/deltat
   delz=1.0_wp/deltaz
   pi2_inv=0.5_wp/pi
   psimax=psi_grid(mpsi)
   psimin=psi_grid(0)

   fdamp0=1.0_wp/tstep*5.0d-4
   if(nonlinear<0.5_wp) fdamp0=0.0_wp
   fdamp0=0.0_wp

   if(irk==1)then
      ! 1st step of Runge-Kutta method
      dtime=0.5_wp*tstep

      ! use (x,y) coordinates near magnetic axis
      !$omp parallel do private(m,xdum,ydum)
      do m=1,mi
         zion0(1:5,m)=zion(1:5,m)
         zion0(7,m)=zion(7,m)
#ifdef _TREAT_MAGNETIC_AXIS
         if(sw(m).eq.1) then
            xdum=zion(1,m)*cos(zion(2,m))
            ydum=zion(1,m)*sin(zion(2,m))
            zion0(1,m)=xdum
            zion0(2,m)=ydum
         end if
#endif
      enddo

      ! 2nd step of Runge-Kutta method
   else
      dtime=tstep
   endif

   !  mu(1:mi)=zion(6,1:mi)**2

   ! gather e_field using 4-point gyro-averaging, sorting in poloidal angle
   ! The following line is an OpenMP directive for loop-level parallelism
   ! on shared memory machines (see http://www.openmp.org).
   !
   ! Ion locations in the 3D field mesh, wzion, wpion, jtion0, jtion1, wtion0, wtion1
   ! are calculated at CHARGE_ION
   !
   !$omp parallel do private(m,e1,e2,e3,wz1,wz0,larmor,ij,wp0,wt00,wt10,wp1,wt01,wt11)
   do m=1,mi
      e1=0.0_wp
      e2=0.0_wp
      e3=0.0_wp
      wz1=wzion(m)
      wz0=1.0_wp-wz1

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


      !MGYOO (10/30/2018; to use analytic profile of ga0(a))
      call Get_MP_alpha_amplitude(zion(1,m), zion(2,m), zion(3,m), istep*tstep, gaE(m), wgaE(1:3,m), gatE(m))



   enddo

   !    prepare equilibrium Er=dPhi/da at particles



   !  ------------------------------------------------------------

   !  apply equilibrium Er in the middle insteady at the beginning

   !  if(istep.lt.2800) then
   if(istep.lt.26500) then
      timef=0.0_wp
   else
      timef=1.0_wp-exp(-(real(istep-26500)/5.0d2)**4)
   end if

   timef=1.0_wp

   Er401_bp=Er401_bp0*timef
   Er_grid=Er_grid0*timef


   !  if(istep.lt.2950) then
   !     timef=0.0_wp
   !  else
   !     timef=1.0_wp-exp(-(real(istep-2950)/5.0d2)**4)
   !  end if
   !  Er401_bp=Er401_bp0*timef
   !  Er_grid=Er_grid0*timef

   !  ------------------------------------------------------------



   dela=1.0_wp/Ydr

   ! if(istep.ge.4000) then
   !$omp parallel do private(m,r,ip,fip)
   do m=1,mi
      r=zion(1,m)-a0
      r=r*dela
      ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
      fip=r-real(ip,wp)
      Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
   enddo
   ! else
   !   Er_particle=0.0
   ! end if

   !  tt0=rtc()
   !  call cpu_time(cput0)

   !!  m=gcmotion(dgr,da,dgq,dgf,zion(4,:),zion(1,:),zion(2,:),mu,wpi(1,:),wpi(2,:),wpi(3,:),Er_particle,nonlinear,mi)

   if(electromagnetic >= 1) then
      !m=gcmotion_em(dgr,da,dgq,dgf,zion,wpi,gaE,wgaE,gatE,Er_particle,nonlinear,mi);
      if(Flag_GC_Solver==0) then
         m=gcmotion_em_MGY(dgr,da,dgq,dgf,zion,wpi,gaE,wgaE,gatE,Er_particle,nonlinear,mi,ziondim);
      elseif(Flag_GC_Solver==1) then
         m=gcmotion_em_MGY_vec_form(dgr,da,dgq,dgf,zion,wpi,gaE,wgaE,gatE,Er_particle,nonlinear,1,mi,ziondim);
      else
         write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
         call flush(stdout)
         stop
      endif
   else
      m=gcmotion(dgr,da,dgq,dgf,zion,wpi,Er_particle,nonlinear,mi,ziondim)
   end if

   !  tt1=rtc()
   !  call cpu_time(cput1)

   !  call n_T_u_group  ! prepare n, T, u_//, dlnn/da, dlnT/da, etc.
   ! moved to collision.f90 -- caution: in this case, BB must be obtained by call esigetb etc.!!!

   !     OUTPUT: dgr=drho/dt;       da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt;
   !             dgf=dzeta/dt;      ghq(:)=da/dt|_EXB;     Ta(:)=dtheta/dt|_EXB;       PP(:)=drho/dt|_EXB;
   !             BB(:)=B;           Ba(:)=dB/da;           Bq(:)=dB/dtheta;
   !             F(:) =F;           Fa(:)=dF/da;           gha(:)=1/|grad a|
   !             no change for zion(1,2,4,:), wpi(1,2,:) and mu;
   !             wpi(3,:)=dphi/dzeta (no change);          Pa(:)=R (major radius)
   !             TT(:)=1/|grad{theta)|; added Aug. 2012 for poloidal momentum flux diagnosis; changed to:
   !             TT(:)=poloidal velocity =vec(v).vec(B_p)/B_p; Oct. 2012

   !
   !        rdum1=maxval(abs(zion(5,1:mi)))
   !        if(rdum1>100.0_wp) print*, 'p1111=',mype, istep, irk, rdum1, mi
   !



   !!!!  if(nonlinear>0.5_wp) then            ! for nonlinear run
   !!!!  if(nonlinear>100.0_wp) then          ! for nonlinear run
   if(wimcompressibility==1) then            ! for implementation of weight imcompressibility
      call weight_field(wfield)      ! calculate weight-field uisng zion0 for irk=1 and zion at half-step for irk=2
      ! used BB obtained from gcmotion, and calculte its own ion locations
      ! which are used to calculate wfield in phase space grids and gether wfield at ions (next)

      !
      ! rdum1=maxval(abs(wfield))
      ! if(rdum1>100.0_wp) then
      !    print*, 'p8888=', mype, istep, irk, rdum1, maxval(abs(wfield(1,1:mgrid,1))) , &
      !                 maxval(abs(wfield(2,1:mgrid,1))) , maxval(abs(wfield(3,1:mgrid,1))) , maxval(abs(wfield(1,1:mgrid,0))) , &
      !                 maxval(abs(wfield(2,1:mgrid,0))) , maxval(abs(wfield(3,1:mgrid,0)))
      ! end if
      !

   end if




   dminv=0.0d0
   dmaxv=0.0d0


   ! update GC position
   ! Another OpenMP parallel loop...
   !$omp parallel do private(m,mu,vpara,energy,kappa,wdrift,weqflow,vdr,wdrive,wdot,xdum,ydum,T_inverse,energy1,B_inverse,u_vpara,vdq,wgradB,wKH, &
   !$omp& omega_t,r,ip,fip,ddT_da_particle,w_anti_relax,muB,instep,outstep,fdamp,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,v2,avw)
   do m=1,mi
      !     dbdz=0.0   for axisymmetric system   !wx

      mu=zion(6,m)*zion(6,m)    ! =v_perp^2/2B
      muB=mu*BB(m)              ! =v_perp^2/2
      T_inverse=1.0_wp/backT(m)
      B_inverse=1.0_wp/BB(m)
      vpara=zion(4,m)*BB(m)     ! BB(m)=B
      omega_t=backu(m)
      backu(m)=backu(m)*F(m)*B_inverse         ! =u_//=I*omega_t/B
      u_vpara=backu(m)/vpara                   ! =(u_//)/(v_//)
      !     energy=0.5_wp*(vpara-backu(m))*(vpara-backu(m))+mu*BB(m)
      energy=0.5_wp*(vpara-backu(m))*(vpara-backu(m))+muB
      energy1=energy-backu(m)*backu(m)*0.5_wp  ! =(v_//-u)^2/2+mu*B-u^2/2 -- see note

      v2=vpara*vpara+muB*2.0_wp        ! v^2

      kappa=(min(umax*umax,energy1*T_inverse)-1.5_wp)*(-dlnTda(m))+(-dlnnda(m))+gh(m) !~-dlnf_SM/da
      !gh=dln<A>/da with A=exp(u^2/2T)
      r=zion(7,m)-s0
      r=r*delr
      ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]; using coordinate s to find particle radial location in grids
      fip=(zion(1,m)-psi_grid(ip))/(psi_grid(ip+1)-psi_grid(ip))
      ddT_da_particle=fip*ddT_da(ip+1)+(1.0_wp-fip)*ddT_da(ip) ! d<delta_Ti>/da at particle's position

      if(ifneoclassical==1) then
         vdr=da(m)              !  full da/at      -- including NC dynamics
         vdq=dgq(m)             !  full dtheta/dt  -- including NC dynamics
         vdgf=dgf(m)            !  full dzeta/dt   -- including NC dynamics
      else
         !   ExB drift in radial direction for w-dot and flux diagnostics
         vdr=ghq(m)             !  da/at|_EXB
         vdq=Ta(m)              !  dtheta/dt|_EXB
         vdgf=dvi_gf(m)
      endif

      wdrive=vdr*kappa

      weqflow=vdr*Ba(m)+vdq*Bq(m)

      ! new term added to maintain constant grad(Ti) drive
      w_anti_relax=(min(umax*umax,energy*T_inverse)-1.5_wp)*ddT_da_particle*T_inverse
      w_anti_relax=vdr*w_anti_relax

      if(zion(1,m).le.lmax.or.zion(1,m).ge.rmin) then
         instep =smoothstep(lmin,lmax,zion(1,m))
         outstep=smoothstep(rmin,rmax,zion(1,m))
         !        fip=smoothstep(0.20_wp,0.25_wp,zion(1,m))*(1.0_wp-smoothstep(0.75_wp,0.80_wp,zion(1,m)))
         fip=instep*(1.0_wp-outstep)
         !  wdrive=wdrive*fip
         w_anti_relax=w_anti_relax*fip
         fdamp=(1.0_wp-fip)*fdamp0
      else
         fip=1.0_wp
         fdamp=0.0_wp
      end if

      if(ifneoclassical==1) then
         wdrift=(dgr(m)*BB(m)+zion(4,m)*weqflow)*(vpara-backu(m))*T_inverse   !  full drho/dt -- including NC dynamics
      else
         wdrift=(PP(m)*BB(m)+zion(4,m)*weqflow)*(vpara-backu(m))*T_inverse   ! see note
      endif

      wgradB=weqflow*T_inverse*mu       ! interesting term propotional to grad_B

      wKH=-(vdr*(duda(m)*F(m)+Fa(m)*omega_t-Ba(m)*backu(m)) &
         -vdq*Bq(m)*backu(m))*B_inverse*vpara*T_inverse       ! major Kelvin-Helmholtz drive


      !
      !      calculate average weight <w> at particles using wieght-field wfield
      !
      !!!!!  if(nonlinear>0.5_wp) then          ! for nonlinear run
      !!!!!  if(nonlinear>100.0_wp) then        ! for nonlinear run
      if(wimcompressibility==1) then          ! for weight imcompressibility
         wz1=wzion(m)
         wz0=1.0_wp-wz1

         ij=jtion0(0,m)
         wp0=1.0_wp-wpion(0,m)

         wt00=1.0_wp-wtion0(0,m)
         e1=wp0*wt00*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
         e2=wp0*wt00*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
         e3=wp0*wt00*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

         ij=ij+1
         wt10=1.0_wp-wt00
         e1=e1+wp0*wt10*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
         e2=e2+wp0*wt10*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
         e3=e3+wp0*wt10*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

         ij=jtion1(0,m)
         wp1=1.0_wp-wp0
         wt01=1.0_wp-wtion1(0,m)
         e1=e1+wp1*wt01*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
         e2=e2+wp1*wt01*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
         e3=e3+wp1*wt01*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

         ij=ij+1
         wt11=1.0_wp-wt01
         e1=e1+wp1*wt11*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
         e2=e2+wp1*wt11*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
         e3=e3+wp1*wt11*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

         avw=e1+e2*vpara+e3*v2


         !
         !      if(abs(avw)>100.0_wp) then
         !         print*, 'p9999=', mype, istep, irk, avw, e1, e2, e3, vpara, v2, maxval(abs(wfield(1,1:mgrid,1))) , &
         !                 maxval(abs(wfield(2,1:mgrid,1))) , maxval(abs(wfield(3,1:mgrid,1))) , maxval(abs(wfield(1,1:mgrid,0))) , &
         !                 maxval(abs(wfield(2,1:mgrid,0))) , maxval(abs(wfield(3,1:mgrid,0)))
         !      end if
         !


      end if


      !     if(mype==100.and.m==mi/2) print*, 'AVW=', avw, zion(5,m), mype
      !     if(mype==0.and.m==mi/4) print*, 'AVW=', avw, zion(5,m), mype
      !     if(mype==200.and.m==mi/8) print*, 'AVW=', avw, zion(5,m), mype

      !     wdot=(zion0(8,m)-nonlinear*zion(5,m))*(wdrive+wdrift+wKH+wgradB+w_anti_relax)

      ! if(wimcompressibility==1) then     ! for weight imcompressibility
      !    !  wdot=(zion0(8,m)-nonlinear*avw)*(wdrive*fip+wdrift+wKH+wgradB)   ! Correct statement for weight_damping
      !    wdot=(zion0(8,m)-nonlinear*avw)*(wdrive+wdrift+wKH+wgradB)*fip   ! Correct statement for weight_damping
      ! else
      !    wdot=(zion0(8,m)-nonlinear*zion(5,m))*(wdrive+wdrift+wKH+wgradB)*fip
      ! end if

      wdot = (wdrive+wdrift+wKH+wgradB)



      !==================================================
      !MGYOO(01/14/2019)
      !Note that f0 is toroidally symmetric so that grad_lnf0_gf=0 and we do not need vd_para_gf or vd_perp_gf terms to calculate weight evolution
      !Caution!! Following routine is only available for neoclassical case which use full velocity.
      if(irk==2)then
         if(ifneoclassical==1) then
            call esiget_Intp_B_vec_cov_and_cntv(zion(1,m), zion(2,m), zion(3,m), istep*tstep, tot_B_vec_cov, tot_B_vec_cntv)
            vd_para_mag = (vdr*tot_B_vec_cov(1) + vdq*tot_B_vec_cov(2) + vdgf*tot_B_vec_cov(3))/BB(m)
            vd_para_a_cntv = vd_para_mag*(tot_B_vec_cntv(1)/BB(m))
            vd_para_gq_cntv = vd_para_mag*(tot_B_vec_cntv(2)/BB(m))
         else
            vd_para_a_cntv = vpara_dvi_a(m)
            vd_para_gq_cntv = vpara_dvi_gq(m)
         endif
         vd_perp_a_cntv = vdr -vd_para_a_cntv
         vd_perp_gq_cntv = vdq -vd_para_gq_cntv

         ! grad_lnf0_a_cov = -kappa - T_inverse*mu*Ba(m) + (duda(m)*F(m)+Fa(m)*omega_t-Ba(m)*backu(m))*B_inverse*vpara*T_inverse
         ! grad_lnf0_gq_cov = -T_inverse*mu*Bq(m) - Bq(m)*backu(m)*B_inverse*vpara*T_inverse

         ! wdot_vpara = -vd_para_a_cntv*grad_lnf0_a_cov -vd_para_gq_cntv*grad_lnf0_gq_cov
         ! wdot_vperp = -vd_perp_a_cntv*grad_lnf0_a_cov -vd_perp_gq_cntv*grad_lnf0_gq_cov
         ! wdot_accel  = wdrift

         !=================================================================
         !MGYOO (01/30/2019)
         ! grad_lnf0 includes (wdrive + wKH) terms
         ! wdot_accel includes (wdrift + wgradB) terms
         grad_lnf0_a_cov = -kappa + (duda(m)*F(m)+Fa(m)*omega_t-Ba(m)*backu(m))*B_inverse*vpara*T_inverse
         grad_lnf0_gq_cov = - Bq(m)*backu(m)*B_inverse*vpara*T_inverse

         wdot_vpara = -vd_para_a_cntv*grad_lnf0_a_cov -vd_para_gq_cntv*grad_lnf0_gq_cov
         wdot_vperp = -vd_perp_a_cntv*grad_lnf0_a_cov -vd_perp_gq_cntv*grad_lnf0_gq_cov
         wdot_accel  = wdrift + wgradB
         !=================================================================

         if(wimcompressibility==1) then
            now_w = avw
         else
            now_w = zion(5,m)
         endif

         wi_vpara(m) = wi_vpara(m) + dtime*( (zion0(8,m)-nonlinear*now_w)*wdot_vpara*fip)
         wi_vperp(m) = wi_vperp(m) + dtime*( (zion0(8,m)-nonlinear*now_w)*wdot_vperp*fip)
         wi_accel(m) = wi_accel(m) + dtime*( (zion0(8,m)-nonlinear*now_w)*wdot_accel*fip)
         wi_fdamp(m) = wi_fdamp(m) + dtime*(-fdamp*zion(5,m) )

      endif
      !=================================================




      !      dminv=min(dminv,(wdrive+wdrift+wKH+wgradB)*tstep)
      !      dmaxv=max(dmaxv,(wdrive+wdrift+wKH+wgradB)*tstep)





      ! if(mype==0.and.m==mi/2) print*, 'GGGGG', wdrive, w_anti_relax

      zion(1,m) = zion0(1,m)+dtime*da(m)
      zion(2,m) = zion0(2,m)+dtime*dgq(m)
      zion(3,m) = zion0(3,m)+dtime*dgf(m)
      zion(4,m) = zion0(4,m)+dtime*dgr(m)

      !  if(zion0(1,m).gt.psi_grid(mpsi-nbound-1))  then
      !     wfix=1.0_wp
      !  else
      !     wfix=1.0_wp
      !  endif

      if(Flag_w_Evol==1) then
         ! From the definition of the transformed weight z=-ln(gaE-w), where gaE=zion0(8,m)
         ! the weight increase becomes delta_w = (gaE-w0)*[1-exp(-F*dt)], where F is the driving force
         ! NOTE This scheme is beneficial for nonlinear case only (dw/dt=(gaE-w)*F).
         ! For linear case ( dw/dt=F), this transform may cause additional numerical errors and accumulation.
         zion(5,m) = zion0(5,m)+(zion0(8,m)-zion0(5,m))*(1.0-exp(-wdot*dtime))
      else
         if(wimcompressibility==1) then     ! for weight imcompressibility
            zion(5,m) = zion0(5,m)+dtime*((zion0(8,m)-nonlinear*avw)*wdot*fip-fdamp*zion(5,m))
         else
            zion(5,m) = zion0(5,m)+dtime*((zion0(8,m)-nonlinear*zion(5,m))*wdot*fip-fdamp*zion(5,m))
         end if
      endif

      !     zion(5,m) = zion0(5,m)+dtime*(wdot-fdamp*zion(5,m))*wfix
      !     zion(7,m) = zion0(7,m)+dtime*sqrt(Tc*T_inverse)*da(m)
      !     zion(7,m) = zion0(7,m)+dtime*(Tc*T_inverse)**0.4_wp*da(m)
      zion(7,m) = zion0(7,m)+dtime*(Tc*T_inverse)**dsda_power*da(m)

      !  if (mype ==0 .AND. MOD(m,(mi/10))==0) then
      !   write(stdout,*) 'I''m= ',mype, ' zion(5,m=',m,')= ',zion(5,m), ' zion0(5,m=',m,')= ',zion0(5,m), ' avw= ',avw
      !   write(stdout,*) 'zion(6,m=',m,')= ',zion(6,m), ' zion0(8,m=',m,')= ',zion0(8,m) , ' wdot= ', wdot
      ! end if


#ifdef _TREAT_MAGNETIC_AXIS
      if(sw(m).eq.1) then    ! convert to (a, theta) coordinates
         xdum=sqrt(zion(1,m)**2+zion(2,m)**2)
         ydum=sign(1.0_wp,zion(2,m))*acos(zion(1,m)/xdum)
         zion(1,m)=xdum
         zion(2,m)=ydum
      end if
#endif

      ! theta and zeta normalize to [0,2*pi), modulo is slower than hand coded procedure
      !     zion(2,m)=zion(2,m)*pi2_inv+100.0 !period of 1
      !     zion(2,m)=2.0*pi*(zion(2,m)-aint(zion(2,m))) ![0,2*pi)
      !     zion(3,m)=zion(3,m)*pi2_inv+100.0
      !     zion(3,m)=2.0*pi*(zion(3,m)-aint(zion(3,m)))

      ! use f90 modulo to normalize theta and zeta to [0,2*pi)
      ! If P > 0, then the function modulo(A,P) determines the positive difference between A and then
      ! next lowest multiple of P. If P < 0, then the function modulo(A,P) determines the negative
      ! difference between A and then next highest multiple of P.
      zion(2,m)=modulo(zion(2,m),twopi)
      zion(3,m)=modulo(zion(3,m),twopi)

      ! store GC information for flux measurements
      wpi(1,m)=ghq(m)    ! at present, wpi(1,*)=da/dt|_EXB
      wpi(2,m)=energy    ! =(v_para-u//)^2/2+v_perp^2/2
      wpi(3,m)=muB       ! =v_perp^2/2, need for calculate distribution function in get_fluxes

      if(irk==2) then
         UnNorm_theta_i(m) = UnNorm_theta_i(m) +dtime*dgq(m)
         UnNorm_zeta_i(m)  = UnNorm_zeta_i(m)  +dtime*dgf(m)
      endif
   enddo


   if(irk==2)then
      ! out of boundary particle
      !$omp parallel do private(m)
      do m=1,mi
         if(zion(1,m) > psimax)then
            !  rdm1=sprng()
            !  rdm1=dsign(1.0d0,(rdm1-0.5d0))
            !  rdm2=sprng()
            !  rdm2=dsign(1.0d0,(rdm2-0.5d0))
            zion(1,m)=zion0(1,m)
            zion(2,m)=2.0_wp*pi-zion0(2,m)
            zion(3,m)=zion0(3,m)
            zion(4,m)=zion0(4,m)
            !  zion(4,m)=zion0(4,m)*rdm1  !  for test
            zion(5,m)=zion0(5,m)  ! this sounds reasonable: energy, momentum and particle conserving boundary
            !  zion(5,m)=zion0(5,m)*rdm2  !  for test
            !  zion(5,m)=-zion0(5,m)  ! remove density (weight) accumulation at edge when ifneoclassical=1
            !  zion(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
            zion(7,m)=zion0(7,m)

         elseif(zion(1,m) < psimin)then
            zion(1,m)=zion0(1,m)
            zion(2,m)=2.0_wp*pi-zion0(2,m)
            zion(3,m)=zion0(3,m)
            zion(4,m)=zion0(4,m)
            zion(5,m)=zion0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
            !  zion(5,m)=-zion0(5,m)  ! remove density (weight) accumulation at edge when ifneoclassical=1
            !  zion(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
            zion(7,m)=zion0(7,m)
         endif
      enddo

      ! ********************************
      !if(mype==2) then
      !!  error_8=1.0d-5
      !!  IFAIL=1
      !!  rdum=D01AHF(0.0d0,zion(1,mi/2),error_8,minpts,acc,dsda,-1,IFAIL)
      !!  print*, '##################', zion(7,mi/2),'  ', rdum
      !   write(stdout,*)'######## gcmotion: wall =',(tt1-tt0),'  cpu =',(cput1-cput0)
      !end if
      ! ********************************

   endif





   !  if(mype==0) print*, 'IONS:       Max Damping=', dmaxv, '   Min Damping=', dminv
   !  if(mype==0) print*, 'IONS:       ', 'Damping=', dmaxv, '   Min Damping=', dminv




   if(idiag==0)then
      ! fluxes diagnose
      !    wpi(1,1:mi)=wpi(1,1:mi)*gha(1:mi)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
      !$omp parallel do private(m)
      do m=1,mi
         wpi(1,m)=wpi(1,m)*gha(m)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
      end do

      rmarker=0.0_wp
      eflux=0.0_wp
      pfluxi=0.0_wp
      dflowi=0.0_wp
      entropyi=0.0_wp
      fluxpsi=0.0_wp
      heatflux=0.0_wp
      En_00=0.0_wp
      En_tur=0.0_wp
      En_tot=0.0_wp
      do m=1,mi
         r=zion0(7,m)         ! old radial location
         ip=max(0,min(mpsi,int((r-s0)*delr+0.5_wp)))    ! ip [0,mpsi]
         vdrenergy=wpi(1,m)*(wpi(2,m)-2.5*vth_grid(ip)**2)*zion0(5,m)
         ! note: wpi(1,m)=dr/dt (EXB part)

         ! radial location in diagnostic bin
         i=max(1,min(mflux,1+int((r-s0)*real(mflux)/(s1-s0))))

         ! marker,energy,particle,momentum fluxes,parallel flows,entropy and kinetic energy
         rmarker(i)=rmarker(i)+zion0(8,m)
         eflux(i)=eflux(i)+vdrenergy
         pfluxi=pfluxi+wpi(1,m)*zion0(5,m)
         dflowi=dflowi+BB(m)*zion0(4,m)*zion0(5,m)
         entropyi=entropyi+zion0(5,m)*zion0(5,m)
      enddo

      if(ZF_SEPARATION==0) then

         !     note: zonal flow component <phi> and phip00=d<phi>/da are calculated in field.f90

         adum=0.0_wp
         adum_00=0.0_wp
         adum_tur=0.0_wp
         adum_tot=0.0_wp

         !$omp parallel private(i,j,ij,E_a,E_t,E_f) reduction(+:adum,adum_00,adum_tur,adum_tot)
         ! Note from OMP compiler: A variable in the REDUCTION clause or any component or element thereof must not have the POINTER
         ! attribute; nor be a pointee, an allocatable object, an assumed-size array, or a THREADLOCAL common block variable.
         !$omp do
         do i=0,mpsi
            do j=1,mtheta(i)
               ij=igrid(i)+j
               adum(i)=adum(i)+dtemper(ij)*(phi(1,ij)-phi00(i))**2  ! surface averaging for delta-phi

               adum_00(i)=adum_00(i)+dtemper(ij)*phip00(i)**2*gij(1,ij)  ! surface average

               E_a=evector(1,1,ij)+Er_grid(i)           ! total radial component   ! d\phi/da
               E_t=evector(2,1,ij)                      ! d\delta\phi/d\theta
               E_f=evector(3,1,ij)                      ! d\delta\phi/d\zeta
               adum_tot(i)=adum_tot(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                  +E_t**2*gij(2,ij)  &
                  +E_f**2*gij(3,ij)  &
                  +2.0_wp*E_a*E_t*gij(4,ij))

               E_a=E_a-phip00(i)-Er_grid(i) ! turbulence component only    ! d\delta\phi/da
               adum_tur(i)=adum_tur(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                  +E_t**2*gij(2,ij)  &
                  +E_f**2*gij(3,ij)  &
                  +2.0_wp*E_a*E_t*gij(4,ij))
            enddo
         enddo
         !$omp end parallel

      else  ! ZF_SEPARATION=1

         adum=0.0_wp
         adum_00=0.0_wp
         adum_tur=0.0_wp
         adum_tot=0.0_wp
         adum_11=0.0_wp

         ! DEBUG
         !        if(mype==100) print*, 'Dirichlet2', phi(1,igrid(0)), phi(1,igrid(0)+mtheta(0)/2),phi(1,igrid(0)+mtheta(0)-1), &
         !                               phi(1,igrid(mpsi)), phi(1,igrid(mpsi)+mtheta(mpsi)/2),phi(1,igrid(mpsi)+mtheta(mpsi)-1)

         !$omp parallel private(i,j,ij,E_a,E_t,E_f) reduction(+:adum,adum_00,adum_tur,adum_tot,adum_11)
         !$omp do
         do i=0,mpsi
            do j=1,mtheta(i)
               ij=igrid(i)+j
               !          fluxpsi(i)=fluxpsi(i)+phi(1,ij)*phi(1,ij)               ! what average?
               adum(i)=adum(i)+dtemper(ij)*phi(1,ij)*phi(1,ij)  ! surface averaging for <delta-phi^2>

               adum_11(i)=adum_11(i)+dtemper(ij)*phi(1,ij)      ! to check surface averaged <delta-phi> which should be zero by defination

               adum_00(i)=adum_00(i)+dtemper(ij)*phip00(i)**2*gij(1,ij)  ! surface average

               E_a=evector(1,1,ij)+Er_grid(i)           ! total radial component   ! d\phi/da
               E_t=evector(2,1,ij)                      ! d\delta\phi/d\theta
               E_f=evector(3,1,ij)                      ! d\delta\phi/d\zeta
               adum_tot(i)=adum_tot(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                  +E_t**2*gij(2,ij)  &
                  +E_f**2*gij(3,ij)  &
                  +2.0*E_a*E_t*gij(4,ij))

               E_a=E_a-phip00(i)-Er_grid(i) ! turbulence component only    ! d\delta\phi/da
               adum_tur(i)=adum_tur(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                  +E_t**2*gij(2,ij)  &
                  +E_f**2*gij(3,ij)  &
                  +2.0_wp*E_a*E_t*gij(4,ij))
            enddo
         enddo
         !$omp end parallel
      end if   !  end if(ZF_SEPARATION==0)

      call MPI_ALLREDUCE(adum_11,flux_11,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      flux_11=flux_11/jaco_sum/Tc   ! =<phi>/Tc -- added in March 2013 to check whether <phi>=0 is obtained

      call MPI_ALLREDUCE(adum,fluxpsi,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      fluxpsi=fluxpsi/jaco_sum/(Tc*Tc)  ! =<(phi/Tc)^2> -- turbulence field intensity
      sumfluxpsi=sumfluxpsi+fluxpsi     ! =int[I(a)dt]

      call MPI_ALLREDUCE(adum_00,En_00,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      En_00=En_00/jaco_sum      ! =<E_00^2>    -- zonal flow energy

      call MPI_ALLREDUCE(adum_tur,En_tur,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      En_tur=En_tur/jaco_sum     ! =<\deltaE^2> -- turbulence energy

      call MPI_ALLREDUCE(adum_tot,En_tot,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      En_tot=En_tot/jaco_sum     ! =<E^2> -- total field energy

      call get_fluxes(wpi,dgf,dgq)

      !   im=mod(istep,100)
      im=mod(istep,ndiag)
      if(im==0) then
         call get_3D_flow(da,dgq,dgf,wpi(2,:))
      end if

   endif



end subroutine pushion



subroutine pushion_RK4(irk_i)
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use field_array
   use diagnosis_array
   use esi_interface
   use function_gg_module
   use Ptl_Src_Loss
   use Diagnostic_Grid

   implicit none
   integer, intent(in):: irk_i
   integer i,im,ip,j,jt,kz,m,ii,j11,j10,j01,j00,jtheta,ierror,larmor,ij,ipjt
   real(wp) rdum,r,grada,vpara,energy,kappa,vdr,wdot,wdrift,&
      vdrenergy,delr,delt(0:mpsi),delz,pi2_inv,adum(0:mpsi),&
      xdum,ydum,xdot,ydot,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,&
      epara,psimax,psimin,pdum,wdrive,ainv,psitmp,thetatmp,&
      zetatmp,rhoi,e1,e2,e3,wpi(3,mi),gaE(mi),wgaE(3,mi),gatE(mi),rdel, &
      dgr(mi),da(mi),dgq(mi),dgf(mi),mu, &
      Er_particle(mi),dela,fip,weqflow,T_inverse,&
      B_inverse,energy1,vdq,wKH,wgradB,omega_t,u_vpara,  &
      v2, avw,   &
      ! following working arrays are defined to use OMP reduction:
   adum_00(0:mpsi),adum_tur(0:mpsi),adum_tot(0:mpsi),adum_11(0:mpsi)
   real(wp) vdgf, vd_para_mag, vd_para_a_cntv, vd_para_gq_cntv, vd_perp_a_cntv, vd_perp_gq_cntv
   real(wp) wdot_vpara, wdot_vperp, wdot_accel, now_w
   real(wp) tot_B_vec_cov(3), tot_B_vec_cntv(3) ! covariant and contravariant components of total B field (B0 + dB) incluidng the magnetic perturbation
   real(wp) grad_lnf0_a_cov, grad_lnf0_gq_cov ! covariant components of grad(ln(f0)) in a and gq axis.
   !(gf component is 0 due to the toroidal axisymmetric assumption )


   integer gcmotion, gcmotion_em, gcmotion_em_MGY, gcmotion_em_MGY_vec_form

   real(wp) tt0,tt1,cput0,cput1,rtc

   real(wp), external :: D01AHF,dsda !,SMOOTHSTEP
   integer minpts,ifail
   real(wp) acc,error_8
   real(wp) E_a,E_t,E_f,timef
   real(wp)  dminv,dmaxv
   real(wp)  wfix     ! =0 in the boundary buffer, and =1 otherwise
   real(wp)  rdm1, rdm2

   real(wp) ddT_da_particle, w_anti_relax, muB, instep,outstep,fdamp,fdamp0, rdum1

   real(wp) wfield(3,mgrid,0:1)   ! weight-field W_0+v_//*W_1+v^2*W_2
   real(wp) RK_cum_coeff, RK_adv_coeff ! Coefficients to cumlate (to get future final value) or advance(to predict next RK step) RK4 terms
   real(wp) FullVelocity_da_dgq_dgf(3,mi) ! ds= a_dot*(Tc/T)**(dsda_power), dgq=gq_dot, dgf=dgf_dot

   !$acc declare create(dgr,da,dgq,dgf,wpi,Er_particle,gaE,wgaE,gatE)

   if(acc_newer_H_zion) then
      !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
      acc_newer_H_zion=.false.
      acc_newer_D_zion=.false.
   endif

   if(istep==istart .and. irk_i==1) then
      time_in_Pushi =0
      timewc_in_Pushi =0
   endif

   delr=1.0_wp/deltar
   delt=1.0_wp/deltat
   delz=1.0_wp/deltaz
   pi2_inv=0.5_wp/pi
   psimax=psi_grid(mpsi)
   psimin=psi_grid(0)

   fdamp0=1.0_wp/tstep*5.0d-4
   if(nonlinear<0.5_wp) fdamp0=0.0_wp
   fdamp0=0.0_wp

   if(irk_i==1)then ! 1st step of Runge-Kutta method
      RK_cum_coeff =(1.0/6.0)
      RK_adv_coeff =0.5
      acc_newer_D_zion=.true. !device side zion will be changed in the following acc routine
      !$acc parallel loop present(zion,zion0,RK_cum_zion)
      do m=1,mi
         ! Note that zion0(8,m) does not mean "mu", this is reserved for the other quantity
         !$acc loop seq
         do i=1,5
            zion0(i,m)=zion(i,m)
         enddo

         zion0(7,m)=zion(7,m)

         !$acc loop seq
         do i=1,size(RK_cum_zion,1)
            RK_cum_zion(i,m)=0.0
         enddo
      enddo
   elseif(irk_i==2) then
      RK_cum_coeff =(1.0/3.0)
      RK_adv_coeff =0.5
   elseif(irk_i==3) then
      RK_cum_coeff =(1.0/3.0)
      RK_adv_coeff =1.0
   elseif(irk_i==4) then
      RK_cum_coeff =(1.0/6.0)
      RK_adv_coeff =1.0
   else
      write(stdout,*) "In pushi_RK4, irk_i=",irk_i,", is weird. Please Check It!!"
      call flush(stdout)
   endif

   ! Followings are not required here, so commented out. (MGYOO: 11/13/2019)
   call ion_locator_at_larmor_points ! This is required to get correct wpi values which would use ion location info. such as jtion, wpion, wtion.
   ! call n_T_u_group  ! prepare n, T, dlnn/da, dlnT/da, etc.


   !  mu(1:mi)=zion(6,1:mi)**2

   ! gather e_field using 4-point gyro-averaging, sorting in poloidal angle
   ! The following line is an OpenMP directive for loop-level parallelism
   ! on shared memory machines (see http://www.openmp.org).
   !
   ! Ion locations in the 3D field mesh, wzion, wpion, jtion0, jtion1, wtion0, wtion1
   ! are calculated at CHARGE_ION
   !
   !$omp parallel do private(m,e1,e2,e3,wz1,wz0,larmor,ij,wp0,wt00,wt10,wp1,wt01,wt11)
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


#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   !MGYOO (10/30/2018; to use analytic profile of ga0(a))
   !$acc parallel loop present(zion,gaE,wgaE,gatE)
   do m=1,mi
      call Get_MP_alpha_amplitude(zion(1,m),zion(2,m),zion(3,m),istep*tstep, gaE(m), wgaE(:,m), gatE(m))
   enddo
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion, Get_MP_alpha_amplitude", tt1-tt0," secs"
   call flush(stdout)
#endif

   !    prepare equilibrium Er=dPhi/da at particles
   !  ------------------------------------------------------------
   !  apply equilibrium Er in the middle insteady at the beginning
   !  if(istep.lt.2800) then
   if(istep.lt.26500) then
      timef=0.0_wp
   else
      timef=1.0_wp-exp(-(real(istep-26500)/5.0d2)**4)
   end if

   timef=1.0_wp

   Er401_bp=Er401_bp0*timef
   Er_grid=Er_grid0*timef
   dela=1.0_wp/Ydr
   ! if(istep.ge.4000) then
   !$omp parallel do private(m,r,ip,fip)
   !$acc parallel loop present(zion,zion0)
   do m=1,mi
      r=zion(1,m)-a0
      r=r*dela
      ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
      fip=r-real(ip,wp)
      Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
   enddo
   !  ------------------------------------------------------------

#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   if(electromagnetic >= 1) then
      !m=gcmotion_em(dgr,da,dgq,dgf,zion,wpi,gaE,wgaE,gatE,Er_particle,nonlinear,mi);
      if(Flag_GC_Solver==0) then
         m=gcmotion_em_MGY(dgr,da,dgq,dgf,zion,wpi,gaE,wgaE,gatE,Er_particle,nonlinear,mi,ziondim);
      elseif(Flag_GC_Solver==1) then
         m=gcmotion_em_MGY_vec_form(dgr,da,dgq,dgf,zion,wpi,gaE,wgaE,gatE,Er_particle,nonlinear,1,mi,ziondim);
      else
         write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
         call flush(stdout)
         stop
      endif
   else
      m=gcmotion(dgr,da,dgq,dgf,zion,wpi,Er_particle,nonlinear,mi,ziondim)
   end if

#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion, ion gc motion takes", tt1-tt0," secs"
   call flush(stdout)
#endif

   !  call n_T_u_group  ! prepare n, T, u_//, dlnn/da, dlnT/da, etc.
   ! moved to collision.f90 -- caution: in this case, BB must be obtained by call esigetb etc.!!!

   !     OUTPUT: dgr=drho/dt;       da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt;
   !             dgf=dzeta/dt;      ghq(:)=da/dt|_EXB;     Ta(:)=dtheta/dt|_EXB;       PP(:)=drho/dt|_EXB;
   !             BB(:)=B;           Ba(:)=dB/da;           Bq(:)=dB/dtheta;
   !             F(:) =F;           Fa(:)=dF/da;           gha(:)=1/|grad a|
   !             no change for zion(1,2,4,:), wpi(1,2,:) and mu;
   !             wpi(3,:)=dphi/dzeta (no change);          Pa(:)=R (major radius)
   !             TT(:)=1/|grad{theta)|; added Aug. 2012 for poloidal momentum flux diagnosis; changed to:
   !             TT(:)=poloidal velocity =vec(v).vec(B_p)/B_p; Oct. 2012

   !
   !        rdum1=maxval(abs(zion(5,1:mi)))
   !        if(rdum1>100.0_wp) print*, 'p1111=',mype, istep, irk_i, rdum1, mi
   !



   !!!!  if(nonlinear>0.5_wp) then            ! for nonlinear run
   !!!!  if(nonlinear>100.0_wp) then          ! for nonlinear run
   if(wimcompressibility==1) then            ! for implementation of weight imcompressibility
      call weight_field(wfield)       ! calculate weight-field uisng zion0 for irk_i=1 and zion at half-step for irk_i=2
      ! used BB obtained from gcmotion, and calculte its own ion locations
      ! which are used to calculate wfield in phase space grids and gether wfield at ions (next)
   end if

   dminv=0.0d0
   dmaxv=0.0d0


#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif

   ! update GC position
   ! Another OpenMP parallel loop...
   !$omp parallel do private(m,mu,vpara,energy,kappa,wdrift,weqflow,vdr,wdrive,wdot,xdum,ydum,T_inverse,energy1,B_inverse,u_vpara,vdq,wgradB,wKH, &
   !$omp& omega_t,r,ip,fip,ddT_da_particle,w_anti_relax,muB,instep,outstep,fdamp,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,v2,avw)
   acc_newer_D_zion=.true. !device side zion will be changed in the following acc routine
   !$acc parallel loop present(zion,zion0,RK_cum_zion)
   do m=1,mi
      !     dbdz=0.0   for axisymmetric system   !wx

      mu=zion(6,m)*zion(6,m)    ! =v_perp^2/2B
      muB=mu*BB(m)              ! =v_perp^2/2
      T_inverse=1.0_wp/backT(m)
      B_inverse=1.0_wp/BB(m)
      vpara=zion(4,m)*BB(m)     ! BB(m)=B
      omega_t=backu(m)
      backu(m)=backu(m)*F(m)*B_inverse         ! =u_//=I*omega_t/B
      u_vpara=backu(m)/vpara                   ! =(u_//)/(v_//)
      !     energy=0.5_wp*(vpara-backu(m))*(vpara-backu(m))+mu*BB(m)
      energy=0.5_wp*(vpara-backu(m))*(vpara-backu(m))+muB
      energy1=energy-backu(m)*backu(m)*0.5_wp  ! =(v_//-u)^2/2+mu*B-u^2/2 -- see note

      v2=vpara*vpara+muB*2.0_wp        ! v^2

      kappa=(min(umax*umax,energy1*T_inverse)-1.5_wp)*(-dlnTda(m))+(-dlnnda(m))+gh(m) !~-dlnf_SM/da
      !gh=dln<A>/da with A=exp(u^2/2T)
      r=zion(7,m)-s0
      r=r*delr
      ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]; using coordinate s to find particle radial location in grids
      fip=(zion(1,m)-psi_grid(ip))/(psi_grid(ip+1)-psi_grid(ip))
      ddT_da_particle=fip*ddT_da(ip+1)+(1.0_wp-fip)*ddT_da(ip) ! d<delta_Ti>/da at particle's position

      if(ifneoclassical==1) then
         vdr=da(m)              !  full da/at      -- including NC dynamics
         vdq=dgq(m)             !  full dtheta/dt  -- including NC dynamics
         vdgf=dgf(m)            !  full dzeta/dt   -- including NC dynamics
      else
         !   ExB drift in radial direction for w-dot and flux diagnostics
         vdr=ghq(m)             !  da/at|_EXB
         vdq=Ta(m)              !  dtheta/dt|_EXB
         vdgf=dvi_gf(m)
      endif

      wdrive=vdr*kappa

      weqflow=vdr*Ba(m)+vdq*Bq(m)

      ! new term added to maintain constant grad(Ti) drive
      w_anti_relax=(min(umax*umax,energy*T_inverse)-1.5_wp)*ddT_da_particle*T_inverse
      w_anti_relax=vdr*w_anti_relax

      if(zion(1,m).le.lmax.or.zion(1,m).ge.rmin) then
         instep =smoothstep_m(lmin,lmax,zion(1,m))
         outstep=smoothstep_m(rmin,rmax,zion(1,m))
         !        fip=smoothstep(0.20_wp,0.25_wp,zion(1,m))*(1.0_wp-smoothstep(0.75_wp,0.80_wp,zion(1,m)))
         fip=instep*(1.0_wp-outstep)
         !  wdrive=wdrive*fip
         w_anti_relax=w_anti_relax*fip
         fdamp=(1.0_wp-fip)*fdamp0
      else
         fip=1.0_wp
         fdamp=0.0_wp
      end if

      if(ifneoclassical==1) then
         wdrift=(dgr(m)*BB(m)+zion(4,m)*weqflow)*(vpara-backu(m))*T_inverse   !  full drho/dt -- including NC dynamics
      else
         wdrift=(PP(m)*BB(m)+zion(4,m)*weqflow)*(vpara-backu(m))*T_inverse   ! see note
      endif

      wgradB=weqflow*T_inverse*mu       ! interesting term propotional to grad_B

      wKH=-(vdr*(duda(m)*F(m)+Fa(m)*omega_t-Ba(m)*backu(m)) &
         -vdq*Bq(m)*backu(m))*B_inverse*vpara*T_inverse       ! major Kelvin-Helmholtz drive


      !
      !      calculate average weight <w> at particles using wieght-field wfield
      !
      !!!!!  if(nonlinear>0.5_wp) then          ! for nonlinear run
      !!!!!  if(nonlinear>100.0_wp) then        ! for nonlinear run
      if(wimcompressibility==1) then          ! for weight imcompressibility
         wz1=wzion(m)
         wz0=1.0_wp-wz1

         ij=jtion0(0,m)
         wp0=1.0_wp-wpion(0,m)

         wt00=1.0_wp-wtion0(0,m)
         e1=wp0*wt00*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
         e2=wp0*wt00*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
         e3=wp0*wt00*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

         ij=ij+1
         wt10=1.0_wp-wt00
         e1=e1+wp0*wt10*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
         e2=e2+wp0*wt10*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
         e3=e3+wp0*wt10*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

         ij=jtion1(0,m)
         wp1=1.0_wp-wp0
         wt01=1.0_wp-wtion1(0,m)
         e1=e1+wp1*wt01*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
         e2=e2+wp1*wt01*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
         e3=e3+wp1*wt01*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

         ij=ij+1
         wt11=1.0_wp-wt01
         e1=e1+wp1*wt11*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
         e2=e2+wp1*wt11*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
         e3=e3+wp1*wt11*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

         avw=e1+e2*vpara+e3*v2
      end if

      if(wimcompressibility==1) then
         now_w = avw
      else
         now_w = zion(5,m)
      endif

      wdot = (wdrive+wdrift+wKH+wgradB)

      ! if(ieee_is_nan(wdot)) then
      !    write(stdout,*) 'In pushiRK4, irk=',irk,', m=',m,', zion=',zion(:,m)
      !    write(stdout,*) 'In pushiRK4, irk=',irk,', m=',m,', zion0=',zion0(:,m)
      !    write(stdout,*) 'In pushiRK4, irk=',irk,', m=',m,', RK_cum_zion=',RK_cum_zion(:,m)
      !    write(stdout,*) 'In pushiRK4, irk=',irk,', m=',m,', now_w=', now_w,', wdot=',wdot,', wdrive=',wdrive,', wdrift=',wdrift,', wKH=',wKH,', wgradB=',wgradB
      !    call flush(stdout)
      ! endif


      !==================================================
      ! MGYOO(03/27/2019)
      ! Note that the following routine should be upgraded to calculate correct RK4
      !==================================================
      !MGYOO(01/14/2019)
      !Note that f0 is toroidally symmetric so that grad_lnf0_gf=0 and we do not need vd_para_gf or vd_perp_gf terms to calculate weight evolution
      !Caution!! Following routine is only available for neoclassical case which use full velocity.
      if(irk_i==4)then
         if(ifneoclassical==1) then
            ! call esiget_Intp_B_vec_cov_and_cntv(zion(1,m), zion(2,m), zion(3,m), istep*tstep, tot_B_vec_cov, tot_B_vec_cntv)
            ! vd_para_mag = (vdr*tot_B_vec_cov(1) + vdq*tot_B_vec_cov(2) + vdgf*tot_B_vec_cov(3))/BB(m)
            ! vd_para_a_cntv = vd_para_mag*(tot_B_vec_cntv(1)/BB(m))
            ! vd_para_gq_cntv = vd_para_mag*(tot_B_vec_cntv(2)/BB(m))
         else
            vd_para_a_cntv = vpara_dvi_a(m)
            vd_para_gq_cntv = vpara_dvi_gq(m)
         endif
         vd_perp_a_cntv = vdr -vd_para_a_cntv
         vd_perp_gq_cntv = vdq -vd_para_gq_cntv


         !=================================================================
         !MGYOO (01/30/2019)
         ! grad_lnf0 includes (wdrive + wKH) terms
         ! wdot_accel includes (wdrift + wgradB) terms
         grad_lnf0_a_cov = -kappa + (duda(m)*F(m)+Fa(m)*omega_t-Ba(m)*backu(m))*B_inverse*vpara*T_inverse
         grad_lnf0_gq_cov = - Bq(m)*backu(m)*B_inverse*vpara*T_inverse

         wdot_vpara = -vd_para_a_cntv*grad_lnf0_a_cov -vd_para_gq_cntv*grad_lnf0_gq_cov
         wdot_vperp = -vd_perp_a_cntv*grad_lnf0_a_cov -vd_perp_gq_cntv*grad_lnf0_gq_cov
         wdot_accel  = wdrift + wgradB
         !=================================================================

         wi_vpara(m) = wi_vpara(m) + tstep*( (zion0(8,m)-nonlinear*now_w)*wdot_vpara*fip)
         wi_vperp(m) = wi_vperp(m) + tstep*( (zion0(8,m)-nonlinear*now_w)*wdot_vperp*fip)
         wi_accel(m) = wi_accel(m) + tstep*( (zion0(8,m)-nonlinear*now_w)*wdot_accel*fip)
         wi_fdamp(m) = wi_fdamp(m) + tstep*(-fdamp*zion(5,m) )

      endif
      !=================================================



      ! Accumulate k-th RK terms (k1,k2,k3,k4) for the future value: f_(n+1) = f_(n) + (1/6)*(k1+2*k2+2*k3+k4)
      RK_cum_zion(1,m) = RK_cum_zion(1,m) +RK_cum_coeff*tstep*da(m)
      RK_cum_zion(2,m) = RK_cum_zion(2,m) +RK_cum_coeff*tstep*dgq(m)
      RK_cum_zion(3,m) = RK_cum_zion(3,m) +RK_cum_coeff*tstep*dgf(m)
      RK_cum_zion(4,m) = RK_cum_zion(4,m) +RK_cum_coeff*tstep*dgr(m)
      RK_cum_zion(7,m) = RK_cum_zion(7,m) +RK_cum_coeff*tstep*(Tc*T_inverse)**dsda_power*da(m)

      if(Flag_w_Evol==1) then
         RK_cum_zion(5,m) = RK_cum_zion(5,m) +(wdot)*(RK_cum_coeff*tstep)
      else
         RK_cum_zion(5,m) = RK_cum_zion(5,m) +RK_cum_coeff*tstep*((zion0(8,m)-nonlinear*now_w)*wdot*fip-fdamp*zion(5,m))
      endif

      ! Advance particle information for the next RK step
      if(irk_i<4) then

         zion(1,m) = zion0(1,m)+RK_adv_coeff*tstep*da(m)
         zion(2,m) = zion0(2,m)+RK_adv_coeff*tstep*dgq(m)
         zion(3,m) = zion0(3,m)+RK_adv_coeff*tstep*dgf(m)
         zion(4,m) = zion0(4,m)+RK_adv_coeff*tstep*dgr(m)
         zion(7,m) = zion0(7,m)+RK_adv_coeff*tstep*(Tc*T_inverse)**dsda_power*da(m)

         if(Flag_w_Evol==1) then
            zion(5,m) = zion0(5,m)+(zion0(8,m)-zion0(5,m))*(1.0-exp(-wdot*RK_adv_coeff*tstep))
         else
            zion(5,m) = zion0(5,m)+RK_adv_coeff*tstep*((zion0(8,m)-nonlinear*now_w)*wdot*fip-fdamp*zion(5,m))
         endif

      else if(irk_i==4) then !irk_i==4
         ! Get the future value at time level (n+1): f_(n+1) = f_(n) + (1/6)*(k1+2*k2+2*k3+k4) = f_(n)+RK_cum_value
         !$acc loop seq
         do i=1,4
            zion(i,m) = zion0(i,m) +RK_cum_zion(i,m)
         enddo
         zion(7,m) = zion0(7,m) +RK_cum_zion(7,m)

         if(Flag_w_Evol==1) then
            ! Get the future value at time level (n+1): using accumulated increase of the transformed weight Z.
            zion(5,m) = zion0(5,m)+(zion0(8,m)-zion0(5,m))*(1.0-exp(-RK_cum_zion(5,m)))
         else
            ! Get the future value at time level (n+1): f_(n+1) = f_(n) + (1/6)*(k1+2*k2+2*k3+k4) = f_(n)+RK_cum_value
            zion(5,m) = zion0(5,m) +RK_cum_zion(5,m)
         endif
      endif



      ! if(ieee_is_nan(zion(5,m))) then
      !    write(stdout,*) 'In pushiRK4, irk=',irk,', m=',m,', zion=',zion(:,m)
      !    write(stdout,*) 'In pushiRK4, irk=',irk,', m=',m,', zion0=',zion0(:,m)
      !    write(stdout,*) 'In pushiRK4, irk=',irk,', m=',m,', RK_cum_zion=',RK_cum_zion(:,m)
      !    write(stdout,*) 'In pushiRK4, irk=',irk,', m=',m,', zion0(8,m)=',zion0(8,m),', now_w=', now_w, ', wdot=',wdot
      !    call flush(stdout)
      ! endif



      ! use f90 modulo to normalize theta and zeta to [0,2*pi)
      ! If P > 0, then the function modulo(A,P) determines the positive difference between A and then
      ! next lowest multiple of P. If P < 0, then the function modulo(A,P) determines the negative
      ! difference between A and then next highest multiple of P.
      zion(2,m)=modulo(zion(2,m),twopi)
      zion(3,m)=modulo(zion(3,m),twopi)

      ! store GC information for flux measurements
      wpi(1,m)=ghq(m)    ! at present, wpi(1,*)=da/dt|_EXB
      wpi(2,m)=energy    ! =(v_para-u//)^2/2+v_perp^2/2
      wpi(3,m)=muB       ! =v_perp^2/2, need for calculate distribution function in get_fluxes

      if(irk_i==4) then
         UnNorm_theta_i(m) = UnNorm_theta_i(m) +RK_cum_zion(2,m)
         UnNorm_zeta_i(m)  = UnNorm_zeta_i(m)  +RK_cum_zion(3,m)
      endif
   enddo ! do m=1,mi



   if(irk_i==4)then
      if(Flag_Measure_Ptl_FootPrint==1) then
      !$acc update self(wpi(:,:mi))
      call Measure_and_Cumulate_Ion_FootPrint_at_BDY1(FP_BDY1,wpi(2,:))
      endif
      call treat_ions_out_of_boundary
      call Cumulate_i_lost_num
   endif

#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion,  ACC 'main do loop in pushion_RK4' takes", tt1-tt0," secs"
   call flush(stdout)
#endif

   if(idiag==0)then

      call cpu_time(cput0)
      tt0=MPI_WTIME()

      !!!!$acc update self(Fa(:mi),PP(:mi),Bq(:mi),Ba(:mi)) ! these are not used for now, but possibly for future diagnostics
      !$acc update self(BB(:mi),TT(:mi),Ta(:mi),ghq(:mi),gha(:mi),Pa(:mi),F(:mi))
      !$acc update self(ii_inv_grad_a(:mi),ii_inv_grad_gq(:mi),ii_inv_grad_gf(:mi))
      !$acc update self(dgr(:mi),da(:mi),dgq(:mi),dgf(:mi),wpi(:,:mi))
      !$acc update self (dvi_gf(:mi),vpara_dvi_a(:mi),vpara_dvi_gq(:mi),vpara_dvi_gf(:mi))
      !$acc update self(UnNorm_theta_i(:mi),UnNorm_zeta_i(:mi))
      !$acc update self(wi_init(:mi),wi_vpara(:mi),wi_vperp(:mi),wi_accel(:mi),wi_fdamp(:mi))
      !$acc update self(backn(:mi),backT(:mi),backu(:mi),dlnnda(:mi),dlnTda(:mi),duda(:mi),gh(:mi))
      !$acc update self(back_f0_i(:mi))

      if(acc_newer_D_zion) then
         !$acc update self(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
         acc_newer_D_zion=.false.
         acc_newer_H_zion=.false.
      endif

      ! fluxes diagnose
      !    wpi(1,1:mi)=wpi(1,1:mi)*gha(1:mi)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
      !$omp parallel do private(m)
      do m=1,mi
         wpi(1,m)=wpi(1,m)*gha(m)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
      end do

      rmarker=0.0_wp
      eflux=0.0_wp
      pfluxi=0.0_wp
      dflowi=0.0_wp
      entropyi=0.0_wp
      fluxpsi=0.0_wp
      heatflux=0.0_wp
      En_00=0.0_wp
      En_tur=0.0_wp
      En_tot=0.0_wp
      !$acc parallel loop present(zion,zion0)
      do m=1,mi
         r=zion0(7,m)         ! old radial location
         ip=max(0,min(mpsi,int((r-s0)*delr+0.5_wp)))    ! ip [0,mpsi]
         vdrenergy=wpi(1,m)*(wpi(2,m)-2.5*vth_grid(ip)**2)*zion0(5,m)
         ! note: wpi(1,m)=dr/dt (EXB part)

         ! radial location in diagnostic bin
         i=max(1,min(mflux,1+int((r-s0)*real(mflux)/(s1-s0))))

         ! marker,energy,particle,momentum fluxes,parallel flows,entropy and kinetic energy
         rmarker(i)=rmarker(i)+zion0(8,m)
         eflux(i)=eflux(i)+vdrenergy
         pfluxi=pfluxi+wpi(1,m)*zion0(5,m)
         dflowi=dflowi+BB(m)*zion0(4,m)*zion0(5,m)
         entropyi=entropyi+zion0(5,m)*zion0(5,m)
      enddo


      if(ZF_SEPARATION==0) then
         !     note: zonal flow component <phi> and phip00=d<phi>/da are calculated in field.f90
         adum=0.0_wp
         adum_00=0.0_wp
         adum_tur=0.0_wp
         adum_tot=0.0_wp

         !$omp parallel private(i,j,ij,E_a,E_t,E_f) reduction(+:adum,adum_00,adum_tur,adum_tot)
         ! Note from OMP compiler: A variable in the REDUCTION clause or any component or element thereof must not have the POINTER
         ! attribute; nor be a pointee, an allocatable object, an assumed-size array, or a THREADLOCAL common block variable.
         !$omp do
         do i=0,mpsi
            do j=1,mtheta(i)
               ij=igrid(i)+j
               adum(i)=adum(i)+dtemper(ij)*(phi(1,ij)-phi00(i))**2  ! surface averaging for delta-phi

               adum_00(i)=adum_00(i)+dtemper(ij)*phip00(i)**2*gij(1,ij)  ! surface average

               E_a=evector(1,1,ij)+Er_grid(i)           ! total radial component   ! d\phi/da
               E_t=evector(2,1,ij)                      ! d\delta\phi/d\theta
               E_f=evector(3,1,ij)                      ! d\delta\phi/d\zeta
               adum_tot(i)=adum_tot(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                  +E_t**2*gij(2,ij)  &
                  +E_f**2*gij(3,ij)  &
                  +2.0_wp*E_a*E_t*gij(4,ij))

               E_a=E_a-phip00(i)-Er_grid(i) ! turbulence component only    ! d\delta\phi/da
               adum_tur(i)=adum_tur(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                  +E_t**2*gij(2,ij)  &
                  +E_f**2*gij(3,ij)  &
                  +2.0_wp*E_a*E_t*gij(4,ij))
            enddo
         enddo
         !$omp end parallel

      else  ! ZF_SEPARATION=1

         adum=0.0_wp
         adum_00=0.0_wp
         adum_tur=0.0_wp
         adum_tot=0.0_wp
         adum_11=0.0_wp

         ! DEBUG
         !        if(mype==100) print*, 'Dirichlet2', phi(1,igrid(0)), phi(1,igrid(0)+mtheta(0)/2),phi(1,igrid(0)+mtheta(0)-1), &
         !                               phi(1,igrid(mpsi)), phi(1,igrid(mpsi)+mtheta(mpsi)/2),phi(1,igrid(mpsi)+mtheta(mpsi)-1)

         !$omp parallel private(i,j,ij,E_a,E_t,E_f) reduction(+:adum,adum_00,adum_tur,adum_tot,adum_11)
         !$omp do
         do i=0,mpsi
            do j=1,mtheta(i)
               ij=igrid(i)+j
               !          fluxpsi(i)=fluxpsi(i)+phi(1,ij)*phi(1,ij)               ! what average?
               adum(i)=adum(i)+dtemper(ij)*phi(1,ij)*phi(1,ij)  ! surface averaging for <delta-phi^2>

               adum_11(i)=adum_11(i)+dtemper(ij)*phi(1,ij)      ! to check surface averaged <delta-phi> which should be zero by defination

               adum_00(i)=adum_00(i)+dtemper(ij)*phip00(i)**2*gij(1,ij)  ! surface average

               E_a=evector(1,1,ij)+Er_grid(i)           ! total radial component   ! d\phi/da
               E_t=evector(2,1,ij)                      ! d\delta\phi/d\theta
               E_f=evector(3,1,ij)                      ! d\delta\phi/d\zeta
               adum_tot(i)=adum_tot(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                  +E_t**2*gij(2,ij)  &
                  +E_f**2*gij(3,ij)  &
                  +2.0*E_a*E_t*gij(4,ij))

               E_a=E_a-phip00(i)-Er_grid(i) ! turbulence component only    ! d\delta\phi/da
               adum_tur(i)=adum_tur(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                  +E_t**2*gij(2,ij)  &
                  +E_f**2*gij(3,ij)  &
                  +2.0_wp*E_a*E_t*gij(4,ij))
            enddo
         enddo
         !$omp end parallel
      end if   !  end if(ZF_SEPARATION==0)

      call MPI_ALLREDUCE(adum_11,flux_11,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      flux_11=flux_11/jaco_sum/Tc   ! =<phi>/Tc -- added in March 2013 to check whether <phi>=0 is obtained

      call MPI_ALLREDUCE(adum,fluxpsi,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      fluxpsi=fluxpsi/jaco_sum/(Tc*Tc)  ! =<(phi/Tc)^2> -- turbulence field intensity
      sumfluxpsi=sumfluxpsi+fluxpsi     ! =int[I(a)dt]

      call MPI_ALLREDUCE(adum_00,En_00,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      En_00=En_00/jaco_sum      ! =<E_00^2>    -- zonal flow energy

      call MPI_ALLREDUCE(adum_tur,En_tur,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      En_tur=En_tur/jaco_sum     ! =<\deltaE^2> -- turbulence energy

      call MPI_ALLREDUCE(adum_tot,En_tot,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      En_tot=En_tot/jaco_sum     ! =<E^2> -- total field energy

      call get_fluxes(wpi,dgf,dgq)

      !   im=mod(istep,100)
      im=mod(istep,ndiag)
      if(im==0) then
         call get_3D_flow(da,dgq,dgf,wpi(2,:))
      end if

      call cpu_time(cput1)
      tt1=MPI_WTIME()
      time_in_Pushi(1) = time_in_Pushi(1) + (cput1-cput0)
      timewc_in_Pushi(1) = timewc_in_Pushi(1) + (tt1-tt0)
   endif

   if(irk_i==1) then
      call cpu_time(cput0)
      tt0=MPI_WTIME()

      call Measure_and_Cumulate_iZonData_on_DG(wpi,da)

      call cpu_time(cput1)
      tt1=MPI_WTIME()
      time_in_Pushi(1) = time_in_Pushi(1) + (cput1-cput0)
      timewc_in_Pushi(1) = timewc_in_Pushi(1) + (tt1-tt0)

#ifdef _MGY_TIMER_
      if(mype==0) then
         write(stdout,*) "Measure_Cumulate_iFlux_a1D takes", tt1-tt0," secs"
         call flush(stdout)
      endif
#endif
  endif

   !Another Diagnostics on DG (Diagnostic Grid)
   if(irk_i==1 .and. mod(istep,msnap)==0) then
      call cpu_time(cput0)
      tt0=MPI_WTIME()

      !!!!$acc update self(Fa(:mi),PP(:mi),Bq(:mi),Ba(:mi)) ! these are not used for now, but possibly for future diagnostics
      !$acc update self(BB(:mi),TT(:mi),Ta(:mi),ghq(:mi),gha(:mi),Pa(:mi),F(:mi))
      !$acc update self(dgr(:mi),da(:mi),dgq(:mi),dgf(:mi),wpi(:,:mi))
      !$acc update self(ii_inv_grad_a(:mi),ii_inv_grad_gq(:mi),ii_inv_grad_gf(:mi))
      !$acc update self (dvi_gf(:mi),vpara_dvi_a(:mi),vpara_dvi_gq(:mi),vpara_dvi_gf(:mi))
      !$acc update self(UnNorm_theta_i(:mi),UnNorm_zeta_i(:mi))
      !$acc update self(wi_init(:mi),wi_vpara(:mi),wi_vperp(:mi),wi_accel(:mi),wi_fdamp(:mi))
      !$acc update self(ilabels(:mi))
      !$acc update self(backn(:mi),backT(:mi),backu(:mi),dlnnda(:mi),dlnTda(:mi),duda(:mi),gh(:mi))
      !$acc update self(back_f0_i(:mi))

      if(acc_newer_D_zion) then
         !$acc update self(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
         acc_newer_D_zion=.false.
         acc_newer_H_zion=.false.
      endif

      do m=1,mi
         FullVelocity_da_dgq_dgf(1,m) = da(m)
         FullVelocity_da_dgq_dgf(2,m) = dgq(m)
         FullVelocity_da_dgq_dgf(3,m) = dgf(m)
      enddo
      call Measure_Ions_on_Diagnostic_Grid(wpi(1:3,:),FullVelocity_da_dgq_dgf)

      call cpu_time(cput1)
      tt1=MPI_WTIME()
      time_in_Pushi(1) = time_in_Pushi(1) + (cput1-cput0)
      timewc_in_Pushi(1) = timewc_in_Pushi(1) + (tt1-tt0)
#ifdef _MGY_TIMER_
      if(mype==0) then
         write(stdout,*) "Measure_Ions_on_Diagnostic_Grid takes", tt1-tt0," secs"
         call flush(stdout)
      endif
#endif
   endif



end subroutine pushion_RK4



subroutine pushion_wEvol3(irk_i,nrk)
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use field_array
   use diagnosis_array
   use esi_interface
   use function_gg_module
   use Ptl_Src_Loss
   use Diagnostic_Grid
   use Structured_Grid
   implicit none
   integer, intent(in):: irk_i,nrk
   integer i,ip,m,ij,ierr,ierror,j
   real(wp) RK_cum_coeff, RK_adv_coeff ! Coefficients to cumulate (to get future final value) or advance(to predict next RK step)
   real(wp) dgr(mi),da(mi),dgq(mi),dgf(mi),gradPhi(3,mi),wpi(3,mi),Er_particle(mi)
   real(wp) gaE(mi),wgaE(3,mi), gatE(mi) ! Magnetic Perturbation Information: alpha, grad_alpha, partial derivative d(gaE)/dt
   real(wp) adum(0:mpsi),adum_00(0:mpsi),adum_tur(0:mpsi),adum_tot(0:mpsi),adum_11(0:mpsi)
   real(wp) tt0,tt1,cput0,cput1 ! For timer
   real(wp) timef,dela,r,fip
   integer gcmotion, gcmotion_em, gcmotion_em_MGY, gcmotion_em_MGY_vec_form
   real(wp) mu,vpara,energy
   real(wp) e1,e2,e3,Ea_00
   real(wp) wz0,wz1,wp0,wt00,wt10,wp1,wt01,wt11
   real(wp) E_a,E_t,E_f,delr
   real(wp) FullVelocity_da_dgq_dgf(3,mi) ! ds= a_dot*(Tc/T)**(dsda_power), dgq=gq_dot, dgf=dgf_dot
   real(wp) vdrenergy,new_f0
   integer larmor, loop

   !$acc declare create(dgr,da,dgq,dgf,gradPhi,wpi,Er_particle,gaE,wgaE,gatE)

   if(acc_newer_H_zion) then
      !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
      acc_newer_H_zion=.false.
      acc_newer_D_zion=.false.
   endif

   if(istep==istart .and. irk_i==1) then
      time_in_Pushi =0
      timewc_in_Pushi =0
   endif


   delr=1.0_wp/deltar


   ! Set RK_cum_coeff and RK_adv_coeff
   if(nrk==4) then
      select case (irk_i)
         case(1)
            RK_adv_coeff = 0.5
            RK_cum_coeff = (1.0/6.0)
         case(2)
            RK_adv_coeff = 0.5
            RK_cum_coeff = (1.0/3.0)
         case(3)
            RK_adv_coeff = 1.0
            RK_cum_coeff = (1.0/3.0)
         case(4)
            RK_adv_coeff = 1.0
            RK_cum_coeff = (1.0/6.0)
         case default
            write(stdout,'((A,I3,A,I3))') "Error!! irk_i=",irk_i,"!! irk_i must be between 1 and nrk=",nrk
            call flush(stdout)
            STOP
      end select
   elseif(nrk==2) then
      select case (irk_i)
         case(1)
            RK_adv_coeff = 0.5
            RK_cum_coeff = 0.0
         case(2)
            RK_cum_coeff = 1.0
            RK_adv_coeff = 1.0
         case default
            write(stdout,'((A,I3,A,I3))') "Error!! irk_i=",irk_i,"!! irk_i must be between 1 and nrk=",nrk
            call flush(stdout)
            STOP
      end select
   endif


   if(irk_i==1)then ! 1st step of Runge-Kutta method
      !$acc parallel loop present(zion,zion0,RK_cum_zion)
      do m=1,mi
         !$acc loop seq
         do i=1,5
            zion0(i,m)=zion(i,m)
         enddo
         zion0(7,m)=zion(7,m)

         !$acc loop seq
         do i=1,size(RK_cum_zion,1)
            RK_cum_zion(i,m)=0.0
         enddo
      enddo
      acc_newer_D_zion=.true. !device side zion will be changed in the following acc routine
   endif

   ! Followings are not required here, so commented out. (MGYOO: 11/13/2019)
   call ion_locator_at_larmor_points ! This is required to get correct wpi values which would use ion location info. such as jtion, wpion, wtion.

   if(irk_i==1) then
      !$acc parallel loop present(zion,back_f0_i)
      do m=1,mi
         zion(8,m) = zion0(6,m)*back_f0_i(m) ! zion(8,m) is temporary old_f0 @ irk=1 step
      enddo
   endif

   ! gather e_field using multi-point gyro-averaging, sorting in poloidal angle
   ! The following line is an OpenMP directive for loop-level parallelism
   ! on shared memory machines (see http://www.openmp.org).

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

      gradPhi(1,m)=e1*inv_larmor_num
      gradPhi(2,m)=e2*inv_larmor_num
      gradPhi(3,m)=e3*inv_larmor_num

   enddo


#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   !MGYOO (10/30/2018; to use analytic profile of ga0(a))
   !$acc parallel loop present(zion,gaE,wgaE,gatE)
   do m=1,mi
      call Get_MP_alpha_amplitude(zion(1,m),zion(2,m),zion(3,m),istep*tstep, gaE(m), wgaE(:,m), gatE(m))
   enddo
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion, Get_MP_alpha_amplitude", tt1-tt0," secs"
   call flush(stdout)
#endif


   timef=1.0_wp
   Er401_bp=Er401_bp0*timef
   Er_grid=Er_grid0*timef
   dela=1.0_wp/Ydr
   !$acc parallel loop present(zion,zion0)
   do m=1,mi
      r=zion(1,m)-a0
      r=r*dela
      ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
      fip=r-real(ip,wp)
      Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
   enddo
   !  ------------------------------------------------------------

#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   if(electromagnetic >= 1) then
      if(Flag_GC_Solver==0) then
         m=gcmotion_em_MGY(dgr,da,dgq,dgf,zion,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,mi,ziondim);
      elseif(Flag_GC_Solver==1) then
         m=gcmotion_em_MGY_vec_form(dgr,da,dgq,dgf,zion,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,1,mi,ziondim);
      else
         write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
         call flush(stdout)
         stop
      endif
   else
      m=gcmotion(dgr,da,dgq,dgf,zion,gradPhi,Er_particle,nonlinear,mi,ziondim)
   end if
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion, ion gc motion takes", tt1-tt0," secs"
   call flush(stdout)
#endif


   !     OUTPUT: dgr=drho/dt;       da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt;
   !             dgf=dzeta/dt;      ghq(:)=da/dt|_EXB;     Ta(:)=dtheta/dt|_EXB;       PP(:)=drho/dt|_EXB;
   !             BB(:)=B;           Ba(:)=dB/da;           Bq(:)=dB/dtheta;
   !             F(:) =F;           Fa(:)=dF/da;           gha(:)=1/|grad a|
   !             no change for zion(1,2,4,:), wpi(1,2,:) and mu;
   !             wpi(3,:)=dphi/dzeta (no change);          Pa(:)=R (major radius)
   !             TT(:)=1/|grad{theta)|; added Aug. 2012 for poloidal momentum flux diagnosis; changed to:
   !             TT(:)=poloidal velocity =vec(v).vec(B_p)/B_p; Oct. 2012




#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif

   !$acc parallel loop present(zion,zion0,RK_cum_zion)
   do m=1,mi
      mu=zion(6,m)**2.0 ! normalized mu of electrons
      vpara = BB(m)*zion(4,m)
      energy = 0.5*vpara**2.0+mu*BB(m)

      ! Accumulate k-th RK terms (k1,k2,k3,k4) for the future value: f_(n+1) = f_(n) + (1/6)*(k1+2*k2+2*k3+k4)
      RK_cum_zion(1,m) = RK_cum_zion(1,m) +RK_cum_coeff*tstep*da(m)
      RK_cum_zion(2,m) = RK_cum_zion(2,m) +RK_cum_coeff*tstep*dgq(m)
      RK_cum_zion(3,m) = RK_cum_zion(3,m) +RK_cum_coeff*tstep*dgf(m)
      RK_cum_zion(4,m) = RK_cum_zion(4,m) +RK_cum_coeff*tstep*dgr(m)
      RK_cum_zion(7,m) = RK_cum_zion(7,m) +RK_cum_coeff*tstep*(Tc/backT(m))**dsda_power*da(m)

      ! Advance particle information for the next RK step
      if(irk_i<nrk) then
         zion(1,m) = zion0(1,m)+RK_adv_coeff*tstep*da(m)
         zion(2,m) = zion0(2,m)+RK_adv_coeff*tstep*dgq(m)
         zion(3,m) = zion0(3,m)+RK_adv_coeff*tstep*dgf(m)
         zion(4,m) = zion0(4,m)+RK_adv_coeff*tstep*dgr(m)
         zion(7,m) = zion0(7,m)+RK_adv_coeff*tstep*(Tc/backT(m))**dsda_power*da(m)
      else if(irk_i==nrk) then !irk_i==4
         !$acc loop seq
         do i=1,4
            zion(i,m) = zion0(i,m) +RK_cum_zion(i,m)
         enddo
         zion(7,m) = zion0(7,m) +RK_cum_zion(7,m)
      endif

      zion(2,m)=modulo(zion(2,m),twopi)
      zion(3,m)=modulo(zion(3,m),twopi)

      ! store GC information for flux measurements
      wpi(1,m)=ghq(m)    ! at present, wpi(1,*)=da/dt|_EXB
      wpi(2,m)=energy    ! =(v_para-u//)^2/2+v_perp^2/2
      wpi(3,m)=mu*BB(m)  ! =v_perp^2/2, need for calculate distribution function in get_fluxes

      if(irk_i==nrk) then
         UnNorm_theta_i(m) = UnNorm_theta_i(m) +RK_cum_zion(2,m)
         UnNorm_zeta_i(m)  = UnNorm_zeta_i(m)  +RK_cum_zion(3,m)
      endif
   enddo ! do m=1,mi
   acc_newer_D_zion=.true. !device side zion will be changed in the following acc routine


   if(irk_i==nrk)then
      if(Flag_Measure_Ptl_FootPrint==1) then
      !$acc update self(wpi(:,:mi))
      call Measure_and_Cumulate_Ion_FootPrint_at_BDY1(FP_BDY1,wpi(2,:))
      endif
      call treat_ions_out_of_boundary
      call Cumulate_i_lost_num
   endif

#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion,  ACC 'main do loop in pushion_RK4' takes", tt1-tt0," secs"
   call flush(stdout)
#endif

   if(idiag==0)then

      call cpu_time(cput0)
      tt0=MPI_WTIME()

      !!!!$acc update self(Fa(:mi),PP(:mi),Bq(:mi),Ba(:mi)) ! these are not used for now, but possibly for future diagnostics
      !$acc update self(BB(:mi),TT(:mi),Ta(:mi),ghq(:mi),gha(:mi),Pa(:mi),F(:mi))
      !$acc update self(ii_inv_grad_a(:mi),ii_inv_grad_gq(:mi),ii_inv_grad_gf(:mi))
      !$acc update self(dgr(:mi),da(:mi),dgq(:mi),dgf(:mi),wpi(:,:mi))
      !$acc update self (dvi_gf(:mi),vpara_dvi_a(:mi),vpara_dvi_gq(:mi),vpara_dvi_gf(:mi))
      !$acc update self(UnNorm_theta_i(:mi),UnNorm_zeta_i(:mi))
      !$acc update self(wi_init(:mi),wi_vpara(:mi),wi_vperp(:mi),wi_accel(:mi),wi_fdamp(:mi))
      !$acc update self(backn(:mi),backT(:mi),backu(:mi),dlnnda(:mi),dlnTda(:mi),duda(:mi),gh(:mi))
      !$acc update self(back_f0_i(:mi))

      if(acc_newer_D_zion) then
         !$acc update self(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
         acc_newer_D_zion=.false.
         acc_newer_H_zion=.false.
      endif

      ! fluxes diagnose
      !    wpi(1,1:mi)=wpi(1,1:mi)*gha(1:mi)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
      !$omp parallel do private(m)
      do m=1,mi
         wpi(1,m)=wpi(1,m)*gha(m)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
      end do

      rmarker=0.0_wp
      eflux=0.0_wp
      pfluxi=0.0_wp
      dflowi=0.0_wp
      entropyi=0.0_wp
      fluxpsi=0.0_wp
      heatflux=0.0_wp
      En_00=0.0_wp
      En_tur=0.0_wp
      En_tot=0.0_wp
      !$acc parallel loop present(zion,zion0)
      do m=1,mi
         r=zion0(7,m)         ! old radial location
         ip=max(0,min(mpsi,int((r-s0)*delr+0.5_wp)))    ! ip [0,mpsi]
         vdrenergy=wpi(1,m)*(wpi(2,m)-2.5*vth_grid(ip)**2)*zion0(5,m)
         ! note: wpi(1,m)=dr/dt (EXB part)

         ! radial location in diagnostic bin
         i=max(1,min(mflux,1+int((r-s0)*real(mflux)/(s1-s0))))

         ! marker,energy,particle,momentum fluxes,parallel flows,entropy and kinetic energy
         rmarker(i)=rmarker(i)+zion0(8,m)
         eflux(i)=eflux(i)+vdrenergy
         pfluxi=pfluxi+wpi(1,m)*zion0(5,m)
         dflowi=dflowi+BB(m)*zion0(4,m)*zion0(5,m)
         entropyi=entropyi+zion0(5,m)*zion0(5,m)
      enddo


      if(ZF_SEPARATION==0) then
         !     note: zonal flow component <phi> and phip00=d<phi>/da are calculated in field.f90
         adum=0.0_wp
         adum_00=0.0_wp
         adum_tur=0.0_wp
         adum_tot=0.0_wp

         !$omp parallel private(i,j,ij,E_a,E_t,E_f) reduction(+:adum,adum_00,adum_tur,adum_tot)
         ! Note from OMP compiler: A variable in the REDUCTION clause or any component or element thereof must not have the POINTER
         ! attribute; nor be a pointee, an allocatable object, an assumed-size array, or a THREADLOCAL common block variable.
         !$omp do
         do i=0,mpsi
            do j=1,mtheta(i)
               ij=igrid(i)+j
               adum(i)=adum(i)+dtemper(ij)*(phi(1,ij)-phi00(i))**2  ! surface averaging for delta-phi

               adum_00(i)=adum_00(i)+dtemper(ij)*phip00(i)**2*gij(1,ij)  ! surface average

               E_a=evector(1,1,ij)+Er_grid(i)           ! total radial component   ! d\phi/da
               E_t=evector(2,1,ij)                      ! d\delta\phi/d\theta
               E_f=evector(3,1,ij)                      ! d\delta\phi/d\zeta
               adum_tot(i)=adum_tot(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                  +E_t**2*gij(2,ij)  &
                  +E_f**2*gij(3,ij)  &
                  +2.0_wp*E_a*E_t*gij(4,ij))

               E_a=E_a-phip00(i)-Er_grid(i) ! turbulence component only    ! d\delta\phi/da
               adum_tur(i)=adum_tur(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                  +E_t**2*gij(2,ij)  &
                  +E_f**2*gij(3,ij)  &
                  +2.0_wp*E_a*E_t*gij(4,ij))
            enddo
         enddo
         !$omp end parallel

      else  ! ZF_SEPARATION=1

         adum=0.0_wp
         adum_00=0.0_wp
         adum_tur=0.0_wp
         adum_tot=0.0_wp
         adum_11=0.0_wp

         ! DEBUG
         !        if(mype==100) print*, 'Dirichlet2', phi(1,igrid(0)), phi(1,igrid(0)+mtheta(0)/2),phi(1,igrid(0)+mtheta(0)-1), &
         !                               phi(1,igrid(mpsi)), phi(1,igrid(mpsi)+mtheta(mpsi)/2),phi(1,igrid(mpsi)+mtheta(mpsi)-1)

         !$omp parallel private(i,j,ij,E_a,E_t,E_f) reduction(+:adum,adum_00,adum_tur,adum_tot,adum_11)
         !$omp do
         do i=0,mpsi
            do j=1,mtheta(i)
               ij=igrid(i)+j
               !          fluxpsi(i)=fluxpsi(i)+phi(1,ij)*phi(1,ij)               ! what average?
               adum(i)=adum(i)+dtemper(ij)*phi(1,ij)*phi(1,ij)  ! surface averaging for <delta-phi^2>

               adum_11(i)=adum_11(i)+dtemper(ij)*phi(1,ij)      ! to check surface averaged <delta-phi> which should be zero by defination

               adum_00(i)=adum_00(i)+dtemper(ij)*phip00(i)**2*gij(1,ij)  ! surface average

               E_a=evector(1,1,ij)+Er_grid(i)           ! total radial component   ! d\phi/da
               E_t=evector(2,1,ij)                      ! d\delta\phi/d\theta
               E_f=evector(3,1,ij)                      ! d\delta\phi/d\zeta
               adum_tot(i)=adum_tot(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                  +E_t**2*gij(2,ij)  &
                  +E_f**2*gij(3,ij)  &
                  +2.0*E_a*E_t*gij(4,ij))

               E_a=E_a-phip00(i)-Er_grid(i) ! turbulence component only    ! d\delta\phi/da
               adum_tur(i)=adum_tur(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                  +E_t**2*gij(2,ij)  &
                  +E_f**2*gij(3,ij)  &
                  +2.0_wp*E_a*E_t*gij(4,ij))
            enddo
         enddo
         !$omp end parallel
      end if   !  end if(ZF_SEPARATION==0)

      call MPI_ALLREDUCE(adum_11,flux_11,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      flux_11=flux_11/jaco_sum/Tc   ! =<phi>/Tc -- added in March 2013 to check whether <phi>=0 is obtained

      call MPI_ALLREDUCE(adum,fluxpsi,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      fluxpsi=fluxpsi/jaco_sum/(Tc*Tc)  ! =<(phi/Tc)^2> -- turbulence field intensity
      sumfluxpsi=sumfluxpsi+fluxpsi     ! =int[I(a)dt]

      call MPI_ALLREDUCE(adum_00,En_00,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      En_00=En_00/jaco_sum      ! =<E_00^2>    -- zonal flow energy

      call MPI_ALLREDUCE(adum_tur,En_tur,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      En_tur=En_tur/jaco_sum     ! =<\deltaE^2> -- turbulence energy

      call MPI_ALLREDUCE(adum_tot,En_tot,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      En_tot=En_tot/jaco_sum     ! =<E^2> -- total field energy

      call get_fluxes(wpi,dgf,dgq)

      if(mod(istep,ndiag)==0) then
         call get_3D_flow(da,dgq,dgf,wpi(2,:))
      end if

      call cpu_time(cput1)
      tt1=MPI_WTIME()
      time_in_Pushi(1) = time_in_Pushi(1) + (cput1-cput0)
      timewc_in_Pushi(1) = timewc_in_Pushi(1) + (tt1-tt0)
   endif

   if(irk_i==1) then
      call cpu_time(cput0)
      tt0=MPI_WTIME()

      call Measure_and_Cumulate_iZonData_on_DG(wpi,da)

      call cpu_time(cput1)
      tt1=MPI_WTIME()
      time_in_Pushi(1) = time_in_Pushi(1) + (cput1-cput0)
      timewc_in_Pushi(1) = timewc_in_Pushi(1) + (tt1-tt0)

#ifdef _MGY_TIMER_
      if(mype==0) then
         write(stdout,*) "Measure_Cumulate_iFlux_a1D takes", tt1-tt0," secs"
         call flush(stdout)
      endif
#endif
  endif

   !Another Diagnostics on DG (Diagnostic Grid)
   if(irk_i==1 .and. mod(istep,msnap)==0) then
      call cpu_time(cput0)
      tt0=MPI_WTIME()

      !!!!$acc update self(Fa(:mi),PP(:mi),Bq(:mi),Ba(:mi)) ! these are not used for now, but possibly for future diagnostics
      !$acc update self(BB(:mi),TT(:mi),Ta(:mi),ghq(:mi),gha(:mi),Pa(:mi),F(:mi))
      !$acc update self(dgr(:mi),da(:mi),dgq(:mi),dgf(:mi),wpi(:,:mi))
      !$acc update self(ii_inv_grad_a(:mi),ii_inv_grad_gq(:mi),ii_inv_grad_gf(:mi))
      !$acc update self (dvi_gf(:mi),vpara_dvi_a(:mi),vpara_dvi_gq(:mi),vpara_dvi_gf(:mi))
      !$acc update self(UnNorm_theta_i(:mi),UnNorm_zeta_i(:mi))
      !$acc update self(wi_init(:mi),wi_vpara(:mi),wi_vperp(:mi),wi_accel(:mi),wi_fdamp(:mi))
      !$acc update self(ilabels(:mi))
      !$acc update self(backn(:mi),backT(:mi),backu(:mi),dlnnda(:mi),dlnTda(:mi),duda(:mi),gh(:mi))
      !$acc update self(back_f0_i(:mi))

      if(acc_newer_D_zion) then
         !$acc update self(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
         acc_newer_D_zion=.false.
         acc_newer_H_zion=.false.
      endif

      do m=1,mi
         FullVelocity_da_dgq_dgf(1,m) = da(m)
         FullVelocity_da_dgq_dgf(2,m) = dgq(m)
         FullVelocity_da_dgq_dgf(3,m) = dgf(m)
      enddo
      call Measure_Ions_on_Diagnostic_Grid(wpi(1:3,:),FullVelocity_da_dgq_dgf)

      call cpu_time(cput1)
      tt1=MPI_WTIME()
      time_in_Pushi(1) = time_in_Pushi(1) + (cput1-cput0)
      timewc_in_Pushi(1) = timewc_in_Pushi(1) + (tt1-tt0)
#ifdef _MGY_TIMER_
      if(mype==0) then
         write(stdout,*) "Measure_Ions_on_Diagnostic_Grid takes", tt1-tt0," secs"
         call flush(stdout)
      endif
#endif
   endif



   ! Update ion particle weight
   call cpu_time(cput0)
   tt0=MPI_WTIME()
   call shifti
   call cpu_time(cput1)
   tt1=MPI_WTIME()
   time_in_Pushi(2) = time_in_Pushi(2) + (cput1-cput0)
   timewc_in_Pushi(2) = timewc_in_Pushi(2) + (tt1-tt0)

   call n_T_u_group(flag_backf0=.true.)

   !$acc parallel loop present(zion,zion0,back_f0_i)
   do m=1,mi
      new_f0 = zion0(6,m)*back_f0_i(m)
      ! zion(5,m) = zion0(5,m) - (new_f0-zion(8,m))
      zion(5,m) = zion0(8,m) - new_f0
      ! zion0(8,m) is conserved
   enddo


end subroutine pushion_wEvol3



!========================================
! MGYOO (03/29/2018) :
!========================================
subroutine ion_locator_at_larmor_points

   use precision
   use global_parameters
   use particle_array
   use field_array
   implicit none

   integer m,i,im,j,k,ip,jt,ii,j11,j10,j01,j00,larmor,ij,ipjt

   real(wp) rdum,tdum,delr,delt(0:mpsi),delz,r,wp1,wp0,wt11,wt10,wt01,wt00,&
      tflr,pi2_inv,psitmp,thetatmp,zetatmp,rhoi,psimax,psimin,pdum,thedum


   integer inner_ip,outter_ip, inner_jt, outter_jt
   real(wp) gyro_phase, intpl_pgyro(larmor_num), intpl_tgyro(larmor_num)

   delr=1.0_wp/deltar
   delt=2.0_wp*pi/deltat
   delz=1.0_wp/deltaz
   pi2_inv=0.5_wp/pi
   !  psimin=psi_grid(0)
   !  psimax=psi_grid(mpsi)
   psimin=s_grid(0)
   psimax=s_grid(mpsi)


   !$omp parallel do private(m,larmor,psitmp,thetatmp,zetatmp,rhoi,r,ip,jt,ipjt,&
   !$omp& pdum,rdum,thedum,ii,wp1,tflr,im,tdum,j00,j01)
   if(acc_newer_H_zion) then
      !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
      acc_newer_H_zion=.false.
      acc_newer_D_zion=.false.
   endif
   !$acc parallel loop present(zion,jtion0,jtion1,wzion,wpion,wtion0,wtion1,zion_s_on_GyroOrbit,zion_gq_on_GyroOrbit) private(intpl_pgyro,intpl_tgyro,thedum)
   do m=1,mi
      !psitmp=zion(1,m)   ! =a
      psitmp=zion(7,m)   ! =s
      thetatmp=zion(2,m)
      zetatmp=zion(3,m)
      rhoi=zion(6,m)

      r=psitmp
      ip=max(0,min(mpsi,int((r-s0)*delr+0.5)))  ! grid size is uniform in s but not a
      !   jt=max(0,min(mtheta(ip),int(thetatmp*pi2_inv*delt(ip)+0.5)))
      thedum=thetatmp-zetatmp*qtinv(ip)    ! corresponding theta on zeta=0 plane
      thedum=thedum*pi2_inv+100.0_wp
      thedum=2.0_wp*pi*(thedum-int(thedum))       ! bring back to range [0,2pi]
      jt=max(0,min(mtheta(ip),int(thedum*pi2_inv*delt(ip)+0.5)))
      ipjt=igrid(ip)+jt     ! note special case of jt=0
      wzion(m)=(zetatmp-zetamin)*delz


      !==========================================
      ! MGYOO (10/24/2018) : to store information of GC position.
      ! GC's lower-left cell-index(jtion) and linear-interpolation-weights (wpion, wtion0, wtion1)
      inner_ip = max(0,min(mpsi-1,int((r-s0)*delr)))
      thedum=zion(2,m)-zetatmp*qtinv(inner_ip)              ! corresponding theta on zeta=0 plane
      thedum=modulo(thedum,2.0_wp*pi)/deltat(inner_ip)  ! bring back to range [0,2pi] and divde by delta_theta
      inner_jt = int(thedum)
      wtion0(0,m) = thedum - int(thedum)
      jtion0(0,m) = igrid(inner_ip) + inner_jt
      wpion(0,m) =  (zion(7,m) - s_grid(inner_ip))/(s_grid(inner_ip+1)-s_grid(inner_ip))

      outter_ip = inner_ip+1
      thedum=zion(2,m)-zetatmp*qtinv(outter_ip)     ! corresponding theta on zeta=0 plane
      thedum=modulo(thedum,2.0_wp*pi)/deltat(outter_ip)  ! bring back to range [0,2pi] and divde by delta_theta
      outter_jt = int(thedum)
      wtion1(0,m) = thedum - int(thedum)
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
         zion_s_on_GyroOrbit(larmor,m)=psitmp+rhoi*intpl_pgyro(larmor)
         rdum=delr*max(0.0_wp,min(s1-s0,pdum-s0))
         ii=max(0,min(mpsi-1,int(rdum)))
         !        wp1=rdum-real(ii)
         !        wpion(larmor,m)=wp1
         wp1=(pdum-s_grid(ii))/(s_grid(ii+1)-s_grid(ii))   !wx  finite rho_i case
         wpion(larmor,m)=wp1                               !wx  finite rho_i case

         ! particle position in theta
         tflr=thetatmp+rhoi*intpl_tgyro(larmor)
         zion_gq_on_GyroOrbit(larmor,m)=modulo(tflr,2.0*pi)

         ! inner flux surface
         im=ii
         tdum=pi2_inv*(tflr-zetatmp*qtinv(im))+10.0_wp
         tdum=(tdum-int(tdum))*delt(im)
         j00=max(0,min(mtheta(im)-1,int(tdum)))
         jtion0(larmor,m)=igrid(im)+j00  ! pay attention to j00=0 case
         wtion0(larmor,m)=tdum-real(j00)

         ! outer flux surface
         im=ii+1
         tdum=pi2_inv*(tflr-zetatmp*qtinv(im))+10.0_wp
         tdum=(tdum-int(tdum))*delt(im)
         j01=max(0,min(mtheta(im)-1,int(tdum)))
         jtion1(larmor,m)=igrid(im)+j01  ! pay attention to j01=0 case
         wtion1(larmor,m)=tdum-real(j01)
      enddo
   enddo

end subroutine ion_locator_at_larmor_points
