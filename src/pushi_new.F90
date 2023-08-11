subroutine pushion_wEvol0or1(irk_i,nrk)
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
   integer i,im,ip,j,kz,m,ii,ierror,larmor,ipjt
   real(wp) rdum,r,grada,vpara,energy,kappa,vdr,wdot,wdrift,&
      vdrenergy,delr,delt(0:mpsi),delz,pi2_inv,adum(0:mpsi),&
      xdum,ydum,xdot,ydot,&
      epara,psimax,psimin,pdum,wdrive,ainv,psitmp,thetatmp,&
      zetatmp,rhoi,wpi(3,mi),gaE(mi),wgaE(3,mi),gatE(mi),rdel, &
      dgr(mi),da(mi),dgq(mi),dgf(mi),mu, &
      Er_particle(mi),dela,fip,weqflow,T_inverse,&
      B_inverse,energy1,vdq,wKH,wgradB,omega_t,u_vpara,  &
      v2 
      ! following working arrays are defined to use OMP reduction:
   real(wp) vdgf, vd_para_mag, vd_para_a_cntv, vd_para_gq_cntv, vd_perp_a_cntv, vd_perp_gq_cntv
   real(wp) wdot_vpara, wdot_vperp, wdot_accel, now_w
   real(wp) tot_B_vec_cov(3), tot_B_vec_cntv(3) ! covariant and contravariant components of total B field (B0 + dB) incluidng the magnetic perturbation
   real(wp) grad_lnf0_a_cov, grad_lnf0_gq_cov ! covariant components of grad(ln(f0)) in a and gq axis.
   !(gf component is 0 due to the toroidal axisymmetric assumption )
   real(wp) gradPhi(3,mi),gaH(mi),gaS(mi),grad_gaH(3,mi),grad_gaS(3,mi)

   integer gcmotion, gcmotion_em, gcmotion_em_MGY, gcmotion_em_MGY_vec_form, gcmotion_em_MGY_vec_form_EM2

   real(wp) tt0,tt1,cput0,cput1,rtc

   real(wp), external :: D01AHF,dsda !,SMOOTHSTEP
   integer minpts,ifail
   real(wp) acc,error_8
   real(wp) E_a,E_t,E_f,timef

   real(wp) ddT_da_particle, muB

   real(wp) RK_cum_coeff, RK_adv_coeff ! Coefficients to cumlate (to get future final value) or advance(to predict next RK step) RK4 terms
   real(wp) FullVelocity_da_dgq_dgf(3,mi) ! ds= a_dot*(Tc/T)**(dsda_power), dgq=gq_dot, dgf=dgf_dot

   !$acc declare create(dgr,da,dgq,dgf,wpi,gradPhi,Er_particle,gaE,wgaE,gatE,gaH,gaS,grad_gaH,grad_gaS)

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



   call Determine_RK_coefficients(RK_adv_coeff, RK_cum_coeff, irk_i, nrk)

   if(irk_i==1)then ! 1st step of Runge-Kutta method
      !$acc parallel loop gang vector present(zion,zion0,RK_cum_zion)
      do m=1,mi
         ! Note that zion0(8,m) does not mean "mu", this is reserved for the other quantity
         zion0(1:5,m)=zion(1:5,m)
         zion0(7,m)=zion(7,m)

         RK_cum_zion(:,m)=0.0
      enddo
      acc_newer_D_zion=.true. !device side zion will be changed in the following acc routine
   endif
   ! Followings are not required here, so commented out. (MGYOO: 11/13/2019)
   call n_T_u_group  ! prepare n, T, dlnn/da, dlnT/da, etc.

   call Cal_ionLarmor_s_gq_gf_using_SG_gyro_over_mu
   if (electromagnetic==2) then
      call Get_gradPhi_gaHS_gradgaHS_at_zion(gradPhi,gaH,gaS,grad_gaH,grad_gaS) 
   else
      call Get_gradPhi_at_zion(gradPhi)
   endif


#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   !MGYOO (10/30/2018; to use analytic profile of ga0(a))
   !$acc parallel loop present(zion,gaE,wgaE,gatE)
   do m=1,mi
      call Get_extMP_alpha_amplitude(zion(1,m),zion(2,m),zion(3,m),istep*tstep, gaE(m), wgaE(:,m), gatE(m))
   enddo
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion, Get_extMP_alpha_amplitude", tt1-tt0," secs"
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
   if(electromagnetic == 1) then
      if(Flag_GC_Solver==0) then
         m=gcmotion_em_MGY(dgr,da,dgq,dgf,zion,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,mi,ziondim);
      elseif(Flag_GC_Solver==1) then
         m=gcmotion_em_MGY_vec_form(dgr,da,dgq,dgf,zion,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,1,mi,ziondim);
      else
         write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
         call flush(stdout)
         stop
      endif
   elseif(electromagnetic==2) then
      m=gcmotion_em_MGY_vec_form_EM2(dgr,da,dgq,dgf,zion,gradPhi,gaE,wgaE,gatE,gaH,gaS,grad_gaH,grad_gaS,Er_particle,nonlinear,1,mi,ziondim)
   else
      m=gcmotion(dgr,da,dgq,dgf,zion,gradPhi,Er_particle,nonlinear,mi,ziondim)
   end if

#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion_wEvol0or1, ion gcmotion takes", tt1-tt0," secs"
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


#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif

   ! update GC position
   acc_newer_D_zion=.true. !device side zion will be changed in the following acc routine
   !$acc parallel loop gang vector present(zion,zion0,RK_cum_zion)
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


      if(ifneoclassical==1) then
         wdrift=(dgr(m)*BB(m)+zion(4,m)*weqflow)*(vpara-backu(m))*T_inverse   !  full drho/dt -- including NC dynamics
      else
         wdrift=(PP(m)*BB(m)+zion(4,m)*weqflow)*(vpara-backu(m))*T_inverse   ! see note
      endif

      wgradB=weqflow*T_inverse*mu       ! interesting term propotional to grad_B

      wKH=-(vdr*(duda(m)*F(m)+Fa(m)*omega_t-Ba(m)*backu(m)) &
         -vdq*Bq(m)*backu(m))*B_inverse*vpara*T_inverse       ! major Kelvin-Helmholtz drive



      now_w = zion(5,m)

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

         wi_vpara(m) = wi_vpara(m) + tstep*( (zion0(8,m)-nonlinear*now_w)*wdot_vpara)
         wi_vperp(m) = wi_vperp(m) + tstep*( (zion0(8,m)-nonlinear*now_w)*wdot_vperp)
         wi_accel(m) = wi_accel(m) + tstep*( (zion0(8,m)-nonlinear*now_w)*wdot_accel)
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
         RK_cum_zion(5,m) = RK_cum_zion(5,m) +RK_cum_coeff*tstep*((zion0(8,m)-nonlinear*now_w)*wdot)
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
            zion(5,m) = zion0(5,m)+RK_adv_coeff*tstep*((zion0(8,m)-nonlinear*now_w)*wdot)
         endif

      else if(irk_i==4) then !irk_i==4
         ! Get the future value at time level (n+1): f_(n+1) = f_(n) + (1/6)*(k1+2*k2+2*k3+k4) = f_(n)+RK_cum_value
         zion(1:4,m) = zion0(1:4,m) +RK_cum_zion(1:4,m)
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
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion_wEvol0or1, main do loop takes", tt1-tt0," secs"
   call flush(stdout)
#endif


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
      !$acc update self(wi_init(:mi),wi_vpara(:mi),wi_vperp(:mi),wi_accel(:mi))
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



end subroutine pushion_wEvol0or1



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
   integer i,ip,m,ierr,ierror,j
   real(wp) RK_cum_coeff, RK_adv_coeff ! Coefficients to cumulate (to get future final value) or advance(to predict next RK step)
   real(wp) dgr(mi),da(mi),dgq(mi),dgf(mi),gradPhi(3,mi),wpi(3,mi),Er_particle(mi),gaH(mi),gaS(mi),grad_gaH(3,mi),grad_gaS(3,mi)
   real(wp) gaE(mi),wgaE(3,mi), gatE(mi) ! Magnetic Perturbation Information: alpha, grad_alpha, partial derivative d(gaE)/dt
   real(wp) tt0,tt1,cput0,cput1 ! For timer
   real(wp) timef,dela,r,fip
   integer gcmotion, gcmotion_em, gcmotion_em_MGY, gcmotion_em_MGY_vec_form, gcmotion_em_MGY_vec_form_EM2
   real(wp) mu,vpara,energy
   real(wp) E_a,E_t,E_f,delr
   real(wp) FullVelocity_da_dgq_dgf(3,mi) ! ds= a_dot*(Tc/T)**(dsda_power), dgq=gq_dot, dgf=dgf_dot
   real(wp) vdrenergy,new_f0
   integer larmor, loop

   !$acc declare create(dgr,da,dgq,dgf,gradPhi,wpi,Er_particle,gaE,wgaE,gatE,gaH,gaS,grad_gaH,grad_gaS)

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

   call Determine_RK_coefficients(RK_adv_coeff, RK_cum_coeff, irk_i, nrk)

   if(irk_i==1)then ! 1st step of Runge-Kutta method
      !$acc parallel loop gang vector present(zion,zion0,RK_cum_zion)
      do m=1,mi
         zion0(1:5,m)=zion(1:5,m)
         zion0(7,m)=zion(7,m)

         RK_cum_zion(:,m)=0.0
      enddo
      acc_newer_D_zion=.true. !device side zion will be changed in the following acc routine
   endif

   ! Followings are not required here, so commented out. (MGYOO: 11/13/2019)

   if(irk_i==1) then
      !$acc parallel loop present(zion,back_f0_i)
      do m=1,mi
         zion(8,m) = zion0(6,m)*back_f0_i(m) ! zion(8,m) is temporary old_f0 @ irk=1 step
      enddo
   endif

   ! gather e_field using multi-point gyro-averaging, sorting in poloidal angle
   ! The following line is an OpenMP directive for loop-level parallelism
   ! on shared memory machines (see http://www.openmp.org).
   call n_T_u_group

   call Cal_ionLarmor_s_gq_gf_using_SG_gyro_over_mu
   if (electromagnetic==2) then
      call Get_gradPhi_gaHS_gradgaHS_at_zion(gradPhi,gaH,gaS,grad_gaH,grad_gaS) 
   else
      call Get_gradPhi_at_zion(gradPhi)
   endif




#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   !MGYOO (10/30/2018; to use analytic profile of ga0(a))
   !$acc parallel loop present(zion,gaE,wgaE,gatE)
   do m=1,mi
      call Get_extMP_alpha_amplitude(zion(1,m),zion(2,m),zion(3,m),istep*tstep, gaE(m), wgaE(:,m), gatE(m))
   enddo
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion, Get_extMP_alpha_amplitude", tt1-tt0," secs"
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
   if(electromagnetic == 1) then
      if(Flag_GC_Solver==0) then
         m=gcmotion_em_MGY(dgr,da,dgq,dgf,zion,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,mi,ziondim);
      elseif(Flag_GC_Solver==1) then
         m=gcmotion_em_MGY_vec_form(dgr,da,dgq,dgf,zion,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,1,mi,ziondim);
      else
         write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
         call flush(stdout)
         stop
      endif
   elseif(electromagnetic==2) then
      m=gcmotion_em_MGY_vec_form_EM2(dgr,da,dgq,dgf,zion,gradPhi,gaE,wgaE,gatE,gaH,gaS,grad_gaH,grad_gaS,Er_particle,nonlinear,1,mi,ziondim)
   else
      m=gcmotion(dgr,da,dgq,dgf,zion,gradPhi,Er_particle,nonlinear,mi,ziondim)
   end if
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion_wEvol3, ion gcmotion takes", tt1-tt0," secs"
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

   !$acc parallel loop gang vector present(zion,zion0,RK_cum_zion)
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
         zion(1:4,m) = zion0(1:4,m) +RK_cum_zion(1:4,m)
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
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion_wEvol3, main do loop takes", tt1-tt0," secs"
   call flush(stdout)
#endif

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
      !$acc update self(wi_init(:mi),wi_vpara(:mi),wi_vperp(:mi),wi_accel(:mi))
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




subroutine pushion_wEvol4(irk_i,nrk)
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
   integer i,ip,m,ierr,ierror,j
   real(wp) RK_cum_coeff, RK_adv_coeff ! Coefficients to cumulate (to get future final value) or advance(to predict next RK step)
   real(wp) dgr(mi),da(mi),dgq(mi),dgf(mi),gradPhi(3,mi),gaH(mi),gaS(mi),grad_gaH(3,mi),grad_gaS(3,mi),wpi(3,mi),Er_particle(mi)
   real(wp) gaE(mi),wgaE(3,mi), gatE(mi) ! Magnetic Perturbation Information: alpha, grad_alpha, partial derivative d(ga)/dt
   real(wp) tt0,tt1,cput0,cput1 ! For timer
   real(wp) timef,dela,r,fip
   integer gcmotion, gcmotion_em, gcmotion_em_MGY, gcmotion_em_MGY_vec_form,gcmotion_em_MGY_vec_form_EM2
   real(wp) mu,vpara,energy
   real(wp) inv_ds
   real(wp) E_a,E_t,E_f,delr
   real(wp) FullVelocity_da_dgq_dgf(3,mi) ! ds= a_dot*(Tc/T)**(dsda_power), dgq=gq_dot, dgf=dgf_dot
   real(wp) vdrenergy,new_f0
   integer larmor, loop

   !$acc declare create(dgr,da,dgq,dgf,gradPhi,gaH,gaS,grad_gaH,grad_gaS,wpi,Er_particle,gaE,wgaE,gatE)

   inv_ds = real(mpsi)/(s1-s0)

   if(acc_newer_H_zion) then
      !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
      acc_newer_H_zion=.false.
      acc_newer_D_zion=.false.
   endif

   if(istep==istart .and. irk_i==1) then
      time_in_Pushi =0
      timewc_in_Pushi =0
   endif


   call Determine_RK_coefficients(RK_adv_coeff, RK_cum_coeff, irk_i, nrk)


   if(irk_i==1)then ! 1st step of Runge-Kutta method
      !$acc parallel loop gang vector present(zion,zion0,RK_cum_zion)
      do m=1,mi
         zion0(1:5,m)=zion(1:5,m)
         zion0(7,m)=zion(7,m)

         RK_cum_zion(:,m)=0.0
      enddo
      acc_newer_D_zion=.true. !device side zion will be changed in the following acc routine
   endif

   if(irk_i==1) then
      !$acc parallel loop present(zion,back_f0_i)
      do m=1,mi
         zion(8,m) = zion0(6,m)*back_f0_i(m) ! zion(8,m) is temporary old_f0 @ irk=1 step
      enddo
   endif

   call n_T_u_group

   call Cal_ionLarmor_s_gq_gf_using_SG_gyro_over_mu
   if (electromagnetic==2) then
      call Get_gradPhi_gaHS_gradgaHS_at_zion(gradPhi,gaH,gaS,grad_gaH,grad_gaS) 
   else
      call Get_gradPhi_at_zion(gradPhi)
   endif



#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   !MGYOO (10/30/2018; to use analytic profile of ga0(a))
   !$acc parallel loop present(zion,gaE,wgaE,gatE)
   do m=1,mi
      call Get_extMP_alpha_amplitude(zion(1,m),zion(2,m),zion(3,m),istep*tstep, gaE(m), wgaE(:,m), gatE(m))
   enddo
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion, Get_extMP_alpha_amplitude", tt1-tt0," secs"
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
   if(electromagnetic == 1) then
      if(Flag_GC_Solver==0) then
         m=gcmotion_em_MGY(dgr,da,dgq,dgf,zion,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,mi,ziondim);
      elseif(Flag_GC_Solver==1) then
         m=gcmotion_em_MGY_vec_form(dgr,da,dgq,dgf,zion,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,1,mi,ziondim);
      else
         write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
         call flush(stdout)
         stop
      endif
   elseif(electromagnetic==2) then
      m=gcmotion_em_MGY_vec_form_EM2(dgr,da,dgq,dgf,zion,gradPhi,gaE,wgaE,gatE,gaH,gaS,grad_gaH,grad_gaS,Er_particle,nonlinear,1,mi,ziondim)
   else
      m=gcmotion(dgr,da,dgq,dgf,zion,gradPhi,Er_particle,nonlinear,mi,ziondim)
   end if
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion_wEvol4, ion gcmotion takes", tt1-tt0," secs"
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

   !$acc parallel loop gang vector present(zion,zion0,RK_cum_zion)
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
         zion(1:4,m) = zion0(1:4,m) +RK_cum_zion(1:4,m)
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
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pushion_wEvol4, main do loop takes", tt1-tt0," secs"
   call flush(stdout)
#endif

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
      !$acc update self(wi_init(:mi),wi_vpara(:mi),wi_vperp(:mi),wi_accel(:mi))
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


end subroutine pushion_wEvol4



subroutine pullback_ion
   use precision
   use global_parameters
   use particle_array
   use Structured_Grid
   use function_gg_module

   implicit none

   integer m
   real(wp) old_gaH_at_zion(mi), old_back_f0_i(mi), old_grad_Phi_at_zion(3,mi)
   !$acc declare create(old_gaH_at_zion,old_back_f0_i)

   call n_T_u_group ! To get back_f0_i

   call Get_SG_scalar_at_zion(SG%gaH,old_gaH_at_zion)
   call Get_SG_vector_at_zion(SG%phi,old_grad_Phi_at_zion)

   !$acc parallel loop present(zion,back_f0_i,old_back_f0_i,old_gaH_at_zion)
   do m=1,mi
      old_back_f0_i(m) = back_f0_i(m)
      zion(4,m) = zion(4,m) - old_gaH_at_zion(m)
   enddo

   call n_T_u_group ! To get new back_f0_i corresponding to the new zion(4,m) (rho_H)

   ! Update particle weight
   !$acc parallel loop present(zion,back_f0_i,old_back_f0_i,old_gaH_at_zion)
   do m=1,mi
      zion(5,m) = zion(5,m) + zion0(6,m)*(old_back_f0_i(m)-back_f0_i(m))
   enddo
end subroutine pullback_ion


!========================================
! MGYOO (03/29/2018) :
!========================================
subroutine ion_locator_at_larmor_points

   use precision
   use global_parameters
   use particle_array
   use field_array
   implicit none

   integer m,i,im,j,k,ip,jt,ii,j11,j10,j01,j00,larmor,ipjt

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
   !$acc parallel loop gang vector present(zion,jtion0,jtion1,wzion,wpion,wtion0,wtion1,ionLarmor_s,ionLarmor_gq) private(intpl_pgyro,intpl_tgyro,thedum)
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
         ionLarmor_s(larmor,m)=psitmp+rhoi*intpl_pgyro(larmor)
         rdum=delr*max(0.0_wp,min(s1-s0,pdum-s0))
         ii=max(0,min(mpsi-1,int(rdum)))
         !        wp1=rdum-real(ii)
         !        wpion(larmor,m)=wp1
         wp1=(pdum-s_grid(ii))/(s_grid(ii+1)-s_grid(ii))   !wx  finite rho_i case
         wpion(larmor,m)=wp1                               !wx  finite rho_i case

         ! particle position in theta
         tflr=thetatmp+rhoi*intpl_tgyro(larmor)
         ionLarmor_gq(larmor,m)=modulo(tflr,2.0*pi)

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
