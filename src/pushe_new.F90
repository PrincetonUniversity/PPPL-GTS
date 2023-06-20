subroutine pusheon_wEvol0or1(solver1, masssolver, irk_e, nrk)
   use precision
   use global_parameters
   use particle_array
   use field_array
   use diagnosis_array
   use esi_interface
   use solver_data_class
   use function_gg_module
   use Ptl_Src_Loss
   use Structured_Grid
   use Diagnostic_Grid
   use MGY_PETSc_Interface_Data 

   implicit none

   type(solver_data) :: solver1, masssolver

   integer, intent(in) :: irk_e,nrk
   integer i,im,ip,j,jt,k,kz,kk,m,ii,ierror,ipjt
   real(wp) rdum,r,grada,energy,kappa,vdr,wdot,wdrift,&
      vdrenergy,delr,delt(0:mpsi),delz,pi2_inv,adum(0:mpsi),&
      xdum,ydum,xdot,ydot,&
      epara,psimax,psimin,pdum,wdrive,ainv,psitmp,&
      zetatmp,gradPhi(3,me),wpi(5,me),gaE(me),wgaE(3,me),gatE(me),&
      dgr(me),da(me),dgq(me),dgf(me),mu, &
      Er_particle(me),dela,fip,weqflow,wgradB,T_inverse,Ea_00,wdphidt,wdphidt2, &
      vpara, v2,  f_rho
      ! following working arrays are defined to use OMP reduction:
   real(wp) Phi_at_zeon(memax),gaHS(me),grad_gaH(3,me),grad_gaS(3,me)
   real(wp) now_wAD, wdot_fullh ! wdot_fullh
   real(wp) tdum,thetatmp
   real(Wp) tmp_old_AD, tmp_new_AD, tmp_old_dADdt, tmp_new_dADdt
   real(wp) vdq, vdgf, vd_para_mag, vd_para_a_cntv, vd_para_gq_cntv, vd_perp_a_cntv, vd_perp_gq_cntv
   real(wp) wdot_vpara, wdot_vperp, wdot_accel, now_w
   real(wp) tot_B_vec_cov(3), tot_B_vec_cntv(3) ! covariant and contravariant components of total B field (B0 + dB) incluidng the magnetic perturbation
   real(wp) grad_lnf0_a_cov, grad_lnf0_gq_cov ! covariant components of grad(ln(f0)) in a and gq axis.
   !(gf component is 0 due to the toroidal axisymmetric assumption )

   integer gcmotion_eon, gcmotion_eon_em, gcmotion_eon_em_MGY, gcmotion_eon_em_MGY_vec_form, gcmotion_eon_em_MGY_vec_form_EM2

   real(wp) tt0,tt1,cput0,cput1,rtc

   real(wp), external :: D01AHF,dsda!,SMOOTHSTEP
   integer minpts,ifail
   real(wp)  acc,error_8
   real(wp)  E_f,timef

   real(wp) muB

   real(wp) RK_cum_coeff, RK_adv_coeff ! Coefficients to cumlate (to get future final value) or advance(to predict next RK step) RK4 terms
   real(wp) FullVelocity_da_dgq_dgf(3,me) ! ds= a_dot*(Tc/T)**(dsda_power), dgq=gq_dot, dgf=dgf_dot
   !$acc declare create(dgr,da,dgq,dgf,gradPhi,wpi,Er_particle,gaE,wgaE,gatE,Phi_at_zeon,gaHS,grad_gaH,grad_gaS)
   real(wp) tmp_zeon5(memax), tmp_zeon8(memax)
   !$acc declare create (tmp_zeon5, tmp_zeon8)
   integer :: loop
   real(wp) prev_phi00(0:mpsi), phi00_convergency

   if(acc_newer_H_zeon) then
      !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
      acc_newer_H_zeon=.false.
      acc_newer_D_zeon=.false.
   endif

   if(istep==istart .and. irk_e==1) then
      time_in_Pushe =0
      timewc_in_Pushe =0
   endif

   delr=1.0_wp/deltar
   delt=1.0_wp/deltat
   delz=1.0_wp/deltaz
   pi2_inv=0.5_wp/pi
   psimax=psi_grid(mpsi)
   psimin=psi_grid(0)

   f_rho=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)

   
   call Determine_RK_coefficients(RK_adv_coeff, RK_cum_coeff, irk_e, nrk)

   if(irk_e==1)then ! 1st step of Runge-Kutta method
      acc_newer_D_zeon=.true. !device side zeon will be changed in the following acc routine
      !$acc parallel loop gang vector present(zeon,zeon0,RK_cum_zeon)
      do m=1,me
         ! Note that zeon0(9,m) does not mean "mu", this is reserved for the other quantity
         ! zeon0(10,m) is the adiabatic term in the weight evolution @ t={n}, x={n}
         ! zeon(10,m) is a kind of temporary pouch to bring wdot term during shifte

         zeon0(1:5,m)=zeon(1:5,m)
         zeon0(7,m)=zeon(7,m)
         zeon0(8,m)=zeon(8,m) ! full-h

         RK_cum_zeon(:,m)=0.0
      enddo
   endif

   call n_T_u_group_e ! call n_T_u here to ensure that backTe is correctly updated

   if(electromagnetic==2) then
      call Get_Phi_gradPhi_gaHS_gradgaHS_at_zeon(Phi_at_zeon, gradPhi, gaHS, grad_gaH, grad_gaS)
   else
      call Get_Phi_and_gradPhi_at_zeon(Phi_at_zeon,gradPhi)
   endif

   ! NOTE: Structure Grid (SG) uses SI units for the Voltage so that it needs to be converted to GTS unit
   !       The SG also uses the flux coordinates for the space on which length unit [L] is defined along each direciton.
   !       Therefore, gaHS [L^(-1)] and other gradient operations [/L] do not need additional conversion.
   !       Only Voltage-related part needs to be converted
   !$acc parallel loop present(Phi_at_zeon,gradPhi)
   do m=1,me
      Phi_at_zeon(m) = Fac_Volt2GTS*Phi_at_zeon(m)
      !$acc loop seq
      do i=1,3
         gradPhi(i,m) = Fac_Volt2GTS*gradPhi(i,m)
      enddo
   enddo 


#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   !$acc parallel loop present(zeon,gaE,wgaE,gatE)
   do m=1,me
      call Get_extMP_alpha_amplitude(zeon(1,m),zeon(2,m),zeon(3,m),istep*tstep, gaE(m), wgaE(:,m), gatE(m))
   enddo
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pusheon, Get_extMP_alpha_amplitude", tt1-tt0," secs"
   call flush(stdout)
#endif

   !    prepare equilibrium Er=dPhi/da at particles
   !  ------------------------------------------------------------
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
   !$acc parallel loop present(zeon,zeon0)
   do m=1,me
      r=zeon(1,m)-a0
      r=r*dela
      ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
      fip=r-real(ip,wp)
      Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
   enddo
   !  ------------------------------------------------------------

#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   if(electromagnetic == 1)then
      if(Flag_GC_Solver==0) then
         m=gcmotion_eon_em_MGY(dgr,da,dgq,dgf,zeon,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,me,zeondim)
      elseif(Flag_GC_Solver==1) then
         m=gcmotion_eon_em_MGY_vec_form(dgr,da,dgq,dgf,zeon,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,1,me,zeondim)
         ! call eon_GC_FullMotion_vec_form(dgr,da,dgq,dgf,gradPhi,Er_particle,gatE)
      else
         write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
         call flush(stdout)
         STOP
      endif
   elseif(electromagnetic==2) then
      m=gcmotion_eon_em_MGY_vec_form_EM2(dgr,da,dgq,dgf,zeon,gradPhi,gaE,wgaE,gatE,gaHS,grad_gaH,grad_gaS,Er_particle,nonlinear,1,me,zeondim)
   else
      m=gcmotion_eon(dgr,da,dgq,dgf,zeon,gradPhi,Er_particle,nonlinear,me,zeondim)
   endif
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pusheon_wEvol0or1, electron gcmotion takes", tt1-tt0, " secs"
   call flush(stdout)
#endif

   !  call n_T_u_group_e  ! prepare n, T, u_//, dlnn/da, dlnT/da, etc at electrons.

   call Ti_at_electrons ! NOTE: array backT is borrowed to store Ti_at_electrons;
   ! WARNING: must make sure there is no memory conflict!!!
   ! backT is used in PUSHI and released afterward

   !     OUTPUT: dgr=drho/dt;        da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt;    dgf=dzeta/dt;
   !             ghq(:)=da/dt|_EXB;  Ta(:)=dtheta/dt|_EXB;  Pa(:)=drho/dt|_EXB;
   !             BB(:)=B;            Ba(:)=dB/da;           Bq(:)=dB/dtheta;
   !             gha(:)=1/|grad a|   F(:) =F
   !             wpi(3,:)=dphi/dzeta (NO change); no change for: zeon(1,2,4,:), wpi(1,2,:) and mu;



   IF(ifelectron==2.and.ifdphi==2) then  ! ifdphi==2 as the similar way as the ifdphi=2 for RK2 of GTS

      ! update GC position

#ifdef _MGY_TIMER_
      tt0=MPI_WTIME()
#endif
      acc_newer_D_zeon=.true.
      !$acc parallel loop gang vector present(zeon,zeon0,RK_cum_zeon,backT,backTe,dlnTdae,dlnndae)
      do m=1,me
         !     dbdz=0.0   for axisymmetric system   !wx
         mu=zeon(6,m)*zeon(6,m)        ! =normalized mu of electrons
         muB=mu*BB(m)                  ! =normalized perpendicular energy of electrons
         T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti
         !     energy=0.5_wp*for_E*(zeon(4,m)*BB(m))**2+mu*BB(m)

         vpara=zeon(4,m)*BB(m)
         v2=vpara**2

         !     energy=0.5_wp*for_E*(zeon(4,m)*BB(m))**2+muB
         energy=0.5_wp*for_E*v2+muB

         vpara=vpara*f_rho                ! =normalized v_//
         v2=vpara*vpara+muB*2.0_wp*m0ms   ! =v^2

         kappa=(min(umax*umax,energy*T_inverse)-1.5_wp)*(-dlnTdae(m))+(-dlnndae(m)) !=-dlnf_M/da

         if(ifneoclassical==1) then
            vdr=da(m)                 !  full da/dt including NC dynamics
            vdq=dgq(m)
            vdgf=dgf(m)
         else
            ! ExB drift in radial direction for w-dot and flux diagnostics
            vdr=ghq(m)             !  da/dt|_EXB
            vdq=Ta(m)
            vdgf=dve_gf(m)
         end if

         wdrive=vdr*kappa
         wgradB=vdr*Ba(m)+vdq*Bq(m)


         !    calculate <phi> at electrons (no FLR effect -- x=R)
         r=zeon(7,m)-s0       ! in radial coordinator s
         r=r*delr
         ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
         fip=r-real(ip,wp)
         Ea_00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons

         if ( Flag_Def_dne_adiab==0 ) then
            E_f=(Phi_at_zeon(m)-real(1.0-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! Adiabatic term (q*dphi_NZ/Te)
         elseif(Flag_Def_dne_adiab==1) then
            ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) for the New Poisson3D solver
            E_f=(Phi_at_zeon(m)+real(ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! Adiabatic term (q*dphi_tot/Te)
         elseif(Flag_Def_dne_adiab==2) then
            E_f=(Phi_at_zeon(m)+real(ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! (E_f=q*dphi_tot/Te). dne_Adiab=(exp(-E_f)-1)*f0
         else
            E_f=-1e100
            write(stdout,*) "Error!! Flag_Def_dne_adiab=",Flag_Def_dne_adiab, " is not defined in pushe_new.F90!!"
            STOP
         endif


         if(irk_e==1) then
            ! MGYOO(06/08/2023) : Save some information @ t{n}, x{n}
            ! zeon(9,m) = f_AD(t{n},x{n})*Z_i, where f_AD=-q*phi/Te*f0=-E_f : adiabatic weight
            ! zeon0(10,m)  = (q*dphi_NZ/Te) @ t{n}, x{n}
            now_wAD = -E_f*back_f0_e(m)*zeon0(6,m) 
            zeon(9,m) = now_wAD 
            zeon0(10,m) = E_f 
         endif

         if(ifneoclassical==1) then
            wdrift=(dgr(m)*BB(m)+zeon(4,m)*wgradB)*zeon(4,m)*BB(m)*for_E*T_inverse   ! see note
         else
            wdrift=(Pa(m)*BB(m)+zeon(4,m)*wgradB)*zeon(4,m)*BB(m)*for_E*T_inverse   ! see note
         end if

         wgradB=wgradB*mu*T_inverse   ! see note: additional term due to v_ExB*grad(f_M)

         wdot=(wdrive+wdrift+wgradB)

         !MGYOO (03/27/2019)
         ! Following evaluation might be wrong
         ! Please check again (dw/dt) equation
         ! MGYOO (12/15/2019):
         ! weight evolution of delta-h scheme:
         ! dw/dt=(1-w)*(-dln(f0)/dt)+(1-w)/(1-E_f)*(d(E_f)/dt)
         ! (1-E_f) is not needed for wdot, it is needed for wdphidt term
         ! wdot=(1.0_wp-E_f)*wdot

          now_w = zeon(5,m)

         !==================================================
         ! MGYOO(03/27/2019)
         ! Note that the following routine should be upgraded to calculate correct RK4
         !==================================================
         !==================================================
         !MGYOO(01/14/2019)
         !Note that f0 is toroidally symmetric so that grad_lnf0_gf=0 and we do not need vd_para_gf or vd_perp_gf terms to calculate weight evolution
         if(irk_e==4) then
            if(ifneoclassical==1) then
               ! !Caution!! Following routine is only available for neoclassical case which use full velocity.
               ! call esiget_Intp_B_vec_cov_and_cntv(zeon(1,m), zeon(2,m), zeon(3,m), istep*tstep, tot_B_vec_cov, tot_B_vec_cntv)
               ! vd_para_mag = (vdr*tot_B_vec_cov(1) + vdq*tot_B_vec_cov(2) + vdgf*tot_B_vec_cov(3))/BB(m)
               ! vd_para_a_cntv = vd_para_mag*(tot_B_vec_cntv(1)/BB(m))
               ! vd_para_gq_cntv = vd_para_mag*(tot_B_vec_cntv(2)/BB(m))
            else ! perturbed components only (ignoring neoclassical parts)
               vd_para_a_cntv = vpara_dve_a(m)
               vd_para_gq_cntv = vpara_dve_gq(m)
            endif
            vd_perp_a_cntv = vdr -vd_para_a_cntv
            vd_perp_gq_cntv = vdq -vd_para_gq_cntv

            !=================================================================
            !MGYOO (01/30/2019)
            ! grad_lnf0 includes (wdrive) terms
            ! wdot_accel includes (wdrift + wgradB) terms
            grad_lnf0_a_cov = -kappa
            grad_lnf0_gq_cov = 0.0

            wdot_vpara = -vd_para_a_cntv*grad_lnf0_a_cov -vd_para_gq_cntv*grad_lnf0_gq_cov
            wdot_vperp = -vd_perp_a_cntv*grad_lnf0_a_cov -vd_perp_gq_cntv*grad_lnf0_gq_cov
            wdot_accel  = wdrift + wgradB
            !=================================================================

            wdot_vpara = (1.0_wp-E_f)*wdot_vpara
            wdot_vperp = (1.0_wp-E_f)*wdot_vperp
            wdot_accel = (1.0_wp-E_f)*wdot_accel

            we_vpara(m) = we_vpara(m) + tstep*( (1.0-nonlinear*now_w)*wdot_vpara)
            we_vperp(m) = we_vperp(m) + tstep*( (1.0-nonlinear*now_w)*wdot_vperp)
            we_accel(m) = we_accel(m) + tstep*( (1.0-nonlinear*now_w)*wdot_accel)
         endif
         !=================================================


         ! Accumulate k-th RK terms (k1,k2,k3,k4) for the future value: f_(n+1) = f_(n) + (1/6)*(k1+2*k2+2*k3+k4)
         RK_cum_zeon(1,m) = RK_cum_zeon(1,m) +RK_cum_coeff*tstep*da(m)
         RK_cum_zeon(2,m) = RK_cum_zeon(2,m) +RK_cum_coeff*tstep*dgq(m)
         RK_cum_zeon(3,m) = RK_cum_zeon(3,m) +RK_cum_coeff*tstep*dgf(m)
         RK_cum_zeon(4,m) = RK_cum_zeon(4,m) +RK_cum_coeff*tstep*dgr(m)
         ! Accumulate weight including dphi/dt term (so the driving terms have (wdot+dphidt))
         RK_cum_zeon(7,m) = RK_cum_zeon(7,m) +RK_cum_coeff*tstep*(Tc/backT(m))**dsda_power*da(m)

         ! Advance particle position for the next RK step
         if(irk_e<4) then

            zeon(1,m) = zeon0(1,m)+RK_adv_coeff*tstep*da(m)
            zeon(2,m) = zeon0(2,m)+RK_adv_coeff*tstep*dgq(m)
            zeon(3,m) = zeon0(3,m)+RK_adv_coeff*tstep*dgf(m)
            zeon(4,m) = zeon0(4,m)+RK_adv_coeff*tstep*dgr(m)
            ! weight advanced including dphi/dt term (so the driving terms have (wdot+dphidt))
            zeon(7,m) = zeon0(7,m) +RK_adv_coeff*tstep*(Tc/backT(m))**dsda_power*da(m)

         else if(irk_e==4) then !irk_i==4
            ! Get the future value at time level (n+1): f_(n+1) = f_(n) + (1/6)*(k1+2*k2+2*k3+k4) = f_(n)+RK_cum_value
            ! zeon(1:5,m) = zeon0(1:5,m) +RK_cum_zeon(1:5,m)

            zeon(1:4,m) = zeon0(1:4,m) +RK_cum_zeon(1:4,m)
            zeon(7,m) = zeon0(7,m) +RK_cum_zeon(7,m)
         endif


         zeon(2,m)=modulo(zeon(2,m),twopi)
         zeon(3,m)=modulo(zeon(3,m),twopi)

         ! store GC information for flux measurements
         wpi(1,m)=ghq(m)    ! at present, wpi(1,*)=da/dt|_EXB
         wpi(2,m)=energy    ! normalized electron energy with u//=0 in Jan.09
         wpi(3,m)=muB       ! normalizd electron perpendicular energy

         if(irk_e==4) then
            UnNorm_theta_e(m) = UnNorm_theta_e(m) +RK_cum_zeon(2,m)
            UnNorm_zeta_e(m)  = UnNorm_zeta_e(m)  +RK_cum_zeon(3,m)
         endif

         ! MGYOO(05/06/2019) : To pass wdot information which are needed in the following for m=1:me routines after shift
         zeon(10,m) = wdot

      enddo


#ifdef _MGY_TIMER_
      tt1=MPI_WTIME()
      if(mype==0) write(stdout,'(a,ES11.3,a)') "In pusheon_wEvol0or1 main do loop takes", tt1-tt0, " secs"
      call flush(stdout)
#endif

      if(irk_e==1) then
         call cpu_time(cput0)
         tt0=MPI_WTIME()

         call Measure_and_Cumulate_eZonData_on_DG(wpi,da)

         call cpu_time(cput1)
         tt1=MPI_WTIME()
         time_in_Pushe(5) = time_in_Pushe(5) + (cput1-cput0)
         timewc_in_Pushe(5) = timewc_in_Pushe(5) + (tt1-tt0)
#ifdef _MGY_TIMER_
         if(mype==0) then
            write(stdout,*) "Measure_Cumulate_eFlux_a1D takes", tt1-tt0," secs"
            call flush(stdout)
         endif
#endif
      endif

      !Another Diagnostics on DG (Diagnostic Grid)
      if(irk_e==1 .and. mod(istep,msnap)==0) then
         call cpu_time(cput0)
         tt0=MPI_WTIME()

         !!!!!!$acc update self (Bq(:me),Ba(:me))  ! these are not used for now, but possibly for future diagnostics
         !$acc update self(BB(:me),Ta(:me),ghq(:me),gha(:me),Pa(:me),F(:me))
         !$acc update self(ee_inv_grad_a(:me),ee_inv_grad_gq(:me),ee_inv_grad_gf(:me))
         !$acc update self(dgr(:me), da(:me), dgq(:me), dgf(:me), wpi(:,:me))
         !$acc update self (dve_gf(:me),vpara_dve_a(:me),vpara_dve_gq(:me),vpara_dve_gf(:me))
         !$acc update self(UnNorm_theta_e(:me),UnNorm_zeta_e(:me))
         !$acc update self(we_init(:me),we_vpara(:me),we_vperp(:me),we_accel(:me),we_adiab(:me))
         !$acc update self(elabels(:me))
         !$acc update self(backTe(:me), dlnTdae(:me), backne(:me), dlnndae(:me))
         !$acc update self(back_f0_e(:me))
         !$acc update self(jteon0(:me),jteon1(:me),wzeon(:me),wpeon(:me),wteon0(:me),wteon1(:me))
         if(acc_newer_D_zeon) then
            !$acc update self(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_D_zeon=.false.
            acc_newer_H_zeon=.false.
         endif

         do m=1,me
            FullVelocity_da_dgq_dgf(1,m) = da(m)
            FullVelocity_da_dgq_dgf(2,m) = dgq(m)
            FullVelocity_da_dgq_dgf(3,m) = dgf(m)
         enddo
         call Measure_Electrons_on_Diagnostic_Grid(wpi(1:3,:),FullVelocity_da_dgq_dgf)

         call cpu_time(cput1)
         tt1=MPI_WTIME()
         time_in_Pushe(5) = time_in_Pushe(5) + (cput1-cput0)
         timewc_in_Pushe(5) = timewc_in_Pushe(5) + (tt1-tt0)
      endif

      call cpu_time(cput0)
      tt0=MPI_WTIME()
      call shifte
      call cpu_time(cput1)
      tt1=MPI_WTIME()
      time_in_Pushe(1) = time_in_Pushe(1) + (cput1-cput0)
      timewc_in_Pushe(1) = timewc_in_Pushe(1) + (tt1-tt0)


      ! call electron_locator
      call n_T_u_group_e

      
      call Get_Phi_at_zeon(Phi_at_zeon)
      !$acc parallel loop present(Phi_at_zeon)
      do m=1,me
         Phi_at_zeon(m) = Fac_Volt2GTS*Phi_at_zeon(m)
      enddo

      acc_newer_D_zeon=.true.
      !$acc parallel loop gang vector present(zeon,zeon0,backTe,tmp_zeon5,tmp_zeon8)
      do m=1,me
         ! calculate d(delta_phi)/dt at particles along trajactory

         tmp_zeon5(m)=zeon(5,m)
         tmp_zeon8(m)=zeon(8,m)

         T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti



         !    calculate <phi> at electrons (no FLR effect -- x=R)
         r=zeon(7,m)-s0       ! in radial coordinator s
         r=r*delr
         ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
         fip=r-real(ip,wp)
         Ea_00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons

         if ( Flag_Def_dne_adiab==0 ) then
            E_f=(Phi_at_zeon(m)-real(1.0-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! Adiabatic term (q*dphi_NZ/Te)
         elseif(Flag_Def_dne_adiab==1) then
            ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) for the New Poisson3D solver
            E_f=(Phi_at_zeon(m)+real(ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! Adiabatic term (q*dphi_tot/Te)
         elseif(Flag_Def_dne_adiab==2) then
            E_f=(Phi_at_zeon(m)+real(ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! (E_f=q*dphi_tot/Te). dne_Adiab=(exp(-E_f)-1)*f0
         else
            E_f=-1e100
            write(stdout,*) "Error!! Flag_Def_dne_adiab=",Flag_Def_dne_adiab, " is not defined in pushe_new.F90!!"
            STOP
         endif

         wdot = zeon(10,m)

         !======================================================
         ! MGYOO(05/10/2019): Here, zeon0(10,m) differs from that of the irk_e=1 stage of pusheon
         ! zeon0(10,m)= (q*dphi_NZ/Te) @ t{n}, x{n}, where {...} indicates the time level
         !======================================================
         ! Here, E_f = (q*dphi_NZ/Te) @ t{n}, x{n+k}, where k is RK_adv_coeff (eg. k=0.5 for irk_e=1,2 | k=1 for irk_e=3)
         ! MGYOO (12/30/2019):
         ! weight evolution of delta-h scheme:
         ! dw/dt=(1-w)*(-dln(f0)/dt)+(1-w)/(1-E_f)*(d(E_f)/dt)
         ! More exact form is the following:
         ! dz/dt= d(-ln(f0))/dt + d(-ln(|1-E_f|))/dt
         ! wdphidt = d(-ln(|1-E_f|))/dt
         wdphidt = -log(abs((1.0-E_f)/(1.0-zeon0(10,m))))/(RK_adv_coeff*tstep)

         if(Flag_Def_dne_adiab==2) then
            wdphidt = (E_f-zeon0(10,m))/(RK_adv_coeff*tstep)
         endif

         ! Note that full-h cannot be solved in the same way of dh because fAD can be zero.
         ! Instead, full-h will be solved directly from dh/dt=-d[f_AD]/dt in the similar way wEvol=3,4 do               
         ! wdot_fullh = -(now_wAD-old_wAD)/dtK
         now_wAD = -E_f*back_f0_e(m)*zeon0(6,m)
         wdot_fullh = int(ifelectron/2)*(-(now_wAD-zeon(9,m))/(RK_adv_coeff*tstep)) ! wdot_fullh=0 when ifelectron<2 

         ! weight advanced including dphi/dt term (so the driving terms have (wdot+dphidt))


         if(Flag_w_Evol==1) then
            ! From the definition of the transformed weight z=-ln(1-w),
            ! the weight increase becomes delta_w = (1-w0)*[1-exp(-F*dt)], where F is the driving force
            ! Note that this scheme is beneficial for nonlinear case only (dw/dt=(1-w)*F).
            ! For linear case ( dw/dt=F), this transform may cause additional numerical errors and accumulation.
            zeon(5,m) = zeon0(5,m) +(zeon0(9,m)-zeon0(5,m))*(1.0-exp(-(wdot+wdphidt)*(RK_adv_coeff*tstep)))
         else
            zeon(5,m) = zeon0(5,m) +RK_adv_coeff*tstep*((zeon0(9,m)-nonlinear*tmp_zeon5(m))*(wdot +wdphidt))
         endif
         zeon(8,m) = zeon0(8,m) +wdot_fullh*(RK_adv_coeff*tstep)
      enddo


      ! Save RK_cum_zeon(5) value to RK_cum_zeon(6) temporary
      !$acc parallel loop present(RK_cum_zeon)
      do m=1,me
         RK_cum_zeon(6,m) = RK_cum_zeon(5,m)
         RK_cum_zeon(9,m) = RK_cum_zeon(8,m)
      enddo

      if(Flag_Poisson==1) then
         do loop=1,MaxLoopNum_Adiab_Weight_Update
            call cpu_time(cput0)
            tt0=MPI_WTIME()
            if(Flag_Poisson3D_on_SG==1) then
               call Measure_Electron_Density_Current_on_SG(0)
            else
               call charge_eon
            endif
            call cpu_time(cput1)
            tt1=MPI_WTIME()
            time_in_Pushe(2) = time_in_Pushe(2) + (cput1-cput0)
            timewc_in_Pushe(2) = timewc_in_Pushe(2) + (tt1-tt0)

            prev_phi00=phi00

            call cpu_time(cput0)
            tt0=MPI_WTIME()
            call poisson( solver1, masssolver )
            call cpu_time(cput1)
            tt1=MPI_WTIME()
            time_in_Pushe(3) = time_in_Pushe(3) + (cput1-cput0)
            timewc_in_Pushe(3) = timewc_in_Pushe(3) + (tt1-tt0)

            phi00_convergency = sum(sqrt((phi00 - prev_phi00)**2.0))/sum(sqrt(prev_phi00**2.0))

            if (mype == 0) then
               open(unit=200903,file='phi00_every_istep_irk.dat', position='append',action='write')
               write(200903,'(I8,I8,I8,(999ES16.8))') istep, irk, loop, phi00_convergency, (E_00*c)/(charg0*1e8)*phi00(:)
               call flush(200903)
               close(200903)
            endif

            if(electromagnetic==2) then
               call cpu_time(cput0)
               tt0=MPI_WTIME()
               call Update_SG_AparaS(irk_e,nrk)
               call Solve_GK_Ampere_3D_on_SG
               call cpu_time(cput1)
               tt1=MPI_WTIME()
               time_in_Pushe(3) = time_in_Pushe(3) + (cput1-cput0)
               timewc_in_Pushe(3) = timewc_in_Pushe(3) + (tt1-tt0)
            endif

            call Get_Phi_at_zeon(Phi_at_zeon)
            !$acc parallel loop present(Phi_at_zeon)
            do m=1,me
               Phi_at_zeon(m) = Fac_Volt2GTS*Phi_at_zeon(m)
            enddo            

            acc_newer_D_zeon=.true.
            !$acc parallel loop gang vector present(zeon,zeon0,backTe,RK_cum_zeon)
            do m=1,me
               ! calculate d(delta_phi)/dt at particles along trajactory

               T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti

               !    calculate <phi> at electrons (no FLR effect -- x=R)
               r=zeon(7,m)-s0       ! in radial coordinator s
               r=r*delr
               ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
               fip=r-real(ip,wp)
               Ea_00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons

               if ( Flag_Def_dne_adiab==0 ) then
                  E_f=(Phi_at_zeon(m)-real(1.0-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! Adiabatic term (q*dphi_NZ/Te)
               elseif(Flag_Def_dne_adiab==1) then
                  ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) for the New Poisson3D solver
                  E_f=(Phi_at_zeon(m)+real(ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! Adiabatic term (q*dphi_tot/Te)
               elseif(Flag_Def_dne_adiab==2) then
                  E_f=(Phi_at_zeon(m)+real(ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! (E_f=q*dphi_tot/Te). dne_Adiab=(exp(-E_f)-1)*f0
               else
                  E_f=-1e100
                  write(stdout,*) "Error!! Flag_Def_dne_adiab=",Flag_Def_dne_adiab, " is not defined in pushe_new.F90!!"
                  STOP
               endif


               wdot = zeon(10,m)

               !======================================================
               ! MGYOO(05/10/2019): Here, zeon0(10,m) differs from that of the irk_e=1 stage of pusheon
               ! zeon0(10,m)= (q*dphi_NZ/Te) @ t{n}, x{n}, where {...} indicates the time level
               !======================================================
               ! Here, E_f = (q*dphi_NZ/Te) @ t{n}, x{n+k}, where k is RK_adv_coeff (eg. k=0.5 for irk_e=1,2 | k=1 for irk_e=3)
               ! MGYOO (12/30/2019):
               ! weight evolution of delta-h scheme:
               ! dw/dt=(1-w)*(-dln(f0)/dt)+(1-w)/(1-E_f)*(d(E_f)/dt)
               ! More exact form is the following:
               ! dz/dt= d(-ln(f0))/dt + d(-ln(|1-E_f|))/dt
               ! wdphidt = d(-ln(|1-E_f|))/dt
               wdphidt = -log(abs((1.0-E_f)/(1.0-zeon0(10,m))))/(RK_adv_coeff*tstep)

               if(Flag_Def_dne_adiab==2) then
                  wdphidt = (E_f-zeon0(10,m))/(RK_adv_coeff*tstep)
               endif

               ! Note that full-h cannot be solved in the same way of dh because fAD can be zero.
               ! Instead, full-h will be solved directly from dh/dt=-d[f_AD]/dt in the similar way wEvol=3,4 do               
               ! wdot_fullh = -(now_wAD-old_wAD)/dtK
               now_wAD = -E_f*back_f0_e(m)*zeon0(6,m)
               wdot_fullh = int(ifelectron/2)*(-(now_wAD-zeon(9,m))/(RK_adv_coeff*tstep))

               ! Advance particle information for the next RK step

               ! Calculate dADdt and evolve weight
               ! Accumulate weight including dphi/dt term (so the driving terms have (wdot+dphidt))
               if(Flag_w_Evol==1) then
                  ! Using the transformed equation, the weight increase becomes delta_w = (1-w0)*[1-exp(-F*dt)], where F is the driving force
                  ! RK_cum_zeon(5,m) = RK_cum_zeon(5,m) +(zeon0(9,m)-zeon0(5,m))*(1.0-exp(-(wdot+wdphidt)*(RK_cum_coeff*tstep)))
                  ! Note that we accumulate the increase of the Transformed wiehgt Z (Z=-ln(1-w); delta_Z=(F*dt)) rather than original weight.
                  RK_cum_zeon(5,m) = RK_cum_zeon(6,m) +(wdot+wdphidt)*(RK_cum_coeff*tstep)
               else
                  RK_cum_zeon(5,m) = RK_cum_zeon(6,m) +RK_cum_coeff*tstep*((zeon0(9,m)-nonlinear*tmp_zeon5(m))*(wdot +wdphidt))
               endif
               RK_cum_zeon(8,m) = RK_cum_zeon(9,m) +(wdot_fullh)*(RK_cum_coeff*tstep)

               ! Advance particle information for the next RK step
               if(irk_e<4) then
                  ! weight advanced including dphi/dt term (so the driving terms have (wdot+dphidt))
                  if(Flag_w_Evol==1) then
                     !Using the transformed equation, the weight increase becomes delta_w = (1-w0)*[1-exp(-F*dt)], where F is the driving force
                     zeon(5,m) = zeon0(5,m) +(zeon0(9,m)-zeon0(5,m))*(1.0-exp(-(wdot+wdphidt)*(RK_adv_coeff*tstep)))
                  else
                     zeon(5,m) = zeon0(5,m) +RK_adv_coeff*tstep*((zeon0(9,m)-nonlinear*tmp_zeon5(m))*(wdot +wdphidt))
                  endif
                  zeon(8,m) = zeon0(8,m) +wdot_fullh*(RK_adv_coeff*tstep)
               else if(irk_e==4) then !irk_i==4
                  if(Flag_w_Evol==1) then
                     ! Get the future value at time level (n+1): using accumulated increase of the transformed weight Z.
                     zeon(5,m) = zeon0(5,m) +(zeon0(9,m)-zeon0(5,m))*(1.0-exp(-RK_cum_zeon(5,m)))
                  else
                     ! Get the future value at time level (n+1): f_(n+1) = f_(n) + (1/6)*(k1+2*k2+2*k3+k4) = f_(n)+RK_cum_value
                     zeon(5,m) = zeon0(5,m) +RK_cum_zeon(5,m)
                  endif
                  zeon(8,m) = zeon0(8,m) +RK_cum_zeon(8,m)
               endif

               !==================================================
               !MGYOO(05/10/2019)
               ! Cal weight evolution due to the adiabatic response (d/dt(q*dphi/Te))
               if(irk_e==4)then
                  if(loop==MaxLoopNum_Adiab_Weight_Update) then
                     we_adiab(m) = we_adiab(m) + tstep*(zeon0(9,m)-nonlinear*tmp_zeon5(m))*wdphidt
                  endif
               endif
            enddo
         enddo
      endif


      if(irk_e==4)then
         if(Flag_Measure_Ptl_FootPrint==1) then
            !$acc update self(wpi(:,:me))
            call Measure_and_Cumulate_Electron_FootPrint_at_BDY1(FP_BDY1,wpi(2,:))
         endif
         call treat_electrons_out_of_boundary
         call Cumulate_e_lost_num

         ! ********************************
         if(mype==2) then
            error_8=1.0e-5_wp
            IFAIL=1
            rdum=D01AHF(0.0d0,zeon(1,me/2),error_8,minpts,acc,dsda,-1,IFAIL)
            print*, 'EEEEEEEEEEEE', zeon(7,me/2),'  ', rdum
            !   write(stdout,*)'######## electron gcmotion: wall =',(tt1-tt0),'  cpu =',(cput1-cput0)
         end if
         ! ********************************
      endif



   ELSE    !IF(ifelectron==1) delta-f scheme

#ifdef _MGY_TIMER_
      tt0=MPI_WTIME()
#endif

      acc_newer_D_zeon=.true.
      !$acc parallel loop gang vector present(zeon,zeon0,RK_cum_zeon,backT,backTe,dlnTdae,dlnndae)
      do m=1,me
         !     dbdz=0.0   for axisymmetric system   !wx
         mu=zeon(6,m)*zeon(6,m)        ! =normalized mu of electrons
         muB=mu*BB(m)                  ! =normalized perpendicular energy of electrons
         T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti
         !     energy=0.5_wp*for_E*(zeon(4,m)*BB(m))**2+mu*BB(m)

         vpara=zeon(4,m)*BB(m)
         v2=vpara**2

         !     energy=0.5_wp*for_E*(zeon(4,m)*BB(m))**2+muB
         energy=0.5_wp*for_E*v2+muB

         vpara=vpara*f_rho                ! =normalized v_//
         v2=vpara*vpara+muB*2.0_wp*m0ms   ! =v^2

         kappa=(min(umax*umax,energy*T_inverse)-1.5_wp)*(-dlnTdae(m))+(-dlnndae(m)) !=-dlnf_M/da

         if(ifneoclassical==1) then
            vdr=da(m)                 !  full da/dt including NC dynamics
            vdq=dgq(m)
            vdgf=dgf(m)
         else
            ! ExB drift in radial direction for w-dot and flux diagnostics
            vdr=ghq(m)             !  da/dt|_EXB
            vdq=Ta(m)
            vdgf=dve_gf(m)
         end if

         wdrive=vdr*kappa
         wgradB=vdr*Ba(m)+vdq*Bq(m)

         if(ifneoclassical==1) then
            wdrift=(dgr(m)*BB(m)+zeon(4,m)*wgradB)*zeon(4,m)*BB(m)*for_E*T_inverse   ! see note
         else
            wdrift=(Pa(m)*BB(m)+zeon(4,m)*wgradB)*zeon(4,m)*BB(m)*for_E*T_inverse   ! see note
         end if

         wgradB=wgradB*mu*T_inverse   ! see note: additional term due to v_ExB*grad(f_M)

         wdot=(wdrive+wdrift+wgradB)

         !==================================================
         ! MGYOO(03/27/2019)
         ! NOTE that the following routine should be upgraded to calculate correct RK4
         !==================================================
         !==================================================
         !MGYOO(01/14/2019)
         !Note that f0 is toroidally symmetric so that grad_lnf0_gf=0 and we do not need vd_para_gf or vd_perp_gf terms to calculate weight evolution
         if(irk_e==4) then
            if(ifneoclassical==1) then
               ! !Caution!! Following routine is only available for neoclassical case which use full velocity.
               ! call esiget_Intp_B_vec_cov_and_cntv(zeon(1,m), zeon(2,m), zeon(3,m), istep*tstep, tot_B_vec_cov, tot_B_vec_cntv)
               ! vd_para_mag = (vdr*tot_B_vec_cov(1) + vdq*tot_B_vec_cov(2) + vdgf*tot_B_vec_cov(3))/BB(m)
               ! vd_para_a_cntv = vd_para_mag*(tot_B_vec_cntv(1)/BB(m))
               ! vd_para_gq_cntv = vd_para_mag*(tot_B_vec_cntv(2)/BB(m))
            else ! perturbed components only (ignoring neoclassical parts)
               vd_para_a_cntv = vpara_dve_a(m)
               vd_para_gq_cntv = vpara_dve_gq(m)
            endif
            vd_perp_a_cntv = vdr -vd_para_a_cntv
            vd_perp_gq_cntv = vdq -vd_para_gq_cntv

            !=================================================================
            !MGYOO (01/30/2019)
            ! grad_lnf0 includes (wdrive) terms
            ! wdot_accel includes (wdrift + wgradB) terms
            grad_lnf0_a_cov = -kappa
            grad_lnf0_gq_cov = 0.0

            wdot_vpara = -vd_para_a_cntv*grad_lnf0_a_cov -vd_para_gq_cntv*grad_lnf0_gq_cov
            wdot_vperp = -vd_perp_a_cntv*grad_lnf0_a_cov -vd_perp_gq_cntv*grad_lnf0_gq_cov
            wdot_accel  = wdrift + wgradB
            !=================================================================

            we_vpara(m) = we_vpara(m) + tstep*( (1.0-nonlinear*zeon(5,m))*wdot_vpara)
            we_vperp(m) = we_vperp(m) + tstep*( (1.0-nonlinear*zeon(5,m))*wdot_vperp)
            we_accel(m) = we_accel(m) + tstep*( (1.0-nonlinear*zeon(5,m))*wdot_accel)
         endif
         !=================================================


         ! Accumulate k-th RK terms (k1,k2,k3,k4) for the future value: f_(n+1) = f_(n) + (1/6)*(k1+2*k2+2*k3+k4)
         RK_cum_zeon(1,m) = RK_cum_zeon(1,m) +RK_cum_coeff*tstep*da(m)
         RK_cum_zeon(2,m) = RK_cum_zeon(2,m) +RK_cum_coeff*tstep*dgq(m)
         RK_cum_zeon(3,m) = RK_cum_zeon(3,m) +RK_cum_coeff*tstep*dgf(m)
         RK_cum_zeon(4,m) = RK_cum_zeon(4,m) +RK_cum_coeff*tstep*dgr(m)
         RK_cum_zeon(7,m) = RK_cum_zeon(7,m) +RK_cum_coeff*tstep*(Tc/backT(m))**dsda_power*da(m)

         if(Flag_w_Evol==1) then
            ! Using the transformed equation, the weight increase becomes delta_w = (1-w0)*[1-exp(-F*dt)], where F is the driving force
            ! RK_cum_zeon(5,m) = RK_cum_zeon(5,m) +(zeon0(9,m)-zeon0(5,m))*(1.0-exp(-(wdot+wdphidt)*(RK_cum_coeff*tstep)))
            ! Note that we accumulate the increase of the Transformed wiehgt Z (Z=-ln(1-w); delta_Z=(F*dt)) rather than original weight.
            RK_cum_zeon(5,m) = RK_cum_zeon(5,m) +(wdot)*(RK_cum_coeff*tstep)
         else
            RK_cum_zeon(5,m) = RK_cum_zeon(5,m) +RK_cum_coeff*tstep*((zeon0(9,m)-nonlinear*zeon(5,m))*(wdot))
         endif
         ! Accumulate weight including dphi/dt term (so the driving terms have (wdot+dphidt))



         ! Advance particle position for the next RK step
         if(irk_e<4) then

            zeon(1,m) = zeon0(1,m)+RK_adv_coeff*tstep*da(m)
            zeon(2,m) = zeon0(2,m)+RK_adv_coeff*tstep*dgq(m)
            zeon(3,m) = zeon0(3,m)+RK_adv_coeff*tstep*dgf(m)
            zeon(4,m) = zeon0(4,m)+RK_adv_coeff*tstep*dgr(m)
            zeon(7,m) = zeon0(7,m) +RK_adv_coeff*tstep*(Tc/backT(m))**dsda_power*da(m)

            ! weight advanced including dphi/dt term (so the driving terms have (wdot+dphidt))
            if(Flag_w_Evol==1) then
               !Using the transformed equation, the weight increase becomes delta_w = (1-w0)*[1-exp(-F*dt)], where F is the driving force
               zeon(5,m) = zeon0(5,m) +(zeon0(9,m)-zeon0(5,m))*(1.0-exp(-(wdot)*(RK_adv_coeff*tstep)))
            else
               zeon(5,m) = zeon0(5,m) +RK_adv_coeff*tstep*((zeon0(9,m)-nonlinear*zeon(5,m))*(wdot))
            endif
         else if(irk_e==4) then !irk_i==4
            ! Get the future value at time level (n+1): f_(n+1) = f_(n) + (1/6)*(k1+2*k2+2*k3+k4) = f_(n)+RK_cum_value
            ! zeon(1:5,m) = zeon0(1:5,m) +RK_cum_zeon(1:5,m)

            zeon(1:4,m) = zeon0(1:4,m) +RK_cum_zeon(1:4,m)
            zeon(7,m) = zeon0(7,m) +RK_cum_zeon(7,m)

            if(Flag_w_Evol==1) then
               ! Get the future value at time level (n+1): using accumulated increase of the transformed weight Z.
               zeon(5,m) = zeon0(5,m) +(zeon0(9,m)-zeon0(5,m))*(1.0-exp(-RK_cum_zeon(5,m)))
            else
               ! Get the future value at time level (n+1): f_(n+1) = f_(n) + (1/6)*(k1+2*k2+2*k3+k4) = f_(n)+RK_cum_value
               zeon(5,m) = zeon0(5,m) +RK_cum_zeon(5,m)
            endif
         endif

         zeon(2,m)=modulo(zeon(2,m),twopi)
         zeon(3,m)=modulo(zeon(3,m),twopi)

         ! store GC information for flux measurements
         wpi(1,m)=ghq(m)    ! at present, wpi(1,*)=da/dt|_EXB
         wpi(2,m)=energy    ! normalized electron energy with u//=0 in Jan.09
         wpi(3,m)=muB       ! normalizd electron perpendicular energy

         if(irk_e==4) then
            UnNorm_theta_e(m) = UnNorm_theta_e(m) +RK_cum_zeon(2,m)
            UnNorm_zeta_e(m)  = UnNorm_zeta_e(m)  +RK_cum_zeon(3,m)
         endif

      enddo

#ifdef _MGY_TIMER_
      tt1=MPI_WTIME()
      if(mype==0) write(stdout,'(a,ES11.3,a)') "In pusheon_wEvol0or1 main do loop takes", tt1-tt0, " secs"
      call flush(stdout)
#endif



      if(irk_e==4)then
         ! out of boundary particle
         !$omp parallel do private(m)
         acc_newer_D_zeon=.true.
         !$acc parallel loop present(zeon,zeon0)
         do m=1,me
            if(zeon(1,m) > psimax)then
               zeon(1,m)=zeon0(1,m)
               zeon(2,m)=2.0_wp*pi-zeon0(2,m)
               zeon(3,m)=zeon0(3,m)
               zeon(4,m)=zeon0(4,m)
               zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
               !   zeon(5,m)=-zeon0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
               !   zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
               zeon(7,m)=zeon0(7,m)

            elseif(zeon(1,m) < psimin)then
               zeon(1,m)=zeon0(1,m)
               zeon(2,m)=2.0_wp*pi-zeon0(2,m)
               zeon(3,m)=zeon0(3,m)
               zeon(4,m)=zeon0(4,m)
               zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
               !   zeon(5,m)=-zeon0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
               !   zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
               zeon(7,m)=zeon0(7,m)
            endif
         enddo

         ! ********************************
         if(mype==2) then
            error_8=1.0e-5_wp
            IFAIL=1
            rdum=D01AHF(0.0d0,zeon(1,me/2),error_8,minpts,acc,dsda,-1,IFAIL)
            print*, 'EEEEEEEEEEEE', zeon(7,me/2),'  ', rdum
            !   write(stdout,*)'######## electron gcmotion: wall =',(tt1-tt0),'  cpu =',(cput1-cput0)
         end if
         ! ********************************

      endif



      if(irk_e==1) then
         call cpu_time(cput0)
         tt0=MPI_WTIME()

         call Measure_and_Cumulate_eZonData_on_DG(wpi,da)

         call cpu_time(cput1)
         tt1=MPI_WTIME()
         time_in_Pushe(5) = time_in_Pushe(5) + (cput1-cput0)
         timewc_in_Pushe(5) = timewc_in_Pushe(5) + (tt1-tt0)
#ifdef _MGY_TIMER_
         if(mype==0) then
            write(stdout,*) "Measure_Cumulate_eFlux_a1D takes", tt1-tt0," secs"
            call flush(stdout)
         endif
#endif
      endif


      !Another Diagnostics on DG (Diagnostic Grid)
      if(irk_e==1 .and. mod(istep,msnap)==0) then
         call cpu_time(cput0)
         tt0=MPI_WTIME()

         !!!!!!$acc update self (Bq(:me),Ba(:me))  ! these are not used for now, but possibly for future diagnostics
         !$acc update self(BB(:me),Ta(:me),ghq(:me),gha(:me),Pa(:me),F(:me))
         !$acc update self(ee_inv_grad_a(:me),ee_inv_grad_gq(:me),ee_inv_grad_gf(:me))
         !$acc update self(dgr(:me), da(:me), dgq(:me), dgf(:me), wpi(:,:me))
         !$acc update self (dve_gf(:me),vpara_dve_a(:me),vpara_dve_gq(:me),vpara_dve_gf(:me))
         !$acc update self(UnNorm_theta_e(:me),UnNorm_zeta_e(:me))
         !$acc update self(we_init(:me),we_vpara(:me),we_vperp(:me),we_accel(:me),we_adiab(:me))
         !$acc update self(elabels(:me))
         !$acc update self(backTe(:me), dlnTdae(:me), backne(:me), dlnndae(:me))
         !$acc update self(back_f0_e(:me))
         !$acc update self(jteon0(:me),jteon1(:me),wzeon(:me),wpeon(:me),wteon0(:me),wteon1(:me))
         if(acc_newer_D_zeon) then
            !$acc update self(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_D_zeon=.false.
            acc_newer_H_zeon=.false.
         endif

         do m=1,me
            FullVelocity_da_dgq_dgf(1,m) = da(m)
            FullVelocity_da_dgq_dgf(2,m) = dgq(m)
            FullVelocity_da_dgq_dgf(3,m) = dgf(m)
         enddo
         call Measure_Electrons_on_Diagnostic_Grid(wpi(1:3,:),FullVelocity_da_dgq_dgf)

         call cpu_time(cput1)
         tt1=MPI_WTIME()
         time_in_Pushe(5) = time_in_Pushe(5) + (cput1-cput0)
         timewc_in_Pushe(5) = timewc_in_Pushe(5) + (tt1-tt0)
#ifdef _MGY_TIMER_
      if(mype==0) then
         write(stdout,*) "Measure_Electrons_on_Diagnostic_Grid takes", tt1-tt0," secs"
         call flush(stdout)
      endif
#endif
      endif



   ENDIF   !IF(ifelectron==2.and.ifdphi==2) then  ! ifdphi==2 as the similar way as the ifdphi=2 for RK2 of GTS

end subroutine pusheon_wEvol0or1


subroutine pusheon_wEvol3(solver1, masssolver, irk_e, nrk)
   use precision
   use global_parameters
   use particle_array
   use mpi_defs
   use field_array
   use esi_interface
   use solver_data_class
   use Structured_Grid
   use Diagnostic_Grid
   use Ptl_Src_Loss
   use function_gg_module
   use MGY_PETSc_Interface_Data 

   implicit none

   type(solver_data) :: solver1, masssolver

   integer, intent(in) :: irk_e, nrk
   integer i,ip,m,ierr
   real(wp) RK_cum_coeff, RK_adv_coeff ! Coefficients to cumulate (to get future final value) or advance(to predict next RK step)
   real(wp) dgr(me),da(me),dgq(me),dgf(me),Phi_at_zeon(memax),gradPhi(3,me),wpi(5,me),Er_particle(me),gaHS(me),grad_gaH(3,me),grad_gaS(3,me)
   real(wp) gaE(me),wgaE(3,me), gatE(me) ! Magnetic Perturbation Information: alpha, grad_alpha, partial derivative d(gaE)/dt
   real(wp) tt0,tt1,cput0,cput1 ! For timer
   real(wp) timef,dela,r,fip
   integer gcmotion_eon, gcmotion_eon_em, gcmotion_eon_em_MGY, gcmotion_eon_em_MGY_vec_form, gcmotion_eon_em_MGY_vec_form_EM2
   real(wp) mu,Te_inverse,vpara,energy
   real(wp) AD ! Adiabatic term: AD=e*phi/Te0  ! Note that e here is the elementary charge which is the same as abs(qe)
   real(wp) Ea_00, inv_ds
   real(wp) new_f0, new_f_AD
   real(wp) FullVelocity_da_dgq_dgf(3,me) ! ds= a_dot*(Tc/T)**(dsda_power), dgq=gq_dot, dgf=dgf_dot
   real(wp) prev_phi00(0:mpsi),phi00_convergency
   integer loop

   !$acc declare create(dgr,da,dgq,dgf,Phi_at_zeon,gradPhi,wpi,Er_particle,gaE,wgaE,gatE,gaHS,grad_gaH,grad_gaS)


   ! if(ifneoclassical/=1 .or. Flag_w_Evol/=3) then
   !    write(stdout,*) "Error!! pusheon_wEvol3 is only available with ifneoclassical=1 .and. Flag_w_Evol=3"
   !    call flush(stdout)
   !    stop
   ! endif


   inv_ds = real(mpsi)/(s1-s0)


   if(acc_newer_H_zeon) then
      !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
      acc_newer_H_zeon=.false.
      acc_newer_D_zeon=.false.
   endif

   if(istep==istart .and. irk_e==1) then
      time_in_Pushe =0
      timewc_in_Pushe =0
   endif

   call Determine_RK_coefficients(RK_adv_coeff, RK_cum_coeff, irk_e, nrk)

   ! Update zeon0 if irk_e==1
   if(irk_e==1) then
      !$acc parallel loop gang vector present(zeon,zeon0,RK_cum_zeon)
      do m=1,me
         wpi(4,m) = zeon(5,m)-zeon0(5,m) ! Change of dh during the last one step
         wpi(5,m) = zeon(8,m)-zeon0(8,m) ! Change of full-h during the last one step

         zeon0(1:5,m)=zeon(1:5,m)
         zeon0(7,m)=zeon(7,m)
         zeon0(8,m)=zeon(8,m) ! full_f or full-h

         RK_cum_zeon(:,m)=0.0
      enddo
      acc_newer_D_zeon=.true. !device side zeon will be changed in the following acc routine
   endif


   call n_T_u_group_e

   if(electromagnetic==2) then
      call Get_Phi_gradPhi_gaHS_gradgaHS_at_zeon(Phi_at_zeon, gradPhi, gaHS, grad_gaH, grad_gaS)
   else
      call Get_Phi_and_gradPhi_at_zeon(Phi_at_zeon,gradPhi)
   endif

   ! NOTE: Structure Grid (SG) uses SI units for the Voltage so that it needs to be converted to GTS unit
   !       The SG also uses the flux coordinates for the space on which length unit [L] is defined along each direciton.
   !       Therefore, gaHS [L^(-1)] and other gradient operations [/L] do not need additional conversion.
   !       Only Voltage-related part needs to be converted
   !$acc parallel loop present(Phi_at_zeon,gradPhi)
   do m=1,me
      Phi_at_zeon(m) = Fac_Volt2GTS*Phi_at_zeon(m)
      !$acc loop seq
      do i=1,3
         gradPhi(i,m) = Fac_Volt2GTS*gradPhi(i,m)
      enddo
   enddo 


#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   !$acc parallel loop present(zeon,gaE,wgaE,gatE)
   do m=1,me
      call Get_extMP_alpha_amplitude(zeon(1,m),zeon(2,m),zeon(3,m),istep*tstep, gaE(m), wgaE(:,m), gatE(m))
   enddo
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pusheon, Get_MP_alpha_amplitude_arr_at_zeon", tt1-tt0," secs"
   call flush(stdout)
#endif


   timef=1.0_wp
   Er401_bp=Er401_bp0*timef
   Er_grid=Er_grid0*timef
   dela=1.0_wp/Ydr
   !$acc parallel loop present(zeon,zeon0)
   do m=1,me
      r=zeon(1,m)-a0
      r=r*dela
      ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
      fip=r-real(ip,wp)
      Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
   enddo
   !  ------------------------------------------------------------


#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   if(electromagnetic == 1)then
      if(Flag_GC_Solver==0) then
         m=gcmotion_eon_em_MGY(dgr,da,dgq,dgf,zeon,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,me,zeondim)
      elseif(Flag_GC_Solver==1) then
         m=gcmotion_eon_em_MGY_vec_form(dgr,da,dgq,dgf,zeon,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,1,me,zeondim)
      else
         write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
         call flush(stdout)
         STOP
      endif
   elseif(electromagnetic==2) then
      m=gcmotion_eon_em_MGY_vec_form_EM2(dgr,da,dgq,dgf,zeon,gradPhi,gaE,wgaE,gatE,gaHS,grad_gaH,grad_gaS,Er_particle,nonlinear,1,me,zeondim)
   else
      m=gcmotion_eon(dgr,da,dgq,dgf,zeon,gradPhi,Er_particle,nonlinear,me,zeondim)
   endif
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pusheon_wEvol3, electron gcmotion takes", tt1-tt0, " secs"
   call flush(stdout)
#endif


   call Ti_at_electrons ! NOTE: array backT is borrowed to store Ti_at_electrons;


   !     OUTPUT: dgr=drho/dt;        da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt;    dgf=dzeta/dt;
   !             ghq(:)=da/dt|_EXB;  Ta(:)=dtheta/dt|_EXB;  Pa(:)=drho/dt|_EXB;
   !             BB(:)=B;            Ba(:)=dB/da;           Bq(:)=dB/dtheta;
   !             gha(:)=1/|grad a|   F(:) =F
   !             wpi(3,:)=dphi/dzeta (NO change); no change for: zeon(1,2,4,:), wpi(1,2,:) and mu;



#ifdef _MGY_TIMER_
      tt0=MPI_WTIME()
#endif

   !$acc parallel loop gang vector present(zeon,zeon0,RK_cum_zeon,backT,backTe,dlnTdae,dlnndae)
   do m=1,me
      Te_inverse=1.0/(backTe(m)*msm0)
      mu=zeon(6,m)**2.0 ! normalized mu of electrons
      vpara = m0ms*csc0*BB(m)*zeon(4,m)
      energy = 0.5*msm0*vpara**2.0+mu*BB(m)


      if(irk_e==1) then
         ! Calculate Adiabatic term at irk_e=1 (so @(t_n,x_n) level)
         ! and save the value intto zeon0(8,m)

         ! calculate <phi> at electrons (no FLR effect -- x=R)
         r=zeon(7,m)-s0       ! in radial coordinator s
         r=r*inv_ds
         ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
         fip=r-real(ip,wp)
         Ea_00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons

         if ( Flag_Def_dne_adiab==0 ) then
            AD=(Phi_at_zeon(m)-real(1.0-ZF_SEPARATION,wp)*Ea_00)*Te_inverse*abs(csc0) ! Adiabatic term (e*dphi_NZ/Te)
         elseif(Flag_Def_dne_adiab==1) then
            ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) for the New Poisson3D solver
            AD=(Phi_at_zeon(m)+real(ZF_SEPARATION,wp)*Ea_00)*Te_inverse*abs(csc0) ! Adiabatic term (e*dphi_tot/Te)
         elseif(Flag_Def_dne_adiab==2) then
            AD=-1.0+exp((Phi_at_zeon(m)+real(ZF_SEPARATION,wp)*Ea_00)*Te_inverse*abs(csc0)) ! (AD=(-1+exp(e*dphi_tot/Te)))
         else
            AD=-1e100
            write(stdout,*) "Error!! Flag_Def_dne_adiab=",Flag_Def_dne_adiab, " is not defined in pushe_new.F90!!"
            STOP
         endif

         ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
         ! zeon0(10,m) is f_0 (@ irk=1)
         zeon(9,m)  = zeon0(6,m)*back_f0_e(m)*AD
         zeon0(10,m) = zeon0(6,m)*back_f0_e(m)
      endif


      ! =============================================================================================
      ! Cumulate (x,v,w) for preparing the final RK push&evolution
      RK_cum_zeon(1,m) = RK_cum_zeon(1,m) +RK_cum_coeff*tstep*da(m)
      RK_cum_zeon(2,m) = RK_cum_zeon(2,m) +RK_cum_coeff*tstep*dgq(m)
      RK_cum_zeon(3,m) = RK_cum_zeon(3,m) +RK_cum_coeff*tstep*dgf(m)
      RK_cum_zeon(4,m) = RK_cum_zeon(4,m) +RK_cum_coeff*tstep*dgr(m)
      RK_cum_zeon(7,m) = RK_cum_zeon(7,m) +RK_cum_coeff*tstep*(Tc/backT(m))**dsda_power*da(m)

      ! Advacne particles (x,v,w) for the next RK step
      if(irk_e<nrk) then
         zeon(1,m) = zeon0(1,m) +RK_adv_coeff*tstep*da(m)
         zeon(2,m) = zeon0(2,m) +RK_adv_coeff*tstep*dgq(m)
         zeon(3,m) = zeon0(3,m) +RK_adv_coeff*tstep*dgf(m)
         zeon(4,m) = zeon0(4,m) +RK_adv_coeff*tstep*dgr(m)
         zeon(7,m) = zeon0(7,m) +RK_adv_coeff*tstep*(Tc/backT(m))**dsda_power*da(m)
      else if(irk_e==nrk) then
         zeon(1:4,m) = zeon0(1:4,m) +RK_cum_zeon(1:4,m)
         zeon(7,m) = zeon0(7,m) +RK_cum_zeon(7,m)
      endif

      zeon(2,m)=modulo(zeon(2,m),twopi)
      zeon(3,m)=modulo(zeon(3,m),twopi)

      ! store GC information for flux measurements
      wpi(1,m)=ghq(m)    ! at present, wpi(1,*)=da/dt|_EXB
      wpi(2,m)=energy    ! normalized electron energy with u//=0 in Jan.09
      wpi(3,m)=mu*BB(m)       ! normalizd electron perpendicular energy


      if(irk_e==nrk) then
         UnNorm_theta_e(m) = UnNorm_theta_e(m) +RK_cum_zeon(2,m)
         UnNorm_zeta_e(m)  = UnNorm_zeta_e(m)  +RK_cum_zeon(3,m)
      endif
      !
      ! =============================================================================================
   enddo
   acc_newer_D_zeon=.true.


#ifdef _MGY_TIMER_
      tt1=MPI_WTIME()
      if(mype==0) write(stdout,'(a,ES11.3,a)') "In pusheon_wEvo3, main do loop takes", tt1-tt0, " secs"
      call flush(stdout)
#endif


   if(irk_e==1) then
      call cpu_time(cput0)
      tt0=MPI_WTIME()

      call Measure_and_Cumulate_eZonData_on_DG(wpi,da)

      call cpu_time(cput1)
      tt1=MPI_WTIME()
      time_in_Pushe(5) = time_in_Pushe(5) + (cput1-cput0)
      timewc_in_Pushe(5) = timewc_in_Pushe(5) + (tt1-tt0)
#ifdef _MGY_TIMER_
      if(mype==0) then
         write(stdout,*) "Measure_Cumulate_eFlux_a1D takes", tt1-tt0," secs"
         call flush(stdout)
      endif
#endif
   endif

   !Another Diagnostics on DG (Diagnostic Grid)
   if(irk_e==1 .and. mod(istep,msnap)==0) then
      call cpu_time(cput0)
      tt0=MPI_WTIME()

      !!!!!!$acc update self (Bq(:me),Ba(:me))  ! these are not used for now, but possibly for future diagnostics
      !$acc update self(BB(:me),Ta(:me),ghq(:me),gha(:me),Pa(:me),F(:me))
      !$acc update self(ee_inv_grad_a(:me),ee_inv_grad_gq(:me),ee_inv_grad_gf(:me))
      !$acc update self(dgr(:me), da(:me), dgq(:me), dgf(:me), wpi(:,:me))
      !$acc update self (dve_gf(:me),vpara_dve_a(:me),vpara_dve_gq(:me),vpara_dve_gf(:me))
      !$acc update self(UnNorm_theta_e(:me),UnNorm_zeta_e(:me))
      !$acc update self(we_init(:me),we_vpara(:me),we_vperp(:me),we_accel(:me),we_adiab(:me))
      !$acc update self(elabels(:me))
      !$acc update self(backTe(:me), dlnTdae(:me), backne(:me), dlnndae(:me))
      !$acc update self(back_f0_e(:me))
      !$acc update self(jteon0(:me),jteon1(:me),wzeon(:me),wpeon(:me),wteon0(:me),wteon1(:me))

      if(acc_newer_D_zeon) then
         !$acc update self(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
         acc_newer_D_zeon=.false.
         acc_newer_H_zeon=.false.
      endif

      do m=1,me
         FullVelocity_da_dgq_dgf(1,m) = da(m)
         FullVelocity_da_dgq_dgf(2,m) = dgq(m)
         FullVelocity_da_dgq_dgf(3,m) = dgf(m)
      enddo
      call Measure_Electrons_on_Diagnostic_Grid(wpi(1:3,:),FullVelocity_da_dgq_dgf)

      call cpu_time(cput1)
      tt1=MPI_WTIME()
      time_in_Pushe(5) = time_in_Pushe(5) + (cput1-cput0)
      timewc_in_Pushe(5) = timewc_in_Pushe(5) + (tt1-tt0)
#ifdef _MGY_TIMER_
      if(mype==0) then
         write(stdout,*) "Measure_Electrons_on_Diagnostic_Grid takes", tt1-tt0," secs"
         call flush(stdout)
      endif
#endif
   endif


   call cpu_time(cput0)
   tt0=MPI_WTIME()
   call shifte
   call cpu_time(cput1)
   tt1=MPI_WTIME()
   time_in_Pushe(1) = time_in_Pushe(1) + (cput1-cput0)
   timewc_in_Pushe(1) = timewc_in_Pushe(1) + (tt1-tt0)

   call n_T_u_group_e


   if(Flag_Poisson==0) then
      !$acc parallel loop present(zeon,zeon0,back_f0_e)
      do m=1,me
         AD=0.0 ! Because phi=0 (no Poisson)

         ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
         ! zeon0(10,m) is f_0 (@ irk=1)
         new_f_AD  = zeon0(6,m)*back_f0_e(m)*AD
         new_f0 = zeon0(6,m)*back_f0_e(m)

         zeon(5,m) = zeon0(5,m) - (new_f0-zeon0(10,m)) - int(ifelectron/2)*(new_f_AD-zeon(9,m)) ! df or dh
         zeon(8,m) = zeon0(8,m) - int(ifelectron/2)*(new_f_AD-zeon(9,m)) ! full-f or full-h
      enddo
      acc_newer_D_zeon=.true. !device side zeon will be changed in the following acc routine
   elseif(Flag_Poisson==1) then
      do loop=1,MaxLoopNum_Adiab_Weight_Update
         call cpu_time(cput0)
         tt0=MPI_WTIME()
            if(Flag_Poisson3D_on_SG==1) then
               call Measure_Electron_Density_Current_on_SG(0)
            else
               call charge_eon
            endif
         call cpu_time(cput1)
         tt1=MPI_WTIME()
         time_in_Pushe(2) = time_in_Pushe(2) + (cput1-cput0)
         timewc_in_Pushe(2) = timewc_in_Pushe(2) + (tt1-tt0)

         prev_phi00=phi00

         call cpu_time(cput0)
         tt0=MPI_WTIME()
         call poisson(solver1, masssolver)
         call cpu_time(cput1)
         tt1=MPI_WTIME()
         time_in_Pushe(3) = time_in_Pushe(3) + (cput1-cput0)
         timewc_in_Pushe(3) = timewc_in_Pushe(3) + (tt1-tt0)

         phi00_convergency = sum(sqrt((phi00 - prev_phi00)**2.0))/sum(sqrt(prev_phi00**2.0))

         if (mype == 0) then
            open(unit=200903,file='phi00_every_istep_irk.dat', position='append',action='write')
            write(200903,'(I8,I8,I8,(999ES16.8))') istep, irk, loop, phi00_convergency, (E_00*c)/(charg0*1e8)*phi00(:)
            call flush(200903)
            close(200903)
         endif

         if(electromagnetic==2) then
            call cpu_time(cput0)
            tt0=MPI_WTIME()
            call Update_SG_AparaS(irk_e,nrk)
            call Solve_GK_Ampere_3D_on_SG
            call cpu_time(cput1)
            tt1=MPI_WTIME()
            time_in_Pushe(3) = time_in_Pushe(3) + (cput1-cput0)
            timewc_in_Pushe(3) = timewc_in_Pushe(3) + (tt1-tt0)
         endif

         call Get_Phi_at_zeon(Phi_at_zeon)
         !$acc parallel loop present(Phi_at_zeon)
         do m=1,me
            Phi_at_zeon(m) = Fac_Volt2GTS*Phi_at_zeon(m)
         enddo 

         !$acc parallel loop gang vector present(zeon,zeon0,RK_cum_zeon,backT,backTe,dlnTdae,dlnndae)
         do m=1,me
            Te_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti

            ! Cal new phi at new position (i.e. phi @ (t_(n+1), x_(n+1)))


            ! calculate <phi> at electrons (no FLR effect -- x=R)
            r=zeon(7,m)-s0       ! in radial coordinator s
            r=r*inv_ds
            ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
            fip=r-real(ip,wp)
            Ea_00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons


            if ( Flag_Def_dne_adiab==0 ) then
               AD=(Phi_at_zeon(m)-real(1.0-ZF_SEPARATION,wp)*Ea_00)*Te_inverse*abs(csc0) ! Adiabatic term (e*dphi_NZ/Te)
            elseif(Flag_Def_dne_adiab==1) then
               ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) for the New Poisson3D solver
               AD=(Phi_at_zeon(m)+real(ZF_SEPARATION,wp)*Ea_00)*Te_inverse*abs(csc0) ! Adiabatic term (e*dphi_tot/Te)
            elseif(Flag_Def_dne_adiab==2) then
               AD=-1.0+exp((Phi_at_zeon(m)+real(ZF_SEPARATION,wp)*Ea_00)*Te_inverse*abs(csc0)) ! (AD=(-1+exp(e*dphi_tot/Te)))
            else
               AD=-1e100
               write(stdout,*) "Error!! Flag_Def_dne_adiab=",Flag_Def_dne_adiab, " is not defined in pushe_new.F90!!"
               STOP
            endif

            ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
            ! zeon0(10,m) is f_0 (@ irk=1)
            new_f_AD  = zeon0(6,m)*back_f0_e(m)*AD
            new_f0 = zeon0(6,m)*back_f0_e(m)

            zeon(5,m) = zeon0(5,m) - (new_f0-zeon0(10,m)) - int(ifelectron/2)*(new_f_AD-zeon(9,m)) ! df or dh
            zeon(8,m) = zeon0(8,m) - int(ifelectron/2)*(new_f_AD-zeon(9,m)) ! full-f or full-h
         enddo
         acc_newer_D_zeon=.true. !device side zeon will be changed in the following acc routine
      enddo
   endif


   if(irk_e==nrk)then
      if(Flag_Measure_Ptl_FootPrint==1) then
         !$acc update self(wpi(:,:me))
         call Measure_and_Cumulate_Electron_FootPrint_at_BDY1(FP_BDY1,wpi(2,:))
      endif
      call treat_electrons_out_of_boundary
      call Cumulate_e_lost_num
   endif


end subroutine pusheon_wEvol3



subroutine pusheon_wEvol4(solver1, masssolver, irk_e, nrk)
   use precision
   use global_parameters
   use particle_array
   use mpi_defs
   use field_array
   use esi_interface
   use solver_data_class
   use Structured_Grid
   use Diagnostic_Grid
   use Ptl_Src_Loss
   use function_gg_module
   use MGY_PETSc_Interface_Data

   implicit none

   type(solver_data) :: solver1, masssolver

   integer, intent(in) :: irk_e, nrk
   integer i,ip,m,ierr
   real(wp) RK_cum_coeff, RK_adv_coeff ! Coefficients to cumulate (to get future final value) or advance(to predict next RK step)
   real(wp) dgr(me),da(me),dgq(me),dgf(me),Phi_at_zeon(memax),gradPhi(3,me),gaHS(me),grad_gaH(3,me),grad_gaS(3,me),wpi(5,me),Er_particle(me)
   real(wp) gaE(me),wgaE(3,me), gatE(me) ! Magnetic Perturbation Information: alpha, grad_alpha, partial derivative d(gaE)/dt
   real(wp) tt0,tt1,cput0,cput1 ! For timer
   real(wp) timef,dela,r,fip
   integer gcmotion_eon, gcmotion_eon_em, gcmotion_eon_em_MGY, gcmotion_eon_em_MGY_vec_form, gcmotion_eon_em_MGY_vec_form_EM2
   real(wp) mu,Te_inverse,vpara,energy
   real(wp) AD ! Adiabatic term: AD=e*phi/Te0  ! Note that e here is the elementary charge which is the same as abs(qe)
   real(wp) inv_ds
   real(wp) new_f0, new_f_AD
   real(wp) FullVelocity_da_dgq_dgf(3,me) ! ds= a_dot*(Tc/T)**(dsda_power), dgq=gq_dot, dgf=dgf_dot
   real(wp) prev_phi00(0:mpsi),phi00_convergency
   integer loop

   !$acc declare create(dgr,da,dgq,dgf,Phi_at_zeon,gradPhi,gaHS,grad_gaH,grad_gaS,wpi,Er_particle,gaE,wgaE,gatE)


   ! if(ifneoclassical/=1 .or. Flag_w_Evol/=3) then
   !    write(stdout,*) "Error!! pusheon_wEvol3 is only available with ifneoclassical=1 .and. Flag_w_Evol=3"
   !    call flush(stdout)
   !    stop
   ! endif


   inv_ds = real(mpsi)/(s1-s0)


   if(acc_newer_H_zeon) then
      !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
      acc_newer_H_zeon=.false.
      acc_newer_D_zeon=.false.
   endif

   if(istep==istart .and. irk_e==1) then
      time_in_Pushe =0
      timewc_in_Pushe =0
   endif

   call Determine_RK_coefficients(RK_adv_coeff, RK_cum_coeff, irk_e, nrk)

   ! Update zeon0 if irk_e==1
   if(irk_e==1) then
      !$acc parallel loop gang vector present(zeon,zeon0,RK_cum_zeon)
      do m=1,me
         wpi(4,m) = zeon(5,m)-zeon0(5,m) ! Change of dh during the last one step
         wpi(5,m) = zeon(8,m)-zeon0(8,m) ! Change of full-h during the last one step
         zeon0(1:5,m)=zeon(1:5,m)
         zeon0(7,m)=zeon(7,m)
         zeon0(8,m)=zeon(8,m) ! full_f or full-h

         RK_cum_zeon(:,m)=0.0
      enddo
      acc_newer_D_zeon=.true. !device side zeon will be changed in the following acc routine
   endif

   call n_T_u_group_e ! call n_T_u here to ensure that backTe is correctly updated
   
   if(electromagnetic==2) then
      call Get_Phi_gradPhi_gaHS_gradgaHS_at_zeon(Phi_at_zeon, gradPhi, gaHS, grad_gaH, grad_gaS)
   else
      call Get_Phi_and_gradPhi_at_zeon(Phi_at_zeon,gradPhi)
   endif

   ! NOTE: Structure Grid (SG) uses SI units for the Voltage so that it needs to be converted to GTS unit
   !       The SG also uses the flux coordinates for the space on which length unit [L] is defined along each direciton.
   !       Therefore, gaHS [L^(-1)] and other gradient operations [/L] do not need additional conversion.
   !       Only Voltage-related part needs to be converted
   !$acc parallel loop present(Phi_at_zeon,gradPhi)
   do m=1,me
      Phi_at_zeon(m) = Fac_Volt2GTS*Phi_at_zeon(m)
      !$acc loop seq
      do i=1,3
         gradPhi(i,m) = Fac_Volt2GTS*gradPhi(i,m)
      enddo
   enddo


#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   !$acc parallel loop present(zeon,gaE,wgaE,gatE)
   do m=1,me
      call Get_extMP_alpha_amplitude(zeon(1,m),zeon(2,m),zeon(3,m),istep*tstep, gaE(m), wgaE(:,m), gatE(m))
   enddo
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pusheon_wEvol4, Get_MP_alpha_amplitude_arr_at_zeon", tt1-tt0," secs"
   call flush(stdout)
#endif



   timef=1.0_wp
   Er401_bp=Er401_bp0*timef
   Er_grid=Er_grid0*timef
   dela=1.0_wp/Ydr
   !$acc parallel loop present(zeon,zeon0)
   do m=1,me
      r=zeon(1,m)-a0
      r=r*dela
      ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
      fip=r-real(ip,wp)
      Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
   enddo
   !  ------------------------------------------------------------



#ifdef _MGY_TIMER_
   tt0=MPI_WTIME()
#endif
   if(electromagnetic == 1)then
      if(Flag_GC_Solver==0) then
         m=gcmotion_eon_em_MGY(dgr,da,dgq,dgf,zeon,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,me,zeondim)
      elseif(Flag_GC_Solver==1) then
         m=gcmotion_eon_em_MGY_vec_form(dgr,da,dgq,dgf,zeon,gradPhi,gaE,wgaE,gatE,Er_particle,nonlinear,1,me,zeondim)
      else
         write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
         call flush(stdout)
         STOP
      endif
   elseif(electromagnetic==2) then
      m=gcmotion_eon_em_MGY_vec_form_EM2(dgr,da,dgq,dgf,zeon,gradPhi,gaE,wgaE,gatE,gaHS,grad_gaH,grad_gaS,Er_particle,nonlinear,1,me,zeondim)
   else
      m=gcmotion_eon(dgr,da,dgq,dgf,zeon,gradPhi,Er_particle,nonlinear,me,zeondim)
   endif
#ifdef _MGY_TIMER_
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,ES11.3,a)') "In pusheon_wEvol4, electron gcmotion takes", tt1-tt0, " secs"
   call flush(stdout)
#endif



   call Ti_at_electrons ! NOTE: array backT is borrowed to store Ti_at_electrons;



   !     OUTPUT: dgr=drho/dt;        da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt;    dgf=dzeta/dt;
   !             ghq(:)=da/dt|_EXB;  Ta(:)=dtheta/dt|_EXB;  Pa(:)=drho/dt|_EXB;
   !             BB(:)=B;            Ba(:)=dB/da;           Bq(:)=dB/dtheta;
   !             gha(:)=1/|grad a|   F(:) =F
   !             wpi(3,:)=dphi/dzeta (NO change); no change for: zeon(1,2,4,:), wpi(1,2,:) and mu;



#ifdef _MGY_TIMER_
      tt0=MPI_WTIME()
#endif

   !$acc parallel loop gang vector present(zeon,zeon0,RK_cum_zeon,backT,backTe,dlnTdae,dlnndae)
   do m=1,me
      Te_inverse=1.0/(backTe(m)*msm0)
      mu=zeon(6,m)**2.0 ! normalized mu of electrons
      vpara = m0ms*csc0*BB(m)*zeon(4,m)
      energy = 0.5*msm0*vpara**2.0+mu*BB(m)


      if(irk_e==1) then
         ! Calculate Adiabatic term at irk_e=1 (so @(t_n,x_n) level)
         ! and save the value intto zeon0(8,m)

         AD = Phi_at_zeon(m)*Te_inverse*abs(csc0)


         ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
         ! zeon0(9,m) is f_0 (@ irk=1)
         zeon(9,m)  = zeon0(6,m)*back_f0_e(m)*AD
         zeon0(9,m) = zeon0(8,m)+int(ifelectron/2)*zeon(9,m) ! #TODO : IMPORTANT! checky why this re-normalization is needed. 


         ! zeon0(9,m) = zeon0(6,m)*back_f0_e(m)

         new_f_AD  = zeon0(6,m)*back_f0_e(m)*AD
         new_f0 = zeon0(6,m)*back_f0_e(m)




         zeon(5,m) = zeon0(9,m) - (new_f0+int(ifelectron/2)*new_f_AD) ! df or dh
         zeon(8,m) = zeon0(9,m) - (int(ifelectron/2)*new_f_AD) ! full-f or full-h


         if(mype==0 .and. m==1) then
            write(stdout,*) 'new_f_AD=', new_f_AD, '; zeon(9,m)=', zeon(9,m),'; zeon0(8)+zeon(9)=', zeon0(8,m)+zeon(9,m)
            write(stdout,*) 'zeon0(8,m)=',zeon0(8,m), '; ori_zeon(8,m)=', zeon0(8,m) - (new_f_AD-zeon(9,m))
            write(stdout,*) 'zeon0(9,m)=',zeon0(9,m), '; new_zeon(8,m)=', zeon0(9,m) - (int(ifelectron/2)*new_f_AD)
         endif

      endif


      ! =============================================================================================
      ! Cumulate (x,v,w) for preparing the final RK push&evolution
      RK_cum_zeon(1,m) = RK_cum_zeon(1,m) +RK_cum_coeff*tstep*da(m)
      RK_cum_zeon(2,m) = RK_cum_zeon(2,m) +RK_cum_coeff*tstep*dgq(m)
      RK_cum_zeon(3,m) = RK_cum_zeon(3,m) +RK_cum_coeff*tstep*dgf(m)
      RK_cum_zeon(4,m) = RK_cum_zeon(4,m) +RK_cum_coeff*tstep*dgr(m)
      RK_cum_zeon(7,m) = RK_cum_zeon(7,m) +RK_cum_coeff*tstep*(Tc/backT(m))**dsda_power*da(m)

      ! Advacne particles (x,v,w) for the next RK step
      if(irk_e<nrk) then
         zeon(1,m) = zeon0(1,m) +RK_adv_coeff*tstep*da(m)
         zeon(2,m) = zeon0(2,m) +RK_adv_coeff*tstep*dgq(m)
         zeon(3,m) = zeon0(3,m) +RK_adv_coeff*tstep*dgf(m)
         zeon(4,m) = zeon0(4,m) +RK_adv_coeff*tstep*dgr(m)
         zeon(7,m) = zeon0(7,m) +RK_adv_coeff*tstep*(Tc/backT(m))**dsda_power*da(m)
      else if(irk_e==nrk) then
         zeon(1:4,m) = zeon0(1:4,m) +RK_cum_zeon(1:4,m)
         zeon(7,m) = zeon0(7,m) +RK_cum_zeon(7,m)
      endif

      zeon(2,m)=modulo(zeon(2,m),twopi)
      zeon(3,m)=modulo(zeon(3,m),twopi)

      ! store GC information for flux measurements
      wpi(1,m)=ghq(m)    ! at present, wpi(1,*)=da/dt|_EXB
      wpi(2,m)=energy    ! normalized electron energy with u//=0 in Jan.09
      wpi(3,m)=mu*BB(m)       ! normalizd electron perpendicular energy


      if(irk_e==nrk) then
         UnNorm_theta_e(m) = UnNorm_theta_e(m) +RK_cum_zeon(2,m)
         UnNorm_zeta_e(m)  = UnNorm_zeta_e(m)  +RK_cum_zeon(3,m)
      endif
      !
      ! =============================================================================================
   enddo
   acc_newer_D_zeon=.true.




#ifdef _MGY_TIMER_
      tt1=MPI_WTIME()
      if(mype==0) write(stdout,'(a,ES11.3,a)') "In pusheon_wEvol4, main do loop takes", tt1-tt0, " secs"
      call flush(stdout)
#endif


   if(irk_e==1) then
      call cpu_time(cput0)
      tt0=MPI_WTIME()


      call Measure_and_Cumulate_eZonData_on_DG(wpi,da)


      call cpu_time(cput1)
      tt1=MPI_WTIME()
      time_in_Pushe(5) = time_in_Pushe(5) + (cput1-cput0)
      timewc_in_Pushe(5) = timewc_in_Pushe(5) + (tt1-tt0)
#ifdef _MGY_TIMER_
      if(mype==0) then
         write(stdout,*) "Measure_Cumulate_eFlux_a1D takes", tt1-tt0," secs"
         call flush(stdout)
      endif
#endif
   endif

   !Another Diagnostics on DG (Diagnostic Grid)
   if(irk_e==1 .and. mod(istep,msnap)==0) then
      call cpu_time(cput0)
      tt0=MPI_WTIME()


      !!!!!!$acc update self (Bq(:me),Ba(:me))  ! these are not used for now, but possibly for future diagnostics
      !$acc update self(BB(:me),Ta(:me),ghq(:me),gha(:me),Pa(:me),F(:me))
      !$acc update self(ee_inv_grad_a(:me),ee_inv_grad_gq(:me),ee_inv_grad_gf(:me))
      !$acc update self(dgr(:me), da(:me), dgq(:me), dgf(:me), wpi(:,:me))
      !$acc update self (dve_gf(:me),vpara_dve_a(:me),vpara_dve_gq(:me),vpara_dve_gf(:me))
      !$acc update self(UnNorm_theta_e(:me),UnNorm_zeta_e(:me))
      !$acc update self(we_init(:me),we_vpara(:me),we_vperp(:me),we_accel(:me),we_adiab(:me))
      !$acc update self(elabels(:me))
      !$acc update self(backTe(:me), dlnTdae(:me), backne(:me), dlnndae(:me))
      !$acc update self(back_f0_e(:me))
      !$acc update self(jteon0(:me),jteon1(:me),wzeon(:me),wpeon(:me),wteon0(:me),wteon1(:me))

      if(acc_newer_D_zeon) then
         !$acc update self(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
         acc_newer_D_zeon=.false.
         acc_newer_H_zeon=.false.
      endif

      do m=1,me
         FullVelocity_da_dgq_dgf(1,m) = da(m)
         FullVelocity_da_dgq_dgf(2,m) = dgq(m)
         FullVelocity_da_dgq_dgf(3,m) = dgf(m)
      enddo


      call Measure_Electrons_on_Diagnostic_Grid(wpi(1:3,:),FullVelocity_da_dgq_dgf)

      call cpu_time(cput1)
      tt1=MPI_WTIME()
      time_in_Pushe(5) = time_in_Pushe(5) + (cput1-cput0)
      timewc_in_Pushe(5) = timewc_in_Pushe(5) + (tt1-tt0)
#ifdef _MGY_TIMER_
      if(mype==0) then
         write(stdout,*) "Measure_Electrons_on_Diagnostic_Grid takes", tt1-tt0," secs"
         call flush(stdout)
      endif
#endif
   endif


   call cpu_time(cput0)
   tt0=MPI_WTIME()
   call shifte
   call cpu_time(cput1)
   tt1=MPI_WTIME()
   time_in_Pushe(1) = time_in_Pushe(1) + (cput1-cput0)
   timewc_in_Pushe(1) = timewc_in_Pushe(1) + (tt1-tt0)

   ! call electron_locator
   call n_T_u_group_e


   ! !$acc update self(wzeon,wteon0,wteon1,jteon0,jteon1)
   ! do m=1,me
   !    if(wzeon(m)<0.0 .or. wzeon(m) >1.0) then
   !      write(0,'((A,I4,A,ES12.4,A,I8,A,10ES12.4))') "11 mype=",mype,"; wzeon(m)=",wzeon(m), "; m=",m,"; zeon(:,m)=",zeon(:,m)
   !    endif
   !    if(jteon0(m)==20050 .or. jteon1(m)==20050)then
   !      write(0,'((A,I4,A,I4,A,I8,A,I4))') "11 mype=",mype,"; jteon0(m)=",jteon0(m), "; m=",m,"; jteon1(m)=",jteon1(m)
   !    endif
   !    if(wteon0(m)<0.0 .or. wteon0(m)>1.0 .or. wteon1(m)<0.0 .or. wteon1(m)>1.0) then
   !      write(0,'((A,I4,A,ES12.4,A,I8,A,10ES12.4))') "11 mype=",mype,"; wteon0(m)=",wteon0(m), "; m=",m,"; wteon1(m)=",wteon1(m)
   !    endif

   ! enddo

   if(Flag_Poisson==0) then
      !$acc parallel loop present(zeon,zeon0,back_f0_e)
      do m=1,me
         AD=0.0 ! Because phi=0 (no Poisson)

         ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
         ! zeon0(9,m) is f_0 (@ irk=1)
         new_f_AD  = zeon0(6,m)*back_f0_e(m)*AD
         new_f0 = zeon0(6,m)*back_f0_e(m)

         ! zeon(5,m) = zeon0(5,m) - (new_f0-zeon0(9,m)) - (new_f_AD-zeon(9,m)) ! df or dh
         ! zeon(8,m) = zeon0(8,m) - (new_f_AD-zeon(9,m)) ! full-f or full-h
         zeon(5,m) = zeon0(9,m) - (new_f0+int(ifelectron/2)*new_f_AD) ! df or dh
         zeon(8,m) = zeon0(9,m) - (int(ifelectron/2)*new_f_AD) ! full-f or full-h
      enddo
      acc_newer_D_zeon=.true. !device side zeon will be changed in the following acc routine
   elseif(Flag_Poisson==1) then
      do loop=1,MaxLoopNum_Adiab_Weight_Update
         call cpu_time(cput0)
         tt0=MPI_WTIME()
            if(Flag_Poisson3D_on_SG==1) then
               call Measure_Electron_Density_Current_on_SG(0)
            else
               call charge_eon
            endif
         call cpu_time(cput1)
         tt1=MPI_WTIME()
         time_in_Pushe(2) = time_in_Pushe(2) + (cput1-cput0)
         timewc_in_Pushe(2) = timewc_in_Pushe(2) + (tt1-tt0)

         prev_phi00=phi00

         call cpu_time(cput0)
         tt0=MPI_WTIME()
         call poisson(solver1, masssolver)
         call cpu_time(cput1)
         tt1=MPI_WTIME()
         time_in_Pushe(3) = time_in_Pushe(3) + (cput1-cput0)
         timewc_in_Pushe(3) = timewc_in_Pushe(3) + (tt1-tt0)

         phi00_convergency = sum(sqrt((phi00 - prev_phi00)**2.0))/sum(sqrt(prev_phi00**2.0))

         if (mype == 0) then
            open(unit=200903,file='phi00_every_istep_irk.dat', position='append',action='write')
            write(200903,'(I8,I8,I8,(999ES16.8))') istep, irk, loop, phi00_convergency, (E_00*c)/(charg0*1e8)*phi00(:)
            call flush(200903)
            close(200903)
         endif


         if(electromagnetic==2) then
            call cpu_time(cput0)
            tt0=MPI_WTIME()
            call Update_SG_AparaS(irk_e,nrk)
            call Solve_GK_Ampere_3D_on_SG
            call cpu_time(cput1)
            tt1=MPI_WTIME()
            time_in_Pushe(3) = time_in_Pushe(3) + (cput1-cput0)
            timewc_in_Pushe(3) = timewc_in_Pushe(3) + (tt1-tt0)
         endif

         call Get_Phi_at_zeon(Phi_at_zeon)
         !$acc parallel loop present(Phi_at_zeon)
         do m=1,me
            Phi_at_zeon(m) = Fac_Volt2GTS*Phi_at_zeon(m)
         enddo


         !$acc parallel loop present(zeon,zeon0,RK_cum_zeon,backT,backTe,dlnTdae,dlnndae)
         do m=1,me
            Te_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti

            AD = Phi_at_zeon(m)*Te_inverse*abs(csc0)

            ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
            ! zeon0(9,m) is f_0 (@ irk=1)
            new_f_AD  = zeon0(6,m)*back_f0_e(m)*AD
            new_f0 = zeon0(6,m)*back_f0_e(m)

            if(mype==0 .and. m==1) then
               write(stdout,*) 'new_f_AD=', new_f_AD, '; zeon(9,m)=', zeon(9,m),'; zeon0(8)+zeon(9)=', zeon0(8,m)+zeon(9,m)
               write(stdout,*) 'zeon0(8,m)=',zeon0(8,m), '; ori_zeon(8,m)=', zeon0(8,m) - (new_f_AD-zeon(9,m))
               write(stdout,*) 'zeon0(9,m)=',zeon0(9,m), '; new_zeon(8,m)=', zeon0(9,m) - (int(ifelectron/2)*new_f_AD)
            endif

            ! zeon(5,m) = zeon0(5,m) - (new_f0-zeon0(9,m)) - (new_f_AD-zeon(9,m)) ! df or dh
            ! zeon(8,m) = zeon0(8,m) - (new_f_AD-zeon(9,m)) ! full-f or full-h
            zeon(5,m) = zeon0(9,m) - (new_f0+int(ifelectron/2)*new_f_AD) ! df or dh
            zeon(8,m) = zeon0(9,m) - (int(ifelectron/2)*new_f_AD) ! full-f or full-h

            ! zeon(5,m) = zeon0(5,m)
            ! zeon(8,m) = zeon0(8,m)
         enddo
         acc_newer_D_zeon=.true. !device side zeon will be changed in the following acc routine
      enddo
   endif


   if(irk_e==nrk)then
      if(Flag_Measure_Ptl_FootPrint==1) then
         !$acc update self(wpi(:,:me))
         call Measure_and_Cumulate_Electron_FootPrint_at_BDY1(FP_BDY1,wpi(2,:))
      endif
      call treat_electrons_out_of_boundary
      call Cumulate_e_lost_num
   endif


end subroutine pusheon_wEvol4



subroutine electron_locator
   use precision
   use global_parameters
   use particle_array
   use field_array
   implicit none

   integer m,im,ii,j01,j00
   real(wp) rdum,tdum,delr,delt(0:mpsi),delz,wz1,wp1,     &
      pi2_inv,psitmp,thetatmp,zetatmp,psimax,psimin


   delr=1.0_wp/deltar
   delt=2.0_wp*pi/deltat
   delz=1.0_wp/deltaz
   pi2_inv=0.5_wp/pi
   psimin=s_grid(0)
   psimax=s_grid(mpsi)

   !$omp parallel do private(m,psitmp,thetatmp,zetatmp,&
   !$omp& rdum,ii,wp1,im,tdum,j00,j01)

   if(acc_newer_H_zeon) then
      !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
      acc_newer_H_zeon=.false.
      acc_newer_D_zeon=.false.
   endif

   !$acc parallel loop present(zeon,jteon0,jteon1,wzeon,wpeon,wteon0,wteon1)
   do m=1,me
      !     psitmp=zeon(1,m)      ! =a
      psitmp=zeon(7,m)       ! =s
      thetatmp=zeon(2,m)
      zetatmp=zeon(3,m)

      wzeon(m)=(zetatmp-zetamin)*delz

      rdum=psitmp
      rdum=delr*max(0.0_wp,min(s1-s0,rdum-s0))
      ii=max(0,min(mpsi-1,int(rdum)))
      wp1=rdum-real(ii,wp)
      wpeon(m)=wp1

      ! inner flux surface
      im=ii
      tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
      tdum=(tdum-int(tdum))*delt(im)
      j00=max(0,min(mtheta(im)-1,int(tdum)))
      jteon0(m)=igrid(im)+j00  ! pay attention to j00=0 case
      wteon0(m)=tdum-real(j00,wp)

      ! outer flux surface
      im=ii+1
      tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
      tdum=(tdum-int(tdum))*delt(im)
      j01=max(0,min(mtheta(im)-1,int(tdum)))
      jteon1(m)=igrid(im)+j01  ! pay attention to j01=0 case
      wteon1(m)=tdum-real(j01,wp)
   enddo

end subroutine electron_locator



subroutine Cal_Current_AdiabaticTerm_at_PrevPosition
   ! This function calculates zeon0(9,m) = e*phi_NZ{t^n,x^(n-1)}/Te0{x^(n-1)}
   ! This value will be needed to calculate central differencing of convective adiabatic driving term
   use precision
   use field_array
   use global_parameters
   use particle_array
   use esi_interface
   use spl_profiles
   use profile_array

   implicit none

   integer m,i,ip,im
   real(wp) r, fip, dela, tempe
   integer ii,j01,j00
   real(wp) rdum,tdum,delr,delt(0:mpsi),delz, &
      pi2_inv,psitmp,thetatmp,zetatmp,psimax,psimin
   real(wp) prev_backTe(me), prev_wzeon(me), prev_wpeon(me), prev_jteon0(me), &
      prev_wteon0(me), prev_jteon1(me),prev_wteon1(me)
   integer ij
   real(wp) wz0,wz1,wp0,wp1,wt00,wt10,wt01,wt11
   real(wp) e1, Ea_00, T_inverse


   ! First, get prev_backTe at zeon0(1,m) position which is the position at time level (n-1)
   do m=1,me
      ! Get Te0 first.
      if(plasma==1) then    ! use experimental ni(r), T(r) and omega(r) profiles
         if(linterp_from_nodes==1) then
            dela=(a1-a0)/real(jr,wp)
            dela=1.0_wp/dela
            do i=1,me
               r=zeon0(1,i)-a0
               r=r*dela
               ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
               fip=r-real(ip,wp)

               prev_backTe(i)=fip*Te_node(ip+1)+(1.0_wp-fip)*Te_node(ip) ! Te at particles
            enddo
         else              ! use spline
            call ezspline_interp(spTe,me,zeon0(1,:),prev_backTe,ier)
         end if
      else                  ! use model profiles
         if(q2.gt.0.0) then  ! first analytic profile
            do i=1,me
               prev_backTe(i)=T_e0*Bt_e*exp(-dTeme*zeon0(1,i)**xt_e)
               prev_backTe(i)=T_e0*At_e+prev_backTe(i)
               prev_backTe(i)=Tc_e*(prev_backTe(i)/Tc_e)**b0     ! for new model Ti profile
            end do
         else      ! second analytic profile
            dela=(a1-a0)/real(jr,wp)
            dela=1.0_wp/dela
            do i=1,me
               r=zeon0(1,i)-a0
               r=r*dela
               ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
               fip=r-real(ip,wp)

               tempe=fip*Te_node(ip+1)+(1.0_wp-fip)*Te_node(ip)
               prev_backTe(i)=T_e0*exp(flow1*tempe)   ! Te at particles
            end do
         end if
      end if
   enddo !m=1,me


   ! Second, calculate prev_wpeon & etc at zeon0 position which is the position at time level (n-1)
   delr=1.0_wp/deltar
   delt=2.0_wp*pi/deltat
   delz=1.0_wp/deltaz
   pi2_inv=0.5_wp/pi
   psimin=s_grid(0)
   psimax=s_grid(mpsi)
   do m=1,me

      psitmp=zeon0(7,m)   ! =s
      thetatmp=zeon0(2,m)
      zetatmp=zeon0(3,m)

      prev_wzeon(m)=(zetatmp-zetamin)*delz

      rdum=psitmp
      rdum=delr*max(0.0_wp,min(s1-s0,rdum-s0))
      ii=max(0,min(mpsi-1,int(rdum)))
      wp1=rdum-real(ii,wp)
      prev_wpeon(m)=wp1

      ! inner flux surface
      im=ii
      tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
      tdum=(tdum-aint(tdum))*delt(im)
      j00=max(0,min(mtheta(im)-1,int(tdum)))
      prev_jteon0(m)=igrid(im)+j00  ! pay attention to j00=0 case
      prev_wteon0(m)=tdum-real(j00,wp)

      ! outer flux surface
      im=ii+1
      tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
      tdum=(tdum-aint(tdum))*delt(im)
      j01=max(0,min(mtheta(im)-1,int(tdum)))
      prev_jteon1(m)=igrid(im)+j01  ! pay attention to j01=0 case
      prev_wteon1(m)=tdum-real(j01,wp)
   enddo

   ! Last, calculate Adiabatic driving terms
   do m=1,me

      T_inverse=1.0_wp/prev_backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti

      wz1=prev_wzeon(m)
      wz0=1.0_wp-wz1

      ij=prev_jteon0(m)
      wp0=1.0_wp-prev_wpeon(m)
      wt00=1.0_wp-prev_wteon0(m)
      e1=   wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=ij+1
      wt10=1.0_wp-wt00
      e1=e1+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=prev_jteon1(m)
      wp1=1.0_wp-wp0
      wt01=1.0_wp-prev_wteon1(m)
      e1=e1+wp1*wt01*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=ij+1
      wt11=1.0_wp-wt01
      e1=e1+wp1*wt11*(wz0*phi(0,ij)+wz1*phi(1,ij)) ! --> phi at electrons

      if(ZF_SEPARATION==0) then
         !    calculate <phi> at electrons (no FLR effect -- x=R)
         r=zeon(7,m)-s0       ! in radial coordinator s
         r=r*delr
         ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
         fip=r-real(ip,wp)
         Ea_00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons
      end if

      zeon0(9,m)=(e1-real(1-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! Adiabatic term (q*dphi_NZ/Te)

   enddo !m=1,me


end subroutine Cal_Current_AdiabaticTerm_at_PrevPosition






subroutine Get_AdiabaticTerm_at_specific_Pos_Time(input_phi, input_phi00, input_zeon_arr, out_val_arr)
   ! This function calculates zeon0(9,m) = e*phi_NZ{t^n,x^(n-1)}/Te0{x^(n-1)}
   ! This value will be needed to calculate central differencing of convective adiabatic driving term
   use precision
   use field_array
   use global_parameters
   use particle_array
   use esi_interface
   use spl_profiles
   use profile_array

   implicit none

   real(wp), intent(in) :: input_phi(0:1,mgrid), input_phi00(0:mpsi), input_zeon_arr(zeondim, me)
   real(wp), intent(out):: out_val_arr(me)

   integer m, i, ip, im
   real(wp) r, fip, dela, tempe
   integer ii,j01,j00
   real(wp) rdum,tdum,delr,delt(0:mpsi),delz, &
      pi2_inv,psitmp,thetatmp,zetatmp,psimax,psimin
   real(wp) this_backTe(me), this_wzeon(me), this_wpeon(me), this_jteon0(me), &
      this_wteon0(me), this_jteon1(me),this_wteon1(me)
   integer ij
   real(wp) wz0,wz1,wp0,wp1,wt00,wt10,wt01,wt11
   real(wp) e1, Ea_00, T_inverse


   ! First, get this_backTe at zeon0(1,m) position which is the position at time level (n-1)
   do m=1,me
      ! Get Te0 first.
      if(plasma==1) then    ! use experimental ni(r), T(r) and omega(r) profiles
         if(linterp_from_nodes==1) then
            dela=(a1-a0)/real(jr,wp)
            dela=1.0_wp/dela
            do i=1,me
               r=zeon0(1,i)-a0
               r=r*dela
               ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
               fip=r-real(ip,wp)

               this_backTe(i)=fip*Te_node(ip+1)+(1.0_wp-fip)*Te_node(ip) ! Te at particles
            enddo
         else              ! use spline
            call ezspline_interp(spTe,me,zeon0(1,:),this_backTe,ier)
         end if
      else                  ! use model profiles
         if(q2.gt.0.0) then  ! first analytic profile
            do i=1,me
               this_backTe(i)=T_e0*Bt_e*exp(-dTeme*zeon0(1,i)**xt_e)
               this_backTe(i)=T_e0*At_e+this_backTe(i)
               this_backTe(i)=Tc_e*(this_backTe(i)/Tc_e)**b0     ! for new model Ti profile
            end do
         else      ! second analytic profile
            dela=(a1-a0)/real(jr,wp)
            dela=1.0_wp/dela
            do i=1,me
               r=zeon0(1,i)-a0
               r=r*dela
               ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
               fip=r-real(ip,wp)

               tempe=fip*Te_node(ip+1)+(1.0_wp-fip)*Te_node(ip)
               this_backTe(i)=T_e0*exp(flow1*tempe)   ! Te at particles
            end do
         end if
      end if
   enddo !m=1,me


   ! Second, calculate this_wpeon & etc at zeon0 position which is the position at time level (n-1)
   delr=1.0_wp/deltar
   delt=2.0_wp*pi/deltat
   delz=1.0_wp/deltaz
   pi2_inv=0.5_wp/pi
   psimin=s_grid(0)
   psimax=s_grid(mpsi)
   do m=1,me

      psitmp=zeon0(7,m)   ! =s
      thetatmp=zeon0(2,m)
      zetatmp=zeon0(3,m)

      this_wzeon(m)=(zetatmp-zetamin)*delz

      rdum=psitmp
      rdum=delr*max(0.0_wp,min(s1-s0,rdum-s0))
      ii=max(0,min(mpsi-1,int(rdum)))
      wp1=rdum-real(ii,wp)
      this_wpeon(m)=wp1

      ! inner flux surface
      im=ii
      tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
      tdum=(tdum-aint(tdum))*delt(im)
      j00=max(0,min(mtheta(im)-1,int(tdum)))
      this_jteon0(m)=igrid(im)+j00  ! pay attention to j00=0 case
      this_wteon0(m)=tdum-real(j00,wp)

      ! outer flux surface
      im=ii+1
      tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
      tdum=(tdum-aint(tdum))*delt(im)
      j01=max(0,min(mtheta(im)-1,int(tdum)))
      this_jteon1(m)=igrid(im)+j01  ! pay attention to j01=0 case
      this_wteon1(m)=tdum-real(j01,wp)
   enddo

   ! Last, calculate Adiabatic driving terms
   do m=1,me

      T_inverse=1.0_wp/this_backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti

      wz1=this_wzeon(m)
      wz0=1.0_wp-wz1

      ij=this_jteon0(m)
      wp0=1.0_wp-this_wpeon(m)
      wt00=1.0_wp-this_wteon0(m)
      e1=   wp0*wt00*(wz0*input_phi(0,ij)+wz1*input_phi(1,ij))

      ij=ij+1
      wt10=1.0_wp-wt00
      e1=e1+wp0*wt10*(wz0*input_phi(0,ij)+wz1*input_phi(1,ij))

      ij=this_jteon1(m)
      wp1=1.0_wp-wp0
      wt01=1.0_wp-this_wteon1(m)
      e1=e1+wp1*wt01*(wz0*input_phi(0,ij)+wz1*input_phi(1,ij))

      ij=ij+1
      wt11=1.0_wp-wt01
      e1=e1+wp1*wt11*(wz0*input_phi(0,ij)+wz1*input_phi(1,ij)) ! --> input_phi at electrons

      if(ZF_SEPARATION==0) then
         !    calculate <input_phi> at electrons (no FLR effect -- x=R)
         r=zeon(7,m)-s0       ! in radial coordinator s
         r=r*delr
         ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
         fip=r-real(ip,wp)
         Ea_00=fip*input_phi00(ip+1)+(1.0_wp-fip)*input_phi00(ip) ! <input_phi> at electrons
      end if

      out_val_arr(m)=(e1-real(1-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! Adiabatic term (q*dphi_NZ/Te)

   enddo !m=1,me


end subroutine Get_AdiabaticTerm_at_specific_Pos_Time



