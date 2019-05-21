subroutine pusheon(solver1, masssolver, irk_e)
  use precision
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface
  use solver_data_class

  implicit none

  type(solver_data) :: solver1, masssolver

!  integer,parameter :: nrbin=5       !wx
  integer i,im,ip,j,jt,k,kz,kk,m,ii,j11,j10,j01,j00,jtheta,ierror,ij,ipjt,irk_e
  real(wp) rdum,dtime,r,grada,energy,kappa,vdr,wdot,wdrift,&
       vdrenergy,delr,delt(0:mpsi),delz,pi2_inv,adum(0:mpsi),&
       xdum,ydum,xdot,ydot,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,&
       epara,psimax,psimin,pdum,wdrive,ainv,psitmp,&
       zetatmp,e1,e2,e3,wpi(3,me),ga(me),wga(3,me),gat(me),rdel, &
       dgr(me),da(me),dgq(me),dgf(me),mu, &
       Er_particle(me),dela,fip,weqflow,wgradB,T_inverse,Ea_00,wdphidt,wdphidt2, & 
       vpara, v2, avw, f_rho,  &
! following working arrays are defined to use OMP reduction:
       adum_00(0:mpsi),adum_tur(0:mpsi),adum_tot(0:mpsi)

  real(wp) vdq, vdgf, vd_para_mag, vd_para_a_cntv, vd_para_gq_cntv, vd_perp_a_cntv, vd_perp_gq_cntv
  real(wp) wdot_vpara, wdot_vperp, wdot_accel, now_w
  real(wp) tot_B_vec_cov(3), tot_B_vec_cntv(3) ! covariant and contravariant components of total B field (B0 + dB) incluidng the magnetic perturbation
  real(wp) grad_lnf0_a_cov, grad_lnf0_gq_cov ! covariant components of grad(ln(f0)) in a and gq axis. 
                                     !(gf component is 0 due to the toroidal axisymmetric assumption )

  real(wp) ga_amplitude  !ga_amplitude: ga(t) is a time-varying amplitude of Magnetic perturbation, only depends on the time in this version.
  real(wp) MP_phase ! MP_phase(x,t) : time-varying phase angle of given Magnetic Perturbation at specific position in radian unit.
  integer gcmotion_eon, gcmotion_eon_em, gcmotion_eon_em_MGY_vec_form

  real(wp) tt0,tt1,cput0,cput1,rtc

  real(wp), external :: D01AHF,dsda,SMOOTHSTEP
  integer minpts,ifail
  real(wp)  acc,error_8
  real(wp)  E_a,E_t,E_f,timef
  real(wp)  dminv,dmaxv
  real(wp)  wfix     ! =0 in the boundary buffer, and =1 otherwise

  real(wp) muB, instep,outstep,fdamp,fdamp0

  real(wp) wfield(3,mgrid,0:1)   ! weight-field W_0+v_//*W_1+v^2*W_2

  delr=1.0_wp/deltar
  delt=1.0_wp/deltat
  delz=1.0_wp/deltaz
  pi2_inv=0.5_wp/pi
  psimax=psi_grid(mpsi)    
  psimin=psi_grid(0)       

  f_rho=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)

  fdamp0=1.0_wp/tstep*5.0d-4
  if(nonlinear<0.5_wp) fdamp0=0.0_wp
  fdamp0=0.0_wp

  if(irk_e==1)then
! 1st step of Runge-Kutta method
     dtime=0.5_wp*tstep/real(nelectron,wp)

! use (x,y) coordinates near magnetic axis
!$omp parallel do private(m,xdum,ydum)
     do m=1,me
        zeon0(1:5,m)=zeon(1:5,m)
        zeon0(7,m)=zeon(7,m)
#ifdef _TREAT_MAGNETIC_AXIS
        if(sw(m).eq.1) then
           xdum=zeon(1,m)*cos(zeon(2,m))
           ydum=zeon(1,m)*sin(zeon(2,m))
           zeon0(1,m)=xdum
           zeon0(2,m)=ydum
        end if
#endif
     enddo

! 2nd step of Runge-Kutta method
  else
     dtime=tstep/real(nelectron,wp)
  endif

  call electron_locator

!  mu(1:me)=zeon(6,1:me)**2

! gather e_field at electron GC positions using linear interpolation.
! No FLR effect for electrons.
! The following line is an OpenMP directive for loop-level parallelism
! on shared memory machines (see http://www.openmp.org).
!$omp parallel do private(m,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11)
  do m=1,me
     wz1=wzeon(m)
     wz0=1.0_wp-wz1

     ij=jteon0(m)
     wp0=1.0_wp-wpeon(m)
     wt00=1.0_wp-wteon0(m)
     e1=wp0*wt00*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
     e2=wp0*wt00*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
     e3=wp0*wt00*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
     
     ij=ij+1
     wt10=1.0_wp-wt00
     e1=e1+wp0*wt10*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
     e2=e2+wp0*wt10*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
     e3=e3+wp0*wt10*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
     
     ij=jteon1(m)
     wp1=1.0_wp-wp0
     wt01=1.0_wp-wteon1(m)
     e1=e1+wp1*wt01*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
     e2=e2+wp1*wt01*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
     e3=e3+wp1*wt01*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
     
     ij=ij+1
     wt11=1.0_wp-wt01
     e1=e1+wp1*wt11*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
     e2=e2+wp1*wt11*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
     e3=e3+wp1*wt11*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

     wpi(1,m)=e1
     wpi(2,m)=e2
     wpi(3,m)=e3

     !MGYOO (10/30/2018; to use analytic profile of ga0(a))
    call Get_MP_alpha_amplitude(zeon(1,m), zeon(2,m), zeon(3,m), istep*tstep, ga(m), wga(1:3,m), gat(m))

  enddo

!    prepare equilibrium Er=dPhi/da at particles 

!  if(istep.lt.2800) then
if(istep.lt.26500) then
   timef=0.0_wp
else
   timef=1.0_wp-exp(-(real(istep-26500)/5.0d2)**4)
end if

   timef=1.0_wp

   Er401_bp=Er401_bp0*timef
   Er_grid=Er_grid0*timef

!  ------------------------------------------------------------



  dela=1.0_wp/Ydr
! if(istep.ge.4000) then
!$omp parallel do private(m,r,ip,fip)
  do m=1,me
     r=zeon(1,m)-a0      
     r=r*dela
     ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
     fip=r-real(ip,wp)
     Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
  enddo
! else
! Er_particle=0.0
! end if

!  tt0=rtc()
! call cpu_time(cput0)

  if(electromagnetic == 1)then
     !m=gcmotion_eon_em(dgr,da,dgq,dgf,zeon,wpi,ga,wga,gat,Er_particle,nonlinear,me,zeondim)
     !m=gcmotion_eon_em_MGY(dgr,da,dgq,dgf,zeon,wpi,ga,wga,gat,Er_particle,nonlinear,me,zeondim)
     m=gcmotion_eon_em_MGY_vec_form(dgr,da,dgq,dgf,zeon,wpi,ga,wga,gat,Er_particle,nonlinear,me,zeondim)
  else
  m=gcmotion_eon(dgr,da,dgq,dgf,zeon,wpi,Er_particle,nonlinear,me,zeondim)
  endif

! tt1=rtc()
! call cpu_time(cput1)

!  call n_T_u_group_e  ! prepare n, T, u_//, dlnn/da, dlnT/da, etc at electrons.
                       ! moved to collision_eon.f90 -- caution: in this case, BB must be obtained by call esigetb etc.!!!

  call Ti_at_electrons ! NOTE: array backT is borrowed to store Ti_at_electrons;
                       ! WARNING: must make sure there is no memory conflict!!!
                       ! backT is used in PUSHI and released afterward

!     OUTPUT: dgr=drho/dt;        da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt;    dgf=dzeta/dt;      
!             ghq(:)=da/dt|_EXB;  Ta(:)=dtheta/dt|_EXB;  Pa(:)=drho/dt|_EXB;
!             BB(:)=B;            Ba(:)=dB/da;           Bq(:)=dB/dtheta;       
!             gha(:)=1/|grad a|   F(:) =F
!             wpi(3,:)=dphi/dzeta (NO change); no change for: zeon(1,2,4,:), wpi(1,2,:) and mu;

!  if(nonlinear>100.0_wp) then            ! for nonlinear run
if(wimcompressibility==1) then            ! for implementation of weight imcompressibility 
  call weight_field_eon(wfield)      ! calculate weight-field uisng zeon0 for irk=1 and zeon at half-step for irk=2
                                     ! used BB obtained from gcmotion, and electron locations are obtain here in pushe 
                                     ! which are used to calculate wfield in phase space grids and gether wfield at electron (next)
end if




dminv=0.0d0
dmaxv=0.0d0



IF(ifelectron==2.and.ifdphi==2) then  ! use direct derivative along trajectory to calculate dphi/dt in delta-h electron case

! update GC position
! Another OpenMP parallel loop...
!$omp parallel do private (m,mu,energy,kappa,wdrift,wgradB,vdr,wdrive,wdot,xdum,ydum,T_inverse,&
!$omp& r,ip,Ea_00,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,wdphidt,wdphidt2,muB,fip,&
!$omp& instep,outstep,fdamp,vpara,v2,avw,E_f)
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

!!!! TEST
     !   dlnTdae(m)=0.0_wp
     !   dlnndae(m)=0.0_wp
!!!! TEST

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


        ! calculate d(delta_phi)/dt at particles along trajactory
        wz1=wzeon(m)
        wz0=1.0_wp-wz1

        ij=jteon0(m)
        wp0=1.0_wp-wpeon(m)
        wt00=1.0_wp-wteon0(m)
        e1=   wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))

        ij=ij+1
        wt10=1.0_wp-wt00
        e1=e1+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))

        ij=jteon1(m)
        wp1=1.0_wp-wp0
        wt01=1.0_wp-wteon1(m)
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

        E_f=(e1-real(1-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0
        if(irk_e==1) zeon0(8,m)=E_f    ! e*delta_phi/Te at electron beginning position



     if(zeon0(1,m).le.lmax.or.zeon0(1,m).ge.rmin) then
        instep =smoothstep(lmin,lmax,zeon(1,m))
        outstep=smoothstep(rmin,rmax,zeon(1,m))
        !  fip=smoothstep(0.20_wp,0.25_wp,zeon(1,m))*(1.0_wp-smoothstep(0.75_wp,0.80_wp,zeon(1,m)))
        fip=instep*(1.0_wp-outstep)
        !  wdrive=wdrive*fip

        fdamp=(1.0_wp-fip)*fdamp0
     else
        fip=1.0_wp
        fdamp=0.0_wp
     end if

     if(ifneoclassical==1) then   
        wdrift=(dgr(m)*BB(m)+zeon(4,m)*wgradB)*zeon(4,m)*BB(m)*for_E*T_inverse   ! see note
     else
        wdrift=(Pa(m)*BB(m)+zeon(4,m)*wgradB)*zeon(4,m)*BB(m)*for_E*T_inverse   ! see note
     end if

     wgradB=wgradB*mu*T_inverse   ! see note: additional term due to v_ExB*grad(f_M)

     wdot=(wdrive+wdrift+wgradB)

     !MGYOO changed this part. (10/31/2018)
     wdot=(1.0_wp-E_f)*wdot



     !!!!!  if(nonlinear>0.5_wp) then     ! for nonlinear run
     !!!!!  if(nonlinear>100.0_wp) then          ! for nonlinear run
     if(wimcompressibility==1) then     ! for weight imcompressibility
     ! calculate average weight <w> at particles using wieght-field wfield
        wz1=wzeon(m)
        wz0=1.0_wp-wz1

        ij=jteon0(m)
        wp0=1.0_wp-wpeon(m)

        wt00=1.0_wp-wteon0(m)
        e1=wp0*wt00*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
        e2=wp0*wt00*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
        e3=wp0*wt00*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

        ij=ij+1
        wt10=1.0_wp-wt00
        e1=e1+wp0*wt10*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
        e2=e2+wp0*wt10*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
        e3=e3+wp0*wt10*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

        ij=jteon1(m)
        wp1=1.0_wp-wp0
        wt01=1.0_wp-wteon1(m)
        e1=e1+wp1*wt01*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
        e2=e2+wp1*wt01*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
        e3=e3+wp1*wt01*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

        ij=ij+1
        wt11=1.0_wp-wt01
        e1=e1+wp1*wt11*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
        e2=e2+wp1*wt11*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
        e3=e3+wp1*wt11*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

        avw=e1+e2*vpara+e3*v2


    !     if(mype==100.and.m==me/2) print*, 'AVW=', avw, zeon(5,m), mype
    !     if(mype==0.and.m==me/4) print*, 'AVW=', avw, zeon(5,m), mype
    !     if(mype==200.and.m==me/8) print*, 'AVW=', avw, zeon(5,m), mype

    end if



    !==================================================
    !MGYOO(01/14/2019)
    !Note that f0 is toroidally symmetric so that grad_lnf0_gf=0 and we do not need vd_para_gf or vd_perp_gf terms to calculate weight evolution
    if(irk_e==2) then
      if(ifneoclassical==1) then
        !Caution!! Following routine is only available for neoclassical case which use full velocity.
        call esiget_Intp_B_vec_cov_and_cntv(zeon(1,m), zeon(2,m), zeon(3,m), istep*tstep, tot_B_vec_cov, tot_B_vec_cntv)
        vd_para_mag = (vdr*tot_B_vec_cov(1) + vdq*tot_B_vec_cov(2) + vdgf*tot_B_vec_cov(3))/BB(m)
        vd_para_a_cntv = vd_para_mag*(tot_B_vec_cntv(1)/BB(m))
        vd_para_gq_cntv = vd_para_mag*(tot_B_vec_cntv(2)/BB(m))
      else ! perturbed components only (ignoring neoclassical parts)
        vd_para_a_cntv = vpara_dve_a(m)
        vd_para_gq_cntv = vpara_dve_gq(m)      
      endif
      vd_perp_a_cntv = vdr -vd_para_a_cntv
      vd_perp_gq_cntv = vdq -vd_para_gq_cntv

      !grad_lnf0_a_cov = -kappa - T_inverse*mu*Ba(m) 
      !grad_lnf0_gq_cov = -T_inverse*mu*Bq(m)

      ! wdot_vpara = -vd_para_a_cntv*grad_lnf0_a_cov -vd_para_gq_cntv*grad_lnf0_gq_cov
      ! wdot_vperp = -vd_perp_a_cntv*grad_lnf0_a_cov -vd_perp_gq_cntv*grad_lnf0_gq_cov
      ! wdot_accel  = wdrift
      
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

      if(wimcompressibility==1) then
        now_w = avw
      else
        now_w = zeon(5,m)
      endif

      we_vpara(m) = we_vpara(m) + dtime*( (1.0-nonlinear*now_w)*wdot_vpara*fip)
      we_vperp(m) = we_vperp(m) + dtime*( (1.0-nonlinear*now_w)*wdot_vperp*fip)
      we_accel(m) = we_accel(m) + dtime*( (1.0-nonlinear*now_w)*wdot_accel*fip)
      we_fdamp(m) = we_fdamp(m) + dtime*(-fdamp*zeon(5,m))
    endif
    !=================================================

     zeon(1,m) = zeon0(1,m)+dtime*da(m)
     zeon(2,m) = zeon0(2,m)+dtime*dgq(m)
     zeon(3,m) = zeon0(3,m)+dtime*dgf(m)
     zeon(4,m) = zeon0(4,m)+dtime*dgr(m)

     if(wimcompressibility==1) then
        zeon(8,m) = avw           ! save average weight avw at electrons
     else
        zeon(8,m) = zeon(5,m)     ! save w 
     end if

!     zeon(5,m) = zeon0(5,m)+dtime*wdot

!     zeon(5,m) = 1.0_wp-(1.0_wp-zeon0(5,m))*exp(-dtime*wdot)

!    two cases below are unified by the replacement using zeon(8,m) in the formula
       !!!!zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*zeon(5,m))*wdot-fdamp*zeon(5,m))*dtime  ! weight advanced w/o including dphi/dt term yet
       !!!!zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*avw)*wdot-fdamp*zeon(5,m))*dtime

     !  if(zeon0(1,m).gt.psi_grid(mpsi-nbound-1))  then 
     !     wfix=1.0_wp
     !  else
     !     wfix=1.0_wp
     !  endif

     !  zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*zeon(8,m))*wdot-fdamp*zeon(5,m))*dtime*wfix  ! weight advanced w/o including dphi/dt term yet
     zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*zeon(8,m))*wdot*fip-fdamp*zeon(5,m))*dtime       ! weight advanced w/o including dphi/dt term yet
                                                                                                ! unified expression April, 2013
    


!     dminv=min(dminv,wdot*tstep)
!     dmaxv=max(dmaxv,wdot*tstep)



!     zeon(7,m) = zeon0(7,m)+dtime*sqrt(Tc/backT(m))*da(m)
!     zeon(7,m) = zeon0(7,m)+dtime*(Tc/backT(m))**0.4_wp*da(m)
     zeon(7,m) = zeon0(7,m)+dtime*(Tc/backT(m))**dsda_power*da(m)

#ifdef _TREAT_MAGNETIC_AXIS
     if(sw(m).eq.1) then    ! convert to (a, theta) coordinates
        xdum=sqrt(zeon(1,m)**2+zeon(2,m)**2)
        ydum=sign(1.0_wp,zeon(2,m))*acos(zeon(1,m)/xdum)
        zeon(1,m)=xdum
        zeon(2,m)=ydum
     end if
#endif

! theta and zeta normalize to [0,2*pi), modulo is slower than hand coded procedure
!     zeon(2,m)=zeon(2,m)*pi2_inv+100.0 !period of 1
!     zeon(2,m)=2.0*pi*(zeon(2,m)-aint(zeon(2,m))) ![0,2*pi)
!     zeon(3,m)=zeon(3,m)*pi2_inv+100.0
!     zeon(3,m)=2.0*pi*(zeon(3,m)-aint(zeon(3,m)))

!     zeon(2,m)=zeon(2,m)*pi2_inv !period of 1
!     zeon(2,m)=10.0+zeon(2,m)-dint(zeon(2,m)) ![0,2*pi)
!     zeon(2,m)=twopi*(zeon(2,m)-dint(zeon(2,m)))
!     zeon(3,m)=zeon(3,m)*pi2_inv
!     zeon(3,m)=10.0+zeon(3,m)-dint(zeon(3,m))
!     zeon(3,m)=twopi*(zeon(3,m)-dint(zeon(3,m)))

! use f90 modulo to normalize theta and zeta to [0,2*pi)
! If P > 0, then the function modulo(A,P) determines the positive difference between A and then
! next lowest multiple of P. If P < 0, then the function modulo(A,P) determines the negative
! difference between A and then next highest multiple of P.
     zeon(2,m)=modulo(zeon(2,m),twopi)
     zeon(3,m)=modulo(zeon(3,m),twopi)


! store GC information for flux measurements
     wpi(1,m)=ghq(m)    ! at present, wpi(1,*)=da/dt|_EXB 
     wpi(2,m)=energy    ! normalized electron energy with u//=0 in Jan.09
     wpi(3,m)=muB       ! normalizd electron perpendicular energy

     if(irk_e==2) then
      UnNorm_theta_e(m) = UnNorm_theta_e(m) +dtime*dgq(m)
      UnNorm_zeta_e(m)  = UnNorm_zeta_e(m)  +dtime*dgf(m)
     endif
  enddo




!  if(mype==0) print*, 'ELECTRONS:  Max Damping=', dmaxv, '   Min Damping=', dminv
!  if(mype==0) print*, 'ELECTRONS:  ', 'Max Damping=', dmaxv, '   Min Damping=', dminv




  !   fluxes diagnose
  if(idiag_e==0.and.irk_e==1) then  ! istep=n*ndiag, irk=1, i=1, and irk_e=1
  !$omp parallel do private(m)
     do m=1,me
        wpi(1,m)=wpi(1,m)*gha(m)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
     end do

     call get_fluxes_e(wpi)

     !   im=mod(istep,100)
     im=mod(istep,ndiag)
     if(im==0) then
        call get_Te_fluctuation(wpi(2,:))
     end if
  endif

   call shifte

 call electron_locator
 call n_T_u_group_e             ! get Te at electron position
 zeon00(8,1:me)=zeon(5,1:me)    ! use for temporary storage -- check again April, 2013
!  E_a=0.0_wp
!  E_f=0.0_wp
!  E_t=0.0_wp

!$omp parallel do private(m,T_inverse,r,ip,fip,Ea_00,instep,outstep,fdamp, &
!$omp& e1,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,wdphidt)
  do m=1,me
     T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti

      ! calculate d(delta_phi)/dt at particles
      wz1=wzeon(m)
      wz0=1.0_wp-wz1

      ij=jteon0(m)
      wp0=1.0_wp-wpeon(m)
      wt00=1.0_wp-wteon0(m)
      e1=   wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=ij+1
      wt10=1.0_wp-wt00
      e1=e1+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=jteon1(m)
      wp1=1.0_wp-wp0
      wt01=1.0_wp-wteon1(m)
      e1=e1+wp1*wt01*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=ij+1
      wt11=1.0_wp-wt01
      e1=e1+wp1*wt11*(wz0*phi(0,ij)+wz1*phi(1,ij)) ! --> phi at electron new positions

      if(ZF_SEPARATION==0) then
        !    calculate <phi> at electrons (no FLR effect -- x=R)
         r=zeon(7,m)-s0       ! in radial coordinator s
         r=r*delr
         ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
         fip=r-real(ip,wp)
         Ea_00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons
      end if

      wdphidt=(e1-real(1-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0   ! e*delta_phi/Te at electron new positions, but with unchanged delta_phi
      wdphidt=(wdphidt-zeon0(8,m))/dtime  ! d(e*delta_phi/Te)/dt purely due to change of electron position

      !  if(zeon0(1,m).gt.psi_grid(mpsi-nbound-1))  then
      !     wfix=1.0_wp
      !  else
      !     wfix=1.0_wp
      !  endif

      if(zeon0(1,m).le.lmax.or.zeon0(1,m).ge.rmin) then
         instep =smoothstep(lmin,lmax,zeon(1,m))
         outstep=smoothstep(rmin,rmax,zeon(1,m))
         fip=instep*(1.0_wp-outstep)
         !  wdrive=wdrive*fip
         fdamp=(1.0_wp-fip)*fdamp0
      else
         fip=1.0_wp
         fdamp=0.0_wp
      end if

      wdphidt=wdphidt

        !  zeon(5,m) = zeon00(8,m)+(zeon0(6,m)-nonlinear*zeon(8,m))*wdphidt*dtime*wfix  ! d(e*delta_phi/Te)/dt term added to weight,  but with delta_phi not advanced
      zeon(5,m) = zeon00(8,m)+(zeon0(6,m)-nonlinear*zeon(8,m))*wdphidt*fip*dtime  ! d(e*delta_phi/Te)/dt term added to weight,  but with delta_phi not advanced

    


      !       E_a=E_a+zeon(5,m)

  end do

  call charge_eon
  call poisson( solver1, masssolver )
  call smooth_PHI

!$omp parallel do private(m,T_inverse,r,ip,fip,Ea_00,instep,outstep,fdamp, &
!$omp& e1,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,wdphidt)
  do m=1,me
      T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti

      ! calculate d(delta_phi)/dt at particles
      wz1=wzeon(m)
      wz0=1.0_wp-wz1

      ij=jteon0(m)
      wp0=1.0_wp-wpeon(m)
      wt00=1.0_wp-wteon0(m)
      e1=   wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=ij+1
      wt10=1.0_wp-wt00
      e1=e1+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=jteon1(m)
      wp1=1.0_wp-wp0
      wt01=1.0_wp-wteon1(m)
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

      wdphidt=(e1-real(1-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0  ! e*delta_phi/Te at electron new positions, but with advanced delta_phi
      wdphidt=(wdphidt-zeon0(8,m))/dtime ! d(e*delta_phi/Te)/dt with correction of time derivative of delta_phi

      !  if(zeon0(1,m).gt.psi_grid(mpsi-nbound-1))  then
      !     wfix=1.0_wp
      !  else
      !     wfix=1.0_wp
      !  endif

      if(zeon0(1,m).le.lmax.or.zeon0(1,m).ge.rmin) then
         instep =smoothstep(lmin,lmax,zeon(1,m))
         outstep=smoothstep(rmin,rmax,zeon(1,m))
         fip=instep*(1.0_wp-outstep)
         !  wdrive=wdrive*fip
         fdamp=(1.0_wp-fip)*fdamp0
      else
         fip=1.0_wp
         fdamp=0.0_wp
      end if

      wdphidt=wdphidt

      !  zeon(5,m) = zeon00(8,m)+(zeon0(6,m)-nonlinear*zeon(8,m))*wdphidt*dtime*wfix ! d(e*delta_phi/Te)/dt term added to weight, with delta_phi advanced
      zeon(5,m) = zeon00(8,m)+(zeon0(6,m)-nonlinear*zeon(8,m))*wdphidt*fip*dtime ! d(e*delta_phi/Te)/dt term added to weight, with delta_phi advanced

      !==================================================
      !MGYOO(01/14/2019)
      ! Cal weight evolution due to the adiabatic response (d/dt(q*dphi/Te))
      if(irk_e==2) then
        we_adiab(m) = we_adiab(m) + (zeon0(6,m)-nonlinear*zeon(8,m))*wdphidt*fip*dtime    
      endif
      !=================================================         



      ! if (mype ==0 .AND. MOD(m,(me/10))==0) then
      !   write(stdout,*) 'my_rank= ',mype, ' zeon(5,m=',m,')= ',zeon(5,m), ' zeon0(5,m=',m,')= ',zeon0(5,m), ' avw= ',avw , ' zeon(8,m=',m,')= ',zeon(8,m), ' zeon00(8,m=',m,')= ',zeon00(8,m)
      !   write(stdout,*) 'zeon(6,m=',m,')= ',zeon(6,m), ' zeon0(6,m=',m,')= ',zeon0(6,m) , ' wdphidt= ', wdphidt
      ! endif

!       E_f=E_f+zeon(5,m)
  end do






!   if(mype==0) print*, 'BBBB=', E_a, '  SSSS=', E_f, '  ZZZZ=', E_t


 if(irk_e==2)then
! out of boundary particle
!$omp parallel do private(m)
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


ELSE    !  for df electrons, adiabatic electrons, and delta-h electrons with dphi/dt equation




! update GC position
! Another OpenMP parallel loop...
!$omp parallel do private(m,mu,energy,kappa,wdrift,wgradB,vdr,wdrive,wdot,xdum,ydum,T_inverse,r,ip,Ea_00, &
!$omp& e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,wdphidt,wdphidt2,muB,fip,instep,&
!$omp& outstep,fdamp,vpara,v2,avw)
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
        vdr=da(m)
        vdq=dgq(m)
        vdgf=dgf(m)
     else
        !  ExB drift in radial direction for w-dot and flux diagnostics
        vdr=ghq(m)             !  da/dt|_EXB
        vdq=Ta(m)
        vdgf=dve_gf(m)
     end if

     wdrive=vdr*kappa
     wgradB=vdr*Ba(m)+vdq*Bq(m)


     if(ifelectron==2) then 

! calculate d(delta_phi)/dt at particles
        wz1=wzeon(m)
        wz0=1.0_wp-wz1

        ij=jteon0(m)
        wp0=1.0_wp-wpeon(m)
        wt00=1.0_wp-wteon0(m)
        e1=   wp0*wt00*(wz0*dphidt(0,ij)+wz1*dphidt(1,ij))

        ij=ij+1
        wt10=1.0_wp-wt00
        e1=e1+wp0*wt10*(wz0*dphidt(0,ij)+wz1*dphidt(1,ij))
  
        ij=jteon1(m)
        wp1=1.0_wp-wp0
        wt01=1.0_wp-wteon1(m)
        e1=e1+wp1*wt01*(wz0*dphidt(0,ij)+wz1*dphidt(1,ij))

        ij=ij+1
        wt11=1.0_wp-wt01
        e1=e1+wp1*wt11*(wz0*dphidt(0,ij)+wz1*dphidt(1,ij)) ! --> dphidt at electrons
  
        wdphidt=e1*csc0*T_inverse    ! when e1 is dphi/dt at electrons -- using FEM solver with order>0
        !  wdphidt=e1*csc0           ! when dphi/dt equation is modified; dphidt=[d(delta_phi)/dt]/Te using old integral solver with order=0

!    calculate d<phi>/da at electrons (no FLR effect -- x=R)
        r=zeon(7,m)-s0       ! in radial coordinator s
        r=r*delr
        ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
        fip=r-real(ip,wp)
        Ea_00=fip*phip00(ip+1)+(1.0_wp-fip)*phip00(ip) ! d<phi>/da at electrons

        wdphidt2=da(m)*(wpi(1,m)-Ea_00)+dgq(m)*wpi(2,m)+dgf(m)*wpi(3,m)
        wdphidt2=wdphidt2*T_inverse*csc0
        wdphidt=wdphidt+wdphidt2          ! d/dt(delta-phi) term in RHS of dhe/dt equation

     else                  ! delta-f electrons
        wdphidt=0.0_wp
     end if

     if(zeon0(1,m).le.lmax.or.zeon0(1,m).ge.rmin) then
        instep =smoothstep(lmin,lmax,zeon(1,m))
        outstep=smoothstep(rmin,rmax,zeon(1,m))
        !  fip=smoothstep(0.20_wp,0.25_wp,zeon(1,m))*(1.0_wp-smoothstep(0.75_wp,0.80_wp,zeon(1,m)))
        fip=instep*(1.0_wp-outstep)
        !  wdrive=wdrive*fip
        !  wdphidt=wdphidt*fip

        fdamp=(1.0_wp-fip)*fdamp0
     else
        fip=1.0_wp
        fdamp=0.0_wp
     end if

     if(ifneoclassical==1) then
        wdrift=(dgr(m)*BB(m)+zeon(4,m)*wgradB)*zeon(4,m)*BB(m)*for_E*T_inverse   ! see note
     else
        wdrift=(Pa(m)*BB(m)+zeon(4,m)*wgradB)*zeon(4,m)*BB(m)*for_E*T_inverse   ! see note
     endif

     wgradB=wgradB*mu*T_inverse   ! see note: additional term due to v_ExB*grad(f_M)

     wdot=(wdrive+wdrift+wgradB+wdphidt)

     wdot=wdot


!!!!!  if(nonlinear>0.5_wp) then     ! for nonlinear run
!!!!!  if(nonlinear>100.0_wp) then        ! for nonlinear run
if(wimcompressibility==1) then     ! for weight imcompressibility
! calculate average weight <w> at particles using wieght-field wfield
     wz1=wzeon(m)
     wz0=1.0_wp-wz1

     ij=jteon0(m)
     wp0=1.0_wp-wpeon(m)

     wt00=1.0_wp-wteon0(m)
     e1=wp0*wt00*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=wp0*wt00*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=wp0*wt00*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt10=1.0_wp-wt00
     e1=e1+wp0*wt10*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp0*wt10*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp0*wt10*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=jteon1(m)
     wp1=1.0_wp-wp0
     wt01=1.0_wp-wteon1(m)
     e1=e1+wp1*wt01*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt01*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt01*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt11=1.0_wp-wt01
     e1=e1+wp1*wt11*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt11*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt11*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     avw=e1+e2*vpara+e3*v2

!     if(mype==100.and.m==me/2) print*, 'AVW=', avw, zeon(5,m), mype
!     if(mype==0.and.m==me/4) print*, 'AVW=', avw, zeon(5,m), mype
!     if(mype==200.and.m==me/8) print*, 'AVW=', avw, zeon(5,m), mype

end if



    !==================================================
    !MGYOO(01/14/2019)
    !Note that f0 is toroidally symmetric so that grad_lnf0_gf=0 and we do not need vd_para_gf or vd_perp_gf terms to calculate weight evolution
    if(irk_e==2) then
      if(ifneoclassical==1) then
        !Caution!! Following routine is only available for neoclassical case which use full velocity.
        call esiget_Intp_B_vec_cov_and_cntv(zeon(1,m), zeon(2,m), zeon(3,m), istep*tstep, tot_B_vec_cov, tot_B_vec_cntv)
        vd_para_mag = (vdr*tot_B_vec_cov(1) + vdq*tot_B_vec_cov(2) + vdgf*tot_B_vec_cov(3))/BB(m)
        vd_para_a_cntv = vd_para_mag*(tot_B_vec_cntv(1)/BB(m))
        vd_para_gq_cntv = vd_para_mag*(tot_B_vec_cntv(2)/BB(m))
      else ! perturbed components only (ignoring neoclassical parts)
        vd_para_a_cntv = vpara_dve_a(m)
        vd_para_gq_cntv = vpara_dve_gq(m)      
      endif
      vd_perp_a_cntv = vdr -vd_para_a_cntv
      vd_perp_gq_cntv = vdq -vd_para_gq_cntv

      ! grad_lnf0_a_cov = -kappa - T_inverse*mu*Ba(m) 
      ! grad_lnf0_gq_cov = -T_inverse*mu*Bq(m)
     
      ! wdot_vpara = -vd_para_a_cntv*grad_lnf0_a_cov -vd_para_gq_cntv*grad_lnf0_gq_cov
      ! wdot_vperp = -vd_perp_a_cntv*grad_lnf0_a_cov -vd_perp_gq_cntv*grad_lnf0_gq_cov
      ! wdot_accel  = wdrift

      !=================================================================
      !MGYOO (01/30/2019)
      ! grad_lnf0 includes (wdrive) term
      ! wdot_accel includes (wdrift + wgradB + wdphidt) terms
      grad_lnf0_a_cov = -kappa 
      grad_lnf0_gq_cov = 0.0

      wdot_vpara = -vd_para_a_cntv*grad_lnf0_a_cov -vd_para_gq_cntv*grad_lnf0_gq_cov
      wdot_vperp = -vd_perp_a_cntv*grad_lnf0_a_cov -vd_perp_gq_cntv*grad_lnf0_gq_cov
      wdot_accel  = wdrift + wgradB + wdphidt
      !=================================================================      

      wdot_vpara = (1.0_wp-E_f)*wdot_vpara
      wdot_vperp = (1.0_wp-E_f)*wdot_vperp
      wdot_accel = (1.0_wp-E_f)*wdot_accel

      if(wimcompressibility==1) then
        now_w = avw
      else
        now_w = zeon(5,m)
      endif

      we_vpara(m) = we_vpara(m) + dtime*( (1.0-nonlinear*now_w)*wdot_vpara*fip)
      we_vperp(m) = we_vperp(m) + dtime*( (1.0-nonlinear*now_w)*wdot_vperp*fip)
      we_accel(m) = we_accel(m) + dtime*( (1.0-nonlinear*now_w)*wdot_accel*fip)
      we_fdamp(m) = we_fdamp(m) + dtime*(-fdamp*zeon(5,m))
    endif
    !=================================================




     zeon(1,m) = zeon0(1,m)+dtime*da(m)
     zeon(2,m) = zeon0(2,m)+dtime*dgq(m)
     zeon(3,m) = zeon0(3,m)+dtime*dgf(m)
     zeon(4,m) = zeon0(4,m)+dtime*dgr(m)

!     zeon(5,m) = zeon0(5,m)+dtime*wdot

!     zeon(5,m) = 1.0_wp-(1.0_wp-zeon0(5,m))*exp(-dtime*wdot)

      if(wimcompressibility==1) then     ! for weight imcompressibility
         zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*avw)*wdot*fip-fdamp*zeon(5,m))*dtime
      else
         zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*zeon(5,m))*wdot*fip-fdamp*zeon(5,m))*dtime
      end if



!     zeon(7,m) = zeon0(7,m)+dtime*sqrt(Tc/backT(m))*da(m)
!     zeon(7,m) = zeon0(7,m)+dtime*(Tc/backT(m))**0.4_wp*da(m)
     zeon(7,m) = zeon0(7,m)+dtime*(Tc/backT(m))**dsda_power*da(m)

#ifdef _TREAT_MAGNETIC_AXIS
     if(sw(m).eq.1) then    ! convert to (a, theta) coordinates
        xdum=sqrt(zeon(1,m)**2+zeon(2,m)**2)
        ydum=sign(1.0_wp,zeon(2,m))*acos(zeon(1,m)/xdum)
        zeon(1,m)=xdum
        zeon(2,m)=ydum
     end if
#endif
     
! theta and zeta normalize to [0,2*pi), modulo is slower than hand coded procedure
!     zeon(2,m)=zeon(2,m)*pi2_inv+100.0 !period of 1
!     zeon(2,m)=2.0*pi*(zeon(2,m)-aint(zeon(2,m))) ![0,2*pi)
!     zeon(3,m)=zeon(3,m)*pi2_inv+100.0
!     zeon(3,m)=2.0*pi*(zeon(3,m)-aint(zeon(3,m)))

!     zeon(2,m)=zeon(2,m)*pi2_inv !period of 1
!     zeon(2,m)=10.0+zeon(2,m)-dint(zeon(2,m)) ![0,2*pi)
!     zeon(2,m)=twopi*(zeon(2,m)-dint(zeon(2,m)))
!     zeon(3,m)=zeon(3,m)*pi2_inv
!     zeon(3,m)=10.0+zeon(3,m)-dint(zeon(3,m))
!     zeon(3,m)=twopi*(zeon(3,m)-dint(zeon(3,m)))

! use f90 modulo to normalize theta and zeta to [0,2*pi)
! If P > 0, then the function modulo(A,P) determines the positive difference between A and then
! next lowest multiple of P. If P < 0, then the function modulo(A,P) determines the negative
! difference between A and then next highest multiple of P.
     zeon(2,m)=modulo(zeon(2,m),twopi)
     zeon(3,m)=modulo(zeon(3,m),twopi)


! store GC information for flux measurements
     wpi(1,m)=ghq(m)    ! at present, wpi(1,*)=da/dt|_EXB 
     wpi(2,m)=energy    ! normalized electron energy with u//=0 in Jan.09
     wpi(3,m)=muB       ! normalizd electron perpendicular energy

     if(irk_e==2) then
      UnNorm_theta_e(m) = UnNorm_theta_e(m) +dtime*dgq(m)
      UnNorm_zeta_e(m)  = UnNorm_zeta_e(m)  +dtime*dgf(m)
     endif     
  enddo

 if(irk_e==2)then
! out of boundary particle
!$omp parallel do private(m)
    do m=1,me
       if(zeon(1,m) > psimax)then
          zeon(1,m)=zeon0(1,m)
          zeon(2,m)=2.0_wp*pi-zeon0(2,m)
          zeon(3,m)=zeon0(3,m)
          zeon(4,m)=zeon0(4,m)
          zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
!          zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
          zeon(7,m)=zeon0(7,m)
          
       elseif(zeon(1,m) < psimin)then
          zeon(1,m)=zeon0(1,m)
          zeon(2,m)=2.0_wp*pi-zeon0(2,m)
          zeon(3,m)=zeon0(3,m)
          zeon(4,m)=zeon0(4,m) 
          zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
!          zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
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

 if(idiag_e==0.and.irk_e==1) then  ! istep=n*ndiag, irk=1, i=1, and irk_e=1
! fluxes diagnose
!$omp parallel do private(m)
    do m=1,me
       wpi(1,m)=wpi(1,m)*gha(m)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
    end do

    call get_fluxes_e(wpi)

!   im=mod(istep,100)
   im=mod(istep,ndiag)
   if(im==0) then
      call get_Te_fluctuation(wpi(2,:))
   end if

 endif

END IF

end subroutine pusheon


subroutine pusheon_RK4(solver1, masssolver, irk_e)
  use precision
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface
  use solver_data_class

  implicit none

  type(solver_data) :: solver1, masssolver
  
  integer, intent(in) :: irk_e
  integer i,im,ip,j,jt,k,kz,kk,m,ii,j11,j10,j01,j00,jtheta,ierror,ij,ipjt
  real(wp) rdum,r,grada,energy,kappa,vdr,wdot,wdrift,&
       vdrenergy,delr,delt(0:mpsi),delz,pi2_inv,adum(0:mpsi),&
       xdum,ydum,xdot,ydot,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,&
       epara,psimax,psimin,pdum,wdrive,ainv,psitmp,&
       zetatmp,e1,e2,e3,wpi(3,me),ga(me),wga(3,me),gat(me),rdel, &
       dgr(me),da(me),dgq(me),dgf(me),mu, &
       Er_particle(me),dela,fip,weqflow,wgradB,T_inverse,Ea_00,wdphidt,wdphidt2, & 
       vpara, v2, avw, f_rho,  &
       ! following working arrays are defined to use OMP reduction:
       adum_00(0:mpsi),adum_tur(0:mpsi),adum_tot(0:mpsi)
  real(wp) tdum,thetatmp
  real(Wp) tmp_old_AD, tmp_new_AD, tmp_old_dADdt, tmp_new_dADdt
  real(wp) vdq, vdgf, vd_para_mag, vd_para_a_cntv, vd_para_gq_cntv, vd_perp_a_cntv, vd_perp_gq_cntv
  real(wp) wdot_vpara, wdot_vperp, wdot_accel, now_w
  real(wp) tot_B_vec_cov(3), tot_B_vec_cntv(3) ! covariant and contravariant components of total B field (B0 + dB) incluidng the magnetic perturbation
  real(wp) grad_lnf0_a_cov, grad_lnf0_gq_cov ! covariant components of grad(ln(f0)) in a and gq axis. 
                                     !(gf component is 0 due to the toroidal axisymmetric assumption )

  real(wp) ga_amplitude  !ga_amplitude: ga(t) is a time-varying amplitude of Magnetic perturbation, only depends on the time in this version.
  real(wp) MP_phase ! MP_phase(x,t) : time-varying phase angle of given Magnetic Perturbation at specific position in radian unit.
  integer gcmotion_eon, gcmotion_eon_em, gcmotion_eon_em_MGY_vec_form

  real(wp) tt0,tt1,cput0,cput1,rtc

  real(wp), external :: D01AHF,dsda,SMOOTHSTEP
  integer minpts,ifail
  real(wp)  acc,error_8
  real(wp)  E_a,E_t,E_f,timef
  real(wp)  dminv,dmaxv
  real(wp)  wfix     ! =0 in the boundary buffer, and =1 otherwise

  real(wp) muB, instep,outstep,fdamp,fdamp0

  real(wp) wfield(3,mgrid,0:1)   ! weight-field W_0+v_//*W_1+v^2*W_2
  real(wp) RK_cum_coeff, RK_adv_coeff ! Coefficients to cumlate (to get future final value) or advance(to predict next RK step) RK4 terms

  delr=1.0_wp/deltar
  delt=1.0_wp/deltat
  delz=1.0_wp/deltaz
  pi2_inv=0.5_wp/pi
  psimax=psi_grid(mpsi)    
  psimin=psi_grid(0)       

  f_rho=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)

  fdamp0=1.0_wp/tstep*5.0d-4
  if(nonlinear<0.5_wp) fdamp0=0.0_wp
  fdamp0=0.0_wp

  if(irk_e==1)then ! 1st step of Runge-Kutta method
    RK_cum_zeon  =0.0
    RK_cum_coeff =(1.0/6.0)
    RK_adv_coeff =0.5
    do m=1,me
      ! Note that zeon0(6,m) does not mean "mu", this is reserved for the other quantity
      ! Please note that zeon0(5,m) and zeon0(8,m) is different quantity.
      ! zeon0(8,m) is the adiabatic drive term in the weight evolution @ t={n}, x={n}
      ! zeon(8,m) is a kind of temporay pouch to bring wdot term during shifte
      zeon0(1:5,m)=zeon(1:5,m)
      zeon0(7,m)=zeon(7,m)
    enddo
  elseif(irk_e==2) then
    RK_cum_coeff =(1.0/3.0)
    RK_adv_coeff =0.5 
  elseif(irk_e==3) then
    RK_cum_coeff =(1.0/3.0)
    RK_adv_coeff =1.0
  elseif(irk_e==4) then
    RK_cum_coeff =(1.0/6.0)
    RK_adv_coeff =1.0 ! adv is not needed for 4th RK step
  else
    write(stdout,*) "In pushe_RK4, irk_e=",irk_e,", is weird. Please Check It!!"
    call flush(stdout)
  endif
  
  call electron_locator
  call n_T_u_group_e             ! get Te at electron position

  !  mu(1:me)=zeon(6,1:me)**2

  ! gather e_field at electron GC positions using linear interpolation.
  ! No FLR effect for electrons.
  ! The following line is an OpenMP directive for loop-level parallelism
  ! on shared memory machines (see http://www.openmp.org).
  !$omp parallel do private(m,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11)
  do m=1,me
     wz1=wzeon(m)
     wz0=1.0_wp-wz1

     ij=jteon0(m)
     wp0=1.0_wp-wpeon(m)
     wt00=1.0_wp-wteon0(m)
     e1=wp0*wt00*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
     e2=wp0*wt00*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
     e3=wp0*wt00*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
     
     ij=ij+1
     wt10=1.0_wp-wt00
     e1=e1+wp0*wt10*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
     e2=e2+wp0*wt10*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
     e3=e3+wp0*wt10*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
     
     ij=jteon1(m)
     wp1=1.0_wp-wp0
     wt01=1.0_wp-wteon1(m)
     e1=e1+wp1*wt01*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
     e2=e2+wp1*wt01*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
     e3=e3+wp1*wt01*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
     
     ij=ij+1
     wt11=1.0_wp-wt01
     e1=e1+wp1*wt11*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
     e2=e2+wp1*wt11*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
     e3=e3+wp1*wt11*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

     wpi(1,m)=e1
     wpi(2,m)=e2
     wpi(3,m)=e3

     ! (MGYOO: 05/01/2019)
     ! For RK4 method, the time here must be updated as well.
     !MGYOO (10/30/2018; to use analytic profile of ga0(a))
    call Get_MP_alpha_amplitude(zeon(1,m), zeon(2,m), zeon(3,m), istep*tstep, ga(m), wga(1:3,m), gat(m))

  enddo

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
  do m=1,me
     r=zeon(1,m)-a0      
     r=r*dela
     ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
     fip=r-real(ip,wp)
     Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
  enddo
  !  ------------------------------------------------------------


  if(electromagnetic == 1)then
    !m=gcmotion_eon_em(dgr,da,dgq,dgf,zeon,wpi,ga,wga,gat,Er_particle,nonlinear,me,zeondim)
    !m=gcmotion_eon_em_MGY(dgr,da,dgq,dgf,zeon,wpi,ga,wga,gat,Er_particle,nonlinear,me,zeondim)
    m=gcmotion_eon_em_MGY_vec_form(dgr,da,dgq,dgf,zeon,wpi,ga,wga,gat,Er_particle,nonlinear,me,zeondim)
  else
    m=gcmotion_eon(dgr,da,dgq,dgf,zeon,wpi,Er_particle,nonlinear,me,zeondim)
  endif


  !  call n_T_u_group_e  ! prepare n, T, u_//, dlnn/da, dlnT/da, etc at electrons.
                       ! moved to collision_eon.f90 -- caution: in this case, BB must be obtained by call esigetb etc.!!!

  call Ti_at_electrons ! NOTE: array backT is borrowed to store Ti_at_electrons;
                       ! WARNING: must make sure there is no memory conflict!!!
                       ! backT is used in PUSHI and released afterward

  !     OUTPUT: dgr=drho/dt;        da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt;    dgf=dzeta/dt;      
  !             ghq(:)=da/dt|_EXB;  Ta(:)=dtheta/dt|_EXB;  Pa(:)=drho/dt|_EXB;
  !             BB(:)=B;            Ba(:)=dB/da;           Bq(:)=dB/dtheta;       
  !             gha(:)=1/|grad a|   F(:) =F
  !             wpi(3,:)=dphi/dzeta (NO change); no change for: zeon(1,2,4,:), wpi(1,2,:) and mu;

  !  if(nonlinear>100.0_wp) then            ! for nonlinear run
  if(wimcompressibility==1) then            ! for implementation of weight imcompressibility 
    call weight_field_eon(wfield)      ! calculate weight-field uisng zeon0 for irk=1 and zeon at half-step for irk=2
                                      ! used BB obtained from gcmotion, and electron locations are obtain here in pushe 
                                      ! which are used to calculate wfield in phase space grids and gether wfield at electron (next)
  end if

  dminv=0.0d0
  dmaxv=0.0d0

  IF(ifelectron==2.and.ifdphi==2) then  ! ifdphi==2 as the similar way as the ifdphi=2 for RK2 of GTS

    ! update GC position
    ! Another OpenMP parallel loop...
    !$omp parallel do private (m,mu,energy,kappa,wdrift,wgradB,vdr,wdrive,wdot,xdum,ydum,T_inverse,&
    !$omp& r,ip,Ea_00,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,wdphidt,wdphidt2,muB,fip,&
    !$omp& instep,outstep,fdamp,vpara,v2,avw,E_f)
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


      ! calculate d(delta_phi)/dt at particles along trajactory
      wz1=wzeon(m)
      wz0=1.0_wp-wz1

      ij=jteon0(m)
      wp0=1.0_wp-wpeon(m)
      wt00=1.0_wp-wteon0(m)
      e1=   wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=ij+1
      wt10=1.0_wp-wt00
      e1=e1+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=jteon1(m)
      wp1=1.0_wp-wp0
      wt01=1.0_wp-wteon1(m)
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

      E_f=(e1-real(1-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! Adiabatic term (q*dphi_NZ/Te)
     

      if(zeon0(1,m).le.lmax.or.zeon0(1,m).ge.rmin) then
        instep =smoothstep(lmin,lmax,zeon(1,m))
        outstep=smoothstep(rmin,rmax,zeon(1,m))
        !  fip=smoothstep(0.20_wp,0.25_wp,zeon(1,m))*(1.0_wp-smoothstep(0.75_wp,0.80_wp,zeon(1,m)))
        fip=instep*(1.0_wp-outstep)
        !  wdrive=wdrive*fip

        fdamp=(1.0_wp-fip)*fdamp0
      else
        fip=1.0_wp
        fdamp=0.0_wp
      end if

      if(irk_e==1) then
        !======================================================
        ! MGYOO(05/10/2019) :
        ! Please note that zeon0(8,m) is the adiabatic drive term
        ! zeon0(8,m)  = (q*dphi_NZ/Te) @ t{n}, x{n}
        zeon0(8,m) = E_f
        !======================================================
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
      wdot=(1.0_wp-E_f)*wdot



      if(wimcompressibility==1) then     ! for weight imcompressibility
      ! calculate average weight <w> at particles using wieght-field wfield
        wz1=wzeon(m)
        wz0=1.0_wp-wz1

        ij=jteon0(m)
        wp0=1.0_wp-wpeon(m)

        wt00=1.0_wp-wteon0(m)
        e1=wp0*wt00*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
        e2=wp0*wt00*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
        e3=wp0*wt00*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

        ij=ij+1
        wt10=1.0_wp-wt00
        e1=e1+wp0*wt10*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
        e2=e2+wp0*wt10*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
        e3=e3+wp0*wt10*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

        ij=jteon1(m)
        wp1=1.0_wp-wp0
        wt01=1.0_wp-wteon1(m)
        e1=e1+wp1*wt01*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
        e2=e2+wp1*wt01*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
        e3=e3+wp1*wt01*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

        ij=ij+1
        wt11=1.0_wp-wt01
        e1=e1+wp1*wt11*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
        e2=e2+wp1*wt11*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
        e3=e3+wp1*wt11*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

        avw=e1+e2*vpara+e3*v2


        !     if(mype==100.and.m==me/2) print*, 'AVW=', avw, zeon(5,m), mype
        !     if(mype==0.and.m==me/4) print*, 'AVW=', avw, zeon(5,m), mype
        !     if(mype==200.and.m==me/8) print*, 'AVW=', avw, zeon(5,m), mype
      end if

      if(wimcompressibility==1) then
        now_w = avw
      else
        now_w = zeon(5,m)
      endif

      !==================================================
      ! MGYOO(03/27/2019)
      ! Note that the following routine should be upgraded to calculate correct RK4
      !==================================================
      !==================================================
      !MGYOO(01/14/2019)
      !Note that f0 is toroidally symmetric so that grad_lnf0_gf=0 and we do not need vd_para_gf or vd_perp_gf terms to calculate weight evolution
      if(irk_e==4) then
        if(ifneoclassical==1) then
          !Caution!! Following routine is only available for neoclassical case which use full velocity.
          call esiget_Intp_B_vec_cov_and_cntv(zeon(1,m), zeon(2,m), zeon(3,m), istep*tstep, tot_B_vec_cov, tot_B_vec_cntv)
          vd_para_mag = (vdr*tot_B_vec_cov(1) + vdq*tot_B_vec_cov(2) + vdgf*tot_B_vec_cov(3))/BB(m)
          vd_para_a_cntv = vd_para_mag*(tot_B_vec_cntv(1)/BB(m))
          vd_para_gq_cntv = vd_para_mag*(tot_B_vec_cntv(2)/BB(m))
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

        we_vpara(m) = we_vpara(m) + tstep*( (1.0-nonlinear*now_w)*wdot_vpara*fip)
        we_vperp(m) = we_vperp(m) + tstep*( (1.0-nonlinear*now_w)*wdot_vperp*fip)
        we_accel(m) = we_accel(m) + tstep*( (1.0-nonlinear*now_w)*wdot_accel*fip)
        we_fdamp(m) = we_fdamp(m) + tstep*(-fdamp*zeon(5,m))
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
        ! zeon(5,m) = zeon0(5,m) +RK_adv_coeff*tstep*((zeon0(6,m)-nonlinear*now_w)*(wdot +zeon(8,m))*fip-fdamp*zeon(5,m))
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
      zeon(8,m) = wdot

    enddo

    !   fluxes diagnose
    if(idiag_e==0) then  ! istep=n*ndiag, irk=1, i=1, and irk_e=1
      !$omp parallel do private(m)
        do m=1,me
            wpi(1,m)=wpi(1,m)*gha(m)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
        end do
  
        call get_fluxes_e(wpi)
  
        !   im=mod(istep,100)
        im=mod(istep,ndiag)
        if(im==0) then
            call get_Te_fluctuation(wpi(2,:))
            ! write(stdout,*) 'irk=',irk,', energy(1:me:100) in keV'
            ! write(stdout,*) (wpi(2,i)*Tem_00,i=1,me,me/20)
            ! call flush(stdout)
        end if
    endif

    call shifte
    call electron_locator
    call n_T_u_group_e
    zeon00(5,1:me)=zeon(5,1:me)    ! use for temporary storage -- check again April, 2013

    do m=1,me
      ! calculate d(delta_phi)/dt at particles along trajactory

      T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti


      wz1=wzeon(m)
      wz0=1.0_wp-wz1
      ij=jteon0(m)
      wp0=1.0_wp-wpeon(m)
      wt00=1.0_wp-wteon0(m)
      e1=   wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))
      ij=ij+1
      wt10=1.0_wp-wt00
      e1=e1+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))
      ij=jteon1(m)
      wp1=1.0_wp-wp0
      wt01=1.0_wp-wteon1(m)
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
      E_f=(e1-real(1-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! Adiabatic term (q*dphi_NZ/Te)

      if(zeon0(1,m).le.lmax.or.zeon0(1,m).ge.rmin) then
        instep =smoothstep(lmin,lmax,zeon(1,m))
        outstep=smoothstep(rmin,rmax,zeon(1,m))
        fip=instep*(1.0_wp-outstep)
        !  wdrive=wdrive*fip
        fdamp=(1.0_wp-fip)*fdamp0
     else
        fip=1.0_wp
        fdamp=0.0_wp
     end if

      wdot = zeon(8,m)

      !======================================================
      ! MGYOO(05/10/2019): Here, zeon0(8,m) differs from that of the irk_e=1 stage of pusheon
      ! zeon0(8,m)= (q*dphi_NZ/Te) @ t{n}, x{n}, where {...} indicates the time level      
      !======================================================
      ! Here, E_f = (q*dphi_NZ/Te) @ t{n}, x{n+k}, where k is RK_adv_coeff (eg. k=0.5 for irk_e=1,2 | k=1 for irk_e=3)
      wdphidt = (E_f-zeon0(8,m))/(RK_adv_coeff*tstep)  
      ! Advance particle information for the next RK step
      ! weight advanced including dphi/dt term (so the driving terms have (wdot+dphidt))
      zeon(5,m) = zeon0(5,m) +RK_adv_coeff*tstep*((zeon0(6,m)-nonlinear*zeon00(5,m))*(wdot +wdphidt)*fip-fdamp*zeon(5,m))
    enddo


    call charge_eon
    call poisson( solver1, masssolver )
    call smooth_PHI

    
    do m=1,me
      ! calculate d(delta_phi)/dt at particles along trajactory

      T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti

      wz1=wzeon(m)
      wz0=1.0_wp-wz1
      ij=jteon0(m)
      wp0=1.0_wp-wpeon(m)
      wt00=1.0_wp-wteon0(m)
      e1=   wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))
      ij=ij+1
      wt10=1.0_wp-wt00
      e1=e1+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))
      ij=jteon1(m)
      wp1=1.0_wp-wp0
      wt01=1.0_wp-wteon1(m)
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
      E_f=(e1-real(1-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0 ! Adiabatic term (q*dphi_NZ/Te)

      if(zeon0(1,m).le.lmax.or.zeon0(1,m).ge.rmin) then
        instep =smoothstep(lmin,lmax,zeon(1,m))
        outstep=smoothstep(rmin,rmax,zeon(1,m))
        fip=instep*(1.0_wp-outstep)
        !  wdrive=wdrive*fip
        fdamp=(1.0_wp-fip)*fdamp0
     else
        fip=1.0_wp
        fdamp=0.0_wp
     end if

      wdot = zeon(8,m)

      !======================================================
      ! MGYOO(05/10/2019): Here, zeon0(8,m) differs from that of the irk_e=1 stage of pusheon
      ! zeon0(8,m)= (q*dphi_NZ/Te) @ t{n}, x{n}, where {...} indicates the time level      
      !======================================================
      ! Here, E_f = (q*dphi_NZ/Te) @ t{n}, x{n+k}, where k is RK_adv_coeff (eg. k=0.5 for irk_e=1,2 | k=1 for irk_e=3)
      wdphidt = (E_f-zeon0(8,m))/(RK_adv_coeff*tstep)  
      ! Advance particle information for the next RK step

      ! Calculate dADdt and evolve weight
      ! Accumulate weight including dphi/dt term (so the driving terms have (wdot+dphidt))     
      RK_cum_zeon(5,m) = RK_cum_zeon(5,m) +RK_cum_coeff*tstep*((zeon0(6,m)-nonlinear*zeon00(5,m))*(wdot +wdphidt)*fip-fdamp*zeon(5,m))

     ! Advance particle information for the next RK step
      if(irk_e<4) then
        ! weight advanced including dphi/dt term (so the driving terms have (wdot+dphidt))
        zeon(5,m) = zeon0(5,m) +RK_adv_coeff*tstep*((zeon0(6,m)-nonlinear*zeon00(5,m))*(wdot +wdphidt)*fip-fdamp*zeon(5,m))
      else if(irk_e==4) then !irk_i==4
        ! Get the future value at time level (n+1): f_(n+1) = f_(n) + (1/6)*(k1+2*k2+2*k3+k4) = f_(n)+RK_cum_value
        zeon(5,m) = zeon0(5,m) +RK_cum_zeon(5,m)
      endif


      !==================================================
      !MGYOO(05/10/2019)
      ! Cal weight evolution due to the adiabatic response (d/dt(q*dphi/Te))
      if(irk_e==4) then
        we_adiab(m) = we_adiab(m) + tstep*(zeon0(6,m)-nonlinear*zeon(8,m))*wdphidt    
      endif
      
    enddo  



    if(irk_e==4)then
    ! out of boundary particle
    !$omp parallel do private(m)
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

  ELSE    !IF(ifelectron==2.and.ifdphi==2) then  ! ifdphi==2 as the similar way as the ifdphi=2 for RK2 of GTS
    write(stdout,*) "In pushe_RK4, wrong flags are observed!! Please Check It!!"
    write(stdout,*) "Pushe_RK4 routine cannot deal with, ifelectron=", ifelectron,", ifdphi=", ifdphi
    call flush(stdout)
  ENDIF   !IF(ifelectron==2.and.ifdphi==2) then  ! ifdphi==2 as the similar way as the ifdphi=2 for RK2 of GTS

end subroutine pusheon_RK4

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
     tdum=(tdum-aint(tdum))*delt(im)
     j00=max(0,min(mtheta(im)-1,int(tdum)))
     jteon0(m)=igrid(im)+j00  ! pay attention to j00=0 case
     wteon0(m)=tdum-real(j00,wp)

! outer flux surface
     im=ii+1
     tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
     tdum=(tdum-aint(tdum))*delt(im)
     j01=max(0,min(mtheta(im)-1,int(tdum)))
     jteon1(m)=igrid(im)+j01  ! pay attention to j01=0 case
     wteon1(m)=tdum-real(j01,wp)
  enddo	

end subroutine electron_locator
