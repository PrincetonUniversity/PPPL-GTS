Subroutine setup

  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface
  use spl_profiles
  use profile_array
  implicit none
  integer i,j,k,ierror,ij,n(num_mode),m(num_mode),mid_theta,ip,jt,indp,indt
  real(wp) tdum,r,q,sint,dtheta_dx,rhoi,b,zdum,delr,delt,wt,a_minor,&
           N_on_ring, Vol_ring

  real(wp) circle,error,psimin,gyroradius_0,eachv,g,bcircle,dpsi,psimax,&
       a_ref,zeta_ref,psiout,pdum,totalv,inv_jaco,munstable,nunstable,gradth,&
       munstable1,nunstable1,inv_grada,Lc,Lci,Lce,Lcn,Lcu,Tcc,backnc,Tau,Tau_e,vthc_e,dela,fip
  integer intgral,ninbin,i_itera

  integer, external :: metric_elements

  CHARACTER(LEN=10) date, time
  CHARACTER(LEN=13) cdum
  CHARACTER(LEN=16) cdum1

! the following definition is for volume integral using NAG 
  integer minpts,ifail
  integer,parameter :: maxpts=1000*3,lenwrk=(3+2)*(1+maxpts/33)
  real(wp) low_lim(3),up_lim(3),wrkstr(lenwrk),acc,error_8,finval,low1(2),up1(2)
  real(wp) rtden_allnodes(nrbin),rtden_global
  real(wp),external :: f_jaco,object_1,object_5,reducedjacobian,FUN1,FUN3,FUN3e,FUN3n,FUN3u, &
                       rho2a_zero,dsda,s2a_zero,FUN_J,smoothstep

  namelist /run_parameters/ micell,misum,mpsi,mthetamax,mzetamax,numberpe,&
       nonlinear,ndiag,tstep,ifelectron,delr,delt,mid_theta,mgrid,mtdiag,&
       a,a0,a1,q0,q1,q2,kappati,kappate,kappan,ulength,utime,gyroradius,&
       irun,nbound,mode00,paranl,iload,icollision,ecollision

  integer esiread,species,get4exb,get4exb12

  real(wp), dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
  real(wp) jacobian
  real(wp) rho_0,rhos1,rhos2,rhos3,magnetic_shear,dqda

  real(wp) grada,rdum,rfac,kappa,dlnT,dlnn,dlnT_e, dlnu, dlnu0, &
           maxx,maxy,maxz,minx,miny,minz

  real(wp), external :: D01AHF,arclength,object_10,temperature,&
                        density,object_9s,dT_da,dn_da,temperature_e,dT_da_e,& 
                        domega_da , object_9,  object_11, u_parallel

  real(wp), dimension(:), allocatable :: workarray_4_ETG
  real(wp), dimension(:), allocatable :: coord_X,coord_Y,coord_Z,&
                                         coord_R,coord_Phi

  integer create_directory

  integer,parameter :: nfft=768
  real(wp), dimension(:,:), allocatable :: jaco4spectrum

  common/stupidpass/pdum
  common/fun1d/Tcc,r

  real(wp) vb(2),vp(nprobe,2)

#ifdef _OPENMP
  integer nthreads,omp_get_num_threads
#endif

! Read run parameters from input file
  call read_input_params

!  kappan=kappati*kappan        ! =R_0*<grad(n)/n>=R_0/L_n
  aelectron=1.0/aelectron      ! aelectron=1.0/1837.0

!  r0=rmajor    !wx  normalized
!  b0=b_axis    !wx  normalized

!  NO. 1 constant
  pi=4.0_wp*atan(1.0_wp)
  twopi=2.0_wp*pi

  dtheta=dtheta*twopi
  tprobe=tprobe*pi

! mass ratio and charge ratio
    msm0=amasse/pmass0
    m0ms=pmass0/amasse
    csc0=charge/charg0
    c0cs=charg0/charge
    for_rho=msm0*c0cs
    for_E=m0ms*csc0**2
    msm0_i=amassi/pmass0
    m0ms_i=pmass0/amassi
    csc0_i=chargi/charg0
    c0cs_i=charg0/chargi

! Set up the particle decomposition within each toroidal domain
  call set_particle_decomp

! make sure mpsi is even
  mpsi=2*(mpsi/2)

! make sure that mzetamax is a multiple of ntoroidal
  mzetamax=ntoroidal*max(1,int(real(mzetamax)/real(ntoroidal)+0.5_wp))

! Domain decomposition in toroidal direction.
  zetamin=2.0_wp*pi*real(toroidal_domain_location)/real(ntoroidal)
  zetamax=2.0_wp*pi*real(toroidal_domain_location+1)/real(ntoroidal)

!  write(0,*)mype,' in set_particle_decomp: mzeta=',mzeta,'  zetamin=',&
!            zetamin,'  zetamax=',zetamax

! While the root process is busy reading the equilibrium and profile files in
! the call to esiread(), we use the last process to create the directories
! where the restart files will go. Their are npartdom such directories but we
! need to use a single process to create them since the top directory "RESTART"
! needs to exist before the other can be created. See the subroutines
! restart_write and restart_read for more details on RESTART files.
  if(mype==numberpe-1)then
     ierror=create_directory('RESTART'//char(0))
     do i=0,npartdom-1
         write(cdum,'("RESTART/D",i4.4)')i
        ierror=create_directory(trim(cdum)//char(0))
     enddo
  endif

! c2.    Reading data file (s) File names: <A-Z,a-z>,<0-9>,<_>,<.>
  k=esiread(B_0,R_0,eq_data_name) ! B_0=B on axis; R_0=major radius (normalized)

  if(k/=0) then
     write(0,*)'bad file for reading'
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

! pass mass ratio, charge ratio to gcmotion_e
  k=species(amasse,pmass0,charge,charg0)

! numerical constant
  mstep=max(2,mstep)
  mstepall=0 ! initailize ??? (ma)
!  msnap=min(msnap,mstep/ndiag)
  idiag1=ipsi
  idiag2=ipsi
  if(nonlinear < 0.5) then   ! linear run
     nonlinear=0.0_wp
     mode00=0                ! exclude zonal flow (in field.90 and smooth.f90) for linear run
     idiag1=1
     idiag2=mpsi
  endif

! allocate memory
  allocate (psi_grid(0:mpsi),plenth(0:mpsi),r_grid(0:mpsi),q_grid(0:mpsi),&
            den_grid(0:mpsi),vth_grid(0:mpsi),vth_inv(0:mpsi),Te_grid(0:mpsi), &
            rho_grid(0:mpsi),s_grid(0:mpsi), &
            eta_grid(0:mpsi),Er_grid(0:mpsi),Er_grid0(0:mpsi),tite(0:mpsi),dadtdz(0:mpsi),ddT_da(0:mpsi))

!reference zeta point
  zeta_ref=0.0_wp

! simulation domains
!  a0=a0*a                  ! inner boundary, unit=R_0
!  a1=a1*a                  ! outer boundary, unit=R_0

  if(mype==0)then
     write(0,*) 'a1=',a1, 'a0=', a0, 'a=',a
  end if

!
!  if plasma=0 and q2<0 --> using another analytical T profile specified by grad(T) profile:
!          R_0/L_T ~ grad(lnT) ~ exp{-[(r-rc)/rw]**q0}
!
!  rc=rc*(a0+a1)
!
  rw=1.0_wp/(rw*(a1-a0))
  flow2=1.0_wp/(flow2*(a1-a0))
  temconst=1.0_wp/(temconst*(a1-a0))
  Bu=1.0_wp/(Bu*(a1-a0))

  if(plasma==0.and.q2.lt.0.0) then 
     allocate(Ti_node(0:jr),Te_node(0:jr),ne_node(0:jr),Iomg_node(0:jr))
     rdum=(a1-a0)/real(jr,wp)
     do i=0,jr
        pdum=a0+rdum*real(i,wp)
        error_8=1.0e-6_wp
        ifail=1
        Ti_node(i)=D01AHF(rc,pdum,error_8,minpts,acc,FUN3,-1,ifail)     !1D integral over r
        ifail=1
        Te_node(i)=D01AHF(rc,pdum,error_8,minpts,acc,FUN3e,-1,ifail)    !1D integral over r
        ifail=1
        ne_node(i)=D01AHF(rc,pdum,error_8,minpts,acc,FUN3n,-1,ifail)    !1D integral over r
        ifail=1
        Iomg_node(i)=D01AHF(rc,pdum,error_8,minpts,acc,FUN3u,-1,ifail)  !1D integral over r

     end do

     call general_flux_average(rc,object_9,q1)  ! q1=<|grad(a)|>
     flow1=-kappate/R_0/q1    ! for Te
     kappan=-kappan/R_0/q1    ! for n
     Au=-Au/R_0/q1            ! for toroidal rotation omega_t
     q1=-kappati/R_0/q1       ! for Ti

  end if

  aaa(1)=1.0_wp
  the(1)=0.0_wp
  call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
  rho_0=sqrt(gFaj(1)/B_0)    ! minor radius defined in (5.35) of HH; da/dr=1/rho_0

! equilibrium unit: length (unit=cm) and time (unit=second) unit
  ulength=R_00                  !old:  ulength=r0
  utime=pmass0*c/(charg0*B_00)  !utime=aion/(9580*qion*B_00)
  uveloci=ulength/utime
  E_00=pmass0*uveloci**2        ! energy unit
  Tem_00=E_00/(1000._wp*trans)  ! temperature unit (unit for input Ti is KeV)


if(mype==0) print*, 'FFFFFFF CCCCCCCC', uveloci*charg0, uveloci*charg0/3.0d5


      if(plasma==1) then           !Read experimental profiles

         call read_profiles(profile_data_name)         ! for JET data output

         Ti_50=Ti_50/Tem_00    ! normalized temperature for which the normalized
                               ! thermal velocity is vth=sqrt(T_i)
         Te_50=Te_50*m0ms/Tem_00  ! normalized temperature corresponding to
                                  ! electron thermal velocity of vth=sqrt(T_e)
         O_50=O_50*utime       ! normalized toroidal rotation angular velocity,
                               ! in radians (rad)
! --------------------
!    O_50=0.0d0
! --------------------

!         if(ETG==1) then  ! exchange Ti and Te so that temperature()-->Te and temperature_e()-->Ti;
!                          ! if using model profiles (plasma=0), Ti=Te is assumed, no exchange is needed
!            allocate (workarray_4_ETG(ndum+2))
!            workarray_4_ETG=Ti_50
!            Ti_50=Te_50
!            Te_50=workarray_4_ETG
!         end if

         call spl_init

!         if(mype==1) call get_profiles
 
! prepare data on nodes for linear interpolation of profiles
         if(linterp_from_nodes==1) call get_profiles_on_nodes

      end if

!  a_ref=0.51         ! radial reference point
!  a_ref=0.65         ! used for NSTX-115821
!  a_ref=0.55_wp         ! used for NSTX-116960
!  a_ref=0.50_wp         ! used for NSTX-112992
!  a_ref=0.60         ! used for NSTX-115821 2nd run
!  a_ref=0.35_wp      ! used for DIIID-125238280 
!  a_ref=0.55_wp         ! esiP64_mono.D3D
!  a_ref=0.62         ! used for NSTX-200155 
  a_ref=0.42_wp      ! used for DIIID-133672-314
  a_ref=0.45_wp      ! used for DIIID-125236_290
  a_ref=0.42_wp      ! used for DIIID-125236_360
  a_ref=0.47_wp      ! used for DIIID-125499_489
  a_ref=0.45_wp      ! used for DIIID-131910
  a_ref=0.56_wp      ! used for DIIID-125499_229
  a_ref=0.50_wp         ! used for NSTX-121012
  a_ref=0.55_wp      ! used for DIIID-122338b02 2nd run

  a_ref=0.55_wp         ! used for DIIID-125229
  a_ref=0.35_wp      ! used for DIIID-129533_405
  a_ref=0.40_wp      ! used for DIIID-129533_405 eq., but with artificial profiles
  a_ref=0.55_wp         ! used for model equilibrium and profiles

!
!  Initially, T_0 and T_e0 are specified in INPUT.d as temperatures at a=a_ref --> change made in Oct. 09
!
  if(plasma==0) then
     if(q2.lt.0.0) then  ! for model profile II
        a_ref=rc
!
!  for model profiles; not used if plasma=1
!  T_0 and T_e0 are temperatures at a=a_ref
!
        T_0=T_0/Tem_00   ! normalized temperature for which the normalized
                      ! thermal velocity(=sqrt(T/m)) is vth=sqrt(T_i)
        T_e0=T_e0/Tem_00 ! electron temperature profile specified too
        T_e0=T_e0*m0ms   ! normalized temperature for which the normalized
                      ! electron thermal velocity(=sqrt(T/m)) is vth=sqrt(T_e)
        Tc=temperature(a_ref,0.0_wp)  ! take T(a=0.5) as reference temperature
        Tc_e=temperature_e(a_ref,0.0_wp)  ! take T(a=0.5) as reference temperature

      !    input u_0 is normalized by thermal velocity [ sqrt(2T/m) ] at r=rc
        u_0=u_0*sqrt(2.0_wp*Tc)   ! this u_0 is u// at r=rc, in unit uveloci now
        Au=Au*(0.1_wp*sqrt(2.0_wp*Tc)/R_0)  ! this roughly sets R_0/L_omg at r=rc ~ Au value originally 
                                            ! specified in INPUT.d with a reference u_// (or omega_t) in 
                                            ! R_0/L_omg taking to be 10% of vth at r=rc 

     else    ! for model profile I
             ! need both T_0 and Tc (T_e0 and Tc_e) to specified temperature profile

        Tc=T_0           ! Initially, T_0 is specified as Ti at a=a_ref, i.e, Tc

        tdum=At+Bt*exp(-dTemi*a_ref**xt)
        T_0=T_0/tdum                              ! --> Now T_0=Ti at a=0 (magnetic axis)
        T_0=T_0/Tem_00   ! normalized temperature for which the normalized
                         ! thermal velocity(=sqrt(T/m)) is vth=sqrt(T_i)
        Tc=Tc/Tem_00

        Tc_e=T_e0        ! Initially, T_e0 is specified as Te at a=a_ref, i.e, Tc_e

        tdum=At_e+Bt_e*exp(-dTeme*a_ref**xt_e)
        T_e0=T_e0/tdum                            ! --> Now T_e0=Te at a=0 (magnetic axis)

        T_e0=T_e0/Tem_00 ! electron temperature profile specified too
        T_e0=T_e0*m0ms   ! normalized temperature for which the normalized
                      ! electron thermal velocity(=sqrt(T/m)) is vth=sqrt(T_e)
        Tc_e=Tc_e/Tem_00
        Tc_e=Tc_e*m0ms

        tdum=An+Bn*exp(-dn*a_ref**xn)  ! Initially, n_0 is specified as ne at a=a_ref, i.e, nc; change made in Nov. 09
        n_0=n_0/tdum                   ! --> Now T_e0=Te at a=0 (magnetic axis)
        !  n_0=4.1d13

      !    input u_0 is normalized by thermal velocity [ sqrt(2T/m) ] at axis
        u_0=u_0*sqrt(2.0_wp*T_0)   ! this u_0 is u// at axis, in unit uveloci now

      ! ******************************
      ! set rigid rotation for q2>0 case
      !  r=domega_da(0.54,0.0,zdum)
      !  u_0=zdum
      !  Au=1.0
      !  Bu=0.0
      !  du=0.0
      ! ******************************

     end if

     u_0=u_0/R_0    ! --> if q2>0, omega_t=omega_0*[Au+Bu*exp(-du*r**xu)] with omega_0=u_0
                    ! --> if q2<0, sepcify d(omega_t)/da=Au*exp{-[(r-rc)/Bu]**xu}  added in 12/09
                    !     u_//=(I*omega_t)/B
  end if

  if(plasma==1) then
     Tc=temperature(a_ref,0.0_wp)  ! take T(a=a_ref) as reference temperature
     Tc_e=temperature_e(a_ref,0.0_wp)
  end if
     
! used for D129533  -- simulation domain [0.05, 0.8]
  lmin=0.10_wp
  lmax=0.15_wp
  rmin=0.70_wp
  rmax=0.75_wp

! used for test
  lmin=0.10_wp
  lmax=0.15_wp
  rmin=0.70_wp
  rmax=0.75_wp

! used for DIIID-125229, and model equilibrium and profiles -- simulation domain [0.1, 0.9]
  lmin=0.20_wp
  lmax=0.25_wp
  rmin=0.75_wp
  rmax=0.80_wp

! used for testing flow in the inner core region due to CTEM in outer core region - simulation domain [0.1, 0.9]
!  lmin=0.15_wp
!  lmax=0.20_wp
!  rmin=0.80_wp
!  rmax=0.85_wp

! used for DIIID-125229,  no boundary effect -- simulation domain [0.1, 0.9]
!  lmin=0.10_wp
!  lmax=0.11_wp
!  rmin=0.89_wp
!  rmax=0.90_wp

  Tcc=Tc                        ! make a copy

  vthc=sqrt(Tc)                 ! normalized reference thermal velocity
  vthc_e=sqrt(Tc_e)             ! normalized reference thermal velocity of electron
  backnc=density(a_ref,0.0_wp)  ! plasma density at a_ref in cm^-3

  Cs=sqrt(Tc_e*msm0)                 ! normalized Cs

  if(mype==1) call get_profiles

! rescale radial coordinator  rho=a*sqrt(Tc/T)
  rho0=a0*sqrt(Tc/temperature(a0,0.0_wp))
  rho1=a1*sqrt(Tc/temperature(a1,0.0_wp))

  error_8=1.0e-5_wp
  IFAIL=1
  s0=D01AHF(0.0_wp,a0,error_8,minpts,acc,dsda,-1,IFAIL)

  error_8=1.0e-5_wp
  IFAIL=1
  s1=D01AHF(0.0_wp,a1,error_8,minpts,acc,dsda,-1,IFAIL)

if(mype==0) write(0,*) 'Tc=', Tc*Tem_00,  'in KeV', backnc, &
            ' a_ref=',a_ref,'  vthc=',vthc,'  s0=',s0,'  s1=',s1

!   vthi=sqrt(temconst*trans/(aion*pmass0)) ! ion thermal velocity=sqrt(T/m)
!   vthi=vthi/uveloci                          ! normalized vthi

!  ion thermal gyroradius=sqrt(T/m)/(eB_0/mc) with B_0=B on axis
!  ion thermal gyroradius in equilibrium unit (for original normalization B_0=1)
!  gyroradius=9.79e5*sqrt(temconst/aion)*utime/ulength

  gyroradius=vthc*aion/qion/B_0   ! normalized ion gyroradius

  call general_flux_average(a_ref,object_9s,inv_grada)  ! inv_grada=<|1/grad_a|>
  Lci=inv_grada/abs(dT_da(a_ref,0.0_wp)/Tc)  ! Lci=<Ti gradient scale length> at a=a_ref, normalized

! for shear-flow-driven turbulence, use
  rdum=domega_da(a_ref,0.0,zdum)
  Lcu=inv_grada/abs(rdum/zdum)  ! Lcu=<omega_t gradient scale length> at a=a_ref, normalized

! for TEM turbulence driven by grad(Te)
  Lce=inv_grada/abs(dT_da_e(a_ref,0.0_wp)/temperature_e(a_ref,0.0_wp))  ! Lce=<Te gradient scale length> at a=a_ref, normalized

! for TEM turbulence driven by grad(ne)
  Lcn=inv_grada/abs(dn_da(a_ref,0.0_wp)/backnc)  ! Lcn=<ne gradient scale length> at a=a_ref, normalized

!  Lc=Lce    ! for D125229
!  Lc=Lce    ! for D129533, but with artificial profiles (grad-Te driven case)
  Lc=Lcn    ! for D129533, but with artificial profiles (grad-n driven case)
!  Lc=Lci    ! for D129533
!  Lc=Lcu    ! for toroidal shear flow driven KH instability and turbulence

if(mype==0) write(0,*) 'R_0/<L_T>=', R_0/Lc   ! 6.9 is often used in cyclone case 

!  tstep=tstep/(vthi*kappati)*R_0  ! normalized delta_t
   
   tstep=tstep*Lc/vthc             ! normalized delta_t

!  kappati=kappati/R_0             ! now -dlnT_i/da=kappati/<|grad_a|>
!  kappate=kappate/R_0             ! now -dlnT_e/da=kappate/<|grad_a|>
!  kappan=kappan/R_0               ! now -dlnn/da=kappan/<|grad_a|>
  
  call get_ave_grada_table        ! get ready for |grad_a| table

if(mype==0) write(0,*) 'time step size=', tstep*utime, gyroradius*ulength, vthc*uveloci

! basic ion-ion collision time, Braginskii definition
  if(icollision.gt.0) then
    ! i-i collision time at reference location (in sec)
      Tau=1.0_wp/4.0_wp/pi*pmass0*(vthc*uveloci)**3.0_wp/(backnc*charg0**2*charg0**2*dlni)*pmass0
    ! taking into account different definition in vthc
      Tau=Tau*(sqrt(2.0_wp))**3
    ! artificially enlarge collision frequency by a factor of coll_enlarge
      Tau=Tau/coll_enlarge
      Tau=Tau/utime     ! normalized collision time at reference location
      dt_coll=tstep/Tau ! step size in terms of reference collision time Tau
      if(mype==0) write(0,*) 'collision is on ', ' with step size=', dt_coll, 'in tau_ii'
    ! now dt_coll*n/v^3 is dt_c/tau for each particle with v normalized; 
    ! collisions are calculated every "icollision" steps
      dt_coll=real(icollision,wp)*dt_coll*(sqrt(2.0_wp)*vthc)**3/backnc 
             ! now dt_coll*n/v^3 is dt_c/tau for each particle with v normalized;
             ! collisions are calculated every icollision steps
  endif


  if(ifelectron/=0) then
     if(ecollision.gt.0) then
    !   e-e collision time at reference location (in sec)
        Tau_e=1.0_wp/4.0_wp/pi*amasse*(vthc_e*uveloci)**3.0_wp/(backnc*charge**2*charge**2*dlni)*amasse

    !   taking into account different definition in vthc
        Tau_e=Tau_e*(sqrt(2.0_wp))**3
    ! artificially enlarge collision frequency by a factor of coll_enlarge
        Tau_e=Tau_e/coll_enlarge
        Tau_e=Tau_e/utime     ! normalized collision time at reference location
        dt_coll_e=tstep/Tau_e ! step size in terms of reference collision time Tau_e
        if(mype==0) write(0,*) 'electron collision is on ', ' with step size=', dt_coll_e, 'in tau_ee'

    ! collisions are calculated every "ecollision" steps
        dt_coll_e=real(ecollision,wp)*dt_coll_e*(sqrt(2.0_wp)*vthc_e)**3/backnc
             ! now dt_coll_e*n/v^3 is dt_c/tau for each particle with v normalized;
             ! collisions are calculated every ecollision steps
     endif
  end if


! allocate memory
  allocate (qtinv(0:mpsi),itran(0:mpsi),mtheta(0:mpsi),deltat(0:mpsi),&
       rtemi(0:mpsi),rden(0:mpsi),igrid(0:mpsi),pmarki(0:mpsi),pmarke(0:mpsi),&
       phi00(0:mpsi),phip00(0:mpsi),fluxpsi(0:mpsi),zonali(0:mpsi),&
       zonale(0:mpsi),jaco_sum(0:mpsi),eigenmode(m_poloidal,num_mode,mpsi), & 
       allnodes(0:mpsi),gaa(0:mpsi),B2_s(0:mpsi),dVda(0:mpsi),En_00(0:mpsi),&
       En_tur(0:mpsi),En_tot(0:mpsi),sn0_g(0:mpsi),dnty_on_bin_g(0:mpsi),&
       dphidt00(0:mpsi))

! grid spacing
!  deltar=(a1-a0)/real(mpsi)
!  deltar=(rho1-rho0)/real(mpsi)

  deltar=(s1-s0)/real(mpsi)  ! uniform in s

  deltaz=(zetamax-zetamin)

!    NOTE:  radial coordinate a=sqrt(pst_t/psi_t,edge)  [0,1]
  psi_grid(0)=a0             
  psi_grid(mpsi)=a1          

  rho_grid(0)=rho0
  rho_grid(mpsi)=rho1

  s_grid(0)=s0
  s_grid(mpsi)=s1

  !if(mype==0) write(0,*) 'i=0', 's_grid(0)=', s0, '  psi_grid(0)=', a0
  !!write(200+mype,*)' i=0', 's_grid(0)=', s0, '  psi_grid(0)=', a0
  do i=1,mpsi-1
!     rho_grid(i)=rho0+deltar*real(i)
!     r=rho_grid(i)

     s_grid(i)=s0+deltar*real(i)
     r=s_grid(i)

     error_8=1.0e-6_wp          !  eps=1.0d-7
     acc=0.0_wp
     IFAIL=1
!        a-grid corresponding to rho-grid
!     call C05ADF(psi_grid(i-1),a1,error_8,acc,rho2a_zero,psi_grid(i),IFAIL) 

!        a-grid corresponding to s-grid
     call C05ADF(psi_grid(i-1),a1,error_8,acc,s2a_zero,psi_grid(i),IFAIL)

     !if(mype==0) write(0,*) 'i=',i, 's_grid(i)=', s_grid(i), '  psi_grid(i)=', psi_grid(i)
     !!write(200+mype,*)' i=',i, 's_grid(i)=', s_grid(i), '  psi_grid(i)=', psi_grid(i)

  enddo
  !if(mype==0) write(0,*) 'i=', mpsi, 's_grid(i)=', s1, '  psi_grid(i)=', a1
  !!write(200+mype,*)' i=', mpsi, 's_grid(i)=', s1, '  psi_grid(i)=', a1

     !!write(200+mype,*)'  calculating plenth(i)'
! field line length on poloidal plane, stored in plenth(i)
  do i=0,mpsi
     aaa(1)=psi_grid(i)
     aaa(2)=psi_grid(i)
     the(1)=0.0_wp
     the(2)=pi
     call esigetrzb(rj,zj,Bj,aaa,the,2)
     r=abs(rj(1)-rj(2))/2.0_wp    ! minor radius
     r_grid(i)=r                  ! minor radius
     den_grid(i)=density(psi_grid(i),0.0_wp)                ! density
     vthi=temperature(psi_grid(i),0.0_wp)

! NOTE!!!!  if all procs write the same thing into the same file, it take terribly long time when using many procs
     if(mype==0) then
         write(10,5) psi_grid(i), Tc*Tc/vthi**2
         call flush(10)
     end if

     vth_grid(i)=sqrt(vthi)      ! thermal velocity sqrt(T/m)
     vth_inv(i)=1.0_wp/vth_grid(i) ! inverse of thermal velocity 

     vthe=temperature_e(psi_grid(i),0.0_wp)
     vthe=vthe*msm0     ! Te with the same unit as Ti
     tite(i)=vthi/vthe  ! Ti/Te
     Te_grid(i)=vthe    ! Te on grid with the same unit as Ti 

     rfac=rw*(aaa(1)-rc)
     rdum=(aaa(1)-a0)/djr
     j=min(jr-1,max(0,int(rdum)))
     if(mype==0) write(0,*)'sos j=', j, jr, size(flux_it),'mmmmmm=',i
     grada=abs((real(j+1,wp)-rdum)*flux_it(j)+(rdum-real(j,wp))*flux_it(j+1))   ! =<|grad_a|>
     kappa=kappati*exp(-rfac**6)  ! =-grad(lnTi)

   ! if(mype==0) write(0,*)'r=',i, rj(1), r, (rj(1)-R_0), kappa, kappa/grada
     if(mype==0) write(0,*)'r=',i, rj(1), r, (rj(1)-R_0), dT_da(aaa(1),0.0)/temperature(aaa(1),0.0_wp)
     if(mype==0) write(0,*)'r=',i, 'vth=', vth_grid(i), 'Ti/Te=', tite(i)
     error=1.0e-4_wp
     intgral=32
     circle=0.0_wp
     delt=2.0_wp*pi/real(intgral,wp)
     do j=1, intgral
        the(1)=delt*real(j-1,wp)  
        the(2)=delt*real(j,wp)  
        call esigetrzb(rj,zj,Bj,aaa,the,2)
        circle=circle+sqrt((rj(2)-rj(1))**2+(zj(2)-zj(1))**2)
     enddo
     bcircle=0.0_wp
     i_itera=0
     do while (abs(circle-bcircle)/circle .gt. error)
        bcircle=circle
        intgral=2*intgral
        circle=0.0_wp
        delt=2.0_wp*pi/real(intgral,wp)
        do j=1, intgral
           the(1)=delt*real(j-1,wp) 
           the(2)=delt*real(j,wp)
           call esigetrzb(rj,zj,Bj,aaa,the,2)
           circle=circle+sqrt((rj(2)-rj(1))**2+(zj(2)-zj(1))**2)
        enddo

        i_itera=i_itera+1
        if(i_itera.gt.12) then
          write(0,*)'check plenth(i) calculation'
          call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
        end if
     enddo
     plenth(i)=circle

     !if(mype==0) write(0,*)'iiiiiii=',i, plenth(i)
     !!write(200+mype,*)'  iiiiiii=',i, plenth(i)


     ! the arc length can also be calculated by int(1/|grad_theta|d_theta):
     ! --> seems problematic 07/06
     pdum=psi_grid(i)
     error_8=1.0e-4_wp
     ifail=1
     ! plenth(i)=D01AHF(0.0d0,twopi,error_8,minpts,acc,arclength,-1,ifail)  !1D integral of NAG
    
     munstable=0.2_wp*plenth(i)/gyroradius/(2.0_wp*pi)
     if(mype==0) write(0,*)'i=',i,'plenth(i)=',plenth(i),'analytic=',&
                            2.0_wp*pi*r,'bcircle=',bcircle,'m=',munstable

     ! calculate dVda=dV/da=int(J*dtheta*dzeta)
     error_8=1.0e-4_wp
     ifail=1
     dVda(i)=D01AHF(0.0_wp,twopi,error_8,minpts,acc,reducedjacobian,-1,ifail)  !1D integral on theta
     dVda(i)=dVda(i)*2.0_wp*pi   ! zeta is ignorable for symmetry

     ! if(mype==0) write(0,*) 'kkkkkkk=',i, plenth(i)

  enddo

  !!!close(200+mype)


!  get q profile information
  if(mype==0) then
     open(501,file='qprofile.dat')
     k=200
     do i=0,k
        aaa(1)=max(1.0e-7_wp,min(1.0_wp,1.0_wp/real(k,wp)*real(i,wp)))
        the(1)=0.0_wp
        call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
        q=-gFaj(1)/gYaj(1) 
        grada=sqrt(rqj(1)**2+zqj(1)**2)/abs(zaj(1)*rqj(1)-zqj(1)*raj(1)) !grad(a) at the=0
! calculate R_0/L_T and R_0/L_n at middle plane:
        dlnT=abs(dT_da(aaa(1),the(1))/temperature(aaa(1),the(1)))*grada*R_0 !=R_0*grad(lnT)
        dlnn=abs(dn_da(aaa(1),the(1))/density(aaa(1),the(1)))*grada*R_0     !=R_0*grad(lnn)
        rdum=domega_da(aaa(1),the(1),zdum)
        if(i==0) tdum=zdum   ! totation at magnetic axis
        dlnu=abs(rdum/zdum)*grada*R_0                                       !=R_0*grad[ln(omega_t)]
        dlnu0=abs(rdum/tdum)*grada*R_0                                       !=R_0*grad[ln(omega_t)] with omega_t at magnetic axis
        write(501,5) aaa(1), q, dlnT, dlnn, dlnT/dlnn, & ! dlnT/dlnn=eta_i
                     temperature(aaa(1),the(1))*Tem_00, dlnT/grada*(rj(1)/R_0), dlnu, rdum, zdum, dlnu0
        call flush(501)
     end do
5    format(12(1x,ES13.6))
  end if


! grid shift associated with fieldline following coordinates
  
!  tdum=plenth(mpsi)/real(mthetamax)
  tdum=0.5_wp*vth_grid(ipsi)/B_0
  k=0
  do i=0,mpsi
!     mtheta(i)=max(2,min(mthetamax,2*int(plenth(i)/tdum/2.0+0.5))) !even # poloidal grid

!     tdum=(0.8*real(mpsi-i)/real(mpsi)+0.2*real(i)/real(mpsi))*vth_grid(i)/B_0 ! for NSTX-115821
!     tdum=(0.6*real(mpsi-i)/real(mpsi)+0.4*real(i)/real(mpsi))*vth_grid(i)/B_0/real(mthetamax)
!     tdum=(0.8*real(mpsi-i)/real(mpsi)+0.2*real(i)/real(mpsi))*vth_grid(i)/B_0/real(mthetamax)   ! for NSTX shots
!     tdum=tdum*0.8       ! for NSTX-116960a05
!     tdum=tdum*0.7       ! for NSTX-121283a01-400
!     tdum=tdum*0.6       ! for NSTX-115821a05 2nd run
     tdum=(0.6_wp*real(mpsi-i,wp)/real(mpsi,wp)+0.4_wp*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)  ! for D3D 122338
     tdum=tdum*1.1_wp       ! for D3D 122338 2nd run
     tdum=tdum*2.0_wp

     tdum=(0.9_wp*real(mpsi-i,wp)/real(mpsi,wp)+0.1_wp*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)  ! for NSTX shots
     tdum=tdum*0.8_wp     ! for NSTX-116960a05
!     tdum=tdum*0.5_wp     ! test higher poloidal resolution for NSTX-116960a05

     tdum=(0.7_wp*real(mpsi-i,wp)/real(mpsi,wp)+0.3_wp*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)  
                          ! for D3D 125238280
     tdum=tdum*1.0_wp     ! for D3D 125238280

!     tdum=(0.6_wp*real(mpsi-i,wp)/real(mpsi,wp)+0.4_wp*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)
!                          ! for esiP64_mono.D3D
!     tdum=tdum*1.2_wp     ! for esiP64_mono.D3D

     tdum=(0.9*real(mpsi-i)/real(mpsi)+0.1*real(i)/real(mpsi))*vth_grid(i)/B_0/real(mthetamax)   ! for NSTX shots
     tdum=tdum*0.6       ! for NSTX-200155

     tdum=(0.6*real(mpsi-i,wp)/real(mpsi,wp)+0.4*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)   ! for model esi_mod.00
     tdum=tdum*1.2       ! for esi_mod.00

     tdum=(0.8*real(mpsi-i,wp)/real(mpsi,wp)+0.2*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)   ! for model esi_mod.00
     tdum=tdum*1.6       ! for esi_mod.00
!     tdum=tdum*2.0       ! for esi_mod.00

     tdum=(0.7_wp*real(mpsi-i,wp)/real(mpsi,wp)+0.3_wp*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)
                          ! for D3D 131910
     tdum=tdum*1.0_wp     ! for D3D 131910

     tdum=(0.6_wp*real(mpsi-i,wp)/real(mpsi,wp)+0.4_wp*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)
                          ! for D3D 125236_360
     tdum=tdum*1.1_wp     ! for D3D 125236_360

     tdum=(0.6_wp*real(mpsi-i,wp)/real(mpsi,wp)+0.4_wp*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)
                          ! for D3D_125236_290
     tdum=tdum*1.0_wp     ! for D3D_125236_290

     tdum=(0.7_wp*real(mpsi-i,wp)/real(mpsi,wp)+0.3_wp*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)
                          ! for D3D_125499_229
     tdum=tdum*1.1_wp     ! for D3D_125499_229

     tdum=(0.6_wp*real(mpsi-i,wp)/real(mpsi,wp)+0.4_wp*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)  ! for D3D 122338
     tdum=tdum*1.1_wp       ! for D3D 122338 2nd run

     tdum=(0.6*real(mpsi-i,wp)/real(mpsi,wp)+0.4*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)   ! for D125229
     tdum=tdum*1.6       ! for D125229
!     tdum=tdum*0.6       ! for D125229

!      tdum=(0.6_wp*real(mpsi-i,wp)/real(mpsi,wp)+0.4_wp*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)
!                           ! for D3D 129533_405
! !     tdum=tdum*1.1_wp     ! for D3D 129533_405
!      tdum=tdum*1.0_wp     ! for D3D 129533_405

!     tdum=(0.6*real(mpsi-i,wp)/real(mpsi,wp)+0.4*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)   ! for model esi_mod.00
!     tdum=tdum*1.2       ! for esi_mod.00 + MODEL_KE31e*

!     tdum=(0.6*real(mpsi-i,wp)/real(mpsi,wp)+0.4*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)   ! for model esi_mod.00
!     tdum=tdum*1.0       ! for esi_mod.00 + MODEL_KE31e*

!     tdum=(0.6*real(mpsi-i,wp)/real(mpsi,wp)+0.4*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)   ! for model esi_mod.00
!     tdum=tdum*1.1       ! for esiP.110421a,b,c

!     tdum=(0.6*real(mpsi-i,wp)/real(mpsi,wp)+0.4*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)   ! for model esi_mod.00
!     tdum=tdum*1.4         ! for esiP.01 -- reversed magnetic shear q=2 at r/a=0.5

!     tdum=(0.8*real(mpsi-i,wp)/real(mpsi,wp)+0.2*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)   ! for esiPq0=1.5B=2
!     tdum=tdum       ! esiPq0=1.5B=2, Ip scan

!     tdum=(0.61*real(mpsi-i,wp)/real(mpsi,wp)+0.39*real(i,wp)/real(mpsi,wp))*vth_grid(i)/B_0/real(mthetamax,wp)   ! for model esi_mod.00
!     tdum=tdum*1.5       ! for esi_mod.00   with T<0.5 kev

     mtheta(i)=max(2,2*int(plenth(i)/tdum/2.0_wp+0.5_wp)) !even # poloidal grid
     deltat(i)=2.0_wp*pi/real(mtheta(i),wp)

     k=max(k,mtheta(i))

     aaa(1)=psi_grid(i)
     the(1)=0.0_wp
     call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
     q=-gFaj(1)/gYaj(1)
     q_grid(i)=q

!  get eta_i function for maintainance of charge conservation
     grada=sqrt(rqj(1)**2+zqj(1)**2)/abs(zaj(1)*rqj(1)-zqj(1)*raj(1)) !grad(a) at the=0
! calculate R_0/L_T and R_0/L_n at middle plane:
     dlnT=abs(dT_da(aaa(1),the(1))/temperature(aaa(1),the(1)))*grada*R_0 !=R_0*grad(lnT)
     dlnn=abs(dn_da(aaa(1),the(1))/density(aaa(1),the(1)))*grada*R_0     !=R_0*grad(lnn)
     dlnT_e=abs(dT_da_e(aaa(1),the(1))/temperature_e(aaa(1),the(1)))*grada*R_0 !=R_0*grad(lnT)
     rdum=domega_da(aaa(1),the(1),zdum)
     dlnu=abs(rdum/zdum)*grada*R_0      ! =R_0*grad[ln(omega_t)]; Note: dlnu may become singular at a radial location (at omega_t=zdum=0)!!!

!     eta_grid(i)=dlnT/dlnn
!     eta_grid(i)=dlnT   ! for ITG
!     eta_grid(i)=dlnn   ! for TEM driven by grad(ne)

!     eta_grid(i)=dlnT_e    ! for TEM driven by grad(Te), such as D125229 
!     eta_grid(i)=dlnT_e    ! for TEM driven by grad(Te)
     eta_grid(i)=dlnn      ! for TEM driven by grad(ne)
!     eta_grid(i)=dlnT      ! for ITG; such as D129533 with KE
!     eta_grid(i)=abs(rdum) ! for KH turbulence driven by toroidal flow shear  ! added 2011

     eta_grid(i)=eta_grid(i)*smoothstep(lmin,lmax,aaa(1))*(1.0_wp-smoothstep(rmin,rmax,aaa(1)))

!    check ZF driving turbulence by changing ZF profile
!     eta_grid(i)=exp(-((psi_grid(i)-0.38)/0.5)**2)

     itran(i)=int(real(mtheta(i),wp)/q+0.5_wp)
     qtinv(i)=real(mtheta(i))/real(itran(i)) !q value for coordinate transformation
 
     gradth=sqrt(zaj(1)**2+raj(1)**2)/abs(zaj(1)*rqj(1)-zqj(1)*raj(1))
     nunstable1=0.2/(q*gradth*gyroradius)
     munstable1=nunstable1*q
     
     call general_flux_average(aaa(1),object_10,gradth)  ! gradth=<|grad_theta|>
     nunstable=0.2_wp/(q*gradth*gyroradius)
     munstable=nunstable*q

     if(mype==0) write(0,*)'i=',i,'mtheta(i)=',mtheta(i),'q=',q,'qtinv(i)=',&
                  qtinv(i),'m=',munstable,munstable1,'n=',nunstable,nunstable1

     if(mype==0.and.i==ipsi) write(0,*) 'gradth.rho_i=', gradth*gyroradius

     qtinv(i)=1.0_wp/qtinv(i) !inverse q to avoid divide operation
     itran(i)=itran(i)-mtheta(i)*(itran(i)/mtheta(i))

!    several surface averages
     call general_flux_average(aaa(1),object_1,gaa(i))   ! gaa=<|grad_a|^2>
     call general_flux_average(aaa(1),object_5,B2_s(i))  ! B2_s=<B^2>

!    compare to direct surface average
     pdum=aaa(1)
     error_8=1.0e-4_wp
     ifail=1
     Fj(9)=D01AHF(0.0_wp,twopi,error_8,minpts,acc,object_1,-1,ifail)  !1D integral on theta

     error_8=1.0e-4_wp
     ifail=1
     Fj(10)=D01AHF(0.0_wp,twopi,error_8,minpts,acc,object_5,-1,ifail)  !1D integral on theta
     Fj(5)=dVda(i)/(2.0_wp*pi)   ! zeta is ignorable for symmetry

     if(mype==0) write(0,*) '<><><>',i, gaa(i), Fj(9)/Fj(5), '   ', B2_s(i), Fj(10)/Fj(5), '   ', dVda(i) 

     if(mype==0) write(20,204) itran(i), q, dVda(i), psi_grid(i), dlnT, dlnn, dlnT_e, eta_grid(i), s_grid(i), & 
                 temperature(aaa(1),the(1)), temperature_e(aaa(1),the(1)), density(aaa(1),the(1)), zdum, qtinv(i), gaa(i), B2_s(i)

204  format(2x, i8, 20(2x, ES13.6))

  enddo

  call flush(20)

  mthetamax=k

if(mype==0) then
  write(0,*)'mthetamax=', mthetamax, ' mtheta(mpsi)=', mtheta(mpsi)
  write(0,*)'compare grid size with local gyroradius:'
  do i=1,mpsi-1
     tdum=deltat(i)*real(2,wp)
     pdum=psi_grid(i)

!     error_8=1.0d-6
!     ifail=1
!     gradth=D01AHF(0.0d0,tdum,error_8,minpts,acc,arclength,-1,ifail)  !1D integral of NAG

! replaced by
     error=1.0e-3_wp
     intgral=4
     circle=0.0_wp
     delt=tdum/real(intgral,wp)
     aaa(1)=pdum
     aaa(2)=pdum
     do j=1, intgral
        the(1)=delt*real(j-1,wp) 
        the(2)=delt*real(j,wp)
        call esigetrzb(rj,zj,Bj,aaa,the,2)
        circle=circle+sqrt((rj(2)-rj(1))**2+(zj(2)-zj(1))**2)
     enddo
     bcircle=0.0_wp
     i_itera=0
     do while (abs(circle-bcircle)/circle .gt. error)
        bcircle=circle
        intgral=2*intgral
        circle=0.0_wp
        delt=tdum/real(intgral,wp)
        do j=1, intgral
           the(1)=delt*real(j-1,wp)
           the(2)=delt*real(j,wp)
           call esigetrzb(rj,zj,Bj,aaa,the,2)
           circle=circle+sqrt((rj(2)-rj(1))**2+(zj(2)-zj(1))**2)
        enddo

        i_itera=i_itera+1
        if(i_itera.gt.12) then
          write(0,*)'check arclength calculation'
          call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
        end if
     enddo
     gradth=circle

     tdum=0.5_wp*gradth/(vth_grid(i)/B_0)                 ! =dl_th/rho_i
     pdum=0.5_wp*(r_grid(i+1)-r_grid(i-1))/(vth_grid(i)/B_0)   ! =dr/rho_i

     write(0,*)'i=',i, '  dr/rho_i=', pdum, '  dl_th/rho_i=', tdum, mtheta(i), plenth(i)
  enddo

  aaa(1)=1.0_wp
  aaa(2)=1.0_wp
  the(1)=0.0_wp
  the(2)=pi
  call esigetrzb(rj,zj,Bj,aaa,the,2)
  r=abs(rj(1)-rj(2))/2.0_wp       ! minor radius
  a_minor=r 
  write(0,*) 'Average a/rho_i=', r/gyroradius, '   Average aspect ratio=', R_0/r

  r=(rj(1)+rj(2))/2.0_wp
 
  write(501,*) 'a-minor=', a_minor, 'R-major=',r

! get q, R_0/L_T etc profile information
  k=200
  do i=0,k
     aaa(1)=max(1.0e-7_wp,min(1.0_wp,1.0_wp/real(k,wp)*real(i,wp)))
     the(1)=0.0_wp
     call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
     q=-gFaj(1)/gYaj(1)
     grada=sqrt(rqj(1)**2+zqj(1)**2)/abs(zaj(1)*rqj(1)-zqj(1)*raj(1))  !grad(a) at the=0
! calculate R_0/L_T and R_0/L_n at middle plane:
     dlnT=abs(dT_da(aaa(1),the(1))/temperature(aaa(1),the(1)))*grada*r !=R_0*grad(lnT)
     dlnn=abs(dn_da(aaa(1),the(1))/density(aaa(1),the(1)))*grada*r
     rdum=domega_da(aaa(1),the(1),zdum)
     dlnu=abs(rdum/zdum)*grada*r                                       !=R_0*grad[ln(omega_t)]
     write(501,5) aaa(1), q, dlnT, dlnn, dlnT/dlnn, & ! dlnT/dlnn=eta_i
                  temperature(aaa(1),the(1))*Tem_00, dlnT/grada*(rj(1)/r),&
                  dlnT/grada/a_minor, dlnu
     call flush(501)
  end do

end if ! if(mype==0)


! un-comment the next two lines to use magnetic coordinate
!  qtinv=0.0
!  itran=0
  mtdiag=max(1,int(real(mtheta(ipsi))/real(2*mzetamax)+0.5_wp))*2*mzetamax
  if(mype==0) write(0,*) 'mtdiag=', mtdiag

! **********************************************
!        mtdiag=1536
        mtdiag=1440
        mtdiag=768*2 ! special FFT requirment consistent with # of tor grid (which is 64)
!        mtdiag=768
! **********************************************

  if(nonlinear>0.5) mtdiag=mtdiag/2 

  if(mype==0) write(0,*) 'mtdiag=', mtdiag

! starting point for a poloidal grid
  igrid(0)=1
  do i=1,mpsi
     igrid(i)=igrid(i-1)+mtheta(i-1)+1 !igrid(i) and igrid(i)+mtheta(i) are identical point
  enddo

! number of grids on a poloidal plane
  mgrid=sum(mtheta+1)
  mi=micell*(mgrid-mpsi-1)/npartdom   !# of particle per processor
! NOTE: misum is now a floating point variable. A 4-byte integer cannot hold
! a value larger than about 2 billion and we always use misum as a real anyway.
  misum=real(mi,wp)*real(numberpe,wp)    !# of total particles 
  me=mi                            !# of electrons per processor
  mesum=misum                      !# of total electrons

!  write(0,*)mype,' micell=',micell,'  mgrid=',mgrid,'  mpsi=',mpsi,'  npartdom=',npartdom,'  mi=',mi

  mimax=mi+50*ceiling(sqrt(real(mi,wp))) !particle array upper bound
  memax=me+50*ceiling(sqrt(real(me,wp))) !particle array upper bound

! temperature and density on the grid, T_i=n_0=1 at mid-radius
  rtemi=1.0_wp
  rden=1.0_wp

! allocate memory
  allocate (pgyro(4,mgrid),tgyro(4,mgrid),markeri(mgrid),densityi(0:1,mgrid),&
       phi(0:1,mgrid),jtp1(2,mgrid),jtp2(2,mgrid),&
       wtp1(2,mgrid),wtp2(2,mgrid),dtemper(mgrid),heatflux(mgrid),&
       phi_tmp(mgrid),gij(4,mgrid),Jcbn(mgrid),pflux_i(3,0:1,mgrid),&
       pflux_e(3,0:1,mgrid),dphidt(0:1,mgrid))

! Added for flow visualization, Feb. 2011
  allocate (on_grid(6,mgrid))

! # of marker per grid, 
  pmarki=0.0_wp

  jaco_sum=0.0_wp

  allocate(jaco4spectrum(0:mpsi,nfft))

!   !$omp parallel do private(i,j,k,ij,zdum,tdum,dpsi)
  do i=0, mpsi
     pdum=psi_grid(i)
     if(i==0) then
          dpsi=(psi_grid(1)-psi_grid(0))*0.5_wp         ! half cell
        else if(i==mpsi) then
          dpsi=(psi_grid(mpsi)-psi_grid(mpsi-1))*0.5_wp ! half cell
        else
          dpsi=(psi_grid(i+1)-psi_grid(i-1))*0.5_wp     ! full cell
     endif

     dadtdz(i)=dpsi*deltat(i)*deltaz

     do j=1,mtheta(i)
        ij=igrid(i)+j
        zdum=zetamin+deltaz
        tdum=real(j,wp)*deltat(i)+zdum*qtinv(i)

        dtemper(ij)=Jacobian(pdum,tdum)  ! dtemper=Jacobian on grid = 1/[grad(a) . grad(theta) x grad(phi)]
        Jcbn(ij)=dtemper(ij)

!    markeri=the volume of a cell x density
!        markeri(ij)=dtemper(ij)*dpsi*deltat(i)  ! note dzeta is even for each grid point
        markeri(ij)=dtemper(ij)*density(pdum,tdum)*dpsi*deltat(i)  ! nouniform loading
        pmarki(i)=pmarki(i)+markeri(ij)

        ip=metric_elements(pdum,tdum,zdum)
        gij(1,ij)=metric(1)         
        gij(2,ij)=metric(2)         
        gij(3,ij)=metric(3)         
        gij(4,ij)=metric(4)         

        dtemper(ij)=dtemper(ij)*deltaz*deltat(i)   ! =J*dtheta*dzeta
        jaco_sum(i)=jaco_sum(i)+dtemper(ij)
     enddo

     do j=1, nfft
        tdum=twopi/real(nfft,wp)*real(j,wp)
        jaco4spectrum(i,j)=Jacobian(pdum,tdum)
     enddo

  enddo

  if(mype==0) then
     open(unit=733,file='Jaco4spectrum.dat',status='replace',RECL=8192)
     do j=1,nfft
        write(733,802) (jaco4spectrum(i,j),i=0,mpsi)
        call flush(733)
     end do
     close(unit=733)
802  format(503(1x,d12.6)) ! for maxium mpsi=500
  end if

  eachv=sum(pmarki)

! use allnodes for storage of sum(pmarki) over nodes
  call MPI_ALLREDUCE(pmarki,allnodes,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)

  call MPI_ALLREDUCE(jaco_sum,fluxpsi,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
  jaco_sum=fluxpsi ! =int(J*dtheta*dzeta); fluxpsi is released after a use

  fluxpsi=1.0_wp    ! initialize to maintain global charge neutrality in terms of turbulence intensity I(r)
  En_tur=1.0_wp    ! initialize to maintain global charge neutrality in terms of turbulence intensity I(r)
  
  do i=0, mpsi
!$omp parallel do private(j,k,ij)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        markeri(ij)=markeri(ij)/eachv
     enddo
!     pmarki is the inverse of total particle # (sum over nodes) at an annulus i
!     pmarki(i)=1.0/(real(numberpe)*real(mi)*(pmarki(i)/eachv))
     markeri(igrid(i))=markeri(igrid(i)+mtheta(i)) ! markeri(:,igrid(i)) is missing in last do-loop
     dtemper(igrid(i))=dtemper(igrid(i)+mtheta(i))
     Jcbn(igrid(i))=Jcbn(igrid(i)+mtheta(i))
     gij(:,igrid(i))=gij(:,igrid(i)+mtheta(i))
  enddo
  
  totalv=sum(allnodes) 

! pmarki is the inverse of total particle # (sum over nodes) at an annulus i
  do i=0, mpsi
     pmarki(i)=allnodes(i)/totalv
  enddo

if(mype==0) then 
   do i=0, mpsi
       write(21,204) i, pmarki(i)
       do j=1,mtheta(i)
          ij=igrid(i)+j
          write(22,204) ij, markeri(ij)
       end do
   end do
end if
call flush(21)
call flush(22)

! Using NAG to calculate 3D integral
  do i=1, nrbin
     low_lim(1)=a0+(a1-a0)/real(nrbin,wp)*real(i-1,wp)
     up_lim(1)=a0+(a1-a0)/real(nrbin,wp)*real(i,wp)
     low_lim(2)=0.0_wp
     up_lim(2)=2.0_wp*pi
     low_lim(3)=zetamin
     up_lim(3)=zetamax
     error_8=1.0e-5_wp
     minpts=0
     ifail=1
     call d01fcf(3,low_lim,up_lim,minpts,maxpts,FUN1, &
                 error_8,acc,lenwrk,wrkstr,finval,ifail)
     rtden(i)=finval     ! rtden(i)=particles of annulus i on my node

     r=(low_lim(1)+up_lim(1))/2.0_wp
     tem_nrbin(i)=temperature(r,0.0_wp)   ! Ti at the center of anulus
     tem_nrbin(i)=1.0_wp/tem_nrbin(i)
  enddo

  call MPI_ALLREDUCE(rtden,rtden_allnodes,nrbin,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
  rtden_global=sum(rtden_allnodes)      ! total volume of simulation domain
  do i=1, nrbin
     rtden(i)=misum*rtden_allnodes(i)/rtden_global
     rtden(i)=1.0_wp/rtden(i)
  enddo

  low_lim(1)=a0
  up_lim(1)=a1
  low_lim(2)=0.0_wp
  up_lim(2)=2.0_wp*pi
  low_lim(3)=0.0_wp
  up_lim(3)=2.0_wp*pi
  error_8=1.0e-5_wp
  minpts=0
  ifail=1
  call d01fcf(3,low_lim,up_lim,minpts,maxpts,FUN1, &
              error_8,acc,lenwrk,wrkstr,finval,ifail)    ! finval=total particle number

  cloud_i=finval
  cloud_e=finval

  error_8=1.0e-5_wp
  minpts=0
  ifail=1
  call d01fcf(3,low_lim,up_lim,minpts,maxpts,f_jaco, &
              error_8,acc,lenwrk,wrkstr,finval_v,ifail)  ! finval_v=volume of simulation torus
  density_ave=finval/finval_v     ! average density
  
if(mype==0) write(0,*) 'average density=', density_ave

! Using NAG to calculate 3D integral
  do i=1, nbin_a
     low_lim(1)=a0+(a1-a0)*real(i-1,wp)/real(nbin_a,wp)
     up_lim(1)=a0+(a1-a0)*real(i,wp)/real(nbin_a,wp)
     low_lim(2)=0.0_wp
     up_lim(2)=2.0_wp*pi
     low_lim(3)=0.0_wp
     up_lim(3)=2.0_wp*pi
     error_8=1.0e-5_wp
     minpts=0
     ifail=1
     call d01fcf(3,low_lim,up_lim,minpts,maxpts,FUN1,error_8,acc,lenwrk,wrkstr,&
                 sn0(i),ifail)    ! sn0(i)=particles of annulus i 

     error_8=1.0e-5_wp
     minpts=0
     ifail=1
     call d01fcf(3,low_lim,up_lim,minpts,maxpts,f_jaco,error_8,acc,lenwrk,&
             wrkstr,dnty_on_bin(i),ifail) ! dnty_on_bin(i)=volume of annulus i

     dnty_on_bin(i)=sn0(i)/dnty_on_bin(i) ! average density of annulus i
     sn0(i)=sn0(i)/finval

     r=(low_lim(1)+up_lim(1))/2.0_wp  ! center of cell
     Tr(i)=temperature(r,0.0_wp)
     dnty_r(i)=density(r,0.0_wp)

     Tr_e(i)=temperature_e(r,0.0_wp)
     Tr_e(i)=Tr_e(i)*msm0     ! Te with the same unit with Ti

     dlnT=dT_da(r,0.0)
     dlnT_e=dT_da_e(r,0.0)
     dlnn=dn_da(r,0.0)
     call general_flux_average(r,object_9,dT_10(i))      ! =<|grad_a|>
     dlnT_e=dT_10(i)*dlnT_e
     dlnT_e=dlnT_e*msm0
     dlnn=dT_10(i)*dlnn
     dT_10(i)=dT_10(i)*dlnT                              ! dT_10(i)==<|grad_a|>*dT/da
     rdum=domega_da(r,0.0,sint)                          ! rdum=d_omega_t/da
     psimax=rdum
     call general_flux_average(r,object_11,dO_10(i))     ! dO_10(i)=<R^2*|grad_a|>
     dO_10(i)=dO_10(i)*rdum                              ! =<R^2*|grad_a|>*d_omega_t/da

!  use 2D integral (due to axisymmetry) for cross check

     low1(1)=low_lim(1)
     up1(1)=up_lim(1)
     low1(2)=low_lim(2)
     up1(2)=up_lim(2)

     error_8=1.0e-6_wp       !  eps=1.0d-7
     IFAIL=1
     minpts=0
     call D01FCF(2,low1,up1,minpts,maxpts,FUN1,error_8,acc, &
                  lenwrk,wrkstr,N_on_ring,IFAIL)           !2D integral of NAG
     error_8=1.0e-6_wp       !  eps=1.0d-7
     IFAIL=1
     minpts=0
     call D01FCF(2,low1,up1,minpts,maxpts,FUN_J,error_8,acc, &
                  lenwrk,wrkstr,Vol_ring,IFAIL)           !2D integral of NAG



     if(mype==0) then
        call general_flux_average(r,object_9s,inv_grada) ! inv_grada=<|1/grad_a|>
        call general_flux_average(r,object_1,grada)      ! grada=<|grad_a|^2>
        dlnT=dT_da(r,0.0_wp)
        rdum=-dlnT*grada*dnty_on_bin(i)*gyroradius**2*vthc/a_minor
        zdum=inv_grada*uveloci*E_00*1.0e-3_wp
        write(9,200) r, rdum, zdum
        write(19,200) r, -1.0/(dnty_r(i)*dO_10(i))/utime,                       &  ! x(this factor) -> chi_phi in units m^2/s 
                         -1.0/(dnty_r(i)*dT_10(i))/(uveloci*E_00*1.0d-3)/utime, &  ! x(this factor) -> chi_i in units m^2/s
                         -1.0/(dnty_r(i)*dlnT_e)/(uveloci*E_00*1.0d-3)/utime,   &  ! x(this factor) -> chi_e in units m^2/s
                         -1.0/(dlnn*1.0d6),                                     &  ! x(this factor) -> D in units m^2/s
                         psimax,sint
200     format(16(2x,ES13.6))
        call flush(9)

        write(0,*) 'DDDSSS', dnty_on_bin(i), N_on_ring/Vol_ring
     end if
  enddo

!  if(mype==0)  call theory

     r=rc
     tdum=temperature(r,0.0_wp)
     psimax=density(r,0.0_wp)

     dlnT=dT_da(r,0.0)
     dlnT_e=dT_da_e(r,0.0)
     dlnn=dn_da(r,0.0)
     call general_flux_average(r,object_9,zdum)      ! =<|grad_a|>
     dlnT_e=zdum*dlnT_e
     dlnT_e=dlnT_e*msm0
     dlnn=zdum*dlnn
     zdum=zdum*dlnT                              ! dT_10(i)==<|grad_a|>*dT/da
     rdum=domega_da(r,0.0,sint)                          ! rdum=d_omega_t/da
     call general_flux_average(r,object_11,q)     ! dO_10(i)=<R^2*|grad_a|>
     q=q*rdum                              ! =<R^2*|grad_a|>*d_omega_t/da

     if(mype==0) then
        write(19,200) r, -1.0/(psimax*q)/utime,                       &  ! x(this factor) -> chi_phi in units m^2/s
                         -1.0/(psimax*zdum)/(uveloci*E_00*1.0d-3)/utime, &  ! x(this factor) -> chi_i in units m^2/s
                         -1.0/(psimax*dlnT_e)/(uveloci*E_00*1.0d-3)/utime,   &  ! x(this factor) -> chi_e in units m^2/s
                         -1.0/(dlnn*1.0d6),     &                             ! x(this factor) -> D in units m^2/s
                         rdum,sint
     end if
     call flush_(19)
     if(mype==0) close(19)




  do i=0, mpsi
     j=max(0,i-1)
     k=min(mpsi,i+1)
     low_lim(1)=(psi_grid(j)+psi_grid(i))*0.5_wp
     up_lim(1)=(psi_grid(i)+psi_grid(k))*0.5_wp
     low_lim(2)=0.0_wp
     up_lim(2)=2.0_wp*pi
     low_lim(3)=0.0_wp
     up_lim(3)=2.0_wp*pi
     error_8=1.0e-5_wp
     minpts=0
     ifail=1
     call d01fcf(3,low_lim,up_lim,minpts,maxpts,FUN1,error_8,acc,lenwrk,wrkstr,&
                 sn0_g(i),ifail)    ! sn0_g(i)=particles of anulus i

     error_8=1.0e-5_wp
     minpts=0
     ifail=1
     call d01fcf(3,low_lim,up_lim,minpts,maxpts,f_jaco,error_8,acc,lenwrk,&
           wrkstr,dnty_on_bin_g(i),ifail) ! dnty_on_bin_g(i)=volume of annulus i

     dnty_on_bin_g(i)=sn0_g(i)/dnty_on_bin_g(i) ! average density of annulus i
     sn0_g(i)=sn0_g(i)/finval

if(mype==0) write(0,*) 'density_bin=', dnty_on_bin_g(i), sn0_g(i)

  enddo

! diagnostic box:  r:     [rprobe-dradius, rprobe+dradius]
!                  theta: [tprobe-dtheta, tprobe+dtheta]
!                  zeta:  [0, 2pi]
  do i=1,nprobe
     low_lim(1)=rprobe(i)-dradius
     up_lim(1)=rprobe(i)+dradius
     low_lim(2)=tprobe(i)-dtheta
     up_lim(2)=tprobe(i)+dtheta
     low_lim(3)=0.0_wp
     up_lim(3)=2.0_wp*pi
     error_8=1.0e-4_wp
     minpts=0
     ifail=1
     call d01fcf(3,low_lim,up_lim,minpts,maxpts,FUN1,error_8,acc,lenwrk,wrkstr,&
                 sn_box(i),ifail)    ! sn_box(i)=particles in diagnostic box i

     error_8=1.0e-4_wp
     minpts=0
     ifail=1
     call d01fcf(3,low_lim,up_lim,minpts,maxpts,f_jaco,error_8,acc,lenwrk,&
           wrkstr,density_in_box(i),ifail) ! density_in_box(i)=volume of box i

     density_in_box(i)=sn_box(i)/density_in_box(i)   ! average density of annulus i
     sn_box(i)=sn_box(i)/finval

     Ti_box(i)=temperature(rprobe(i),tprobe(i))

     Te_box(i)=temperature_e(rprobe(i),tprobe(i))
     Te_box(i)=Te_box(i)*msm0     ! Te with the same unit with Ti
  enddo

! 4-point gyro-averaging for sqrt(mu)=gyroradius on grid of magnetic coordinates
! rhoi=(rho_th*sqrt(2/b))*(sqrt(mu)/rho_th), rho_th=gyroradius
! dtheta/delta_x=1/(r*(1+r*cos(theta))), delta_x=poloidal length increasement
!!!!        !$omp parallel do private(i,j,ij,zdum,pdum,tdum,q,b,dtheta_dx,rhoi)
  do i=0,mpsi
     zdum=zetamin+deltaz     ! zeta local to me 
     !   zdum=0.0_wp         ! old before 2013
     pdum=psi_grid(i)

     do j=0,mtheta(i)
        ij=igrid(i)+j
        tdum=deltat(i)*real(j,wp)+zdum*qtinv(i)     ! local theta of grid point
        tdum=tdum/(2.0_wp*pi)
        tdum=2.0_wp*pi*(tdum-aint(tdum))

        aaa(1)=pdum
        the(1)=tdum
        call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
        q=-gFaj(1)/gYaj(1)
        inv_jaco=1.0_wp/jacobian(pdum,tdum)

!          in metric, 1=grad_a^2,2=grad_theta^2,3=grad_zeta^2
!          4=grad_a.grad_theta,5=grad_a.grad_zeta,6=grad_theta.grad_zeta

        k=metric_elements(pdum,tdum,zdum) 

!        dtheta_dx=1.0/(r*(1.0+r*cos(tdum)))

!        dtheta_dx=1.0/sqrt(metric(3)*Bj(1)**2/gFaj(1)**2-inv_jaco**2) 

!        for some ESC equilibrium, (grad_phiXB)^2=metric(3)*Bj(1)**2/gFaj(1)**2-inv_jaco**2<0 
!        dtheta_dx=1.0/sqrt(abs(metric(3)*Bj(1)**2/gFaj(1)**2-inv_jaco**2))

! first two points perpendicular to field line on poloidal surface            
! first two points (1 and 2) are in grad_zetaxB direction with rho_i unchanged 

        rhoi=sqrt(2.0_wp/Bj(1))     ! rho=rhoi*sqrt(mu)

!        pgyro(1,ij)=-rhoi
!        pgyro(2,ij)=rhoi

!        pgyro(1,ij)=-((metric(5)*metric(5)-metric(1)*metric(3))/q+  &
!                      (metric(1)*metric(6)-metric(5)*metric(4)))*dtheta_dx*rhoi   !wx gg

        pgyro(1,ij)=-sqrt(metric(1))*rhoi           ! for tokamak; da/dl
        rdum=dsda(pdum)                             ! = dsda
        pgyro(1,ij)=pgyro(1,ij)*rdum                ! ds/dl=(ds/da)*(da/dl)

        pgyro(2,ij)=-pgyro(1,ij)       !wx gg

! non-orthorgonality between psi and theta: tgyro=-rhoi*dtheta_dx*r*sin(tdum)
!        tgyro(1,ij)=rhoi*dtheta_dx*r*sin(tdum)
!        tgyro(2,ij)=-rhoi*dtheta_dx*r*sin(tdum)

!        tgyro(1,ij)=-((metric(5)*metric(6)-metric(3)*metric(4))/q+  &
!                      (metric(4)*metric(6)-metric(2)*metric(5)))*dtheta_dx*rhoi   !wx gg

        tgyro(1,ij)=-metric(4)/sqrt(metric(1))*rhoi   ! for tokamak

        tgyro(2,ij)=-tgyro(1,ij)

! the other two points (3 and 4) are in (grad_zetaxB)xgrad_zeta direction 
! with rhoi enlarged by a factor JB*|grad_zeta| 
!        pgyro(3,ij)=-metric(5)*dtheta_dx*Bj(1)*rhoi/gFaj(1)
        pgyro(3,ij)=-metric(5)*q/sqrt(metric(1))/metric(3)*Bj(1)*rhoi/gFaj(1)
        pgyro(3,ij)=pgyro(3,ij)*rdum             ! ds/dl=(ds/da)*(da/dl)
!        pgyro(3,ij)=pgyro(3,ij)*vthc/vth_grid(i) ! ds/dl=(ds/da)*(da/dl)

        pgyro(4,ij)=-pgyro(3,ij)

!        tgyro(3,ij)=(metric(3)/q-metric(6))*dtheta_dx*Bj(1)*rhoi/gFaj(1)
        tgyro(3,ij)=Bj(1)*rhoi/gFaj(1)/sqrt(metric(1))

        tgyro(4,ij)=-tgyro(3,ij)
     enddo
! to be safe do the following
     pgyro(:,igrid(i))=pgyro(:,igrid(i)+mtheta(i)) !(:,igrid(i)) is missing in last do-loop
     tgyro(:,igrid(i))=tgyro(:,igrid(i)+mtheta(i)) !(:,igrid(i)) is missing in last do-loop
  enddo

!  for zero gyroradius limit, set tgyro,pgyro=0 !wx
!  pgyro=0.0
!  tgyro=0.0

! open file for standard output, record program starting time
  if(mype == 0)then
     if(stdout /= 6 .and. stdout /= 0)open(stdout,file='stdout.out',status='replace')
     call date_and_time(date,time)
     write(stdout,*) 'Program starts at DATE=', date, 'TIME=', time
     delr=(r_grid(ipsi)-r_grid(ipsi-1))/gyroradius     ! delr=deltar/gyroradius    ! ???
     delt=deltat(ipsi)*r_grid(ipsi)/gyroradius         
                         ! delt=deltat(mpsi/2)*(a0+deltar*real(mpsi/2))/gyroradius     ! ??? 
     mid_theta=mtheta(ipsi)
     write(stdout,run_parameters)

#ifdef _OPENMP
!$omp parallel private(nthreads)
     nthreads=omp_get_num_threads()  !Get the number of threads if using OMP
!$omp single
     write(stdout,'(/,"===================================")')
     write(stdout,*)' Number of OpenMP threads = ',nthreads
     write(stdout,'("===================================",/)')
!$omp end single nowait
!$omp end parallel
#else
     write(stdout,'(/,"===================================")')
     write(stdout,*)' Run without OpenMP threads'
     write(stdout,'("===================================",/)')
#endif
  endif

! initiate radial interpolation for grid
  zdum=zetamin+deltaz
!$omp parallel do private(i,ip,j,indp,indt,ij,tdum,jt,wt)
!  do i=1,mpsi-1     ! changed to do i=0,mpsi for radial smoothing
  do i=0,mpsi        ! incoporating boundary condition dphi/da=0
     do ip=1,2
        indp=min(mpsi,i+ip)
        indt=max(0,i-ip)
        do j=1,mtheta(i)
           ij=igrid(i)+j
! upward
           tdum=(real(j,wp)*deltat(i)+zdum*(qtinv(i)-qtinv(indp)))/deltat(indp)
           jt=floor(tdum)
           wt=tdum-real(jt,wp)
           jt=mod(jt+mtheta(indp),mtheta(indp))
           if(ip==1)then
              wtp1(1,ij)=wt
              jtp1(1,ij)=igrid(indp)+jt
           else
              wtp2(1,ij)=wt
              jtp2(1,ij)=igrid(indp)+jt
           endif
! downward
               
           tdum=(real(j,wp)*deltat(i)+zdum*(qtinv(i)-qtinv(indt)))/deltat(indt)
           jt=floor(tdum)
           wt=tdum-real(jt,wp)
           jt=mod(jt+mtheta(indt),mtheta(indt))
           if(ip==1)then
              wtp1(2,ij)=wt
              jtp1(2,ij)=igrid(indt)+jt
           else
              wtp2(2,ij)=wt
              jtp2(2,ij)=igrid(indt)+jt
           endif
        enddo
     enddo
  enddo
! note added: w(j)tp1(2)(:,igrid(i),:) has no valid data

! Import from GTC-Neo or estimate equilibrium Er. If from GTC_Neo, the file
! name is stored in the character variable "Er_data_name", which can be found
! in the global_parameters module 
 ! call Er_provider(Er_from_GTC_Neo,n_node) 
  call Er_provider

!  call poisson_initial

! allocate memory
  allocate (F(mimax),Fa(mimax),gFa(mimax),gFaa(mimax),gYa(mimax),gYaa(mimax),TT(mimax), &
            Ta(mimax),PP(mimax),Pa(mimax), rr(mimax),ra(mimax),rq(mimax),z(mimax), &
            za(mimax),zq(mimax),BB(mimax),Ba(mimax),Bq(mimax),gh(mimax),gha(mimax), &
            ghq(mimax),sw(mimax),backn(mimax),backu(mimax),backT(mimax),&
            dlnnda(mimax),dlnTda(mimax),backne(mimax),backue(mimax),&
            backTe(mimax),dlnndae(mimax),dlnTdae(mimax),duda(mimax))

! c1.     Making link with C-ESI
!  1. Supplying addresses to C-code (not a, gq)
      call esilink2c(F,Fa,gFa,gFaa,gYa,gYaa,TT,Ta,PP,Pa, &
                   rr,ra,rq,z,za,zq,BB,Ba,Bq,gh,gha,ghq,sw)

  call get_table_4_dlnAda    ! get ready for dln<A>/da and make A table;
                             ! where A=exp(u^2/2T). NOTE: must be put after
                             ! call to esilink2c to have F, Fa, ... defined





! *******************************************************************

! Post-processor calculation of zonal flow ExB shearing rate 
! need a separate directory with ZF data copied 
!  i=30265    ! number of lines in ZF data files
!  call get_ZF_shear(i)
!  call mpi_finalize(ierror)
!  STOP

! *******************************************************************





! allocate memory
  allocate (zion(7,mimax),zion0(7,mimax),jtion0(4,mimax),jtion1(4,mimax),&
       wzion(mimax),wpion(4,mimax),wtion0(4,mimax),wtion1(4,mimax),enion(mimax))
  if(ifelectron/=0) then
     allocate (zeon(7,memax),zeon0(7,memax),zeon00(7,memax),jteon0(memax),&
               jteon1(memax),wzeon(memax),wpeon(memax),&
               wteon0(memax),wteon1(memax),markere(mgrid),&
             !   densitye(0:1,mgrid),cij(3,mgrid))
               densitye(0:1,mgrid),cij(5,mgrid))

if(mype==-100) then  ! before 2012

     dela=1.0_wp/Ydr
     do i=0, mpsi
        pdum=psi_grid(i)

        r=pdum-a0
        r=r*dela
        ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
        fip=r-real(ip,wp)

        vthe=temperature_e(pdum,0.0_wp)
        vthe=vthe*msm0     ! Te with the same unit with Ti

        do j=1,mtheta(i)
           ij=igrid(i)+j
           zdum=zetamin+deltaz
           tdum=real(j,wp)*deltat(i)+zdum*qtinv(i)
           cij(1,ij)=pdum        ! when input --> a
           cij(2,ij)=tdum        ! when input --> theta
           cij(3,ij)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! --> dphi_0/da at grid points
           cij(3,ij)=cij(3,ij)/vthe     ! when input --> (dphi_0/da)/Te at grid points
        enddo

        !       igrid(i)=igrid(i)+mtheta(i)
        cij(1,igrid(i))=cij(1,igrid(i)+mtheta(i))
        cij(2,igrid(i))=cij(2,igrid(i)+mtheta(i))
        cij(3,igrid(i))=cij(3,igrid(i)+mtheta(i))
     enddo

     !          a      theta    Er_grid     n        as input
     !          |        |        |         |        
     i=get4exb(cij(1,:),cij(2,:),cij(3,:),mgrid)
     !          |        |        |         |        
     !         cij1     cij2     cij3       n        as output

end if

     do i=0, mpsi
        pdum=psi_grid(i)
        do j=1,mtheta(i)
           ij=igrid(i)+j
           zdum=zetamin+deltaz
           tdum=real(j,wp)*deltat(i)+zdum*qtinv(i)
           cij(1,ij)=pdum        ! when input --> a
           cij(2,ij)=tdum        ! when input --> theta
        enddo

        !       igrid(i)=igrid(i)+mtheta(i)
        cij(1,igrid(i))=cij(1,igrid(i)+mtheta(i))
        cij(2,igrid(i))=cij(2,igrid(i)+mtheta(i))
     enddo

            !     a      theta     null      null     null      n        as input
            !     |        |        |         |        |        |
     i=get4exb12(cij(1,:),cij(2,:),cij(3,:),cij(4,:),cij(5,:),mgrid)
            !     |        |        |         |        |
            !    cij1     cij2     cij3      cij4     cij5               as output

  endif

! calculate grid point location in Cartisian Coordinates
if(irun==0) then
  allocate(coord_X(mgrid),coord_Y(mgrid),coord_Z(mgrid),&
           coord_R(mgrid),coord_Phi(mgrid))
  do i=0,mpsi
     pdum=psi_grid(i)

     do j=0,mtheta(i)
        ij=igrid(i)+j
        zdum=zetamin+deltaz
        tdum=real(j,wp)*deltat(i)+zdum*qtinv(i)

        aaa(1)=pdum
        the(1)=tdum
        call esigetrzb(rj,zj,Bj,aaa,the,1)

        coord_X(ij)=rj(1)*cos(zdum)   ! --> x
        coord_Y(ij)=rj(1)*sin(zdum)   ! --> y
        coord_Z(ij)=zj(1)             ! --> z
        coord_R(ij)=rj(1)             ! --> R
        coord_Phi(ij)=zdum            ! --> phi
     end do
     !evector(:,1,igrid(i))=evector(:,1,igrid(i)+mtheta(i))
  end do

  if(myrank_partd==0) then   ! myrank_partd==1 will write x,y,z
      write(cdum1,'("GRID_COORD.",i5.5)')myrank_toroidal
      open(unit=995,file=cdum1,form='unformatted')
      write(995)coord_X
      write(995)coord_Y
      write(995)coord_Z
      write(995)coord_R
      write(995)coord_Phi
      close(995)
  endif

  deallocate(coord_X,coord_Y,coord_Z,coord_R,coord_Phi)
end if

! Now that we have freed  up some memory by deallocating the coordinate
! arrays coord_{X,Y,Z,r,phi} (when irun=0) we can allocate evector()
  allocate(evector(3,0:1,mgrid))

! Added for flow visualization, Feb. 2011
  do i=0,mpsi
     pdum=psi_grid(i)

     do j=0,mtheta(i)
        ij=igrid(i)+j
        zdum=zetamin+deltaz
        tdum=real(j,wp)*deltat(i)+zdum*qtinv(i)

        aaa(1)=pdum
        the(1)=tdum
        call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)

        on_grid(1,ij)=raj(1)
        on_grid(2,ij)=rqj(1)
        on_grid(3,ij)=zaj(1)
        on_grid(4,ij)=zqj(1)
        on_grid(5,ij)=rj(1)
        on_grid(6,ij)=abs(zaj(1)*rqj(1)-zqj(1)*raj(1))/sqrt(zaj(1)**2+raj(1)**2) ! =1/|grad(theta)|

     end do
     on_grid(:,igrid(i))=on_grid(:,igrid(i)+mtheta(i))
  end do

     do k=1,nprobe
        pdum=rprobe(k)
        tdum=tprobe(k)
        vthprobe(k)=temperature(pdum,tdum)
        vthprobe(k)=sqrt(2.0d0*vthprobe(k)) 
        uflow(k)=u_parallel(pdum,tdum)
        uflow(k)=uflow(k)/vthprobe(k)    ! normalized by local thermal velocity

        vthprobe_e(k)=temperature_e(pdum,tdum)
        vthprobe_e(k)=sqrt(2.0d0*vthprobe_e(k))

! calculate trap-passing boundary
        aaa(1)=rprobe(k)
        aaa(2)=rprobe(k)
        the(1)=tprobe(k)
        the(2)=pi
        call esigetrzb(rj,zj,Bj,aaa,the,2)
        vb(1)=0.0_wp
        vb(2)=4.0_wp
        vp(k,1)=0.0
        vp(k,2)=vb(2)/sqrt(Bj(2)/Bj(1)-1.0_wp)
     end do
     if(mype==0) then
        write(11,55) vb(1), (vp(k,1),k=1,nprobe)
        write(11,55) vb(2), (vp(k,2),k=1,nprobe)
        call flush(11)
     end if
55   format(30(1x,ES13.6))
     close(11)


!  add info needed for kr spectrum calculation:
! a_ref=a_ref
pdum=temperature(a_ref,0.0_wp)   ! Ti
rdum=temperature_e(a_ref,0.0_wp)
rdum=rdum*msm0                   ! Te with the same unit as Ti
aaa(1)=a_ref
aaa(2)=a_ref
aaa(3)=a_ref
the(1)=0.0_wp
the(2)=0.5_wp*pi
the(3)=pi
call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,3)
grada=sqrt(rqj(1)**2+zqj(1)**2)/abs(zaj(1)*rqj(1)-zqj(1)*raj(1)) !grad(a) at the=0
tdum= sqrt(rqj(2)**2+zqj(2)**2)/abs(zaj(2)*rqj(2)-zqj(2)*raj(2)) !grad(a) at the=pi/2
zdum= sqrt(rqj(3)**2+zqj(3)**2)/abs(zaj(3)*rqj(3)-zqj(3)*raj(3)) !grad(a) at the=pi
! normalized rho_s using Te
rhos1=sqrt(rdum)/Bj(1)    ! normalized rho_s at the=0
rhos2=sqrt(rdum)/Bj(2)    ! normalized rho_s at the=pi/2
rhos3=sqrt(rdum)/Bj(3)    ! normalized rho_s at the=pi

if(mype==0) then
   print*, ' reference point a_ref', a_ref
   print*, ' grid size deltar in terms of s', deltar 
   print*, '|grad(a)| (a_ref,    0)', grada, '  ', rhos1
   print*, '|grad(a)| (a_ref, pi/2)', tdum, '  ', rhos2
   print*, '|grad(a)| (a_ref,   pi)', zdum, '  ', rhos3
end if

!  add info needed for k_theta spectrum calculation:
! a_ref=0.6d0                      ! radial location for  k_theta spectrum calculation
! a_ref=0.5733537d0                ! radial location for  k_theta spectrum calculation
 a_ref=0.59516d0                ! radial location for  k_theta spectrum calculation
pdum=temperature(a_ref,0.0_wp)   ! Ti
rdum=temperature_e(a_ref,0.0_wp)
rdum=rdum*msm0                   ! Te with the same unit as Ti
aaa(1)=a_ref
aaa(2)=a_ref
aaa(3)=a_ref
the(1)=0.0_wp
the(2)=0.5_wp*pi
the(3)=pi
call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,3)
grada=zaj(1)*rqj(1)-zqj(1)*raj(1)
tdum= zaj(2)*rqj(2)-zqj(2)*raj(2)
zdum= zaj(3)*rqj(3)-zqj(3)*raj(3)

grada=sqrt((zaj(1)**2+raj(1)**2)/grada**2) ! |grad(theta)| at the=0
tdum= sqrt((zaj(2)**2+raj(2)**2)/tdum**2)  ! |grad(theta)| at the=pi/2
zdum= sqrt((zaj(3)**2+raj(3)**2)/zdum**2)  ! |grad(theta)| at the=pi
! normalized rho_s using Te
rhos1=sqrt(pdum)/Bj(1)*sqrt(rdum/pdum)  ! normalized rho_s at the=0
rhos2=sqrt(pdum)/Bj(2)*sqrt(rdum/pdum)  ! normalized rho_s at the=pi/2
rhos3=sqrt(pdum)/Bj(3)*sqrt(rdum/pdum)  ! normalized rho_s at the=pi

if(mype==0) then
   print*, '(|(grad(theta)|,    0)', grada, '  ', rhos1
   print*, '(|(grad(theta)|, pi/2)', tdum, '  ', rhos2
   print*, '(|(grad(theta)|,   pi)', zdum, '  ', rhos3
end if

!    generate (R,Z) coordinates for contour plot of n=0 mode
if(mype==0) then
   k=768
   tdum=twopi/real(k)
   open(unit=633,file='RZ_grid.dat')
   do i=0,mpsi
      aaa(1)=psi_grid(i)
      do j=1,k
         the(1)=tdum*real(j)
         call esigetrzb(rj,zj,Bj,aaa,the,1)
         write(633,100) i, j, rj(1), zj(1)
      end do
      the(1)=tdum
      call esigetrzb(rj,zj,Bj,aaa,the,1)  
      write(633,100) i, (k+1), rj(1), zj(1)
   end do
100  format(2x, i4, 2x, i4, 2x, d12.6, 2x, d12.6)
   close(633)
end if


! write down local parameters for flux tube simulations; May 09
  if(mype==0) then
     aaa(1)=a_ref
     the(1)=0.0_wp
     call esiget2dfunctions(aaa,the,1)

     q=-gFa(1)/gYa(1)
!     magnetic_shear=a_ref*(gFaa(1)/gFa(1)-gYaa(1)/gYa(1)) ! s=r/q*(dq/dr)
     magnetic_shear=(gFaa(1)/gFa(1)-gYaa(1)/gYa(1))        ! s=r0/q*(dq/dr); some use r0 instead of r

     vthi=temperature(a_ref,0.0_wp)
     vthe=temperature_e(a_ref,0.0_wp)
     vthe=vthe*msm0     ! Te with the same unit as Ti
     tdum=vthe/vthi     ! Ti/Te

     call general_flux_average(a_ref,object_9,grada)        ! =<|grad_a|>
     call general_flux_average(a_ref,object_9s,inv_grada)   ! inv_grada=<|1/grad_a|>

     dlnT=abs(dT_da(aaa(1),the(1))/temperature(aaa(1),the(1)))*grada*R_0       !=R_0*grad(lnT)
     dlnn=abs(dn_da(aaa(1),the(1))/density(aaa(1),the(1)))*grada*R_0           !=R_0*grad(lnn)
     dlnT_e=abs(dT_da_e(aaa(1),the(1))/temperature_e(aaa(1),the(1)))*grada*R_0 !=R_0*grad(lnT)
     rdum=domega_da(aaa(1),the(1),zdum)
     dlnu=abs(rdum/zdum)*grada*R_0                                             !=R_0*grad[ln(omega_t)]

     write(501,*) 'm_e/m_i=', amasse/amassi
     write(501,*) 'Te/Ti=', tdum
     write(501,*) 'q=', q
     write(501,*) 'magnetic shear s=', magnetic_shear
     write(501,*) 'a/R_0=', a_minor/R_0
     write(501,*) 'R_0/L_Ti=', dlnT, R_0/Lci
     write(501,*) 'R_0/L_Te=', dlnT_e, R_0/Lce
     write(501,*) 'R_0/L_n=', dlnn, R_0/Lcn
     write(501,*) 'R_0/L_omega_t=', dlnu, R_0/Lcu
     write(501,*) 'B_0=', B_0, 'B(r=a_ref, theta=0)=', BB(1), 'B_0/B=', B_0/BB(1)
     write(501,*) 'Ti=', vthi*Tem_00, ' Kev'
     write(501,*) 'Te=', vthe*Tem_00, ' Kev'
     write(501,*) 'ne=', density(aaa(1),the(1)), 'cm^3'
     write(501,*) 'conversion factor for calculated total intrinsic torque to units Nm/m^3', &
                  1.0d10/21.5317*pmass0*R_00*uveloci/((Lc/vthc)*5.5)*0.1/utime  ! for MODEL_KE94 -- see note

     call flush(501)

!  get q profile information
     k=200
     do i=0,k
        aaa(1)=max(1.0e-7_wp,min(1.0_wp,1.0_wp/real(k,wp)*real(i,wp)))
        the(1)=0.0_wp
        call esiget2dfunctions(aaa,the,1)

        q=-gFa(1)/gYa(1)

        dqda=-gFaa(1)/gYa(1)+gFa(1)*gYaa(1)/gYa(1)**2
        magnetic_shear=aaa(1)*(gFaa(1)/gFa(1)-gYaa(1)/gYa(1))       ! s=r/q*(dq/dr)

        grada=sqrt(rq(1)**2+zq(1)**2)/abs(za(1)*rq(1)-zq(1)*ra(1)) !grad(a) at the=0
! calculate R_0/L_T and R_0/L_n at middle plane:
        dlnT=abs(dT_da(aaa(1),the(1))/temperature(aaa(1),the(1)))*grada*R_0 !=R_0*grad(lnT)
        dlnn=abs(dn_da(aaa(1),the(1))/density(aaa(1),the(1)))*grada*R_0     !=R_0*grad(lnn)
        rdum=domega_da(aaa(1),the(1),zdum)
        dlnu=abs(rdum/zdum)*grada*R_0                                       !=R_0*grad[ln(omega_t)]
        write(501,5) aaa(1), q, dqda, magnetic_shear, dlnT, dlnn, dlnT/dlnn, & ! dlnT/dlnn=eta_i
                     temperature(aaa(1),the(1))*Tem_00, dlnT/grada*(rr(1)/R_0), dlnu
        call flush(501)
     end do

     close(501)
  end if


  w2_ion0=1.0_wp
  w2_eon0=1.0_wp


end subroutine setup

!=============================================================================

  Subroutine read_input_params

!=============================================================================

  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface
  use spl_profiles
  use profile_array

  implicit none

  integer,parameter :: n_integers=47, n_reals=58, n_chars=3
  integer :: integer_params(n_integers)
  character(25),dimension(n_chars) :: char_params
! only the first processor will open files...
  real*8       :: real_params(n_reals)
  integer ierror,n(num_mode),m(num_mode)

  namelist/ read1 /  eq_data_name,       &   ! MHD equilibrium file name
                     irun,               &   ! =0, 1st run; >0 continueing run
                     mstep,              &   ! number of last step
                     msnap,              &   ! step interval for snapshot, PHI record and checkpoint
                     ndiag,              &   ! calculate fluxes, tur intensity, etc, every ndiag-step
                     nonlinear,          &   ! =1.0, nonlinear run, =0.0, linear run
                     paranl,             &   ! =1, parallel nonlinearity is on (always on so far) 
                     mode00,             &   ! =1, include (0,0) mode, =0, exclude (0,0) mode
                     maintnT,            &   ! =0, no T profile restore; =1, T profile restore (be careful!)
                     tstep,              &   ! time step size ((in unit L_T/V_th)
                     micell,             &   ! number of particles per cell
                     mpsi,               &   ! number of radial grid points
                     ipsi,               &   ! radial grid point where mode spectrum is calculated (normally, most unstable location)
                     mthetamax,          &   ! (standard) poloidal grid points are increased by a factor of mthetamax
                     mzetamax,           &   ! number of poloidal planes (toroidal grid points)
                     npartdom,           &   ! number of MPI tasks between which the particles are divided within a toroidal domain
                     ifelectron,         &   ! =0 or 1 or 2; 0: adiabatic electrons; 1: full delta-f electrons; 2: delta-h (non-adiabatic) electrons 
                     a,                  &   ! normalized minor radius, always =1, not really neceassary
                     a0,                 &   ! inner boundary (r/a) of radial domain 
                     a1,                 &   ! outer boundary (r/a) of radial domain
                     q0,                 &   ! abandoned; --> when q2<0 and plasms=0, use to specify grad(T) ~ -exp{-[(r-rc)/rw]**q0} 10/10
                     q1,                 &   ! abandoned
                     q2,                 &   ! abandoned  --> use another analytic T profile by specifying grad(T) when q2<0 and plasms=0 10/10
                     rc,                 &   ! abandoned
                     rw,                 &   ! abandoned
                     aion,               &   
                     qion,               &
                     aelectron,          &
                     qelectron,          &
                     kappati,            &   ! abandoned
                     kappate,            &   ! abandoned
                     kappan,             &   ! abandoned
                     flow0,              &   ! abandoned
                     flow1,              &   ! abandoned
                     flow2,              &   ! abandoned
                     r0,                 &   ! abandoned  --> use to specify the 1st analytic Ti profile made in Oct. 09
                     b0,                 &   ! abandoned  --> use to specify the 1st analytic Te profile made in Oct. 09
                     temconst,           &   ! abandoned
                     edensity0,          &   ! abandoned
                     n_0,                &   ! model density profile: 
                     dn,                 &   ! n=n_0*[An+Bn*exp(-dn*r**xn)]
                     xn,                 &
                     An,                 &
                     Bn,                 &
                     T_0,                &   ! model ion temperature profile:
                     dTemi,              &   ! Ti=T_0*[At+Bt*exp(-dTemi*r**xt)]
                     xt,                 &
                     At,                 &   
                     Bt,                 &
                     T_e0,               &   ! model electron temperature profile:
                     dTeme,              &   ! Te=T_e0*[At_e+Bt_e*exp(-dTeme*r**xt_e)]
                     xt_e,               &
                     At_e,               &   
                     Bt_e,               &
                     u_0,                &   ! model paralle flow profile:
                     du,                 &   ! u_//=u_0*[Au+Bu*exp(-du*r**xu)]*(B_0/B) with 1/B dependence
                     xu,                 &
                     Au,                 &
                     Bu,                 &
                     stdout,             &
                     nbound,             &
                     umax,               &
                     iload,              &   ! not really useful
                     icollision,         &   ! =0: no collision; =i>0: calculate collision every i steps
                     ecollision,         &   ! =0: no electron collision; =i>0: calculate electron collision every i steps
                     coll_enlarge,       &   ! enlarge collision fequency by a factor of coll_enlarge; coll_enlarge=1.0: real frequency
                     nss,                &   ! use smaller collision step-size (1/nss) for particles with small velocity
                     n,                  &   ! (m,n): specified modes for diagnosis
                     m,                  &
                     jr,                 &   ! number of radial nodes for plasma profiles and Er(?) (used for linear interpolation)
                     plasma,             &   ! =0, use specified model profiles; >0, experimental profiles
                     linterp_from_nodes, &   ! =0, use EZ spline (slow); >0, put on radial nodes first and then use linear interpolation
                     profile_data_name,  &   ! plasma profiles from TRANSP data base
                     Er_from_GTC_Neo,    &   ! =0, no equilibrium Er; >0, include equilibrium Er
                     n_node,             &   ! number of radial nodes for Er data
                     Er_data_name,       &   ! name of Er data file
                     movie_3D_h5,        &   ! =0, 3D-movie data off; >0, write 3D-movie data every movie_3D_h5 step
                     ZF_SEPARATION,      &   ! =1, solving zonal flow separately; =0, solve total potential (phi_00 + delta-phi)
                     order                   ! =0, original solver; =1, linear finite elements; =2 higher order elements

  namelist/ read2 /  gama0,              &   ! one of two coefficients used in solving Poisson iteratively
                     its,                &   ! number of iteration used for superLU interation in poisson solver
                     its_jacobi,         &   ! used to be for Jacobi iteration in poisson solver -> used control parameter for charge density conservation 05/12
                     nelectron,          &   ! number of sub-step  for pushing electrons; must take 1 or an even number due to 2nd RK used
                     nv,                 &   ! number of grids in vz and vper 
                     delta_v,            &   ! grid size in velocity space 
                     dradius,            &   ! box size in radial direction 
                     dtheta,             &   ! box size in theta direction 
                     rprobe,             &   ! radial locations of probes
                     tprobe,             &   ! poloidal locations of probes
                     gama1                   ! one of two coefficients used in solving Poisson iteratively


  if (mype==0) then
     open(unit=5,file='INPUT.d')
     read(unit=5,nml=read1)
     read(unit=5,nml=read2)
     close(5)

     char_params(1)(1:25) = eq_data_name
     char_params(2)(1:25) = profile_data_name
     char_params(3)(1:25) = Er_data_name
     !write(0,*)'char_params(1)=',char_params(1)
     integer_params(1)=irun
     integer_params(2)=mstep
     integer_params(3)=msnap
     integer_params(4)=ndiag
     integer_params(5)=mode00
     integer_params(6)=maintnT
     integer_params(7)=micell
     integer_params(8)=mpsi
     integer_params(9)=ipsi
     integer_params(10)=mthetamax
     integer_params(11)=mzetamax
     integer_params(12)=npartdom
     integer_params(13)=ifelectron
     integer_params(14)=stdout
     integer_params(15)=nbound
     integer_params(16)=iload
     integer_params(17)=icollision
     integer_params(18)=nss
     integer_params(19)=jr
     integer_params(20)=plasma
     integer_params(21)=linterp_from_nodes
     integer_params(22)=Er_from_GTC_Neo
     integer_params(23)=n_node
     integer_params(24)=movie_3D_h5
     integer_params(25)=ZF_SEPARATION
     integer_params(26)=its
     integer_params(27)=its_jacobi
     integer_params(28)=nelectron
     integer_params(29:36)=n
     integer_params(37:44)=m
     integer_params(45)=nv
     integer_params(46)=ecollision
     integer_params(47)= order

     real_params(1)=nonlinear
     real_params(2)=paranl
     real_params(3)=tstep
     real_params(4)=a
     real_params(5)=a0
     real_params(6)=a1
     real_params(7)=q0
     real_params(8)=q1
     real_params(9)=q2
     real_params(10)=rc
     real_params(11)=rw
     real_params(12)=aion
     real_params(13)=qion
     real_params(14)=aelectron
     real_params(15)=qelectron
     real_params(16)=kappati
     real_params(17)=kappate
     real_params(18)=kappan
     real_params(19)=flow0
     real_params(20)=flow1
     real_params(21)=flow2
     real_params(22)=r0
     real_params(23)=b0
     real_params(24)=temconst
     real_params(25)=edensity0
     real_params(26)=n_0
     real_params(27)=dn
     real_params(28)=xn
     real_params(29)=An
     real_params(30)=Bn
     real_params(31)=T_0
     real_params(32)=dTemi
     real_params(33)=xt
     real_params(34)=At
     real_params(35)=Bt
     real_params(36)=T_e0
     real_params(37)=dTeme
     real_params(38)=xt_e
     real_params(39)=At_e
     real_params(40)=Bt_e
     real_params(41)=umax
     real_params(42)=u_0
     real_params(43)=du
     real_params(44)=xu
     real_params(45)=Au
     real_params(46)=Bu
     real_params(47)=coll_enlarge
     real_params(48)=gama0
     real_params(49)=delta_v
     real_params(50)=dradius
     real_params(51)=dtheta
     real_params(52:54)=rprobe
     real_params(55:57)=tprobe
     real_params(58)=gama1

     write(0,*) 'equillibrium=', eq_data_name
     write(0,*) 'irun=', irun
     write(0,*) 'mstep=', mstep
     write(0,*) 'nonlinear=', nonlinear
     write(0,*) 'paranl=', paranl
     write(0,*) 'tstep=', tstep
     write(0,*) 'micell=', micell
     write(0,*) 'mpsi=', mpsi
     write(0,*) 'ipsi=', ipsi
     write(0,*) 'mthetamax=', mthetamax
     write(0,*) 'mzetamax=', mzetamax
     write(0,*) 'kappati=', kappati
     write(0,*) 'kappan=', kappan
     write(0,*) 'B0=', b0
     write(0,*) 'temconst=', temconst
     write(0,*) 'n_0=',n_0, ' dn=',dn, ' xn=',xn, ' An=', An, ' Bn=', Bn
     write(0,*) 'T_0=',T_0, ' dTemi=',dTemi, ' xt=',xt, ' At=', At, ' Bt=', Bt
     write(0,*) 'T_e0=',T_e0, ' dTeme=',dTeme, ' xt_e=',xt_e, ' At_e=', At_e, ' Bt_e=', Bt_e
     write(0,*) 'n=', n
     write(0,*) 'm=', m
     write(0,*) 'jr=', jr, 'print end'

     write(0,*) 'gama0=', gama0, 'its=', its,' nelectron=',nelectron     !cj...added....sept...28...2006
     write(0,*) 'finite element order=', order
  endif 

  call mpi_bcast(char_params,25*n_chars,MPI_CHARACTER,0,MPI_COMM_WORLD,ierror)
!  call mpi_bcast(eq_data_name,25,MPI_CHARACTER,0,MPI_COMM_WORLD,ierror)
!  call mpi_bcast(profile_data_name,25,MPI_CHARACTER,0,MPI_COMM_WORLD,ierror)
!  call mpi_bcast(Er_data_name,25,MPI_CHARACTER,0,MPI_COMM_WORLD,ierror)

  call mpi_bcast(integer_params,n_integers,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  call mpi_bcast(real_params,n_reals,MPI_Rsize,0,MPI_COMM_WORLD,ierror)

  if (mype/=0) then

     eq_data_name(1:25)      =char_params(1)(1:25)
     profile_data_name(1:25) =char_params(2)(1:25)
     Er_data_name(1:25)      =char_params(3)(1:25)

!     write(0,*)'eq_data_name=',eq_data_name

     irun              =integer_params(1)
     mstep             =integer_params(2)
     msnap             =integer_params(3)
     ndiag             =integer_params(4)
     mode00            =integer_params(5)
     maintnT           =integer_params(6)
     micell            =integer_params(7)
     mpsi              =integer_params(8)
     ipsi              =integer_params(9)
     mthetamax         =integer_params(10)
     mzetamax          =integer_params(11)
     npartdom          =integer_params(12)
     ifelectron        =integer_params(13)
     stdout            =integer_params(14)
     nbound            =integer_params(15)
     iload             =integer_params(16)
     icollision        =integer_params(17)
     nss               =integer_params(18)
     jr                =integer_params(19)
     plasma            =integer_params(20)
     linterp_from_nodes=integer_params(21)
     Er_from_GTC_Neo   =integer_params(22)
     n_node            =integer_params(23)
     movie_3D_h5       =integer_params(24)
     ZF_SEPARATION     =integer_params(25)
     its               =integer_params(26)
     its_jacobi        =integer_params(27)
     nelectron         =integer_params(28)
     n                 =integer_params(29:36)
     m                 =integer_params(37:44)
     nv                =integer_params(45)
     ecollision        =integer_params(46)
     order             =integer_params(47)

     nonlinear         =real_params(1)
     paranl            =real_params(2)
     tstep             =real_params(3)
     a                 =real_params(4)
     a0                =real_params(5)
     a1                =real_params(6)
     q0                =real_params(7)
     q1                =real_params(8)
     q2                =real_params(9)
     rc                =real_params(10)
     rw                =real_params(11)
     aion              =real_params(12)
     qion              =real_params(13)
     aelectron         =real_params(14)
     qelectron         =real_params(15)
     kappati           =real_params(16)
     kappate           =real_params(17)
     kappan            =real_params(18)
     flow0             =real_params(19)
     flow1             =real_params(20)
     flow2             =real_params(21)
     r0                =real_params(22)
     b0                =real_params(23)
     temconst          =real_params(24)
     edensity0         =real_params(25)
     n_0               =real_params(26)
     dn                =real_params(27)
     xn                =real_params(28)
     An                =real_params(29)
     Bn                =real_params(30)
     T_0               =real_params(31)
     dTemi             =real_params(32)
     xt                =real_params(33)
     At                =real_params(34)
     Bt                =real_params(35)
     T_e0              =real_params(36)
     dTeme             =real_params(37)
     xt_e              =real_params(38)
     At_e              =real_params(39)
     Bt_e              =real_params(40)
     umax              =real_params(41)
     u_0               =real_params(42)
     du                =real_params(43)
     xu                =real_params(44)
     Au                =real_params(45)
     Bu                =real_params(46)
     coll_enlarge      =real_params(47)
     gama0             =real_params(48)
     delta_v           =real_params(49)
     dradius           =real_params(50)
     dtheta            =real_params(51)
     rprobe            =real_params(52:54)
     tprobe            =real_params(55:57)
     gama1             =real_params(58)
     
     
  end if

  nmode=n
  mmode=m

  end subroutine read_input_params


!=============================================================================

    Subroutine set_particle_decomp

!=============================================================================

  use mpi_defs
  use global_parameters
  use particle_decomp

  implicit none
  integer  :: i,j,k,pe_number,mtest,ierror

! ----- First we verify the consistency of ntoroidal and npartdom -------
! The number of toroidal domains (ntoroidal) times the number of particle
! "domains" (npartdom) needs to be equal to the number of processor "numberpe".
! numberpe cannot be changed since it is given on the command line.

! numberpe must be a multiple of npartdom so change npartdom accordingly
  do while (mod(numberpe,npartdom) /= 0)
     npartdom=npartdom-1
     if(npartdom==1)exit
  enddo
  ntoroidal=numberpe/npartdom
  if(mype==0)then
    write(stdout,*)'*******************************************************'
    write(stdout,*)'  Using npartdom=',npartdom,' and ntoroidal=',ntoroidal
    write(stdout,*)'*******************************************************'
    write(stdout,*)
  endif

! We now give each PE (task) a unique domain identified by 2 numbers: the
! particle and toroidal domain numbers.
!    particle_domain_location = rank of the particle domain holding mype
!    toroidal_domain_location = rank of the toroidal domain holding mype
! 
! On the IBM SP, the MPI tasks are distributed in an orderly fashion to each
! node unless the LoadLeveler instruction "#@ blocking = unlimited" is used.
! On Seaborg for example, the first 16 tasks (mype=0-15) will be assigned to
! the first node that has been allocated to the job, then the next 16
! (mype=16-31) will be assigned to the second node, etc. When not using the
! OpenMP, we want the particle domains to sit on the same node because
! communication is more intensive. To achieve this, successive PE numbers are
! assigned to the particle domains first.
! It is easy to achieve this ordering by simply using mype/npartdom for
! the toroidal domain and mod(mype,npartdom) for the particle domain.
!
  particle_domain_location=mod(mype,npartdom)
  toroidal_domain_location=mype/npartdom

!  write(0,*)'mype=',mype,"  particle_domain_location =",&
!            particle_domain_location,' toroidal_domain_location =',&
!            toroidal_domain_location,' pi=',pi

! ---- Create particle domain communicator and toroidal communicator -----
! We now need to create a new communicator which will include only the
! processes located in the same toroidal domain. The particles inside
! each toroidal domain are divided equally between "npartdom" processes.
! Each one of these processes will do a charge deposition on a copy of
! the same grid, requiring a toroidal-domain-wide reduction after the
! deposition. The new communicator will allow the reduction to be done
! only between those processes located in the same toroidal domain.
!
! We also need to create a purely toroidal communicator so that the
! particles with the same particle domain id can exchange with their
! toroidal neighbors.
!
! Form 2 subcommunicators: one that includes all the processes located in
! the same toroidal domain (partd_comm), and one that includes all the
! processes part of the same particle domain (toroidal_comm).
! Here is how to create a new communicator from an old one by using
! the MPI call "MPI_COMM_SPLIT()".
! All the processes passing the same value of "color" will be placed in
! the same communicator. The "rank_in_new_comm" value will be used to
! set the rank of that process on the communicator.
!  call MPI_COMM_SPLIT(old_comm,color,rank_in_new_comm,new_comm,ierror)

! particle domain communicator (for communications between the particle
! domains WITHIN the same toroidal domain)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,toroidal_domain_location,&
                      particle_domain_location,partd_comm,ierror)

! toroidal communicator (for communications BETWEEN toroidal domains of same
! particle domain number)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,particle_domain_location,&
                      toroidal_domain_location,toroidal_comm,ierror)

  call mpi_comm_size(partd_comm,nproc_partd,ierror)
  call mpi_comm_rank(partd_comm,myrank_partd,ierror)

  call mpi_comm_size(toroidal_comm,nproc_toroidal,ierror)
  call mpi_comm_rank(toroidal_comm,myrank_toroidal,ierror)

!  write(0,*)'mype=',mype,'  nproc_toroidal=',nproc_toroidal,&
!       ' myrank_toroidal=',myrank_toroidal,'  nproc_partd=',nproc_partd,&
!       ' myrank_partd=',myrank_partd

  if(nproc_partd/=npartdom)then
    write(0,*)'*** nproc_partd=',nproc_partd,' NOT EQUAL to npartdom=',npartdom
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  if(nproc_toroidal/=ntoroidal)then
    write(0,*)'*** nproc_toroidal=',nproc_toroidal,' NOT EQUAL to ntoroidal=',&
              ntoroidal
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

! We now find the toroidal neighbors of the current toroidal domain and
! store that information in 2 easily accessible variables. This information
! is needed several times inside the code, such as when particles cross
! the domain boundaries. We will use the toroidal communicator to do these
! transfers so we don't need to worry about the value of myrank_partd.
! We have periodic boundary conditions in the toroidal direction so the
! neighbor to the left of myrank_toroidal=0 is (ntoroidal-1).

  left_pe=mod(myrank_toroidal-1+ntoroidal,ntoroidal)
  right_pe=mod(myrank_toroidal+1,ntoroidal)

end subroutine set_particle_decomp

