module mpi_defs
  use mpi
end module mpi_defs

module precision
   use ieee_arithmetic
   use mpi_defs, only:MPI_DOUBLE_PRECISION,MPI_REAL,MPI_DOUBLE_COMPLEX
   integer, parameter :: mpi_Rsize=MPI_DOUBLE_PRECISION,&
      mpi_complex_size=MPI_DOUBLE_COMPLEX

   integer, parameter :: doubleprec=selected_real_kind(12),singleprec=selected_real_kind(6),int8=selected_int_kind(10), wp=doubleprec
end module precision

module global_parameters
   use precision
   integer,parameter :: ihistory=22,snapout=33
   integer micell,mi,mimax,me,ini_me,ini_mi,ini_me_npartd, ini_mi_npartd, memax,mgrid,mpsi,mthetamax,& ! ,mzeta
      mzetamax,istep,ndiag,istart,movie_3D_h5,ntracer,msnap,mstep,mstepall,&
      stdout,mype,numberpe,mode00,nbound,irun,iload,irk,&
      idiag,ncycle,mtdiag,idiag1,idiag2,ipsi,maintnT,&
      nelectron,ntracer_e,ntracer00,me00,idiag_e,ierrr
   integer :: nbound_inner=-100, nbound_outer=-100
   real(wp) Gt_start, Gt_end ! for measuring the elpased time for something
   real(wp) misum,mesum  !Total number of ions and electrons can be >2 billion
   real(wp) nonlinear,a0,a1,a,q0,q1,q2,pi,twopi,inv_twopi,tstep,kappati,kappate,&
      kappan,flow0,flow1,flow2,ulength,utime,E_00,gyroradius,deltar,&
      deltaz,zetamax,zetamin,umax,rc,rw,qion,qelectron,aion,aelectron,&
      uveloci,Tem_00,vthi,vthe,vthc,density_ave,finval_v,rho0,rho1,dt_coll,dt_coll_e, &
      s0,s1,B_0,R_0,Cs   ! B_0=normalized on-axis B, R_0=normalized major radius
   !$acc declare create(a0,a1)
   !$acc declare create(s0,mzetamax)

   real(wp) MM_Cs, MM_Lc_over_vthci, MM_LcTi, MM_LcTe, MM_Lcn, MM_Lcu, MM_vthci, MM_vthce, MM_R_axis, MM_B_axis ! Variables for MeasureMent
   integer,dimension(:),allocatable :: mtheta
   real(wp),dimension(:),allocatable :: deltat

   real(wp),parameter :: B_00=1.0e4_wp,  &  ! code unit for B, in unit gauss
      R_00=100._wp,  &  ! code unit for length, in cm
      c=2.9979e10_wp, &
      pmass0=3.3452e-24_wp, charg0=4.8032e-10_wp,  &  ! Deuturium
      amassi=3.3452e-24_wp, chargi=4.8032e-10_wp,  &  ! Deuturium
      amasse=9.1095e-28_wp, charge=-4.8032e-10_wp, &
      trans=1.6022e-12_wp, &  ! energy (erg) of 1eV
      dlni=16.0_wp, &
      !dsda_power=0.2_wp
   dsda_power =0.0
   !  pmass0=1.6734e-24_wp, charg0=4.8032e-10_wp,  &   ! hydrogen
   !  amassi=1.6734e-24_wp, chargi=4.8032e-10_wp,  &   ! hydrogen
   !  amasse=1.6734e-24_wp, charge=-4.8032e-10_wp, &  ! test mass ratio variation: mi/me=100

   real(wp) :: eNmNr, eNrNm, iNmNr, iNrNm ! (electron & ion)'s (Nm/Np, Np/Nm), where Nm and Np are total number of markers and real particles

   ! When particle weight exceeds threshold value, it would not contribute to charge or current deposition
   ! Default values are  w_min=-50 (df/f=-0.9803892168), w_max=0.9803921568 (df/f=50)
   ! These values can be controlled in INPUT.d file
   integer :: Flag_w_Threshold = 0 ! Default=0 (OFF), can be turned on (=1;ON) via INPUT.d, (=2: Abandon abnormal particles)
   integer :: Flag_v_Threshold = 0 ! Default=0 (OFF), can be turned on (=1;ON: Abandon abnormal particles)
   real(wp) :: w_Thold_min = -50.0, w_Thold_max = 1.0

   integer jr

   real(wp)  T_0,dTemi,xt,At,Bt,Tc,T_e0,dTeme,xt_e,At_e,Bt_e,Tc_e
   real(wp)  n_0,dn,xn,An,Bn
   real(wp)  du,xu,Au,Bu,u_0

   real(wp),dimension(:),allocatable :: Yr, coty, Er401_pc, Er401_bp, Er401_bp0
   real(wp) Ydr
   integer je
   integer seed

   real(wp),dimension(:),allocatable :: tite
   integer plasma,linterp_from_nodes,icollision,nss,ZF_SEPARATION,ifelectron,ecollision,  &
      ifdphi,zeondim,ziondim

   real(wp) :: adiab=0.0 !(Default=0.0) In the case of delta-f electrons, this is used for poisson equation for better convergence

   real(wp)  msm0,m0ms,csc0,c0cs,for_rho,for_E,msm0_i,m0ms_i,csc0_i,c0cs_i
   real(wp)  ntracer_energy_ion0, ntracer_energy_electron0
   real(wp)  cloud_i, cloud_e
   real(wp),dimension(:),allocatable :: dadtdz

   ! The following input parameters are read from INPUT.d and used in setup
   ! to calculate derived quantities.
   integer Er_from_GTC_Neo,n_node
   real(wp) r0,b0,edensity0,temconst,coll_enlarge
   character(128) :: Er_data_name,profile_data_name,eq_data_name

   ! the following vaiables and arrays are used for diagnosis of distribution function
   integer,parameter :: nprobe=3
   real(wp)  rprobe(nprobe),tprobe(nprobe),uflow(nprobe),vthprobe(nprobe),vthprobe_e(nprobe),  &
      sn_box(nprobe),density_in_box(nprobe), &
      Ti_box(nprobe),Te_box(nprobe)
   integer   nv
   real(wp)  dradius, dtheta, delta_v

   ! Diagnostic Grid, PVS information
   integer :: num_of_PVS=0
   integer :: PVS_max_v=4, PVS_nv=20
   real(wp):: PVS_a0(100), PVS_a1(100), PVS_gq0(100), PVS_gq1(100) ! at zetamax (temporal storage here)


   ! The following input parameters are specified in setup for radial boundary layers with damping
   real(wp) lmin, lmax, rmin, rmax

   ! The following  parameters are used to sepcify weight damping implementation
   !     wimcompressibility=0 --> no damping
   !     wimcompressibility=1 --> use the algorithm to assure imcompressibility in weight
   !     wimcompressibility=2 --> apply weight damping with damping rate specified by eta_w=weight_damping_rate/tstep
   integer,parameter  :: wimcompressibility=0       ! switch  added April, 2013 for implementing weight imcompressibility
   !     integer  :: wimcompressibility=1       ! switch  added April, 2013 for implementing weight imcompressibility
   real(wp),parameter :: weight_damping_rate=2.0d-7 ! =eta_w*tstep

   ! finite element order
   ! order = 0 original solver
   ! order = 1 linear finite element
   ! order = 2 2nd order finite element
   integer ::  order=0

   ! boundary type
   ! bc_type = 0 neumann boundary conditions at both a=a0 and a=a1
   ! bc_type = 1 dirichlet boundary condition at a=a0, neumann boundary conditions at a=a1
   ! bc_type = 2 dirichlet boundary conditions at both a=a0 and a=a1
   ! bc_type = 3 neumann boundary condition at a=a0, dirichlet boundary conditions at a=a1
   integer   bc_type

   ! The following  parameters are used to sepcify the inclusion of neoclassical physics
   !     ifneoclassical=0 --> no neoclassical include
   !     ifneoclassical=1 --> include neoclassical
   ! integer,parameter  :: ifneoclassical=1
   integer:: ifneoclassical  !--ZHIXIN: this is passed from INPUT.d

   integer,dimension(:),allocatable :: mappingL2G, mappingG2L

   integer :: electromagnetic, Flag_Streaming=1 ! Default Value: Flag_Streaming=1(ON)
   integer :: num_of_EMP
   type External_Magnetic_Perturbation
      integer :: gam, gan
      real(wp) :: ga0, ga_rc, ga_rw, ga_pp, gagw, gagf0, dga0dt
   end type
   type(External_Magnetic_Perturbation), dimension(:), allocatable :: EMP
   !$acc declare create(electromagnetic, num_of_EMP)
   !$acc declare create(EMP)

   integer :: num_of_FLabel_pcws
   real(wp),dimension(:),allocatable :: FLabel_pcws_a, FLabel_pcws_y ! piecewise_radial_base_function for FLabel

   ! real(wp) :: ga0, ga_rc, ga_rw, ga_pp, gagw, gagf0, dga0dt;
   ! integer :: gam, gan;

   !==============================================================================
   !MGYOO (10/24/2018) to control accuracy of gyro-average on Larmor orbit
   integer :: larmor_num = 4 ! how many points are used to do something on gyro orbit
   real(wp) inv_larmor_num ! (=1.0/larmor_num) For faster calculation
   !==============================================================================

   integer :: Flag_RK4=0 ! Default value is 0 (OFF), it can be controlled by INPUT.d file
   integer :: Flag_Poisson=1 ! Default=1 (ON), if this is 0 then does not solve poisson equation so becomes the single particle simulation.
   integer :: Flag_Poisson3D_on_SG=0 ! Default=0 (OFF; using original equation), can be set =1 (new 3D Poisson Solver on SG)
   integer :: Flag_Def_dne_adiab=0 ! Default=0 (original definition; using nonzonal phi), can be set =1 (using full phi with Te0), =2 (using full phi with FLA_Te0)
   integer :: Flag_w_Evol=0  ! Default=0 (OFF;original equation), can be set =1(Transformed Weight Evolution Equation) by INPUT.d
   integer :: Flag_GC_Solver=1 ! Default=1 (Vector form GC equation), can be set =0(Matrix form GC equation)
   integer :: Flag_Load_and_WeightDef=0 ! Default=0 (original), can be set =1 (Uniform Loading and double-weight definition for df and f_tot)
   integer :: Flag_Quiet_Start=1 ! Default=1 (Quiet Start), can be set =0 (Non-Quiet Start)
   integer :: Flag_Ptl_BC0_type=0 ! Default=0 (Original reflecting particle boundary), can be set =1 (Delete particles when hit the boundary)
   integer :: Flag_Ptl_BC1_type=0 ! Default=0 (Original reflecting particle boundary), can be set =1 (Delete particles when hit the boundary)
   real(wp) :: Ptl_BC0_val=1.0, Ptl_BC1_val=0.0 ! Default value of Particle Flux Generating at both edge (0.0=loss boundary; 1.0=flux from original n0+dne_adiab)
   integer :: Flag_Update_Pol_Shielding=0 ! Default=0 (Not updating Polarization Shielding over time)
   integer :: Flag_DK_Ion=0 ! Default=0 (GyroKinetic Ion), can be set =1 (DriftKinetic Ion; zero Gyro-radius)
   integer :: Flag_Measure_Ptl_FootPrint=0 ! Default=0(Not Meausre), can be set =1 (Do Measure)
   integer :: Flag_noExB=0 ! Default=0 (ExB is considered), can be set =1 (noExB)
   integer :: Flag_Recycling=0 ! Default=0 (no Recycling), can be set = 1 (Recycling of Ion at outer wall)
   integer :: Flag_Weight_dfh=0 ! Default=0 (delta-f), can be set = 1(full-f or full-h), or =2(hybrid two-weight; using both delta and full weight)
   integer :: Ptl_Src_nstep=50 ! Generate src from boundaires at every Ptl_Src_nstep (default=50 steps)

   integer :: MaxLoopNum_Adiab_Weight_Update=1 ! Default=1 (number of loop for updating adiabatic weight evolution)

   logical :: acc_newer_H_zion=.true., acc_newer_D_zion=.true., acc_newer_H_zeon=.true., acc_newer_D_zeon=.true.

   ! To measure the elapsed times of the following routines in the Pusheon
   ! (1) shifte, (2) charge_eon, (3) poisson, (4) Smooth_Phi (5) Diagnostics
   real(wp) time_in_Pushe(5), timewc_in_Pushe(5)
   ! To measure the elapsed times of the following routines in the Pushion
   ! (1) Diagnostics, (2) shifti (in the case of pushion_wEvol3)
   real(wp) time_in_Pushi(2), timewc_in_Pushi(2)


end module global_parameters

module FieldLineFollowing

   use precision
   type Poincare_Point_Information
      integer Fid  ! Fid (Field ID) for each field line
      real(wp) target_normalized_gf  ! target_zeta in [0,...,2*pi]
      real(wp) s, a, gq, gf
      real(wp) R, Z
      real(wp) dl, inv_B, B
      real(wp) ini_FLabel ! FLabel_Function will be used to estimate modified flux-surface labels
   end type

   integer :: FLF_max_tor_num = 400, gPPI_max_tor_num = 200

   integer :: num_of_SampledFieldLines
   integer, dimension(:), allocatable :: num_of_gPPI_for_each_FieldLine
   ! Global PPI ==> gPPI
   ! gPPI(gPPI_max_tor_rotation, 1:num_of_SampledFieldLines)
   type(Poincare_Point_Information), dimension(:,:), allocatable:: gPPI

end module FieldLineFollowing
module particle_array
   use precision
   integer,dimension(:),allocatable :: jteon0,jteon1
   integer,dimension(:,:),allocatable :: jtion0,jtion1
   real(wp),dimension(:),allocatable :: wzion,wzeon,wpeon,wteon0,wteon1
   real(wp),dimension(:,:),allocatable :: wpion,wtion0,wtion1
   !$acc declare create(jtion0,jtion1,wzion,wpion,wtion0,wtion1)
   !$acc declare create(jteon0,jteon1,wzeon,wpeon,wteon0,wteon1)
   real(wp),dimension(:,:),allocatable :: zion_s_on_GyroOrbit, zion_gq_on_GyroOrbit
   !$acc declare create(zion_s_on_GyroOrbit,zion_gq_on_GyroOrbit)
   ! Here we declare zion,zion0,zeon and zeon0 as array pointers
   ! so that we can pass them as dummy arguments to a subroutine and be allowed
   ! to reallocate them inside that subroutine if needed.
   real(wp),dimension(:,:), allocatable :: zion,zion0,zeon,zeon0,zeon00
   !$acc declare create(zion,zion0,zeon,zeon0)

   real(wp), dimension(:,:),allocatable :: RK_cum_zion, RK_cum_zeon
   !$acc declare create(RK_cum_zion,RK_cum_zeon)
   ! ============================
   ! MGYOO (07/31/2018)
   ! Information of Multiple-Tracers
   real(wp), dimension(:), allocatable :: elabels, ilabels
   !$acc declare create(elabels,ilabels)
   real(wp) :: elabel_of_ntracer_e, ilabel_of_ntracer
   integer :: num_of_etracers=200, num_of_itracers=200, nstep_for_tracer_record=5 ! Default values
   !real(wp), dimension(:), allocatable :: etracers_labels(num_of_etracers), itracers_labels(num_of_itracers)
   real(wp), dimension(:), allocatable :: etracers_labels, itracers_labels
   ! Unnormalized zeta&theta (not normalized into [0,2*pi))
   real(wp), dimension(:), allocatable :: UnNorm_zeta_e, UnNorm_theta_e, UnNorm_zeta_i, UnNorm_theta_i  ! These are not perfectly consistent when the particle touches the boundaries
   !$acc declare create(UnNorm_zeta_e, UnNorm_theta_e, UnNorm_zeta_i, UnNorm_theta_i)
   ! ===========================

   ! ============================
   ! MGYOO (07/31/2018)
   ! Additional information about weight component ratio
   real(wp), dimension(:), allocatable :: wi_init, wi_vpara, wi_vperp, wi_accel, wi_fdamp ! for ion weight
   !$acc declare create(wi_init, wi_vpara, wi_vperp, wi_accel, wi_fdamp)
   real(wp), dimension(:), allocatable :: we_init, we_vpara, we_vperp, we_accel, we_adiab, we_fdamp ! for electron weight
   !$acc declare create(we_init, we_vpara, we_vperp, we_accel, we_adiab, we_fdamp)
   ! ============================


   real(wp),dimension(:),allocatable :: enion
   real(wp),dimension(:),allocatable :: backn,backu,backT,dlnnda,dlnTda,&
      backne,backue,backTe,dlnndae,&
      dlnTdae,duda
   real(wp),dimension(:),allocatable :: back_f0_i, back_f0_e
   !$acc declare create(backn,backu,backT,dlnnda,dlnTda,backne,backue,backTe,dlnndae,dlnTdae,duda)
   !$acc declare create(back_f0_i,back_f0_e)
end module particle_array

module field_array
   use precision
   integer,dimension(:),allocatable :: itran,igrid
   integer,dimension(:,:),allocatable :: jtp1,jtp2
   real(wp),dimension(:),allocatable :: phi00,phip00,rtemi,rteme,rden,qtinv,&
      pmarki,pmarke,zonali,zonale,psi_grid,plenth,allnodes,r_grid,q_grid,&
      den_grid,vth_grid,vth_inv,gaa,B2_s,dVda,rho_grid,s_grid,eta_grid,&
      Er_grid,Er_grid0,dphidt00,ddT_da,Te_grid,gsgaa,&
      dni_grid,dne_grid
   real(wp),dimension(:),allocatable :: zonal_dTe, zonal_dTi
   real(wp),dimension(:,:),allocatable :: tot_dne,tot_dni, FLA_Denominator, FLA_dphi, FLA_dTe, FLA_dTi, FLA_dne, FLA_dni, FLA_ne_adiabatic
   real(wp),dimension(:,:),allocatable :: FLA_psi0, FLA_Te0, FLA_Ti0, FLA_n0, FLA_FLabel
   real(wp),dimension(:,:),allocatable :: FLabel_Function
   real(wp),dimension(:,:),allocatable :: Weighted_Vol
   real(wp),dimension(:),allocatable :: Node_simple_Vol
   real(wp),dimension(:,:),allocatable :: phi,densityi,densitye,pgyro,tgyro,dphidt
   real(wp),dimension(:),allocatable :: test_src ! test charge src in ratio ((ni-ne)/n0)
   real(wp),dimension(:,:),allocatable :: old_phi, tmp_old_phi ! phi at previous time step
   real(wp),dimension(:),allocatable :: old_phi00, tmp_old_phi00 ! phi00 at previous time step
   real(wp),dimension(:),allocatable :: cum_phi00, cum_phi00_SQ, cum_zonali, cum_zonali_SQ, cum_zonale, cum_zonale_SQ
   real(wp),dimension(:),allocatable :: cum_zonal_dTe, cum_zonal_dTe_SQ, cum_zonal_dTi, cum_zonal_dTi_SQ
   real(wp),dimension(:,:),allocatable :: cum_full_dphi, cum_full_dphi_SQ, cum_full_dne, cum_full_dne_SQ, cum_full_dni, cum_full_dni_SQ
   real(wp),dimension(:,:),allocatable :: cum_full_dTe, cum_full_dTe_SQ, cum_full_dTi, cum_full_dTi_SQ

   real(wp),dimension(:),allocatable :: tphi, bphi, left1, left2   ! for testing solver

   real(wp),dimension(:,:),allocatable :: dni_init, dni_vpara, dni_vperp, dni_accel, dni_fdamp
   real(wp),dimension(:,:),allocatable :: dne_init, dne_vpara, dne_vperp, dne_accel, dne_adiab, dne_fdamp
   real(wp),dimension(:,:),allocatable :: dTi_init, dTi_vpara, dTi_vperp, dTi_accel, dTi_fdamp
   real(wp),dimension(:,:),allocatable :: dTe_init, dTe_vpara, dTe_vperp, dTe_accel, dTe_adiab, dTe_fdamp
   ! Added for flow visualization, Feb. 2011
   real(wp),dimension(:,:),allocatable :: on_grid, flow_grid

   real(wp),dimension(:),allocatable :: markeri,markere,dtemper,heatflux,&
      phi_tmp,Jcbn,gsij,div1,div2,     &
      dRda, dRdt, dZda, Dzdt, gbeta
   real(wp),dimension(:,:,:),allocatable :: evector,pflux_i,pflux_e
   real(wp),dimension(:,:),allocatable :: wtp1,wtp2,gij,cij

   !integer,dimension(:,:),allocatable :: nindex
   !integer,dimension(:,:,:),allocatable :: indexp
   !real(wp),dimension(:,:,:),allocatable :: ring
   real(wp) gama0, gama1
   integer(wp) its,its_jacobi
   real(wp) :: its_phip00_tolerance=1e-5 ! tolerance for inner iteraions of zonal and nonzonal possion solve

   real,dimension(:),allocatable :: save_flux_ave

   real(wp) Er_a0, Er_a1


end module field_array

module diagnosis_array
   use precision
   integer,parameter :: mflux=5,num_mode=8,m_poloidal=9,nrbin=5,nbin_a=10
   integer nmode(num_mode),mmode(num_mode)
   real(wp) efluxi,efluxe,pfluxi,pfluxe,ddeni,ddene,dflowi,dflowe,entropyi,&
      entropye,efield,eradial,eflux(mflux),rmarker(mflux),&
      amp_mode(2,num_mode,2),rdtem(nrbin),rtden(nrbin),rdtem_all(nrbin),&
      rdtem0(nrbin),tem_nrbin(nrbin)   !wx
   real(wp) sn0(nbin_a),dnty_on_bin(nbin_a),Tr(nbin_a),Tr_e(nbin_a), &
      dnty_r(nbin_a),dO_10(nbin_a),dT_10(nbin_a)
   real(wp),dimension(:),allocatable :: fluxpsi,jaco_sum,En_00,En_tur,En_tot,sn0_g,dnty_on_bin_g
   real(wp),dimension(:),allocatable :: sumfluxpsi,flux_11
   real(wp),dimension(:,:,:),allocatable :: eigenmode
   integer,dimension(:),allocatable :: nf,nfL,nf_e,nfL_e
   integer my_nf, my_nfL, my_ZL(2),nfile
   integer,parameter :: nn=17, ncolm=20




   real(wp) w2_ion0(nbin_a), w2_ion1(nbin_a), w2_eon0(nbin_a), w2_eon1(nbin_a)


end module diagnosis_array




module esi_interface
   use precision
   real(wp), dimension(:), allocatable :: F,Fa,gFa,gFaa,gYa,gYaa,TT,Ta,PP,Pa, &
      rr,ra,rq,z,za,zq,BB,Ba,Bq,gh,gha,ghq,sw
   !$acc declare create(F,Fa,TT,Ta,PP,Pa,BB,Ba,Bq,gh,gha,ghq)

   ! gradient of (a,gq,gf) at electron and ion particle positino
   real(wp), dimension(:), allocatable :: ee_inv_grad_a, ee_inv_grad_gq, ee_inv_grad_gf
   real(wp), dimension(:), allocatable :: ii_inv_grad_a, ii_inv_grad_gq, ii_inv_grad_gf
   !$acc declare create(ee_inv_grad_a, ee_inv_grad_gq, ee_inv_grad_gf)
   !$acc declare create(ii_inv_grad_a, ii_inv_grad_gq, ii_inv_grad_gf)

   real(wp), dimension(:), allocatable :: dve_gf, dvi_gf  ! gf(zeta) component of delta_v (due to plasma potential and magnetic perturbation)
   !$acc declare create(dve_gf,dvi_gf)
   real(wp), dimension(:), allocatable :: vpara_dvi_a, vpara_dvi_gq,vpara_dvi_gf  ! delta components (perturbed part) of ion v_para
   !$acc declare create (vpara_dvi_a, vpara_dvi_gq,vpara_dvi_gf)
   real(wp), dimension(:), allocatable :: vpara_dve_a, vpara_dve_gq,vpara_dve_gf  ! delta components (perturbed part) of electron v_para
   !$acc declare create (vpara_dve_a, vpara_dve_gq,vpara_dve_gf)
   real(wp)  metric(6)
   real(wp)  djr
   real(wp), dimension(:), allocatable :: flux_it,dlnA_da,flux_A
end module esi_interface


!module linear_system
!  use precision
!   integer superlu_iopt
!   integer superlu_nprocs
!   integer superlu_n
!   integer superlu_nnz
!   integer superlu_nrhs
!   real(wp),dimension(:,:),allocatable :: superlu_csr
!   integer,dimension(:,:),allocatable :: superlu_colind
!   integer,dimension(:,:),allocatable :: superlu_rowptr
!   real(wp),dimension(:,:),allocatable :: superlu_csc
!   integer,dimension(:,:),allocatable :: superlu_rowind
!   integer,dimension(:,:),allocatable :: superlu_colptr
!   real(wp),dimension(:,:),allocatable :: superlu_rhs
!   real(wp),dimension(:,:),allocatable :: one
!   integer*8 superlu_factors  ! 64 bit
!   integer superlu_info
!end module linear_system

module particle_decomp
   integer  :: ntoroidal,npartdom
   integer  :: partd_comm,nproc_partd,myrank_partd
   integer  :: toroidal_comm,nproc_toroidal,myrank_toroidal
   integer  :: left_pe,right_pe
   integer  :: toroidal_domain_location,particle_domain_location
end module particle_decomp


!ZHIXIN
module control_parameters
   use precision
   real(wp) ratio_dthe0,ratio_dthe1,a0_dthe,a1_dthe,a_ref_in
   real(wp) ratio_Er401,ratio_Er_grid,w_initial_in
   integer flag_Lc,iratio_dthe,irewrite
   !-- switch for recording corresponding variables
   integer flag_DEN,flag_FLOW,flag_PHI,flag_RESTART,flag_GRID
   real(wp) flip_rotation,ain,aout   ! WX added Nov. 2014

end module control_parameters
