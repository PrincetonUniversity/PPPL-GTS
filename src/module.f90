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
   integer micell,mi,mimax,me,ini_me, ini_mi, memax,mgrid,mpsi,mthetamax,& ! ,mzeta
      mzetamax,istep,ndiag,istart,movie_3D_h5,ntracer,msnap,mstep,mstepall,&
      stdout,mype,numberpe,mode00,nbound,irun,iload,irk,&
      idiag,ncycle,mtdiag,idiag1,idiag2,ipsi,maintnT,&
      nelectron,ntracer_e,ntracer00,me00,idiag_e
   real(wp) misum,mesum  !Total number of ions and electrons can be >2 billion
   real(wp) nonlinear,a0,a1,a,q0,q1,q2,pi,twopi,inv_twopi,tstep,kappati,kappate,&
      kappan,flow0,flow1,flow2,ulength,utime,E_00,gyroradius,deltar,&
      deltaz,zetamax,zetamin,umax,rc,rw,qion,qelectron,aion,aelectron,&
      uveloci,Tem_00,vthi,vthe,vthc,density_ave,finval_v,rho0,rho1,dt_coll,dt_coll_e, &
      s0,s1,B_0,R_0,Cs   ! B_0=normalized on-axis B, R_0=normalized major radius
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
   integer :: Flag_Ptl_BC_type=0 ! Default=0 (Original reflecting particle boundary), can be set =1 (Delete particles when hit the boundary)

   logical :: acc_newer_H_zion=.true., acc_newer_D_zion=.true., acc_newer_H_zeon=.true., acc_newer_D_zeon=.true.

   ! To measure the elapsed times of the following routines in the Pusheon
   ! (1) shifte, (2) charge_eon, (3) poisson, (4) Smooth_Phi (5) Diagnostics
   real(wp) time_in_Pushe(5), timewc_in_Pushe(5)
   ! To measure the elapsed times of the following routines in the Pushion
   ! (1) Diagnostics
   real(wp) time_in_Pushi(1), timewc_in_Pushi(1)

end module global_parameters

module FieldLineFollowing

   use precision
   type Poincare_Point_Information
      integer Fid  ! Fid (Field ID) for each field line
      real(wp) target_normalized_gf  ! target_zeta in [0,...,2*pi]
      real(wp) s, a, gq, gf
      real(wp) R, Z
      real(wp) dl, inv_B
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
   !$acc declare create (phi)
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


module Diagnostic_Grid
   use precision
   implicit none

   type Grid_Info
      integer :: gq_jump=3, a_jump=1
      integer :: dim_a, dim_gq
      integer :: num_eSD, num_eVD ! number of eScaData and eVecData (electrons)
      integer :: num_iSD, num_iVD ! number of iScaData and iVecData (ions)
      real(wp) zetamin, zetamax ! zeta value of this plane
      real(wp) dgq, inv_dgq ! dgq is constant since it is uniform along theta direction
      real(wp) inv_ds
      real(wp), allocatable, dimension(:) :: a_1D, ori_aid, s_1D, gq_1D, n0_1D, Te0_1D, Ti0_1D
      real(wp), allocatable, dimension(:,:) :: a_2D, gq_2D, n0_2D, Te0_2D, Ti0_2D
      real(wp), allocatable, dimension(:,:) :: Cell_dV, inv_Cell_dV ! Cell Volume and its inverse(Weighted Volume, but now it is approximated)
      real(wp), allocatable, dimension(:,:) :: Cell_dS_a, inv_Cell_dS_a, Cell_dS_gq, inv_Cell_dS_gq, Cell_dS_gf, inv_Cell_dS_gf ! Node_Surface in each direction
      real(wp), allocatable, dimension(:,:,:,:) :: eScaData, iScaData!
      real(wp), allocatable, dimension(:,:,:,:,:) :: eVecData,iVecData ! Flux..
      real(wp), allocatable, dimension(:,:) :: dne_adiab ! Electron adiabatic density for diagnostic purpose
      real(wp), allocatable, dimension(:,:,:) :: FLabel_Function
      real(wp), allocatable, dimension(:,:,:) :: FLA_Denominator, FLA_FLabel, FLA_psi0
      real(wp), allocatable, dimension(:,:) :: FLA_SM_Denominator, FLA_SM_FLabel, FLA_SM_psi0
      real(wp), allocatable, dimension(:,:) :: FLA_SM_FLabel_Left, FLA_SM_FLabel_Right

      ! Hybrid (full-f + df) density and temperature
      real(wp), allocatable, dimension(:,:) :: w_dne, w_dni ! wdne=abs(snap.DG.dne_tot)./(abs(snap.DG.ne_tot)+abs(snap.DG.dne_tot));
      real(wp), allocatable, dimension(:,:) :: ne_tot_hybrid, ni_hybrid, Te_hybrid, Ti_hybrid ! (n_hybrid = w_dne.*snap.DG.ne_tot+(1.0-w_dne).*(snap.DG.dne_tot+snap.n0_2D))


      real(wp), allocatable, dimension(:,:,:) :: mag_grad_FLA_SM_FLabel
      real(wp), allocatable, dimension(:,:,:,:) :: grad_FLA_SM_FLabel, uvec_grad_FL_cov

      ! Magnetic Field Information (ESI variables)
      real(wp), allocatable, dimension(:,:) :: F,dFda,dgFda,dgYda,R,dRda,dRdgq,Z,dZda,dZdgq,BB,dBBda,dBBdgq,P,T
      real(wp), allocatable, dimension(:,:) :: Jacobian, esi_N, esi_K, g_aa_cov, g_aq_cov, g_qq_cov, g_ff_cov
      real(wp), allocatable, dimension(:,:) :: D, g_aa_cntv, g_aq_cntv, g_qq_cntv, g_ff_cntv
      real(wp), allocatable, dimension(:,:,:) :: B0_cov3, B0_cntv3, j0_cov3, j0_cntv3

      ! Metrics and some vectors regarding s coordinates
      real(wp), allocatable, dimension(:,:) :: g_ss_cntv, g_sq_cntv, g_ss_cov, g_sq_cov, Jacob_sqf
      real(wp), allocatable, dimension(:,:,:) :: B0_cntv3_sqf

   end type


   type Probe_Velocity_Space
      real(wp) a, gq, gf ! center position of the probe
      real(wp) L_a, L_gq, L_gf ! length of the probe in each direction

      real(wp) n0, Te0, Ti0
      real(wp) inv_n0, inv_Te0, inv_Ti0
      real(wp) max_norm_v ! maximum normalized v (norm_v=v/v_th; v_th=sqrt(T/m))
      integer :: nv ! number of velocity (vpara=(-nv:nv), vperp=(0:nv))
      real(wp) dv, inv_dv ! each cell size
      real(wp) :: Cell_dA ! Area of each cell (dvpara*dvperp)
      integer :: dim_vpara, dim_vperp

      real(wp), allocatable, dimension(:) :: vpara_1D, vperp_1D  ! vpara_1D(-nv:nv), vperp_1D(0:nv)

      integer :: num_of_fi, num_of_fe
      real(wp), allocatable, dimension(:,:,:) :: fe, fi
   end type


   type(Grid_Info) :: DG
   type(Probe_Velocity_Space), allocatable, dimension(:) :: PVS
   integer :: nPVS
end module Diagnostic_Grid

module DG_operations

   use Diagnostic_Grid

   contains
      subroutine Smooth_DG_Data_2D(num_smooth, target_Data2D, Reference_Data2D)

         use precision
         use global_parameters
         use Diagnostic_Grid
         implicit none

         integer, intent(in) :: num_smooth
         real(wp),intent(inout) :: target_Data2D(DG%dim_gq, DG%dim_a)
         real(wp),optional:: Reference_Data2D(DG%dim_gq, DG%dim_a)
         real(wp) :: VV(DG%dim_gq, DG%dim_a)
         integer :: smid, aid, gqid, Na,Ngq, vid,k
         real(wp) :: SM_Data(0:DG%dim_gq+1, 0:DG%dim_a+1)  ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
         real(wp) :: SM_gq_Coeff(-1:1,0:DG%dim_gq+1, 0:DG%dim_a+1), SM_a_Coeff(-1:1,0:DG%dim_gq+1, 0:DG%dim_a+1)

         Na  = DG%dim_a
         Ngq = DG%dim_gq

         VV=1.0
         if(present(Reference_Data2D)) VV=Reference_Data2D

         ! Get Smooth Coeffs which conserves Volume Weighted Integration
         call Get_Conserving_Smooth_Coeffs_on_DG(VV, SM_a_Coeff, SM_gq_Coeff)

         SM_Data = 0.0
         SM_Data(1:Ngq,1:Na) = target_Data2D
         do smid=1,num_smooth
            ! Set ghost cell periodic in theta direction
            SM_Data(0,:)     = SM_Data(Ngq,:)
            SM_Data(Ngq+1,:) = SM_Data(1,:)
            ! Smooth in theta direction
            SM_Data(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*SM_Data(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*SM_Data(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*SM_Data(2:Ngq+1,:)
            ! Smooth in radial direction
            SM_Data(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*SM_Data(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*SM_Data(:,1:Na) + SM_a_Coeff(1,:,1:Na)*SM_Data(:,2:Na+1)
         enddo
         target_Data2D = SM_Data(1:Ngq,1:Na)

      end subroutine Smooth_DG_Data_2D
end module DG_operations


module esi_interface
   use precision
   real(wp), dimension(:), allocatable :: F,Fa,gFa,gFaa,gYa,gYaa,TT,Ta,PP,Pa, &
      rr,ra,rq,z,za,zq,BB,Ba,Bq,gh,gha,ghq,sw
   !$acc declare create(F,Fa,TT,Ta,PP,Pa,BB,Ba,Bq,gh,gha,ghq)
   
   ! gradient of (a,gq,gf) at electron and ion particle positino
   real(wp), dimension(:), allocatable :: ee_grad_a, ee_grad_gq, ee_grad_gf
   real(wp), dimension(:), allocatable :: ii_grad_a, ii_grad_gq, ii_grad_gf
   !$acc declare create(ee_grad_a, ee_grad_gq, ee_grad_gf)  
   !$acc declare create(ii_grad_a, ii_grad_gq, ii_grad_gf)  

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
