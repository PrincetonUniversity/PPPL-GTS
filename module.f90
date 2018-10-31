module mpi_defs
  include 'mpif.h'
end module mpi_defs

module precision
  use mpi_defs, only:MPI_DOUBLE_PRECISION,MPI_REAL,MPI_DOUBLE_COMPLEX
  integer, parameter :: mpi_Rsize=MPI_DOUBLE_PRECISION,&
                        mpi_complex_size=MPI_DOUBLE_COMPLEX

  integer, parameter :: doubleprec=selected_real_kind(12),&
                       &singleprec=selected_real_kind(6),&
                       &int8=selected_int_kind(10),&  !Can hold int >10 billion
                       &wp=doubleprec   ! working precision
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
  real(wp) nonlinear,a0,a1,a,q0,q1,q2,pi,twopi,tstep,kappati,kappate,&
           kappan,flow0,flow1,flow2,ulength,utime,E_00,gyroradius,deltar,&
           deltaz,zetamax,zetamin,umax,rc,rw,qion,qelectron,aion,aelectron,&
           uveloci,Tem_00,vthi,vthe,vthc,density_ave,finval_v,rho0,rho1,dt_coll,dt_coll_e, &
           s0,s1,B_0,R_0,Cs   ! B_0=normalized on-axis B, R_0=normalized major radius 
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
                        dsda_power=0.2_wp
                      !  pmass0=1.6734e-24_wp, charg0=4.8032e-10_wp,  &   ! hydrogen
                      !  amassi=1.6734e-24_wp, chargi=4.8032e-10_wp,  &   ! hydrogen
                      !  amasse=1.6734e-24_wp, charge=-4.8032e-10_wp, &  ! test mass ratio variation: mi/me=100

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
          ifdphi,zeondim

  real(wp)  msm0,m0ms,csc0,c0cs,for_rho,for_E,msm0_i,m0ms_i,csc0_i,c0cs_i
  real(wp)  ntracer_energy_ion0, ntracer_energy_electron0
  real(wp)  cloud_i, cloud_e
  real(wp),dimension(:),allocatable :: dadtdz

! The following input parameters are read from INPUT.d and used in setup
! to calculate derived quantities.
  integer Er_from_GTC_Neo,n_node
  real(wp) r0,b0,edensity0,temconst,coll_enlarge
  character(25) :: Er_data_name,profile_data_name,eq_data_name

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
  integer,parameter  :: wimcompressibility=1       ! switch  added April, 2013 for implementing weight imcompressibility
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

  integer :: electromagnetic
  real(wp) :: ga0, ga_rc, ga_rw, ga_pp, gagw, gagf0, dga0dt;
  integer :: gam, gan;

  !==============================================================================
  !MGYOO (10/24/2018) to control accuracy of gyro-average on Larmor orbit
  integer :: larmor_num = 16 ! how many points are used to do something on gyro orbit
  real(wp) inv_larmor_num ! (=1.0/larmor_num) For faster calculation
  !==============================================================================

end module global_parameters

module particle_array
  use precision
  integer,dimension(:),allocatable :: jteon0,jteon1
  integer,dimension(:,:),allocatable :: jtion0,jtion1
  real(wp),dimension(:),allocatable :: wzion,wzeon,wpeon,wteon0,wteon1
  real(wp),dimension(:,:),allocatable :: wpion,wtion0,wtion1
! Here we declare zion,zion0,zeon and zeon0 as array pointers
! so that we can pass them as dummy arguments to a subroutine and be allowed
! to reallocate them inside that subroutine if needed.
  real(wp),pointer,dimension(:,:) :: zion,zion0,zeon,zeon0,zeon00

  ! ============================
  ! MGYOO (07/31/2018)
  ! Information of Multiple-Tracers
  real(wp), dimension(:), allocatable :: elabels, ilabels
  integer :: num_of_etracers = 200, num_of_itracers =200, nstep_for_tracer_record = 5
  !real(wp), dimension(:), allocatable :: etracers_labels(num_of_etracers), itracers_labels(num_of_itracers)
  real(wp), dimension(:), allocatable :: etracers_labels, itracers_labels
  ! ===========================

  real(wp),dimension(:),allocatable :: enion
  real(wp),dimension(:),allocatable :: backn,backu,backT,dlnnda,dlnTda,&
                                       backne,backue,backTe,dlnndae,&
                                       dlnTdae,duda
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
  real(wp),dimension(:,:),allocatable :: phi,densityi,densitye,pgyro,tgyro,dphidt
  real(wp),dimension(:),allocatable :: tphi, bphi, left1, left2   ! for testing solver

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
  integer,parameter :: nn=13, ncolm=20




  real(wp) w2_ion0(nbin_a), w2_ion1(nbin_a), w2_eon0(nbin_a), w2_eon1(nbin_a)


end module diagnosis_array

module esi_interface
  use precision
   real(wp), dimension(:), allocatable :: F,Fa,gFa,gFaa,gYa,gYaa,TT,Ta,PP,Pa, &
          rr,ra,rq,z,za,zq,BB,Ba,Bq,gh,gha,ghq,sw
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

