&read1           ! short description, first given is the "more relevant" case
  eq_data_name='esiP.00_cyclon', ! equilibrium file name

  !general major model used: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nonlinear=1.0d0, ! nonlinear for 1.0, linear for 0.0
  maintnT=0,     ! don't restore T profile if =0, do if =1 (caution!)
  ifelectron=2,  ! electron model: 0=adiabatic, 1=full delta-f, 2=delta-h (non-adiabatic)
  ifdphi=2,      ! usage of dphi/dt equation, related to ifelectron
  mode00=1,      ! include (0,0) mode if =1, not if =0
  ZF_SEPARATION=1, ! =1, solving zonal flow separately; =0, solve total potential (phi_00 + delta-phi)
  order=1        ! =0, original solver; =1, linear finite elements; =2 higher order elements (for ZonalFlows only!?)

  ! particle physics: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  aion=1.0d0,      ! ion mass [normalized by Deuterium mass 3.3452e-27 kg]
  qion=1.0d0,      ! ion charge [in units of abs(elementary charge)]
  aelectron=3674.0d0, !inverse electron mass [normalized by Deuterium mass]
  qelectron=-1.0d0,! electron charge [in units of abs(elementary charge)] (usually =-1)


  ! Numerics: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  irun=0,        ! restart run if >0
  mstep=100000,   ! number of timesteps
  tstep=0.0020d0,   ! time step size [in unit L_T/v_th]
  micell=50,   ! number of markers per cell for first ion species
  mpsi=100,      ! number of radial grid points
  ipsi=55,       ! radial grid point, where mode spectrum is calculated (usually most unstable region)
  mthetamax=1,   ! factor by which (standard) poloidal grid points are increased
  mzetamax=16,    ! number of toroidal grid points (=poloidal planes)
  npartdom=12,   ! number of MPI tasks between which the particles are divided within a toroidal domain
                 ! mzetamax*npartdom is the number of total cores for the simulation!
  nbound=5,      ! number of boundary radial cells. ?? some arrays (that span the radius) don't go beyond mpsi-nbound !???
  umax=4.0d0, ! cutoff velocity [in terms of v_th] during the initial particle loading
     
  a=1.0d0,       ! normalized minor radius (not really necessary)
  a0=0.1d0,      ! inner boundary of radial domain [in r/a]
  a1=0.9d0,      ! outer boundary of radial domain [in r/a]

  iload=0,       ! different loading schemes, but so far, only =0 implemented!!

  ! Output settings: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  stdout=6,      ! file pointer for job output file (why?)
  msnap=5,      ! step interval for writing out PHI record
  ndiag=2,       ! step interval for fluxes, tur intensity, energy
  movie_3D_h5=0, ! =0, 3D-movie data off; >0, write 3D-movie data every movie_3D_h5 step


  ! Collisions: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  icollision=0,  ! ion collisions are included if >0. Collisions are calculated every 'icollision'-step
  ecollision=0,  ! electron collisions are included if >0. Collisions are calculated every 'ecollision'-step
  coll_enlarge=1.0d0, ! collision frequency [in real frequency]: enlarge collision fequency by a factor of 'coll_enlarge'
  nss=10,        ! use smaller collision step-size (1/nss) for particles with small velocity


  ! Modes to diagnose: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  n=1,2,2,2,3,3,4,4       ! specified modes (m,n) for diagnosis: toroidal mode number n
  m=6,7,8,9,10,12,14,16   ! specified modes (m,n) for diagnosis: poloidal mode number m



  ! Equilibrium profiles: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  jr=1000,       ! number of radial nodes for plasma profiles and Er(?) (used for linear interpolation)
  linterp_from_nodes=1,    ! for getting profiles on radial nodes. =0, use EZ spline (slow); >0, put on radial nodes first and then use linear interpolation     

  plasma=0,        ! profile model for Ti, Te, n, u: 0=analytic profiles, 1=numerical profiles from file (below)
                   ! Ti=ion temperature profile, Te=electron temperature profile, Tn=density profile, u=parallel flow profile
  q2=-2.184d0,     ! defines formula for profiles of Ti, Te, n, omega; 
                   ! if <0 formula is: R0/L = -kappa*exp{-[(r-rc)/rw]**alpha}          ! r is radial coord.
  rc=0.5d0,        ! parameter rc in profile formula of Ti, Te, n, omega (see q2), if 'plasma=0' and q2<0

  ! electron related:
  T_e0=5.0d0,     ! parameter T_e0  for Te profile (see q2), if 'plasma=0'
  kappate=6.92d0,   ! parameter kappa for Te profile formula (see q2), if 'plasma=0' and q2<0
  flow0=6.0d0,     ! parameter alpha for Te profile formula (see q2), if 'plasma=0' and q2<0
  flow2=0.20d0,    ! parameter rw in Te profile formula (see q2), if 'plasma=0' and q2<0

  ! 1st ion species related:
  !distfun_type=0,   ! decides upon ion specie's distribution function (default= 0=Maxwellian)
  T_0    =5.0d0,    ! parameter T_0   for Ti profile (see q2), if 'plasma=0'
  q0     =6.0d0,     ! parameter alpha in Ti profile formula (see q2), if 'plasma=0' and q2<0
  rw     =0.20d0,    ! parameter rw in Ti profile formula (see q2), if 'plasma=0' and q2<0
  kappati=6.92d0,     ! parameter kappa for Ti profile formula (see q2), if 'plasma=0' and q2<0
  n_0    =1.6d13,    ! parameter n_0 for n profile (see q2), if 'plasma=0'
  kappan =2.2222d0,     ! parameter kappa for n  profile formula (see q2), if 'plasma=0' and q2<0 (and q>0?)
  temconst=0.2d0,   ! array-parameter rw for ni profile formula (see q2), if 'plasma==0' and q2<0
  edensity0=6.0d0,   ! array-parameter alpha for density formula (see q2), if 'plasma==0' and q2<0

  u_0=0.0d0,     ! parameter u_0 for u profile (see q2), if 'plasma=0'
  du=4.0d0,      ! parameter du  for u profile (see q2), if 'plasma=0'
  xu=6.0d0,      ! parameter xu  for u profile (see q2), if 'plasma=0'
  Au=0.0d0,      ! parameter Au  for u profile (see q2), if 'plasma=0'
  Bu=0.28d0,     ! parameter Bu  for u profile (see q2), if 'plasma=0'


  ! Equilibrium radial electric field: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Er_from_GTC_Neo=0,        ! include equilibrium radial electrif field? 0=no equilibrium Er; >0 = include equilibrium Er
  n_node=10,                ! number of radial nodes for Er data

  ! Analytic profile for Magnetic Perturbation 
  ! alpha_amplitude(r) = r*(ga0+dgadt*t)*exp(-abs(((r-rc)/rw))**pp)
  ! alpha_phase(theta,zeta) = (gan*zeta - gam*theta - gagw*t + gagf0)
  ! alpha(r,theta,zeta) = alpha_amplitude(r) * cos(alpha_phase)
  electromagnetic=1         ! include electromagnetic perturbation?
  Streaming_Flag=1          ! (include=1) or (exclude=0) perturbed streaming term and mirror force
  num_of_EMP=2
  ga0=1e-4, 1e-4                ! Initial magnitude of alpha at r=rc
  ga_rc=2*0.5,                ! rc of alpha
  ga_rw=2*0.3,                ! rw of alpha
  ga_pp=2*6.0,                  ! pp of alpha
  gam=3,5                    ! m mode number of alpha
  gan=2,4                   ! n mode number of alpha
  gagw=2*0.00000,             ! frequency (omega) of alpha, in the code normalization unit (1/Omega_ci)
  gagf0=2*3.14,               ! Initial phase of alpha
  dga0dt=2*0.0		              ! Linear growth rate of ga0, [in unit L_T/v_th]
/

&read2 ! related to a "diagnostic box" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  gama0=0.21d0,             ! one of two coefficients used in solving Poisson iteratively
  its=5,                    ! number of iteration used for superLU interation in poisson solver
  its_jacobi=1,             ! used to be for Jacobi iteration in poisson solver 
                            ! -> used control parameter for charge density conservation 05/12
  nelectron=1,              ! number of sub-step for pushing electrons; must take 1 or an even number due to 2nd RK used
  nv=16,                    ! number of grids in vz and vper 
  delta_v=0.2d0,                 ! grid size in velocity space 
  dradius=0.025d0,               ! box size in radial direction
  dtheta=0.02d0,                 ! box size in theta direction
  rprobe=0.50d0, 0.50d0, 0.50d0  ! radial locations of probes
  tprobe=0.00d0, 0.50d0, 1.00d0  ! poloidal locations of probes
  gama1=0.09d0                   ! one of two coefficients used in solving Poisson iteratively   
/

&var_control ! Lu ZhiXin's extra parameters: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! physical effects/model !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ifneoclassical=0,   ! switch on neoclassical effect if=1 (?)
  ratio_Er401=0.0d0,  ! modify equilibrium field: multiply Er401 by thisfactor
  ratio_Er_grid=0.d0, ! modify equilibrium field: multiply Er_grid by thisfactor
  flip_rotation=0.0d0,

  ! control radial grid: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lmin=-0.02d0,       ! buffer regime from lmin to lmax at inner boundary
  lmax=-0.04d0,       ! buffer regime from lmin to lmax at inner boundary
  rmin=-1.0d0,        ! buffer regime from rmin to rmax at outer boundary
  rmax=-1.0d0,        ! buffer regime from rmin to rmax at outer boundary
                      ! -1 sets buffer regime such that rmin=psi_grid(mpsi-nbound), rmax=psi_grid(mpsi-nbound/2)
                      ! otherwise, one can set e.g. if a1=0.9: rmax=0.9, rmin=0.87
  ! control poloidal grid: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  iratio_dthe=10,     ! number of closest surfaces to the edge, for which theta grid size is affected by ratio_dthe:
                      ! determine ploidal grid size: the length between ploidal grids is:
  ratio_dthe0=1.50d0, ! in inner (bulk) regime = ratio_dthe1*rho_i
  ratio_dthe1=3.00d0, ! in outer (edge) regime = ratio_dthe0*rho_i
                      ! determine poloidal grid size: length between poloidal grids is linear from surface 0 to edge:
  a0_dthe=0.60d0,     ! linear coefficient to control inner grid
  a1_dthe=0.40d0,     ! linear coefficient to control outer grid

  ! profiles: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  flag_Lc=2,          ! characteristic length: defined by 1=density, 2=ion temperature, 3=electron temperature, 4=flow profile
  a_ref_in=0.50d0,    ! reference radius

  ! particles/markers: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  w_initial_in=1.0d-4,! initial weight of nonlinear run; this value will be multiplied by 1e-16 if linear run

  ! output: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  flag_DEN=1,         ! records density in DEN.* if =1
  flag_FLOW=1,        ! records flow in FLOW.* if =1
  flag_PHI=1,         ! records el. perturbation in PHI.* if =1
  flag_RESTART=1,     ! producese RESTART.* if =1
  flag_GRID=1,        ! records grid info in GRID_COORD.* if =1

  irewrite=100,       ! number of steps to write restart file (?)

  ain=0.20d0,         ! gives start of interval for the "work array" 'worka' to be printed out (??)
  aout=0.80d0         ! gives end of interval for the "work array" 'worka' to be printed out (??)
/
