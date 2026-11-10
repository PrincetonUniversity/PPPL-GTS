&read1  
  eq_data_name='esiP.00_cyclon',

  irun=0,
  mstep=80000,
  msnap=50,
  ndiag=1,
  nonlinear=0.0,
  ifdphi=2,
  mode00=0,
  maintnT=0,

  tstep=0.010,
  micell=50,
  mpsi=40,
  ipsi=20,
  mthetamax=1,
  mzetamax=32,
  npartdom=24,
  ifelectron=2,

  fltr=0.0,
     
  a=1.0,
  a0=0.1,
  a1=0.9,

  q0=6.0,
  q1=0.0,
  q2=-2.184,
  rc=0.5,
  rw=0.3,

  aion=1.0,
  qion=1.0,
  aelectron=1837.0,
  qelectron=-1.0,

  kappati=6.96,
  kappate=6.96,
  kappan=2.2,
  flow0=6.0,
  flow1=0.0,
  flow2=0.3,

  r0=93.2,
  b0=19100.0,
  temconst=0.3,
  edensity0=6.0,

  n_0=1.6d13,
  dn=1.0d0,
  xn=2.0d0,
  An=0.36d0,
  Bn=0.64d0,

  T_0=5.25d0,
  dTemi=6.0d0,
  xt=3.0d0,
  At=0.3d0,
  Bt=0.7d0,

  T_e0=5.25d0,
  dTeme=6.0d0,
  xt_e=3.0d0,
  At_e=0.3d0,
  Bt_e=0.7d0,

  u_0=0.0d0,
  du=4.0d0,
  xu=6.0d0,
  Au=0.0d0,
  Bu=0.28d0,

  stdout=6,
  nbound=5,
  umax=5.0,
  iload=0,

  icollision=0,
  ecollision=0,
  coll_enlarge=50.,
  nss=1,

  n=15,15,15,15,15,15,15,15
  m=6,7,8,9,10,12,14,16

  jr=1000,

  plasma=0,
  linterp_from_nodes=1,
  profile_data_name='prof129533w05-405',

  Er_from_GTC_Neo=0,
  n_node=10,
  Er_data_name='Er_eq.dat',

  movie_3D_h5=0,
  ZF_SEPARATION=1,
  order=1
/
&read2
  gama0=0.,
  its=1,
  its_jacobi=1,
  nelectron=1,
  nv=16,
  delta_v=0.2,
  dradius=0.05,
  dtheta=0.1,
  rprobe=0.50, 0.55, 0.60
  tprobe=0.00, 1.57, 3.14
  gama1=0.0
/
&var_control
  lmin=0.02d0,
  lmax=0.03d0,
  rmin=0.997d0,
  rmax=0.998d0,
  iratio_dthe=30,
  ratio_dthe0=1.0d0,
  ratio_dthe1=1.0d0,
  a0_dthe=0.5d0,
  a1_dthe=0.5d0,
  flag_Lc=2,
  a_ref_in=0.5d0,

  flag_DEN=0,
  flag_FLOW=0,
  flag_PHI=1,
  flag_RESTART=1,
  flag_GRID=1,

  irewrite=2000,
  ifneoclassical=0,
  ratio_Er401=0.d0,
  ratio_Er_grid=0.d0,
  w_initial_in=1.0d-10,
  flip_rotation=0.0d0,
  ain=0.20d0,
  aout=0.80d0,
  pade=1
/
