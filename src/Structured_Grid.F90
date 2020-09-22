module Structured_Grid
   use precision
   use particle_array
   use global_parameters
   implicit none

   type Grid_Info

      integer :: gq_jump=1,s_jump=1
      integer :: dim_s,dim_a,dim_gq
      real(wp) zetamin,zetamax,dgf,inv_dgf ! zeta value of this plane
      real(wp) dgq,inv_dgq ! dgq is constant since it is uniform along theta direction
      real(wp) ds,inv_ds   ! ds is constant since this is the structured grid on (s-gq-gf) coordinate.

      real(wp),allocatable,dimension(:) :: a_1D,ori_sid,s_1D,dsda_1D,gq_1D
      real(wp),allocatable,dimension(:) :: safety_q_1D,n0_1D,Te0_1D,Te0_eV_1D,Ti0_1D ! some equilibrium 1D profile information
      real(wp),allocatable,dimension(:,:) :: a_2D,s_2D,dsda_2D,gq_2D
      real(wp),allocatable,dimension(:,:) :: Cell_dV,inv_Cell_dV ! Cell Volume and its inverse
      real(wp),allocatable,dimension(:,:) :: Cell_dS_s,inv_Cell_dS_s,Cell_dS_gq,inv_Cell_dS_gq,Cell_dS_gf,inv_Cell_dS_gf ! Node_Surface in each direction

      ! Magnetic Field Information (ESI variables)
      real(wp),allocatable,dimension(:,:) :: F,dFda,dgFda,dgYda,R,dRda,dRdgq,Z,dZda,dZdgq,BB,dBBda,dBBdgq,P,T
      real(wp),allocatable,dimension(:,:) :: Jacobian,esi_N,esi_K,g_aa_cov,g_aq_cov,g_qq_cov,g_ff_cov
      real(wp),allocatable,dimension(:,:) :: D,g_aa_cntv,g_aq_cntv,g_qq_cntv,g_ff_cntv
      real(wp),allocatable,dimension(:,:,:) :: B0_cov3,B0_cntv3,j0_cov3,j0_cntv3

      ! Metrics and some vectors regarding s coordinates
      real(wp),allocatable,dimension(:,:) :: g_ss_cntv,g_sq_cntv,g_ss_cov,g_sq_cov,Jacob_sqf
      real(wp),allocatable,dimension(:,:,:) :: B0_cntv3_sqf

      ! Gyro-Orbit positioning information
      ! real(wp),allocatable,dimension(:,:,:,:) :: gyro_over_sqrt_mu ! Dimension: (3,larmor_num,dim_gq,dim_s); First dimension is for (s,gq,gf)

      ! Plasma information
      real(wp),allocatable, dimension(:,:,:) :: dni,dTi ! dimension: (dim_gq, dim_s, 0:1)
      real(wp),allocatable, dimension(:,:,:) :: dne,dTe
      real(wp),allocatable, dimension(:,:) :: dni_SM, dne_SM ! Smoothed Densities for Poisson3D equation and updating Polarization Shielding Components
      real(wp),allocatable, dimension(:,:,:) :: ne,ni ! These are total electron density and ion density including equilibrium density (So, f=f0+df)

      real(wp), allocatable, dimension(:,:,:) :: phi ! dimension (dim_gq, dim_s,0:2) including neighboring left and right poloidal planes
      real(wp), allocatable, dimension(:,:,:,:) :: grad_phi_cov
      ! Test 2D phi calculated by solving 3D Poission equation
      real(wp), allocatable, dimension(:,:) :: test_phi_2D, test_src_2D, test_res_2D
      real(wp), allocatable, dimension(:,:,:) :: test_phi_3D, test_src_3D, test_res_3D

   end type

   type(Grid_Info) :: SG ! Structure Grid

   ! Gyro-Orbit positioning information
   real(wp),allocatable,dimension(:,:,:,:) :: gyro_over_sqrt_mu ! Dimension: (3,larmor_num,dim_gq,dim_s); First dimension is for (s,gq,gf)
   !$acc declare create(gyro_over_sqrt_mu)

contains

   subroutine Initialize_Structured_Grid
      use field_array
      use esi_interface
      implicit none

      integer i,j,k,sid,gqid,ij,Ns,Nq
      real(wp), external :: dsda
      real(wp) tmp, this_ds
      integer larmor,ss,qq
      real(wp) gyro_phase,coeff,dR_du(3,2)

      SG%dim_s=(mpsi+1)/SG%s_jump
      SG%dim_a=SG%dim_s
      SG%dim_gq=mtheta(ipsi)/SG%gq_jump+1

      SG%dgq=2.0*pi/(SG%dim_gq-1)
      SG%inv_dgq=1.0/SG%dgq

      SG%ds=((s1-s0)/real(mpsi))*SG%s_jump
      SG%inv_ds=1.0/SG%ds

      SG%zetamax=zetamax
      SG%zetamin=zetamin
      SG%dgf = zetamax-zetamin
      SG%inv_dgf = 1.0/SG%dgf

      Ns=SG%dim_s
      Nq=SG%dim_gq

      allocate(SG%gq_1D(Nq))
      allocate(SG%a_1D(Ns),SG%ori_sid(Ns),SG%s_1D(Ns),SG%n0_1D(Ns),SG%Te0_1D(Ns),SG%Te0_eV_1D(Ns),SG%Ti0_1D(Ns),SG%dsda_1D(Ns),SG%safety_q_1D(Ns))
      allocate(SG%a_2D(Nq,Ns),SG%dsda_2D(Nq,Ns),SG%gq_2D(Nq,Ns),SG%Cell_dV(Nq,Ns),SG%inv_Cell_dV(Nq,Ns))
      allocate(SG%Cell_dS_s(Nq,Ns),SG%inv_Cell_dS_s(Nq,Ns),SG%Cell_dS_gq(Nq,Ns),SG%inv_Cell_dS_gq(Nq,Ns),SG%Cell_dS_gf(Nq,Ns),SG%inv_Cell_dS_gf(Nq,Ns))

      k=1
      do i=0,mpsi,SG%s_jump
         SG%ori_sid(k)=i
         SG%a_1D(k)=psi_grid(i)
         SG%a_2D(:,k)=psi_grid(i)
         SG%s_1D(k)=s_grid(i)
         SG%n0_1D(k)=den_grid(i)
         SG%Te0_1D(k)=Te_grid(i)
         SG%Te0_eV_1D(k)=1e3*Tem_00*Te_grid(i)
         SG%Ti0_1D(k)=vth_grid(i)**2.0
         SG%safety_q_1D(k)=q_grid(i)
         k=k+1
      enddo

      do i=1,Nq
         SG%gq_1D(i)=(2.0*pi)/(Nq-1)*(i-1)
         SG%gq_2D(i,:)=SG%gq_1D(i)
      enddo

      do sid=1,SG%dim_s
         SG%dsda_1D(sid)=dsda(SG%a_1D(sid))
         SG%dsda_2D(:,sid)=SG%dsda_1D(sid)
      enddo

      allocate(SG%Jacobian(Nq,Ns),SG%F(Nq,Ns),SG%dFda(Nq,Ns),SG%dgFda(Nq,Ns),SG%dgYda(Nq,Ns),SG%R(Nq,Ns),SG%dRda(Nq,Ns),SG%dRdgq(Nq,Ns),SG%Z(Nq,Ns),SG%dZda(Nq,Ns),SG%dZdgq(Nq,Ns),SG%BB(Nq,Ns),SG%dBBda(Nq,Ns),SG%dBBdgq(Nq,Ns),SG%P(Nq,Ns),SG%T(Nq,Ns))
      allocate(SG%esi_N(Nq,Ns),SG%esi_K(Nq,Ns),SG%g_aa_cov(Nq,Ns),SG%g_aq_cov(Nq,Ns),SG%g_qq_cov(Nq,Ns),SG%g_ff_cov(Nq,Ns))
      allocate(SG%D(Nq,Ns),SG%g_aa_cntv(Nq,Ns),SG%g_aq_cntv(Nq,Ns),SG%g_qq_cntv(Nq,Ns),SG%g_ff_cntv(Nq,Ns))
      allocate(SG%B0_cov3(3,Nq,Ns))
      allocate(SG%B0_cntv3(3,Nq,Ns),SG%j0_cov3(3,Nq,Ns),SG%j0_cntv3(3,Nq,Ns))

      allocate(SG%g_ss_cntv(Nq,Ns),SG%g_sq_cntv(Nq,Ns),SG%g_ss_cov(Nq,Ns),SG%g_sq_cov(Nq,Ns),SG%Jacob_sqf(Nq,Ns))
      allocate(SG%B0_cntv3_sqf(3,Nq,Ns))

      allocate(SG%dne(Nq,Ns,0:1),SG%dni(Nq,Ns,0:1))
      allocate(SG%dne_SM(Nq,Ns),SG%dni_SM(Nq,Ns))
      allocate(SG%dTe(Nq,Ns,0:1),SG%dTi(Nq,Ns,0:1))
      allocate(SG%ne(Nq,Ns,0:1),SG%ni(Nq,Ns,0:1))

      call esiget2dthese_more(SG%F,SG%dFda,SG%dgFda,SG%dgYda,SG%R,SG%dRda,SG%dRdgq,SG%Z,SG%dZda,SG%dZdgq,SG%BB,SG%dBBda,SG%dBBdgq,SG%P,SG%T,SG%a_2D,SG%gq_2D,Nq*Ns)

      SG%Jacobian=SG%R*(SG%dZda*SG%dRdgq-SG%dZdgq*SG%dRda)
      SG%esi_N=(SG%dRda*SG%dRdgq+SG%dZda*SG%dZdgq)/SG%Jacobian
      SG%esi_K=(SG%dRdgq*SG%dRdgq+SG%dZdgq*SG%dZdgq)/SG%Jacobian

      SG%B0_cov3(1,:,:)=-SG%esi_N*SG%dgYda
      SG%B0_cov3(2,:,:)=-SG%esi_K*SG%dgYda
      SG%B0_cov3(3,:,:)=SG%F

      SG%B0_cntv3(1,:,:)=0.0
      SG%B0_cntv3(2,:,:)=-SG%dgYda/SG%Jacobian
      SG%B0_cntv3(3,:,:)=SG%F/(SG%R**2)

      SG%j0_cntv3(1,:,:)=0.0
      SG%j0_cntv3(2,:,:)=-SG%dFda/SG%Jacobian
      SG%j0_cntv3(3,:,:)=SG%P+SG%T/(SG%R**2)

      SG%g_aa_cov=SG%dRda*SG%dRda+SG%dZda*SG%dZda
      SG%g_aq_cov=SG%dRda*SG%dRdgq+SG%dZda*SG%dZdgq
      SG%g_qq_cov=SG%dRdgq*SG%dRdgq+SG%dZdgq*SG%dZdgq
      SG%g_ff_cov=SG%R**2

      SG%D=SG%dZda*SG%dRdgq-SG%dZdgq*SG%dRda

      SG%g_aa_cntv=(SG%dRdgq**2+SG%dZdgq**2)/SG%D**2   ! SG%grada^2
      SG%g_qq_cntv=(SG%dZda**2+SG%dRda**2)/SG%D**2   ! SG%gradtheta^2
      SG%g_ff_cntv=1.0/SG%R**2                 ! SG%gradphi^2
      SG%g_aq_cntv=(-SG%dRda*SG%dRdgq-SG%dZda*SG%dZdgq)/SG%D**2  ! SG%gradagradtheta

      ! Metrics regarding s coordinate
      do sid=1,SG%dim_s
         tmp=dsda(SG%a_1D(sid))
         SG%g_ss_cov(:,sid)=SG%g_aa_cov(:,sid)/(tmp**2.0)
         SG%g_sq_cov(:,sid)=SG%g_aq_cov(:,sid)/tmp

         SG%g_ss_cntv(:,sid)=SG%g_aa_cntv(:,sid)*(tmp**2.0)
         SG%g_sq_cntv(:,sid)=SG%g_aq_cntv(:,sid)*tmp

         SG%Jacob_sqf(:,sid)=SG%Jacobian(:,sid)/tmp

         SG%B0_cntv3_sqf(1,:,sid)=SG%B0_cntv3(1,:,sid)*tmp
      enddo
      SG%B0_cntv3_sqf(2:3,:,:)=SG%B0_cntv3(2:3,:,:)

      SG%j0_cov3(1,:,:)=SG%g_aa_cov*SG%j0_cntv3(1,:,:)+SG%g_aq_cov*SG%j0_cntv3(2,:,:)
      SG%j0_cov3(2,:,:)=SG%g_aq_cov*SG%j0_cntv3(1,:,:)+SG%g_qq_cov*SG%j0_cntv3(2,:,:)
      SG%j0_cov3(3,:,:)=SG%g_ff_cov*SG%j0_cntv3(3,:,:)



      ! Calculate gyro-radius information
      ! gyro-orbit position = GC_position + (gyro_over_sqrt_mu)*sqrt(mu)
      allocate(gyro_over_sqrt_mu(3,larmor_num,SG%dim_gq,SG%dim_s))
      do ss=1,SG%dim_s
         do qq=1,SG%dim_gq
            do larmor=1,larmor_num
               gyro_phase = (larmor-1)*2.0*pi*inv_larmor_num

               ! Note that dgq_du2 does not need (1/cos_alpha^2) correction unlike the original GTS.
               ! This is because, in the case of Structured Grid (SG),
               ! we first calculate the exact GyroOrbit positions at the Gyro-Plane using gyro_over_sqrt_mu
               ! and then project them to the nearby poloidal planes along the magnetic field, which is done at Measure_Ions_on_SG
               dR_du(1,1) = sqrt(SG%g_ss_cntv(qq,ss)) !ds_du1
               dR_du(1,2) = 0.0 ! ds_du2
               dR_du(2,1) = SG%g_sq_cntv(qq,ss)/sqrt(SG%g_ss_cntv(qq,ss)) ! dgq_du1
               dR_du(2,2) = SG%B0_cov3(3,qq,ss)/(SG%Jacob_sqf(qq,ss)*SG%BB(qq,ss)*sqrt(SG%g_ss_cntv(qq,ss))) ! dgq_du2
               dR_du(3,1) = 0.0 ! dgf_du1
               dR_du(3,2) = -SG%B0_cov3(2,qq,ss)/(SG%Jacob_sqf(qq,ss)*SG%BB(qq,ss)*sqrt(SG%g_ss_cntv(qq,ss))) ! dgf_du2

               coeff = sqrt(2.0/SG%BB(qq,ss))
               do i=1,3
                  gyro_over_sqrt_mu(i,larmor,qq,ss) = coeff*(cos(gyro_phase)*dR_du(i,1) + sin(gyro_phase)*dR_du(i,2))
               enddo
            enddo
         enddo
      enddo

      if(Flag_DK_Ion==1) then
         gyro_over_sqrt_mu=0.0
      endif

      !$acc update device(gyro_over_sqrt_mu)

      ! NOTE : Get Cell volume like dV=(Jacobian_sqf)*ds*dtheta*dgf (Calculate on (s,gq,gf) grid)
      !        Following Cell_dV is an approximation
      !        For higher accuracy, a Weighted Volume by considering helical particle shape and varying da is required

      do sid=1,Ns
         if(sid==1 .or. sid==Ns) then
            this_ds=0.5*SG%ds
         else
            this_ds=SG%ds
         endif

         do gqid=1,Nq
            SG%Cell_dV(gqid,sid)=SG%Jacob_sqf(gqid,sid)*this_ds*SG%dgq*SG%dgf
            SG%Cell_dS_s(gqid,sid)=sqrt(SG%g_aa_cntv(gqid,sid))*SG%Jacob_sqf(gqid,sid)*SG%dgq*SG%dgf
            SG%Cell_dS_gq(gqid,sid)=sqrt(SG%g_qq_cntv(gqid,sid))*SG%Jacob_sqf(gqid,sid)*this_ds*SG%dgf
            SG%Cell_dS_gf(gqid,sid)=sqrt(SG%g_ff_cntv(gqid,sid))*SG%Jacob_sqf(gqid,sid)*this_ds*SG%dgq
         enddo
      enddo
      SG%inv_Cell_dV=1.0/SG%Cell_dV
      SG%inv_Cell_dS_s=1.0/SG%Cell_dS_s
      SG%inv_Cell_dS_gq=1.0/SG%Cell_dS_gq
      SG%inv_Cell_dS_gf=1.0/SG%Cell_dS_gf


      ! allocate phi
      allocate(SG%phi(Nq,Ns,0:2)) ! 0: left plane (@zetamin), 1: current plane (@zetamax), 2: right plane (@zeta=zetamax+delta_zeta)
      SG%phi=0.0
      allocate(SG%grad_phi_cov(3,Nq,Ns,0:1))

      allocate(SG%test_phi_2D(Nq,Ns),SG%test_src_2D(Nq,Ns),SG%test_res_2D(Nq,Ns))
      allocate(SG%test_phi_3D(Nq,Ns,mzetamax),SG%test_src_3D(Nq,Ns,mzetamax),SG%test_res_3D(Nq,Ns,mzetamax))


   end subroutine Initialize_Structured_Grid

   subroutine Measure_Electron_Density_on_SG(flag_which_target)
      use precision
      use global_parameters,only:me
      use particle_array,only:zeon,zeon0
      use particle_decomp
      use function_gg_module
      use mpi_defs

      implicit none

      integer, intent(in), optional :: flag_which_target
      integer :: flagT ! flag which target will be measured
      real(wp) ratio(4), ws0,ws1,wgq0,wgq1,wgf0,wgf1
      real(wp) P_s, P_gq, P_gf, P_safety_q
      real(wp) target_gq, target_gf
      real(wp) tmp_var(SG%dim_gq,SG%dim_s,0:1)
      integer m, vid, sid, gqid, ss(4), qq(4), k, length, istatus(MPI_STATUS_SIZE),ierror

      real(wp), allocatable, dimension(:,:,:,:) :: eData, tmp_eData
      real(wp), allocatable, dimension(:,:) :: ptl_data
      integer num_eD ! num of eData (First Dimension)
      real(wp) tt0,tt1,tmp


      !$acc declare create(ptl_data)


      ! Determine Flag target data
      ! flagT=0 (Default; only measures a term from particle weight zeon(5,m))
      ! flagT=1 (Diagnostic puprpose; measures various physical quantities)
      flagT=0 ! Default value
      if(present(flag_which_target)) flagT = flag_which_target


#ifdef _MGY_TIMER_
      tt0=MPI_WTIME()
#endif
      if(acc_newer_H_zeon) then
         !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
         acc_newer_H_zeon=.false.
         acc_newer_D_zeon=.false.
      endif

      select case(flagT)
         case(0) ! Default; only zeon(5,m)
            num_eD=1
         case(1) ! Diagnostic purpose
            num_eD=2
         case default
            write(stdout,*) "Error!! flagT=",flagT," in Measure_Electron_Density_on_SG"
            call flush(stdout)
            STOP
      end select


      allocate(eData(num_eD,SG%dim_gq,SG%dim_s,0:1),tmp_eData(num_eD,SG%dim_gq,SG%dim_s,0:1))
      allocate(ptl_data(num_eD,me))
      eData = 0.0



      select case(flagT)
         case(0)
            !$acc parallel loop present(zeon)
            do m=1,me
               ptl_data(1,m) = zeon(5,m)
            enddo
         case(1)
            !$acc parallel loop present(zeon,zeon0,ptl_data,backTe)
            do m=1,me
               ptl_data(1,m) = zeon(5,m)
               ptl_data(2,m) = zeon0(8,m)
            enddo
      end select


      ! Mulitplying coeff in the case of Flag_Load_and_WeightDef==0
      if(Flag_Load_and_WeightDef==0) then
         !$acc parallel loop present(ptl_data)
         do m=1,me
            ptl_data(:,m) = eNrNm*ptl_data(:,m)
         enddo
      endif

      !$acc parallel loop gang vector present(zeon,zeon0,ptl_data) private(ratio,qq,ss)
      do m=1,me
         P_s = zeon(7,m)
         P_gq = modulo(zeon(2,m),2*pi)
         P_gf = modulo(zeon(3,m),2*pi)

         wgf1 = abs(P_gf-SG%zetamin)*SG%inv_dgf
         wgf0 = (1.0-wgf1)

         sid = max(1,min(int((P_s-s0)*SG%inv_ds) +1,SG%dim_s-1)) ! +1 is due to Fortran's default index starting from 1
         ws1 = (P_s-SG%s_1D(sid))*SG%inv_ds
         ws0 = (1.0-ws1)

         P_safety_q = ws0*SG%safety_q_1D(sid) + ws1*SG%safety_q_1D(sid+1)

         ! Scatter particle to zetamin plane first
         target_gf = zetamin
         target_gq = modulo(P_gq + (target_gf-P_gf)/P_safety_q, 2.0*pi)
         gqid  = target_gq*SG%inv_dgq+1
         wgq1 = (target_gq-SG%gq_1D(gqid))*SG%inv_dgq
         wgq0 = 1.0-wgq1

         ratio(1)=wgf0*ws0*wgq0; qq(1)=gqid;    ss(1)=sid;
         ratio(2)=wgf0*ws0*wgq1; qq(2)=gqid+1;  ss(2)=sid;
         ratio(3)=wgf0*ws1*wgq0; qq(3)=gqid;    ss(3)=sid+1;
         ratio(4)=wgf0*ws1*wgq1; qq(4)=gqid+1;  ss(4)=sid+1;
         !$acc loop seq
         do k=1,4
            do vid=1,num_eD
               !$acc atomic update
               eData(vid,qq(k),ss(k),0) = eData(vid,qq(k),ss(k),0) + ratio(k)*ptl_data(vid,m)
            enddo
         enddo


         ! Scatter particle to zetamax plane
         target_gf = zetamax
         target_gq = modulo(P_gq + (target_gf-P_gf)/P_safety_q, 2.0*pi)
         gqid  = target_gq*SG%inv_dgq+1
         wgq1 = (target_gq-SG%gq_1D(gqid))*SG%inv_dgq
         wgq0 = 1.0-wgq1


         ratio(1)=wgf1*ws0*wgq0; qq(1)=gqid;    ss(1)=sid;
         ratio(2)=wgf1*ws0*wgq1; qq(2)=gqid+1;  ss(2)=sid;
         ratio(3)=wgf1*ws1*wgq0; qq(3)=gqid;    ss(3)=sid+1;
         ratio(4)=wgf1*ws1*wgq1; qq(4)=gqid+1;  ss(4)=sid+1;
         !$acc loop seq
         do k=1,4
            do vid=1,num_eD
               !$acc atomic update
               eData(vid,qq(k),ss(k),1) = eData(vid,qq(k),ss(k),1) + ratio(k)*ptl_data(vid,m)
            enddo
         enddo

      enddo


      if(npartdom>1) then
         tmp_eData = eData
         length  = num_eD*SG%dim_gq*SG%dim_s*2
         call MPI_ALLREDUCE(tmp_eData, eData, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      endif

      ! Match poloidal end cell
      ! (gq=0 and gq=2pi is the same position physically
      ! although they have different gqid)
      do sid=1,SG%dim_s
         eData(:,1,sid,:) = eData(:,1,sid,:) + eData(:,SG%dim_gq,sid,:)
         eData(:,SG%dim_gq,sid,:) = eData(:,1,sid,:)
      enddo

      ! toroidal communication
      length = num_eD*SG%dim_gq*SG%dim_s
      call MPI_SENDRECV(eData(:,:,:,0),length,mpi_Rsize,left_pe,1, tmp_eData(:,:,:,1),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
      eData(:,:,:,1) = eData(:,:,:,1) + tmp_eData(:,:,:,1)
      call MPI_SENDRECV(eData(:,:,:,1),length,mpi_Rsize,right_pe,1, eData(:,:,:,0),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

      ! Normalize the results
      do vid=1,num_eD
         eData(vid,:,:,0) = eData(vid,:,:,0)*SG%inv_Cell_dV
         eData(vid,:,:,1) = eData(vid,:,:,1)*SG%inv_Cell_dV
      enddo

      ! Artificial damping near the both inner and outer boundary
      do sid=1,1+nbound_inner
         tmp = SmootherStep(SG%s_1D(sid),SG%s_1D(1),0.0,SG%s_1D(1+nbound_inner),1.0)
         eData(:,:,sid,:) = tmp*eData(:,:,sid,:)
      enddo
      do sid=SG%dim_s-nbound_outer,SG%dim_s
         tmp = SmootherStep(SG%s_1D(sid),SG%s_1D(SG%dim_s-nbound_outer),1.0,SG%s_1D(SG%dim_s),0.0)
         eData(:,:,sid,:) = tmp*eData(:,:,sid,:)
      enddo

      select case(flagT)
         case(0) ! Default; only zeon(5,m)
            SG%dne = eData(1,:,:,:)
         case(1) ! Diagnostic purpose
            SG%dne = eData(1,:,:,:)
            SG%ne  = eData(2,:,:,:)
      end select

      deallocate(eData,tmp_eData)
      deallocate(ptl_data)

#ifdef _MGY_TIMER_
      tt1=MPI_WTIME()
      if(mype==0) write(stdout,*) "Measure_Electron_Density_on_SG takes", tt1-tt0, " secs"
      call flush(stdout)
#endif

   end subroutine Measure_Electron_Density_on_SG

   subroutine Measure_Ion_Density_on_SG(flag_which_target)
      use precision
      use global_parameters,only:mi
      use particle_array,only:zion,zion0
      use particle_decomp
      use function_gg_module
      use mpi_defs

      implicit none

      integer, intent(in), optional :: flag_which_target
      integer :: flagT ! flag which target will be measured
      real(wp) ratio(4), ws0,ws1,wgq0,wgq1,wgf0,wgf1
      real(wp) P_s, P_gq, P_gf, P_safety_q
      real(wp) target_gq, target_gf
      real(wp) tmp_var(SG%dim_gq,SG%dim_s,0:1)
      integer m, k, sid, gqid, vid, ss(4), qq(4), length, istatus(MPI_STATUS_SIZE),ierror, larmor
      real(wp) tt0, tt1, tmp
      real(wp), allocatable, dimension(:,:,:,:) :: iData, tmp_iData
      real(wp), allocatable, dimension(:,:) :: ptl_data
      integer num_iD ! num_of iData (First dimension)
      real(wp) :: intpl_gyro(3,larmor_num)

      !$acc declare create(ptl_data)

      ! Determine Flag target data
      ! flagT=0 (Default; only measures a term from particle weight zeon(5,m))
      ! flagT=1 (Diagnostic puprpose; measures various physical quantities)
      flagT=0 ! Default value
      if(present(flag_which_target)) flagT = flag_which_target


#ifdef _MGY_TIMER_
      tt0=MPI_WTIME()
#endif

      if(acc_newer_H_zion) then
         !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
         acc_newer_D_zion=.false.
         acc_newer_H_zion=.false.
      endif


     select case(flagT)
         case(0) ! Default; only zion(5,m)
            num_iD=1
         case(1) ! Diagnostic purpose
            num_iD=2
         case default
            write(stdout,*) "Error!! flagT=",flagT," in Measure_Ion_Density_on_SG"
            call flush(stdout)
            STOP
      end select


      allocate(iData(num_iD,SG%dim_gq,SG%dim_s,0:1), tmp_iData(num_iD,SG%dim_gq,SG%dim_s,0:1))
      allocate(ptl_data(num_iD,mi))
      iData=0.0


      select case(flagT)
         case(0)
            !$acc parallel loop present(zion)
            do m=1,mi
               ptl_data(1,m) = zion(5,m)
            enddo
         case(1)
            !$acc parallel loop present(zion,zion0,ptl_data)
            do m=1,mi
               ptl_data(1,m) = zion(5,m)
               ptl_data(2,m) = zion0(8,m)
            enddo
      end select


      ! Mulitplying coeff in the case of Flag_Load_and_WeightDef==0
      if(Flag_Load_and_WeightDef==0) then
         !$acc parallel loop present(ptl_data)
         do m=1,mi
            ptl_data(:,m) = iNrNm*ptl_data(:,m)
         enddo
      endif

      !$acc parallel loop gang vector present(zion,zion0,zion_s_on_GyroOrbit,zion_gq_on_GyroOrbit,gyro_over_sqrt_mu,ptl_data) private(intpl_gyro,ratio,qq,ss)
      do m=1,mi

         ! Measure the intpl_gyro at the Guiding Center position
         P_s = zion(7,m)
         P_gq = modulo(zion(2,m),2*pi)
         P_gf = modulo(zion(3,m),2*pi)

         sid = max(1,min(int((P_s-s0)*SG%inv_ds) +1,SG%dim_s-1)) ! +1 is due to Fortran's default index starting from 1
         ws1 = (P_s-SG%s_1D(sid))*SG%inv_ds
         ws0 = (1.0-ws1)
         gqid = P_gq*SG%inv_dgq+1
         wgq1 = (P_gq-SG%gq_1D(gqid))*SG%inv_dgq
         wgq0 = 1.0-wgq1
         !$acc loop seq
         do larmor=1,larmor_num
            intpl_gyro(:,larmor) = ws0*(wgq0*gyro_over_sqrt_mu(:,larmor,gqid,sid) + wgq1*gyro_over_sqrt_mu(:,larmor,gqid+1,sid)) &
                                  +ws1*(wgq0*gyro_over_sqrt_mu(:,larmor,gqid,sid+1) + wgq1*gyro_over_sqrt_mu(:,larmor,gqid+1,sid+1))
         enddo


         ! Using the intpl_gyro, calculate Larmor orbit positions and project them to the nearby poloidal planes
         !$acc loop seq
         do larmor=1,larmor_num

            zion_s_on_GyroOrbit(larmor,m) = zion(7,m) + zion(6,m)*intpl_gyro(1,larmor)
            zion_gq_on_GyroOrbit(larmor,m) =  modulo(zion(2,m) + zion(6,m)*intpl_gyro(2,larmor), 2.0*pi)

            P_s = max(s0,min(s1,zion_s_on_GyroOrbit(larmor,m)))
            P_gq = zion_gq_on_GyroOrbit(larmor,m)
            P_gf = modulo(zion(3,m) + zion(6,m)*intpl_gyro(3,larmor), 2.0*pi)
            P_gf = max(zetamin,min(zetamax,P_gf))

            wgf1 = abs(P_gf-SG%zetamin)*SG%inv_dgf
            wgf0 = (1.0-wgf1)

            sid = max(1,min(int((P_s-s0)*SG%inv_ds) +1,SG%dim_s-1)) ! +1 is due to Fortran's default index starting from 1
            ws1 = (P_s-SG%s_1D(sid))*SG%inv_ds
            ws0 = (1.0-ws1)

            P_safety_q = ws0*SG%safety_q_1D(sid) + ws1*SG%safety_q_1D(sid+1)

            ! Scatter particle to zetamin plane first
            target_gf = zetamin
            target_gq = modulo(P_gq + (target_gf-P_gf)/P_safety_q, 2.0*pi)
            gqid  = target_gq*SG%inv_dgq+1
            wgq1 = (target_gq-SG%gq_1D(gqid))*SG%inv_dgq
            wgq0 = 1.0-wgq1


            ratio(1)=inv_larmor_num*wgf0*ws0*wgq0; qq(1)=gqid;    ss(1)=sid;
            ratio(2)=inv_larmor_num*wgf0*ws0*wgq1; qq(2)=gqid+1;  ss(2)=sid;
            ratio(3)=inv_larmor_num*wgf0*ws1*wgq0; qq(3)=gqid;    ss(3)=sid+1;
            ratio(4)=inv_larmor_num*wgf0*ws1*wgq1; qq(4)=gqid+1;  ss(4)=sid+1;
            !$acc loop seq
            do k=1,4
               do vid=1,num_iD
                  !$acc atomic update
                  iData(vid,qq(k),ss(k),0) = iData(vid,qq(k),ss(k),0) + ratio(k)*ptl_data(vid,m)
               enddo
            enddo


            ! Scatter particle to zetamax plane
            target_gf = zetamax
            target_gq = modulo(P_gq + (target_gf-P_gf)/P_safety_q, 2.0*pi)
            gqid  = target_gq*SG%inv_dgq+1
            wgq1 = (target_gq-SG%gq_1D(gqid))*SG%inv_dgq
            wgq0 = 1.0-wgq1

            ratio(1)=inv_larmor_num*wgf1*ws0*wgq0; qq(1)=gqid;    ss(1)=sid;
            ratio(2)=inv_larmor_num*wgf1*ws0*wgq1; qq(2)=gqid+1;  ss(2)=sid;
            ratio(3)=inv_larmor_num*wgf1*ws1*wgq0; qq(3)=gqid;    ss(3)=sid+1;
            ratio(4)=inv_larmor_num*wgf1*ws1*wgq1; qq(4)=gqid+1;  ss(4)=sid+1;
            !$acc loop seq
            do k=1,4
               do vid=1,num_iD
                  !$acc atomic update
                  iData(vid,qq(k),ss(k),1) = iData(vid,qq(k),ss(k),1) + ratio(k)*ptl_data(vid,m)
               enddo
            enddo
         enddo
      enddo

      if(npartdom>1) then
         tmp_iData = iData
         length = num_iD*SG%dim_gq*SG%dim_s*2
         call MPI_ALLREDUCE(tmp_iData, iData, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      endif

      ! Match poloidal end cell
      ! (gq=0 and gq=2pi is the same position physically
      ! although they have different gqid)
      do sid=1,SG%dim_s
         iData(:,1,sid,:) = iData(:,1,sid,:) + iData(:,SG%dim_gq,sid,:)
         iData(:,SG%dim_gq,sid,:) = iData(:,1,sid,:)
      enddo

      ! toroidal communication
      length = num_iD*SG%dim_gq*SG%dim_s
      call MPI_SENDRECV(iData(:,:,:,0),length,mpi_Rsize,left_pe,1, tmp_iData(:,:,:,1),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
      iData(:,:,:,1) = iData(:,:,:,1) + tmp_iData(:,:,:,1)
      call MPI_SENDRECV(iData(:,:,:,1),length,mpi_Rsize,right_pe,1, iData(:,:,:,0),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

      ! Normalize the results
      do vid=1,num_iD
         iData(vid,:,:,0) = iData(vid,:,:,0)*SG%inv_Cell_dV
         iData(vid,:,:,1) = iData(vid,:,:,1)*SG%inv_Cell_dV
      enddo


      ! Artificial damping near the both inner and outer boundary
      do sid=1,1+nbound_inner
         tmp = SmootherStep(SG%s_1D(sid),SG%s_1D(1),0.0,SG%s_1D(1+nbound_inner),1.0)
         iData(:,:,sid,:) = tmp*iData(:,:,sid,:)
      enddo
      do sid=SG%dim_s-nbound_outer,SG%dim_s
         tmp = SmootherStep(SG%s_1D(sid),SG%s_1D(SG%dim_s-nbound_outer),1.0,SG%s_1D(SG%dim_s),0.0)
         iData(:,:,sid,:) = tmp*iData(:,:,sid,:)
      enddo



      select case(flagT)
         case(0) ! Default; only zion(5,m)
            SG%dni = iData(1,:,:,:)
         case(1) ! Diagnostic purpose
            SG%dni = iData(1,:,:,:)
            SG%ni  = iData(2,:,:,:)
      end select

      deallocate(iData,tmp_iData)
      deallocate(ptl_data)

#ifdef _MGY_TIMER_
      tt1=MPI_WTIME()
      if(mype==0) write(stdout,*) "Measure_Ion_Density_on_SG takes", tt1-tt0, " secs"
      call flush(stdout)
#endif
   end subroutine Measure_Ion_Density_on_SG

   subroutine Smooth_SG_Data_2D(num_smooth_gq_s, target_Data2D, Reference_Data2D)

      use precision
      use global_parameters
      implicit none

      integer, intent(in) :: num_smooth_gq_s(:)
      real(wp),intent(inout) :: target_Data2D(SG%dim_gq, SG%dim_s)
      real(wp),optional:: Reference_Data2D(SG%dim_gq, SG%dim_s)
      real(wp) :: VV(SG%dim_gq, SG%dim_s)
      integer :: smid, aid, gqid, Ns,Nq, vid,k
      real(wp) :: SM_Data(0:SG%dim_gq+1, 0:SG%dim_s+1)  ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
      real(wp) :: SM_gq_Coeff(-1:1,0:SG%dim_gq+1, 0:SG%dim_s+1), SM_s_Coeff(-1:1,0:SG%dim_gq+1, 0:SG%dim_s+1)

      Ns = SG%dim_s
      Nq = SG%dim_gq

      VV=1.0
      if(present(Reference_Data2D)) VV=Reference_Data2D

      ! Get Smooth Coeffs which conserves Volume Weighted Integration
      call Get_Conserving_Smooth_Coeffs_on_SG(VV, SM_s_Coeff, SM_gq_Coeff)


      if(size(num_smooth_gq_s)==1) then
         SM_Data = 0.0
         SM_Data(1:Nq,1:Ns) = target_Data2D
         do smid=1,num_smooth_gq_s(1)
            ! Set ghost cell periodic in theta direction
            SM_Data(0,:)    = SM_Data(Nq-1,:)
            SM_Data(Nq+1,:) = SM_Data(2,:)
            ! Smooth in theta direction
            SM_Data(1:Nq,:) = SM_gq_Coeff(-1,1:Nq,:)*SM_Data(0:Nq-1,:) + SM_gq_Coeff(0,1:Nq,:)*SM_Data(1:Nq,:) + SM_gq_Coeff(1,1:Nq,:)*SM_Data(2:Nq+1,:)
            ! Smooth in radial direction
            SM_Data(:,1:Ns) = SM_s_Coeff(-1,:,1:Ns)*SM_Data(:,0:Ns-1) + SM_s_Coeff(0,:,1:Ns)*SM_Data(:,1:Ns) + SM_s_Coeff(1,:,1:Ns)*SM_Data(:,2:Ns+1)
         enddo
         target_Data2D = SM_Data(1:Nq,1:Ns)
      elseif(size(num_smooth_gq_s)==2) then
         SM_Data = 0.0
         SM_Data(1:Nq,1:Ns) = target_Data2D
         do smid=1,num_smooth_gq_s(1)
            ! Set ghost cell periodic in theta direction
            SM_Data(0,:)    = SM_Data(Nq-1,:)
            SM_Data(Nq+1,:) = SM_Data(2,:)
            ! Smooth in theta direction
            SM_Data(1:Nq,:) = SM_gq_Coeff(-1,1:Nq,:)*SM_Data(0:Nq-1,:) + SM_gq_Coeff(0,1:Nq,:)*SM_Data(1:Nq,:) + SM_gq_Coeff(1,1:Nq,:)*SM_Data(2:Nq+1,:)
         enddo
         do smid=1,num_smooth_gq_s(2)
            ! Smooth in radial direction
            SM_Data(:,1:Ns) = SM_s_Coeff(-1,:,1:Ns)*SM_Data(:,0:Ns-1) + SM_s_Coeff(0,:,1:Ns)*SM_Data(:,1:Ns) + SM_s_Coeff(1,:,1:Ns)*SM_Data(:,2:Ns+1)
         enddo
         target_Data2D = SM_Data(1:Nq,1:Ns)
      else
         write(stdout,*) "Error! Check num_smooth_gq_s = ",num_smooth_gq_s(:)
         stop
      endif

   end subroutine Smooth_SG_Data_2D


   subroutine Get_Conserving_Smooth_Coeffs_on_SG(in_VOLorAREA, SM_s_Coeff, SM_gq_Coeff)

      use precision
      use global_parameters
      use field_array
      implicit none

      integer :: Ns,Nq
      real(wp), intent(in) :: in_VOLorAREA(SG%dim_gq, SG%dim_s) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
      real(wp) :: VV(0:SG%dim_gq+1, 0:SG%dim_s+1) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
      real(wp), intent(out):: SM_gq_Coeff(-1:1,0:SG%dim_gq+1, 0:SG%dim_s+1), SM_s_Coeff(-1:1,0:SG%dim_gq+1, 0:SG%dim_s+1)

      Ns = SG%dim_s
      Nq = SG%dim_gq

      VV = 0
      VV(1:Nq,1:Ns) = in_VOLorAREA

      ! Note that theta direction is periodic
      VV(0,:) = VV(Nq-1,:)
      VV(Nq+1,:) = VV(2,:)

      SM_gq_Coeff(-1,1:Nq,:) = 0.5*VV(0:Nq-1,:)/(VV(0:Nq-1,:) + VV(1:Nq,:)) ! Coeff_L = 0.5*VL/(VL+VC)
      SM_gq_Coeff( 0,1:Nq,:) = 0.5*VV(1:Nq,:)*(VV(0:Nq-1,:)+2.0*VV(1:Nq,:)+VV(2:Nq+1,:))/((VV(0:Nq-1,:)+VV(1:Nq,:))*(VV(1:Nq,:)+VV(2:Nq+1,:))) ! 0.5*VC*(VL+2*VC+VR)/((VL+VC)*(VC+VR))
      SM_gq_Coeff( 1,1:Nq,:) = 0.5*VV(2:Nq+1,:)/(VV(2:Nq+1,:) + VV(1:Nq,:)) ! Coeff_R = 0.5*VR/(VR+VC)

      SM_s_Coeff(-1,:,1:Ns) = 0.5*VV(:,0:Ns-1)/(VV(:,0:Ns-1) + VV(:,1:Ns)) ! Coeff_L = 0.5*VL/(VL+VC)
      SM_s_Coeff( 0,:,1:Ns) = 0.5*VV(:,1:Ns)*(VV(:,0:Ns-1)+2.0*VV(:,1:Ns)+VV(:,2:Ns+1))/((VV(:,0:Ns-1)+VV(:,1:Ns))*(VV(:,1:Ns)+VV(:,2:Ns+1))) ! 0.5*VC*(VL+2*VC+VR)/((VL+VC)*(VC+VR))
      SM_s_Coeff( 1,:,1:Ns) = 0.5*VV(:,2:Ns+1)/(VV(:,2:Ns+1) + VV(:,1:Ns)) ! Coeff_R = 0.5*VR/(VR+VC)

      ! To remove NsN values at the raidal ghost cells
      SM_gq_Coeff(:,:,0)   = 0
      SM_gq_Coeff(:,:,Ns+1)= 0
      SM_s_Coeff(:,:,0)    = 0
      SM_s_Coeff(:,:,Ns+1) = 0

   end subroutine Get_Conserving_Smooth_Coeffs_on_SG

   subroutine Interpolate_OriGridData_into_SG(oriData_1D, out_IntpData_2D)

      use precision
      use global_parameters
      use field_array
      implicit none

      real(wp), dimension(mgrid), intent(in) :: oriData_1D
      real(wp), dimension(SG%dim_gq,SG%dim_s), intent(out) :: out_IntpData_2D

      integer i,j,ij,jt, ori_sid

      real(wp) ori_inv_dgq(0:mpsi), gqdum, wt

      ori_inv_dgq=1.0/deltat

      do i=1,SG%dim_s
         do j=1,SG%dim_gq

            ori_sid = SG%ori_sid(i)

            gqdum = modulo(SG%gq_1D(j) - qtinv(ori_sid)*zetamax, 2.0*pi) * ori_inv_dgq(ori_sid)
            jt=int(gqdum);  wt=gqdum-jt
            ij=igrid(ori_sid)+jt
            out_IntpData_2D(j,i) =(1.0-wt)*oriData_1D(ij) + wt*oriData_1D(ij+1)
         enddo
      enddo

      ! Periodic in theta direction (to make sure they are exactly the same)
      out_IntpData_2D(SG%dim_gq,:) = out_IntpData_2D(1,:)

   end subroutine Interpolate_OriGridData_into_SG

   subroutine Interpolate_SG_into_GTSmgrid(SG_Data_2D, GTS_mgrid_Data1D)
      use precision
      use field_array
      use global_parameters
      use particle_decomp
      implicit none
      real(wp), dimension(SG%dim_gq,SG%dim_s), intent(in) :: SG_Data_2D
      real(wp), dimension(mgrid), intent(out) :: GTS_mgrid_Data1D
      integer i,j,ij,sid, gqid
      real(wp) tmp_gq, tmp_s, ws0, ws1, wq0, wq1

      do i=0,mpsi
         tmp_s  = s_grid(i)
         sid = max(1,min(int((tmp_s-s0)*SG%inv_ds)+1,SG%dim_s-1))  ! +1 is due to Fortran's default index starting from 1
         ws1 = (tmp_s-SG%s_1D(sid))*SG%inv_ds
         ws0 = 1.0 - ws1
         do j=0,mtheta(i)
            ij = igrid(i)+j


            tmp_gq = modulo(j*deltat(i) + qtinv(i)*zetamax,2.0*pi)
            gqid= tmp_gq*SG%inv_dgq +1  !+1 is due to Fortran's default index starting from 1
            wq1 = (tmp_gq-SG%gq_1D(gqid))*SG%inv_dgq
            wq0 = 1.0 - wq1

            GTS_mgrid_Data1D(ij) = ws0*(wq0*SG_Data_2D(gqid,sid)+wq1*SG_Data_2D(gqid+1,sid)) + ws1*(wq0*SG_Data_2D(gqid,sid+1)+wq1*SG_Data_2D(gqid+1,sid+1))
         enddo
         ! Periodic in theta direction (to make sure they are exactly the same)
         GTS_mgrid_Data1D(igrid(i)) = GTS_mgrid_Data1D(igrid(i)+mtheta(i))
      enddo

   end subroutine Interpolate_SG_into_GTSmgrid


   subroutine Cal_gradPhi_on_SG
      use precision
      use global_parameters
      use field_array
      use mpi_defs
      use particle_decomp
      implicit none

      integer sid,gqid,istatus(MPI_STATUS_SIZE),ierr

      ! Calculate grad_phi_cov at zetamax (so, the last index is 1)
      ! radial direction
      SG%grad_phi_cov(1,:,1,1) = SG%inv_ds*(SG%phi(:,2,1)-SG%phi(:,1,1))
      do sid=2,SG%dim_s-1
         SG%grad_phi_cov(1,:,sid,1) = 0.5*SG%inv_ds*(SG%phi(:,sid+1,1)-SG%phi(:,sid-1,1))
      enddo
      SG%grad_phi_cov(1,:,SG%dim_s,1) = SG%inv_ds*(SG%phi(:,SG%dim_s,1)-SG%phi(:,SG%dim_s-1,1))

      ! theta direction with the periodic condition
      SG%grad_phi_cov(2,1,:,1) = SG%inv_dgq*(SG%phi(2,:,1)-SG%phi(SG%dim_gq-1,:,1))
      do gqid=2,SG%dim_gq-1
         SG%grad_phi_cov(2,gqid,:,1) = 0.5*SG%inv_dgq*(SG%phi(gqid+1,:,1)-SG%phi(gqid-1,:,1))
      enddo
      SG%grad_phi_cov(2,SG%dim_gq,:,1) = SG%grad_phi_cov(2,1,:,1)

      ! toroidal direction
      SG%grad_phi_cov(3,:,:,1) = 0.5*SG%inv_dgf*(SG%phi(:,:,2)-SG%phi(:,:,0))


      ! Communicate in toroidal direction so that get grad_phi_cov at zetamin (grad_phi_cov(:,:,:,0)) at from the left poloidal plane
      call MPI_SENDRECV(SG%grad_phi_cov(:,:,:,1),3*SG%dim_gq*SG%dim_s,mpi_Rsize,right_pe,0, &
                        SG%grad_phi_cov(:,:,:,0),3*SG%dim_gq*SG%dim_s,mpi_Rsize,left_pe,0,toroidal_comm,istatus,ierr)

   end subroutine Cal_gradPhi_on_SG


end module Structured_Grid
