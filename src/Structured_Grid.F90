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
      real(wp),allocatable,dimension(:) :: safety_q_1D,n0_1D,Te0_1D,Ti0_1D ! some equilibrium 1D profile information
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

      ! Plasma information
      real(wp),allocatable, dimension(:,:,:) :: dni,dTi ! dimension: (dim_gq, dim_s, 0:1)
      real(wp),allocatable, dimension(:,:,:) :: dne,dTe

      ! Test 2D phi calculated by solving 3D Poission equation
      real(wp), allocatable, dimension(:,:) :: test_phi_2D, test_src_2D, test_res_2D
      real(wp), allocatable, dimension(:,:,:) :: test_phi_3D, test_src_3D, test_res_3D

   end type

   type(Grid_Info) :: SG ! Structure Grid

contains

   subroutine Initialize_Structured_Grid
      use field_array
      use esi_interface
      implicit none

      integer i,j,k,sid,gqid,ij,Ns,Nq
      real(wp), external :: dsda
      real(wp) tmp, this_ds

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
      allocate(SG%a_1D(Ns),SG%ori_sid(Ns),SG%s_1D(Ns),SG%n0_1D(Ns),SG%Te0_1D(Ns),SG%Ti0_1D(Ns),SG%dsda_1D(Ns),SG%safety_q_1D(Ns))
      allocate(SG%a_2D(Nq,Ns),SG%dsda_2D(Nq,Ns),SG%gq_2D(Nq,Ns),SG%Cell_dV(Nq,Ns),SG%inv_Cell_dV(Nq,Ns))
      allocate(SG%Cell_dS_s(Nq,Ns),SG%inv_Cell_dS_s(Nq,Ns),SG%Cell_dS_gq(Nq,Ns),SG%inv_Cell_dS_gq(Nq,Ns),SG%Cell_dS_gf(Nq,Ns),SG%inv_Cell_dS_gf(Nq,Ns))

      k=1
      do i=0,mpsi,SG%s_jump
         SG%ori_sid(k)=i
         SG%a_1D(k)=psi_grid(i)
         SG%a_2D(:,k)=psi_grid(i)
         SG%s_1D(k)=s_grid(i)
         SG%Te0_1D(k)=Te_grid(i)
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
      allocate(SG%dTe(Nq,Ns,0:1),SG%dTi(Nq,Ns,0:1))

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


      allocate(SG%test_phi_2D(Nq,Ns),SG%test_src_2D(Nq,Ns),SG%test_res_2D(Nq,Ns))
      allocate(SG%test_phi_3D(Nq,Ns,mzetamax),SG%test_src_3D(Nq,Ns,mzetamax),SG%test_res_3D(Nq,Ns,mzetamax))


   end subroutine Initialize_Structured_Grid

   subroutine Measure_Electron_Density_on_SG
      use precision
      use global_parameters,only:me
      use particle_array
      use particle_decomp
      use mpi_defs

      implicit none

      real(wp) ratio, ws0,ws1,wgq0,wgq1,wgf0,wgf1
      real(wp) P_s, P_gq, P_gf, P_safety_q, P_weight
      real(wp) target_gq, target_gf
      real(wp) tmp_dne(SG%dim_gq,SG%dim_s,0:1)
      integer m, sid, gqid, length, istatus(MPI_STATUS_SIZE),ierror

      real(wp) tt0,tt1

#ifdef _MGY_TIMER_
      tt0=MPI_WTIME()
#endif
      SG%dne = 0.0

      do m=1,me
         P_s = zeon(7,m)
         P_gq = modulo(zeon(2,m),2*pi)
         P_gf = modulo(zeon(3,m),2*pi)
         P_weight = zeon(5,m)

         wgf1 = abs(P_gf-SG%zetamin)*SG%inv_dgf
         wgf0 = (1.0-wgf1)

         sid = max(0,min(int((P_s-s0)*SG%inv_ds) +1,SG%dim_s-1)) ! +1 is due to Fortran's default index starting from 1
         ws1 = (P_s-SG%s_1D(sid))*SG%inv_ds
         ws0 = (1.0-ws0)

         P_safety_q = ws0*SG%safety_q_1D(sid) + ws1*SG%safety_q_1D(sid+1)

         ! Scatter particle to zetamin plane first
         target_gf = zetamin
         target_gq = modulo(P_gq + (target_gf-P_gf)/P_safety_q, 2.0*pi)
         gqid  = target_gq*SG%inv_dgq+1
         wgq1 = (target_gq-SG%gq_1D(gqid))*SG%inv_dgq
         wgq0 = 1.0-wgq1

         SG%dne(gqid  ,sid  ,0)   =  SG%dne(gqid  ,sid  ,0) + wgf0*ws0*wgq0*P_weight
         SG%dne(gqid+1,sid  ,0)   =  SG%dne(gqid+1,sid  ,0) + wgf0*ws0*wgq1*P_weight
         SG%dne(gqid  ,sid+1,0)   =  SG%dne(gqid  ,sid+1,0) + wgf0*ws1*wgq0*P_weight
         SG%dne(gqid+1,sid+1,0)   =  SG%dne(gqid+1,sid+1,0) + wgf0*ws1*wgq1*P_weight

         ! Scatter particle to zetamax plane
         target_gf = zetamax
         target_gq = modulo(P_gq + (target_gf-P_gf)/P_safety_q, 2.0*pi)
         gqid  = target_gq*SG%inv_dgq+1
         wgq1 = (target_gq-SG%gq_1D(gqid))*SG%inv_dgq
         wgq0 = 1.0-wgq1

         SG%dne(gqid  ,sid  ,1)   =  SG%dne(gqid  ,sid  ,1) + wgf1*ws0*wgq0*P_weight
         SG%dne(gqid+1,sid  ,1)   =  SG%dne(gqid+1,sid  ,1) + wgf1*ws0*wgq1*P_weight
         SG%dne(gqid  ,sid+1,1)   =  SG%dne(gqid  ,sid+1,1) + wgf1*ws1*wgq0*P_weight
         SG%dne(gqid+1,sid+1,1)   =  SG%dne(gqid+1,sid+1,1) + wgf1*ws1*wgq1*P_weight

      enddo

      if(npartdom>1) then
         tmp_dne = SG%dne
         length  = SG%dim_gq*SG%dim_s*2
         call MPI_ALLREDUCE(tmp_dne, SG%dne, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      endif

      ! Match poloidal end cell
      ! (gq=0 and gq=2pi is the same position physically
      ! although they have different gqid)
      do sid=1,SG%dim_s
         SG%dne(1,sid,:) = SG%dne(1,sid,:) + SG%dne(SG%dim_gq,sid,:)
         SG%dne(SG%dim_gq,sid,:) = SG%dne(1,sid,:)
      enddo

      ! toroidal communication
      length = SG%dim_gq*SG%dim_s
      call MPI_SENDRECV(SG%dne(:,:,0),length,mpi_Rsize,left_pe,1, tmp_dne(:,:,1),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
      SG%dne(:,:,1) = SG%dne(:,:,1) + tmp_dne(:,:,1)
      call MPI_SENDRECV(SG%dne(:,:,1),length,mpi_Rsize,right_pe,1, tmp_dne(:,:,0),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)
      SG%dne(:,:,0) = SG%dne(:,:,0) + tmp_dne(:,:,0)

      if(Flag_Load_and_WeightDef==0) then
         SG%dne(:,:,0) = eNrNm*SG%dne(:,:,0)*SG%inv_Cell_dV
         SG%dne(:,:,1) = eNrNm*SG%dne(:,:,1)*SG%inv_Cell_dV
      else
         SG%dne(:,:,0) = SG%dne(:,:,0)*SG%inv_Cell_dV
         SG%dne(:,:,1) = SG%dne(:,:,1)*SG%inv_Cell_dV
      endif

#ifdef _MGY_TIMER_
      tt1=MPI_WTIME()
      if(mype==0) write(stdout,*) "Measure_Electron_Density_on_SG takes", tt1-tt0, " secs"
      call flush(stdout)
#endif

   end subroutine Measure_Electron_Density_on_SG

   subroutine Measure_Ion_Density_on_SG
      use precision
      use global_parameters,only:mi
      use particle_array
      use particle_decomp
      use mpi_defs

      implicit none

      real(wp) ratio, ws0,ws1,wgq0,wgq1,wgf0,wgf1
      real(wp) P_s, P_gq, P_gf, P_safety_q, P_weight
      real(wp) target_gq, target_gf
      real(wp) tmp_dni(SG%dim_gq,SG%dim_s,0:1)
      integer m, sid, gqid, length, istatus(MPI_STATUS_SIZE),ierror
      real(wp) tt0, tt1

#ifdef _MGY_TIMER_
      tt0=MPI_WTIME()
#endif
      SG%dni = 0.0

      do m=1,mi
         P_s = zion(7,m)
         P_gq = modulo(zion(2,m),2*pi)
         P_gf = modulo(zion(3,m),2*pi)
         P_weight = zion(5,m)

         wgf1 = abs(P_gf-SG%zetamin)*SG%inv_dgf
         wgf0 = (1.0-wgf1)

         sid = max(0,min(int((P_s-s0)*SG%inv_ds) +1,SG%dim_s-1)) ! +1 is due to Fortran's default index starting from 1
         ws1 = (P_s-SG%s_1D(sid))*SG%inv_ds
         ws0 = (1.0-ws0)

         P_safety_q = ws0*SG%safety_q_1D(sid) + ws1*SG%safety_q_1D(sid+1)

         ! Scatter particle to zetamin plane first
         target_gf = zetamin
         target_gq = modulo(P_gq + (target_gf-P_gf)/P_safety_q, 2.0*pi)
         gqid  = target_gq*SG%inv_dgq+1
         wgq1 = (target_gq-SG%gq_1D(gqid))*SG%inv_dgq
         wgq0 = 1.0-wgq1

         SG%dni(gqid  ,sid  ,0)   =  SG%dni(gqid  ,sid  ,0) + wgf0*ws0*wgq0*P_weight
         SG%dni(gqid+1,sid  ,0)   =  SG%dni(gqid+1,sid  ,0) + wgf0*ws0*wgq1*P_weight
         SG%dni(gqid  ,sid+1,0)   =  SG%dni(gqid  ,sid+1,0) + wgf0*ws1*wgq0*P_weight
         SG%dni(gqid+1,sid+1,0)   =  SG%dni(gqid+1,sid+1,0) + wgf0*ws1*wgq1*P_weight

         ! Scatter particle to zetamax plane
         target_gf = zetamax
         target_gq = modulo(P_gq + (target_gf-P_gf)/P_safety_q, 2.0*pi)
         gqid  = target_gq*SG%inv_dgq+1
         wgq1 = (target_gq-SG%gq_1D(gqid))*SG%inv_dgq
         wgq0 = 1.0-wgq1

         SG%dni(gqid  ,sid  ,1)   =  SG%dni(gqid  ,sid  ,1) + wgf1*ws0*wgq0*P_weight
         SG%dni(gqid+1,sid  ,1)   =  SG%dni(gqid+1,sid  ,1) + wgf1*ws0*wgq1*P_weight
         SG%dni(gqid  ,sid+1,1)   =  SG%dni(gqid  ,sid+1,1) + wgf1*ws1*wgq0*P_weight
         SG%dni(gqid+1,sid+1,1)   =  SG%dni(gqid+1,sid+1,1) + wgf1*ws1*wgq1*P_weight

      enddo

      if(npartdom>1) then
         tmp_dni = SG%dni
         length  = SG%dim_gq*SG%dim_s*2
         call MPI_ALLREDUCE(tmp_dni, SG%dni, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      endif

      ! Match poloidal end cell
      ! (gq=0 and gq=2pi is the same position physically
      ! although they have different gqid)
      do sid=1,SG%dim_s
         SG%dni(1,sid,:) = SG%dni(1,sid,:) + SG%dni(SG%dim_gq,sid,:)
         SG%dni(SG%dim_gq,sid,:) = SG%dni(1,sid,:)
      enddo

      ! toroidal communication
      length = SG%dim_gq*SG%dim_s
      call MPI_SENDRECV(SG%dni(:,:,0),length,mpi_Rsize,left_pe,1, tmp_dni(:,:,1),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
      SG%dni(:,:,1) = SG%dni(:,:,1) + tmp_dni(:,:,1)
      call MPI_SENDRECV(SG%dni(:,:,1),length,mpi_Rsize,right_pe,1, tmp_dni(:,:,0),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)
      SG%dni(:,:,0) = SG%dni(:,:,0) + tmp_dni(:,:,0)

      if(Flag_Load_and_WeightDef==0) then
         SG%dni(:,:,0) = iNrNm*SG%dni(:,:,0)*SG%inv_Cell_dV
         SG%dni(:,:,1) = iNrNm*SG%dni(:,:,1)*SG%inv_Cell_dV
      else
         SG%dni(:,:,0) = SG%dni(:,:,0)*SG%inv_Cell_dV
         SG%dni(:,:,1) = SG%dni(:,:,1)*SG%inv_Cell_dV
      endif

#ifdef _MGY_TIMER_
      tt1=MPI_WTIME()
      if(mype==0) write(stdout,*) "Measure_Ion_Density_on_SG takes", tt1-tt0, " secs"
      call flush(stdout)
#endif
   end subroutine Measure_Ion_Density_on_SG

   subroutine Smooth_SG_Data_2D(num_smooth, target_Data2D, Reference_Data2D)

      use precision
      use global_parameters
      implicit none

      integer, intent(in) :: num_smooth
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

      SM_Data = 0.0
      SM_Data(1:Nq,1:Ns) = target_Data2D
      do smid=1,num_smooth
         ! Set ghost cell periodic in theta direction
         SM_Data(0,:)    = SM_Data(Nq,:)
         SM_Data(Nq+1,:) = SM_Data(1,:)
         ! Smooth in theta direction
         SM_Data(1:Nq,:) = SM_gq_Coeff(-1,1:Nq,:)*SM_Data(0:Nq-1,:) + SM_gq_Coeff(0,1:Nq,:)*SM_Data(1:Nq,:) + SM_gq_Coeff(1,1:Nq,:)*SM_Data(2:Nq+1,:)
         ! Smooth in radial direction
         SM_Data(:,1:Ns) = SM_s_Coeff(-1,:,1:Ns)*SM_Data(:,0:Ns-1) + SM_s_Coeff(0,:,1:Ns)*SM_Data(:,1:Ns) + SM_s_Coeff(1,:,1:Ns)*SM_Data(:,2:Ns+1)
      enddo
      target_Data2D = SM_Data(1:Nq,1:Ns)

   end subroutine Smooth_SG_Data_2D


   subroutine Get_Conserving_Smooth_Coeffs_on_SG(in_VOLorAREA, SM_s_Coeff, SM_gq_Coeff)

      use precision
      use global_parameters
      use field_array
      use Diagnostic_Grid
      implicit none

      integer :: Ns,Nq
      real(wp), intent(in) :: in_VOLorAREA(SG%dim_gq, SG%dim_a) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
      real(wp) :: VV(0:SG%dim_gq+1, 0:SG%dim_a+1) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
      real(wp), intent(out):: SM_gq_Coeff(-1:1,0:SG%dim_gq+1, 0:SG%dim_a+1), SM_s_Coeff(-1:1,0:SG%dim_gq+1, 0:SG%dim_a+1)

      Ns = SG%dim_a
      Nq = SG%dim_gq

      VV = 0
      VV(1:Nq,1:Ns) = in_VOLorAREA

      ! Note that theta direction is periodic
      VV(0,:) = VV(Nq,:)
      VV(Nq+1,:) = VV(1,:)

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
      use Diagnostic_Grid
      implicit none

      real(wp), dimension(mgrid), intent(in) :: oriData_1D
      real(wp), dimension(SG%dim_gq,SG%dim_a), intent(out) :: out_IntpData_2D

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

   end subroutine Interpolate_OriGridData_into_SG

end module Structured_Grid
