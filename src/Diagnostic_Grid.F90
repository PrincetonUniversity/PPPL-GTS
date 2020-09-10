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
      integer, allocatable, dimension(:) :: ori_aid
      real(wp), allocatable, dimension(:) :: a_1D, s_1D, gq_1D, n0_1D, Te0_1D, Ti0_1D
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

      integer :: num_eZD_1D, num_iZD_1D ! number of eZonData_1D and iZonData_1D
      real(wp), allocatable, dimension(:,:) :: eZonData_1D, iZonData_1D !  Zonal Data (particle, Energy, etc)
      real(wp), allocatable, dimension(:,:) :: cum_eZonData_1D, cum_iZonData_1D ! cumulative Zonal Data (particle, Energy, etc)
      real(wp), allocatable, dimension(:) :: Anulus_dV, inv_Anulus_dV

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
      real(wp), allocatable, dimension(:,:) :: fe0, fi0 ! background distribution function
   end type


   type(Grid_Info) :: DG
   type(Probe_Velocity_Space), allocatable, dimension(:) :: PVS
   integer :: nPVS


   contains

      subroutine Init_Diagnostic_Grid

         use global_parameters
         use field_array
         use precision
         use particle_decomp
         use mpi_defs
         use particle_array
         use function_gg_module
         use esi_interface
         implicit none

         integer i,j,k,aid,gqid,ij,Na,Nq
         integer pbid ! probe id
         real(wp), external :: temperature, temperature_e, density, Cal_FLabel_Function, dsda

         ! approximated delta_a and delta_gf for approximated calculation of the weighted cell volume
         real(wp) approx_da, approx_dgf, tmp
         real(wp), allocatable, dimension(:,:,:) :: tmp_DG_Data1, tmp_DG_Data2
         integer ierror, istatus(MPI_STATUS_SIZE)
         real(wp) uvec_B_cov(3), uvec_N_cov(3), uvec_BxN_cov(3)
         integer :: mi0,mi1,me0,me1

         DG%dim_a  = (mpsi+1)/DG%a_jump
         DG%dim_gq = mtheta(ipsi)/DG%gq_jump +1
         DG%num_eSD = 6 +6 +6 +1 +4 !23
         DG%num_eVD = 5

         DG%num_iSD = 8 +6 +6 +1 +4 ! 2 additional dni_Larmor, dTi_Larmor by considering larmor orbit of ions
         DG%num_iVD = 5+4 ! 4 additional energy, particle flows (total and only ExB)

         DG%num_eZD_1D = 8 ! Particle and Energy fluxes (Total; ExB; streaming) + Particle Density + Energy Density
         DG%num_iZD_1D = 8 ! Particle and Energy fluxes (Total; ExB; streaming) + Particle Density + Energy Density

         DG%dgq   = 2.0*pi/(DG%dim_gq-1)
         DG%inv_dgq = 1.0/DG%dgq
         DG%inv_ds  = (real(mpsi)/(s1-s0))/DG%a_jump

         DG%zetamax = zetamax
         DG%zetamin = zetamin

         Na=DG%dim_a
         Nq=DG%dim_gq

         allocate(DG%gq_1D(Nq))
         allocate(DG%a_1D(Na), DG%ori_aid(Na), DG%s_1D(Na), DG%n0_1D(Na), DG%Te0_1D(Na), DG%Ti0_1D(Na))
         allocate(DG%a_2D(Nq,Na))
         allocate(DG%n0_2D(Nq,Na),DG%Te0_2D(Nq,Na),DG%Ti0_2D(Nq,Na))
         allocate(DG%gq_2D(Nq,Na), DG%Cell_dV(Nq,Na), DG%inv_Cell_dV(Nq,Na), DG%dne_adiab(Nq,Na))
         allocate(DG%Cell_dS_a(Nq,Na), DG%inv_Cell_dS_a(Nq,Na), DG%Cell_dS_gq(Nq,Na), DG%inv_Cell_dS_gq(Nq,Na), DG%Cell_dS_gf(Nq,Na), DG%inv_Cell_dS_gf(Nq,Na))
         allocate(DG%FLabel_Function(Nq,Na,0:1))
         allocate(DG%FLA_Denominator(Nq,Na,0:1), DG%FLA_FLabel(Nq,Na,0:1), DG%FLA_psi0(Nq,Na,0:1))
         allocate(DG%FLA_SM_Denominator(Nq,Na), DG%FLA_SM_FLabel(Nq,Na), DG%FLA_SM_psi0(Nq,Na))
         allocate(DG%FLA_SM_FLabel_Left(Nq,Na), DG%FLA_SM_FLabel_Right(Nq,Na))
         allocate(DG%grad_FLA_SM_FLabel(3,Nq,Na,0:1),DG%uvec_grad_FL_cov(3,Nq,Na,0:1), DG%mag_grad_FLA_SM_FLabel(Nq,Na,0:1))

         allocate(DG%eScaData(Nq,Na,DG%num_eSD,0:1), DG%eVecData(3,Nq,Na,DG%num_eVD,0:1) )
         allocate(DG%iScaData(Nq,Na,DG%num_iSD,0:1), DG%iVecData(3,Nq,Na,DG%num_iVD,0:1) )

         allocate(DG%w_dne(Nq,Na), DG%w_dni(Nq,Na))
         allocate(DG%ne_tot_hybrid(Nq,Na), DG%ni_hybrid(Nq,Na), DG%Te_hybrid(Nq,Na), DG%Ti_hybrid(Nq,Na))

         allocate(DG%eZonData_1D(DG%num_eZD_1D,Na), DG%iZonData_1D(DG%num_iZD_1D,Na)) ! Zonal Data (particle, Energy, etc)
         allocate(DG%cum_eZonData_1D(DG%num_eZD_1D,Na), DG%cum_iZonData_1D(DG%num_iZD_1D,Na)) ! cumulative Zonal Data (particle, Energy, etc)
         allocate(DG%Anulus_dV(Na), DG%inv_Anulus_dV(Na))

         k=1
         do i=0,mpsi,DG%a_jump
            DG%ori_aid(k) = i
            DG%a_1D(k) = psi_grid(i)
            DG%a_2D(:,k) = psi_grid(i)
            DG%s_1D(k) = s_grid(i)

            DG%n0_1D(k)  = den_grid(i)
            DG%Te0_1D(k) = Te_grid(i)
            DG%Ti0_1D(k) = vth_grid(i)**2.0

            DG%n0_2D(:,k)= DG%n0_1D(k)
            DG%Te0_2D(:,k)= DG%Te0_1D(k)
            DG%Ti0_2D(:,k)= DG%Ti0_1D(k)
            k=k+1
         enddo

         do i=1,Nq
            DG%gq_1D(i) = (2.0*pi)/(Nq-1)*(i-1)
            DG%gq_2D(i,:) = DG%gq_1D(i)
         enddo

         ! Calculate FLabel Function and add FLA_psi0
         do j=1,Nq
            do i=1,Na
               DG%FLabel_Function(j,i,0) = Cal_FLabel_Function(DG%a_2D(j,i), DG%gq_2D(j,i), zetamin)
               DG%FLabel_Function(j,i,1) = Cal_FLabel_Function(DG%a_2D(j,i), DG%gq_2D(j,i), zetamax)
            enddo
         enddo

         !Set and Calculate Magnetic Information (ESI variables)
         allocate(DG%Jacobian(Nq,Na), DG%F(Nq,Na),DG%dFda(Nq,Na),DG%dgFda(Nq,Na),DG%dgYda(Nq,Na),DG%R(Nq,Na),DG%dRda(Nq,Na),DG%dRdgq(Nq,Na),DG%Z(Nq,Na),DG%dZda(Nq,Na),DG%dZdgq(Nq,Na),DG%BB(Nq,Na),DG%dBBda(Nq,Na),DG%dBBdgq(Nq,Na),DG%P(Nq,Na),DG%T(Nq,Na))
         allocate(DG%esi_N(Nq,Na),DG%esi_K(Nq,Na),DG%g_aa_cov(Nq,Na),DG%g_aq_cov(Nq,Na),DG%g_qq_cov(Nq,Na),DG%g_ff_cov(Nq,Na))
         allocate(DG%D(Nq,Na),DG%g_aa_cntv(Nq,Na),DG%g_aq_cntv(Nq,Na),DG%g_qq_cntv(Nq,Na),DG%g_ff_cntv(Nq,Na))
         allocate(DG%B0_cov3(3,Nq,Na))
         allocate(DG%B0_cntv3(3,Nq,Na),DG%j0_cov3(3,Nq,Na),DG%j0_cntv3(3,Nq,Na))

         allocate(DG%g_ss_cntv(Nq,Na),DG%g_sq_cntv(Nq,Na),DG%g_ss_cov(Nq,Na),DG%g_sq_cov(Nq,Na),DG%Jacob_sqf(Nq,Na))
         allocate(DG%B0_cntv3_sqf(3,Nq,Na))

         call esiget2dthese_more(DG%F,DG%dFda,DG%dgFda,DG%dgYda,DG%R,DG%dRda,DG%dRdgq,DG%Z,DG%dZda,DG%dZdgq,DG%BB,DG%dBBda,DG%dBBdgq,DG%P,DG%T,DG%a_2D,DG%gq_2D, Nq*Na)

         DG%Jacobian = DG%R*(DG%dZda*DG%dRdgq-DG%dZdgq*DG%dRda)
         DG%esi_N = (DG%dRda*DG%dRdgq+DG%dZda*DG%dZdgq)/DG%Jacobian
         DG%esi_K = (DG%dRdgq*DG%dRdgq+DG%dZdgq*DG%dZdgq)/DG%Jacobian

         DG%B0_cov3(1,:,:) = -DG%esi_N*DG%dgYda
         DG%B0_cov3(2,:,:) = -DG%esi_K*DG%dgYda
         DG%B0_cov3(3,:,:) = DG%F

         DG%B0_cntv3(1,:,:) = 0.0
         DG%B0_cntv3(2,:,:) = -DG%dgYda/DG%Jacobian
         DG%B0_cntv3(3,:,:) = DG%F/(DG%R**2)

         DG%j0_cntv3(1,:,:) = 0.0
         DG%j0_cntv3(2,:,:) = -DG%dFda/DG%Jacobian
         DG%j0_cntv3(3,:,:) = DG%P + DG%T/(DG%R**2)

         DG%g_aa_cov = DG%dRda*DG%dRda+DG%dZda*DG%dZda
         DG%g_aq_cov = DG%dRda*DG%dRdgq+DG%dZda*DG%dZdgq
         DG%g_qq_cov = DG%dRdgq*DG%dRdgq+DG%dZdgq*DG%dZdgq
         DG%g_ff_cov = DG%R**2

         DG%D=DG%dZda*DG%dRdgq-DG%dZdgq*DG%dRda

         DG%g_aa_cntv=(DG%dRdgq**2+DG%dZdgq**2)/DG%D**2   ! DG%grada^2
         DG%g_qq_cntv=(DG%dZda**2+DG%dRda**2)/DG%D**2   ! DG%gradtheta^2
         DG%g_ff_cntv=1.0/DG%R**2                 ! DG%gradphi^2
         DG%g_aq_cntv=(-DG%dRda*DG%dRdgq-DG%dZda*DG%dZdgq)/DG%D**2  ! DG%gradagradtheta


         ! Metrics regarding s coordinate
         do aid=1,DG%dim_a
            tmp=dsda(DG%a_1D(aid))
            DG%g_ss_cov(:,aid)   = DG%g_aa_cov(:,aid)/(tmp**2.0)
            DG%g_sq_cov(:,aid)   = DG%g_aq_cov(:,aid)/tmp

            DG%g_ss_cntv(:,aid)  = DG%g_aa_cntv(:,aid)*(tmp**2.0)
            DG%g_sq_cntv(:,aid)  = DG%g_aq_cntv(:,aid)*tmp

            DG%Jacob_sqf(:,aid) = DG%Jacobian(:,aid)/tmp

            DG%B0_cntv3_sqf(1,:,aid) = DG%B0_cntv3(1,:,aid)*tmp
         enddo
         DG%B0_cntv3_sqf(2:3,:,:)=DG%B0_cntv3(2:3,:,:)

         DG%j0_cov3(1,:,:) = DG%g_aa_cov*DG%j0_cntv3(1,:,:) + DG%g_aq_cov*DG%j0_cntv3(2,:,:)
         DG%j0_cov3(2,:,:) = DG%g_aq_cov*DG%j0_cntv3(1,:,:) + DG%g_qq_cov*DG%j0_cntv3(2,:,:)
         DG%j0_cov3(3,:,:) = DG%g_ff_cov*DG%j0_cntv3(3,:,:)


         ! NOTE : Get Cell volume like dV=(Jacobian)*dpsi*dtheta*dgf
         !        Following Cell_dV is an approximation
         !        For higher accuracy, a Weighted Volume by considering hellical particle shape and varying da is required
         approx_dgf=zetamax-zetamin
         do aid=1,Na
            if(aid==1) then
               approx_da=0.5*(DG%a_1D(2)-DG%a_1D(1))    ! half cell
            else if(aid==Na) then
               approx_da=0.5*(DG%a_1D(Na)-DG%a_1D(Na-1)) ! half cell
            else
               approx_da=0.5*(DG%a_1D(aid+1)-DG%a_1D(aid-1))  ! full cell
            endif

            do gqid=1,Nq
               DG%Cell_dV(gqid,aid)=DG%Jacobian(gqid,aid)*approx_da*DG%dgq*approx_dgf
               DG%Cell_dS_a(gqid,aid) = sqrt(DG%g_aa_cntv(gqid,aid))*DG%Jacobian(gqid,aid)*DG%dgq*approx_dgf
               DG%Cell_dS_gq(gqid,aid) = sqrt(DG%g_qq_cntv(gqid,aid))*DG%Jacobian(gqid,aid)*approx_da*approx_dgf
               DG%Cell_dS_gf(gqid,aid) = sqrt(DG%g_ff_cntv(gqid,aid))*DG%Jacobian(gqid,aid)*approx_da*DG%dgq
            enddo
         enddo
         DG%inv_Cell_dV=1.0/DG%Cell_dV
         DG%inv_Cell_dS_a  =1.0/DG%Cell_dS_a
         DG%inv_Cell_dS_gq =1.0/DG%Cell_dS_gq
         DG%inv_Cell_dS_gf =1.0/DG%Cell_dS_gf

         DG%Anulus_dV = mzetamax*sum(DG%Cell_dV(1:DG%dim_gq-1,:), 1)
         DG%inv_Anulus_dV = 1.0/DG%Anulus_dV

         ! Set PVS (ProbeVelocitySpace) parameters

         nPVS = nprobe
         allocate(PVS(0:nPVS)) ! 0-index PVS will be used to measure particles in a whole poloidal section.

         ! Set location and size of probes
         PVS(0)%a  = 0.5*(a0+a1)
         PVS(0)%gq = 0
         PVS(0)%gf = 0.5*(zetamax+zetamin)
         PVS(0)%L_a  = a1-a0  !To deal with a whole poloidal section
         PVS(0)%L_gq = 2.0*pi
         PVS(0)%L_gf = zetamax-zetamin

         do pbid=1,nPVS
            PVS(pbid)%a  = rprobe(pbid)
            PVS(pbid)%gq = tprobe(pbid)
            PVS(pbid)%gf = 0.5*(zetamax+zetamin)
            PVS(pbid)%L_a  = dradius
            PVS(pbid)%L_gq = dtheta
            PVS(pbid)%L_gf = zetamax-zetamin
         enddo

         ! Set other information
         do pbid=0,nPVS
            PVS(pbid)%n0  = density(PVS(pbid)%a,0.0)
            PVS(pbid)%Ti0 = temperature(PVS(pbid)%a,PVS(pbid)%gq)
            PVS(pbid)%Te0 = msm0*temperature_e(PVS(pbid)%a,PVS(pbid)%gq)

            PVS(pbid)%inv_n0  = 1.0/PVS(pbid)%n0
            PVS(pbid)%inv_Ti0 = 1.0/PVS(pbid)%Ti0
            PVS(pbid)%inv_Te0 = 1.0/PVS(pbid)%Te0

            PVS(pbid)%nv  = nv
            PVS(pbid)%dv  = delta_v
            PVS(pbid)%inv_dv = 1.0/PVS(pbid)%dv
            PVS(pbid)%max_norm_v = PVS(pbid)%nv*PVS(pbid)%dv ! or just use umax parameters (ex. 1.2*umax)

            PVS(pbid)%dim_vpara = 2*PVS(pbid)%nv+1
            PVS(pbid)%dim_vperp = PVS(pbid)%nv+1

            PVS(pbid)%num_of_fi  = 5
            PVS(pbid)%num_of_fe  = 5

            allocate(PVS(pbid)%vpara_1D(-PVS(pbid)%nv:PVS(pbid)%nv))
            allocate(PVS(pbid)%vperp_1D(0:PVS(pbid)%nv))
            allocate(PVS(pbid)%fi(-PVS(pbid)%nv:PVS(pbid)%nv, 0:PVS(pbid)%nv, PVS(pbid)%num_of_fi))
            allocate(PVS(pbid)%fe(-PVS(pbid)%nv:PVS(pbid)%nv, 0:PVS(pbid)%nv, PVS(pbid)%num_of_fe))
            allocate(PVS(pbid)%fi0(-PVS(pbid)%nv:PVS(pbid)%nv, 0:PVS(pbid)%nv))
            allocate(PVS(pbid)%fe0(-PVS(pbid)%nv:PVS(pbid)%nv, 0:PVS(pbid)%nv))

            do i=0,PVS(pbid)%nv
               PVS(pbid)%vperp_1D(i) = i*PVS(pbid)%dv
            enddo
            do i=-PVS(pbid)%nv,PVS(pbid)%nv
               PVS(pbid)%vpara_1D(i) = i*PVS(pbid)%dv
            enddo
         enddo


         ! Calculate fe0 and fi0 using temporary particles and n_T_u_group subroutines
         do pbid=0,nPVS

            ! Get back_f0_i first
            mi0=mi+1
            mi1=mi+(2*PVS(pbid)%nv+1)*(PVS(pbid)%nv+1)
            zion(1,mi0:mi1) = PVS(pbid)%a
            zion(2,mi0:mi1) = PVS(pbid)%gq
            zion(3,mi0:mi1) = PVS(pbid)%gf
            !$acc update device(zion(:,mi0:mi1))
            call n_T_u_group(mi0,mi1)
            !$acc update self(backT(mi0:mi1))
            call esigetbf_partial(zion,mi0,mi1,size(zion,1))
            k=0
            do j=0,PVS(pbid)%nv
               do i=-PVS(pbid)%nv,PVS(pbid)%nv
                  ! Note that PVS%v_para, %v_perp are the normalized velocities
                  zion(4,mi0+k) = PVS(pbid)%vpara_1D(i)*sqrt(backT(mi0+k))/BB(mi0+k)
                  zion(6,mi0+k) = PVS(pbid)%vperp_1D(j)*sqrt(backT(mi0+k))*sqrt(1.0/(2.0*BB(mi0+k)))
                  k=k+1
               enddo
            enddo
            !$acc update device(zion(:,mi0:mi1))
            call n_T_u_group(mi0,mi1,flag_backf0=.true.) ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))
            !$acc update self(back_f0_i(mi0:mi1))


            ! Calculate electron background distribution function
            me0=me+1
            me1=me+(2*PVS(pbid)%nv+1)*(PVS(pbid)%nv+1)
            zeon(1,me0:me1) = PVS(pbid)%a
            zeon(2,me0:me1) = PVS(pbid)%gq
            zeon(3,me0:me1) = PVS(pbid)%gf
            !$acc update device(zeon(:,me0:me1))
            call n_T_u_group_e(me0,me1)
            !$acc update self(backTe(me0:me1))
            zeon(4,me0:me1) = zion(4,mi0:mi1)*sqrt(msm0)*c0cs*sqrt(msm0*backTe(me0:me1)/backT(mi0:mi1))
            zeon(6,me0:me1) = zion(6,mi0:mi1)*sqrt(msm0*backTe(me0:me1)/backT(mi0:mi1))
            !$acc update device(zeon(:,me0:me1))
            call n_T_u_group_e(me0,me1,flag_backf0=.true.) ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))
            !$acc update self(back_f0_e(me0:me1))


            ! Substitute particle's back_f0 to PVS f0 grid value
            k=0
            do j=0,PVS(pbid)%nv
               do i=-PVS(pbid)%nv,PVS(pbid)%nv
                  PVS(pbid)%fi0(i,j) = back_f0_i(mi0+k)
                  PVS(pbid)%fe0(i,j) = back_f0_e(me0+k)
                  k=k+1
               enddo
            enddo
         enddo



         if(myrank_partd==nproc_partd-1) then

            call Cal_FLA_Denominator_on_DG_using_gPPI

            allocate(tmp_DG_Data1(Nq, Na,0:1),tmp_DG_Data2(Nq, Na,0:1))

            DG%FLA_SM_Denominator=DG%FLA_Denominator(:,:,1)
            call Smooth_DG_Data_2D(3, DG%FLA_SM_Denominator)

            call Do_FieldLineAverage_on_DG_by_using_gPPI([DG%a_2D,DG%a_2D], DG%FLA_psi0)
            DG%FLA_SM_psi0 = DG%FLA_psi0(:,:,1)
            call Smooth_DG_Data_2D(3, DG%FLA_SM_psi0)

            ! Cal FLabel using repetitive FieldLineAverage
            call Cal_FLA_FLabel_on_DG_by_using_gPPI_iniFLabel(DG%FLA_FLabel)
            tmp_DG_Data1 = DG%FLA_FLabel
            do i=1,4
               call Smooth_DG_Data_2D(2, tmp_DG_Data1(:,:,0))
               call Smooth_DG_Data_2D(2, tmp_DG_Data1(:,:,1))
               call Do_FieldLineAverage_on_DG_by_using_gPPI(tmp_DG_Data1, tmp_DG_Data2)
               tmp_DG_Data1=tmp_DG_Data2
            enddo
            ! DG%FLA_FLabel = tmp_DG_Data1
            DG%FLA_SM_FLabel = tmp_DG_Data1(:,:,1)
            call Smooth_DG_Data_2D(1, DG%FLA_SM_FLabel)

            deallocate(tmp_DG_Data1,tmp_DG_Data2)
         endif

         ! Broadcast calculated FLA data
         call MPI_BCAST(DG%FLA_Denominator, 2*Na*Nq,mpi_Rsize, nproc_partd-1,partd_comm,ierror)
         call MPI_BCAST(DG%FLA_FLabel,      2*Na*Nq,mpi_Rsize, nproc_partd-1,partd_comm,ierror)
         call MPI_BCAST(DG%FLA_psi0,        2*Na*Nq,mpi_Rsize, nproc_partd-1,partd_comm,ierror)
         call MPI_BCAST(DG%FLA_SM_Denominator, Na*Nq,mpi_Rsize, nproc_partd-1,partd_comm,ierror)
         call MPI_BCAST(DG%FLA_SM_FLabel,      Na*Nq,mpi_Rsize, nproc_partd-1,partd_comm,ierror)
         call MPI_BCAST(DG%FLA_SM_psi0,        Na*Nq,mpi_Rsize, nproc_partd-1,partd_comm,ierror)

         ! Communicate btw toroidal_comm to get Left and Right FLA_SM_FLabel
         call MPI_SENDRECV(DG%FLA_SM_FLabel, Nq*Na, mpi_Rsize, right_pe, myrank_toroidal, &
         DG%FLA_SM_FLabel_Left, Nq*Na, mpi_Rsize, left_pe, left_pe, toroidal_comm, istatus, ierror)
         call MPI_SENDRECV(DG%FLA_SM_FLabel, Nq*Na, mpi_Rsize, left_pe, myrank_toroidal, &
         DG%FLA_SM_FLabel_Right, Nq*Na, mpi_Rsize, right_pe, right_pe, toroidal_comm, istatus, ierror)


         !Calculate normal vector (gradient) of FLA_SM_FLabel
         Na=Na
         Nq=Nq
         ! d(FL)/da
         DG%grad_FLA_SM_FLabel(1,:,1,1) = (DG%FLA_SM_FLabel(:,2)-DG%FLA_SM_FLabel(:,1))/(DG%a_1D(2)-DG%a_1D(1))
         DG%grad_FLA_SM_FLabel(1,:,2:Na-1,1) = (DG%FLA_SM_FLabel(:,3:Na)-DG%FLA_SM_FLabel(:,1:Na-2))/(DG%a_2D(:,3:Na)-DG%a_2D(:,1:Na-2))
         DG%grad_FLA_SM_FLabel(1,:,Na,1) = (DG%FLA_SM_FLabel(:,Na)-DG%FLA_SM_FLabel(:,Na-1))/(DG%a_1D(Na)-DG%a_1D(Na-1))

         ! d(FL)/dgq with periodicity
         DG%grad_FLA_SM_FLabel(2,1,:,1) = (DG%FLA_SM_FLabel(2,:)-DG%FLA_SM_FLabel(Nq,:))/(2.0*DG%dgq)
         DG%grad_FLA_SM_FLabel(2,2:Nq-1,:,1) = (DG%FLA_SM_FLabel(3:Nq,:)-DG%FLA_SM_FLabel(1:Nq-2,:))/(2.0*DG%dgq)
         DG%grad_FLA_SM_FLabel(2,Nq,:,1) = DG%grad_FLA_SM_FLabel(2,1,:,1)

         ! d(FL)/dgf with periodicity
         DG%grad_FLA_SM_FLabel(3,:,:,1) = (DG%FLA_SM_FLabel_Right - DG%FLA_SM_FLabel_Left)/(2.0*(zetamax-zetamin))

         call MPI_SENDRECV(DG%grad_FLA_SM_FLabel(:,:,:,1), 3*Nq*Na, mpi_Rsize, right_pe, myrank_toroidal, &
                           DG%grad_FLA_SM_FLabel(:,:,:,0), 3*Nq*Na, mpi_Rsize, left_pe, left_pe, toroidal_comm, istatus, ierror)



         ! Calculate magnitude of grad(FL)
         do k=0,1
            DG%mag_grad_FLA_SM_FLabel(:,:,k) = sqrt(  DG%g_aa_cntv*(DG%grad_FLA_SM_FLabel(1,:,:,k)**2.0) &
                                                      +2.0*DG%g_aq_cntv*DG%grad_FLA_SM_FLabel(1,:,:,k)*DG%grad_FLA_SM_FLabel(2,:,:,k) &
                                                      +DG%g_qq_cntv*(DG%grad_FLA_SM_FLabel(2,:,:,k)**2.0)  &
                                                      +DG%g_ff_cntv*(DG%grad_FLA_SM_FLabel(3,:,:,k)**2.0))
         enddo

         do k=1,3
            DG%uvec_grad_FL_cov(k,:,:,:) = DG%grad_FLA_SM_FLabel(k,:,:,:)/(DG%mag_grad_FLA_SM_FLabel)
         enddo

      end subroutine Init_Diagnostic_Grid


      subroutine Cal_unit_vectors_tangential_normal_binormal_to_FluxSurface(in_a, in_gq, in_gf, in_s, uvec_B_cov, uvec_B_cntv, uvec_N_cov, uvec_N_cntv, uvec_BxN_cov, uvec_BxN_cntv)
         use precision
         use global_parameters
         use field_array
         use esi_interface

         implicit none
         real(wp),intent(in) :: in_a, in_gq, in_gf, in_s

         real(wp) target_gq ! projected gq at nearby poloidal plane following the B0 field
         integer :: ori_sid !sid of original full grid
         real(wp) ori_inv_ds ! inv_ds of original full grid
         real(wp) ff, ws0,ws1, wgq0,wgq1, wgf,ratio ! ratios along each dimension
         real(wp) target_zeta, safety_q
         integer :: j,k, sid, gqid, qq, ss

         real(wp):: F_,dFda_,dgFda_,dgYda_,R_,dRda_,dRdgq_,Z_,dZda_,dZdgq_,BB_,dBBda_,dBBdgq_,P_,T_
         real(wp) Jacobian_, D_
         real(wp):: g_aa_cov_, g_aq_cov_, g_qq_cov_, g_ff_cov_, g_aa_cntv_, g_aq_cntv_, g_qq_cntv_, g_ff_cntv_

         real(wp) tot_B_vec_cov(3), tot_B_vec_cntv(3)
         real(wp) vec_BxN_cntv(3), vec_BxN_cov(3) !Note tha theses are not the unit vectors
         real(wp), intent(out) :: uvec_B_cov(3), uvec_N_cov(3), uvec_BxN_cov(3) ! unit vectors in total B direction, Normal (gradFL), and binormal
         real(wp), intent(out) :: uvec_B_cntv(3), uvec_N_cntv(3), uvec_BxN_cntv(3)

         ori_inv_ds = real(mpsi)/(s1-s0)

         !Calculate Interpolation of uvec_N_cov which uses the precalculated DG%uvec_grad_FL_cov on DG nodes
         !Get uvec_B_cov(3), uvec_N_cov(3), uvec_BxN_cov(3)
         call esiget_Intp_B_vec_cov_and_cntv(in_a, in_gq, in_gf, istep*tstep, tot_B_vec_cov, tot_B_vec_cntv)
         uvec_B_cov   = tot_B_vec_cov/sqrt(dot_product(tot_B_vec_cov,tot_B_vec_cntv))

         call esiget2dthese_more(F_,dFda_,dgFda_,dgYda_,R_,dRda_,dRdgq_,Z_,dZda_,dZdgq_,BB_,dBBda_,dBBdgq_,P_,T_,in_a,in_gq, 1)
         Jacobian_ = R_*(dZda_*dRdgq_-dZdgq_*dRda_)
         D_=dZda_*dRdgq_-dZdgq_*dRda_

         g_aa_cov_   = dRda_*dRda_+dZda_*dZda_;
         g_aq_cov_  = dRda_*dRdgq_+dZda_*dZdgq_;
         g_qq_cov_ = dRdgq_*dRdgq_+dZdgq_*dZdgq_;
         g_ff_cov_ = R_**2

         g_aa_cntv_=(dRdgq_**2+dZdgq_**2)/D_**2   ! grada^2
         g_qq_cntv_=(dZda_**2+dRda_**2)/D_**2   ! gradtheta^2
         g_ff_cntv_=1.0_wp/R_**2                 ! gradphi^2
         g_aq_cntv_=(-dRda_*dRdgq_-dZda_*dZdgq_)/D_**2  ! gradagradtheta


         ori_sid = (in_s-s0)*ori_inv_ds ! s index of original full grid
         ff  = (in_s-s_grid(ori_sid))*ori_inv_ds
         safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)
         sid = max(1,min(int((in_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
         ws1 = (in_s-DG%s_1D(sid))*DG%inv_ds
         ws0 = 1.0-ws1

         uvec_N_cov = 0
         do k=0,1
            if(k==0) target_zeta =zetamin
            if(k==1) target_zeta =zetamax

            wgf = (1.0-abs(in_gf-target_zeta)/deltaz)

            target_gq = modulo(in_gq + (target_zeta-in_gf)/safety_q, 2.0*pi)
            gqid = (target_gq)*DG%inv_dgq +1
            wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
            wgq0 = 1.0-wgq1
            do j=1,4
               select case(j)
                  case(1)
                     qq=gqid; ss=sid
                     ratio=ws0*wgq0*wgf
                  case(2)
                     qq=gqid+1; ss=sid
                     ratio=ws0*wgq1*wgf
                  case(3)
                     qq=gqid; ss=sid+1
                     ratio=ws1*wgq0*wgf
                  case(4)
                     qq=gqid+1; ss=sid+1
                     ratio=ws1*wgq1*wgf
               end select
               uvec_N_cov = uvec_N_cov+ratio*DG%uvec_grad_FL_cov(:,qq,ss,k)
            enddo
         enddo

         ! Remove a possible small error part of uvec_N_cov which is parallel to the magnetic field
         uvec_B_cntv = tot_B_vec_cntv/sqrt(dot_product(tot_B_vec_cov,tot_B_vec_cntv))
         uvec_N_cov= uvec_N_cov - dot_product(uvec_N_cov,uvec_B_cntv)*uvec_B_cov

         uvec_N_cov = uvec_N_cov/sqrt( g_aa_cntv_*uvec_N_cov(1)**2.0 + 2.0*g_aq_cntv_*uvec_N_cov(1)*uvec_N_cov(2) + g_qq_cntv_*uvec_N_cov(2)**2.0 + g_ff_cntv_*uvec_N_cov(3)**2.0)
         uvec_N_cntv(1) = g_aa_cntv_*uvec_N_cov(1) + g_aq_cntv_*uvec_N_cov(2)
         uvec_N_cntv(2) = g_aq_cntv_*uvec_N_cov(1) + g_qq_cntv_*uvec_N_cov(2)
         uvec_N_cntv(3) = g_ff_cntv_*uvec_N_cov(3)

         vec_BxN_cntv(1) = (1.0/Jacobian_)*(uvec_B_cov(2)*uvec_N_cov(3)-uvec_B_cov(3)*uvec_N_cov(2))
         vec_BxN_cntv(2) = (1.0/Jacobian_)*(uvec_B_cov(3)*uvec_N_cov(1)-uvec_B_cov(1)*uvec_N_cov(3))
         vec_BxN_cntv(3) = (1.0/Jacobian_)*(uvec_B_cov(1)*uvec_N_cov(2)-uvec_B_cov(2)*uvec_N_cov(1))

         vec_BxN_cov(1) = g_aa_cov_*vec_BxN_cntv(1) + g_aq_cov_*vec_BxN_cntv(2)
         vec_BxN_cov(2) = g_aq_cov_*vec_BxN_cntv(1) + g_qq_cov_*vec_BxN_cntv(2)
         vec_BxN_cov(3) = g_ff_cov_*vec_BxN_cntv(3)

         uvec_BxN_cov   = vec_BxN_cov/sqrt(dot_product(vec_BxN_cov,vec_BxN_cntv))
         uvec_BxN_cntv  = vec_BxN_cntv/sqrt(dot_product(vec_BxN_cov,vec_BxN_cntv))

      end subroutine Cal_unit_vectors_tangential_normal_binormal_to_FluxSurface


      subroutine Measure_Electrons_on_Diagnostic_Grid(in_P_data, Full_dadgqdgf)

         use precision
         use global_parameters
         use field_array
         use mpi_defs
         use particle_array
         use particle_decomp
         use esi_interface

         implicit none

         real(wp), intent(in) :: in_P_data(3,me) ! additional particle data from input
         real(wp), intent(in) :: Full_dadgqdgf(3,me)
         real(wp) energy(me), KE_para(me), KE_perp(me)
         !P varaiables are at particle position
         real(wp) P_s, P_safety_q !safety factor
         real(wp) P_a, P_gq, P_gf

         real(wp) delta_f(me), f_tot(me) ! delta_f and f_tot (f_tot=f0+delta_f)

         real(wp) target_gq ! projected gq at nearby poloidal plane following the B0 field
         integer :: ori_sid !sid of original full grid
         real(wp) ori_inv_ds,ori_inv_dgf, P_inv_ori_dgq ! inv_ds of original full grid
         real(wp) ff, ws0,ws1, wgq0,wgq1, wgf,ratio ! ratios along each dimension
         real(wp) target_zeta
         integer :: i,m, j,k, aid, sid, vid, gqid, qq, ss
         integer ::pbid,vxid,vyid,fid
         real(wp):: norm_vpara, norm_vperp, rvx,rvy ! norm_vpara=vpara/vth, norm_vperp=vperp/vth, where vth=sqrt(T/m)
         real(wp),allocatable, dimension(:,:,:) :: tmp_PVS_fe
         integer :: length, istatus(MPI_STATUS_SIZE),ierror

         real(wp) tmp_Scalar_Data(DG%dim_gq,DG%dim_a,DG%num_eSD,0:1), tmp_Vector_Data(3,DG%dim_gq,DG%dim_a,DG%num_eVD,0:1)
         real(wp) :: intpl_dphi(DG%dim_gq, DG%dim_a), inv_tot_ne(DG%dim_gq, DG%dim_a)
         real(wp) :: intpl_FLA_Te0(DG%dim_gq,DG%dim_a)
         real(wp) :: uvec_N_cov(3), uvec_B_cov(3), uvec_BxN_cov(3) ! unit vectors in total B direction, Normal (gradFL), and binormal
         real(wp) :: uvec_N_cntv(3), uvec_B_cntv(3), uvec_BxN_cntv(3) ! unit vectors in total B direction, Normal (gradFL), and binormal


         if(acc_newer_D_zeon) then
            !$acc update self(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_H_zeon=.false.
            acc_newer_D_zeon=.false.
         endif

         DG%eScaData=0
         DG%eVecData=0
         do pbid=0,nPVS
            PVS(pbid)%fe = 0
         enddo

         ori_inv_ds = real(mpsi)/(s1-s0)
         ori_inv_dgf= 1.0/(zetamax-zetamin)

         energy(:) = in_P_data(2,:)
         KE_para = 0.5*for_E*(zeon(4,1:me)*BB(1:me))**2.0
         KE_perp = zeon(6,1:me)*zeon(6,1:me)*BB(1:me)

         ! Calculate dne_adiabatic by interpolating dphi on original grids
         call Interpolate_OriGridData_into_DG(phi(1,:), intpl_dphi)
         call Interpolate_OriGridData_into_DG(FLA_Te0(1,:), intpl_FLA_Te0)
         if (Flag_Def_dne_adiab==0) then
            do aid=1,DG%dim_a
               i = DG%ori_aid(aid)
               DG%dne_adiab(:,aid)=den_grid(i)*abs(csc0)*(intpl_dphi(:,aid)-real(1.0-ZF_SEPARATION)*phi00(i))/Te_grid(i)
            enddo
         elseif(Flag_Def_dne_adiab==1) then
            do aid=1,DG%dim_a
               i = DG%ori_aid(aid)
               DG%dne_adiab(:,aid)=den_grid(i)*abs(csc0)*(intpl_dphi(:,aid)+real(ZF_SEPARATION)*phi00(i))/Te_grid(i)
            enddo
         elseif(Flag_Def_dne_adiab==2) then
            do aid=1,DG%dim_a
               do gqid=1,DG%dim_gq
                  i = DG%ori_aid(aid)
                  DG%dne_adiab(gqid,aid)=den_grid(i)*abs(csc0)*(intpl_dphi(gqid,aid)+real(ZF_SEPARATION)*phi00(i))/intpl_FLA_Te0(gqid,aid)
               enddo
            enddo
         elseif(Flag_Def_dne_adiab==3) then
            do aid=1,DG%dim_a
               i = DG%ori_aid(aid)
               DG%dne_adiab(:,aid)=den_grid(i)*(-1.0+exp(abs(csc0)*(intpl_dphi(:,aid)+real(ZF_SEPARATION)*phi00(i))/Te_grid(i)))
            enddo
         endif


         if(Flag_Load_and_WeightDef==0)then
            delta_f(:) = zeon0(5,1:me)*eNrNm
            f_tot(:)   = eNrNm
         else
            delta_f(:) = zeon0(5,1:me)
            f_tot(:)   = zeon0(6,1:me)
         endif

         do m=1,me

            if (Flag_w_Threshold==1 .and. (zeon0(5,m)<w_Thold_min .or. zeon0(5,m)>w_Thold_max)) cycle

            P_a  = zeon0(1,m)
            P_gq = modulo(zeon0(2,m),2*pi)
            P_gf = modulo(zeon0(3,m),2*pi)
            P_s  = zeon0(7,m)

            ori_sid = (P_s-s0)*ori_inv_ds ! s index of original full grid
            ff  = (P_s-s_grid(ori_sid))*ori_inv_ds
            P_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)
            P_inv_ori_dgq = 1.0/(ff*deltat(ori_sid+1)+(1.0-ff)*deltat(ori_sid))

            sid = max(1,min(int((P_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
            ws1 = (P_s-DG%s_1D(sid))*DG%inv_ds
            ws0 = 1.0-ws1

            call Cal_unit_vectors_tangential_normal_binormal_to_FluxSurface(P_a,P_gq,P_gf,P_s, uvec_B_cov, uvec_B_cntv, uvec_N_cov, uvec_N_cntv, uvec_BxN_cov, uvec_BxN_cntv)

            ! Scatter data into neighboring DG nodes
            do k=0,1
               if(k==0) target_zeta =zetamin
               if(k==1) target_zeta =zetamax

               wgf = (1.0-abs(P_gf-target_zeta)/deltaz)

               target_gq = modulo(P_gq + (target_zeta-P_gf)/P_safety_q, 2.0*pi)
               gqid = (target_gq)*DG%inv_dgq +1
               wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
               wgq0 = 1.0-wgq1

               do j=1,4
                  select case(j)
                     case(1)
                        qq=gqid; ss=sid
                        ratio=ws0*wgq0*wgf
                     case(2)
                        qq=gqid+1; ss=sid
                        ratio=ws0*wgq1*wgf
                     case(3)
                        qq=gqid; ss=sid+1
                        ratio=ws1*wgq0*wgf
                     case(4)
                        qq=gqid+1; ss=sid+1
                        ratio=ws1*wgq1*wgf
                  end select

                  vid=1; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m) ! Density
                  vid=2; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*(2.0/3.0*energy(m)-backTe(m)*msm0)*delta_f(m) !Temperature definition 1
                  vid=3; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*(2.0/3.0*energy(m))*delta_f(m) !
                  vid=4; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio ! Marker Distribution
                  vid=5; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
                  vid=6; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*f_tot(m)*(2.0/3.0*energy(m)) ! Integrating f_tot*energy, so Te_tot

                  vid=20; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*(2.0*KE_para(m))*delta_f(m) !
                  vid=21; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*(KE_perp(m))*delta_f(m) !
                  vid=22; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*f_tot(m)*(2.0*KE_para(m)) ! Te_tot_para
                  vid=23; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*f_tot(m)*(KE_perp(m)) ! Te_tot_perp

                  if(elabels(m)==-1) then
                     vid=19; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio ! Negative df Markers from Walls
                  endif

                  ! Volumetric delta-f Partcle Flux (=df*v_tot) (#/(m^2*s)) along tangential, normal, binormal directions, respectively
                  vid=7; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_B_cov)
                  vid=8; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_N_cov)
                  vid=9; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_BxN_cov)

                  vid=10; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_B_cntv(1))*uvec_B_cov)
                  vid=11; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_N_cntv(1))*uvec_N_cov)
                  vid=12; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_BxN_cntv(1))*uvec_BxN_cov)

                  ! Volumetric delta-f Energy Flux (=df*KE*v_tot) (#/(m^2*s)) along tangential, normal, binormal directions, respectively
                  vid=13; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_B_cov)
                  vid=14; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_N_cov)
                  vid=15; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_BxN_cov)

                  vid=16; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_B_cntv(1))*uvec_B_cov)
                  vid=17; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_N_cntv(1))*uvec_N_cov)
                  vid=18; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_BxN_cntv(1))*uvec_BxN_cov)


                  vid=1 ! Turbulent Particle Flux
                  ! DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*(delta_f(m)*ghq(m)) ! Radial Flux
                  ! DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*(delta_f(m)*Ta(m)) ! theta Flux
                  ! DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*(delta_f(m)*dve_gf(m)) ! zeta Flux
                  DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*delta_f(m)*Full_dadgqdgf(1,m)*ee_inv_grad_a(m) ! Radial Flux
                  DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*delta_f(m)*Full_dadgqdgf(2,m)*ee_inv_grad_gq(m) ! theta Flux
                  DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*delta_f(m)*Full_dadgqdgf(3,m)*ee_inv_grad_gf(m) ! zeta Flux


                  vid=2 ! Turbulent Particle Flux due to ExB only
                  DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*(delta_f(m)*(ghq(m)-vpara_dve_a(m)))*ee_inv_grad_a(m) ! Radial Flux
                  DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*(delta_f(m)*(Ta(m)-vpara_dve_gq(m)))*ee_inv_grad_gq(m) ! theta Flux
                  DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*(delta_f(m)*(dve_gf(m)-vpara_dve_gf(m)))*ee_inv_grad_gf(m)! zeta Flux

                  vid=3 ! Turbulent Energy Flux
                  ! DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*ghq(m)) ! Radial Flux
                  ! DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Ta(m)) ! theta Flux
                  ! DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*dve_gf(m)) ! zeta Flux
                  DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Full_dadgqdgf(1,m))*ee_inv_grad_a(m) ! Radial Flux
                  DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Full_dadgqdgf(2,m))*ee_inv_grad_gq(m) ! theta Flux
                  DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Full_dadgqdgf(3,m))*ee_inv_grad_gf(m) ! zeta Flux

                  vid=4 ! Turbulent Energy Flux due to ExB only
                  DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*(ghq(m)-vpara_dve_a(m)))*ee_inv_grad_a(m) ! Radial Flux
                  DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*(Ta(m)-vpara_dve_gq(m)))*ee_inv_grad_gq(m) ! theta Flux
                  DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*(dve_gf(m)-vpara_dve_gf(m)))*ee_inv_grad_gf(m)! zeta Flux

                  vid=5 ! To measure CFL condition (v*dt/dx)
                  DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*abs(tstep*((Tc/backT(m))**dsda_power*Full_dadgqdgf(1,m))*ori_inv_ds)
                  DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*abs(tstep*(Full_dadgqdgf(2,m)-Full_dadgqdgf(3,m)/P_safety_q)*P_inv_ori_dgq) !time derivative of target_gq=gq-gf/safety_q
                  DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*abs(tstep*Full_dadgqdgf(3,m)*ori_inv_dgf)
               enddo
            enddo


            ! Velocity Space Measurement
            do pbid=0,nPVS
               target_gq = modulo(P_gq + (PVS(pbid)%gf-P_gf)/P_safety_q, 2.0*pi)

               if( abs(P_gf-PVS(pbid)%gf)<0.5*PVS(pbid)%L_gf &
               .and. abs(P_a-PVS(pbid)%a) <0.5*PVS(pbid)%L_a &
               .and. min(abs(target_gq-PVS(pbid)%gq), 2.0*pi-abs(target_gq-PVS(pbid)%gq)) < 0.5*PVS(pbid)%L_gq) then

                  norm_vpara = zeon0(4,m)*csc0*BB(m)*sqrt(m0ms*PVS(pbid)%inv_Te0) ! norm_vpara = vpara/sqrt(T/m) = rho_para*(eB/m)/sqrt(T/m) = rho_para*(e*B)/sqrt(m*T)
                  norm_vperp = zeon(6,m)*sqrt(2.0*BB(m))*sqrt(PVS(pbid)%inv_Te0) ! norm_vperp = sqrt(mu)*sqrt(2B/m)/sqrt(T/m) = sqrt(mu)*sqrt(2B)/(sqrt(T))

                  vxid = floor(norm_vpara*PVS(pbid)%inv_dv)
                  vyid = floor(norm_vperp*PVS(pbid)%inv_dv)
                  rvx  = (norm_vpara-PVS(pbid)%vpara_1D(vxid))*PVS(pbid)%inv_dv
                  rvy  = (norm_vperp-PVS(pbid)%vperp_1D(vyid))*PVS(pbid)%inv_dv

                  ! Do not calculate any contribution from "too" outside of our probe boudnary
                  ! "too" outside means that even over the ghost cells around the probe boundary
                  if(vxid< -PVS(pbid)%nv-1) then
                     continue
                  elseif(vxid>=PVS(pbid)%nv+1) then
                     continue
                  elseif(vyid>=PVS(pbid)%nv+1) then
                     continue
                  endif

                  if(vxid< -PVS(pbid)%nv) then
                     vxid= -PVS(pbid)%nv
                     rvx=0.0
                  elseif (vxid>= PVS(pbid)%nv) then
                     vxid=  PVS(pbid)%nv-1
                     rvx=1.0
                  endif

                  if(vyid>= PVS(pbid)%nv) then
                     vyid=  PVS(pbid)%nv-1
                     rvy=1.0
                  endif

                  do j=1,4
                     select case(j)
                        case(1)
                           qq=vxid; ss=vyid
                           ratio=(1.0-rvx)*(1.0-rvy)
                        case(2)
                           qq=vxid+1; ss=vyid
                           ratio=(rvx)*(1.0-rvy)
                        case(3)
                           qq=vxid; ss=vyid+1
                           ratio=(1.0-rvx)*(rvy)
                        case(4)
                           qq=vxid+1; ss=vyid+1
                           ratio=(rvx)*(rvy)
                     end select

                     fid=1; PVS(pbid)%fe(qq,ss,fid) = PVS(pbid)%fe(qq,ss,fid) + ratio*delta_f(m) ! Density
                     fid=2; PVS(pbid)%fe(qq,ss,fid) = PVS(pbid)%fe(qq,ss,fid) + ratio*delta_f(m)*(2.0/3.0*energy(m)-backT(m)) ! Temperature definition
                     fid=3; PVS(pbid)%fe(qq,ss,fid) = PVS(pbid)%fe(qq,ss,fid) + ratio ! Marker Density
                     fid=4; PVS(pbid)%fe(qq,ss,fid) = PVS(pbid)%fe(qq,ss,fid) + ratio*f_tot(m) ! f_tot
                     fid=5; PVS(pbid)%fe(qq,ss,fid) = PVS(pbid)%fe(qq,ss,fid) + ratio*f_tot(m)*(2.0/3.0*energy(m)) ! f_tot*energy
                  enddo
               endif
            enddo

         enddo


         if(npartdom>1) then
            tmp_Scalar_Data = DG%eScaData
            tmp_Vector_Data = DG%eVecData

            length = DG%dim_gq*DG%dim_a*DG%num_eSD*2
            call MPI_ALLREDUCE(tmp_Scalar_Data, DG%eScaData, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)

            length = 3*DG%dim_gq*DG%dim_a*DG%num_eVD*2
            call MPI_ALLREDUCE(tmp_Vector_Data, DG%eVecData, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)

            allocate(tmp_PVS_fe, mold=PVS(0)%fe)
            do pbid=0,nPVS
               tmp_PVS_fe = PVS(pbid)%fe
               length = PVS(pbid)%dim_vpara*PVS(pbid)%dim_vperp*PVS(pbid)%num_of_fe
               call MPI_ALLREDUCE(tmp_PVS_fe, PVS(pbid)%fe, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
            enddo
            deallocate(tmp_PVS_fe)
         endif


         ! Match poloidal end cell
         ! (gq=0 and gq=2pi is the same position physically
         ! although they have different gqid)
         do sid=1,DG%dim_a
            DG%eScaData(1,sid,:,:) = DG%eScaData(1,sid,:,:) + DG%eScaData(DG%dim_gq,sid,:,:)
            DG%eScaData(DG%dim_gq,sid,:,:) = DG%eScaData(1,sid,:,:)

            DG%eVecData(:,1,sid,:,:) = DG%eVecData(:,1,sid,:,:) + DG%eVecData(:,DG%dim_gq,sid,:,:)
            DG%eVecData(:,DG%dim_gq,sid,:,:) = DG%eVecData(:,1,sid,:,:)
         enddo


         ! toroidal end cell
         length = DG%dim_gq*DG%dim_a*DG%num_eSD
         call MPI_SENDRECV(DG%eScaData(:,:,:,0),length,mpi_Rsize,left_pe,1, tmp_Scalar_Data(:,:,:,1),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         DG%eScaData(:,:,:,1) = DG%eScaData(:,:,:,1) + tmp_Scalar_Data(:,:,:,1)



         length = 3*DG%dim_gq*DG%dim_a*DG%num_eVD
         call MPI_SENDRECV(DG%eVecData(:,:,:,:,0),length,mpi_Rsize,left_pe,1, tmp_Vector_Data(:,:,:,:,1),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         DG%eVecData(:,:,:,:,1) = DG%eVecData(:,:,:,:,1) + tmp_Vector_Data(:,:,:,:,1)


         ! Do normalization and post-processing
         ! Divide by Cell Volume to get density quantities
         do vid=1,DG%num_eSD
            DG%eScaData(:,:,vid,1) = DG%eScaData(:,:,vid,1)*DG%inv_Cell_dV
         enddo

         ! Calculate hybrid particle and energy density
         DG%w_dne(:,:) = abs(DG%eScaData(:,:,1,1))/(abs(DG%eScaData(:,:,1,1))+abs(DG%eScaData(:,:,5,1))) ! wdne=abs(snap.DG.dne_tot)./(abs(snap.DG.ne_tot)+abs(snap.DG.dne_tot));
         DG%ne_tot_hybrid = DG%w_dne*DG%eScaData(:,:,5,1) + (1.0-DG%w_dne)*(DG%eScaData(:,:,1,1)+DG%n0_2D + int(ifelectron/2)*DG%dne_adiab) ! (n_hybrid = w_dne.*snap.DG.ne_tot+(1.0-w_dne).*(snap.DG.dne_tot+snap.n0_2D))
         DG%Te_hybrid  = DG%w_dne*DG%eScaData(:,:,6,1) + (1.0-DG%w_dne)*(DG%eScaData(:,:,3,1) +DG%Te0_2D*(DG%n0_2D+int(ifelectron/2)*DG%dne_adiab)) ! (hybrid (2/3)*energy_density, then will be divded by n_hybrid to get T_hybrid)
         DG%Te_hybrid  = DG%Te_hybrid / DG%ne_tot_hybrid

         ! Calculate inverse of tot density (including dne_adiabatic in the case of delta-h electron (ifelectron==2))
         inv_tot_ne = 1.0/(DG%n0_2D + DG%eScaData(:,:,1,1) + int(ifelectron/2)*DG%dne_adiab) ! smoothed total density using df or dh
         ! Convert Energy to Temperature
         DG%eScaData(:,:,2,1) = DG%eScaData(:,:,2,1)/DG%n0_2D ! Convert Energy to Temperature
         ! Consistent Temperature by considering density fluctuation in the denominator
         ! delta_T = ((2/3)*dEnergy - T0*dn)/(n0+dn_h+dne_adiab)], where dEnergy=Integral[(df or dh)*energy(m)]+1.5*T0*dne_adiab
         DG%eScaData(:,:,3,1) = (DG%eScaData(:,:,3,1)-DG%Te0_2D*DG%eScaData(:,:,1,1))*inv_tot_ne !delta_T = (2/3*dEnergy - T0*dn)/(n0+dn_h+dne_adiab)]
         DG%eScaData(:,:,20,1) = (DG%eScaData(:,:,20,1)-DG%Te0_2D*DG%eScaData(:,:,1,1))*inv_tot_ne ! delta_T_para
         DG%eScaData(:,:,21,1) = (DG%eScaData(:,:,21,1)-DG%Te0_2D*DG%eScaData(:,:,1,1))*inv_tot_ne ! delta_T_perp
         DG%eScaData(:,:,6,1) = DG%eScaData(:,:,6,1)/DG%eScaData(:,:,5,1)! Te_tot = Energy_tot/Density_tot
         DG%eScaData(:,:,22,1) = DG%eScaData(:,:,22,1)/DG%eScaData(:,:,5,1)! Te_tot_para
         DG%eScaData(:,:,23,1) = DG%eScaData(:,:,23,1)/DG%eScaData(:,:,5,1)! Te_tot_perp


         ! Get Flux Density by dividing Area in each direction
         do vid=1,4
            do j=1,3
               DG%eVecData(j,:,:,vid,1) = DG%eVecData(j,:,:,vid,1)*DG%inv_Cell_dV  ! Flux Density
            enddo
         enddo

         ! To Normalize CFL condtion values, divide by number of markers (4th iSca)
         vid=5
         do j=1,3
            DG%eVecData(j,:,:,vid,1) = DG%eVecData(j,:,:,vid,1)/(DG%eScaData(:,:,4,1)*DG%Cell_dV)
         enddo

      end subroutine Measure_Electrons_on_Diagnostic_Grid



      subroutine Measure_Ions_on_Diagnostic_Grid(in_P_data, Full_dadgqdgf)

         use precision
         use global_parameters
         use field_array
         use mpi_defs
         use particle_array
         use particle_decomp
         use esi_interface

         implicit none

         real(wp), intent(in) :: in_P_data(3,mi) ! additional particle data from input
         real(wp), intent(in) :: Full_dadgqdgf(3,mi)
         real(wp) energy(mi), KE_para(mi), KE_perp(mi)
         !P varaiables are at particle position
         real(wp) P_s, P_safety_q !safety factor
         real(wp) P_a, P_gq, P_gf

         real(wp) delta_f(mi), f_tot(mi) ! delta_f and f_tot (f_tot=f0+delta_f)

         real(wp) target_gq ! projected gq at nearby poloidal plane following the B0 field
         integer :: ori_sid !sid of original full grid
         real(wp) ori_inv_ds,ori_inv_dgf, P_inv_ori_dgq ! inv_ds of original full grid
         real(wp) ff, ws0,ws1, wgq0,wgq1, wgf,ratio ! ratios along each dimension
         real(wp) target_zeta
         integer :: i,m, j,k, aid, sid, vid, gqid, qq, ss
         integer ::pbid,vxid,vyid,fid
         real(wp):: norm_vpara, norm_vperp, rvx,rvy ! norm_vpara=vpara/vth, norm_vperp=vperp/vth, where vth=sqrt(T/m)
         real(wp),allocatable, dimension(:,:,:) :: tmp_PVS_fi
         integer :: length, istatus(MPI_STATUS_SIZE),ierror
         integer :: larmor
         real(wp) tmp_Scalar_Data(DG%dim_gq,DG%dim_a,DG%num_iSD,0:1), tmp_Vector_Data(3,DG%dim_gq,DG%dim_a,DG%num_iVD,0:1)
         real(wp) ::inv_tot_ni(DG%dim_gq, DG%dim_a)

         real(wp) :: uvec_N_cov(3), uvec_B_cov(3), uvec_BxN_cov(3) ! unit vectors in total B direction, Normal (gradFL), and binormal
         real(wp) :: uvec_N_cntv(3), uvec_B_cntv(3), uvec_BxN_cntv(3) ! unit vectors in total B direction, Normal (gradFL), and binormal

         if(acc_newer_D_zion) then
            !$acc update self(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
            acc_newer_D_zion=.false.
            acc_newer_H_zion=.false.
         endif

         DG%iScaData=0
         DG%iVecData=0
         do pbid=0,nPVS
            PVS(pbid)%fi = 0
         enddo


         ori_inv_ds = real(mpsi)/(s1-s0)
         ori_inv_dgf= 1.0/(zetamax-zetamin)

         energy(:) = in_P_data(2,:)
         KE_para = 0.5*(zion(4,1:mi)*BB(1:mi)-backu(1:mi))**2.0
         KE_perp = zion(6,1:mi)*zion(6,1:mi)*BB(1:mi)

         if(Flag_Load_and_WeightDef==0)then
            delta_f(:) = zion0(5,1:mi)*iNrNm
            f_tot(:)   = iNrNm
         else
            delta_f(:) = zion0(5,1:mi)
            f_tot(:)   = zion0(6,1:mi)
         endif

         do m=1,mi

            if (Flag_w_Threshold==1 .and. (zion0(5,m)<w_Thold_min .or. zion0(5,m)>w_Thold_max)) cycle

            P_a  = zion0(1,m)
            P_gq = modulo(zion0(2,m),2*pi)
            P_gf = modulo(zion0(3,m),2*pi)
            P_s  = zion0(7,m)

            ori_sid = (P_s-s0)*ori_inv_ds ! s index of original full grid
            ff  = (P_s-s_grid(ori_sid))*ori_inv_ds
            P_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)
            P_inv_ori_dgq = 1.0/(ff*deltat(ori_sid+1)+(1.0-ff)*deltat(ori_sid))

            sid = max(1,min(int((P_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
            ws1 = (P_s-DG%s_1D(sid))*DG%inv_ds
            ws0 = 1.0-ws1

            call Cal_unit_vectors_tangential_normal_binormal_to_FluxSurface(P_a,P_gq,P_gf,P_s, uvec_B_cov, uvec_B_cntv, uvec_N_cov, uvec_N_cntv, uvec_BxN_cov, uvec_BxN_cntv)

            ! Scatter data into neighboring DG nodes
            do k=0,1
               if(k==0) target_zeta =zetamin
               if(k==1) target_zeta =zetamax

               wgf = (1.0-abs(P_gf-target_zeta)/deltaz)

               target_gq = modulo(P_gq + (target_zeta-P_gf)/P_safety_q, 2.0*pi)
               gqid = (target_gq)*DG%inv_dgq +1
               wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
               wgq0 = 1.0-wgq1

               do j=1,4
                  select case(j)
                     case(1)
                        qq=gqid; ss=sid
                        ratio=ws0*wgq0*wgf
                     case(2)
                        qq=gqid+1; ss=sid
                        ratio=ws0*wgq1*wgf
                     case(3)
                        qq=gqid; ss=sid+1
                        ratio=ws1*wgq0*wgf
                     case(4)
                        qq=gqid+1; ss=sid+1
                        ratio=ws1*wgq1*wgf
                  end select

                  vid=1; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m) ! Density
                  vid=2; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*(2.0/3.0*energy(m)-backT(m))*delta_f(m) !Temperature definition 1
                  vid=3; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*(2.0/3.0*energy(m))*delta_f(m) !
                  vid=4; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio ! Marker Distribution
                  vid=5; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
                  vid=6; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*f_tot(m)*(2.0/3.0*energy(m)) ! Integrating f_tot*energy, so Te_tot

                  vid=22; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*(2.0*KE_para(m))*delta_f(m) !
                  vid=23; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*(KE_perp(m))*delta_f(m) !
                  vid=24; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*f_tot(m)*(2.0*KE_para(m)) ! Ti_tot_para
                  vid=25; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*f_tot(m)*(KE_perp(m)) ! Ti_tot_perp

                  if(ilabels(m)==-1) then
                     vid=21; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio !  Negative df Markers from Walls
                  endif

                  ! Volumetric delta-f Partcle Flux (=df*v_tot) (#/(m^2*s)) along tangential, normal, binormal directions, respectively
                  vid=9; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_B_cov)
                  vid=10;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_N_cov)
                  vid=11;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_BxN_cov)

                  vid=12;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_B_cntv(1))*uvec_B_cov)
                  vid=13;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_N_cntv(1))*uvec_N_cov)
                  vid=14;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_BxN_cntv(1))*uvec_BxN_cov)

                  ! Volumetric delta-f Energy Flux (=df*KE*v_tot) (#/(m^2*s)) along tangential, normal, binormal directions, respectively
                  vid=15;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_B_cov)
                  vid=16;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_N_cov)
                  vid=17;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_BxN_cov)

                  vid=18;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_B_cntv(1))*uvec_B_cov)
                  vid=19;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_N_cntv(1))*uvec_N_cov)
                  vid=20;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_BxN_cntv(1))*uvec_BxN_cov)

                  vid=1 ! Turbulent Particle Flux
                  ! DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*(delta_f(m)*ghq(m)) ! Radial Flux
                  ! DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*(delta_f(m)*Ta(m)) ! theta Flux
                  ! DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*(delta_f(m)*dvi_gf(m)) ! zeta Flux
                  DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*delta_f(m)*Full_dadgqdgf(1,m)*ii_inv_grad_a(m) ! Radial Flux
                  DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*delta_f(m)*Full_dadgqdgf(2,m)*ii_inv_grad_gq(m) ! theta Flux
                  DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*delta_f(m)*Full_dadgqdgf(3,m)*ii_inv_grad_gf(m) ! zeta Flux


                  vid=2 ! Turbulent Particle Flux due to ExB only
                  DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*(delta_f(m)*(ghq(m)-vpara_dvi_a(m)))*ii_inv_grad_a(m) ! Radial Flux
                  DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*(delta_f(m)*(Ta(m)-vpara_dvi_gq(m)))*ii_inv_grad_gq(m) ! theta Flux
                  DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*(delta_f(m)*(dvi_gf(m)-vpara_dvi_gf(m)))*ii_inv_grad_gf(m)! zeta Flux

                  vid=3 ! Turbulent Energy Flux
                  ! DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*ghq(m)) ! Radial Flux
                  ! DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Ta(m)) ! theta Flux
                  ! DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*dvi_gf(m)) ! zeta Flux
                  DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Full_dadgqdgf(1,m))*ii_inv_grad_a(m) ! Radial Flux
                  DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Full_dadgqdgf(2,m))*ii_inv_grad_gq(m) ! theta Flux
                  DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Full_dadgqdgf(3,m))*ii_inv_grad_gf(m) ! zeta Flux



                  vid=4 ! Turbulent Energy Flux due to ExB only
                  DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*(ghq(m)-vpara_dvi_a(m)))*ii_inv_grad_a(m) ! Radial Flux
                  DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*(Ta(m)-vpara_dvi_gq(m)))*ii_inv_grad_gq(m) ! theta Flux
                  DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*(dvi_gf(m)-vpara_dvi_gf(m)))*ii_inv_grad_gf(m)! zeta Flux

                  vid=5 ! To measure CFL condition (v*dt/dx)
                  DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*abs(tstep*((Tc/backT(m))**dsda_power*Full_dadgqdgf(1,m))*ori_inv_ds)
                  DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*abs(tstep*(Full_dadgqdgf(2,m)-Full_dadgqdgf(3,m)/P_safety_q)*P_inv_ori_dgq) !time derivative of target_gq=gq-gf/safety_q
                  DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*abs(tstep*Full_dadgqdgf(3,m)*ori_inv_dgf)
               enddo
            enddo

            ! Velocity Space Measurement
            do pbid=0,nPVS
               target_gq = modulo(P_gq + (PVS(pbid)%gf-P_gf)/P_safety_q, 2.0*pi)

               if( abs(P_gf-PVS(pbid)%gf)<0.5*PVS(pbid)%L_gf &
               .and. abs(P_a-PVS(pbid)%a) <0.5*PVS(pbid)%L_a &
               .and. min(abs(target_gq-PVS(pbid)%gq), 2.0*pi-abs(target_gq-PVS(pbid)%gq)) < 0.5*PVS(pbid)%L_gq) then

                  norm_vpara = zion0(4,m)*BB(m)*sqrt(PVS(pbid)%inv_Ti0) ! norm_vpara = vpara/sqrt(T/m)
                  norm_vperp = zion(6,m)*sqrt(2.0*BB(m))*sqrt(PVS(pbid)%inv_Ti0) ! norm_vperp = sqrt(mu)*sqrt(2B/m)/sqrt(T/m)

                  vxid = floor(norm_vpara*PVS(pbid)%inv_dv)
                  vyid = floor(norm_vperp*PVS(pbid)%inv_dv)
                  rvx  = (norm_vpara-PVS(pbid)%vpara_1D(vxid))*PVS(pbid)%inv_dv
                  rvy  = (norm_vperp-PVS(pbid)%vperp_1D(vyid))*PVS(pbid)%inv_dv

                  ! Do not calculate any contribution from "too" outside of our probe boudnary
                  ! "too" outside means that even over the ghost cells around the probe boundary
                  if(vxid< -PVS(pbid)%nv-1) then
                     continue
                  elseif(vxid>=PVS(pbid)%nv+1) then
                     continue
                  elseif(vyid>=PVS(pbid)%nv+1) then
                     continue
                  endif

                  if(vxid< -PVS(pbid)%nv) then
                     vxid= -PVS(pbid)%nv
                     rvx=0.0
                  elseif (vxid>= PVS(pbid)%nv) then
                     vxid=  PVS(pbid)%nv-1
                     rvx=1.0
                  endif

                  if(vyid>= PVS(pbid)%nv) then
                     vyid=  PVS(pbid)%nv-1
                     rvy=1.0
                  endif

                  do j=1,4
                     select case(j)
                        case(1)
                           qq=vxid; ss=vyid
                           ratio=(1.0-rvx)*(1.0-rvy)
                        case(2)
                           qq=vxid+1; ss=vyid
                           ratio=(rvx)*(1.0-rvy)
                        case(3)
                           qq=vxid; ss=vyid+1
                           ratio=(1.0-rvx)*(rvy)
                        case(4)
                           qq=vxid+1; ss=vyid+1
                           ratio=(rvx)*(rvy)
                     end select

                     fid=1; PVS(pbid)%fi(qq,ss,fid) = PVS(pbid)%fi(qq,ss,fid) + ratio*delta_f(m) ! Density
                     fid=2; PVS(pbid)%fi(qq,ss,fid) = PVS(pbid)%fi(qq,ss,fid) + ratio*delta_f(m)*(2.0/3.0*energy(m)-backT(m)) ! Temperature definition
                     fid=3; PVS(pbid)%fi(qq,ss,fid) = PVS(pbid)%fi(qq,ss,fid) + ratio ! Marker Density
                     fid=4; PVS(pbid)%fi(qq,ss,fid) = PVS(pbid)%fi(qq,ss,fid) + ratio*f_tot(m) ! f_tot
                     fid=5; PVS(pbid)%fi(qq,ss,fid) = PVS(pbid)%fi(qq,ss,fid) + ratio*f_tot(m)*(2.0/3.0*energy(m)) ! f_tot*energy
                  enddo
               endif
            enddo
         enddo

         !$acc update self (zion_s_on_GyroOrbit(:,mi),zion_gq_on_GyroOrbit(:,mi))
         do m=1,mi
            if (Flag_w_Threshold==1 .and. (zion0(5,m)<w_Thold_min .or. zion0(5,m)>w_Thold_max)) cycle

            do larmor=1,larmor_num
               P_a  = zion0(1,m)
               P_gq = modulo(zion_gq_on_GyroOrbit(larmor,m),2*pi)
               P_gf = modulo(zion0(3,m),2*pi)
               P_s  = zion_s_on_GyroOrbit(larmor,m)

               ! ori_sid = (P_s-s0)*ori_inv_ds ! s index of original full grid
               ori_sid = max(0,min(int((P_s-s0)*ori_inv_ds),mpsi-1)) ! To avoid segmentation fault
               ff  = (P_s-s_grid(ori_sid))*ori_inv_ds
               P_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)
               P_inv_ori_dgq = 1.0/(ff*deltat(ori_sid+1)+(1.0-ff)*deltat(ori_sid))

               sid = max(1,min(int((P_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
               ws1 = (P_s-DG%s_1D(sid))*DG%inv_ds
               ws0 = 1.0-ws1

               do k=0,1
                  if(k==0) target_zeta =zetamin
                  if(k==1) target_zeta =zetamax

                  wgf = (1.0-abs(P_gf-target_zeta)/deltaz)

                  target_gq = modulo(P_gq + (target_zeta-P_gf)/P_safety_q, 2.0*pi)
                  gqid = (target_gq)*DG%inv_dgq +1
                  wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
                  wgq0 = 1.0-wgq1

                  do j=1,4
                     select case(j)
                        case(1)
                           qq=gqid; ss=sid
                           ratio=ws0*wgq0*wgf
                        case(2)
                           qq=gqid+1; ss=sid
                           ratio=ws0*wgq1*wgf
                        case(3)
                           qq=gqid; ss=sid+1
                           ratio=ws1*wgq0*wgf
                        case(4)
                           qq=gqid+1; ss=sid+1
                           ratio=ws1*wgq1*wgf
                     end select

                     ! Scalar Data
                     vid=7; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*inv_larmor_num*delta_f(m) ! Density
                     vid=8; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*inv_larmor_num*(2.0/3.0*energy(m))*delta_f(m) !Temperature definition 1

                     ! Vector Data
                     vid=6 ! Turbulent Particle Flux
                     DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*inv_larmor_num*(delta_f(m)*ghq(m))*ii_inv_grad_a(m) ! Radial Flux
                     DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*inv_larmor_num*(delta_f(m)*Ta(m))*ii_inv_grad_gq(m) ! theta Flux
                     DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*inv_larmor_num*(delta_f(m)*dvi_gf(m))*ii_inv_grad_gf(m) ! zeta Flux

                     vid=7 ! Turbulent Particle Flux due to ExB only
                     DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*inv_larmor_num*(delta_f(m)*(ghq(m)-vpara_dvi_a(m)))*ii_inv_grad_a(m) ! Radial Flux
                     DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*inv_larmor_num*(delta_f(m)*(Ta(m)-vpara_dvi_gq(m)))*ii_inv_grad_gq(m) ! theta Flux
                     DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*inv_larmor_num*(delta_f(m)*(dvi_gf(m)-vpara_dvi_gf(m)) )*ii_inv_grad_gf(m)! zeta Flux

                     vid=8 ! Turbulent Energy Flux
                     DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*energy(m)*inv_larmor_num*(delta_f(m)*ghq(m))*ii_inv_grad_a(m) ! Radial Flux
                     DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*energy(m)*inv_larmor_num*(delta_f(m)*Ta(m))*ii_inv_grad_gq(m) ! theta Flux
                     DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*energy(m)*inv_larmor_num*(delta_f(m)*dvi_gf(m))*ii_inv_grad_gf(m) ! zeta Flux

                     vid=9 ! Turbulent Energy Flux due to ExB only
                     DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*energy(m)*inv_larmor_num*(delta_f(m)*(ghq(m)-vpara_dvi_a(m)))*ii_inv_grad_a(m) ! Radial Flux
                     DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*energy(m)*inv_larmor_num*(delta_f(m)*(Ta(m)-vpara_dvi_gq(m)))*ii_inv_grad_gq(m) ! theta Flux
                     DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*energy(m)*inv_larmor_num*(delta_f(m)*(dvi_gf(m)-vpara_dvi_gf(m)) )*ii_inv_grad_gf(m)! zeta Flux
                  enddo
               enddo
            enddo
         enddo

         if(npartdom>1) then
            tmp_Scalar_Data = DG%iScaData
            tmp_Vector_Data = DG%iVecData

            length = DG%dim_gq*DG%dim_a*DG%num_iSD*2
            call MPI_ALLREDUCE(tmp_Scalar_Data, DG%iScaData, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)

            length = 3*DG%dim_gq*DG%dim_a*DG%num_iVD*2
            call MPI_ALLREDUCE(tmp_Vector_Data, DG%iVecData, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)

            allocate(tmp_PVS_fi, mold=PVS(0)%fi)
            do pbid=0,nPVS
               tmp_PVS_fi = PVS(pbid)%fi
               length = PVS(pbid)%dim_vpara*PVS(pbid)%dim_vperp*PVS(pbid)%num_of_fi
               call MPI_ALLREDUCE(tmp_PVS_fi, PVS(pbid)%fi, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
            enddo
            deallocate(tmp_PVS_fi)
         endif


         ! Match poloidal end cell
         ! (gq=0 and gq=2pi is the same position physically
         ! although they have different gqid)
         do sid=1,DG%dim_a
            DG%iScaData(1,sid,:,:) = DG%iScaData(1,sid,:,:) + DG%iScaData(DG%dim_gq,sid,:,:)
            DG%iScaData(DG%dim_gq,sid,:,:) = DG%iScaData(1,sid,:,:)

            DG%iVecData(:,1,sid,:,:) = DG%iVecData(:,1,sid,:,:) + DG%iVecData(:,DG%dim_gq,sid,:,:)
            DG%iVecData(:,DG%dim_gq,sid,:,:) = DG%iVecData(:,1,sid,:,:)
         enddo


         ! toroidal end cell
         length = DG%dim_gq*DG%dim_a*DG%num_iSD
         call MPI_SENDRECV(DG%iScaData(:,:,:,0),length,mpi_Rsize,left_pe,1, tmp_Scalar_Data(:,:,:,1),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         DG%iScaData(:,:,:,1) = DG%iScaData(:,:,:,1) + tmp_Scalar_Data(:,:,:,1)



         length = 3*DG%dim_gq*DG%dim_a*DG%num_iVD
         call MPI_SENDRECV(DG%iVecData(:,:,:,:,0),length,mpi_Rsize,left_pe,1, tmp_Vector_Data(:,:,:,:,1),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         DG%iVecData(:,:,:,:,1) = DG%iVecData(:,:,:,:,1) + tmp_Vector_Data(:,:,:,:,1)


         ! Do normalization and post-processing
         ! Divide by Cell Volume to get density quantities
         do vid=1,DG%num_iSD
            DG%iScaData(:,:,vid,1) = DG%iScaData(:,:,vid,1)*DG%inv_Cell_dV
         enddo


         ! Calculate hybrid particle and energy density
         DG%w_dni      = abs(DG%iScaData(:,:,1,1))/(abs(DG%iScaData(:,:,1,1))+abs(DG%iScaData(:,:,5,1))) ! wdni=abs(snap.DG.dni)./(abs(snap.DG.ni)+abs(snap.DG.dni));
         DG%ni_hybrid  = DG%w_dni*DG%iScaData(:,:,5,1) + (1.0-DG%w_dni)*(DG%iScaData(:,:,1,1)+DG%n0_2D) ! (n_hybrid = w_dni.*snap.DG.ni+(1.0-w_dni).*(snap.DG.dni+snap.n0_2D))
         DG%Ti_hybrid  = DG%w_dni*DG%iScaData(:,:,6,1) + (1.0-DG%w_dni)*(DG%iScaData(:,:,3,1)+DG%n0_2D*DG%Ti0_2D) ! (hybrid energy_density, then will be divded by n_hybrid to get T_hybrid)
         DG%Ti_hybrid  = DG%Ti_hybrid / DG%ni_hybrid

         ! Calculate Temperatures
         DG%iScaData(:,:,2,1) = DG%iScaData(:,:,2,1)/DG%n0_2D ! Convert Energy to Temperature
         ! Consistent Temperature by considering density fluctuation in the denominator
         DG%iScaData(:,:,3,1) = (DG%iScaData(:,:,3,1)-DG%Ti0_2D*DG%iScaData(:,:,1,1))/(DG%n0_2D+DG%iScaData(:,:,1,1)) ! dTi_consst using consistent varying ion density
         DG%iScaData(:,:,22,1) = (DG%iScaData(:,:,22,1)-DG%Ti0_2D*DG%iScaData(:,:,1,1))/(DG%n0_2D+DG%iScaData(:,:,1,1)) ! dTi_para_consst
         DG%iScaData(:,:,23,1) = (DG%iScaData(:,:,23,1)-DG%Ti0_2D*DG%iScaData(:,:,1,1))/(DG%n0_2D+DG%iScaData(:,:,1,1)) ! dTi_perp_consst
         DG%iScaData(:,:,8,1) = (DG%iScaData(:,:,8,1)-DG%Ti0_2D*DG%iScaData(:,:,7,1))/(DG%n0_2D+DG%iScaData(:,:,7,1)) ! dTi_consst using consistent ion density on Larmor points
         DG%iScaData(:,:,6,1) = DG%iScaData(:,:,6,1)/(DG%iScaData(:,:,5,1))! Ti_tot using f_tot on GC
         DG%iScaData(:,:,24,1) = DG%iScaData(:,:,24,1)/(DG%iScaData(:,:,5,1))! Ti_tot_para
         DG%iScaData(:,:,25,1) = DG%iScaData(:,:,25,1)/(DG%iScaData(:,:,5,1))! Ti_tot_perp

         ! Get Flux Density by dividing Cell Volume
         do vid=1,4
            do j=1,3
               DG%iVecData(j,:,:,vid,1) = DG%iVecData(j,:,:,vid,1)*DG%inv_Cell_dV  ! Flux Density
            enddo
         enddo

         ! To Normalize CFL condtion values, divide by number of markers (4th iSca)
         vid=5
         do j=1,3
            DG%iVecData(j,:,:,vid,1) = DG%iVecData(j,:,:,vid,1)/(DG%iScaData(:,:,4,1)*DG%Cell_dV)
         enddo

         ! Get Flux Density (using the positions on Larmor Radius) by dividing Cell Volume
         do vid=6,9
            do j=1,3
               DG%iVecData(j,:,:,vid,1) = DG%iVecData(j,:,:,vid,1)*DG%inv_Cell_dV  ! Flux Density
            enddo
         enddo

      end subroutine Measure_Ions_on_Diagnostic_Grid

      subroutine Measure_and_Cumulate_eZonData_on_DG(in_P_data, Full_da)
         use precision
         use global_parameters
         use field_array
         use mpi_defs
         use particle_array
         use particle_decomp
         use esi_interface

         implicit none

         real(wp), intent(in) :: in_P_data(3,me) ! additional particle data from input
         real(wp), intent(in) :: Full_da(me)
         real(wp) ws0,ws1,tmp,P_s
         real(wp) eData1D(DG%num_eZD_1D,DG%dim_a)
         real(wp), allocatable, dimension(:,:) ::ptl_data
         integer :: sid,vid,m,size,ierr

         !$acc declare create(ptl_data)

         if(acc_newer_H_zeon) then
            !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_H_zeon=.false.
            acc_newer_D_zeon=.false.
         endif


         allocate(ptl_data(DG%num_eZD_1D,me))
         eData1D=0.0
         ! in_P_data(2,:)= energy(:)

         !$acc parallel loop present(zeon0,ptl_data,in_P_data,Full_da,ee_inv_grad_a,ghq,vpara_dve_a)
         do m=1,me
            ! Fluxes by Full Velocity
            ptl_data(1,m) = zeon0(5,m)*Full_da(m)*ee_inv_grad_a(m)  ! Particle Flux
            ptl_data(2,m) = ptl_data(1,m)*in_P_data(2,m) ! Energy Flux

            ! Fluxes by ExB contribution
            ptl_data(3,m) = zeon0(5,m)*(ghq(m)-vpara_dve_a(m))*ee_inv_grad_a(m) ! Particle Flux by ExB only
            ptl_data(4,m) = ptl_data(3,m)*in_P_data(2,m) ! Energy Flux by ExB only

            ! Fluxes by Streaming contribution
            ptl_data(5,m) = zeon0(5,m)*(vpara_dve_a(m))*ee_inv_grad_a(m) ! Particle Flux by streaming only
            ptl_data(6,m) = ptl_data(5,m)*in_P_data(2,m) ! Energy Flux by streaming only

            ! Particle Density and Energy Density
            ptl_data(7,m) = zeon0(5,m) ! Particle Density
            ptl_data(8,m) = zeon0(5,m)*in_P_data(2,m) ! Energy Density
         enddo
         if(Flag_Load_and_WeightDef==0)then
            !$acc parallel loop present(ptl_data)
            do m=1,me
               ptl_data(:,m) = eNrNm*ptl_data(:,m)
            enddo
         endif


         !$acc parallel loop gang vector present(zeon,zeon0,ptl_data)
         do m=1,me
            P_s  = zeon0(7,m)
            sid = max(1,min(int((P_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
            ws1 = (P_s-DG%s_1D(sid))*DG%inv_ds
            ws0 = 1.0-ws1

            !$acc loop seq
            do vid=1,DG%num_eZD_1D
               !$acc atomic update
               eData1D(vid,sid  ) = eData1D(vid,sid  ) + ws0*ptl_data(vid,m)
               !$acc atomic update
               eData1D(vid,sid+1) = eData1D(vid,sid+1) + ws1*ptl_data(vid,m)
            enddo
         enddo


         ! Cumulate Zonal Data
         DG%eZonData_1D = eData1D
         DG%cum_eZonData_1D = DG%cum_eZonData_1D + eData1D


         ! Communicate between Ranks and Normalize the data, when they are needed
         ! Not every time step, to reduce computational cost
         if(mod(istep,msnap)==0) then

            eData1D = DG%eZonData_1D
            size = DG%dim_a*DG%num_eZD_1D
            call MPI_ALLREDUCE(eData1D, DG%eZonData_1D, size, mpi_Rsize, MPI_SUM, MPI_COMM_WORLD, ierr)

            eData1D = DG%cum_eZonData_1D
            size = DG%dim_a*DG%num_eZD_1D
            call MPI_ALLREDUCE(eData1D, DG%cum_eZonData_1D, size, mpi_Rsize, MPI_SUM, MPI_COMM_WORLD, ierr)

            ! Normalize Data by dividing Anulus volume
            do vid=1,DG%num_eZD_1D
               DG%eZonData_1D(vid,:) = DG%eZonData_1D(vid,:)*DG%inv_Anulus_dV
               DG%cum_eZonData_1D(vid,:) = DG%cum_eZonData_1D(vid,:)*DG%inv_Anulus_dV
            enddo
         endif

         deallocate(ptl_data)

      end subroutine Measure_and_Cumulate_eZonData_on_DG

      subroutine Measure_and_Cumulate_iZonData_on_DG(in_P_data, Full_da)
         use precision
         use global_parameters
         use field_array
         use mpi_defs
         use particle_array
         use particle_decomp
         use esi_interface

         implicit none

         real(wp), intent(in) :: in_P_data(3,mi) ! additional particle data from input
         real(wp), intent(in) :: Full_da(mi)
         real(wp) ws0,ws1,tmp,P_s
         real(wp) iData1D(DG%num_iZD_1D,DG%dim_a)
         real(wp), allocatable, dimension(:,:) ::ptl_data
         integer :: sid,vid,m,size,ierr

         !$acc declare create(ptl_data)

         if(acc_newer_H_zion) then
            !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
            acc_newer_H_zion=.false.
            acc_newer_D_zion=.false.
         endif


         allocate(ptl_data(DG%num_iZD_1D,mi))
         iData1D=0.0

         ! in_P_data(2,:)= energy(:)

         !$acc parallel loop present(zion0,ptl_data,in_P_data,Full_da,ii_inv_grad_a,ghq,vpara_dvi_a)
         do m=1,mi
            ! Fluxes by Full Velocity
            ptl_data(1,m) = zion0(5,m)*Full_da(m)*ii_inv_grad_a(m)  ! Particle Flux
            ptl_data(2,m) = ptl_data(1,m)*in_P_data(2,m) ! Energy Flux

            ! Fluxes by ExB contribution
            ptl_data(3,m) = zion0(5,m)*(ghq(m)-vpara_dvi_a(m))*ii_inv_grad_a(m) ! Particle Flux by ExB only
            ptl_data(4,m) = ptl_data(3,m)*in_P_data(2,m) ! Energy Flux by ExB only

            ! Fluxes by Streaming contribution
            ptl_data(5,m) = zion0(5,m)*(vpara_dvi_a(m))*ii_inv_grad_a(m) ! Particle Flux by streaming only
            ptl_data(6,m) = ptl_data(5,m)*in_P_data(2,m) ! Energy Flux by streaming only

            ! Particle Density and Energy Density
            ptl_data(7,m) = zion0(5,m) ! Particle Density
            ptl_data(8,m) = zion0(5,m)*in_P_data(2,m) ! Energy Density
         enddo
         if(Flag_Load_and_WeightDef==0)then
            !$acc parallel loop present(ptl_data)
            do m=1,me
               ptl_data(:,m) = iNrNm*ptl_data(:,m)
            enddo
         endif

         !$acc parallel loop gang vector present(zion0,ptl_data)
         do m=1,mi
            P_s  = zion0(7,m)
            sid = max(1,min(int((P_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
            ws1 = (P_s-DG%s_1D(sid))*DG%inv_ds
            ws0 = 1.0-ws1

            !$acc loop seq
            do vid=1,DG%num_iZD_1D
               !$acc atomic update
               iData1D(vid,sid  ) = iData1D(vid,sid  ) + ws0*ptl_data(vid,m)
               !$acc atomic update
               iData1D(vid,sid+1) = iData1D(vid,sid+1) + ws1*ptl_data(vid,m)
            enddo
         enddo


         ! Cumulate Zonal Data
         DG%iZonData_1D = iData1D
         DG%cum_iZonData_1D = DG%cum_iZonData_1D + iData1D


         ! Communicate between Ranks and Normalize the data, when they are needed
         ! Not every time step, to reduce computational cost
         if(mod(istep,msnap)==0) then

            iData1D = DG%iZonData_1D
            size = DG%dim_a*DG%num_iZD_1D
            call MPI_ALLREDUCE(iData1D, DG%iZonData_1D, size, mpi_Rsize, MPI_SUM, MPI_COMM_WORLD, ierr)

            iData1D = DG%cum_iZonData_1D
            size = DG%dim_a*DG%num_iZD_1D
            call MPI_ALLREDUCE(iData1D, DG%cum_iZonData_1D, size, mpi_Rsize, MPI_SUM, MPI_COMM_WORLD, ierr)

            ! Normalize Data by dividing Cell Volume
            do vid=1,DG%num_iZD_1D
               DG%iZonData_1D(vid,:) = DG%iZonData_1D(vid,:)*DG%inv_Anulus_dV
               DG%cum_iZonData_1D(vid,:) = DG%cum_iZonData_1D(vid,:)*DG%inv_Anulus_dV
            enddo
         endif

         deallocate(ptl_data)

      end subroutine Measure_and_Cumulate_iZonData_on_DG

      subroutine SetZeros_cum_Radial_FluxData_a1D_on_DG
         implicit none

         DG%cum_eZonData_1D = 0.0
         DG%cum_iZonData_1D = 0.0

      end subroutine

      subroutine Interpolate_OriGridData_into_DG(oriData_1D, out_IntpData_2D)

         use precision
         use global_parameters
         use field_array
         implicit none

         real(wp), dimension(mgrid), intent(in) :: oriData_1D
         real(wp), dimension(DG%dim_gq,DG%dim_a), intent(out) :: out_IntpData_2D

         integer i,j,ij,jt, ori_aid

         real(wp) ori_inv_dgq(0:mpsi), gqdum, wt

         ori_inv_dgq=1.0/deltat

         do i=1,DG%dim_a
            do j=1,DG%dim_gq

               ori_aid = DG%ori_aid(i)

               gqdum = modulo(DG%gq_1D(j) - qtinv(ori_aid)*zetamax, 2.0*pi) * ori_inv_dgq(ori_aid)
               jt=int(gqdum);  wt=gqdum-jt
               ij=igrid(ori_aid)+jt
               out_IntpData_2D(j,i) =(1.0-wt)*oriData_1D(ij) + wt*oriData_1D(ij+1)
            enddo
         enddo

         ! Periodic in theta direction (to make sure they are exactly the same)
         out_IntpData_2D(DG%dim_gq,:) = out_IntpData_2D(1,:)

      end subroutine Interpolate_OriGridData_into_DG



      subroutine Smooth_DG_Scalar_Vector_Data(num_smooth)

         use precision
         use global_parameters
         implicit none

         integer, intent(in) :: num_smooth
         integer :: smid, aid, gqid, Na,Ngq, vid,k
         real(wp) :: VV(0:DG%dim_gq+1, 0:DG%dim_a+1),  target_Data2D(0:DG%dim_gq+1, 0:DG%dim_a+1) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
         real(wp) :: SM_gq_Coeff(-1:1,0:DG%dim_gq+1, 0:DG%dim_a+1), SM_a_Coeff(-1:1,0:DG%dim_gq+1, 0:DG%dim_a+1)
         real(wp) :: tot_ne(DG%dim_gq, DG%dim_a)

         Na  = DG%dim_a
         Ngq = DG%dim_gq


         target_Data2D=0

         ! Get Smooth Coeffs which conserves Volume Weighted Integration
         call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dV, SM_a_Coeff, SM_gq_Coeff)
         ! =======================================
         ! Smooth Electron Scalar Data by conserving Volume Weighted Integration
         ! =======================================

         ! Convert temperature to EnergyDensity before smoothing
         tot_ne = DG%n0_2D + DG%eScaData(:,:,1,1) + int(ifelectron/2)*DG%dne_adiab ! using df or dh
         DG%eScaData(:,:,2,1)=DG%eScaData(:,:,2,1)*DG%n0_2D ! dTe which is nonconsistent because of using n0 density
         DG%eScaData(:,:,3,1)=DG%eScaData(:,:,3,1)*tot_ne+DG%Te0_2D*DG%eScaData(:,:,1,1) ! dTe_consst using consistent varying electron density
         DG%eScaData(:,:,20,1)=DG%eScaData(:,:,20,1)*tot_ne+DG%Te0_2D*DG%eScaData(:,:,1,1) ! dTe_consst_para
         DG%eScaData(:,:,21,1)=DG%eScaData(:,:,21,1)*tot_ne+DG%Te0_2D*DG%eScaData(:,:,1,1) ! dTe_consst_perp
         DG%eScaData(:,:,6,1)=DG%eScaData(:,:,6,1)*DG%eScaData(:,:,5,1) ! using full-f, Energy_tot=Te_tot*ne_tot
         DG%eScaData(:,:,22,1)=DG%eScaData(:,:,22,1)*DG%eScaData(:,:,5,1) ! using full-f, Te_tot_para
         DG%eScaData(:,:,23,1)=DG%eScaData(:,:,23,1)*DG%eScaData(:,:,5,1) ! using full-f, Te_tot_perp

         ! In the case of electrons, dne_adiab should be smoothed as well.
         target_Data2D(1:Ngq,1:Na) = DG%dne_adiab
         do smid=1,num_smooth
            ! Set ghost cell periodic in theta direction
            target_Data2D(0,:)     = target_Data2D(Ngq-1,:)
            target_Data2D(Ngq+1,:) = target_Data2D(2,:)
            ! Smooth in theta direction
            target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
            ! Smooth in radial direction
            target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
         enddo
         DG%dne_adiab = target_Data2D(1:Ngq,1:Na)

         do vid=1,DG%num_eSD
            target_Data2D(1:Ngq,1:Na) = DG%eScaData(:,:,vid,1)
            do smid=1,num_smooth
               ! Set ghost cell periodic in theta direction
               target_Data2D(0,:)     = target_Data2D(Ngq-1,:)
               target_Data2D(Ngq+1,:) = target_Data2D(2,:)
               ! Smooth in theta direction
               target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
               ! Smooth in radial direction
               target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
            enddo
            DG%eScaData(:,:,vid,1) = target_Data2D(1:Ngq,1:Na)
         enddo

         ! Convert smoothed energy density to the temperature using smoothed plasma density
         tot_ne = DG%n0_2D + DG%eScaData(:,:,1,1) + int(ifelectron/2)*DG%dne_adiab ! smoothed total density using df or dh
         DG%eScaData(:,:,2,1)=DG%eScaData(:,:,2,1)/DG%n0_2D ! dTe which is nonconsistent because of using n0 density
         DG%eScaData(:,:,3,1)=(DG%eScaData(:,:,3,1)-DG%Te0_2D*DG%eScaData(:,:,1,1))/tot_ne !delta_T = (2/3*dEnergy - T0*dn)/(n0+dn_h+dne_adiab)]
         DG%eScaData(:,:,20,1)=(DG%eScaData(:,:,20,1)-DG%Te0_2D*DG%eScaData(:,:,1,1))/tot_ne !delta_T_consst_para
         DG%eScaData(:,:,21,1)=(DG%eScaData(:,:,21,1)-DG%Te0_2D*DG%eScaData(:,:,1,1))/tot_ne !delta_T_consst_perp
         DG%eScaData(:,:,6,1)=DG%eScaData(:,:,6,1)/DG%eScaData(:,:,5,1) ! using full-f, Energy_tot=Te_tot*ne_tot
         DG%eScaData(:,:,22,1)=DG%eScaData(:,:,22,1)/DG%eScaData(:,:,5,1) ! using full-f, Te_tot_para
         DG%eScaData(:,:,23,1)=DG%eScaData(:,:,23,1)/DG%eScaData(:,:,5,1) ! using full-f, Te_tot_perp

         ! =======================================
         ! Smooth Ion Scalar Data by conserving Volume Weighted Integration
         ! =======================================

         ! Convert temperature to EnergyDensity before smoothing
         DG%iScaData(:,:,2,1)=DG%iScaData(:,:,2,1)*DG%n0_2D ! dTi which is nonconsistent because of using n0 density
         DG%iScaData(:,:,3,1)=DG%iScaData(:,:,3,1)*(DG%n0_2D+DG%iScaData(:,:,1,1)) +DG%Ti0_2D*DG%iScaData(:,:,1,1)  ! dTi_consst using consistent varying ion density
         DG%iScaData(:,:,22,1)=DG%iScaData(:,:,22,1)*(DG%n0_2D+DG%iScaData(:,:,1,1)) +DG%Ti0_2D*DG%iScaData(:,:,1,1)  ! dTi_consst_para
         DG%iScaData(:,:,23,1)=DG%iScaData(:,:,23,1)*(DG%n0_2D+DG%iScaData(:,:,1,1)) +DG%Ti0_2D*DG%iScaData(:,:,1,1)  ! dTi_consst_perp
         DG%iScaData(:,:,8,1)=DG%iScaData(:,:,8,1)*(DG%n0_2D+DG%iScaData(:,:,7,1)) +DG%Ti0_2D*DG%iScaData(:,:,7,1)  ! dTi_consst using consistent ion density on Larmor points
         DG%iScaData(:,:,6,1)=DG%iScaData(:,:,6,1)*DG%iScaData(:,:,5,1) ! using full-f, Energy_tot=Ti_tot*ni_tot
         DG%iScaData(:,:,24,1)=DG%iScaData(:,:,24,1)*DG%iScaData(:,:,5,1) ! using full-f, Ti_tot_para
         DG%iScaData(:,:,25,1)=DG%iScaData(:,:,25,1)*DG%iScaData(:,:,5,1) ! using full-f, Ti_tot_perp


         do vid=1,DG%num_iSD
            target_Data2D(1:Ngq,1:Na) = DG%iScaData(:,:,vid,1)
            do smid=1,num_smooth
               ! Set ghost cell periodic in theta direction
               target_Data2D(0,:)     = target_Data2D(Ngq-1,:)
               target_Data2D(Ngq+1,:) = target_Data2D(2,:)
               ! Smooth in theta direction
               target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
               ! Smooth in radial direction
               target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
            enddo
            DG%iScaData(:,:,vid,1) = target_Data2D(1:Ngq,1:Na)
         enddo

         ! Convert smoothed energy density to the temperature using smoothed plasma density
         DG%iScaData(:,:,2,1)=DG%iScaData(:,:,2,1)/DG%n0_2D ! dTi which is nonconsistent because of using n0 density
         DG%iScaData(:,:,3,1)=(DG%iScaData(:,:,3,1)-DG%Ti0_2D*DG%iScaData(:,:,1,1))/(DG%n0_2D+DG%iScaData(:,:,1,1)) ! dTi_consst using consistent varying ion density
         DG%iScaData(:,:,22,1)=(DG%iScaData(:,:,22,1)-DG%Ti0_2D*DG%iScaData(:,:,1,1))/(DG%n0_2D+DG%iScaData(:,:,1,1)) ! dTi_consst_para
         DG%iScaData(:,:,23,1)=(DG%iScaData(:,:,23,1)-DG%Ti0_2D*DG%iScaData(:,:,1,1))/(DG%n0_2D+DG%iScaData(:,:,1,1)) ! dTi_consst_perp
         DG%iScaData(:,:,8,1)=(DG%iScaData(:,:,8,1)-DG%Ti0_2D*DG%iScaData(:,:,7,1))/(DG%n0_2D+DG%iScaData(:,:,7,1)) ! dTi_consst using consistent ion density on Larmor points
         DG%iScaData(:,:,6,1)=DG%iScaData(:,:,6,1)/DG%iScaData(:,:,5,1) ! using full-f, Energy_tot=Ti_tot*ni_tot
         DG%iScaData(:,:,24,1)=DG%iScaData(:,:,24,1)/DG%iScaData(:,:,5,1) ! using full-f, Ti_tot_para
         DG%iScaData(:,:,25,1)=DG%iScaData(:,:,25,1)/DG%iScaData(:,:,5,1) ! using full-f, Ti_tot_perp

         ! =======================================
         ! Smooth Electron Vector Data by conserving Surface Area Weighted Integration in each direction
         ! =======================================
         do vid=1,DG%num_eVD
            do k=1,3
               select case(k)
                  case(1)
                     call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dS_a, SM_a_Coeff, SM_gq_Coeff)
                  case(2)
                     call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dS_gq, SM_a_Coeff, SM_gq_Coeff)
                  case(3)
                     call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dS_gf, SM_a_Coeff, SM_gq_Coeff)
               end select

               target_Data2D(1:Ngq,1:Na) = DG%eVecData(k,:,:,vid,1)
               do smid=1,num_smooth
                  ! Set ghost cell periodic in theta direction
                  target_Data2D(0,:)     = target_Data2D(Ngq-1,:)
                  target_Data2D(Ngq+1,:) = target_Data2D(2,:)
                  ! Smooth in theta direction
                  target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
                  ! Smooth in radial direction
                  target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
               enddo
               DG%eVecData(k,:,:,vid,1) = target_Data2D(1:Ngq,1:Na)
            enddo
         enddo

         ! =======================================
         ! Smooth Ion Vector Data by conserving Surface Area Weighted Integration in each direction
         ! =======================================
         do vid=1,DG%num_iVD
            do k=1,3
               select case(k)
               case(1)
                  call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dS_a, SM_a_Coeff, SM_gq_Coeff)
               case(2)
                  call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dS_gq, SM_a_Coeff, SM_gq_Coeff)
               case(3)
                  call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dS_gf, SM_a_Coeff, SM_gq_Coeff)
               end select

               target_Data2D(1:Ngq,1:Na) = DG%iVecData(k,:,:,vid,1)
               do smid=1,num_smooth
                  ! Set ghost cell periodic in theta direction
                  target_Data2D(0,:)     = target_Data2D(Ngq-1,:)
                  target_Data2D(Ngq+1,:) = target_Data2D(2,:)
                  ! Smooth in theta direction
                  target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
                  ! Smooth in radial direction
                  target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
               enddo
               DG%iVecData(k,:,:,vid,1) = target_Data2D(1:Ngq,1:Na)
            enddo
         enddo


         ! =======================================
         ! Smooth Electron's hybrid density and temperature
         ! =======================================


         ! Convert electron temperautre to Energy_Density first.
         DG%Te_hybrid  = DG%Te_hybrid * DG%ne_tot_hybrid

         !Smooth Electron's hybrid density
         target_Data2D(1:Ngq,1:Na) = DG%ne_tot_hybrid
         do smid=1,num_smooth
            ! Set ghost cell periodic in theta direction
            target_Data2D(0,:)     = target_Data2D(Ngq-1,:)
            target_Data2D(Ngq+1,:) = target_Data2D(2,:)
            ! Smooth in theta direction
            target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
            ! Smooth in radial direction
            target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
         enddo
         DG%ne_tot_hybrid = target_Data2D(1:Ngq,1:Na)

         !Smooth Electron's hybrid Temperature
         target_Data2D(1:Ngq,1:Na) = DG%Te_hybrid
         do smid=1,num_smooth
            ! Set ghost cell periodic in theta direction
            target_Data2D(0,:)     = target_Data2D(Ngq-1,:)
            target_Data2D(Ngq+1,:) = target_Data2D(2,:)
            ! Smooth in theta direction
            target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
            ! Smooth in radial direction
            target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
         enddo
         DG%Te_hybrid = target_Data2D(1:Ngq,1:Na) ! This is smoothed Energy Density
         DG%Te_hybrid = DG%Te_hybrid / DG%ne_tot_hybrid ! Now this is the smoothed temperature


         ! =======================================
         ! Smooth Ion's hybrid density and temperature
         ! =======================================

         ! Convert Ion temperature to Energy_Density first.
         DG%Ti_hybrid  = DG%Ti_hybrid * DG%ni_hybrid
         !Smooth Ion's hybrid density
         target_Data2D(1:Ngq,1:Na) = DG%ni_hybrid
         do smid=1,num_smooth
            ! Set ghost cell periodic in theta direction
            target_Data2D(0,:)     = target_Data2D(Ngq-1,:)
            target_Data2D(Ngq+1,:) = target_Data2D(2,:)
            ! Smooth in theta direction
            target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
            ! Smooth in radial direction
            target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
         enddo
         DG%ni_hybrid = target_Data2D(1:Ngq,1:Na)

         !Smooth Ion's hybrid Temperature
         target_Data2D(1:Ngq,1:Na) = DG%Ti_hybrid
         do smid=1,num_smooth
            ! Set ghost cell periodic in theta direction
            target_Data2D(0,:)     = target_Data2D(Ngq-1,:)
            target_Data2D(Ngq+1,:) = target_Data2D(2,:)
            ! Smooth in theta direction
            target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
            ! Smooth in radial direction
            target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
         enddo
         DG%Ti_hybrid = target_Data2D(1:Ngq,1:Na) ! This is smoothed Energy Density
         DG%Ti_hybrid = DG%Ti_hybrid / DG%ni_hybrid ! Now this is the smoothed temperature




      end subroutine Smooth_DG_Scalar_Vector_Data


      subroutine Get_Conserving_Smooth_Coeffs_on_DG(in_VOLorAREA, SM_a_Coeff, SM_gq_Coeff)

         use precision
         use global_parameters
         use field_array
         implicit none

         integer :: Na,Ngq
         real(wp), intent(in) :: in_VOLorAREA(DG%dim_gq, DG%dim_a) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
         real(wp) :: VV(0:DG%dim_gq+1, 0:DG%dim_a+1) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
         real(wp), intent(out):: SM_gq_Coeff(-1:1,0:DG%dim_gq+1, 0:DG%dim_a+1), SM_a_Coeff(-1:1,0:DG%dim_gq+1, 0:DG%dim_a+1)

         Na  = DG%dim_a
         Ngq = DG%dim_gq

         VV = 0
         VV(1:Ngq,1:Na) = in_VOLorAREA

         ! Note that theta direction is periodic
         VV(0,:) = VV(Ngq-1,:)
         VV(Ngq+1,:) = VV(2,:)

         SM_gq_Coeff(-1,1:Ngq,:) = 0.5*VV(0:Ngq-1,:)/(VV(0:Ngq-1,:) + VV(1:Ngq,:)) ! Coeff_L = 0.5*VL/(VL+VC)
         SM_gq_Coeff( 0,1:Ngq,:) = 0.5*VV(1:Ngq,:)*(VV(0:Ngq-1,:)+2.0*VV(1:Ngq,:)+VV(2:Ngq+1,:))/((VV(0:Ngq-1,:)+VV(1:Ngq,:))*(VV(1:Ngq,:)+VV(2:Ngq+1,:))) ! 0.5*VC*(VL+2*VC+VR)/((VL+VC)*(VC+VR))
         SM_gq_Coeff( 1,1:Ngq,:) = 0.5*VV(2:Ngq+1,:)/(VV(2:Ngq+1,:) + VV(1:Ngq,:)) ! Coeff_R = 0.5*VR/(VR+VC)

         SM_a_Coeff(-1,:,1:Na) = 0.5*VV(:,0:Na-1)/(VV(:,0:Na-1) + VV(:,1:Na)) ! Coeff_L = 0.5*VL/(VL+VC)
         SM_a_Coeff( 0,:,1:Na) = 0.5*VV(:,1:Na)*(VV(:,0:Na-1)+2.0*VV(:,1:Na)+VV(:,2:Na+1))/((VV(:,0:Na-1)+VV(:,1:Na))*(VV(:,1:Na)+VV(:,2:Na+1))) ! 0.5*VC*(VL+2*VC+VR)/((VL+VC)*(VC+VR))
         SM_a_Coeff( 1,:,1:Na) = 0.5*VV(:,2:Na+1)/(VV(:,2:Na+1) + VV(:,1:Na)) ! Coeff_R = 0.5*VR/(VR+VC)

         ! To remove NaN values at the raidal ghost cells
         SM_gq_Coeff(:,:,0)   = 0
         SM_gq_Coeff(:,:,Na+1)= 0
         SM_a_Coeff(:,:,0)    = 0
         SM_a_Coeff(:,:,Na+1) = 0

      end subroutine Get_Conserving_Smooth_Coeffs_on_DG

      subroutine Smooth_DG_Data_2D(num_smooth, target_Data2D, Reference_Data2D)

         use precision
         use global_parameters
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
            SM_Data(0,:)     = SM_Data(Ngq-1,:)
            SM_Data(Ngq+1,:) = SM_Data(2,:)
            ! Smooth in theta direction
            SM_Data(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*SM_Data(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*SM_Data(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*SM_Data(2:Ngq+1,:)
            ! Smooth in radial direction
            SM_Data(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*SM_Data(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*SM_Data(:,1:Na) + SM_a_Coeff(1,:,1:Na)*SM_Data(:,2:Na+1)
         enddo
         target_Data2D = SM_Data(1:Ngq,1:Na)

      end subroutine Smooth_DG_Data_2D

      subroutine Do_FieldLineAverage_from_oriGrid_to_DG(oriG_TargetData, DG_FLA_ResultData)
         use mpi_defs
         use precision
         use global_parameters
         use particle_decomp
         use FieldLineFollowing
         use field_array
         implicit none

         real(wp), intent(in) :: oriG_TargetData(0:1,mgrid)
         real(wp), intent(out):: DG_FLA_ResultData(DG%dim_gq, DG%dim_a,0:1)
         real(wp) tmp_DG_Data(DG%dim_gq, DG%dim_a)

         real(wp) each_FLA_val(num_of_SampledFieldLines), each_FLA_denominator(num_of_SampledFieldLines), intpl_val ! FLA=Field Line Average
         real(wp) this_gq, this_s, this_gf, sdum, gqdum,inv_ds, inv_dgq(0:mpsi), inv_gf
         real(wp) ws0, ws1, wgq00, wgq01, wgq10, wgq11, wgf0, wgf1
         integer Fid, Pid ! (Fid=Field ID), (Pid=gPPI ID)
         integer inner_sid, outter_sid, inner_gqid, outter_gqid, inner_ij, outter_ij
         integer ierror

         real(wp) recvr(mgrid)
         integer istatus(MPI_STATUS_SIZE)

         real(wp) target_gq, target_zeta
         real(wp) wgq0,wgq1,wgf,ratio
         integer ori_sid, j,k, sid, gqid, qq, ss, length
         real(wp) ff, this_safety_q



         if (num_of_SampledFieldLines ==0) then
            write(*,*) "Error! Please construct gPPI first before using FieldLineAverage!!"
            STOP
         endif

         if(myrank_toroidal==0) then
            write(*,*) "Do_FieldLineAverage_by_using_gPPI"
         endif
         inv_ds = real(mpsi)/(s1-s0)
         inv_dgq = 1.0/deltat
         inv_gf = 1.0/deltaz

         ! Get FieldLineAverages of Target_data
         do Fid=1,num_of_SampledFieldLines
            each_FLA_val(Fid)=0
            each_FLA_denominator(Fid)=0
            do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

               this_s = gPPI(Pid,Fid)%s
               this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
               this_gf = gPPI(Pid,Fid)%target_normalized_gf

               sdum = (this_s-s0)*inv_ds
               inner_sid = int(sdum)
               ws1  = sdum - inner_sid
               ws0  = 1.0 - ws1
               gqdum = inv_dgq(inner_sid)*modulo(this_gq - qtinv(inner_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
               inner_gqid = int(gqdum)
               wgq01 = gqdum - inner_gqid
               wgq00 = 1.0 - wgq01
               inner_ij = igrid(inner_sid) + inner_gqid

               outter_sid = inner_sid+1
               gqdum = inv_dgq(outter_sid)*modulo(this_gq - qtinv(outter_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
               outter_gqid = int(gqdum)
               wgq11 = gqdum - outter_gqid
               wgq10 = 1.0 - wgq11
               outter_ij = igrid(outter_sid) + outter_gqid


               wgf1 = (this_gf-zetamin)*inv_gf
               wgf0 = 1.0 - wgf1

               intpl_val = wgf0*(ws0*(wgq00*oriG_TargetData(0,inner_ij) + wgq01*oriG_TargetData(0,inner_ij+1)) + ws1*(wgq10*oriG_TargetData(0,outter_ij)+wgq11*oriG_TargetData(0,outter_ij+1))) &
                  +wgf1*(ws0*(wgq00*oriG_TargetData(1,inner_ij) + wgq01*oriG_TargetData(1,inner_ij+1)) + ws1*(wgq10*oriG_TargetData(1,outter_ij)+wgq11*oriG_TargetData(1,outter_ij+1)))

               each_FLA_val(Fid) = each_FLA_val(Fid) + intpl_val*gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
               each_FLA_denominator(Fid) = each_FLA_denominator(Fid) + gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B

            enddo
         enddo

         ! Do Allreduce across toroidal_comm to get a sum of FLA
         call MPI_ALLREDUCE(each_FLA_val,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
         each_FLA_val = recvr(1:num_of_SampledFieldLines)
         call MPI_ALLREDUCE(each_FLA_denominator,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
         each_FLA_denominator = recvr(1:num_of_SampledFieldLines)

         each_FLA_val =  each_FLA_val/ each_FLA_denominator ! Take average



         ! Scatter Averaged data to Diagnostic Grid (DG) nodes
         DG_FLA_ResultData = 0.0
         do Fid=1,num_of_SampledFieldLines
            do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

               this_s = gPPI(Pid,Fid)%s
               this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
               this_gf = gPPI(Pid,Fid)%target_normalized_gf

               ori_sid = (this_s-s0)*inv_ds ! s index of original full grid
               ff  = (this_s-s_grid(ori_sid))*inv_ds
               this_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)

               sid = max(1,min(int((this_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
               ws1 = (this_s-DG%s_1D(sid))*DG%inv_ds
               ws0 = 1.0-ws1

               do k=0,1
                  if(k==0) target_zeta =zetamin
                  if(k==1) target_zeta =zetamax

                  wgf = (1.0-abs(this_gf-target_zeta)/deltaz)

                  target_gq = modulo(this_gq + (target_zeta-this_gf)/this_safety_q, 2.0*pi)
                  gqid = (target_gq)*DG%inv_dgq +1
                  wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
                  wgq0 = 1.0-wgq1

                  do j=1,4
                     select case(j)
                        case(1)
                           qq=gqid; ss=sid
                           ratio=ws0*wgq0*wgf
                        case(2)
                           qq=gqid+1; ss=sid
                           ratio=ws0*wgq1*wgf
                        case(3)
                           qq=gqid; ss=sid+1
                           ratio=ws1*wgq0*wgf
                        case(4)
                           qq=gqid+1; ss=sid+1
                           ratio=ws1*wgq1*wgf
                     end select
                     DG_FLA_ResultData(qq,ss,k)  = DG_FLA_ResultData(qq,ss,k)  + ratio*each_FLA_val(Fid)
                  enddo
               enddo
            enddo
         enddo



         ! Match poloidal end cell
         ! (gq=0 and gq=2pi is the same position physically
         ! although they have different gqid)
         do sid=1,DG%dim_a
            DG_FLA_ResultData(1,sid,:) = DG_FLA_ResultData(1,sid,:) + DG_FLA_ResultData(DG%dim_gq,sid,:)
            DG_FLA_ResultData(DG%dim_gq,sid,:) = DG_FLA_ResultData(1,sid,:)
         enddo

         ! Communicate DG_FLA_ResultData between nodes
         ! toroidal end cell
         length = DG%dim_gq*DG%dim_a
         call MPI_SENDRECV(DG_FLA_ResultData(:,:,0),length,mpi_Rsize,left_pe,1, tmp_DG_Data,length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         DG_FLA_ResultData(:,:,1) = DG_FLA_ResultData(:,:,1) + tmp_DG_Data
         call MPI_SENDRECV(DG_FLA_ResultData(:,:,1),length,mpi_Rsize,right_pe,1, DG_FLA_ResultData(:,:,0),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

         ! Last part is to divide FLA_ResultData by DG%FLA_Denominator to get correct line-averaged FLA_ResultData
         DG_FLA_ResultData = DG_FLA_ResultData/DG%FLA_Denominator

      end subroutine Do_FieldLineAverage_from_oriGrid_to_DG

      subroutine Do_FieldLineAverage_on_DG_by_using_gPPI(DG_TargetData, DG_FLA_ResultData)
         use mpi_defs
         use precision
         use global_parameters
         use particle_decomp
         use FieldLineFollowing
         use field_array
         implicit none

         real(wp), intent(in) :: DG_TargetData(DG%dim_gq, DG%dim_a,0:1)
         real(wp), intent(out):: DG_FLA_ResultData(DG%dim_gq, DG%dim_a,0:1)
         real(wp) tmp_DG_Data(DG%dim_gq, DG%dim_a)

         real(wp) each_FLA_val(num_of_SampledFieldLines), each_FLA_denominator(num_of_SampledFieldLines), intpl_val ! FLA=Field Line Average
         real(wp) this_gq, this_s, this_gf, sdum, gqdum,ori_inv_ds
         real(wp) ws0, ws1, wgq00, wgq01, wgq10, wgq11, wgf0, wgf1
         integer Fid, Pid ! (Fid=Field ID), (Pid=gPPI ID)
         integer inner_sid, outter_sid, inner_gqid, outter_gqid, inner_ij, outter_ij
         integer ierror

         real(wp) recvr(mgrid)
         integer istatus(MPI_STATUS_SIZE)

         real(wp) target_gq, target_zeta
         real(wp) wgq0,wgq1,wgf,ratio
         integer ori_sid, j,k, sid, gqid, qq, ss, length
         real(wp) ff, this_safety_q



         if (num_of_SampledFieldLines ==0) then
            write(*,*) "Error! Please construct gPPI first before using FieldLineAverage!!"
            STOP
         endif

         if(myrank_toroidal==0) then
            write(*,*) "Do_FieldLineAverage_by_using_gPPI"
         endif
         ori_inv_ds = real(mpsi)/(s1-s0)



         ! Scatter Averaged data to Diagnostic Grid (DG) nodes
         DG_FLA_ResultData = 0.0
         each_FLA_val = 0.0
         each_FLA_denominator = 0.0
         do Fid=1,num_of_SampledFieldLines
            do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

               this_s = gPPI(Pid,Fid)%s
               this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
               this_gf = gPPI(Pid,Fid)%target_normalized_gf

               ori_sid = (this_s-s0)*ori_inv_ds ! s index of original full grid
               ff  = (this_s-s_grid(ori_sid))*ori_inv_ds
               this_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)

               sid = max(1,min(int((this_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
               ws1 = (this_s-DG%s_1D(sid))*DG%inv_ds
               ws0 = 1.0-ws1

               intpl_val=0
               do k=0,1
                  if(k==0) target_zeta =zetamin
                  if(k==1) target_zeta =zetamax

                  wgf = (1.0-abs(this_gf-target_zeta)/deltaz)

                  target_gq = modulo(this_gq + (target_zeta-this_gf)/this_safety_q, 2.0*pi)
                  gqid = (target_gq)*DG%inv_dgq +1
                  wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
                  wgq0 = 1.0-wgq1

                  do j=1,4
                     select case(j)
                        case(1)
                           qq=gqid; ss=sid
                           ratio=ws0*wgq0*wgf
                        case(2)
                           qq=gqid+1; ss=sid
                           ratio=ws0*wgq1*wgf
                        case(3)
                           qq=gqid; ss=sid+1
                           ratio=ws1*wgq0*wgf
                        case(4)
                           qq=gqid+1; ss=sid+1
                           ratio=ws1*wgq1*wgf
                     end select
                     intpl_val = intpl_val + ratio*DG_TargetData(qq,ss,k)
                  enddo
               enddo
               each_FLA_val(Fid) = each_FLA_val(Fid) + intpl_val*gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
               each_FLA_denominator(Fid) = each_FLA_denominator(Fid) + gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
            enddo
         enddo

         ! Do Allreduce across toroidal_comm to get a sum of FLA
         call MPI_ALLREDUCE(each_FLA_val,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
         each_FLA_val = recvr(1:num_of_SampledFieldLines)
         call MPI_ALLREDUCE(each_FLA_denominator,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
         each_FLA_denominator = recvr(1:num_of_SampledFieldLines)

         each_FLA_val =  each_FLA_val/ each_FLA_denominator ! Take average



         ! Scatter Averaged data to Diagnostic Grid (DG) nodes
         DG_FLA_ResultData = 0.0
         do Fid=1,num_of_SampledFieldLines
            do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

               this_s = gPPI(Pid,Fid)%s
               this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
               this_gf = gPPI(Pid,Fid)%target_normalized_gf

               ori_sid = (this_s-s0)*ori_inv_ds ! s index of original full grid
               ff  = (this_s-s_grid(ori_sid))*ori_inv_ds
               this_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)

               sid = max(1,min(int((this_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
               ws1 = (this_s-DG%s_1D(sid))*DG%inv_ds
               ws0 = 1.0-ws1

               do k=0,1
                  if(k==0) target_zeta =zetamin
                  if(k==1) target_zeta =zetamax

                  wgf = (1.0-abs(this_gf-target_zeta)/deltaz)

                  target_gq = modulo(this_gq + (target_zeta-this_gf)/this_safety_q, 2.0*pi)
                  gqid = (target_gq)*DG%inv_dgq +1
                  wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
                  wgq0 = 1.0-wgq1

                  do j=1,4
                     select case(j)
                        case(1)
                           qq=gqid; ss=sid
                           ratio=ws0*wgq0*wgf
                        case(2)
                           qq=gqid+1; ss=sid
                           ratio=ws0*wgq1*wgf
                        case(3)
                           qq=gqid; ss=sid+1
                           ratio=ws1*wgq0*wgf
                        case(4)
                           qq=gqid+1; ss=sid+1
                           ratio=ws1*wgq1*wgf
                     end select
                     DG_FLA_ResultData(qq,ss,k)  = DG_FLA_ResultData(qq,ss,k)  + ratio*each_FLA_val(Fid)
                  enddo
               enddo
            enddo
         enddo


         ! Match poloidal end cell
         ! (gq=0 and gq=2pi is the same position physically
         ! although they have different gqid)
         do sid=1,DG%dim_a
            DG_FLA_ResultData(1,sid,:) = DG_FLA_ResultData(1,sid,:) + DG_FLA_ResultData(DG%dim_gq,sid,:)
            DG_FLA_ResultData(DG%dim_gq,sid,:) = DG_FLA_ResultData(1,sid,:)
         enddo

         ! Communicate DG_FLA_ResultData between nodes
         ! toroidal end cell
         length = DG%dim_gq*DG%dim_a
         call MPI_SENDRECV(DG_FLA_ResultData(:,:,0),length,mpi_Rsize,left_pe,1, tmp_DG_Data,length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         DG_FLA_ResultData(:,:,1) = DG_FLA_ResultData(:,:,1) + tmp_DG_Data
         call MPI_SENDRECV(DG_FLA_ResultData(:,:,1),length,mpi_Rsize,right_pe,1, DG_FLA_ResultData(:,:,0),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

         ! Last part is to divide FLA_ResultData by DG%FLA_Denominator to get correct line-averaged FLA_ResultData
         DG_FLA_ResultData = DG_FLA_ResultData/DG%FLA_Denominator

      end subroutine Do_FieldLineAverage_on_DG_by_using_gPPI

      subroutine Cal_FLA_FLabel_on_DG_by_using_gPPI_iniFLabel(DG_FLA_ResultData)
         use mpi_defs
         use precision
         use global_parameters
         use particle_decomp
         use FieldLineFollowing
         use field_array
         implicit none

         real(wp), intent(out):: DG_FLA_ResultData(DG%dim_gq, DG%dim_a,0:1)
         real(wp) tmp_DG_Data(DG%dim_gq, DG%dim_a)

         real(wp) each_FLA_val(num_of_SampledFieldLines), each_FLA_denominator(num_of_SampledFieldLines), intpl_val ! FLA=Field Line Average
         real(wp) this_gq, this_s, this_gf, sdum, gqdum,ori_inv_ds
         real(wp) ws0, ws1, wgq00, wgq01, wgq10, wgq11, wgf0, wgf1
         integer Fid, Pid ! (Fid=Field ID), (Pid=gPPI ID)
         integer inner_sid, outter_sid, inner_gqid, outter_gqid, inner_ij, outter_ij
         integer ierror

         real(wp) recvr(mgrid)
         integer istatus(MPI_STATUS_SIZE)

         real(wp) target_gq, target_zeta
         real(wp) wgq0,wgq1,wgf,ratio
         integer ori_sid, j,k, sid, gqid, qq, ss, length
         real(wp) ff, this_safety_q



         if (num_of_SampledFieldLines ==0) then
            write(*,*) "Error! Please construct gPPI first before using FieldLineAverage!!"
            STOP
         endif

         if(myrank_toroidal==0) then
            write(*,*) "Do_FieldLineAverage_by_using_gPPI"
         endif
         ori_inv_ds = real(mpsi)/(s1-s0)

         ! Scatter Averaged data to Diagnostic Grid (DG) nodes
         DG_FLA_ResultData = 0.0
         each_FLA_val = 0.0
         each_FLA_denominator = 0.0
         do Fid=1,num_of_SampledFieldLines
            do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)
               each_FLA_val(Fid) = each_FLA_val(Fid) + gPPI(Pid,Fid)%ini_FLabel*gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
               each_FLA_denominator(Fid) = each_FLA_denominator(Fid) + gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
            enddo
         enddo

         ! Do Allreduce across toroidal_comm to get a sum of FLA
         call MPI_ALLREDUCE(each_FLA_val,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
         each_FLA_val = recvr(1:num_of_SampledFieldLines)
         call MPI_ALLREDUCE(each_FLA_denominator,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
         each_FLA_denominator = recvr(1:num_of_SampledFieldLines)

         each_FLA_val =  each_FLA_val/ each_FLA_denominator ! Take average



         ! Scatter Averaged data to Diagnostic Grid (DG) nodes
         DG_FLA_ResultData = 0.0
         do Fid=1,num_of_SampledFieldLines
            do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

               this_s = gPPI(Pid,Fid)%s
               this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
               this_gf = gPPI(Pid,Fid)%target_normalized_gf

               ori_sid = (this_s-s0)*ori_inv_ds ! s index of original full grid
               ff  = (this_s-s_grid(ori_sid))*ori_inv_ds
               this_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)

               sid = max(1,min(int((this_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
               ws1 = (this_s-DG%s_1D(sid))*DG%inv_ds
               ws0 = 1.0-ws1

               do k=0,1
                  if(k==0) target_zeta =zetamin
                  if(k==1) target_zeta =zetamax

                  wgf = (1.0-abs(this_gf-target_zeta)/deltaz)

                  target_gq = modulo(this_gq + (target_zeta-this_gf)/this_safety_q, 2.0*pi)
                  gqid = (target_gq)*DG%inv_dgq +1
                  wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
                  wgq0 = 1.0-wgq1

                  do j=1,4
                     select case(j)
                        case(1)
                           qq=gqid; ss=sid
                           ratio=ws0*wgq0*wgf
                        case(2)
                           qq=gqid+1; ss=sid
                           ratio=ws0*wgq1*wgf
                        case(3)
                           qq=gqid; ss=sid+1
                           ratio=ws1*wgq0*wgf
                        case(4)
                           qq=gqid+1; ss=sid+1
                           ratio=ws1*wgq1*wgf
                     end select
                     DG_FLA_ResultData(qq,ss,k)  = DG_FLA_ResultData(qq,ss,k)  + ratio*each_FLA_val(Fid)
                  enddo
               enddo
            enddo
         enddo


         ! Match poloidal end cell
         ! (gq=0 and gq=2pi is the same position physically
         ! although they have different gqid)
         do sid=1,DG%dim_a
            DG_FLA_ResultData(1,sid,:) = DG_FLA_ResultData(1,sid,:) + DG_FLA_ResultData(DG%dim_gq,sid,:)
            DG_FLA_ResultData(DG%dim_gq,sid,:) = DG_FLA_ResultData(1,sid,:)
         enddo

         ! Communicate DG_FLA_ResultData between nodes
         ! toroidal end cell
         length = DG%dim_gq*DG%dim_a
         call MPI_SENDRECV(DG_FLA_ResultData(:,:,0),length,mpi_Rsize,left_pe,1, tmp_DG_Data,length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         DG_FLA_ResultData(:,:,1) = DG_FLA_ResultData(:,:,1) + tmp_DG_Data
         call MPI_SENDRECV(DG_FLA_ResultData(:,:,1),length,mpi_Rsize,right_pe,1, DG_FLA_ResultData(:,:,0),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

         ! Last part is to divide FLA_ResultData by DG%FLA_Denominator to get correct line-averaged FLA_ResultData
         DG_FLA_ResultData = DG_FLA_ResultData/DG%FLA_Denominator

      end subroutine Cal_FLA_FLabel_on_DG_by_using_gPPI_iniFLabel

      subroutine Cal_FLA_Denominator_on_DG_using_gPPI
         use mpi_defs
         use precision
         use global_parameters
         use particle_decomp
         use FieldLineFollowing
         use field_array

         implicit none

         real(wp) tmp_Data(DG%dim_gq,DG%dim_a)
         real(wp) this_gq, this_s, this_gf, ori_inv_ds, this_safety_q
         real(wp) target_gq, target_zeta
         real(wp) ff, ws0,ws1, wgq0,wgq1, wgf,ratio ! ratios along each dimension
         integer Fid, Pid ! (Fid=Field ID), (Pid=gPPI ID)
         integer :: ori_sid, j,k, sid, gqid, qq, ss, length
         integer istatus(MPI_STATUS_SIZE),ierror


         if (num_of_SampledFieldLines ==0) then
            write(*,*) "Error! Please construct gPPI first before using FieldLineAverage!!"
            STOP
         endif

         ori_inv_ds = real(mpsi)/(s1-s0)

         ! Scatter Averaged data to grid nodes
         DG%FLA_Denominator = 0.0
         do Fid=1,num_of_SampledFieldLines
            do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

               this_s = gPPI(Pid,Fid)%s
               this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
               this_gf = gPPI(Pid,Fid)%target_normalized_gf

               ori_sid = (this_s-s0)*ori_inv_ds ! s index of original full grid
               ff  = (this_s-s_grid(ori_sid))*ori_inv_ds
               this_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)


               sid = max(1,min(int((this_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
               ws1 = (this_s-DG%s_1D(sid))*DG%inv_ds
               ws0 = 1.0-ws1

               do k=0,1
                  if(k==0) target_zeta =zetamin
                  if(k==1) target_zeta =zetamax

                  wgf = (1.0-abs(this_gf-target_zeta)/deltaz)

                  target_gq = modulo(this_gq + (target_zeta-this_gf)/this_safety_q, 2.0*pi)
                  gqid = (target_gq)*DG%inv_dgq +1
                  wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
                  wgq0 = 1.0-wgq1

                  do j=1,4
                     select case(j)
                        case(1)
                           qq=gqid; ss=sid
                           ratio=ws0*wgq0*wgf
                        case(2)
                           qq=gqid+1; ss=sid
                           ratio=ws0*wgq1*wgf
                        case(3)
                           qq=gqid; ss=sid+1
                           ratio=ws1*wgq0*wgf
                        case(4)
                           qq=gqid+1; ss=sid+1
                           ratio=ws1*wgq1*wgf
                     end select
                     DG%FLA_Denominator(qq,ss,k) = DG%FLA_Denominator(qq,ss,k) + ratio
                  enddo
               enddo
            enddo
         enddo

         ! Match poloidal end cell
         ! (gq=0 and gq=2pi is the same position physically
         ! although they have different gqid)
         do sid=1,DG%dim_a
            DG%FLA_Denominator(1,sid,:) = DG%FLA_Denominator(1,sid,:) + DG%FLA_Denominator(DG%dim_gq,sid,:)
            DG%FLA_Denominator(DG%dim_gq,sid,:) = DG%FLA_Denominator(1,sid,:)
         enddo

         ! Communicate DG%FLA_FLabel between nodes
         ! toroidal end cell
         length = DG%dim_gq*DG%dim_a
         call MPI_SENDRECV(DG%FLA_Denominator(:,:,0),length,mpi_Rsize,left_pe,1, tmp_Data,length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         DG%FLA_Denominator(:,:,1) = DG%FLA_Denominator(:,:,1) + tmp_Data
         call MPI_SENDRECV(DG%FLA_Denominator(:,:,1),length,mpi_Rsize,right_pe,1, DG%FLA_Denominator(:,:,0),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

      end subroutine Cal_FLA_Denominator_on_DG_using_gPPI

end module Diagnostic_Grid
