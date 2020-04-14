!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!            GYROKINETIC TOKAMAK SIMULATION (GTS) CODE                       !
!                                                                            !
!   1. This code simulates micro-turbulence and transport                    !
!      in toroidal (tokamak) experiment;                                     !
!   2. It is based on a generalized gyrokinetic particle simulation model;   !
!   3. It is developed using GTC architecture                             !
!                                                                            !
!              Princeton Plasma Physics Laboratory                           !
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program GTS
   use precision
   !2015-0ct-24 already included in "solver_data_class": use mpi_defs
   use global_parameters
   use particle_array
   use field_array
   use diagnosis_array
   use solver_data_class
   use fem_data_class
   use particle_decomp
   use Ptl_Src_Loss

#ifdef _ACCEL
   use accel_lib
#endif _ACCEL
   implicit none

   ! to use parallel random number generator
#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

   SPRNG_POINTER stream
   character buffer1(MAX_PACKED_LENGTH)

   integer i,ierror,rng_type
   integer,parameter :: nc=18
   real(doubleprec) program_start_t0, program_start_t0wc, t0,t0wc,time(nc),timewc(nc),dt,dtwc
   real(doubleprec) tmp_time(nc), tmp_timewc(nc), loop_start_time,loop_start_timewc, mainloop_time, mainloop_timewc
   character(len=16) ic(nc)
   type(solver_data) :: solver1,solver_dphi_dt,masssolver
   type(fem_data) :: fem1,femt,femm
   character(len=40) cdum
   integer DT_int(8)
   character(len=10) DT_str(3)
   integer:: nrk !number of RK step (nrk=2 for RK2, nrk=4 for RK4)
   integer:: ngpus_per_node !number of gpus per node
   integer:: device_type,my_gpu_device_num

   interface
      subroutine shift(phase,phase0,mp,mpmax)
         real,pointer,dimension(:,:) :: phase,phase0
         integer,intent(INOUT) :: mp,mpmax
      end subroutine shift
   end interface

   call petsc_init()

   ! total # of PE and rank of PE
   call mpi_comm_size(mpi_comm_world,numberpe,ierror)
   call mpi_comm_rank(mpi_comm_world,mype,ierror)


   if(mype==0) then
      call Print_Git_Info()
   endif

#ifdef _ACCEL
   device_type=acc_get_device_type()
   ngpus_per_node=acc_get_num_devices(device_type)
   my_gpu_device_num = acc_get_device_num(device_type)
   write(stdout,'(3(a,I5))') "Before set_device, mype=",mype,"; ngpus_per_node",ngpus_per_node, "; my_gpu_device_num=",my_gpu_device_num
   call flush(stdout)
   call acc_set_device_num(mod(mype,ngpus_per_node),device_type)
   my_gpu_device_num = acc_get_device_num(device_type)
   write(stdout,'(3(a,I5))') "After set_device,  mype=",mype,"; ngpus_per_node",ngpus_per_node, "; my_gpu_device_num=",my_gpu_device_num
   call flush(stdout)
#endif _ACCEL

   time=0.0_wp
   t0=0.0_wp
   timewc=0.0_wp
   t0wc=0.0_wp

   call timer(t0,dt,t0wc,dtwc)
   program_start_t0=t0
   program_start_t0wc=t0wc

   ! Set name of each clock
   ic(1)='push_ion'
   ic(2)='shifti'
   ic(3)='charge_ion'
   ic(4)='poisson'
   ic(5)='smooth'
   ic(6)='field'
   ic(7)='init'
   ic(8)='collision_ion'
   ic(9)='collision_eon'
   ic(10)='push_eon'
   ic(11)='shifte'
   ic(12)='charge_eon'
   ic(13)='check_tracers'
   ic(14)='snapshot'
   ic(15)='diagnostics'
   ic(16)='Main_Loop'
   ic(17)='others'
   ic(nc)='TOTAL'

   if(irun==0) then      !  SPRNG initialize
      seed=985456376     !  random number seed advised in sprng 1.0 docum.
      rng_type = SPRNG_LFG
      stream=init_sprng(rng_type,seed,SPRNG_DEFAULT)
      istart=1
   end if

   ! input parameters, setup equilibrium, allocate memory
   call setup

   !call Cal_Weighted_Neighboring_Volume_on_each_GridNode
   !call Cal_Weighted_Neighboring_Volume_on_each_GridNode_by_MCI
   if(myrank_partd==nproc_partd-1) then
      call Do_FieldLineFollowing_for_Poincare_Plot
      ! call Construct_gPPI_for_FieldLineAverage
      call Construct_gPPI_for_FieldLineAverage_Parallelized
   endif
   ! solver
   ! which solver to use
   ! order = 0 original solver
   ! order = 1 linear finite element
   ! order = 2 2nd order finite element
   ! which bc_type to use
   ! bc_type=0 neumann bc at a0, neumann bc at a1
   ! bc_type=1 dirichlet bc at a0, neumann bc at a1
   ! bc_type=2 dirichlet bc at a0, dirichlet bc at a1
   ! bc_type=3 neumann bc at a0, dirichlet bc at a1
   if( order .le. 0 ) then
      solver1%order=0
      solver_dphi_dt%order=0
      if(mype == 0) write (0, *) "	integral equation order=", order
      ! mat type 1, bc type 1
      solver1%mat%bc_type = 2
      call new_gk_solver( solver1, '', 1, 2 )
      ! mat type 2, bc type 2
      solver_dphi_dt%mat%bc_type = 2
      if(ifelectron==2.and.ifdphi==1) call new_gk_solver( solver_dphi_dt, 'dphidt_', 2, 2 )
   else          !  for order > 0
      if(mype == 0) write (0, *) "	finite element order=", order

      solver1%order = order
      solver1%mat%mat_type = 1
      solver1%mat%bc_type = 0
      fem1%order = order
      call new_fem_data( fem1 )
#ifdef CJTEST
      call new_fem_solver( fem1, solver1, '', solver1%mat%mat_type, 2 )
      call testsolve( solver1, masssolver )
      call delete_fem_data( fem1 )
      call destroy_solver( solver1 )
#else
      call new_fem_solver( fem1, solver1, '', solver1%mat%mat_type, solver1%mat%bc_type )
#endif

      masssolver%order = 0
      masssolver%mat%mat_type = -1
      masssolver%mat%bc_type = -1
      femm%order = 1
      call new_fem_data( femm )
      call new_fem_masssolver( femm, masssolver )
#ifdef CJTEST
      call testmasssolve( masssolver )
      call delete_fem_data( femm )
      call destroy_solver( masssolver )
#endif

      if(ifelectron==2.and.ifdphi==1) then
         solver_dphi_dt%order = order
         solver_dphi_dt%mat%mat_type = 2
         solver_dphi_dt%mat%bc_type = 2
         femt%order = order
         call new_fem_data( femt )
         call new_fem_solver( femt,solver_dphi_dt, 'dphidt_', solver_dphi_dt%mat%mat_type, solver_dphi_dt%mat%bc_type )
#ifdef CJTEST
         call testsolve( solver_dphi_dt, masssolver )
         call delete_fem_data( femt )
         call destroy_solver( solver_dphi_dt )
#endif
      endif

#ifdef CJTEST
      call PetscFinalize( ierror )
      stop
#endif
   endif

   ! initialize particle position and velocity
   if(Flag_Load_and_WeightDef==0) then
      call load
   else
      call load_uniform_MGY()
   endif
   if(mype==0) then
      call check_eq
      !     call view_distribution(0)
   end if

   if(mype==0) then
      call system("mkdir -p " // 'snapshots_hdf5/AllPlanes')
      call system("mkdir -p " // 'snapshots_hdf5')
   endif


   istep=0
   !  if(mype==0) call view_orbit

   ! initial field
   irk=2
   idiag=1
   call charge_ion

   if(ifelectron/=0) call charge_eon

   if(Flag_Poisson==1) call poisson( solver1, masssolver  )
   !  if(Flag_Poisson==1) call poisson( solver_dphi_dt )

   call smooth(0)

   call smooth(1)

   call field(masssolver)

   if(ifelectron==2.and.ifdphi==1) then
      call current_ion
      call current_eon
      call solve_dphidt(solver_dphi_dt,masssolver)    ! --> d(delta-phi)/dt
   end if

   !  call check_temperature

   call open_files

   ! ********************
   call Init_Diagnostic_Grid
   if(irun==0) then
      istep=0
      ! call get_Te_fluctuation_DirectlyFromParticles
      ! call get_Ti_fluctuation_DirectlyFromParticles
      call SetZeros_Cumulative_Variables
      call Cumulate_Some_Variables_for_Measure
      call snapshot_hdf5
      call snapshot_Magnetic_Information
      call check_tracer_e
      call check_tracer_i
      call snapshot     ! snapshot at t=0
   end if

   if(irun==-100) then
      mstep=(istart-1)+mstep   ! for continous debug runs with irun=-100
   end if

   !  phi_tmp=0.0
   ! ********************

   call MPI_BARRIER(MPI_COMM_WORLD,ierror)

   call timer(t0,dt,t0wc,dtwc)
   time(7)=time(7)+dt
   timewc(7)=timewc(7)+dtwc
   loop_start_time=t0
   loop_start_timewc=t0wc

   call system_mem_usage

   call date_and_time(DT_str(1),DT_str(2),DT_str(3),DT_int)
   if(mype==0) then
      ! write(cdum,'("elapsed_Time.dat_",A,".",A)')trim(DT_str(1)(3:8)),trim(DT_str(2)(1:4))
      ! open(unit=190515,file=cdum,action='write')
      open(unit=190515,file="elapsed_Time.dat",action="write",position="append")
   endif

   if(Flag_RK4==1) then
      nrk=4 !RK4
   else
      nrk=2 !RK2
   endif

   call SetZeros_Cumulative_Variables
   !           MAIN TIME LOOP ----- MAIN TIME LOOP ----- MAIN TIME LOOP ----- MAIN TIME LOOP
   do istep=istart,mstep
      if(mype==0) then
         write(stdout,'(a,I8)')'istep=',istep
         call flush(stdout)
      endif
      do irk=1,nrk

         ! idiag=0: do time history diagnosis
         idiag=mod(irk-1,nrk)+mod(istep,ndiag)

         ! calculate test like-particle collisions
         call collision  ! collision is on when icollision>0 and calculated every icollision-step at irk=2

         call timer(t0,dt,t0wc,dtwc)
         time(8)=time(8)+dt
         timewc(8)=timewc(8)+dtwc

         ! push ion
         if(Flag_RK4==1) then
            call pushion_RK4(irk)
         else
            call pushion
         endif
         call timer(t0,dt,t0wc,dtwc)
         time(1)=time(1)+dt
         timewc(1)=timewc(1)+dtwc


         if(istep==mstep) then
            if(mype==0) then
               call view_distribution(irk-1)
            endif
         endif
         call timer(t0,dt,t0wc,dtwc)
         time(17)=time(17)+dt
         timewc(17)=timewc(17)+dtwc
         ! redistribute ion across PEs
         call shifti
         call timer(t0,dt,t0wc,dtwc)
         time(2)=time(2)+dt
         timewc(2)=timewc(2)+dtwc

         ! Calculate ion weight damping
         if(irk==nrk.and.wimcompressibility==2) call weight_damping_ion

         ! Calculate ion perturbed density "delta-ni" and ion locations in the 3D field
         ! mesh, which are ued in CURRENT_EON for calculating gamma_i as well as
         ! in PUSHI for calculating electric force on ions
         call charge_ion
         call timer(t0,dt,t0wc,dtwc)
         time(3)=time(3)+dt
         timewc(3)=timewc(3)+dtwc


         if(ifelectron/=0) then    ! include kinetic electrons

            ! NOTE: nelectron must take 1 or even number due to the use of 2nd order Runge-Kutta
            if(nelectron==1) then  ! electrons have the same step size as ions
               idiag_e=mod(irk-1,nrk)+mod(istep,ndiag)
               ! if idiag_e=0, then istep=n*ndiag, irk=1, calculate fluxes
               ! via calling get_fluexs_e in pushe

               !  March, 2012: PUSHEON uses more a different way (more precise in form -- see note) to calculate RHS of delta-fe equation,
               !  but give almost identical results for particle weight and ererything else according to a short test (20 steps)

               ! calculate electron collisions
               call collision_eon  ! collision is on when icollision>0 and calculated every icollision-step at irk=2
               call timer(t0,dt,t0wc,dtwc)
               time(9)=time(9)+dt
               timewc(9)=timewc(9)+dtwc

               !     CALL PUSHE(irk)
               if(Flag_RK4==1) then
                  call pusheon_RK4(solver1, masssolver, irk)
               else
                  call pusheon(solver1, masssolver, irk)
               endif
               call timer(t0,dt,t0wc,dtwc)
               time(10)=time(10)+dt
               timewc(10)=timewc(10)+dtwc

               call shifte
               call timer(t0,dt,t0wc,dtwc)
               time(11)=time(11)+dt
               timewc(11)=timewc(11)+dtwc

            else
               ! Note that Following Subelectron routines are only valid for RK2 scheme
               ! 1 ion step = nelectron(>1) electron steps
               ! originally, 1 ion step = 2*nelectron electron steps
               if(irk==1) then
                  me00=me
                  ntracer00=ntracer_e
                  zeon00=0.0_wp
                  !$omp parallel do private(i)
                  do i=1,me
                     zeon00(:,i)=zeon(:,i)
                  end do
               end if

               if(irk==2) then
                  me=me00
                  ntracer_e=ntracer00
                  zeon=0.0_wp
                  !$omp parallel do private(i)
                  do i=1,me
                     zeon(:,i)=zeon00(:,i)
                  end do
               end if

               do i=1,(nelectron/2)*irk
                  idiag_e=mod(istep,ndiag)+mod(irk-1,2)+mod(i-1,(nelectron/2)*irk)
                  ! if idiag_e=0, then istep=n*ndiag, irk=1, and i=1, calculate fluxes via calling
                  ! get_fluexs_e in pushe
                  ! push electron -- the 1st step of 2nd order Runge-Kutta
                  !    CALL PUSHE(1)
                  call pusheon(solver1, masssolver, 1)
                  call timer(t0,dt,t0wc,dtwc)
                  time(10)=time(10)+dt
                  timewc(10)=timewc(10)+dtwc

                  call shifte
                  call timer(t0,dt,t0wc,dtwc)
                  time(11)=time(11)+dt
                  timewc(11)=timewc(11)+dtwc

                  ! push electron -- the 2nd step of 2nd order Runge-Kutta
                  !    CALL PUSHE(2)
                  call pusheon(solver1, masssolver, 2)
                  call timer(t0,dt,t0wc,dtwc)
                  time(10)=time(10)+dt
                  timewc(10)=timewc(10)+dtwc

                  call shifte
                  call timer(t0,dt,t0wc,dtwc)
                  time(11)=time(11)+dt
                  timewc(11)=timewc(11)+dtwc

               end do
            end if    ! end if(nelectron==1)

            if(irk==nrk.and.wimcompressibility==2) call weight_damping_eon

            ! Calculate perturbed density and current of nonadiabatic electrons "delta-ne^h"
            ! and electron locations in the 3D field mesh,  which are needed in CURRENT_EON
            ! for calculating gamma_e
            call charge_eon
            call timer(t0,dt,t0wc,dtwc)
            time(12)=time(12)+dt
            timewc(12)=timewc(12)+dtwc

         endif    ! end if(ifelectron/=0)

         ! solve GK Poisson equation using adiabatic electron
         if(Flag_Poisson==1) call poisson( solver1, masssolver )
         !if(Flag_Poisson==1) call poisson( solver_dphi_dt )

         call timer(t0,dt,t0wc,dtwc)
         time(4)=time(4)+dt
         timewc(4)=timewc(4)+dtwc

         ! smoothing,field diagnosis and finite difference for e_field
         call smooth(0)

         call timer(t0,dt,t0wc,dtwc)
         time(5)=time(5)+dt
         timewc(5)=timewc(5)+dtwc

         idiag=mod(irk,nrk)+mod(istep+1,ndiag)

         call smooth(1)  ! 1. solve zonal flow for nonlinear run;
         ! 2. smooth the complete phi;
         ! 3. fft and filter

         call timer(t0,dt,t0wc,dtwc)
         time(5)=time(5)+dt
         timewc(5)=timewc(5)+dtwc

         ! field
         call field(masssolver)

         call timer(t0,dt,t0wc,dtwc)
         time(6)=time(6)+dt
         timewc(6)=timewc(6)+dtwc

         if(ifelectron==2.and.ifdphi==1) then
            call current_ion
            call current_eon
            call solve_dphidt(solver_dphi_dt,masssolver)    ! --> d(delta-phi)/dt
         end if

      enddo  ! enddo irk=1,nrk

      call Print_Reference_ntracer_Information

      ! profile snapshots, write particle information to restart file
      !  phi_tmp(:)=phi_tmp(:)+phi(1,:) ! a time average
      call Cumulate_Some_Variables_for_Measure

      if(Flag_Ptl_BC_type==1) then
         call Inject_electrons_from_Boundaries(tstep*utime)
         call Inject_ions_from_Boundaries(tstep*utime)
      endif
      ! call load_electrons_In_or_Out_Flux(100,tstep*utime)
      ! call load_electrons_In_or_Out_Flux(10,tstep)

      if(mod(istep,msnap) .eq. 0) then
         !$acc update self(elabels(:me))
         !$acc update self(UnNorm_theta_e(:me),UnNorm_zeta_e(:me))
         !$acc update self(we_vpara(:me), we_vperp(:me), we_accel(:me), we_fdamp(:me))
         if(acc_newer_D_zeon) then
            !$acc update self(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_D_zeon=.false.
            acc_newer_H_zeon=.false.
         endif
         if(acc_newer_D_zion) then
            !$acc update self(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
            acc_newer_D_zion=.false.
            acc_newer_H_zion=.false.
         endif
         !  phi_tmp=phi_tmp/real(msnap,wp)
         call snapshot_hdf5
         call snapshot ! snapshot must be called after snapshot_hdf5, because snapshot changed densitye values inside
         call SetZeros_Cumulative_Variables
         call system_mem_usage
         !  phi_tmp=0.0_wp
      end if
      call timer(t0,dt,t0wc,dtwc)
      time(14)=time(14)+dt
      timewc(14)=timewc(14)+dtwc

      if(mod(istep,nstep_for_tracer_record)==0) then
         !$acc update self(elabels(:me))
         !$acc update self(UnNorm_theta_e(:me),UnNorm_zeta_e(:me))
         !$acc update self(we_vpara(:me), we_vperp(:me), we_accel(:me), we_fdamp(:me))
         if(acc_newer_D_zeon) then
            !$acc update self(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_D_zeon=.false.
            acc_newer_H_zeon=.false.
         endif
         !$acc update self(ilabels(:mi))
         !$acc update self(UnNorm_theta_i(:mi),UnNorm_zeta_i(:mi))
         !$acc update self(wi_vpara(:mi), wi_vperp(:mi), wi_accel(:mi), wi_fdamp(:mi))
         if(acc_newer_D_zion) then
            !$acc update self(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
            acc_newer_D_zion=.false.
            acc_newer_H_zion=.false.
         endif

         call check_tracer_e
         call check_tracer_i
         if(Flag_Load_and_WeightDef==0) then
            call check_abnormal_weight_electrons
            call check_abnormal_weight_ions
         endif
         if(Flag_w_Threshold==2) call Abandon_Abnormal_Particles
      endif
      call timer(t0,dt,t0wc,dtwc)
      time(13)=time(13)+dt
      timewc(13)=timewc(13)+dtwc



      call timer(t0,dt,t0wc,dtwc)
      time(nc)=t0-program_start_t0
      timewc(nc)=t0wc-program_start_t0wc
      time(16)=t0-loop_start_time
      timewc(16)=t0wc-loop_start_timewc
      ! Print Elapsed Times
      if(mype==0 .and. mod(istep,1)==0) then

         tmp_time    = time
         tmp_timewc  = timewc

         tmp_time(1)  = tmp_time(1)  - sum(time_in_Pushi) ! pushion
         tmp_time(4)  = tmp_time(4)  + time_in_Pushe(3)   ! Poisson
         tmp_time(5)  = tmp_time(5)  + time_in_Pushe(4)   ! smooth
         tmp_time(10) = tmp_time(10) - sum(time_in_Pushe) ! pusheon
         tmp_time(11) = tmp_time(11) + time_in_Pushe(1)   ! shifte
         tmp_time(12) = tmp_time(12) + time_in_Pushe(2)   ! charge_eon
         tmp_time(15) = tmp_time(15) + time_in_Pushe(5) +time_in_Pushi(1) ! diagnostics in push routines

         tmp_timewc(1)  = tmp_timewc(1)  - sum(timewc_in_Pushi) ! pushion
         tmp_timewc(4)  = tmp_timewc(4)  + timewc_in_Pushe(3)   ! Poisson
         tmp_timewc(5)  = tmp_timewc(5)  + timewc_in_Pushe(4)   ! smooth
         tmp_timewc(10) = tmp_timewc(10) - sum(timewc_in_Pushe) ! pusheon
         tmp_timewc(11) = tmp_timewc(11) + timewc_in_Pushe(1)   ! shifte
         tmp_timewc(12) = tmp_timewc(12) + timewc_in_Pushe(2)   ! charge_eon
         tmp_timewc(15) = tmp_timewc(15) + timewc_in_Pushe(5) +timewc_in_Pushi(1) ! diagnostics in push routines

         write(190515,*) 'istep=', istep
         write(190515,'(5a16)')' ROUTINES       ',' CPU TIME(sec)',' WALL CLOCK TIME', '#/Total (%)', '#/MainLoop (%)'
         write(190515,'(5a16)')'==========      ',' =============',' ===============',' ===============',' ==============='
         do i=1,nc
            write(190515,'(a16,2f16.6, f16.4)',advance="no")ic(i),tmp_time(i),tmp_timewc(i), tmp_timewc(i)/tmp_timewc(nc)*100
            if(i/=7 .and. i<nc-1) then
               write(190515,'(f16.4)') tmp_timewc(i)/tmp_timewc(16)*100
            else
               write(190515,*) ' '
            endif
         enddo
         write(190515,*) ' '
         write(190515,*) ' '
         write(190515,*) ' '
         call flush(190515)
      endif
   enddo      !    MAIN TIME LOOP ----- MAIN TIME LOOP ----- MAIN TIME LOOP ----- MAIN TIME LOOP

   if(mype==0) close(190515)

   if(mype==0) write(stdout,*) 'end of main time loop'

   if(mype==0) then
      !     call view_distribution(1)
   end if

   call timer(t0,dt,t0wc,dtwc)
   time(nc)=t0-program_start_t0
   timewc(nc)=t0wc-program_start_t0wc
   time(16)=t0-loop_start_time
   timewc(16)=t0wc-loop_start_timewc

   if(mype==0)then

      tmp_time    = time
      tmp_timewc  = timewc

      tmp_time(1)  = tmp_time(1)  - sum(time_in_Pushi) ! pushion
      tmp_time(4)  = tmp_time(4)  + time_in_Pushe(3)   ! Poisson
      tmp_time(5)  = tmp_time(5)  + time_in_Pushe(4)   ! smooth
      tmp_time(10) = tmp_time(10) - sum(time_in_Pushe) ! pusheon
      tmp_time(11) = tmp_time(11) + time_in_Pushe(1)   ! shifte
      tmp_time(12) = tmp_time(12) + time_in_Pushe(2)   ! charge_eon
      tmp_time(15) = tmp_time(15) + time_in_Pushe(5) +time_in_Pushi(1) ! diagnostics in push routines

      tmp_timewc(1)  = tmp_timewc(1)  - sum(timewc_in_Pushi) ! pushion
      tmp_timewc(4)  = tmp_timewc(4)  + timewc_in_Pushe(3)   ! Poisson
      tmp_timewc(5)  = tmp_timewc(5)  + timewc_in_Pushe(4)   ! smooth
      tmp_timewc(10) = tmp_timewc(10) - sum(timewc_in_Pushe) ! pusheon
      tmp_timewc(11) = tmp_timewc(11) + timewc_in_Pushe(1)   ! shifte
      tmp_timewc(12) = tmp_timewc(12) + timewc_in_Pushe(2)   ! charge_eon
      tmp_timewc(15) = tmp_timewc(15) + timewc_in_Pushe(5) +timewc_in_Pushi(1) ! diagnostics in push routines

      write(stdout,'(5a16)')' ROUTINES       ',' CPU TIME(sec)',' WALL CLOCK TIME', 'REL_WC_TIME (%)', '#/MainLoop (%)'
      write(stdout,'(5a16)')'==========      ',' =============',' ===============',' ===============',' ==============='
      do i=1,nc
         write(stdout,'(a16,2f16.6, f16.4)',advance="no")ic(i),tmp_time(i),tmp_timewc(i), tmp_timewc(i)/tmp_timewc(nc)*100
         if(i/=7 .and. i<nc-1) then
            write(stdout,'(f16.4)') tmp_timewc(i)/tmp_timewc(16)*100
         else
            write(stdout,*) ' '
         endif
      enddo
   endif

   if( order .le. 0 ) then
      if(ifelectron==2.and.ifdphi==1) call destroy_solver(solver_dphi_dt)
      call destroy_solver(solver1)
   else
      call delete_fem_data( fem1 )
      call destroy_solver( solver1 )
      call delete_fem_data( femm )
      call destroy_solver( masssolver )
      if(ifelectron==2.and.ifdphi==1) then
         call delete_fem_data( femt )
         call destroy_solver( solver_dphi_dt )
      endif
   endif

   !free all allocated arrays
   !  call GTS_finalize(ierror)

   ! MPI finalize
   !  call mpi_finalize(ierror)
   call PetscFinalize(ierror)

end program GTS

!=========================================
subroutine timer(t0,dt,t0wc,dtwc)
   !=========================================
   use precision
   use mpi_defs
   implicit none
   real(doubleprec) t0,dt,t0wc,dtwc
   real(doubleprec) t1,t1wc

   ! Get cpu usage time since the beginning of the run and subtract value
   ! from the previous call
   call cpu_time(t1)
   dt=t1-t0
   t0=t1

   ! Get wall clock time and subtract value from the previous call
   t1wc=MPI_WTIME()
   dtwc=t1wc-t0wc
   t0wc=t1wc

end subroutine timer


subroutine system_mem_usage

   use mpi_defs
   use global_parameters
   implicit none

   integer :: val_VmRSS, val_VmSize, val_VmPeak

   character(len=200):: filename=' '
   character(len=80) :: line
   character(len=8)  :: pid_char=' '
   integer :: pid,getpid,rid,ierror
   logical :: ifxst

   integer :: int_arr(8), gathered_int_arr(8,numberpe)

   ! return negative number if not found
   val_VmRSS=-1
   val_VmSize=-1
   val_VmPeak=-1

   !--- get process ID
   pid=getpid()
   write(pid_char,'(I8)') pid
   filename='/proc/'//trim(adjustl(pid_char))//'/status'


   !--- read system file
   inquire (file=filename,exist=ifxst)
   if (.not.ifxst) then
      write (*,*) 'system file does not exist'
      return
   endif

   open(unit=1234, file=filename, action='read')
   do
      read (1234,'(a)') line
      if (line(1:7).eq.'VmSize:') then
         read (line(8:),*) val_VmSize
      elseif (line(1:7).eq.'VmPeak:') then
         read (line(8:),*) val_VmPeak
      elseif (line(1:6).eq.'VmRSS:') then
         read (line(7:),*) val_VmRSS
         exit
      endif
   enddo
   close(1234)

   !--- Print out mem. usage information
   if(mype==0) then
      write(stdout,'(3(a,I8),3(a,I12,a))')   'istep=', istep, ', mype=',mype,', pid=',pid, '; VmRSS=', val_VmRSS, ' kB;', ' VmSize=', val_VmSize, ' kB; ',' VmPeak=', val_VmPeak, ' kB; '
      call flush(stdout)
   endif

   int_arr(1:8)=[istep, mype, pid, mi, me, val_VmRSS, val_VmSize, val_VmPeak]

   call MPI_GATHER(int_arr,8,MPI_INTEGER, gathered_int_arr, 8, MPI_INTEGER, 0,MPI_COMM_WORLD,ierror)

   if(mype==0) then
      open(unit=190923,file='gts_mem_usage_every_rank.dat',position='append',action='write')
      ! write(190923,*) istep, mype, pid, mi, me, val_VmRSS, val_VmSize, val_VmPeak
      do rid=1,numberpe
         write(190923, "(8(I12,1x))") gathered_int_arr(:,rid)
      enddo
      call flush(190923)
      close(190923)

      open(unit=190924,file='gts_mem_usage_total_sum.dat',position='append',action='write')
      ! write(190924,*) istep, mi, me, val_VmRSS, val_VmSize, val_VmPeak
      write(190924,"(6(I12,1x))") gathered_int_arr(1,1), sum(gathered_int_arr(4:8,:),dim=2)
      call flush(190924)
      close(190924)
   endif

end subroutine
