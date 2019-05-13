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
  implicit none
  
! to use parallel random number generator
#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

 SPRNG_POINTER stream
 character buffer1(MAX_PACKED_LENGTH)

  integer i,ierror,rng_type
  integer,parameter :: nc=13
  real(doubleprec) t0,t0wc,time(nc),timewc(nc),dt,dtwc,loop_time
  character(len=16) ic(nc)
  type(solver_data) :: solver1,solver_dphi_dt,masssolver
  type(fem_data) :: fem1,femt,femm

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

  time=0.0_wp
  t0=0.0_wp
  timewc=0.0_wp
  t0wc=0.0_wp

  call timer(t0,dt,t0wc,dtwc)
  time(nc)=t0
  timewc(nc)=t0wc

  if(irun==0) then      !  SPRNG initialize
     seed=985456376     !  random number seed advised in sprng 1.0 docum.
     rng_type = SPRNG_LFG
     stream=init_sprng(rng_type,seed,SPRNG_DEFAULT)
     istart=1
  end if

  ! input parameters, setup equilibrium, allocate memory 
  call setup
  call snapshot_Magnetic_Information

  !call Cal_Weighted_Neighboring_Volume_on_each_GridNode
  !call Cal_Weighted_Neighboring_Volume_on_each_GridNode_by_MCI
  if(myrank_partd==nproc_partd-1) then    
    call Do_FieldLineFollowing_for_Poincare_Plot
    call Construct_gPPI_for_FieldLineAverage
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
  call load
  
  if(mype==0) then 
     call check_eq
!     call view_distribution(0)
  end if

  istep=0
!  if(mype==0) call view_orbit

! initial field
  irk=2
  idiag=1
  call charge_ion

  if(ifelectron/=0) call charge_eon

  call poisson( solver1, masssolver  )
!  call poisson( solver_dphi_dt )

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
if(irun==0) then 
   istep=0
   call SetZeros_Cumulative_Variables
   call Cumulate_Some_Variables_for_Measure   
   call snapshot     ! snapshot at t=0
   call snapshot_hdf5
   call check_tracer_e
   call check_tracer_i 
end if 

if(irun==-100) then
   mstep=(istart-1)+400   ! for continous debug runs with irun=-100
end if

!  phi_tmp=0.0
! ********************

  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  call timer(t0,dt,t0wc,dtwc)
  time(7)=time(7)+dt
  timewc(7)=timewc(7)+dtwc
  loop_time=t0wc

!           MAIN TIME LOOP ----- MAIN TIME LOOP ----- MAIN TIME LOOP ----- MAIN TIME LOOP 
  if (RK4_Flag<=0) then
    do istep=istart,mstep
      do irk=1,2

    ! idiag=0: do time history diagnosis
          idiag=mod(irk+1,2)+mod(istep,ndiag)

    ! calculate test like-particle collisions
          call collision  ! collision is on when icollision>0 and calculated every icollision-step at irk=2

          call timer(t0,dt,t0wc,dtwc)
          time(8)=time(8)+dt
          timewc(8)=timewc(8)+dtwc

    ! push ion
        !  call pushi
          call pushion

          call timer(t0,dt,t0wc,dtwc)
          time(1)=time(1)+dt
          timewc(1)=timewc(1)+dtwc
            
    !        if(idiag==0)then
    !           CALL DIAGNOSIS
    !           CALL OUTPUT
    !           call timer(t0,dt,t0wc,dtwc)
    !        endif

          if(istep==mstep) then
            if(mype==0) then
              call view_distribution(irk-1)
            endif
          endif
          
    ! redistribute ion across PEs
          call shifti

          call timer(t0,dt,t0wc,dtwc)
          time(2)=time(2)+dt
          timewc(2)=timewc(2)+dtwc

    ! Calculate ion weight damping 
          if(irk==2.and.wimcompressibility==2) call weight_damping_ion

          if(irk==2) call check_tracer_i
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
                idiag_e=mod(irk+1,2)+mod(istep,ndiag)
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
                call pusheon(solver1, masssolver, irk)
                call timer(t0,dt,t0wc,dtwc)
                time(10)=time(10)+dt
                timewc(10)=timewc(10)+dtwc

                call shifte
                call timer(t0,dt,t0wc,dtwc)
                time(11)=time(11)+dt
                timewc(11)=timewc(11)+dtwc

                if(irk==2) call check_tracer_e
            else                   ! 1 ion step = nelectron(>1) electron steps
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
                  idiag_e=mod(istep,ndiag)+mod(irk+1,2)+mod(i-1,(nelectron/2)*irk)
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

                  if(irk==2) call check_tracer_e
                end do
            end if    ! end if(nelectron==1)

            if(irk==2.and.wimcompressibility==2) call weight_damping_eon

    ! Calculate perturbed density and current of nonadiabatic electrons "delta-ne^h"
    ! and electron locations in the 3D field mesh,  which are needed in CURRENT_EON
    ! for calculating gamma_e
            call charge_eon
            call timer(t0,dt,t0wc,dtwc)
            time(12)=time(12)+dt
            timewc(12)=timewc(12)+dtwc

          endif    ! end if(ifelectron/=0)

    ! solve GK Poisson equation using adiabatic electron
          call poisson( solver1, masssolver )
    !        call poisson( solver_dphi_dt )

          call timer(t0,dt,t0wc,dtwc)
          time(4)=time(4)+dt
          timewc(4)=timewc(4)+dtwc

    ! smoothing,field diagnosis and finite difference for e_field
          call smooth(0)

          call timer(t0,dt,t0wc,dtwc)
          time(5)=time(5)+dt
          timewc(5)=timewc(5)+dtwc

          idiag=mod(irk,2)+mod(istep+1,ndiag)

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

      enddo  ! enddo irk=1,2

    ! profile snapshots, write particle information to restart file
          !  phi_tmp(:)=phi_tmp(:)+phi(1,:) ! a time average
          call Cumulate_Some_Variables_for_Measure

          if(mod(istep,msnap) .eq. 0) then 
            !  phi_tmp=phi_tmp/real(msnap,wp)
            call snapshot
            call snapshot_hdf5
            call SetZeros_Cumulative_Variables
            !  phi_tmp=0.0_wp
            call check_abnormal_weight_electrons
            call check_abnormal_weight_ions
          end if
          

    !        if(movie_3D_h5>0) then                   ! write potential 3D-data
    !           if(mod(istep,movie_3D_h5)==0) then    ! every movie_3D_h5 step
    !              call output3d
    !           end if
    !        end if
                
    ! if(mype==0) CALL view_orbit

      !if(mype==0) print*,' MAIN_TIME_LOOP',  istep
      if(ntracer>0)write(0,*)'istep=',istep,'  tracer=',zion(1:6,ntracer)
      if(ntracer_e>0)write(0,*)'istep=',istep,'  tracer_e=',zeon(1:6,ntracer_e)

    !  call check_temperature

      call timer(t0,dt,t0wc,dtwc) 

    enddo      !    MAIN TIME LOOP ----- MAIN TIME LOOP ----- MAIN TIME LOOP ----- MAIN TIME LOOP
  
  
  
  else if (RK4_Flag==1 .and. nelectron==1) then! When RK4_Flag>=1
    do istep=istart,mstep
      do irk=1,4
        ! calculate test like-particle collisions
        call collision  ! collision is on when icollision>0 and calculated every icollision-step at irk=2 (This might be a typo, collision is calculated when irk==1)
        call timer(t0,dt,t0wc,dtwc)
        time(8)=time(8)+dt
        timewc(8)=timewc(8)+dtwc

        ! Solve GK equation of motion using Runge-Kutta 4th order method
        idiag=mod(irk+3,4)+mod(istep,ndiag) ! idaig=0 when irk==4 and istep is a multiple of ndiag
        call pushion_RK4(irk)
        call timer(t0,dt,t0wc,dtwc)
        time(1)=time(1)+dt
        timewc(1)=timewc(1)+dtwc

        ! redistribute ion across PEs
        call shifti 
        call timer(t0,dt,t0wc,dtwc)
        time(2)=time(2)+dt
        timewc(2)=timewc(2)+dtwc

        ! Calculate ion weight damping 
        if(wimcompressibility==2) call weight_damping_ion
  
        call check_tracer_i
        ! Calculate ion perturbed density "delta-ni" and ion locations in the 3D field
        ! mesh, which are ued in CURRENT_EON for calculating gamma_i as well as
        ! in PUSHI for calculating electric force on ions
        call charge_ion
        call timer(t0,dt,t0wc,dtwc)
        time(3)=time(3)+dt
        timewc(3)=timewc(3)+dtwc
  
  
        if(ifelectron/=0) then    ! include kinetic electrons 
          ! NOTE: nelectron must take 1 or even number due to the use of 2nd order Runge-Kutta

          ! calculate electron collisions
          call collision_eon  ! collision is on when icollision>0 and calculated every icollision-step at irk=2 (This might be a typo, collision is calculated when irk==1)
          call timer(t0,dt,t0wc,dtwc)
          time(9)=time(9)+dt
          timewc(9)=timewc(9)+dtwc

          idiag_e=mod(irk+3,4)+mod(istep,ndiag) ! idaig=0 when irk==4 and istep is a multiple of ndiag
          call pusheon_RK4(solver1, masssolver, irk)
          call timer(t0,dt,t0wc,dtwc)
          time(10)=time(10)+dt
          timewc(10)=timewc(10)+dtwc

          call shifte
          call timer(t0,dt,t0wc,dtwc)
          time(11)=time(11)+dt
          timewc(11)=timewc(11)+dtwc

          call check_tracer_e
            
        endif ! end if(ifelectron/=0)
  
        if(wimcompressibility==2) call weight_damping_eon
  
        ! Calculate perturbed density and current of nonadiabatic electrons "delta-ne^h"
        ! and electron locations in the 3D field mesh,  which are needed in CURRENT_EON
        ! for calculating gamma_e
        call charge_eon
        call timer(t0,dt,t0wc,dtwc)
        time(12)=time(12)+dt
        timewc(12)=timewc(12)+dtwc
  
        ! solve GK Poisson equation using adiabatic electron
        call poisson( solver1, masssolver )
        !call poisson( solver_dphi_dt ) 
        call timer(t0,dt,t0wc,dtwc)
        time(4)=time(4)+dt
        timewc(4)=timewc(4)+dtwc

        ! smoothing,field diagnosis and finite difference for e_field
        call smooth(0)
        call timer(t0,dt,t0wc,dtwc)
        time(5)=time(5)+dt
        timewc(5)=timewc(5)+dtwc

        !idiag=mod(irk,2)+mod(istep+1,ndiag)

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
        endif
      enddo
 
      ! profile snapshots, write particle information to restart file
      !  phi_tmp(:)=phi_tmp(:)+phi(1,:) ! a time average
      call Cumulate_Some_Variables_for_Measure
 
      if(mod(istep,msnap) .eq. 0) then 
        !  phi_tmp=phi_tmp/real(msnap,wp)
        call snapshot
        call snapshot_hdf5
        call SetZeros_Cumulative_Variables
        !  phi_tmp=0.0_wp
        call check_abnormal_weight_electrons
        call check_abnormal_weight_ions
      end if
      
      if(mype==0) then
        write(stdout,*)'istep=',istep
        call flush(stdout)
      endif

      CALL DIAGNOSIS

      ! Find a specific tracer and print Tracer Information
      ! This routine may need to use multi-tracers.. etracers_labels, itracers_labels
      ! For example, zeon(1:6, etracer_labels(etracer_num/2))
       call timer(t0,dt,t0wc,dtwc) 
 
   enddo      !    MAIN TIME LOOP ----- MAIN TIME LOOP ----- MAIN TIME LOOP ----- MAIN TIME LOOP

  else
    write(stdout,*) 'RK4_Flag and associate settings might be wrong. Check the flags!!'
    call flush(stdout)
  endif


  loop_time=t0wc-loop_time

  if(mype==0) write(stdout,*) 'end of main time loop'

  if(mype==0) then
!     call view_distribution(1)
  end if

  call timer(t0,dt,t0wc,dtwc)
  time(nc)=t0-time(nc)
  timewc(nc)=t0wc-timewc(nc)
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
  ic(13)='Total'
  if(mype==0)then
     write(stdout,'(3a16)')' ROUTINES       ',' CPU TIME(sec)',' WALL CLOCK TIME'
     write(stdout,'(3a16)')'==========      ',' =============',' ==============='
     do i=1,nc
        write(stdout,'(a16,2f16.6)')ic(i),time(i),timewc(i)
     enddo
     write(stdout,'(/,"MAIN LOOP TIME(SEC):",f15.6,/)')loop_time
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

