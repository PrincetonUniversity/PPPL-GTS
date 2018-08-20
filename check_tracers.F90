subroutine check_tracer_e
  use precision
  use mpi_defs
  use global_parameters
  use particle_array
  use esi_interface
  use particle_decomp
  implicit none

  logical tf
  integer ierror,istatus(MPI_STATUS_SIZE)
  real etracer,tracer(7), tracer0(7), tracer_at_0(7)
  integer m,i, kk
  integer my_mpi_id, ierr
  real(wp), dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
  integer my_num_of_tracers, my_tracer_indices(max(me,mi))
  real normalized_zeta, normalized_zeta0, acute_zeta, acute_zeta0, ratio_for_tracer, ratio_for_tracer0

  character(len=30) file_name


  if(size(zeon,1) /= size(zeon0,1)) then
    write(*,*) "size(zeon,1)=", size(zeon,1), "size(zeon0,1)=", size(zeon0,1), " is differnet!!!!!"
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif



! send tracer information to PE=0
  tracer=0.0_wp
  tracer0=0.0_wp
  if(ntracer_e>0) then
     tracer(1:7)=zeon(1:7,ntracer_e)
     tracer0(1:7)=zeon0(1:7,ntracer_e)
     my_mpi_id=mype
     if(mype.ne.0) then 
        call MPI_SEND(tracer,7,mpi_Rsize,0,1,MPI_COMM_WORLD,ierror)
        call MPI_SEND(tracer0,7,mpi_Rsize,0,2,MPI_COMM_WORLD,ierror)
        call MPI_SEND(my_mpi_id,1,mpi_integer,0,1,MPI_COMM_WORLD,ierror)
     end if
  endif
  if(mype==0) then
     if(ntracer_e==0) then 
        call MPI_RECV(tracer,7,mpi_Rsize,MPI_ANY_SOURCE,1,&
                                    MPI_COMM_WORLD,istatus,ierror)
        call MPI_RECV(tracer0,7,mpi_Rsize,MPI_ANY_SOURCE,2,&
                                    MPI_COMM_WORLD,istatus,ierror)                                    
        call MPI_RECV(my_mpi_id,1,mpi_integer,MPI_ANY_SOURCE,1,&
                                    MPI_COMM_WORLD,istatus,ierror)
     end if

! tracer particle energy and angular momentum
     aaa(1)=tracer(1)
     the(1)=tracer(2)
     call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)

     etracer=0.5_wp*(tracer(4)*Bj(1))**2*m0ms*csc0**2+tracer(6)**2*Bj(1)  ! +epotential

! record tracer information
     write(89,102) istep,(rj(1)-R_0), zj(1), tracer(2), tracer(3), &
                    ((etracer-ntracer_energy_electron0)/ntracer_energy_electron0), etracer, my_mpi_id
     call flush(89)

! For Poincar plot of tracer particles
     open(unit=899, file="tracer_Poincare_electron.dat", position='append',action='write')
     if ( ABS(tracer(3) - tracer0(3)) > 5.0 ) then
  ! record tracer information      
      write(899,*) istep, (rj(1)-R_0), zj(1), tracer(2), tracer(3), etracer, my_mpi_id
      call flush(899)
    endif  
    close(899)

  endif

  ! call system('mkdir Poincare_Plot')
  ! write(file_name,'("Poincare_Plot/electron_",i3.3,".dat")') mype

  ! open(unit=12345, file=file_name, position='append',action='write')

  ! do m=1,me

  !   if ( ABS(zeon(3,m) - zeon0(3,m)) > 6.0 ) then

  !     tracer(1:7)=zeon(1:7,m)
  !     aaa(1)=tracer(1)
  !     the(1)=tracer(2)
  !     call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)

  !     etracer=0.5_wp*(tracer(4)*Bj(1))**2*m0ms*csc0**2+tracer(6)**2*Bj(1)  ! +epotential

  ! ! record tracer information
  !     write(12345,*) istep, (rj(1)-R_0), zj(1), tracer(2), tracer(3), 0.0, etracer, mype
  !     ! write(12345,*) istep, (rj(1)-R_0), zj(1), mype
  !   end if
  ! enddo
  ! call flush(12345)
  ! close(12345)

102 format(1x,i7,6(1x,ES14.6),1x,i6)



  ! check tracer
  my_num_of_tracers = 0

  do m=1,me
    do i=1,num_of_etracers
      if (elabels(m) == etracers_labels(i)) then
        my_num_of_tracers = my_num_of_tracers + 1      
        my_tracer_indices(my_num_of_tracers) = m      
      endif
    enddo
  enddo

  do i=0,numberpe-1
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (my_num_of_tracers > 0 .AND. i == mype) then    

      open(unit=1515, file="real_tracers_electron.dat", position='append',action='write')
      open(unit=1516, file="real_tracers_electron_Poincare.dat", position='append',action='write')

      do m=1,my_num_of_tracers

        tracer(:) = zeon(:,my_tracer_indices(m))
        tracer0(:) = zeon0(:,my_tracer_indices(m))   
    
        aaa(1)=tracer(1)
        the(1)=tracer(2)
        call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
    
        etracer=0.5_wp*(tracer(4)*Bj(1))**2*m0ms*csc0**2+tracer(6)**2*Bj(1)  ! +epotential
        ! record tracer information
        if (MOD(istep, nstep_for_tracer_record) == 0) then
          ! write(1515,*) istep, (rj(1)-R_0), zj(1), tracer(2), tracer(3), etracer, elabels(my_tracer_indices(m)),mype
          write(1515,*) istep, rj(1), zj(1), tracer(1:6), etracer, elabels(my_tracer_indices(m)),mype

          call flush(1515)
        endif

        if ( ABS(tracer(3) - tracer0(3)) > 5.0 ) then
          
          normalized_zeta =  tracer(3) - floor(tracer(3)/(2.0*pi)) * 2.0*pi
          normalized_zeta0 =  tracer0(3) - floor(tracer0(3)/(2.0*pi)) * 2.0*pi

          acute_zeta = min(normalized_zeta, 2*pi-normalized_zeta)
          acute_zeta0 = min(normalized_zeta0, 2*pi-normalized_zeta0)

          ratio_for_tracer = acute_zeta0/ (acute_zeta+ acute_zeta0)
          ratio_for_tracer0 = acute_zeta/ (acute_zeta+ acute_zeta0)
          
          do kk=1,size(zeon,1)
              tracer_at_0(kk) = ratio_for_tracer*tracer(kk) + ratio_for_tracer0*tracer0(kk)
          enddo

          ! new angle (theta and zeta) should be calculated using two angles whose difference is smaller than pi
          ! Angle of tracer0 is modified to satisfy the above condition.
          ! A = tracer angle, B = tracer0;   B' = B + 2*n*pi where n = floor((A-B+pi)/(2*pi));  
          do kk=2,3
            tracer_at_0(kk) = ratio_for_tracer*tracer(kk) + ratio_for_tracer0*  ( tracer0(kk) + 2.0*pi* floor( (tracer(kk) - tracer0(kk) + pi )/(2.0*pi)  ))            
            tracer_at_0(kk) = tracer_at_0(kk) - floor(tracer_at_0(kk)/(2.0*pi)) * 2.0*pi ! normalization
          end do



          aaa(1)=tracer_at_0(1)
          the(1)=tracer_at_0(2)
          call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
          etracer=0.5_wp*(tracer_at_0(4)*Bj(1))**2*m0ms*csc0**2+tracer_at_0(6)**2*Bj(1)  ! +epotential


          !write(1516,*) istep, (rj(1)-R_0), zj(1), tracer_at_0(2), tracer_at_0(3), etracer, elabels(my_tracer_indices(m)),mype
          write(1516,*) istep, rj(1), zj(1), tracer_at_0(1:6), etracer, elabels(my_tracer_indices(m)),mype
          call flush(1516)
        endif 
   
      enddo  
 
      close(1515)
      close(1516)    
    endif
  enddo



end subroutine check_tracer_e


subroutine check_tracer_i
  use precision
  use mpi_defs
  use global_parameters
  use particle_array
  use esi_interface
  use particle_decomp
  implicit none
  integer m,i,kk
  logical tf
  integer my_mpi_id
  integer ierror,istatus(MPI_STATUS_SIZE)
  real etracer,tracer(7),tracer0(7), tracer_at_0(7)
  real normalized_zeta, normalized_zeta0, acute_zeta, acute_zeta0, ratio_for_tracer, ratio_for_tracer0




  real, dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
  character(len=30) file_name
  
  integer my_num_of_tracers, my_tracer_indices(mi)

  
  
  ! send tracer information to PE=0
  tracer=0.0_wp

  if(ntracer>0) then
     tracer(1:7)=zion(1:7,ntracer)     
     tracer0(1:7)=zion0(1:7,ntracer)     
     my_mpi_id=mype
     if(mype.ne.0) then 
        call MPI_SEND(tracer,7,mpi_Rsize,0,1,MPI_COMM_WORLD,ierror)
        call MPI_SEND(tracer0,7,mpi_Rsize,0,2,MPI_COMM_WORLD,ierror)
        call MPI_SEND(my_mpi_id,1,mpi_integer,0,1,MPI_COMM_WORLD,ierror)
     end if
  endif
  if(mype==0) then
     if(ntracer==0) then 
        call MPI_RECV(tracer,7,mpi_Rsize,MPI_ANY_SOURCE,1,&
                                  MPI_COMM_WORLD,istatus,ierror)
        call MPI_RECV(tracer0,7,mpi_Rsize,MPI_ANY_SOURCE,2,&
                                  MPI_COMM_WORLD,istatus,ierror)                                  
        call MPI_RECV(my_mpi_id,1,mpi_integer,MPI_ANY_SOURCE,1,&
                                    MPI_COMM_WORLD,istatus,ierror)
     end if
! tracer particle energy and angular momentum
     aaa(1)=tracer(1)
     the(1)=tracer(2)
     call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)

     etracer=0.5_wp*(tracer(4)*Bj(1))**2+tracer(6)**2*Bj(1)  ! +epotential

! record tracer information
     write(189,102) istep, (rj(1)-R_0), zj(1), tracer(2), tracer(3), &
                      ((etracer-ntracer_energy_ion0)/ntracer_energy_ion0), etracer, my_mpi_id
     call flush(189)




! For Poincar plot of tracer particles
     open(unit=1899, file="tracer_Poincare_ion.dat", position='append',action='write')
     if ( ABS(tracer(3) - tracer0(3)) > 6.0 ) then
  ! record tracer information
      write(1899,*) istep, (rj(1)-R_0), zj(1), tracer(2), tracer(3), etracer, my_mpi_id
      call flush(1899)
    endif  
     close(1899)


  endif


  ! call system('mkdir Poincare_Plot')
  ! write(file_name,'("Poincare_Plot/ion_",i3.3,".dat")')mype
  ! open(unit=12346, file=file_name, position='append',action='write')
  ! do m=1,mi

  !   if ( ABS(zion(3,m) - zion0(3,m)) > 5.0 ) then
  !     tracer(1:7)=zion(1:7,m)
  !     aaa(1)=tracer(1)
  !     the(1)=tracer(2)
  !     call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)

  !     etracer=0.5_wp*(tracer(4)*Bj(1))**2*m0ms*csc0**2+tracer(6)**2*Bj(1)  ! +epotential

  ! ! record tracer information
  !     write(12346,*) istep, (rj(1)-R_0), zj(1), tracer(2), tracer(3), 0.0, etracer, mype
  !     ! write(12346,*) istep, (rj(1)-R_0), zj(1), mype
  !   endif 
  ! enddo
  ! call flush(12346)
  ! close(12346)


102 format(1x,i7,6(1x,ES14.6),1x,i6)


    ! check tracer
  my_num_of_tracers = 0

  do m=1,mi
    do i=1,num_of_itracers
      if (ilabels(m) == itracers_labels(i)) then

        my_num_of_tracers = my_num_of_tracers + 1      
        my_tracer_indices(my_num_of_tracers) = m      
      endif
    enddo
  enddo

  do i=0,numberpe-1
    call MPI_Barrier(MPI_COMM_WORLD, ierror)
    if (i == mype) then    

      open(unit=1515, file="real_tracers_ion.dat", position='append',action='write')
      open(unit=1516, file="real_tracers_ion_Poincare.dat", position='append',action='write')

      do m=1,my_num_of_tracers

        tracer(:) = zion(:,my_tracer_indices(m))
        tracer0(:) = zion0(:,my_tracer_indices(m))   
    
        aaa(1)=tracer(1)
        the(1)=tracer(2)
        call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
    
        etracer=0.5_wp*(tracer(4)*Bj(1))**2+tracer(6)**2*Bj(1)  ! +epotential

        ! record tracer information      
        if (MOD(istep, nstep_for_tracer_record) == 0) then
          !write(1515,*) istep, (rj(1)-R_0), zj(1), tracer(2), tracer(3), etracer, ilabels(my_tracer_indices(m)),mype
          write(1515,*) istep, rj(1), zj(1), tracer(1:6), etracer, ilabels(my_tracer_indices(m)),mype
          call flush(1515)
        endif
        if ( ABS(tracer(3) - tracer0(3)) > 5.0 ) then
          
          normalized_zeta =  tracer(3) - floor(tracer(3)/(2.0*pi)) * 2.0*pi
          normalized_zeta0 =  tracer0(3) - floor(tracer0(3)/(2.0*pi)) * 2.0*pi

          acute_zeta = min(normalized_zeta, 2*pi-normalized_zeta)
          acute_zeta0 = min(normalized_zeta0, 2*pi-normalized_zeta0)

          ratio_for_tracer = acute_zeta0/ (acute_zeta+ acute_zeta0)
          ratio_for_tracer0 = acute_zeta/ (acute_zeta+ acute_zeta0)
          
          do kk=1,size(zion,1)
              tracer_at_0(kk) = ratio_for_tracer*tracer(kk) + ratio_for_tracer0*tracer0(kk)
          enddo

          ! new angle (theta and zeta) should be calculated using two angles whose difference is smaller than pi
          ! Angle of tracer0 is modified to satisfy the above condition.
          ! A = tracer angle, B = tracer0;   B' = B + 2*n*pi where n = floor((A-B+pi)/(2*pi));  
          do kk=2,3
            tracer_at_0(kk) = ratio_for_tracer*tracer(kk) + ratio_for_tracer0*  ( tracer0(kk) + 2.0*pi* floor( (tracer(kk) - tracer0(kk) + pi )/(2.0*pi)  ))
          end do

          aaa(1)=tracer_at_0(1)
          the(1)=tracer_at_0(2)
          call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
          etracer=0.5_wp*(tracer_at_0(4)*Bj(1))**2+tracer_at_0(6)**2*Bj(1)  ! +epotential


          !write(1516,*) istep, (rj(1)-R_0), zj(1), tracer_at_0(2), tracer_at_0(3), etracer, ilabels(my_tracer_indices(m)),mype
          write(1516,*) istep, rj(1), zj(1), tracer_at_0(1:6), etracer, ilabels(my_tracer_indices(m)),mype
          call flush(1516)
        endif        
      enddo  
      
      close(1515)
      close(1516)    
    endif
  enddo


end subroutine check_tracer_i




! subroutine Am_I_tracer_e(m, tf) 
!   use global_parameters
!   use particle_decomp
!   implicit none
!   integer, intent(in) :: m
!   logical, intent(in) :: tf  
!   integer ori_me_per_processor
!   integer :: i,t_max 
!   integer, DIMENSION(:), allocatable :: tracer_e_list

!   t_max = 20


!   allocate( tracer_e_list(t_max) )
  
!   do i=1,t_max
!     tracer_e_list(i) = 1e5*((real(i)/tmax) * ori_me_per_processor) + 0
!   enddo

!   do i=1,t_max
!      if(elabels(m) == tracer_e_list(i)) then
!         tf = .true.
!         return
!     endif
!   enddo

!   tf = .false.
!   return 

!   deallocate(tracer_e_list)
  
! end subroutine Am_I_tracer_e


! subroutine Am_I_tracer_i(m, tf)
!   use global_parameters
!   use particle_decomp
!   implicit none
!   integer, intent(in) :: m
!   logical, intent(in) :: tf
!   integer ori_mi_per_processor
!   integer :: i,t_max
!   integer, DIMENSION(:), allocatable :: tracer_i_list



!   t_max = 20


!   allocate( tracer_i_list(t_max) )
  
!   do i=1,t_max
!     tracer_i_list(i) = 1e5*((real(i)/tmax) * ori_mi_per_processor) + 0
!   enddo

!   do i=1,t_max
!      if(ilabels(m) == tracer_i_list(i)) then
!         tf = .true.
!         return
!     endif
!   enddo

!   tf = .false.
!   return 

!   deallocate(tracer_i_list)
  
! end subroutine Am_I_tracer_i