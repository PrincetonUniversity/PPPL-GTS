subroutine check_tracer_e
  use precision
  use mpi_defs
  use global_parameters
  use particle_array
  use esi_interface
  implicit none

  integer ierror,istatus(MPI_STATUS_SIZE)
  real etracer,tracer(7)
  integer my_mpi_id
 
  real(wp), dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj

! send tracer information to PE=0
  tracer=0.0_wp

  if(ntracer_e>0) then
     tracer(1:7)=zeon(1:7,ntracer_e)
     my_mpi_id=mype
     if(mype.ne.0) then 
        call MPI_SEND(tracer,7,mpi_Rsize,0,1,MPI_COMM_WORLD,ierror)
        call MPI_SEND(my_mpi_id,1,mpi_integer,0,1,MPI_COMM_WORLD,ierror)
     end if
  endif
  if(mype==0) then
     if(ntracer_e==0) then 
        call MPI_RECV(tracer,7,mpi_Rsize,MPI_ANY_SOURCE,1,&
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
                     abs(1.0_wp-etracer/etracer_e0), etracer, my_mpi_id
     call flush(89)
  endif

102 format(1x,i7,6(1x,ES14.6),1x,i6)

end subroutine check_tracer_e


subroutine check_tracer_i
  use precision
  use mpi_defs
  use global_parameters
  use particle_array
  use esi_interface
  implicit none

  integer ierror,istatus(MPI_STATUS_SIZE)
  real etracer,tracer(7)
  integer my_mpi_id

  real, dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj

! send tracer information to PE=0
  tracer=0.0_wp

  if(ntracer>0) then
     tracer(1:7)=zion(1:7,ntracer)
     my_mpi_id=mype
     if(mype.ne.0) then 
        call MPI_SEND(tracer,7,mpi_Rsize,0,1,MPI_COMM_WORLD,ierror)
        call MPI_SEND(my_mpi_id,1,mpi_integer,0,1,MPI_COMM_WORLD,ierror)
     end if
  endif
  if(mype==0) then
     if(ntracer==0) then 
        call MPI_RECV(tracer,7,mpi_Rsize,MPI_ANY_SOURCE,1,&
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
                      abs(1.0_wp-etracer/etracer_i0), etracer, my_mpi_id
     call flush(189)
  endif

102 format(1x,i7,6(1x,ES14.6),1x,i6)

end subroutine check_tracer_i


subroutine check_abnormal_weight_electrons

   use precision
   use mpi_defs
   use global_parameters
   use particle_array
   use particle_decomp
   implicit none

   integer ierr,istatus(MPI_STATUS_SIZE)
   integer m,i, ALL_PE_me
   real(wp) MY_PE_weight_sum, ALL_PE_weight_sum
   real(wp) MY_PE_abnormal_weight_sum, ALL_PE_abnormal_weight_sum
   real(wp), dimension(5) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
   integer MY_PTL_indices(memax)
   real(wp), dimension(:,:), allocatable:: send_buf_PTL, recv_buf_PTL
   integer each_PTL_length, num_of_Send_PTL, num_of_Recv_PTL
   integer recv_PTL_counts(numberpe), recv_PTL_displs(numberpe)
   integer MPI_PTL_Type ! PTL = TRacer Information

   real(wp) weight, energy
   character(len=30) file_name

   real(wp) w_Thold_min, w_Thold_max

   each_PTL_length = 19 ! including we_adiab + we_fdamp

   w_Thold_min = -50.0
   w_Thold_max = 50.0

   MY_PE_weight_sum=0
   MY_PE_abnormal_weight_sum=0
   num_of_Send_PTL=0
   do m=1,me
      weight = zeon(5,m)
      MY_PE_weight_sum = MY_PE_weight_sum +weight
      if(weight<w_Thold_min .or. weight>w_Thold_max) then
         num_of_Send_PTL = num_of_Send_PTL +1
         MY_PTL_indices(num_of_Send_PTL) = m
         MY_PE_abnormal_weight_sum = MY_PE_abnormal_weight_sum + weight
      endif
   enddo

   call MPI_REDUCE(me, ALL_PE_me, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
   call MPI_REDUCE(MY_PE_weight_sum, ALL_PE_weight_sum, 1, mpi_Rsize, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
   call MPI_REDUCE(MY_PE_abnormal_weight_sum, ALL_PE_abnormal_weight_sum, 1, mpi_Rsize, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
   call MPI_REDUCE(num_of_Send_PTL, num_of_Recv_PTL, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
   call MPI_GATHER(num_of_Send_PTL, 1, MPI_INTEGER, recv_PTL_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   allocate( send_buf_PTL(each_PTL_length, num_of_Send_PTL))

   if(mype==0) then
      recv_PTL_displs(1) = 0
      do i=2, numberpe
         recv_PTL_displs(i) = recv_PTL_displs(i-1) + recv_PTL_counts(i-1)
      enddo
      allocate(recv_buf_PTL(each_PTL_length, num_of_Recv_PTL))
   endif

   102 format(999(1x,ES16.8))

   call MPI_Type_contiguous(each_PTL_length, MPI_Rsize, MPI_PTL_Type, ierr);
   call MPI_Type_commit(MPI_PTL_Type, ierr);

   ! check electrons
   do i=1,num_of_Send_PTL

      m = MY_PTL_indices(i)

      aaa(1)=zeon(1,m)
      the(1)=zeon(2,m)
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      energy=0.5_wp*(zeon(4,m)*Bj(1))**2*m0ms*csc0**2+zeon(6,m)**2*Bj(1)  ! +epotential

      send_buf_PTL(1,i) = istep
      send_buf_PTL(2,i) = rj(1)
      send_buf_PTL(3,i) = zj(1)
      send_buf_PTL(4:10,i) = zeon(:,m)
      send_buf_PTL(11,i) = energy
      !  send_buf_PTL(12,i) = elabels(m)
      send_buf_PTL(12,i) = 0
      send_buf_PTL(13,i) = mype

      send_buf_PTL(14,i) = 0.0
      send_buf_PTL(15,i) = 0.0
      send_buf_PTL(16,i) = 0.0
      send_buf_PTL(17,i) = 0.0
      send_buf_PTL(18,i) = 0.0
      send_buf_PTL(19,i) = 0.0

      !  send_buf_PTL(14,i) = we_init(m)
      !  send_buf_PTL(15,i) = we_vpara(m)
      !  send_buf_PTL(16,i) = we_vperp(m)
      !  send_buf_PTL(17,i) = we_accel(m)
      !  send_buf_PTL(18,i) = we_adiab(m)
      !  send_buf_PTL(19,i) = we_fdamp(m)
   enddo

   ! Gather tracers_information to Root Rank.
   call MPI_Gatherv(send_buf_PTL, num_of_Send_PTL, MPI_PTL_Type, recv_buf_PTL, recv_PTL_counts, recv_PTL_displs, MPI_PTL_Type, 0, MPI_COMM_WORLD, ierr)

   if (mype == 0 ) then
      open(unit=190321, file="Abnormal_Electrons_Num.dat", position='append',action='write')
      ! write(190321,"(3(1x,I12),2(1x,ES16.8))") istep, ALL_PE_me, num_of_Recv_PTL, ALL_PE_weight_sum, ALL_PE_abnormal_weight_sum
      write(190321,102) real(istep), real(ALL_PE_me), real(num_of_Recv_PTL), ALL_PE_weight_sum, ALL_PE_abnormal_weight_sum
      call flush(190321)
      close(190321)
      if(num_of_Recv_PTL>0) then
         open(unit=190320, file="Abnormal_Electrons_Weights.dat", position='append',action='write')
         do i=1, num_of_Recv_PTL
            write(190320,102) recv_buf_PTL(:,i)
         enddo
         call flush(190320)
         close(190320)
      endif
   endif

   deallocate( send_buf_PTL)
   if(mype==0) then
      deallocate( recv_buf_PTL)
   endif

end subroutine check_abnormal_weight_electrons


subroutine check_abnormal_weight_ions

   use precision
   use mpi_defs
   use global_parameters
   use particle_array
   use particle_decomp
   implicit none

   integer ierr,istatus(MPI_STATUS_SIZE)
   integer m,i, ALL_PE_mi
   real(wp) MY_PE_weight_sum, ALL_PE_weight_sum
   real(wp) MY_PE_abnormal_weight_sum, ALL_PE_abnormal_weight_sum
   real(wp), dimension(5) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
   integer MY_PTL_indices(mimax)
   real(wp), dimension(:,:), allocatable:: send_buf_PTL, recv_buf_PTL
   integer each_PTL_length, num_of_Send_PTL, num_of_Recv_PTL
   integer recv_PTL_counts(numberpe), recv_PTL_displs(numberpe)
   integer MPI_PTL_Type ! PTL = TRacer Information

   real(wp) w_Thold_min, w_Thold_max
   real(wp) weight, energy
   character(len=30) file_name


   w_Thold_min = -50.0
   w_Thold_max = 50.0

   each_PTL_length = 18 !

   MY_PE_weight_sum=0
   MY_PE_abnormal_weight_sum=0
   num_of_Send_PTL=0
   do m=1,mi
      weight = zion(5,m)
      MY_PE_weight_sum = MY_PE_weight_sum +weight
      if(weight<w_Thold_min .or. weight>w_Thold_max) then
         num_of_Send_PTL = num_of_Send_PTL +1
         MY_PTL_indices(num_of_Send_PTL) = m
         MY_PE_abnormal_weight_sum = MY_PE_abnormal_weight_sum + weight
      endif
   enddo

   call MPI_REDUCE(mi, ALL_PE_mi, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
   call MPI_REDUCE(MY_PE_weight_sum, ALL_PE_weight_sum, 1, mpi_Rsize, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
   call MPI_REDUCE(MY_PE_abnormal_weight_sum, ALL_PE_abnormal_weight_sum, 1, mpi_Rsize, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
   call MPI_REDUCE(num_of_Send_PTL, num_of_Recv_PTL, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
   call MPI_GATHER(num_of_Send_PTL, 1, MPI_INTEGER, recv_PTL_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   allocate( send_buf_PTL(each_PTL_length, num_of_Send_PTL))

   if(mype==0) then
      recv_PTL_displs(1) = 0
      do i=2, numberpe
         recv_PTL_displs(i) = recv_PTL_displs(i-1) + recv_PTL_counts(i-1)
      enddo
      allocate(recv_buf_PTL(each_PTL_length, num_of_Recv_PTL))
   endif



   102 format(999(1x,ES16.8))

   call MPI_Type_contiguous(each_PTL_length, MPI_Rsize, MPI_PTL_Type, ierr);
   call MPI_Type_commit(MPI_PTL_Type, ierr);

   ! check electrons
   do i=1,num_of_Send_PTL

      m = MY_PTL_indices(i)

      aaa(1)=zion(1,m)
      the(1)=zion(2,m)
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      energy=0.5_wp*(zion(4,m)*Bj(1))**2+zion(6,m)**2*Bj(1)  ! +epotential

      send_buf_PTL(1,i) = istep
      send_buf_PTL(2,i) = rj(1)
      send_buf_PTL(3,i) = zj(1)
      send_buf_PTL(4:10,i) = zion(:,m)
      send_buf_PTL(11,i) = energy
      !  send_buf_PTL(12,i) = ilabels(m)
      send_buf_PTL(12,i) = 0
      send_buf_PTL(13,i) = mype

      send_buf_PTL(14,i) = 0.0
      send_buf_PTL(15,i) = 0.0
      send_buf_PTL(16,i) = 0.0
      send_buf_PTL(17,i) = 0.0
      send_buf_PTL(18,i) = 0.0

      !  send_buf_PTL(14,i) = wi_init(m)
      !  send_buf_PTL(15,i) = wi_vpara(m)
      !  send_buf_PTL(16,i) = wi_vperp(m)
      !  send_buf_PTL(17,i) = wi_accel(m)
      !  send_buf_PTL(18,i) = wi_fdamp(m)
   enddo

   ! Gather tracers_information to Root Rank.
   call MPI_Gatherv(send_buf_PTL, num_of_Send_PTL, MPI_PTL_Type, recv_buf_PTL, recv_PTL_counts, recv_PTL_displs, MPI_PTL_Type, 0, MPI_COMM_WORLD, ierr)

   if (mype == 0) then
      open(unit=190321, file="Abnormal_Ions_Num.dat", position='append',action='write')
      ! write(190321,"(3(1x,I12),2(1x,ES16.8))") istep, ALL_PE_mi, num_of_Recv_PTL, ALL_PE_weight_sum, ALL_PE_abnormal_weight_sum
      write(190321,102) real(istep), real(ALL_PE_mi), real(num_of_Recv_PTL), ALL_PE_weight_sum, ALL_PE_abnormal_weight_sum
      call flush(190321)
      close(190321)
      if(num_of_Recv_PTL>0) then
         open(unit=190320, file="Abnormal_Ions_Weights.dat", position='append',action='write')
         do i=1, num_of_Recv_PTL
            write(190320,102) recv_buf_PTL(:,i)
         enddo
         call flush(190320)
         close(190320)
      endif
   endif

   deallocate( send_buf_PTL)
   if(mype==0) then
      deallocate( recv_buf_PTL)
   endif

end subroutine check_abnormal_weight_ions
