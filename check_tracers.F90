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

