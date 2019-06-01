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
   real tracer_energy, tracer(7), tracer0(7), tracer_at_ZeroPhi(7)
   integer m,i, kk
   integer my_mpi_id, ierr
   real(wp), dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
   integer my_tracer_indices(num_of_etracers)
   integer num_of_MY_tracers, recv_num_of_TRI, num_of_MY_tracers_Poincare, recv_num_of_TRI_Poincare
   !real send_buf_TRI(13, num_of_etracers), recv_buf_TRI(13, num_of_etracers), send_buf_TRI_Poincare(13, num_of_etracers), recv_buf_TRI_Poincare(13, num_of_etracers)
   real(wp), dimension(:,:), allocatable:: send_buf_TRI, recv_buf_TRI, send_buf_TRI_Poincare, recv_buf_TRI_Poincare
   integer each_TRI_length
   integer recv_TRI_counts(nproc_toroidal), recv_TRI_displs(nproc_toroidal), recv_TRI_Poincare_counts(nproc_toroidal), recv_TRI_Poincare_displs(nproc_toroidal)
   integer MPI_TRI_Type ! TRI = TRacer Information
   real normalized_zeta, normalized_zeta0, acute_zeta, acute_zeta0, ratio_for_tracer, ratio_for_tracer0

   character(len=30) file_name


   if(size(zeon,1) /= size(zeon0,1)) then
      write(*,*) "size(zeon,1)=", size(zeon,1), "size(zeon0,1)=", size(zeon0,1), " is differnet!!!!!"
      call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
   endif


   each_TRI_length = 21 ! including we_adiab + we_fdamp + 2 Unnormalized Angles

   allocate( send_buf_TRI(each_TRI_length, num_of_etracers), recv_buf_TRI(each_TRI_length, num_of_etracers), send_buf_TRI_Poincare(each_TRI_length, num_of_etracers), recv_buf_TRI_Poincare(each_TRI_length, num_of_etracers))

   ! 102 format(1x,i7,6(1x,ES16.8),1x,i6)
   !102 format(13(1x,ES16.8))
   102 format(999(1x,ES16.8))

   !call MPI_Type_contiguous(13, MPI_Rsize, MPI_TRI_Type, ierr);
   call MPI_Type_contiguous(each_TRI_length, MPI_Rsize, MPI_TRI_Type, ierr);
   call MPI_Type_commit(MPI_TRI_Type, ierr);

   ! check tracer
   num_of_MY_tracers = 0
   do m=1,me
      do i=1,num_of_etracers
         if (elabels(m) == etracers_labels(i)) then
            num_of_MY_tracers = num_of_MY_tracers + 1
            my_tracer_indices(num_of_MY_tracers) = m

            ! Print out a reference tracer information
            if(elabels(m) == etracers_labels((num_of_etracers/2))) then
               aaa(1)=zeon(1,m)
               the(1)=zeon(2,m)
               call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
               tracer_energy=0.5_wp*(zeon(4,m)*Bj(1))**2*m0ms*csc0**2+zeon(6,m)**2*Bj(1)
               write(stdout,*) 'istep=',istep,', KE=',tracer_energy*Tem_00, 'keV',', etracer(',num_of_etracers/2,') =', zeon(:,m)
               call flush(stdout)
            endif
         endif
      enddo
   enddo


   ! Real_tracers_electron.dat output routine
   if (MOD(istep, nstep_for_tracer_record) == 0) then

      call MPI_GATHER(num_of_MY_tracers, 1, MPI_INTEGER, recv_TRI_counts, 1, MPI_INTEGER, 0, toroidal_comm, ierror)

      if(mype == 0) then
         recv_TRI_displs(1) = 0
         recv_num_of_TRI = recv_TRI_counts(1)
         do i=2, nproc_toroidal
            recv_TRI_displs(i) = recv_TRI_displs(i-1) + recv_TRI_counts(i-1)
            recv_num_of_TRI = recv_num_of_TRI + recv_TRI_counts(i)
         enddo
      endif

      ! Save TRacers Information to send_buf_TRI
      do m=1, num_of_MY_tracers
         tracer(:) = zeon(:,my_tracer_indices(m))
         tracer0(:) = zeon0(:,my_tracer_indices(m))
         aaa(1)=tracer(1)
         the(1)=tracer(2)
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
         tracer_energy=0.5_wp*(tracer(4)*Bj(1))**2*m0ms*csc0**2+tracer(6)**2*Bj(1)  ! +epotential

         send_buf_TRI(1,m) = istep
         send_buf_TRI(2,m) = rj(1)
         send_buf_TRI(3,m) = zj(1)
         send_buf_TRI(4:10,m) = tracer(1:7)
         send_buf_TRI(11,m) = tracer_energy
         send_buf_TRI(12,m) = elabels(my_tracer_indices(m))
         send_buf_TRI(13,m) = mype
         send_buf_TRI(14,m) = we_init(my_tracer_indices(m))
         send_buf_TRI(15,m) = we_vpara(my_tracer_indices(m))
         send_buf_TRI(16,m) = we_vperp(my_tracer_indices(m))
         send_buf_TRI(17,m) = we_accel(my_tracer_indices(m))
         send_buf_TRI(18,m) = we_adiab(my_tracer_indices(m))
         send_buf_TRI(19,m) = we_fdamp(my_tracer_indices(m))
         send_buf_TRI(20,m) = UnNorm_theta_e(my_tracer_indices(m))
         send_buf_TRI(21,m) = UnNorm_zeta_e(my_tracer_indices(m))
      enddo

      ! Gather tracers_information to Root Rank.
      call MPI_Gatherv(send_buf_TRI, num_of_MY_tracers, MPI_TRI_Type, recv_buf_TRI, recv_TRI_counts, recv_TRI_displs, MPI_TRI_Type, 0, toroidal_comm, ierror)

      if (mype == 0) then
         open(unit=1515, file="multi_tracers_electron.dat", position='append',action='write')
         do i=1, recv_num_of_TRI
            !write(1515,*) recv_buf_TRI(:,i)
            write(1515,102) recv_buf_TRI(:,i)
         enddo
         call flush(1515)
         close(1515)
      endif
   endif


   ! Check particles passing a specific poloidal plane (phi=0) for a Poincare plot
   ! Real_tracers_electron_Poincare.dat output routine
   num_of_MY_tracers_Poincare = 0
   recv_num_of_TRI_Poincare = 0
   do m=1, num_of_MY_tracers
      tracer(:) = zeon(:,my_tracer_indices(m))
      tracer0(:) = zeon0(:,my_tracer_indices(m))

      if ( ABS(tracer(3) - tracer0(3)) > 5.0 ) then

         normalized_zeta =  tracer(3) - floor(tracer(3)/(2.0*pi)) * 2.0*pi
         normalized_zeta0 =  tracer0(3) - floor(tracer0(3)/(2.0*pi)) * 2.0*pi
         acute_zeta = min(normalized_zeta, 2*pi-normalized_zeta)
         acute_zeta0 = min(normalized_zeta0, 2*pi-normalized_zeta0)
         ratio_for_tracer = acute_zeta0/ (acute_zeta+ acute_zeta0)
         ratio_for_tracer0 = acute_zeta/ (acute_zeta+ acute_zeta0)

         do kk=1,size(zeon,1)
            tracer_at_ZeroPhi(kk) = ratio_for_tracer*tracer(kk) + ratio_for_tracer0*tracer0(kk)
         enddo

         ! tracer(6) = sqrt(magnetic moment) is constant during simulation except collisional process.
         ! Therefore, tracer(6) is always enough to deal with magnetic moment mu = tracer(6)*tracer(6).
         ! In this sense, tracer0(6) is not used for sqrt(mu), but it is set by value 1 (in loading routine) for other purpose.
         ! We need to recover tracer_at_ZeroPhi(6) using tracer(6) only, not including tracer0(6) which is not sqrt(mu).
         tracer_at_ZeroPhi(6) = tracer(6)

         ! new angle (theta and zeta) should be calculated using two angles whose difference is smaller than pi
         ! Angle of tracer0 is modified to satisfy the above condition.
         ! A = tracer angle, B = tracer0;   B' = B + 2*n*pi where n = floor((A-B+pi)/(2*pi));
         do kk=2,3
            tracer_at_ZeroPhi(kk) = ratio_for_tracer*tracer(kk) + ratio_for_tracer0*  ( tracer0(kk) + 2.0*pi* floor( (tracer(kk) - tracer0(kk) + pi )/(2.0*pi)  ))
            tracer_at_ZeroPhi(kk) = tracer_at_ZeroPhi(kk) - floor(tracer_at_ZeroPhi(kk)/(2.0*pi)) * 2.0*pi ! normalization
         end do

         aaa(1)=tracer_at_ZeroPhi(1)
         the(1)=tracer_at_ZeroPhi(2)
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
         tracer_energy=0.5_wp*(tracer_at_ZeroPhi(4)*Bj(1))**2*m0ms*csc0**2+tracer_at_ZeroPhi(6)**2*Bj(1)  ! +epotential

         ! Save tracers_Poincare Information to send_buf_TRI_Poincare
         num_of_MY_tracers_Poincare = num_of_MY_tracers_Poincare + 1

         send_buf_TRI_Poincare(1,num_of_MY_tracers_Poincare) = istep
         send_buf_TRI_Poincare(2,num_of_MY_tracers_Poincare) = rj(1)
         send_buf_TRI_Poincare(3,num_of_MY_tracers_Poincare) = zj(1)
         send_buf_TRI_Poincare(4:10,num_of_MY_tracers_Poincare) = tracer_at_ZeroPhi(1:7)
         send_buf_TRI_Poincare(11,num_of_MY_tracers_Poincare) = tracer_energy
         send_buf_TRI_Poincare(12,num_of_MY_tracers_Poincare) = elabels(my_tracer_indices(m))
         send_buf_TRI_Poincare(13,num_of_MY_tracers_Poincare) = mype
         send_buf_TRI_Poincare(14,num_of_MY_tracers_Poincare) = we_init(my_tracer_indices(m))
         send_buf_TRI_Poincare(15,num_of_MY_tracers_Poincare) = we_vpara(my_tracer_indices(m))
         send_buf_TRI_Poincare(16,num_of_MY_tracers_Poincare) = we_vperp(my_tracer_indices(m))
         send_buf_TRI_Poincare(17,num_of_MY_tracers_Poincare) = we_accel(my_tracer_indices(m))
         send_buf_TRI_Poincare(18,num_of_MY_tracers_Poincare) = we_adiab(my_tracer_indices(m))
         send_buf_TRI_Poincare(19,num_of_MY_tracers_Poincare) = we_fdamp(my_tracer_indices(m))
         send_buf_TRI_Poincare(20,num_of_MY_tracers_Poincare) = UnNorm_theta_e(my_tracer_indices(m))
         send_buf_TRI_Poincare(21,num_of_MY_tracers_Poincare) = UnNorm_zeta_e(my_tracer_indices(m))
      endif
   enddo

   call MPI_GATHER(num_of_MY_tracers_Poincare, 1, MPI_INTEGER, recv_TRI_Poincare_counts, 1, MPI_INTEGER, 0, toroidal_comm, ierror)

   if(mype == 0) then
      recv_TRI_Poincare_displs(1) = 0
      recv_num_of_TRI_Poincare = recv_TRI_Poincare_counts(1)
      do i=2, nproc_toroidal
         recv_TRI_Poincare_displs(i) = recv_TRI_Poincare_displs(i-1) + recv_TRI_Poincare_counts(i-1)
         recv_num_of_TRI_Poincare = recv_num_of_TRI_Poincare + recv_TRI_Poincare_counts(i)
      enddo
   endif

   call MPI_Gatherv(send_buf_TRI_Poincare, num_of_MY_tracers_Poincare, MPI_TRI_Type, recv_buf_TRI_Poincare, recv_TRI_Poincare_counts, recv_TRI_Poincare_displs, MPI_TRI_Type, 0, toroidal_comm, ierror)

   if (mype == 0) then
      open(unit=1516, file="multi_tracers_electron_Poincare.dat", position='append',action='write')
      do i=1, recv_num_of_TRI_Poincare
         !write(1516,*) recv_buf_TRI_Poincare(:,i)
         write(1516,102) recv_buf_TRI_Poincare(:,i)
      enddo
      call flush(1516)
      close(1516)
   endif


   deallocate( send_buf_TRI, recv_buf_TRI, send_buf_TRI_Poincare, recv_buf_TRI_Poincare)

end subroutine check_tracer_e





subroutine check_tracer_i
   use precision
   use mpi_defs
   use global_parameters
   use particle_array
   use esi_interface
   use particle_decomp
   implicit none

   logical tf
   integer ierror,istatus(MPI_STATUS_SIZE)
   real tracer_energy, tracer(7), tracer0(7), tracer_at_ZeroPhi(7)
   integer m,i, kk
   integer my_mpi_id, ierr
   real(wp), dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
   integer my_tracer_indices(num_of_itracers)
   integer num_of_MY_tracers, recv_num_of_TRI, num_of_MY_tracers_Poincare, recv_num_of_TRI_Poincare
   !real send_buf_TRI(13, num_of_itracers), recv_buf_TRI(13, num_of_itracers), send_buf_TRI_Poincare(13, num_of_itracers), recv_buf_TRI_Poincare(13, num_of_itracers)
   real(wp), dimension(:,:), allocatable:: send_buf_TRI, recv_buf_TRI, send_buf_TRI_Poincare, recv_buf_TRI_Poincare
   integer each_TRI_length
   integer recv_TRI_counts(nproc_toroidal), recv_TRI_displs(nproc_toroidal), recv_TRI_Poincare_counts(nproc_toroidal), recv_TRI_Poincare_displs(nproc_toroidal)
   integer MPI_TRI_Type ! TRI = TRacer Information
   real normalized_zeta, normalized_zeta0, acute_zeta, acute_zeta0, ratio_for_tracer, ratio_for_tracer0

   character(len=30) file_name


   if(size(zion,1) /= size(zion0,1)) then
      write(*,*) "size(zion,1)=", size(zion,1), "size(zion0,1)=", size(zion0,1), " is differnet!!!!!"
      call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
   endif

   each_TRI_length = 20 ! + 2 Unnormalized Angles
   allocate( send_buf_TRI(each_TRI_length, num_of_itracers), recv_buf_TRI(each_TRI_length, num_of_itracers), send_buf_TRI_Poincare(each_TRI_length, num_of_itracers), recv_buf_TRI_Poincare(each_TRI_length, num_of_itracers))


   ! 102 format(1x,i7,6(1x,ES16.8),1x,i6)
   ! 102 format(13(1x,ES16.8))
   102 format(999(1x,ES16.8))
   !call MPI_Type_contiguous(13, MPI_Rsize, MPI_TRI_Type, ierr);
   call MPI_Type_contiguous(each_TRI_length, MPI_Rsize, MPI_TRI_Type, ierr);
   call MPI_Type_commit(MPI_TRI_Type, ierr);

   ! check tracer
   num_of_MY_tracers = 0
   do m=1,mi
      do i=1,num_of_itracers
         if (ilabels(m) == itracers_labels(i)) then
            num_of_MY_tracers = num_of_MY_tracers + 1
            my_tracer_indices(num_of_MY_tracers) = m

            ! Print out a reference tracer information
            if(ilabels(m) == itracers_labels((num_of_itracers/2))) then
               aaa(1)=zion(1,m)
               the(1)=zion(2,m)
               call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
               tracer_energy=0.5_wp*(zion(4,m)*Bj(1))**2+zion(6,m)**2*Bj(1)  ! Ion energy

               write(stdout,*) 'istep=', istep,', KE=',tracer_energy*Tem_00, 'keV',', itracer(',num_of_itracers/2,') =', zion(:,m)
               call flush(stdout)
            endif
         endif
      enddo
   enddo


   ! Real_tracers_ion.dat output routine
   if (MOD(istep, nstep_for_tracer_record) == 0) then

      call MPI_GATHER(num_of_MY_tracers, 1, MPI_INTEGER, recv_TRI_counts, 1, MPI_INTEGER, 0, toroidal_comm, ierror)

      if(mype == 0) then
         recv_TRI_displs(1) = 0
         recv_num_of_TRI = recv_TRI_counts(1)
         do i=2, nproc_toroidal
            recv_TRI_displs(i) = recv_TRI_displs(i-1) + recv_TRI_counts(i-1)
            recv_num_of_TRI = recv_num_of_TRI + recv_TRI_counts(i)
         enddo
      endif

      ! Save TRacers Information to send_buf_TRI
      do m=1, num_of_MY_tracers
         tracer(:) = zion(:,my_tracer_indices(m))
         tracer0(:) = zion0(:,my_tracer_indices(m))
         aaa(1)=tracer(1)
         the(1)=tracer(2)
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
         tracer_energy=0.5_wp*(tracer(4)*Bj(1))**2+tracer(6)**2*Bj(1)  ! Ion energy
         send_buf_TRI(1,m) = istep
         send_buf_TRI(2,m) = rj(1)
         send_buf_TRI(3,m) = zj(1)
         send_buf_TRI(4:10,m) = tracer(1:7)
         send_buf_TRI(11,m) = tracer_energy
         send_buf_TRI(12,m) = ilabels(my_tracer_indices(m))
         send_buf_TRI(13,m) = mype
         send_buf_TRI(14,m) = wi_init(my_tracer_indices(m))
         send_buf_TRI(15,m) = wi_vpara(my_tracer_indices(m))
         send_buf_TRI(16,m) = wi_vperp(my_tracer_indices(m))
         send_buf_TRI(17,m) = wi_accel(my_tracer_indices(m))
         send_buf_TRI(18,m) = wi_fdamp(my_tracer_indices(m))
         send_buf_TRI(19,m) = UnNorm_theta_i(my_tracer_indices(m))
         send_buf_TRI(20,m) = UnNorm_zeta_i(my_tracer_indices(m))
      enddo

      ! Gather tracers_information to Root Rank.
      call MPI_Gatherv(send_buf_TRI, num_of_MY_tracers, MPI_TRI_Type, recv_buf_TRI, recv_TRI_counts, recv_TRI_displs, MPI_TRI_Type, 0, toroidal_comm, ierror)

      if (mype == 0) then
         open(unit=1515, file="multi_tracers_ion.dat", position='append',action='write')
         do i=1, recv_num_of_TRI
            !write(1515,*) recv_buf_TRI(:,i)
            write(1515,102) recv_buf_TRI(:,i)
         enddo
         call flush(1515)
         close(1515)
      endif
   endif


   ! Check particles passing a specific poloidal plane (phi=0) for a Poincare plot
   ! Real_tracers_ion_Poincare.dat output routine
   num_of_MY_tracers_Poincare = 0
   recv_num_of_TRI_Poincare = 0
   do m=1, num_of_MY_tracers
      tracer(:) = zion(:,my_tracer_indices(m))
      tracer0(:) = zion0(:,my_tracer_indices(m))

      if ( ABS(tracer(3) - tracer0(3)) > 5.0 ) then

         normalized_zeta =  tracer(3) - floor(tracer(3)/(2.0*pi)) * 2.0*pi
         normalized_zeta0 =  tracer0(3) - floor(tracer0(3)/(2.0*pi)) * 2.0*pi
         acute_zeta = min(normalized_zeta, 2*pi-normalized_zeta)
         acute_zeta0 = min(normalized_zeta0, 2*pi-normalized_zeta0)
         ratio_for_tracer = acute_zeta0/ (acute_zeta+ acute_zeta0)
         ratio_for_tracer0 = acute_zeta/ (acute_zeta+ acute_zeta0)

         do kk=1,size(zion,1)
            tracer_at_ZeroPhi(kk) = ratio_for_tracer*tracer(kk) + ratio_for_tracer0*tracer0(kk)
         enddo

         ! tracer(6) = sqrt(magnetic moment) is constant during simulation except collisional process.
         ! Therefore, tracer(6) is always enough to deal with magnetic moment mu = tracer(6)*tracer(6).
         ! In this sense, tracer0(6) is not used for sqrt(mu), but it is set by value 1 (in loading routine) for other purpose.
         ! We need to recover tracer_at_ZeroPhi(6) using tracer(6) only, not including tracer0(6) which is not sqrt(mu).
         tracer_at_ZeroPhi(6) = tracer(6)

         ! new angle (theta and zeta) should be calculated using two angles whose difference is smaller than pi
         ! Angle of tracer0 is modified to satisfy the above condition.
         ! A = tracer angle, B = tracer0;   B' = B + 2*n*pi where n = floor((A-B+pi)/(2*pi));
         do kk=2,3
            tracer_at_ZeroPhi(kk) = ratio_for_tracer*tracer(kk) + ratio_for_tracer0*  ( tracer0(kk) + 2.0*pi* floor( (tracer(kk) - tracer0(kk) + pi )/(2.0*pi)  ))
            tracer_at_ZeroPhi(kk) = tracer_at_ZeroPhi(kk) - floor(tracer_at_ZeroPhi(kk)/(2.0*pi)) * 2.0*pi ! normalization
         end do

         aaa(1)=tracer_at_ZeroPhi(1)
         the(1)=tracer_at_ZeroPhi(2)
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
         tracer_energy=0.5_wp*(tracer_at_ZeroPhi(4)*Bj(1))**2+tracer_at_ZeroPhi(6)**2*Bj(1)  ! Ion energy
         ! Save tracers_Poincare Information to send_buf_TRI_Poincare
         num_of_MY_tracers_Poincare = num_of_MY_tracers_Poincare + 1

         send_buf_TRI_Poincare(1,num_of_MY_tracers_Poincare) = istep
         send_buf_TRI_Poincare(2,num_of_MY_tracers_Poincare) = rj(1)
         send_buf_TRI_Poincare(3,num_of_MY_tracers_Poincare) = zj(1)
         send_buf_TRI_Poincare(4:10,num_of_MY_tracers_Poincare) = tracer_at_ZeroPhi(1:7)
         send_buf_TRI_Poincare(11,num_of_MY_tracers_Poincare) = tracer_energy
         send_buf_TRI_Poincare(12,num_of_MY_tracers_Poincare) = ilabels(my_tracer_indices(m))
         send_buf_TRI_Poincare(13,num_of_MY_tracers_Poincare) = mype
         send_buf_TRI_Poincare(14,num_of_MY_tracers_Poincare) = wi_init(my_tracer_indices(m))
         send_buf_TRI_Poincare(15,num_of_MY_tracers_Poincare) = wi_vpara(my_tracer_indices(m))
         send_buf_TRI_Poincare(16,num_of_MY_tracers_Poincare) = wi_vperp(my_tracer_indices(m))
         send_buf_TRI_Poincare(17,num_of_MY_tracers_Poincare) = wi_accel(my_tracer_indices(m))
         send_buf_TRI_Poincare(18,num_of_MY_tracers_Poincare) = wi_fdamp(my_tracer_indices(m))
         send_buf_TRI_Poincare(19,num_of_MY_tracers_Poincare) = UnNorm_theta_i(my_tracer_indices(m))
         send_buf_TRI_Poincare(20,num_of_MY_tracers_Poincare) = UnNorm_zeta_i(my_tracer_indices(m))
      endif
   enddo

   call MPI_GATHER(num_of_MY_tracers_Poincare, 1, MPI_INTEGER, recv_TRI_Poincare_counts, 1, MPI_INTEGER, 0, toroidal_comm, ierror)

   if(mype == 0) then
      recv_TRI_Poincare_displs(1) = 0
      recv_num_of_TRI_Poincare = recv_TRI_Poincare_counts(1)
      do i=2, nproc_toroidal
         recv_TRI_Poincare_displs(i) = recv_TRI_Poincare_displs(i-1) + recv_TRI_Poincare_counts(i-1)
         recv_num_of_TRI_Poincare = recv_num_of_TRI_Poincare + recv_TRI_Poincare_counts(i)
      enddo
   endif

   call MPI_Gatherv(send_buf_TRI_Poincare, num_of_MY_tracers_Poincare, MPI_TRI_Type, recv_buf_TRI_Poincare, recv_TRI_Poincare_counts, recv_TRI_Poincare_displs, MPI_TRI_Type, 0, toroidal_comm, ierror)

   if (mype == 0) then
      open(unit=1516, file="multi_tracers_ion_Poincare.dat", position='append',action='write')
      do i=1, recv_num_of_TRI_Poincare
         !write(1516,*) recv_buf_TRI_Poincare(:,i)
         write(1516,102) recv_buf_TRI_Poincare(:,i)
      enddo
      call flush(1516)
      close(1516)
   endif

   deallocate( send_buf_TRI, recv_buf_TRI, send_buf_TRI_Poincare, recv_buf_TRI_Poincare)

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


   each_PTL_length = 19 ! including we_adiab + we_fdamp


   MY_PE_weight_sum=0
   MY_PE_abnormal_weight_sum=0
   num_of_Send_PTL=0
   do m=1,me
      weight = zeon(5,m)
      MY_PE_weight_sum = MY_PE_weight_sum +weight
      if(weight<min_thold_weight .or. weight>max_thold_weight) then
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
      send_buf_PTL(12,i) = elabels(m)
      send_buf_PTL(13,i) = mype
      send_buf_PTL(14,i) = we_init(m)
      send_buf_PTL(15,i) = we_vpara(m)
      send_buf_PTL(16,i) = we_vperp(m)
      send_buf_PTL(17,i) = we_accel(m)
      send_buf_PTL(18,i) = we_adiab(m)
      send_buf_PTL(19,i) = we_fdamp(m)
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

   real(wp) weight, energy
   character(len=30) file_name


   each_PTL_length = 18 !

   MY_PE_weight_sum=0
   MY_PE_abnormal_weight_sum=0
   num_of_Send_PTL=0
   do m=1,mi
      weight = zion(5,m)
      MY_PE_weight_sum = MY_PE_weight_sum +weight
      if(weight<min_thold_weight .or. weight>max_thold_weight) then
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
      send_buf_PTL(12,i) = ilabels(m)
      send_buf_PTL(13,i) = mype
      send_buf_PTL(14,i) = wi_init(m)
      send_buf_PTL(15,i) = wi_vpara(m)
      send_buf_PTL(16,i) = wi_vperp(m)
      send_buf_PTL(17,i) = wi_accel(m)
      send_buf_PTL(18,i) = wi_fdamp(m)
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
