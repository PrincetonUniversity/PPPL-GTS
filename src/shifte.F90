subroutine shifte
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use esi_interface

   implicit none

   integer i,k,tid, m, mtop,lasth, dest_tor_rank, dim_zeon, dim_zeon0, dim_RK_cum_zeon, dim_Ptl_package,ierr
   integer tor_send_counts(ntoroidal), tor_send_displs(ntoroidal), tot_send_num
   integer tor_recv_counts(ntoroidal), tor_recv_displs(ntoroidal), tot_recv_num
   integer send_mid_arr(memax), send_dest_arr(memax)

   real(wp) normalized_zeta, Delta_zeta
   real(wp),dimension(:,:),allocatable :: tor_send_Ptls_data, tor_recv_Ptls_data

   integer MPI_Ptl_type

   character(len=30) file_name

   ! write(file_name,'("shifte_out_",i2.2,".dat")')mype
   ! open(unit=4321, file=file_name, position='append', action='write')

   ! do m=1,me
   !   if(isnan(zeon(3,m)) .or. isnan(RK_cum_zeon(3,m))) then
   !     write(stdout,*) 'Before Shifte, Im=',mype,', m=',m,', zeon=', zeon(:,m)
   !     write(stdout,*) 'Before Shifte, Im=',mype,', m=',m,', zeon0=', zeon0(:,m)
   !     write(stdout,*) 'Before Shifte, Im=',mype,', m=',m,', RK_cum_zeon=', RK_cum_zeon(:,m)
   !     call flush(stdout)
   !   endif
   ! enddo

   ! write(4321,*) "shift_e start! I'm=",mype,"mytorrank",myrank_toroidal,"ntoroidal",ntoroidal, "toroidal_comm", toroidal_comm
   ! write(4321,*) "zetamin", zetamin, "zetamax", zetamax
   ! call flush(4321)



   dim_zeon = size(zeon,1)
   dim_zeon0 = size(zeon0,1)
   if(RK4_Flag==1) then
      dim_RK_cum_zeon = size(RK_cum_zeon,1)
   else
      dim_RK_cum_zeon = 0
   endif
   !dim_Ptl_package = dim_zeon + dim_zeon0 + 1  ! zeon + zeon0 + elabels
   ! write(4321,*) "dim_zeon=", dim_zeon, "dim_zeon0=",dim_zeon0
   ! call flush(4321)

   dim_Ptl_package = dim_zeon + dim_zeon0 + dim_RK_cum_zeon + 9  ! zeon + zeon0 +cum_RK4 +elabels + 7 values about we evolution +2 Unnormalized angles


   call MPI_Type_contiguous( dim_Ptl_package, MPI_Rsize, MPI_Ptl_type, ierr);
   call MPI_Type_commit(MPI_Ptl_type, ierr);


   Delta_zeta = 2.0*pi/ntoroidal


   ! write(4321,*) "Delta_zeta =", Delta_zeta
   ! call flush(4321)



   ! find send_counts and related information
   tor_send_counts = 0
   tor_send_displs = 0
   tot_send_num = 0
   send_mid_arr = 0
   send_dest_arr = 0
   do m=1, me
      ! normalized_zeta =  zeon(3,m) - floor(zeon(3,m)/(2.0*pi)) * 2.0*pi
      normalized_zeta =  modulo(zeon(3,m),2.0*pi)
      dest_tor_rank = floor(normalized_zeta / Delta_zeta)

      if(dest_tor_rank<0 .OR. dest_tor_rank>ntoroidal) then
         write(file_name,'("shifte_out_",i2.2,".dat")')mype
         open(unit=4321, file=file_name, position='append', action='write')
         write(4321,*) "Error!! dest_tor_rank=", dest_tor_rank, 'zeta',zeon(3,m), 'norm_zeta', normalized_zeta,'zeta0', zeon0(3,m)
         call flush(4321)
         close(4321)
      endif

      if(dest_tor_rank /= myrank_toroidal) then
         tor_send_counts(1+dest_tor_rank) = tor_send_counts(1+dest_tor_rank) + 1
         tot_send_num = tot_send_num + 1
         send_mid_arr(tot_send_num) = m
         send_dest_arr(tot_send_num) = dest_tor_rank
      endif
   enddo



   do tid=2, ntoroidal
      tor_send_displs(tid) = tor_send_displs(tid-1) + tor_send_counts(tid-1)
   enddo


   ! establish recv_counts using send_counts
   tor_recv_counts = 0
   tor_recv_displs = 0
   tot_recv_num = 0
   call MPI_Alltoall(tor_send_counts, 1, MPI_INTEGER, tor_recv_counts, 1, MPI_INTEGER, toroidal_comm, ierr);
   do tid=1, ntoroidal
      tot_recv_num = tot_recv_num + tor_recv_counts(tid)
   enddo
   do tid=2, ntoroidal
      tor_recv_displs(tid) = tor_recv_displs(tid-1) + tor_recv_counts(tid-1)
   enddo


   ! write(4321,*) mype, myrank_toroidal, "tot_send_num =", tot_send_num, "tot_recv_num=", tot_recv_num
   ! write(4321,*) mype, myrank_toroidal, "tor_send_counts =", tor_send_counts
   ! write(4321,*) mype, myrank_toroidal, "tor_send_displs =", tor_send_displs
   ! write(4321,*) mype, myrank_toroidal, "tor_recv_counts =", tor_recv_counts
   ! write(4321,*) mype, myrank_toroidal, "tor_recv_displs =", tor_recv_displs
   ! call flush(4321)


   !allocate required memories
   allocate(tor_send_Ptls_data(dim_Ptl_package, tot_send_num))
   allocate(tor_recv_Ptls_data(dim_Ptl_package, tot_recv_num))


   ! write(4321,*) "mype=",mype,", myrank_toroidal=",myrank_toroidal, "tot_send_num",tot_send_num, "tot_recv_num", tot_recv_num
   ! call flush(4321)



   ! 여기서 뿌리기전에 target_tor_rank에 맞게 send_mid_arr 을 sorting 미리 해줘야한다.

   ! sort send_dest_arr ascending order; send_mid_arr is also sorted by the ascending order of send_dest_arr
   call quicksort2(send_dest_arr, send_mid_arr, 1, tot_send_num)

   ! copy all information of sending particles into tor_send_Ptls
   do i=1, tot_send_num
      k=0
      tor_send_Ptls_data(k+1:k+dim_zeon, i) =  zeon(:,send_mid_arr(i))
      k=k+dim_zeon
      tor_send_Ptls_data(k+1:k+dim_zeon0, i) =  zeon0(:,send_mid_arr(i))
      k=k+dim_zeon0
      if(RK4_Flag==1) then
         tor_send_Ptls_data(k+1:k+dim_RK_cum_zeon, i) =  RK_cum_zeon(:,send_mid_arr(i))
         k=k+dim_RK_cum_zeon
      endif
      tor_send_Ptls_data(k+1, i) =  elabels(send_mid_arr(i))
      tor_send_Ptls_data(k+2, i) =  we_init(send_mid_arr(i))
      tor_send_Ptls_data(k+3, i) =  we_vpara(send_mid_arr(i))
      tor_send_Ptls_data(k+4, i) =  we_vperp(send_mid_arr(i))
      tor_send_Ptls_data(k+5, i) =  we_accel(send_mid_arr(i))
      tor_send_Ptls_data(k+6, i) =  we_adiab(send_mid_arr(i))
      tor_send_Ptls_data(k+7, i) =  we_fdamp(send_mid_arr(i))
      tor_send_Ptls_data(k+8, i) =  UnNorm_theta_e(send_mid_arr(i))
      tor_send_Ptls_data(k+9, i) =  UnNorm_zeta_e(send_mid_arr(i))
   enddo


   ! AlltoAllv
   call MPI_Alltoallv(tor_send_Ptls_data, tor_send_counts, tor_send_displs, MPI_Ptl_type, tor_recv_Ptls_data, tor_recv_counts, tor_recv_displs, MPI_Ptl_type, toroidal_comm, ierr)





   ! fill holes
   call quicksort(send_mid_arr, 1, tot_send_num)
   mtop = me
   me = me - tot_send_num ! remained me after sending
   lasth = tot_send_num
   do i=1, tot_send_num
      m = send_mid_arr(i)

      if(m>me) exit
      do while(mtop == send_mid_arr(lasth))
         mtop=mtop-1
         lasth=lasth-1
      enddo

      zeon(:, m)      = zeon(:, mtop)
      zeon0(:, m)     = zeon0(:, mtop)
      elabels(m)      = elabels(mtop)
      we_init(m)      = we_init(mtop)
      we_vpara(m)     = we_vpara(mtop)
      we_vperp(m)     = we_vperp(mtop)
      we_accel(m)     = we_accel(mtop)
      we_adiab(m)     = we_adiab(mtop)
      we_fdamp(m)     = we_fdamp(mtop)

      UnNorm_theta_e(m) =UnNorm_theta_e(mtop)
      UnNorm_zeta_e(m) =UnNorm_zeta_e(mtop)

      if(RK4_Flag==1) RK_cum_zeon(:,m)= RK_cum_zeon(:,mtop)

      mtop = mtop -1
   enddo


   ! ! append recv particles
   ! zeon(:, me+1:me+tot_recv_num) = tor_recv_Ptls_data(1:dim_zeon, :)
   ! zeon0(:, me+1:me+tot_recv_num) = tor_recv_Ptls_data(dim_zeon+1:dim_zeon+dim_zeon0, :)
   ! elabels(me+1:me+tot_recv_num) =  tor_recv_Ptls_data(dim_zeon+dim_zeon0+1, :)

   ! append recv particles
   do i=1, tot_recv_num
      k=0
      zeon(:,me+i)    = tor_recv_Ptls_data(k+1:k+dim_zeon, i)
      k=k+dim_zeon
      zeon0(:,me+i)   = tor_recv_Ptls_data(k+1:k+dim_zeon0, i)
      k=k+dim_zeon0
      if(RK4_Flag==1) then
         RK_cum_zeon(:,me+i) = tor_recv_Ptls_data(k+1:k+dim_RK_cum_zeon, i)
         k=k+dim_RK_cum_zeon
      endif
      elabels(me+i)   =  tor_recv_Ptls_data(k+1, i)
      we_init(me+i)   =  tor_recv_Ptls_data(k+2, i)
      we_vpara(me+i)  =  tor_recv_Ptls_data(k+3, i)
      we_vperp(me+i)  =  tor_recv_Ptls_data(k+4, i)
      we_accel(me+i)  =  tor_recv_Ptls_data(k+5, i)
      we_adiab(me+i)  =  tor_recv_Ptls_data(k+6, i)
      we_fdamp(me+i)  =  tor_recv_Ptls_data(k+7, i)
      UnNorm_theta_e(me+i)= tor_recv_Ptls_data(k+8, i)
      UnNorm_zeta_e(me+i) = tor_recv_Ptls_data(k+9, i)
   enddo

   me = me + tot_recv_num

   deallocate(tor_send_Ptls_data, tor_recv_Ptls_data)


   ! close(4321)

   ! do m=1,me
   !   if(isnan(zeon(3,m)) .or. isnan(RK_cum_zeon(3,m))) then
   !     write(stdout,*) 'After Shifte, Im=',mype,', m=',m,', zeon=', zeon(:,m)
   !     write(stdout,*) 'After Shifte, Im=',mype,', m=',m,', zeon0=', zeon0(:,m)
   !     write(stdout,*) 'After Shifte, Im=',mype,', m=',m,', RK_cum_zeon=', RK_cum_zeon(:,m)
   !     call flush(stdout)
   !   endif
   ! enddo


end subroutine shifte






recursive subroutine quicksort(a, first, last)
   ! quicksort.f -*-f90-*-
   ! Author: t-nissie
   ! License: GPLv3
   ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
   implicit none
   integer  a(*), x, t
   integer first, last
   integer i, j

   x = a( (first+last) / 2 )
   i = first
   j = last
   do
      do while (a(i) < x)
         i=i+1
      end do
      do while (x < a(j))
         j=j-1
      end do
      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t

      i=i+1
      j=j-1
   end do
   if (first < i-1) call quicksort(a,first, i-1)
   if (j+1 < last)  call quicksort(a,j+1, last)
end subroutine quicksort


recursive subroutine quicksort2(a, b, first, last)
   ! quicksort.f -*-f90-*-
   ! Author: t-nissie
   ! License: GPLv3
   ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
   implicit none
   integer  a(*),b(*), x, t
   integer first, last
   integer i, j

   x = a( (first+last) / 2 )
   i = first
   j = last
   do
      do while (a(i) < x)
         i=i+1
      end do
      do while (x < a(j))
         j=j-1
      end do
      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      t = b(i);  b(i) = b(j);  b(j) = t
      i=i+1
      j=j-1
   end do
   if (first < i-1) call quicksort2(a, b,first, i-1)
   if (j+1 < last)  call quicksort2(a, b,j+1, last)
end subroutine quicksort2
