subroutine shifti
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use esi_interface

   implicit none

   integer i,k,tid, m, mtop,lasth, dim_zion, dim_zion0, dim_RK_cum_zion, dim_Ptl_package, ierr,istatus(MPI_STATUS_SIZE)
   integer send_L_count,send_R_count, send_num, send_mid_arr(mimax), left_mid_arr(mimax), right_mid_arr(mimax)
   integer recv_L_count,recv_R_count, recv_num
   integer, dimension(:),allocatable :: fill_mid_arr

   integer m0, sid, rid, idest, isource, loop,tor_dir, tot_send_num

   real(wp) normalized_zeta
   real(wp),dimension(:,:),allocatable :: send_Ptls, recv_Ptls

   integer MPI_Ptl_type

   !!!! $acc declare create (send_mid_arr,left_mid_arr,right_mid_arr)

   if(acc_newer_D_zion) then
      !$acc update self (zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
      acc_newer_D_zion=.false.
      acc_newer_H_zion=.false.
   endif


   dim_zion = size(zion,1)
   dim_zion0 = size(zion0,1)
   if(Flag_RK4==1) then
      dim_RK_cum_zion = size(RK_cum_zion,1)
   else
      dim_RK_cum_zion =0
   endif

   dim_Ptl_package = dim_zion + dim_zion0 + dim_RK_cum_zion + 8  ! zeon + zeon0 + RK_cum_zion+ ilables + 5 values about wi evolution +2 Unnormalized Angles


   call MPI_Type_contiguous( dim_Ptl_package, MPI_Rsize, MPI_Ptl_type, ierr);
   call MPI_Type_commit(MPI_Ptl_type, ierr);

   m0=1
   loop=0
   do while(1)
      loop = loop +1
      if(loop>ntoroidal)then
         write(0,*)'endless ION shift loop at PE=',mype
         call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      endif
      ! find send_counts and related information
      ! send_counts = 0
      send_L_count=0
      send_R_count=0
      send_num = 0

      write(stdout, '(5(a,I8))') "istep=",istep,"loop=",loop,"mype=",mype, "m0=",m0,"mi=",mi
      call flush(stdout)


      !$acc parallel loop seq reduction(+:send_num,send_L_count,send_R_count)
      do m=m0,mi
         normalized_zeta =  modulo(zion(3,m),2.0*pi)

         if( (normalized_zeta-zetamax)*(normalized_zeta-zetamin) >0 ) then
            send_num=send_num+1
            send_mid_arr(send_num) = m

            if( modulo(zion(3,m)-zetamax, 2.0*pi)<pi) then
               ! send to right rank
               send_R_count = send_R_count +1
               right_mid_arr(send_R_count) = m
            else
               ! send to left rank
               send_L_count = send_L_count +1
               left_mid_arr(send_L_count) = m
            endif
         endif
      enddo


      !!!!$acc update self (send_mid_arr,left_mid_arr,right_mid_arr)
      call quicksort_integer(send_mid_arr,1,send_num)
      call quicksort_integer(right_mid_arr,1,send_R_count)
      call quicksort_integer(left_mid_arr,1,send_L_count)
      !!!!$acc update device (send_mid_arr,left_mid_arr,right_mid_arr)

      write(stdout, '(a,I8,a,2I8)') "mype=",mype,"; SHIFTI!! send_mid_arr(1,send_num)=",send_mid_arr(1),send_mid_arr(send_num)
      call flush(stdout)


      write(stdout, '(2(a,I8),a,2I8)') "mype=",mype,"; SHIFTI!! send_num=",send_num,", send_counts=",send_L_count, send_R_count
      write(stdout, '(a,I8,a,6I7))') "mype=",mype,"; send(1:2), right(1:2), left(1:2)=",send_mid_arr(1:2),right_mid_arr(1:2),left_mid_arr(1:2)
      call flush(stdout)



      call MPI_ALLREDUCE(send_num,tot_send_num,1,MPI_INTEGER,MPI_SUM,toroidal_comm,ierr)
      if(tot_send_num==0) then
         if(loop>1) then
            !acc_newer_H_zion=.true.
         endif
         write(stdout, '(5(a,I8))') "SHIFTI Finished!! istep=",istep,"loop=",loop,"mype=",mype, "m0=",m0,"mi=",mi
         call flush(stdout)

         !$acc update self(zion0(:,:mi),RK_cum_zion(:,:mi))
         return
      endif

      call MPI_SENDRECV(send_R_count,1,MPI_INTEGER,right_pe,0, recv_L_count,1,MPI_INTEGER,left_pe ,0,toroidal_comm,istatus,ierr)
      call MPI_SENDRECV(send_L_count,1,MPI_INTEGER,left_pe ,0, recv_R_count,1,MPI_INTEGER,right_pe,0,toroidal_comm,istatus,ierr)

      recv_num=recv_L_count+recv_R_count

      write(stdout, '(2(a,I8),a,2I8)') "mype=",mype,"; SHIFTI!! recv_num=",recv_num,", recv_counts=",recv_L_count,recv_R_count
      call flush(stdout)

      allocate(send_Ptls(dim_Ptl_package, send_num), recv_Ptls(dim_Ptl_package, recv_num))

      ! copy all information of sending particles into tor_send_Ptls

      !$acc parallel loop seq copyin(left_mid_arr(:send_L_count),right_mid_arr(:send_R_count))
      do i=1, send_num

         if (i<=send_L_count) then
            m=left_mid_arr(i)
         else
            m=right_mid_arr(i-send_L_count)
         endif

         k=0
         send_Ptls(k+1:k+dim_zion, i) =  zion(:,m)
         k=k+dim_zion
         send_Ptls(k+1:k+dim_zion0, i) =  zion0(:,m)
         k=k+dim_zion0
         if(Flag_RK4==1) then
            send_Ptls(k+1:k+dim_RK_cum_zion, i) =  RK_cum_zion(:,m)
            k=k+dim_RK_cum_zion
         endif
         send_Ptls(k+1, i) =  ilabels(m)
         send_Ptls(k+2, i) =  wi_init(m)
         send_Ptls(k+3, i) =  wi_vpara(m)
         send_Ptls(k+4, i) =  wi_vperp(m)
         send_Ptls(k+5, i) =  wi_accel(m)
         send_Ptls(k+6, i) =  wi_fdamp(m)
         send_Ptls(k+7, i) =  UnNorm_theta_i(m)
         send_Ptls(k+8, i) =  UnNorm_zeta_i(m)
      enddo


      ! Send to Left
      call MPI_SENDRECV(send_Ptls(:,1:send_L_count),send_L_count,MPI_Ptl_type,left_pe,0, recv_Ptls(:,1+recv_L_count:recv_num),recv_R_count,MPI_Ptl_type,right_pe ,0,toroidal_comm,istatus,ierr)
      ! Send to Right
      call MPI_SENDRECV(send_Ptls(:,1+send_L_count:send_num),send_R_count,MPI_Ptl_type,right_pe,0, recv_Ptls(:,1:recv_L_count),recv_L_count,MPI_Ptl_type,left_pe ,0,toroidal_comm,istatus,ierr)

      write(stdout, '(a,I8,a,5F8.4))') "mype=",mype,"; send_Ptls_right(1:5,1)=",send_Ptls(1:5,1)
      write(stdout, '(a,I8,a,5F8.4))') "mype=",mype,";  recv_Ptls_left(1:5,1)=",recv_Ptls(1:5,1+recv_L_count)
      call flush(stdout)
      write(stdout, '(a,I8,a,5F8.4))') "mype=",mype,"; send_Ptls_left(1:5,end)=",send_Ptls(1:5,send_num)
      write(stdout, '(a,I8,a,5F8.4))') "mype=",mype,"; recv_Ptls_right(1:5,end)=",recv_Ptls(1:5,recv_L_count)
      call flush(stdout)




      !!!$acc parallel loop seq copy(mtop,lasth) copyin(mi) reduction(+:lasth,mtop) private(send_mid_arr)
      !!$acc kernels copy(mtop,lasth,mi)




      ! fill holes
      mtop = mi
      mi = mi - send_num ! remained mi after sending
      lasth = send_num

      allocate(fill_mid_arr(send_num))
      do i=1, send_num
         m = send_mid_arr(i)
         if(m>mi) exit
         do while(mtop == send_mid_arr(lasth))
            mtop=mtop-1
            lasth=lasth-1
         enddo
         fill_mid_arr(i) = mtop
         mtop = mtop -1
      enddo


      write(stdout, '((a,I8),2(a,2(I8,1x)))') "mype=",mype,"; SHIFTI!! send_mid_arr=",send_mid_arr(1),send_mid_arr(send_num),"fill_mid_arr=",fill_mid_arr(1), fill_mid_arr(send_num)
      write(stdout, '(4(a,I8))') "mype=",mype,"; SHIFTI!! mtop=",mtop,", mi=",mi,"  lasth=",lasth
      call flush(stdout)

      !$acc parallel loop seq
      do i=1, send_num
         m = send_mid_arr(i)
         if(m>mi) exit
         mtop= fill_mid_arr(i)

         zion(:, m)      = zion(:, mtop)
         zion0(:, m)     = zion0(:, mtop)
         ilabels(m)      = ilabels(mtop)
         wi_init(m)      = wi_init(mtop)
         wi_vpara(m)     = wi_vpara(mtop)
         wi_vperp(m)     = wi_vperp(mtop)
         wi_accel(m)     = wi_accel(mtop)
         wi_fdamp(m)     = wi_fdamp(mtop)

         UnNorm_theta_i(m) =UnNorm_theta_i(mtop)
         UnNorm_zeta_i(m) =UnNorm_zeta_i(mtop)
         if(Flag_RK4==1) RK_cum_zion(:,m)= RK_cum_zion(:,mtop)
      enddo

      ! write(stdout, '(4(a,I8),a,2I8)') "mype=",mype,"; SHIFTI!! mtop=",mtop,", mi=",mi,"  lasth=",lasth,  "; send_mid_arr(1,lasth)=",send_mid_arr(1), send_mid_arr(lasth)
      ! call flush(stdout)
      ! append recv particles
      !$acc parallel loop seq
      do i=1, recv_num
         k=0
         zion(:,mi+i)    = recv_Ptls(k+1:k+dim_zion, i)
         k=k+dim_zion
         zion0(:,mi+i)   = recv_Ptls(k+1:k+dim_zion0, i)
         k=k+dim_zion0
         if(Flag_RK4==1) then
            RK_cum_zion(:,mi+i) = recv_Ptls(k+1:k+dim_RK_cum_zion, i)
            k=k+dim_RK_cum_zion
         endif
         ilabels(mi+i)   =  recv_Ptls(k+1, i)
         wi_init(mi+i)   =  recv_Ptls(k+2, i)
         wi_vpara(mi+i)  =  recv_Ptls(k+3, i)
         wi_vperp(mi+i)  =  recv_Ptls(k+4, i)
         wi_accel(mi+i)  =  recv_Ptls(k+5, i)
         wi_fdamp(mi+i)  =  recv_Ptls(k+6, i)
         UnNorm_theta_i(mi+i)= recv_Ptls(k+7, i)
         UnNorm_zeta_i(mi+i) = recv_Ptls(k+8, i)
      enddo
      mi = mi + recv_num

      if(allocated(send_Ptls)) deallocate(send_Ptls)
      if(allocated(recv_Ptls)) deallocate(recv_Ptls)
      if(allocated(fill_mid_arr)) deallocate(fill_mid_arr)


      m0=mi-recv_num+1  ! new starting index of next loop
   enddo

end subroutine shifti
