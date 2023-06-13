subroutine shifti
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use esi_interface

   implicit none

   integer i,j,k,tid, m, mtop,lasth, dim_zion, dim_zion0, dim_RK_cum_zion, dim_Ptl_package, ierr,istatus(MPI_STATUS_SIZE)
   integer send_L_count,send_R_count, send_num, send_mid_arr(mimax), left_mid_arr(mimax), right_mid_arr(mimax)
   integer recv_L_count,recv_R_count, recv_num
   integer, dimension(:),allocatable :: fill_mid_arr
   integer fill_num
   real(wp) zion_zeta_arr(mimax)

   integer m0, sid, rid, idest, isource, loop,tor_dir, tot_send_num

   real(wp) normalized_zeta
   real(wp),dimension(:,:),allocatable :: send_Ptls, recv_Ptls

   integer request_send_L, request_send_R, request_recv_L, request_recv_R
   integer MPI_Ptl_type

   !$acc declare create (zion_zeta_arr)

   if(acc_newer_H_zion) then
      !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
      acc_newer_H_zion=.false.
      acc_newer_D_zion=.false.
   endif

   dim_zion = size(zion,1)
   dim_zion0 = size(zion0,1)
   dim_RK_cum_zion = size(RK_cum_zion,1)


   dim_Ptl_package = dim_zion + dim_zion0 + dim_RK_cum_zion + 7  ! zeon + zeon0 + RK_cum_zion+ ilables + 4 values about wi evolution +2 Unnormalized Angles


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


      !$acc parallel loop
      do m=m0,mi
         zion_zeta_arr(m)=zion(3,m)
      enddo
      !$acc update self(zion_zeta_arr(m0:mi))

      do m=m0,mi
         normalized_zeta =  modulo(zion_zeta_arr(m),2.0*pi)

         if( (normalized_zeta-zetamax)*(normalized_zeta-zetamin) >0 ) then
            send_num=send_num+1
            send_mid_arr(send_num) = m

            if( modulo(zion_zeta_arr(m)-zetamax, 2.0*pi)<pi) then
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


      call MPI_ALLREDUCE(send_num,tot_send_num,1,MPI_INTEGER,MPI_SUM,toroidal_comm,ierr)
      if(tot_send_num==0) then
         if(loop>1) then
            acc_newer_D_zion=.true.
            acc_newer_H_zion=.false.

            ntracer=0
            !$acc parallel loop present(ilabels) copy(ntracer) reduction(max:ntracer)
            do m=1,mi
               if(ilabels(m)==ilabel_of_ntracer) then
                  ntracer=m
               endif
            enddo
         endif

#ifdef _MGY_DEBUG_
         if(mype==0) then
            write(stdout,'(3(A,I8))') "In shifti, istep=",istep," irk=",irk, ", loop=",loop
            call flush(stdout)
         endif
#endif

         return
      endif

      call MPI_SENDRECV(send_R_count,1,MPI_INTEGER,right_pe,0, recv_L_count,1,MPI_INTEGER,left_pe ,0,toroidal_comm,istatus,ierr)
      call MPI_SENDRECV(send_L_count,1,MPI_INTEGER,left_pe ,0, recv_R_count,1,MPI_INTEGER,right_pe,0,toroidal_comm,istatus,ierr)

      recv_num=recv_L_count+recv_R_count

      allocate(send_Ptls(dim_Ptl_package, send_num), recv_Ptls(dim_Ptl_package, recv_num))

      ! copy all information of sending particles into tor_send_Ptls

      !$acc parallel loop copyin(left_mid_arr(:send_L_count),right_mid_arr(:send_R_count))
      do i=1, send_num

         if (i<=send_L_count) then
            m=left_mid_arr(i)
         else
            m=right_mid_arr(i-send_L_count)
         endif

         k=0
         !$acc loop seq
         do j=1,dim_zion
            send_Ptls(k+j, i) =  zion(j,m)
         enddo
         k=k+dim_zion

         !$acc loop seq
         do j=1,dim_zion0
            send_Ptls(k+j, i) =  zion0(j,m)
         enddo
         k=k+dim_zion0

         !$acc loop seq
         do j=1,dim_RK_cum_zion
            send_Ptls(k+j, i) =  RK_cum_zion(j,m)
         enddo
         k=k+dim_RK_cum_zion
         send_Ptls(k+1, i) =  ilabels(m)
         send_Ptls(k+2, i) =  wi_init(m)
         send_Ptls(k+3, i) =  wi_vpara(m)
         send_Ptls(k+4, i) =  wi_vperp(m)
         send_Ptls(k+5, i) =  wi_accel(m)
         send_Ptls(k+6, i) =  UnNorm_theta_i(m)
         send_Ptls(k+7, i) =  UnNorm_zeta_i(m)
      enddo


      ! Send to Left and Right 
      if(send_L_count > 0) call MPI_ISEND(send_Ptls(:,1:send_L_count), send_L_count, MPI_Ptl_type, left_pe, 0, toroidal_comm, request_send_L, ierr)
      if(send_R_count > 0) call MPI_ISEND(send_Ptls(:,1+send_L_count:send_num), send_R_count, MPI_Ptl_type, right_pe, 0, toroidal_comm, request_send_R, ierr)

      ! Recv from Left and Right
      if(recv_L_count > 0) call MPI_IRECV(recv_Ptls(:,1:recv_L_count), recv_L_count, MPI_Ptl_type, left_pe, 0, toroidal_comm, request_recv_L, ierr)
      if(recv_R_count > 0) call MPI_IRECV(recv_Ptls(:,1+recv_L_count:recv_num), recv_R_count, MPI_Ptl_type, right_pe, 0, toroidal_comm, request_recv_R, ierr)
     
      ! Wait all non-blocking communications are finished
      if(send_L_count > 0) call MPI_WAIT(request_send_L, istatus, ierr)
      if(send_R_count > 0) call MPI_WAIT(request_send_R, istatus, ierr)
      if(recv_L_count > 0) call MPI_WAIT(request_recv_L, istatus, ierr)
      if(recv_R_count > 0) call MPI_WAIT(request_recv_R, istatus, ierr)



      ! fill holes
      mtop = mi
      mi = mi - send_num ! remained mi after sending
      lasth = send_num

      allocate(fill_mid_arr(send_num))
      fill_num=0
      do i=1, send_num
         m = send_mid_arr(i)
         if(m>mi) exit
         do while(mtop == send_mid_arr(lasth))
            mtop=mtop-1
            lasth=lasth-1
         enddo
         fill_mid_arr(i) = mtop
         mtop = mtop -1
         fill_num=i
      enddo


      !$acc parallel loop
      do i=1, fill_num
         m = send_mid_arr(i)
         mtop= fill_mid_arr(i)
         !$acc loop seq
         do j=1,dim_zion
            zion(j, m)      = zion(j, mtop)
         enddo
         !$acc loop seq
         do j=1,dim_zion0
            zion0(j, m)     = zion0(j, mtop)
         enddo
         ilabels(m)      = ilabels(mtop)
         wi_init(m)      = wi_init(mtop)
         wi_vpara(m)     = wi_vpara(mtop)
         wi_vperp(m)     = wi_vperp(mtop)
         wi_accel(m)     = wi_accel(mtop)

         UnNorm_theta_i(m) =UnNorm_theta_i(mtop)
         UnNorm_zeta_i(m) =UnNorm_zeta_i(mtop)
         !$acc loop seq
         do j=1,dim_RK_cum_zion
            RK_cum_zion(j,m)= RK_cum_zion(j,mtop)
         enddo
      enddo


      ! append recv particles
      !$acc parallel loop
      do i=1, recv_num
         k=0
         !$acc loop seq
         do j=1,dim_zion
            zion(j,mi+i)    = recv_Ptls(k+j, i)
         enddo
         k=k+dim_zion
         !$acc loop seq
         do j=1,dim_zion0
            zion0(j,mi+i)   = recv_Ptls(k+j, i)
         enddo
         k=k+dim_zion0
         !$acc loop seq
         do j=1,dim_RK_cum_zion
            RK_cum_zion(j,mi+i) = recv_Ptls(k+j, i)
         enddo
         k=k+dim_RK_cum_zion
         ilabels(mi+i)   =  recv_Ptls(k+1, i)
         wi_init(mi+i)   =  recv_Ptls(k+2, i)
         wi_vpara(mi+i)  =  recv_Ptls(k+3, i)
         wi_vperp(mi+i)  =  recv_Ptls(k+4, i)
         wi_accel(mi+i)  =  recv_Ptls(k+5, i)
         UnNorm_theta_i(mi+i)= recv_Ptls(k+6, i)
         UnNorm_zeta_i(mi+i) = recv_Ptls(k+7, i)
      enddo
      mi = mi + recv_num

      if(allocated(send_Ptls)) deallocate(send_Ptls)
      if(allocated(recv_Ptls)) deallocate(recv_Ptls)
      if(allocated(fill_mid_arr)) deallocate(fill_mid_arr)


      m0=mi-recv_num+1  ! new starting index of next loop
   enddo

end subroutine shifti
