subroutine shifte
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use esi_interface

   implicit none

   integer i,k,tid, m, mtop,lasth, dim_zeon, dim_zeon0, dim_RK_cum_zeon, dim_Ptl_package, ierr,istatus(MPI_STATUS_SIZE)
   integer send_counts(0:1), send_num, send_mid_arr(memax), left_mid_arr(memax), right_mid_arr(memax)
   integer recv_counts(0:1), recv_num


   integer m0, sid, rid, idest, isource, loop,tor_dir, tot_send_num

   real(wp) normalized_zeta
   real(wp),dimension(:,:),allocatable :: tor_send_Ptls_data, tor_recv_Ptls_data

   integer MPI_Ptl_type



   dim_zeon = size(zeon,1)
   dim_zeon0 = size(zeon0,1)
   if(Flag_RK4==1) then
      dim_RK_cum_zeon = size(RK_cum_zeon,1)
   else
      dim_RK_cum_zeon =0
   endif

   dim_Ptl_package = dim_zeon + dim_zeon0 + dim_RK_cum_zeon + 9  ! zeon + zeon0 + RK_cum_zeon+ elabels + 6 values about we evolution +2 Unnormalized Angles


   call MPI_Type_contiguous( dim_Ptl_package, MPI_Rsize, MPI_Ptl_type, ierr);
   call MPI_Type_commit(MPI_Ptl_type, ierr);

   m0=1
   loop=0
   do while(1)
      loop = loop +1
      if(loop>ntoroidal)then
         write(0,*)'endless EON shift loop at PE=',mype
         call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      endif
      ! find send_counts and related information
      send_counts = 0
      send_num = 0
      do m=m0,me
         normalized_zeta =  modulo(zeon(3,m),2.0*pi)

         if( (normalized_zeta-zetamax)*(normalized_zeta-zetamin) >0 ) then
            send_num=send_num+1
            send_mid_arr(send_num) = m

            if( modulo(zeon(3,m)-zetamax, 2.0*pi)<pi) then
               ! send to right rank
               send_counts(1) = send_counts(1) +1
               right_mid_arr(send_counts(1)) = m
            else
               ! send to left rank
               send_counts(0) = send_counts(0) +1
               left_mid_arr(send_counts(0)) = m
            endif
         endif
      enddo


      call MPI_ALLREDUCE(send_num,tot_send_num,1,MPI_INTEGER,MPI_SUM,toroidal_comm,ierr)
      if(tot_send_num==0) then
         return
      endif

      write(stdout,*) "I'm=",mype,"SHIFTE loop=",loop,m0, "send_counts", send_counts
      call flush(stdout)

      call MPI_SENDRECV(send_counts(1),1,MPI_INTEGER,right_pe,0, recv_counts(0),1,MPI_INTEGER,left_pe ,0,toroidal_comm,istatus,ierr)
      call MPI_SENDRECV(send_counts(0),1,MPI_INTEGER,left_pe ,0, recv_counts(1),1,MPI_INTEGER,right_pe,0,toroidal_comm,istatus,ierr)

      recv_num=recv_counts(0)+recv_counts(1)


      allocate(tor_send_Ptls_data(dim_Ptl_package, send_num), tor_recv_Ptls_data(dim_Ptl_package, recv_num))

      ! copy all information of sending particles into tor_send_Ptls
      do i=1, send_num

         if (i<=send_counts(0)) then
            m=left_mid_arr(i)
         else
            m=right_mid_arr(i-send_counts(0))
         endif

         k=0
         tor_send_Ptls_data(k+1:k+dim_zeon, i) =  zeon(:,m)
         k=k+dim_zeon
         tor_send_Ptls_data(k+1:k+dim_zeon0, i) =  zeon0(:,m)
         k=k+dim_zeon0
         if(Flag_RK4==1) then
            tor_send_Ptls_data(k+1:k+dim_RK_cum_zeon, i) =  RK_cum_zeon(:,m)
            k=k+dim_RK_cum_zeon
         endif
         tor_send_Ptls_data(k+1, i) =  elabels(m)
         tor_send_Ptls_data(k+2, i) =  we_init(m)
         tor_send_Ptls_data(k+3, i) =  we_vpara(m)
         tor_send_Ptls_data(k+4, i) =  we_vperp(m)
         tor_send_Ptls_data(k+5, i) =  we_accel(m)
         tor_send_Ptls_data(k+6, i) =  we_adiab(m)
         tor_send_Ptls_data(k+7, i) =  we_fdamp(m)
         tor_send_Ptls_data(k+8, i) =  UnNorm_theta_e(m)
         tor_send_Ptls_data(k+9, i) =  UnNorm_zeta_e(m)
      enddo

      ! Send to Left
      call MPI_SENDRECV(tor_send_Ptls_data(:,1:send_counts(0)),send_counts(0),MPI_Ptl_type,left_pe,0, tor_recv_Ptls_data(:,1+recv_counts(0):recv_num),recv_counts(1),MPI_Ptl_type,right_pe ,0,toroidal_comm,istatus,ierr)
      ! Send to Right
      call MPI_SENDRECV(tor_send_Ptls_data(:,1+send_counts(0):send_num),send_counts(1),MPI_Ptl_type,right_pe,0, tor_recv_Ptls_data(:,1:recv_counts(0)),recv_counts(0),MPI_Ptl_type,left_pe ,0,toroidal_comm,istatus,ierr)


      ! fill holes
      mtop = me
      me = me - send_num ! remained me after sending
      lasth = send_num
      do i=1, send_num
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
         if(Flag_RK4==1) RK_cum_zeon(:,m)= RK_cum_zeon(:,mtop)

         mtop = mtop -1
      enddo

      ! append recv particles
      do i=1, recv_num
         k=0
         zeon(:,me+i)    = tor_recv_Ptls_data(k+1:k+dim_zeon, i)
         k=k+dim_zeon
         zeon0(:,me+i)   = tor_recv_Ptls_data(k+1:k+dim_zeon0, i)
         k=k+dim_zeon0
         if(Flag_RK4==1) then
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
      me = me + recv_num

      if(allocated(tor_send_Ptls_data)) deallocate(tor_send_Ptls_data)
      if(allocated(tor_recv_Ptls_data)) deallocate(tor_recv_Ptls_data)



      m0=me-recv_num+1  ! new starting index of next loop
   enddo

end subroutine shifte
