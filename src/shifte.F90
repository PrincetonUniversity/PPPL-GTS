subroutine shifte
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use esi_interface

   implicit none

   integer i,j,k,tid, m, mtop,lasth, dim_zeon, dim_zeon0, dim_RK_cum_zeon, dim_Ptl_package, ierr,istatus(MPI_STATUS_SIZE)
   integer send_L_count,send_R_count, send_num, send_mid_arr(memax), left_mid_arr(memax), right_mid_arr(memax)
   integer recv_L_count,recv_R_count, recv_num
   integer, dimension(:),allocatable :: fill_mid_arr
   integer fill_num
   real(wp) zeon_zeta_arr(memax)

   integer m0, sid, rid, idest, isource, loop,tor_dir, tot_send_num

   real(wp) normalized_zeta
   real(wp),dimension(:,:),allocatable :: send_Ptls, recv_Ptls

   integer request_send_L, request_send_R, request_recv_L, request_recv_R
   integer MPI_Ptl_type

   !$acc declare create (zeon_zeta_arr)


   if(acc_newer_H_zeon) then
      !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
      acc_newer_H_zeon=.false.
      acc_newer_D_zeon=.false.
   endif

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
      ! if(loop>ntoroidal)then
      !    !$acc update self(zeon)
      !    write(0,*) 'mype=',mype,'; zeon=',zeon()

      !    write(0,*)'endless EON shift loop at PE=',mype
      !    call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      ! endif
      ! find send_counts and related information
      send_L_count=0
      send_R_count=0
      send_num = 0

      !$acc parallel loop
      do m=m0,me
         zeon_zeta_arr(m)=zeon(3,m)
      enddo
      !$acc update self(zeon_zeta_arr(m0:me))

      do m=m0,me
         normalized_zeta =  modulo(zeon_zeta_arr(m),2.0*pi)

         if( (normalized_zeta-zetamax)*(normalized_zeta-zetamin) >0 ) then

            if(loop>ntoroidal) then
               write(0,'(2(A,I4),A,I8, 2(A,ES12.4))') 'LOOP=',loop, ',mype=',mype,'; m=',m,', zeon_zeta_arr(m)=',zeon_zeta_arr(m),'; normal_zeta=',normalized_zeta
            endif

            send_num=send_num+1
            send_mid_arr(send_num) = m

            if( modulo(zeon_zeta_arr(m)-zetamax, 2.0*pi)<pi) then
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
      if(loop>ntoroidal)then
         !$acc update self(zeon)
         write(0,'(2(A,I4),A,2I8,A,2ES12.4,2(A,I8))') 'LOOP=',loop, ',mype=',mype,'; m0,me=',m0,me,'; zetamin,max=',zetamin,zetamax,'; send_L_count=', send_L_count, '; send_R_count=',send_R_count
         write(0,'(2(A,I4),3(A,I8))') 'LOOP=',loop,',mype=',mype,'; S1=',send_mid_arr(1),'; R1=',right_mid_arr(1),'; L1=',left_mid_arr(1)
         write(0,'(2(A,I4),A,ES12.4,A,10ES12.4)') 'LOOP=',loop,',mype=',mype,'; zeon_zeta(L1)=',zeon_zeta_arr(left_mid_arr(1)),'; zeon(L1)=',zeon(:,left_mid_arr(1))
         write(0,'(2(A,I4),A,ES12.4,A,10ES12.4)') 'LOOP=',loop,',mype=',mype,'; zeon_zeta(R1)=',zeon_zeta_arr(right_mid_arr(1)),'; zeon(R1)=',zeon(:,right_mid_arr(1))
         write(0,'(2(A,I4))')'endless EON shift loop at PE=',mype ,'; LOOP=',loop
         if(loop>2*ntoroidal) then
            call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
         endif
      endif

      if(tot_send_num==0) then
         if(loop>1) then
            acc_newer_D_zeon=.true.
            acc_newer_H_zeon=.false.

            ntracer_e=0
            !$acc parallel loop present(elabels) copy(ntracer_e) reduction(max:ntracer_e)
            do m=1,me
               if(elabels(m)==elabel_of_ntracer_e) then
                  ntracer_e=m
               endif
            enddo
         endif

#ifdef _MGY_DEBUG_
         if(mype==0) then
            write(stdout,'(3(A,I8))') "In shifte, istep=",istep," irk=",irk, ", loop=",loop
            call flush(stdout)
         endif
#endif

         return
      endif

      call MPI_SENDRECV(send_R_count,1,MPI_INTEGER,right_pe,0, recv_L_count,1,MPI_INTEGER,left_pe ,0,toroidal_comm,istatus,ierr)
      call MPI_SENDRECV(send_L_count,1,MPI_INTEGER,left_pe ,0, recv_R_count,1,MPI_INTEGER,right_pe,0,toroidal_comm,istatus,ierr)

      recv_num=recv_L_count+recv_R_count


      allocate(send_Ptls(dim_Ptl_package, send_num), recv_Ptls(dim_Ptl_package, recv_num))
      send_Ptls=0.0
      recv_Ptls=0.0

      !$acc enter data create(send_Ptls, recv_Ptls)


      ! copy all information of sending particles into tor_send_Ptls
      !$acc parallel loop present(send_Ptls) copyin(left_mid_arr(:send_L_count),right_mid_arr(:send_R_count))
      do i=1, send_num

         if (i<=send_L_count) then
            m=left_mid_arr(i)
         else
            m=right_mid_arr(i-send_L_count)
         endif

         k=0
         !$acc loop seq
         do j=1,dim_zeon
            send_Ptls(k+j, i) =  zeon(j,m)
         enddo
         k=k+dim_zeon
         !$acc loop seq
         do j=1,dim_zeon0
            send_Ptls(k+j, i) =  zeon0(j,m)
         enddo
         k=k+dim_zeon0
         if(Flag_RK4==1) then
            !$acc loop seq
            do j=1,dim_RK_cum_zeon
               send_Ptls(k+j, i) =  RK_cum_zeon(j,m)
            enddo
            k=k+dim_RK_cum_zeon
         endif
         send_Ptls(k+1, i) =  elabels(m)
         send_Ptls(k+2, i) =  we_init(m)
         send_Ptls(k+3, i) =  we_vpara(m)
         send_Ptls(k+4, i) =  we_vperp(m)
         send_Ptls(k+5, i) =  we_accel(m)
         send_Ptls(k+6, i) =  we_adiab(m)
         send_Ptls(k+7, i) =  we_fdamp(m)
         send_Ptls(k+8, i) =  UnNorm_theta_e(m)
         send_Ptls(k+9, i) =  UnNorm_zeta_e(m)
      enddo

      if(send_num>0) then
         !$acc update self(send_Ptls(:,:send_num))
      endif

#ifdef _MGY_DEBUG_
      if(loop>ntoroidal)then
         write(0,'(2(A,I4),2(A,3I8))') 'LOOP=',loop, ',mype=',mype,'; send_num,L,R=', send_num, send_L_count, send_R_count, '; recv_num,L,R=', recv_num, recv_L_count, recv_R_count
      endif
#endif

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


      if(recv_num>0) then
         !$acc update device(recv_Ptls(:,:recv_num))
      endif

#ifdef _MGY_DEBUG_
      if(loop>ntoroidal)then
         if(send_num>0) then
            write(0,'(2(A,I4),A,5ES12.4)') 'LOOP=',loop, ',mype=',mype,'; send_Ptls(3,1:5)=', send_Ptls(3,1:5)
         endif
         if(recv_num>0) then
            write(0,'(2(A,I4),A,5ES12.4)') 'LOOP=',loop, ',mype=',mype,'; recv_Ptls(3,1:5)=', recv_Ptls(3,1:5)
         endif
      endif
#endif

      ! fill holes
      mtop = me
      me = me - send_num ! remained me after sending
      lasth = send_num

      allocate(fill_mid_arr(send_num))
      fill_num=0
      do i=1, send_num
         m = send_mid_arr(i)
         if(m>me) exit
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
         do j=1,dim_zeon
            zeon(j, m)   = zeon(j, mtop)
         enddo
         !$acc loop seq
         do j=1,dim_zeon0
            zeon0(j, m)  = zeon0(j, mtop)
         enddo
         elabels(m)      = elabels(mtop)
         we_init(m)      = we_init(mtop)
         we_vpara(m)     = we_vpara(mtop)
         we_vperp(m)     = we_vperp(mtop)
         we_accel(m)     = we_accel(mtop)
         we_adiab(m)     = we_adiab(mtop)
         we_fdamp(m)     = we_fdamp(mtop)

         UnNorm_theta_e(m) =UnNorm_theta_e(mtop)
         UnNorm_zeta_e(m) =UnNorm_zeta_e(mtop)
         if(Flag_RK4==1) then
            !$acc loop seq
            do j=1,dim_RK_cum_zeon
               RK_cum_zeon(j,m)= RK_cum_zeon(j,mtop)
            enddo
         endif
      enddo

#ifdef _MGY_DEBUG_
      if(loop>ntoroidal)then
         if(recv_num>0) then
            write(0,'(2(A,I4),A,5ES12.4)') 'AA LOOP=',loop, ',mype=',mype,'; recv_Ptls(3,1:5)=', recv_Ptls(3,1:5)
         endif
      endif
#endif

      ! append recv particles
      !$acc parallel loop present(recv_Ptls)
      do i=1, recv_num
         k=0
         !$acc loop seq
         do j=1,dim_zeon
            zeon(j,me+i)    = recv_Ptls(k+j, i)
         enddo
         k=k+dim_zeon
         !$acc loop seq
         do j=1,dim_zeon0
            zeon0(j,me+i)   = recv_Ptls(k+j, i)
         enddo
         k=k+dim_zeon0
         if(Flag_RK4==1) then
            !$acc loop seq
            do j=1,dim_RK_cum_zeon
               RK_cum_zeon(j,me+i) = recv_Ptls(k+j, i)
            enddo
            k=k+dim_RK_cum_zeon
         endif
         elabels(me+i)   =  recv_Ptls(k+1, i)
         we_init(me+i)   =  recv_Ptls(k+2, i)
         we_vpara(me+i)  =  recv_Ptls(k+3, i)
         we_vperp(me+i)  =  recv_Ptls(k+4, i)
         we_accel(me+i)  =  recv_Ptls(k+5, i)
         we_adiab(me+i)  =  recv_Ptls(k+6, i)
         we_fdamp(me+i)  =  recv_Ptls(k+7, i)
         UnNorm_theta_e(me+i)= recv_Ptls(k+8, i)
         UnNorm_zeta_e(me+i) = recv_Ptls(k+9, i)
      enddo
      me = me + recv_num
      m0=me-recv_num+1  ! new starting index of next loop


#ifdef _MGY_DEBUG_
      if(loop>ntoroidal)then
         if(recv_num>0) then
            write(0,'(3(A,I4),A,5ES12.4)') 'BB LOOP=',loop, ',mype=',mype,', recv_num=',recv_num,'; recv_Ptls(3,1:5)=', recv_Ptls(3,1:5)
         endif

         !$acc update self(zeon)
         write(0,'(2(A,I4),A,2I8,A,12ES12.4)') 'BB LOOP=',loop, ',mype=',mype,'; m0,me=',m0,me,'; zeon(3,m0:m0+5)=',zeon(3,m0:m0+5)
      endif
#endif

      !$acc exit data delete(send_Ptls,recv_Ptls)

      if(allocated(send_Ptls)) deallocate(send_Ptls)
      if(allocated(recv_Ptls)) deallocate(recv_Ptls)
      if(allocated(fill_mid_arr)) deallocate(fill_mid_arr)


      ! m0=me-recv_num+1  ! new starting index of next loop
   enddo

end subroutine shifte
