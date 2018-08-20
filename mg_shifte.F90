subroutine shifte
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use esi_interface

  implicit none
    
  integer i,tid, m, mtop,lasth, dest_tor_rank, dim_zeon, dim_zeon0, ierr
  integer tor_send_counts(ntoroidal), tor_send_displs(ntoroidal), tot_send_num
  integer tor_recv_counts(ntoroidal), tor_recv_displs(ntoroidal), tot_recv_num
  integer send_mid_arr(memax), send_dest_arr(memax)

  real(wp) normalized_zeta, Delta_zeta
  real(wp),dimension(:,:),allocatable :: tor_send_Ptls_data, tor_recv_Ptls_data

  integer MPI_Ptl_type
  
  character(len=30) file_name

  ! write(file_name,'("shifte_out_",i2.2,".dat")')mype
  ! open(unit=4321, file=file_name, position='append', action='write')

  
  
  ! write(4321,*) "shift_e start! I'm=",mype,"mytorrank",myrank_toroidal,"ntoroidal",ntoroidal, "toroidal_comm", toroidal_comm
  ! write(4321,*) "zetamin", zetamin, "zetamax", zetamax
  ! call flush(4321)
 


  dim_zeon = size(zeon,1)
  dim_zeon0 = size(zeon0,1)

  ! write(4321,*) "dim_zeon=", dim_zeon, "dim_zeon0=",dim_zeon0
  ! call flush(4321)

  
  call MPI_Type_contiguous( (dim_zeon + dim_zeon0 + 1), MPI_Rsize, MPI_Ptl_type, ierr);
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
    normalized_zeta =  zeon(3,m) - floor(zeon(3,m)/(2.0*pi)) * 2.0*pi
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
  allocate(tor_send_Ptls_data(dim_zeon + dim_zeon0 + 1, tot_send_num))
  allocate(tor_recv_Ptls_data(dim_zeon + dim_zeon0 + 1, tot_recv_num))


  ! write(4321,*) "mype=",mype,", myrank_toroidal=",myrank_toroidal, "tot_send_num",tot_send_num, "tot_recv_num", tot_recv_num
  ! call flush(4321)
  


  ! 여기서 뿌리기전에 target_tor_rank에 맞게 send_mid_arr 을 sorting 미리 해줘야한다.
  
  ! sort send_dest_arr ascending order; send_mid_arr is also sorted by the ascending order of send_dest_arr
   call quicksort2(send_dest_arr, send_mid_arr, 1, tot_send_num)

  ! copy all information of sending particles into tor_send_Ptls
  do i=1, tot_send_num
    tor_send_Ptls_data(1:dim_zeon, i) =  zeon(:,send_mid_arr(i))
    tor_send_Ptls_data(dim_zeon+1:dim_zeon+dim_zeon0, i) =  zeon0(:,send_mid_arr(i))
    tor_send_Ptls_data(dim_zeon+dim_zeon0+1, i) =  elabels(send_mid_arr(i))
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

    zeon(:, m) = zeon(:, mtop)
    zeon0(:, m) = zeon0(:, mtop)
    elabels(m) = elabels(mtop)
    mtop = mtop -1
  enddo


  ! append recv particles
  zeon(:, me+1:me+tot_recv_num) = tor_recv_Ptls_data(1:dim_zeon, :)
  zeon0(:, me+1:me+tot_recv_num) = tor_recv_Ptls_data(dim_zeon+1:dim_zeon+dim_zeon0, :)  
  elabels(me+1:me+tot_recv_num) =  tor_recv_Ptls_data(dim_zeon+dim_zeon0+1, :)  
  
  me = me + tot_recv_num  

  deallocate(tor_send_Ptls_data, tor_recv_Ptls_data)
  
  
  ! close(4321)


end subroutine shifte






recursive subroutine quicksort(a, first, last)
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