subroutine shifti
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use esi_interface

  implicit none
    
  integer i,tid, m, mtop,lasth, dest_tor_rank, dim_zion, dim_zion0, dim_Ptl_package, ierr
  integer tor_send_counts(ntoroidal), tor_send_displs(ntoroidal), tot_send_num
  integer tor_recv_counts(ntoroidal), tor_recv_displs(ntoroidal), tot_recv_num
  integer send_mid_arr(mimax), send_dest_arr(mimax)

  real(wp) normalized_zeta, Delta_zeta
  real(wp),dimension(:,:),allocatable :: tor_send_Ptls_data, tor_recv_Ptls_data

  integer MPI_Ptl_type
  
  character(len=30) file_name



  ! write(file_name,'("shifti_out_",i2.2,".dat")')mype
  ! open(unit=4321, file=file_name, position='append', action='write')

  
  
  ! write(4321,*) "shift_i start! I'm=",mype,"mytorrank",myrank_toroidal,"ntoroidal",ntoroidal, "toroidal_comm", toroidal_comm
  ! write(4321,*) "zetamin", zetamin, "zetamax", zetamax
  ! call flush(4321)
 


  dim_zion = size(zion,1)
  dim_zion0 = size(zion0,1)

  !dim_Ptl_package = dim_zion + dim_zion0 + 1  ! zeon + zeon0 + ilables
  dim_Ptl_package = dim_zion + dim_zion0 + 5  ! zeon + zeon0 + ilables + 4 values about wi evolution
  ! write(4321,*) "dim_zion=", dim_zion, "dim_zion0=",dim_zion0
  ! call flush(4321)

  
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
  do m=1, mi 
    normalized_zeta =  zion(3,m) - floor(zion(3,m)/(2.0*pi)) * 2.0*pi
    dest_tor_rank = floor(normalized_zeta / Delta_zeta)
    
    if(dest_tor_rank<0 .OR. dest_tor_rank>ntoroidal) then
      
      write(file_name,'("shifti_out_",i2.2,".dat")')mype
      open(unit=4321, file=file_name, position='append', action='write')
  
      write(4321,*) "Error!! dest_tor_rank=", dest_tor_rank, 'zeta',zion(3,m), 'norm_zeta', normalized_zeta,'zeta0', zion0(3,m)
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
    tor_send_Ptls_data(1:dim_zion, i) =  zion(:,send_mid_arr(i))
    tor_send_Ptls_data(dim_zion+1:dim_zion+dim_zion0, i) =  zion0(:,send_mid_arr(i))
    tor_send_Ptls_data(dim_zion+dim_zion0+1, i) =  ilabels(send_mid_arr(i))  
    tor_send_Ptls_data(dim_zion+dim_zion0+2, i) =  wi_init(send_mid_arr(i))  
    tor_send_Ptls_data(dim_zion+dim_zion0+3, i) =  wi_vpara(send_mid_arr(i))
    tor_send_Ptls_data(dim_zion+dim_zion0+4, i) =  wi_vperp(send_mid_arr(i))  
    tor_send_Ptls_data(dim_zion+dim_zion0+5, i) =  wi_accel(send_mid_arr(i))    
  enddo
  

  ! AlltoAllv
  call MPI_Alltoallv(tor_send_Ptls_data, tor_send_counts, tor_send_displs, MPI_Ptl_type, tor_recv_Ptls_data, tor_recv_counts, tor_recv_displs, MPI_Ptl_type, toroidal_comm, ierr)




  
  ! fill holes 
  call quicksort(send_mid_arr, 1, tot_send_num)
  mtop = mi
  mi = mi - tot_send_num ! remained mi after sending
  lasth = tot_send_num
  do i=1, tot_send_num
    m = send_mid_arr(i)

    if(m>mi) exit
    do while(mtop == send_mid_arr(lasth))
        mtop=mtop-1
        lasth=lasth-1
    enddo

    zion(:, m) = zion(:, mtop)
    zion0(:, m) = zion0(:, mtop)
    ilabels(m) = ilabels(mtop)
    wi_init(m)     = wi_init(mtop)
    wi_vpara(m)   = wi_vpara(mtop)
    wi_vperp(m)  = wi_vperp(mtop)
    wi_accel(m)= wi_accel(mtop)

    mtop = mtop -1
  enddo


  ! append recv particles
  ! zion(:, mi+1:mi+tot_recv_num) = tor_recv_Ptls_data(1:dim_zion, :)
  ! zion0(:, mi+1:mi+tot_recv_num) = tor_recv_Ptls_data(dim_zion+1:dim_zion+dim_zion0, :)
  ! ilabels(mi+1:mi+tot_recv_num) =  tor_recv_Ptls_data(dim_zion+dim_zion0+1, :)   
  ! wi_init(mi+1:mi+tot_recv_num) =  tor_recv_Ptls_data(dim_zion+dim_zion0+2, :) 
  ! wi_vpara(mi+1:mi+tot_recv_num) =  tor_recv_Ptls_data(dim_zion+dim_zion0+3, :) 
  ! wi_vperp(mi+1:mi+tot_recv_num) =  tor_recv_Ptls_data(dim_zion+dim_zion0+4, :) 
  ! wi_accel(mi+1:mi+tot_recv_num) =  tor_recv_Ptls_data(dim_zion+dim_zion0+5, :) 
 
  ! append recv particles
  do i=1, tot_recv_num
    zion(:,mi+i) = tor_recv_Ptls_data(1:dim_zion, i)
    zion0(:,mi+i) = tor_recv_Ptls_data(dim_zion+1:dim_zion+dim_zion0, i)
    ilabels(mi+i) =  tor_recv_Ptls_data(dim_zion+dim_zion0+1, i)   
    wi_init(mi+i) =  tor_recv_Ptls_data(dim_zion+dim_zion0+2, i) 
    wi_vpara(mi+i) =  tor_recv_Ptls_data(dim_zion+dim_zion0+3, i) 
    wi_vperp(mi+i) =  tor_recv_Ptls_data(dim_zion+dim_zion0+4, i) 
    wi_accel(mi+i) =  tor_recv_Ptls_data(dim_zion+dim_zion0+5, i) 
  enddo

  mi = mi + tot_recv_num  

  deallocate(tor_send_Ptls_data, tor_recv_Ptls_data)
  
  
  !close(4321)


end subroutine shifti


