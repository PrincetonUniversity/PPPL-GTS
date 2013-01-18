subroutine shifte
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use esi_interface

  implicit none
  
  integer i,m,msendleft(2),msendright(2),mrecvleft(2),mrecvright(2),mtop,m0,msend,mrecv, &
       idest,isource,isendtag,irecvtag,nzeon,iright(memax),ileft(memax),kzeon(memax), &
       isendcount,irecvcount,istatus(MPI_STATUS_SIZE),ierror,iteration,lasth
  real(wp),dimension(:,:),allocatable :: recvleft,recvright,sendleft,sendright
  real(wp) zetaright,zetaleft,pi_inv
  character(len=8) cdum

#ifdef _OPENMP
  integer msleft(32,0:15),msright(32,0:15)
  integer nthreads,gnthreads,iam,delm,mbeg,mend,omp_get_num_threads,&
       omp_get_thread_num
#endif

  nzeon=14    ! 7 zeon + 7 zeon0
  pi_inv=1.0_wp/pi
  m0=1
  iteration=0

100 iteration=iteration+1
  if(iteration>numberpe)then
     write(0,*)'endless ELECTRON sorting loop at PE=',mype, ' mrecv=', mrecv
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  msend=0
  msendright=0
  msendleft=0

  if(m0 <= me)then

!  !$omp parallel do private(m)
!        do m=1,me
!           kzeon(m)=0
!        enddo

  kzeon=0

#ifdef _OPENMP
! This section of code (down to #else) is included by the preprocessor when
! the compilation is performed with OpenMP support. We must then use a few
! temporary arrays and add some work distribution code.
! First, we initialize the shared (and temporary) arrays msleft and msright
! to zero.
  msleft=0
  msright=0

! Then we start the parallel region with !$omp parallel

!$omp parallel private(nthreads,iam,delm,i,mbeg,mend,m,zetaright,zetaleft) &
!$omp& shared(gnthreads)
  nthreads=omp_get_num_threads()    !Get the number of threads ready to work
  iam=omp_get_thread_num()          !Get my thread number (position)
  delm=(me-m0+1)/nthreads       !Calculate the number of steps per thread
  i=mod((me-m0+1),nthreads)
!$omp single              !Put nthread in global variable for later use.
  gnthreads=nthreads      !nthread is the same for all threads so only one
!$omp end single nowait   !of them needs to copy the value in gnthreads

! We now distribute the work between the threads. The loop over the particles
! is distributed equally (as much as possible) between them.
  mbeg=m0+min(iam,i)*(delm+1)+max(0,(iam-i))*delm
  mend=mbeg+delm+(min((iam+1),i)/(iam+1))-1

! label particle to be moved
     do m=mbeg,mend
        zetaright=min(2.0_wp*pi,zeon(3,m))-zetamax
        zetaleft=zeon(3,m)-zetamin

        if( zetaright*zetaleft > 0 )then
           zetaright=zetaright*0.5_wp*pi_inv
           zetaright=zetaright-real(floor(zetaright),wp)
           msright(3,iam)=msright(3,iam)+1
           kzeon(mbeg+msright(3,iam)-1)=m

           if( zetaright < 0.5_wp )then
! particle to move right
              msright(1,iam)=msright(1,iam)+1
              iright(mbeg+msright(1,iam)-1)=m
! keep track of tracer
              if( m == ntracer_e )then
                 msright(2,iam)=msright(1,iam)
                 ntracer_e=0
              endif

! particle to move left
           else
              msleft(1,iam)=msleft(1,iam)+1
              ileft(mbeg+msleft(1,iam)-1)=m
              if( m == ntracer_e )then
                 msleft(2,iam)=msleft(1,iam)
                 ntracer_e=0
              endif
           endif
        endif
     enddo
! End of the OpenMP parallel region
!$omp end parallel

! Now that we are out of the parallel region we need to gather and rearrange
! the results of the multi-thread calculation. We need to end up with the
! same arrays as for the sequential (single-threaded) calculation.
     do m=0,gnthreads-1
        delm=(me-m0+1)/gnthreads
        i=mod((me-m0+1),gnthreads)
        mbeg=m0+min(m,i)*(delm+1)+max(0,(m-i))*delm
        if( msleft(2,m) /= 0 )msendleft(2)=msendleft(1)+msleft(2,m)
        do i=1,msleft(1,m)
           ileft(msendleft(1)+i)=ileft(mbeg+i-1)
        enddo
        msendleft(1)=msendleft(1)+msleft(1,m)
        if( msright(2,m) /= 0 )msendright(2)=msendright(1)+msright(2,m)
        do i=1,msright(1,m)
           iright(msendright(1)+i)=iright(mbeg+i-1)
        enddo
        msendright(1)=msendright(1)+msright(1,m)
        do i=1,msright(3,m)
           kzeon(msend+i)=kzeon(mbeg+i-1)
        enddo
        msend=msend+msright(3,m)
     enddo

#else
!  This section of code replaces the section above when the compilation does
!  NOT include the OpenMP support option. Temporary arrays msleft and msright
!  are not needed as well as the extra code for thread work distribution.

     do m=m0,me
        zetaright=min(2.0_wp*pi,zeon(3,m))-zetamax
        zetaleft=zeon(3,m)-zetamin

        if( zetaright*zetaleft > 0 )then
           zetaright=zetaright*0.5_wp*pi_inv
           zetaright=zetaright-real(floor(zetaright),wp)
           msend=msend+1
           kzeon(msend)=m

           if( zetaright < 0.5_wp )then
! # of particle to move right
              msendright(1)=msendright(1)+1
              iright(msendright(1))=m
! keep track of tracer
              if( m == ntracer_e )then
                 msendright(2)=msendright(1)
                 ntracer_e=0
              endif

! # of particle to move left
           else
              msendleft(1)=msendleft(1)+1
              ileft(msendleft(1))=m
              if( m == ntracer_e )then
                 msendleft(2)=msendleft(1)
                 ntracer_e=0
              endif
           endif
        endif
     enddo

#endif

  endif

  if (msend /= (msendleft(1)+msendright(1))) then
     write(0,*)'mype=',mype,'  msend NOT equal to msendleft+msendright'
     msend=msendleft(1)+msendright(1)
  endif

  if(iteration>1)then
! total # of particles to be shifted
     mrecv=0

     call MPI_ALLREDUCE(msend,mrecv,1,MPI_INTEGER,MPI_SUM,toroidal_comm,ierror)

! no particle to be shifted, return
     if ( mrecv == 0 ) then
!        write(*,*)istep,irk,mype,me,m0,iteration

        return
     endif
  endif

! an extra space to prevent zero size when msendright(1)=msendleft(1)=0
  allocate(sendright(nzeon,max(1,msendright(1))),sendleft(nzeon,max(1,msendleft(1))))

! pack particle to move right
!$omp parallel do private(m)
     do m=1,msendright(1)
        sendright(1:nzeon/2,m)=zeon(1:nzeon/2,iright(m))
        sendright((nzeon/2+1):nzeon,m)=zeon0(1:nzeon/2,iright(m))
     enddo

! pack particle to move left
!$omp parallel do private(m)
     do m=1,msendleft(1)    
        sendleft(1:nzeon/2,m)=zeon(1:nzeon/2,ileft(m))
        sendleft((nzeon/2+1):nzeon,m)=zeon0(1:nzeon/2,ileft(m))
     enddo

     mtop=me

! # of particles remain on local PE
     me=me-msendleft(1)-msendright(1)
! fill the hole
     lasth=msend
     do i=1,msend
        m=kzeon(i)
        if(m>me) exit  !Break out of the DO loop
        do while(mtop == kzeon(lasth))
           mtop=mtop-1
           lasth=lasth-1
        enddo

        if(mtop==ntracer_e) ntracer_e=m

        zeon(1:nzeon/2,m)=zeon(1:nzeon/2,mtop)
        zeon0(1:nzeon/2,m)=zeon0(1:nzeon/2,mtop)
        mtop=mtop-1
!        if (mtop <= me) exit  !Break out of the DO loop
     enddo

! send # of particle to move right
  mrecvleft=0
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(msendright,2,MPI_INTEGER,idest,isendtag,&
       mrecvleft,2,MPI_INTEGER,isource,irecvtag,toroidal_comm,istatus,ierror)

  allocate(recvleft(nzeon,max(1,mrecvleft(1))))
 
! send particle to right and receive from left
  recvleft=0.0_wp
  isendcount=max(1,msendright(1))*nzeon
  irecvcount=max(1,mrecvleft(1))*nzeon
  call MPI_SENDRECV(sendright,isendcount,mpi_Rsize,idest,isendtag,recvleft,&
       irecvcount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
  
! send # of particle to move left
  mrecvright=0
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(msendleft,2,MPI_INTEGER,idest,isendtag,&
       mrecvright,2,MPI_INTEGER,isource,irecvtag,toroidal_comm,istatus,ierror)
  
  allocate(recvright(nzeon,max(1,mrecvright(1))))

! send particle to left and receive from right
  recvright=0.0_wp
  isendcount=max(1,msendleft(1))*nzeon
  irecvcount=max(1,mrecvright(1))*nzeon
  call MPI_SENDRECV(sendleft,isendcount,mpi_Rsize,idest,isendtag,recvright,&
       irecvcount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
  
! tracer particle
  if( mrecvleft(2) > 0 )then
     ntracer_e=me+mrecvleft(2)
  elseif( mrecvright(2) > 0 )then
     ntracer_e=me+mrecvleft(1)+mrecvright(2)
  endif

! need extra particle array
  if(me+mrecvleft(1)+mrecvright(1) > memax)then
     write(0,*)"need bigger particle array",mype,memax,mrecvleft(1),mrecvright(1),me+mrecvleft(1)+mrecvright(1)
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
! open disk file
!     if(mype < 10)then
!        write(cdum,'("TEMP.00",i1)')mype
!     elseif(mype < 100)then
!        write(cdum,'("TEMP.0",i2)')mype
!     else
!        write(cdum,'("TEMP.",i3)')mype
!     endif
!     open(111,file=cdum,status='replace',form='unformatted')
     
! record particle information to file
!     write(111)zeon(1:nzeon/2,1:me)
!     write(111)zeon0(1:nzeon/2,1:me)
     
! make bigger array
!     deallocate(zeon,zeon0)
!     memax=2*(me+mrecvleft(1)+mrecvright(1))-memax
!     allocate(zeon(nzeon/2,memax),zeon0(nzeon/2,memax))
     
! read back particle information
!     rewind(111)
!     read(111)zeon(1:nzeon/2,1:me)
!     read(111)zeon0(1:nzeon/2,1:me)
!     close(111)
  endif

! unpack particle, particle moved from left
!$omp parallel do private(m)
  do m=1,mrecvleft(1)
     zeon(1:nzeon/2,m+me)=recvleft(1:nzeon/2,m)
     zeon0(1:nzeon/2,m+me)=recvleft((nzeon/2+1):nzeon,m)
  enddo

! particle moved from right
!$omp parallel do private(m)
  do m=1,mrecvright(1)
     zeon(1:nzeon/2,m+me+mrecvleft(1))=recvright(1:nzeon/2,m)
     zeon0(1:nzeon/2,m+me+mrecvleft(1))=recvright((nzeon/2+1):nzeon,m)
  enddo

  me=me+mrecvleft(1)+mrecvright(1)  ! current partcle number= original - sent + received
  
  deallocate(sendleft,sendright,recvleft,recvright)
  m0=me-mrecvright(1)-mrecvleft(1)+1  ! starting index of received particles
  goto 100
  
end subroutine shifte


