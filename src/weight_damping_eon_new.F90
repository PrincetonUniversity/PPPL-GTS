subroutine weight_field_eon(wfield)
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use esi_interface

  implicit none

  integer,parameter ::     lwork=100

  integer m,i,j,k,jt,ii,ierror,ij,info,  &
       icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
  real(wp) rdum(9),wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,vz,mu,enew,f_rho 
  real(wp) send(3,mgrid),recv(3,mgrid)

  real(wp) dnitmp(9,mgrid,0:1), dampf(9,mgrid,0:1), wfield(3,mgrid,0:1)
  real(wp) aav(3,3),bv(3),work(lwork)
  integer ipiv(3)

  real(doubleprec) t0wc,t1wc,t2wc,t3wc

  t0wc=MPI_WTIME()

!  *****************************************
!  calculate electron weight field
!

!  electron_locator is called in advance to give: wzeon, wteon0, wteon1, jteon0, jteon1, 
!  in pushe.F90 where wieght_field_eon.F90 is called
!  BE SURE no memory conflict (!!!)
!  call electron_locator         ! don't need to call electron_locator here

  f_rho=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)

  dnitmp=0.0_wp

#ifdef _OPENMP
!$omp parallel private(dampf)
  dampf=0.0_wp   ! Set array elements to zero
#endif
!$omp do private(m,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij,rdum,vz,mu,enew)
do m=1,me
     wz1=wzeon(m)
     wz0=1.0-wz1

        wp1=wpeon(m)
        wp0=1.0_wp-wp1

        wt10=wp0*wteon0(m)
        wt00=wp0-wt10

        wt11=wp1*wteon1(m)
        wt01=wp1-wt11
  
        vz=zeon(4,m)*BB(m)*f_rho           ! =v_//; BB is obtained from call gcmotion in pushi
        mu=zeon(6,m)*zeon(6,m)             ! =v_perp^2/2B
        enew=vz*vz+mu*BB(m)*2.0_wp*m0ms    ! =v^2

        rdum(1)=1.0_wp
        rdum(2)=vz
        rdum(3)=enew
        rdum(4)=vz*vz
        rdum(5)=vz*enew
        rdum(6)=enew*enew
        rdum(7)=zeon(5,m)
        rdum(8)=vz*zeon(5,m)
        rdum(9)=enew*zeon(5,m)

#ifdef _OPENMP
        ij=jteon0(m)
        dampf(1:9,ij,0) =dampf(1:9,ij,0) + wz0*wt00*rdum(1:9)
        dampf(1:9,ij,1) =dampf(1:9,ij,1) + wz1*wt00*rdum(1:9)

        ij=ij+1
        dampf(1:9,ij,0) =dampf(1:9,ij,0) + wz0*wt10*rdum(1:9)
        dampf(1:9,ij,1) =dampf(1:9,ij,1) + wz1*wt10*rdum(1:9)

        ij=jteon1(m)
        dampf(1:9,ij,0) =dampf(1:9,ij,0) + wz0*wt01*rdum(1:9)
        dampf(1:9,ij,1) =dampf(1:9,ij,1) + wz1*wt01*rdum(1:9)

        ij=ij+1
        dampf(1:9,ij,0) =dampf(1:9,ij,0) + wz0*wt11*rdum(1:9)
        dampf(1:9,ij,1) =dampf(1:9,ij,1) + wz1*wt11*rdum(1:9)
#else
        ij=jteon0(m)
        dnitmp(1:9,ij,0) =dnitmp(1:9,ij,0) + wz0*wt00*rdum(1:9)
        dnitmp(1:9,ij,1) =dnitmp(1:9,ij,1) + wz1*wt00*rdum(1:9)

        ij=ij+1
        dnitmp(1:9,ij,0) =dnitmp(1:9,ij,0) + wz0*wt10*rdum(1:9)
        dnitmp(1:9,ij,1) =dnitmp(1:9,ij,1) + wz1*wt10*rdum(1:9)

        ij=jteon1(m)
        dnitmp(1:9,ij,0) =dnitmp(1:9,ij,0) + wz0*wt01*rdum(1:9)
        dnitmp(1:9,ij,1) =dnitmp(1:9,ij,1) + wz1*wt01*rdum(1:9)

        ij=ij+1
        dnitmp(1:9,ij,0) =dnitmp(1:9,ij,0) + wz0*wt11*rdum(1:9)
        dnitmp(1:9,ij,1) =dnitmp(1:9,ij,1) + wz1*wt11*rdum(1:9)
#endif

enddo    ! enddo  m=1,me

#ifdef _OPENMP
!$omp critical
    dnitmp(1:9,1:mgrid,0:1) = dnitmp(1:9,1:mgrid,0:1) + dampf(1:9,1:mgrid,0:1)
!$omp end critical
#endif
!$omp end parallel

  t1wc=MPI_WTIME()

! If we have a particle decomposition on the toroidal domains, do a reduce
! operation to add up all the contributions to charge density on the grid
  if(npartdom>1)then
   !$omp parallel do private(ij)
     dampf=0._wp
     call MPI_ALLREDUCE(dnitmp,dampf,(9*mgrid*2),mpi_Rsize,&
                        MPI_SUM,partd_comm,ierror)
  endif

! poloidal end cell
  do j=0,mpsi
     dampf(:,igrid(j)+mtheta(j),:)=dampf(:,igrid(j)+mtheta(j),:)+dampf(:,igrid(j),:)
     dampf(:,igrid(j),:)=dampf(:,igrid(j)+mtheta(j),:)
  enddo

!do i=1,9  ! corresponding 9 coefficients using do-loop to save memory

! toroidal end cell
! Use dnitmp as MPI send {dnitmp(:,:,0)} and receive {dnitmp(:,:,1)} buffers
   dnitmp(:,:,0)=dampf(:,:,0)  ! send buffer
   dnitmp(:,:,1)=0.0_wp        ! recv buffer
   icount=9*mgrid
   idest=left_pe
   isource=right_pe
   isendtag=myrank_toroidal
   irecvtag=isource

! send pflux_e to left and receive from right
   call MPI_SENDRECV(dnitmp(:,:,0),icount,mpi_Rsize,idest,isendtag,dnitmp(:,:,1),&
                   icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

  if(myrank_toroidal == ntoroidal-1) then
! B.C. at zeta=2*pi is shifted
     do j=0,mpsi
        ii=igrid(j)
        jt=mtheta(j)
        dampf(:,ii+1:ii+jt,1)=dampf(:,ii+1:ii+jt,1)+&
              cshift(dnitmp(:,ii+1:ii+jt,1),itran(j),2) 
     enddo 
  else
! B.C. at zeta<2*pi is continuous
     dampf(:,:,1)=dampf(:,:,1)+dnitmp(:,:,1)
  endif

!end do   ! do i=1,9

!  NOTE: dampf(:,:,1) has valid data, but no valid data for dampf(:,:,0) and dampf(:,igrid(j),:)

! identical grid points:   igrid(j)=igrid(j)+mtheta(j)
  do j=0,mpsi
     dampf(:,igrid(j),:)=dampf(:,igrid(j)+mtheta(j),:)
  enddo

  t2wc=MPI_WTIME()

  do ij=1,mgrid
     aav(1,1)=dampf(1,ij,1)
     aav(2,1)=dampf(2,ij,1)
     aav(2,2)=dampf(4,ij,1)
     aav(3,1)=dampf(3,ij,1)
     aav(3,2)=dampf(5,ij,1)
     aav(3,3)=dampf(6,ij,1)
     bv(1)=dampf(7,ij,1)
     bv(2)=dampf(8,ij,1)
     bv(3)=dampf(9,ij,1)

!     if(mype==100) print*, 'Weight_Damp=', aav(1,1), bv(1), bv(2), bv(3), ij

     call dsytrf('L',3,aav,3,ipiv,work,lwork,info)
     call dsytrs('L',3,1,aav,3,ipiv,bv,3,info)

!     if(info < 0) then 
!        if(myrank_partd==1) print*, 'INFO=', info, 'ij=', ij, myrank_toroidal
!     end if

     wfield(1,ij,1)=bv(1)   
     wfield(2,ij,1)=bv(2)   
     wfield(3,ij,1)=bv(3)   
  end do

! toroidal BC: send to right and receive from left
!  do j=1,3
     send=wfield(:,:,1)
     recv(:,:)=0.0_wp
     icount=mgrid*3
     idest=right_pe
     isource=left_pe
     isendtag=myrank_toroidal
     irecvtag=isource
     call MPI_SENDRECV(send,icount,mpi_Rsize,idest,isendtag,recv,icount, &
                   mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

     if(myrank_toroidal==0)then
        do i=0,mpsi
           ii=igrid(i)
           jt=mtheta(i)
           wfield(1:3,ii+1:ii+jt,0)=cshift(recv(1:3,ii+1:ii+jt),-itran(i),2)
           !wfield(2,ii+1:ii+jt,0)=cshift(recv(2,ii+1:ii+jt),-itran(i),2)
           !wfield(3,ii+1:ii+jt,0)=cshift(recv(3,ii+1:ii+jt),-itran(i),2)
        enddo
     else
        wfield(:,:,0)=recv
     endif
!  end do

! now both weight-field  wfield(:,:,1) and  wfield(:,:,0) have valid data

  t3wc=MPI_WTIME()

  if(mype==0)write(0,*)istep,'  in weight_field_eon: t1=',(t1wc-t0wc),'  t2=',&
     (t2wc-t1wc),'  t3=',(t3wc-t2wc),'  tot=',(t3wc-t0wc)

end subroutine weight_field_eon


