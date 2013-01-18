subroutine weight_field(wfield)
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use esi_interface

  use diagnosis_array

  implicit none

  integer,parameter ::     lwork=100

  integer m,i,j,k,jt,ii,im,ierror,ij,j00,j01,info,  &
       icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
  real(wp) rdum(9),wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,vz,mu,enew,rdum0,rdum1 
  real(wp) send(3,mgrid),recv(3,mgrid)

  real(wp) dnitmp(9,mgrid,0:1), dampf(9,mgrid,0:1), wfield(3,mgrid,0:1)
  real(wp) aav(3,3),bv(3),work(lwork)
  integer ipiv(3)

  integer displacement(0:nproc_partd-1),n_points(0:nproc_partd-1)
  integer npt_tot,my_npts

  real(doubleprec) t0wc,t1wc,t2wc,t3wc,t4wc
  real(wp) tdum,delr,delt(0:mpsi),delz,pi2_inv,psitmp,thetatmp,zetatmp

  t0wc=MPI_WTIME()

!  *****************************************
!  calculate ion weight field
!
!  BORROW array wzion, wtion0, wtion1, jtion0, jtion1
!  BE SURE no memory conflict (!!!)
!

  delr=1.0_wp/deltar
  delt=2.0_wp*pi/deltat
  delz=1.0_wp/deltaz
  pi2_inv=0.5_wp/pi

! ion_locator is called to give ion locations 
!!!!!!!  call ion_locator
! No need to call the ion locator. We now do this calculation in the loop below.

!$omp parallel do private(m,psitmp,thetatmp,zetatmp,vz,mu,enew,rdum,wz1,wz0,rdum0, &
!$omp& ii,wp1,wp0,im,tdum,j00,ij,wt10,wt00,j01,wt11,wt01) reduction(+:dnitmp)
do m=1,mi
     psitmp=zion(7,m)       ! =s
     thetatmp=zion(2,m)
     zetatmp=zion(3,m)

     vz=zion(4,m)*BB(m)          ! =v_//; BB is obtained from call gcmotion in pushi
     mu=zion(6,m)*zion(6,m)       ! =v_perp^2/2B
     enew=vz*vz+mu*BB(m)*2.0_wp   ! =v^2

     rdum(1)=1.0_wp
     rdum(2)=vz
     rdum(3)=enew
     rdum(4)=vz*vz
     rdum(5)=vz*enew
     rdum(6)=enew*enew
     rdum(7)=zion(5,m)
     rdum(8)=vz*zion(5,m)
     rdum(9)=enew*zion(5,m)

     wz1=(zetatmp-zetamin)*delz
     wz0=1.0-wz1

     rdum0=delr*max(0.0_wp,min(s1-s0,psitmp-s0))
     ii=max(0,min(mpsi-1,int(rdum0)))
     wp1=rdum0-real(ii,wp)
     wp0=1.0_wp-wp1

! inner flux surface
     im=ii
     tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
     tdum=(tdum-aint(tdum))*delt(im)
     j00=max(0,min(mtheta(im)-1,int(tdum)))
     ij=igrid(im)+j00  ! pay attention to j00=0 case
     wt10=wp0*(tdum-real(j00,wp))
     wt00=wp0-wt10

     dnitmp(1:9,ij,0) =dnitmp(1:9,ij,0) + wz0*wt00*rdum(1:9)
     dnitmp(1:9,ij,1) =dnitmp(1:9,ij,1) + wz1*wt00*rdum(1:9)

     ij=ij+1
     dnitmp(1:9,ij,0) =dnitmp(1:9,ij,0) + wz0*wt10*rdum(1:9)
     dnitmp(1:9,ij,1) =dnitmp(1:9,ij,1) + wz1*wt10*rdum(1:9)

! outer flux surface
     im=ii+1
     tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
     tdum=(tdum-aint(tdum))*delt(im)
     j01=max(0,min(mtheta(im)-1,int(tdum)))
     ij=igrid(im)+j01  ! pay attention to j01=0 case
     wt11=wp1*(tdum-real(j01,wp))
     wt01=wp1-wt11
  
     dnitmp(1:9,ij,0) =dnitmp(1:9,ij,0) + wz0*wt01*rdum(1:9)
     dnitmp(1:9,ij,1) =dnitmp(1:9,ij,1) + wz1*wt01*rdum(1:9)

     ij=ij+1
     dnitmp(1:9,ij,0) =dnitmp(1:9,ij,0) + wz0*wt11*rdum(1:9)
     dnitmp(1:9,ij,1) =dnitmp(1:9,ij,1) + wz1*wt11*rdum(1:9)

enddo    ! enddo  m=1,mi

  t1wc=MPI_WTIME()

! If we have a particle decomposition on the toroidal domains, do a reduce
! operation to add up all the contributions to charge density on the grid
  if(npartdom>1)then
     dampf=0._wp
     call MPI_ALLREDUCE(dnitmp,dampf,(9*mgrid*2),mpi_Rsize,&
                        MPI_SUM,partd_comm,ierror)
  endif
! poloidal end cell
  do j=0,mpsi
     dampf(:,igrid(j)+mtheta(j),:)=dampf(:,igrid(j)+mtheta(j),:)+dampf(:,igrid(j),:)
     dampf(:,igrid(j),:)=dampf(:,igrid(j)+mtheta(j),:)
  enddo

  t2wc=MPI_WTIME()


!do i=1,9  ! corresponding 9 coefficients using do-loop to save memory

! toroidal end cell
! Use dnitmp as MPI send {dnitmp(:,:,0)} and receive {dnitmp(:,:,1)} buffers
   dnitmp(:,:,0)=dampf(:,:,0)  ! send buffer
   dnitmp(:,:,1)=0.0_wp        ! receive buffer
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

!  NOTE: dampf(:,:,1) has valid data, but no valid data for dampf(:,:,0) and dampf(igrid(j),:,:)

! identical grid points:   igrid(j)=igrid(j)+mtheta(j)
  do j=0,mpsi
     dampf(:,igrid(j),1)=dampf(:,igrid(j)+mtheta(j),1)
  enddo






  w2_ion0=abs(w2_ion1-w2_ion0)/w2_ion0
  rdum0=maxval(w2_ion0)
  w2_ion0=w2_ion1



  t3wc=MPI_WTIME()

! We now split the work of the following loop between the nproc_partd processes that
! occupy the same toroidal domain and share the same planes. We divide the mgrid points
! as equally as possible. Since we mgrid might not be a multiple of "nproc_partd", the
! number of processes in a toroidal domain, we need to know the number of points and
! the displacement in the final array for all the processes in the domain. For each
! point we send 3 quantities.

  npt_tot=0
  do i=0,nproc_partd-1
     n_points(i)=mgrid/nproc_partd
     if(i < mod(mgrid,nproc_partd))n_points(i)=n_points(i)+1
     npt_tot=npt_tot+n_points(i)
  enddo

  if(npt_tot /= mgrid)then
     write(0,*)mype,'  in weight_damping_eon_new2: npt_tot /= mgrid , npt_tot=',npt_tot,'  mgrid=',mgrid
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  displacement(0)=0
  do i=1,nproc_partd-1
     displacement(i)=displacement(i-1)+n_points(i-1)
  enddo

  send=0.0_wp

! My starting point on the plane
  ii=displacement(myrank_partd)

  do ij=ii+1,ii+n_points(myrank_partd)
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

!     wfield(ij,1,1)=bv(1)   
!     wfield(ij,2,1)=bv(2)   
!     wfield(ij,3,1)=bv(3)   
     send(1,ij-ii)=bv(1)
     send(2,ij-ii)=bv(2)
     send(3,ij-ii)=bv(3)
  end do

! We now call MPI_GATHERV to collect the contributions from all the processes
! to make a full plane
! First, we need to reajust the number of elements sent and the displacement to
! take into account the 3 values per point being sent

  do i=0,nproc_partd-1
     n_points(i)=3*n_points(i)
     displacement(i)=3*displacement(i)
  enddo

  my_npts=n_points(myrank_partd)
  call MPI_ALLGATHERV(send,my_npts,mpi_Rsize,wfield(1,1,1),n_points,&
                      displacement,mpi_Rsize,partd_comm,ierror)




wfield(:,:,0)=0.0_wp    !!! ATTENTION: without this, cray compiler default option seems to casue problem for following MPI_SENDRECV
!
!  rdum0=maxval(abs(wfield))
!  if(rdum0>100.0_wp) then
!     print*, 'v7777=', mype, istep, irk, rdum0,  &
!                  maxval(abs(wfield(1,:,1))) , maxval(abs(wfield(2,:,1))) , maxval(abs(wfield(3,:,1))) ,  &
!                  maxval(abs(wfield(1,:,0))) , maxval(abs(wfield(2,:,0))) , maxval(abs(wfield(3,:,0)))
!  end if
!




! toroidal BC: send to right and receive from left
     send=wfield(:,:,1)
     recv=0.0_wp
     icount=3*mgrid
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
        enddo
     else
        wfield(:,:,0)=recv
     endif





! 
! rdum0=maxval(abs(wfield))
! if(rdum0>100.0_wp) then
!    print*, 'w7777=', mype, istep, irk, rdum0,  &
!                 maxval(abs(wfield(1,:,1))) , maxval(abs(wfield(2,:,1))) , maxval(abs(wfield(3,:,1))) ,  &
!                 maxval(abs(wfield(1,:,0))) , maxval(abs(wfield(2,:,0))) , maxval(abs(wfield(3,:,0)))
! end if  
! 





! now both weight-field  wfield(:,:,1) and  wfield(:,:,0) have valid data

  t4wc=MPI_WTIME()

  if(mype==0)write(0,*)istep,'  weight_field_ion: t1=',(t1wc-t0wc),'  t2=',&
     (t2wc-t1wc),'  t3=',(t3wc-t2wc),'  t4=',(t4wc-t3wc),'  tot=',(t4wc-t0wc)

end subroutine weight_field


