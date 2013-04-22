subroutine weight_damping_ion
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use esi_interface

  implicit none

  integer m,ij
  real(wp) mu,muB,vpara,v2,e1,e2,e3,wz1,wz0,wp0,wt00,wt10,wp1,wt01,wt11,avw,  &
           eta_w

  real(wp) wfield(3,mgrid,0:1)

  !  eta_w=2.0d-7/tstep
  eta_w=weight_damping_rate/tstep

  call esigetb(zion,mi,7)   ! get B-->BB
  call weight_field(wfield)
  

  !  calculate average weight <w> at ions using wieght-field wfield and then weight damping term

!$omp parallel do private(m,mu,muB,vpara,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,v2,avw)
  do m=1,mi
     mu=zion(6,m)*zion(6,m)           ! =v_perp^2/2B
     muB=mu*BB(m)                     ! =v_perp^2/2
     vpara=zion(4,m)*BB(m)
     v2=vpara*vpara+muB*2.0_wp        ! v^2

     wz1=wzion(m)
     wz0=1.0_wp-wz1

     ij=jtion0(1,m)
     wp0=1.0_wp-wpion(1,m)

     wt00=1.0_wp-wtion0(1,m)
     e1=wp0*wt00*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=wp0*wt00*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=wp0*wt00*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt10=1.0_wp-wt00
     e1=e1+wp0*wt10*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp0*wt10*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp0*wt10*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=jtion1(1,m)
     wp1=1.0_wp-wp0
     wt01=1.0_wp-wtion1(1,m)
     e1=e1+wp1*wt01*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt01*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt01*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt11=1.0_wp-wt01
     e1=e1+wp1*wt11*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt11*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt11*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     avw=e1+e2*vpara+e3*v2

     zion(5,m) = avw+(zion(5,m)-avw)*exp(-eta_w*tstep)   ! dw/dt=-eta_w*(w-<w>)
! 
!      if(abs(avw)>100.0_wp) then
!         print*, 'p9999=', mype, istep, irk, avw, e1, e2, e3, vpara, v2, maxval(abs(wfield(1,1:mgrid,1))) , &
!                 maxval(abs(wfield(2,1:mgrid,1))) , maxval(abs(wfield(3,1:mgrid,1))) , maxval(abs(wfield(1,1:mgrid,0))) , & 
!                 maxval(abs(wfield(2,1:mgrid,0))) , maxval(abs(wfield(3,1:mgrid,0)))
!      end if
! 


  end do


end subroutine weight_damping_ion






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

  integer m,i,j,k,jt,ii,ierror,ij,info,  &
       icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
  real(wp) rdum(9),wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,vz,mu,enew,rdum0,rdum1 
  real(wp) send(3,mgrid),recv(3,mgrid)

  real(wp) dnitmp(9,mgrid,0:1), dampf(9,mgrid,0:1), wfield(3,mgrid,0:1)
  real(wp) aav(3,3),bv(3),work(lwork)
  integer ipiv(3)

  integer displacement(0:nproc_partd-1),n_points(0:nproc_partd-1)
  integer npt_tot,my_npts

  real(doubleprec) t0wc,t1wc,t2wc,t3wc,t4wc

!   t0wc=MPI_WTIME()

!  *****************************************
!  calculate ion weight field
!
!  BORROW array wzion, wtion0, wtion1, jtion0, jtion1
!  BE SURE no memory conflict (!!!)
!

! ion_locator is called to give ion locations 
  call ion_locator

  dnitmp=0.0_wp

#ifdef _OPENMP
!$omp parallel private(dampf)
  dampf=0.0_wp   ! Set array elements to zero
#endif
!$omp do private(m,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij,rdum,vz,mu,enew)
do m=1,mi
     wz1=wzion(m)
     wz0=1.0-wz1

        wp1=wpion(1,m)
        wp0=1.0_wp-wp1

        wt10=wp0*wtion0(1,m)
        wt00=wp0-wt10

        wt11=wp1*wtion1(1,m)
        wt01=wp1-wt11
  
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

#ifdef _OPENMP
        ij=jtion0(1,m)
        dampf(1:9,ij,0) =dampf(1:9,ij,0) + wz0*wt00*rdum(1:9)
        dampf(1:9,ij,1) =dampf(1:9,ij,1) + wz1*wt00*rdum(1:9)

        ij=ij+1
        dampf(1:9,ij,0) =dampf(1:9,ij,0) + wz0*wt10*rdum(1:9)
        dampf(1:9,ij,1) =dampf(1:9,ij,1) + wz1*wt10*rdum(1:9)

        ij=jtion1(1,m)
        dampf(1:9,ij,0) =dampf(1:9,ij,0) + wz0*wt01*rdum(1:9)
        dampf(1:9,ij,1) =dampf(1:9,ij,1) + wz1*wt01*rdum(1:9)

        ij=ij+1
        dampf(1:9,ij,0) =dampf(1:9,ij,0) + wz0*wt11*rdum(1:9)
        dampf(1:9,ij,1) =dampf(1:9,ij,1) + wz1*wt11*rdum(1:9)
#else
        ij=jtion0(1,m)
        dnitmp(1:9,ij,0) =dnitmp(1:9,ij,0) + wz0*wt00*rdum(1:9)
        dnitmp(1:9,ij,1) =dnitmp(1:9,ij,1) + wz1*wt00*rdum(1:9)

        ij=ij+1
        dnitmp(1:9,ij,0) =dnitmp(1:9,ij,0) + wz0*wt10*rdum(1:9)
        dnitmp(1:9,ij,1) =dnitmp(1:9,ij,1) + wz1*wt10*rdum(1:9)

        ij=jtion1(1,m)
        dnitmp(1:9,ij,0) =dnitmp(1:9,ij,0) + wz0*wt01*rdum(1:9)
        dnitmp(1:9,ij,1) =dnitmp(1:9,ij,1) + wz1*wt01*rdum(1:9)

        ij=ij+1
        dnitmp(1:9,ij,0) =dnitmp(1:9,ij,0) + wz0*wt11*rdum(1:9)
        dnitmp(1:9,ij,1) =dnitmp(1:9,ij,1) + wz1*wt11*rdum(1:9)
#endif

enddo    ! enddo  m=1,mi




! 
! if(istep>74) then
!   do ij=1,mgrid
! 
!      rdum0=maxval(abs(dnitmp(7:9,ij,0)))
!      if(rdum0 > 0.001_wp) then
!         print*, 'AAAAA=', mype, istep, irk, ij, dnitmp(7,ij,0), dnitmp(8,ij,0), dnitmp(9,ij,0)
!         rdum1=maxval(abs(zion(5,1:mi)))
!         print*, 'aaaaa=',mype, istep, irk, rdum1, mi
!      end if
! 
!      rdum0=maxval(abs(dnitmp(7:9,ij,1)))
!      if(rdum0 > 0.001_wp) then
!         print*, 'BBBBB=', mype, istep, irk, ij, dnitmp(7,ij,1), dnitmp(8,ij,1), dnitmp(9,ij,1)
!         rdum1=maxval(abs(zion(5,1:mi)))
!         print*, 'bbbbb=',mype, istep, irk, rdum1, mi
!      end if
! 
!   end do
! end if
! 




#ifdef _OPENMP
!$omp critical
  dnitmp(1:9,1:mgrid,0:1) = dnitmp(1:9,1:mgrid,0:1) + dampf(1:9,1:mgrid,0:1)
!$omp end critical
#endif
!$omp end parallel

!   t1wc=MPI_WTIME()

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

!   t2wc=MPI_WTIME()





! 
!  if(myrank_partd==0) then
!    do ij=1,mgrid
! 
!      rdum0=maxval(abs(dampf(7:9,ij,0)))
!      if(rdum0 > 0.01_wp) then
!         print*, 'XXXXX',  istep, irk, ij, dampf(7,ij,0), dampf(8,ij,0), dampf(9,ij,0)
!      end if
! 
!      rdum0=maxval(abs(dampf(7:9,ij,1)))
!      if(rdum0 > 0.01_wp) then
!         print*, 'YYYYY',  istep, irk, ij, dampf(7,ij,1), dampf(8,ij,1), dampf(9,ij,1)
!      end if
! 
!   end do
!  end if
!
!
!  do ij=1,mgrid
!
!     rdum0=maxval(abs(dampf(1:9,ij,0)))
!     if(rdum0 > 1000.0_wp) then
!        print*, 'xxxxx',  istep, irk, ij,rdum0
!     end if
!
!     rdum0=maxval(abs(dampf(1:9,ij,1)))
!     if(rdum0 > 1000.0_wp) then
!        print*, 'yyyyy',  istep, irk, ij, rdum0
!     end if
!
!  end do
!





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
     dampf(:,igrid(j),:)=dampf(:,igrid(j)+mtheta(j),:)
  enddo






  w2_ion0=abs(w2_ion1-w2_ion0)/w2_ion0
  rdum0=maxval(w2_ion0)
  w2_ion0=w2_ion1

!   do ij=1,mgrid
!
!! ---- DEBUG ----- ---- DEBUG ----- ---- DEBUG -----
!   if(myrank_partd==0) then
! !
! !     if(rdum0 > 0.3_wp) then
!         write(677, 200) ij, dampf(1,ij,1), dampf(2,ij,1), dampf(4,ij,1), dampf(3,ij,1), dampf(5,ij,1), dampf(6,ij,1), dampf(7,ij,1), dampf(8,ij,1), dampf(9,ij,1)
!         call flush(677)
! !     end if
!   end if
! 200 format(1x,i7,12(2x,d12.6))
!! ---- DEBUG ----- ---- DEBUG ----- ---- DEBUG -----
!
!   end do






!   t3wc=MPI_WTIME()

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

!   t4wc=MPI_WTIME()

!   if(mype==0)write(0,*)istep,'  weight_field_ion: t1=',(t1wc-t0wc),'  t2=',&
!      (t2wc-t1wc),'  t3=',(t3wc-t2wc),'  t4=',(t4wc-t3wc),'  tot=',(t4wc-t0wc)

end subroutine weight_field


