subroutine weight_field_eon(wfield)
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
  real(wp) rdum,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,vz,mu,enew,f_rho 
  real(wp) send(mgrid),recv(mgrid)

  real(wp) dnitmp(mgrid,9,0:1), dampf(mgrid,9,0:1), wfield(mgrid,3,0:1)
  real(wp) aav(3,3),bv(3),work(lwork)
  integer ipiv(3)

!  *****************************************
!  calculate electron weight field
!

!  electron_locator is called in advance to give: wzeon, wteon0, wteon1, jteon0, jteon1, 
!  in pushe.F90 where wieght_field_eon.F90 is called
!  BE SURE no memory conflict (!!!)
!  call electron_locator         ! don't need to call electron_locator here

  f_rho=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)

  dampf=0.0_wp

#ifdef _OPENMP
!$omp parallel private(dnitmp)
   dnitmp=0.0_wp   ! Set array elements to zero
#endif
!$omp do private(m,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,j,ij,rdum,vz,mu,enew)
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
 
#ifdef _OPENMP
! Use thread-private temporary array dnitmp to store the results
        do j=1,9
           if(j==1) rdum=1.0_wp 
           if(j==2) rdum=vz 
           if(j==3) rdum=enew
           if(j==4) rdum=vz*vz 
           if(j==5) rdum=vz*enew 
           if(j==6) rdum=enew*enew 
           if(j==7) rdum=zeon(5,m) 
           if(j==8) rdum=vz*zeon(5,m) 
           if(j==9) rdum=enew*zeon(5,m)

           ij=jteon0(m)
           dnitmp(ij,j,0) = dnitmp(ij,j,0) + wz0*wt00*rdum
           dnitmp(ij,j,1) = dnitmp(ij,j,1) + wz1*wt00*rdum

           ij=ij+1
           dnitmp(ij,j,0) = dnitmp(ij,j,0) + wz0*wt10*rdum
           dnitmp(ij,j,1) = dnitmp(ij,j,1) + wz1*wt10*rdum

           ij=jteon1(m)
           dnitmp(ij,j,0) = dnitmp(ij,j,0) + wz0*wt01*rdum
           dnitmp(ij,j,1) = dnitmp(ij,j,1) + wz1*wt01*rdum

           ij=ij+1
           dnitmp(ij,j,0) = dnitmp(ij,j,0) + wz0*wt11*rdum
           dnitmp(ij,j,1) = dnitmp(ij,j,1) + wz1*wt11*rdum
        end do
#else
! If no loop-level parallelism, use original algorithm (write directly
! into array "dampf(*,*,*)".
        do j=1,9
           if(j==1) rdum=1.0_wp
           if(j==2) rdum=vz
           if(j==3) rdum=enew
           if(j==4) rdum=vz*vz
           if(j==5) rdum=vz*enew
           if(j==6) rdum=enew*enew
           if(j==7) rdum=zeon(5,m)
           if(j==8) rdum=vz*zeon(5,m)
           if(j==9) rdum=enew*zeon(5,m)

           ij=jteon0(m)
           dampf(ij,j,0) = dampf(ij,j,0) + wz0*wt00*rdum
           dampf(ij,j,1) = dampf(ij,j,1) + wz1*wt00*rdum

           ij=ij+1
           dampf(ij,j,0) = dampf(ij,j,0) + wz0*wt10*rdum
           dampf(ij,j,1) = dampf(ij,j,1) + wz1*wt10*rdum

           ij=jteon1(m)
           dampf(ij,j,0) = dampf(ij,j,0) + wz0*wt01*rdum
           dampf(ij,j,1) = dampf(ij,j,1) + wz1*wt01*rdum

           ij=ij+1
           dampf(ij,j,0) = dampf(ij,j,0) + wz0*wt11*rdum
           dampf(ij,j,1) = dampf(ij,j,1) + wz1*wt11*rdum
        end do
#endif

enddo    ! enddo  m=1,me

#ifdef _OPENMP
!$omp critical
  do ij=1,mgrid
     dampf(ij,1:9,0:1) = dampf(ij,1:9,0:1) + dnitmp(ij,1:9,0:1)
  enddo
!$omp end critical
#endif
!$omp end parallel

! If we have a particle decomposition on the toroidal domains, do a reduce
! operation to add up all the contributions to charge density on the grid
  if(npartdom>1)then
   !$omp parallel do private(ij)
     do ij=1,mgrid
        dnitmp(ij,1:9,0:1)=dampf(ij,1:9,0:1)
     enddo
     dampf=0._wp
     call MPI_ALLREDUCE(dnitmp,dampf,(mgrid*9*2),mpi_Rsize,&
                        MPI_SUM,partd_comm,ierror)
  endif

! poloidal end cell
  do j=0,mpsi
     dampf(igrid(j)+mtheta(j),:,:)=dampf(igrid(j)+mtheta(j),:,:)+dampf(igrid(j),:,:)
     dampf(igrid(j),:,:)=dampf(igrid(j)+mtheta(j),:,:)
  enddo


do i=1,9  ! corresponding 9 coefficients using do-loop to save memory

! toroidal end cell
   send=dampf(:,i,0)
   recv=0.0_wp
   icount=mgrid
   idest=left_pe
   isource=right_pe
   isendtag=myrank_toroidal
   irecvtag=isource

! send pflux_e to left and receive from right
   call MPI_SENDRECV(send,icount,mpi_Rsize,idest,isendtag,recv,icount,&
                     mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

  if(myrank_toroidal == ntoroidal-1) then
! B.C. at zeta=2*pi is shifted
     do j=0,mpsi
        ii=igrid(j)
        jt=mtheta(j)
        dampf(ii+1:ii+jt,i,1)=dampf(ii+1:ii+jt,i,1)+&
              cshift(recv(ii+1:ii+jt),itran(j)) 
     enddo 
  else
! B.C. at zeta<2*pi is continuous
     dampf(:,i,1)=dampf(:,i,1)+recv(:) 
  endif

end do   ! do i=1,9

!  NOTE: dampf(:,:,1) has valid data, but no valid data for dampf(:,:,0) and dampf(igrid(j),:,:)

! identical grid points:   igrid(j)=igrid(j)+mtheta(j)
  do j=0,mpsi
     dampf(igrid(j),:,:)=dampf(igrid(j)+mtheta(j),:,:)
  enddo


  w2_eon0=abs(w2_eon1-w2_eon0)/w2_eon0
  rdum=maxval(w2_eon0)
  w2_eon0=w2_eon1
  
!$omp parallel do private(ij,aav,bv,info,ipiv,work)
  do ij=1,mgrid
     aav=0.0_wp
     bv=0.0_wp
  
     aav(1,1)=dampf(ij,1,1)
     aav(2,1)=dampf(ij,2,1)
     aav(2,2)=dampf(ij,4,1)
     aav(3,1)=dampf(ij,3,1)
     aav(3,2)=dampf(ij,5,1)
     aav(3,3)=dampf(ij,6,1)
     bv(1)=dampf(ij,7,1)
     bv(2)=dampf(ij,8,1)
     bv(3)=dampf(ij,9,1)

!! ---- DEBUG ----- ---- DEBUG ----- ---- DEBUG -----
!  if(myrank_partd==0) then
!     if(rdum > 0.3_wp) then
!        write(678, 200) ij, aav(1,1), aav(2,1), aav(2,2), aav(3,1), aav(3,2), aav(3,3), bv(1), bv(2), bv(3)
!        call flush(678)
!     end if
!  end if
!200 format(1x,i7,12(2x,d12.6))
!! ---- DEBUG ----- ---- DEBUG ----- ---- DEBUG -----

!     if(mype==100) print*, 'Weight_Damp=', aav(1,1), bv(1), bv(2), bv(3), ij

     info=0
     call dsytrf('L',3,aav,3,ipiv,work,lwork,info)
     info=0
     call dsytrs('L',3,1,aav,3,ipiv,bv,3,info)

!      if(info < 0) then 
!         print*, 'electron INFO=', info, 'ij=', ij, myrank_toroidal
!      end if

     wfield(ij,1,1)=bv(1)   
     wfield(ij,2,1)=bv(2)   
     wfield(ij,3,1)=bv(3)   
     wfield(ij,1,0)=0.0_wp
     wfield(ij,2,0)=0.0_wp   
     wfield(ij,3,0)=0.0_wp   
  end do

! toroidal BC: send to right and receive from left
  do j=1,3
     send=wfield(:,j,1)
     recv=0.0_wp
     icount=mgrid
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
           wfield(ii+1:ii+jt,j,0)=cshift(recv(ii+1:ii+jt),-itran(i))
        enddo
     else
        wfield(:,j,0)=recv
     endif
  end do

! now both weight-field  wfield(:,:,1) and  wfield(:,:,0) have valid data

end subroutine weight_field_eon


