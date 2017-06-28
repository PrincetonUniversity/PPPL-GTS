subroutine current_ion
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface
  implicit none

  integer m,i,j,k,ip,jt,kk,ii,ierror,larmor,ij,&
       icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
  real(wp) weight,rdum,r,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00
  real(wp) sendl(mgrid),recvr(mgrid)

  real(wp) e1,e2,e3,wpi(3,mi),dgr(mi),da(mi),dgq(mi), &
       dgf(mi),Er_particle(mi),dela,fip
  real(wp) dnitmp(0:1,mgrid),buff(0:1,mgrid)

  integer gcmotion

!  **************************************************************
!  calculate ion particle flux
 
! ion locations [wzion, wpion, wtion0, jtion0, wtion1, jtion1] are calculated in CHARGE_ION 

     pflux_i=0.0_wp

! gather e_field using 4-point gyro-averaging, sorting in poloidal angle
! The following line is an OpenMP directive for loop-level parallelism
! on shared memory machines (see http://www.openmp.org).
!$omp parallel do private(m,e1,e2,e3,wz1,wz0,larmor,ij,wp0,wt00,wt10,wp1,wt01,wt11)
  do m=1,mi
     e1=0.0_wp
     e2=0.0_wp
     e3=0.0_wp
     wz1=wzion(m)
     wz0=1.0_wp-wz1

     do larmor=1,4

        ij=jtion0(larmor,m)
        wp0=1.0_wp-wpion(larmor,m)
        wt00=1.0_wp-wtion0(larmor,m)
        e1=e1+wp0*wt00*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
        e2=e2+wp0*wt00*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
        e3=e3+wp0*wt00*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
        
        ij=ij+1
        wt10=1.0_wp-wt00
        e1=e1+wp0*wt10*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
        e2=e2+wp0*wt10*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
        e3=e3+wp0*wt10*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
        
        ij=jtion1(larmor,m)
        wp1=1.0_wp-wp0
        wt01=1.0_wp-wtion1(larmor,m)
        e1=e1+wp1*wt01*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
        e2=e2+wp1*wt01*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
        e3=e3+wp1*wt01*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
        
        ij=ij+1
        wt11=1.0_wp-wt01
        e1=e1+wp1*wt11*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
        e2=e2+wp1*wt11*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
        e3=e3+wp1*wt11*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

     enddo

     wpi(1,m)=0.25_wp*e1
     wpi(2,m)=0.25_wp*e2
     wpi(3,m)=0.25_wp*e3

  enddo

!    prepare equilibrium Er=dPhi/da at particles 
  dela=1.0_wp/Ydr
!$omp parallel do private(m,r,ip,fip)
  do m=1,mi
     r=zion(1,m)-a0      
     r=r*dela
     ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
     fip=r-real(ip,wp)
     Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particles
  enddo

  m=gcmotion(dgr,da,dgq,dgf,zion,wpi,Er_particle,nonlinear,mi)  ! should use nonlinear=1.0 here for linear run? --> calculate current

!     OUTPUT: dgr=drho/dt;       da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt; 
!             dgf=dzeta/dt;      ghq(:)=da/dt|_EXB;     wpi(3,:)=dphi/dzeta, NO CHANGE; 
!             BB(:)=B;           no change for zion(1,2,4,:), wpi(1,2,:) 
!             gha(:)=1/|grad a|

 do i=1,3  ! corresponding to a, theta, zeta; using do-loop to save memory

#ifdef USE_OPENMP
!$omp parallel private(dnitmp)
  dnitmp=0.   ! Set array elements to zero
#endif

!$omp do private(m,larmor,weight,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij,rdum)
  do m=1,mi
     weight=zion(5,m)

     wz1=weight*wzion(m)
     wz0=weight-wz1

     if(i==1) then
        rdum=da(m)
     else if(i==2) then
        rdum=dgq(m)
     else
        rdum=dgf(m)
     end if

     wz1=wz1*rdum
     wz0=wz0*rdum
     
     do larmor=1,4
        wp1=wpion(larmor,m)
        wp0=1.0_wp-wp1

        wt10=wp0*wtion0(larmor,m)
        wt00=wp0-wt10

        wt11=wp1*wtion1(larmor,m)
        wt01=wp1-wt11

#ifdef USE_OPENMP
! Use thread-private temporary array dnitmp to store the results
        ij=jtion0(larmor,m)

        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt00
        dnitmp(1,ij)   = dnitmp(1,ij)   + wz1*wt00


        ij=ij+1
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt10
        dnitmp(1,ij)   = dnitmp(1,ij)   + wz1*wt10

        ij=jtion1(larmor,m)
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt01
        dnitmp(1,ij)   = dnitmp(1,ij)   + wz1*wt01

        ij=ij+1
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt11
        dnitmp(1,ij)   = dnitmp(1,ij)   + wz1*wt11
#else
! If no loop-level parallelism, use original algorithm
! (write directly into array "pflux_i()".
        ij=jtion0(larmor,m)
        pflux_i(i,0,ij) = pflux_i(i,0,ij) + wz0*wt00
        pflux_i(i,1,ij)   = pflux_i(i,1,ij)   + wz1*wt00

        ij=ij+1
        pflux_i(i,0,ij) = pflux_i(i,0,ij) + wz0*wt10
        pflux_i(i,1,ij)   = pflux_i(i,1,ij)   + wz1*wt10

        ij=jtion1(larmor,m)
        pflux_i(i,0,ij) = pflux_i(i,0,ij) + wz0*wt01
        pflux_i(i,1,ij)   = pflux_i(i,1,ij)   + wz1*wt01

        ij=ij+1
        pflux_i(i,0,ij) = pflux_i(i,0,ij) + wz0*wt11
        pflux_i(i,1,ij)   = pflux_i(i,1,ij)   + wz1*wt11
#endif

     enddo
  enddo

#ifdef USE_OPENMP
!$omp critical
  do ij=1,mgrid
     pflux_i(i,0,ij) = pflux_i(i,0,ij) + dnitmp(0,ij)
     pflux_i(i,1,ij) = pflux_i(i,1,ij) + dnitmp(1,ij)
  enddo
!$omp end critical
!$omp end parallel   ! match "!$omp parallel private(dnitmp)"
#endif

! If we have a particle decomposition on the toroidal domains, do a reduce
! operation to add up all the contributions to charge density on the grid
  if(npartdom>1)then
   !$omp parallel do private(ij)
    do ij=1,mgrid
       dnitmp(0,ij)=pflux_i(i,0,ij)
       dnitmp(1,ij)=pflux_i(i,1,ij)
    enddo
    buff=0._wp
    call MPI_ALLREDUCE(dnitmp,buff,(mgrid*2),mpi_Rsize, MPI_SUM,&
                       partd_comm,ierror)
    do ij=1,mgrid
       pflux_i(i,0,ij)=buff(0,ij)
       pflux_i(i,1,ij)=buff(1,ij)
    enddo
  endif

! poloidal end cell
  do j=0,mpsi
     pflux_i(i,:,igrid(j)+mtheta(j))=pflux_i(i,:,igrid(j)+mtheta(j))+pflux_i(i,:,igrid(j))
     pflux_i(i,:,igrid(j))=pflux_i(i,:,igrid(j)+mtheta(j))
  enddo

! toroidal end cell
  sendl=pflux_i(i,0,:)
  recvr=0.0_wp
  icount=mgrid
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource

! send pflux_i to left and receive from right
  call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,&
       recvr,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

  if(myrank_toroidal == ntoroidal-1) then
! B.C. at zeta=2*pi is shifted
     do j=0,mpsi
        ii=igrid(j)
        jt=mtheta(j)
        pflux_i(i,1,ii+1:ii+jt)=pflux_i(i,1,ii+1:ii+jt)+&
             cshift(recvr(ii+1:ii+jt),itran(j))
     enddo
  else
! B.C. at zeta<2*pi is continuous
     pflux_i(i,1,:)=pflux_i(i,1,:)+recvr(:)
  endif

 end do   ! do i=1,3

! zero out particle flux in radial boundary cell
!  do i=0,nbound-1
!     pflux_i(:,:,igrid(i):igrid(i)+mtheta(i))=&
!        pflux_i(:,:,igrid(i):igrid(i)+mtheta(i))*float(i)/float(nbound)
!     pflux_i(:,:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))=&
!        pflux_i(:,:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))*float(i)/float(nbound)
!  enddo

!  ! get Gamma_i/n_0 for using the following approximation later:
!  !     d/da(Gamma)/n_0=d/da(Gamma/n_0) (consistent with GK ordering)
!  !     so as to save storage for n_0 at grid points
!  !$omp parallel do private(ij)
!    do ij=1,mgrid
!       pflux_i(1,1,ij)=pflux_i(1,1,ij)*markeri(ij)*0.25  ! Gamma_e/n_0
!       pflux_i(2,1,ij)=pflux_i(2,1,ij)*markeri(ij)*0.25  ! Gamma_e/n_0
!       pflux_i(3,1,ij)=pflux_i(3,1,ij)*markeri(ij)*0.25  ! Gamma_e/n_0
!  !                factor 0.25=1/4 comes from 4-point average
!    enddo

!$omp parallel do private(i,j,ij)
  do i=0,mpsi
     do j=1,mtheta(i)
        ij=igrid(i)+j
        pflux_i(1,1,ij)=0.25_wp*pflux_i(1,1,ij)*cloud_i/dadtdz(i)  ! Gamma_e*J
        pflux_i(2,1,ij)=0.25_wp*pflux_i(2,1,ij)*cloud_i/dadtdz(i)  ! Gamma_e*J
        pflux_i(3,1,ij)=0.25_wp*pflux_i(3,1,ij)*cloud_i/dadtdz(i)  ! Gamma_e*J
!                factor 0.25=1/4 comes from 4-point average
     enddo
  enddo

! identical grid points:   igrid(j)=igrid(j)+mtheta(j)
  do j=0,mpsi
     pflux_i(:,:,igrid(j))=pflux_i(:,:,igrid(j)+mtheta(j))
  enddo

! NOTE: data on pflux_i(:,0,:) is not valid data, but we need it duirng gathering process

end subroutine current_ion
