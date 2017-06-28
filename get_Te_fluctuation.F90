subroutine get_Te_fluctuation(wpi2)
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface, only: BB

  implicit none

  integer m,i,j,k,jt,ii,ierror,ij,&
          icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
  real(wp) weight,rdum,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00
  real(wp) sendl(mgrid),recvr(mgrid)
  real(wp) dnitmp(0:1,mgrid),buff(0:1,mgrid)

  real(wp) upara, f_rho

  real(wp) wpi2(me)

  heatflux=0.0_wp   ! borrow array heatflux(mgrid) to calculate Te fluctuations
  buff=0.0_wp

#ifdef USE_OPENMP
!$omp parallel private(dnitmp)
  dnitmp=0._wp   ! Set array elements to zero
#endif
!$omp do private(m,weight,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij)
  do m=1,me
     !  weight=zeon0(5,m)*wpi2(m)    ! calculation using zeon0 at irk_e=1 determined by idiag_e=0
     !                               ! wzeon, wpeon, wteon0, wteon1 etc are calculated in pusheon
     weight=zeon0(5,m)*(wpi2(m)*2.0_wp/3.0_wp-backTe(m)*msm0)    ! calculation using zeon0 at irk_e=1 determined by idiag_e=0
                                                                 ! wzeon, wpeon, wteon0, wteon1 etc are calculated in pusheon

     wz1=weight*wzeon(m)
     wz0=weight-wz1     

        wp1=wpeon(m)
        wp0=1.0-wp1

        wt10=wp0*wteon0(m)
        wt00=wp0-wt10

        wt11=wp1*wteon1(m)
        wt01=wp1-wt11

#ifdef USE_OPENMP
! Use thread-private temporary array dnitmp to store the results
        ij=jteon0(m)
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt00
        dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt00
        
        ij=ij+1
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt10
        dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt10

        ij=jteon1(m)
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt01
        dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt01
        
        ij=ij+1
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt11
        dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt11
#else
! If no loop-level parallelism, use original algorithm (write directly
! into array "buff()".
        ij=jteon0(m)
        buff(0,ij) = buff(0,ij) + wz0*wt00
        buff(1,ij) = buff(1,ij) + wz1*wt00

        ij=ij+1
        buff(0,ij) = buff(0,ij) + wz0*wt10
        buff(1,ij) = buff(1,ij) + wz1*wt10

        ij=jteon1(m)
        buff(0,ij) = buff(0,ij) + wz0*wt01
        buff(1,ij) = buff(1,ij) + wz1*wt01

        ij=ij+1
        buff(0,ij) = buff(0,ij) + wz0*wt11
        buff(1,ij) = buff(1,ij) + wz1*wt11
#endif

  enddo

#ifdef USE_OPENMP
! For OpenMP, we now write the results accumulated in each thread-private
! array dnitmp() back into the shared array buff(). The loop is enclosed
! in a critical section so that one thread at a time updates buff(). 
!$omp critical
  do ij=1,mgrid
     buff(0,ij) = buff(0,ij) + dnitmp(0,ij)
     buff(1,ij) = buff(1,ij) + dnitmp(1,ij)
  enddo
!$omp end critical
!$omp end parallel
#endif

! If we have a particle decomposition on the toroidal domains, do a reduce
! operation to add up all the contributions to temperature on the grid
  if(npartdom>1)then
   !$omp parallel do private(ij)
    do ij=1,mgrid
       dnitmp(0,ij)=buff(0,ij)
       dnitmp(1,ij)=buff(1,ij)
    enddo
    buff=0._wp
    call MPI_ALLREDUCE(dnitmp,buff,(mgrid*2),mpi_Rsize,&
                       MPI_SUM,partd_comm,ierror)
  endif

! poloidal end cell
  do i=0,mpsi
     buff(:,igrid(i)+mtheta(i))=buff(:,igrid(i)+mtheta(i))+buff(:,igrid(i))
     buff(:,igrid(i))=buff(:,igrid(i)+mtheta(i))  
  enddo

! toroidal end cell
  sendl=buff(0,:)
  recvr=0.0_wp
  icount=mgrid
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource

! send Te to left and receive from right
  call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,recvr,icount,&
                    mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
     
  if(myrank_toroidal == ntoroidal-1)then
! B.C. at zeta=2*pi is shifted
     do i=0,mpsi
        ii=igrid(i)
        jt=mtheta(i)
        buff(1,ii+1:ii+jt)=buff(1,ii+1:ii+jt)+&
             cshift(recvr(ii+1:ii+jt),itran(i))
     enddo
  else
! B.C. at zeta<2*pi is continuous
     buff(1,:)=buff(1,:)+recvr(:)
  endif

!  NOTE: buff(1,:) has valid data, but no valid data for buff(0,:)

!$omp parallel do private(ij)
  do ij=1,mgrid
     buff(1,ij)=buff(1,ij)*markere(ij)  ! Te fluctuations on grid
  enddo

  do j=0,mpsi
     buff(1,igrid(j))=buff(1,igrid(j)+mtheta(j))

     rdum=1.0/Te_grid(j)         ! = 1/Te_0 with Te_0 normalized with same unit as Ti
     !  rdum=rdum*2.0_wp/3.0_wp
     !$omp parallel do private(i,ij)
     do i=0,mtheta(j)
        ij=igrid(j)+i
        buff(1,ij)=buff(1,ij)*rdum    ! = delta-Te/Te_0 on grid   !  added Oct. 2012
     end do
  enddo

  heatflux(1:mgrid)= buff(1,1:mgrid)   ! use heatflux(1:mgrid) to store Te fluctuations and write to output DENs at SNAPSHOT



!!!!!!!!!!!!!!!!!!!!! Calculate electron current fluctuations on grids !!!!!!!!!!!!!!!!!!!!!!!!!!!!


  jflux=0.0_wp   ! jflux(mgrid) to store B*j_e// fluctuations on grids
  buff=0.0_wp

  f_rho=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)

#ifdef USE_OPENMP
!$omp parallel private(dnitmp)
  dnitmp=0._wp   ! Set array elements to zero
#endif
!$omp do private(m,upara,weight,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij)
  do m=1,me
     upara=zeon0(4,m)*BB(m)*f_rho   ! normalized v_//
     weight=zeon0(5,m)*upara*BB(m)  ! electron bootstrap current = v_//*B

     wz1=weight*wzeon(m)
     wz0=weight-wz1     

        wp1=wpeon(m)
        wp0=1.0-wp1

        wt10=wp0*wteon0(m)
        wt00=wp0-wt10

        wt11=wp1*wteon1(m)
        wt01=wp1-wt11

#ifdef USE_OPENMP
! Use thread-private temporary array dnitmp to store the results
        ij=jteon0(m)
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt00
        dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt00
        
        ij=ij+1
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt10
        dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt10

        ij=jteon1(m)
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt01
        dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt01
        
        ij=ij+1
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt11
        dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt11
#else
! If no loop-level parallelism, use original algorithm (write directly
! into array "buff()".
        ij=jteon0(m)
        buff(0,ij) = buff(0,ij) + wz0*wt00
        buff(1,ij) = buff(1,ij) + wz1*wt00

        ij=ij+1
        buff(0,ij) = buff(0,ij) + wz0*wt10
        buff(1,ij) = buff(1,ij) + wz1*wt10

        ij=jteon1(m)
        buff(0,ij) = buff(0,ij) + wz0*wt01
        buff(1,ij) = buff(1,ij) + wz1*wt01

        ij=ij+1
        buff(0,ij) = buff(0,ij) + wz0*wt11
        buff(1,ij) = buff(1,ij) + wz1*wt11
#endif

  enddo

#ifdef USE_OPENMP
! For OpenMP, we now write the results accumulated in each thread-private
! array dnitmp() back into the shared array buff(). The loop is enclosed
! in a critical section so that one thread at a time updates buff(). 
!$omp critical
  do ij=1,mgrid
     buff(0,ij) = buff(0,ij) + dnitmp(0,ij)
     buff(1,ij) = buff(1,ij) + dnitmp(1,ij)
  enddo
!$omp end critical
!$omp end parallel
#endif

! If we have a particle decomposition on the toroidal domains, do a reduce
! operation to add up all the contributions to temperature on the grid
  if(npartdom>1)then
   !$omp parallel do private(ij)
    do ij=1,mgrid
       dnitmp(0,ij)=buff(0,ij)
       dnitmp(1,ij)=buff(1,ij)
    enddo
    buff=0._wp
    call MPI_ALLREDUCE(dnitmp,buff,(mgrid*2),mpi_Rsize,&
                       MPI_SUM,partd_comm,ierror)
  endif

! poloidal end cell
  do i=0,mpsi
     buff(:,igrid(i)+mtheta(i))=buff(:,igrid(i)+mtheta(i))+buff(:,igrid(i))
     buff(:,igrid(i))=buff(:,igrid(i)+mtheta(i))  
  enddo

! toroidal end cell
  sendl=buff(0,:)
  recvr=0.0_wp
  icount=mgrid
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource

! send Te to left and receive from right
  call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,recvr,icount,&
                    mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
     
  if(myrank_toroidal == ntoroidal-1)then
! B.C. at zeta=2*pi is shifted
     do i=0,mpsi
        ii=igrid(i)
        jt=mtheta(i)
        buff(1,ii+1:ii+jt)=buff(1,ii+1:ii+jt)+&
             cshift(recvr(ii+1:ii+jt),itran(i))
     enddo
  else
! B.C. at zeta<2*pi is continuous
     buff(1,:)=buff(1,:)+recvr(:)
  endif

!  NOTE: buff(1,:) has valid data, but no valid data for buff(0,:)

!$omp parallel do private(ij)
  do ij=1,mgrid
     buff(1,ij)=buff(1,ij)*markere(ij)  !  ~ (B*j_e//)/n_0 
  enddo

  do j=0,mpsi
     buff(1,igrid(j))=buff(1,igrid(j)+mtheta(j))

     rdum=dnty_on_bin_g(j)         
     !$omp parallel do private(i,ij)
     do i=0,mtheta(j)
        ij=igrid(j)+i
        buff(1,ij)=buff(1,ij)*rdum    ! ~ B*j_e// fluctuations on grid   !  added  2017
     end do
  enddo

  jflux(1:mgrid)= buff(1,1:mgrid)   ! use ijflux(1:mgrid) to store B*j_e// fluctuations and write to output DENs at SNAPSHOT


end subroutine get_Te_fluctuation
