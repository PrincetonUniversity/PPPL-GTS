subroutine get_Te_fluctuation(wpi2)
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use field_array
   use diagnosis_array

   implicit none

   integer m,i,j,k,jt,ii,ierror,ij,&
      icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
   real(wp) weight,rdum,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00
   real(wp) sendl(mgrid),recvr(mgrid)
   real(wp) dnitmp(0:1,mgrid),buff(0:1,mgrid)
   real(wp) mype_zonal_dTe(0:mpsi) ! zonal component of dTe in my cpu part
   real(wp) wpi2(me)

   !call electron_locator

   heatflux=0.0_wp   ! borrow array heatflux(mgrid) to calculate Te fluctuations
   buff=0.0_wp
   dnitmp=0.0_wp
#ifdef _OPENMP
   !$omp parallel private(dnitmp)
   dnitmp=0._wp   ! Set array elements to zero
#endif
   !$omp do private(m,weight,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij)
   do m=1,me
      !  weight=zeon0(5,m)*wpi2(m)    ! calculation using zeon0 at irk_e=1 determined by idiag_e=0
      !                               ! wzeon, wpeon, wteon0, wteon1 etc are calculated in pusheon
      weight=zeon0(5,m)*(wpi2(m)*2.0_wp/3.0_wp-backTe(m)*msm0)    ! calculation using zeon0 at irk_e=1 determined by idiag_e=0
      ! wzeon, wpeon, wteon0, wteon1 etc are calculated in pusheon

      ! if(zeon0(5,m)==0 .or. wpi2(m)==0 .or. zeon0(5,m)<min_thold_weight .or. zeon0(5,m)>max_thold_weight) then
      !   write(stdout,*) 'm=',m, 'zeon0(5,m)=',zeon0(5,m), 'wpi2(m)=',wpi2(m)
      !   write(stdout,*) 'zeon0(:,m)=',zeon0(:,m)
      !   write(stdout,*) 'zeon(:,m)=',zeon(:,m)
      !   call flush(stdout)
      ! endif

      if (zeon0(5,m)<min_thold_weight .or. zeon0(5,m)>max_thold_weight) cycle

      wz1=weight*wzeon(m)
      wz0=weight-wz1

      wp1=wpeon(m)
      wp0=1.0-wp1

      wt10=wp0*wteon0(m)
      wt00=wp0-wt10

      wt11=wp1*wteon1(m)
      wt01=wp1-wt11

#ifdef _OPENMP
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

#ifdef _OPENMP
   ! For OpenMP, we now write the results accumulated in each thread-private
   ! array dnitmp() back into the shared array buff(). The loop is enclosed
   ! in a critical section so that one thread at a time updates buff().
   !$omp critical
   do ij=1,mgrid
      buff(0,ij) = buff(0,ij) + dnitmp(0,ij)
      buff(1,ij) = buff(1,ij) + dnitmp(1,ij)
   enddo
   !$omp end critical
#endif
   !$omp end parallel

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

   ! !==============================================================================
   ! ! MGYOO (01/09/2019)
   ! ! Convert buff to get a correct form of dTe by considering dne.
   ! ! buff = Sigma(weight*(2/3*energy-T0))*markere = (1/n0)*Sigma(weight*(2/3*energy-T0))/node_Vol
   ! ! dTe = (1/(n0+dn))*Sigma(weight*(2/3*energy-T0))/node_Vol = n0/(n0+dn)*buff
   ! do i=0,mpsi
   !   do j=0,mtheta(i)
   !     ij = igrid(i) + j
   !     buff(1,ij) = den_grid(i)/(den_grid(i)+tot_dne(1,ij))*buff(1,ij)
   !   enddo
   ! enddo
   ! !==============================================================================


   !====================================================
   ! MGYOO (11/27/2018) : implementation of zonal and non-zonal terms
   ! zonal_dTe is the zonal component of delta_Te
   !====================================================
   mype_zonal_dTe = 0.0
   do i=0,mpsi
      do j=1,mtheta(i) ! j=0 is ignored to avoid over-estimation at j=0 or j=mtheta(i) which is theta=0 or 2*pi
         ij=igrid(i)+j
         mype_zonal_dTe(i) = mype_zonal_dTe(i) + dtemper(ij)*buff(1,ij) !
      enddo
   enddo

   call MPI_ALLREDUCE(mype_zonal_dTe, zonal_dTe, mpsi+1, mpi_Rsize, MPI_SUM, toroidal_comm, ierror)
   zonal_dTe = zonal_dTe/jaco_sum

   !====================================================================
   ! MGYOO(12/03/2018) : To measure FLA(Field Line Average) component
   ! This FLA might be important in the case of external magnetic perturbation
   if((myrank_partd==nproc_partd-1) .and. (mod(istep,msnap) .eq. 0)) then
      call Do_FieldLineAverage_by_using_gPPI(buff, FLA_dTe)
   endif
   !===================================================================

   ! IF ZE_SEPARATION=0: heatflux is (1/Te0)*(full delta_Te (zonal+non-zonal)): i.e. heatflux = (delta_Te)/Te0
   ! IF ZE_SEPARATION=1: heatflux is now calculated from non-zonal delta_Te : i.e. heatflux = (delta_Te-zonal_dTe)/Te0
   do j=0,mpsi
      rdum=1.0/Te_grid(j)         ! = 1/Te_0 with Te_0 normalized with same unit as Ti
      !$omp parallel do private(i,ij)
      do i=0,mtheta(j)
         ij=igrid(j)+i
         buff(1,ij)=(buff(1,ij)-(ZF_SEPARATION)*zonal_dTe(j))*rdum    ! = delta-Te/Te_0 on grid   !  added Oct. 2012


         ! if( abs(buff(1,ij))<1e-4 ) then
         !   write(stdout,*) 'get_Te2, mype=',mype,',m=',m,'ij=',ij,'weight=',weight,'buff(:,ij)=',buff(:,ij)
         !   call flush(stdout)
         ! endif
      end do
      buff(1,igrid(j))=buff(1,igrid(j)+mtheta(j))
   enddo

   heatflux(1:mgrid)= buff(1,1:mgrid)   ! use heatflux(1:mgrid) to store Te fluctuations and write to output DENs at SNAPSHOT

end subroutine get_Te_fluctuation
