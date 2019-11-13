subroutine current_eon
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use field_array
   use diagnosis_array
   use esi_interface
   implicit none

   integer m,i,j,k,ip,jt,kk,ii,ierror,ij,&
      icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
   real(wp) weight,rdum,r,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00
   real(wp) sendl(mgrid),recvr(mgrid)

   real(wp) e1,e2,e3,wpi(3,me), &
      dgr(me),da(me),dgq(me),dgf(me), &
      Er_particle(me),dela,fip

   integer gcmotion_eon

   real(wp) dni(0:mpsi),dne(0:mpsi),tdni(0:mpsi),tdne(0:mpsi),sne,sni,tsne,tsni
   real(wp) dnitmp(0:1,mgrid),buff(0:1,mgrid)


   !  *****************************************
   !  calculate electron particle flux

   pflux_e=0.0_wp

   ! electron_locator is called in CHARGE_EON to give electron locations

   ! gather e_field at electron GC positions using linear interpolation. No FLR effect for electrons
   ! The following line is an OpenMP directive for loop-level parallelism
   ! on shared memory machines (see http://www.openmp.org).
   !$omp parallel do private(m,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11)
   !$acc parallel loop present(zeon,jteon0,jteon1,wzeon,wpeon,wteon0,wteon1)
   do m=1,me
      wz1=wzeon(m)
      wz0=1.0_wp-wz1

      ij=jteon0(m)
      wp0=1.0_wp-wpeon(m)
      wt00=1.0_wp-wteon0(m)
      e1=wp0*wt00*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
      e2=wp0*wt00*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
      e3=wp0*wt00*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

      ij=ij+1
      wt10=1.0_wp-wt00
      e1=e1+wp0*wt10*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
      e2=e2+wp0*wt10*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
      e3=e3+wp0*wt10*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

      ij=jteon1(m)
      wp1=1.0_wp-wp0
      wt01=1.0_wp-wteon1(m)
      e1=e1+wp1*wt01*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
      e2=e2+wp1*wt01*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
      e3=e3+wp1*wt01*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

      ij=ij+1
      wt11=1.0_wp-wt01
      e1=e1+wp1*wt11*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
      e2=e2+wp1*wt11*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
      e3=e3+wp1*wt11*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

      wpi(1,m)=e1
      wpi(2,m)=e2
      wpi(3,m)=e3
   enddo

   !    prepare equilibrium Er=dPhi/da at particles
   dela=1.0_wp/Ydr
   !$omp parallel do private(m,r,ip,fip)
   do m=1,me
      r=zeon(1,m)-a0
      r=r*dela
      ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
      fip=r-real(ip,wp)
      Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particles
   enddo

   !    m=gcmotion_e(dgr,da,dgq,dgf,zeon,wpi,Er_particle,nonlinear,me,zeondim)  ! should use nonlinear=1.0 here for linear run? --> calculate current
   m=gcmotion_eon(dgr,da,dgq,dgf,zeon,wpi,Er_particle,nonlinear,me,zeondim)  ! should use nonlinear=1.0 here for linear run? --> calculate current

   !     OUTPUT: dgr=drho/dt;       da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt;
   !             dgf=dzeta/dt;      ghq(:)=da/dt|_EXB;     wpi(3,:)=dphi/dzeta, NO CHANGE;
   !             BB(:)=B;           no change for zeon(1,2,4,:), wpi(1,2,:)
   !             gha(:)=1/|grad a|  F(:) =F

   do i=1,3  ! corresponding to a, theta, zeta; using do-loop to save memory

#ifdef _OPENMP
      !$omp parallel private(dnitmp)
      dnitmp=0.   ! Set array elements to zero
#endif
      !$omp do private(m,weight,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij,rdum)
      !$acc parallel loop present(zeon,jteon0,jteon1,wzeon,wpeon,wteon0,wteon1)
      do m=1,me
         weight=zeon(5,m)

         wz1=weight*wzeon(m)
         wz0=weight-wz1

         wp1=wpeon(m)
         wp0=1.0_wp-wp1

         wt10=wp0*wteon0(m)
         wt00=wp0-wt10

         wt11=wp1*wteon1(m)
         wt01=wp1-wt11

         if(i==1) then
            rdum=da(m)
         else if(i==2) then
            rdum=dgq(m)
         else
            rdum=dgf(m)
         end if

#ifdef _OPENMP
         ! Use thread-private temporary array dnitmp to store the results
         ij=jteon0(m)
         dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt00*rdum
         dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt00*rdum

         ij=ij+1
         dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt10*rdum
         dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt10*rdum

         ij=jteon1(m)
         dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt01*rdum
         dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt01*rdum

         ij=ij+1
         dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt11*rdum
         dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt11*rdum
#else
         ! If no loop-level parallelism, use original algorithm (write directly
         ! into array "pflux_e(*,*,*)".
         ij=jteon0(m)
         pflux_e(i,0,ij) = pflux_e(i,0,ij) + wz0*wt00*rdum
         pflux_e(i,1,ij) = pflux_e(i,1,ij) + wz1*wt00*rdum

         ij=ij+1
         pflux_e(i,0,ij) = pflux_e(i,0,ij) + wz0*wt10*rdum
         pflux_e(i,1,ij) = pflux_e(i,1,ij) + wz1*wt10*rdum

         ij=jteon1(m)
         pflux_e(i,0,ij) = pflux_e(i,0,ij) + wz0*wt01*rdum
         pflux_e(i,1,ij) = pflux_e(i,1,ij) + wz1*wt01*rdum

         ij=ij+1
         pflux_e(i,0,ij) = pflux_e(i,0,ij) + wz0*wt11*rdum
         pflux_e(i,1,ij) = pflux_e(i,1,ij) + wz1*wt11*rdum
#endif

      enddo

#ifdef _OPENMP
      !$omp critical
      do ij=1,mgrid
         pflux_e(i,0,ij) = pflux_e(i,0,ij) + dnitmp(0,ij)
         pflux_e(i,1,ij) = pflux_e(i,1,ij) + dnitmp(1,ij)
      enddo
      !$omp end critical
#endif
      !$omp end parallel

      ! If we have a particle decomposition on the toroidal domains, do a reduce
      ! operation to add up all the contributions to charge density on the grid
      if(npartdom>1)then
         !$omp parallel do private(ij)
         do ij=1,mgrid
            dnitmp(0,ij)=pflux_e(i,0,ij)
            dnitmp(1,ij)=pflux_e(i,1,ij)
         enddo
         buff=0._wp
         call MPI_ALLREDUCE(dnitmp,buff,(mgrid*2),mpi_Rsize,&
            MPI_SUM,partd_comm,ierror)
         do ij=1,mgrid
            pflux_e(i,0,ij)=buff(0,ij)
            pflux_e(i,1,ij)=buff(1,ij)
         enddo
      endif

      ! poloidal end cell
      do j=0,mpsi
         pflux_e(i,:,igrid(j)+mtheta(j))=pflux_e(i,:,igrid(j)+mtheta(j))+pflux_e(i,:,igrid(j))
         pflux_e(i,:,igrid(j))=pflux_e(i,:,igrid(j)+mtheta(j))
      enddo

      ! toroidal end cell
      sendl=pflux_e(i,0,:)
      recvr=0.0_wp
      icount=mgrid
      idest=left_pe
      isource=right_pe
      isendtag=myrank_toroidal
      irecvtag=isource

      ! send pflux_e to left and receive from right
      call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,recvr,icount,&
         mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

      if(myrank_toroidal == ntoroidal-1) then
         ! B.C. at zeta=2*pi is shifted
         do j=0,mpsi
            ii=igrid(j)
            jt=mtheta(j)
            pflux_e(i,1,ii+1:ii+jt)=pflux_e(i,1,ii+1:ii+jt)+&
               cshift(recvr(ii+1:ii+jt),itran(j))
         enddo
      else
         ! B.C. at zeta<2*pi is continuous
         pflux_e(i,1,:)=pflux_e(i,1,:)+recvr(:)
      endif

   end do   ! do i=1,3

   !  NOTE: pflux_e(:,1,:) has valid data, but no valid data for pflux_e(:,0,:) and pflux_e(:,:,igrid(j))

   ! zero out particle flux in radial boundary cell
   !  do i=0,nbound-1
   !     pflux_e(:,:,igrid(i):igrid(i)+mtheta(i))=&
   !        pflux_e(:,:,igrid(i):igrid(i)+mtheta(i))*float(i)/float(nbound)
   !     pflux_e(:,:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))=&
   !        pflux_e(:,:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))*float(i)/float(nbound)
   !  enddo


   !!!!!!!!!! TEST !!!!!!!!!!!
   dni=0.0_wp
   dne=0.0_wp
   sni=0.0_wp
   sne=0.0_wp

   do i=0,mpsi
      do j=1,mtheta(i)
         ij=igrid(i)+j
         dne(i)=dne(i)+pflux_e(3,1,ij)**2  ! surface averaging
         sne=sne+pflux_e(3,1,ij)**2
         dni(i)=dni(i)+(pflux_i(3,1,ij)*dadtdz(i)/cloud_i)**2 ! surface averaging
         sni=sni+(pflux_i(3,1,ij)*dadtdz(i)/cloud_i)**2
      enddo
   enddo

   call MPI_ALLREDUCE(dne,tdne,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
   tdne=tdne/real(mzetamax*mtheta,wp)   ! average pflux_e^2

   call MPI_ALLREDUCE(dni,tdni,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
   tdni=tdni/real(mzetamax*mtheta,wp)   ! average pflux_e^2

   call MPI_ALLREDUCE(sne,tsne,1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
   call MPI_ALLREDUCE(sni,tsni,1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
   tsne=tsne/real(mzetamax*(mgrid-mpsi-1),wp)
   tsni=tsni/real(mzetamax*(mgrid-mpsi-1),wp)

   if(mype==1) then

      do i=0,mpsi
         write(7,100) i,  tdne(i), tsne, tdni(i), tsni
         100 format(2x, 4hpsi=, I3, 3x, 10h<pflux_e>=,ES14.8, 3x, ES14.8, 3x, 10h<pflux_i>=,ES14.8,  3x, ES14.8)
      end do
      call flush(7)

   end if

   !!!!!!!!!!!END TEST!!!!!!!!!!!!!!!!!!!!!!

   !  ! get Gamma_e/n_0 for using the following approximation later:
   !  !     d/da(Gamma)/n_0=d/da(Gamma/n_0) (consistent with GK ordering)
   !  !     so as to save storage of n_0 at grid points
   !
   !  !$omp parallel do private(ij)
   !    do ij=1,mgrid
   !       pflux_e(1,1,ij)=pflux_e(1,1,ij)*markere(ij)  ! Gamma_e/n_0
   !       pflux_e(2,1,ij)=pflux_e(2,1,ij)*markere(ij)  ! Gamma_e/n_0
   !       pflux_e(3,1,ij)=pflux_e(3,1,ij)*markere(ij)  ! Gamma_e/n_0
   !    enddo

   !$omp parallel do private(i,j,ij)
   do i=0,mpsi
      do j=1,mtheta(i)
         ij=igrid(i)+j
         pflux_e(1,1,ij)=pflux_e(1,1,ij)*cloud_e/dadtdz(i)  ! Gamma_e*J
         pflux_e(2,1,ij)=pflux_e(2,1,ij)*cloud_e/dadtdz(i)  ! Gamma_e*J
         pflux_e(3,1,ij)=pflux_e(3,1,ij)*cloud_e/dadtdz(i)  ! Gamma_e*J
      enddo
   enddo

   ! identical grid points:   igrid(j)=igrid(j)+mtheta(j)
   do j=0,mpsi
      pflux_e(:,:,igrid(j))=pflux_e(:,:,igrid(j)+mtheta(j))
   enddo

   ! NOTE: data on pflux_e(:,0,:) is not valid data, but we need it during gathering process

end subroutine current_eon
