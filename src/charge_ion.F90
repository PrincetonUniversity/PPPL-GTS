subroutine charge_ion
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use field_array
   use diagnosis_array
   use esi_interface
   implicit none

   integer m,i,im,j,k,ip,jt,ii,j11,j10,j01,j00,ierror,larmor,ij,ipjt,&
      icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
   real(wp) weight,wemark,rdum,tdum,delr,delt(0:mpsi),delz,r,wz1,wz0,wp1,&
      wp0,wt11,wt10,wt01,wt00,wdum(0:mpsi),tflr,damping,pi2_inv,&
      psitmp,thetatmp,zetatmp,rhoi,deltheta, &
      psimax,psimin,pdum,total_w,wf,wn,thedum, zdum            !wx
   real(wp) sendl(mgrid),recvr(mgrid)
   real(wp) dnitmp(0:1,mgrid)

   real(wp), external :: SMOOTHSTEP

   integer inner_ip,outter_ip, inner_jt, outter_jt
   real(wp) gyro_phase, intpl_pgyro(larmor_num), intpl_tgyro(larmor_num)

#ifdef _OPENMP
   !! We need the following temporary array only when using OpenMP
   !  real,dimension(:,:),allocatable :: dnitmp
#endif

   call ion_locator_at_larmor_points

   densityi=0.0_wp

#ifdef _OPENMP
   ! The following lines are OpenMP directives for loop-level parallelism
   ! on shared memory machines (see http://www.openmp.org).
   ! When compiling with the OpenMP option (-qsmp=omp for the XLF compiler on
   ! AIX and -mp for the MIPSpro compiler on IRIX), the preprocessor variable
   ! _OPENMP is automatically defined and we use it here to add pieces of codes
   ! needed to avoid contentions between threads. In the following parallel
   ! loop each thread has its own private array "dnitmp()" which is used to
   ! prevent the thread from writing into the shared array "densityi()"
   !
   !$omp parallel private(dnitmp)
   dnitmp=0.   ! Set array elements to zero
#endif
   !$omp do private(m,larmor,weight,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij)
   if(acc_newer_H_zion) then
      !$acc update device(zion(:,:mi))
      acc_newer_H_zion=.false.
      acc_newer_D_zion=.false.
   endif
   !$acc parallel loop gang vector present(zion)
   do m=1,mi
      weight=zion(5,m)

      if (Flag_w_Threshold==1 .and. (weight<w_Thold_min .or. weight>w_Thold_max)) cycle

      wz1=weight*wzion(m)
      wz0=weight-wz1
      !$acc loop seq
      do larmor=1, larmor_num
         wp1=wpion(larmor,m)
         wp0=1.0_wp-wp1

         wt10=wp0*wtion0(larmor,m)
         wt00=wp0-wt10

         wt11=wp1*wtion1(larmor,m)
         wt01=wp1-wt11

#ifdef _OPENMP
         ! Use thread-private temporary array dnitmp to store the results
         ij=jtion0(larmor,m)
         dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt00
         dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt00

         ij=ij+1
         dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt10
         dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt10

         ij=jtion1(larmor,m)
         dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt01
         dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt01

         ij=ij+1
         dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt11
         dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt11
#else
         ! If no loop-level parallelism, use original algorithm (write directly
         ! into array "densityi()".
         ij=jtion0(larmor,m)
         !$acc atomic update
         densityi(0,ij) = densityi(0,ij) + wz0*wt00
         !$acc atomic update
         densityi(1,ij) = densityi(1,ij) + wz1*wt00

         ij=ij+1
         !$acc atomic update
         densityi(0,ij) = densityi(0,ij) + wz0*wt10
         !$acc atomic update
         densityi(1,ij) = densityi(1,ij) + wz1*wt10

         ij=jtion1(larmor,m)
         !$acc atomic update
         densityi(0,ij) = densityi(0,ij) + wz0*wt01
         !$acc atomic update
         densityi(1,ij) = densityi(1,ij) + wz1*wt01

         ij=ij+1
         !$acc atomic update
         densityi(0,ij) = densityi(0,ij) + wz0*wt11
         !$acc atomic update
         densityi(1,ij) = densityi(1,ij) + wz1*wt11
#endif

         ! if(ieee_is_nan(densityi(0,ij)) .or. ieee_is_nan(densityi(1,ij))) then
         !    write(stdout,*) 'charge_ion, mype=',mype,',m=',m,'ij=',ij,'weight=',weight,'densityi(:,ij)=',densityi(:,ij)
         !    call flush(stdout)
         ! endif

      enddo
   enddo



#ifdef _OPENMP
   ! For OpenMP, we now write the results accumulated in each thread-private
   ! array dnitmp() back into the shared array densityi(). The loop is enclosed
   ! in a critical section so that one thread at a time updates densityi().
   !$omp critical
   do ij=1,mgrid
      densityi(0,ij) = densityi(0,ij) + dnitmp(0,ij)
      densityi(1,ij) = densityi(1,ij) + dnitmp(1,ij)
   enddo
   !$omp end critical
#endif
   !$omp end parallel

   ! If we have a particle decomposition on the toroidal domains, do a reduce
   ! operation to add up all the contributions to charge density on the grid
   if(npartdom>1)then
      !$omp parallel do private(ij)
      do ij=1,mgrid
         dnitmp(0,ij)=densityi(0,ij)
         dnitmp(1,ij)=densityi(1,ij)
         densityi(0,ij)=0.
         densityi(1,ij)=0.
      enddo
      call MPI_ALLREDUCE(dnitmp,densityi,(mgrid*2),mpi_Rsize, MPI_SUM,&
         partd_comm,ierror)
   endif

   do ij=1,mgrid
      if(ieee_is_nan(densityi(0,ij)) .or. ieee_is_nan(densityi(1,ij))) then
         write(stdout,*) 'charge_ion2, mype=',mype,',m=',m,'ij=',ij,'weight=',weight,'densityi(:,ij)=',densityi(:,ij)
         call flush(stdout)
      endif
   enddo

   ! poloidal end cell; need valid data for densityi(:,igrid(i))
   do i=0,mpsi
      densityi(:,igrid(i)+mtheta(i))=densityi(:,igrid(i)+mtheta(i))+densityi(:,igrid(i))
      densityi(:,igrid(i))=densityi(:,igrid(i)+mtheta(i))
   enddo

   do ij=1,mgrid
      if(ieee_is_nan(densityi(0,ij)) .or. ieee_is_nan(densityi(1,ij))) then
         write(stdout,*) 'charge_ion3, mype=',mype,',m=',m,'ij=',ij,'weight=',weight,'densityi(:,ij)=',densityi(:,ij)
         call flush(stdout)
      endif
   enddo

   ! toroidal end cell
   sendl=densityi(0,:)
   recvr=0.0_wp
   icount=mgrid
   idest=left_pe
   isource=right_pe
   isendtag=myrank_toroidal
   irecvtag=isource

   ! send densityi to left and receive from right
   call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,recvr,icount,&
      mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

   if(myrank_toroidal == ntoroidal-1)then
      ! B.C. at zeta=2*pi is shifted
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         densityi(1,ii+1:ii+jt)=densityi(1,ii+1:ii+jt)+&
            cshift(recvr(ii+1:ii+jt),itran(i))
         !        densityi(mzeta,ii)=densityi(mzeta,ii+jt)  !wx not necessary
      enddo
   else
      ! B.C. at zeta<2*pi is continuous
      densityi(1,:)=densityi(1,:)+recvr
   endif

   do ij=1,mgrid
      if(ieee_is_nan(densityi(0,ij)) .or. ieee_is_nan(densityi(1,ij))) then
         write(stdout,*) 'charge_ion4, mype=',mype,',m=',m,'ij=',ij,'weight=',weight,'densityi(:,ij)=',densityi(:,ij)
         call flush(stdout)
      endif
   enddo

   if(ZF_SEPARATION==0) then  !  solve total potential

      ! #ifdef _NO_DENSITY_CONSERV

      if(its_jacobi==0) then     ! not enforcing density conservation; use existing, unused parameter its_jacobi 05/2012

         ! if enforcing density conservation is NOT necessary, then use
         ! ************************************************
         !$omp parallel do private(i,j,ij)
         do i=0,mpsi
            do j=1,mtheta(i)
               ij=igrid(i)+j
               densityi(1,ij)=(inv_larmor_num)*densityi(1,ij)*markeri(ij)  ! =dn/n_0
            enddo
            densityi(1,igrid(i))=densityi(1,igrid(i)+mtheta(i)) ! validate densityi(1,igrid(i))
         enddo
         ! ************************************************
         ! #else

      else    ! enforcing density conservation

         ! if enforcing density conservation, then use
         ! ++++++++++++++++++++++++++++++++++++++++++++++++

         !!!!!! ATTENTION: the following do when using cray compiler with default option may produce rdum=NaN !!!!!!!!

         !    rdum=0.0_wp
         !  !$omp parallel private(i,j,ij) reduction(+:rdum)
         !  !$omp do
         !    do i=0,mpsi
         !       do j=1,mtheta(i)
         !          ij=igrid(i)+j
         !          rdum=rdum+(inv_larmor_num)*densityi(1,ij)
         !          densityi(1,ij)=(inv_larmor_num)*densityi(1,ij)*markeri(ij)  ! =dn/n_0
         !       enddo
         !    enddo
         !  !$omp end parallel

         ! replaced by:

         wdum=0.0_wp
         !$omp parallel private(i,j,ij) reduction(+:wdum)
         !$omp do
         do i=0,mpsi
            do j=1,mtheta(i)
               ij=j+igrid(i)
               wdum(i)=wdum(i)+(inv_larmor_num)*densityi(1,ij)
               densityi(1,ij)=(inv_larmor_num)*densityi(1,ij)*markeri(ij)  ! =dn/n_0
            enddo
         enddo
         !$omp end parallel

         rdum=sum(wdum)

         total_w=0.0_wp
         call MPI_ALLREDUCE(rdum,total_w,1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
         ! total_w=total particle weight (appro to int(delta-n d^3x) )

         ! select weight function for maintaining charge neutrality, see note
         !  wdum=fluxpsi   ! in terms turbulence intensity I(a)

         !  wdum(1:mpsi-1)=sqrt(fluxpsi(1:mpsi-1))   ! in terms sqrt[I(a)] -- mostly used
         !  wdum(0)=0.0_wp
         !  wdum(mpsi)=0.0_wp

         !  wdum(1:mpsi-1)=sqrt(En_tur(1:mpsi-1))   ! in terms of sqrt[turbulence energy]
         !  wdum(0)=0.0
         !  wdum(mpsi)=0.0

         wdum=eta_grid**2 ! in terms eta_i

         !  wdum=1.0  ! uniform compensation used previously: may generate artificial ZF
         ! in linearly stable zones!!!

         ! calculate int[dVda*I(a)*da]
         wf=0.0_wp
         do i=1,mpsi
            wf=wf+0.5_wp*(psi_grid(i)-psi_grid(i-1))*(dVda(i)*wdum(i)+dVda(i-1)*wdum(i-1)) ! =int[I(a)d^3x]
         end do

         !$omp parallel do private(i,j,ij)
         do i=0,mpsi
            do j=1,mtheta(i)
               ij=igrid(i)+j
               densityi(1,ij)=densityi(1,ij)-(density_ave/den_grid(i))*&
                  (total_w/misum)*(wdum(i)*finval_v/wf)  ! enforce density conservation int(delta-n d^3x)=0 (see note)
            enddo
            densityi(1,igrid(i))=densityi(1,igrid(i)+mtheta(i))  ! validate densityi(1,igrid(i)
         enddo
         ! NOTE: enforcing density conservation modifies only <delta-n>, and no change for (delta-n)-<delta-n>
         ! ++++++++++++++++++++++++++++++++++++++++++++++++
         ! #endif

      endif    ! end if(its_jacobi==0)

      !============================================================
      ! MGYOO (01/08/2019): To diagnose zonal component for later snapshot,
      ! This zonal component is only used in the snapshot routine to separate zonal and non-zonal components
      !============================================================
      if(mod(istep,msnap) .eq. 0) then
         ! flux surface average and normalization
         wdum=0.0
         !$omp parallel private(i,j,ij) reduction(+:wdum)
         !$omp do
         do i=0,mpsi
            do j=1,mtheta(i) ! j=0 is ignored to avoid over-estimation at j=0 or j=mtheta(i) which is theta=0 or 2*pi
               ij=igrid(i)+j
               wdum(i)=wdum(i)+dtemper(ij)*densityi(1,ij)  ! new flux-surface avarage --> zonali; used in Jan. 2015
            enddo
         enddo
         !$omp end parallel
         call MPI_ALLREDUCE(wdum,zonali,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
         zonali=zonali/jaco_sum    ! =dn_00/n_0; new flux-surface avarage --> zonali; used in Jan. 2015
      endif


      !====================================================================
      ! MGYOO(12/03/2018) : To measure FLA(Field Line Average) component
      ! This FLA might be important in the case of external magnetic perturbation
      ! MGYOO(12/13/2018) : tot_dni is calculated rather than (dh/f or df/f)
      if((myrank_partd==nproc_partd-1) .and. (mod(istep,msnap) .eq. 0)) then
         do i=0,mpsi
            do j=0,mtheta(i)
               ij = igrid(i)+j
               tot_dni(:,ij) = densityi(:,ij)*den_grid(i)  ! Here, (delta-f) = (densityi*den_grid), because densityi=(delta-f)/n0
            enddo
         enddo
         call Do_FieldLineAverage_by_using_gPPI(tot_dni, FLA_dni)
      endif
      !===================================================================

      ! zero out charge in boundary cell
      !  do i=0,nbound-1
      !     densityi(:,igrid(i):igrid(i)+mtheta(i))=   &
      !        densityi(:,igrid(i):igrid(i)+mtheta(i))*real(i,wp)/real(nbound,wp)
      !     densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))=   &
      !        densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))*real(i,wp)/real(nbound,wp)
      !  enddo
      do i=0,nbound
         rdum=smoothstep(0.0_wp, real(nbound,wp),real(i,wp))

         ! pdum=abs(dni_grid(nbound)/dni_grid(i))
         ! zdum=abs(dni_grid(mpsi-nbound)/dni_grid(mpsi-i))
         pdum=1.0_wp
         zdum=1.0_wp

         densityi(:,igrid(i):igrid(i)+mtheta(i))=   &
            densityi(:,igrid(i):igrid(i)+mtheta(i))*rdum*pdum
         densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))=   &
            densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))*rdum*zdum
      enddo

   else  ! solve zonal flows separately when ZF_SEPARATION==1

      ! zero out charge in boundary cell
      !  do i=0,nbound-1
      !     densityi(:,igrid(i):igrid(i)+mtheta(i))=   &
      !        densityi(:,igrid(i):igrid(i)+mtheta(i))*real(i,wp)/real(nbound,wp)
      !     densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))=   &
      !        densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))*real(i,wp)/real(nbound,wp)
      !  enddo

      ! flux surface average and normalization
      wdum=0.0
      !$omp parallel private(i,j,ij) reduction(+:wdum)
      !$omp do
      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid(i)+j
            !  wdum(i)=wdum(i)+(inv_larmor_num)*densityi(1,ij)      ! original flux-surface avarage --> zonali; replaced in Jan. 2015
            densityi(1,ij)=(inv_larmor_num)*densityi(1,ij)*markeri(ij)  ! =dn/n_0
            wdum(i)=wdum(i)+dtemper(ij)*densityi(1,ij)  ! new flux-surface avarage --> zonali; used in Jan. 2015
         enddo
         densityi(1,igrid(i))=densityi(1,igrid(i)+mtheta(i))
      enddo
      !$omp end parallel

      ! global sum of phi00, broadcast to every PE
      call MPI_ALLREDUCE(wdum,zonali,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      !  zonali=zonali*pmarki   ! =dn_00/n_0  -- original flux-surface avarage --> zonali; replaced in Jan. 2015
      !  total_w=sum(wdum)      ! total particle weight (appro to int(delta-n d^3x) )
      zonali=zonali/jaco_sum    ! =dn_00/n_0; new flux-surface avarage --> zonali; used in Jan. 2015

      do i=0,mpsi
         if(ieee_is_nan(zonali(i))) then
            write(stdout,*) 'charge_ion,hello0, mype=',mype,',zonali=',zonali
            call flush(stdout)
            exit
         endif
      enddo

      !====================================================================
      ! MGYOO(01/08/2019) : To measure FLA(Field Line Average) component
      ! This FLA might be important in the case of external magnetic perturbation
      ! MGYOO(01/08/2019) : tot_dni is calculated rather than (dh/f0 or df/f0)
      if((myrank_partd==nproc_partd-1) .and. (mod(istep,msnap) .eq. 0)) then
         do i=0,mpsi
            do j=0,mtheta(i)
               ij = igrid(i)+j
               tot_dni(:,ij) = densityi(:,ij)*den_grid(i)  ! Here, (delta-f) = (densityi*den_grid), because densityi=(delta-f)/n0
            enddo
         enddo
         call Do_FieldLineAverage_by_using_gPPI(tot_dni, FLA_dni)
      endif
      !===================================================================


      ! densityi subtracted (0,0) mode
      !$omp parallel do private(i,j,ij)
      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid(i)+j
            densityi(1,ij)=densityi(1,ij)-zonali(i)  ! tilde{dn}/n_0
         enddo
         densityi(1,igrid(i))=densityi(1,igrid(i)+mtheta(i)) !wx validate densityi(1,igrid(i))
      enddo

      if(its_jacobi==0) then

         do i=0,nbound
            rdum=smoothstep(0.0_wp, real(nbound,wp),real(i,wp))

            densityi(:,igrid(i):igrid(i)+mtheta(i))=   &
               densityi(:,igrid(i):igrid(i)+mtheta(i))*rdum
            densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))=   &
               densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))*rdum

            zonali(i)=zonali(i)*rdum
            zonali(mpsi-i)=zonali(mpsi-i)*rdum
         enddo

      else     !  its_jacobi=1 --> enforcing enforce charge conservation for zonal flow mode

         ! enforce charge conservation for zonal flow mode
         ! wx note: why should charge conservation be enforced? Should it be automatically satisfied ???
         ! int(d^3x A)=int(<A>dV)=int(<A>dV/dpsi *dpsi; dV/dpsi=int(J dtheta*dzeta);
         ! dzeta are even; allnodes(i)=sum_jk(J_ijk*dpsi_i*dtheta_j);
         ! assuming n_0 is uniform
         rdum=0.0
         tdum=0.0
         do i=1,mpsi-1
            rdum=rdum+allnodes(i)                           !wx
            tdum=tdum+allnodes(i)*zonali(i)                 !wx
         enddo
         tdum=tdum/rdum
         ddeni=tdum !for diagnostic, recorded only in history.out at diagnosis.f90

         ! select weight function for maintaining charge neutrality, see note
         !  wdum=fluxpsi   ! in terms turbulence intensity I(a)

         !   wdum(1:mpsi-1)=sqrt(fluxpsi(1:mpsi-1))   ! in terms sqrt[I(a)]
         !   wdum(0)=0.0_wp
         !   wdum(mpsi)=0.0_wp


         !  added Jan. 2015
         if(irun==0.and.istep < 2*ndiag) then
            wdum=1.0  ! this is because initially sumfluxpsi=0 which cannot be used as weighting function
         else
            wdum=sqrt(sumfluxpsi)   ! in terms of sqrt{int[I(a)dt]}
         end if


         !  wdum(1:mpsi-1)=sqrt(En_tur(1:mpsi-1))   ! in terms of sqrt[turbulence energy]
         !  wdum(0)=0.0
         !  wdum(mpsi)=0.0

         !  wdum=eta_grid**2 ! in terms eta_i

         !  wdum=1.0  ! uniform compensation used previously: may generate artificial ZF
         ! in linearly stable zones!!!

         !  wdum=eta_grid  ! check ZF driving turbulence

         do i=0,nbound
            rdum=smoothstep(0.0_wp, real(nbound,wp),real(i,wp))
            wdum(i)=wdum(i)*rdum
            wdum(mpsi-i)=wdum(mpsi-i)*rdum
         enddo


         ! calculate int[dVda*I(a)*da]
         wf=0.0_wp
         wn=0.0_wp
         do i=1,mpsi

            ! if(order==1) then  ! FEM solver
            !    wf=wf+0.5*(psi_grid(i)-psi_grid(i-1))*(dVda(i)*wdum(i)+dVda(i-1)*wdum(i-1))
            ! else               !  old solver

            if(order==1) then  ! FEM solver, added Jan. 2015 See note -- charge conservation consistent with E_zf=0 at boundary
               wf=wf+0.5*(psi_grid(i)-psi_grid(i-1))*(wdum(i)+wdum(i-1))
               wn=wn+0.5_wp*(psi_grid(i)-psi_grid(i-1))*(den_grid(i-1)*dVda(i-1)*zonali(i-1)  &
                  +den_grid(i)*dVda(i)*zonali(i))
            else               !  old solver
               wf=wf+0.5_wp*(psi_grid(i)-psi_grid(i-1))*(dVda(i)*wdum(i)*B2_s(i)  &
                  +dVda(i-1)*wdum(i-1)*B2_s(i-1))
               wn=wn+0.5_wp*(psi_grid(i)-psi_grid(i-1))*(dVda(i)*zonali(i)*B2_s(i)  &
                  +dVda(i-1)*zonali(i-1)*B2_s(i-1))
            end if
         end do

         do i=1,mpsi-1

            !     zonali(i)=zonali(i)-tdum
            !   Q: may change ZF behavior in linear stable zone--> drive a ZF in stable zone ????
            !     zonali(i)=zonali(i)-(density_ave/den_grid(i))*(total_w/misum) ! see note

            ! if(order==1) then     ! FEM solver
            !    zonali(i)=zonali(i)-(density_ave/den_grid(i))*(total_w/misum)  &
            !              *(wdum(i)*finval_v/wf)            ! see note
            ! else                  !  old solver

            if(order==1) then     ! FEM solver, added Jan. 2015 See note -- charge conservation consistent with E_zf=0 at boundary
               !  zonali(i)=zonali(i)-wdum(i)*(wn/wf)/(den_grid(i)*dVda(i))
               wdum(i)=wdum(i)*(wn/wf)/(den_grid(i)*dVda(i))
            else                  !  old solver
               !  zonali(i)=zonali(i)-wdum(i)*(wn/wf)            ! see note
               wdum(i)=wdum(i)*(wn/wf)            ! see note
            endif

         enddo

         if(mype==0.and.idiag==0) then
            wdum(0)=0.0
            wdum(mpsi)=0.0
            write(800,100) istep, ( wdum(i),i=0,mpsi), ( zonali(i),i=0,mpsi), ( sumfluxpsi(i),i=0,mpsi)
            100  format(1x,i6, 606(1x,d12.6)) ! for maxium mpsi=200
            call flush(800)
         end if

         zonali(1:mpsi-1)=zonali(1:mpsi-1)-wdum(1:mpsi-1)

         do i=0,nbound
            rdum=smoothstep(0.0_wp, real(nbound,wp),real(i,wp))

            densityi(:,igrid(i):igrid(i)+mtheta(i))=   &
               densityi(:,igrid(i):igrid(i)+mtheta(i))*rdum
            densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))=   &
               densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))*rdum

            zonali(i)=zonali(i)*rdum
            zonali(mpsi-i)=zonali(mpsi-i)*rdum
         enddo

      end if   ! end (its_jacobi==1)

      do i=0,mpsi
         if(ieee_is_nan(zonali(i))) then
            write(stdout,*) 'charge_ion,hello1, mype=',mype,',zonali=',zonali
            call flush(stdout)
            exit
         endif
      enddo

   end if

   ! ----- add artificial zonal flow damping to identify momentum pinch and residual stress 09/2008 ---------
   !
   !                                     zonali=0.0_wp
   !                                   zonali=zonali/10.0_wp
   !
   ! --------------------------------------------------------------------------------------------------------


   ! NOTE: data on densityi(0,:) is not valid data  --- 2005

   ! if(mype==0) then
   ! do i=1,mgrid
   ! write(889,*) i, densityi(1,i)
   ! end do
   ! end if

end subroutine charge_ion
