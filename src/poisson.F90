subroutine poisson( solver, masssolver )
   use precision
   !2015-oct-24 already defined in "solver_data_class": use mpi_defs
   use global_parameters
   use particle_decomp
   use field_array
   use esi_interface
   use particle_array
   use diagnosis_array
   !use linear_system
   use solver_data_class
   use Structured_Grid
   use MGY_PETSc_Interface_Data
   implicit none
   !
   type(solver_data) :: solver, masssolver
   ! local vars
   integer i,ii,jt,it,ij,j,jj,k,iteration,mring,mindex,ierr,it2,loop
   integer istatus(MPI_STATUS_SIZE)
   !  integer,dimension(:,:),allocatable :: nindex
   !  integer,dimension(:,:,:),allocatable :: indexp
   !  real,dimension(:,:,:),allocatable :: ring

   real(wp) tmp,rms ! gamma0,

   real(wp), allocatable :: ori_SG_phi(:,:)
   real(wp) RMS_deviations, RMS_old_phi
   real(wp) adum(0:mpsi),phi_flux_ave(0:mpsi),adum0(0:mpsi),adum_11(0:mpsi)
   real(wp) old_phip00(0:mpsi), phip00_convergency
   !real(wp) ptilde(mgrid),dentmp(mzeta,mgrid),titer(mzeta,mgrid)
   real(wp) dentmp(mgrid),xx(mgrid) !,titer(mgrid)
   real(wp) add_phi(mgrid), yy(mgrid), bb2(mgrid)
   real(wp) recvl(mgrid)
   real(wp), external :: SMOOTHSTEP
   real(wp) instep,outstep

   real(wp) zdum,tdum,tdum0,tdum1
   real(wp), external :: temperature_e
   real(wp) tmp_Te0_inverse

   real(wp) dpsi_inv,e_psi2,e_psi1,dpsi2_inv,dpsi1,dpsi2,dpsi1_inv
   real(wp) inv_ds
   real(wp), external :: dsda

   real(doubleprec),save :: t0=0.0_wp
   real(doubleprec) t1,rtc

   !gamma0=0.8: max. resolution for k=0.5. for k=[0,0.8], err=0.08%
   ! gamma0=0.75
   ! iteration=10
   ! gamma0=gama0
   ! iteration=its

   !  zero gyroradius limit
   !  iteration=1        !wx

   ! initialize poisson solver
   if(istep==0)then
      if(irun==0) phi=0.0_wp     ! for follow-up runs, phi is read from restart
      if(irun==0) phip00=0.0_wp  ! for follow-up runs, phi is read from restart
      !gamma=gamma0

      !  initialize terms for solving <phi> and delta-phi interatively using FEM solver -- see notes
      if( solver%order > 0 ) then
         allocate ( div1(mgrid),div2(mgrid) )
         div1(:)=Jcbn(:)*gsij(:)*gij(1,:)
         div2(:)=Jcbn(:)*gsij(:)*gij(4,:)

         call find_dxdadxdt( masssolver, mgrid, div1, xx, yy, ierr )
         div1=xx/Jcbn

         call find_dxdadxdt( masssolver, mgrid, div2, xx, yy, ierr )
         div2=yy/Jcbn

         do i=0, mpsi
            do j=0,mtheta(i)
               ij=igrid(i)+j
               div1(ij)=(div1(ij)+div2(ij))/den_grid(i)   ! = -beta*div[sum(gs)*grad(a)]
            end do
            div1(igrid(i))=div1(igrid(i)+mtheta(i))
         enddo

         div2(:)=gsij(:)*gij(1,:)
         do i=0, mpsi
            do j=0,mtheta(i)
               ij=igrid(i)+j
               div2(ij)=div2(ij)/den_grid(i)   ! = -beta*[sum(gs)*gaa)]
            end do
            div2(igrid(i))=div2(igrid(i)+mtheta(i))
         enddo
      endif
   endif

   !  tmp=1.0/(tite+1.0-gamma) --> replaced by titer(ij)
   !  zero gyroradius limit
   !  tmp=1.0/tite           !wx

   ! Note added Feb. 2012 -- control parameter for electron model:
   !      ifelectron=0 --> adiabatic electrons
   !      ifelectron=1 --> full delta-f electrons
   !      ifelectron=2 --> delta-h electrons (non-adiabatic electrons)
   old_phi=phi
   if(ifelectron==1) then  ! full delta-f electron case
      ! this part works for both ZF_SEPARATION=0 and 1

      ! for ZF_SEPARATION=0, namely solving total potential:
      ! densityi=n_i/n_0 and densitye=n_e/n_0 are total density fluctuations
      ! with n_e being the full (ad. + non-ad.) electron perturbation

      ! for ZF_SEPARATION=1, namely solving ZF separately
      ! densityi=(n_i-<n_i>)/n_0 and densitye=(n_e-<n_e>)/n_0
      ! with <n_i> and <n_e> are flux-surface-averaged fluctuations

      !$omp parallel do private(i,j,ij)
      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid(i)+j
            dentmp(ij)=densityi(1,ij)-densitye(1,ij)
            !... finite element needs not this modification
            if( solver%order <= 0 ) then
               dentmp(ij)=dentmp(ij)*vth_grid(i)**2
            else
               if(ZF_SEPARATION==1) dentmp(ij)=dentmp(ij)+(zonali(i)-zonale(i))  ! total density fluctuations used for FEM solver with ZF_SEPARATION=1
            endif
         enddo
         dentmp(igrid(i))=dentmp(igrid(i)+mtheta(i))  ! which is missing in last do-loop
      enddo


      if(Flag_Poisson3D_on_SG==0) then
         do i=0,mpsi
            tmp_Te0_inverse = 1.0/temperature_e(psi_grid(i),0.0_wp)
            do j=1,mtheta(i)
               ij=igrid(i)+j
               if(adiab > 0.0_wp) then
                  dentmp(ij) = dentmp(ij) + adiab * abs(csc0) * m0ms * tmp_Te0_inverse * phi(1, ij)
               else
                  dentmp(ij) = dentmp(ij) + 0.01/(Cs*Cs) * phi(1,ij)
               endif
            enddo
            dentmp(igrid(i))=dentmp(igrid(i)+mtheta(i))  ! which is missing in last do-loop
         enddo
      else
         ! This adding (artificial constant or adiabatic-like) part will be conducted in Solve_test_Poisson3D_onSG routine in mgy_pestc.F90
         ! In the case of the comparison test case (when Flag_Poisson3D_on_SG=2), Poisson3D_on_SG will also add an artificial term.
         ! However, for a realistic delat-f simulation using SG (Structured Grid) will not use any adding artificial diagonal term.
      endif

      if(Flag_Poisson3D_on_SG>=1) then

         ! Calculate dn
         if(Flag_Weight_dfh==0) then
            SG%dne_SM = SG%dne(:,:,1)
            SG%dni_SM = SG%dni(:,:,1)
         elseif(Flag_Weight_dfh==1) then
            do i=1,SG%dim_s
               SG%dne_SM(:,i) = SG%ne_by_fh(:,i,1) -SG%n0_1D(i)
               SG%dni_SM(:,i) = SG%ni_by_ff(:,i,1) -SG%n0_1D(i)
            enddo
         endif

         ! Then, smooth them
         call Smooth_SG_Data_2D([2,0],SG%dne_SM,SG%Cell_dV)
         call Smooth_SG_Data_2D([2,0],SG%dni_SM,SG%Cell_dV)

         ! Solve Poisson3D Equation on SG using SG%dni and SG%dne
         call Solve_Poisson3D_on_SG
         call Cal_gradPhi_on_SG

         ! Note that potential on SG (StructuredGrid) uses SI unit (Volt)
         ! Therefore, a unit conversion is needed afther the interpolation
         ! call Interpolate_SG_into_GTSmgrid(SG%test_phi_2D,phi(1,:))
         call Interpolate_SG_into_GTSmgrid(SG%phi(:,:,1),phi(1,:))
         phi(1,:) = ((charg0*1e8)/(E_00*c))*phi(1,:) ! Convert from V unit to GTS normalized unit
         ! Note that phi(0,:) will be communicated in smooth_phi or field subroutine..
         ! Also, phi00 will be calculated in smooth_phi in the case of ZE_SEPRATION=0

         ! Send phi(1,:) to right poloidal plane
         call MPI_SENDRECV(phi(1,:),mgrid,mpi_Rsize,right_pe,0,recvl,mgrid,mpi_Rsize,left_pe,0,toroidal_comm,istatus,ierr)
         if(myrank_toroidal==0)then
            do i=0,mpsi
               ii=igrid(i)
               jt=mtheta(i)
               phi(0,ii+1:ii+jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
            enddo
         else
            phi(0,:)=recvl
         endif

         ! poloidal BC
         do i=0,mpsi
            phi(:,igrid(i)) = phi(:,igrid(i)+mtheta(i))
         enddo

         ! flux surface average and normalization
         adum=0.0_wp
         do i=0,mpsi
            do j=1,mtheta(i)  !j=0 is ignored to avoid over-estimation at j=0 or j=mtheta(i) which is theta=0 or 2*pi
               ij=igrid(i)+j
               adum(i)=adum(i)+dtemper(ij)*phi(1,ij)   ! surface averaging
            enddo
         enddo
         phi00=0.0_wp
         call MPI_ALLREDUCE(adum,phi00,mpsi+1,mpi_Rsize, MPI_SUM,toroidal_comm,ierr)
         phi00=phi00/jaco_sum     !phi00=flux-surface avarage

         ! direct numerical derivative to get d<phi>/da
         inv_ds = real(mpsi)/(s1-s0)
         do i=1,mpsi-1
            phip00(i) = dsda(psi_grid(i))*0.5*inv_ds*(phi00(i+1)-phi00(i-1))
         enddo
         ! 3-point forward (or backward) difference at the boundaries
         phip00(0) = dsda(psi_grid(0))*0.5*inv_ds*(-3.0*phi00(0)+4.0*phi00(1)-phi00(2))
         phip00(mpsi) = dsda(psi_grid(mpsi))*0.5*(-inv_ds)*(-3.0*phi00(mpsi)+4.0*phi00(mpsi-1)-phi00(mpsi-2))

         if(ZF_SEPARATION==1) then
            ! Subtract zonal_phi so get nonzonal phi
            do i=0,mpsi
               do j=0,mtheta(i)
                  ij=igrid(i)+j
                  phi(:,ij) = phi(:,ij) - phi00(i)
               enddo
            enddo
         endif

         return
      endif


      if(ZF_SEPARATION==1) then

         if(mype==0 .and. mod(istep,msnap)==0 ) open(unit=191209,file='phip00_its.dat',position='append',action='write')
         phip00_convergency=1.0 ! any initial value to go inside do while
         it=0

         do while(phip00_convergency>its_phip00_tolerance .and. it<its)
            it=it+1

            call solve( solver, masssolver, mgrid, dentmp, xx, ierr )  !  --> delta-phi for given phip00=d<phi>/da
            do i=0, mpsi
               xx(igrid(i))=xx(igrid(i)+mtheta(i))  !validate xx(1,igrid(i)) for safe
            enddo

            call find_dxdadxdt( masssolver, mgrid, xx, add_phi, yy, ierr )

            add_phi(:)=(add_phi(:)*gij(1,:)+yy(:)*gij(4,:))*gsij(:)   ! --> sum(gs)*grad_perp(delta-phi).grad(a)
            adum=0.0_wp
            do i=0,mpsi
               do j=1,mtheta(i)
                  ij=igrid(i)+j
                  adum(i)=adum(i)+dtemper(ij)*add_phi(ij)   ! surface averaging
               enddo
            enddo

            phi_flux_ave=0.0_wp
            call MPI_ALLREDUCE(adum,phi_flux_ave,mpsi+1,mpi_Rsize, MPI_SUM,toroidal_comm,ierr)
            phi_flux_ave=phi_flux_ave/jaco_sum    !  =<sum(gs)*grad_perp(delta-phi).grad(a)>

            old_phip00 = phip00

            phip00=0.0_wp
            do i=1,mpsi
               phip00(i)=phip00(i-1)+0.5_wp*(psi_grid(i)-psi_grid(i-1))*  &
                  (-den_grid(i-1)*dVda(i-1)*(zonali(i-1)-zonale(i-1))  &
                  -den_grid(i)*dVda(i)*(zonali(i)-zonale(i)) )
            enddo

            phip00=phip00/dVda/gsgaa
            phip00=phip00 -(phi_flux_ave*dVda-phi_flux_ave(0)*dVda(0))/dVda/gsgaa   ! --> phip00=dphi_00/da with phip00(0)=0

            phip00=phip00+Er_a0        ! --> phip00=dphi_00/da with phip00(0)=Er_a0 which is Er_NC when ifneoclassical=1

            phip00_convergency = sum(sqrt((phip00 - old_phip00)**2.0))/sum(sqrt(old_phip00**2.0))
            if(mype==0 .and. mod(istep,msnap)==0 ) write(191209,"(2(I8,1x),200(ES16.8,1x))") istep, it, phip00_convergency, phip00(0:mpsi)
         enddo
         if(mype==0 .and. mod(istep,msnap)==0 ) then
            call flush(191209)
            close(191209)
         endif

         phi(1,:) = xx

         ! (0,0) mode potential store in phi00
         phi00=0.0_wp   ! phi00(0)=0.0 is chosen
         do i=1,mpsi
            phi00(i)=phi00(i-1)+0.5_wp*(psi_grid(i)-psi_grid(i-1))*(phip00(i-1)+phip00(i))
         enddo

      else  ! ZF_SEPARATION=0 --> solve total potential

         ! TODO : convergence check is not implemented yet.
         adum0=phi00
         do it=1,its
            !   calculate <phi> using phi of last time step for the 1st iteration

            ! if(mype==-100) then      !  keep the following coding in the record
            adum=0.0_wp
            rms=0.0_wp
            phi_flux_ave=0.0_wp
            !$omp parallel private(i,j,ij) reduction(+:adum,rms)
            !$omp do
            do i=0,mpsi
               do j=1,mtheta(i)
                  ij=igrid(i)+j
                  adum(i)=adum(i)+dtemper(ij)*phi(1,ij)  ! surface averaging
                  rms=rms+phi(1,ij)**2
               enddo
            enddo
            !$omp end parallel
            call MPI_ALLREDUCE(adum,phi_flux_ave,mpsi+1,mpi_Rsize, MPI_SUM,toroidal_comm,ierr)
            phi_flux_ave=phi_flux_ave/jaco_sum    !  =<phi>

            rtc=phi_flux_ave(0)

            ! re-level phi and <phi> by subtracting a constant for original integral solver
            if(solver%order <= 0 )  phi_flux_ave=phi_flux_ave - phi_flux_ave(0)

            t1=sum((phi_flux_ave-adum0)**2)/sum(adum0**2)  ! calculated for checking convergence
            adum0=phi_flux_ave

            phi(1,:)=phi(1,:)-phi_flux_ave(0)
            if(solver%order <= 0 ) then   ! use original integral solver
               phi_flux_ave=phi_flux_ave*tite        !  =<phi>*(Ti/Te)
            else                          ! use FEM solver
               phi_flux_ave=phi_flux_ave*abs(csc0)/Te_grid   ! =(e/e0)*(m0/ms)*<phi>/T_e with <phi>, T_e normalized
            endif

            if(mype==0)   print*, 'mype=', mype, 'its=', it, '  <phi>(',ipsi,')=',  phi_flux_ave(ipsi), 'rms=', rms,  '  <phi>(0)=', rtc, &
               'relative error in <phi>^2=', t1

            add_phi=0.0_wp
            xx=0.0_wp
            do i=0,mpsi
               do j=1,mtheta(i)
                  ij=igrid(i)+j
                  add_phi(ij)=dentmp(ij) + phi_flux_ave(i)
               end do
               add_phi(igrid(i))=add_phi(igrid(i)+mtheta(i))
            end do

            call solve( solver, masssolver, mgrid, add_phi, xx, ierr )
            phi(1,:) = xx

            do i=0, mpsi
               phi(1,igrid(i))=phi(1,igrid(i)+mtheta(i))  !validate phi(1,igrid(i)) for safe
            enddo

         enddo  ! end of iteration do it=1,its
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   else     ! for both delta-h (ifelectron=2) and adiabatic (ifelectron=0) electrons

      ! for ZF_SEPARATION=0, namely solving total potential:
      ! densityi=n_i/n_0 is total ion fluctuation
      ! densitye=n_h/n_0 is non-adiabatic electron density fluctuation

      ! for ZF_SEPARATION=1, namely solving ZF separately
      ! densityi=(n_i-<n_i>)/n_0 and densitye=(n_h-<n_h>)/n_0
      ! with <n_i> and <n_h> are flux-surface-averaged fluctuations

      !$omp parallel do private(i,j,ij)
      do i=0, mpsi
         do j=1,mtheta(i)
            ij=igrid(i)+j
            dentmp(ij)=densityi(1,ij)
            if(ifelectron==2) dentmp(ij)=dentmp(ij)-densitye(1,ij)   ! with non-adiabatic electrons
            if( solver%order <= 0 )  then
               dentmp(ij)=dentmp(ij)*vth_grid(i)**2
            else
               if(ZF_SEPARATION==1) dentmp(ij)=dentmp(ij)+(zonali(i)-zonale(i))  ! total density fluctuations used for FEM solver with ZF_SEPARATION=1
            endif
         enddo
      enddo

      do i=0, mpsi
         dentmp(igrid(i))=dentmp(igrid(i)+mtheta(i))  ! which is missing in last do-loop
      enddo

      if(Flag_Poisson3D_on_SG==0) then
         if(Flag_Def_dne_adiab==1) then
            do i=0, mpsi
               do j=0,mtheta(i)
                  ij=igrid(i)+j
                  dentmp(ij)=dentmp(ij) - abs(csc0)*phi00(i)/Te_grid(i)
                enddo
            enddo
         endif
      else
          if(Flag_Def_dne_adiab==0) then
            do i=0, mpsi
               do j=0,mtheta(i)
                  ij=igrid(i)+j
                  dentmp(ij)=dentmp(ij) + abs(csc0)*phi00(i)/Te_grid(i)
                enddo
            enddo
         endif
      endif


      if(Flag_Poisson3D_on_SG>=1) then

         ! Calculate dn
         if(Flag_Weight_dfh==0) then
            SG%dne_SM = SG%dne(:,:,1)
            SG%dni_SM = SG%dni(:,:,1)
         elseif(Flag_Weight_dfh==1) then
            do i=1,SG%dim_s
               SG%dne_SM(:,i) = SG%ne_by_fh(:,i,1) -SG%n0_1D(i)
               SG%dni_SM(:,i) = SG%ni_by_ff(:,i,1) -SG%n0_1D(i)
            enddo
         endif

         ! and smooth them
         call Smooth_SG_Data_2D([2,0],SG%dne_SM,SG%Cell_dV)
         call Smooth_SG_Data_2D([2,0],SG%dni_SM,SG%Cell_dV)

         if(Flag_Def_dne_adiab==3) then
            call Solve_Nonlinear_Poisson3D_on_SG()
         else
            call Solve_Poisson3D_on_SG
         endif


         ! call Measure_Electron_Density_on_SG
         ! call Measure_Ion_Density_on_SG
         ! ! call Smooth_SG_Data_2D([3,0],SG%dni(:,:,1),SG%Cell_dV)
         ! ! call Smooth_SG_Data_2D([3,0],SG%dne(:,:,1),SG%Cell_dV)
         ! call Solve_Poisson3D_on_SG
         call Cal_gradPhi_on_SG
         ! Note that potential on SG (StructuredGrid) uses SI unit (Volt)
         ! Therefore, a unit conversion is needed afther the interpolation
         ! call Interpolate_SG_into_GTSmgrid(SG%test_phi_2D,phi(1,:))
         call Interpolate_SG_into_GTSmgrid(SG%phi(:,:,1),phi(1,:))
         phi(1,:) = ((charg0*1e8)/(E_00*c))*phi(1,:)

         ! Send phi(1,:) to right poloidal plane
         call MPI_SENDRECV(phi(1,:),mgrid,mpi_Rsize,right_pe,0,recvl,mgrid,mpi_Rsize,left_pe,0,toroidal_comm,istatus,ierr)
         if(myrank_toroidal==0)then
            do i=0,mpsi
               ii=igrid(i)
               jt=mtheta(i)
               phi(0,ii+1:ii+jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
            enddo
         else
            phi(0,:)=recvl
         endif

         ! poloidal BC
         do i=0,mpsi
            phi(:,igrid(i)) = phi(:,igrid(i)+mtheta(i))
         enddo

         ! flux surface average and normalization
         adum=0.0_wp
         do i=0,mpsi
            do j=1,mtheta(i)  !j=0 is ignored to avoid over-estimation at j=0 or j=mtheta(i) which is theta=0 or 2*pi
               ij=igrid(i)+j
               adum(i)=adum(i)+dtemper(ij)*phi(1,ij)   ! surface averaging
            enddo
         enddo
         phi00=0.0_wp
         call MPI_ALLREDUCE(adum,phi00,mpsi+1,mpi_Rsize, MPI_SUM,toroidal_comm,ierr)
         phi00=phi00/jaco_sum     !phi00=flux-surface avarage

         ! direct numerical derivative to get d<phi>/da
         inv_ds = real(mpsi)/(s1-s0)
         do i=1,mpsi-1
            phip00(i) = dsda(psi_grid(i))*0.5*inv_ds*(phi00(i+1)-phi00(i-1))
         enddo
         ! 3-point forward (or backward) difference at the boundaries
         phip00(0) = dsda(psi_grid(0))*0.5*inv_ds*(-3.0*phi00(0)+4.0*phi00(1)-phi00(2))
         phip00(mpsi) = dsda(psi_grid(mpsi))*0.5*(-inv_ds)*(-3.0*phi00(mpsi)+4.0*phi00(mpsi-1)-phi00(mpsi-2))

         if(ZF_SEPARATION==1) then
            ! Subtract zonal_phi so get nonzonal phi
            do i=0,mpsi
               do j=0,mtheta(i)
                  ij=igrid(i)+j
                  phi(:,ij) = phi(:,ij) - phi00(i)
               enddo
            enddo
         endif

         return
      endif

      if(ZF_SEPARATION==1) then

         if( solver%order <= 0 ) then    ! old solver
            call solve( solver, masssolver, mgrid, dentmp, xx, ierr )

            phi(1,:) = xx

            ! radial boundary  Note wx: setting phi=0 at both coundaries could enforce a constraint on
            !                           long wave length modes ---- a constraint on large scale Er (???)
            !    3 cases:
            ! 1. phi=0 at both inner and outer boundaries
            ! 2. phi=0 at inner boundary and float phi at outer boundary
            ! 3. phi is float at both inner and outer boundaries
            !    for above three cases, electric field E=grad\phi=0 is enforced at both boundaries

            !        phi(1,igrid(0):igrid(0)+mtheta(0))=0.0
            !          potential at one of boundaries should be float
            !        phi(1,igrid(mpsi):igrid(mpsi)+mtheta(mpsi))=0.0

            do i=0, mpsi
               phi(1,igrid(i))=phi(1,igrid(i)+mtheta(i))  !validate phi(1,igrid(i)) for safe
            enddo

         else             ! FEM solver solves delta-phi and d<phi>/da interatively

            !   adum_11=0.0_wp
            ! do it=1,its

            if(mype==0 .and. mod(istep,msnap)==0 ) open(unit=191209,file='phip00_its.dat',position='append',action='write')
            phip00_convergency=1.0 ! any initial value to go inside do while
            it=0

            do while(phip00_convergency>its_phip00_tolerance .and. it<its)
               it=it+1
               !  phip00=0.0_wp    ! for debug

               if(ZF_SEPARATION==-100) then   !  direct numerical derivative to calculate -beta*div[sum(gs)*d<phi>/da*grad(a)]

                  ! direct numerical derivative to get d/da(d<phi>/da)
                  !$omp parallel do private(i,dpsi1,dpsi2,dpsi1_inv,dpsi2_inv,dpsi_inv,e_psi2,e_psi1)
                  do i=1,mpsi-1
                     dpsi1=psi_grid(i)-psi_grid(i-1)
                     dpsi2=psi_grid(i+1)-psi_grid(i)
                     dpsi1_inv=1.0_wp/dpsi1
                     dpsi2_inv=1.0_wp/dpsi2
                     dpsi_inv=1.0_wp/(dpsi1+dpsi2)
                     e_psi2=(phip00(i+1)-phip00(i))*dpsi2_inv
                     e_psi1=(phip00(i)-phip00(i-1))*dpsi1_inv
                     adum0(i)=(e_psi1*dpsi2+e_psi2*dpsi1)*dpsi_inv    !  =d/da(d<phi>/da)
                  enddo
                  !   boundary conditions d<phi>/da=0 (not phi=0 at both boundaries)
                  adum0(0)=0.0_wp
                  adum0(mpsi)=0.0_wp

                  !omp parallel do private(i,j,ij)
                  do i=0, mpsi
                     do j=0,mtheta(i)
                        ij=igrid(i)+j
                        yy(ij)=div1(ij)*phip00(i)+div2(ij)*adum0(i)   ! = -beta*div[sum(gs)*d<phi>/da*grad(a)]
                     end do
                     yy(igrid(i))=yy(igrid(i)+mtheta(i))
                  enddo

                  bb2(:)=dentmp(:) + yy(:)    ! = b - beta*div[sum(gs)*d<phi>/da*grad(a)]

                  call solve( solver, masssolver, mgrid, bb2, xx, ierr )  !  --> delta-phi for given phip00=d<phi>/da

               else      ! use FEM inside the solver for adding -beta*div[sum(gs)*d<phi>/da*grad(a)]

                  call solve( solver, masssolver, mgrid, dentmp, xx, ierr )  !  --> delta-phi for given phip00=d<phi>/da


                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !                  !  replaced by the following to ensure <delta-phi>=0
                  !
                  !    !omp parallel do private(i,j,ij)
                  !                  do i=0, mpsi
                  !                     do j=0,mtheta(i)
                  !                        ij=igrid(i)+j
                  !                        bb2(ij)=dentmp(ij)+adum_11(i)   ! = add alpha*<delta-phi> to RHS of delta-phi equation solved by FEM
                  !                     end do
                  !                     bb2(igrid(i))=bb2(igrid(i)+mtheta(i))
                  !                  enddo
                  !
                  !                  call solve( solver, masssolver, mgrid, bb2, xx, ierr )  !  --> delta-phi for given phip00=d<phi>/da
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


               end if

               !  if(mype==0)   print*, 'it=', it, '   istep=', istep, '   irk=', irk, '   d<phi>/da=', phip00(mpsi/2),  '   delta-phi=', xx(mgrid/2)

               do i=0, mpsi
                  xx(igrid(i))=xx(igrid(i)+mtheta(i))  !validate xx(1,igrid(i)) for safe
               enddo


               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !                  adum_11=0.0_wp
               !                  !$omp parallel private(i,j,ij) reduction(+:adum_11)
               !                  !$omp do
               !                  do i=0,mpsi
               !                     do j=1,mtheta(i)
               !                        ij=igrid(i)+j
               !                        adum_11(i)=adum_11(i)+dtemper(ij)*xx(ij)  ! surface averaged <delta-phi>
               !                     enddo
               !                  enddo
               !                  !$omp end parallel
               !
               !                  phi_flux_ave=0.0_wp
               !                  call MPI_ALLREDUCE(adum_11,phi_flux_ave,mpsi+1,mpi_Rsize, MPI_SUM,toroidal_comm,ierr)
               !                  phi_flux_ave=phi_flux_ave/jaco_sum       ! =<delta-phi>
               !                  adum_11=phi_flux_ave*abs(csc0)/Te_grid   ! =(e/e0)*(m0/ms)*<delta-phi>/T_e with <phi>, T_e normalized, which should be zero by defination
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               !  call find_dxda( masssolver, mgrid, xx, add_phi, ierr )
               !  call find_dxdt( masssolver, mgrid, add_phi, yy, ierr )

               call find_dxdadxdt( masssolver, mgrid, xx, add_phi, yy, ierr )

               add_phi(:)=(add_phi(:)*gij(1,:)+yy(:)*gij(4,:))*gsij(:)   ! --> sum(gs)*grad_perp(delta-phi).grad(a)
               adum=0.0_wp
               !$omp parallel private(i,j,ij) reduction(+:adum)
               !$omp do
               do i=0,mpsi
                  do j=1,mtheta(i)
                     ij=igrid(i)+j
                     adum(i)=adum(i)+dtemper(ij)*add_phi(ij)   ! surface averaging
                  enddo
               enddo
               !$omp end parallel

               phi_flux_ave=0.0_wp
               call MPI_ALLREDUCE(adum,phi_flux_ave,mpsi+1,mpi_Rsize, MPI_SUM,toroidal_comm,ierr)
               phi_flux_ave=phi_flux_ave/jaco_sum    !  =<sum(gs)*grad_perp(delta-phi).grad(a)>

               old_phip00 = phip00

               phip00=0.0_wp
               do i=1,mpsi
                  phip00(i)=phip00(i-1)+0.5_wp*(psi_grid(i)-psi_grid(i-1))*  &
                     (-den_grid(i-1)*dVda(i-1)*(zonali(i-1)-zonale(i-1))  &
                     -den_grid(i)*dVda(i)*(zonali(i)-zonale(i)) )
               enddo

               phip00=phip00/dVda/gsgaa
               phip00=phip00 -(phi_flux_ave*dVda-phi_flux_ave(0)*dVda(0))/dVda/gsgaa   ! --> phip00=dphi_00/da with phip00(0)=0

               phip00=phip00+Er_a0        ! --> phip00=dphi_00/da with phip00(0)=Er_a0 which is Er_NC when ifneoclassical=1

               phip00_convergency = sum(sqrt((phip00 - old_phip00)**2.0))/sum(sqrt(old_phip00**2.0))
               if(mype==0 .and. mod(istep,msnap)==0 ) write(191209,"(2(I8,1x),200(ES16.8,1x))") istep, it, phip00_convergency, phip00(0:mpsi)
            end do     !  end iteration for FEM
            if(mype==0 .and. mod(istep,msnap)==0 ) then
               call flush(191209)
               close(191209)
            endif

            phi(1,:) = xx

            ! (0,0) mode potential store in phi00
            phi00=0.0_wp   ! phi00(0)=0.0 is chosen
            do i=1,mpsi
               phi00(i)=phi00(i-1)+0.5_wp*(psi_grid(i)-psi_grid(i-1))*(phip00(i-1)+phip00(i))
            enddo

         end if

      else      ! ZF_SEPARATION=0 --> solve total potential

         adum0=phi00
         do it=1,its
            !   calculate <phi> using phi of last time step for the 1st iteration

            ! if(mype==-100) then      !  keep the following coding in the record
            adum=0.0_wp
            rms=0.0_wp
            phi_flux_ave=0.0_wp
            !$omp parallel private(i,j,ij) reduction(+:adum,rms)
            !$omp do
            do i=0,mpsi
               do j=1,mtheta(i)
                  ij=igrid(i)+j
                  adum(i)=adum(i)+dtemper(ij)*phi(1,ij)  ! surface averaging
                  rms=rms+phi(1,ij)**2
               enddo
            enddo
            !$omp end parallel
            call MPI_ALLREDUCE(adum,phi_flux_ave,mpsi+1,mpi_Rsize, MPI_SUM,toroidal_comm,ierr)
            phi_flux_ave=phi_flux_ave/jaco_sum    !  =<phi>

            rtc=phi_flux_ave(0)

            ! re-level phi and <phi> by subtracting a constant for original integral solver
            if(solver%order <= 0 )  phi_flux_ave=phi_flux_ave - phi_flux_ave(0)

            t1=sum((phi_flux_ave-adum0)**2)/sum(adum0**2)  ! calculated for checking convergence
            adum0=phi_flux_ave

            phi(1,:)=phi(1,:)-phi_flux_ave(0)
            if(solver%order <= 0 ) then   ! use original integral solver
               phi_flux_ave=phi_flux_ave*tite        !  =<phi>*(Ti/Te)
            else                          ! use FEM solver
               phi_flux_ave=phi_flux_ave*abs(csc0)/Te_grid   ! =(e/e0)*(m0/ms)*<phi>/T_e with <phi>, T_e normalized
            endif
            ! endif


            if(ZF_SEPARATION==-100) then

               ! Re-level phi first (there is an arbitrary constant in the potential)

               !   calculate <phi> at inner boundary a=a0
               rms=0.0_wp
               !$omp parallel private(j,ij) reduction(+:rms)
               !$omp do
               do j=1,mtheta(0)
                  ij=igrid(0)+j
                  rms=rms+dtemper(ij)*phi(1,ij)  ! surface averaging
               enddo
               !$omp end parallel
               call MPI_ALLREDUCE(rms,tmp,1,mpi_Rsize,MPI_SUM,toroidal_comm,ierr)
               tmp=tmp/jaco_sum(0)    !  =<phi> at inner boundary a=a0

               phi=phi-tmp  ! Re-level potential phi by a constant to make <phi>(a0)=0 at inner boundary

               adum=0.0_wp
               rms=0.0_wp
               phi_flux_ave=0.0_wp
               !$omp parallel private(i,j,ij) reduction(+:adum,rms)
               !$omp do
               do i=0,mpsi
                  do j=1,mtheta(i)
                     ij=igrid(i)+j
                     adum(i)=adum(i)+dtemper(ij)*phi(1,ij)  ! surface averaging
                     rms=rms+phi(1,ij)**2
                  enddo
               enddo
               !$omp end parallel
               call MPI_ALLREDUCE(adum,phi_flux_ave,mpsi+1,mpi_Rsize, MPI_SUM,toroidal_comm,ierr)
               phi_flux_ave=phi_flux_ave/jaco_sum    !  =<phi>
               adum=adum/jaco_sum

               rtc=phi_flux_ave(0)

               !    move part of flux surface average <phi> to left
               phi_flux_ave=phi_flux_ave - adum
               phi_flux_ave=phi_flux_ave*tite        !  =<phi>*(Ti/Te)
            endif


            if(mype==0)   print*, 'mype=', mype, 'its=', it, '  <phi>(',ipsi,')=',  phi_flux_ave(ipsi), 'rms=', rms,  '  <phi>(0)=', rtc, &
               'relative error in <phi>^2=', t1
            !   if(mype==180) print*, 'mype=', mype, 'its=', it, '  <phi>(',ipsi,')=',  phi_flux_ave(ipsi), 'rms=', rms,  '  <phi>(0)=', rtc

            add_phi=0.0_wp
            xx=0.0_wp
            do i=0,mpsi
               do j=1,mtheta(i)
                  ij=igrid(i)+j
                  add_phi(ij)=dentmp(ij) + phi_flux_ave(i)
               end do
               add_phi(igrid(i))=add_phi(igrid(i)+mtheta(i))
            end do

            call solve( solver, masssolver, mgrid, add_phi, xx, ierr )
            phi(1,:) = xx

            do i=0, mpsi
               phi(1,igrid(i))=phi(1,igrid(i)+mtheta(i))  !validate phi(1,igrid(i)) for safe
            enddo

         end do  ! end of iteration do it=1,its

      end if   ! end for ZF_SEPARATION

   end if   ! end for ifelectron

   ! debug
   do i=0,mpsi
      if(phi(1,igrid(i)).ne.phi(1,igrid(i)+mtheta(i)))then
         !validate phi(1,igrid(i)) for safe
         stop 'x != x'
      endif
   enddo

   if(mype==-10)then
      write(300+numberpe,*)phi(1,:)
      !stop 'DEBUG'
   endif

   !
   !   if(mype==0.and.istep==10) then
   !      do i=1,mgrid
   !         write(18, 100) i, phi(1,i), add_phi(i), dentmp(i)
   !      end do
   !   end if
   ! 100 format(1x,i6, 3(2x,d12.6))
   !     call flush(18)
   !

end subroutine poisson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine poisson_initial( mat )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine poisson_initial( mat )
   use precision
   use mpi_defs
   use global_parameters
   use field_array
   use diagnosis_array, only: jaco_sum
   use esi_interface
   use phit_matrix_class
   implicit none
   !
   type(phit_matrix_data) :: mat
   ! locals
   integer mring,mindex,ierr

   integer i,ii,ij,ij0,ipjt,j,jj,jt,j1,j0,jm,k,kr,kp,nt,np,minn,maxn
   real(wp) vring(3),fring(3),rgrid,tgrid,ddelr,ddelt,wght,r,rs,t,rdum,wr,&
      wt1,wt0,tdum,zdum,pi2_inv,delr,delt(0:mpsi),psimin,psimax, vring_r(3)

   real(wp), dimension(10) :: aaa,the,rj,zj,Bj

   real(wp) tmp,rms

   ! number of gyro-ring
   mring=2

   ! number of summation: maximum is 32*mring+1
   !  mindex=32*mring+1
   !  mindex=45
   mindex=58
   !  mindex=37 !for a=700,1000

   !allocate(indexp(mindex,mgrid,mzeta),ring(mindex,mgrid,mzeta),nindex(mgrid,mzeta))
   !  allocate(save_flux_ave(0:mpsi))

   if(mring==1)then
      ! one ring, velocity in unit of v_th: vring=v_perp/v_th
      vring(1)=sqrt(2.0_wp)
      fring(1)=1.0_wp
   elseif(mring==2)then
      ! two rings good for up to k_perp=1.5
      vring(1)=0.9129713024553_wp
      vring(2)=2.233935334042_wp
      fring(1)=0.7193896325719_wp
      fring(2)=0.2806103674281_wp
   else
      ! three rings: exact(<0.8%) for up to k_perp=1.5
      vring(1)=0.388479356715_wp
      vring(2)=1.414213562373_wp
      vring(3)=2.647840808818_wp
      fring(1)=0.3043424333839_wp
      fring(2)=0.5833550690524_wp
      fring(3)=0.1123024975637_wp
   endif

   !  vring=vring*vthi    ! for uniform Ti

   psimin=s_grid(0)     ! working with radial coordinator s
   psimax=s_grid(mpsi)
   pi2_inv=0.5_wp/pi
   delr=1.0_wp/deltar
   delt=2.0_wp*pi/deltat
   !nindex=0
   !ring=0.0_wp
   !indexp=1
   do k=1,1 !mzeta
      zdum=zetamin+deltaz*real(k)
      !  !$omp parallel do private(i,j,kr,kp,np,nt,ij0,rgrid,tgrid,jt,aaa,the, &
      !  !$omp&                    rj,zj,Bj,ipjt,ddelr,ddelt,wght,&
      !  !$omp&                    r,t,rdum,ii,wr,tdum,j1,wt1,j0,wt0,ij,rs)
      do i=0,mpsi
         do j=1,mring
            vring_r(j)=vring(j)*vth_grid(i)
         enddo
         do j=1,mtheta(i)
            !do j=0,mtheta(i)
            ij0=igrid(i)+j
#ifndef NO_PETSC
            rs = 0.25_wp
            call add_gk_value( mat, ij0, ij0, rs )
#else
            ! 1st point is the original grid point
            nindex(ij0,k)=1
            indexp(1,ij0,k)=ij0
            ring(1,ij0,k)=0.25_wp
#endif
            ! position of grid points
            rgrid=s_grid(i)
            tgrid=deltat(i)*real(j)+zdum*qtinv(i)
            tgrid=tgrid*pi2_inv
            tgrid=2.0_wp*pi*(tgrid-aint(tgrid))
            !  jt=max(0,min(mtheta(i),int(pi2_inv*delt(i)*tgrid+0.5_wp)))

            ! B-field
            aaa(1)=psi_grid(i)  ! equilibrium interface using a as radial coordinate
            the(1)=tgrid
            call esigetrzb(rj,zj,Bj,aaa,the,1)
            !  ipjt=igrid(i)+jt
            ipjt=ij0            ! 2013

            do kr=1,mring

               ! FLR from grid point and weight of larmor_num-points for each ring
               do kp=1,larmor_num
                  ddelr=pgyro(kp,ipjt)
                  ddelt=tgyro(kp,ipjt)
                  wght=(inv_larmor_num)*fring(kr)

                  ! position for each point with rho_i=2.0*vring
                  r=max(psimin,min(psimax,rgrid+ddelr*2.0_wp*vring_r(kr)*sqrt(0.5_wp/Bj(1)))) ! r is s
                  t=tgrid+ddelt*2.0_wp*vring_r(kr)*sqrt(0.5_wp/Bj(1))

                  ! linear interpolation
                  rdum=delr*max(0.0_wp,min(s1-s0,r-s0))
                  ii=max(0,min(mpsi-1,int(rdum)))
                  wr=rdum-real(ii)
                  if(wr>0.95_wp)wr=1.0_wp
                  if(wr<0.05_wp)wr=0.0_wp

                  ! outer flux surface
                  tdum=t-zdum*qtinv(ii+1)
                  tdum=tdum*pi2_inv+10.0_wp
                  tdum=delt(ii+1)*(tdum-aint(tdum))
                  j1=max(0,min(mtheta(ii+1)-1,int(tdum)))
                  wt1=tdum-real(j1)
                  if(wt1>0.95_wp)wt1=1.0_wp
                  if(wt1<0.05_wp)wt1=0.0_wp

                  ! inner flux surface
                  tdum=t-zdum*qtinv(ii)
                  tdum=tdum*pi2_inv+10.0_wp
                  tdum=delt(ii)*(tdum-aint(tdum))
                  j0=max(0,min(mtheta(ii)-1,int(tdum)))
                  wt0=tdum-real(j0)
                  if(wt0>0.95_wp)wt0=1.0_wp
                  if(wt0<0.05_wp)wt0=0.0_wp

                  ! index and weight of each point
                  do np=1,4
                     if(np==1)then
                        ij=igrid(ii+1)+j1+1
                        rs=wght*wr*wt1
                     elseif(np==2)then
                        if(j1==0)j1=mtheta(ii+1)
                        ij=igrid(ii+1)+j1
                        rs=wght*wr*(1.0-wt1)
                     elseif(np==3)then
                        ij=igrid(ii)+j0+1
                        rs=wght*(1.0_wp-wr)*wt0
                     else
                        if(j0==0)j0=mtheta(ii)
                        ij=igrid(ii)+j0
                        rs=wght*(1.0_wp-wr)*(1.0_wp-wt0)
                     endif
#ifndef NO_PETSC
                     call add_gk_value( mat, ij0, ij, rs )
#else
                     ! insignificant point replaced by the original grid point
                     if(rs<0.001_wp)then
                        ring(1,ij0,k)=ring(1,ij0,k)+rs
                        goto 100
                     endif
                     do nt=1,nindex(ij0,k)
                        ! redundant point
                        if(ij==indexp(nt,ij0,k))then
                           ring(nt,ij0,k)=ring(nt,ij0,k)+rs
                           goto 100
                        endif
                     enddo
                     ! new point
                     nindex(ij0,k)=nindex(ij0,k)+1
                     nt=nindex(ij0,k)
                     indexp(nt,ij0,k)=ij
                     ring(  nt,ij0,k)=rs
                     100                 continue
#endif
                  enddo  !end of 4-point interpolation loop
               enddo     !end of 8-point-per-ring loop
            enddo        !end of ring loop
         enddo           !end of poloidal loop
      enddo              !end of radial loop
   enddo                 !end of toroidal loop


   if(ZF_SEPARATION==-100) then
      if(ifelectron/=1.and.ZF_SEPARATION==0) then
         !    move part of flux surface average <phi> to left
         do i=0,mpsi
            do j=1,mtheta(i)
               ij0=igrid(i)+j
               do k=1,mtheta(i)
                  ij=igrid(i)+k
                  rs=dtemper(ij)/jaco_sum(i)*tite(i)
                  call add_gk_value( mat, ij0, ij, rs )
               enddo
            enddo
         enddo
      endif
   endif


   ! Note: no valid data for nindex, ring and indexp at igrid(i)

   ! check array size
   !  if(maxval(nindex)>mindex)then
   !     write(*,*)'Poisson error',mype,maxval(nindex),mindex
   !     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
   !  endif

   !  rdum=0.0_wp
   !  tdum=0.0_wp
   !  zdum=1.0_wp
   !  do k=1,mzeta
   !     do i=1,mgrid
   !        rdum=sum(ring(1:nindex(i,k),i,k))
   !        tdum=max(tdum,rdum)
   !        zdum=min(zdum,rdum)
   !     enddo
   !  enddo
   !  if(mype==0)write(stdout,*)'poisson solver=',maxval(nindex),minval(nindex),&
   !       tdum,zdum,mgrid,sum(nindex)

end subroutine poisson_initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
