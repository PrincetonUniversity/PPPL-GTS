subroutine snapshot
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use field_array
   use diagnosis_array
   use esi_interface
   !ZHIXIN control output of PHI
   use  control_parameters

   implicit none

   integer,parameter :: mbin_psi=5, mbin_u=32
   integer i,ibin,ip,iu,j,jj,jm,jt,k,kz,m,nsnap,ierror,icount,ij
   real(wp) Emax,eright,uright,uleft,pdum,tdum,b,upara,energy,pitch,q,dela,&
            delu,g_inv,qdum,r,delr,delz,phisn(mtheta(ipsi)/2),&
            eachpsi(mtheta(ipsi)/2*mzetamax),wt,delt(0:mpsi),&
            weight,adum(0:mpsi),bdum(mbin_u,mbin_psi)
   real(wp),dimension(0:mpsi) :: marker,fflows,dflows,ftem,dtem
   real(wp),dimension(mbin_u,mbin_psi) :: ubin,dubin,pbin,dpbin,ebin,debin
   real(wp) ubin_max,pbin_max,ebin_max
   !  real(wp),dimension(mzeta,0:mpsi,0:mthetamax) :: ftotal,tpara,tperp

   character(len=13) cdum,date,time

   real(wp),dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj


   ! write particle data for re-run
   !  call restart_write
   !ZHIXIN use flag irewrite to output RESTART
   ! if(mod(istep,40) .eq. 0) call restart_write
   if(mod(istep,irewrite) .eq. 0) call restart_write


   ! number of poloidal grid
   jm=mtheta(ipsi)/2

   ! zion space grids: parallel velocity, kinetic energy, pitch angle
   if(ifelectron==0) then
      uright=umax/sqrt(aion)
   else
      uright=umax/sqrt(aelectron)
   end if
   uleft=-uright
   dela=real(mbin_psi)/(a1-a0)
   delu=1.0/(uright-uleft)
   Emax=umax*umax/4.0
   eright=1.0/Emax !checked,same for e and i with T_e ~ T_i ?why Emax=umax^2/4.0
   g_inv=1.0/vthc
   delr=1.0/deltar
   delt=1.0/deltat
   delz=1.0/deltaz

   ! radial and 3-d volume data
   ubin=0.0
   pbin=0.0
   ebin=0.0
   dubin=0.0
   dpbin=0.0
   debin=0.0
   marker=0.0
   fflows=0.0
   dflows=0.0
   ftem=0.0
   dtem=0.0

   if(ifelectron==0) then

      !     call esigetb(zion(1,:),zion(2,:),mi)   ! get B
      call esigetb(zion,mi,7)   ! get B-->BB

      do m=1,mi
         weight=zion0(6,m)
         r=zion(1,m)
         ibin=max(1,min(mbin_psi,1+int((r-a0)*dela)))
         ip=max(0,min(mpsi,int((zion(7,m)-s0)*delr+0.5)))
                            ! zion(2,*) and zion(3,*) are within [0,2*pi)
         jt=max(0,min(mtheta(ip),int(zion(2,m)*delt(ip)+0.5)))

         upara=zion(4,m)*BB(m)*qion/aion*vth_inv(ip)
         energy=max(1.0e-20,0.5*aion*upara*upara+&
              zion(6,m)*zion(6,m)*BB(m)*vth_inv(ip)*vth_inv(ip))
         pitch=upara/sqrt(2.0*energy/aion)

         !       iu=max(1,min(mbin_u,1+int(real(mbin_u)*(upara-uleft)*delu)))
         iu=1+int(real(mbin_u)*(upara-uleft)*delu)
         if(iu>=1.and.iu<=mbin_u) then
            ubin(iu,ibin)=ubin(iu,ibin)+weight
            dubin(iu,ibin)=dubin(iu,ibin)+zion(5,m)
         end if

         !       iu=max(1,min(mbin_u,1+int(real(mbin_u)*energy*eright)))
         iu=1+int(real(mbin_u)*energy*eright)
         if(iu>=1.and.iu<=mbin_u) then
            ebin(iu,ibin)=ebin(iu,ibin)+weight
            debin(iu,ibin)=debin(iu,ibin)+zion(5,m)
         end if

         iu=max(1,min(mbin_u,1+int(real(mbin_u)*0.5*(pitch+1.0))))
         pbin(iu,ibin)=pbin(iu,ibin)+weight
         dpbin(iu,ibin)=dpbin(iu,ibin)+zion(5,m)

         ! radial profile of marker, v_para, and temperature
         marker(ip)=marker(ip)+weight
         fflows(ip)=fflows(ip)+upara*weight
         dflows(ip)=dflows(ip)+zion(5,m)*upara
         ftem(ip)=ftem(ip)+energy*weight
         dtem(ip)=dtem(ip)+zion(5,m)*energy

      ! parallel and perpendicular pressure perturbation
      !     ftotal(ip,jt)=ftotal(ip,jt)+weight
      !     tpara(ip,jt)=tpara(ip,jt)+zion(5,m)*0.5*upara*upara
      !     tperp(ip,jt)=tperp(ip,jt)+zion(5,m)*energy
      enddo

   else

      !     call esigetb(zeon(1,:),zeon(2,:),me)   ! get B
      call esigetb(zeon,me,zeondim)   ! get B

      do m=1,me
         weight=zeon0(6,m)

         r=zeon(1,m)       ! r=sqrt(2.0*zeon(1,m))
         ibin=max(1,min(mbin_psi,1+int((r-a0)*dela)))
         ip=max(0,min(mpsi,int((r-a0)*delr+0.5)))
         jt=max(0,min(mtheta(ip),int(zeon(2,m)*delt(ip)+0.5)))

         upara=zeon(4,m)*BB(m)*qelectron/aelectron*g_inv
         energy=max(1.0e-20,0.5*aelectron*upara*upara+  &
              zeon(6,m)*zeon(6,m)*BB(m)*g_inv*g_inv)
         pitch=upara/sqrt(2.0*energy/aelectron)

         iu=max(1,min(mbin_u,1+int(real(mbin_u)*(upara-uleft)*delu)))
         ubin(iu,ibin)=ubin(iu,ibin)+weight
         dubin(iu,ibin)=dubin(iu,ibin)+zeon(5,m)

         iu=max(1,min(mbin_u,1+int(real(mbin_u)*energy*eright)))
         ebin(iu,ibin)=ebin(iu,ibin)+weight
         debin(iu,ibin)=debin(iu,ibin)+zeon(5,m)

         iu=max(1,min(mbin_u,1+int(real(mbin_u)*0.5*(pitch+1.0))))
         pbin(iu,ibin)=pbin(iu,ibin)+weight
         dpbin(iu,ibin)=dpbin(iu,ibin)+zeon(5,m)

         ! radial profile of marker, v_para, and temperature
         marker(ip)=marker(ip)+weight
         fflows(ip)=fflows(ip)+upara*weight
         dflows(ip)=dflows(ip)+zeon(5,m)*upara
         ftem(ip)=ftem(ip)+energy*weight
         dtem(ip)=dtem(ip)+zeon(5,m)*energy

      ! parallel and perpendicular pressure perturbation
      !     ftotal(ip,jt)=ftotal(ip,jt)+weight
      !     tpara(ip,jt)=tpara(ip,jt)+zeon(5,m)*0.5*upara*upara
      !     tperp(ip,jt)=tperp(ip,jt)+zeon(5,m)*energy
      enddo

   endif

   !  do i=0,mpsi
   !     tpara(:,i,mtheta(i))=tpara(:,i,mtheta(i))+tpara(:,i,0)
   !     tperp(:,i,mtheta(i))=tperp(:,i,mtheta(i))+tperp(:,i,0)
   !     ftotal(:,i,mtheta(i))=ftotal(:,i,mtheta(i))+ftotal(:,i,0)
   !  enddo

   ! do i=0,mpsi
   !    ftotal(:,i,1:mtheta(i))=max(1.0e-6,ftotal(:,i,1:mtheta(i)))
   !     tpara(:,i,1:mtheta(i))=tpara(:,i,1:mtheta(i))/ftotal(:,i,1:mtheta(i))
   !     tperp(:,i,1:mtheta(i))=tperp(:,i,1:mtheta(i))/ftotal(:,i,1:mtheta(i))
   !     ftotal(:,i,1:mtheta(i))=ftotal(:,i,1:mtheta(i))*markeri(:,i,1:mtheta(i))
   !     ftotal(:,i,0)=ftotal(:,i,mtheta(i))
   !  enddo

   icount=mbin_psi*mbin_u
   call MPI_REDUCE(ubin,bdum,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
   if(mype==0)ubin=bdum
   call MPI_REDUCE(dubin,bdum,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
   if(mype==0)dubin=bdum
   call MPI_REDUCE(pbin,bdum,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
   if(mype==0)pbin=bdum
   call MPI_REDUCE(dpbin,bdum,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
   if(mype==0)dpbin=bdum
   call MPI_REDUCE(ebin,bdum,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
   if(mype==0)ebin=bdum
   call MPI_REDUCE(debin,bdum,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
   if(mype==0)debin=bdum
   call MPI_REDUCE(marker,adum,mpsi+1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
   if(mype==0)marker=adum
   call MPI_REDUCE(fflows,adum,mpsi+1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
   if(mype==0)fflows=adum
   call MPI_REDUCE(dflows,adum,mpsi+1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
   if(mype==0)dflows=adum
   call MPI_REDUCE(ftem,adum,mpsi+1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
   if(mype==0)ftem=adum
   call MPI_REDUCE(dtem,adum,mpsi+1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
   if(mype==0)dtem=adum

   where(marker .lt. 1.0)marker=1.0
   fflows=fflows/marker
   dflows=dflows/marker
   ftem=ftem/marker
   dtem=dtem/marker
   marker(1:mpsi)=marker(1:mpsi)*pmarki(1:mpsi)

   ! gather flux surface profile
   eachpsi=0.0
   do j=1,jm
      phisn(j)=phi(1,igrid(ipsi)+2*j)-real(1-ZF_SEPARATION,wp)*phi00(ipsi)  ! --> delta-phi
   enddo
   call MPI_GATHER(phisn,jm,mpi_Rsize,eachpsi,&
        jm,mpi_Rsize,0,toroidal_comm,ierror)

   if(mype == 0)then !mype=0 must also have myrank_toroidal=0
      ! normalization
      ubin_max=maxval(ubin)
      pbin_max=maxval(pbin)
      ebin_max=maxval(ebin)
      dubin=dubin/ubin_max
      dpbin=dpbin/pbin_max
      debin=debin/ebin_max

      ! open snapshot output file
      nsnap=istep/msnap
      write(cdum,'("snap",i5.5,".out")')nsnap

      open(snapout,file=cdum,status='replace')

      ! time step and grid dimension,write every other point in poloidal direction
      write(snapout,102) real(istep)*tstep, q_grid(ipsi)
      write(snapout,101) mbin_u,mbin_psi,mpsi/2,jm,mzetamax,mpsi/4   ! mpsi/2 specifies
                     ! how many number of data are written out; don't confuse with ipsi!

      ! write velocity space grids: V_para, pitch, energy to data file
      ! number of quantities
      i=9
      write(snapout,101)i
      do i=1,mbin_u
         write(snapout,102)uleft+(uright-uleft)*real(i-1)/real(mbin_u-1)
      enddo
      do i=1,mbin_u
         write(snapout,102)-1.0+2.0*real(i-1)/real(mbin_u-1)
      enddo
      do i=1,mbin_u
         write(snapout,102)1.0/eright*real(i-1)/real(mbin_u-1)
      enddo
      ! velocity space distribution function
      do j=1,mbin_psi
         write(snapout,102)(ubin(i,j),i=1,mbin_u)
         write(snapout,102)(dubin(i,j),i=1,mbin_u)
         write(snapout,102)(pbin(i,j),i=1,mbin_u)
         write(snapout,102)(dpbin(i,j),i=1,mbin_u)
         write(snapout,102)(ebin(i,j),i=1,mbin_u)
         write(snapout,102)(debin(i,j),i=1,mbin_u)
      enddo

      ! radial grid
      i=8
      write(snapout,101)i
      do i=1,mpsi,2
         write(snapout,102)  psi_grid(i)  !  (a0+deltar*real(i))/gyroradius  ????
      enddo
      write(snapout,102)(marker(i),i=1,mpsi,2)
      !     write(snapout,102)(phi00(i)/vthc**2,i=1,mpsi)
      do i=1,mpsi,2
         r=psi_grid(i)       !   r=(a0+deltar*real(i))/a  ????
         !        r=r_grid(i)        !   use this minor radiu
         q=q_grid(i)
         write(snapout,102)q/r
      enddo
      write(snapout,102)(phip00(i)/gyroradius,i=1,mpsi,2)   ! ???
      write(snapout,102)(fflows(i),i=1,mpsi,2)
      write(snapout,102)(dflows(i),i=1,mpsi,2)
      write(snapout,102)(ftem(i),i=1,mpsi,2)
      write(snapout,102)(dtem(i),i=1,mpsi,2)

      ! poloidal cross section: marker,density,potential,parallel flow and temperature
      i=3
      write(snapout,101)i
      do j=0,jm
         jj=j
         if(j>jm/2)jj=j-jm
         do i=1,mpsi,2
            pdum=psi_grid(i)  !wx
            tdum=pi+2.0*pi*real(jj)/real(jm)
            aaa(1)=pdum
            the(1)=tdum
            call esigetrzb(rj,zj,Bj,aaa,the,1)
            write(snapout,102) rj(1)
         enddo
      enddo
      do j=0,jm
         jj=j
         if(j>jm/2)jj=j-jm
         do i=1,mpsi,2
            pdum=psi_grid(i)      !wx
            tdum=pi+2.0*pi*real(jj)/real(jm)
            aaa(1)=pdum
            the(1)=tdum
            call esigetrzb(rj,zj,Bj,aaa,the,1)
            write(snapout,102) zj(1)
         enddo
      enddo

      do j=0,jm
         jj=j
         if(j>jm/2)jj=j-jm
         do i=1,mpsi,2
            tdum=pi+2.0*pi*real(jj)/real(jm)
            jt=max(1,min(mtheta(i),1+int(tdum/deltat(i))))
            wt=tdum/deltat(i)-real(jt-1)
            write(snapout,102)(wt*(phi(0,igrid(i)+jt)-real(1-ZF_SEPARATION)*phi00(i))+  &
                              (1.0-wt)*(phi(0,igrid(i)+jt-1)-real(1-ZF_SEPARATION)*phi00(i)))/  &
                               vthc**2  ! vth_grid(i)**2  ! should replace by local Ti??

            !           write(snapout,102)wt*densityi(0,igrid(i)+jt)+(1.0-wt)*densityi(0,igrid(i)+jt-1)

         enddo
      enddo

      ! flux surface profile for marker, density, potential, parallel flow and temperature
      i=1
      write(snapout,101)i
      do j=1,jm
         do kz=1,ntoroidal
            write(snapout,102) eachpsi((kz-1)*jm+(j-1)+1)/ &
                  vth_grid(ipsi)**2  ! eachpsi is on ipsi with local Ti=vth**2
         enddo
      enddo

      ! S.Ethier 7/5/2007
      ! eigenmode structure
      eigenmode=0._wp
      write(snapout,101)num_mode,m_poloidal,nmode
      do i=1,mpsi,2
         do kz=1,num_mode
            do j=1,m_poloidal
               write(snapout,102)eigenmode(j,kz,i)
            enddo
         enddo
      enddo

      ! close snapshot file
      close(snapout)
   endif

   ! record program end time
   if(istep .eq. mstep)then
      call date_and_time(date,time)
      if(mype == 0)write(stdout,*) 'Program ends at DATE=', date, 'TIME=', time
   endif

   101 format(i6)
   !  102 format(e10.4)


   !  105 format(e10.4,4x,6hflux2x,4x,i4)
   !  106 format(e10.4,4x,6hflux2z,4x,i4)
   !  107 format(e10.4,2x,i4,2x,4hpsi=,e10.4,2x,4hthe=,e10.4,2x,2hR=,e10.4,2x,2hz=,e10.4)

   102 format(d16.10)

   105 format(d16.10,4x,6hflux2x,4x,i4)
   106 format(e16.10,4x,6hflux2z,4x,i4)
   107 format(e16.10,2x,i4,2x,4hpsi=,d16.10,2x,4hthe=,d16.10,2x,2hR=,d16.10,2x,2hz=,d16.10)

   ! REcord Potential. Each myrank_toroidal has its own file (see open_file).
   if(movie_3D_h5==0) then  ! no 3D movie data created
      !ZHIXIN use flag to output PHI
      !     if(istep>0 .and. myrank_partd==0) then
      if(istep>0 .and. myrank_partd==0.and.flag_PHI>0) then
         write(777) istep
         write(777) phi(1,1:mgrid)  ! Note(5/08): phi contains zonal flow when ZF_SEPARATION=0 (solving total potential)
         !        write(777) phi_tmp(1,:)   ! time averaged phi
         call flush(777)  ! otherwise, if job crashes PHI data is not completed, which cause problem for reading PHI
      end if
   end if

   ! Record Density. Each myrank_toroidal has its own file (see open_file).
   !ZHIXIN use flag to output DEN
   !  if(istep>0 .and. myrank_partd==1) then
   if(istep>0 .and. myrank_partd==1.and.flag_DEN>0) then
      write(888) istep
      if(ifelectron==0) then
         write(888) densityi(1,1:mgrid),phi_tmp(1:mgrid)
      else

         !   Total electron density fluctuations = adiabatic + nonadiabatic
         !$omp parallel do private(i,j,ij)
         do i=0,mpsi
            do j=1,mtheta(i)
               ij=igrid(i)+j
               densitye(1,ij)=densitye(1,ij)+phi(1,ij)/Te_grid(i)*abs(csc0)  ! total electron density fluctuation dn_e/n_0
            enddo
            densitye(1,igrid(i))=densitye(1,igrid(i)+mtheta(i))  ! validate densitye(:,igrid(i)
         enddo

         write(888) densityi(1,1:mgrid),phi_tmp(1:mgrid),densitye(1,1:mgrid),heatflux(1:mgrid)   ! delta-ne output added Oct. 2012
                             !  phi_tmp:  Ti fluctuations
                             !  heatflux: Te fluctuations

   !          print*, 'mype= ', mype, 'istep=', istep, phi_tmp(mgrid/2), heatflux(mgrid/2), densitye(1,mgrid/2), densityi(1,mgrid/2)

      end if
      call flush(888)
   end if

   !   if(mype==0) print*, 'Write DENs + PHIs at ', 'istep=', istep,  'irk=', irk

 end subroutine snapshot



! module snapsnap

!   implicit none

!   contains

!     subroutine Interpolate_Data_into_Coarser_Grids(in_Data_Arr1D, gq_jump, a_jump, unit_coeff, out_intpl_Arr1D)

!       use precision
!       use global_parameters
!       use field_array
!       !implicit none

!       integer, intent(in) :: gq_jump, a_jump
!       real(wp), intent(in) :: unit_coeff
!       real(wp), dimension(:), intent(inout) :: in_Data_Arr1D
!       real(wp), dimension(:), intent(inout) :: out_intpl_Arr1D
!       integer i,j,ij,k,jm,jt

!       real(wp) inv_dgq(0:mpsi), gqdum, wt

!       inv_dgq=1.0/deltat

!       jm=mtheta(ipsi)/gq_jump

!       k=1
!       do i=0,mpsi,a_jump
!         do j=0,jm
!           gqdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
!           jt=int(gqdum);  wt=gqdum-jt
!           out_intpl_Arr1D(k) =(1.0-wt)*in_Data_Arr1D(igrid(i)+jt) + wt*in_Data_Arr1D(igrid(i)+jt+1)
!           k=k+1
!         enddo
!       enddo
!       out_intpl_Arr1D = unit_coeff*out_intpl_Arr1D

!     end subroutine Interpolate_Data_into_Coarser_Grids
! end module snapsnap


subroutine Interpolate_Data_into_Coarser_Grids(in_Data_Arr1D, gq_jump, a_jump, out_intpl_Arr1D)

   use precision
   use global_parameters
   use field_array
   implicit none

   integer, intent(in) :: gq_jump, a_jump
   real(wp), dimension(mgrid), intent(inout) :: in_Data_Arr1D
   real(wp), dimension(((mpsi+1)/a_jump)*(mtheta(ipsi)/gq_jump+1)), intent(out) :: out_intpl_Arr1D
   integer i,j,ij,k,jt,jm

   real(wp) inv_dgq(0:mpsi), gqdum, wt

   inv_dgq=1.0/deltat

   jm=mtheta(ipsi)/gq_jump

   k=1
   do i=0,mpsi,a_jump
      do j=0,jm
         gqdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
         jt=int(gqdum);  wt=gqdum-jt
         out_intpl_Arr1D(k) =(1.0-wt)*in_Data_Arr1D(igrid(i)+jt) + wt*in_Data_Arr1D(igrid(i)+jt+1)
         k=k+1
      enddo
   enddo

end subroutine Interpolate_Data_into_Coarser_Grids


subroutine snapshot_hdf5
   use precision
   use global_parameters
   use control_parameters
   use particle_decomp
   use field_array
   use diagnosis_array
   use mpi_defs
   use hdf5
   use h5lt
   use Diagnostic_Grid
   use Structured_Grid
   use MGY_PETSc_Interface_Data,only:CT,Metric_T,Polar_T
   use Ptl_Src_Loss
   implicit none

   integer, parameter :: gq_jump=3, a_jump=1
   integer :: msnap_for_All
   integer :: dim_a, dim_gq
   real(wp) :: gqdum, wt, inv_dgq(0:mpsi)
   ! integer(HSIZE_T), dimension(2) :: h5dims=(/dim_a, dim_gq/))
   integer(HSIZE_T) :: h5dims(5)
   integer(HID_T) :: file_id, gid_PARAMS, gid_FLA, gid_DCOMP, gid_BIN, gid_DG, gid_SG, gid_SG_SM, gid_DG_SM,gid_DG_FLA,gid_DG_FLA_SM, gid_PVS, gid, gid_FP, gid_FP_cum
   integer :: h5err
   logical :: snap_flag=.false.
   ! real(wp) :: intp_val_Arr(((mpsi+1)/a_jump)*(jm+1))  ! 1D representation for 2D data(0:a_jump:mpsi, 0:jm)
   real(wp) :: intp_val, intp_val2
   real(wp), dimension(:), allocatable ::intp_val_Arr, test_dphi_3D  ! 1D representation for 1D data(dims(1)) , 2D data(dims(1), dims(2)) or 3D data(dims(1), dims(2), dims(3))
   integer :: i,j,ij,k,jt,jm,istatus(MPI_STATUS_SIZE),ierror, pbid
   real(wp) full_dphi(0:1,mgrid), intpl_dphi_NZ(0:1), dne_adiabatic(0:1,mgrid)
   real(wp) full_dni(0:1,mgrid), full_dne(0:1,mgrid), FLA_full_dni(0:1,mgrid), FLA_full_dne(0:1,mgrid)
   real(wp) inv_msnap,ndiag_over_msnap, mean_full_dni(0:1,mgrid), mean_full_dne(0:1,mgrid), mean_full_dphi(0:1,mgrid), mean_dne_adiab(0:1)
   real(wp) mean_zonal_dne(0:mpsi), mean_zonal_dni(0:mpsi), mean_zonal_dphi(0:mpsi)
   real(wp) mean_zonal_dTe(0:mpsi), mean_zonal_dTi(0:mpsi), mean_full_dTe(0:1,mgrid), mean_full_dTi(0:1,mgrid)

   real(wp) mean_zonal_dphi_SQ(0:mpsi), mean_zonal_dni_SQ(0:mpsi), mean_zonal_dne_SQ(0:mpsi), mean_zonal_dTi_SQ(0:mpsi), mean_zonal_dTe_SQ(0:mpsi)
   real(wp) mean_full_dphi_SQ(0:1,mgrid), mean_full_dni_SQ(0:1,mgrid), mean_full_dne_SQ(0:1,mgrid), mean_full_dTe_SQ(0:1,mgrid), mean_full_dTi_SQ(0:1,mgrid)

   real(wp) Unit_Conv_Coeff ! Converstion coefficient from CGS to SI
   real(wp) Device_simple_Vol, EE_tot, EE_Z,EE_NZ, EE_FLA, EE_NFLA, recv_0D, tmp_FLA_dphi_nonNaN(0:1,mgrid) ! Electric Energy (EE) with total, Z(Zonal), NZ(NonZonal), FLA and NFLA components of dphi

   real(wp) coarser_dgq, coarser_dpsi1, coarser_dpsi2
   real(wp), dimension(:), allocatable :: intp_a_1D, recvl, recvr
   real(wp), dimension(:,:), allocatable :: intp_mean_full_dphi, right_intp_mean_full_dphi, left_intp_mean_full_dphi
   real(wp), dimension(:,:), allocatable :: intp_FLA_dphi, right_intp_FLA_dphi, left_intp_FLA_dphi
   real(wp), dimension(:,:,:), allocatable :: intp_mean_dE_cov3, right_mean_dE_cov3, left_mean_dE_cov3, intp_FLA_dE_cov3, right_FLA_dE_cov3, left_FLA_dE_cov3
   real(wp) dE_cov3(3,0:1,mgrid)  ! (-evector)=-grad(Phi) , This is instantaneous electric fields which might have many noises.
   real(wp), dimension(:), allocatable :: intp_R_2D, intp_Z_2D, intp_x_2D, intp_y_2D
   real(wp),dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
   character(LEN=128) :: file_path, file_path2
   character(LEN=128) :: cdum

   real(wp) DG_FLA_Data(DG%dim_gq,DG%dim_a,0:1)
   real(wp), dimension(:,:,:),allocatable :: DG_FLA_eScaData, DG_FLA_iScaData
   real(wp) tmp_vec_mgrid_size(mgrid),tmp_Te0_inverse
   real(wp), external::density,temperature_e

   if((myrank_partd==nproc_partd-1) .and. (mod(istep,msnap) .eq. 0)) then
      allocate(DG_FLA_eScaData(DG%dim_gq,DG%dim_a,7:18), DG_FLA_iScaData(DG%dim_gq,DG%dim_a,9:20))
      ! Calculate some DG FLA quantities
      do i=7,18
         call Do_FieldLineAverage_on_DG_by_using_gPPI((DG%eScaData(:,:,i,:)), DG_FLA_Data)
         DG_FLA_eScaData(:,:,i) = DG_FLA_Data(:,:,1)
      enddo
      do i=9,20
         call Do_FieldLineAverage_on_DG_by_using_gPPI((DG%iScaData(:,:,i,:)), DG_FLA_Data)
         DG_FLA_iScaData(:,:,i) = DG_FLA_Data(:,:,1)
      enddo
   endif

   msnap_for_All = msnap*50


   if(istep==0) then
      inv_msnap=1.0
      ndiag_over_msnap=1.0
   else
      inv_msnap=(1.0/msnap)
      ndiag_over_msnap=(real(ndiag)/real(msnap)) ! The quantities cumulated at every ndiag step have (ndiag/msnap) times accumulated data at each msnap step.
   endif


   inv_dgq=1.0/deltat
   jm=mtheta(ipsi)/gq_jump
   dim_a = (mpsi+1)/a_jump
   dim_gq = jm+1

   coarser_dgq = 2.0*pi/(dim_gq-1)
   ! Write Restart files
   !if(mod(istep,irewrite) .eq. 0) call restart_write


   ! Calculate some data before the record
   call Cal_dni_dTi_from_each_ion_weight_component
   call Cal_dne_dTe_from_each_electron_weight_component

   mean_zonal_dphi = inv_msnap*cum_phi00
   mean_zonal_dni  = inv_msnap*den_grid*cum_zonali
   mean_zonal_dne  = inv_msnap*den_grid*cum_zonale
   mean_zonal_dTe  = ndiag_over_msnap*cum_zonal_dTe
   mean_zonal_dTi  = ndiag_over_msnap*cum_zonal_dTi

   mean_zonal_dphi_SQ = inv_msnap*cum_phi00_SQ
   mean_zonal_dni_SQ  = inv_msnap*den_grid*cum_zonali_SQ
   mean_zonal_dne_SQ  = inv_msnap*den_grid*cum_zonale_SQ
   mean_zonal_dTe_SQ  = ndiag_over_msnap*cum_zonal_dTe_SQ
   mean_zonal_dTi_SQ  = ndiag_over_msnap*cum_zonal_dTi_SQ

   mean_full_dphi  = inv_msnap*cum_full_dphi
   mean_full_dne   = inv_msnap*cum_full_dne
   mean_full_dni   = inv_msnap*cum_full_dni
   mean_full_dTe   = ndiag_over_msnap*cum_full_dTe
   mean_full_dTi   = ndiag_over_msnap*cum_full_dTi

   mean_full_dphi_SQ  = inv_msnap*cum_full_dphi_SQ
   mean_full_dne_SQ   = inv_msnap*cum_full_dne_SQ
   mean_full_dni_SQ   = inv_msnap*cum_full_dni_SQ
   mean_full_dTe_SQ   = ndiag_over_msnap*cum_full_dTe_SQ
   mean_full_dTi_SQ  = ndiag_over_msnap*cum_full_dTi_SQ

   do i=0,mpsi
      do j=0,mtheta(i)
         ij=igrid(i) + j
         full_dphi(:,ij) = phi(:,ij) +real(ZF_SEPARATION)*phi00(i)
         if (Flag_Def_dne_adiab==0) then
            intpl_dphi_NZ = phi(:,ij) -real(1.0-ZF_SEPARATION)*phi00(i)
            dne_adiabatic(:,ij) = den_grid(i)*abs(csc0)*intpl_dphi_NZ/Te_grid(i)
         elseif(Flag_Def_dne_adiab==1) then !Flag_Def_dne_adiab==1
            ! adiabatic electron is using full_dphi instead of nonzonal phi
            dne_adiabatic(:,ij) = den_grid(i)*abs(csc0)*(phi(:,ij)+real(ZF_SEPARATION)*phi00(i))/Te_grid(i)
         elseif(Flag_Def_dne_adiab==2) then
            dne_adiabatic(:,ij) = den_grid(i)*abs(csc0)*(phi(:,ij)+real(ZF_SEPARATION)*phi00(i))/FLA_Te0(:,ij)
         end if

         full_dne(:,ij) =  den_grid(i)*(densitye(:,ij)+real(ZF_SEPARATION)*zonale(i)) + int(ifelectron/2)*dne_adiabatic(:,ij)
         full_dni(:,ij) =  den_grid(i)*(densityi(:,ij)+real(ZF_SEPARATION)*zonali(i))
      enddo
   enddo

   !==============================================
   ! Record Information at zeta=0 poloidal plane
   !==============================================
   if((myrank_partd==nproc_partd-1) .and. (mod(istep,msnap) .eq. 0)) then
      ! call Do_FieldLineAverage_by_using_gPPI(full_dphi, FLA_dphi)
      ! call Do_FieldLineAverage_by_using_gPPI(dne_adiabatic, FLA_ne_adiabatic)
      ! call Do_FieldLineAverage_by_using_gPPI(full_dne, FLA_full_dne)
      ! call Do_FieldLineAverage_by_using_gPPI(full_dni, FLA_full_dni)
      call Do_FieldLineAverage_by_using_gPPI(mean_full_dphi, FLA_dphi)
      call Do_FieldLineAverage_by_using_gPPI(dne_adiabatic, FLA_ne_adiabatic)
      call Do_FieldLineAverage_by_using_gPPI(mean_full_dne, FLA_full_dne)
      call Do_FieldLineAverage_by_using_gPPI(mean_full_dni, FLA_full_dni)
      call Do_FieldLineAverage_by_using_gPPI(mean_full_dTe, FLA_dTe)
      call Do_FieldLineAverage_by_using_gPPI(mean_full_dTi, FLA_dTi)
   endif


   !=================================
   ! Calculate electrostatic energy (EE) in the system
   !================================
   tmp_FLA_dphi_nonNaN = FLA_dphi
   where(FLA_Denominator==0.0) tmp_FLA_dphi_nonNaN=0.0 ! Remove NaN elements in FLA_dphi to get nonNaN EE values for FLA and NFLA components

   ! Calculate EE on each poloidal plane

   Device_simple_Vol = sum(Node_simple_Vol)
   EE_tot  = sum(Node_simple_Vol*(mean_full_dni(1,:) +csc0*mean_full_dne(1,:))*mean_full_dphi(1,:))
   EE_FLA  = sum(Node_simple_Vol*(mean_full_dni(1,:) +csc0*mean_full_dne(1,:))*(tmp_FLA_dphi_nonNaN(1,:)))
   EE_NFLA = sum(Node_simple_Vol*(mean_full_dni(1,:) +csc0*mean_full_dne(1,:))*(mean_full_dphi(1,:)-tmp_FLA_dphi_nonNaN(1,:)))
   EE_Z  = 0.0
   EE_NZ = 0.0
   do i=0,mpsi
      do j=0,mtheta(i)
         ij=igrid(i) +j
         EE_Z  = EE_Z  +Node_simple_Vol(ij)*(mean_full_dni(1,ij)+csc0*mean_full_dne(1,ij))*mean_zonal_dphi(i)
         EE_NZ = EE_NZ +Node_simple_Vol(ij)*(mean_full_dni(1,ij)+csc0*mean_full_dne(1,ij))*(mean_full_dphi(1,ij)-mean_zonal_dphi(i))
      enddo
   enddo
   ! Take a summation over toroidal_comm to get overall EE energy
   call MPI_ALLREDUCE(Device_simple_Vol,recv_0D, 1, mpi_Rsize, MPI_SUM, toroidal_comm, ierror)
   Device_simple_Vol = recv_0D
   call MPI_ALLREDUCE(EE_tot,recv_0D, 1, mpi_Rsize, MPI_SUM, toroidal_comm, ierror)
   EE_tot = recv_0D
   call MPI_ALLREDUCE(EE_FLA,recv_0D, 1, mpi_Rsize, MPI_SUM, toroidal_comm, ierror)
   EE_FLA = recv_0D
   call MPI_ALLREDUCE(EE_NFLA,recv_0D, 1, mpi_Rsize, MPI_SUM, toroidal_comm, ierror)
   EE_NFLA = recv_0D
   call MPI_ALLREDUCE(EE_Z,recv_0D, 1, mpi_Rsize, MPI_SUM, toroidal_comm, ierror)
   EE_Z = recv_0D
   call MPI_ALLREDUCE(EE_NZ,recv_0D, 1, mpi_Rsize, MPI_SUM, toroidal_comm, ierror)
   EE_NZ = recv_0D


   !=================================
   ! Calculate Electric fields & mean electric fields
   !================================
   dE_cov3 = -evector ! This is an instanteneous

   allocate(intp_a_1D(dim_a), recvl(dim_gq*dim_a), recvr(dim_gq*dim_a), intp_val_Arr(dim_gq*dim_a))
   allocate(intp_mean_full_dphi(dim_gq, dim_a), right_intp_mean_full_dphi(dim_gq, dim_a), left_intp_mean_full_dphi(dim_gq, dim_a))
   allocate(intp_FLA_dphi(dim_gq, dim_a), right_intp_FLA_dphi(dim_gq, dim_a), left_intp_FLA_dphi(dim_gq, dim_a))
   allocate(intp_mean_dE_cov3(3,dim_gq,dim_a), right_mean_dE_cov3(3,dim_gq,dim_a),left_mean_dE_cov3(3,dim_gq,dim_a))
   allocate(intp_FLA_dE_cov3(3,dim_gq,dim_a), right_FLA_dE_cov3(3,dim_gq,dim_a), left_FLA_dE_cov3(3,dim_gq,dim_a))


   k=1
   do i=0,mpsi,a_jump
      intp_a_1D(k) = psi_grid(i)
      k=k+1
   enddo


   ! Cal interpolated phi
   call Interpolate_Data_into_Coarser_Grids(mean_full_dphi(1,:),gq_jump,a_jump,intp_val_Arr)
   ! send phi to right and receive from left
   call MPI_SENDRECV(intp_val_Arr, dim_gq*dim_a, mpi_Rsize, right_pe, myrank_toroidal, &
      recvl, dim_gq*dim_a, mpi_Rsize, left_pe, left_pe, toroidal_comm, istatus, ierror)
   ! send phi to left and receive from right
   call MPI_SENDRECV(intp_val_Arr, dim_gq*dim_a, mpi_Rsize, left_pe, myrank_toroidal, &
      recvr, dim_gq*dim_a, mpi_Rsize, right_pe, right_pe, toroidal_comm, istatus, ierror)
   intp_mean_full_dphi = reshape(intp_val_Arr,(/dim_gq, dim_a/))
   left_intp_mean_full_dphi = reshape(recvl,(/dim_gq, dim_a/))
   right_intp_mean_full_dphi = reshape(recvr,(/dim_gq, dim_a/))

   ! Cal interpolated FLA phi
   call Interpolate_Data_into_Coarser_Grids(FLA_dphi(1,:),gq_jump,a_jump,intp_val_Arr)
   ! send phi to right and receive from left
   call MPI_SENDRECV(intp_val_Arr, dim_gq*dim_a, mpi_Rsize, right_pe, myrank_toroidal, &
      recvl, dim_gq*dim_a, mpi_Rsize, left_pe, left_pe, toroidal_comm, istatus, ierror)
   ! send phi to left and receive from right
   call MPI_SENDRECV(intp_val_Arr, dim_gq*dim_a, mpi_Rsize, left_pe, myrank_toroidal, &
      recvr, dim_gq*dim_a, mpi_Rsize, right_pe, right_pe, toroidal_comm, istatus, ierror)
   intp_FLA_dphi = reshape(intp_val_Arr,(/dim_gq, dim_a/))
   left_intp_FLA_dphi = reshape(recvl,(/dim_gq, dim_a/))
   right_intp_FLA_dphi = reshape(recvr,(/dim_gq, dim_a/))

   deallocate(recvl, recvr, intp_val_Arr)


   ! Radial electric fields (delta_a may not be constant; so it should be calculated at each )
   intp_mean_dE_cov3(1,:,1) = -(intp_mean_full_dphi(:,2)-intp_mean_full_dphi(:,1))/(intp_a_1D(2)-intp_a_1D(1))
   do i=2,dim_a-1
      coarser_dpsi1 = intp_a_1D(i)-intp_a_1D(i-1)
      coarser_dpsi2 = intp_a_1D(i+1)-intp_a_1D(i)
      intp_mean_dE_cov3(1,:,i) = -(coarser_dpsi2*(intp_mean_full_dphi(:,i)-intp_mean_full_dphi(:,i-1))/coarser_dpsi1 + coarser_dpsi1*(intp_mean_full_dphi(:,i+1)-intp_mean_full_dphi(:,i))/coarser_dpsi2)/(coarser_dpsi1+coarser_dpsi2)
   enddo
   intp_mean_dE_cov3(1,:,dim_a) = -(intp_mean_full_dphi(:,dim_a)-intp_mean_full_dphi(:,dim_a-1))/(intp_a_1D(dim_a)-intp_a_1D(dim_a-1))
   ! theta electric fields
   intp_mean_dE_cov3(2,1,:)    = -(intp_mean_full_dphi(dim_gq-1,:)-intp_mean_full_dphi(2,:))/(2.0*coarser_dgq)
   do j=2,dim_gq-1
      intp_mean_dE_cov3(2,j,:)  = -(intp_mean_full_dphi(j+1,:)-intp_mean_full_dphi(j-1,:))/(2.0*coarser_dgq)
   enddo
   intp_mean_dE_cov3(2,dim_gq,:) = intp_mean_dE_cov3(2,1,:)
   ! toroidal electric fields
   intp_mean_dE_cov3(3,:,:) = -(right_intp_mean_full_dphi-left_intp_mean_full_dphi)/(2.0*deltaz)


   ! Radial electric fields (delta_a may not be constant; so it should be calculated at each )
   intp_FLA_dE_cov3(1,:,1) = -(intp_FLA_dphi(:,2)-intp_FLA_dphi(:,1))/(intp_a_1D(2)-intp_a_1D(1))
   do i=2,dim_a-1
      coarser_dpsi1 = intp_a_1D(i)-intp_a_1D(i-1)
      coarser_dpsi2 = intp_a_1D(i+1)-intp_a_1D(i)
      intp_FLA_dE_cov3(1,:,i) = -(coarser_dpsi2*(intp_FLA_dphi(:,i)-intp_FLA_dphi(:,i-1))/coarser_dpsi1 + coarser_dpsi1*(intp_FLA_dphi(:,i+1)-intp_FLA_dphi(:,i))/coarser_dpsi2)/(coarser_dpsi1+coarser_dpsi2)
   enddo
   intp_FLA_dE_cov3(1,:,dim_a) = -(intp_FLA_dphi(:,dim_a)-intp_FLA_dphi(:,dim_a-1))/(intp_a_1D(dim_a)-intp_a_1D(dim_a-1))
   ! theta electric fields
   intp_FLA_dE_cov3(2,1,:)    = -(intp_FLA_dphi(dim_gq-1,:)-intp_FLA_dphi(2,:))/(2.0*coarser_dgq)
   do j=2,dim_gq-1
      intp_FLA_dE_cov3(2,j,:)  = -(intp_FLA_dphi(j+1,:)-intp_FLA_dphi(j-1,:))/(2.0*coarser_dgq)
   enddo
   intp_FLA_dE_cov3(2,dim_gq,:) = intp_FLA_dE_cov3(2,1,:)
   ! toroidal electric fields
   intp_FLA_dE_cov3(3,:,:) = -(right_intp_FLA_dphi-left_intp_FLA_dphi)/(2.0*deltaz)

   ! Send&Recv intp_mean_dE_cov3 b.t.w. adjacent polodial planes
   call MPI_SENDRECV(intp_mean_dE_cov3, 3*dim_gq*dim_a, mpi_Rsize, right_pe, myrank_toroidal, &
      left_mean_dE_cov3, 3*dim_gq*dim_a, mpi_Rsize, left_pe, left_pe, toroidal_comm, istatus, ierror)
   call MPI_SENDRECV(intp_mean_dE_cov3, 3*dim_gq*dim_a, mpi_Rsize, left_pe, myrank_toroidal, &
      right_mean_dE_cov3, 3*dim_gq*dim_a, mpi_Rsize, right_pe, right_pe, toroidal_comm, istatus, ierror)
   ! Send&Recv intp_FLA_dE_cov3 b.t.w. adjacent polodial planes
   call MPI_SENDRECV(intp_FLA_dE_cov3, 3*dim_gq*dim_a, mpi_Rsize, right_pe, myrank_toroidal, &
      left_FLA_dE_cov3, 3*dim_gq*dim_a, mpi_Rsize, left_pe, left_pe, toroidal_comm, istatus, ierror)
   call MPI_SENDRECV(intp_FLA_dE_cov3, 3*dim_gq*dim_a, mpi_Rsize, left_pe, myrank_toroidal, &
      right_FLA_dE_cov3, 3*dim_gq*dim_a, mpi_Rsize, right_pe, right_pe, toroidal_comm, istatus, ierror)


   ! Gather mean_full_dphi into test_dphi_3D to record full 3D potential
   if(myrank_partd==nproc_partd-1) then
      h5dims(1:2)=(/dim_gq, dim_a/)
      allocate(intp_val_Arr(dim_a*dim_gq))
      allocate(test_dphi_3D(dim_a*dim_gq*mzetamax))

      ! (temporal) Mean full_dphi
      call Interpolate_Data_into_Coarser_Grids(mean_full_dphi(1,:),gq_jump,a_jump,intp_val_Arr)
      call MPI_ALLGATHER(intp_val_Arr,dim_a*dim_gq,mpi_Rsize,test_dphi_3D,dim_a*dim_gq,mpi_Rsize,toroidal_comm,h5err)
      deallocate(intp_val_Arr)
   endif




   ! ==================================
   ! Record some data on sparse grid points
   ! ==================================
   snap_flag=.false.
   if(mod(istep,msnap_for_All)==0 .and. myrank_partd==nproc_partd-1) then
      call system("mkdir -p " // 'snapshots_hdf5/AllPlanes')
      write(file_path, '("snapshots_hdf5/AllPlanes/P",i3.3,"_snap",i5.5,".h5")') myrank_toroidal, istep/msnap
      snap_flag=.true.
   endif
   if(mype == numberpe-1)then !mype=0 must also have myrank_toroidal=0
      call system("mkdir -p " // 'snapshots_hdf5')
      write(file_path, '("snapshots_hdf5/snap",i5.5,".h5")') istep/msnap
      !open(snapout,file=file_path,status='replace')
      snap_flag=.true.
   endif

   if(snap_flag) then
      call h5open_f(h5err)
      call h5fcreate_f(file_path, H5F_ACC_TRUNC_F, file_id, h5err)

      call h5gcreate_f(file_id, "PARAMS", gid_PARAMS, h5err)
      call h5gcreate_f(file_id, "FLA", gid_FLA, h5err)
      call h5gcreate_f(file_id, "DCOMP", gid_DCOMP, h5err)
      call h5gcreate_f(file_id, "BIN", gid_BIN, h5err)
      call h5gcreate_f(file_id, "SG", gid_SG, h5err) ! Diagnostic Grid
      call h5gcreate_f(file_id, "SG/SM", gid_SG_SM, h5err) ! Diagnostic Grid
      call h5gcreate_f(file_id, "DG", gid_DG, h5err) ! Diagnostic Grid
      call h5gcreate_f(file_id, "DG/SM", gid_DG_SM, h5err) ! Diagnostic Grid
      call h5gcreate_f(file_id, "DG/FLA", gid_DG_FLA, h5err) ! Diagnostic Grid
      call h5gcreate_f(file_id, "DG/FLA/SM", gid_DG_FLA_SM, h5err) ! Diagnostic Grid
      call h5gcreate_f(file_id, "FP", gid_FP, h5err) ! Diagnostic Grid
      call h5gcreate_f(file_id, "FP/cum", gid_FP_cum, h5err) ! Diagnostic Grid

      !
      h5dims(1) = 1
      call h5ltmake_dataset_f(file_id, "time_s",1,h5dims,H5T_NATIVE_DOUBLE, real(istep)*tstep*utime,h5err)
      call h5ltset_attribute_string_f(file_id, "time_s", "unit", "sec", h5err)
      call h5ltmake_dataset_f(file_id, "norm_time",1,h5dims,H5T_NATIVE_DOUBLE, real(istep)*tstep/MM_Lc_over_vthci,h5err)
      call h5ltset_attribute_string_f(file_id, "norm_time", "unit", "L_T/v_th", h5err)


      ! PARAMS: parameters about the simulation and diagnostic
      h5dims(1) = 1
      call h5ltmake_dataset_f(file_id, "PARAMS/zetamax",1,h5dims,H5T_NATIVE_DOUBLE, zetamax,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/zetamax", "unit", "rad", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/utime_s",1,h5dims,H5T_NATIVE_DOUBLE, utime,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/utime_s", "unit", "sec", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/dt_s",1,h5dims,H5T_NATIVE_DOUBLE, tstep*utime,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/dt_s", "unit", "sec", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/norm_dt",1,h5dims,H5T_NATIVE_DOUBLE, tstep/MM_Lc_over_vthci,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/norm_dt", "unit", "L_T/v_th", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/mpsi",1,h5dims,H5T_NATIVE_INTEGER, mpsi/a_jump, h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/mtheta",1,h5dims,H5T_NATIVE_INTEGER, jm, h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/mzetamax",1,h5dims,H5T_NATIVE_INTEGER, mzetamax, h5err)

      call h5ltmake_dataset_f(file_id, "PARAMS/Lc_over_vthci",1,h5dims,H5T_NATIVE_DOUBLE, MM_Lc_over_vthci*utime,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/Lc_over_vthci", "unit", "sec", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/Cs",1,h5dims,H5T_NATIVE_DOUBLE, MM_Cs/utime,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/Cs", "unit", "m/s", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/vthci",1,h5dims,H5T_NATIVE_DOUBLE, MM_vthci/utime,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/vthci", "unit", "m/s", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/vthce",1,h5dims,H5T_NATIVE_DOUBLE, MM_vthce/utime,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/vthce", "unit", "m/s", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/LcTi",1,h5dims,H5T_NATIVE_DOUBLE, MM_LcTi,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/LcTi", "unit", "m", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/LcTe",1,h5dims,H5T_NATIVE_DOUBLE, MM_LcTe,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/LcTe", "unit", "m", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/Lcn",1,h5dims,H5T_NATIVE_DOUBLE, MM_Lcn,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/Lcn", "unit", "m", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/LcU",1,h5dims,H5T_NATIVE_DOUBLE, MM_LcU,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/LcU", "unit", "m", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/axis_R",1,h5dims,H5T_NATIVE_DOUBLE, MM_R_axis,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/axis_R", "unit", "m", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/axis_B",1,h5dims,H5T_NATIVE_DOUBLE, MM_B_axis,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/axis_B", "unit", "T", h5err)
      call h5ltmake_dataset_f(file_id, "PARAMS/Device_Simple_Vol",1,h5dims,H5T_NATIVE_DOUBLE, Device_simple_Vol,h5err)
      call h5ltset_attribute_string_f(file_id, "PARAMS/Device_Simple_Vol", "unit", "m^3", h5err)




      ! 0D Data
      h5dims(1) = 1
      Unit_Conv_Coeff = (charg0*10/c)*1e6*((E_00*c)/(charg0*1e8))  ! e*n*phi =  (Charge)*(Density)*(Volatege) Conversion
      call h5ltmake_dataset_f(file_id, "EE_tot",1,h5dims,H5T_NATIVE_DOUBLE, Unit_Conv_Coeff*EE_tot, h5err)
      call h5ltset_attribute_string_f(file_id, "EE_tot", "unit", "J", h5err)
      call h5ltmake_dataset_f(file_id, "EE_Z",1,h5dims,H5T_NATIVE_DOUBLE, Unit_Conv_Coeff*EE_Z, h5err)
      call h5ltset_attribute_string_f(file_id, "EE_Z", "unit", "J", h5err)
      call h5ltmake_dataset_f(file_id, "EE_NZ",1,h5dims,H5T_NATIVE_DOUBLE, Unit_Conv_Coeff*EE_NZ, h5err)
      call h5ltset_attribute_string_f(file_id, "EE_NZ", "unit", "J", h5err)
      call h5ltmake_dataset_f(file_id, "EE_FLA",1,h5dims,H5T_NATIVE_DOUBLE, Unit_Conv_Coeff*EE_FLA, h5err)
      call h5ltset_attribute_string_f(file_id, "EE_FLA", "unit", "J", h5err)
      call h5ltmake_dataset_f(file_id, "EE_NFLA",1,h5dims,H5T_NATIVE_DOUBLE, Unit_Conv_Coeff*EE_NFLA, h5err)
      call h5ltset_attribute_string_f(file_id, "EE_NFLA", "unit", "J", h5err)


      ! EQ: About the equilibrium information
      ! integer(HSIZE_T), dimension(2) :: h5dims=(/dim_a, dim_gq/))
      h5dims(1) = dim_a
      call h5ltmake_dataset_f(file_id, "a_1D",1,h5dims,H5T_NATIVE_DOUBLE, psi_grid(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "a_1D", "unit", "norm_psi_t", h5err)
      call h5ltmake_dataset_f(file_id, "n0_1D",1,h5dims,H5T_NATIVE_DOUBLE, 1e6*den_grid(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "n0_1D", "unit", "#/m^3", h5err)
      call h5ltmake_dataset_f(file_id, "Ti0_1D",1,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(vth_grid(0:mpsi:a_jump)**2), h5err)
      call h5ltset_attribute_string_f(file_id, "Ti0_1D", "unit", "eV", h5err)
      call h5ltmake_dataset_f(file_id, "Te0_1D",1,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*Te_grid(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "Te0_1D", "unit", "eV", h5err)
      call h5ltmake_dataset_f(file_id, "zonal_dphi",1,h5dims,H5T_NATIVE_DOUBLE, ((E_00*c)/(charg0*1e8))*phi00(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "zonal_dphi", "unit", "V", h5err)
      call h5ltmake_dataset_f(file_id, "zonal_ddphidr",1,h5dims,H5T_NATIVE_DOUBLE, ((E_00*c)/(charg0*1e8))*phip00(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "zonal_ddphidr", "unit", "V/(norm_psi_t)", h5err)
      call h5ltmake_dataset_f(file_id, "zonal_dni",1,h5dims,H5T_NATIVE_DOUBLE, 1e6*den_grid(0:mpsi:a_jump)*zonali(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "zonal_dni", "unit", "#/m^3", h5err)
      call h5ltmake_dataset_f(file_id, "zonal_dne",1,h5dims,H5T_NATIVE_DOUBLE, 1e6*den_grid(0:mpsi:a_jump)*zonale(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "zonal_dne", "unit", "#/m^3", h5err)
      call h5ltmake_dataset_f(file_id, "zonal_dTi",1,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*zonal_dTi(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "zonal_dTi", "unit", "eV", h5err)
      call h5ltmake_dataset_f(file_id, "zonal_dTe",1,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*zonal_dTe(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "zonal_dTe", "unit", "eV", h5err)

      ! 1D mean zonal quantities
      call h5ltmake_dataset_f(file_id, "mean_zonal_dphi",1,h5dims,H5T_NATIVE_DOUBLE, ((E_00*c)/(charg0*1e8))*mean_zonal_dphi(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "mean_zonal_dphi", "unit", "V", h5err)
      call h5ltmake_dataset_f(file_id, "mean_zonal_dni",1,h5dims,H5T_NATIVE_DOUBLE, 1e6*mean_zonal_dni(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "mean_zonal_dni", "unit", "#/m^3", h5err)
      call h5ltmake_dataset_f(file_id, "mean_zonal_dne",1,h5dims,H5T_NATIVE_DOUBLE, 1e6*mean_zonal_dne(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "mean_zonal_dne", "unit", "#/m^3", h5err)
      call h5ltmake_dataset_f(file_id, "mean_zonal_dTe",1,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*mean_zonal_dTe(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "mean_zonal_dTe", "unit", "eV", h5err)
      call h5ltmake_dataset_f(file_id, "mean_zonal_dTi",1,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*mean_zonal_dTi(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "mean_zonal_dTi", "unit", "eV", h5err)

      call h5ltmake_dataset_f(file_id, "mean_zonal_dphi_SQ",1,h5dims,H5T_NATIVE_DOUBLE, ((E_00*c)/(charg0*1e8))**2.0*mean_zonal_dphi_SQ(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "mean_zonal_dphi_SQ", "unit", "V", h5err)
      call h5ltmake_dataset_f(file_id, "mean_zonal_dni_SQ",1,h5dims,H5T_NATIVE_DOUBLE, (1e6)**2.0*mean_zonal_dni_SQ(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "mean_zonal_dni_SQ", "unit", "#/m^3", h5err)
      call h5ltmake_dataset_f(file_id, "mean_zonal_dne_SQ",1,h5dims,H5T_NATIVE_DOUBLE, (1e6)**2.0*mean_zonal_dne_SQ(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "mean_zonal_dne_SQ", "unit", "#/m^3", h5err)
      call h5ltmake_dataset_f(file_id, "mean_zonal_dTi_SQ",1,h5dims,H5T_NATIVE_DOUBLE, (1e3)**2.0*Tem_00*mean_zonal_dTi_SQ(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "mean_zonal_dTi_SQ", "unit", "eV", h5err)
      call h5ltmake_dataset_f(file_id, "mean_zonal_dTe_SQ",1,h5dims,H5T_NATIVE_DOUBLE, (1e3)**2.0*Tem_00*mean_zonal_dTe_SQ(0:mpsi:a_jump), h5err)
      call h5ltset_attribute_string_f(file_id, "mean_zonal_dTe_SQ", "unit", "eV", h5err)

      ! Some information about diagnostic BINs
      h5dims(1) = 1
      call h5ltmake_dataset_f(file_id, "BIN/nbin",1,h5dims,H5T_NATIVE_INTEGER, nbin_a, h5err)
      h5dims(1) = nbin_a
      call h5ltmake_dataset_f(file_id, "BIN/sn0",1,h5dims,H5T_NATIVE_DOUBLE, sn0, h5err)
      call h5ltmake_dataset_f(file_id, "BIN/dnty_on_bin",1,h5dims,H5T_NATIVE_DOUBLE, dnty_on_bin, h5err)
      call h5ltmake_dataset_f(file_id, "BIN/Tr",1,h5dims,H5T_NATIVE_DOUBLE, Tr, h5err)
      call h5ltmake_dataset_f(file_id, "BIN/Tr_e",1,h5dims,H5T_NATIVE_DOUBLE, Tr_e, h5err)
      call h5ltmake_dataset_f(file_id, "BIN/dnty_r",1,h5dims,H5T_NATIVE_DOUBLE, dnty_r, h5err)
      call h5ltmake_dataset_f(file_id, "BIN/dO_10",1,h5dims,H5T_NATIVE_DOUBLE, dO_10, h5err)
      call h5ltmake_dataset_f(file_id, "BIN/dT_10",1,h5dims,H5T_NATIVE_DOUBLE, dT_10, h5err)

      h5dims(1)=1
      call h5ltmake_dataset_f(file_id, "num_of_FLabel_pcws",1,h5dims,H5T_NATIVE_INTEGER, num_of_FLabel_pcws,h5err)
      h5dims(1)=num_of_FLabel_pcws
      call h5ltmake_dataset_f(file_id, "FLabel_pcws_a",1,h5dims,H5T_NATIVE_DOUBLE, FLabel_pcws_a,h5err)
      call h5ltmake_dataset_f(file_id, "FLabel_pcws_y",1,h5dims,H5T_NATIVE_DOUBLE, FLabel_pcws_y,h5err)

      ! =======================
      ! 2D data
      ! =======================
      h5dims(1:2)=(/dim_gq, dim_a/)
      allocate(intp_val_Arr(dim_a*dim_gq))
      allocate(intp_R_2D(dim_a*dim_gq), intp_Z_2D(dim_a*dim_gq), intp_x_2D(dim_a*dim_gq), intp_y_2D(dim_a*dim_gq))

      k=1
      do i=0,mpsi,a_jump
         do j=0,jm
            aaa(1)=psi_grid(i)
            the(1)=2.0*pi*real(j)/real(jm)
            call esigetrzb(rj,zj,Bj,aaa,the,1)
            intp_R_2D(k) = rj(1)
            intp_Z_2D(k) = zj(1)
            intp_x_2D(k) = rj(1)*cos(zetamax)
            intp_y_2D(k) = rj(1)*sin(zetamax)
            k=k+1
         enddo
      enddo

      ! R_2D
      call h5ltmake_dataset_f(file_id, "R_2D",2,h5dims,H5T_NATIVE_DOUBLE, intp_R_2D,h5err)
      call h5ltset_attribute_string_f(file_id, "R_2D", "unit", "m", h5err)
      ! Z_2D
      call h5ltmake_dataset_f(file_id, "Z_2D",2,h5dims,H5T_NATIVE_DOUBLE, intp_Z_2D,h5err)
      call h5ltset_attribute_string_f(file_id, "Z_2D", "unit", "m", h5err)
      ! x_2D
      call h5ltmake_dataset_f(file_id, "x_2D",2,h5dims,H5T_NATIVE_DOUBLE, intp_x_2D,h5err)
      call h5ltset_attribute_string_f(file_id, "x_2D", "unit", "m", h5err)
      ! y_2D
      call h5ltmake_dataset_f(file_id, "y_2D",2,h5dims,H5T_NATIVE_DOUBLE, intp_y_2D,h5err)
      call h5ltset_attribute_string_f(file_id, "y_2D", "unit", "m", h5err)

      deallocate(intp_R_2D, intp_Z_2D, intp_x_2D, intp_y_2D)


      ! a_2D
      k=1
      do i=0,mpsi,a_jump
         do j=0,jm
            intp_val_Arr(k) = psi_grid(i)
            k=k+1
         enddo
      enddo
      call h5ltmake_dataset_f(file_id, "a_2D",2,h5dims,H5T_NATIVE_DOUBLE, intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "a_2D", "unit", "norm_psi_t", h5err)

      ! gq_2D
      k=1
      do i=0,mpsi,a_jump
         do j=0,jm
            intp_val_Arr(k) = 2.0*pi*real(j)/real(jm)
            k=k+1
         enddo
      enddo
      call h5ltmake_dataset_f(file_id, "gq_2D",2,h5dims,H5T_NATIVE_DOUBLE, intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "gq_2D", "unit", "radians", h5err)


      call Interpolate_Data_into_Coarser_Grids(Node_simple_Vol,gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "ori_Node_simple_Vol",2,h5dims,H5T_NATIVE_DOUBLE, intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "ori_Node_simple_Vol", "unit", "m^3", h5err)

      call Interpolate_Data_into_Coarser_Grids(FLabel_Function(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLabel_Function",2,h5dims,H5T_NATIVE_DOUBLE, intp_val_Arr,h5err)




      ! dphi_2D
      k=1
      do i=0,mpsi,a_jump
         do j=0,jm
            gqdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(gqdum);  wt=gqdum-jt
            intp_val = ((1.0-wt)*phi(1,igrid(i)+jt)+  wt*phi(1,igrid(i)+jt+1)) +real(ZF_SEPARATION)*phi00(i)
            intp_val_Arr(k) = ((E_00*c)/(charg0*1e8))*intp_val
            k=k+1
         enddo
      enddo
      call h5ltmake_dataset_f(file_id, "dphi",2,h5dims,H5T_NATIVE_DOUBLE, intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "dphi", "unit", "V", h5err)

      call Interpolate_Data_into_Coarser_Grids(phi(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "test_dphi",2,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "test_dphi", "unit", "V", h5err)


      ! ion density
      ! from tot_dni with no enforced charge-conservation
      call Interpolate_Data_into_Coarser_Grids(tot_dni(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "dni_from_tot_dni",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "dni_from_tot_dni", "unit", "#/m^3", h5err)


      ! from densityi directly
      ! If(its_jacobi=1). this is the result enforced by charge-conservation
      call Interpolate_Data_into_Coarser_Grids(full_dni(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "full_dni",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "full_dni", "unit", "#/m^3", h5err)

      ! ion temperature fluctuation
      k=1
      do i=0,mpsi,a_jump
         do j=0,jm
            gqdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(gqdum);  wt=gqdum-jt
            intp_val = real(ZF_SEPARATION)*(zonal_dTi(i)) + (vth_grid(i)**2)*((1.0-wt)*phi_tmp(igrid(i)+jt)+  wt*phi_tmp(igrid(i)+jt+1))
            intp_val_Arr(k) = 1e3*Tem_00*intp_val
            k=k+1
         enddo
      enddo
      call h5ltmake_dataset_f(file_id, "dTi",2,h5dims,H5T_NATIVE_DOUBLE, intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "dTi", "unit", "eV", h5err)

      ! electron density
      ! from tot_dne with no enforced charge-conservation
      call Interpolate_Data_into_Coarser_Grids(tot_dne(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "dne_from_tot_dne",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "dne_from_tot_dne", "unit", "#/m^3", h5err)


      ! from densitye directly
      ! If(its_jacobi=1). this is the result enforced by charge-conservation
      k=1
      do i=0,mpsi,a_jump
         do j=0,jm
            gqdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(gqdum);  wt=gqdum-jt
            intp_val = (1.0-wt)*full_dne(1,igrid(i)+jt) + wt*full_dne(1,igrid(i)+jt+1)
            intp_val_Arr(k) = 1e6*intp_val
            k=k+1
         enddo
      enddo
      call Interpolate_Data_into_Coarser_Grids(full_dne(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "full_dne",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "full_dne", "unit", "#/m^3", h5err)

      ! Adiabatic density calculated from nonzonal dphi
      call Interpolate_Data_into_Coarser_Grids(dne_adiabatic(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "dne_adiabatic",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "dne_adiabatic", "unit", "#/m^3", h5err)

      ! Electron temperature fluctuation
      k=1
      do i=0,mpsi,a_jump
         do j=0,jm
            gqdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(gqdum);  wt=gqdum-jt
            intp_val = Te_grid(i)*(((1.0-wt)*heatflux(igrid(i)+jt)+  wt*heatflux(igrid(i)+jt+1))) +real(ZF_SEPARATION)*(zonal_dTe(i))
            intp_val_Arr(k) = 1e3*Tem_00*(intp_val)
            k=k+1
         enddo
      enddo
      call h5ltmake_dataset_f(file_id, "dTe",2,h5dims,H5T_NATIVE_DOUBLE, intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "dTe", "unit", "eV", h5err)


      ! (temporal) Mean full_dni
      call Interpolate_Data_into_Coarser_Grids(mean_full_dni(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "mean_full_dni",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dni", "unit", "#/m^3", h5err)

      ! (temporal) Mean full_dne
      call Interpolate_Data_into_Coarser_Grids(mean_full_dne(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "mean_full_dne",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dne", "unit", "#/m^3", h5err)

      ! (temporal) Mean full_dphi
      call Interpolate_Data_into_Coarser_Grids(mean_full_dphi(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "mean_full_dphi",2,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dphi", "unit", "V", h5err)

      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid(i)+j
            tmp_vec_mgrid_size(ij) = test_src(ij)*density(psi_grid(i),0)
         enddo
      enddo

      call Interpolate_Data_into_Coarser_Grids(tmp_vec_mgrid_size,gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "test_src",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "test_src", "unit", "/#^3", h5err)

      if(ifelectron==1) then
         tmp_vec_mgrid_size=0.0
         do i=0,mpsi
            tmp_Te0_inverse = 1.0/temperature_e(psi_grid(i),0.0_wp)
            do j=1,mtheta(i)
               ij=igrid(i)+j
               if(adiab > 0.0_wp) then
                  tmp_vec_mgrid_size(ij) = adiab * abs(csc0) * m0ms * tmp_Te0_inverse * (phi(1, ij)-old_phi(1,ij))
               else
                  tmp_vec_mgrid_size(ij) = 0.01/(Cs*Cs) * (phi(1,ij)-old_phi(1,ij))
               endif
               tmp_vec_mgrid_size(ij)=tmp_vec_mgrid_size(ij)*density(psi_grid(i),0)
            enddo
         enddo
         do i=0, mpsi
            tmp_vec_mgrid_size(igrid(i))=tmp_vec_mgrid_size(igrid(i)+mtheta(i))  ! which is missing in last do-loop
         enddo

         call Interpolate_Data_into_Coarser_Grids(old_phi(1,:),gq_jump,a_jump,intp_val_Arr)
         call h5ltmake_dataset_f(file_id, "old_dphi",2,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*intp_val_Arr,h5err)
         call h5ltset_attribute_string_f(file_id, "old_dphi", "unit", "V", h5err)

         call Interpolate_Data_into_Coarser_Grids(tmp_vec_mgrid_size,gq_jump,a_jump,intp_val_Arr)
         call h5ltmake_dataset_f(file_id, "missing_src",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
         call h5ltset_attribute_string_f(file_id, "missing_src", "unit", "/#^3", h5err)
      endif

      h5dims(1:3)=(/dim_gq, dim_a, mzetamax/)
      call h5ltmake_dataset_f(file_id, "mean_full_dphi_3D",3,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*test_dphi_3D,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dphi_3D", "unit", "V", h5err)

      ! Mean full_dphi of Right poloidal plane
      call h5ltmake_dataset_f(file_id, "mean_full_dphi_Right",2,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*right_intp_mean_full_dphi,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dphi_Right", "unit", "V", h5err)

      call h5ltmake_dataset_f(file_id, "mean_full_dphi_Left",2,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*left_intp_mean_full_dphi,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dphi_Left", "unit", "V", h5err)


      ! (temporal) Mean dTe
      call Interpolate_Data_into_Coarser_Grids(mean_full_dTe(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "mean_full_dTe",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dTe", "unit", "eV", h5err)

      ! (temporal) Mean dTi
      call Interpolate_Data_into_Coarser_Grids(mean_full_dTi(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "mean_full_dTi",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dTi", "unit", "eV", h5err)


      ! (temporal) Mean of full dphi_SQ
      call Interpolate_Data_into_Coarser_Grids(mean_full_dphi_SQ(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "mean_full_dphi_SQ",2,h5dims,H5T_NATIVE_DOUBLE, ((E_00*c)/(charg0*1e8))**2.0*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dphi_SQ", "unit", "V", h5err)

      ! (temporal) Mean of full dni_SQ
      call Interpolate_Data_into_Coarser_Grids(mean_full_dni_SQ(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "mean_full_dni_SQ",2,h5dims,H5T_NATIVE_DOUBLE, (1e6)**2.0*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dni_SQ", "unit", "#/m^3", h5err)

      ! (temporal) Mean of full dne_SQ
      call Interpolate_Data_into_Coarser_Grids(mean_full_dne_SQ(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "mean_full_dne_SQ",2,h5dims,H5T_NATIVE_DOUBLE, (1e6)**2.0*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dne_SQ", "unit", "#/m^3", h5err)

      ! (temporal) Mean of full dTi_SQ
      call Interpolate_Data_into_Coarser_Grids(mean_full_dTi_SQ(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "mean_full_dTi_SQ",2,h5dims,H5T_NATIVE_DOUBLE,  (1e3*Tem_00)**2.0*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dTi_SQ", "unit", "eV", h5err)

      ! (temporal) Mean of full dTe_SQ
      call Interpolate_Data_into_Coarser_Grids(mean_full_dTe_SQ(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "mean_full_dTe_SQ",2,h5dims,H5T_NATIVE_DOUBLE, (1e3*Tem_00)**2.0*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_full_dTe_SQ", "unit", "eV", h5err)


      ! FLA FLabel (Field Fid)
      call Interpolate_Data_into_Coarser_Grids(FLA_FLabel(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/FLabel",2,h5dims,H5T_NATIVE_DOUBLE, intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/FLabel", "unit", "", h5err)

      ! FLA psi0
      call Interpolate_Data_into_Coarser_Grids(FLA_psi0(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/psi0",2,h5dims,H5T_NATIVE_DOUBLE, intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/psi0", "unit", "norm_psi_t", h5err)

      ! FLA Ti0
      call Interpolate_Data_into_Coarser_Grids(FLA_Ti0(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/Ti0",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/Ti0", "unit", "eV", h5err)

      ! FLA Te0
      call Interpolate_Data_into_Coarser_Grids(FLA_Te0(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/Te0",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/Te0", "unit", "eV", h5err)

      ! FLA Te0_Left
      call Interpolate_Data_into_Coarser_Grids(FLA_Te0(0,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/Te0_Left",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/Te0_Left", "unit", "eV", h5err)


      ! FLA n0
      call Interpolate_Data_into_Coarser_Grids(FLA_n0(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/n0",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/n0", "unit", "#/m^3", h5err)

      ! FLA dphi
      call Interpolate_Data_into_Coarser_Grids(FLA_dphi(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/dphi",2,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/dphi", "unit", "V", h5err)

      ! FLA dphi of Right poloidal plane
      call h5ltmake_dataset_f(file_id, "FLA/dphi_Right",2,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*right_intp_FLA_dphi,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/dphi_Right", "unit", "V", h5err)
      ! FLA dphi of Left poloidal plane
      call h5ltmake_dataset_f(file_id, "FLA/dphi_Left",2,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*left_intp_FLA_dphi,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/dphi_Left", "unit", "V", h5err)


      ! FLA dTi
      call Interpolate_Data_into_Coarser_Grids(FLA_dTi(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/dTi",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/dTi", "unit", "eV", h5err)

      ! FLA dni
      call Interpolate_Data_into_Coarser_Grids(FLA_dni(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/dni_from_tot_dni",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/dni_from_tot_dni", "unit", "#/m^3", h5err)

      ! FLA dni
      call Interpolate_Data_into_Coarser_Grids(FLA_full_dni(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/full_dni",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/full_dni", "unit", "#/m^3", h5err)

      ! FLA dTe
      call Interpolate_Data_into_Coarser_Grids(FLA_dTe(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/dTe",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/dTe", "unit", "eV", h5err)

      ! FLA_dne
      call Interpolate_Data_into_Coarser_Grids(FLA_dne(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/dne_from_tot_dne",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/dne_from_tot_dne", "unit", "#/m^3", h5err)

      ! FLA dne
      call Interpolate_Data_into_Coarser_Grids(FLA_full_dne(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/full_dne",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/full_dne", "unit", "#/m^3", h5err)


      ! FLA_ne_adiabatic
      call Interpolate_Data_into_Coarser_Grids(FLA_dne(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "FLA/dne_adiabatic",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/dne_adiabatic", "unit", "#/m^3", h5err)



      ! Decomposition
      ! dni each components
      call Interpolate_Data_into_Coarser_Grids(dni_init(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dni_init",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dni_init", "unit", "#/m^3", h5err)

      call Interpolate_Data_into_Coarser_Grids(dni_vpara(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dni_vpara",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dni_vpara", "unit", "#/m^3", h5err)

      call Interpolate_Data_into_Coarser_Grids(dni_vperp(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dni_vperp",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dni_vperp", "unit", "#/m^3", h5err)

      call Interpolate_Data_into_Coarser_Grids(dni_accel(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dni_accel",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dni_accel", "unit", "#/m^3", h5err)

      call Interpolate_Data_into_Coarser_Grids(dni_fdamp(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dni_fdamp",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dni_fdamp", "unit", "#/m^3", h5err)

      ! dTi each components
      call Interpolate_Data_into_Coarser_Grids(dTi_init(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dTi_init",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dTi_init", "unit", "eV", h5err)

      call Interpolate_Data_into_Coarser_Grids(dTi_vpara(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dTi_vpara",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dTi_vpara", "unit", "eV", h5err)

      call Interpolate_Data_into_Coarser_Grids(dTi_vperp(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dTi_vperp",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dTi_vperp", "unit", "eV", h5err)

      call Interpolate_Data_into_Coarser_Grids(dTi_fdamp(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dTi_fdamp",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dTi_fdamp", "unit", "eV", h5err)

      ! dne each components
      call Interpolate_Data_into_Coarser_Grids(dne_init(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dne_init",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dne_init", "unit", "#/m^3", h5err)

      call Interpolate_Data_into_Coarser_Grids(dne_vpara(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dne_vpara",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dne_vpara", "unit", "#/m^3", h5err)

      call Interpolate_Data_into_Coarser_Grids(dne_vperp(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dne_vperp",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dne_vperp", "unit", "#/m^3", h5err)

      call Interpolate_Data_into_Coarser_Grids(dne_accel(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dne_accel",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dne_accel", "unit", "#/m^3", h5err)

      call Interpolate_Data_into_Coarser_Grids(dne_adiab(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dne_dADdt",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dne_dADdt", "unit", "#/m^3", h5err)

      call Interpolate_Data_into_Coarser_Grids(dne_fdamp(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dne_fdamp",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dne_fdamp", "unit", "#/m^3", h5err)

      ! dTe each components
      call Interpolate_Data_into_Coarser_Grids(dTe_init(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dTe_init",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dTe_init", "unit", "eV", h5err)

      call Interpolate_Data_into_Coarser_Grids(dTe_vpara(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dTe_vpara",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dTe_vpara", "unit", "eV", h5err)

      call Interpolate_Data_into_Coarser_Grids(dTe_vperp(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dTe_vperp",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dTe_vperp", "unit", "eV", h5err)

      call Interpolate_Data_into_Coarser_Grids(dTe_accel(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dTe_accel",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dTe_accel", "unit", "eV", h5err)

      call Interpolate_Data_into_Coarser_Grids(dTe_adiab(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dTe_dADdt",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dTe_dADdt", "unit", "eV", h5err)

      call Interpolate_Data_into_Coarser_Grids(dTe_fdamp(1,:),gq_jump,a_jump,intp_val_Arr)
      call h5ltmake_dataset_f(file_id, "DCOMP/dTe_fdamp",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "DCOMP/dTe_fdamp", "unit", "eV", h5err)

      ! 2D data record finished
      deallocate(intp_val_Arr)




      !===========================
      !3D data  (vector data on 2D grids)
      !==========================

      h5dims(1:3)=(/3,dim_gq,dim_a/)
      allocate(intp_val_Arr(3*dim_gq*dim_a))
      call Interpolate_Data_into_Coarser_Grids(dE_cov3(1,1,:),gq_jump,a_jump,intp_val_Arr(1:3*dim_gq*dim_a:3))
      call Interpolate_Data_into_Coarser_Grids(dE_cov3(2,1,:),gq_jump,a_jump,intp_val_Arr(2:3*dim_gq*dim_a:3))
      call Interpolate_Data_into_Coarser_Grids(dE_cov3(3,1,:),gq_jump,a_jump,intp_val_Arr(3:3*dim_gq*dim_a:3))
      call h5ltmake_dataset_f(file_id, "dE_cov",3,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*intp_val_Arr,h5err)
      call h5ltset_attribute_string_f(file_id, "dE_cov", "unit", "V/L", h5err)
      call h5ltset_attribute_string_f(file_id, "dE_cov", "description", "This is -grad(phi) in PEST coordinate. The length unit is parametric length in PEST coordinate", h5err)

      call h5ltmake_dataset_f(file_id, "mean_dE_cov",3,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*intp_mean_dE_cov3,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_dE_cov", "unit", "V/L", h5err)
      call h5ltset_attribute_string_f(file_id, "mean_dE_cov", "description", "This is -grad(phi) in PEST coordinate. The length unit is parametric length in PEST coordinate", h5err)
      call h5ltmake_dataset_f(file_id, "mean_dE_cov_Right",3,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*right_mean_dE_cov3,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_dE_cov_Right", "unit", "V/L", h5err)
      call h5ltmake_dataset_f(file_id, "mean_dE_cov_Left",3,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*left_mean_dE_cov3,h5err)
      call h5ltset_attribute_string_f(file_id, "mean_dE_cov_Left", "unit", "V/L", h5err)

      call h5ltmake_dataset_f(file_id, "FLA/dE_cov",3,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*intp_FLA_dE_cov3,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/dE_cov", "unit", "V/L", h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/dE_cov", "description", "This is -grad(phi) in PEST coordinate. The length unit is parametric length in PEST coordinate", h5err)
      call h5ltmake_dataset_f(file_id, "FLA/dE_cov_Right",3,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*right_FLA_dE_cov3,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/dE_cov_Right", "unit", "V/L", h5err)
      call h5ltmake_dataset_f(file_id, "FLA/dE_cov_Left",3,h5dims,H5T_NATIVE_DOUBLE, (E_00*c)/(charg0*1e8)*left_FLA_dE_cov3,h5err)
      call h5ltset_attribute_string_f(file_id, "FLA/dE_cov_Left", "unit", "V/L", h5err)


      deallocate(intp_val_Arr)


      ! Print SG (Diagnostic Grid) Data
      ! Note that dimension has to use SG%dim_* rather than dim_*
      ! because the dimension of snapshot_hdf5 and SG might not be the same
      h5dims(1)=1
      call h5ltmake_dataset_f(gid_SG, "ds",1,h5dims,H5T_NATIVE_DOUBLE, SG%ds,h5err)
      call h5ltmake_dataset_f(gid_SG, "dgq",1,h5dims,H5T_NATIVE_DOUBLE, SG%dgq,h5err)

      h5dims(1)=SG%dim_s
      call h5ltmake_dataset_f(gid_SG, "a_1D",1,h5dims,H5T_NATIVE_DOUBLE, SG%a_1D,h5err)
      call h5ltmake_dataset_f(gid_SG, "s_1D",1,h5dims,H5T_NATIVE_DOUBLE, SG%s_1D,h5err)
      call h5ltmake_dataset_f(gid_SG, "safety_q_1D",1,h5dims,H5T_NATIVE_DOUBLE, SG%safety_q_1D,h5err)

      h5dims(1)=SG%dim_gq
      call h5ltmake_dataset_f(gid_SG, "gq_1D",1,h5dims,H5T_NATIVE_DOUBLE, SG%gq_1D,h5err)

      h5dims(1:2)=(/SG%dim_gq,SG%dim_s/)
      call h5ltmake_dataset_f(gid_SG, "a_2D",2,h5dims,H5T_NATIVE_DOUBLE, SG%a_2D,h5err)
      call h5ltmake_dataset_f(gid_SG, "gq_2D",2,h5dims,H5T_NATIVE_DOUBLE, SG%gq_2D,h5err)
      call h5ltmake_dataset_f(gid_SG, "Cell_dV",2,h5dims,H5T_NATIVE_DOUBLE, SG%Cell_dV,h5err)
      call h5ltmake_dataset_f(gid_SG, "Cell_dS_s",2,h5dims,H5T_NATIVE_DOUBLE, SG%Cell_dS_s,h5err)
      call h5ltmake_dataset_f(gid_SG, "Cell_dS_gq",2,h5dims,H5T_NATIVE_DOUBLE, SG%Cell_dS_gq,h5err)
      call h5ltmake_dataset_f(gid_SG, "Cell_dS_gf",2,h5dims,H5T_NATIVE_DOUBLE, SG%Cell_dS_gf,h5err)

      call h5ltmake_dataset_f(gid_SG, "dni",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*SG%dni(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_SG, "dne",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*SG%dne(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_SG, "ni",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*SG%ni(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_SG, "ne",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*SG%ne(:,:,1),h5err)

      call Smooth_SG_Data_2D(3,SG%dni(:,:,1),SG%Cell_dV)
      call Smooth_SG_Data_2D(3,SG%dne(:,:,1),SG%Cell_dV)
      call Smooth_SG_Data_2D(3,SG%ni(:,:,1),SG%Cell_dV)
      call Smooth_SG_Data_2D(3,SG%ne(:,:,1),SG%Cell_dV)
      call h5ltmake_dataset_f(gid_SG_SM, "dni",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*SG%dni(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_SG_SM, "dne",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*SG%dne(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_SG_SM, "ni",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*SG%ni(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_SG_SM, "ne",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*SG%ne(:,:,1),h5err)

      call h5ltmake_dataset_f(gid_SG, "phi",2,h5dims,H5T_NATIVE_DOUBLE, SG%phi(:,:,1),h5err)
      h5dims(1:3)=(/3,SG%dim_gq,SG%dim_s/)
      call h5ltmake_dataset_f(gid_SG, "grad_phi_cov",3,h5dims,H5T_NATIVE_DOUBLE, SG%grad_phi_cov(:,:,:,1),h5err)


      h5dims(1:2)=(/SG%dim_gq,SG%dim_s/)
      ! Followings are for the debugging purpose
      allocate(intp_val_Arr(SG%dim_gq*SG%dim_s))
      call Interpolate_SG_into_GTSmgrid(1e6*SG%dne(:,:,1),tmp_vec_mgrid_size)
      call Interpolate_OriGridData_into_SG(tmp_vec_mgrid_size,intp_val_Arr)
      call h5ltmake_dataset_f(gid_SG, "intp_SM_dne",2,h5dims,H5T_NATIVE_DOUBLE, intp_val_Arr,h5err)

      call Interpolate_SG_into_GTSmgrid((1.0/1.602176634e-19)*SG%test_src_2D,tmp_vec_mgrid_size)
      call Interpolate_OriGridData_into_SG(tmp_vec_mgrid_size,intp_val_Arr)
      call h5ltmake_dataset_f(gid_SG, "intp_test_src_2D",2,h5dims,H5T_NATIVE_DOUBLE, intp_val_Arr,h5err)
      deallocate(intp_val_Arr)

      call h5ltmake_dataset_f(gid_SG, "test_phi_2D",2,h5dims,H5T_NATIVE_DOUBLE, SG%test_phi_2D,h5err)
      call h5ltmake_dataset_f(gid_SG, "test_src_2D",2,h5dims,H5T_NATIVE_DOUBLE, (1.0/1.602176634e-19)*SG%test_src_2D,h5err)
      call h5ltmake_dataset_f(gid_SG, "test_res_2D",2,h5dims,H5T_NATIVE_DOUBLE, (1.0/1.602176634e-19)*SG%test_res_2D,h5err)

      !Following 3D data are for temporary diagnostic purpose
      h5dims(1:3)=(/SG%dim_gq,SG%dim_s,mzetamax/)
      call h5ltmake_dataset_f(gid_SG, "test_phi_3D",3,h5dims,H5T_NATIVE_DOUBLE, SG%test_phi_3D,h5err)
      call h5ltmake_dataset_f(gid_SG, "test_src_3D",3,h5dims,H5T_NATIVE_DOUBLE, (1.0/1.602176634e-19)*SG%test_src_3D,h5err)
      call h5ltmake_dataset_f(gid_SG, "test_res_3D",3,h5dims,H5T_NATIVE_DOUBLE, (1.0/1.602176634e-19)*SG%test_res_3D,h5err)

      call Smooth_SG_Data_2D(3,SG%test_phi_2D,SG%Cell_dV)
      call Smooth_SG_Data_2D(3,SG%test_src_2D,SG%Cell_dV)
      call h5ltmake_dataset_f(gid_SG_SM, "test_phi_2D",2,h5dims,H5T_NATIVE_DOUBLE, SG%test_phi_2D,h5err)
      call h5ltmake_dataset_f(gid_SG_SM, "test_src_2D",2,h5dims,H5T_NATIVE_DOUBLE, (1.0/1.602176634e-19)*SG%test_src_2D,h5err)


      h5dims(1:2)=3
      h5dims(3)=SG%dim_gq+1
      h5dims(4)=SG%dim_a+2
      h5dims(5)=3
      call h5ltmake_dataset_f(gid_SG, "test_CT",5,h5dims,H5T_NATIVE_DOUBLE, CT,h5err)
      call h5ltmake_dataset_f(gid_SG, "test_Metric_T",4,h5dims,H5T_NATIVE_DOUBLE, Metric_T,h5err)
      call h5ltmake_dataset_f(gid_SG, "test_Polar_T",4,h5dims,H5T_NATIVE_DOUBLE, Polar_T,h5err)


      ! Print FP (FootPrint) Data
      h5dims(1)=FP_BDY1%dim_gq
      h5dims(2)=FP_BDY1%dim_gf

      call h5ltmake_dataset_f(gid_FP, "R_2D",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%R_2D,h5err)
      call h5ltmake_dataset_f(gid_FP, "Z_2D",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%Z_2D,h5err)
      call h5ltmake_dataset_f(gid_FP, "a_2D",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%a_2D,h5err)
      call h5ltmake_dataset_f(gid_FP, "gq_2D",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%gq_2D,h5err)
      call h5ltmake_dataset_f(gid_FP, "gf_2D",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%gf_2D,h5err)
      call h5ltmake_dataset_f(gid_FP, "BB_2D",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%BB_2D,h5err)
      call h5ltmake_dataset_f(gid_FP, "eFlux_dP",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%eScaData(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_FP, "eFlux_dT",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%eScaData(:,:,2),h5err)
      call h5ltmake_dataset_f(gid_FP, "eFlux_dE",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%eScaData(:,:,3),h5err)
      call h5ltmake_dataset_f(gid_FP, "eFlux_Marker",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%eScaData(:,:,4),h5err)
      call h5ltmake_dataset_f(gid_FP, "eFlux_P",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%eScaData(:,:,5),h5err)
      call h5ltmake_dataset_f(gid_FP, "eFlux_E",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%eScaData(:,:,6),h5err)
      call h5ltmake_dataset_f(gid_FP, "iFlux_dP",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%iScaData(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_FP, "iFlux_dT",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%iScaData(:,:,2),h5err)
      call h5ltmake_dataset_f(gid_FP, "iFlux_dE",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%iScaData(:,:,3),h5err)
      call h5ltmake_dataset_f(gid_FP, "iFlux_Marker",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%iScaData(:,:,4),h5err)
      call h5ltmake_dataset_f(gid_FP, "iFlux_P",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%iScaData(:,:,5),h5err)
      call h5ltmake_dataset_f(gid_FP, "iFlux_E",2,h5dims,H5T_NATIVE_DOUBLE, FP_BDY1%iScaData(:,:,6),h5err)



      call h5ltmake_dataset_f(gid_FP_cum, "R_2D",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%R_2D,h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "Z_2D",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%Z_2D,h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "a_2D",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%a_2D,h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "gq_2D",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%gq_2D,h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "gf_2D",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%gf_2D,h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "BB_2D",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%BB_2D,h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "eFlux_dP",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%eScaData(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "eFlux_dT",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%eScaData(:,:,2),h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "eFlux_dE",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%eScaData(:,:,3),h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "eFlux_Marker",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%eScaData(:,:,4),h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "eFlux_P",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%eScaData(:,:,5),h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "eFlux_E",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%eScaData(:,:,6),h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "iFlux_dP",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%iScaData(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "iFlux_dT",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%iScaData(:,:,2),h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "iFlux_dE",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%iScaData(:,:,3),h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "iFlux_Marker",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%iScaData(:,:,4),h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "iFlux_P",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%iScaData(:,:,5),h5err)
      call h5ltmake_dataset_f(gid_FP_cum, "iFlux_E",2,h5dims,H5T_NATIVE_DOUBLE, cum_FP_BDY1%iScaData(:,:,6),h5err)


      ! Print DG (Diagnostic Grid) Data
      ! Note that dimension has to use DG%dim_* rather than dim_*
      ! because the dimension of snapshot_hdf5 and DG might not be the same
      h5dims(1)=DG%dim_a
      call h5ltmake_dataset_f(gid_DG, "a_1D",1,h5dims,H5T_NATIVE_DOUBLE, DG%a_1D,h5err)
      call h5ltmake_dataset_f(gid_DG, "s_1D",1,h5dims,H5T_NATIVE_DOUBLE, DG%s_1D,h5err)
      call h5ltmake_dataset_f(gid_DG, "Anulus_dV",1,h5dims,H5T_NATIVE_DOUBLE, DG%Anulus_dV,h5err)

      call h5ltmake_dataset_f(gid_DG, "eFlux1D_P",1,h5dims,H5T_NATIVE_DOUBLE, DG%eFluxData_a1D(:,1)*uveloci*1e4,h5err)
      call h5ltmake_dataset_f(gid_DG, "eFlux1D_E",1,h5dims,H5T_NATIVE_DOUBLE, DG%eFluxData_a1D(:,2)*uveloci*E_00*1e-3,h5err)
      call h5ltmake_dataset_f(gid_DG, "eFlux1D_P_ExB",1,h5dims,H5T_NATIVE_DOUBLE, DG%eFluxData_a1D(:,3)*uveloci*1e4,h5err)
      call h5ltmake_dataset_f(gid_DG, "eFlux1D_E_ExB",1,h5dims,H5T_NATIVE_DOUBLE, DG%eFluxData_a1D(:,4)*uveloci*E_00*1e-3,h5err)
      call h5ltmake_dataset_f(gid_DG, "eFlux1D_P_streaming",1,h5dims,H5T_NATIVE_DOUBLE, DG%eFluxData_a1D(:,5)*uveloci*1e4,h5err)
      call h5ltmake_dataset_f(gid_DG, "eFlux1D_E_streaming",1,h5dims,H5T_NATIVE_DOUBLE, DG%eFluxData_a1D(:,6)*uveloci*E_00*1e-3,h5err)
      call h5ltmake_dataset_f(gid_DG, "mean_eFlux1D_P",1,h5dims,H5T_NATIVE_DOUBLE, inv_msnap*DG%cum_eFluxData_a1D(:,1)*uveloci*1e4,h5err)
      call h5ltmake_dataset_f(gid_DG, "mean_eFlux1D_E",1,h5dims,H5T_NATIVE_DOUBLE, inv_msnap*DG%cum_eFluxData_a1D(:,2)*uveloci*E_00*1e-3,h5err)
      call h5ltmake_dataset_f(gid_DG, "mean_eFlux1D_P_ExB",1,h5dims,H5T_NATIVE_DOUBLE, inv_msnap*DG%cum_eFluxData_a1D(:,3)*uveloci*1e4,h5err)
      call h5ltmake_dataset_f(gid_DG, "mean_eFlux1D_E_ExB",1,h5dims,H5T_NATIVE_DOUBLE, inv_msnap*DG%cum_eFluxData_a1D(:,4)*uveloci*E_00*1e-3,h5err)
      call h5ltmake_dataset_f(gid_DG, "mean_eFlux1D_P_streaming",1,h5dims,H5T_NATIVE_DOUBLE, inv_msnap*DG%cum_eFluxData_a1D(:,5)*uveloci*1e4,h5err)
      call h5ltmake_dataset_f(gid_DG, "mean_eFlux1D_E_streaming",1,h5dims,H5T_NATIVE_DOUBLE, inv_msnap*DG%cum_eFluxData_a1D(:,6)*uveloci*E_00*1e-3,h5err)

      call h5ltmake_dataset_f(gid_DG, "iFlux1D_P",1,h5dims,H5T_NATIVE_DOUBLE, DG%iFluxData_a1D(:,1)*uveloci*1e4,h5err)
      call h5ltmake_dataset_f(gid_DG, "iFlux1D_E",1,h5dims,H5T_NATIVE_DOUBLE, DG%iFluxData_a1D(:,2)*uveloci*E_00*1e-3,h5err)
      call h5ltmake_dataset_f(gid_DG, "iFlux1D_P_ExB",1,h5dims,H5T_NATIVE_DOUBLE, DG%iFluxData_a1D(:,3)*uveloci*1e4,h5err)
      call h5ltmake_dataset_f(gid_DG, "iFlux1D_E_ExB",1,h5dims,H5T_NATIVE_DOUBLE, DG%iFluxData_a1D(:,4)*uveloci*E_00*1e-3,h5err)
      call h5ltmake_dataset_f(gid_DG, "iFlux1D_P_streaming",1,h5dims,H5T_NATIVE_DOUBLE, DG%iFluxData_a1D(:,5)*uveloci*1e4,h5err)
      call h5ltmake_dataset_f(gid_DG, "iFlux1D_E_streaming",1,h5dims,H5T_NATIVE_DOUBLE, DG%iFluxData_a1D(:,6)*uveloci*E_00*1e-3,h5err)
      call h5ltmake_dataset_f(gid_DG, "mean_iFlux1D_P",1,h5dims,H5T_NATIVE_DOUBLE, inv_msnap*DG%cum_iFluxData_a1D(:,1)*uveloci*1e4,h5err)
      call h5ltmake_dataset_f(gid_DG, "mean_iFlux1D_E",1,h5dims,H5T_NATIVE_DOUBLE, inv_msnap*DG%cum_iFluxData_a1D(:,2)*uveloci*E_00*1e-3,h5err)
      call h5ltmake_dataset_f(gid_DG, "mean_iFlux1D_P_ExB",1,h5dims,H5T_NATIVE_DOUBLE, inv_msnap*DG%cum_iFluxData_a1D(:,3)*uveloci*1e4,h5err)
      call h5ltmake_dataset_f(gid_DG, "mean_iFlux1D_E_ExB",1,h5dims,H5T_NATIVE_DOUBLE, inv_msnap*DG%cum_iFluxData_a1D(:,4)*uveloci*E_00*1e-3,h5err)
      call h5ltmake_dataset_f(gid_DG, "mean_iFlux1D_P_streaming",1,h5dims,H5T_NATIVE_DOUBLE, inv_msnap*DG%cum_iFluxData_a1D(:,5)*uveloci*1e4,h5err)
      call h5ltmake_dataset_f(gid_DG, "mean_iFlux1D_E_streaming",1,h5dims,H5T_NATIVE_DOUBLE, inv_msnap*DG%cum_iFluxData_a1D(:,6)*uveloci*E_00*1e-3,h5err)


      h5dims(1)=DG%dim_gq
      call h5ltmake_dataset_f(gid_DG, "gq_1D",1,h5dims,H5T_NATIVE_DOUBLE, DG%gq_1D,h5err)

      h5dims(1:2)=(/DG%dim_gq,DG%dim_a/)
      call h5ltmake_dataset_f(gid_DG, "a_2D",2,h5dims,H5T_NATIVE_DOUBLE, DG%a_2D,h5err)
      call h5ltmake_dataset_f(gid_DG, "gq_2D",2,h5dims,H5T_NATIVE_DOUBLE, DG%gq_2D,h5err)
      call h5ltmake_dataset_f(gid_DG, "Cell_dV",2,h5dims,H5T_NATIVE_DOUBLE, DG%Cell_dV,h5err)
      call h5ltmake_dataset_f(gid_DG, "Cell_dS_a",2,h5dims,H5T_NATIVE_DOUBLE, DG%Cell_dS_a,h5err)
      call h5ltmake_dataset_f(gid_DG, "Cell_dS_gq",2,h5dims,H5T_NATIVE_DOUBLE, DG%Cell_dS_gq,h5err)
      call h5ltmake_dataset_f(gid_DG, "Cell_dS_gf",2,h5dims,H5T_NATIVE_DOUBLE, DG%Cell_dS_gf,h5err)
      call h5ltmake_dataset_f(gid_DG, "FLabel_Function",2,h5dims,H5T_NATIVE_DOUBLE, DG%FLabel_Function(:,:,1),h5err)

      call h5ltmake_dataset_f(gid_DG_FLA, "Denominator",2,h5dims,H5T_NATIVE_DOUBLE, DG%FLA_Denominator(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "FLabel",2,h5dims,H5T_NATIVE_DOUBLE, DG%FLA_FLabel(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "psi0",2,h5dims,H5T_NATIVE_DOUBLE, DG%FLA_psi0(:,:,1),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "psi0_Left",2,h5dims,H5T_NATIVE_DOUBLE, DG%FLA_psi0(:,:,0),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "SM/Denominator",2,h5dims,H5T_NATIVE_DOUBLE, DG%FLA_SM_Denominator,h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "SM/FLabel",2,h5dims,H5T_NATIVE_DOUBLE, DG%FLA_SM_FLabel,h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "SM/FLabel_Left",2,h5dims,H5T_NATIVE_DOUBLE, DG%FLA_SM_FLabel_Left,h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "SM/FLabel_Right",2,h5dims,H5T_NATIVE_DOUBLE, DG%FLA_SM_FLabel_Right,h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "SM/psi0",2,h5dims,H5T_NATIVE_DOUBLE, DG%FLA_SM_psi0,h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "SM/mag_grad_FL",2,h5dims,H5T_NATIVE_DOUBLE, DG%mag_grad_FLA_SM_FLabel(:,:,1),h5err)


      call h5ltmake_dataset_f(gid_DG_FLA, "eFlux_p_B",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*DG_FLA_eScaData(:,:,7),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "eFlux_p_N",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*DG_FLA_eScaData(:,:,8),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "eFlux_p_BxN",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*DG_FLA_eScaData(:,:,9),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "eFlux_p_B_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*DG_FLA_eScaData(:,:,10),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "eFlux_p_N_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*DG_FLA_eScaData(:,:,11),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "eFlux_p_BxN_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*DG_FLA_eScaData(:,:,12),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "eFlux_E_B",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*DG_FLA_eScaData(:,:,13),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "eFlux_E_N",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*DG_FLA_eScaData(:,:,14),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "eFlux_E_BxN",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*DG_FLA_eScaData(:,:,15),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "eFlux_E_B_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*DG_FLA_eScaData(:,:,16),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "eFlux_E_N_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*DG_FLA_eScaData(:,:,17),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "eFlux_E_BxN_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*DG_FLA_eScaData(:,:,18),h5err)

      call h5ltmake_dataset_f(gid_DG_FLA, "iFlux_p_B",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*DG_FLA_iScaData(:,:,9),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "iFlux_p_N",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*DG_FLA_iScaData(:,:,10),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "iFlux_p_BxN",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*DG_FLA_iScaData(:,:,11),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "iFlux_p_B_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*DG_FLA_iScaData(:,:,12),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "iFlux_p_N_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*DG_FLA_iScaData(:,:,13),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "iFlux_p_BxN_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*DG_FLA_iScaData(:,:,14),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "iFlux_E_B",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*DG_FLA_iScaData(:,:,15),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "iFlux_E_N",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*DG_FLA_iScaData(:,:,16),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "iFlux_E_BxN",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*DG_FLA_iScaData(:,:,17),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "iFlux_E_B_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*DG_FLA_iScaData(:,:,18),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "iFlux_E_N_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*DG_FLA_iScaData(:,:,19),h5err)
      call h5ltmake_dataset_f(gid_DG_FLA, "iFlux_E_BxN_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*DG_FLA_iScaData(:,:,20),h5err)



      h5dims(1:3)=(/3,DG%dim_gq,DG%dim_a/)
      call h5ltmake_dataset_f(gid_DG_FLA, "SM/grad_FL",3,h5dims,H5T_NATIVE_DOUBLE, DG%grad_FLA_SM_FLabel(:,:,:,1),h5err)
      h5dims(1:2)=(/DG%dim_gq,DG%dim_a/)

      do i=1,2
         select case(i)
            case(1)
               ! Write GD variables
               gid=gid_DG
            case(2)
               ! Smooth GD variables and write them into /GD/SM group
               call Smooth_DG_Scalar_Vector_Data(2) ! smooth 3 times
               gid=gid_DG_SM
         end select

         h5dims(1:2)=(/DG%dim_gq,DG%dim_a/)
         call h5ltmake_dataset_f(gid, "dne_tot",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*(DG%eScaData(:,:,1,1)+int(ifelectron/2)*DG%dne_adiab),h5err)
         call h5ltset_attribute_string_f(gid, "dne_tot", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid, "dne_adiab",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*(DG%dne_adiab),h5err)
         call h5ltset_attribute_string_f(gid, "dne_adiab", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid, "dTe",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%eScaData(:,:,2,1)),h5err)
         call h5ltset_attribute_string_f(gid, "dTe", "unit", "eV", h5err)
         call h5ltmake_dataset_f(gid, "dTe_consst",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%eScaData(:,:,3,1)),h5err)
         call h5ltset_attribute_string_f(gid, "dTe_consst", "unit", "eV", h5err)
         call h5ltmake_dataset_f(gid, "dTe_consst_para",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%eScaData(:,:,20,1)),h5err)
         call h5ltmake_dataset_f(gid, "dTe_consst_perp",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%eScaData(:,:,21,1)),h5err)
         call h5ltmake_dataset_f(gid, "ne_marker",2,h5dims,H5T_NATIVE_DOUBLE, DG%eScaData(:,:,4,1),h5err)
         call h5ltset_attribute_string_f(gid, "ne_marker", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid, "ne_marker_from_walls",2,h5dims,H5T_NATIVE_DOUBLE, DG%eScaData(:,:,19,1),h5err)
         call h5ltset_attribute_string_f(gid, "ne_marker_from_walls", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid, "ne_tot",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*(DG%eScaData(:,:,5,1)),h5err)
         call h5ltset_attribute_string_f(gid, "ne_tot", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid, "Te",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%eScaData(:,:,6,1)),h5err)
         call h5ltset_attribute_string_f(gid, "Te", "unit", "eV", h5err)
         call h5ltmake_dataset_f(gid, "Te_para",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%eScaData(:,:,22,1)),h5err)
         call h5ltmake_dataset_f(gid, "Te_perp",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%eScaData(:,:,23,1)),h5err)
         call h5ltmake_dataset_f(gid, "w_dne",2,h5dims,H5T_NATIVE_DOUBLE, DG%w_dne,h5err)
         call h5ltmake_dataset_f(gid, "ne_tot_hybrid",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*(DG%ne_tot_hybrid),h5err)
         call h5ltset_attribute_string_f(gid, "ne_tot_hybrid", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid, "Te_hybrid",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%Te_hybrid),h5err)
         call h5ltset_attribute_string_f(gid, "Te_hybrid", "unit", "eV", h5err)
         call h5ltmake_dataset_f(gid, "eFlux_p_B",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%eScaData(:,:,7,1)),h5err)
         call h5ltmake_dataset_f(gid, "eFlux_p_N",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%eScaData(:,:,8,1)),h5err)
         call h5ltmake_dataset_f(gid, "eFlux_p_BxN",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%eScaData(:,:,9,1)),h5err)
         call h5ltmake_dataset_f(gid, "eFlux_p_B_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%eScaData(:,:,10,1)),h5err)
         call h5ltmake_dataset_f(gid, "eFlux_p_N_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%eScaData(:,:,11,1)),h5err)
         call h5ltmake_dataset_f(gid, "eFlux_p_BxN_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%eScaData(:,:,12,1)),h5err)
         call h5ltmake_dataset_f(gid, "eFlux_E_B",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%eScaData(:,:,13,1)),h5err)
         call h5ltmake_dataset_f(gid, "eFlux_E_N",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%eScaData(:,:,14,1)),h5err)
         call h5ltmake_dataset_f(gid, "eFlux_E_BxN",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%eScaData(:,:,15,1)),h5err)
         call h5ltmake_dataset_f(gid, "eFlux_E_B_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%eScaData(:,:,16,1)),h5err)
         call h5ltmake_dataset_f(gid, "eFlux_E_N_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%eScaData(:,:,17,1)),h5err)
         call h5ltmake_dataset_f(gid, "eFlux_E_BxN_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%eScaData(:,:,18,1)),h5err)


         call h5ltmake_dataset_f(gid, "dni",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*(DG%iScaData(:,:,1,1)),h5err)
         call h5ltset_attribute_string_f(gid, "dni", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid, "dTi",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%iScaData(:,:,2,1)),h5err)
         call h5ltset_attribute_string_f(gid, "dTi", "unit", "eV", h5err)
         call h5ltmake_dataset_f(gid, "dTi_consst",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%iScaData(:,:,3,1)),h5err)
         call h5ltset_attribute_string_f(gid, "dTi_consst", "unit", "eV", h5err)
         call h5ltmake_dataset_f(gid, "dTi_consst_para",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%iScaData(:,:,22,1)),h5err)
         call h5ltmake_dataset_f(gid, "dTi_consst_perp",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%iScaData(:,:,23,1)),h5err)
         call h5ltmake_dataset_f(gid, "ni_marker",2,h5dims,H5T_NATIVE_DOUBLE, DG%iScaData(:,:,4,1),h5err)
         call h5ltset_attribute_string_f(gid, "ni_marker", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid, "ni_marker_from_walls",2,h5dims,H5T_NATIVE_DOUBLE, DG%iScaData(:,:,21,1),h5err)
         call h5ltset_attribute_string_f(gid, "ni_marker_from_walls", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid, "ni",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*(DG%iScaData(:,:,5,1)),h5err)
         call h5ltset_attribute_string_f(gid, "ni", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid, "Ti",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%iScaData(:,:,6,1)),h5err)
         call h5ltmake_dataset_f(gid, "Ti_para",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%iScaData(:,:,24,1)),h5err)
         call h5ltmake_dataset_f(gid, "Ti_perp",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%iScaData(:,:,25,1)),h5err)
         call h5ltmake_dataset_f(gid, "w_dni",2,h5dims,H5T_NATIVE_DOUBLE, DG%w_dni,h5err)
         call h5ltmake_dataset_f(gid, "ni_hybrid",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*(DG%ni_hybrid),h5err)
         call h5ltset_attribute_string_f(gid, "ni_hybrid", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid, "Ti_hybrid",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%Ti_hybrid),h5err)
         call h5ltset_attribute_string_f(gid, "Ti_hybrid", "unit", "eV", h5err)
         call h5ltset_attribute_string_f(gid, "Ti", "unit", "eV", h5err)
         call h5ltmake_dataset_f(gid, "dni_Larmor",2,h5dims,H5T_NATIVE_DOUBLE, 1e6*(DG%iScaData(:,:,7,1)),h5err)
         call h5ltset_attribute_string_f(gid, "dni_Larmor", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid, "dTi_Larmor_consst",2,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*(DG%iScaData(:,:,8,1)),h5err)
         call h5ltset_attribute_string_f(gid, "dTi_Larmor_consst", "unit", "eV", h5err)
         call h5ltmake_dataset_f(gid, "iFlux_p_B",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%iScaData(:,:,9,1)),h5err)
         call h5ltmake_dataset_f(gid, "iFlux_p_N",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%iScaData(:,:,10,1)),h5err)
         call h5ltmake_dataset_f(gid, "iFlux_p_BxN",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%iScaData(:,:,11,1)),h5err)
         call h5ltmake_dataset_f(gid, "iFlux_p_B_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%iScaData(:,:,12,1)),h5err)
         call h5ltmake_dataset_f(gid, "iFlux_p_N_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%iScaData(:,:,13,1)),h5err)
         call h5ltmake_dataset_f(gid, "iFlux_p_BxN_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%iScaData(:,:,14,1)),h5err)
         call h5ltmake_dataset_f(gid, "iFlux_E_B",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%iScaData(:,:,15,1)),h5err)
         call h5ltmake_dataset_f(gid, "iFlux_E_N",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%iScaData(:,:,16,1)),h5err)
         call h5ltmake_dataset_f(gid, "iFlux_E_BxN",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%iScaData(:,:,17,1)),h5err)
         call h5ltmake_dataset_f(gid, "iFlux_E_B_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%iScaData(:,:,18,1)),h5err)
         call h5ltmake_dataset_f(gid, "iFlux_E_N_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%iScaData(:,:,19,1)),h5err)
         call h5ltmake_dataset_f(gid, "iFlux_E_BxN_outwards",2,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%iScaData(:,:,20,1)),h5err)


         h5dims(1:3)=(/3,DG%dim_gq,DG%dim_a/)
         call h5ltmake_dataset_f(gid, "eFlux_p",3,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%eVecData(:,:,:,1,1)),h5err)
         call h5ltset_attribute_string_f(gid, "eFlux_p", "unit", "1/(m^2*s)", h5err)
         call h5ltmake_dataset_f(gid, "eFlux_p_ExB",3,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%eVecData(:,:,:,2,1)),h5err)
         call h5ltset_attribute_string_f(gid, "eFlux_p_ExB", "unit", "1/(m^2*s)", h5err)
         call h5ltmake_dataset_f(gid, "eFlux_E",3,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%eVecData(:,:,:,3,1)),h5err)
         call h5ltset_attribute_string_f(gid, "eFlux_E", "unit", "J/(m^2*s)", h5err)
         call h5ltmake_dataset_f(gid, "eFlux_E_ExB",3,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%eVecData(:,:,:,4,1)),h5err)
         call h5ltset_attribute_string_f(gid, "eFlux_E_ExB", "unit", "J/(m^2*s)", h5err)
         call h5ltmake_dataset_f(gid, "eAvgCFL",3,h5dims,H5T_NATIVE_DOUBLE, (DG%eVecData(:,:,:,5,1)),h5err)

         call h5ltmake_dataset_f(gid, "iFlux_p",3,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%iVecData(:,:,:,1,1)),h5err)
         call h5ltset_attribute_string_f(gid, "iFlux_p", "unit", "1/(m^2*s)", h5err)
         call h5ltmake_dataset_f(gid, "iFlux_p_ExB",3,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%iVecData(:,:,:,2,1)),h5err)
         call h5ltset_attribute_string_f(gid, "iFlux_p_ExB", "unit", "1/(m^2*s)", h5err)
         call h5ltmake_dataset_f(gid, "iFlux_E",3,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%iVecData(:,:,:,3,1)),h5err)
         call h5ltset_attribute_string_f(gid, "iFlux_E", "unit", "J/(m^2*s)", h5err)
         call h5ltmake_dataset_f(gid, "iFlux_E_ExB",3,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%iVecData(:,:,:,4,1)),h5err)
         call h5ltset_attribute_string_f(gid, "iFlux_E_ExB", "unit", "J/(m^2*s)", h5err)
         call h5ltmake_dataset_f(gid, "iAvgCFL",3,h5dims,H5T_NATIVE_DOUBLE, (DG%iVecData(:,:,:,5,1)),h5err)

         call h5ltmake_dataset_f(gid, "iFlux_p_Larmor",3,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%iVecData(:,:,:,6,1)),h5err)
         call h5ltset_attribute_string_f(gid, "iFlux_p_Larmor", "unit", "1/(m^2*s)", h5err)
         call h5ltmake_dataset_f(gid, "iFlux_p_ExB_Larmor",3,h5dims,H5T_NATIVE_DOUBLE, uveloci*1e4*(DG%iVecData(:,:,:,7,1)),h5err)
         call h5ltset_attribute_string_f(gid, "iFlux_p_ExB_Larmor", "unit", "1/(m^2*s)", h5err)
         call h5ltmake_dataset_f(gid, "iFlux_E_Larmor",3,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%iVecData(:,:,:,8,1)),h5err)
         call h5ltset_attribute_string_f(gid, "iFlux_E_Larmor", "unit", "J/(m^2*s)", h5err)
         call h5ltmake_dataset_f(gid, "iFlux_E_ExB_Larmor",3,h5dims,H5T_NATIVE_DOUBLE, uveloci*E_00*1e-3*(DG%iVecData(:,:,:,9,1)),h5err)
      enddo


      ! Write PVS information
      h5dims(1) = 1
      call h5ltmake_dataset_f(file_id, "nPVS",1,h5dims, H5T_NATIVE_INTEGER, nPVS,h5err)
      do pbid=0,nPVS

         write(cdum,'("PVS",i2.2)') pbid
         call h5gcreate_f(file_id, cdum, gid_PVS, h5err) ! PVS(ProbeVelocitySpace)


         h5dims(1)=1
         call h5ltmake_dataset_f(gid_PVS, "a",1,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%a,h5err)
         call h5ltmake_dataset_f(gid_PVS, "gq",1,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%gq,h5err)
         call h5ltmake_dataset_f(gid_PVS, "gf",1,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%gf,h5err)
         call h5ltmake_dataset_f(gid_PVS, "L_a",1,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%L_a,h5err)
         call h5ltmake_dataset_f(gid_PVS, "L_gq",1,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%L_gq,h5err)
         call h5ltmake_dataset_f(gid_PVS, "L_gf",1,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%L_gf,h5err)

         call h5ltmake_dataset_f(gid_PVS, "n0",1,h5dims,H5T_NATIVE_DOUBLE, 1e6*PVS(pbid)%n0,h5err)
         call h5ltset_attribute_string_f(gid_PVS, "n0", "unit", "#/m^3", h5err)
         call h5ltmake_dataset_f(gid_PVS, "Te0",1,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*PVS(pbid)%Te0,h5err)
         call h5ltset_attribute_string_f(gid_PVS, "Te0", "unit", "eV", h5err)
         call h5ltmake_dataset_f(gid_PVS, "Ti0",1,h5dims,H5T_NATIVE_DOUBLE, 1e3*Tem_00*PVS(pbid)%Ti0,h5err)
         call h5ltset_attribute_string_f(gid_PVS, "Ti0", "unit", "eV", h5err)

         h5dims(1)=PVS(pbid)%dim_vpara
         call h5ltmake_dataset_f(gid_PVS, "vpara_1D",1,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%vpara_1D,h5err)
         call h5ltset_attribute_string_f(gid_PVS, "vpara_1D", "unit", "v/vth", h5err)
         h5dims(1)=PVS(pbid)%dim_vperp
         call h5ltmake_dataset_f(gid_PVS, "vperp_1D",1,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%vperp_1D,h5err)
         call h5ltset_attribute_string_f(gid_PVS, "vperp_1D", "unit", "v/vth", h5err)

         h5dims(1:2)=(/PVS(pbid)%dim_vpara,PVS(pbid)%dim_vperp/)
         call h5ltmake_dataset_f(gid_PVS, "dfi",2,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%fi(:,:,1),h5err)
         call h5ltmake_dataset_f(gid_PVS, "dfi_dTi",2,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%fi(:,:,2),h5err)
         call h5ltmake_dataset_f(gid_PVS, "fi_marker",2,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%fi(:,:,3),h5err)
         call h5ltmake_dataset_f(gid_PVS, "fi",2,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%fi(:,:,4),h5err)
         call h5ltmake_dataset_f(gid_PVS, "fi_Ti",2,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%fi(:,:,5),h5err)

         call h5ltmake_dataset_f(gid_PVS, "dfe",2,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%fe(:,:,1),h5err)
         call h5ltmake_dataset_f(gid_PVS, "dfe_dTe",2,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%fe(:,:,2),h5err)
         call h5ltmake_dataset_f(gid_PVS, "fe_marker",2,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%fe(:,:,3),h5err)
         call h5ltmake_dataset_f(gid_PVS, "fe",2,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%fe(:,:,4),h5err)
         call h5ltmake_dataset_f(gid_PVS, "fe_Te",2,h5dims,H5T_NATIVE_DOUBLE, PVS(pbid)%fe(:,:,5),h5err)

         call h5gclose_f(gid_PVS, h5err)
      enddo

      call h5gclose_f(gid_PARAMS, h5err)
      call h5gclose_f(gid_FLA, h5err)
      call h5gclose_f(gid_DCOMP, h5err)
      call h5gclose_f(gid_BIN, h5err)
      call h5gclose_f(gid_SG, h5err)
      call h5gclose_f(gid_SG_SM, h5err)
      call h5gclose_f(gid_DG, h5err)
      call h5gclose_f(gid_DG_SM, h5err)
      call h5gclose_f(gid_DG_FLA, h5err)
      call h5gclose_f(gid_DG_FLA_SM, h5err)
      call h5gclose_f(gid_FP, h5err)
      call h5gclose_f(gid_FP_cum, h5err)

      call h5fclose_f(file_id, h5err)
      call h5close_f(h5err)

      !close(snapout)
   endif ! if(mype == numberpe-1)then !mype=0 must also have myrank_toroidal=0


   ! Copy snapshot data of last poloidal plane into AllPlanes if istep==msnap_for_All
   if(mod(istep,msnap_for_All)==0 .and. mype == numberpe-1) then
      write(file_path, '("snapshots_hdf5/snap",i5.5,".h5")') istep/msnap
      write(file_path2,'("snapshots_hdf5/AllPlanes/P",i3.3,"_snap",i5.5,".h5")') myrank_toroidal, istep/msnap
      call system("cp " // file_path // file_path2)
   endif


   deallocate(intp_a_1D)
   deallocate(intp_mean_full_dphi, right_intp_mean_full_dphi, left_intp_mean_full_dphi)
   deallocate(intp_FLA_dphi, right_intp_FLA_dphi, left_intp_FLA_dphi)
   deallocate(intp_mean_dE_cov3, right_mean_dE_cov3, left_mean_dE_cov3)
   deallocate(intp_FLA_dE_cov3, right_FLA_dE_cov3, left_FLA_dE_cov3)

   if(myrank_partd==nproc_partd-1) then
      deallocate(test_dphi_3D)
   endif

   if((myrank_partd==nproc_partd-1) .and. (mod(istep,msnap) .eq. 0)) then
      deallocate(DG_FLA_eScaData, DG_FLA_iScaData)
   endif



   ! ! ===========================
   ! ! Following part is for 3D full data
   ! ! ===========================
   ! ! Adiitional Information
   ! ! REcord Potential. Each myrank_toroidal has its own file (see open_file).
   ! if(movie_3D_h5==0) then  ! no 3D movie data created
   !   !ZHIXIN use flag to output PHI
   !   !     if(istep>0 .and. myrank_partd==0) then
   !   if(istep>0 .and. myrank_partd==0.and.flag_PHI>0) then
   !     write(777) istep
   !     write(777) phi(1,1:mgrid)  ! Note(5/08): phi contains zonal flow when ZF_SEPARATION=0 (solving total potential)
   !     call flush(777)  ! otherwise, if job crashes PHI data is not completed, which cause problem for reading PHI
   !   endif
   ! endif

   ! ! Record Density. Each myrank_toroidal has its own file (see open_file).
   ! !ZHIXIN use flag to output DEN
   ! if(istep>0 .and. myrank_partd==1.and.flag_DEN>0) then
   !   write(888) istep
   !   if(ifelectron==0) then
   !     write(888) densityi(1,1:mgrid),phi_tmp(1:mgrid)
   !   else
   !     !   Total electron density fluctuations = adiabatic + nonadiabatic
   !     !$omp parallel do private(i,j,ij)
   !     do i=0,mpsi
   !       do j=1,mtheta(i)
   !         ij=igrid(i)+j
   !         densitye(1,ij)=densitye(1,ij)+phi(1,ij)/Te_grid(i)*abs(csc0)  ! total electron density fluctuation dn_e/n_0
   !       enddo
   !       densitye(1,igrid(i))=densitye(1,igrid(i)+mtheta(i))  ! validate densitye(:,igrid(i)
   !     enddo

   !     write(888) densityi(1,1:mgrid),phi_tmp(1:mgrid),densitye(1,1:mgrid),heatflux(1:mgrid)   ! delta-ne output added Oct. 2012
   !                             !  phi_tmp:  Ti fluctuations
   !                             !  heatflux: Te fluctuations
   !     !          print*, 'mype= ', mype, 'istep=', istep, phi_tmp(mgrid/2), heatflux(mgrid/2), densitye(1,mgrid/2), densityi(1,mgrid/2)
   !   endif
   !   call flush(888)
   ! endif




end subroutine snapshot_hdf5


subroutine snapshot_Magnetic_Information

   use precision
   use global_parameters
   use control_parameters
   use particle_decomp
   use field_array
   use hdf5
   use h5lt
   use Structured_Grid
   implicit none

   integer, parameter :: gq_jump=3, a_jump=1
   integer :: dim_a, dim_gq
   real(wp) :: intp_val, inv_dgq(0:mpsi)
   integer :: i,j,ij,k,jt,jm

   real(wp), dimension(:), allocatable :: a_,s_, gq_,F_,dFda_,dgFda_,dgYda_,R_,dRda_,dRdgq_,Z_,dZda_,dZdgq_,BB_,dBBda_,dBBdgq_,P_,T_
   real(wp), dimension(:), allocatable ::Jacobian_, esi_N_, esi_K_, g_aa_cov_, g_aq_cov_, g_qq_cov_, g_ff_cov_  ! 1D representation for 1D data(dims(1)) , 2D data(dims(1), dims(2)) or 3D data(dims(1), dims(2), dims(3))
   real(wp), dimension(:), allocatable :: D_, g_aa_cntv_, g_aq_cntv_, g_qq_cntv_, g_ff_cntv_
   real(wp), dimension(:,:), allocatable :: B0_cov3_, B0_cntv3_, j0_cov3_, j0_cntv3_
   real(wp), dimension(:), allocatable :: ini_ga_, ini_dga0dt_ ! initial magnetic perturbation
   real(wp), dimension(:,:), allocatable :: ini_dB_cov3_, ini_dB_cntv3_, ini_grad_ga_cov3_! initial magnetic perturbation at t=0
   real(wp), dimension(:,:), allocatable :: uvec_B_cov, uvec_N_cov, uvec_BxN_cov
   real(wp), dimension(:,:), allocatable :: uvec_B_cntv, uvec_N_cntv, uvec_BxN_cntv

   character(LEN=128) :: file_path
   integer(HSIZE_T) :: h5dims(3)
   integer(HID_T) :: file_id, gid_MP, gid_SG! gid_MP: group id for MP (Magnetic Perturbation)
   integer :: h5err
   integer Naq

   inv_dgq=1.0/deltat
   jm=mtheta(ipsi)/gq_jump
   dim_a = (mpsi+1)/a_jump
   dim_gq = jm+1

   Naq=dim_a*dim_gq

   if(mype == numberpe-1)then !mype=0 must also have myrank_toroidal=0


      allocate(a_(Naq))
      allocate(s_(Naq),gq_(Naq),F_(Naq),dFda_(Naq),dgFda_(Naq),dgYda_(Naq),R_(Naq),dRda_(Naq),dRdgq_(Naq),Z_(Naq),dZda_(Naq),dZdgq_(Naq),BB_(Naq),dBBda_(Naq),dBBdgq_(Naq),P_(Naq),T_(Naq))
      allocate(Jacobian_(Naq), esi_N_(Naq),esi_K_(Naq),g_aa_cov_(Naq),g_aq_cov_(Naq),g_qq_cov_(Naq),g_ff_cov_(Naq))
      allocate(D_(Naq), g_aa_cntv_(Naq), g_aq_cntv_(Naq), g_qq_cntv_(Naq), g_ff_cntv_(Naq))
      allocate(B0_cov3_(3,Naq))
      allocate(B0_cntv3_(3,Naq),j0_cov3_(3,Naq),j0_cntv3_(3,Naq))
      allocate(uvec_B_cov(3,Naq), uvec_N_cov(3,Naq), uvec_BxN_cov(3,Naq))
      allocate(uvec_B_cntv(3,Naq), uvec_N_cntv(3,Naq), uvec_BxN_cntv(3,Naq))

      if(electromagnetic==1) then
         allocate(ini_ga_(Naq),ini_dga0dt_(Naq))
         allocate(ini_dB_cov3_(3,Naq),ini_dB_cntv3_(3,Naq), ini_grad_ga_cov3_(3,Naq))
      endif

      k=1
      do i=0,mpsi,a_jump
         do j=0,jm
            a_(k)=psi_grid(i)
            s_(k)=s_grid(i)
            gq_(k)=2.0*pi*real(j)/real(jm)
            k=k+1
         enddo
      enddo
      call esiget2dthese_more(F_,dFda_,dgFda_,dgYda_,R_,dRda_,dRdgq_,Z_,dZda_,dZdgq_,BB_,dBBda_,dBBdgq_,P_,T_,a_,gq_, dim_gq*dim_a)

      Jacobian_ = R_*(dZda_*dRdgq_-dZdgq_*dRda_)
      esi_N_ = (dRda_*dRdgq_+dZda_*dZdgq_)/Jacobian_
      esi_K_ = (dRdgq_*dRdgq_+dZdgq_*dZdgq_)/Jacobian_

      B0_cov3_(1,:) = -esi_N_*dgYda_
      B0_cov3_(2,:) = -esi_K_*dgYda_
      B0_cov3_(3,:) = F_

      B0_cntv3_(1,:) = 0.0
      B0_cntv3_(2,:) = -dgYda_/Jacobian_
      B0_cntv3_(3,:) = F_/(R_**2)

      j0_cntv3_(1,:) = 0.0
      j0_cntv3_(2,:) = -dFda_/Jacobian_
      j0_cntv3_(3,:) = P_ + T_/(R_**2)

      g_aa_cov_   = dRda_*dRda_+dZda_*dZda_;
      g_aq_cov_  = dRda_*dRdgq_+dZda_*dZdgq_;
      g_qq_cov_ = dRdgq_*dRdgq_+dZdgq_*dZdgq_;
      g_ff_cov_ = R_**2

      D_=dZda_*dRdgq_-dZdgq_*dRda_

      g_aa_cntv_=(dRdgq_**2+dZdgq_**2)/D_**2   ! grada^2
      g_qq_cntv_=(dZda_**2+dRda_**2)/D_**2   ! gradtheta^2
      g_ff_cntv_=1.0_wp/R_**2                 ! gradphi^2
      g_aq_cntv_=(-dRda_*dRdgq_-dZda_*dZdgq_)/D_**2  ! gradagradtheta

      j0_cov3_(1,:) = g_aa_cov_*j0_cntv3_(1,:) + g_aq_cov_*j0_cntv3_(2,:)
      j0_cov3_(2,:) = g_aq_cov_*j0_cntv3_(1,:) + g_qq_cov_*j0_cntv3_(2,:)
      j0_cov3_(3,:) = g_ff_cov_*j0_cntv3_(3,:)

      do k=1,size(a_)
         call Cal_unit_vectors_tangential_normal_binormal_to_FluxSurface(a_(k), gq_(k), zetamax,s_(k),uvec_B_cov(:,k),uvec_B_cntv(:,k),uvec_N_cov(:,k),uvec_N_cntv(:,k),uvec_BxN_cov(:,k),uvec_BxN_cntv(:,k))
      enddo

      if(electromagnetic==1)then
         ! Get initial magnetic perturbation information
         do k=1,size(a_)
            call Get_MP_alpha_amplitude(a_(k), gq_(k),modulo(zetamax,2.0*pi),0.0, ini_ga_(k), ini_grad_ga_cov3_(1:3,k), ini_dga0dt_(k))
         enddo

         ini_dB_cntv3_(1,:) = ini_ga_*j0_cntv3_(1,:) + (ini_grad_ga_cov3_(2,:)*B0_cov3_(3,:) - ini_grad_ga_cov3_(3,:)*B0_cov3_(2,:))/Jacobian_
         ini_dB_cntv3_(2,:) = ini_ga_*j0_cntv3_(2,:) + (ini_grad_ga_cov3_(3,:)*B0_cov3_(1,:) - ini_grad_ga_cov3_(1,:)*B0_cov3_(3,:))/Jacobian_
         ini_dB_cntv3_(3,:) = ini_ga_*j0_cntv3_(3,:) + (ini_grad_ga_cov3_(1,:)*B0_cov3_(2,:) - ini_grad_ga_cov3_(2,:)*B0_cov3_(1,:))/Jacobian_

         ini_dB_cov3_(1,:) = g_aa_cov_*ini_dB_cntv3_(1,:) + g_aq_cov_*ini_dB_cntv3_(2,:)
         ini_dB_cov3_(2,:) = g_aq_cov_*ini_dB_cntv3_(1,:) + g_qq_cov_*ini_dB_cntv3_(2,:)
         ini_dB_cov3_(3,:) = g_ff_cov_*ini_dB_cntv3_(3,:)
      endif


      ! Record Information using hdf5 format
      call system("mkdir -p " // 'snapshots_hdf5')
      write(file_path, '("SnapMagneticInfo.h5")')

      call h5open_f(h5err)
      call h5fcreate_f(file_path, H5F_ACC_TRUNC_F, file_id, h5err)

      ! 0D data
      h5dims(1) = 1
      call h5ltmake_dataset_f(file_id, "mpsi",1,h5dims,H5T_NATIVE_INTEGER, mpsi/a_jump, h5err)
      call h5ltmake_dataset_f(file_id, "mtheta",1,h5dims,H5T_NATIVE_INTEGER, jm, h5err)
      call h5ltmake_dataset_f(file_id, "zeta",1,h5dims,H5T_NATIVE_DOUBLE, zetamax, h5err)

      ! 2D data
      h5dims(1:2)=(/dim_gq, dim_a/)
      call h5ltmake_dataset_f(file_id, "a",2,h5dims,H5T_NATIVE_DOUBLE, a_,h5err)
      call h5ltmake_dataset_f(file_id, "gq",2,h5dims,H5T_NATIVE_DOUBLE, gq_,h5err)
      call h5ltmake_dataset_f(file_id, "F",2,h5dims,H5T_NATIVE_DOUBLE, F_,h5err)
      call h5ltmake_dataset_f(file_id, "dFda",2,h5dims,H5T_NATIVE_DOUBLE, dFda_,h5err)
      call h5ltmake_dataset_f(file_id, "dgFda",2,h5dims,H5T_NATIVE_DOUBLE, dgFda_,h5err)
      call h5ltmake_dataset_f(file_id, "dgYda",2,h5dims,H5T_NATIVE_DOUBLE, dgYda_,h5err)
      call h5ltmake_dataset_f(file_id, "R",2,h5dims,H5T_NATIVE_DOUBLE, R_,h5err)
      call h5ltmake_dataset_f(file_id, "dRda",2,h5dims,H5T_NATIVE_DOUBLE, dRda_,h5err)
      call h5ltmake_dataset_f(file_id, "dRdgq",2,h5dims,H5T_NATIVE_DOUBLE, dRdgq_,h5err)
      call h5ltmake_dataset_f(file_id, "Z",2,h5dims,H5T_NATIVE_DOUBLE, Z_,h5err)
      call h5ltmake_dataset_f(file_id, "dZda",2,h5dims,H5T_NATIVE_DOUBLE, dZda_,h5err)
      call h5ltmake_dataset_f(file_id, "dZdgq",2,h5dims,H5T_NATIVE_DOUBLE, dZdgq_,h5err)
      call h5ltmake_dataset_f(file_id, "BB",2,h5dims,H5T_NATIVE_DOUBLE, BB_,h5err)
      call h5ltmake_dataset_f(file_id, "dBBda",2,h5dims,H5T_NATIVE_DOUBLE, dBBda_,h5err)
      call h5ltmake_dataset_f(file_id, "dBBdgq",2,h5dims,H5T_NATIVE_DOUBLE, dBBdgq_,h5err)
      call h5ltmake_dataset_f(file_id, "P",2,h5dims,H5T_NATIVE_DOUBLE, P_,h5err)
      call h5ltmake_dataset_f(file_id, "T",2,h5dims,H5T_NATIVE_DOUBLE, T_,h5err)
      call h5ltmake_dataset_f(file_id, "Jacobian",2,h5dims,H5T_NATIVE_DOUBLE, Jacobian_,h5err)
      call h5ltmake_dataset_f(file_id, "esi_N",2,h5dims,H5T_NATIVE_DOUBLE, esi_N_,h5err)
      call h5ltmake_dataset_f(file_id, "esi_K",2,h5dims,H5T_NATIVE_DOUBLE, esi_K_,h5err)
      call h5ltmake_dataset_f(file_id, "g_aa_cov",2,h5dims,H5T_NATIVE_DOUBLE, g_aa_cov_,h5err)
      call h5ltmake_dataset_f(file_id, "g_aq_cov",2,h5dims,H5T_NATIVE_DOUBLE, g_aq_cov_,h5err)
      call h5ltmake_dataset_f(file_id, "g_qq_cov",2,h5dims,H5T_NATIVE_DOUBLE, g_qq_cov_,h5err)
      call h5ltmake_dataset_f(file_id, "g_ff_cov",2,h5dims,H5T_NATIVE_DOUBLE, g_ff_cov_,h5err)
      call h5ltmake_dataset_f(file_id, "g_aa_cntv",2,h5dims,H5T_NATIVE_DOUBLE, g_aa_cntv_,h5err)
      call h5ltmake_dataset_f(file_id, "g_aq_cntv",2,h5dims,H5T_NATIVE_DOUBLE, g_aq_cntv_,h5err)
      call h5ltmake_dataset_f(file_id, "g_qq_cntv",2,h5dims,H5T_NATIVE_DOUBLE, g_qq_cntv_,h5err)
      call h5ltmake_dataset_f(file_id, "g_ff_cntv",2,h5dims,H5T_NATIVE_DOUBLE, g_ff_cntv_,h5err)
      call h5ltmake_dataset_f(file_id, "Jaco_over_R",2,h5dims,H5T_NATIVE_DOUBLE, D_,h5err)


      call h5ltset_attribute_string_f(file_id, "a", "unit", "norm_psi_t", h5err)
      call h5ltset_attribute_string_f(file_id, "gq", "unit", "radians", h5err)
      call h5ltset_attribute_string_f(file_id, "R", "unit", "m", h5err)
      call h5ltset_attribute_string_f(file_id, "Z", "unit", "m", h5err)
      call h5ltset_attribute_string_f(file_id, "BB", "unit", "T", h5err)
      call h5ltset_attribute_string_f(file_id, "BB", "description", "magnitude of equilibrium B field (B0)", h5err)

      ! 3D data (vectors on 2D grids)
      h5dims(1:3)=(/3, dim_gq, dim_a/)
      call h5ltmake_dataset_f(file_id, "B0_cov",3,h5dims,H5T_NATIVE_DOUBLE, B0_cov3_,h5err)
      call h5ltmake_dataset_f(file_id, "B0_cntv",3,h5dims,H5T_NATIVE_DOUBLE, B0_cntv3_,h5err)
      call h5ltmake_dataset_f(file_id, "j0_cov",3,h5dims,H5T_NATIVE_DOUBLE, j0_cov3_,h5err)
      call h5ltmake_dataset_f(file_id, "j0_cntv",3,h5dims,H5T_NATIVE_DOUBLE, j0_cntv3_,h5err)
      call h5ltmake_dataset_f(file_id, "uvec_B_cov",3,h5dims,H5T_NATIVE_DOUBLE, uvec_B_cov,h5err)
      call h5ltmake_dataset_f(file_id, "uvec_N_cov",3,h5dims,H5T_NATIVE_DOUBLE, uvec_N_cov,h5err)
      call h5ltmake_dataset_f(file_id, "uvec_BxN_cov",3,h5dims,H5T_NATIVE_DOUBLE, uvec_BxN_cov,h5err)
      call h5ltmake_dataset_f(file_id, "uvec_B_cntv",3,h5dims,H5T_NATIVE_DOUBLE, uvec_B_cntv,h5err)
      call h5ltmake_dataset_f(file_id, "uvec_N_cntv",3,h5dims,H5T_NATIVE_DOUBLE, uvec_N_cntv,h5err)
      call h5ltmake_dataset_f(file_id, "uvec_BxN_cntv",3,h5dims,H5T_NATIVE_DOUBLE, uvec_BxN_cntv,h5err)



      if(electromagnetic==1) then
         call h5gcreate_f(file_id, "iniMP", gid_MP, h5err)
         ! 2D data
         h5dims(1:2)=(/dim_gq, dim_a/)
         call h5ltmake_dataset_f(file_id, "iniMP/ga",2,h5dims,H5T_NATIVE_DOUBLE, ini_ga_,h5err)
         call h5ltmake_dataset_f(file_id, "iniMP/dga0dt",2,h5dims,H5T_NATIVE_DOUBLE, ini_dga0dt_,h5err)

         ! 3D data (vectors on 2D grids)
         h5dims(1:3)=(/3, dim_gq, dim_a/)
         call h5ltmake_dataset_f(file_id, "iniMP/dB_cov",3,h5dims,H5T_NATIVE_DOUBLE, ini_dB_cov3_,h5err)
         call h5ltmake_dataset_f(file_id, "iniMP/dB_cntv",3,h5dims,H5T_NATIVE_DOUBLE, ini_dB_cntv3_,h5err)
         call h5ltmake_dataset_f(file_id, "iniMP/grad_ga_cov",3,h5dims,H5T_NATIVE_DOUBLE, ini_grad_ga_cov3_,h5err)


         call h5gclose_f(gid_MP, h5err)
      endif


      if(Flag_Poisson3D_on_SG>=1) then
         call h5gcreate_f(file_id, "SG", gid_SG, h5err)

         h5dims(1:2)=(/SG%dim_gq, SG%dim_s/)
         call h5ltmake_dataset_f(gid_SG, "a_2D",2,h5dims,H5T_NATIVE_DOUBLE, SG%a_2D,h5err)
         call h5ltmake_dataset_f(gid_SG, "gq_2D",2,h5dims,H5T_NATIVE_DOUBLE, SG%gq_2D,h5err)
         call h5ltmake_dataset_f(gid_SG, "F",2,h5dims,H5T_NATIVE_DOUBLE, SG%F,h5err)
         call h5ltmake_dataset_f(gid_SG, "dFda",2,h5dims,H5T_NATIVE_DOUBLE, SG%dFda,h5err)
         call h5ltmake_dataset_f(gid_SG, "dgFda",2,h5dims,H5T_NATIVE_DOUBLE, SG%dgFda,h5err)
         call h5ltmake_dataset_f(gid_SG, "dgYda",2,h5dims,H5T_NATIVE_DOUBLE, SG%dgYda,h5err)
         call h5ltmake_dataset_f(gid_SG, "R",2,h5dims,H5T_NATIVE_DOUBLE, SG%R,h5err)
         call h5ltmake_dataset_f(gid_SG, "dRda",2,h5dims,H5T_NATIVE_DOUBLE, SG%dRda,h5err)
         call h5ltmake_dataset_f(gid_SG, "dRdgq",2,h5dims,H5T_NATIVE_DOUBLE, SG%dRdgq,h5err)
         call h5ltmake_dataset_f(gid_SG, "Z",2,h5dims,H5T_NATIVE_DOUBLE, SG%Z,h5err)
         call h5ltmake_dataset_f(gid_SG, "dZda",2,h5dims,H5T_NATIVE_DOUBLE, SG%dZda,h5err)
         call h5ltmake_dataset_f(gid_SG, "dZdgq",2,h5dims,H5T_NATIVE_DOUBLE, SG%dZdgq,h5err)
         call h5ltmake_dataset_f(gid_SG, "BB",2,h5dims,H5T_NATIVE_DOUBLE, SG%BB,h5err)
         call h5ltmake_dataset_f(gid_SG, "dBBda",2,h5dims,H5T_NATIVE_DOUBLE, SG%dBBda,h5err)
         call h5ltmake_dataset_f(gid_SG, "dBBdgq",2,h5dims,H5T_NATIVE_DOUBLE, SG%dBBdgq,h5err)
         call h5ltmake_dataset_f(gid_SG, "P",2,h5dims,H5T_NATIVE_DOUBLE, SG%P,h5err)
         call h5ltmake_dataset_f(gid_SG, "T",2,h5dims,H5T_NATIVE_DOUBLE, SG%T,h5err)
         call h5ltmake_dataset_f(gid_SG, "Jacobian",2,h5dims,H5T_NATIVE_DOUBLE, SG%Jacobian,h5err)
         call h5ltmake_dataset_f(gid_SG, "esi_N",2,h5dims,H5T_NATIVE_DOUBLE, SG%esi_N,h5err)
         call h5ltmake_dataset_f(gid_SG, "esi_K",2,h5dims,H5T_NATIVE_DOUBLE, SG%esi_K,h5err)
         call h5ltmake_dataset_f(gid_SG, "g_aa_cov",2,h5dims,H5T_NATIVE_DOUBLE, SG%g_aa_cov,h5err)
         call h5ltmake_dataset_f(gid_SG, "g_aq_cov",2,h5dims,H5T_NATIVE_DOUBLE, SG%g_aq_cov,h5err)
         call h5ltmake_dataset_f(gid_SG, "g_qq_cov",2,h5dims,H5T_NATIVE_DOUBLE, SG%g_qq_cov,h5err)
         call h5ltmake_dataset_f(gid_SG, "g_ff_cov",2,h5dims,H5T_NATIVE_DOUBLE, SG%g_ff_cov,h5err)
         call h5ltmake_dataset_f(gid_SG, "g_aa_cntv",2,h5dims,H5T_NATIVE_DOUBLE, SG%g_aa_cntv,h5err)
         call h5ltmake_dataset_f(gid_SG, "g_aq_cntv",2,h5dims,H5T_NATIVE_DOUBLE, SG%g_aq_cntv,h5err)
         call h5ltmake_dataset_f(gid_SG, "g_qq_cntv",2,h5dims,H5T_NATIVE_DOUBLE, SG%g_qq_cntv,h5err)
         call h5ltmake_dataset_f(gid_SG, "g_ff_cntv",2,h5dims,H5T_NATIVE_DOUBLE, SG%g_ff_cntv,h5err)
         call h5ltmake_dataset_f(gid_SG, "Jaco_over_R",2,h5dims,H5T_NATIVE_DOUBLE, SG%D,h5err)

         call h5ltmake_dataset_f(gid_SG, "Jacob_sqf",2,h5dims,H5T_NATIVE_DOUBLE, SG%Jacob_sqf,h5err)
         call h5ltmake_dataset_f(gid_SG, "g_ss_cov",2,h5dims,H5T_NATIVE_DOUBLE, SG%g_ss_cov,h5err)
         call h5ltmake_dataset_f(gid_SG, "g_sq_cov",2,h5dims,H5T_NATIVE_DOUBLE, SG%g_sq_cov,h5err)
         call h5ltmake_dataset_f(gid_SG, "g_ss_cntv",2,h5dims,H5T_NATIVE_DOUBLE, SG%g_ss_cntv,h5err)
         call h5ltmake_dataset_f(gid_SG, "g_sq_cntv",2,h5dims,H5T_NATIVE_DOUBLE, SG%g_sq_cntv,h5err)

         call h5ltmake_dataset_f(gid_SG, "Cell_dV",2,h5dims,H5T_NATIVE_DOUBLE, SG%Cell_dV,h5err)
         call h5ltmake_dataset_f(gid_SG, "Cell_dS_s",2,h5dims,H5T_NATIVE_DOUBLE, SG%Cell_dS_s,h5err)
         call h5ltmake_dataset_f(gid_SG, "Cell_dS_gq",2,h5dims,H5T_NATIVE_DOUBLE, SG%Cell_dS_gq,h5err)
         call h5ltmake_dataset_f(gid_SG, "Cell_dS_gf",2,h5dims,H5T_NATIVE_DOUBLE, SG%Cell_dS_gf,h5err)

         call h5gclose_f(gid_SG, h5err)
      end if


      call h5fclose_f(file_id, h5err)
      call h5close_f(h5err)


      deallocate(a_,s_,gq_,F_,dFda_,dgFda_,dgYda_,R_,dRda_,dRdgq_,Z_,dZda_,dZdgq_,BB_,dBBda_,dBBdgq_,P_,T_)
      deallocate(Jacobian_, esi_N_, esi_K_, g_aa_cov_, g_aq_cov_, g_qq_cov_, g_ff_cov_)
      deallocate(D_, g_aa_cntv_, g_aq_cntv_, g_qq_cntv_, g_ff_cntv_)
      deallocate(B0_cov3_, B0_cntv3_, j0_cov3_, j0_cntv3_)
      deallocate(uvec_B_cov, uvec_N_cov, uvec_BxN_cov)
      deallocate(uvec_B_cntv, uvec_N_cntv, uvec_BxN_cntv)
      if(electromagnetic==1)then
         deallocate(ini_ga_,ini_dga0dt_)
         deallocate(ini_dB_cov3_,ini_dB_cntv3_, ini_grad_ga_cov3_)
      endif

   endif


end subroutine snapshot_Magnetic_Information


subroutine Print_FLA_Denominator


   use global_parameters
   use precision
   use global_parameters
   use particle_decomp
   use FieldLineFollowing
   use field_array

   implicit none

   integer i, j,ij
   character(len=32) chdum

   write(chdum,'("FLA_Denominator_",i5.5,".dat")') myrank_toroidal
   open(181206,file=chdum,status='replace',action='write')
   do i=0,mpsi
      do j=0, mtheta(i)
         ij=igrid(i) + j
         write(181206, "((ES16.8,1x))") FLA_Denominator(1,ij)
      enddo
   enddo

   close(181206)
end subroutine Print_FLA_Denominator


subroutine Cal_dni_dTi_from_each_ion_weight_component

   use precision
   use mpi_defs
   use global_parameters
   use particle_array
   use field_array
   use particle_decomp
   use esi_interface

   implicit none

   integer m,i,j,k,jt,ii,ierror,ij, istatus(MPI_STATUS_SIZE)
   integer larmor

   real(wp) wa0,wa1,inner_wgq0, inner_wgq1, outter_wgq0, outter_wgq1, wgf0, wgf1
   real(wp) tmp1, tmp2
   real(wp) weight, vpara, energy

   real(wp) recv_array(0:1,mgrid)

   dni_init=0.0; dni_vpara=0.0; dni_vperp=0.0; dni_accel=0.0; dni_fdamp=0.0;
   dTi_init=0.0; dTi_vpara=0.0; dTi_vperp=0.0; dTi_accel=0.0; dTi_fdamp=0.0;

   !$acc parallel loop present(zion,jtion0,jtion1,wzion,wpion,wtion0,wtion1)
   do m=1,mi

      if (Flag_w_Threshold==1 .and. (zion(5,m)<w_Thold_min .or. zion(5,m)>w_Thold_max)) cycle

      wgf1=wzion(m)
      wgf0=1.0-wgf1

      !! Caution!! Is backu(m) here omega or real velocity?? Check it! Refer pushion routine
      energy = 0.5*(zion(4,m)*BB(m)-backu(m))**2.0 + ((zion(6,m)**2.0)*BB(m))

      !$acc loop seq
      do larmor=1, larmor_num
         wa1=wpion(larmor,m)
         wa0=1.0-wa1
         inner_wgq1=wtion0(larmor,m)
         inner_wgq0=1.0-inner_wgq1
         outter_wgq1=wtion1(larmor,m)
         outter_wgq0=1.0-outter_wgq1

         !inner radius, smaller theta index
         ij=jtion0(larmor,m)
         tmp1 = wgf0*wa0*inner_wgq0
         dni_init(0,ij)  = dni_init(0,ij)  + tmp1*wi_init(m)
         dni_vpara(0,ij) = dni_vpara(0,ij) + tmp1*wi_vpara(m)
         dni_vperp(0,ij) = dni_vperp(0,ij) + tmp1*wi_vperp(m)
         dni_accel(0,ij) = dni_accel(0,ij) + tmp1*wi_accel(m)
         dni_fdamp(0,ij) = dni_fdamp(0,ij) + tmp1*wi_fdamp(m)
         tmp2 = tmp1*((2.0/3.0)*energy-backT(m))
         dTi_init(0,ij)  = dTi_init(0,ij)  + tmp2*wi_init(m)
         dTi_vpara(0,ij) = dTi_vpara(0,ij) + tmp2*wi_vpara(m)
         dTi_vperp(0,ij) = dTi_vperp(0,ij) + tmp2*wi_vperp(m)
         dTi_accel(0,ij) = dTi_accel(0,ij) + tmp2*wi_accel(m)
         dTi_fdamp(0,ij) = dTi_fdamp(0,ij) + tmp2*wi_fdamp(m)
         tmp1 = wgf1*wa0*inner_wgq0
         dni_init(1,ij)  = dni_init(1,ij)  + tmp1*wi_init(m)
         dni_vpara(1,ij) = dni_vpara(1,ij) + tmp1*wi_vpara(m)
         dni_vperp(1,ij) = dni_vperp(1,ij) + tmp1*wi_vperp(m)
         dni_accel(1,ij) = dni_accel(1,ij) + tmp1*wi_accel(m)
         dni_fdamp(1,ij) = dni_fdamp(1,ij) + tmp1*wi_fdamp(m)
         tmp2 = tmp1*((2.0/3.0)*energy-backT(m))
         dTi_init(1,ij)  = dTi_init(1,ij)  + tmp2*wi_init(m)
         dTi_vpara(1,ij) = dTi_vpara(1,ij) + tmp2*wi_vpara(m)
         dTi_vperp(1,ij) = dTi_vperp(1,ij) + tmp2*wi_vperp(m)
         dTi_accel(1,ij) = dTi_accel(1,ij) + tmp2*wi_accel(m)
         dTi_fdamp(1,ij) = dTi_fdamp(1,ij) + tmp2*wi_fdamp(m)

         !inner radius, larger theta index
         ij=ij+1
         tmp1 = wgf0*wa0*inner_wgq1
         dni_init(0,ij)  = dni_init(0,ij)  + tmp1*wi_init(m)
         dni_vpara(0,ij) = dni_vpara(0,ij) + tmp1*wi_vpara(m)
         dni_vperp(0,ij) = dni_vperp(0,ij) + tmp1*wi_vperp(m)
         dni_accel(0,ij) = dni_accel(0,ij) + tmp1*wi_accel(m)
         dni_fdamp(0,ij) = dni_fdamp(0,ij) + tmp1*wi_fdamp(m)
         tmp2 = tmp1*((2.0/3.0)*energy-backT(m))
         dTi_init(0,ij)  = dTi_init(0,ij)  + tmp2*wi_init(m)
         dTi_vpara(0,ij) = dTi_vpara(0,ij) + tmp2*wi_vpara(m)
         dTi_vperp(0,ij) = dTi_vperp(0,ij) + tmp2*wi_vperp(m)
         dTi_accel(0,ij) = dTi_accel(0,ij) + tmp2*wi_accel(m)
         dTi_fdamp(0,ij) = dTi_fdamp(0,ij) + tmp2*wi_fdamp(m)
         tmp1 = wgf1*wa0*inner_wgq1
         dni_init(1,ij)  = dni_init(1,ij)  + tmp1*wi_init(m)
         dni_vpara(1,ij) = dni_vpara(1,ij) + tmp1*wi_vpara(m)
         dni_vperp(1,ij) = dni_vperp(1,ij) + tmp1*wi_vperp(m)
         dni_accel(1,ij) = dni_accel(1,ij) + tmp1*wi_accel(m)
         dni_fdamp(1,ij) = dni_fdamp(1,ij) + tmp1*wi_fdamp(m)
         tmp2 = tmp1*((2.0/3.0)*energy-backT(m))
         dTi_init(1,ij)  = dTi_init(1,ij)  + tmp2*wi_init(m)
         dTi_vpara(1,ij) = dTi_vpara(1,ij) + tmp2*wi_vpara(m)
         dTi_vperp(1,ij) = dTi_vperp(1,ij) + tmp2*wi_vperp(m)
         dTi_accel(1,ij) = dTi_accel(1,ij) + tmp2*wi_accel(m)
         dTi_fdamp(1,ij) = dTi_fdamp(1,ij) + tmp2*wi_fdamp(m)

         !outter radius, smaller theta index
         ij=jtion1(larmor,m)
         tmp1 = wgf0*wa1*outter_wgq0
         dni_init(0,ij)  = dni_init(0,ij)  + tmp1*wi_init(m)
         dni_vpara(0,ij) = dni_vpara(0,ij) + tmp1*wi_vpara(m)
         dni_vperp(0,ij) = dni_vperp(0,ij) + tmp1*wi_vperp(m)
         dni_accel(0,ij) = dni_accel(0,ij) + tmp1*wi_accel(m)
         dni_fdamp(0,ij) = dni_fdamp(0,ij) + tmp1*wi_fdamp(m)
         tmp2 = tmp1*((2.0/3.0)*energy-backT(m))
         dTi_init(0,ij)  = dTi_init(0,ij)  + tmp2*wi_init(m)
         dTi_vpara(0,ij) = dTi_vpara(0,ij) + tmp2*wi_vpara(m)
         dTi_vperp(0,ij) = dTi_vperp(0,ij) + tmp2*wi_vperp(m)
         dTi_accel(0,ij) = dTi_accel(0,ij) + tmp2*wi_accel(m)
         dTi_fdamp(0,ij) = dTi_fdamp(0,ij) + tmp2*wi_fdamp(m)
         tmp1 = wgf1*wa1*outter_wgq0
         dni_init(1,ij)  = dni_init(1,ij)  + tmp1*wi_init(m)
         dni_vpara(1,ij) = dni_vpara(1,ij) + tmp1*wi_vpara(m)
         dni_vperp(1,ij) = dni_vperp(1,ij) + tmp1*wi_vperp(m)
         dni_accel(1,ij) = dni_accel(1,ij) + tmp1*wi_accel(m)
         dni_fdamp(1,ij) = dni_fdamp(1,ij) + tmp1*wi_fdamp(m)
         tmp2 = tmp1*((2.0/3.0)*energy-backT(m))
         dTi_init(1,ij)  = dTi_init(1,ij)  + tmp2*wi_init(m)
         dTi_vpara(1,ij) = dTi_vpara(1,ij) + tmp2*wi_vpara(m)
         dTi_vperp(1,ij) = dTi_vperp(1,ij) + tmp2*wi_vperp(m)
         dTi_accel(1,ij) = dTi_accel(1,ij) + tmp2*wi_accel(m)
         dTi_fdamp(1,ij) = dTi_fdamp(1,ij) + tmp2*wi_fdamp(m)

         !outter radius, larger theta index
         ij=ij+1
         tmp1 = wgf0*wa1*outter_wgq1
         dni_init(0,ij)  = dni_init(0,ij)  + tmp1*wi_init(m)
         dni_vpara(0,ij) = dni_vpara(0,ij) + tmp1*wi_vpara(m)
         dni_vperp(0,ij) = dni_vperp(0,ij) + tmp1*wi_vperp(m)
         dni_accel(0,ij) = dni_accel(0,ij) + tmp1*wi_accel(m)
         dni_fdamp(0,ij) = dni_fdamp(0,ij) + tmp1*wi_fdamp(m)
         tmp2 = tmp1*((2.0/3.0)*energy-backT(m))
         dTi_init(0,ij)  = dTi_init(0,ij)  + tmp2*wi_init(m)
         dTi_vpara(0,ij) = dTi_vpara(0,ij) + tmp2*wi_vpara(m)
         dTi_vperp(0,ij) = dTi_vperp(0,ij) + tmp2*wi_vperp(m)
         dTi_accel(0,ij) = dTi_accel(0,ij) + tmp2*wi_accel(m)
         dTi_fdamp(0,ij) = dTi_fdamp(0,ij) + tmp2*wi_fdamp(m)
         tmp1 = wgf1*wa1*outter_wgq1
         dni_init(1,ij)  = dni_init(1,ij)  + tmp1*wi_init(m)
         dni_vpara(1,ij) = dni_vpara(1,ij) + tmp1*wi_vpara(m)
         dni_vperp(1,ij) = dni_vperp(1,ij) + tmp1*wi_vperp(m)
         dni_accel(1,ij) = dni_accel(1,ij) + tmp1*wi_accel(m)
         dni_fdamp(1,ij) = dni_fdamp(1,ij) + tmp1*wi_fdamp(m)
         tmp2 = tmp1*((2.0/3.0)*energy-backT(m))
         dTi_init(1,ij)  = dTi_init(1,ij)  + tmp2*wi_init(m)
         dTi_vpara(1,ij) = dTi_vpara(1,ij) + tmp2*wi_vpara(m)
         dTi_vperp(1,ij) = dTi_vperp(1,ij) + tmp2*wi_vperp(m)
         dTi_accel(1,ij) = dTi_accel(1,ij) + tmp2*wi_accel(m)
         dTi_fdamp(1,ij) = dTi_fdamp(1,ij) + tmp2*wi_fdamp(m)
      enddo
   enddo


   if(npartdom>1) then
      call MPI_ALLREDUCE(dni_init, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dni_init=recv_array
      call MPI_ALLREDUCE(dni_vpara, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dni_vpara=recv_array
      call MPI_ALLREDUCE(dni_vperp, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dni_vperp=recv_array
      call MPI_ALLREDUCE(dni_accel, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dni_accel=recv_array
      call MPI_ALLREDUCE(dni_fdamp, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dni_fdamp=recv_array
      call MPI_ALLREDUCE(dTi_init, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dTi_init=recv_array
      call MPI_ALLREDUCE(dTi_vpara, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dTi_vpara=recv_array
      call MPI_ALLREDUCE(dTi_vperp, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dTi_vperp=recv_array
      call MPI_ALLREDUCE(dTi_accel, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dTi_accel=recv_array
      call MPI_ALLREDUCE(dTi_fdamp, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dTi_fdamp=recv_array
   endif

   ! poloidal end cell; need valid data for densityi(:,igrid(i))
   do i=0,mpsi

      dni_init(:,igrid(i)+mtheta(i))=dni_init(:,igrid(i)+mtheta(i))+dni_init(:,igrid(i))
      dni_init(:,igrid(i))=dni_init(:,igrid(i)+mtheta(i))
      dni_vpara(:,igrid(i)+mtheta(i))=dni_vpara(:,igrid(i)+mtheta(i))+dni_vpara(:,igrid(i))
      dni_vpara(:,igrid(i))=dni_vpara(:,igrid(i)+mtheta(i))
      dni_vperp(:,igrid(i)+mtheta(i))=dni_vperp(:,igrid(i)+mtheta(i))+dni_vperp(:,igrid(i))
      dni_vperp(:,igrid(i))=dni_vperp(:,igrid(i)+mtheta(i))
      dni_accel(:,igrid(i)+mtheta(i))=dni_accel(:,igrid(i)+mtheta(i))+dni_accel(:,igrid(i))
      dni_accel(:,igrid(i))=dni_accel(:,igrid(i)+mtheta(i))
      dni_fdamp(:,igrid(i)+mtheta(i))=dni_fdamp(:,igrid(i)+mtheta(i))+dni_fdamp(:,igrid(i))
      dni_fdamp(:,igrid(i))=dni_fdamp(:,igrid(i)+mtheta(i))

      dTi_init(:,igrid(i)+mtheta(i))=dTi_init(:,igrid(i)+mtheta(i))+dTi_init(:,igrid(i))
      dTi_init(:,igrid(i))=dTi_init(:,igrid(i)+mtheta(i))
      dTi_vpara(:,igrid(i)+mtheta(i))=dTi_vpara(:,igrid(i)+mtheta(i))+dTi_vpara(:,igrid(i))
      dTi_vpara(:,igrid(i))=dTi_vpara(:,igrid(i)+mtheta(i))
      dTi_vperp(:,igrid(i)+mtheta(i))=dTi_vperp(:,igrid(i)+mtheta(i))+dTi_vperp(:,igrid(i))
      dTi_vperp(:,igrid(i))=dTi_vperp(:,igrid(i)+mtheta(i))
      dTi_accel(:,igrid(i)+mtheta(i))=dTi_accel(:,igrid(i)+mtheta(i))+dTi_accel(:,igrid(i))
      dTi_accel(:,igrid(i))=dTi_accel(:,igrid(i)+mtheta(i))
      dTi_fdamp(:,igrid(i)+mtheta(i))=dTi_fdamp(:,igrid(i)+mtheta(i))+dTi_fdamp(:,igrid(i))
      dTi_fdamp(:,igrid(i))=dTi_fdamp(:,igrid(i)+mtheta(i))
   enddo


   ! send densityi to left and receive from right
   call MPI_SENDRECV(dni_init(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dni_init(1,ii+1:ii+jt)=dni_init(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dni_init(1,:)=dni_init(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dni_vpara(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dni_vpara(1,ii+1:ii+jt)=dni_vpara(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dni_vpara(1,:)=dni_vpara(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dni_vperp(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dni_vperp(1,ii+1:ii+jt)=dni_vperp(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dni_vperp(1,:)=dni_vperp(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dni_accel(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dni_accel(1,ii+1:ii+jt)=dni_accel(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dni_accel(1,:)=dni_accel(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dni_fdamp(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dni_fdamp(1,ii+1:ii+jt)=dni_fdamp(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dni_fdamp(1,:)=dni_fdamp(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dTi_init(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dTi_init(1,ii+1:ii+jt)=dTi_init(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dTi_init(1,:)=dTi_init(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dTi_vpara(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dTi_vpara(1,ii+1:ii+jt)=dTi_vpara(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dTi_vpara(1,:)=dTi_vpara(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dTi_vperp(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dTi_vperp(1,ii+1:ii+jt)=dTi_vperp(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dTi_vperp(1,:)=dTi_vperp(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dTi_accel(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dTi_accel(1,ii+1:ii+jt)=dTi_accel(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dTi_accel(1,:)=dTi_accel(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dTi_fdamp(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dTi_fdamp(1,ii+1:ii+jt)=dTi_fdamp(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dTi_fdamp(1,:)=dTi_fdamp(1,:)+recv_array(1,:)
   endif

   !===============================================================================================
   ! Caution. This routine does not enforce charge conservation regardless of its_jacobi parameter
   !===============================================================================================
   tmp1=1.0/real(larmor_num)
   do i=0,mpsi
      do j=1,mtheta(i)
         ij=igrid(i)+j
         dni_init(1,ij) =tmp1*dni_init(1,ij)*markeri(ij)
         dni_vpara(1,ij)=tmp1*dni_vpara(1,ij)*markeri(ij)
         dni_vperp(1,ij)=tmp1*dni_vperp(1,ij)*markeri(ij)
         dni_accel(1,ij)=tmp1*dni_accel(1,ij)*markeri(ij)
         dni_fdamp(1,ij)=tmp1*dni_fdamp(1,ij)*markeri(ij)

         dTi_init(1,ij) =tmp1*dTi_init(1,ij)*markeri(ij)
         dTi_vpara(1,ij)=tmp1*dTi_vpara(1,ij)*markeri(ij)
         dTi_vperp(1,ij)=tmp1*dTi_vperp(1,ij)*markeri(ij)
         dTi_accel(1,ij)=tmp1*dTi_accel(1,ij)*markeri(ij)
         dTi_fdamp(1,ij)=tmp1*dTi_fdamp(1,ij)*markeri(ij)
      enddo
      dni_init(1,igrid(i))=dni_init(1,igrid(i)+mtheta(i))
      dni_vpara(1,igrid(i))=dni_vpara(1,igrid(i)+mtheta(i))
      dni_vperp(1,igrid(i))=dni_vperp(1,igrid(i)+mtheta(i))
      dni_accel(1,igrid(i))=dni_accel(1,igrid(i)+mtheta(i))
      dni_fdamp(1,igrid(i))=dni_fdamp(1,igrid(i)+mtheta(i))

      dTi_init(1,igrid(i))=dTi_init(1,igrid(i)+mtheta(i))
      dTi_vpara(1,igrid(i))=dTi_vpara(1,igrid(i)+mtheta(i))
      dTi_vperp(1,igrid(i))=dTi_vperp(1,igrid(i)+mtheta(i))
      dTi_accel(1,igrid(i))=dTi_accel(1,igrid(i)+mtheta(i))
      dTi_fdamp(1,igrid(i))=dTi_fdamp(1,igrid(i)+mtheta(i))
   enddo

   !==============================================================================
   ! MGYOO (01/16/2019)
   ! Get real dne by multiplying ne0
   do i=0,mpsi
      do j=0,mtheta(i)
         ij = igrid(i) + j
         dni_init(1,ij)  = den_grid(i)*dni_init(1,ij)
         dni_vpara(1,ij) = den_grid(i)*dni_vpara(1,ij)
         dni_vperp(1,ij) = den_grid(i)*dni_vperp(1,ij)
         dni_accel(1,ij) = den_grid(i)*dni_accel(1,ij)
         dni_fdamp(1,ij) = den_grid(i)*dni_fdamp(1,ij)
      enddo
   enddo
   ! ! Convert dTe to get a correct form of dTe by considering dne.
   ! ! buff = Sigma(weight*(2/3*energy-T0))*markere = (1/n0)*Sigma(weight*(2/3*energy-T0))/node_Vol
   ! ! dTe = (1/(n0+dn))*Sigma(weight*(2/3*energy-T0))/node_Vol = n0/(n0+dn)*buff
   ! do i=0,mpsi
   !   do j=0,mtheta(i)
   !     ij = igrid(i) + j
   !     tmp1 =  den_grid(i)/(den_grid(i)+tot_dni(1,ij))
   !     dTi_init(1,ij)  = tmp1*dTi_init(1,ij)
   !     dTi_vpara(1,ij) = tmp1*dTi_vpara(1,ij)
   !     dTi_vperp(1,ij) = tmp1*dTi_vperp(1,ij)
   !     dTi_accel(1,ij) = tmp1*dTi_accel(1,ij)
   !     dTi_fdamp(1,ij) = tmp1*dTi_fdamp(1,ij)
   !   enddo
   ! enddo
   !==============================================================================

end subroutine Cal_dni_dTi_from_each_ion_weight_component




subroutine Cal_dne_dTe_from_each_electron_weight_component

   use precision
   use mpi_defs
   use global_parameters
   use particle_array
   use field_array
   use particle_decomp
   use esi_interface

   implicit none

   integer m,i,j,k,jt,ii,ierror,ij, istatus(MPI_STATUS_SIZE)

   real(wp) wa0,wa1,inner_wgq0, inner_wgq1, outter_wgq0, outter_wgq1, wgf0, wgf1
   real(wp) tmp1, tmp2
   real(wp) weight, vpara, energy

   real(wp) recv_array(0:1,mgrid)


   dne_init=0.0; dne_vpara=0.0; dne_vperp=0.0; dne_accel=0.0; dne_adiab=0.0; dne_fdamp=0.0;
   dTe_init=0.0; dTe_vpara=0.0; dTe_vperp=0.0; dTe_accel=0.0; dTe_adiab=0.0; dTe_fdamp=0.0;

   !$acc parallel loop present(zeon,jteon0,jteon1,wzeon,wpeon,wteon0,wteon1)
   do m=1,me

      if (Flag_w_Threshold==1 .and. (zeon(5,m)<w_Thold_min .or. zeon(5,m)>w_Thold_max)) cycle

      wgf1=wzeon(m)
      wgf0=1.0-wgf1

      energy = 0.5*for_E*(zeon(4,m)*BB(m))**2.0 + ((zeon(6,m)**2.0)*BB(m))


      wa1=wpeon(m)
      wa0=1.0-wa1
      inner_wgq1=wteon0(m)
      inner_wgq0=1.0-inner_wgq1
      outter_wgq1=wteon1(m)
      outter_wgq0=1.0-outter_wgq1

      !inner radius, smaller theta index
      ij=jteon0(m)
      tmp1 = wgf0*wa0*inner_wgq0
      dne_init(0,ij)  = dne_init(0,ij)  + tmp1*we_init(m)
      dne_vpara(0,ij) = dne_vpara(0,ij) + tmp1*we_vpara(m)
      dne_vperp(0,ij) = dne_vperp(0,ij) + tmp1*we_vperp(m)
      dne_accel(0,ij) = dne_accel(0,ij) + tmp1*we_accel(m)
      dne_adiab(0,ij) = dne_adiab(0,ij) + tmp1*we_adiab(m)
      dne_fdamp(0,ij) = dne_fdamp(0,ij) + tmp1*we_fdamp(m)
      tmp2 = tmp1*((2.0/3.0)*energy-backTe(m)*msm0)
      dTe_init(0,ij)  = dTe_init(0,ij)  + tmp2*we_init(m)
      dTe_vpara(0,ij) = dTe_vpara(0,ij) + tmp2*we_vpara(m)
      dTe_vperp(0,ij) = dTe_vperp(0,ij) + tmp2*we_vperp(m)
      dTe_accel(0,ij) = dTe_accel(0,ij) + tmp2*we_accel(m)
      dTe_adiab(0,ij) = dTe_adiab(0,ij) + tmp2*we_adiab(m)
      dTe_fdamp(0,ij) = dTe_fdamp(0,ij) + tmp2*we_fdamp(m)
      tmp1 = wgf1*wa0*inner_wgq0
      dne_init(1,ij)  = dne_init(1,ij)  + tmp1*we_init(m)
      dne_vpara(1,ij) = dne_vpara(1,ij) + tmp1*we_vpara(m)
      dne_vperp(1,ij) = dne_vperp(1,ij) + tmp1*we_vperp(m)
      dne_accel(1,ij) = dne_accel(1,ij) + tmp1*we_accel(m)
      dne_adiab(1,ij) = dne_adiab(1,ij) + tmp1*we_adiab(m)
      dne_fdamp(1,ij) = dne_fdamp(1,ij) + tmp1*we_fdamp(m)
      tmp2 = tmp1*((2.0/3.0)*energy-backTe(m)*msm0)
      dTe_init(1,ij)  = dTe_init(1,ij)  + tmp2*we_init(m)
      dTe_vpara(1,ij) = dTe_vpara(1,ij) + tmp2*we_vpara(m)
      dTe_vperp(1,ij) = dTe_vperp(1,ij) + tmp2*we_vperp(m)
      dTe_accel(1,ij) = dTe_accel(1,ij) + tmp2*we_accel(m)
      dTe_adiab(1,ij) = dTe_adiab(1,ij) + tmp2*we_adiab(m)
      dTe_fdamp(1,ij) = dTe_fdamp(1,ij) + tmp2*we_fdamp(m)

      !inner radius, larger theta index
      ij=ij+1
      tmp1 = wgf0*wa0*inner_wgq1
      dne_init(0,ij)  = dne_init(0,ij)  + tmp1*we_init(m)
      dne_vpara(0,ij) = dne_vpara(0,ij) + tmp1*we_vpara(m)
      dne_vperp(0,ij) = dne_vperp(0,ij) + tmp1*we_vperp(m)
      dne_accel(0,ij) = dne_accel(0,ij) + tmp1*we_accel(m)
      dne_adiab(0,ij) = dne_adiab(0,ij) + tmp1*we_adiab(m)
      dne_fdamp(0,ij) = dne_fdamp(0,ij) + tmp1*we_fdamp(m)
      tmp2 = tmp1*((2.0/3.0)*energy-backTe(m)*msm0)
      dTe_init(0,ij)  = dTe_init(0,ij)  + tmp2*we_init(m)
      dTe_vpara(0,ij) = dTe_vpara(0,ij) + tmp2*we_vpara(m)
      dTe_vperp(0,ij) = dTe_vperp(0,ij) + tmp2*we_vperp(m)
      dTe_accel(0,ij) = dTe_accel(0,ij) + tmp2*we_accel(m)
      dTe_adiab(0,ij) = dTe_adiab(0,ij) + tmp2*we_adiab(m)
      dTe_fdamp(0,ij) = dTe_fdamp(0,ij) + tmp2*we_fdamp(m)
      tmp1 = wgf1*wa0*inner_wgq1
      dne_init(1,ij)  = dne_init(1,ij)  + tmp1*we_init(m)
      dne_vpara(1,ij) = dne_vpara(1,ij) + tmp1*we_vpara(m)
      dne_vperp(1,ij) = dne_vperp(1,ij) + tmp1*we_vperp(m)
      dne_accel(1,ij) = dne_accel(1,ij) + tmp1*we_accel(m)
      dne_adiab(1,ij) = dne_adiab(1,ij) + tmp1*we_adiab(m)
      dne_fdamp(1,ij) = dne_fdamp(1,ij) + tmp1*we_fdamp(m)
      tmp2 = tmp1*((2.0/3.0)*energy-backTe(m)*msm0)
      dTe_init(1,ij)  = dTe_init(1,ij)  + tmp2*we_init(m)
      dTe_vpara(1,ij) = dTe_vpara(1,ij) + tmp2*we_vpara(m)
      dTe_vperp(1,ij) = dTe_vperp(1,ij) + tmp2*we_vperp(m)
      dTe_accel(1,ij) = dTe_accel(1,ij) + tmp2*we_accel(m)
      dTe_adiab(1,ij) = dTe_adiab(1,ij) + tmp2*we_adiab(m)
      dTe_fdamp(1,ij) = dTe_fdamp(1,ij) + tmp2*we_fdamp(m)

      !outter radius, smaller theta index
      ij=jteon1(m)
      tmp1 = wgf0*wa1*outter_wgq0
      dne_init(0,ij)  = dne_init(0,ij)  + tmp1*we_init(m)
      dne_vpara(0,ij) = dne_vpara(0,ij) + tmp1*we_vpara(m)
      dne_vperp(0,ij) = dne_vperp(0,ij) + tmp1*we_vperp(m)
      dne_accel(0,ij) = dne_accel(0,ij) + tmp1*we_accel(m)
      dne_adiab(0,ij) = dne_adiab(0,ij) + tmp1*we_adiab(m)
      dne_fdamp(0,ij) = dne_fdamp(0,ij) + tmp1*we_fdamp(m)
      tmp2 = tmp1*((2.0/3.0)*energy-backTe(m)*msm0)
      dTe_init(0,ij)  = dTe_init(0,ij)  + tmp2*we_init(m)
      dTe_vpara(0,ij) = dTe_vpara(0,ij) + tmp2*we_vpara(m)
      dTe_vperp(0,ij) = dTe_vperp(0,ij) + tmp2*we_vperp(m)
      dTe_accel(0,ij) = dTe_accel(0,ij) + tmp2*we_accel(m)
      dTe_adiab(0,ij) = dTe_adiab(0,ij) + tmp2*we_adiab(m)
      dTe_fdamp(0,ij) = dTe_fdamp(0,ij) + tmp2*we_fdamp(m)
      tmp1 = wgf1*wa1*outter_wgq0
      dne_init(1,ij)  = dne_init(1,ij)  + tmp1*we_init(m)
      dne_vpara(1,ij) = dne_vpara(1,ij) + tmp1*we_vpara(m)
      dne_vperp(1,ij) = dne_vperp(1,ij) + tmp1*we_vperp(m)
      dne_accel(1,ij) = dne_accel(1,ij) + tmp1*we_accel(m)
      dne_adiab(1,ij) = dne_adiab(1,ij) + tmp1*we_adiab(m)
      dne_fdamp(1,ij) = dne_fdamp(1,ij) + tmp1*we_fdamp(m)
      tmp2 = tmp1*((2.0/3.0)*energy-backTe(m)*msm0)
      dTe_init(1,ij)  = dTe_init(1,ij)  + tmp2*we_init(m)
      dTe_vpara(1,ij) = dTe_vpara(1,ij) + tmp2*we_vpara(m)
      dTe_vperp(1,ij) = dTe_vperp(1,ij) + tmp2*we_vperp(m)
      dTe_accel(1,ij) = dTe_accel(1,ij) + tmp2*we_accel(m)
      dTe_adiab(1,ij) = dTe_adiab(1,ij) + tmp2*we_adiab(m)
      dTe_fdamp(1,ij) = dTe_fdamp(1,ij) + tmp2*we_fdamp(m)

      !outter radius, larger theta index
      ij=ij+1
      tmp1 = wgf0*wa1*outter_wgq1
      dne_init(0,ij)  = dne_init(0,ij)  + tmp1*we_init(m)
      dne_vpara(0,ij) = dne_vpara(0,ij) + tmp1*we_vpara(m)
      dne_vperp(0,ij) = dne_vperp(0,ij) + tmp1*we_vperp(m)
      dne_accel(0,ij) = dne_accel(0,ij) + tmp1*we_accel(m)
      dne_adiab(0,ij) = dne_adiab(0,ij) + tmp1*we_adiab(m)
      dne_fdamp(0,ij) = dne_fdamp(0,ij) + tmp1*we_fdamp(m)
      tmp2 = tmp1*((2.0/3.0)*energy-backTe(m)*msm0)
      dTe_init(0,ij)  = dTe_init(0,ij)  + tmp2*we_init(m)
      dTe_vpara(0,ij) = dTe_vpara(0,ij) + tmp2*we_vpara(m)
      dTe_vperp(0,ij) = dTe_vperp(0,ij) + tmp2*we_vperp(m)
      dTe_accel(0,ij) = dTe_accel(0,ij) + tmp2*we_accel(m)
      dTe_adiab(0,ij) = dTe_adiab(0,ij) + tmp2*we_adiab(m)
      dTe_fdamp(0,ij) = dTe_fdamp(0,ij) + tmp2*we_fdamp(m)
      tmp1 = wgf1*wa1*outter_wgq1
      dne_init(1,ij)  = dne_init(1,ij)  + tmp1*we_init(m)
      dne_vpara(1,ij) = dne_vpara(1,ij) + tmp1*we_vpara(m)
      dne_vperp(1,ij) = dne_vperp(1,ij) + tmp1*we_vperp(m)
      dne_accel(1,ij) = dne_accel(1,ij) + tmp1*we_accel(m)
      dne_adiab(1,ij) = dne_adiab(1,ij) + tmp1*we_adiab(m)
      dne_fdamp(1,ij) = dne_fdamp(1,ij) + tmp1*we_fdamp(m)
      tmp2 = tmp1*((2.0/3.0)*energy-backTe(m)*msm0)
      dTe_init(1,ij)  = dTe_init(1,ij)  + tmp2*we_init(m)
      dTe_vpara(1,ij) = dTe_vpara(1,ij) + tmp2*we_vpara(m)
      dTe_vperp(1,ij) = dTe_vperp(1,ij) + tmp2*we_vperp(m)
      dTe_accel(1,ij) = dTe_accel(1,ij) + tmp2*we_accel(m)
      dTe_adiab(1,ij) = dTe_adiab(1,ij) + tmp2*we_adiab(m)
      dTe_fdamp(1,ij) = dTe_fdamp(1,ij) + tmp2*we_fdamp(m)

   enddo


   if(npartdom>1) then
      call MPI_ALLREDUCE(dne_init, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dne_init=recv_array
      call MPI_ALLREDUCE(dne_vpara, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dne_vpara=recv_array
      call MPI_ALLREDUCE(dne_vperp, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dne_vperp=recv_array
      call MPI_ALLREDUCE(dne_accel, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dne_accel=recv_array
      call MPI_ALLREDUCE(dne_adiab, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dne_adiab=recv_array
      call MPI_ALLREDUCE(dne_fdamp, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dne_fdamp=recv_array
      call MPI_ALLREDUCE(dTe_init, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dTe_init=recv_array
      call MPI_ALLREDUCE(dTe_vpara, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dTe_vpara=recv_array
      call MPI_ALLREDUCE(dTe_vperp, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dTe_vperp=recv_array
      call MPI_ALLREDUCE(dTe_accel, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dTe_accel=recv_array
      call MPI_ALLREDUCE(dTe_adiab, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dTe_adiab=recv_array
      call MPI_ALLREDUCE(dTe_fdamp, recv_array, 2*mgrid, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      dTe_fdamp=recv_array
   endif

   ! poloidal end cell; need valid data for densityi(:,igrid(i))
   do i=0,mpsi

      dne_init(:,igrid(i)+mtheta(i))=dne_init(:,igrid(i)+mtheta(i))+dne_init(:,igrid(i))
      dne_init(:,igrid(i))=dne_init(:,igrid(i)+mtheta(i))
      dne_vpara(:,igrid(i)+mtheta(i))=dne_vpara(:,igrid(i)+mtheta(i))+dne_vpara(:,igrid(i))
      dne_vpara(:,igrid(i))=dne_vpara(:,igrid(i)+mtheta(i))
      dne_vperp(:,igrid(i)+mtheta(i))=dne_vperp(:,igrid(i)+mtheta(i))+dne_vperp(:,igrid(i))
      dne_vperp(:,igrid(i))=dne_vperp(:,igrid(i)+mtheta(i))
      dne_accel(:,igrid(i)+mtheta(i))=dne_accel(:,igrid(i)+mtheta(i))+dne_accel(:,igrid(i))
      dne_accel(:,igrid(i))=dne_accel(:,igrid(i)+mtheta(i))
      dne_adiab(:,igrid(i)+mtheta(i))=dne_adiab(:,igrid(i)+mtheta(i))+dne_adiab(:,igrid(i))
      dne_adiab(:,igrid(i))=dne_adiab(:,igrid(i)+mtheta(i))
      dne_fdamp(:,igrid(i)+mtheta(i))=dne_fdamp(:,igrid(i)+mtheta(i))+dne_fdamp(:,igrid(i))
      dne_fdamp(:,igrid(i))=dne_fdamp(:,igrid(i)+mtheta(i))

      dTe_init(:,igrid(i)+mtheta(i))=dTe_init(:,igrid(i)+mtheta(i))+dTe_init(:,igrid(i))
      dTe_init(:,igrid(i))=dTe_init(:,igrid(i)+mtheta(i))
      dTe_vpara(:,igrid(i)+mtheta(i))=dTe_vpara(:,igrid(i)+mtheta(i))+dTe_vpara(:,igrid(i))
      dTe_vpara(:,igrid(i))=dTe_vpara(:,igrid(i)+mtheta(i))
      dTe_vperp(:,igrid(i)+mtheta(i))=dTe_vperp(:,igrid(i)+mtheta(i))+dTe_vperp(:,igrid(i))
      dTe_vperp(:,igrid(i))=dTe_vperp(:,igrid(i)+mtheta(i))
      dTe_accel(:,igrid(i)+mtheta(i))=dTe_accel(:,igrid(i)+mtheta(i))+dTe_accel(:,igrid(i))
      dTe_accel(:,igrid(i))=dTe_accel(:,igrid(i)+mtheta(i))
      dTe_adiab(:,igrid(i)+mtheta(i))=dTe_adiab(:,igrid(i)+mtheta(i))+dTe_adiab(:,igrid(i))
      dTe_adiab(:,igrid(i))=dTe_adiab(:,igrid(i)+mtheta(i))
      dTe_fdamp(:,igrid(i)+mtheta(i))=dTe_fdamp(:,igrid(i)+mtheta(i))+dTe_fdamp(:,igrid(i))
      dTe_fdamp(:,igrid(i))=dTe_fdamp(:,igrid(i)+mtheta(i))
   enddo


   ! send densityi to left and receive from right
   call MPI_SENDRECV(dne_init(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dne_init(1,ii+1:ii+jt)=dne_init(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dne_init(1,:)=dne_init(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dne_vpara(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dne_vpara(1,ii+1:ii+jt)=dne_vpara(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dne_vpara(1,:)=dne_vpara(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dne_vperp(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dne_vperp(1,ii+1:ii+jt)=dne_vperp(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dne_vperp(1,:)=dne_vperp(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dne_accel(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dne_accel(1,ii+1:ii+jt)=dne_accel(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dne_accel(1,:)=dne_accel(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dne_adiab(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dne_adiab(1,ii+1:ii+jt)=dne_adiab(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dne_adiab(1,:)=dne_adiab(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dne_fdamp(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dne_fdamp(1,ii+1:ii+jt)=dne_fdamp(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dne_fdamp(1,:)=dne_fdamp(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dTe_init(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dTe_init(1,ii+1:ii+jt)=dTe_init(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dTe_init(1,:)=dTe_init(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dTe_vpara(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dTe_vpara(1,ii+1:ii+jt)=dTe_vpara(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dTe_vpara(1,:)=dTe_vpara(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dTe_vperp(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dTe_vperp(1,ii+1:ii+jt)=dTe_vperp(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dTe_vperp(1,:)=dTe_vperp(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dTe_accel(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dTe_accel(1,ii+1:ii+jt)=dTe_accel(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dTe_accel(1,:)=dTe_accel(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dTe_adiab(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dTe_adiab(1,ii+1:ii+jt)=dTe_adiab(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dTe_adiab(1,:)=dTe_adiab(1,:)+recv_array(1,:)
   endif

   call MPI_SENDRECV(dTe_fdamp(0,:),mgrid,mpi_Rsize,left_pe,1,recv_array(1,:),mgrid,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      do i=0,mpsi
         ii=igrid(i); jt=mtheta(i)
         dTe_fdamp(1,ii+1:ii+jt)=dTe_fdamp(1,ii+1:ii+jt)+ cshift(recv_array(1,ii+1:ii+jt),itran(i))
      enddo
   else
      dTe_fdamp(1,:)=dTe_fdamp(1,:)+recv_array(1,:)
   endif


   !===============================================================================================
   ! Caution. This routine does not enforce charge conservation regardless of its_jacobi parameter
   !===============================================================================================
   do i=0,mpsi
      do j=1,mtheta(i)
         ij=igrid(i)+j
         dne_init(1,ij) =dne_init(1,ij)*markere(ij)
         dne_vpara(1,ij)=dne_vpara(1,ij)*markere(ij)
         dne_vperp(1,ij)=dne_vperp(1,ij)*markere(ij)
         dne_accel(1,ij)=dne_accel(1,ij)*markere(ij)
         dne_adiab(1,ij)=dne_adiab(1,ij)*markere(ij)
         dne_fdamp(1,ij)=dne_fdamp(1,ij)*markere(ij)

         dTe_init(1,ij) =dTe_init(1,ij)*markere(ij)
         dTe_vpara(1,ij)=dTe_vpara(1,ij)*markere(ij)
         dTe_vperp(1,ij)=dTe_vperp(1,ij)*markere(ij)
         dTe_accel(1,ij)=dTe_accel(1,ij)*markere(ij)
         dTe_adiab(1,ij)=dTe_adiab(1,ij)*markere(ij)
         dTe_fdamp(1,ij)=dTe_fdamp(1,ij)*markere(ij)
      enddo
      dne_init(1,igrid(i))=dne_init(1,igrid(i)+mtheta(i))
      dne_vpara(1,igrid(i))=dne_vpara(1,igrid(i)+mtheta(i))
      dne_vperp(1,igrid(i))=dne_vperp(1,igrid(i)+mtheta(i))
      dne_accel(1,igrid(i))=dne_accel(1,igrid(i)+mtheta(i))
      dne_adiab(1,igrid(i))=dne_adiab(1,igrid(i)+mtheta(i))
      dne_fdamp(1,igrid(i))=dne_fdamp(1,igrid(i)+mtheta(i))

      dTe_init(1,igrid(i))=dTe_init(1,igrid(i)+mtheta(i))
      dTe_vpara(1,igrid(i))=dTe_vpara(1,igrid(i)+mtheta(i))
      dTe_vperp(1,igrid(i))=dTe_vperp(1,igrid(i)+mtheta(i))
      dTe_accel(1,igrid(i))=dTe_accel(1,igrid(i)+mtheta(i))
      dTe_adiab(1,igrid(i))=dTe_adiab(1,igrid(i)+mtheta(i))
      dTe_fdamp(1,igrid(i))=dTe_fdamp(1,igrid(i)+mtheta(i))
   enddo


   !==============================================================================
   ! MGYOO (01/16/2019)
   ! Get real dne by multiplying ne0
   do i=0,mpsi
      do j=0,mtheta(i)
         ij = igrid(i) + j
         dne_init(1,ij)  = den_grid(i)*dne_init(1,ij)
         dne_vpara(1,ij) = den_grid(i)*dne_vpara(1,ij)
         dne_vperp(1,ij) = den_grid(i)*dne_vperp(1,ij)
         dne_accel(1,ij) = den_grid(i)*dne_accel(1,ij)
         dne_adiab(1,ij) = den_grid(i)*dne_adiab(1,ij)
         dne_fdamp(1,ij) = den_grid(i)*dne_fdamp(1,ij)
      enddo
   enddo

   ! ! Convert dTe to get a correct form of dTe by considering dne.
   ! ! buff = Sigma(weight*(2/3*energy-T0))*markere = (1/n0)*Sigma(weight*(2/3*energy-T0))/node_Vol
   ! ! dTe = (1/(n0+dn))*Sigma(weight*(2/3*energy-T0))/node_Vol = n0/(n0+dn)*buff
   ! do i=0,mpsi
   !   do j=0,mtheta(i)
   !     ij = igrid(i) + j
   !     tmp1 =  den_grid(i)/(den_grid(i)+tot_dne(1,ij))
   !     dTe_init(1,ij)  = tmp1*dTe_init(1,ij)
   !     dTe_vpara(1,ij) = tmp1*dTe_vpara(1,ij)
   !     dTe_vperp(1,ij) = tmp1*dTe_vperp(1,ij)
   !     dTe_accel(1,ij) = tmp1*dTe_accel(1,ij)
   !     dTe_adiab(1,ij) = tmp1*dTe_adiab(1,ij)
   !     dTe_fdamp(1,ij) = tmp1*dTe_fdamp(1,ij)
   !   enddo
   ! enddo
   ! !==============================================================================


end subroutine Cal_dne_dTe_from_each_electron_weight_component
