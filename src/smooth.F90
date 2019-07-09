subroutine smooth(iflag)

   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use field_array
   use diagnosis_array
   use esi_interface
   implicit none
   integer iflag,i,j,k,ii,ij,ip,jt,jte,kz,ismooth,mz,mzmax,mzbig,jpe,indp,indt,indp1,indt1,&
      meachtheta,icount,ierror,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
   real(wp) sendl(mgrid),sendr(mgrid),recvl(mgrid),recvr(mgrid),phism(mgrid),&
      pleft(mthetamax),pright(mthetamax),phitmp(0:1,mgrid),filter(mtdiag/2+1),wt,r,&
      dt,wz,zdum,tdum,pi2_inv,den00(0:mpsi),phiflux(mtdiag/ntoroidal,mtdiag,idiag1:idiag2),&
      eachzeta((idiag2-idiag1+1)*mtdiag*mtdiag/ntoroidal/ntoroidal),&
      allzeta((idiag2-idiag1+1)*mtdiag*mtdiag/ntoroidal),xz(mtdiag)
   complex(wp) y_eigen(mtdiag/ntoroidal*num_mode),yz(mtdiag/2+1),yt(mtdiag*num_mode),&
      ye(mtdiag),yp(ipsi+1)
   ! the dummy array size for thread-safe FFT assume data array size <16384
   real(wp) :: aux1f(25000),aux1b(25000),aux2f(20000),aux2b(20000)
   real(wp) :: aux3f(1),aux3b(1)
   real(wp) :: tmp_zonal_term(0:mpsi)




   ! wx: solve zonal flow for general geometry (July, 04):
   !     assume <tilde{phi}>=tilde{<phi>} -- probabaly not justified for gg

   if(iflag==1)then
      if(ZF_SEPARATION==1.and.order <= 0) then    ! solve zonal flows separately
         phip00=zonali                           ! =<n_i>/n_0
         if(ifelectron/=0) phip00=phip00-zonale  ! =(<n_i>-<n_e>)/n_0 with non-adiabatic electrons
         ! (-0.0625 0.25 0.625 0.25 -0.0625) radial smoothing of (0,0) mode density phip00
         do ismooth=1,1
            den00(0)=phip00(0)
            den00(mpsi)=phip00(mpsi)
            den00(1)=0.5_wp*phip00(1)+0.25_wp*(phip00(0)+phip00(2))
            den00(mpsi-1)=0.5_wp*phip00(mpsi-1)+0.25_wp*(phip00(mpsi-2)+phip00(mpsi))
            den00(2:mpsi-2)=phip00(0:mpsi-4)+phip00(4:mpsi)
            den00(2:mpsi-2)=0.625_wp*phip00(2:mpsi-2)+0.25_wp*(phip00(1:mpsi-3)+&
               phip00(3:mpsi-1))-0.0625_wp*den00(2:mpsi-2)
            phip00=den00
         enddo

         ! solve zonal flow for gg: phip00=dphi_00/da, dphi_00/da(a0)=0. Trapezoid rule
         den00=phip00
         phip00=0.0_wp
         do i=1,mpsi
            phip00(i)=phip00(i-1)+0.5_wp*(psi_grid(i)-psi_grid(i-1))*  &
               (dVda(i-1)*B2_s(i-1)*den00(i-1)+dVda(i)*B2_s(i)*den00(i))
         enddo

         do i=0,mpsi
            phip00(i)=-phip00(i)/dVda(i)/gaa(i)
         enddo

         ! add FLR contribution using Pade approximation: b*<phi>=(1+b)*<n>
         phi00=den00*vth_grid**2   ! =T*(dn_00/n_0)
         do i=1,mpsi-1
            phip00(i)=phip00(i) + (phi00(i+1)-phi00(i-1))/(psi_grid(i+1)-psi_grid(i-1))  ! =dphi_00/da
         enddo
         phip00(mpsi)=phip00(mpsi) ! taking (d/da)(Ti*S)=0 at right boundary

         ! enforce boundary condition
         phip00(0)=0.0_wp
         !     phip00(mpsi)=0.0   ! only one boundary condition is needed to solve dphi_00/da

         ! zero gyroradius limit OR linear simulation
         !     phip00=0.0     !wx      no zonal flow

         ! (0,0) mode potential store in phi00
         phi00=0.0_wp   ! phi00(0)=0.0 is chosen
         do i=1,mpsi
            phi00(i)=phi00(i-1)+0.5_wp*(psi_grid(i)-psi_grid(i-1))*(phip00(i-1)+phip00(i))
         enddo

         ! cancel following do loop           !wx unreasonable enforcement
         ! enforce <phi>=0 at radial boundary
         !     do i=1,mpsi
         !        phi00(i)=phi00(i)-phi00(mpsi)*real(i)/real(mpsi)
         !     enddo

         if(mode00==0) phip00=0.0_wp  ! exclude zonal flow
      endif    ! end ZF_SEPARATION==1
   endif     ! end iflag=1






   if(iflag==1) then
      !$omp parallel do private(i)
      do i=1,mgrid
         phitmp(1,i)=phi(1,i)  !phitmp=phi  ! smooth phi for numerical derivative -->E field
      enddo
   else                                    !   if(iflag==0), smooth n_i or n_e just for output; used nowhere
      if(ifelectron==0) then               !if(ifelectron==0), phitmp=densityi
         !$omp parallel do private(i)
         ! MGYOO (06/11/2019): For consistent smoothing including zonal components
         if(ZF_SEPARATION==1) then
            do i=0,mpsi
               do j=0,mtheta(i)
                  ij=igrid(i)+j
                  phitmp(1,ij)=densityi(1,ij)+zonali(i)  ! tilde{dn_i}/n_0
               enddo
            enddo
         else
            do i=1,mgrid
               phitmp(1,i)=densityi(1,i)  !just for output; used nowhere
            enddo
         endif
      else                                 !if(ifelectron==1 or 2), phitmp=densitye
         ! MGYOO (06/11/2019): For consistent smoothing including zonal components
         if(ZF_SEPARATION==1) then
            do i=0,mpsi
               do j=0,mtheta(i)
                  ij=igrid(i)+j
                  phitmp(1,ij)=densitye(1,ij)+zonale(i)  ! tilde{dn_e}/n_0
               enddo
            enddo
         else
            !$omp parallel do private(i)
            do i=1,mgrid
               phitmp(1,i)=densitye(1,i)  !just for output; used nowhere
            enddo
         endif
      end if
   endif   ! end if(iflag==1)

   ismooth=1
   do ip=1,ismooth

      ! radial smoothing
      do i=1,mpsi-1
         phitmp(:,igrid(i))=phitmp(:,igrid(i)+mtheta(i))
      enddo
      phism=0.0
      !$omp parallel do private(i,j,ij)
      do i=2,mpsi-2
         do j=1,mtheta(i)
            ij=igrid(i)+j
            !
            ! need to take into account uneven radial grid size dpsi when doing radial smoothing? wx
            !
            phism(ij)=0.25_wp*((1.0_wp-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
               wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1)+&
               (1.0_wp-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
               wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))-&
               0.0625_wp*((1.0_wp-wtp2(1,ij))*phitmp(1,jtp2(1,ij))+&
               wtp2(1,ij)*phitmp(1,jtp2(1,ij)+1)+&
               (1.0_wp-wtp2(2,ij))*phitmp(1,jtp2(2,ij))+&
               wtp2(2,ij)*phitmp(1,jtp2(2,ij)+1))

         enddo
      enddo

      !
      !       radial smoothing I
      !
      !       i=1
      !  !$omp parallel do private(j,ij)
      !          do j=1,mtheta(i)
      !             ij=igrid(i)+j
      !             phism(ij)=0.25*((1.0-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
      !                  wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1)+&
      !                  (1.0-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
      !                  wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))
      !          enddo
      !
      !       i=mpsi-1
      !  !$omp parallel do private(j,ij)
      !          do j=1,mtheta(i)
      !             ij=igrid(i)+j
      !             phism(ij)=0.25*((1.0-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
      !                  wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1)+&
      !                  (1.0-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
      !                  wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))
      !          enddo
      !
      !  !$omp parallel do private(i)
      !       do i=igrid(2)+1, igrid(mpsi-2)+mtheta(mpsi-2)
      !          phitmp(1,i)=0.625*phitmp(1,i)+phism(i) ! 5-point smoothing
      !       end do
      !
      !  !$omp parallel do private(i)
      !       do i=igrid(1)+1,igrid(1)+mtheta(1)
      !          phitmp(1,i)=0.5*phitmp(1,i)+phism(i)   ! 3-point smoothing
      !       end do
      !
      !  !$omp parallel do private(i)
      !       do i=igrid(mpsi-1)+1,igrid(mpsi-1)+mtheta(mpsi-1)
      !          phitmp(1,i)=0.5*phitmp(1,i)+phism(i)   ! 3-point smoothing
      !       end do

      !
      !       radial smoothing II
      !
      i=1
      !$omp parallel do private(j,ij)
      do j=1,mtheta(i)
         ij=igrid(i)+j
         phism(ij)=0.25_wp*((1.0_wp-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
            wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1)+&
            (1.0_wp-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
            wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))-&
            0.0625_wp*((1.0_wp-wtp2(1,ij))*phitmp(1,jtp2(1,ij))+&
            wtp2(1,ij)*phitmp(1,jtp2(1,ij)+1)+&
            phitmp(1,ij))
      enddo

      i=mpsi-1
      !$omp parallel do private(j,ij)
      do j=1,mtheta(i)
         ij=igrid(i)+j
         phism(ij)=0.25_wp*((1.0_wp-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
            wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1)+&
            (1.0_wp-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
            wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))-&
            0.0625_wp*(phitmp(1,ij)+&
            (1.0_wp-wtp2(2,ij))*phitmp(1,jtp2(2,ij))+&
            wtp2(2,ij)*phitmp(1,jtp2(2,ij)+1))
      enddo

      i=0
      !$omp parallel do private(j,ij)
      do j=1,mtheta(i)
         ij=igrid(i)+j
         phism(ij)=0.5_wp*((1.0_wp-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
            wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1))-&
            0.125_wp*((1.0-wtp2(1,ij))*phitmp(1,jtp2(1,ij))+&
            wtp2(1,ij)*phitmp(1,jtp2(1,ij)+1))
      enddo

      i=mpsi
      !$omp parallel do private(j,ij)
      do j=1,mtheta(i)
         ij=igrid(i)+j
         phism(ij)=0.5_wp*((1.0_wp-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
            wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))-&
            0.125_wp*((1.0_wp-wtp2(2,ij))*phitmp(1,jtp2(2,ij))+&
            wtp2(2,ij)*phitmp(1,jtp2(2,ij)+1))
      enddo

      !$omp parallel do private(i)
      do i=1,mgrid
         phitmp(1,i)=0.625_wp*phitmp(1,i)+phism(i) ! 5-point smoothing
      end do


      do i=0,mpsi
         phitmp(1,igrid(i))=phitmp(1,igrid(i)+mtheta(i))
      enddo

      ! poloidal smoothing (-0.0625 0.25 0.625 0.25 -0.0625)
      !$omp parallel do private(i,ii,jt,pright)
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         pright(1:jt)=phitmp(1,ii+1:ii+jt)
         phitmp(1,ii+1:ii+jt)=0.625_wp*pright(1:jt)+&
            0.25_wp*(cshift(pright(1:jt),-1)+cshift(pright(1:jt),1))-&
            0.0625_wp*(cshift(pright(1:jt),-2)+cshift(pright(1:jt),2))
      enddo

      do i=0,mpsi
         phitmp(:,igrid(i))=phitmp(:,igrid(i)+mtheta(i))
      enddo

      ! parallel smoothing
      ! send phi to right and receive from left
      sendr=phitmp(1,:)
      recvl=0.0_wp
      icount=mgrid
      idest=right_pe
      isource=left_pe
      isendtag=myrank_toroidal
      irecvtag=isource
      call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
         recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

      ! send phi to left and receive from right
      sendl=phitmp(1,:)
      recvr=0.0_wp
      idest=left_pe
      isource=right_pe
      isendtag=myrank_toroidal
      irecvtag=isource
      call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,&
         recvr,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

      !$omp parallel do private(i,ii,jt,j,ij,pleft,pright)
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         if(myrank_toroidal==0)then !down-shift for zeta=0
            pleft(1:jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
            pright(1:jt)=recvr(ii+1:ii+jt)
         elseif(myrank_toroidal==ntoroidal-1)then !up-shift for zeta=2*pi
            pright(1:jt)=cshift(recvr(ii+1:ii+jt),itran(i))
            pleft(1:jt)=recvl(ii+1:ii+jt)
         else
            pleft(1:jt)=recvl(ii+1:ii+jt)
            pright(1:jt)=recvr(ii+1:ii+jt)
         endif

         do j=1,mtheta(i)
            ij=igrid(i)+j
            phitmp(1,ij)=0.5_wp*phitmp(1,ij)+0.25_wp*(pleft(j)+pright(j))
         enddo
      enddo   ! end do i=0,mpsi
   enddo      ! end do ip=1,ismooth

   ! toroidal BC: send phi to right and receive from left
   sendr=phitmp(1,:)
   recvl=0.0_wp
   icount=mgrid
   idest=right_pe
   isource=left_pe
   isendtag=myrank_toroidal
   irecvtag=isource
   call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
      recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

   if(myrank_toroidal==0)then
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         phitmp(0,ii+1:ii+jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
      enddo
   else
      phitmp(0,:)=recvl
   endif

   ! poloidal BC
   do i=0,mpsi
      phitmp(:,igrid(i))=phitmp(:,igrid(i)+mtheta(i))
   enddo

   ! radial boundary: 3 cases
   ! 1. phi=0 at both inner and outer boundaries
   ! 2. phi=0 at inner boundary and float phi at outer boundary
   ! 3. phi is float at both inner and outer boundaries
   !    for above three cases, electric field E=grad\phi=0 is enforced at both boundaries

   !    phitmp(:,igrid(0):igrid(0)+mtheta(0))=0.0
   !    potential at right boundary is float
   !    phitmp(:,igrid(mpsi):igrid(mpsi)+mtheta(mpsi))=0.0

   !adiabatic electron response
   if(iflag==1) then
      !$omp parallel do private(i)
      do i=1,mgrid
         phi(:,i)=phitmp(:,i)    ! now phi(0,:) also has valid data
         ! phi(:,igird(i)) also has valid data
      enddo
   else  ! if(iflag==0)

      if(ifelectron==0) then
         !$omp parallel do private(i)
         do i=1,mgrid
            densityi(:,i)=phitmp(:,i) !just for output; used nowhere
         enddo

         ! MGYOO (06/11/2019): For consistent smoothing including zonal components
         tmp_zonal_term=0.0
         do i=0,mpsi
            do j=1,mtheta(i)
               ij=igrid(i)+j
               tmp_zonal_term(i)=tmp_zonal_term(i)+dtemper(ij)*densityi(1,ij)
            enddo
         enddo
         call MPI_ALLREDUCE(tmp_zonal_term,zonali,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
         zonali=zonali/jaco_sum

         if(ZF_SEPARATION==1) then
            ! Subtract zonal term to get correct non-zonal density
            do i=0,mpsi
               do j=1,mtheta(i)
                  ij=igrid(i)+j
                  densityi(1,ij)=densityi(1,ij)-zonali(i)
               enddo
               densityi(1,igrid(i))=densityi(1,igrid(i)+mtheta(i))
            enddo
         endif
      else    ! if(ifelectron==1 or 2)
         !$omp parallel do private(i)
         do i=1,mgrid
            densitye(:,i)=phitmp(:,i) !just for output; used nowhere
         enddo

         ! MGYOO (06/11/2019): For consistent smoothing including zonal components
         tmp_zonal_term=0.0
         do i=0,mpsi
            do j=1,mtheta(i)
               ij=igrid(i)+j
               tmp_zonal_term(i)=tmp_zonal_term(i)+dtemper(ij)*densitye(1,ij)
            enddo
         enddo
         call MPI_ALLREDUCE(tmp_zonal_term,zonale,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
         zonale=zonale/jaco_sum

         if(ZF_SEPARATION==1) then
            ! Subtract zonal term to get correct non-zonal density
            do i=0,mpsi
               do j=1,mtheta(i)
                  ij=igrid(i)+j
                  densitye(1,ij)=densitye(1,ij)-zonale(i)
               enddo
               densitye(1,igrid(i))=densitye(1,igrid(i)+mtheta(i))
            enddo
         endif
      end if
   endif




   ! interpolate on a flux surface from fieldline coordinates to magnetic coordinates
   ! use mtdiag for both poloidal and toroidal grid points
   !    if(iflag==1 .and. (nonlinear<0.5 .or. idiag==0))then


#ifdef USE_OLD_DIAGNOSTICS
   go to 1000

   if(mype==0) print*, 'INSIDE SMOOTH(1) 111'

   ! diagnostic
   if(idiag==0 .and. iflag==1)then

      if(mype==0) print*, 'INSIDE SMOOTH(1) 222', 'idiag=', idiag, 'iflag=', iflag

      ! shear flow amplitude
      eradial=sqrt(sum(phip00(1:mpsi)**2)/real(mpsi,wp))/vthc    ! ???

      if(mype==0) print*, 'INSIDE SMOOTH(1) 333', 'idiag=', idiag, 'iflag=', iflag

      ! RMS of potential
      efield=0.0_wp
      do i=1,mpsi
         efield=efield+sum(phi(1,igrid(i)+1:igrid(i)+mtheta(i))**2)
      enddo
      efield=sqrt(efield/real(sum(mtheta(1:mpsi)),wp))/(vthc**2)  !  ???

      if(mype==0) print*, 'INSIDE SMOOTH(1) 6aaa7'

      ! record zonal flow mode history
      call fftr1d(1,mpsi,1.0_wp,phip00(1:mpsi),yp,1)

      if(mype==0) print*, 'INSIDE SMOOTH(1) 6bbb7'

      do i=1,num_mode
         amp_mode(1,i,1)=real(yp(i+1))/(real(mpsi,wp)*vthc)  ! ???
         amp_mode(2,i,1)=aimag(yp(i+1))/(real(mpsi,wp)*vthc) ! ???
      enddo

      if(mype==0) print*, 'INSIDE SMOOTH(1) 6ccc7'

      ! dominant (n,m) mode history data for flux surface at i=ipsi
      icount=meachtheta*num_mode
      call MPI_GATHER(y_eigen,icount,mpi_complex_size,yt,icount,&
         mpi_complex_size,0,toroidal_comm,ierror)

      if(mype==0) print*, 'INSIDE SMOOTH(1) 6ddd7'

      if(myrank_toroidal == 0)then
         do kz=1,num_mode

            do i=0,ntoroidal-1
               do j=1,meachtheta
                  ye(j+i*meachtheta)=yt(j+(kz-1)*meachtheta+i*icount)
               enddo
            enddo

            if(mype==0) print*, 'INSIDE SMOOTH(1) do-loop 1'

            call fftc1d(1,mtdiag,1.0_wp,ye)

            if(mype==0) print*, 'INSIDE SMOOTH(1) do-loop 2'

            amp_mode(1,kz,2)=real(ye(mtdiag-mmode(kz)+1))/ &
               (real(mzmax*mtdiag,wp)*vth_grid(ipsi)**2)
            amp_mode(2,kz,2)=aimag(ye(mtdiag-mmode(kz)+1))/ &
               (real(mzmax*mtdiag,wp)*vth_grid(ipsi)**2)

            if(mype==0) print*, 'INSIDE SMOOTH(1) do-loop 3'

         enddo
      endif
   endif   ! end if(idiag==0 .and. iflag==1)

   if(mype==0) print*, 'INSIDE SMOOTH(1) 777 end'

   1000 continue
#endif

end subroutine smooth







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine smooth_phi

   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use field_array
   use diagnosis_array
   implicit none

   integer i,j,k,ii,ij,ip,jt,ismooth, &
      icount,ierror,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
   real(wp) sendl(mgrid),sendr(mgrid),recvl(mgrid),recvr(mgrid),phism(mgrid), &
      pleft(mthetamax),pright(mthetamax),phitmp(0:1,mgrid), &
      adum(0:mpsi)
   real(wp) dum, tdum

   !$omp parallel do private(i)
   do i=1,mgrid
      phitmp(1,i)=phi(1,i)  !phitmp=phi  ! smooth phi for numerical derivative -->E field
   enddo

   ismooth=1
   do ip=1,ismooth

      ! radial smoothing
      do i=1,mpsi-1
         phitmp(:,igrid(i))=phitmp(:,igrid(i)+mtheta(i))
      enddo
      phism=0.0
      !$omp parallel do private(i,j,ij)
      do i=2,mpsi-2
         do j=1,mtheta(i)
            ij=igrid(i)+j
            !
            ! need to take into account uneven radial grid size dpsi when doing radial smoothing? wx
            !
            phism(ij)=0.25_wp*((1.0_wp-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
               wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1)+&
               (1.0_wp-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
               wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))-&
               0.0625_wp*((1.0_wp-wtp2(1,ij))*phitmp(1,jtp2(1,ij))+&
               wtp2(1,ij)*phitmp(1,jtp2(1,ij)+1)+&
               (1.0_wp-wtp2(2,ij))*phitmp(1,jtp2(2,ij))+&
               wtp2(2,ij)*phitmp(1,jtp2(2,ij)+1))

         enddo
      enddo

      !
      !       radial smoothing I
      !
      !       i=1
      !  !$omp parallel do private(j,ij)
      !          do j=1,mtheta(i)
      !             ij=igrid(i)+j
      !             phism(ij)=0.25*((1.0-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
      !                  wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1)+&
      !                  (1.0-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
      !                  wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))
      !          enddo
      !
      !       i=mpsi-1
      !  !$omp parallel do private(j,ij)
      !          do j=1,mtheta(i)
      !             ij=igrid(i)+j
      !             phism(ij)=0.25*((1.0-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
      !                  wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1)+&
      !                  (1.0-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
      !                  wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))
      !          enddo
      !
      !  !$omp parallel do private(i)
      !       do i=igrid(2)+1, igrid(mpsi-2)+mtheta(mpsi-2)
      !          phitmp(1,i)=0.625*phitmp(1,i)+phism(i) ! 5-point smoothing
      !       end do
      !
      !  !$omp parallel do private(i)
      !       do i=igrid(1)+1,igrid(1)+mtheta(1)
      !          phitmp(1,i)=0.5*phitmp(1,i)+phism(i)   ! 3-point smoothing
      !       end do
      !
      !  !$omp parallel do private(i)
      !       do i=igrid(mpsi-1)+1,igrid(mpsi-1)+mtheta(mpsi-1)
      !          phitmp(1,i)=0.5*phitmp(1,i)+phism(i)   ! 3-point smoothing
      !       end do

      !
      !       radial smoothing II
      !
      i=1
      !$omp parallel do private(j,ij)
      do j=1,mtheta(i)
         ij=igrid(i)+j
         phism(ij)=0.25_wp*((1.0_wp-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
            wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1)+&
            (1.0_wp-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
            wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))-&
            0.0625_wp*((1.0_wp-wtp2(1,ij))*phitmp(1,jtp2(1,ij))+&
            wtp2(1,ij)*phitmp(1,jtp2(1,ij)+1)+&
            phitmp(1,ij))
      enddo

      i=mpsi-1
      !$omp parallel do private(j,ij)
      do j=1,mtheta(i)
         ij=igrid(i)+j
         phism(ij)=0.25_wp*((1.0_wp-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
            wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1)+&
            (1.0_wp-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
            wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))-&
            0.0625_wp*(phitmp(1,ij)+&
            (1.0_wp-wtp2(2,ij))*phitmp(1,jtp2(2,ij))+&
            wtp2(2,ij)*phitmp(1,jtp2(2,ij)+1))
      enddo

      i=0
      !$omp parallel do private(j,ij)
      do j=1,mtheta(i)
         ij=igrid(i)+j
         phism(ij)=0.5_wp*((1.0_wp-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
            wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1))-&
            0.125_wp*((1.0-wtp2(1,ij))*phitmp(1,jtp2(1,ij))+&
            wtp2(1,ij)*phitmp(1,jtp2(1,ij)+1))
      enddo

      i=mpsi
      !$omp parallel do private(j,ij)
      do j=1,mtheta(i)
         ij=igrid(i)+j
         phism(ij)=0.5_wp*((1.0_wp-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
            wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))-&
            0.125_wp*((1.0_wp-wtp2(2,ij))*phitmp(1,jtp2(2,ij))+&
            wtp2(2,ij)*phitmp(1,jtp2(2,ij)+1))
      enddo

      !$omp parallel do private(i)
      do i=1,mgrid
         phitmp(1,i)=0.625_wp*phitmp(1,i)+phism(i) ! 5-point smoothing
      end do


      do i=0,mpsi
         phitmp(1,igrid(i))=phitmp(1,igrid(i)+mtheta(i))
      enddo

      ! poloidal smoothing (-0.0625 0.25 0.625 0.25 -0.0625)
      !$omp parallel do private(i,ii,jt,pright)
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         pright(1:jt)=phitmp(1,ii+1:ii+jt)
         phitmp(1,ii+1:ii+jt)=0.625_wp*pright(1:jt)+&
            0.25_wp*(cshift(pright(1:jt),-1)+cshift(pright(1:jt),1))-&
            0.0625_wp*(cshift(pright(1:jt),-2)+cshift(pright(1:jt),2))
      enddo

      do i=0,mpsi
         phitmp(:,igrid(i))=phitmp(:,igrid(i)+mtheta(i))
      enddo

      ! parallel smoothing
      ! send phi to right and receive from left
      sendr=phitmp(1,:)
      recvl=0.0_wp
      icount=mgrid
      idest=right_pe
      isource=left_pe
      isendtag=myrank_toroidal
      irecvtag=isource
      call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
         recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

      ! send phi to left and receive from right
      sendl=phitmp(1,:)
      recvr=0.0_wp
      idest=left_pe
      isource=right_pe
      isendtag=myrank_toroidal
      irecvtag=isource
      call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,&
         recvr,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

      !$omp parallel do private(i,ii,jt,j,ij,pleft,pright)
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         if(myrank_toroidal==0)then !down-shift for zeta=0
            pleft(1:jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
            pright(1:jt)=recvr(ii+1:ii+jt)
         elseif(myrank_toroidal==ntoroidal-1)then !up-shift for zeta=2*pi
            pright(1:jt)=cshift(recvr(ii+1:ii+jt),itran(i))
            pleft(1:jt)=recvl(ii+1:ii+jt)
         else
            pleft(1:jt)=recvl(ii+1:ii+jt)
            pright(1:jt)=recvr(ii+1:ii+jt)
         endif

         do j=1,mtheta(i)
            ij=igrid(i)+j
            phitmp(1,ij)=0.5_wp*phitmp(1,ij)+0.25_wp*(pleft(j)+pright(j))
         enddo
      enddo   ! end do i=0,mpsi
   enddo      ! end do ip=1,ismooth

   ! toroidal BC: send phi to right and receive from left
   sendr=phitmp(1,:)
   recvl=0.0_wp
   icount=mgrid
   idest=right_pe
   isource=left_pe
   isendtag=myrank_toroidal
   irecvtag=isource
   call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
      recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

   if(myrank_toroidal==0)then
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         phitmp(0,ii+1:ii+jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
      enddo
   else
      phitmp(0,:)=recvl
   endif

   ! poloidal BC
   do i=0,mpsi
      phitmp(:,igrid(i))=phitmp(:,igrid(i)+mtheta(i))
   enddo

   ! radial boundary: 3 cases
   ! 1. phi=0 at both inner and outer boundaries
   ! 2. phi=0 at inner boundary and float phi at outer boundary
   ! 3. phi is float at both inner and outer boundaries
   !    for above three cases, electric field E=grad\phi=0 is enforced at both boundaries

   !    phitmp(:,igrid(0):igrid(0)+mtheta(0))=0.0
   !    potential at right boundary is float
   !    phitmp(:,igrid(mpsi):igrid(mpsi)+mtheta(mpsi))=0.0

   !$omp parallel do private(i)
   do i=1,mgrid
      phi(:,i)=phitmp(:,i)    ! now phi(0,:) also has valid data
      ! phi(:,igird(i)) also has valid data
   enddo




   if(ZF_SEPARATION==0) then  ! Calculate zonal flow component of phi, <phi>


      ! Note: phi and phi00 are calculated at particle coordinates (x), but not guiding center coordinates (R)

      ! There is arbitrary constant in the potential;
      ! Re-level potential phi to make <phi>(a0)=0 at inner boundary 12/07

      !   calculate <phi> at inner boundary a=a0
      dum=0.0_wp
      !$omp parallel private(j,ij) reduction(+:dum)
      !$omp do
      do j=1,mtheta(0)
         ij=igrid(0)+j
         dum=dum+dtemper(ij)*phi(1,ij)  ! surface averaging
      enddo
      !$omp end parallel
      call MPI_ALLREDUCE(dum,tdum,1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      tdum=tdum/jaco_sum(0)    !  =<phi> at inner boundary a=a0

      phi=phi-tdum  ! Re-level potential phi by a constant to make <phi>(a0)=0 at inner boundary

      do i=0, mpsi
         phi(:,igrid(i))=phi(:,igrid(i)+mtheta(i))  !validate phi(1,igrid(i)) for safe
      enddo

      ! now phi(0,:) also has valid data;  phi(:,igird(i)) also has valid data


      adum=0.0
      !$omp parallel private(i,j,ij) reduction(+:adum)
      !$omp do
      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid(i)+j
            adum(i)=adum(i)+dtemper(ij)*phi(1,ij)  ! surface averaging
         enddo
      enddo
      !$omp end parallel
      call MPI_ALLREDUCE(adum,phi00,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
      phi00=phi00/jaco_sum    !  =<phi>


   end if  ! end ZF_SEPARATION==0


end subroutine smooth_phi
