subroutine field(masssolver)
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use field_array
  use diagnosis_array
  use esi_interface

  use solver_data_class

  implicit none

  type(solver_data) :: masssolver

  integer i,ii,ij,j,k,icount,idest,isource,isendtag,irecvtag,ierror,ip,jt,&
       istatus(MPI_STATUS_SIZE)
  real(wp) difft(0:mpsi),diffz,r,pleft(mthetamax),pright(mthetamax),&
       sendl(mgrid),recvr(mgrid),sendr(mgrid),recvl(mgrid),sendrs(3,mgrid),&
       recvls(3,mgrid),q,delq,adum(0:mpsi)

  real(wp) dpsi_inv,e_psi2,e_psi1,dpsi2_inv,dpsi1,dpsi2,dpsi1_inv,dum,tdum  !wx

  real(wp) instep,outstep,fdamp,fdamp0,fip,dtime
  real(wp), external :: SMOOTHSTEP 




! Note: evector,phip00,phi00 etc are calculated at particle coordinates (x), but not guiding center coordinates (R)

!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(order <= 0.and.ZF_SEPARATION==0) then

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

! !$omp parallel do private(ij)
!      do ij=1,mgrid
!         phi(:,ij)=phi(:,ij)-tdum ! Re-level potential phi by a constant to make <phi>(a0)=0 at inner boundary
!      enddo

     phi=phi-tdum  ! Re-level potential phi by a constant to make <phi>(a0)=0 at inner boundary

     do i=0, mpsi
        phi(:,igrid(i))=phi(:,igrid(i)+mtheta(i))  !validate phi(1,igrid(i)) for safe
     enddo

! now phi(0,:) also has valid data;  phi(:,igird(i)) also has valid data

  end if
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! finite difference for e-field in equilibrium unit
  difft=0.5_wp/deltat
  diffz=0.5_wp/deltaz
!$omp parallel do private(i,j,k)
  do i=1,mgrid
     do j=0,1
        do k=1,3
           evector(k,j,i)=0.0_wp
        enddo
     enddo
  enddo



if(order <= 0) then   ! using direct numerical derivative to claculate d_phi/d_a and d_phi/d_theta when using integral Poisson solver

! d_phi/d_a,    d_phi/d_psi   
!$omp parallel do private(i,j,ij,dpsi1,dpsi2,dpsi1_inv,dpsi2_inv,dpsi_inv,e_psi2,e_psi1)
  do i=1,mpsi-1
     dpsi1=psi_grid(i)-psi_grid(i-1)
     dpsi2=psi_grid(i+1)-psi_grid(i)
     dpsi1_inv=1.0_wp/dpsi1
     dpsi2_inv=1.0_wp/dpsi2
     dpsi_inv=1.0_wp/(dpsi1+dpsi2)
     do j=1,mtheta(i)
        ij=igrid(i)+j

        e_psi2=((1.0_wp-wtp1(1,ij))*phi(1,jtp1(1,ij))+&
               wtp1(1,ij)*phi(1,jtp1(1,ij)+1)-phi(1,ij))*dpsi2_inv

        e_psi1=(phi(1,ij)-((1.0_wp-wtp1(2,ij))*phi(1,jtp1(2,ij))+&
               wtp1(2,ij)*phi(1,jtp1(2,ij)+1)))*dpsi1_inv

!           evector(1,1,ij)=(e_psi1+e_psi2)*0.5
        evector(1,1,ij)=(e_psi1*dpsi2+e_psi2*dpsi1)*dpsi_inv  ! note: evector(1,*,*)=d\delta-phi/da, ZF_SEPARATION=1
                                                              !       evector(1,*,*)=dphi/da,        ZF_SEPARATION=0
     enddo
  enddo
! at radial boundaries i=0(a0) and mpsi(a1) evector(1,*,*)=0 

! d_phi/d_theta
!$omp parallel do private(i,j,ij,jt)
  do i=1,mpsi-1
     do j=1,mtheta(i)
        ij=igrid(i)+j
        jt=mtheta(i)*(j/mtheta(i))
        evector(2,1,ij)=difft(i)*(phi(1,ij+1-jt)-phi(1,ij-1))
     enddo
  enddo

else                  ! using FEM to d_phi/d_a and d_phi/d_theta when using FEM solver
  
  sendl(:)=phi(1,:)
  ierror=0

  !  call find_dxdadxdt( masssolver, mgrid, phi(1,:), evector(1,1,:), evector(2,1,:), ierror )
  call find_dxdadxdt( masssolver, mgrid, sendl, recvr, recvl, ierror )

!  if(mype==0) then
!     if(istep==1.and.irk==1) then 
!     do ij=1,mgrid
!        write(16, 206) ij,evector(1,1,ij),evector(2,1,ij),recvr(ij),recvl(ij)
!     end do
!     end if
!     call flush_(16)
!  206 format(2x, i6, 4(2x,ES13.6))
!  end if
   

  evector(1,1,:)=recvr(:)
  evector(2,1,:)=recvl(:)

end if
  
! send phi to right and receive from left
  sendr=phi(1,:)
  recvl=0.0_wp
  icount=mgrid
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
       recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
  
! send phi to left and receive from right
  sendl=phi(1,:)
  recvr=0.0_wp
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,&
       recvr,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
! unpack phi_boundary and calculate E_zeta at boundaries, mzeta=1
!$omp parallel do private(i,j,ii,jt,ij,pleft,pright)
  do i=1,mpsi-1
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

! d_phi/d_zeta at a=const. and theta-zeta/q_N=const. 
     do j=1,mtheta(i)
        ij=igrid(i)+j
        evector(3,1,ij)=(pright(j)-pleft(j))*diffz
     enddo
  enddo

! adjust the difference between safety factor q and qtinv for fieldline coordinate
!omp parallel do private(i,j,delq,ij)
  do i=1,mpsi-1
     delq=qtinv(i)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        evector(3,1,ij)=evector(3,1,ij)-delq*evector(2,1,ij) ! --> now evector(3...=d_phi/d_zeta
                     ! if delq=(1.0/q_grid(i)-qtinv(i)), evector(3...=b.grad_phi/b.grad_zeta
     enddo
  enddo
  
  if(ZF_SEPARATION==0) then   ! calculate zonal flow component <phi> and phip00=d<phi>/da

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

! direct numerical derivative to get d<phi>/da
!$omp parallel do private(i,dpsi1,dpsi2,dpsi1_inv,dpsi2_inv,dpsi_inv,e_psi2,e_psi1)
    do i=1,mpsi-1
       dpsi1=psi_grid(i)-psi_grid(i-1)
       dpsi2=psi_grid(i+1)-psi_grid(i)
       dpsi1_inv=1.0_wp/dpsi1
       dpsi2_inv=1.0_wp/dpsi2
       dpsi_inv=1.0_wp/(dpsi1+dpsi2)
       e_psi2=(phi00(i+1)-phi00(i))*dpsi2_inv
       e_psi1=(phi00(i)-phi00(i-1))*dpsi1_inv
       phip00(i)=(e_psi1*dpsi2+e_psi2*dpsi1)*dpsi_inv
    enddo
! boundary conditions d<phi>/da=0 (not phi=0 at both boundaries)
    phip00(0)=0.0_wp
    !  phip00(mpsi)=0.0_wp
    if(ifneoclassical==0) then
       phip00(mpsi)=0.0_wp
    else        ! d<phi>/da /= 0 when NC physics included
       phip00(mpsi)=(phi00(mpsi)-phi00(mpsi-1))/(psi_grid(mpsi)-psi_grid(mpsi-1))
       if(order <= 0) evector(1,1,igrid(mpsi):(igrid(mpsi)+mtheta(mpsi)))=phip00(mpsi)   ! also set dphi/dtheta=dphi/dzeta=0 at mpsi
    endif

  end if     ! end if(ZF_SEPARATION==0)

! add (0,0) mode, d phi/da
  if(mode00==1)then
    if(ZF_SEPARATION==1) then   ! add ZF component when ZF is solved separately
!omp parallel do private(i,j,ii,jt)
     do i=1,mpsi-1

! zero gyroradius limit OR linear     ????
!        phip00(i)=0.0             !wx    no zonal flow

        ii=igrid(i)
        jt=mtheta(i)
        do j=ii+1,ii+jt
           evector(1,1,j)=evector(1,1,j)+phip00(i)   ! --> dphi/da (total)
        enddo
     enddo
    endif
  else                          ! mode00==0 (exclude ZF case)
    if(ZF_SEPARATION==0) then   ! substract ZF component when total potential is solved
!omp parallel do private(i,j,ii,jt)
     do i=1,mpsi-1
        ii=igrid(i)
        jt=mtheta(i)
        do j=ii+1,ii+jt
           evector(1,1,j)=evector(1,1,j)-phip00(i)   ! --> d\delta-phi/da (total)
        enddo
     enddo
    endif  ! end if(ZF_SEPARATION==0)
  endif    ! end if(mode00==1)
  
! toroidal end point, pack end point (k=mzeta) data 
! send E to right and receive from left
!$omp parallel do private(i)
  do i=1,mgrid
     sendrs(1,i)=evector(1,1,i)
     sendrs(2,i)=evector(2,1,i)
     sendrs(3,i)=evector(3,1,i)
     recvls(1,i)=0.0
     recvls(2,i)=0.0
     recvls(3,i)=0.0
  enddo
  icount=3*mgrid
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendrs,icount,mpi_Rsize,idest,isendtag,&
       recvls,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

! unpack end point data for k=0
  if(myrank_toroidal==0)then
!$omp parallel do private(i,ii,jt)
     do i=1,mpsi-1
        ii=igrid(i)
        jt=mtheta(i)
        evector(1,0,ii+1:ii+jt)=cshift(recvls(1,ii+1:ii+jt),-itran(i))
        evector(2,0,ii+1:ii+jt)=cshift(recvls(2,ii+1:ii+jt),-itran(i))
        evector(3,0,ii+1:ii+jt)=cshift(recvls(3,ii+1:ii+jt),-itran(i))
     enddo
  else
!$omp parallel do private(i)
     do i=1,mgrid
        evector(1,0,i)=recvls(1,i)
        evector(2,0,i)=recvls(2,i)
        evector(3,0,i)=recvls(3,i)
     enddo
  endif
  
! poloidal end point
!omp parallel do private(i,j)
  do i=1,mpsi-1
     evector(1,0,igrid(i))=evector(1,0,igrid(i)+mtheta(i))
     evector(2,0,igrid(i))=evector(2,0,igrid(i)+mtheta(i))
     evector(3,0,igrid(i))=evector(3,0,igrid(i)+mtheta(i))

     evector(1,1,igrid(i))=evector(1,1,igrid(i)+mtheta(i))
     evector(2,1,igrid(i))=evector(2,1,igrid(i)+mtheta(i))
     evector(3,1,igrid(i))=evector(3,1,igrid(i)+mtheta(i))
  enddo
! now evector(:,0,:) and evector(:,:,igrid(i)) have valid data 
! evector(1:3,*,*)=0 at radial boundaries






!  ! ******************************* added on may 5, 09
!  !  fdamp0=1.0_wp/tstep*5.0d-4
!    fdamp0=4.0_wp/tstep
!    dtime=float(irk)/2.0_wp*tstep
!    do i=1,mpsi-1
!       instep =smoothstep(lmin,lmax,psi_grid(i))
!       outstep=smoothstep(rmin,rmax,psi_grid(i))
!       fip=instep*(1.0_wp-outstep)
!       fdamp=(1.0_wp-fip)*fdamp0
!       ii=igrid(i)
!       jt=mtheta(i)
!       evector(:,:,ii:ii+jt)=evector(:,:,ii:ii+jt)*exp(-dtime*fdamp)
!       phip00(i)=phip00(i)*exp(-dtime*fdamp)
!    end do
!  ! **********************************




     

!$omp parallel do private(i,j,k)
!   do i=igrid(0),igrid(0)+mtheta(0)
!      do j=0,1
!         do k=1,3
!            if(evector(k,j,i).ne.0.0_wp)write(0,*)'E is not zero at inner boundary', mype, k, j, i, evector(k,j,i)
!         enddo
!      enddo
!   enddo

!$omp parallel do private(i,j,k)
!   do i=igrid(mpsi),igrid(mpsi)+mtheta(mpsi)
!      do j=0,1
!         do k=1,3
!            if(evector(k,j,i).ne.0.0_wp)write(0,*)'E is not zero at outer boundary', mype, k, j, i, evector(k,j,i)
!         enddo
!      enddo
!   enddo
  
end subroutine field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

