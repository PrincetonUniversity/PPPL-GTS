subroutine solve_dphidt( solver, masssolver )
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use field_array
  use diagnosis_array
  use esi_interface
  use solver_data_class
  implicit none
  !
  type(solver_data) :: solver 
  type(solver_data) :: masssolver 
  !
  integer i,j,k,ii,jj,ij,ip,jt,ismooth,it,n, &
       icount,ierr,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
  real(wp) sendl(mgrid),sendr(mgrid),recvl(mgrid),recvr(mgrid),phism(mgrid),&
       pleft(mthetamax),pright(mthetamax),phitmp(0:1,mgrid)
 
  real(wp) tmp,rms
  
  real(wp) adum(0:mpsi), bdum(0:mpsi), tbdum(0:mpsi), rms_b, tdum_i, tdum_e
  real(wp) dni(0:mpsi),dne(0:mpsi),tdni(0:mpsi),tdne(0:mpsi),sne,sni,tsne,tsni
  real(wp) bb_rhs(mgrid),xx(mgrid) 

  real(wp), external :: SMOOTHSTEP
  real(wp) instep,outstep

  real(wp) fstep1, fstep2, fstep3
  fstep1=1.00_wp
  fstep2=1.0_wp
  fstep3=1.0_wp

  do jj=1,3

!$omp parallel do private(ij)
     do ij=1,mgrid
!          phitmp(1,ij)=(fstep3*pflux_e(jj,1,ij)-pflux_i(jj,1,ij))*Jcbn(ij) ! =(Gamma_e-Gamma_i)*J
        phitmp(1,ij)=(fstep3*pflux_e(jj,1,ij)-pflux_i(jj,1,ij))          ! =(Gamma_e-Gamma_i)*J
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
!         radial smoothing I
!
!            i=1
!  !$omp parallel do private(j,ij)
!          do j=1,mtheta(i)
!             ij=igrid(i)+j
!             phism(ij)=0.25_wp*((1.0_wp-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
!                  wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1)+&
!                  (1.0_wp-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
!                  wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))
!          enddo
!  
!          i=mpsi-1
!  !$omp parallel do private(j,ij)
!          do j=1,mtheta(i)
!             ij=igrid(i)+j
!             phism(ij)=0.25_wp*((1.0-wtp1(1,ij))*phitmp(1,jtp1(1,ij))+&
!                  wtp1(1,ij)*phitmp(1,jtp1(1,ij)+1)+&
!                  (1.0_wp-wtp1(2,ij))*phitmp(1,jtp1(2,ij))+&
!                  wtp1(2,ij)*phitmp(1,jtp1(2,ij)+1))
!          enddo
!  
!  !$omp parallel do private(i)
!          do i=igrid(2)+1, igrid(mpsi-2)+mtheta(mpsi-2)
!             phitmp(1,i)=0.625_wp*phitmp(1,i)+phism(i) ! 5-point smoothing
!          end do
!  
!  !$omp parallel do private(i)
!          do i=igrid(1)+1,igrid(1)+mtheta(1)
!             phitmp(1,i)=0.5_wp*phitmp(1,i)+phism(i)   ! 3-point smoothing
!          end do
!  
!  !$omp parallel do private(i)
!          do i=igrid(mpsi-1)+1,igrid(mpsi-1)+mtheta(mpsi-1)
!             phitmp(1,i)=0.5_wp*phitmp(1,i)+phism(i)   ! 3-point smoothing
!          end do

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
           phitmp(:,igrid(i))=phitmp(:,igrid(i)+mtheta(i))
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
! send to right and receive from left
        sendr=phitmp(1,:)
        recvl=0.0_wp
        icount=mgrid
        idest=right_pe
        isource=left_pe
        isendtag=myrank_toroidal
        irecvtag=isource
        call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,recvl,icount,&
                      mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierr)
     
! send to left and receive from right
        sendl=phitmp(1,:)
        recvr=0.0_wp
        idest=left_pe
        isource=right_pe
        isendtag=myrank_toroidal
        irecvtag=isource
        call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,recvr,icount,&
                      mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierr)
    
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

! *************** validate phitmp(0,*) and phitmp(*,igrid(i)) ******
! *************** may not needed ***********************************

! toroidal BC: send to right and receive from left
     sendr=phitmp(1,:)
     recvl=0.0_wp
     icount=mgrid
     idest=right_pe
     isource=left_pe
     isendtag=myrank_toroidal
     irecvtag=isource
     call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,recvl,icount,&
                   mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierr)
     
     if(myrank_toroidal==0)then
        do i=0,mpsi
           ii=igrid(i)
           jt=mtheta(i)
           phitmp(0,ii+1:ii+jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
        enddo
     else
        phitmp(0,:)=recvl
     endif

! ******************************************************************

! poloidal BC
     do i=0,mpsi
        phitmp(:,igrid(i))=phitmp(:,igrid(i)+mtheta(i))
     end do

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
        pflux_e(jj,:,i)=phitmp(:,i)    ! pflux_e=(Gamma_e-Gamma_i)*J/n_0; 
                                       ! now pflux_e(:,0,:) also has valid data;
                                       ! pflux_e(:,:,igird(i)) also has valid data
     enddo

  enddo   ! end do jj=1,3

  call divergence_j  ! --> pflux_i=dvi[(Gamma_e-Gamma_i)*J]/(n_0*J)  -- May, 2014

!$omp parallel do private(i,j,ij,sne,tsne,tsni)
  do i=0, mpsi
     sne=(phip00(i)+Er_grid(i))/Te_grid(i)*csc0
     do j=1,mtheta(i)
        ij=igrid(i)+j

!         pflux_e(1,1,ij)=fstep1*(pflux_i(1,1,ij)+pflux_i(2,1,ij)+pflux_i(3,1,ij))+   &
!                         fstep2*(cij(1,ij)*(evector(1,1,ij)-phip00(i))+           &
!                         cij(2,ij)*evector(2,1,ij)+cij(3,ij)*evector(3,1,ij))   
!                         ! --> This is the right hand side of dphi/dt equation


!  replaced by (2012): 
        tsne=sne+cij(4,ij)
        tsni=cij(3,ij)*cij(5,ij)*(evector(1,1,ij)-phip00(i))- &
             cij(3,ij)*tsne*evector(2,1,ij)+ &
             (cij(2,ij)*tsne-cij(1,ij)*cij(5,ij))*evector(3,1,ij)

        pflux_e(1,1,ij)=fstep1*(pflux_i(1,1,ij)+pflux_i(2,1,ij)+pflux_i(3,1,ij))+   &
                        fstep2*tsni
                   ! --> This is the right hand side of normalized dphi/dt equation which is [dphi/dt-tilde(dphi/dt)]/Ti = ... 
                   ! --> also the RHS of normalized dphi/dt equation of PDE form which is: 
                   !     beta*nabla_perp[sum(gs)*nabla_perp(dphi/dt)] = ... (see note)

!        pflux_e(1,1,ij)=pflux_e(1,1,ij)*vth_grid(i)*vth_grid(i)  ! multiply by Ti

        if(solver%order <= 0) then 
           pflux_e(1,1,ij)=pflux_e(1,1,ij)*tite(i)  ! multiply by Ti/Te --> RHS of following equation
                                                    ! Eq.--> [dphi/dt-tilde(dphi/dt)]/Te = ...
        endif

     enddo
  enddo

  do i=0, mpsi
     pflux_e(1,1,igrid(i))=pflux_e(1,1,igrid(i)+mtheta(i))  ! which is missing in last do-loop
  enddo

! ! initialize dphidt which is needed to estimate tilde(dphidt) using Jacobian interation method
!   if(istep==0)then
!      if(irun==0) then 
!         dphidt=0.0_wp    ! for follow-up runs, read dphidt from restart
!         dphidt00=0.0_wp  ! for follow-up runs, read dphidt00 from restart
!      endif
!   endif


! Solve  

if(solver%order <= 0) then   ! use integral solver -- never worked!!
! Note: dphidt=(dphi/dt)/Te
  dphidt(1,:)=0.0_wp  ! addition to RHS with dphidt=0 for 1st iteration

  do it=1,its
     bb_rhs=0.0_wp
     xx=0.0_wp

     do i=0,mpsi
        outstep=smoothstep(rmin,rmax,psi_grid(i))
       !       instep=smoothstep(lmin,lmax,psi_grid(i))
       !       instep=instep*(1.0d0-outstep)
       !       outstep=1.0d0-instep
        ! tmp = (0.21 + 0.09*outstep)     ! full electron (no split) --> tmp=0.01+0.09=0.1  for r > rmax  -- see solver_petsc.f90
        tmp = (gama0 + gama1*outstep)     ! full electron (no split) --> tmp=0.01+0.09=0.1  for r > rmax  -- see solver_petsc.f90
        do j=1,mtheta(i)
           ij=igrid(i)+j
           bb_rhs(ij)=pflux_e(1,1,ij)+tmp*dphidt(1,ij)
        end do
        bb_rhs(igrid(i))=bb_rhs(igrid(i)+mtheta(i))
     end do

     call solve( solver, masssolver, mgrid, bb_rhs, xx, ierr )
     dphidt(1,:) = xx       ! =(dphi/dt)/Te

     ! radial boundary: dphidt=0 at both inner and outer boundaries
     dphidt(1,igrid(0):igrid(0)+mtheta(0))=0.0_wp
     ! dphidt(1,igrid(mpsi):igrid(mpsi)+mtheta(mpsi))=0.0_wp

     do i=0, mpsi
        dphidt(1,igrid(i))=dphidt(1,igrid(i)+mtheta(i))  !validate dphidt(1,igrid(i)) for safe
     enddo

     if(mype==0)   print*, 'dphidt solver', 'its=', it, dphidt(1,1), dphidt(1,mgrid/3), dphidt(1,mgrid) 

  end do

else     ! using FEM solver

  bb_rhs=0.0_wp
  xx=0.0_wp
  bb_rhs(:)=pflux_e(1,1,:)

!!  bb_rhs=bb_rhs*1.0d-30

!!  if(mype==1) print*, 'HHHHEEEEEERRRRR'

  call solve( solver, masssolver, mgrid, bb_rhs, xx, ierr )
  dphidt(1,:) = xx       ! =dphi/dt

  do i=0, mpsi
     dphidt(1,igrid(i))=dphidt(1,igrid(i)+mtheta(i))  !validate dphidt(1,igrid(i)) for safe
  enddo

endif 







  adum=0.0_wp
  rms=0.0_wp
  bdum=0.0_wp
  rms_b=0.0_wp
  dphidt00=0.0_wp

!$omp parallel private(i,j,ij) reduction(+:adum,rms,bdum,rms_b)
!$omp do
  do i=0,mpsi
     do j=1,mtheta(i)
        ij=igrid(i)+j
        adum(i)=adum(i)+dtemper(ij)*dphidt(1,ij)  ! surface averaging
        rms=rms+dphidt(1,ij)**2
        bdum(i)=bdum(i)+dtemper(ij)*bb_rhs(ij)  ! surface averaging
        rms_b=rms_b+bb_rhs(ij)**2
     enddo
  enddo
!$omp end parallel

  call MPI_ALLREDUCE(adum,dphidt00,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierr)
  dphidt00=dphidt00/jaco_sum    !  =<dphi/dt>

  call MPI_ALLREDUCE(bdum,tbdum,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierr)
  tbdum=tbdum/jaco_sum    !  =<dphi/dt>

  if(mype==0) print*, 'petsc:  <dphi/dt>(',ipsi,')=',  dphidt00(ipsi), 'rms=', rms

! at this point, dphidt=d(phi)/dt; d(delta_phi)/dt=dphidt-dphidt00
!$omp parallel do private(i,j,ij)
  do i=0, mpsi
     do j=1,mtheta(i)
        ij=igrid(i)+j
        dphidt(1,ij)=dphidt(1,ij)-dphidt00(i)  ! --> [d(delta_phi)/dt]/Te when order=0
                                               ! --> d(delta_phi)/dt when order>0
     enddo
  enddo

! validate dphidt(0,:)
  
! toroidal BC: send to right and receive from left
  sendr=dphidt(1,:)
  recvl=0.0_wp
  icount=mgrid
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
       recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierr)

  if(myrank_toroidal==0) then
     do i=0,mpsi
        ii=igrid(i)
        jt=mtheta(i)
        dphidt(0,ii+1:ii+jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
     enddo
  else
     dphidt(0,:)=recvl
  endif

! poloidal BC
  do i=0,mpsi
     dphidt(:,igrid(i))=dphidt(:,igrid(i)+mtheta(i))
  end do

! output: dphidt=d(delta_phi)/dt for order>0; both dphidt(0,:) and dphidt(:,igrid(i)) have valid data
! NOTE: with the modified RHS, dphidt represents [d(delta_phi)/dt]/Te (normalized) when order=0

end subroutine solve_dphidt





subroutine divergence_j
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use field_array
  use diagnosis_array
  use esi_interface
  implicit none

  integer i,ii,ij,j,k,icount,idest,isource,isendtag,irecvtag,ierr,jt,&
          istatus(MPI_STATUS_SIZE)
  real(wp) difft(0:mpsi),diffz,r,pleft(mthetamax),pright(mthetamax),&
       sendl(mgrid),recvr(mgrid),sendr(mgrid),recvl(mgrid),delq

  real(wp) dpsi_inv,e_psi2,e_psi1,dpsi2_inv,dpsi1,dpsi2,dpsi1_inv     !wx

!
! Note: This subroutine calculates divergence-j at grid points
!       using finite difference:  
!       INOUT:  pflux_e(1:3, 0:1, mgrid)=(Gamma_e-Gamma_i)*J
!       OUTPUT: pflux_i(1:3, 0:1, mgrid)=div[(Gamma_e-Gamma_i)*J]/(n_0*J)
!       pflux_i(1,*,*)=dpflux_e(1,*,*)/da/(n_0*J); 
!       pflux_i(2,*,*)=dpflux_e(2,*,*)/dtheta/(n_0*J);
!       pflux_i(3,*,*)=dpflux_e(3,*,*)/dzeta/(n_0*J);
!       not involving pflux_e(:,0,:) and pflux_i(:,0,:).
!       We set pflux_i(:,:,:)=0 at the radial boundries, which may be not true, 
!       but pflux_i at the radial boundries are not used elsewhere. 
!

  difft=0.5_wp/deltat
  diffz=0.5_wp/deltaz

  pflux_i=0.0_wp

! d_pflux_e(1,*,*)/d_a
!$omp parallel do private(i,j,ij,dpsi1,dpsi2,dpsi1_inv,dpsi2_inv,dpsi_inv,e_psi2,e_psi1)
  do i=1,mpsi-1
     dpsi1=psi_grid(i)-psi_grid(i-1)
     dpsi2=psi_grid(i+1)-psi_grid(i)
     dpsi1_inv=1.0_wp/dpsi1
     dpsi2_inv=1.0_wp/dpsi2
     dpsi_inv=1.0_wp/(dpsi1+dpsi2)
     do j=1,mtheta(i)
        ij=igrid(i)+j

        e_psi2=((1.0-wtp1(1,ij))*pflux_e(1,1,jtp1(1,ij))+&
              wtp1(1,ij)*pflux_e(1,1,jtp1(1,ij)+1)-pflux_e(1,1,ij))*dpsi2_inv

        e_psi1=(pflux_e(1,1,ij)-((1.0_wp-wtp1(2,ij))*pflux_e(1,1,jtp1(2,ij))+&
              wtp1(2,ij)*pflux_e(1,1,jtp1(2,ij)+1)))*dpsi1_inv

!           pflux_i(1,1,ij)=(e_psi1+e_psi2)*0.5_wp
        pflux_i(1,1,ij)=(e_psi1*dpsi2+e_psi2*dpsi1)*dpsi_inv  ! note: pflux_i(1,*,*)=d/da[(Gamma_e-Gamma_i)^a*J/n_0]
     enddo
  enddo
! at radial boundaries i=0(a0) and mpsi(a1) pflux_i(1,*,*)=0. This is not true! 

! d_pflux_e(2,*,*)/d_theta and d_pflux_e(3,*,*)/d_theta
!$omp parallel do private(i,j,ij,jt)
  do i=1,mpsi-1
     do j=1,mtheta(i)
        ij=igrid(i)+j
        jt=mtheta(i)*(j/mtheta(i))
        pflux_i(2,1,ij)=difft(i)*(pflux_e(2,1,ij+1-jt)-pflux_e(2,1,ij-1))

        pflux_i(2,0,ij)=difft(i)*(pflux_e(3,1,ij+1-jt)-pflux_e(3,1,ij-1))
     enddo
  enddo
  
! send pflux_e(3,:,:) to right and receive from left
  sendr=pflux_e(3,1,:)
  recvl=0.0_wp
  icount=mgrid
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
       recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierr)
  
! send pflux_e(3,:,:) to left and receive from right
  sendl=pflux_e(3,1,:)
  recvr=0.0_wp
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,&
       recvr,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierr)
! unpack pflux_e(3,:,:)_boundary and calculate d_pflux_e(3,:,:)/d_zeta at boundaries, mzeta=1
!$omp parallel do private(i,j,ii,jt,ij,pleft,pright)
  do i=1,mpsi-1
     ii=igrid(i)
     jt=mtheta(i)
     if(myrank_toroidal==0)then !down-shift for zeta=0
        pleft(1:jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
        pright(1:jt)=recvr(ii+1:ii+jt)
     elseif(myrank_toroidal==numberpe-1)then !up-shift for zeta=2*pi
        pright(1:jt)=cshift(recvr(ii+1:ii+jt),itran(i))
        pleft(1:jt)=recvl(ii+1:ii+jt)
     else
        pleft(1:jt)=recvl(ii+1:ii+jt)
        pright(1:jt)=recvr(ii+1:ii+jt)
     endif

! d_pflux_e(3,:,:)/d_zeta
     do j=1,mtheta(i)
        ij=igrid(i)+j
        pflux_i(3,1,ij)=(pright(j)-pleft(j))*diffz
     enddo
  enddo

! adjust the difference between safety factor q and qtinv for fieldline coordinate
!omp parallel do private(i,j,delq,ij)
  do i=1,mpsi-1
     delq=qtinv(i)

     do j=1,mtheta(i)
        ij=igrid(i)+j
        !   pflux_i(3,1,ij)=pflux_i(3,1,ij)-delq*pflux_i(2,1,ij) !    WRONG!!!     pflux_i(3...=dpflux_e(3,:,:)/dzeta| with a,theta=const.
        pflux_i(3,1,ij)=pflux_i(3,1,ij)-delq*pflux_i(2,0,ij) ! pflux_i(3...=dpflux_e(3,:,:)/dzeta| with a,theta=const.
     enddo
  enddo

!  !omp parallel do private(ij)
!    do ij=1,mgrid
!       pflux_i(1,1,ij)=pflux_i(1,1,ij)/Jcbn(ij)
!       pflux_i(2,1,ij)=pflux_i(2,1,ij)/Jcbn(ij)
!       pflux_i(3,1,ij)=pflux_i(3,1,ij)/Jcbn(ij)
!    enddo

!$omp parallel do private(i,j,ij)
  do i=0,mpsi
     do j=1,mtheta(i)
        ij=igrid(i)+j
        pflux_i(1,1,ij)=pflux_i(1,1,ij)/(Jcbn(ij)*den_grid(i))  ! (1/j*n_0)*d/da[J*(Gamma_e-Gamma_i)^a]
        pflux_i(2,1,ij)=pflux_i(2,1,ij)/(Jcbn(ij)*den_grid(i))  ! (1/j*n_0)*d/theta[J*(Gamma_e-Gamma_i)^theta]
        pflux_i(3,1,ij)=pflux_i(3,1,ij)/(Jcbn(ij)*den_grid(i))  ! (1/j*n_0)*d/dzeta[J*(Gamma_e-Gamma_i)^zeta]
     enddo
  enddo
  
! poloidal end point
!omp parallel do private(i)
  do i=1,mpsi-1
     pflux_i(1,1,igrid(i))=pflux_i(1,1,igrid(i)+mtheta(i))
     pflux_i(2,1,igrid(i))=pflux_i(2,1,igrid(i)+mtheta(i))
     pflux_i(3,1,igrid(i))=pflux_i(3,1,igrid(i)+mtheta(i))
  enddo
! now pflux_i(:,:,igrid(i)) have valid data 
!     pflux_i(:,k=0,:) has NO valid data
! pflux_i(1:3,*,*)=0 at radial boundaries

!     calculate pflux_i at radial boundaries i=0(a0) and mpsi(a1) here if nessary in the future

end subroutine divergence_j

