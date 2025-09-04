!============================================================
!  This file includes operators of field
!  The purpose of these operators is to implement electromagnetic formulae
!  zx

!============================================================
subroutine field_calc0(masssolver,phi,phi00,phip00,evector) !--to do...
!Purpose: 
!Input: 
!Output: 
!-------------------------------------------------------------
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use field_array,only:igrid,itran,qtinv,dtemper,psi_grid,wtp1,jtp1
  use diagnosis_array
  use esi_interface

  use solver_data_class

  implicit none

  real(wp) phi(0:1,mgrid),phi00(0:mpsi),phip00(0:mpsi),evector(1:3,0:1,mgrid)
  type(solver_data) :: masssolver

end subroutine field_calc0 
!============================================================
subroutine smooth_scalar0(phi,phi00) !--to do...
!Purpose: 
!Input: 
!Output: 
!-------------------------------------------------------------
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use field_array,only:igrid,itran,wtp1,wtp2,jtp1,jtp2,dtemper
  use diagnosis_array
  implicit none

  real(wp) phi(0:1,mgrid),phi00(0:mpsi)

end subroutine smooth_scalar0
!============================================================
subroutine fd_fetch_scalar(myphi,myrecv,mysend,idirect)
!Purpose: fetch myphi(1,:) from right or left; can be used for smooth and toroidal derivative
!Input: myphi(0:1,1:mgrid);
!       idirect=-1,1: direction of data flow, fetch myphi(1,:) from  right, left and put to myrecv
!Output: myrecv, with A. cshift operation and B. poloidal B.C. applied
!Workspace: mysend
!Convention: (1,:) is source, ghost (0,:) is destination 
!-------------------------------------------------------------
 use precision
 use mpi_defs
 use global_parameters
 use particle_decomp
 use field_array,only:itran,igrid
 
 implicit none
 integer idirect
 real(wp),dimension(0:1,mgrid)::myphi
 real(wp),dimension(mgrid)::myrecv,mysend

 integer icount,idest,isource,isendtag,irecvtag,ierror,istatus(MPI_STATUS_SIZE)
 integer i,ii,jt

 mysend(:)=myphi(1,:)
 myrecv=0.0_wp
 icount=mgrid

 select case(idirect)
  case(1) !fetch from left, myphi(1,:)-->myrecv
   idest=right_pe
   isource=left_pe
  case(-1) !fetch from right, myrecv<--myphi(1,:)
   idest=left_pe
   isource=right_pe
 end select 

 isendtag=myrank_toroidal
 irecvtag=isource
!--what's purpose of isendtag and irecvtag??
 call  MPI_SENDRECV(mysend,icount,mpi_Rsize,idest,isendtag,&
       myrecv,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

!--cshift of myrank=0 or ntoroidal-1
 do i=0,mpsi
  ii=igrid(i)
  jt=mtheta(i)
  if(idirect==1.and.myrank_toroidal==0)then !down-shift for zeta=0
    myrecv(ii+1:ii+jt)=cshift(myrecv(ii+1:ii+jt),-itran(i))
  elseif(idirect==-1.and.myrank_toroidal==ntoroidal-1)then !up-shift for zeta=2*pi
    myrecv(ii+1:ii+jt)=cshift(myrecv(ii+1:ii+jt),itran(i))
  else
  endif
 end do

!--poloidal B.C.
 do i=0,mpsi
  myrecv(igrid(i))=myrecv(igrid(i)+mtheta(i))
 end do
end subroutine fd_fetch_scalar
!============================================================
subroutine fd_fetch_vector(myphi,myrecv,mysend,ndim,idirect)
!Purpose: fetch myphi(1,:) from right or left; can be used for smooth and toroidal derivative
!Input: myphi(0:1,1:mgrid);
!       idirect=-1,1: direction of data flow, fetch myphi(:,1,:) from  right, left and put to myrecv
!Output: myrecv, with A. cshift operation and B. poloidal B.C. applied
!Workspace: mysend
!Convention: (1,:) is source, ghost (0,:) is destination 
!-------------------------------------------------------------
 use precision
 use mpi_defs
 use global_parameters
 use particle_decomp
 use field_array,only:itran,igrid
 
 implicit none
 integer idirect,ndim
 real(wp),dimension(ndim,0:1,mgrid)::myphi
 real(wp),dimension(ndim,mgrid)::myrecv,mysend

 integer icount,idest,isource,isendtag,irecvtag,ierror,istatus(MPI_STATUS_SIZE)
 integer i,ii,jt

 mysend(:,:)=myphi(:,1,:)
 myrecv=0.0_wp
 icount=mgrid*ndim

 select case(idirect)
  case(1) !fetch from left, myphi(1,:)-->myrecv
   idest=right_pe
   isource=left_pe
  case(-1) !fetch from right, myrecv<--myphi(1,:)
   idest=left_pe
   isource=right_pe
 end select 

 isendtag=myrank_toroidal
 irecvtag=isource
!--what's purpose of isendtag and irecvtag??
 call  MPI_SENDRECV(mysend,icount,mpi_Rsize,idest,isendtag,&
       myrecv,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

!--cshift of myrank=0 or ntoroidal-1
 do i=0,mpsi
  ii=igrid(i)
  jt=mtheta(i)
   if(idirect==1.and.myrank_toroidal==0)then !down-shift for zeta=0
     myrecv(:,ii+1:ii+jt)=cshift(myrecv(:,ii+1:ii+jt),-itran(i),2)
   elseif(idirect==-1.and.myrank_toroidal==ntoroidal-1)then !up-shift for zeta=2*pi
     myrecv(:,ii+1:ii+jt)=cshift(myrecv(:,ii+1:ii+jt),itran(i),2)
   else
   endif
 end do

!--poloidal B.C.
 do i=0,mpsi
  myrecv(:,igrid(i))=myrecv(:,igrid(i)+mtheta(i))
 end do
end subroutine fd_fetch_vector

!============================================================
subroutine fd_smooth_vector(myphi,ndim,ismooth,iskip)
!Purpose: smooth myphi
!Input: myphi
!Output: smoothed myphi
!Procedure: 
!	I.   phitmp(1,:)=myphi(1,:)
!	II.  phitmp(1,:)<--phism(:)<--radially smoothed phitmp(1,:)
!	III. phitmp(1,:)<--poloidally smoothed phitmp(1,:)
!	IV.  phitmp(1,:)<--recvl,recvr(mgrid)<--phitmp(1,:); sync phitmp(0:1,:)
!	V.   myphi(:,:)<--phitmp(:,:)
!	Originally, VI.   Re-level myphi00(0:mpsi), myphi(0:1,mgrid); not here but in FD_RELEVEL_SCALAR 
!Called suroutine: FD_SYNC_SCALAR, FD_FETCH_SCALAR
!Calling subroutines: EM version SOLVE_DPHIDT only
!-------------------------------------------------------------
 use precision
 use mpi_defs
 use global_parameters
 use particle_decomp
 use field_array,only:igrid,wtp1,wtp2,jtp1,jtp2
 
 implicit none
 integer ismooth,ndim,iskip
 !!real(wp),dimension(ndim,0:1,mgrid)::myphi
 real(wp) myphi(ndim,0:1,mgrid)

 integer i,j,jt,ii,ij,ip
 real(wp) sendl(ndim,mgrid),sendr(ndim,mgrid),recvl(ndim,mgrid),recvr(ndim,mgrid),&
          phism(ndim,mgrid),pright(ndim,mthetamax),phitmp(ndim,0:1,mgrid)

!$omp parallel do private(i)
     do i=1,mgrid
        phitmp(:,:,i)=myphi(:,:,i)  !phitmp=phi  ! smooth phi for numerical derivative -->E field
     enddo
 do ip=1,ismooth

!---- radial smoothing
!     do i=0,mpsi  !--different than original gts where i=1,mpsi-1
!     do i=1,mpsi-1  !--same as original gts where i=1,mpsi-1
     do i=iskip,mpsi-iskip  
        phitmp(:,:,igrid(i))=phitmp(:,:,igrid(i)+mtheta(i))
     enddo

     phism=0.0
!$omp parallel do private(i,j,ij)
     do i=2,mpsi-2
        do j=1,mtheta(i)
           ij=igrid(i)+j
!
! need to take into account uneven radial grid size dpsi when doing radial smoothing? wx
!
              phism(:,ij)=0.25_wp*((1.0_wp-wtp1(1,ij))*phitmp(:,1,jtp1(1,ij))+&
                   wtp1(1,ij)*phitmp(:,1,jtp1(1,ij)+1)+&
                   (1.0_wp-wtp1(2,ij))*phitmp(:,1,jtp1(2,ij))+&
                   wtp1(2,ij)*phitmp(:,1,jtp1(2,ij)+1))-&
                   0.0625_wp*((1.0_wp-wtp2(1,ij))*phitmp(:,1,jtp2(1,ij))+&
                   wtp2(1,ij)*phitmp(:,1,jtp2(1,ij)+1)+&
                   (1.0_wp-wtp2(2,ij))*phitmp(:,1,jtp2(2,ij))+&
                   wtp2(2,ij)*phitmp(:,1,jtp2(2,ij)+1))

        enddo
     enddo

     i=1
!$omp parallel do private(j,ij)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        phism(:,ij)=0.25_wp*((1.0_wp-wtp1(1,ij))*phitmp(:,1,jtp1(1,ij))+&
             wtp1(1,ij)*phitmp(:,1,jtp1(1,ij)+1)+&
             (1.0_wp-wtp1(2,ij))*phitmp(:,1,jtp1(2,ij))+&
             wtp1(2,ij)*phitmp(:,1,jtp1(2,ij)+1))-&
             0.0625_wp*((1.0_wp-wtp2(1,ij))*phitmp(:,1,jtp2(1,ij))+&
             wtp2(1,ij)*phitmp(:,1,jtp2(1,ij)+1)+&
             phitmp(:,1,ij))

     enddo

     i=mpsi-1
!$omp parallel do private(j,ij)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        phism(:,ij)=0.25_wp*((1.0_wp-wtp1(1,ij))*phitmp(:,1,jtp1(1,ij))+&
             wtp1(1,ij)*phitmp(:,1,jtp1(1,ij)+1)+&
             (1.0_wp-wtp1(2,ij))*phitmp(:,1,jtp1(2,ij))+&
             wtp1(2,ij)*phitmp(:,1,jtp1(2,ij)+1))-&
             0.0625_wp*((1.0_wp-wtp2(2,ij))*phitmp(:,1,jtp2(2,ij))+&
             wtp2(2,ij)*phitmp(:,1,jtp2(2,ij)+1)+&
             phitmp(:,1,ij))
!modified to the same form as i=1; otherwise cause different results between subroutine calling and direct calling, zx
!!             0.0625_wp*(phitmp(:,1,ij)+&
!!             (1.0_wp-wtp2(2,ij))*phitmp(:,1,jtp2(2,ij))+&
!!             wtp2(2,ij)*phitmp(:,1,jtp2(2,ij)+1))

     enddo

     i=0
!$omp parallel do private(j,ij)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        phism(:,ij)=0.5_wp*((1.0_wp-wtp1(1,ij))*phitmp(:,1,jtp1(1,ij))+&
             wtp1(1,ij)*phitmp(:,1,jtp1(1,ij)+1))-&
             0.125_wp*((1.0-wtp2(1,ij))*phitmp(:,1,jtp2(1,ij))+&
             wtp2(1,ij)*phitmp(:,1,jtp2(1,ij)+1))

     enddo

     i=mpsi
!$omp parallel do private(j,ij)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        phism(:,ij)=0.5_wp*((1.0_wp-wtp1(2,ij))*phitmp(:,1,jtp1(2,ij))+&
             wtp1(2,ij)*phitmp(:,1,jtp1(2,ij)+1))-&
             0.125_wp*((1.0_wp-wtp2(2,ij))*phitmp(:,1,jtp2(2,ij))+&
             wtp2(2,ij)*phitmp(:,1,jtp2(2,ij)+1))

     enddo

!$omp parallel do private(i)
     do i=1,mgrid
        phitmp(:,1,i)=0.625_wp*phitmp(:,1,i)+phism(:,i) ! 5-point smoothing
     end do
!endif!false

     do i=0,mpsi
!        phitmp(:,1,igrid(i))=phitmp(:,1,igrid(i)+mtheta(i)) !same as original 
        phitmp(:,:,igrid(i))=phitmp(:,:,igrid(i)+mtheta(i))  !same as Edward
     enddo

!---- poloidal smoothing (-0.0625 0.25 0.625 0.25 -0.0625)
!$omp parallel do private(i,ii,jt,pright)
     do i=0,mpsi
        ii=igrid(i)
        jt=mtheta(i)
        pright(:,1:jt)=phitmp(:,1,ii+1:ii+jt)
        phitmp(:,1,ii+1:ii+jt)=0.625_wp*pright(:,1:jt)+&   
               0.25_wp*(cshift(pright(:,1:jt),-1,2)+cshift(pright(:,1:jt),1,2))-&
             0.0625_wp*(cshift(pright(:,1:jt),-2,2)+cshift(pright(:,1:jt),2,2))
     enddo

     do i=0,mpsi
           phitmp(:,:,igrid(i))=phitmp(:,:,igrid(i)+mtheta(i))
     enddo

!---- parallel smoothing. idirect=1,-1 can be used without variable idirect
  call  fd_fetch_vector(phitmp,recvl,sendr,ndim,1)
  call  fd_fetch_vector(phitmp,recvr,sendl,ndim,-1)

  phitmp(:,1,:)=0.5_wp*phitmp(:,1,:)+0.25_wp*(recvl(:,:)+recvr(:,:))

 end do ! end do ip=1,ismooth

  call fd_sync_vector(phitmp,ndim,iskip) !--iskip=0 to cmp with Edward solve_dphidt
  myphi(:,:,:)=phitmp(:,:,:)

!--re-level here not done, see original SMOOTH_PHI 
  !--need more to do
end subroutine fd_smooth_vector
!============================================================
subroutine fd_smooth_scalar(myphi,ismooth)
!Purpose: smooth myphi
!Input: myphi
!Output: smoothed myphi
!Procedure: 
!	I.   phitmp(1,:)=myphi(1,:)
!	II.  phitmp(1,:)<--phism(:)<--radially smoothed phitmp(1,:)
!	III. phitmp(1,:)<--poloidally smoothed phitmp(1,:)
!	IV.  phitmp(1,:)<--recvl,recvr(mgrid)<--phitmp(1,:); sync phitmp(0:1,:)
!	V.   myphi(:,:)<--phitmp(:,:)
!	Originally, VI.   Re-level myphi00(0:mpsi), myphi(0:1,mgrid); not here but in FD_RELEVEL_SCALAR 
!Called suroutine: FD_SYNC_SCALAR, FD_FETCH_SCALAR
!-------------------------------------------------------------
 use precision
 use global_parameters
 use particle_decomp
 use field_array,only:igrid,wtp1,wtp2,jtp1,jtp2
 
 implicit none
 integer ismooth
 real(wp),dimension(0:1,mgrid)::myphi

 integer i,j,jt,ii,ij,ip
 real(wp) sendl(mgrid),sendr(mgrid),recvl(mgrid),recvr(mgrid),&
          phism(mgrid),pright(mthetamax),phitmp(0:1,mgrid)

!$omp parallel do private(i)
     do i=1,mgrid
        phitmp(1,i)=myphi(1,i)  !phitmp=phi  ! smooth phi for numerical derivative -->E field
     enddo

 do ip=1,ismooth

!---- radial smoothing
!     do i=0,mpsi  !--different than original gts where i=1,mpsi-1
     do i=1,mpsi-1  !--same as original gts where i=1,mpsi-1
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

!---- poloidal smoothing (-0.0625 0.25 0.625 0.25 -0.0625)
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

!---- parallel smoothing. idirect=1,-1 can be used without variable idirect??
  call  fd_fetch_scalar(phitmp,recvl,sendr,1)
  call  fd_fetch_scalar(phitmp,recvr,sendl,-1)

!  phitmp(1,:)=0.5_wp*(recvl(:)+recvr(:))
  phitmp(1,:)=0.5_wp*phitmp(1,:)+0.25_wp*(recvl(:)+recvr(:))

 end do ! end do ip=1,ismooth

 call fd_sync_scalar(phitmp)

! poloidal BC
!!  do i=0,mpsi
!!     phitmp(:,igrid(i))=phitmp(:,igrid(i)+mtheta(i))
!!  enddo

  myphi(:,:)=phitmp(:,:)

!--re-level here not done, see original SMOOTH_PHI 
  !--need more to do
end subroutine fd_smooth_scalar
!============================================================
subroutine fd_calc_deriv_par(myphi,myderiv,masssolver)
!Purpose: calculate derivative of myphi
!Input: myphi
!Output: myderiv=parallel derivative of myphi
!Called subroutines: FD_FETCH_SCALAR,FD_SYNC_SCALAR
!-------------------------------------------------------------
 use precision
 use global_parameters
! use particle_decomp
! use field_array,only:igrid,qtinv,psi_grid,wtp1,jtp1

 use solver_data_class

 implicit none

 type(solver_data) :: masssolver

 real(wp),dimension(0:1,mgrid)::myphi,myderiv

end subroutine fd_calc_deriv_par
!============================================================
subroutine fd_calc_deriv(myphi,myevector,masssolver)
!Purpose: calculate derivative of myphi
!Input: myphi
!Output: myevector=(d/da,d/dtheta,d/dzeta)myphi
!Procedure:
!	I. calculation of d/da,d/dtheta using direct or FEM method depending on order
!	II.calclation of d/dzeta: first parallel derivative and then d/dzeta 
!Called subroutines: FD_FETCH_SCALAR,FD_SYNC_VECTOR
!-------------------------------------------------------------
 use precision
 use mpi_defs
 use global_parameters
 use particle_decomp
 use field_array,only:igrid,qtinv,psi_grid,wtp1,jtp1,itran
 
 use solver_data_class

 implicit none

 type(solver_data) :: masssolver

 real(wp),dimension(0:1,mgrid)::myphi
 real(wp),dimension(3,0:1,mgrid)::myevector

 integer i,j,k,ij,ierror,jt,ii
 integer icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
 real(wp) difft(0:mpsi),diffz,delq,&
      sendl(mgrid),sendr(mgrid),recvl(mgrid),recvr(mgrid),pleft(mthetamax),pright(mthetamax), sendrs(3, mgrid), recvls(3, mgrid)
 real(wp) dpsi_inv,e_psi2,e_psi1,dpsi2_inv,dpsi1,dpsi2,dpsi1_inv !,dum,tdum  !wx
 real(wp) ftmp(0:1, mgrid);

! finite difference for e-field in equilibrium unit
  difft=0.5_wp/deltat
  diffz=0.5_wp/deltaz
!$omp parallel do private(i,j,k)
  do i=1,mgrid
     do j=0,1
        do k=1,3
           myevector(k,j,i)=0.0_wp
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

        e_psi2=((1.0_wp-wtp1(1,ij))*myphi(1,jtp1(1,ij))+&
               wtp1(1,ij)*myphi(1,jtp1(1,ij)+1)-myphi(1,ij))*dpsi2_inv

        e_psi1=(myphi(1,ij)-((1.0_wp-wtp1(2,ij))*myphi(1,jtp1(2,ij))+&
               wtp1(2,ij)*myphi(1,jtp1(2,ij)+1)))*dpsi1_inv

!           myevector(1,1,ij)=(e_psi1+e_psi2)*0.5
        myevector(1,1,ij)=(e_psi1*dpsi2+e_psi2*dpsi1)*dpsi_inv  ! note: evector(1,*,*)=d\delta-phi/da, ZF_SEPARATION=1
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
        myevector(2,1,ij)=difft(i)*(myphi(1,ij+1-jt)-myphi(1,ij-1))
     enddo
  enddo

else                  ! using FEM to d_phi/d_a and d_phi/d_theta when using FEM solver

  sendl(:)=myphi(1,:)
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


  myevector(1,1,:)=recvr(:)
  myevector(2,1,:)=recvl(:)

end if
 
!---- parallel derivative. idirect=1,-1 
  call  fd_fetch_scalar(myphi,recvl,sendr,1)
  call  fd_fetch_scalar(myphi,recvr,sendl,-1)

!!  myevector(3,1,:)=(recvr(:)-recvl(:))*diffz   !different than original where (1,igrid(i)) not included!!
 do i=1,mpsi-1
 ! d_phi/d_zeta at a=const. and theta-zeta/q_N=const.
     do j=1,mtheta(i)
        ij=igrid(i)+j
        myevector(3,1,ij)=(recvr(ij)-recvl(ij))*diffz
     enddo
 enddo

! adjust the difference between safety factor q and qtinv for fieldline coordinate
!omp parallel do private(i,j,delq,ij)
  do i=1,mpsi-1
     delq=qtinv(i)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        myevector(3,1,ij)=myevector(3,1,ij)-delq*myevector(2,1,ij) ! --> now evector(3...=d_phi/d_zeta
                     ! if delq=(1.0/q_grid(i)-qtinv(i)), evector(3...=b.grad_phi/b.grad_zeta
     enddo
  enddo

  ! call smooth_fft_vector(myevector);

! toroidal end point, pack end point (k=mzeta) data 
! send E to right and receive from left
!$omp parallel do private(i)
  do i=1,mgrid
     sendrs(1,i)=myevector(1,1,i)
     sendrs(2,i)=myevector(2,1,i)
     sendrs(3,i)=myevector(3,1,i)
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
        myevector(1,0,ii+1:ii+jt)=cshift(recvls(1,ii+1:ii+jt),-itran(i))
        myevector(2,0,ii+1:ii+jt)=cshift(recvls(2,ii+1:ii+jt),-itran(i))
        myevector(3,0,ii+1:ii+jt)=cshift(recvls(3,ii+1:ii+jt),-itran(i))
     enddo
  else
!$omp parallel do private(i)
     do i=1,mgrid
        myevector(1,0,i)=recvls(1,i)
        myevector(2,0,i)=recvls(2,i)
        myevector(3,0,i)=recvls(3,i)
     enddo
  endif

  
  
! poloidal end point
!omp parallel do private(i,j)
  do i=1,mpsi-1
     myevector(1,0,igrid(i))=myevector(1,0,igrid(i)+mtheta(i))
     myevector(2,0,igrid(i))=myevector(2,0,igrid(i)+mtheta(i))
     myevector(3,0,igrid(i))=myevector(3,0,igrid(i)+mtheta(i))

     myevector(1,1,igrid(i))=myevector(1,1,igrid(i)+mtheta(i))
     myevector(2,1,igrid(i))=myevector(2,1,igrid(i)+mtheta(i))
     myevector(3,1,igrid(i))=myevector(3,1,igrid(i)+mtheta(i))
  enddo
! now evector(:,0,:) and evector(:,:,igrid(i)) have valid data 
! evector(1:3,*,*)=0 at radial boundaries


end subroutine fd_calc_deriv 
!============================================================
subroutine fd_sync_scalar(myphi)
!Purpose: send to right and receive from left, to sync field in neighbouring poloidal planes
!Input: myphi(0:1,1:mgrid);
!Procedure: A. send (1,:) to right (0,:), i.e., (1,:)-->(0,:); 
!           B. For myrank_toroidal==0, (1,:) from left is cshifted by -itran
!Convention: (1,:) is source, (0,:) is destination (ghost)
!Note: fd_sync_scalar is less used as FD_SYNC_VECTOR
!-------------------------------------------------------------
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use field_array,only:itran,igrid

   implicit none
   real(wp),dimension(0:1,mgrid)::myphi
   real(wp) recvls(mgrid), sendrs(mgrid)
   integer i,ij,ii,jt

   integer icount,idest,isource,isendtag,irecvtag,ierror,istatus(MPI_STATUS_SIZE)

   sendrs(:)=myphi(1,:)
   recvls=0.0_wp

   icount=mgrid
   idest=right_pe
   isource=left_pe
   isendtag=myrank_toroidal
   irecvtag=isource
   call MPI_SENDRECV(sendrs,icount,mpi_Rsize,idest,isendtag,&
      recvls,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

   if(myrank_toroidal==0) then
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         myphi(0,ii+1:ii+jt)=cshift(recvls(ii+1:ii+jt),-itran(i)) !original SMOOTH_PHI has pleft, pright in cshift.??
      end do
   else
      do i=1,mgrid
         myphi(0,i)=recvls(i)
      end do
   end if
   !--poloidal B.C.
   do i=0,mpsi
      myphi(0,igrid(i))=myphi(0,igrid(i)+mtheta(i))
      myphi(1,igrid(i))=myphi(1,igrid(i)+mtheta(i))
   end do
end subroutine fd_sync_scalar

!============================================================
subroutine fd_set_vector(myphi)
!Purpose: set myphi with testing value 
!Input: myphi(0:1,1:mgrid)
!Output: myphi<--fijc+myrank_partd*1d-4+myrank_toroidal*1d-7; +0.1 for myphy(1,:)
!-------------------------------------------------------------
 use precision
 use mpi_defs
 use global_parameters
 use particle_decomp
 use field_array !might be not used

 implicit none
 real(wp),dimension(0:1,mgrid)::myphi
 integer fic,fjc,fijc

end subroutine fd_set_vector
!============================================================
subroutine fd_sync_vector(myevector,ndim,iskip)
!Purpose: send to right and receive from left, to sync field in neighbouring poloidal planes
!Input: myphi(1:3,0:1,1:mgrid);
!	iskip: # of skipped grid at inner/outer boundary when doing sync and poloidal bc check
!  convention: send (:,1,:) to right (:,0,:), i.e., (:,1,:)-->(:,0,:)
!-------------------------------------------------------------
 use precision
 use mpi_defs
 use global_parameters
 use particle_decomp
 use field_array,only:itran,igrid
 
 implicit none
 integer ndim,iskip
 real(wp),dimension(ndim,0:1,mgrid)::myevector
 real(wp) recvls(ndim,mgrid),sendrs(ndim,mgrid)
 integer icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE),ierror,&
         i,ij,ii,jt

! toroidal end point, pack end point (k=mzeta) data
! send E to right and receive from left
!! $omp parallel do private(i)
!!  do i=1,mgrid
!!     sendrs(1,i)=myevector(1,1,i)
!!     sendrs(2,i)=myevector(2,1,i)
!!     sendrs(3,i)=myevector(3,1,i)
!!     recvls(1,i)=0.0
!!     recvls(2,i)=0.0
!!     recvls(3,i)=0.0
!!  enddo
  sendrs(:,:)=myevector(:,1,:)

  icount=ndim*mgrid
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendrs,icount,mpi_Rsize,idest,isendtag,&
       recvls,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
! unpack end point data for k=0
  if(myrank_toroidal==0)then
!$omp parallel do private(i,ii,jt)
!!     do i=1,mpsi-1 !--same as original gts but confusing 
     do i=iskip,mpsi-iskip !zx
        ii=igrid(i)
        jt=mtheta(i)
!!        myevector(1,0,ii+1:ii+jt)=cshift(recvls(1,ii+1:ii+jt),-itran(i))
!!        myevector(2,0,ii+1:ii+jt)=cshift(recvls(2,ii+1:ii+jt),-itran(i))
!!        myevector(3,0,ii+1:ii+jt)=cshift(recvls(3,ii+1:ii+jt),-itran(i))
        myevector(:,0,ii+1:ii+jt)=cshift(recvls(:,ii+1:ii+jt),-itran(i),2)
     enddo
  else
!! $omp parallel do private(i)
!!     do i=1,mgrid
!!        myevector(1,0,i)=recvls(1,i)
!!        myevector(2,0,i)=recvls(2,i)
!!        myevector(3,0,i)=recvls(3,i)
!!     enddo
        myevector(:,0,:)=recvls(:,:)
  endif

! poloidal end point
!omp parallel do private(i,j)
!notice: i=0,mpsi not treated!!!
!!  do i=1,mpsi-1 !--same as original gts but confusing 
  do i=iskip,mpsi-iskip !zx
!!     myevector(1,0,igrid(i))=myevector(1,0,igrid(i)+mtheta(i))
!!     myevector(2,0,igrid(i))=myevector(2,0,igrid(i)+mtheta(i))
!!     myevector(3,0,igrid(i))=myevector(3,0,igrid(i)+mtheta(i))
     myevector(:,:,igrid(i))=myevector(:,:,igrid(i)+mtheta(i))

!!     myevector(1,1,igrid(i))=myevector(1,1,igrid(i)+mtheta(i))
!!     myevector(2,1,igrid(i))=myevector(2,1,igrid(i)+mtheta(i))
!!     myevector(3,1,igrid(i))=myevector(3,1,igrid(i)+mtheta(i))
  enddo
! now evector(:,0,:) and evector(:,:,igrid(i)) have valid data
! evector(1:3,*,*)=0 at radial boundaries 
end subroutine fd_sync_vector

!============================================================
subroutine fd_add_scalar(myphi)
!-------------------------------------------------------------
 use precision
 use mpi_defs
 use global_parameters
 use particle_decomp
 use field_array !might be not used

end subroutine fd_add_scalar
!============================================================
!!subroutine fd_relevel_scalar(myphi,myphi00)
subroutine fd_relevel_scalar(myphi)
!Purpose: re-level myphi, calculate re-leveled myphi00 for e.g. ZF_SEPARATION==0, order==0 (different in SMOOTH_PHI, POISSON, FIELD)
!Input: myphi, myphi00
!Output: re-levelled myphi, myphi00
!Procedure: 
!	I.  calc. phi00(0), then re-level: myphi(1,:)<--myphi(1,:)-phi00(0)
!	Not included but put to FD_SURF_AVG: II. calc. phi00(:) using re-leveled myphi(:,1)
!to do: need verify 
!-------------------------------------------------------------
 use mpi_defs
 use precision
 use global_parameters
 use particle_decomp
 use field_array
 use  diagnosis_array,only:jaco_sum

 implicit none
 real(wp),dimension(0:1,mgrid)::myphi
!! real(wp),dimension(0:mpsi)::myphi00
 integer i,j,ij,ierror
 real(wp) dum, tdum
!! real(wp) adum(0:mpsi)
! Calculate zonal flow component of phi, <phi>
! Note: phi and phi00 are calculated at particle coordinates (x), but not guiding center coordinates (R)

! There is arbitrary constant in the potential;
! Re-level potential phi to make <phi>(a0)=0 at inner boundary 12/07

!   calculate <phi> at inner boundary a=a0
     dum=0.0_wp
!$omp parallel private(j,ij) reduction(+:dum)
!$omp do
     do j=1,mtheta(0)
        ij=igrid(0)+j
        dum=dum+dtemper(ij)*myphi(1,ij)  ! surface averaging
     enddo
!$omp end parallel
     call MPI_ALLREDUCE(dum,tdum,1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
     tdum=tdum/jaco_sum(0)    !  =<phi> at inner boundary a=a0

     myphi=myphi-tdum  ! Re-level potential phi by a constant to make <phi>(a0)=0 at inner boundary

     do i=0, mpsi
        myphi(:,igrid(i))=myphi(:,igrid(i)+mtheta(i))  !validate phi(1,igrid(i)) for safe
     enddo

! now phi(0,:) also has valid data;  phi(:,igird(i)) also has valid data

!!    adum=0.0
!!!$omp parallel private(i,j,ij) reduction(+:adum)
!!!$omp do
!!    do i=0,mpsi
!!       do j=1,mtheta(i)
!!          ij=igrid(i)+j
!!          adum(i)=adum(i)+dtemper(ij)*myphi(1,ij)  ! surface averaging
!!       enddo
!!    enddo
!! $omp end parallel
!!    call MPI_ALLREDUCE(adum,myphi00,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
!!    myphi00=myphi00/jaco_sum    !  =<phi>

end subroutine fd_relevel_scalar
!============================================================
subroutine fd_surf_avg(myphi,myphi00)
!Purpose: calculate flux surface average (not rigorous?!!) of myphi
!Input: myphi
!Output: myphi00<--flux surface average of myphi
!Used variable (dependence): dtemper, jaco_sum
!-------------------------------------------------------------
 use precision
 use mpi_defs
 use global_parameters
 use particle_decomp
 use field_array,only:dtemper,igrid
 use diagnosis_array,only:jaco_sum
 
 implicit none
 real(wp),dimension(0:1,mgrid)::myphi
 real(wp),dimension(0:mpsi)::myphi00
 integer i,j,ij,ierror
 real(wp) adum(0:mpsi)

! if(mype==0) write(*,*)'--field surface average --' 
 myphi00(:)=0.0_wp
! --following copied from SMOOTH_PHI 
     adum=0.0
!$omp parallel private(i,j,ij) reduction(+:adum)
!$omp do
    do i=0,mpsi
       do j=1,mtheta(i)
          ij=igrid(i)+j
          adum(i)=adum(i)+dtemper(ij)*myphi(1,ij)  ! surface averaging
       enddo
    enddo
!$omp end parallel
    call MPI_ALLREDUCE(adum,myphi00,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
    myphi00=myphi00/jaco_sum    !  =<phi>
 
end subroutine fd_surf_avg

!============================================================
subroutine fd_calc_deriv00(myphi00,myphip00)
!Purpose: direct numerical derivative to get d<phi>/da
!Input: myphip
!Output: myphip00
!Procedure:
!	I.  calc. phip00(1:mpsi-1) from phi00(0:mpsi)
!	II. set phip00(0)=0; set phip00(mpsi) to 0 or from phi00(mpsi-1:mpsi) for IFNEOCLASSICAL==0/1
!-------------------------------------------------------------
 use precision
 use global_parameters
 use field_array,only:psi_grid
 implicit none
 real(wp),dimension(0:mpsi)::myphi00,myphip00
 integer i
 real(wp) dpsi1,dpsi2,dpsi1_inv,dpsi2_inv,dpsi_inv,e_psi1,e_psi2

! direct numerical derivative to get d<phi>/da
!$omp parallel do private(i,dpsi1,dpsi2,dpsi1_inv,dpsi2_inv,dpsi_inv,e_psi2,e_psi1)
    do i=1,mpsi-1
       dpsi1=psi_grid(i)-psi_grid(i-1)
       dpsi2=psi_grid(i+1)-psi_grid(i)
       dpsi1_inv=1.0_wp/dpsi1
       dpsi2_inv=1.0_wp/dpsi2
       dpsi_inv=1.0_wp/(dpsi1+dpsi2)
       e_psi2=(myphi00(i+1)-myphi00(i))*dpsi2_inv
       e_psi1=(myphi00(i)-myphi00(i-1))*dpsi1_inv
       myphip00(i)=(e_psi1*dpsi2+e_psi2*dpsi1)*dpsi_inv
    enddo
! boundary conditions d<phi>/da=0 (not phi=0 at both boundaries)
    myphip00(0)=0.0_wp
    !  myphip00(mpsi)=0.0_wp
    if(ifneoclassical==0) then
       myphip00(mpsi)=0.0_wp
    else        ! d<phi>/da /= 0 when NC physics included
       myphip00(mpsi)=(myphi00(mpsi)-myphi00(mpsi-1))/(psi_grid(mpsi)-psi_grid(mpsi-1))
!!!!       if(order <= 0) evector(1,1,igrid(mpsi):(igrid(mpsi)+mtheta(mpsi)))=phip00(mpsi)   ! also set dphi/dtheta=dphi/dzeta=0 at mpsi; zx move this outside since it's not relavant to DERIV00
    endif

end subroutine fd_calc_deriv00
!============================================================
subroutine fd_dump_scalar1d(myphi00,ivar)
!Purpose: dump (write) myphi to files 
!Input: myphi00(0:mpsi)
!Output: dump(ivar)_(myrank_partd)_(myrank_toroidal) files 
!-------------------------------------------------------------
  use  precision,only:wp
  use global_parameters,only:mpsi
  use particle_decomp
!  use particle_array
!  use field_array
!  use esi_interface
!  use diagnosis_array

  implicit none
  integer ivar
  real(wp),dimension(0:mpsi)::myphi00
  character(len=30) cdum
!  integer fic,fjc,fijc,ni
!  real(wp) showdata(0:1,5,0:mpsi)

!--set file name and write data 
  write(cdum,'("dump",i2.2,"_",i4.4,"_",i3.3)')ivar,myrank_partd,myrank_toroidal
  open(222,file=cdum,status='replace')
  write(222,110)myphi00(:)
  close(222)

110 format(10(1x,F13.6))
end subroutine fd_dump_scalar1d
!============================================================
subroutine fd_dump_scalar2d(myphi,ivar)
!Purpose: dump (write) myphi to files 
!Input: myphi(0:1,1:mgrid)
!Output: dump(myrank_partd)_(myrank_toroidal) files 
!-------------------------------------------------------------
  use global_parameters
  use particle_decomp
!  use particle_array
  use field_array
  use esi_interface
  use diagnosis_array

  implicit none
  integer ivar
  real(wp),dimension(0:1,mgrid)::myphi
  character(len=30) cdum
  integer fic,fjc,fijc,ni
  real(wp) showdata(0:1,5,0:mpsi)
!--copy myphi to reduced data set showdata
  do fic=0,mpsi
   fijc=igrid(fic)
   showdata(:,1,fic)=myphi(:,fijc)
   showdata(:,2,fic)=myphi(:,fijc+1)

   showdata(:,3,fic)=myphi(:,fijc+int(mtheta(fic)/2))

   fijc=fijc+mtheta(fic)
   showdata(:,4,fic)=myphi(:,fijc-1)
   showdata(:,5,fic)=myphi(:,fijc)
  end do
!--set file name and write reduced/original data set
  write(cdum,'("dump",i2.2,"_",i4.4,"_",i3.3)')ivar,myrank_partd,myrank_toroidal
  open(222,file=cdum,status='replace')
  do fic=0,mpsi
   write(222,110)showdata(0,:,fic)
   write(222,110)showdata(1,:,fic)
  end do
 write(222,*)
 write(222,110)myphi(0,:)
 write(222,110)myphi(1,:)
  close(222)

110 format(10(1x,F13.6))
end subroutine fd_dump_scalar2d

!============================================================
subroutine fd_set_scalar(myphi)
!Purpose: set myphi with testing value 
!Input: myphi(0:1,1:mgrid)
!Output: myphi<--fijc+myrank_partd*1d-4+myrank_toroidal*1d-7; +0.1 for myphy(1,:)
!-------------------------------------------------------------
 use precision
 use mpi_defs
 use global_parameters
 use particle_decomp
 use field_array,only:igrid 

 implicit none
 real(wp),dimension(0:1,mgrid)::myphi
 integer fic,fjc,fijc

 do fic=0,mpsi
  fijc=igrid(fic)
  do fjc=0,mtheta(fic)
   myphi(:,fijc)=fijc*1.d-5+myrank_partd*1.d3+myrank_toroidal*1.d0
   myphi(1,fijc)=myphi(1,fijc)+1.d5
   fijc=fijc+1
  end do
 end do
!--poloidal B.C.
 do fic=0,mpsi
  myphi(:,igrid(fic))=myphi(:,igrid(fic)+mtheta(fic))
 end do
end subroutine fd_set_scalar

!============================================================
subroutine fd_g2i_vector(evector,wpi)
   !Purpose:
   !Input:
   !Output:
   !Procedure:
   !-------------------------------------------------------------
   use precision
   use global_parameters,only:mi,mgrid
   use particle_array,only:wzion,wpion,wtion0,wtion1,jtion0,jtion1
   ! use field_array

   implicit none
   real(wp) evector(3,0:1,mgrid),wpi(3,mi)

   integer m,ij, larmor
   real(wp) wz0,wz1,wp0,wp1,wt00,wt10,wt01,wt11,e1,e2,e3

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

end subroutine fd_g2i_vector

!============================================================
subroutine fd_g2i_scalar(y_g, y_p)
   !Purpose:
   !Input:
   !Output:
   !Procedure:
   !-------------------------------------------------------------
   use precision
   use global_parameters,only:mi,mgrid
   use particle_array,only:wzion,wpion,wtion0,wtion1,jtion0,jtion1
   ! use field_array

   implicit none
   real(wp) y_g(0:1,mgrid),y_p(mi)

   integer m,ij, larmor
   real(wp) wz0,wz1,wp0,wp1,wt00,wt10,wt01,wt11,e1

!$omp parallel do private(m,e1,wz1,wz0,larmor,ij,wp0,wt00,wt10,wp1,wt01,wt11)
   do m=1,mi
      e1=0.0_wp
      wz1=wzion(m)
      wz0=1.0_wp-wz1

      do larmor=1,4

         ij=jtion0(larmor,m)
         wp0=1.0_wp-wpion(larmor,m)
         wt00=1.0_wp-wtion0(larmor,m)
         e1=e1+wp0*wt00*(wz0*y_g(0,ij)+wz1*y_g(1,ij))

         ij=ij+1
         wt10=1.0_wp-wt00
         e1=e1+wp0*wt10*(wz0*y_g(0,ij)+wz1*y_g(1,ij))

         ij=jtion1(larmor,m)
         wp1=1.0_wp-wp0
         wt01=1.0_wp-wtion1(larmor,m)
         e1=e1+wp1*wt01*(wz0*y_g(0,ij)+wz1*y_g(1,ij))

         ij=ij+1
         wt11=1.0_wp-wt01
         e1=e1+wp1*wt11*(wz0*y_g(0,ij)+wz1*y_g(1,ij))

      enddo

      y_p(m)=0.25_wp*e1

   enddo

end subroutine fd_g2i_scalar

!============================================================
 subroutine fd_g2e_vector(evector,wpi)
!Purpose:
!Input:
!Output:
!Procedure:
!-------------------------------------------------------------
 use precision
 use global_parameters,only:me,mgrid
 use particle_array,only:wzeon,wpeon,wteon0,wteon1,jteon0,jteon1
! use field_array

  implicit none
  real(wp) evector(3,0:1,mgrid),wpi(3,me)

  integer m,ij
  real(wp) wz0,wz1,wp0,wp1,wt00,wt10,wt01,wt11,e1,e2,e3

!--make sure electron_locator is called before: (me) wpeon,jteon0,wteon0,jteon1,wteon1
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
 end subroutine fd_g2e_vector

!============================================================
 subroutine fd_g2e_scalar(y_g,y_p)
!Purpose:
!Input: y_g: scalar on grid 
!Output: y_p scalar on particle
!Procedure:
!-------------------------------------------------------------
 use precision
 use global_parameters,only:me,mgrid
 use particle_array,only:wzeon,wpeon,wteon0,wteon1,jteon0,jteon1
! use field_array

  implicit none
  real(wp) y_g(0:1,mgrid),y_p(me)

  integer m,ij
  real(wp) wz0,wz1,wp0,wp1,wt00,wt10,wt01,wt11,e1

  call electron_locator
!--make sure electron_locator is called before: (me) wpeon,jteon0,wteon0,jteon1,wteon1
  do m=1,me
     wz1=wzeon(m)
     wz0=1.0_wp-wz1

     ij=jteon0(m)
     wp0=1.0_wp-wpeon(m)
     wt00=1.0_wp-wteon0(m)
     e1=wp0*wt00*(wz0*y_g(0,ij)+wz1*y_g(1,ij))

     ij=ij+1
     wt10=1.0_wp-wt00
     e1=e1+wp0*wt10*(wz0*y_g(0,ij)+wz1*y_g(1,ij))

     ij=jteon1(m)
     wp1=1.0_wp-wp0
     wt01=1.0_wp-wteon1(m)
     e1=e1+wp1*wt01*(wz0*y_g(0,ij)+wz1*y_g(1,ij))

     ij=ij+1
     wt11=1.0_wp-wt01
     e1=e1+wp1*wt11*(wz0*y_g(0,ij)+wz1*y_g(1,ij))

     y_p(m)=e1
  enddo
 end subroutine fd_g2e_scalar

