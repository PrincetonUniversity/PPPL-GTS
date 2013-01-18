subroutine get_3D_flow(da,dgq,dgf,wpi2)
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

  real(wp) da(mi),dgq(mi),dgf(mi),wpi2(mi)

  real(wp) dnitmp(0:1,mgrid),buff(0:1,mgrid)

!  *****************************************
!  calculate ion flow delta-Vi
!
!  BORROW array pflux_e, wzion, wtion0, wtion1, jtion0, jtion1
!  BE SURE no memory conflict (!!!)
!

! ion_locator is called to give ion locations 
  call ion_locator0

     pflux_e=0.0_wp
     pflux_i=0.0_wp

do i=1,4  ! corresponding to a, theta, zeta; using do-loop to save memory

#ifdef _OPENMP
!$omp parallel private(dnitmp)
  dnitmp=0.   ! Set array elements to zero
#endif
!$omp do private(m,weight,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij,rdum)
  do m=1,mi
     weight=zion0(5,m)

     wz1=weight*wzion(m)
     wz0=weight-wz1

        wp1=wpion(1,m)
        wp0=1.0_wp-wp1

        wt10=wp0*wtion0(1,m)
        wt00=wp0-wt10

        wt11=wp1*wtion1(1,m)
        wt01=wp1-wt11
  
        if(i==1) then 
           rdum=da(m)
        else if(i==2) then 
           rdum=dgq(m)
        else if(i==3) then
           rdum=dgf(m)
        else
           rdum=wpi2(m)  ! added Oct. 2012 for calculating delta-Ti on grid
        end if

#ifdef _OPENMP
! Use thread-private temporary array dnitmp to store the results
        ij=jtion0(1,m)
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt00*rdum
        dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt00*rdum

        ij=ij+1
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt10*rdum
        dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt10*rdum

        ij=jtion1(1,m)
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt01*rdum
        dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt01*rdum

        ij=ij+1
        dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt11*rdum
        dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt11*rdum
#else
! If no loop-level parallelism, use original algorithm (write directly
! into array "pflux_e(*,*,*)".
        ij=jtion0(1,m)
        if(i < 4) then
           pflux_e(i,0,ij) = pflux_e(i,0,ij) + wz0*wt00*rdum
           pflux_e(i,1,ij) = pflux_e(i,1,ij) + wz1*wt00*rdum

           ij=ij+1
           pflux_e(i,0,ij) = pflux_e(i,0,ij) + wz0*wt10*rdum
           pflux_e(i,1,ij) = pflux_e(i,1,ij) + wz1*wt10*rdum
        else
           pflux_i(1,0,ij) = pflux_i(1,0,ij) + wz0*wt00*rdum
           pflux_i(1,1,ij) = pflux_i(1,1,ij) + wz1*wt00*rdum

           ij=ij+1
           pflux_i(1,0,ij) = pflux_i(1,0,ij) + wz0*wt10*rdum
           pflux_i(1,1,ij) = pflux_i(1,1,ij) + wz1*wt10*rdum
        end if


        ij=jtion1(1,m)
        if(i < 4) then
           pflux_e(i,0,ij) = pflux_e(i,0,ij) + wz0*wt01*rdum
           pflux_e(i,1,ij) = pflux_e(i,1,ij) + wz1*wt01*rdum

           ij=ij+1
           pflux_e(i,0,ij) = pflux_e(i,0,ij) + wz0*wt11*rdum
           pflux_e(i,1,ij) = pflux_e(i,1,ij) + wz1*wt11*rdum
        else
           pflux_i(1,0,ij) = pflux_i(1,0,ij) + wz0*wt01*rdum
           pflux_i(1,1,ij) = pflux_i(1,1,ij) + wz1*wt01*rdum

           ij=ij+1
           pflux_i(1,0,ij) = pflux_i(1,0,ij) + wz0*wt11*rdum
           pflux_i(1,1,ij) = pflux_i(1,1,ij) + wz1*wt11*rdum
        end if
#endif

  enddo

#ifdef _OPENMP
!$omp critical
  do ij=1,mgrid
     if(i < 4) then
        pflux_e(i,0,ij) = pflux_e(i,0,ij) + dnitmp(0,ij)
        pflux_e(i,1,ij) = pflux_e(i,1,ij) + dnitmp(1,ij)
     else
        pflux_i(1,0,ij) = pflux_i(1,0,ij) + dnitmp(0,ij)
        pflux_i(1,1,ij) = pflux_i(1,1,ij) + dnitmp(1,ij)
     end if
  enddo
!$omp end critical
#endif
!$omp end parallel

! If we have a particle decomposition on the toroidal domains, do a reduce
! operation to add up all the contributions to charge density on the grid
  if(npartdom>1)then
   !$omp parallel do private(ij)
    do ij=1,mgrid
       if(i < 4) then
          dnitmp(0,ij)=pflux_e(i,0,ij)
          dnitmp(1,ij)=pflux_e(i,1,ij)
       else
          dnitmp(0,ij)=pflux_i(1,0,ij)
          dnitmp(1,ij)=pflux_i(1,1,ij)
       end if
    enddo
    buff=0._wp
    call MPI_ALLREDUCE(dnitmp,buff,(mgrid*2),mpi_Rsize,&
                       MPI_SUM,partd_comm,ierror)
    do ij=1,mgrid
       if(i < 4) then
          pflux_e(i,0,ij)=buff(0,ij)
          pflux_e(i,1,ij)=buff(1,ij)
       else
          pflux_i(1,0,ij)=buff(0,ij)
          pflux_i(1,1,ij)=buff(1,ij)
       end if
    enddo
  endif

! poloidal end cell
  do j=0,mpsi
     if(i < 4) then
        pflux_e(i,:,igrid(j)+mtheta(j))=pflux_e(i,:,igrid(j)+mtheta(j))+pflux_e(i,:,igrid(j))
        pflux_e(i,:,igrid(j))=pflux_e(i,:,igrid(j)+mtheta(j))
     else
        pflux_i(1,:,igrid(j)+mtheta(j))=pflux_i(1,:,igrid(j)+mtheta(j))+pflux_i(1,:,igrid(j))
        pflux_i(1,:,igrid(j))=pflux_i(1,:,igrid(j)+mtheta(j))
     end if
  enddo

! toroidal end cell
  if(i < 4) then
     sendl=pflux_e(i,0,:)
  else
     sendl=pflux_i(1,0,:)
  end if
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
        if(i < 4) then
           pflux_e(i,1,ii+1:ii+jt)=pflux_e(i,1,ii+1:ii+jt)+&
                cshift(recvr(ii+1:ii+jt),itran(j))
        else
           pflux_i(1,1,ii+1:ii+jt)=pflux_i(1,1,ii+1:ii+jt)+&
                cshift(recvr(ii+1:ii+jt),itran(j))
        end if
     enddo
  else
! B.C. at zeta<2*pi is continuous
     if(i < 4) then
        pflux_e(i,1,:)=pflux_e(i,1,:)+recvr(:)
     else
        pflux_i(1,1,:)=pflux_i(1,1,:)+recvr(:)
     end if
  endif

end do   ! do i=1,4

!  NOTE: pflux_e(:,1,:) has valid data, but no valid data for pflux_e(:,0,:) and pflux_e(:,:,igrid(j))

!$omp parallel do private(ij,rdum)
  do ij=1,mgrid
     pflux_e(1,1,ij)=pflux_e(1,1,ij)*markeri(ij)  ! da/dt on grid
     pflux_e(2,1,ij)=pflux_e(2,1,ij)*markeri(ij)  ! dtheta/dt on grid
     pflux_e(3,1,ij)=pflux_e(3,1,ij)*markeri(ij)  ! dphi/dt on grid

     pflux_e(1,0,ij)=on_grid(1,ij)*pflux_e(1,1,ij)+on_grid(2,ij)*pflux_e(2,1,ij)  ! delta-Vi in R direction
     pflux_e(2,0,ij)=on_grid(5,ij)*pflux_e(3,1,ij)                                ! delta-Vi in phi (toroidal) direction
     pflux_e(3,0,ij)=on_grid(3,ij)*pflux_e(1,1,ij)+on_grid(4,ij)*pflux_e(2,1,ij)  ! delta-Vi in Z direction

     pflux_i(1,1,ij)=pflux_i(1,1,ij)*markeri(ij)  ! Ti fluctuations on grid

!      pflux_e(2,1,ij)=pflux_e(2,1,ij)*on_grid(6,ij)                                ! delta-Vi in theta (poloidal) direction; added Aug. 2012
     rdum=sqrt(gij(1,ij))   ! =|grad(a)|
     pflux_e(2,1,ij)=(pflux_e(2,1,ij)*rdum+pflux_e(1,1,ij)*(-gij(4,ij))/rdum)*Jcbn(ij)/on_grid(5,ij)  ! delta-Vi in poloidal (B_p) direction; corrected Oct. 2012 - see note
  enddo

! identical grid points:   igrid(j)=igrid(j)+mtheta(j)
  do j=0,mpsi
     pflux_e(:,:,igrid(j))=pflux_e(:,:,igrid(j)+mtheta(j))
     pflux_i(1,:,igrid(j))=pflux_i(1,:,igrid(j)+mtheta(j))
     
     rdum=vth_inv(j)*vth_inv(j)         ! = 1/Ti_0 with Ti_0 normalized 
     rdum=rdum*2.0_wp/3.0_wp
     !$omp parallel do private(i,ij)
     do i=0,mtheta(j)
        ij=igrid(j)+i
        pflux_i(1,1,ij)=pflux_i(1,1,ij)*rdum    ! = delta-Ti/Ti_0 on grid   !  added Oct. 2012
     end do
  enddo

! Record Density. Each myrank_toroidal has its own file (see open_file).
  if(istep>0 .and. myrank_partd==2) then
     write(808) istep
     write(808) pflux_e(1,0,1:mgrid), pflux_e(2,0,1:mgrid), pflux_e(3,0,1:mgrid), pflux_e(2,1,1:mgrid)   !   , pflux_i(1,1,1:mgrid)
!                   |                     |                     |                     |                  !        |
!                   |                     |                     |                     |                  !        |
!                   V_R                   V_phi                 V_Z                   V_p                !    delta-Ti/Ti_0

     call flush(808)
  end if

  phi_tmp(1:mgrid)=pflux_i(1,1,1:mgrid)    ! use phi_tmp(1:mgrid) to store Ti fluctuations and write to output DENs at SNAPSHOT


!   if(mype==0) print*, 'Write FLOWs + delta-Ti at ', 'istep=', istep,  'irk=', irk

end subroutine get_3D_flow






subroutine ion_locator0  ! Added to calculate delta-vi and delta Ti on grid
                         ! using zion0 rather than zion
  use precision
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer m,im,ii,j01,j00
  real(wp) rdum,tdum,delr,delt(0:mpsi),delz,wp1,     &
       pi2_inv,psitmp,thetatmp,zetatmp,psimax,psimin

! Borrow array wzion, wpion, jtion0, wtion0, jtion1, wtion1;
! No memory conflict since ion_locator is called right after PUSHI,
! and these array are oroginally calculated in CHARGE_ION, used
! in CHARGE_ION, CURRENT_ION and PUSHI, and is released after PUSHI

  delr=1.0_wp/deltar
  delt=2.0_wp*pi/deltat
  delz=1.0_wp/deltaz
  pi2_inv=0.5_wp/pi
  psimin=s_grid(0)
  psimax=s_grid(mpsi)

!$omp parallel do private(m,psitmp,thetatmp,zetatmp, &
!$omp& rdum,ii,wp1,im,tdum,j00,j01)
  do m=1,mi
!     psitmp=zion0(1,m)      ! =a
     psitmp=zion0(7,m)       ! =s
     thetatmp=zion0(2,m)
     zetatmp=zion0(3,m)

     wzion(m)=(zetatmp-zetamin)*delz

     rdum=psitmp
     rdum=delr*max(0.0_wp,min(s1-s0,rdum-s0))
     ii=max(0,min(mpsi-1,int(rdum)))
     wp1=rdum-real(ii,wp)
     wpion(1,m)=wp1

! inner flux surface
     im=ii
     tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
     tdum=(tdum-aint(tdum))*delt(im)
     j00=max(0,min(mtheta(im)-1,int(tdum)))
     jtion0(1,m)=igrid(im)+j00  ! pay attention to j00=0 case
     wtion0(1,m)=tdum-real(j00,wp)

! outer flux surface
     im=ii+1
     tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
     tdum=(tdum-aint(tdum))*delt(im)
     j01=max(0,min(mtheta(im)-1,int(tdum)))
     jtion1(1,m)=igrid(im)+j01  ! pay attention to j01=0 case
     wtion1(1,m)=tdum-real(j01,wp)
  enddo

end subroutine ion_locator0



subroutine ion_locator  ! Added to calculate delta-vi and delta Ti on grid
                        ! using zion rather than zion0
  use precision
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer m,im,ii,j01,j00
  real(wp) rdum,tdum,delr,delt(0:mpsi),delz,wp1,     &
       pi2_inv,psitmp,thetatmp,zetatmp,psimax,psimin

! Borrow array wzion, wpion, jtion0, wtion0, jtion1, wtion1;
! No memory conflict since ion_locator is called right after PUSHI,
! and these array are oroginally calculated in CHARGE_ION, used
! in CHARGE_ION, CURRENT_ION and PUSHI, and is released after PUSHI

  delr=1.0_wp/deltar
  delt=2.0_wp*pi/deltat
  delz=1.0_wp/deltaz
  pi2_inv=0.5_wp/pi
  psimin=s_grid(0)
  psimax=s_grid(mpsi)

!$omp parallel do private(m,psitmp,thetatmp,zetatmp, &
!$omp& rdum,ii,wp1,im,tdum,j00,j01)
  do m=1,mi
!     psitmp=zion(1,m)      ! =a
     psitmp=zion(7,m)       ! =s
     thetatmp=zion(2,m)
     zetatmp=zion(3,m)

     wzion(m)=(zetatmp-zetamin)*delz

     rdum=psitmp
     rdum=delr*max(0.0_wp,min(s1-s0,rdum-s0))
     ii=max(0,min(mpsi-1,int(rdum)))
     wp1=rdum-real(ii,wp)
     wpion(1,m)=wp1

! inner flux surface
     im=ii
     tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
     tdum=(tdum-aint(tdum))*delt(im)
     j00=max(0,min(mtheta(im)-1,int(tdum)))
     jtion0(1,m)=igrid(im)+j00  ! pay attention to j00=0 case
     wtion0(1,m)=tdum-real(j00,wp)

! outer flux surface
     im=ii+1
     tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
     tdum=(tdum-aint(tdum))*delt(im)
     j01=max(0,min(mtheta(im)-1,int(tdum)))
     jtion1(1,m)=igrid(im)+j01  ! pay attention to j01=0 case
     wtion1(1,m)=tdum-real(j01,wp)
  enddo

end subroutine ion_locator
