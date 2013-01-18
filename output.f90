Subroutine output
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use field_array
  implicit none
  include 'netcdf.inc'

! # of grid in R and Z direction of R-Z plane
  integer,parameter :: nr=300,nz=360
  integer i,j,k,ip,jt,kz,kk,nfile,npoint,icount,ierror,istatus,dataid,ncid,&
          dimid(2),kp,km,jtp,jtm
  real(wp) eachdata(mgrid),alldata(mgrid,mzetamax),psizeta(mtdiag,5:mpsi/2-4),&
           delt(0:mpsi),delz,tdum,zdum,rdum,r1,z1,r2,z2,t,r,wt,wtp,wtm,delr,&
           pi2_inv,poloidal(nr,nz)
  character(len=12) cdum
  character(len=4) ddum(3)

  write(ddum(1),'("dpot")')
  write(ddum(2),'("dtem")')
  write(ddum(3),'("heat")')

  do kk=1,3

! gather data to PE0
     if (kk==1)then
        eachdata(:)=phi(1,:)/(gyroradius*gyroradius)
     elseif(kk==2)then
        eachdata=dtemper
     else
        eachdata=heatflux
     endif

   ! Only one PE per toroidal domain needs to send its data. We choose the
   ! PE with rank 0 within the domain, i.e. myrank_partd=0
     if(myrank_partd==0)then
        icount=mgrid
        call MPI_GATHER(eachdata,icount,mpi_Rsize,alldata,icount,&
                        mpi_Rsize,0,toroidal_comm,ierror)
     end if

     if(myrank_partd==0 .and. myrank_toroidal==0)then  ! Equivalent to mype=0

! cell center value
        if(kk==1)then
           do i=1,mgrid
              alldata(i,1:mzetamax)=0.5*(alldata(i,1:mzetamax)+&
                   cshift(alldata(i,1:mzetamax),-1))
           enddo
        endif

! 2D data with field-line averaged around theta=0
        npoint=mzetamax/10
!        npoint=0
        delt=2.0*pi/deltat
        delz=1.0/real(mtdiag)
        do j=1,mtdiag
           zdum=delz*real(j)
           kz=max(1,min(mzetamax,1+int(zdum*real(mzetamax))))
           do i=5,mpsi/2-4
              ip=2*i-1
              tdum=10.0-zdum*qtinv(ip)
              tdum=(tdum-aint(tdum))*delt(ip)+0.5
              jt=max(0,min(mtheta(ip),int(tdum)))
              psizeta(j,i)=alldata(igrid(ip)+jt,kz)
              do k=1,npoint
                 kp=kz+k
                 jtp=jt
                 if(kp>mzetamax)then
                    kp=kp-mzetamax
                    jtp=jtp+itran(ip)
                    if(jtp>mtheta(ip))jtp=jtp-mtheta(ip)
                 endif
                 km=kz-k
                 jtm=jt
                 if(km<1)then
                    km=km+mzetamax
                    jtm=jtm-itran(ip)
                    if(jtm<0)jtm=jtm+mtheta(ip)
                 endif
                 psizeta(j,i)=psizeta(j,i)+alldata(igrid(ip)+jtp,kp)+&
                      alldata(igrid(ip)+jtm,km)
              enddo
           enddo
        enddo
        psizeta=psizeta/real(2*npoint+1)

! smooth temperature and heatflux
        if(kk>1)then
           do i=5,mpsi/2-4
              psizeta(:,i)=0.5*psizeta(:,i)+0.25*(cshift(psizeta(:,i),1)+&
                   cshift(psizeta(:,i),-1))
           enddo
           do j=1,mtdiag
              psizeta(j,:)=0.5*psizeta(j,:)+0.25*(cshift(psizeta(j,:),1)+&
                   cshift(psizeta(j,:),-1))
           enddo
        endif

! write to data file in NetCDF format
        nfile=(mstepall+istep)/ndiag
        if(nfile<10)then
           write(cdum,'(a4,"000",i1,".ncd")')ddum(kk),nfile
        elseif(nfile<100)then
           write(cdum,'(a4,"00",i2,".ncd")')ddum(kk),nfile
        elseif(nfile<1000)then
           write(cdum,'(a4,"0",i3,".ncd")')ddum(kk),nfile
        else
           write(cdum,'(a4,i4,".ncd")')ddum(kk),nfile
        endif
        
! open netcdf data file
        istatus=nf_create(cdum,nf_clobber,ncid)
! define data array dimension
        istatus=nf_def_dim(ncid,'toroidal_grid',mtdiag,dimid(1))
        istatus=nf_def_dim(ncid,'radial_grid',mpsi/2-8,dimid(2))
! define data array id
        istatus=nf_def_var(ncid,ddum(kk),nf_real,2,dimid,dataid)
! end of netcdf definition
        istatus=nf_enddef(ncid)
! write data
        istatus=nf_put_var_real(ncid,dataid,psizeta)
        istatus=nf_close(ncid)        

     endif
  enddo

  return

! write 2D potential at poloidal plane at zeta=0 for Nazikian
  if(mype==0)then
     r1=0.2*a   ! starting R coordinate measured from magnetic axis  ! ????
     z1=-0.36*a   ! starting Z coordinate measured from the mid-plane
     pi2_inv=0.5/pi
     delr=1.0/deltar
     poloidal=0.0
     do i=1,nr       ! ????
        do j=1,nz
           r2=r1+gyroradius*real(i-1)   ! grid size=gyroradius
           z2=z1+gyroradius*real(j-1)
           r=sqrt(r2*r2+z2*z2)
           t=atan(z2/r2)
           t=pi2_inv*(t-r*sin(t))+10.0
           tdum=t-aint(t)

           ip=int((r-a0)*delr)
           wt=(r-a0)*delr-real(ip)

           jtm=int(tdum*delt(ip))
           wtm=tdum*delt(ip)-real(jtm)

           jtp=int(tdum*delt(ip+1))
           wtp=tdum*delt(ip+1)-real(jtp)

           poloidal(i,j)=wt*(wtp*phi(0,igrid(ip+1)+jtp+1)+&
                       (1.0-wtp)*phi(0,igrid(ip+1)+jtp))+&
                   (1.0-wt)*(wtm*phi(0,igrid(ip)+jtm+1)+&
                       (1.0-wtm)*phi(0,igrid(ip)+jtm))
        enddo
     enddo

! smooth
     do i=1,nr
        poloidal(i,:)=0.5*poloidal(i,:)+&
             0.25*(cshift(poloidal(i,:),1)+cshift(poloidal(i,:),-1))
     enddo
     do j=1,nz
        poloidal(:,j)=0.5*poloidal(:,j)+&
             0.25*(cshift(poloidal(:,j),1)+cshift(poloidal(:,j),-1))
     enddo

! writeout data in NetCDF
        if(nfile<10)then
           write(cdum,'(a4,"000",i1,".ncd")')'deln',nfile
        elseif(nfile<100)then
           write(cdum,'(a4,"00",i2,".ncd")')'deln',nfile
        elseif(nfile<1000)then
           write(cdum,'(a4,"0",i3,".ncd")')'deln',nfile
        else
           write(cdum,'(a4,i4,".ncd")')'deln',nfile
        endif
        
! open netcdf data file
        istatus=nf_create(cdum,nf_clobber,ncid)
! define data array dimension
        istatus=nf_def_dim(ncid,'R',nr,dimid(1))
        istatus=nf_def_dim(ncid,'Z',nz,dimid(2))
! define data array id
        istatus=nf_def_var(ncid,'deln',nf_real,2,dimid,dataid)
! end of netcdf definition
        istatus=nf_enddef(ncid)
! write data
        istatus=nf_put_var_real(ncid,dataid,poloidal)
        istatus=nf_close(ncid)        
  endif

end subroutine output
