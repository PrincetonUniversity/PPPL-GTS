subroutine snapshot
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface
  
  implicit none

  integer,parameter :: mbin_psi=5, mbin_u=32
  integer i,ibin,ip,iu,j,jj,jm,jt,k,kz,m,nsnap,ierror,icount
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
  if(mod(istep,40) .eq. 0) call restart_write

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
102 format(e10.4)
  
105 format(e10.4,4x,6hflux2x,4x,i4)
106 format(e10.4,4x,6hflux2z,4x,i4)
107 format(e10.4,2x,i4,2x,4hpsi=,e10.4,2x,4hthe=,e10.4,2x,2hR=,e10.4,2x,2hz=,e10.4)


! REcord Potential. Each myrank_toroidal has its own file (see open_file).
  if(movie_3D_h5==0) then  ! no 3D movie data created
     if(istep>0 .and. myrank_partd==0) then
        write(777) istep
        write(777) phi(1,1:mgrid)  ! Note(5/08): phi contains zonal flow when ZF_SEPARATION=0 (solving total potential)
!        write(777) phi_tmp(1,:)   ! time averaged phi
        call flush(777)  ! otherwise, if job crashes PHI data is not completed, which cause problem for reading PHI 
     end if
  end if

! Record Density. Each myrank_toroidal has its own file (see open_file).
  if(istep>0 .and. myrank_partd==1) then
     write(888) istep
     if(ifelectron==0) then
        write(888) densityi(1,1:mgrid),phi_tmp(1:mgrid)
     else 
        write(888) densityi(1,1:mgrid),phi_tmp(1:mgrid),densitye(1,1:mgrid),heatflux(1:mgrid)   ! delta-ne output added Oct. 2012
                            !  phi_tmp:  Ti fluctuations
                            !  heatflux: Te fluctuations

!          print*, 'mype= ', mype, 'istep=', istep, phi_tmp(mgrid/2), heatflux(mgrid/2), densitye(1,mgrid/2), densityi(1,mgrid/2)

     end if
     call flush(888)
  end if

!   if(mype==0) print*, 'Write DENs + PHIs at ', 'istep=', istep,  'irk=', irk

end subroutine snapshot
