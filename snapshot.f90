!=======================================================
! This is snapshot ver.5 (01/24/2019)
!=======================================================

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
           eachpsi(mtheta(ipsi)/2*mzetamax),wt,inv_dgq(0:mpsi),&
           weight,adum(0:mpsi),bdum(mbin_u,mbin_psi)
  real(wp),dimension(0:mpsi) :: marker,fflows,dflows,ftem,dtem
  real(wp),dimension(mbin_u,mbin_psi) :: ubin,dubin,pbin,dpbin,ebin,debin
  real(wp) ubin_max,pbin_max,ebin_max
!  real(wp),dimension(mzeta,0:mpsi,0:mthetamax) :: ftotal,tpara,tperp

  real(wp) intpl_dphi, intpl_dTe, intpl_dTi, intpl_dni, intpl_dne, dn_adiabatic, intpl_dphi_NZ
  real(wp) tmp_full_dphi(0:1,mgrid), tmp_ne_adiabatic(0:1,mgrid), intp_val  
  integer gq_jump, a_jump
  character(len=32) cdum,date,time

  real(wp),dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj

! write particle data for re-run
!  call restart_write
!ZHIXIN use flag irewrite to output RESTART
! if(mod(istep,40) .eq. 0) call restart_write
  if(mod(istep,irewrite) .eq. 0) call restart_write

    gq_jump = 2
    a_jump = 1

  call Cal_dni_dTi_from_each_ion_weight_component
  call Cal_dne_dTe_from_each_electron_weight_component

! number of poloidal grid
  jm=mtheta(ipsi)/gq_jump

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
  inv_dgq=1.0/deltat
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

  !====================================================================
  ! MGYOO(01/09/2019) : To measure FLA(Field Line Average) component 
  ! This FLA might be important in the case of external magnetic perturbation

  do i=0,mpsi
    do j=0,mtheta(i)
      ij=igrid(i) + j
      tmp_full_dphi(:,ij) = phi(:,ij) +real(ZF_SEPARATION)*phi00(i)        
    enddo
  enddo

  if((myrank_partd==nproc_partd-1) .and. (mod(istep,msnap) .eq. 0)) then 
    call Do_FieldLineAverage_by_using_gPPI(tmp_full_dphi, FLA_dphi)    
    call Do_FieldLineAverage_by_using_gPPI(tmp_ne_adiabatic, FLA_ne_adiabatic)
  endif
  !=========================================================================

  if(ifelectron==0) then

!     call esigetb(zion(1,:),zion(2,:),mi)   ! get B
     call esigetb(zion,mi,7)   ! get B-->BB

     do m=1,mi 
        weight=zion0(6,m)
        r=zion(1,m)        
        ibin=max(1,min(mbin_psi,1+int((r-a0)*dela)))
        ip=max(0,min(mpsi,int((zion(7,m)-s0)*delr+0.5)))
                           ! zion(2,*) and zion(3,*) are within [0,2*pi)
        jt=max(0,min(mtheta(ip),int(zion(2,m)*inv_dgq(ip)+0.5)))
 
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
        jt=max(0,min(mtheta(ip),int(zeon(2,m)*inv_dgq(ip)+0.5)))
 
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
  call MPI_REDUCE(ubin,bdum,icount,mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierror)
  if(mype==numberpe-1)ubin=bdum
  call MPI_REDUCE(dubin,bdum,icount,mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierror)
  if(mype==numberpe-1)dubin=bdum
  call MPI_REDUCE(pbin,bdum,icount,mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierror)
  if(mype==numberpe-1)pbin=bdum
  call MPI_REDUCE(dpbin,bdum,icount,mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierror)
  if(mype==numberpe-1)dpbin=bdum
  call MPI_REDUCE(ebin,bdum,icount,mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierror)
  if(mype==numberpe-1)ebin=bdum
  call MPI_REDUCE(debin,bdum,icount,mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierror)
  if(mype==numberpe-1)debin=bdum
  call MPI_REDUCE(marker,adum,mpsi+1,mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierror)
  if(mype==numberpe-1)marker=adum
  call MPI_REDUCE(fflows,adum,mpsi+1,mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierror)
  if(mype==numberpe-1)fflows=adum
  call MPI_REDUCE(dflows,adum,mpsi+1,mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierror)
  if(mype==numberpe-1)dflows=adum
  call MPI_REDUCE(ftem,adum,mpsi+1,mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierror)
  if(mype==numberpe-1)ftem=adum
  call MPI_REDUCE(dtem,adum,mpsi+1,mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierror)
  if(mype==numberpe-1)dtem=adum

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
  
  if(mype == numberpe-1)then !mype=0 must also have myrank_toroidal=0
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
     write(snapout,101) mbin_u,mbin_psi,mpsi/a_jump,jm,mzetamax,mpsi/4   ! mpsi/2 specifies  
                    ! how many number of data are written out; don't confuse with ipsi! 
     write(snapout,102)  1e3*Tem_00*(vth_grid(ipsi)**2), 1e3*Tem_00*Te_grid(ipsi) ! Ti and Te at reference position in eV unit


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
     do i=0,mpsi,a_jump
        write(snapout,102)  psi_grid(i)  !  (a0+deltar*real(i))/gyroradius  ????
     enddo
     write(snapout,102)(marker(i),i=0,mpsi,a_jump)
!     write(snapout,102)(phi00(i)/vthc**2,i=1,mpsi)
     do i=0,mpsi,a_jump
        r=psi_grid(i)       !   r=(a0+deltar*real(i))/a  ????
!        r=r_grid(i)        !   use this minor radiu      
        q=q_grid(i) 
        write(snapout,102)q/r
     enddo
     write(snapout,102)(phip00(i)/gyroradius,i=0,mpsi,a_jump)   ! ???
     write(snapout,102)(fflows(i),i=0,mpsi,a_jump)
     write(snapout,102)(dflows(i),i=0,mpsi,a_jump)
     write(snapout,102)(ftem(i),i=0,mpsi,a_jump)
     write(snapout,102)(dtem(i),i=0,mpsi,a_jump)



     !==========================================================
     ! MGYOO (01/09/2019): From here, most outputs are re-scaled to be SI unit.
     ! EX) density in (#/m^3), temperature in (eV), potential in (V)
     !=========================================================
     ! Equilibrium density and temperatures
     write(snapout,101) 1
     write(snapout,102)(1e6*den_grid(i),i=0,mpsi,a_jump)    ! density in (#/m^3) unit  
     write(snapout,101) 1
     write(snapout,102)(1e3*Tem_00*(vth_grid(i)**2),i=0,mpsi,a_jump)      ! Ti in (eV) unit
     write(snapout,101) 1
     write(snapout,102)(1e3*Tem_00*Te_grid(i),i=0,mpsi,a_jump)          ! Te in same unit (normalized unit is Tem_00)

     ! Zonal component of potential, density and temperatures
     write(snapout,101) 1
     write(snapout,102)(((E_00*c)/(charg0*1e8))*phi00(i),i=0,mpsi,a_jump)   ! potential in (V) unit
     write(snapout,101) 1
     write(snapout,102)(((E_00*c)/(charg0*1e8))*phip00(i),i=0,mpsi,a_jump)  ! potential derivative in (V/ psi) unit
     write(snapout,101) 1
     write(snapout,102)(1e6*den_grid(i)*zonali(i),i=0,mpsi,a_jump)   ! zonal ion density in (#/m^3) unit  

     write(snapout,101) 1
     write(snapout,102) (1e6*den_grid(i)*zonale(i),i=0,mpsi,a_jump)   ! zonal electron density in (#/m^3) unit  
     

     write(snapout,101) 1
     write(snapout,102)(1e3*Tem_00*zonal_dTi(i),i=0,mpsi,a_jump)  ! zonal_dTi is (eV) unit
     write(snapout,101) 1
     write(snapout,102)(1e3*Tem_00*zonal_dTe(i),i=0,mpsi,a_jump)   !zonal_dTe is (eV) unit

! poloidal cross section: marker,density,potential,parallel flow and temperature
     write(snapout,101) 1
     do i=0,mpsi,a_jump
        do j=0,jm
           pdum=psi_grid(i)  !wx
           tdum=2.0*pi*real(j)/real(jm)
           aaa(1)=pdum
           the(1)=tdum
           call esigetrzb(rj,zj,Bj,aaa,the,1)  
           write(snapout,102) rj(1)
        enddo
     enddo

     write(snapout,101) 1
     do i=0,mpsi,a_jump
        do j=0,jm
           pdum=psi_grid(i)      !wx
           tdum=2.0*pi*real(j)/real(jm)
           aaa(1)=pdum
           the(1)=tdum
           call esigetrzb(rj,zj,Bj,aaa,the,1)
           write(snapout,102) zj(1) 
        enddo
     enddo

     write(snapout,101) 1
     do i=0,mpsi,a_jump
        do j=0,jm
           pdum=psi_grid(i)    
           tdum=2.0*pi*real(j)/real(jm)
           write(snapout,102) pdum
        enddo
     enddo 
     
     write(snapout,101) 1
     do i=0,mpsi,a_jump
        do j=0,jm
           pdum=psi_grid(i)    
           tdum=2.0*pi*real(j)/real(jm)
           write(snapout,102) tdum
        enddo
     enddo
     
     
     !==============================================
     ! MGYOO (11/13/2018) : record potential, density, temperature fluctuation at zeta=0 plane
     !==============================================
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
          tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
          jt=int(tdum);  wt=tdum-jt
          intpl_dphi = ((1.0-wt)*phi(1,igrid(i)+jt)+  wt*phi(1,igrid(i)+jt+1)) +real(ZF_SEPARATION)*phi00(i)
          write(snapout,102) ((E_00*c)/(charg0*1e8))*intpl_dphi ! total delta-phi (including zonal phi) in (V) unit
      enddo
    enddo

    ! ion density
    ! from tot_dni with no enforced charge-conservation
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intpl_dni = ((1.0-wt)*tot_dni(1,igrid(i)+jt) + wt*tot_dni(1,igrid(i)+jt+1)) 
            write(snapout,102) 1e6*intpl_dni ! ion total perturbed density (dni) in (#/m^3) unit            
        enddo
    enddo

    ! CAUTION!! Following routine is just for test
    ! from densityi directly
    ! if(its_jacobi=1) This is the result enforced by charge-conservation
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intpl_dni = den_grid(i)*((1.0-wt)*densityi(1,igrid(i)+jt) + wt*densityi(1,igrid(i)+jt+1) + real(ZF_SEPARATION)*zonali(i)) 
            write(snapout,102) 1e6*intpl_dni ! ion total perturbed density (dni) in (#/m^3) unit            
        enddo
    enddo

    ! ion temperature fluctuation 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
          tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
          jt=int(tdum);  wt=tdum-jt
          intpl_dTi = (vth_grid(i)**2)*((1.0-wt)*phi_tmp(igrid(i)+jt)+  wt*phi_tmp(igrid(i)+jt+1)) +real(ZF_SEPARATION)*(zonal_dTi(i))
          write(snapout,102) 1e3*Tem_00*intpl_dTi ! dTi in (eV) unit
        enddo
    enddo

    ! electron density
    ! from tot_dne with no enforced charge-conservation    
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intpl_dne = ((1.0-wt)*tot_dne(1,igrid(i)+jt) + wt*tot_dne(1,igrid(i)+jt+1)) 
            write(snapout,102) 1e6*intpl_dne ! electron total perturbed density (dne) in (#/m^3) unit            
        enddo
    enddo

    ! CAUTION!! Following routine is just for test
    ! from densitye directly
    ! if(its_jacobi=1) This is the result enforced by charge-conservation
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intpl_dne = den_grid(i)*((1.0-wt)*densitye(1,igrid(i)+jt) + wt*densitye(1,igrid(i)+jt+1) + real(ZF_SEPARATION)*zonale(i)) 
            write(snapout,102) 1e6*intpl_dne ! electron total perturbed density (dne) in (#/m^3) unit            
        enddo
    enddo

    ! electron temperature fluctuation 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
          tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
          jt=int(tdum);  wt=tdum-jt
          intpl_dTe = Te_grid(i)*(((1.0-wt)*heatflux(igrid(i)+jt)+  wt*heatflux(igrid(i)+jt+1))) +real(ZF_SEPARATION)*(zonal_dTe(i))
          write(snapout,102) 1e3*Tem_00*intpl_dTe  ! dTi in (eV) unit
        enddo
    enddo



  !====================================================================
  ! MGYOO(12/03/2018) : To measure FLA(Field Line Average) component 
    ! FLA_psi0 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intp_val = ((1.0-wt)*FLA_psi0(1,igrid(i)+jt)+  wt*FLA_psi0(1,igrid(i)+jt+1))
            write(snapout,102) intp_val
        enddo
    enddo

    ! FLA_Ti0 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intp_val = ((1.0-wt)*FLA_Ti0(1,igrid(i)+jt)+  wt*FLA_Ti0(1,igrid(i)+jt+1))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo
    
     ! FLA_Te0 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intp_val = ((1.0-wt)*FLA_Te0(1,igrid(i)+jt)+  wt*FLA_Te0(1,igrid(i)+jt+1))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo   

    ! FLA_n0 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intp_val = ((1.0-wt)*FLA_n0(1,igrid(i)+jt)+  wt*FLA_n0(1,igrid(i)+jt+1))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo  
    
    ! FLA_dphi 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intp_val = ((1.0-wt)*FLA_dphi(1,igrid(i)+jt)+  wt*FLA_dphi(1,igrid(i)+jt+1))
            write(snapout,102) ((E_00*c)/(charg0*1e8))*intp_val  ! potential in (V) unit
        enddo
    enddo    

    ! FLA_dTi 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intp_val = ((1.0-wt)*FLA_dTi(1,igrid(i)+jt)+  wt*FLA_dTi(1,igrid(i)+jt+1))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo

    ! FLA_dni 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intp_val = ((1.0-wt)*FLA_dni(1,igrid(i)+jt)+  wt*FLA_dni(1,igrid(i)+jt+1))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo    

    ! FLA_dTe
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intp_val = ((1.0-wt)*FLA_dTe(1,igrid(i)+jt)+  wt*FLA_dTe(1,igrid(i)+jt+1))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo

    ! FLA_dne 
    write(snapout,101) 1   
    do i=0,mpsi,a_jump
      do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intp_val = ((1.0-wt)*(FLA_dne(1,igrid(i)+jt))+  wt*(FLA_dne(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
      enddo
    enddo

    ! FLA_ne_adiabatic 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt
            intp_val = ((1.0-wt)*(FLA_ne_adiabatic(1,igrid(i)+jt))+  wt*(FLA_ne_adiabatic(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo

    !==============================================================================
    ! MGYOO (01/16/2019)
    ! From here, dn and dT from each weight component for ions and electrons

    ! dni each components 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dni_init(1,igrid(i)+jt))+  wt*(dni_init(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dni_vpara(1,igrid(i)+jt))+  wt*(dni_vpara(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dni_vperp(1,igrid(i)+jt))+  wt*(dni_vperp(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo  
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dni_accel(1,igrid(i)+jt))+  wt*(dni_accel(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo  
    
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dni_fdamp(1,igrid(i)+jt))+  wt*(dni_fdamp(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo  
    
    ! dTi each components
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dTi_init(1,igrid(i)+jt))+  wt*(dTi_init(1,igrid(i)+jt+1)))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dTi_vpara(1,igrid(i)+jt))+  wt*(dTi_vpara(1,igrid(i)+jt+1)))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dTi_vperp(1,igrid(i)+jt))+  wt*(dTi_vperp(1,igrid(i)+jt+1)))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dTi_accel(1,igrid(i)+jt))+  wt*(dTi_accel(1,igrid(i)+jt+1)))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo
    
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dTi_fdamp(1,igrid(i)+jt))+  wt*(dTi_fdamp(1,igrid(i)+jt+1)))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo


    ! dne each components 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dne_init(1,igrid(i)+jt))+  wt*(dne_init(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dne_vpara(1,igrid(i)+jt))+  wt*(dne_vpara(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dne_vperp(1,igrid(i)+jt))+  wt*(dne_vperp(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo  
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dne_accel(1,igrid(i)+jt))+  wt*(dne_accel(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo  
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dne_adiab(1,igrid(i)+jt))+  wt*(dne_adiab(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo    
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dne_fdamp(1,igrid(i)+jt))+  wt*(dne_fdamp(1,igrid(i)+jt+1)))
            write(snapout,102) 1e6*intp_val ! density in (#/m^3) unit
        enddo
    enddo    
    ! dTe each components
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dTe_init(1,igrid(i)+jt))+  wt*(dTe_init(1,igrid(i)+jt+1)))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dTe_vpara(1,igrid(i)+jt))+  wt*(dTe_vpara(1,igrid(i)+jt+1)))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo 
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dTe_vperp(1,igrid(i)+jt))+  wt*(dTe_vperp(1,igrid(i)+jt+1)))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dTe_accel(1,igrid(i)+jt))+  wt*(dTe_accel(1,igrid(i)+jt+1)))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo    
    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dTe_adiab(1,igrid(i)+jt))+  wt*(dTe_adiab(1,igrid(i)+jt+1)))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo 

    write(snapout,101) 1
    do i=0,mpsi,a_jump
        do j=0,jm
            tdum = modulo(2.0*pi*real(j)/real(jm) - qtinv(i)*zetamax, 2.0*pi) * inv_dgq(i)
            jt=int(tdum);  wt=tdum-jt                        
            intp_val = ((1.0-wt)*(dTe_fdamp(1,igrid(i)+jt))+  wt*(dTe_fdamp(1,igrid(i)+jt+1)))
            write(snapout,102) 1e3*Tem_00*intp_val ! temperature in (eV) unit
        enddo
    enddo 
  ! MGYOO (01/16/2019)
    !==============================================================================



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
     do i=0,mpsi,a_jump
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
     if(mype == numberpe-1)write(stdout,*) 'Program ends at DATE=', date, 'TIME=', time
  endif
  
101 format(i6)
!  102 format(e10.4)

  
!  105 format(e10.4,4x,6hflux2x,4x,i4)
!  106 format(e10.4,4x,6hflux2z,4x,i4)
!  107 format(e10.4,2x,i4,2x,4hpsi=,e10.4,2x,4hthe=,e10.4,2x,2hR=,e10.4,2x,2hz=,e10.4)

102 format(ES16.8)

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

  if(istep==0 .and. myrank_partd==nproc_partd-1) then
    call Print_FLA_Denominator
  endif

!   if(mype==0) print*, 'Write DENs + PHIs at ', 'istep=', istep,  'irk=', irk

end subroutine snapshot


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
      write(181206, "((ES14.6,1x))") FLA_Denominator(1,ij)
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

  do m=1,mi

    wgf1=wzion(m)
    wgf0=1.0-wgf1

    !! Caution!! Is backu(m) here omega or real velocity?? Check it! Refer pushion routine
    energy = 0.5*(zion(4,m)*BB(m)-backu(m))**2.0 + ((zion(6,m)**2.0)*BB(m)) 

    
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
  ! Convert dTe to get a correct form of dTe by considering dne.
  ! buff = Sigma(weight*(2/3*energy-T0))*markere = (1/n0)*Sigma(weight*(2/3*energy-T0))/node_Vol 
  ! dTe = (1/(n0+dn))*Sigma(weight*(2/3*energy-T0))/node_Vol = n0/(n0+dn)*buff
  do i=0,mpsi
    do j=0,mtheta(i)
      ij = igrid(i) + j
      tmp1 =  den_grid(i)/(den_grid(i)+tot_dni(1,ij))
      dTi_init(1,ij)  = tmp1*dTi_init(1,ij)
      dTi_vpara(1,ij) = tmp1*dTi_vpara(1,ij)
      dTi_vperp(1,ij) = tmp1*dTi_vperp(1,ij)
      dTi_accel(1,ij) = tmp1*dTi_accel(1,ij)    
      dTi_fdamp(1,ij) = tmp1*dTi_fdamp(1,ij)  
    enddo
  enddo
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


  do m=1,me

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

  ! Convert dTe to get a correct form of dTe by considering dne.
  ! buff = Sigma(weight*(2/3*energy-T0))*markere = (1/n0)*Sigma(weight*(2/3*energy-T0))/node_Vol 
  ! dTe = (1/(n0+dn))*Sigma(weight*(2/3*energy-T0))/node_Vol = n0/(n0+dn)*buff
  do i=0,mpsi
    do j=0,mtheta(i)
      ij = igrid(i) + j
      tmp1 =  den_grid(i)/(den_grid(i)+tot_dne(1,ij))
      dTe_init(1,ij)  = tmp1*dTe_init(1,ij)
      dTe_vpara(1,ij) = tmp1*dTe_vpara(1,ij)
      dTe_vperp(1,ij) = tmp1*dTe_vperp(1,ij)
      dTe_accel(1,ij) = tmp1*dTe_accel(1,ij)    
      dTe_adiab(1,ij) = tmp1*dTe_adiab(1,ij)    
      dTe_fdamp(1,ij) = tmp1*dTe_fdamp(1,ij)    
    enddo
  enddo
  !==============================================================================


end subroutine Cal_dne_dTe_from_each_electron_weight_component
