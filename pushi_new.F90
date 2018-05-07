subroutine pushi
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface

  implicit none

!  integer,parameter :: nrbin=5       !wx
  integer i,im,ip,j,jt,kz,m,ii,j11,j10,j01,j00,jtheta,ierror,larmor,ij,ipjt
  real(wp) rdum,dtime,r,grada,vpara,energy,kappa,vdr,wdot,wdrift,&
       vdrenergy,delr,delt(0:mpsi),delz,pi2_inv,adum(0:mpsi),&
       xdum,ydum,xdot,ydot,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,&
       epara,psimax,psimin,pdum,wdrive,ainv,psitmp,thetatmp,&
       zetatmp,rhoi,e1,e2,e3,wpi(3,mi),rdel, &
       dgr(mi),da(mi),dgq(mi),dgf(mi),mu, &
       Er_particle(mi),dela,fip,weqflow,T_inverse,& 
       B_inverse,energy1,vdq,wKH,wgradB,omega_t,u_vpara,  &
       v2, avw,   &
! following working arrays are defined to use OMP reduction:
       adum_00(0:mpsi),adum_tur(0:mpsi),adum_tot(0:mpsi) 

!  integer ip_index(mimax)   !wx

  integer gcmotion

  real(wp) tt0,tt1,cput0,cput1,rtc

  real(wp), external :: D01AHF,dsda,SMOOTHSTEP
  integer minpts,ifail
  real(wp) acc,error_8
  real(wp) E_a,E_t,E_f,timef

  real(wp) ddT_da_particle, w_anti_relax, muB, instep,outstep,fdamp,fdamp0

  real(wp) wfield(3,mgrid,0:1)   ! weight-field W_0+v_//*W_1+v^2*W_2

  delr=1.0_wp/deltar
  delt=1.0_wp/deltat
  delz=1.0_wp/deltaz
  pi2_inv=0.5_wp/pi
  psimax=psi_grid(mpsi)    
  psimin=psi_grid(0)       

  fdamp0=1.0_wp/tstep*5.0d-4
  if(nonlinear<0.5_wp) fdamp0=0.0_wp
!  fdamp0=0.0_wp

  if(irk==1)then
! 1st step of Runge-Kutta method
     dtime=0.5_wp*tstep

! use (x,y) coordinates near magnetic axis
!$omp parallel do private(m,xdum,ydum)
     do m=1,mi
        zion0(1:5,m)=zion(1:5,m)
        zion0(7,m)=zion(7,m)
#ifdef _TREAT_MAGNETIC_AXIS
        if(sw(m).eq.1) then
           xdum=zion(1,m)*cos(zion(2,m))
           ydum=zion(1,m)*sin(zion(2,m))
           zion0(1,m)=xdum
           zion0(2,m)=ydum
        end if
#endif
     enddo

! 2nd step of Runge-Kutta method
  else
     dtime=tstep
  endif

!  mu(1:mi)=zion(6,1:mi)**2

! gather e_field using 4-point gyro-averaging, sorting in poloidal angle
! The following line is an OpenMP directive for loop-level parallelism
! on shared memory machines (see http://www.openmp.org).
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

!    prepare equilibrium Er=dPhi/da at particles 



!  ------------------------------------------------------------

!  apply equilibrium Er in the middle insteady at the beginning

!  if(istep.lt.2800) then
!  if(istep.lt.2950) then
!     timef=0.0_wp
!  else
!     timef=1.0_wp-exp(-(real(istep-2950)/5.0d2)**4)
!  end if
!  Er401_bp=Er401_bp0*timef
!  Er_grid=Er_grid0*timef

!  ------------------------------------------------------------



  dela=1.0_wp/Ydr

! if(istep.ge.4000) then
!$omp parallel do private(m,r,ip,fip)
  do m=1,mi
     r=zion(1,m)-a0      
     r=r*dela
     ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
     fip=r-real(ip,wp)
     Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
  enddo
! else
!   Er_particle=0.0
! end if

!  tt0=rtc()
!  call cpu_time(cput0)

!!  m=gcmotion(dgr,da,dgq,dgf,zion(4,:),zion(1,:),zion(2,:),mu,wpi(1,:),wpi(2,:),wpi(3,:),Er_particle,nonlinear,mi)

  m=gcmotion(dgr,da,dgq,dgf,zion,wpi,Er_particle,nonlinear,mi)

!  tt1=rtc()
!  call cpu_time(cput1)

!  call n_T_u_group  ! prepare n, T, u_//, dlnn/da, dlnT/da, etc. 
                     ! moved to collision.f90 -- caution: in this case, BB must be obtained by call esigetb etc.!!!

!     OUTPUT: dgr=drho/dt;       da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt; 
!             dgf=dzeta/dt;      ghq(:)=da/dt|_EXB;     Ta(:)=dtheta/dt|_EXB; 
!             BB(:)=B;           Ba(:)=dB/da;           Bq(:)=dB/dtheta;  
!             F(:) =F;           Fa(:)=dF/da;           gha(:)=1/|grad a|
!             no change for zion(1,2,4,:), wpi(1,2,:) and mu;
!             wpi(3,:)=dphi/dzeta (no change);          Pa(:)=R (major radius)

if(nonlinear>0.5_wp) then        ! for nonlinear run
  call weight_field(wfield)      ! calculate weight-field uisng zion0 for irk=1 and zion at half-step for irk=2
                                 ! used BB is obtained from gcmotion, and calculte its own ion locations  
                                 ! which are used to calculate wfield in phase space grids and gether wfield at ions (next)
end if

! update GC position
! Another OpenMP parallel loop...
!$omp parallel do private(m,mu,vpara,energy,kappa,wdrift,weqflow,vdr,wdrive,wdot,xdum,ydum,T_inverse,energy1,B_inverse,u_vpara,vdq,wgradB,wKH, &
!$omp& omega_t,r,ip,fip,ddT_da_particle,w_anti_relax,muB,instep,outstep,fdamp,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,v2,avw)
  do m=1,mi
!     dbdz=0.0   for axisymmetric system   !wx

     mu=zion(6,m)*zion(6,m)    ! =v_perp^2/2B
     muB=mu*BB(m)              ! =v_perp^2/2
     T_inverse=1.0_wp/backT(m)
     B_inverse=1.0_wp/BB(m)
     vpara=zion(4,m)*BB(m)     ! BB(m)=B
     omega_t=backu(m)
     backu(m)=backu(m)*F(m)*B_inverse         ! =u_//=I*omega_t/B
     u_vpara=backu(m)/vpara                   ! =(u_//)/(v_//)
!     energy=0.5_wp*(vpara-backu(m))*(vpara-backu(m))+mu*BB(m)
     energy=0.5_wp*(vpara-backu(m))*(vpara-backu(m))+muB
     energy1=energy-backu(m)*backu(m)*0.5_wp  ! =(v_//-u)^2/2+mu*B-u^2/2 -- see note

     v2=vpara*vpara+muB*2.0_wp        ! v^2

     kappa=(min(umax*umax,energy1*T_inverse)-1.5_wp)*(-dlnTda(m))+(-dlnnda(m))+gh(m) !~-dlnf_SM/da
                                                                                  !gh=dln<A>/da with A=exp(u^2/2T)
     r=zion(7,m)-s0
     r=r*delr
     ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]; using coordinate s to find particle radial location in grids
     fip=(zion(1,m)-psi_grid(ip))/(psi_grid(ip+1)-psi_grid(ip))
     ddT_da_particle=fip*ddT_da(ip+1)+(1.0_wp-fip)*ddT_da(ip) ! d<delta_Ti>/da at particle's position

! ExB drift in radial direction for w-dot and flux diagnostics
     vdr=ghq(m)             !  da/at|_EXB
     vdq=Ta(m)              !  dtheta/dt|_EXB

     wdrive=vdr*kappa

! new term added to maintain constant grad(Ti) drive
     w_anti_relax=(min(umax*umax,energy*T_inverse)-1.5_wp)*ddT_da_particle*T_inverse
     w_anti_relax=vdr*w_anti_relax

     if(zion(1,m).le.lmax.or.zion(1,m).ge.rmin) then
        instep =smoothstep(lmin,lmax,zion(1,m))
        outstep=smoothstep(rmin,rmax,zion(1,m))
!        fip=smoothstep(0.20_wp,0.25_wp,zion(1,m))*(1.0_wp-smoothstep(0.75_wp,0.80_wp,zion(1,m)))
        fip=instep*(1.0_wp-outstep)
        wdrive=wdrive*fip
        w_anti_relax=w_anti_relax*fip
        fdamp=(1.0_wp-fip)*fdamp0
     else
        fdamp=0.0_wp
     end if

     weqflow=-vdr*Er_particle(m)*T_inverse*(1.0_wp-u_vpara)   ! additional term due to E_0XB drift 

     wdrift=-(da(m)*wpi(1,m)+dgq(m)*wpi(2,m)+dgf(m)*wpi(3,m))*T_inverse*(1.0_wp-u_vpara)  !wx  ??? x,y

!     wKH=-(vdr*(duda(m)/omega_t+Fa(m)/F(m)-Ba(m)*B_inverse)     &
!           -vdq*Bq(m)*B_inverse)*backu(m)*vpara*T_inverse       ! major Kelvin-Helmholtz drive
!        the above expression for wKH cause problem when omega_t=0. so replaced by:
     wKH=-(vdr*(duda(m)*F(m)+Fa(m)*omega_t-Ba(m)*backu(m)) &
           -vdq*Bq(m)*backu(m))*B_inverse*vpara*T_inverse       ! major Kelvin-Helmholtz drive

     wgradB=(vdr*Ba(m)+vdq*Bq(m))*u_vpara*T_inverse*mu   ! interesting term propotional to grad_B


if(nonlinear>0.5_wp) then     ! for nonlinear run
! calculate average weight <w> at particles using wieght-field wfield
     wz1=wzion(m)
     wz0=1.0_wp-wz1

     ij=jtion0(1,m)
     wp0=1.0_wp-wpion(1,m)

     wt00=1.0_wp-wtion0(1,m)
     e1=wp0*wt00*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=wp0*wt00*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=wp0*wt00*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt10=1.0_wp-wt00
     e1=e1+wp0*wt10*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp0*wt10*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp0*wt10*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=jtion1(1,m)
     wp1=1.0_wp-wp0
     wt01=1.0_wp-wtion1(1,m)
     e1=e1+wp1*wt01*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt01*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt01*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt11=1.0_wp-wt01
     e1=e1+wp1*wt11*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt11*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt11*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     avw=e1+e2*vpara+e3*v2
end if


!     wdot=(zion0(6,m)-nonlinear*zion(5,m))*(wdrive+wdrift+weqflow+wKH+wgradB+w_anti_relax)

!     wdot=(zion0(6,m)-nonlinear*zion(5,m))*(wdrive+wdrift+weqflow+wKH+wgradB)

     wdot=(zion0(6,m)-nonlinear*avw)*(wdrive+wdrift+weqflow+wKH+wgradB)

! if(mype==0.and.m==mi/2) print*, 'GGGGG', wdrive, w_anti_relax

     zion(1,m) = zion0(1,m)+dtime*da(m)
     zion(2,m) = zion0(2,m)+dtime*dgq(m)
     zion(3,m) = zion0(3,m)+dtime*dgf(m)
     zion(4,m) = zion0(4,m)+dtime*dgr(m)
     zion(5,m) = zion0(5,m)+dtime*(wdot-fdamp*zion(5,m))
!     zion(7,m) = zion0(7,m)+dtime*sqrt(Tc*T_inverse)*da(m)
!     zion(7,m) = zion0(7,m)+dtime*(Tc*T_inverse)**0.4_wp*da(m)
     zion(7,m) = zion0(7,m)+dtime*(Tc*T_inverse)**dsda_power*da(m)

#ifdef _TREAT_MAGNETIC_AXIS
     if(sw(m).eq.1) then    ! convert to (a, theta) coordinates
        xdum=sqrt(zion(1,m)**2+zion(2,m)**2)
        ydum=sign(1.0_wp,zion(2,m))*acos(zion(1,m)/xdum)
        zion(1,m)=xdum
        zion(2,m)=ydum
     end if
#endif
     
! theta and zeta normalize to [0,2*pi), modulo is slower than hand coded procedure
!     zion(2,m)=zion(2,m)*pi2_inv+100.0 !period of 1
!     zion(2,m)=2.0*pi*(zion(2,m)-aint(zion(2,m))) ![0,2*pi)
!     zion(3,m)=zion(3,m)*pi2_inv+100.0
!     zion(3,m)=2.0*pi*(zion(3,m)-aint(zion(3,m)))

! use f90 modulo to normalize theta and zeta to [0,2*pi)
! If P > 0, then the function modulo(A,P) determines the positive difference between A and then
! next lowest multiple of P. If P < 0, then the function modulo(A,P) determines the negative
! difference between A and then next highest multiple of P.
     zion(2,m)=modulo(zion(2,m),twopi)
     zion(3,m)=modulo(zion(3,m),twopi)

! store GC information for flux measurements
     wpi(1,m)=vdr       ! at present, wpi(1,*)=da/dt|_EXB 
     wpi(2,m)=energy    ! =(v_para-u//)^2/2+v_perp^2/2
     wpi(3,m)=muB       ! =v_perp^2/2, need for calculate distribution function in get_fluxes

  enddo
      
 if(irk==2)then
! out of boundary particle
!$omp parallel do private(m)
    do m=1,mi
       if(zion(1,m) > psimax)then
          zion(1,m)=zion0(1,m)
          zion(2,m)=2.0_wp*pi-zion0(2,m)
          zion(3,m)=zion0(3,m)
          zion(4,m)=zion0(4,m)
          zion(5,m)=zion0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
!          zion(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
          zion(7,m)=zion0(7,m)
          
       elseif(zion(1,m) < psimin)then
          zion(1,m)=zion0(1,m)
          zion(2,m)=2.0_wp*pi-zion0(2,m)
          zion(3,m)=zion0(3,m)
          zion(4,m)=zion0(4,m) 
          zion(5,m)=zion0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
!          zion(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
          zion(7,m)=zion0(7,m)
       endif
    enddo







! ********************************
!if(mype==2) then
!!  error_8=1.0d-5
!!  IFAIL=1
!!  rdum=D01AHF(0.0d0,zion(1,mi/2),error_8,minpts,acc,dsda,-1,IFAIL)
!!  print*, '##################', zion(7,mi/2),'  ', rdum
!   write(stdout,*)'######## gcmotion: wall =',(tt1-tt0),'  cpu =',(cput1-cput0)
!end if
! ********************************

 endif


 if(idiag==0)then
! fluxes diagnose
!    wpi(1,1:mi)=wpi(1,1:mi)*gha(1:mi)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
!$omp parallel do private(m)
    do m=1,mi
       wpi(1,m)=wpi(1,m)*gha(m)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
    end do

    rmarker=0.0_wp
    eflux=0.0_wp
    pfluxi=0.0_wp
    dflowi=0.0_wp
    entropyi=0.0_wp
    fluxpsi=0.0_wp
    heatflux=0.0_wp
    En_00=0.0_wp
    En_tur=0.0_wp
    En_tot=0.0_wp
    do m=1,mi
       r=zion0(7,m)         ! old radial location
       ip=max(0,min(mpsi,int((r-s0)*delr+0.5_wp)))    ! ip [0,mpsi]
       vdrenergy=wpi(1,m)*(wpi(2,m)-2.5*vth_grid(ip)**2)*zion0(5,m) 
                ! note: wpi(1,m)=dr/dt (EXB part) 

! radial location in diagnostic bin
       i=max(1,min(mflux,1+int((r-s0)*real(mflux)/(s1-s0))))

! marker,energy,particle,momentum fluxes,parallel flows,entropy and kinetic energy
       rmarker(i)=rmarker(i)+zion0(6,m)
       eflux(i)=eflux(i)+vdrenergy
       pfluxi=pfluxi+wpi(1,m)*zion0(5,m)
       dflowi=dflowi+BB(m)*zion0(4,m)*zion0(5,m)
       entropyi=entropyi+zion0(5,m)*zion0(5,m)
    enddo

   if(ZF_SEPARATION==0) then   

!     note: zonal flow component <phi> and phip00=d<phi>/da are calculated in field.f90

    adum=0.0_wp
    adum_00=0.0_wp
    adum_tur=0.0_wp
    adum_tot=0.0_wp

!$omp parallel private(i,j,ij,E_a,E_t,E_f) reduction(+:adum,adum_00,adum_tur,adum_tot)
! Note from OMP compiler: A variable in the REDUCTION clause or any component or element thereof must not have the POINTER 
! attribute; nor be a pointee, an allocatable object, an assumed-size array, or a THREADLOCAL common block variable.
!$omp do
    do i=0,mpsi
       do j=1,mtheta(i)
          ij=igrid(i)+j
          adum(i)=adum(i)+dtemper(ij)*(phi(1,ij)-phi00(i))**2  ! surface averaging for delta-phi

          adum_00(i)=adum_00(i)+dtemper(ij)*phip00(i)**2*gij(1,ij)  ! surface average

          E_a=evector(1,1,ij)+Er_grid(i)           ! total radial component   ! d\phi/da
          E_t=evector(2,1,ij)                      ! d\delta\phi/d\theta
          E_f=evector(3,1,ij)                      ! d\delta\phi/d\zeta
          adum_tot(i)=adum_tot(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                                              +E_t**2*gij(2,ij)  &
                                              +E_f**2*gij(3,ij)  &
                                              +2.0_wp*E_a*E_t*gij(4,ij))

          E_a=E_a-phip00(i)-Er_grid(i) ! turbulence component only    ! d\delta\phi/da
          adum_tur(i)=adum_tur(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                                              +E_t**2*gij(2,ij)  &
                                              +E_f**2*gij(3,ij)  &
                                              +2.0_wp*E_a*E_t*gij(4,ij))
       enddo
    enddo
!$omp end parallel

   else  ! ZF_SEPARATION=1

    adum=0.0_wp
    adum_00=0.0_wp
    adum_tur=0.0_wp
    adum_tot=0.0_wp

! DEBUG 
!        if(mype==100) print*, 'Dirichlet2', phi(1,igrid(0)), phi(1,igrid(0)+mtheta(0)/2),phi(1,igrid(0)+mtheta(0)-1), &
!                               phi(1,igrid(mpsi)), phi(1,igrid(mpsi)+mtheta(mpsi)/2),phi(1,igrid(mpsi)+mtheta(mpsi)-1)

!$omp parallel private(i,j,ij,E_a,E_t,E_f) reduction(+:adum,adum_00,adum_tur,adum_tot)
!$omp do
    do i=0,mpsi
       do j=1,mtheta(i)
          ij=igrid(i)+j
!          fluxpsi(i)=fluxpsi(i)+phi(1,ij)*phi(1,ij)               ! what average?
          adum(i)=adum(i)+dtemper(ij)*phi(1,ij)*phi(1,ij)  ! surface averaging for delta-phi

          adum_00(i)=adum_00(i)+dtemper(ij)*phip00(i)**2*gij(1,ij)  ! surface average

          E_a=evector(1,1,ij)+Er_grid(i)           ! total radial component   ! d\phi/da
          E_t=evector(2,1,ij)                      ! d\delta\phi/d\theta
          E_f=evector(3,1,ij)                      ! d\delta\phi/d\zeta
          adum_tot(i)=adum_tot(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                                              +E_t**2*gij(2,ij)  &
                                              +E_f**2*gij(3,ij)  &
                                              +2.0*E_a*E_t*gij(4,ij))

          E_a=E_a-phip00(i)-Er_grid(i) ! turbulence component only    ! d\delta\phi/da
          adum_tur(i)=adum_tur(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                                              +E_t**2*gij(2,ij)  &
                                              +E_f**2*gij(3,ij)  &
                                              +2.0_wp*E_a*E_t*gij(4,ij))
       enddo
    enddo
!$omp end parallel
   end if   !  end if(ZF_SEPARATION==0)
 
    call MPI_ALLREDUCE(adum,fluxpsi,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
    fluxpsi=fluxpsi/jaco_sum/(Tc*Tc)  ! =<(phi/Tc)^2> -- turbulence field intensity

    call MPI_ALLREDUCE(adum_00,En_00,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
    En_00=En_00/jaco_sum      ! =<E_00^2>    -- zonal flow energy

    call MPI_ALLREDUCE(adum_tur,En_tur,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
    En_tur=En_tur/jaco_sum     ! =<\deltaE^2> -- turbulence energy

    call MPI_ALLREDUCE(adum_tot,En_tot,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
    En_tot=En_tot/jaco_sum     ! =<E^2> -- total field energy

!   call get_fluxes(wpi,dgf)

!   im=mod(istep,100)
!   if(im==0) then
!      call get_3D_flow(da,dgq,dgf)
!   end if




   call get_fluxes(wpi,dgf,dgq)

!   im=mod(istep,100)
   im=mod(istep,msnap)
   if(im==0) then
      call get_3D_flow(da,dgq,dgf,wpi(2,:))
   end if


 endif


end subroutine pushi




!
!   NEW PUSHI  -- 2012 Spring
!

subroutine pushion
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface

  implicit none

!  #define SIMPLE_SPRNG
!  #define USE_MPI
!  #include "sprng_f.h"

!  integer,parameter :: nrbin=5       !wx
  integer i,im,ip,j,jt,kz,m,ii,j11,j10,j01,j00,jtheta,ierror,larmor,ij,ipjt
  real(wp) rdum,dtime,r,grada,vpara,energy,kappa,vdr,wdot,wdrift,&
       vdrenergy,delr,delt(0:mpsi),delz,pi2_inv,adum(0:mpsi),&
       xdum,ydum,xdot,ydot,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,&
       epara,psimax,psimin,pdum,wdrive,ainv,psitmp,thetatmp,&
       zetatmp,rhoi,e1,e2,e3,wpi(3,mi),rdel, &
       dgr(mi),da(mi),dgq(mi),dgf(mi),mu, &
       Er_particle(mi),dela,fip,weqflow,T_inverse,& 
       B_inverse,energy1,vdq,wKH,wgradB,omega_t,u_vpara,  &
       v2, avw,   &
! following working arrays are defined to use OMP reduction:
       adum_00(0:mpsi),adum_tur(0:mpsi),adum_tot(0:mpsi),adum_11(0:mpsi) 

!  integer ip_index(mimax)   !wx

  integer gcmotion

  real(wp) tt0,tt1,cput0,cput1,rtc

  real(wp), external :: D01AHF,dsda,SMOOTHSTEP
  integer minpts,ifail
  real(wp) acc,error_8
  real(wp) E_a,E_t,E_f,timef
  real(wp)  dminv,dmaxv
  real(wp)  wfix     ! =0 in the boundary buffer, and =1 otherwise
  real(wp)  rdm1, rdm2

  real(wp) ddT_da_particle, w_anti_relax, muB, instep,outstep,fdamp,fdamp0, rdum1

  real(wp) wfield(3,mgrid,0:1)   ! weight-field W_0+v_//*W_1+v^2*W_2

  delr=1.0_wp/deltar
  delt=1.0_wp/deltat
  delz=1.0_wp/deltaz
  pi2_inv=0.5_wp/pi
  psimax=psi_grid(mpsi)    
  psimin=psi_grid(0)       

  fdamp0=1.0_wp/tstep*5.0d-4
  if(nonlinear<0.5_wp) fdamp0=0.0_wp
  fdamp0=0.0_wp

  if(irk==1)then
! 1st step of Runge-Kutta method
     dtime=0.5_wp*tstep

! use (x,y) coordinates near magnetic axis
!$omp parallel do private(m,xdum,ydum)
     do m=1,mi
        zion0(1:5,m)=zion(1:5,m)
        zion0(7,m)=zion(7,m)
#ifdef _TREAT_MAGNETIC_AXIS
        if(sw(m).eq.1) then
           xdum=zion(1,m)*cos(zion(2,m))
           ydum=zion(1,m)*sin(zion(2,m))
           zion0(1,m)=xdum
           zion0(2,m)=ydum
        end if
#endif
     enddo

! 2nd step of Runge-Kutta method
  else
     dtime=tstep
  endif

!  mu(1:mi)=zion(6,1:mi)**2

! gather e_field using 4-point gyro-averaging, sorting in poloidal angle
! The following line is an OpenMP directive for loop-level parallelism
! on shared memory machines (see http://www.openmp.org).
!
! Ion locations in the 3D field mesh, wzion, wpion, jtion0, jtion1, wtion0, wtion1 
! are calculated at CHARGE_ION
!
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

!    prepare equilibrium Er=dPhi/da at particles 



!  ------------------------------------------------------------

!  apply equilibrium Er in the middle insteady at the beginning

!  if(istep.lt.2800) then
if(istep.lt.26500) then
   timef=0.0_wp
else
   timef=1.0_wp-exp(-(real(istep-26500)/5.0d2)**4)
end if

  timef=1.0_wp

  Er401_bp=Er401_bp0*timef
  Er_grid=Er_grid0*timef


!  if(istep.lt.2950) then
!     timef=0.0_wp
!  else
!     timef=1.0_wp-exp(-(real(istep-2950)/5.0d2)**4)
!  end if
!  Er401_bp=Er401_bp0*timef
!  Er_grid=Er_grid0*timef

!  ------------------------------------------------------------



  dela=1.0_wp/Ydr

! if(istep.ge.4000) then
!$omp parallel do private(m,r,ip,fip)
  do m=1,mi
     r=zion(1,m)-a0      
     r=r*dela
     ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
     fip=r-real(ip,wp)
     Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
  enddo
! else
!   Er_particle=0.0
! end if

!  tt0=rtc()
!  call cpu_time(cput0)

!!  m=gcmotion(dgr,da,dgq,dgf,zion(4,:),zion(1,:),zion(2,:),mu,wpi(1,:),wpi(2,:),wpi(3,:),Er_particle,nonlinear,mi)

  m=gcmotion(dgr,da,dgq,dgf,zion,wpi,Er_particle,nonlinear,mi)

!  tt1=rtc()
!  call cpu_time(cput1)

!  call n_T_u_group  ! prepare n, T, u_//, dlnn/da, dlnT/da, etc. 
                     ! moved to collision.f90 -- caution: in this case, BB must be obtained by call esigetb etc.!!!

!     OUTPUT: dgr=drho/dt;       da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt; 
!             dgf=dzeta/dt;      ghq(:)=da/dt|_EXB;     Ta(:)=dtheta/dt|_EXB;       PP(:)=drho/dt|_EXB;
!             BB(:)=B;           Ba(:)=dB/da;           Bq(:)=dB/dtheta;  
!             F(:) =F;           Fa(:)=dF/da;           gha(:)=1/|grad a|
!             no change for zion(1,2,4,:), wpi(1,2,:) and mu;
!             wpi(3,:)=dphi/dzeta (no change);          Pa(:)=R (major radius)
!             TT(:)=1/|grad{theta)|; added Aug. 2012 for poloidal momentum flux diagnosis; changed to:
!             TT(:)=poloidal velocity =vec(v).vec(B_p)/B_p; Oct. 2012

!
!        rdum1=maxval(abs(zion(5,1:mi)))
!        if(rdum1>100.0_wp) print*, 'p1111=',mype, istep, irk, rdum1, mi
!



!!!!  if(nonlinear>0.5_wp) then            ! for nonlinear run
!!!!  if(nonlinear>100.0_wp) then          ! for nonlinear run
if(wimcompressibility==1) then            ! for implementation of weight imcompressibility
   call weight_field(wfield)      ! calculate weight-field uisng zion0 for irk=1 and zion at half-step for irk=2
                                  ! used BB obtained from gcmotion, and calculte its own ion locations  
                                  ! which are used to calculate wfield in phase space grids and gether wfield at ions (next)

! 
! rdum1=maxval(abs(wfield))
! if(rdum1>100.0_wp) then
!    print*, 'p8888=', mype, istep, irk, rdum1, maxval(abs(wfield(1,1:mgrid,1))) , &
!                 maxval(abs(wfield(2,1:mgrid,1))) , maxval(abs(wfield(3,1:mgrid,1))) , maxval(abs(wfield(1,1:mgrid,0))) , & 
!                 maxval(abs(wfield(2,1:mgrid,0))) , maxval(abs(wfield(3,1:mgrid,0)))
! end if
! 

end if




dminv=0.0d0
dmaxv=0.0d0


! update GC position
! Another OpenMP parallel loop...
!$omp parallel do private(m,mu,vpara,energy,kappa,wdrift,weqflow,vdr,wdrive,wdot,xdum,ydum,T_inverse,energy1,B_inverse,u_vpara,vdq,wgradB,wKH, &
!$omp& omega_t,r,ip,fip,ddT_da_particle,w_anti_relax,muB,instep,outstep,fdamp,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,v2,avw)
  do m=1,mi
!     dbdz=0.0   for axisymmetric system   !wx

     mu=zion(6,m)*zion(6,m)    ! =v_perp^2/2B
     muB=mu*BB(m)              ! =v_perp^2/2
     T_inverse=1.0_wp/backT(m)
     B_inverse=1.0_wp/BB(m)
     vpara=zion(4,m)*BB(m)     ! BB(m)=B
     omega_t=backu(m)
     backu(m)=backu(m)*F(m)*B_inverse         ! =u_//=I*omega_t/B
     u_vpara=backu(m)/vpara                   ! =(u_//)/(v_//)
!     energy=0.5_wp*(vpara-backu(m))*(vpara-backu(m))+mu*BB(m)
     energy=0.5_wp*(vpara-backu(m))*(vpara-backu(m))+muB
     energy1=energy-backu(m)*backu(m)*0.5_wp  ! =(v_//-u)^2/2+mu*B-u^2/2 -- see note

     v2=vpara*vpara+muB*2.0_wp        ! v^2

     kappa=(min(umax*umax,energy1*T_inverse)-1.5_wp)*(-dlnTda(m))+(-dlnnda(m))+gh(m) !~-dlnf_SM/da
                                                                                  !gh=dln<A>/da with A=exp(u^2/2T)
     r=zion(7,m)-s0
     r=r*delr
     ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]; using coordinate s to find particle radial location in grids
     fip=(zion(1,m)-psi_grid(ip))/(psi_grid(ip+1)-psi_grid(ip))
     ddT_da_particle=fip*ddT_da(ip+1)+(1.0_wp-fip)*ddT_da(ip) ! d<delta_Ti>/da at particle's position

     if(ifneoclassical==1) then 
        vdr=da(m)              !  full da/at      -- including NC dynamics
        vdq=dgq(m)             !  full dtheta/dt  -- including NC dynamics
     else
        !   ExB drift in radial direction for w-dot and flux diagnostics
        vdr=ghq(m)             !  da/at|_EXB
        vdq=Ta(m)              !  dtheta/dt|_EXB
     endif

     wdrive=vdr*kappa

     weqflow=vdr*Ba(m)+vdq*Bq(m)

! new term added to maintain constant grad(Ti) drive
     w_anti_relax=(min(umax*umax,energy*T_inverse)-1.5_wp)*ddT_da_particle*T_inverse
     w_anti_relax=vdr*w_anti_relax

     if(zion(1,m).le.lmax.or.zion(1,m).ge.rmin) then
        instep =smoothstep(lmin,lmax,zion(1,m))
        outstep=smoothstep(rmin,rmax,zion(1,m))
!        fip=smoothstep(0.20_wp,0.25_wp,zion(1,m))*(1.0_wp-smoothstep(0.75_wp,0.80_wp,zion(1,m)))
        fip=instep*(1.0_wp-outstep)
        !  wdrive=wdrive*fip
        w_anti_relax=w_anti_relax*fip
        fdamp=(1.0_wp-fip)*fdamp0
     else
        fip=1.0_wp
        fdamp=0.0_wp
     end if

     if(ifneoclassical==1) then
        wdrift=(dgr(m)*BB(m)+zion(4,m)*weqflow)*(vpara-backu(m))*T_inverse   !  full drho/dt -- including NC dynamics
     else
        wdrift=(PP(m)*BB(m)+zion(4,m)*weqflow)*(vpara-backu(m))*T_inverse   ! see note
     endif

     wgradB=weqflow*T_inverse*mu       ! interesting term propotional to grad_B

     wKH=-(vdr*(duda(m)*F(m)+Fa(m)*omega_t-Ba(m)*backu(m)) &
           -vdq*Bq(m)*backu(m))*B_inverse*vpara*T_inverse       ! major Kelvin-Helmholtz drive


!
!      calculate average weight <w> at particles using wieght-field wfield
!
!!!!!  if(nonlinear>0.5_wp) then          ! for nonlinear run
!!!!!  if(nonlinear>100.0_wp) then        ! for nonlinear run
  if(wimcompressibility==1) then          ! for weight imcompressibility
     wz1=wzion(m)
     wz0=1.0_wp-wz1

     ij=jtion0(1,m)
     wp0=1.0_wp-wpion(1,m)

     wt00=1.0_wp-wtion0(1,m)
     e1=wp0*wt00*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=wp0*wt00*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=wp0*wt00*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt10=1.0_wp-wt00
     e1=e1+wp0*wt10*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp0*wt10*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp0*wt10*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=jtion1(1,m)
     wp1=1.0_wp-wp0
     wt01=1.0_wp-wtion1(1,m)
     e1=e1+wp1*wt01*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt01*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt01*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt11=1.0_wp-wt01
     e1=e1+wp1*wt11*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt11*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt11*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     avw=e1+e2*vpara+e3*v2


! 
!      if(abs(avw)>100.0_wp) then
!         print*, 'p9999=', mype, istep, irk, avw, e1, e2, e3, vpara, v2, maxval(abs(wfield(1,1:mgrid,1))) , &
!                 maxval(abs(wfield(2,1:mgrid,1))) , maxval(abs(wfield(3,1:mgrid,1))) , maxval(abs(wfield(1,1:mgrid,0))) , & 
!                 maxval(abs(wfield(2,1:mgrid,0))) , maxval(abs(wfield(3,1:mgrid,0)))
!      end if
! 


  end if


!     if(mype==100.and.m==mi/2) print*, 'AVW=', avw, zion(5,m), mype
!     if(mype==0.and.m==mi/4) print*, 'AVW=', avw, zion(5,m), mype
!     if(mype==200.and.m==mi/8) print*, 'AVW=', avw, zion(5,m), mype

!     wdot=(zion0(6,m)-nonlinear*zion(5,m))*(wdrive+wdrift+wKH+wgradB+w_anti_relax)

      if(wimcompressibility==1) then     ! for weight imcompressibility
         !  wdot=(zion0(6,m)-nonlinear*avw)*(wdrive*fip+wdrift+wKH+wgradB)   ! Correct statement for weight_damping
         wdot=(zion0(6,m)-nonlinear*avw)*(wdrive+wdrift+wKH+wgradB)*fip   ! Correct statement for weight_damping
      else
         wdot=(zion0(6,m)-nonlinear*zion(5,m))*(wdrive+wdrift+wKH+wgradB)*fip
      end if




      
!      dminv=min(dminv,(wdrive+wdrift+wKH+wgradB)*tstep)
!      dmaxv=max(dmaxv,(wdrive+wdrift+wKH+wgradB)*tstep)





! if(mype==0.and.m==mi/2) print*, 'GGGGG', wdrive, w_anti_relax

     zion(1,m) = zion0(1,m)+dtime*da(m)
     zion(2,m) = zion0(2,m)+dtime*dgq(m)
     zion(3,m) = zion0(3,m)+dtime*dgf(m)
     zion(4,m) = zion0(4,m)+dtime*dgr(m)

     !  if(zion0(1,m).gt.psi_grid(mpsi-nbound-1))  then
     !     wfix=1.0_wp
     !  else
     !     wfix=1.0_wp
     !  endif

     zion(5,m) = zion0(5,m)+dtime*(wdot-fdamp*zion(5,m))
!     zion(5,m) = zion0(5,m)+dtime*(wdot-fdamp*zion(5,m))*wfix
!     zion(7,m) = zion0(7,m)+dtime*sqrt(Tc*T_inverse)*da(m)
!     zion(7,m) = zion0(7,m)+dtime*(Tc*T_inverse)**0.4_wp*da(m)
     zion(7,m) = zion0(7,m)+dtime*(Tc*T_inverse)**dsda_power*da(m)

#ifdef _TREAT_MAGNETIC_AXIS
     if(sw(m).eq.1) then    ! convert to (a, theta) coordinates
        xdum=sqrt(zion(1,m)**2+zion(2,m)**2)
        ydum=sign(1.0_wp,zion(2,m))*acos(zion(1,m)/xdum)
        zion(1,m)=xdum
        zion(2,m)=ydum
     end if
#endif
     
! theta and zeta normalize to [0,2*pi), modulo is slower than hand coded procedure
!     zion(2,m)=zion(2,m)*pi2_inv+100.0 !period of 1
!     zion(2,m)=2.0*pi*(zion(2,m)-aint(zion(2,m))) ![0,2*pi)
!     zion(3,m)=zion(3,m)*pi2_inv+100.0
!     zion(3,m)=2.0*pi*(zion(3,m)-aint(zion(3,m)))

! use f90 modulo to normalize theta and zeta to [0,2*pi)
! If P > 0, then the function modulo(A,P) determines the positive difference between A and then
! next lowest multiple of P. If P < 0, then the function modulo(A,P) determines the negative
! difference between A and then next highest multiple of P.
     zion(2,m)=modulo(zion(2,m),twopi)
     zion(3,m)=modulo(zion(3,m),twopi)

! store GC information for flux measurements
     wpi(1,m)=ghq(m)    ! at present, wpi(1,*)=da/dt|_EXB 
     wpi(2,m)=energy    ! =(v_para-u//)^2/2+v_perp^2/2
     wpi(3,m)=muB       ! =v_perp^2/2, need for calculate distribution function in get_fluxes

  enddo
      

 if(irk==2)then
! out of boundary particle
!$omp parallel do private(m)
    do m=1,mi
       if(zion(1,m) > psimax)then
          !  rdm1=sprng()
          !  rdm1=dsign(1.0d0,(rdm1-0.5d0))
          !  rdm2=sprng()
          !  rdm2=dsign(1.0d0,(rdm2-0.5d0))
          zion(1,m)=zion0(1,m)
          zion(2,m)=2.0_wp*pi-zion0(2,m)
          zion(3,m)=zion0(3,m)
          zion(4,m)=zion0(4,m)
          !  zion(4,m)=zion0(4,m)*rdm1  !  for test 
          zion(5,m)=zion0(5,m)  ! this sounds reasonable: energy, momentum and particle conserving boundary
          !  zion(5,m)=zion0(5,m)*rdm2  !  for test 
          !  zion(5,m)=-zion0(5,m)  ! remove density (weight) accumulation at edge when ifneoclassical=1
          !  zion(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
          zion(7,m)=zion0(7,m)
          
       elseif(zion(1,m) < psimin)then
          zion(1,m)=zion0(1,m)
          zion(2,m)=2.0_wp*pi-zion0(2,m)
          zion(3,m)=zion0(3,m)
          zion(4,m)=zion0(4,m) 
          zion(5,m)=zion0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
          !  zion(5,m)=-zion0(5,m)  ! remove density (weight) accumulation at edge when ifneoclassical=1
          !  zion(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
          zion(7,m)=zion0(7,m)
       endif
    enddo

! ********************************
!if(mype==2) then
!!  error_8=1.0d-5
!!  IFAIL=1
!!  rdum=D01AHF(0.0d0,zion(1,mi/2),error_8,minpts,acc,dsda,-1,IFAIL)
!!  print*, '##################', zion(7,mi/2),'  ', rdum
!   write(stdout,*)'######## gcmotion: wall =',(tt1-tt0),'  cpu =',(cput1-cput0)
!end if
! ********************************

 endif





!  if(mype==0) print*, 'IONS:       Max Damping=', dmaxv, '   Min Damping=', dminv
!  if(mype==0) print*, 'IONS:       ', 'Damping=', dmaxv, '   Min Damping=', dminv




 if(idiag==0)then
! fluxes diagnose
!    wpi(1,1:mi)=wpi(1,1:mi)*gha(1:mi)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
!$omp parallel do private(m)
    do m=1,mi
       wpi(1,m)=wpi(1,m)*gha(m)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
    end do

    rmarker=0.0_wp
    eflux=0.0_wp
    pfluxi=0.0_wp
    dflowi=0.0_wp
    entropyi=0.0_wp
    fluxpsi=0.0_wp
    heatflux=0.0_wp
    En_00=0.0_wp
    En_tur=0.0_wp
    En_tot=0.0_wp
    do m=1,mi
       r=zion0(7,m)         ! old radial location
       ip=max(0,min(mpsi,int((r-s0)*delr+0.5_wp)))    ! ip [0,mpsi]
       vdrenergy=wpi(1,m)*(wpi(2,m)-2.5*vth_grid(ip)**2)*zion0(5,m) 
                ! note: wpi(1,m)=dr/dt (EXB part) 

! radial location in diagnostic bin
       i=max(1,min(mflux,1+int((r-s0)*real(mflux)/(s1-s0))))

! marker,energy,particle,momentum fluxes,parallel flows,entropy and kinetic energy
       rmarker(i)=rmarker(i)+zion0(6,m)
       eflux(i)=eflux(i)+vdrenergy
       pfluxi=pfluxi+wpi(1,m)*zion0(5,m)
       dflowi=dflowi+BB(m)*zion0(4,m)*zion0(5,m)
       entropyi=entropyi+zion0(5,m)*zion0(5,m)
    enddo

   if(ZF_SEPARATION==0) then   

!     note: zonal flow component <phi> and phip00=d<phi>/da are calculated in field.f90

    adum=0.0_wp
    adum_00=0.0_wp
    adum_tur=0.0_wp
    adum_tot=0.0_wp

!$omp parallel private(i,j,ij,E_a,E_t,E_f) reduction(+:adum,adum_00,adum_tur,adum_tot)
! Note from OMP compiler: A variable in the REDUCTION clause or any component or element thereof must not have the POINTER 
! attribute; nor be a pointee, an allocatable object, an assumed-size array, or a THREADLOCAL common block variable.
!$omp do
    do i=0,mpsi
       do j=1,mtheta(i)
          ij=igrid(i)+j
          adum(i)=adum(i)+dtemper(ij)*(phi(1,ij)-phi00(i))**2  ! surface averaging for delta-phi

          adum_00(i)=adum_00(i)+dtemper(ij)*phip00(i)**2*gij(1,ij)  ! surface average

          E_a=evector(1,1,ij)+Er_grid(i)           ! total radial component   ! d\phi/da
          E_t=evector(2,1,ij)                      ! d\delta\phi/d\theta
          E_f=evector(3,1,ij)                      ! d\delta\phi/d\zeta
          adum_tot(i)=adum_tot(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                                              +E_t**2*gij(2,ij)  &
                                              +E_f**2*gij(3,ij)  &
                                              +2.0_wp*E_a*E_t*gij(4,ij))

          E_a=E_a-phip00(i)-Er_grid(i) ! turbulence component only    ! d\delta\phi/da
          adum_tur(i)=adum_tur(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                                              +E_t**2*gij(2,ij)  &
                                              +E_f**2*gij(3,ij)  &
                                              +2.0_wp*E_a*E_t*gij(4,ij))
       enddo
    enddo
!$omp end parallel

   else  ! ZF_SEPARATION=1

    adum=0.0_wp
    adum_00=0.0_wp
    adum_tur=0.0_wp
    adum_tot=0.0_wp
    adum_11=0.0_wp

! DEBUG 
!        if(mype==100) print*, 'Dirichlet2', phi(1,igrid(0)), phi(1,igrid(0)+mtheta(0)/2),phi(1,igrid(0)+mtheta(0)-1), &
!                               phi(1,igrid(mpsi)), phi(1,igrid(mpsi)+mtheta(mpsi)/2),phi(1,igrid(mpsi)+mtheta(mpsi)-1)

!$omp parallel private(i,j,ij,E_a,E_t,E_f) reduction(+:adum,adum_00,adum_tur,adum_tot,adum_11)
!$omp do
    do i=0,mpsi
       do j=1,mtheta(i)
          ij=igrid(i)+j
!          fluxpsi(i)=fluxpsi(i)+phi(1,ij)*phi(1,ij)               ! what average?
          adum(i)=adum(i)+dtemper(ij)*phi(1,ij)*phi(1,ij)  ! surface averaging for <delta-phi^2>

          adum_11(i)=adum_11(i)+dtemper(ij)*phi(1,ij)      ! to check surface averaged <delta-phi> which should be zero by defination

          adum_00(i)=adum_00(i)+dtemper(ij)*phip00(i)**2*gij(1,ij)  ! surface average

          E_a=evector(1,1,ij)+Er_grid(i)           ! total radial component   ! d\phi/da
          E_t=evector(2,1,ij)                      ! d\delta\phi/d\theta
          E_f=evector(3,1,ij)                      ! d\delta\phi/d\zeta
          adum_tot(i)=adum_tot(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                                              +E_t**2*gij(2,ij)  &
                                              +E_f**2*gij(3,ij)  &
                                              +2.0*E_a*E_t*gij(4,ij))

          E_a=E_a-phip00(i)-Er_grid(i) ! turbulence component only    ! d\delta\phi/da
          adum_tur(i)=adum_tur(i)+dtemper(ij)*(E_a**2*gij(1,ij)  &  ! surface average
                                              +E_t**2*gij(2,ij)  &
                                              +E_f**2*gij(3,ij)  &
                                              +2.0_wp*E_a*E_t*gij(4,ij))
       enddo
    enddo
!$omp end parallel
   end if   !  end if(ZF_SEPARATION==0)

    call MPI_ALLREDUCE(adum_11,flux_11,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
    flux_11=flux_11/jaco_sum/Tc   ! =<phi>/Tc -- added in March 2013 to check whether <phi>=0 is obtained
 
    call MPI_ALLREDUCE(adum,fluxpsi,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
    fluxpsi=fluxpsi/jaco_sum/(Tc*Tc)  ! =<(phi/Tc)^2> -- turbulence field intensity
    sumfluxpsi=sumfluxpsi+fluxpsi     ! =int[I(a)dt]

    call MPI_ALLREDUCE(adum_00,En_00,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
    En_00=En_00/jaco_sum      ! =<E_00^2>    -- zonal flow energy

    call MPI_ALLREDUCE(adum_tur,En_tur,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
    En_tur=En_tur/jaco_sum     ! =<\deltaE^2> -- turbulence energy

    call MPI_ALLREDUCE(adum_tot,En_tot,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
    En_tot=En_tot/jaco_sum     ! =<E^2> -- total field energy

   call get_fluxes(wpi,dgf,dgq)

!   im=mod(istep,100)
   im=mod(istep,msnap)
   if(im==0) then
      call get_3D_flow(da,dgq,dgf,wpi(2,:))
   end if

 endif


end subroutine pushion
