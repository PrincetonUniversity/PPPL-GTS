subroutine pushe(irk_e)
  use precision
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface

  implicit none

!  integer,parameter :: nrbin=5       !wx
  integer i,im,ip,j,jt,k,kz,kk,m,ii,j11,j10,j01,j00,jtheta,ierror,ij,ipjt,irk_e
  real(wp) rdum,dtime,r,grada,energy,kappa,vdr,wdot,wdrift,&
       vdrenergy,delr,delt(0:mpsi),delz,pi2_inv,adum(0:mpsi),&
       xdum,ydum,xdot,ydot,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,&
       epara,psimax,psimin,pdum,wdrive,ainv,psitmp,&
       zetatmp,e1,e2,e3,wpi(3,me),rdel, &
       dgr(me),da(me),dgq(me),dgf(me),mu, &
       Er_particle(me),dela,fip,weqflow,T_inverse,wdphidt, & 
       vpara, v2, avw, f_rho,  &
! following working arrays are defined to use OMP reduction:
       adum_00(0:mpsi),adum_tur(0:mpsi),adum_tot(0:mpsi)

  integer gcmotion_e

  real(wp) tt0,tt1,cput0,cput1,rtc

  real(wp), external :: D01AHF,dsda,SMOOTHSTEP
  integer minpts,ifail
  real(wp)  acc,error_8
  real(wp)  E_a,E_t,E_f

  real(wp) muB, instep,outstep,fdamp,fdamp0

  real(wp) wfield(3,mgrid,0:1)   ! weight-field W_0+v_//*W_1+v^2*W_2

  delr=1.0_wp/deltar
  delt=1.0_wp/deltat
  delz=1.0_wp/deltaz
  pi2_inv=0.5_wp/pi
  psimax=psi_grid(mpsi)    
  psimin=psi_grid(0)       

  f_rho=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)

  fdamp0=1.0_wp/tstep*5.0d-4
  if(nonlinear<0.5_wp) fdamp0=0.0_wp
!  fdamp0=0.0_wp

  if(irk_e==1)then
! 1st step of Runge-Kutta method
     dtime=0.5_wp*tstep/real(nelectron,wp)

! use (x,y) coordinates near magnetic axis
!$omp parallel do private(m,xdum,ydum)
     do m=1,me
        zeon0(1:5,m)=zeon(1:5,m)
        zeon0(7,m)=zeon(7,m)
#ifdef _TREAT_MAGNETIC_AXIS
        if(sw(m).eq.1) then
           xdum=zeon(1,m)*cos(zeon(2,m))
           ydum=zeon(1,m)*sin(zeon(2,m))
           zeon0(1,m)=xdum
           zeon0(2,m)=ydum
        end if
#endif
     enddo

! 2nd step of Runge-Kutta method
  else
     dtime=tstep/real(nelectron,wp)
  endif

  call electron_locator

!  mu(1:me)=zeon(6,1:me)**2

! gather e_field at electron GC positions using linear interpolation.
! No FLR effect for electrons.
! The following line is an OpenMP directive for loop-level parallelism
! on shared memory machines (see http://www.openmp.org).
!$omp parallel do private(m,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11)
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

!    prepare equilibrium Er=dPhi/da at particles 
  dela=1.0_wp/Ydr
! if(istep.ge.4000) then
!$omp parallel do private(m,r,ip,fip)
  do m=1,me
     r=zeon(1,m)-a0      
     r=r*dela
     ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
     fip=r-real(ip,wp)
     Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
  enddo
! else
! Er_particle=0.0
! end if

!  tt0=rtc()
! call cpu_time(cput0)

  m=gcmotion_e(dgr,da,dgq,dgf,zeon,wpi,Er_particle,nonlinear,me,zeondim)

!  tt1=rtc()
!  call cpu_time(cput1)

!  call n_T_u_group_e  ! prepare n, T, u_//, dlnn/da, dlnT/da, etc at electrons.
                       ! moved to collision_eon.f90 -- caution: in this case, BB must be obtained by call esigetb etc.!!!

  call Ti_at_electrons ! NOTE: array backT is borrowed to store Ti_at_electrons;
                       ! WARNING: must make sure there is no memory conflict!!!
                       ! backT is used in PUSHI and released afterward

!     OUTPUT: dgr=drho/dt;       da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt; 
!             dgf=dzeta/dt;      ghq(:)=da/dt|_EXB;     wpi(3,:)=dphi/dzeta (NO change); 
!             BB(:)=B;           no change for zeon(1,2,4,:), wpi(1,2,:) and mu;
!             gha(:)=1/|grad a|  F(:) =F

if(nonlinear>0.5_wp) then            ! for nonlinear run
  call weight_field_eon(wfield)      ! calculate weight-field uisng zeon0 for irk=1 and zeon at half-step for irk=2
                                     ! used BB obtained from gcmotion, and electron locations is obtained here in pushe  
                                     ! which are used to calculate wfield in phase space grids and gether wfield at electrons (next)
end if

! update GC position
! Another OpenMP parallel loop...
!$omp parallel do private(m,mu,energy,kappa,wdrift,weqflow,vdr,wdrive,wdot,xdum,ydum,T_inverse,r, &
!$omp& e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,wdphidt,muB,fip,instep,outstep,fdamp,&
!$omp& vpara,v2,avw)
  do m=1,me
!     dbdz=0.0   for axisymmetric system   !wx

     mu=zeon(6,m)*zeon(6,m)        ! =normalized mu of electrons
     muB=mu*BB(m)                  ! =normalized perpendicular energy of electrons
     T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti
!     energy=0.5_wp*for_E*(zeon(4,m)*BB(m))**2+mu*BB(m)

     vpara=zeon(4,m)*BB(m)
     v2=vpara**2

!     energy=0.5_wp*for_E*(zeon(4,m)*BB(m))**2+muB
     energy=0.5_wp*for_E*v2+muB

     vpara=vpara*f_rho                ! =normalized v_//
     v2=vpara*vpara+muB*2.0_wp*m0ms   ! =v^2

     kappa=(min(umax*umax,energy*T_inverse)-1.5_wp)*(-dlnTdae(m))+(-dlnndae(m)) !=-dlnf_M/da

! ExB drift in radial direction for w-dot and flux diagnostics
     vdr=ghq(m)             !  da/dt|_EXB
     wdrive=vdr*kappa
     weqflow=-vdr*Er_particle(m)*csc0*T_inverse   ! additional term due to E_0XB drift 

     if(zeon(1,m).le.lmax.or.zeon(1,m).ge.rmin) then
        instep =smoothstep(lmin,lmax,zeon(1,m))
        outstep=smoothstep(rmin,rmax,zeon(1,m))
!        fip=smoothstep(0.20_wp,0.25_wp,zeon(1,m))*(1.0_wp-smoothstep(0.75_wp,0.80_wp,zeon(1,m)))
        fip=instep*(1.0_wp-outstep)
        wdrive=wdrive*fip
        fdamp=(1.0_wp-fip)*fdamp0
     else
        fdamp=0.0_wp
     end if

     wdrift=-(da(m)*wpi(1,m)+dgq(m)*wpi(2,m)+dgf(m)*wpi(3,m))*csc0*T_inverse

     wdot=(wdrive+wdrift+weqflow)


if(nonlinear>0.5_wp) then     ! for nonlinear run
! calculate average weight <w> at particles using wieght-field wfield
     wz1=wzeon(m)
     wz0=1.0_wp-wz1

     ij=jteon0(m)
     wp0=1.0_wp-wpeon(m)

     wt00=1.0_wp-wteon0(m)
     e1=wp0*wt00*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=wp0*wt00*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=wp0*wt00*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt10=1.0_wp-wt00
     e1=e1+wp0*wt10*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp0*wt10*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp0*wt10*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=jteon1(m)
     wp1=1.0_wp-wp0
     wt01=1.0_wp-wteon1(m)
     e1=e1+wp1*wt01*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt01*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt01*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt11=1.0_wp-wt01
     e1=e1+wp1*wt11*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt11*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt11*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     avw=e1+e2*vpara+e3*v2

!     if(mype==100.and.m==me/2) print*, 'AVW=', avw, zeon(5,m), mype
!     if(mype==0.and.m==me/4) print*, 'AVW=', avw, zeon(5,m), mype
!     if(mype==200.and.m==me/8) print*, 'AVW=', avw, zeon(5,m), mype
end if


     zeon(1,m) = zeon0(1,m)+dtime*da(m)
     zeon(2,m) = zeon0(2,m)+dtime*dgq(m)
     zeon(3,m) = zeon0(3,m)+dtime*dgf(m)
     zeon(4,m) = zeon0(4,m)+dtime*dgr(m)

!     zeon(5,m) = zeon0(5,m)+dtime*wdot

!     zeon(5,m) = 1.0_wp-(1.0_wp-zeon0(5,m))*exp(-dtime*wdot)

!     zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*zeon(5,m))*wdot-fdamp*zeon(5,m))*dtime

     zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*avw)*wdot-fdamp*zeon(5,m))*dtime

!     zeon(7,m) = zeon0(7,m)+dtime*sqrt(Tc/backT(m))*da(m)
!     zeon(7,m) = zeon0(7,m)+dtime*(Tc/backT(m))**0.4_wp*da(m)
     zeon(7,m) = zeon0(7,m)+dtime*(Tc/backT(m))**dsda_power*da(m)

#ifdef _TREAT_MAGNETIC_AXIS
     if(sw(m).eq.1) then    ! convert to (a, theta) coordinates
        xdum=sqrt(zeon(1,m)**2+zeon(2,m)**2)
        ydum=sign(1.0_wp,zeon(2,m))*acos(zeon(1,m)/xdum)
        zeon(1,m)=xdum
        zeon(2,m)=ydum
     end if
#endif
     
! theta and zeta normalize to [0,2*pi), modulo is slower than hand coded procedure
!     zeon(2,m)=zeon(2,m)*pi2_inv+100.0 !period of 1
!     zeon(2,m)=2.0*pi*(zeon(2,m)-aint(zeon(2,m))) ![0,2*pi)
!     zeon(3,m)=zeon(3,m)*pi2_inv+100.0
!     zeon(3,m)=2.0*pi*(zeon(3,m)-aint(zeon(3,m)))

!     zeon(2,m)=zeon(2,m)*pi2_inv !period of 1
!     zeon(2,m)=10.0+zeon(2,m)-dint(zeon(2,m)) ![0,2*pi)
!     zeon(2,m)=twopi*(zeon(2,m)-dint(zeon(2,m)))
!     zeon(3,m)=zeon(3,m)*pi2_inv
!     zeon(3,m)=10.0+zeon(3,m)-dint(zeon(3,m))
!     zeon(3,m)=twopi*(zeon(3,m)-dint(zeon(3,m)))

! use f90 modulo to normalize theta and zeta to [0,2*pi)
! If P > 0, then the function modulo(A,P) determines the positive difference between A and then
! next lowest multiple of P. If P < 0, then the function modulo(A,P) determines the negative
! difference between A and then next highest multiple of P.
     zeon(2,m)=modulo(zeon(2,m),twopi)
     zeon(3,m)=modulo(zeon(3,m),twopi)


! store GC information for flux measurements
     wpi(1,m)=vdr       ! at present, wpi(1,*)=da/dt|_EXB 
     wpi(2,m)=energy    ! normalized electron energy with u//=0 in Jan.09
     wpi(3,m)=muB       ! normalizd electron perpendicular energy

  enddo

 if(irk_e==2)then
! out of boundary particle
!$omp parallel do private(m)
    do m=1,me
       if(zeon(1,m) > psimax)then
          zeon(1,m)=zeon0(1,m)
          zeon(2,m)=2.0_wp*pi-zeon0(2,m)
          zeon(3,m)=zeon0(3,m)
          zeon(4,m)=zeon0(4,m)
          zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
!          zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
          zeon(7,m)=zeon0(7,m)
          
       elseif(zeon(1,m) < psimin)then
          zeon(1,m)=zeon0(1,m)
          zeon(2,m)=2.0_wp*pi-zeon0(2,m)
          zeon(3,m)=zeon0(3,m)
          zeon(4,m)=zeon0(4,m) 
          zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
!          zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
          zeon(7,m)=zeon0(7,m)
       endif
    enddo

! ********************************
if(mype==2) then
  error_8=1.0e-5_wp
  IFAIL=1
  rdum=D01AHF(0.0d0,zeon(1,me/2),error_8,minpts,acc,dsda,-1,IFAIL)
  print*, 'EEEEEEEEEEEE', zeon(7,me/2),'  ', rdum
!   write(stdout,*)'######## electron gcmotion: wall =',(tt1-tt0),'  cpu =',(cput1-cput0)
end if
! ********************************

 endif

 if(idiag_e==0.and.irk_e==1) then  ! istep=n*ndiag, irk=1, i=1, and irk_e=1
! fluxes diagnose
!$omp parallel do private(m)
    do m=1,me
       wpi(1,m)=wpi(1,m)*gha(m)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
    end do

    call get_fluxes_e(wpi)
 endif

end subroutine pushe







!   ********************************************************************************************************************
!   ********************************************************************************************************************
!   ********************************************************************************************************************







subroutine pusheon(solver1, masssolver, irk_e)
  use precision
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface
  use solver_data_class

  implicit none

  type(solver_data) :: solver1, masssolver

!  integer,parameter :: nrbin=5       !wx
  integer i,im,ip,j,jt,k,kz,kk,m,ii,j11,j10,j01,j00,jtheta,ierror,ij,ipjt,irk_e
  real(wp) rdum,dtime,r,grada,energy,kappa,vdr,wdot,wdrift,&
       vdrenergy,delr,delt(0:mpsi),delz,pi2_inv,adum(0:mpsi),&
       xdum,ydum,xdot,ydot,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,&
       epara,psimax,psimin,pdum,wdrive,ainv,psitmp,&
       zetatmp,e1,e2,e3,wpi(3,me),rdel, &
       dgr(me),da(me),dgq(me),dgf(me),mu, &
       Er_particle(me),dela,fip,weqflow,T_inverse,Ea_00,wdphidt,wdphidt2, & 
       vpara, v2, avw, f_rho,  &
! following working arrays are defined to use OMP reduction:
       adum_00(0:mpsi),adum_tur(0:mpsi),adum_tot(0:mpsi)

  integer gcmotion_eon

  real(wp) tt0,tt1,cput0,cput1,rtc

  real(wp), external :: D01AHF,dsda,SMOOTHSTEP
  integer minpts,ifail
  real(wp)  acc,error_8
  real(wp)  E_a,E_t,E_f
  real(wp)  dminv,dmaxv
  real(wp)  wfix     ! =0 in the boundary buffer, and =1 otherwise

  real(wp) muB, instep,outstep,fdamp,fdamp0

  real(wp) wfield(3,mgrid,0:1)   ! weight-field W_0+v_//*W_1+v^2*W_2

  delr=1.0_wp/deltar
  delt=1.0_wp/deltat
  delz=1.0_wp/deltaz
  pi2_inv=0.5_wp/pi
  psimax=psi_grid(mpsi)    
  psimin=psi_grid(0)       

  f_rho=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)

  fdamp0=1.0_wp/tstep*5.0d-4
  if(nonlinear<0.5_wp) fdamp0=0.0_wp
  fdamp0=0.0_wp

  if(irk_e==1)then
! 1st step of Runge-Kutta method
     dtime=0.5_wp*tstep/real(nelectron,wp)

! use (x,y) coordinates near magnetic axis
!$omp parallel do private(m,xdum,ydum)
     do m=1,me
        zeon0(1:5,m)=zeon(1:5,m)
        zeon0(7,m)=zeon(7,m)
#ifdef _TREAT_MAGNETIC_AXIS
        if(sw(m).eq.1) then
           xdum=zeon(1,m)*cos(zeon(2,m))
           ydum=zeon(1,m)*sin(zeon(2,m))
           zeon0(1,m)=xdum
           zeon0(2,m)=ydum
        end if
#endif
     enddo

! 2nd step of Runge-Kutta method
  else
     dtime=tstep/real(nelectron,wp)
  endif

  call electron_locator

!  mu(1:me)=zeon(6,1:me)**2

! gather e_field at electron GC positions using linear interpolation.
! No FLR effect for electrons.
! The following line is an OpenMP directive for loop-level parallelism
! on shared memory machines (see http://www.openmp.org).
!$omp parallel do private(m,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11)
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

!    prepare equilibrium Er=dPhi/da at particles 
  dela=1.0_wp/Ydr
! if(istep.ge.4000) then
!$omp parallel do private(m,r,ip,fip)
  do m=1,me
     r=zeon(1,m)-a0      
     r=r*dela
     ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
     fip=r-real(ip,wp)
     Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
  enddo
! else
! Er_particle=0.0
! end if

!  tt0=rtc()
! call cpu_time(cput0)

  m=gcmotion_eon(dgr,da,dgq,dgf,zeon,wpi,Er_particle,nonlinear,me,zeondim)

! tt1=rtc()
! call cpu_time(cput1)

!  call n_T_u_group_e  ! prepare n, T, u_//, dlnn/da, dlnT/da, etc at electrons.
                       ! moved to collision_eon.f90 -- caution: in this case, BB must be obtained by call esigetb etc.!!!

  call Ti_at_electrons ! NOTE: array backT is borrowed to store Ti_at_electrons;
                       ! WARNING: must make sure there is no memory conflict!!!
                       ! backT is used in PUSHI and released afterward

!     OUTPUT: dgr=drho/dt;        da=da/dt or dx/dt;     dgq=dtheta/dt or dy/dt;    dgf=dzeta/dt;      
!             ghq(:)=da/dt|_EXB;  Ta(:)=dtheta/dt|_EXB;  Pa(:)=drho/dt|_EXB;
!             BB(:)=B;            Ba(:)=dB/da;           Bq(:)=dB/dtheta;       
!             gha(:)=1/|grad a|   F(:) =F
!             wpi(3,:)=dphi/dzeta (NO change); no change for: zeon(1,2,4,:), wpi(1,2,:) and mu;

!  if(nonlinear>100.0_wp) then            ! for nonlinear run
if(wimcompressibility==1) then            ! for implementation of weight imcompressibility 
  call weight_field_eon(wfield)      ! calculate weight-field uisng zeon0 for irk=1 and zeon at half-step for irk=2
                                     ! used BB obtained from gcmotion, and electron locations are obtain here in pushe 
                                     ! which are used to calculate wfield in phase space grids and gether wfield at electron (next)
end if




dminv=0.0d0
dmaxv=0.0d0



IF(ifelectron==2.and.ifdphi==2) then  ! use direct derivative along trajectory to calculate dphi/dt in delta-h electron case

! update GC position
! Another OpenMP parallel loop...
!$omp parallel do private (m,mu,energy,kappa,wdrift,weqflow,vdr,wdrive,wdot,xdum,ydum,T_inverse,&
!$omp& r,ip,Ea_00,e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,wdphidt,wdphidt2,muB,fip,&
!$omp& instep,outstep,fdamp,vpara,v2,avw,E_f)
  do m=1,me
!     dbdz=0.0   for axisymmetric system   !wx

     mu=zeon(6,m)*zeon(6,m)        ! =normalized mu of electrons
     muB=mu*BB(m)                  ! =normalized perpendicular energy of electrons
     T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti
!     energy=0.5_wp*for_E*(zeon(4,m)*BB(m))**2+mu*BB(m)

     vpara=zeon(4,m)*BB(m)
     v2=vpara**2

!     energy=0.5_wp*for_E*(zeon(4,m)*BB(m))**2+muB
     energy=0.5_wp*for_E*v2+muB

     vpara=vpara*f_rho                ! =normalized v_//
     v2=vpara*vpara+muB*2.0_wp*m0ms   ! =v^2

!!!! TEST
     !   dlnTdae(m)=0.0_wp
     !   dlnndae(m)=0.0_wp
!!!! TEST

     kappa=(min(umax*umax,energy*T_inverse)-1.5_wp)*(-dlnTdae(m))+(-dlnndae(m)) !=-dlnf_M/da

     if(ifneoclassical==1) then
        vdr=da(m)                 !  full da/dt including NC dynamics
        weqflow=vdr*Ba(m)+dgq(m)*Bq(m)
     else
        ! ExB drift in radial direction for w-dot and flux diagnostics
        vdr=ghq(m)             !  da/dt|_EXB
        weqflow=vdr*Ba(m)+Ta(m)*Bq(m)
     end if

     wdrive=vdr*kappa


! calculate d(delta_phi)/dt at particles along trajactory
        wz1=wzeon(m)
        wz0=1.0_wp-wz1

        ij=jteon0(m)
        wp0=1.0_wp-wpeon(m)
        wt00=1.0_wp-wteon0(m)
        e1=   wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))

        ij=ij+1
        wt10=1.0_wp-wt00
        e1=e1+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))

        ij=jteon1(m)
        wp1=1.0_wp-wp0
        wt01=1.0_wp-wteon1(m)
        e1=e1+wp1*wt01*(wz0*phi(0,ij)+wz1*phi(1,ij))

        ij=ij+1
        wt11=1.0_wp-wt01
        e1=e1+wp1*wt11*(wz0*phi(0,ij)+wz1*phi(1,ij)) ! --> phi at electrons

        if(ZF_SEPARATION==0) then
!    calculate <phi> at electrons (no FLR effect -- x=R)
           r=zeon(7,m)-s0       ! in radial coordinator s
           r=r*delr
           ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
           fip=r-real(ip,wp)
           Ea_00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons
        end if

        E_f=(e1-real(1-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0
        if(irk_e==1) zeon0(8,m)=E_f    ! e*delta_phi/Te at electron beginning position



     if(zeon0(1,m).le.lmax.or.zeon0(1,m).ge.rmin) then
        instep =smoothstep(lmin,lmax,zeon(1,m))
        outstep=smoothstep(rmin,rmax,zeon(1,m))
        !  fip=smoothstep(0.20_wp,0.25_wp,zeon(1,m))*(1.0_wp-smoothstep(0.75_wp,0.80_wp,zeon(1,m)))
        fip=instep*(1.0_wp-outstep)
        !  wdrive=wdrive*fip

        fdamp=(1.0_wp-fip)*fdamp0
     else
        fip=1.0_wp
        fdamp=0.0_wp
     end if

     if(ifneoclassical==1) then   
        wdrift=(dgr(m)*BB(m)+zeon(4,m)*weqflow)*zeon(4,m)*BB(m)*for_E*T_inverse   ! see note
     else
        wdrift=(Pa(m)*BB(m)+zeon(4,m)*weqflow)*zeon(4,m)*BB(m)*for_E*T_inverse   ! see note
     end if

     weqflow=weqflow*mu*T_inverse   ! see note: additional term due to v_ExB*grad(f_M)

     wdot=(wdrive+wdrift+weqflow)

     wdot=(1.0_wp-E_f)*wdot*fip   ! (1-e*delta_phi/Te)*(...) term -- see note


!!!!!  if(nonlinear>0.5_wp) then     ! for nonlinear run
!!!!!  if(nonlinear>100.0_wp) then          ! for nonlinear run
if(wimcompressibility==1) then     ! for weight imcompressibility
! calculate average weight <w> at particles using wieght-field wfield
     wz1=wzeon(m)
     wz0=1.0_wp-wz1

     ij=jteon0(m)
     wp0=1.0_wp-wpeon(m)

     wt00=1.0_wp-wteon0(m)
     e1=wp0*wt00*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=wp0*wt00*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=wp0*wt00*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt10=1.0_wp-wt00
     e1=e1+wp0*wt10*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp0*wt10*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp0*wt10*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=jteon1(m)
     wp1=1.0_wp-wp0
     wt01=1.0_wp-wteon1(m)
     e1=e1+wp1*wt01*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt01*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt01*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt11=1.0_wp-wt01
     e1=e1+wp1*wt11*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt11*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt11*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     avw=e1+e2*vpara+e3*v2

!     if(mype==100.and.m==me/2) print*, 'AVW=', avw, zeon(5,m), mype
!     if(mype==0.and.m==me/4) print*, 'AVW=', avw, zeon(5,m), mype
!     if(mype==200.and.m==me/8) print*, 'AVW=', avw, zeon(5,m), mype

end if


     zeon(1,m) = zeon0(1,m)+dtime*da(m)
     zeon(2,m) = zeon0(2,m)+dtime*dgq(m)
     zeon(3,m) = zeon0(3,m)+dtime*dgf(m)
     zeon(4,m) = zeon0(4,m)+dtime*dgr(m)

     if(wimcompressibility==1) then
        zeon(8,m) = avw           ! save average weight avw at electrons
     else
        zeon(8,m) = zeon(5,m)     ! save w 
     end if

!     zeon(5,m) = zeon0(5,m)+dtime*wdot

!     zeon(5,m) = 1.0_wp-(1.0_wp-zeon0(5,m))*exp(-dtime*wdot)

!    two cases below are unified by the replacement using zeon(8,m) in the formula
       !!!!zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*zeon(5,m))*wdot-fdamp*zeon(5,m))*dtime  ! weight advanced w/o including dphi/dt term yet
       !!!!zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*avw)*wdot-fdamp*zeon(5,m))*dtime

     !  if(zeon0(1,m).gt.psi_grid(mpsi-nbound-1))  then 
     !     wfix=1.0_wp
     !  else
     !     wfix=1.0_wp
     !  endif

     !  zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*zeon(8,m))*wdot-fdamp*zeon(5,m))*dtime*wfix  ! weight advanced w/o including dphi/dt term yet
     zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*zeon(8,m))*wdot-fdamp*zeon(5,m))*dtime       ! weight advanced w/o including dphi/dt term yet
                                                                                                ! unified expression April, 2013




!     dminv=min(dminv,wdot*tstep)
!     dmaxv=max(dmaxv,wdot*tstep)



!     zeon(7,m) = zeon0(7,m)+dtime*sqrt(Tc/backT(m))*da(m)
!     zeon(7,m) = zeon0(7,m)+dtime*(Tc/backT(m))**0.4_wp*da(m)
     zeon(7,m) = zeon0(7,m)+dtime*(Tc/backT(m))**dsda_power*da(m)

#ifdef _TREAT_MAGNETIC_AXIS
     if(sw(m).eq.1) then    ! convert to (a, theta) coordinates
        xdum=sqrt(zeon(1,m)**2+zeon(2,m)**2)
        ydum=sign(1.0_wp,zeon(2,m))*acos(zeon(1,m)/xdum)
        zeon(1,m)=xdum
        zeon(2,m)=ydum
     end if
#endif

! theta and zeta normalize to [0,2*pi), modulo is slower than hand coded procedure
!     zeon(2,m)=zeon(2,m)*pi2_inv+100.0 !period of 1
!     zeon(2,m)=2.0*pi*(zeon(2,m)-aint(zeon(2,m))) ![0,2*pi)
!     zeon(3,m)=zeon(3,m)*pi2_inv+100.0
!     zeon(3,m)=2.0*pi*(zeon(3,m)-aint(zeon(3,m)))

!     zeon(2,m)=zeon(2,m)*pi2_inv !period of 1
!     zeon(2,m)=10.0+zeon(2,m)-dint(zeon(2,m)) ![0,2*pi)
!     zeon(2,m)=twopi*(zeon(2,m)-dint(zeon(2,m)))
!     zeon(3,m)=zeon(3,m)*pi2_inv
!     zeon(3,m)=10.0+zeon(3,m)-dint(zeon(3,m))
!     zeon(3,m)=twopi*(zeon(3,m)-dint(zeon(3,m)))

! use f90 modulo to normalize theta and zeta to [0,2*pi)
! If P > 0, then the function modulo(A,P) determines the positive difference between A and then
! next lowest multiple of P. If P < 0, then the function modulo(A,P) determines the negative
! difference between A and then next highest multiple of P.
     zeon(2,m)=modulo(zeon(2,m),twopi)
     zeon(3,m)=modulo(zeon(3,m),twopi)


! store GC information for flux measurements
     wpi(1,m)=ghq(m)    ! at present, wpi(1,*)=da/dt|_EXB 
     wpi(2,m)=energy    ! normalized electron energy with u//=0 in Jan.09
     wpi(3,m)=muB       ! normalizd electron perpendicular energy

  enddo




!  if(mype==0) print*, 'ELECTRONS:  Max Damping=', dmaxv, '   Min Damping=', dminv
!  if(mype==0) print*, 'ELECTRONS:  ', 'Max Damping=', dmaxv, '   Min Damping=', dminv




  !   fluxes diagnose
  if(idiag_e==0.and.irk_e==1) then  ! istep=n*ndiag, irk=1, i=1, and irk_e=1
!$omp parallel do private(m)
     do m=1,me
        wpi(1,m)=wpi(1,m)*gha(m)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
     end do

     call get_fluxes_e(wpi)

     !   im=mod(istep,100)
     im=mod(istep,msnap)
     if(im==0) then
        call get_Te_fluctuation(wpi(2,:))
     end if
  endif

   CALL SHIFTE

 call electron_locator
 call n_T_u_group_e             ! get Te at electron position
 zeon00(8,1:me)=zeon(5,1:me)    ! use for temporary storage -- check again April, 2013

!  E_a=0.0_wp
!  E_f=0.0_wp
!  E_t=0.0_wp

!$omp parallel do private(m,T_inverse,r,ip,fip,Ea_00, &
!$omp& e1,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,wdphidt)
  do m=1,me
     T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti

! calculate d(delta_phi)/dt at particles
      wz1=wzeon(m)
      wz0=1.0_wp-wz1

      ij=jteon0(m)
      wp0=1.0_wp-wpeon(m)
      wt00=1.0_wp-wteon0(m)
      e1=   wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=ij+1
      wt10=1.0_wp-wt00
      e1=e1+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=jteon1(m)
      wp1=1.0_wp-wp0
      wt01=1.0_wp-wteon1(m)
      e1=e1+wp1*wt01*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=ij+1
      wt11=1.0_wp-wt01
      e1=e1+wp1*wt11*(wz0*phi(0,ij)+wz1*phi(1,ij)) ! --> phi at electron new positions

      if(ZF_SEPARATION==0) then
!    calculate <phi> at electrons (no FLR effect -- x=R)
         r=zeon(7,m)-s0       ! in radial coordinator s
         r=r*delr
         ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
         fip=r-real(ip,wp)
         Ea_00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons
      end if

      wdphidt=(e1-real(1-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0   ! e*delta_phi/Te at electron new positions, but with unchanged delta_phi
      wdphidt=(wdphidt-zeon0(8,m))/dtime  ! d(e*delta_phi/Te)/dt purely due to change of electron position

      !  if(zeon0(1,m).gt.psi_grid(mpsi-nbound-1))  then
      !     wfix=1.0_wp
      !  else
      !     wfix=1.0_wp
      !  endif

      if(zeon0(1,m).le.lmax.or.zeon0(1,m).ge.rmin) then
         instep =smoothstep(lmin,lmax,zeon(1,m))
         outstep=smoothstep(rmin,rmax,zeon(1,m))
         fip=instep*(1.0_wp-outstep)
         !  wdrive=wdrive*fip
         fdamp=(1.0_wp-fip)*fdamp0
      else
         fip=1.0_wp
         fdamp=0.0_wp
      end if

      wdphidt=wdphidt*fip

      !  zeon(5,m) = zeon00(8,m)+(zeon0(6,m)-nonlinear*zeon(8,m))*wdphidt*dtime*wfix  ! d(e*delta_phi/Te)/dt term added to weight,  but with delta_phi not advanced
      zeon(5,m) = zeon00(8,m)+(zeon0(6,m)-nonlinear*zeon(8,m))*wdphidt*dtime  ! d(e*delta_phi/Te)/dt term added to weight,  but with delta_phi not advanced

!       E_a=E_a+zeon(5,m)

  end do

  CALL CHARGE_EON
  CALL POISSON( solver1, masssolver )
  CALL SMOOTH_PHI

!$omp parallel do private(m,T_inverse,r,ip,fip,Ea_00, &
!$omp& e1,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,wdphidt)
  do m=1,me
     T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti

! calculate d(delta_phi)/dt at particles
      wz1=wzeon(m)
      wz0=1.0_wp-wz1

      ij=jteon0(m)
      wp0=1.0_wp-wpeon(m)
      wt00=1.0_wp-wteon0(m)
      e1=   wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=ij+1
      wt10=1.0_wp-wt00
      e1=e1+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=jteon1(m)
      wp1=1.0_wp-wp0
      wt01=1.0_wp-wteon1(m)
      e1=e1+wp1*wt01*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=ij+1
      wt11=1.0_wp-wt01
      e1=e1+wp1*wt11*(wz0*phi(0,ij)+wz1*phi(1,ij)) ! --> phi at electrons

      if(ZF_SEPARATION==0) then
!    calculate <phi> at electrons (no FLR effect -- x=R)
         r=zeon(7,m)-s0       ! in radial coordinator s
         r=r*delr
         ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
         fip=r-real(ip,wp)
         Ea_00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons
      end if

      wdphidt=(e1-real(1-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0  ! e*delta_phi/Te at electron new positions, but with advanced delta_phi
      wdphidt=(wdphidt-zeon0(8,m))/dtime ! d(e*delta_phi/Te)/dt with correction of time derivative of delta_phi

      !  if(zeon0(1,m).gt.psi_grid(mpsi-nbound-1))  then
      !     wfix=1.0_wp
      !  else
      !     wfix=1.0_wp
      !  endif

      if(zeon0(1,m).le.lmax.or.zeon0(1,m).ge.rmin) then
         instep =smoothstep(lmin,lmax,zeon(1,m))
         outstep=smoothstep(rmin,rmax,zeon(1,m))
         fip=instep*(1.0_wp-outstep)
         !  wdrive=wdrive*fip
         fdamp=(1.0_wp-fip)*fdamp0
      else
         fip=1.0_wp
         fdamp=0.0_wp
      end if

      wdphidt=wdphidt*fip

      !  zeon(5,m) = zeon00(8,m)+(zeon0(6,m)-nonlinear*zeon(8,m))*wdphidt*dtime*wfix ! d(e*delta_phi/Te)/dt term added to weight, with delta_phi advanced
      zeon(5,m) = zeon00(8,m)+(zeon0(6,m)-nonlinear*zeon(8,m))*wdphidt*dtime ! d(e*delta_phi/Te)/dt term added to weight, with delta_phi advanced

!       E_f=E_f+zeon(5,m)
  end do





if(ifelectron==-100) then    ! need to do one more time ?

  CALL CHARGE_EON
  CALL POISSON( solver1, masssolver )
  CALL SMOOTH_PHI

!$omp parallel do private(m,T_inverse,r,ip,fip,Ea_00, &
!$omp& e1,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,wdphidt)
  do m=1,me
     T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti

! calculate d(delta_phi)/dt at particles
      wz1=wzeon(m)
      wz0=1.0_wp-wz1

      ij=jteon0(m)
      wp0=1.0_wp-wpeon(m)
      wt00=1.0_wp-wteon0(m)
      e1=   wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=ij+1
      wt10=1.0_wp-wt00
      e1=e1+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=jteon1(m)
      wp1=1.0_wp-wp0
      wt01=1.0_wp-wteon1(m)
      e1=e1+wp1*wt01*(wz0*phi(0,ij)+wz1*phi(1,ij))

      ij=ij+1
      wt11=1.0_wp-wt01
      e1=e1+wp1*wt11*(wz0*phi(0,ij)+wz1*phi(1,ij)) ! --> phi at electrons

      if(ZF_SEPARATION==0) then
!    calculate <phi> at electrons (no FLR effect -- x=R)
         r=zeon(7,m)-s0       ! in radial coordinator s
         r=r*delr
         ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
         fip=r-real(ip,wp)
         Ea_00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons
      end if

      wdphidt=(e1-real(1-ZF_SEPARATION,wp)*Ea_00)*T_inverse*csc0  ! e*delta_phi/Te at electron new positions, with one more interation for advancing delta_phi
      wdphidt=(wdphidt-zeon0(8,m))/dtime ! d(e*delta_phi/Te)/dt with one more interation for advancing delta_phi 

      zeon(5,m) = zeon00(8,m)+(zeon0(6,m)-nonlinear*zeon(8,m))*wdphidt*dtime  ! d(e*delta_phi/Te)/dt term added to weight 

!       E_t=E_t+zeon(5,m)
  end do

end if




!   if(mype==0) print*, 'BBBB=', E_a, '  SSSS=', E_f, '  ZZZZ=', E_t


 if(irk_e==2)then
! out of boundary particle
!$omp parallel do private(m)
    do m=1,me
       if(zeon(1,m) > psimax)then
          zeon(1,m)=zeon0(1,m)
          zeon(2,m)=2.0_wp*pi-zeon0(2,m)
          zeon(3,m)=zeon0(3,m)
          zeon(4,m)=zeon0(4,m)
          zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
          !   zeon(5,m)=-zeon0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
          !   zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
          zeon(7,m)=zeon0(7,m)

       elseif(zeon(1,m) < psimin)then
          zeon(1,m)=zeon0(1,m)
          zeon(2,m)=2.0_wp*pi-zeon0(2,m)
          zeon(3,m)=zeon0(3,m)
          zeon(4,m)=zeon0(4,m)
          zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
          !   zeon(5,m)=-zeon0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
          !   zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
          zeon(7,m)=zeon0(7,m)
       endif
    enddo

! ********************************
if(mype==2) then
  error_8=1.0e-5_wp
  IFAIL=1
  rdum=D01AHF(0.0d0,zeon(1,me/2),error_8,minpts,acc,dsda,-1,IFAIL)
  print*, 'EEEEEEEEEEEE', zeon(7,me/2),'  ', rdum
!   write(stdout,*)'######## electron gcmotion: wall =',(tt1-tt0),'  cpu =',(cput1-cput0)
end if
! ********************************

 endif


ELSE    !  for df electrons, adiabatic electrons, and delta-h electrons with dphi/dt equation




! update GC position
! Another OpenMP parallel loop...
!$omp parallel do private(m,mu,energy,kappa,wdrift,weqflow,vdr,wdrive,wdot,xdum,ydum,T_inverse,r,ip,Ea_00, &
!$omp& e1,e2,e3,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11,wdphidt,wdphidt2,muB,fip,instep,&
!$omp& outstep,fdamp,vpara,v2,avw)
  do m=1,me
!     dbdz=0.0   for axisymmetric system   !wx

     mu=zeon(6,m)*zeon(6,m)        ! =normalized mu of electrons
     muB=mu*BB(m)                  ! =normalized perpendicular energy of electrons
     T_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti
!     energy=0.5_wp*for_E*(zeon(4,m)*BB(m))**2+mu*BB(m)

     vpara=zeon(4,m)*BB(m)
     v2=vpara**2

!     energy=0.5_wp*for_E*(zeon(4,m)*BB(m))**2+muB
     energy=0.5_wp*for_E*v2+muB

     vpara=vpara*f_rho                ! =normalized v_//
     v2=vpara*vpara+muB*2.0_wp*m0ms   ! =v^2

     kappa=(min(umax*umax,energy*T_inverse)-1.5_wp)*(-dlnTdae(m))+(-dlnndae(m)) !=-dlnf_M/da

     if(ifneoclassical==1) then
        vdr=da(m)
        weqflow=vdr*Ba(m)+dgq(m)*Bq(m)
     else
        !  ExB drift in radial direction for w-dot and flux diagnostics
        vdr=ghq(m)             !  da/dt|_EXB
        weqflow=vdr*Ba(m)+Ta(m)*Bq(m)
     end if

     wdrive=vdr*kappa


     if(ifelectron==2) then 

! calculate d(delta_phi)/dt at particles
        wz1=wzeon(m)
        wz0=1.0_wp-wz1

        ij=jteon0(m)
        wp0=1.0_wp-wpeon(m)
        wt00=1.0_wp-wteon0(m)
        e1=   wp0*wt00*(wz0*dphidt(0,ij)+wz1*dphidt(1,ij))

        ij=ij+1
        wt10=1.0_wp-wt00
        e1=e1+wp0*wt10*(wz0*dphidt(0,ij)+wz1*dphidt(1,ij))
  
        ij=jteon1(m)
        wp1=1.0_wp-wp0
        wt01=1.0_wp-wteon1(m)
        e1=e1+wp1*wt01*(wz0*dphidt(0,ij)+wz1*dphidt(1,ij))

        ij=ij+1
        wt11=1.0_wp-wt01
        e1=e1+wp1*wt11*(wz0*dphidt(0,ij)+wz1*dphidt(1,ij)) ! --> dphidt at electrons
  
        wdphidt=e1*csc0*T_inverse    ! when e1 is dphi/dt at electrons -- using FEM solver with order>0
        !  wdphidt=e1*csc0           ! when dphi/dt equation is modified; dphidt=[d(delta_phi)/dt]/Te using old integral solver with order=0

!    calculate d<phi>/da at electrons (no FLR effect -- x=R)
        r=zeon(7,m)-s0       ! in radial coordinator s
        r=r*delr
        ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
        fip=r-real(ip,wp)
        Ea_00=fip*phip00(ip+1)+(1.0_wp-fip)*phip00(ip) ! d<phi>/da at electrons

        wdphidt2=da(m)*(wpi(1,m)-Ea_00)+dgq(m)*wpi(2,m)+dgf(m)*wpi(3,m)
        wdphidt2=wdphidt2*T_inverse*csc0
        wdphidt=wdphidt+wdphidt2          ! d/dt(delta-phi) term in RHS of dhe/dt equation

     else                  ! delta-f electrons
        wdphidt=0.0_wp
     end if

     if(zeon0(1,m).le.lmax.or.zeon0(1,m).ge.rmin) then
        instep =smoothstep(lmin,lmax,zeon(1,m))
        outstep=smoothstep(rmin,rmax,zeon(1,m))
        !  fip=smoothstep(0.20_wp,0.25_wp,zeon(1,m))*(1.0_wp-smoothstep(0.75_wp,0.80_wp,zeon(1,m)))
        fip=instep*(1.0_wp-outstep)
        !  wdrive=wdrive*fip
        !  wdphidt=wdphidt*fip

        fdamp=(1.0_wp-fip)*fdamp0
     else
        fip=1.0_wp
        fdamp=0.0_wp
     end if

     if(ifneoclassical==1) then
        wdrift=(dgr(m)*BB(m)+zeon(4,m)*weqflow)*zeon(4,m)*BB(m)*for_E*T_inverse   ! see note
     else
        wdrift=(Pa(m)*BB(m)+zeon(4,m)*weqflow)*zeon(4,m)*BB(m)*for_E*T_inverse   ! see note
     endif

     weqflow=weqflow*mu*T_inverse   ! see note: additional term due to v_ExB*grad(f_M)

     wdot=(wdrive+wdrift+weqflow+wdphidt)

     wdot=wdot*fip


!!!!!  if(nonlinear>0.5_wp) then     ! for nonlinear run
!!!!!  if(nonlinear>100.0_wp) then        ! for nonlinear run
if(wimcompressibility==1) then     ! for weight imcompressibility
! calculate average weight <w> at particles using wieght-field wfield
     wz1=wzeon(m)
     wz0=1.0_wp-wz1

     ij=jteon0(m)
     wp0=1.0_wp-wpeon(m)

     wt00=1.0_wp-wteon0(m)
     e1=wp0*wt00*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=wp0*wt00*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=wp0*wt00*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt10=1.0_wp-wt00
     e1=e1+wp0*wt10*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp0*wt10*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp0*wt10*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=jteon1(m)
     wp1=1.0_wp-wp0
     wt01=1.0_wp-wteon1(m)
     e1=e1+wp1*wt01*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt01*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt01*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     ij=ij+1
     wt11=1.0_wp-wt01
     e1=e1+wp1*wt11*(wz0*wfield(1,ij,0)+wz1*wfield(1,ij,1))
     e2=e2+wp1*wt11*(wz0*wfield(2,ij,0)+wz1*wfield(2,ij,1))
     e3=e3+wp1*wt11*(wz0*wfield(3,ij,0)+wz1*wfield(3,ij,1))

     avw=e1+e2*vpara+e3*v2

!     if(mype==100.and.m==me/2) print*, 'AVW=', avw, zeon(5,m), mype
!     if(mype==0.and.m==me/4) print*, 'AVW=', avw, zeon(5,m), mype
!     if(mype==200.and.m==me/8) print*, 'AVW=', avw, zeon(5,m), mype

end if


     zeon(1,m) = zeon0(1,m)+dtime*da(m)
     zeon(2,m) = zeon0(2,m)+dtime*dgq(m)
     zeon(3,m) = zeon0(3,m)+dtime*dgf(m)
     zeon(4,m) = zeon0(4,m)+dtime*dgr(m)

!     zeon(5,m) = zeon0(5,m)+dtime*wdot

!     zeon(5,m) = 1.0_wp-(1.0_wp-zeon0(5,m))*exp(-dtime*wdot)

      if(wimcompressibility==1) then     ! for weight imcompressibility
         zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*avw)*wdot-fdamp*zeon(5,m))*dtime
      else
         zeon(5,m) = zeon0(5,m)+((zeon0(6,m)-nonlinear*zeon(5,m))*wdot-fdamp*zeon(5,m))*dtime
      end if

!     zeon(7,m) = zeon0(7,m)+dtime*sqrt(Tc/backT(m))*da(m)
!     zeon(7,m) = zeon0(7,m)+dtime*(Tc/backT(m))**0.4_wp*da(m)
     zeon(7,m) = zeon0(7,m)+dtime*(Tc/backT(m))**dsda_power*da(m)

#ifdef _TREAT_MAGNETIC_AXIS
     if(sw(m).eq.1) then    ! convert to (a, theta) coordinates
        xdum=sqrt(zeon(1,m)**2+zeon(2,m)**2)
        ydum=sign(1.0_wp,zeon(2,m))*acos(zeon(1,m)/xdum)
        zeon(1,m)=xdum
        zeon(2,m)=ydum
     end if
#endif
     
! theta and zeta normalize to [0,2*pi), modulo is slower than hand coded procedure
!     zeon(2,m)=zeon(2,m)*pi2_inv+100.0 !period of 1
!     zeon(2,m)=2.0*pi*(zeon(2,m)-aint(zeon(2,m))) ![0,2*pi)
!     zeon(3,m)=zeon(3,m)*pi2_inv+100.0
!     zeon(3,m)=2.0*pi*(zeon(3,m)-aint(zeon(3,m)))

!     zeon(2,m)=zeon(2,m)*pi2_inv !period of 1
!     zeon(2,m)=10.0+zeon(2,m)-dint(zeon(2,m)) ![0,2*pi)
!     zeon(2,m)=twopi*(zeon(2,m)-dint(zeon(2,m)))
!     zeon(3,m)=zeon(3,m)*pi2_inv
!     zeon(3,m)=10.0+zeon(3,m)-dint(zeon(3,m))
!     zeon(3,m)=twopi*(zeon(3,m)-dint(zeon(3,m)))

! use f90 modulo to normalize theta and zeta to [0,2*pi)
! If P > 0, then the function modulo(A,P) determines the positive difference between A and then
! next lowest multiple of P. If P < 0, then the function modulo(A,P) determines the negative
! difference between A and then next highest multiple of P.
     zeon(2,m)=modulo(zeon(2,m),twopi)
     zeon(3,m)=modulo(zeon(3,m),twopi)


! store GC information for flux measurements
     wpi(1,m)=ghq(m)    ! at present, wpi(1,*)=da/dt|_EXB 
     wpi(2,m)=energy    ! normalized electron energy with u//=0 in Jan.09
     wpi(3,m)=muB       ! normalizd electron perpendicular energy

  enddo

 if(irk_e==2)then
! out of boundary particle
!$omp parallel do private(m)
    do m=1,me
       if(zeon(1,m) > psimax)then
          zeon(1,m)=zeon0(1,m)
          zeon(2,m)=2.0_wp*pi-zeon0(2,m)
          zeon(3,m)=zeon0(3,m)
          zeon(4,m)=zeon0(4,m)
          zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
!          zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
          zeon(7,m)=zeon0(7,m)
          
       elseif(zeon(1,m) < psimin)then
          zeon(1,m)=zeon0(1,m)
          zeon(2,m)=2.0_wp*pi-zeon0(2,m)
          zeon(3,m)=zeon0(3,m)
          zeon(4,m)=zeon0(4,m) 
          zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
!          zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
          zeon(7,m)=zeon0(7,m)
       endif
    enddo

! ********************************
if(mype==2) then
  error_8=1.0e-5_wp
  IFAIL=1
  rdum=D01AHF(0.0d0,zeon(1,me/2),error_8,minpts,acc,dsda,-1,IFAIL)
  print*, 'EEEEEEEEEEEE', zeon(7,me/2),'  ', rdum
!   write(stdout,*)'######## electron gcmotion: wall =',(tt1-tt0),'  cpu =',(cput1-cput0)
end if
! ********************************

 endif

 if(idiag_e==0.and.irk_e==1) then  ! istep=n*ndiag, irk=1, i=1, and irk_e=1
! fluxes diagnose
!$omp parallel do private(m)
    do m=1,me
       wpi(1,m)=wpi(1,m)*gha(m)       ! now wpi(1,*)=dr/dt|_EXB=da/dt|_EXB *(1/|grad a|)
    end do

    call get_fluxes_e(wpi)

!   im=mod(istep,100)
   im=mod(istep,msnap)
   if(im==0) then
      call get_Te_fluctuation(wpi(2,:))
   end if

 endif

END IF

end subroutine pusheon




subroutine electron_locator
  use precision
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer m,im,ii,j01,j00
  real(wp) rdum,tdum,delr,delt(0:mpsi),delz,wz1,wp1,     &
       pi2_inv,psitmp,thetatmp,zetatmp,psimax,psimin


  delr=1.0_wp/deltar 
  delt=2.0_wp*pi/deltat
  delz=1.0_wp/deltaz
  pi2_inv=0.5_wp/pi
  psimin=s_grid(0)
  psimax=s_grid(mpsi)

!$omp parallel do private(m,psitmp,thetatmp,zetatmp,&
!$omp& rdum,ii,wp1,im,tdum,j00,j01)
  do m=1,me
!     psitmp=zeon(1,m)      ! =a  
     psitmp=zeon(7,m)       ! =s
     thetatmp=zeon(2,m)
     zetatmp=zeon(3,m)

     wzeon(m)=(zetatmp-zetamin)*delz

     rdum=psitmp
     rdum=delr*max(0.0_wp,min(s1-s0,rdum-s0))
     ii=max(0,min(mpsi-1,int(rdum)))
     wp1=rdum-real(ii,wp)
     wpeon(m)=wp1

! inner flux surface
     im=ii
     tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
     tdum=(tdum-aint(tdum))*delt(im)
     j00=max(0,min(mtheta(im)-1,int(tdum)))
     jteon0(m)=igrid(im)+j00  ! pay attention to j00=0 case
     wteon0(m)=tdum-real(j00,wp)

! outer flux surface
     im=ii+1
     tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
     tdum=(tdum-aint(tdum))*delt(im)
     j01=max(0,min(mtheta(im)-1,int(tdum)))
     jteon1(m)=igrid(im)+j01  ! pay attention to j01=0 case
     wteon1(m)=tdum-real(j01,wp)
  enddo	

end subroutine electron_locator
