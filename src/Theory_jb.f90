module data_pass
  real*8 dlamda
  integer nbad
  real*8, dimension(:), allocatable ::  bad
end module data_pass




subroutine theory
 
!      Evaluation of neoclassical theory results of fluxes etc
 
  use global_parameters
  use diagnosis_array
  use particle_array
  use field_array
  use esi_interface

  use data_pass

  implicit none
 
  real*8, dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the
  real*8  dlnT, dlnT_e, dlnn, vthr, vthri, ep, qq, B0B2IN, &
          density_ei, T_elec, T_ion, B_s, B_max, f_t, adum, up, lamda1,       &
          ajavbs,ajavbs_e,tempxe,tempxi,pvale,pvali,ppvale,ppvali,xnuse,xnusi,tempxpe,tempxpi, &
          rm2avx,ftrapx,epsx,rmajx,zeff,qx,mu0, &
          dr, r,  R2, eps, acc, ierr


  real*8, external :: temperature, temperature_e, density, dT_da, dT_da_e, dn_da, &
                      trap_fraction, D01AHF, trap,   R2_inverse, object_5
  integer i,   ifail,minpts, nout,iboot

  common/stupidpass/adum
  common/sendtotrap/B_s,lamda1
      
      nbad=10000
      allocate(bad(0:nbad))

!       dr=(a1-a0)/real(nbin_a,wp)
!       do 19 i=1, nbin_a
!       dr=(a1-a0)/real(200)
      do 19 i=0, mpsi
!          r=a0+dr*(dfloat(i-1)+0.5d0)
         r=psi_grid(i)

         density_ei=density(r,0.0_wp)
         T_elec=temperature_e(r,0.0_wp)
         T_ion=temperature(r,0.0_wp)
         vthr=sqrt(2.0d0*T_elec)  ! normalized elec thermal velocity
         vthri=sqrt(2.0d0*T_ion)  ! normalized ion thermal velocity

         aaa(1)=r
         the(1)=0.0d0
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
   
         R2=rj(1)
         qq=-gFaj(1)/gYaj(1)
         the(1)=pi
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
         ep=abs(R2-rj(1))/R_0/2.0d0                     ! inverse aspect ratio
         
         epsx=(abs(R2-rj(1))/2.0d0)/((R2+rj(1))/2.0d0)  ! inverse aspect ratio
         rmajx=(R2+rj(1))/2.0d0*(R_00/1.0d2)            ! major radius in meter

         B_max=Bj(1)

         dlnT_e=dT_da_e(r,0.0d0)/T_elec
         dlnn=dn_da(r,0.0d0)/density_ei
         dlnT=dT_da(r,0.0d0)/T_ion

         call general_flux_average(r,object_5,B0B2IN) ! <B^2>
         B_s=sqrt(B0B2IN)

         adum=r
         up=B_s/B_max

         dlamda=up/real(nbad)
         do IFAIL=1,nbad
            lamda1=dlamda*dfloat(IFAIL) 
            call general_flux_average(adum,trap,bad(IFAIL))  ! bad=<sqrt(1-lamda/h)>
         end do
         bad(0)=1.0d0
   
         eps=1.0d-5      !  eps=1.0d-7
         IFAIL=1
         f_t=D01AHF(0.0d0,up,eps,minpts,acc,trap_fraction,-1,IFAIL)
         f_t=1.0d0-3.0d0/4.0d0*f_t   ! trapped particle fraction

         if(mype==0) print*,'trapped particle fraction=', f_t, '1.46*sqrt(ep)=', 1.46d0*sqrt(ep), IFAIL

!     prepare input for subroutine bootstrap
         ftrapx=f_t
         call general_flux_average(adum,R2_inverse,rm2avx)  ! rm2avx=<1/R^2>
         rm2avx=rm2avx/(R_00/1.0d2)**2  ! convert to unit of 1/meter^2
         mu0=4.0d0*pi*1.0d-7  ! permeability
         iboot=2
         nout=80
         zeff=1.0d0
         ftrapx=f_t
         qx=qq
         pvale=density_ei*1.0d6*T_elec*E_00*msm0/1.0d7 ! in unit Newton/meter^2
            ! energy associated with normalized energy is T_elec*(m/m_0) in unit E_00
         pvali=density_ei*1.0d6*T_ion*E_00*msm0_i/1.0d7 ! in unit Newton/meter^2
         pvale=pvale*mu0   ! electron temperature time mu0
         pvali=pvali*mu0   ! ion temperature time mu0
         ppvale=pvale*(dlnT_e+dlnn)/(-gYaj(1)*(B_00/1.0d4)*(R_00/1.0d2)**2)  !dp_e/dpsi
         ppvali=pvali*(dlnT + dlnn)/(-gYaj(1)*(B_00/1.0d4)*(R_00/1.0d2)**2)  !dp_i/dpsi
         tempxe=T_elec*msm0*Tem_00    ! T_e in kev
         tempxi=T_ion*msm0_i*Tem_00   ! T_i in kev
         tempxpe=tempxe*dlnT_e/(-gYaj(1)*(B_00/1.0d4)*(R_00/1.0d2)**2)  ! dT_e/dpsi
         tempxpi=tempxi*dlnT/(-gYaj(1)*(B_00/1.0d4)*(R_00/1.0d2)**2) ! dT_i/dpsi
!     estimate j_b in terms of Hirshman's formula         
         call bootstrap(ajavbs, ajavbs_e, iboot, ierr, nout,                 &
                        zeff,rm2avx,epsx,rmajx,qx,ftrapx,          &
                        pvale,ppvale,tempxe,tempxpe,               &
                        pvali,ppvali,tempxi,tempxpi,xnuse,xnusi)
         ajavbs=ajavbs/mu0
         ajavbs=ajavbs*Fj(1)*rm2avx*utime/(charg0/3.0d9)/1.0d6   ! in unit e0*v_00*B_00

!      electron BS current added 2017
         ajavbs_e=ajavbs_e/mu0
         ajavbs_e=ajavbs_e*Fj(1)*rm2avx*utime/(charg0/3.0d9)/1.0d6   ! in unit e0*v_00*B_00

!         write(21,1299) r, ajavbs, pvale, pvali, ppvale, ppvali, tempxe, tempxi, tempxpe, tempxpi
         write(21,1299) r, ajavbs, ajavbs_e, xnusi, xnuse
1299     format(10(2x,ES13.6))

19    continue
!         call flush_(21)
         close(21)
      
end subroutine theory



real*8 FUNCTION trap_fraction(lamda)

!  provide function  lamda/<sqrt(1-lamda/h)>

      use data_pass

      implicit  none
      real*8 lamda,xx, ave_trap
      integer i

      xx=lamda/dlamda
      i=int(xx)
      ave_trap=(dfloat(i+1)-xx)*bad(i)+(xx-dfloat(i))*bad(i+1)
                        ! ave_trap=<sqrt(1-lamda/h)>
      trap_fraction=lamda/ave_trap

      RETURN
   end function trap_fraction


real*8 FUNCTION trap(thedum)

!  provide function  J*sqrt(1-lamda/h)

      use esi_interface

      use data_pass

      implicit  none

      real*8, external:: reducedjacobian
      real*8 adum,thedum,udum,vth2dum,aaa(2),the(2),grad_a2
      real*8, dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      real*8 B_s,lamda1
      common/stupidpass/adum
      common/sendtotrap/B_s,lamda1

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      trap=reducedjacobian(thedum)*sqrt(1.0d0-lamda1/(B_s/Bj(1)))

   end function trap


real*8 FUNCTION R2_inverse(thedum)

!  provide function  J*(1/R^2)

      use esi_interface

      implicit  none

      real*8, external:: reducedjacobian
      real*8 adum,thedum,udum,vth2dum,aaa(2),the(2),grad_a2
      real*8, dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      R2_inverse=reducedjacobian(thedum)*(1.0d0/rj(1)**2)

   end function R2_inverse




!
!    analytical estimate of bootstrap current in term of Hirshman's formula
!                       by Charles Kessel   June 2003
!
      subroutine bootstrap(ajavbs, ajavbs_e, iboot, ierr, nout, &
                           zeff,rm2avx,epsx,rmajx,qx,ftrapx,          &
                           pvale,ppvale,tempxe,tempxpe,               &
                           pvali,ppvali,tempxi,tempxpi,xnuse,xnusi)
!
!     DESCRIPTION OF RETURNED PARAMETERS:
!     ajavbs   surface averaged parallel current due to bootstrap
!              divided by surface averaged B dot grad phi
!              (divide by mu0=4*pi*1.e-7 to get amperes/meter)
!
!              ajavbs = <mu_0*j_//*B>/<B.gard(phi)>  ! noted by wx
!
!     ajavbs_e electron BS current  added 2017
!   
!     ierr     =0 for normal calculation
!              =1 for error exit
!
!     DESCRIPTION OF INPUT PARAMETERS:
!     iboot    =1 for collisionless bootstrap calculation
!              =2 to include collisonal corrections
!     nout     logical unit for error output
!
!     zeff     effective charge
!     rm2avx   surface average of 1/R**2
!     epsx     inverse aspect ratio of this flux surface
!     rmajx    major radius of this flux surface (m)
!     qx       safety factor
!     ftrapx   trapped particle fraction
!     pvale    electron pressure times mu0
!     ppvale   derivative of electron pressure wrt psi  (PF/radian)
!     pvali    ion pressure times mu0
!     ppvali   derivative of ion pressure wrt psi       (PF/radian)
!     tempxe   electron temperature (kev)
!     tempxpe  derivative of electron pressure wrt psi  (PF/radian)
!     tempxi   ion temperature      (kev)
!     tempxpi  derivative of ion temperature wrt psi    (PF/radian)
!
!

      use global_parameters, only: coll_enlarge

      implicit real*8(a-h,o-z)

      data amu0/1.256637061e-6/
      data pi/3.1415926535/
      data echarge/1.6022e-19/
      data emass/9.1095e-31/
      data pmass/1.6726e-27/
      data xeps0/8.854e-12/
 
      ierr = 0
      if(tempxe.le.0 .or. tempxi.le.0) go to 101
      if(ftrapx.le.0 .or. ftrapx.ge.1) go to 102
      if(pvale .le.0 .or. pvali .le.0) go to 103
      if(rm2avx .le. 0) go to 104
      if(zeff .lt. 1.0) go to 105
!     if(iabs(iboot).ne.1 .or. iabs(iboot).ne.2
!    + .or. iabs(iboot) .ne. 4) go to 106
      if(iabs(iboot) .lt. 1 .or. iabs(iboot) .gt. 4) go to 106
      zave = (pvale/tempxe)/(pvali/tempxi)
      xf=ftrapx/(1.0-ftrapx)
      alfi = 1.172/(1.0+0.462*xf)
      denom = 1.414*zeff + zeff**2 + xf*(0.754+2.657*zeff+2.*zeff**2)    &
                                + xf**2*(0.348+1.243*zeff+   zeff**2)
      anum1 = xf*(0.754 + 2.21*zeff + zeff**2    &
            + xf*(0.348 + 1.24*zeff + zeff**2) )
      anum2 = xf*(0.884 + 2.07*zeff)
      pepope=ppvale/pvale
      pipopi=ppvali/pvali
      tepote=tempxpe/tempxe
      tipoti=tempxpi/tempxi
      if(abs(iboot) .eq. 1) then
      ajavbs=-(pvale/rm2avx)*(anum1*(pepope+(pvali/pvale)         &
                *(pipopi - alfi*tipoti) ) - anum2*tepote) / denom
      return
      endif
 
      if(zeff .le. 2.0) then
      a13=1.02-0.13*(zeff-1.0)
      b13=1.07-0.45*(zeff-1.0)
      c13=1.07-0.38*(zeff-1.0)
      a23=0.57-0.05*(zeff-1.0)
      b23=0.61-0.27*(zeff-1.0)
      c23=0.61-0.23*(zeff-1.0)
      else
      a13=0.89-0.05*(zeff-2.0)
      b13=0.62-0.03*(zeff-2.0)
      c13=0.69-0.09*(zeff-2.0)
      a23=0.52-0.02*(zeff-2.0)
      b23=0.34-0.005*(zeff-2.0)
      c23=0.38-0.05*(zeff-2.0)
      endif
      xmi=1.0
      xclog=16.0
      xdense=pvale/(amu0*tempxe*echarge*1.e3)
      xdensi=pvali/(amu0*tempxi*echarge*1.e3)
      xvthe=sqrt(2.0*echarge*1.e3*tempxe/emass)
      xvthi=sqrt(2.0*echarge*1.e3*tempxi/(xmi*pmass))
      xtaue=(12.0*pi**1.5*xeps0**2*(emass)**2*xvthe**3)/  &
            (4.0*xdense*zeff*(echarge)**4*xclog)
      xtaui=1.414*(12.0*pi**1.5*xeps0**2*(xmi*pmass)**2*  &
            xvthi**3)/(4.0*xdensi*(zave*echarge)**4*xclog)

      !  electron collision frequency is scaled by a factor of coll_enlarge -- 2016
      xtaue=xtaue/coll_enlarge

      xnuse=1.414*(rmajx*qx)/(xtaue*xvthe*epsx**1.5)
      xnusi=1.414*(rmajx*qx)/(xtaui*xvthi*epsx**1.5)
      xl31=anum1/denom
      xl32=anum2/denom
      col1=(1.0/(1.0+a13*sqrt(xnuse)+b13*xnuse))*  &
           (1.0/(1.0+c13*xnuse*epsx**1.5))
      col2=(1.0/(1.0+a23*sqrt(xnuse)+b23*xnuse))*  &
           (1.0/(1.0+c23*xnuse*epsx**1.5))
      alfhh=((alfi-0.35*sqrt(xnusi))/(1.0+0.7*sqrt(xnusi))  &
            -2.1*xnusi**2*epsx**3)*(1.0/(1.0+xnusi**2*epsx**3))*   &
            (1.0/(1.0+xnuse**2*epsx**3))
      if(iabs(iboot) .eq. 2) then
      ajavbs=-(pvale/rm2avx)*(xl31*col1*(pepope+(pvali/pvale)   &
             *(pipopi-alfhh*tipoti))-(2.5*xl31*col1-(2.5*xl31-xl32)*col2)*tepote)
      ajavbs_e=-(pvale/rm2avx)*(xl31*col1*(pepope)   &
             -(2.5*xl31*col1-(2.5*xl31-xl32)*col2)*tepote)
      return
      endif
 
 
 
      ft33=ftrapx/(1.+(0.55-0.1*ftrapx)*sqrt(xnuse)+0.45*   &
           (((1.-ftrapx)*xnuse)/zeff**1.5))
 
      ft31=ftrapx/(1.+(1.-0.1*ftrapx)*sqrt(xnuse)+0.5*      &
           (((1.-ftrapx)*xnuse)/zeff))
 
      ft32ee=ftrapx/(1.+0.26*(1.-ftrapx)*sqrt(xnuse)+0.18*  &
             (((1.-0.37*ftrapx)*xnuse)/sqrt(zeff)))
 
      ft32ei=ftrapx/(1.+(1.+0.6*ftrapx)*sqrt(xnuse)+0.85*   &
             (((1.-0.37*ftrapx)*xnuse)*(1.+zeff)))
 
      ft34=ftrapx/(1.+(1.-0.1*ftrapx)*sqrt(xnuse)+0.51*     &
           (((1.-0.5*ftrapx)*xnuse)/zeff))
 
      alp0=-(1.17*(1.-ftrapx))/(1.-0.22*ftrapx-0.19*ftrapx**2)
 
      alp=(((alp0+0.25*(1.-ftrapx**2)*sqrt(xnusi))/         &
          (1.+0.5*sqrt(xnusi)))                             &
          - 0.315*xnusi**2*ftrapx**6)*(1./(1.+0.15*xnusi**2*ftrapx**6))
 
      xl31=(1.+(1.4/(1.+zeff)))*ft31-(1.9/(1.+zeff))*ft31**2   &
           +(0.3/(1.+zeff))*ft31**3+(0.2/(1.+zeff))*ft31**4
 
      fl32ee=(((0.05+0.62*zeff)/(zeff*(1.+0.44*zeff)))*      &
             (ft32ee-ft32ee**4)                              &
             +(1./(1.+0.22*zeff))*(ft32ee**2-ft32ee**4-1.2*(ft32ee**3  &
             -ft32ee**4))                                    &
             +(1.2/(1.+0.5*zeff))*ft32ee**4)
 
      fl32ei=((-(0.56+1.93*zeff)/(zeff*(1.+0.44*zeff)))*(ft32ei  &
             -ft32ei**4)   &
             +(4.95/(1.+2.48*zeff))*(ft32ei**2-ft32ei**4-0.55*(ft32ei**3  &
             -ft32ei**4))  &
             -(1.2/(1.+0.5*zeff))*ft32ei**4)
 
      xl32=fl32ee+fl32ei
 
      xl34=(1.+(1.4/(1.+zeff)))*ft34-(1.9/(1.+zeff))*ft34**2   &
           +(0.3/(1.+zeff))*ft34**3+(0.2/(1.+zeff))*ft34**4
 
      ajavbs=-(pvale/rm2avx)*(xl31*(pepope+(pvali/pvale)*pipopi)+  &
             xl32*tepote+xl34*alp*tipoti)
      ajavbs_e=-(pvale/rm2avx)*(xl31*(pepope)+  &
             xl32*tepote)
      return
 
!.....error exit
  101 continue
      write(nout,1101) tempxe,tempxi
 1101 format(" Error in bootstrap, tempxe,tempxi=",1p2e12.4)
      go to 200
  102 continue
      write(nout,1102) ftrapx
 1102 format(" Error in bootstrap, ftrapx=",1pe12.4)
      go to 200
  103 continue
      write(nout,1103) pvale,pvali
 1103 format(" Error in bootstrap, pvale,pvali =",1p2e12.4)
      go to 200
  104 continue
      write(nout,1104) rm2avx
 1104 format(" Error in bootstrap, rm2avx=",1pe12.4)
      go to 200
  105 continue
      write(nout,1105) zeff
 1105 format(" Error in bootstrap, zeff=",1pe12.4)
      go to 200
  106 continue
      write(nout,1106) iboot
 1106 format(" Error in bootstrap, iboot=",i5)
 
  200 ierr = 1
      return
      end

