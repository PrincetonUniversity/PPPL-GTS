function metric_elements(pdum,tdum,zdum)
  use precision
  use esi_interface

  implicit none
  integer metric_elements
  real(wp)  pdum,tdum,zdum,D
  real(wp), dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
  
  aaa(1)=pdum
  the(1)=tdum
  
  call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
  D=zaj(1)*rqj(1)-zqj(1)*raj(1)
  
  metric(1)=(rqj(1)**2+zqj(1)**2)/D**2   ! grada^2
  metric(2)=(zaj(1)**2+raj(1)**2)/D**2   ! gradtheta^2
  metric(3)=1.0_wp/rj(1)**2                 ! gradphi^2
  metric(4)=(-raj(1)*rqj(1)-zaj(1)*zqj(1))/D**2  ! gradagradtheta
  metric(5)=0.0_wp                                  ! gradagradphi
  metric(6)=0.0_wp                                  ! gradthetagradphi

  metric_elements=0

  return

end function metric_elements


function f_jaco(ndim,xyz)
  use precision
  implicit none

  real(wp) f_jaco
  integer ndim
  real(wp) xyz(ndim)
  real(wp),external :: jacobian

     f_jaco=jacobian(xyz(1),xyz(2))

end function f_jaco



FUNCTION FUN1(ndim,xy)
      use precision
!
!      FUN1=<ni>*J =n_0*(An+Bn*exp(-dn*(a)**xn))*J
!
!      use esi_interface

      implicit  none
      real(wp) FUN1
      integer ndim
      real(wp) xy(ndim), density, jacobian
!      aaa(10),the(10)

!      aaa(1)=xy(1)
!      the(1)=xy(2)
!      call esiget2dfunctions(aaa,the,1)

!      FUN1=r(1)*(za(1)*rq(1)-zq(1)*ra(1))*density(xy(1),xy(2))
      FUN1=jacobian(xy(1),xy(2))*density(xy(1),xy(2))

      RETURN
end function FUN1 



FUNCTION FUN2(ndim,xy)
      use precision
!
!      FUN2=ni*J ( not <ni>*J )
!
!      use esi_interface

      implicit  none
      real(wp) FUN2
      integer ndim
      real(wp) xy(ndim),density,jacobian,u_parallel,temperature,surface_average,vary

      vary=exp(u_parallel(xy(1),xy(2))**2/(2.0_wp*temperature(xy(1),xy(2)))) &
           /surface_average(xy(1))
      FUN2=jacobian(xy(1),xy(2))*density(xy(1),xy(2))*vary

      RETURN
end function FUN2



FUNCTION FUN_J(ndim,xy)
      use precision
!
!      FUN_J=Jacobian for 2D integral 
!

      implicit  none
      real(wp) FUN_J
      integer ndim
      real(wp) xy(ndim), jacobian

      FUN_J=jacobian(xy(1),xy(2))

      RETURN
end function FUN_J




FUNCTION Jacobian(adum,thedum)
!
!   Jacobian=[grad(a) x grad(theta) . grad(phi)]^(-1)
!
      use precision
      use esi_interface

      implicit none
      real(wp) Jacobian
      real(wp) adum,thedum,aaa(10),the(10)
      real(wp), dimension(10) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj

      aaa(1)=adum
      the(1)=thedum
!      call esiget2dfunctions(aaa,the,1)  !replaced by:
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)

      Jacobian=rj(1)*(zaj(1)*rqj(1)-zqj(1)*raj(1))

      RETURN
end function Jacobian



FUNCTION reducedjacobian(thedum)
      use precision
      implicit  none
      real(wp) reducedjacobian
      real(wp) adum,thedum,jacobian
      common/stupidpass/adum

      reducedjacobian=Jacobian(adum,thedum)

      RETURN
end function reducedjacobian



FUNCTION fun_zero(thedum)
      use precision      
      implicit  none
      real(wp) fun_zero
      real(wp) adum,thedum,acc,eps,the_0,term_2,D01AHF
      real(wp), external:: reducedjacobian
      integer minpts,IFAIL
      common/for_the_search/the_0,term_2

      eps=1.0e-5_wp
      IFAIL=1
      fun_zero=D01AHF(the_0,thedum,eps,minpts,acc,reducedjacobian,-1,IFAIL)*  & 
               term_2-1.0_wp
      RETURN
end function fun_zero


FUNCTION distribution_para(vz)
      use precision
      implicit  none
      real(wp) distribution_para
      real(wp) vz

      distribution_para=exp(-vz*vz)

      RETURN
end function distribution_para


FUNCTION fun_para_zero(vz)
      use precision
      implicit  none
      real(wp) fun_para_zero
      real(wp) vz,acc,eps,low_para,S_para
      real(wp), external:: distribution_para,D01AHF
      integer minpts,IFAIL
      common/badpass/low_para,S_para

      eps=1.0e-5_wp
      IFAIL=1
      fun_para_zero=D01AHF(low_para,vz,eps,minpts,acc,distribution_para,-1,&
            IFAIL)*S_para-1.0_wp
      RETURN
end function fun_para_zero


FUNCTION fun_perp_zero(v_perp)
      use precision
      implicit  none
      real(wp) fun_perp_zero
      real(wp) v_perp,low_perp,S_perp
      common/badpass1/low_perp,S_perp

      fun_perp_zero=(exp(-low_perp*low_perp)-exp(-v_perp*v_perp))*S_perp-1.0_wp
      RETURN
end function fun_perp_zero


FUNCTION rho2a_zero(adum)
      use precision
      implicit  none
      real(wp) rho2a_zero
      real(wp) adum,Tcc,r
      real(wp), external:: temperature
      common/fun1d/Tcc,r

      rho2a_zero=adum**2*Tcc/temperature(adum,0.0_wp)-r**2

      RETURN
end function rho2a_zero


FUNCTION dsda(adum)
      use precision
      use global_parameters,only: dsda_power
      implicit  none
      real(wp) dsda
      real(wp) adum,Tcc,r
      real(wp), external:: temperature
      common/fun1d/Tcc,r

!      dsda=sqrt(Tcc/temperature(adum,0.0))
!      dsda=(Tcc/temperature(adum,0.0_wp))**0.4_wp
      dsda=(Tcc/temperature(adum,0.0_wp))**dsda_power

      RETURN
end function dsda


FUNCTION s2a_zero(adum)
      use precision
      implicit  none
      real(wp) s2a_zero
      real(wp) adum,Tcc,r,error_8,acc
      real(wp), external:: D01AHF,dsda
      integer minpts,IFAIL
      common/fun1d/Tcc,r

      error_8=1.0e-5_wp
      IFAIL=1

      s2a_zero=D01AHF(0.0_wp,adum,error_8,minpts,acc,dsda,-1,IFAIL)-r

      RETURN
end function s2a_zero


FUNCTION surface_average(adum)
      use precision
!
!    calculate <exp(u^2/v_th^2)>
!
      use global_parameters

      implicit  none
      real(wp) surface_average

      real(wp) adum,acc,eps,D01AHF,adum_share
      real(wp), external:: reducedjacobian, density_vary
      integer minpts,IFAIL
      real(wp) j_average,jA_average
      common/stupidpass/adum_share

      adum_share=adum
      eps=1.0e-5_wp
      IFAIL=1
      j_average=D01AHF(0.0_wp,2.0_wp*pi,eps,minpts,acc,reducedjacobian,-1,IFAIL)
      eps=1.0e-5_wp
      IFAIL=1
      jA_average=D01AHF(0.0_wp,2.0_wp*pi,eps,minpts,acc,density_vary,-1,IFAIL)
      surface_average=jA_average/j_average

      RETURN
end function surface_average 



FUNCTION density_vary(thedum) 
      use precision

!  provide function    density_vary(aa,theta)=Jacobian*exp(u^2/v_th^2) 

      implicit  none
      real(wp) density_vary

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum
      common/stupidpass/adum
      
      udum=u_parallel(adum,thedum)
      vth2dum=2.0_wp*temperature(adum,thedum)  ! v_th^2
      density_vary=reducedjacobian(thedum)*exp(udum*udum/vth2dum)

      return
end function density_vary



FUNCTION object_1(thedum)
      use precision

!  provide function  J*(grad_a)^2 

      use esi_interface

      implicit  none
      real(wp) object_1

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),grad_a2
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum
     
      aaa(1)=adum
      the(1)=thedum
!      call esiget2dfunctions(aaa,the,1)  ! replaced by:
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      grad_a2=(rqj(1)*rqj(1)+zqj(1)*zqj(1))/(zaj(1)*rqj(1)-zqj(1)*raj(1))**2
      object_1=reducedjacobian(thedum)*grad_a2

end function object_1


FUNCTION object_2(thedum)
      use precision

!  provide function   J*exp(u^2/v_th^2)/B^2 *(grad_a)^2

      use esi_interface

      implicit  none
      real(wp) object_2

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),grad_a2
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum
    
      udum=u_parallel(adum,thedum)
      vth2dum=2.0_wp*temperature(adum,thedum)  ! v_th^2

      aaa(1)=adum
      the(1)=thedum
!      call esiget2dfunctions(aaa,the,1)   ! replaced by:
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      grad_a2=(rqj(1)*rqj(1)+zqj(1)*zqj(1))/(zaj(1)*rqj(1)-zqj(1)*raj(1))**2

      object_2=reducedjacobian(thedum)*exp(udum*udum/vth2dum)/(Bj(1)*Bj(1))*grad_a2

end function object_2


FUNCTION object_3(thedum)
      use precision

!  provide function  J*u_//

      use esi_interface

      implicit  none
      real(wp) object_3

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),grad_a2
      common/stupidpass/adum

      object_3=reducedjacobian(thedum)*u_parallel(adum,thedum)

end function object_3


FUNCTION object_4(thedum)
      use precision

!  provide function  J*(B_0^2/B^2)

      use esi_interface

      implicit  none
      real(wp) object_4

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),grad_a2
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      object_4=reducedjacobian(thedum)*(1.0_wp/Bj(1))**2

end function object_4


FUNCTION object_5(thedum)
      use precision

!  provide function  J*(B^2/B_0^2)

      use esi_interface

      implicit  none
      real(wp) object_5

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),grad_a2
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      object_5=reducedjacobian(thedum)*Bj(1)**2

end function object_5


FUNCTION object_6(thedum)
      use precision

!  provide function  J*(B_t/B^2)*|grad_a|

      use esi_interface

      implicit  none
      real(wp) object_6

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),grad_a
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum
    
      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      grad_a=sqrt(rqj(1)**2+zqj(1)**2)/abs(zaj(1)*rqj(1)-zqj(1)*raj(1))
      object_6=reducedjacobian(thedum)*(Fj(1)/rj(1))/Bj(1)**2*grad_a

end function object_6


FUNCTION object_7(thedum)
      use precision

!  provide function  J*(B_p/B)*u_//

      use esi_interface

      implicit  none
      real(wp) object_7

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),grad_a
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum
    
      udum=u_parallel(adum,thedum)

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      grad_a=sqrt(rqj(1)**2+zqj(1)**2)/abs(zaj(1)*rqj(1)-zqj(1)*raj(1))
      object_7=reducedjacobian(thedum)*udum*(grad_a/rj(1)*(-gYaj(1)))/Bj(1)

end function object_7


FUNCTION object_8(thedum)
      use precision

!  provide function  J*(B_t/B^2)

      use esi_interface

      implicit  none
      real(wp) object_8

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),grad_a2
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      object_8=reducedjacobian(thedum)*(Fj(1)/rj(1))/Bj(1)**2

end function object_8


FUNCTION object_9(thedum)
      use precision

!  provide function  J*|grad_a|

      use esi_interface

      implicit  none
      real(wp) object_9

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),grad_a
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      grad_a=sqrt(rqj(1)**2+zqj(1)**2)/abs(zaj(1)*rqj(1)-zqj(1)*raj(1))
      object_9=reducedjacobian(thedum)*grad_a

end function object_9


FUNCTION object_9s(thedum)
      use precision

!  provide function  J/|grad_a|

      use esi_interface

      implicit  none
      real(wp) object_9s

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),grad_a
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      grad_a=sqrt(rqj(1)**2+zqj(1)**2)/abs(zaj(1)*rqj(1)-zqj(1)*raj(1))
      object_9s=reducedjacobian(thedum)/grad_a

end function object_9s


FUNCTION object_10(thedum)
      use precision

!  provide function  J*|grad_theta|

      use esi_interface

      implicit  none
      real(wp) object_10

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),grad_th
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      grad_th=sqrt(zaj(1)**2+raj(1)**2)/abs(zaj(1)*rqj(1)-zqj(1)*raj(1))
      object_10=reducedjacobian(thedum)*grad_th

end function object_10


real*8 FUNCTION object_11(thedum)

!  provide function  J*R^2*|grad_a|

      use esi_interface

      implicit  none

      real*8, external:: reducedjacobian
      real*8 adum,thedum,aaa(2),the(2),grad_a
      real*8, dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      grad_a=sqrt(rqj(1)**2+zqj(1)**2)/abs(zaj(1)*rqj(1)-zqj(1)*raj(1))
      object_11=reducedjacobian(thedum)*rj(1)**2*grad_a

   end function object_11


FUNCTION object_u0(thedum)
      use precision

!  provide function  J*exp(u^2/2T) for calculating <A>=<exp(u^2/2T)>

      use esi_interface

      implicit  none
      real(wp) object_u0, tmp, tmp2

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),Ti,ui
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum    ! here adum is reveived, and then sent to reducedjacobian

      Ti=temperature(adum,thedum)
      ui=u_parallel(adum,thedum)
!cj      object_u0=exp(ui*ui/Ti*0.5_wp)
!cj      object_u0=object_u0*reducedjacobian(thedum)
      tmp=exp(ui*ui/Ti*0.5_wp)
      tmp2=tmp*reducedjacobian(thedum)
      object_u0=tmp2

end function object_u0


FUNCTION object_u1(thedum)
      use precision

!  provide function  J*(d/da)[exp(u^2/2T)] for calculating <dA/da>

      use esi_interface

      implicit  none
      real(wp) object_u1, tmp, tmp2, tmp3

      real(wp), external:: u_parallel, temperature, reducedjacobian,dT_da,domega_da
      real(wp)  adum,thedum,udum,vth2dum,aaa(2),the(2),Ti,ui,pupa,dTda,omega_t,doda
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum    ! here adum is reveived, and then sent to reducedjacobian

      aaa(1)=adum
      the(1)=thedum

      Ti=temperature(adum,thedum)
      ui=u_parallel(adum,thedum)
      dTda=dT_da(adum,thedum)

      call esiget2dfunctions(aaa,the,1)

      doda=domega_da(adum,thedum,omega_t)  ! =d(omega_t)/da
!      pupa=ui*(doda/omega_t+Fa(1)/F(1)-Ba(1)/BB(1))  ! =d(u_//)/da NOTE: cause problem when omega_t=0
      pupa=F(1)/BB(1)*doda+omega_t/BB(1)*Fa(1)-ui*Ba(1)/BB(1)   ! =d(u_//)/da

!cj      object_u1=exp(ui*ui/Ti*0.5_wp)
!cj      object_u1=object_u1*(ui/Ti*pupa-ui*ui/Ti*0.5_wp*dTda/Ti)
!cj      object_u1=object_u1*reducedjacobian(thedum)
      tmp=exp(ui*ui/Ti*0.5_wp)
      tmp2=tmp*(ui/Ti*pupa-ui*ui/Ti*0.5_wp*dTda/Ti)
      tmp3=tmp2*reducedjacobian(thedum)

      object_u1=tmp3

end function object_u1


FUNCTION object_u2(thedum)
      use precision

!  provide function  J*[dR/da*exp(u^2/2T)/R] for calculating <A*(dR/da)/R>

      use esi_interface

      implicit  none
      real(wp) object_u2, tmp, tmp2, tmp3

      real(wp), external:: u_parallel, temperature, reducedjacobian
      real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2),Ti,ui
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum    ! here adum is reveived, and then sent to reducedjacobian

      Ti=temperature(adum,thedum)
      ui=u_parallel(adum,thedum)

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
!cj      object_u2=exp(ui*ui/Ti*0.5_wp)
!cj      object_u2=object_u2*raj(1)/rj(1)
!cj      object_u2=object_u2*reducedjacobian(thedum)
      tmp=exp(ui*ui/Ti*0.5_wp)
      tmp2=tmp*raj(1)/rj(1)
      tmp3=tmp2*reducedjacobian(thedum)
      object_u2=tmp3

end function object_u2


FUNCTION object_u3(thedum)
      use precision

!  provide function  J*[(dR/da)/R] for calculating <(dR/da)/R>

      use esi_interface

      implicit  none
      real(wp) object_u3, tmp, tmp2

      real(wp), external:: reducedjacobian
      real(wp) adum,thedum,aaa(2),the(2)
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum    ! here adum is reveived, and then sent to reducedjacobian

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
!cj      object_u3=raj(1)/rj(1)
!cj      object_u3=object_u3*reducedjacobian(thedum)
      tmp=raj(1)/rj(1)
      tmp2=tmp*reducedjacobian(thedum)

      object_u3=tmp2

end function object_u3


FUNCTION object_12(thedum)
      use precision

!  provide function  J*gs*(grad_a)^2  for flux-surface_average

      use esi_interface

      implicit  none
      real(wp) object_12

      real(wp), external:: gs, reducedjacobian
      real(wp) adum,thedum,aaa(2),the(2),grad_a2
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
      common/stupidpass/adum

      aaa(1)=adum
      the(1)=thedum

!      call esiget2dfunctions(aaa,the,1)  ! replaced by:
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)

      grad_a2=(rqj(1)*rqj(1)+zqj(1)*zqj(1))/(zaj(1)*rqj(1)-zqj(1)*raj(1))**2

      object_12=reducedjacobian(thedum)*grad_a2*gs(adum,thedum)

end function object_12


subroutine get_table_4_dlnAda
!
!   calculate dln<A>/da and <A> and create tables for them; A=exp(u^2/2T)
!
!   d<A>/da=<dA/da>+2[<A*(dR/da)/R>-<A><(dR/da)/R>]
!   IMPORTANT NOTE: This subroutine assumes axisymmetry coordinates and equilibrium
!                   in which Jacobian \prop R^2  is used to obtain above expression!!!
  use precision
  use global_parameters
  use particle_array
  use esi_interface

  implicit none

  real(wp), external:: object_u0,object_u1,object_u2,object_u3
  integer i
  real(wp)  ajr,ave0,ave1,ave2,ave3

    allocate(dlnA_da(0:jr),flux_A(0:jr))
    djr=(a1-a0)/real(jr,wp)

    do i=0,jr
      ajr=min(a1,max(a0,a0+djr*real(i,wp)))
      call general_flux_average(ajr,object_u0,ave0) ! ave0=<A>
      call general_flux_average(ajr,object_u1,ave1) ! ave1=<dA/da>
      call general_flux_average(ajr,object_u2,ave2) ! ave2=<A*(dR/da)/R>
      call general_flux_average(ajr,object_u3,ave3) ! ave3=<(dR/da)/R>
      dlnA_da(i)=ave1+2.0_wp*(ave2-ave0*ave3)
      dlnA_da(i)=dlnA_da(i)/ave0
      flux_A(i)=ave0

!      if(mype==0) print*, i, 'get_table_4_dlnAda', ave0, ave1, ave2, ave3, dlnA_da(i)
    end do

end subroutine get_table_4_dlnAda


subroutine general_flux_average(adum,object_fun,averaged_result)
!
!    calculate average_result=<object_fun/J>
!
      use precision
      use global_parameters

      implicit  none

      real(wp) adum,acc,eps,D01AHF,adum_share,averaged_result
      real(wp), external:: reducedjacobian
      integer minpts,IFAIL
      real(wp) j_average,jA_average
      common/stupidpass/adum_share  ! here adum_share is sent out,
                                    ! and received by reducedjacobian and object_fun
      interface
        FUNCTION object_fun(thedum)
        use precision
        use esi_interface
        implicit  none
        real(wp) object_fun
        real(wp), external:: u_parallel, temperature, reducedjacobian
        real(wp) adum,thedum,udum,vth2dum,aaa(2),the(2)
        common/stupidpass/adum
        end function object_fun
      end interface

      adum_share=adum
      eps=1.0e-5_wp
      IFAIL=1
      j_average=D01AHF(0.0_wp,2.0_wp*pi,eps,minpts,acc,reducedjacobian,-1,IFAIL)
      eps=1.0e-5_wp
      IFAIL=1
      jA_average=D01AHF(0.0_wp,2.0_wp*pi,eps,minpts,acc,object_fun,-1,IFAIL)
      averaged_result=jA_average/j_average

end subroutine general_flux_average


Subroutine get_ave_grada_table
  use precision
  use global_parameters
  use particle_array
  use esi_interface

  implicit none

  real(wp), external:: object_9
  integer i
  real(wp)  ajr

    allocate(flux_it(0:jr))
    djr=(a1-a0)/real(jr,wp)
    do i=0,jr
      ajr=min(a1,max(a0,a0+djr*real(i,wp)))
      call general_flux_average(ajr,object_9,flux_it(i))
    end do

end subroutine get_ave_grada_table


FUNCTION arclength(thedum)
      use precision

!  provide function  1/|grad_theta|

      use esi_interface

      implicit  none
      real(wp) arclength

      real(wp) adum,thedum
      real(wp), dimension(2) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the
      common/stupidpass/adum

      aaa(1)=adum
      the(1)=thedum
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)

      arclength=abs(zaj(1)*rqj(1)-zqj(1)*raj(1))/sqrt(zaj(1)**2+raj(1)**2)   ! =1/|grad_theta|

      RETURN
end function arclength


FUNCTION density(adum,thedum)
      use precision
!
!      density=<ni>=n_0*(An+Bn*exp(-dn*(a)**xn))
!
      use global_parameters
      use spl_profiles
      use profile_array

      implicit  none
      real(wp) density, tmp, tmp2
      real(wp) adum,thedum

      integer i, ip
      real(wp)  r, fip, dela

      dela=(a1-a0)/real(jr,wp)
      dela=1.0_wp/dela

      if(plasma==1) then    ! use experimental n(r) profile
         call ezspline_interp(spn,adum,density,ier)

      else               ! use model n(r) profile
         if(q2.gt.0.0) then
            density=n_0*(An+Bn*dexp(-dn*adum**xn))
         else
            r=adum-a0
            r=r*dela
            ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
            fip=r-real(ip,wp)
            tmp=fip*ne_node(ip+1)+(1.0_wp-fip)*ne_node(ip)     
            tmp2=n_0*exp(kappan*tmp)      ! ne
            density=tmp2
!cj            density=fip*ne_node(ip+1)+(1.0_wp-fip)*ne_node(ip)     
!cj            density=n_0*exp(kappan*density)      ! ne
         end if

      end if

      return
end function density


FUNCTION dn_da(adum,thedum)
      use precision
!
!      dn_da=d<ni>/da
!
      use global_parameters
      use spl_profiles
      use profile_array

      implicit  none
      real(wp) dn_da,densi, tmp, tmp2
      real(wp) adum,thedum

      integer i, ip
      real(wp)  r, fip, dela, x

      dela=(a1-a0)/real(jr,wp)
      dela=1.0_wp/dela

      if(plasma==1) then    ! use experimental n(r) profile
         call ezspline_derivative(spn,1,adum,dn_da,ier)

      else               ! use model n(r) profile
         if(q2.gt.0.0) then
            dn_da=n_0*Bn*dexp(-dn*adum**xn)*(-dn*xn*adum**(xn-1.0_wp))
         else
            r=adum-a0
            r=r*dela
            ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
            fip=r-real(ip,wp)
            densi=fip*ne_node(ip+1)+(1.0_wp-fip)*ne_node(ip)     
            densi=n_0*exp(kappan*densi)     ! ne
            x=temconst*(adum-rc)
!cj            dn_da=kappan*exp(-x**edensity0) ! dlnn/da
!cj            dn_da=densi*dn_da               ! dn/da
            tmp=kappan*exp(-ABS(x)**edensity0) ! dlnn/da
            tmp2=densi*tmp               ! dn/da

            dn_da=tmp2
         end if
      end if

      return
end function dn_da


FUNCTION temperature(adum,thedum)
      use precision
!
!      Ti=Tc*[T_0*(At+Bt*exp(-dTemi*(a)**xt))/Tc]**r0  change made in Oct. 09
!
      use global_parameters
      use spl_profiles
      use profile_array

      implicit  none
      real(wp) temperature, tmp, tmp2
      real(wp) adum,thedum

      integer i, ip
      real(wp)  r, fip, dela

      dela=(a1-a0)/real(jr,wp)
      dela=1.0_wp/dela

      if(plasma==1) then    ! use experimental Ti(r) profile
           call ezspline_interp(spTi,adum,temperature,ier)

         else               ! use model Ti(r) profile
           if(q2.gt.0.0) then
              temperature=Tc*(T_0*(At+Bt*exp(-dTemi*adum**xt))/Tc)**r0 ! model Ti profile
           else
              r=adum-a0
              r=r*dela
              ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
              fip=r-real(ip,wp)
!cj              temperature=fip*Ti_node(ip+1)+(1.0_wp-fip)*Ti_node(ip)     ! Ti at particles
!cj              temperature=T_0*exp(q1*temperature)
              tmp=fip*Ti_node(ip+1)+(1.0_wp-fip)*Ti_node(ip)     ! Ti at particles
              tmp2=T_0*exp(q1*tmp)

              temperature=tmp2
           end if
      end if

      return
end function temperature


FUNCTION temperature_e(adum,thedum)
      use precision
!
!      Te=Tc_e*[T_e0*(At_e+Bt_e*exp(-dTeme*(a)**xt_e))/Tc_e]**b0     change made in Oct. 09
!
      use global_parameters
      use spl_profiles
      use profile_array

      implicit  none
      real(wp) temperature_e, tmp, tmp2
      real(wp) adum,thedum

      integer i, ip
      real(wp)  r, fip, dela

      dela=(a1-a0)/real(jr,wp)
      dela=1.0_wp/dela

      if(plasma==1) then    ! use experimental Te(r) profile
           call ezspline_interp(spTe,adum,temperature_e,ier)

      else               ! use model Te(r) profile
        if(q2.gt.0.0) then 
           tmp=T_e0*(At_e+Bt_e*exp(-dTeme*adum**xt_e)) ! model Te profil, the same as Ti
           tmp2=Tc_e*(tmp/Tc_e)**b0

           temperature_e=tmp2
!cj           temperature_e=T_e0*(At_e+Bt_e*exp(-dTeme*adum**xt_e)) ! model Te profil, the same as Ti
!cj           temperature_e=Tc_e*(temperature_e/Tc_e)**b0
        else
           r=adum-a0
           r=r*dela
           ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
           fip=r-real(ip,wp)
           tmp=fip*Te_node(ip+1)+(1.0_wp-fip)*Te_node(ip)     
           tmp2=T_e0*exp(flow1*tmp)    ! =Te

           temperature_e=tmp2
!cj           temperature_e=fip*Te_node(ip+1)+(1.0_wp-fip)*Te_node(ip)     
!cj           temperature_e=T_e0*exp(flow1*temperature_e)    ! =Te
        end if
      end if

      return
end function temperature_e


subroutine Get_MP_alpha_amplitude(in_a, in_theta, in_zeta, in_time, out_ga, out_grad_ga, out_partial_dga0dt) bind (C, name="get_mp_alpha_amplitude_")

      use iso_c_binding, only : c_double
      use global_parameters
      implicit  none
      integer k
      real(wp) MP_phase, ga_amplitude
      real(wp) in_a, in_theta, in_zeta, in_time
      real(wp) out_ga, out_grad_ga(3), out_partial_dga0dt

      if(electromagnetic==0) then ! when electromagnetic is turned off.
        out_ga = 0
        out_grad_ga(1:3) = 0
        out_partial_dga0dt = 0
      else
        out_ga = 0
        out_grad_ga = 0
        out_partial_dga0dt = 0
        do k=1,num_of_EMP
          MP_phase = EMP(k)%gan*in_zeta- EMP(k)%gam*in_theta - EMP(k)%gagw*in_time + EMP(k)%gagf0           
          ga_amplitude = (in_a/EMP(k)%ga_rc)*(EMP(k)%ga0 + EMP(k)%dga0dt*in_time) * exp(-abs((in_a-EMP(k)%ga_rc)/EMP(k)%ga_rw)**EMP(k)%ga_pp)
                    
          out_ga = out_ga + ga_amplitude* cos(MP_phase);

          out_grad_ga(1) = out_grad_ga(1) + ga_amplitude* cos(MP_phase)*(1.0/in_a  + (-sign(EMP(k)%ga_pp/EMP(k)%ga_rw, in_a-EMP(k)%ga_rc)*abs((in_a-EMP(k)%ga_rc)/EMP(k)%ga_rw)**(EMP(k)%ga_pp-1)))
          out_grad_ga(2) = out_grad_ga(2) + EMP(k)%gam*ga_amplitude*sin(MP_phase)
          out_grad_ga(3) = out_grad_ga(3) - EMP(k)%gan*ga_amplitude*sin(MP_phase)
          out_partial_dga0dt = out_partial_dga0dt + EMP(k)%dga0dt*(in_a/EMP(k)%ga_rc)*exp(-abs((in_a-EMP(k)%ga_rc)/EMP(k)%ga_rw)**EMP(k)%ga_pp)*cos(MP_phase) + EMP(k)%gagw*ga_amplitude*sin(MP_phase)
        enddo
      endif 
      
end subroutine Get_MP_alpha_amplitude


FUNCTION dT_da(adum,thedum)
      use precision
!
!      dT_da=dTi/da
!
      use global_parameters
      use spl_profiles
      use profile_array

      implicit  none
      real(wp) dT_da,tempi, tmp, tmp2
      real(wp) adum,thedum

      integer i, ip
      real(wp)  r, fip, dela, x

      dela=(a1-a0)/real(jr,wp)
      dela=1.0_wp/dela


      if(plasma==1) then    ! use experimental Ti(r) profile
         call ezspline_derivative(spTi,1,adum,dT_da,ier)

      else               ! use model Ti(r) profile
         if(q2.gt.0.0) then
            tempi=T_0*(At+Bt*exp(-dTemi*adum**xt))     ! original Ti profile
!cj            dT_da=T_0*Bt*exp(-dTemi*adum**xt)*(-dTemi*xt*adum**(xt-1.0_wp)) ! for original Ti profile
!cj            dT_da=r0*(tempi/Tc)**(r0-1.0d0)*dT_da      ! for new model Ti profile
            tmp=T_0*Bt*exp(-dTemi*adum**xt)*(-dTemi*xt*adum**(xt-1.0_wp)) ! for original Ti profile
            tmp2=r0*(tempi/Tc)**(r0-1.0d0)*tmp      ! for new model Ti profile

            dT_da=tmp2
         else
            r=adum-a0
            r=r*dela
            ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
            fip=r-real(ip,wp)
            tempi=fip*Ti_node(ip+1)+(1.0_wp-fip)*Ti_node(ip)     ! Ti at particles
            tempi=T_0*exp(q1*tempi)
            x=rw*(adum-rc)
!cj            dT_da=q1*exp(-x**q0)
!cj            dT_da=tempi*dT_da
            tmp=q1*exp(-ABS(x)**q0)
            tmp2=tempi*tmp

            dT_da=tmp2
         end if
      end if

      return
end function dT_da


FUNCTION dT_da_e(adum,thedum)
      use precision
!
!      dT_da_e=dTe/da
!
      use global_parameters
      use spl_profiles
      use profile_array

      implicit  none
      real(wp) dT_da_e, tempe, tmp, tmp2
      real(wp) adum,thedum

      integer i, ip
      real(wp)  r, fip, dela, x

      dela=(a1-a0)/real(jr,wp)
      dela=1.0_wp/dela

      if(plasma==1) then    ! use experimental Te(r) profile
         call ezspline_derivative(spTe,1,adum,dT_da_e,ier)

      else               ! use model Te(r) profile
         if(q2.gt.0.0) then
            tempe=T_e0*(At_e+Bt_e*exp(-dTeme*adum**xt_e))   ! original Te profile
!cj            dT_da_e=T_e0*Bt_e*exp(-dTeme*adum**xt_e)*(-dTeme*xt_e*adum**(xt_e-1.0_wp))  ! for original Te profile
!cj            dT_da_e=b0*(tempe/Tc_e)**(b0-1.0d0)*dT_da_e     ! for new model Te profile
            tmp=T_e0*Bt_e*exp(-dTeme*adum**xt_e)*(-dTeme*xt_e*adum**(xt_e-1.0_wp))  ! for original Te profile
            tmp2=b0*(tempe/Tc_e)**(b0-1.0d0)*tmp     ! for new model Te profile

            dT_da_e=tmp2
         else
            r=adum-a0
            r=r*dela
            ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
            fip=r-real(ip,wp)
            tempe=fip*Te_node(ip+1)+(1.0_wp-fip)*Te_node(ip)     ! Ti at particles
            tempe=T_e0*exp(flow1*tempe)
            x=flow2*(adum-rc)
!cj            dT_da_e=flow1*exp(-x**flow0)
!cj            dT_da_e=tempe*dT_da_e
            tmp=flow1*exp(-ABS(x)**flow0)
            tmp2=tempe*tmp

            dT_da_e=tmp2
         end if
      end if

      return
end function dT_da_e


FUNCTION domega_da(adum,thedum,omega_t)
      use precision
!
!      domega_da=d(omega_t)/da
!
      use global_parameters
      use spl_profiles
      use profile_array

      implicit  none
      real(wp) domega_da
      real(wp) adum,tempe,thedum,omega_t

      integer i, ip
      real(wp)  r, fip, dela, x

      dela=(a1-a0)/real(jr,wp)
      dela=1.0_wp/dela

      if(plasma==1) then    ! use experimental omega_t(r) profile
         call ezspline_derivative(spO,1,adum,domega_da,ier)
         call ezspline_interp(spO,adum,omega_t,ier)

      else               ! use model omega_t(r) profile
         if(q2.gt.0.0) then  ! model profile I
            domega_da=u_0*Bu*exp(-du*adum**xu)*(-du*xu*adum**(xu-1.0_wp)) 
            omega_t=u_0*(Au+Bu*exp(-du*adum**xu))
         else                ! model profile II
            r=adum-a0
            r=r*dela
            ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
            fip=r-real(ip,wp)
            tempe=fip*Iomg_node(ip+1)+(1.0_wp-fip)*Iomg_node(ip)     
            omega_t=u_0+Au*tempe         ! omega_t
            x=Bu*(adum-rc)
            domega_da=Au*exp(-ABS(x)**xu)
         end if
      end if

      return
end function domega_da


FUNCTION u_parallel(adum,thedum)
      use precision
!
!      model profile: omega_t=u_0*(Au+Bu*exp(-du*(a)**xu));
!      u_//=(I*omega_t)/B with 1/B dependence
!
      use global_parameters
      use esi_interface
      use spl_profiles
      use profile_array

      implicit  none
      real(wp) u_parallel, tmp, tmp2
      real(wp) adum,thedum,aaa(10),the(10),rj(10),zz(10),Bj(10)

      real(wp), dimension(10) :: Fj,Faj,gFaj,gYaj,raj,rqj,zaj,zqj

      integer i, ip
      real(wp)  r, fip, dela, x, tempe

      dela=(a1-a0)/real(jr,wp)
      dela=1.0_wp/dela


      aaa(1)=adum
      the(1)=thedum
      if(plasma==1) then    ! use experimental omega(r) profile
!cj         call ezspline_interp(spO,adum,u_parallel,ier)
         call ezspline_interp(spO,adum,tmp,ier)
      else                  ! use model omega_t profile
         if(q2.gt.0.0) then  ! model profile I
!cj            u_parallel=u_0*(Au+Bu*exp(-du*adum**xu))  ! this is omega_t (12/07)
            tmp=u_0*(Au+Bu*exp(-du*adum**xu))  ! this is omega_t (12/07)
         else                ! model profile II
            r=adum-a0
            r=r*dela
            ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
            fip=r-real(ip,wp)
            tempe=fip*Iomg_node(ip+1)+(1.0_wp-fip)*Iomg_node(ip)
!cj            u_parallel=u_0+Au*tempe                   ! this is omega_t
            tmp=u_0+Au*tempe                   ! this is omega_t
         end if 
      end if

      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zz,zaj,zqj,Bj,aaa,the,1)
      tmp2=Fj(1)*tmp/Bj(1)  ! =u_//=(I/B)*omega

!           call esigetrzb(rj,zz,Bj,aaa,the,1)
!           u_parallel=u_parallel*(B_0/Bj(1))

      u_parallel=tmp2

      return
end function u_parallel


subroutine n_T_u_group
!
!      Note: output includes backu=omega_t and duda=d(omega_t)/da;
!            then u_//=I*backu/B
!      To calculate u_// ~ 1/B and so on, esiget2dfunctions or others should be called to get B, I 
!
!   model profiles
!      ni=n_0*(An+Bn*exp(-dn*(a)**xn))  !flux surface averaged ni profile
!      Ti=T_0*(At+Bt*exp(-dTem*(a)**xt))
!      omega_t=u_0*(Au+Bu*exp(-du*(a)**xu))
!      borrow array gh(:) which is defined in esilink2c
!      gh(:)=dln<A>/da at particles with A=exp(u^2/2T)
!
      use precision
      use global_parameters
      use particle_array
      use esi_interface
      use spl_profiles
      use profile_array

      implicit  none
      integer i, ip
      real(wp)  r, fip, dela, tempi

      dela=(a1-a0)/real(jr,wp)
      dela=1.0_wp/dela

      if(plasma==1) then    ! use experimental ni(r), T(r) and omega(r) profiles
         if(linterp_from_nodes==1) then 

!$omp parallel do private(i,r,ip,fip)
            do i=1,mi
               r=zion(1,i)-a0
               r=r*dela
               ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
               fip=r-real(ip,wp)

               backn(i)=fip*ne_node(ip+1)+(1.0_wp-fip)*ne_node(ip)     ! density at particles
               backT(i)=fip*Ti_node(ip+1)+(1.0_wp-fip)*Ti_node(ip)     ! Ti at particles
               backu(i)=fip*Iomg_node(ip+1)+(1.0_wp-fip)*Iomg_node(ip) ! omega_t at particles

               dlnnda(i)=fip*dne_node(ip+1)+(1.0_wp-fip)*dne_node(ip)   ! dn/da  at particles
               dlnTda(i)=fip*dTi_node(ip+1)+(1.0_wp-fip)*dTi_node(ip)   ! dTi/da at particles
               duda(i)=fip*dIomg_node(ip+1)+(1.0_wp-fip)*dIomg_node(ip) ! d(omega_t)/da at particles
               gh(i)=fip*dlnA_da(ip+1)+(1.0_wp-fip)*dlnA_da(ip)         ! dln<A>/da at particles with
                                                                     ! A=exp(u^2/2T)
            enddo

         else              ! use spline 
            call ezspline_interp(spn,mi,zion(1,:),backn,ier)
            call ezspline_interp(spTi,mi,zion(1,:),backT,ier)
            call ezspline_interp(spO,mi,zion(1,:),backu,ier)
            call ezspline_derivative(spn,1,mi,zion(1,:),dlnnda,ier)
            call ezspline_derivative(spTi,1,mi,zion(1,:),dlnTda,ier)
            call ezspline_derivative(spO,1,mi,zion(1,:),duda,ier)
!$omp parallel do private(i,r,ip,fip)
            do i=1,mi
               r=zion(1,i)-a0
               r=r*dela
               ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
               fip=r-real(ip,wp)
               gh(i)=fip*dlnA_da(ip+1)+(1.0_wp-fip)*dlnA_da(ip)         ! dln<A>/da at particles with
                                                                     ! A=exp(u^2/2T)
            enddo

         end if

!$omp parallel do private(i)
         do i=1,mi
            dlnnda(i)=dlnnda(i)/backn(i)      !=dlnn/da
            dlnTda(i)=dlnTda(i)/backT(i)      !=dlnT/da
         end do

      else                  ! use model profiles

         if(q2.gt.0.0) then  ! first analytic profile
!$omp parallel do private(i,r,ip,fip)
            do i=1,mi
               backn(i)=n_0*Bn*exp(-dn*zion(1,i)**xn)
               backT(i)=T_0*Bt*exp(-dTemi*zion(1,i)**xt)
               backu(i)=u_0*Bu*exp(-du*zion(1,i)**xu)
   
               dlnnda(i)=backn(i)*(-dn*xn*zion(1,i)**(xn-1.0_wp))       !=dn/da
               dlnTda(i)=backT(i)*(-dTemi*xt*zion(1,i)**(xt-1.0_wp))    !=dT/da
               duda(i)  =backu(i)*(-du*xu*zion(1,i)**(xu-1.0_wp))       !=d[omega_t(a)]/da 

               backn(i)=n_0*An+backn(i)
               backT(i)=T_0*At+backT(i)
               backu(i)=u_0*Au+backu(i)        ! then u_//=I*backu/B
                                            ! backu=u_0*(Au+Bu*exp(-du*aa**xu))

               dlnnda(i)=dlnnda(i)/backn(i)      !=dlnn/da
               dlnTda(i)=dlnTda(i)/backT(i)      !=dlnT/da
   
               backT(i)=Tc*(backT(i)/Tc)**r0     ! for new model Ti profile
               dlnTda(i)=dlnTda(i)*r0            ! for new model Ti profile

               r=zion(1,i)-a0
               r=r*dela
               ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
               fip=r-real(ip,wp)
               gh(i)=fip*dlnA_da(ip+1)+(1.0_wp-fip)*dlnA_da(ip)         ! dln<A>/da at particles with
                                                                  ! A=exp(u^2/2T)
            end do
         else      ! second analytic profile
!$omp parallel do private(i,r,ip,fip,tempi)
            do i=1,mi
               r=zion(1,i)-a0
               r=r*dela
               ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
               fip=r-real(ip,wp)
               gh(i)=fip*dlnA_da(ip+1)+(1.0_wp-fip)*dlnA_da(ip)         ! dln<A>/da at particles with

               tempi=fip*Ti_node(ip+1)+(1.0_wp-fip)*Ti_node(ip)     
               backT(i)=T_0*exp(q1*tempi)                              ! Ti at particles

               tempi=rw*(zion(1,i)-rc)
               dlnTda(i)=q1*exp(-ABS(tempi)**q0)                            ! =dlnT/da

               tempi=fip*ne_node(ip+1)+(1.0_wp-fip)*ne_node(ip)
               backn(i)=n_0*exp(kappan*tempi)                              ! ne at particles

               tempi=temconst*(zion(1,i)-rc)
               dlnnda(i)=kappan*exp(-ABS(tempi)**edensity0)                     ! =dlnn/da

               tempi=fip*Iomg_node(ip+1)+(1.0_wp-fip)*Iomg_node(ip)     
               backu(i)=u_0+Au*tempi                                   ! this is omega_t; then u_//=I*backu/B

               tempi=Bu*(zion(1,i)-rc)
               duda(i)  =Au*exp(-ABS(tempi)**xu)       ! =d[omega_t(a)]/da
            end do
         end if
      end if
 
!     neglect rotation
!       backu=0.0    
!       duda=0.0

end subroutine n_T_u_group


subroutine n_T_u_group_e
!
! get n,Te,ue for electrons
!
!      esiget2dfunctions or others should be called first to get B for this subroutine, 
!      which is needed to calculate u_// ~ 1/B
!   model profiles
!      ne=n_0*(An+Bn*exp(-dn*(a)**xn))  !flux surface averaged ni profile
!      Te=T_0*(At+Bt*exp(-dTem*(a)**xt))
!      u_//=u_0*(Au+Bu*exp(-du*(a)**xu))*(B_0/B) with 1/B dependence
!
      use precision
      use global_parameters
      use particle_array
      use esi_interface
      use spl_profiles
      use profile_array

      implicit  none
      integer i, ip
      real(wp)  r, fip, dela, tempe

      if(plasma==1) then    ! use experimental ni(r), T(r) and omega(r) profiles
         if(linterp_from_nodes==1) then
            dela=(a1-a0)/real(jr,wp)
            dela=1.0_wp/dela

!$omp parallel do private(i,r,ip,fip)
            do i=1,me
               r=zeon(1,i)-a0
               r=r*dela
               ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
               fip=r-real(ip,wp)

               backne(i)=fip*ne_node(ip+1)+(1.0_wp-fip)*ne_node(ip) ! density at particles
               backTe(i)=fip*Te_node(ip+1)+(1.0_wp-fip)*Te_node(ip) ! Te at particles

               dlnndae(i)=fip*dne_node(ip+1)+(1.0_wp-fip)*dne_node(ip) ! dn/da  at particles
               dlnTdae(i)=fip*dTe_node(ip+1)+(1.0_wp-fip)*dTe_node(ip) ! dTe/da at particles
            enddo

         else              ! use spline
            call ezspline_interp(spn,me,zeon(1,:),backne,ier)
            call ezspline_interp(spTe,me,zeon(1,:),backTe,ier)
!            call ezspline_interp(spO,me,zeon(1,:),backue,ier)
!            backu=backu*F/B
            call ezspline_derivative(spn,1,me,zeon(1,:),dlnndae,ier)
            call ezspline_derivative(spTe,1,me,zeon(1,:),dlnTdae,ier)
         end if

!$omp parallel do private(i)
         do i=1,me
            dlnndae(i)=dlnndae(i)/backne(i)      !=dlnn/da
            dlnTdae(i)=dlnTdae(i)/backTe(i)      !=dlnT/da
         end do

      else                  ! use model profiles

        if(q2.gt.0.0) then  ! first analytic profile
!$omp parallel do private(i)
           do i=1,me
              backne(i)=n_0*Bn*exp(-dn*zeon(1,i)**xn)
              backTe(i)=T_e0*Bt_e*exp(-dTeme*zeon(1,i)**xt_e)

              dlnndae(i)=backne(i)*(-dn*xn*zeon(1,i)**(xn-1.0_wp))       !=dn/da
              dlnTdae(i)=backTe(i)*(-dTeme*xt_e*zeon(1,i)**(xt_e-1.0_wp))    !=dT/da

              backne(i)=n_0*An+backne(i)
              backTe(i)=T_e0*At_e+backTe(i)

              dlnndae(i)=dlnndae(i)/backne(i)      !=dlnn/da
              dlnTdae(i)=dlnTdae(i)/backTe(i)      !=dlnT/da

              backTe(i)=Tc_e*(backTe(i)/Tc_e)**b0     ! for new model Ti profile
              dlnTdae(i)=dlnTdae(i)*b0                ! for new model Ti profile

!               backu=u_0*(Au+Bu*exp(-du*aa**xu))*(B_0/B)
           end do
        else      ! second analytic profile
           dela=(a1-a0)/real(jr,wp)
           dela=1.0_wp/dela
!$omp parallel do private(i,r,ip,fip,tempe)
           do i=1,me
              r=zeon(1,i)-a0
              r=r*dela
              ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
              fip=r-real(ip,wp)

              tempe=fip*Te_node(ip+1)+(1.0_wp-fip)*Te_node(ip)
              backTe(i)=T_e0*exp(flow1*tempe)                              ! Te at particles

              tempe=flow2*(zeon(1,i)-rc)
              dlnTdae(i)=flow1*exp(-ABS(tempe)**flow0)                            ! =dlnT/da

              tempe=fip*ne_node(ip+1)+(1.0_wp-fip)*ne_node(ip)
              backne(i)=n_0*exp(kappan*tempe)                              ! ne at particles

              tempe=temconst*(zeon(1,i)-rc)
              dlnndae(i)=kappan*exp(-ABS(tempe)**edensity0)                     ! =dlnn/da

!               backu=u_0*(Au+Bu*exp(-du*aa**xu))*(B_0/B)
           end do
        end if

      end if

       backue=0._wp    ! neglect rotation

end subroutine n_T_u_group_e


subroutine Do_FieldLineFollowing_for_Poincare_Plot
  
    use precision
    use global_parameters
    
    type Poincare_Point_Information
        integer Label  ! Label for each field line        
        real(wp) a, gq, gf
        real(wp) R, Z
    end type
    
    real(wp) out_ga, out_grad_ga(3), out_partial_dg0dt
    real(wp) B_vec_cntv(3)

    real(wp) dgf ! toroidal step for field-line-following integration
    real(wp) a_RK(4), gq_RK(4)

    real(wp) now_a, now_gq, now_gf, specific_time
    real(wp) prev_a, prev_gq, prev_gf
    integer :: maximum_tor_rotation = 100
    integer tor_dir, aid, num_a, gqid, num_gq, PPI_id
    type(Poincare_Point_Information), dimension(:,:,:), allocatable:: PPI
    real(wp) intp_weight, dummy_real

    specific_time = 0.0 ! time is constant because we want the field Poincare Plot at specific time t=t0
    
    dgf = 2.0*pi/1000.0 ! toroidal step for field-line-following integration
    num_a = 50
    num_gq = 7

    allocate( PPI(2*maximum_tor_rotation, num_gq, num_a))
    !PPI(:,:,:) = 0

    do aid=1,num_a
        do gqid = 1, num_gq

            PPI_id = 0
            do tor_dir=-1,1,2
                dgf = (tor_dir)*2.0*pi/1000.0

                now_a = a0 + (aid-1)*(a1-a0)/(num_a-1.0)
                now_gq = (gqid-1)*(2.0*pi)/(num_gq) ! To avoid overlap between gq=0 and gq=2pi points
                now_gf = 0.0        

                prev_a = now_a
                prev_gq = now_gq
                prev_gf = now_gf

                do while(.true.)

                    !if( modulo(abs(now_gf), 2*pi) < dgf) then ! When B-line-integrated point just passed gf=0 poloidal plane.
                    if( floor(now_gf/(2.0*pi)) /= floor(prev_gf/(2.0*pi)) ) then
                        PPI_id = PPI_id + 1

                        PPI(PPI_id, gqid, aid)%Label = (aid-1)*num_gq + (gqid-1)            
                        !intp_weight = modulo(abs(now_gf), 2.0*pi) / dgf
                        intp_weight = abs(min(modulo(now_gf, 2.0*pi), 2.0*pi - modulo(now_gf, 2.0*pi))/dgf) ! interpolation weight for previous point

                        PPI(PPI_id, gqid, aid)%a  = intp_weight*prev_a + (1.0-intp_weight)*now_a
                        PPI(PPI_id, gqid, aid)%gq = intp_weight*prev_gq + (1.0-intp_weight)*now_gq
                        PPI(PPI_id, gqid, aid)%gf = intp_weight*prev_gf + (1.0-intp_weight)*now_gf

                        call esigetrzb(PPI(PPI_id, gqid, aid)%R, PPI(PPI_id, gqid, aid)%Z, dummy_real, PPI(PPI_id, gqid, aid)%a, modulo(PPI(PPI_id, gqid, aid)%gq, 2*pi), 1)
                    endif

                    if(abs(now_gf)>maximum_tor_rotation*(2*pi)) then
                        exit      
                    endif

                    prev_a = now_a
                    prev_gq = now_gq
                    prev_gf = now_gf

                    ! Runge-Kutta 4th order
                    call esiget_Intp_B_vec_contravariant(now_a, now_gq, now_gf, specific_time, B_vec_cntv)
                    a_RK(1)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
                    gq_RK(1)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
                    call esiget_Intp_B_vec_contravariant(now_a+0.5*a_RK(1), now_gq+0.5*gq_RK(1), now_gf+0.5*dgf, specific_time, B_vec_cntv)
                    a_RK(2)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
                    gq_RK(2)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
                    call esiget_Intp_B_vec_contravariant(now_a+0.5*a_RK(2), now_gq+0.5*gq_RK(2), now_gf+0.5*dgf, specific_time, B_vec_cntv)
                    a_RK(3)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
                    gq_RK(3)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
                    call esiget_Intp_B_vec_contravariant(now_a+a_RK(3), now_gq+gq_RK(3), now_gf+dgf, specific_time, B_vec_cntv)
                    a_RK(4)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
                    gq_RK(4)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
                    
                    now_a = now_a + (a_RK(1)+2.0*a_RK(2)+2.0*a_RK(3)+a_RK(4))/6.0
                    now_gq = now_gq + (gq_RK(1)+2.0*gq_RK(2)+2.0*gq_RK(3)+gq_RK(4))/6.0
                    now_gf = now_gf + dgf                

                end do
            end do
        end do 
    end do


    open(181125,file='FLF_Poincare.dat',status='replace',action='write')
    do aid=1,num_a
        do gqid = 1, num_gq
            do tor_dir=1,2*maximum_tor_rotation
                write(181125, "(I4,1x,5(ES14.6,1x))") PPI(tor_dir,gqid,aid)%Label, PPI(tor_dir,gqid,aid)%a, PPI(tor_dir,gqid,aid)%gq, PPI(tor_dir,gqid,aid)%gf, PPI(tor_dir,gqid,aid)%R, PPI(tor_dir,gqid,aid)%Z
            enddo
        enddo
    enddo
    flush(181125)
    close(181125)

    write(*,*) "Do_FieldLineFollowing_for_Poincare_Plot OK"


    deallocate(PPI)

end subroutine Do_FieldLineFollowing_for_Poincare_Plot


subroutine Ti_at_electrons
!
!   model profiles
!      Ti=T_0*(At+Bt*exp(-dTem*(a)**xt))
!
!   NOTE: borrow array backT to store Ti at electrons -- make sure no memory conflict!!!
!
      use precision
      use global_parameters
      use particle_array
      use spl_profiles
      use profile_array

      implicit  none
      integer i, ip
      real(wp)  r, fip, dela, tempi

      if(plasma==1) then    ! use experimental ni(r), T(r) and omega(r) profiles
         if(linterp_from_nodes==1) then
            dela=(a1-a0)/real(jr,wp)
            dela=1.0_wp/dela

!$omp parallel do private(i,r,ip,fip)
            do i=1,me
               r=zeon(1,i)-a0
               r=r*dela
               ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
               fip=r-real(ip,wp)
               backT(i)=fip*Ti_node(ip+1)+(1.0_wp-fip)*Ti_node(ip) ! Ti at particles
            enddo

         else              ! use spline
            call ezspline_interp(spTi,me,zeon(1,:),backT,ier)
         end if

      else                  ! use model profiles

         if(q2.gt.0.0) then  ! first analytic profile
!$omp parallel do private(i)
            do i=1,me
               backT(i)=T_0*Bt*exp(-dTemi*zeon(1,i)**xt)
               backT(i)=Tc*(backT(i)/Tc)**r0
            end do
         else
            dela=(a1-a0)/real(jr,wp)
            dela=1.0_wp/dela
!$omp parallel do private(i,r,ip,fip,tempi)
            do i=1,me
               r=zeon(1,i)-a0
               r=r*dela
               ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
               fip=r-real(ip,wp)
               tempi=fip*Ti_node(ip+1)+(1.0_wp-fip)*Ti_node(ip)
               backT(i)=T_0*exp(q1*tempi)                              ! Ti at particles
            end do
         end if

      end if

end subroutine Ti_at_electrons


FUNCTION SMOOTHSTEP(xmin,xmax,x)
      use precision
      use global_parameters,only: pi
!
!      SMOOTHSTEP generates a smooth step function using 
!      SMOOTHSTEP=1/2-cos(pi*x)/2 
!      SMOOTHSTEP=0 if x<xmin; =1 if x>xmax; otherwise 0<SMOOTHSTEP<1
!      with xmax > xmin
!

      implicit  none
      real(wp) SMOOTHSTEP
      real(wp) xmin,xmax,x,y

      y=(x-xmin)/(xmax-xmin)
      y=max(0.0_wp,min(1.0_wp,y))
      SMOOTHSTEP=0.5_wp-0.5_wp*cos(pi*y)

      RETURN
end function SMOOTHSTEP


FUNCTION FUN3(adum)
      use precision
!
!    provide function FUN3=exp{-[(r-rc)/rw]**q0} for calculating A=int(FUN3*dr) 
!
      use global_parameters, only: rc, rw, q0

      implicit  none
      real(wp) FUN3, x, adum

      x=ABS(rw*(adum-rc))
      FUN3=exp(-ABS(x)**q0)

      RETURN
end function FUN3


FUNCTION FUN3e(adum)
      use precision
!
!    electron version of FUN3 for sepcifing T_e
!
!    provide function FUN3e=exp{-[(r-rc)/flow2]**flow0} for calculating A=int(FUN3e*dr)
!    with rw -> flow2, q0 -> flow0
!
      use global_parameters, only: rc, flow2, flow0

      implicit  none
      real(wp) FUN3e, x, adum

      x=ABS(flow2*(adum-rc))
      FUN3e=exp(-ABS(x)**flow0)

      RETURN
end function FUN3e


FUNCTION FUN3n(adum)
      use precision
!      
!    density version of FUN3 for sepcifing ne
!
!    provide function FUN3n=exp{-[(r-rc)/temconst]**edensity0} for calculating A=int(FUN3n*dr)
!    with rw -> temconst, q0 -> edensity0
!
      use global_parameters, only: rc, edensity0, temconst

      implicit  none
      real(wp) FUN3n, x, adum

      x=ABS(temconst*(adum-rc))
      FUN3n=exp(-ABS(x)**edensity0)

      RETURN
end function FUN3n


FUNCTION FUN3u(adum)
      use precision
!
!    density version of FUN3 for sepcifing omega_t
!
!    provide function FUN3u=exp{-[(r-rc)/Bu]**xu} for calculating A=int(FUN3u*dr)
!    with rw -> Bu, q0 -> xu
!
      use global_parameters, only: rc, xu, Bu

      implicit  none
      real(wp) FUN3u, x, adum

      x=ABS(Bu*(adum-rc))
      FUN3u=exp(-ABS(x)**xu)

      RETURN
end function FUN3u



FUNCTION gs(adum,thedum)
      use precision
!
!      gs=sum_over_species[Zs/B^2*n_s0*(e0/es)*(ms/m0)]
!
      use global_parameters
      use esi_interface
      use spl_profiles
      use profile_array

      implicit  none
      real(wp) gs, tmp
      real(wp) adum,thedum,aaa(10),the(10),rj(10),zz(10),Bj(10)

      real(wp), dimension(10) :: Fj,Faj,gFaj,gYaj,raj,rqj,zaj,zqj

      integer ns, i 
      real(wp)  Zs

      real(wp), external :: density

      ns=1        ! one ion species
      Zs=1.0_wp   ! ion charge

      aaa(1)=adum
      the(1)=thedum

      tmp=0.0_wp
      do i=1,ns
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zz,zaj,zqj,Bj,aaa,the,1)
         tmp=tmp+Zs/(Bj(1)*Bj(1))*density(adum,thedum)*c0cs_i*msm0_i
      end do
      gs=tmp

      return
end function gs
