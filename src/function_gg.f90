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


subroutine Get_MP_alpha_amplitude_arr_at_zion(in_time, out_ga, out_grad_ga, out_partial_dga0dt)

   use iso_c_binding, only : c_double
   use global_parameters
   use particle_array
   implicit  none
   integer k
   real(wp) MP_phase, ga_amplitude
   real(wp),intent(in) :: in_time
   real(wp),intent(out) :: out_ga(mi), out_grad_ga(3,mi), out_partial_dga0dt(mi)
   integer m


   if(electromagnetic==0) then ! when electromagnetic is turned off.
      out_ga = 0
      out_grad_ga(1:3,:) = 0
      out_partial_dga0dt = 0
   else
      !$acc parallel loop present(zion)
      do m=1,mi
         out_ga(m) = 0
         !$acc loop seq
         do k=1,3
            out_grad_ga(k,m) = 0
         enddo
         out_partial_dga0dt(m) = 0
         !$acc loop seq
         do k=1,num_of_EMP
            MP_phase = EMP(k)%gan*zion(3,m)- EMP(k)%gam*zion(2,m) - EMP(k)%gagw*in_time + EMP(k)%gagf0
            ga_amplitude = (zion(1,m)/EMP(k)%ga_rc)*(EMP(k)%ga0 + EMP(k)%dga0dt*in_time) * exp(-abs((zion(1,m)-EMP(k)%ga_rc)/EMP(k)%ga_rw)**EMP(k)%ga_pp)

            out_ga(m) = out_ga(m) + ga_amplitude* cos(MP_phase);

            out_grad_ga(1,m) = out_grad_ga(1,m) + ga_amplitude* cos(MP_phase)*(1.0/zion(1,m)  + (-sign(EMP(k)%ga_pp/EMP(k)%ga_rw, zion(1,m)-EMP(k)%ga_rc)*abs((zion(1,m)-EMP(k)%ga_rc)/EMP(k)%ga_rw)**(EMP(k)%ga_pp-1)))
            out_grad_ga(2,m) = out_grad_ga(2,m) + EMP(k)%gam*ga_amplitude*sin(MP_phase)
            out_grad_ga(3,m) = out_grad_ga(3,m) - EMP(k)%gan*ga_amplitude*sin(MP_phase)
            out_partial_dga0dt(m) = out_partial_dga0dt(m) + EMP(k)%dga0dt*(zion(1,m)/EMP(k)%ga_rc)*exp(-abs((zion(1,m)-EMP(k)%ga_rc)/EMP(k)%ga_rw)**EMP(k)%ga_pp)*cos(MP_phase) + EMP(k)%gagw*ga_amplitude*sin(MP_phase)
         enddo
      enddo
   endif
end subroutine Get_MP_alpha_amplitude_arr_at_zion

subroutine Get_MP_alpha_amplitude_arr_at_zeon(in_time, out_ga, out_grad_ga, out_partial_dga0dt)

   use iso_c_binding, only : c_double
   use global_parameters
   use particle_array
   implicit  none
   integer k
   real(wp) MP_phase, ga_amplitude
   real(wp),intent(in) :: in_time
   real(wp),intent(out) :: out_ga(me), out_grad_ga(3,me), out_partial_dga0dt(me)
   integer m


   if(electromagnetic==0) then ! when electromagnetic is turned off.
      out_ga = 0
      out_grad_ga(1:3,:) = 0
      out_partial_dga0dt = 0
   else
      !$acc parallel loop present(zeon)
      do m=1,me
         out_ga(m) = 0
         !$acc loop seq
         do k=1,3
            out_grad_ga(k,m) = 0
         enddo
         out_partial_dga0dt(m) = 0
         !$acc loop seq
         do k=1,num_of_EMP
            MP_phase = EMP(k)%gan*zeon(3,m)- EMP(k)%gam*zeon(2,m) - EMP(k)%gagw*in_time + EMP(k)%gagf0
            ga_amplitude = (zeon(1,m)/EMP(k)%ga_rc)*(EMP(k)%ga0 + EMP(k)%dga0dt*in_time) * exp(-abs((zeon(1,m)-EMP(k)%ga_rc)/EMP(k)%ga_rw)**EMP(k)%ga_pp)

            out_ga(m) = out_ga(m) + ga_amplitude* cos(MP_phase);

            out_grad_ga(1,m) = out_grad_ga(1,m) + ga_amplitude* cos(MP_phase)*(1.0/zeon(1,m)  + (-sign(EMP(k)%ga_pp/EMP(k)%ga_rw, zeon(1,m)-EMP(k)%ga_rc)*abs((zeon(1,m)-EMP(k)%ga_rc)/EMP(k)%ga_rw)**(EMP(k)%ga_pp-1)))
            out_grad_ga(2,m) = out_grad_ga(2,m) + EMP(k)%gam*ga_amplitude*sin(MP_phase)
            out_grad_ga(3,m) = out_grad_ga(3,m) - EMP(k)%gan*ga_amplitude*sin(MP_phase)
            out_partial_dga0dt(m) = out_partial_dga0dt(m) + EMP(k)%dga0dt*(zeon(1,m)/EMP(k)%ga_rc)*exp(-abs((zeon(1,m)-EMP(k)%ga_rc)/EMP(k)%ga_rw)**EMP(k)%ga_pp)*cos(MP_phase) + EMP(k)%gagw*ga_amplitude*sin(MP_phase)
         enddo
      enddo
   endif
end subroutine Get_MP_alpha_amplitude_arr_at_zeon

subroutine Generate_FLabel_piecewise_radial_base_function

   use precision
   use global_parameters

   implicit none

   integer i,k
   real(wp) resonant_q(num_of_EMP), resonant_a(num_of_EMP), island_w(num_of_EMP)
   real(wp):: Find_a_at_specific_safety_q
   real(wp) a_ori(400), a_tilde(400),da, ff
   integer dim_a, aid
   integer flag_island_inside,prev_flag
   real(wp) tmp_pcws_a(200), tmp_pcws_y(200)
   integer :: pcws_id=0

   if(electromagnetic/=1 .or. num_of_EMP==0) then
      return
   endif

   do k=1,num_of_EMP
      resonant_q(k) = real(EMP(k)%gam)/real(EMP(k)%gan)
      resonant_a(k) = Find_a_at_specific_safety_q(resonant_q(k))
      island_w(k) = 0.15
   enddo

   dim_a = size(a_ori)
   da = (a1-a0)/real(dim_a-1)
   a_ori = (/(a0+da*real(i-1),i=1,dim_a)/)

   a_tilde(1)=a0

   pcws_id=pcws_id+1
   tmp_pcws_a(pcws_id)=a0-1e-5 ! Subtract small number at left edge for later safe interpolation
   tmp_pcws_y(pcws_id)=a0-1e-5 ! Subtract small number at left edge for later safe interpolation

   do i=2,dim_a
      flag_island_inside = 0
      do k=1,num_of_EMP
         if(abs(a_ori(i)-resonant_a(k)) < 0.5*island_w(k)) then
            flag_island_inside=1
         endif
      enddo

      if(prev_flag+flag_island_inside==1) then ! only one inside island
         if(flag_island_inside==1) then
            a_tilde(i) =a_tilde(i-1) + minval(abs(resonant_a(:)-0.5*island_w(:)-a_ori(i-1)))

            pcws_id=pcws_id+1
            tmp_pcws_a(pcws_id) = a_ori(i-1) + minval(abs(resonant_a(:)-0.5*island_w(:)-a_ori(i-1)))
            tmp_pcws_y(pcws_id) = a_tilde(i)
         else
            a_tilde(i) =a_tilde(i-1) + minval(abs(a_ori(i)- (resonant_a(:)+0.5*island_w(:))))

            pcws_id=pcws_id+1
            tmp_pcws_a(pcws_id) = a_ori(i)-minval(abs(a_ori(i)- (resonant_a(:)+0.5*island_w(:))))
            tmp_pcws_y(pcws_id) = a_tilde(i-1)
         endif
      elseif(prev_flag+flag_island_inside==2) then !both inside island
         a_tilde(i) =a_tilde(i-1)
      else ! both outside island
         a_tilde(i) =a_tilde(i-1)+da
      endif
      prev_flag = flag_island_inside
   enddo
   pcws_id=pcws_id+1
   tmp_pcws_a(pcws_id) = a_ori(dim_a)  +1e-5 ! Add small number at right edge for later safe interpolation
   tmp_pcws_y(pcws_id) = a_tilde(dim_a)+1e-5 ! Add small number at right edge for later safe interpolation


   num_of_FLabel_pcws = pcws_id
   write(stdout,*) "In pcws, num_of_FLabel_pcws=",num_of_FLabel_pcws
   call flush(stdout)
   ! Allocate piecewise elements and copy from temporary things
   allocate(FLabel_pcws_a(num_of_FLabel_pcws), FLabel_pcws_y(num_of_FLabel_pcws))
   FLabel_pcws_a=tmp_pcws_a(1:num_of_FLabel_pcws)
   FLabel_pcws_y=tmp_pcws_y(1:num_of_FLabel_pcws)

end subroutine Generate_FLabel_piecewise_radial_base_function

! function Cal_FLabel_Function(in_a,in_gq,in_gf) result(out_val)

!    use precision
!    use global_parameters

!    implicit none

!    real(wp), intent(in) :: in_a,in_gq,in_gf
!    integer i,k
!    real(wp) resonant_q(num_of_EMP), resonant_a(num_of_EMP), island_w(num_of_EMP)
!    real(wp):: Find_a_at_specific_safety_q, out_val
!    real(wp) MP_phase
!    integer target_mode_n,target_mode_m,num_of_target_safety_q
!    real(wp) target_safety_q(100)
!    real(wp) a_ori(200), a_tilde(200),da, ff
!    integer dim_a, aid
!    integer flag_island_inside,prev_flag
!    logical ex


!    do k=1,num_of_EMP
!       resonant_q(k) = real(EMP(k)%gam)/real(EMP(k)%gan)
!       resonant_a(k) = Find_a_at_specific_safety_q(resonant_q(k))
!       island_w(k) = 0.1
!    enddo

!    dim_a = size(a_ori)
!    da = (a1-a0)/real(dim_a-1)
!    a_ori = (/(a0+da*real(i-1),i=1,dim_a)/)

!    a_tilde(1)=a0
!    do i=2,dim_a
!       flag_island_inside = 0
!       do k=1,num_of_EMP
!          if(abs(a_ori(i)-resonant_a(k)) < 0.5*island_w(k)) then
!             flag_island_inside=1
!          endif
!       enddo

!       if(prev_flag+flag_island_inside==1) then ! only one inside island
!          ! a_tilde(i) =a_tilde(i-1) + min([abs(resonant_a(:)-0.5*island_w(:)-a_ori(i)), abs(a_ori(i)- (resonant_a(:)+0.5*island_w(:)))])
!          if(flag_island_inside==1) then
!             a_tilde(i) =a_tilde(i-1) + minval(abs(resonant_a(:)-0.5*island_w(:)-a_ori(i-1)))
!          else
!             a_tilde(i) =a_tilde(i-1) + minval(abs(a_ori(i)- (resonant_a(:)+0.5*island_w(:))))
!          endif
!       elseif(prev_flag+flag_island_inside==2) then !both inside island
!          a_tilde(i) =a_tilde(i-1)
!       else ! both outside island
!          a_tilde(i) =a_tilde(i-1)+da
!       endif
!       prev_flag = flag_island_inside
!    enddo

!    ! inquire(file='da.dat',exist=ex)
!    ! if(.not.ex) then
!    !    open(12321,file='da.dat',action='write')
!    !    do i=1,dim_a
!    !       write(12321, *) a_tilde(i)
!    !    enddo
!    !    call flush(12321)
!    !    close(12321)
!    ! endif

!    ! Get linear interpolation of a_tilde as a background FLabel value
!    aid=max(1,min(int((in_a-a0)/da)+1,dim_a-1))
!    ff = (in_a-(a0+(aid-1)*da))/da
!    out_val = (1.0-ff)*a_tilde(aid) + ff*a_tilde(aid+1)

!    ! Add quadratic function values inside island
!    do k=1,num_of_EMP
!       if(abs(in_a-resonant_a(k)) < 0.5*island_w(k)) then
!          MP_phase = EMP(k)%gan*in_gf- EMP(k)%gam*in_gq + EMP(k)%gagf0 !- EMP(k)%gagw*in_time
!          ! out_val= out_val + 2*(1.0-cos(MP_phase))*((in_a-resonant_a(k))**2.0- 0.5*island_w(k)**2.0)
!          out_val= out_val + (2.0-cos(MP_phase))*(abs(in_a-resonant_a(k))- 0.5*island_w(k))
!       endif
!       ! out_val= out_val + 0.2*(1.0-cos(MP_phase))*exp(-(abs(in_a-resonant_a(k))/(0.5*island_w)))**pp
!    enddo



!    ! ! Put small Peaky FLabel values at rational q surfaces (regardless of EMP)
!    ! ! This is for better representation of possible resonant structure due to nonlinear interactions of EMPs.
!    ! k=1
!    ! do target_mode_n = 1,3
!    !    do target_mode_m = 1,10
!    !       if(target_mode_m <= 3*target_mode_n) then ! search for the safety_q<=3
!    !          target_safety_q(k) = real(target_mode_m)/real(target_mode_n)
!    !          k=k+1
!    !       endif
!    !    enddo
!    ! enddo
!    ! call find_unique_sorted_array_and_its_num_of_elements(target_safety_q, 1, k-1, num_of_target_safety_q)

!    ! do k=1,num_of_target_safety_q
!    !    resonant_q = target_safety_q(k)
!    !    resonant_a = Find_a_at_specific_safety_q(resonant_q)
!    !    pp = 4.0
!    !    out_val= out_val + 0.1*exp(-(abs(in_a-resonant_a)/0.05))**pp
!    ! enddo

! end function Cal_FLabel_Function

function Cal_FLabel_Function(in_a,in_gq,in_gf) result(out_val)

   use precision
   use global_parameters

   implicit none

   real(wp), intent(in) :: in_a,in_gq,in_gf
   integer i,k,pwid
   real(wp) resonant_q(num_of_EMP), resonant_a(num_of_EMP), island_w(num_of_EMP)
   real(wp),external:: Find_a_at_specific_safety_q
   real(wp) out_val
   real(wp) MP_phase
   real(wp) ff
   real(wp) aa(2),yy(2)

   if(.not.allocated(FLabel_pcws_a)) then
      call Generate_FLabel_piecewise_radial_base_function
   endif

   if(electromagnetic/=1 .or. num_of_EMP==0) then
      out_val = in_a
      return
   endif


   do k=1,num_of_EMP
      resonant_q(k) = real(EMP(k)%gam)/real(EMP(k)%gan)
      resonant_a(k) = Find_a_at_specific_safety_q(resonant_q(k))
      island_w(k) = 0.15
   enddo


   ! Get linear interpolation of piecewise base function as a background FLabel value
   out_val = 0
   do pwid=2,num_of_FLabel_pcws
      if( (in_a-FLabel_pcws_a(pwid-1))*(in_a-FLabel_pcws_a(pwid)) <= 0) then
         aa(1:2) = FLabel_pcws_a(pwid-1:pwid)
         yy(1:2) = FLabel_pcws_y(pwid-1:pwid)
   out_val = (yy(2)-yy(1))/(aa(2)-aa(1))*(in_a-aa(1)) +yy(1)
         exit
      endif
   enddo

   ! Add shaped values inside island
   do k=1,num_of_EMP
      if(abs(in_a-resonant_a(k)) < 0.5*island_w(k)) then
         MP_phase = EMP(k)%gan*in_gf- EMP(k)%gam*in_gq + EMP(k)%gagf0 !- EMP(k)%gagw*in_time
         out_val= out_val + (1.0-cos(MP_phase))*(abs(in_a-resonant_a(k))- 0.5*island_w(k))
         ! out_val= out_val + (1.0-cos(MP_phase))*((in_a-resonant_a(k))**2.0-(0.5*island_w(k))**2.0)
      endif
   enddo
end function Cal_FLabel_Function


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
   use mpi_defs

   implicit  none
   integer i, ip
   real(wp)  r, fip, dela, tempi
   real(wp) back_flux_A(mi), back_U_para, A_over_meanA, E_over_T
   real(wp) vpara, vperp
   real(wp) time0, time1

#ifdef _MGY_TIMER_
   time0 = MPI_WTIME()
#endif


   dela=(a1-a0)/real(jr,wp)
   dela=1.0_wp/dela

   if(acc_newer_H_zion) then
      !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
      acc_newer_H_zion=.false.
      acc_newer_D_zion=.false.
   endif
   if(plasma==1) then    ! use experimental ni(r), T(r) and omega(r) profiles
      if(linterp_from_nodes==1) then

         !$omp parallel do private(i,r,ip,fip)
         !$acc parallel loop present(zion,backn,backT,backu,dlnnda,dlnTda,duda,gh)
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
            back_flux_A(i)  =(fip*flux_A(ip+1)+(1.0_wp-fip)*flux_A(ip))
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
            back_flux_A(i)  =(fip*flux_A(ip+1)+(1.0_wp-fip)*flux_A(ip))
            ! A=exp(u^2/2T)
         enddo

      end if

      !$omp parallel do private(i)
      !$acc parallel loop present(dlnnda,backn,dlnTda,backT)
      do i=1,mi
         dlnnda(i)=dlnnda(i)/backn(i)      ! =dlnn/da
         dlnTda(i)=dlnTda(i)/backT(i)      ! =dlnT/da
      end do

   else                  ! use model profiles

      if(q2.gt.0.0) then  ! first analytic profile
         !$omp parallel do private(i,r,ip,fip)
         !$acc parallel loop present(zion,backn,backT,backu,dlnnda,dlnTda,duda,gh)
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

            back_flux_A(i)  =(fip*flux_A(ip+1)+(1.0_wp-fip)*flux_A(ip))
         end do
      else      ! second analytic profile
         !$omp parallel do private(i,r,ip,fip,tempi)
         !$acc parallel loop present(zion,backn,backT,backu,dlnnda,dlnTda,duda,gh)
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

            back_flux_A(i)  =(fip*flux_A(ip+1)+(1.0_wp-fip)*flux_A(ip))

         end do
      end if
   end if

   !$acc update self(backn(:mi),backT(:mi),backu(:mi),dlnnda(:mi),dlnTda(:mi),duda(:mi),gh(:mi))


   !     neglect rotation
   !       backu=0.0
   !       duda=0.0

   ! Calculate Shifted Maxwellian f0 for ions.
   call esigetbf(zion, mi, size(zion,1))

   !$acc parallel loop present(zion,backu,backT,back_f0_i,backn)
   do i=1, mi
      vpara = BB(i)*zion(4,i)  ! vpara = (eB/m)*rho_para
      vperp = sqrt(2.0*BB(i))*zion(6,i) ! vperp=sqrt(2B/m)*sqrt(mu)

      back_U_para = F(i)/BB(i)*backu(i)
      A_over_meanA = exp(0.5*(back_U_para**2.0)/backT(i))/(back_flux_A(i))
      E_over_T = (0.5*((vpara-back_U_para)**2.0+vperp**2.0))/backT(i)

      back_f0_i(i) = backn(i)*A_over_meanA*sqrt(inv_twopi)*((1.0/backT(i))**1.5)*vperp*exp(-E_over_T)
   enddo
   !$acc update self(back_f0_i(:mi))

#ifdef _MGY_TIMER_
   time1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "I'm=",mype, ", n_T_u_group takes", time1-time0, " secs"
   call flush(stdout)
#endif


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
   use mpi_defs

   implicit  none
   integer i, ip
   real(wp)  r, fip, dela, tempe
   real(wp) inv_Te_norm ! inverse of normalized Te which has same unit of Ti (Te_norm = Te_real/(mo*v0^2))
   real(wp) E_over_T
   real(wp) vpara, vperp
   real(wp) back_flux_A(me), omega_ui(me),backTi(me), back_Ui_para, A_over_meanA ! This is about given ion's rotation which determines theta dependency of density
   real(wp) time0, time1


#ifdef _MGY_TIMER_
   time0= MPI_WTIME()
#endif

   if(acc_newer_H_zeon) then
      !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
      acc_newer_H_zeon=.false.
      acc_newer_D_zeon=.false.
   endif

   if(plasma==1) then    ! use experimental ni(r), T(r) and omega(r) profiles
      if(linterp_from_nodes==1) then
         dela=(a1-a0)/real(jr,wp)
         dela=1.0_wp/dela

         !$omp parallel do private(i,r,ip,fip)
         !$acc parallel loop present(zeon,backTe,dlnTdae,backne,dlnndae)
         do i=1,me
            r=zeon(1,i)-a0
            r=r*dela
            ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
            fip=r-real(ip,wp)

            backne(i)=fip*ne_node(ip+1)+(1.0_wp-fip)*ne_node(ip) ! density at particles
            backTe(i)=fip*Te_node(ip+1)+(1.0_wp-fip)*Te_node(ip) ! Te at particles
            backTi(i)=fip*Ti_node(ip+1)+(1.0_wp-fip)*Ti_node(ip) ! Ti at particles

            dlnndae(i)=fip*dne_node(ip+1)+(1.0_wp-fip)*dne_node(ip) ! dn/da  at particles
            dlnTdae(i)=fip*dTe_node(ip+1)+(1.0_wp-fip)*dTe_node(ip) ! dTe/da at particles

            omega_ui(i)=fip*Iomg_node(ip+1)+(1.0_wp-fip)*Iomg_node(ip) ! omega_t at particles
            back_flux_A(i)  =(fip*flux_A(ip+1)+(1.0_wp-fip)*flux_A(ip))
         enddo

      else              ! use spline
         call ezspline_interp(spn,me,zeon(1,:),backne,ier)
         call ezspline_interp(spTe,me,zeon(1,:),backTe,ier)
         call ezspline_interp(spTi,mi,zeon(1,:),backTi,ier)
         call ezspline_interp(spO,me,zeon(1,:),omega_ui,ier)
         !            backu=backu*F/B
         call ezspline_derivative(spn,1,me,zeon(1,:),dlnndae,ier)
         call ezspline_derivative(spTe,1,me,zeon(1,:),dlnTdae,ier)

         !$acc parallel loop present(zeon)
         do i=1,me
            r=zeon(1,i)-a0
            r=r*dela
            ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
            fip=r-real(ip,wp)
            back_flux_A(i)  =(fip*flux_A(ip+1)+(1.0_wp-fip)*flux_A(ip)) ! A=exp(u^2/2T)
         enddo
      end if

      !$omp parallel do private(i)
      !$acc parallel loop present(backne,backTe,dlnndae,dlnTdae)
      do i=1,me
         dlnndae(i)=dlnndae(i)/backne(i)      !=dlnn/da
         dlnTdae(i)=dlnTdae(i)/backTe(i)      !=dlnT/da
      end do

   else                  ! use model profiles

      if(q2.gt.0.0) then  ! first analytic profile
         !$omp parallel do private(i)
         !$acc parallel loop present(zeon,backTe,dlnTdae,backne,dlnndae)
         do i=1,me
            backne(i)=n_0*Bn*exp(-dn*zeon(1,i)**xn)
            backTe(i)=T_e0*Bt_e*exp(-dTeme*zeon(1,i)**xt_e)
            backTi(i)=T_0*Bt*exp(-dTemi*zeon(1,i)**xt)

            omega_ui(i)=u_0*Bu*exp(-du*zeon(1,i)**xu)

            dlnndae(i)=backne(i)*(-dn*xn*zeon(1,i)**(xn-1.0_wp))       !=dn/da
            dlnTdae(i)=backTe(i)*(-dTeme*xt_e*zeon(1,i)**(xt_e-1.0_wp))    !=dT/da

            backne(i)=n_0*An+backne(i)
            backTe(i)=T_e0*At_e+backTe(i)
            backTi(i)=T_0*At+backTi(i)
            omega_ui(i)=u_0*Au+omega_ui(i)        ! then u_//=I*omega_ui/B


            dlnndae(i)=dlnndae(i)/backne(i)      !=dlnn/da
            dlnTdae(i)=dlnTdae(i)/backTe(i)      !=dlnT/da

            backTe(i)=Tc_e*(backTe(i)/Tc_e)**b0     ! for new model Ti profile
            backTi(i)=Tc*(backTi(i)/Tc)**r0     ! for new model Ti profile

            dlnTdae(i)=dlnTdae(i)*b0                ! for new model Ti profile

            r=zeon(1,i)-a0
            r=r*dela
            ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
            fip=r-real(ip,wp)
            back_flux_A(i)  =(fip*flux_A(ip+1)+(1.0_wp-fip)*flux_A(ip))
            !               backu=u_0*(Au+Bu*exp(-du*aa**xu))*(B_0/B)
         end do
      else      ! second analytic profile
         dela=(a1-a0)/real(jr,wp)
         dela=1.0_wp/dela
         !$omp parallel do private(i,r,ip,fip,tempe)
         !$acc parallel loop present(zeon,backTe,dlnTdae,backne,dlnndae)
         do i=1,me
            r=zeon(1,i)-a0
            r=r*dela
            ip=max(0,min(jr-1,int(r)))      ! ip [0,jr-1]
            fip=r-real(ip,wp)

            tempe=fip*Te_node(ip+1)+(1.0_wp-fip)*Te_node(ip)
            backTe(i)=T_e0*exp(flow1*tempe)                              ! Te at particles

            tempe=fip*Ti_node(ip+1)+(1.0_wp-fip)*Ti_node(ip)
            backTi(i)=T_0*exp(q1*tempe)                              ! Ti at particles


            tempe=flow2*(zeon(1,i)-rc)
            dlnTdae(i)=flow1*exp(-ABS(tempe)**flow0)                            ! =dlnT/da

            tempe=fip*ne_node(ip+1)+(1.0_wp-fip)*ne_node(ip)
            backne(i)=n_0*exp(kappan*tempe)                              ! ne at particles

            tempe=temconst*(zeon(1,i)-rc)
            dlnndae(i)=kappan*exp(-ABS(tempe)**edensity0)                     ! =dlnn/da


            tempe=fip*Iomg_node(ip+1)+(1.0_wp-fip)*Iomg_node(ip)
            omega_ui(i)=u_0+Au*tempe                                   ! this is omega_t; then u_//=I*backu/B

            back_flux_A(i)  =(fip*flux_A(ip+1)+(1.0_wp-fip)*flux_A(ip))
            !               backu=u_0*(Au+Bu*exp(-du*aa**xu))*(B_0/B)
         end do
      end if

   end if

   !$acc update self(backTe(:me), dlnTdae(:me), backne(:me), dlnndae(:me))


   backue=0._wp    ! neglect rotation


   ! Calculate Maxwellian f0 for electrons.
   call esigetb(zeon, me, size(zeon,1))

   !$acc parallel loop present(zeon,backne,backTe,back_f0_e)
   do i=1, me
      ! Following part is to calculate theta dependency of density which needs ion rotation information
      back_Ui_para = F(i)/BB(i)*omega_ui(i)
      A_over_meanA = exp(0.5*(back_Ui_para**2.0)/backTi(i))/(back_flux_A(i))

      inv_Te_norm = 1.0/(backTe(i)*msm0)  ! normalized-Te=Te/E_0, the same unit as normalized Ti
      vpara = m0ms*csc0*BB(i)*zeon(4,i)  ! vpara = (eB/m)*rho_para
      vperp = sqrt(2.0*BB(i)*m0ms)*zeon(6,i) ! vperp=sqrt(2B/m)*sqrt(mu)

      E_over_T= (0.5*msm0*(vpara**2.0+vperp**2.0))*inv_Te_norm

      back_f0_e(i) = backne(i)*A_over_meanA*sqrt(inv_twopi)*((msm0*inv_Te_norm)**1.5)*vperp*exp(-E_over_T)
   enddo
   !$acc update self(back_f0_e(:me))

#ifdef _MGY_TIMER_
   time1=MPI_WTIME()
   if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "I'm=",mype, ", n_T_u_group_e takes", time1-time0, " secs"
   call flush(stdout)
#endif

end subroutine n_T_u_group_e



subroutine Cumulate_Some_Variables_for_Measure

   use global_parameters
   use field_array
   implicit none

   integer i,j,ij
   real(wp) full_dphi(0:1), full_dni(0:1), full_dne(0:1), adiab_dne_n0(0:1), full_dTe, full_dTi




   cum_phi00       = cum_phi00    +phi00
   cum_phi00_SQ    = cum_phi00_SQ +phi00**2.0

   cum_zonali      = cum_zonali    +zonali
   cum_zonali_SQ   = cum_zonali_SQ +zonali**2.0

   cum_zonale    = cum_zonale    +zonale
   cum_zonale_SQ = cum_zonale_SQ +zonale**2.0


   do i=0,mpsi
      do j=0,mtheta(i)
         ij=igrid(i)+j

         full_dphi     = phi(:,ij)+(ZF_SEPARATION)*phi00(i)
         full_dni      = den_grid(i)*(densityi(:,ij) +(ZF_SEPARATION)*zonali(i))
         adiab_dne_n0  = abs(csc0)*(full_dphi-(ZF_SEPARATION)*phi00(i))/Te_grid(i)
         full_dne      = den_grid(i)*(densitye(:,ij) +(ZF_SEPARATION)*zonale(i) +int(ifelectron/2)*adiab_dne_n0)

         cum_full_dphi(:,ij)    = cum_full_dphi(:,ij)    +full_dphi
         cum_full_dphi_SQ(:,ij) = cum_full_dphi_SQ(:,ij) +full_dphi**2.0

         cum_full_dni(:,ij)    = cum_full_dni(:,ij)    +full_dni
         cum_full_dni_SQ(:,ij) = cum_full_dni_SQ(:,ij) +full_dni**2.0

         cum_full_dne(:,ij)    = cum_full_dne(:,ij)    +full_dne
         cum_full_dne_SQ(:,ij) = cum_full_dne_SQ(:,ij) +full_dne**2.0
      enddo
   enddo


   if(mod(istep,ndiag)==0) then
      cum_zonal_dTe    = cum_zonal_dTe    +zonal_dTe
      cum_zonal_dTe_SQ = cum_zonal_dTe_SQ +zonal_dTe**2.0
      cum_zonal_dTi     = cum_zonal_dTi   +zonal_dTi
      cum_zonal_dTi_SQ  = cum_zonal_dTi_SQ+zonal_dTi**2.0

      do i=0,mpsi
         do j=0,mtheta(i)
            ij=igrid(i)+j

            full_dTe = Te_grid(i)*heatflux(ij) +(ZF_SEPARATION)*zonal_dTe(i)
            full_dTi = (vth_grid(i)**2.0)*phi_tmp(ij) +(ZF_SEPARATION)*zonal_dTi(i)

            cum_full_dTe(:,ij)    =cum_full_dTe(:,ij)   +full_dTe
            cum_full_dTe_SQ(:,ij) =cum_full_dTe_SQ(:,ij)+full_dTe**2.0

            cum_full_dTi(:,ij)    =cum_full_dTi(:,ij)   +full_dTi
            cum_full_dTi_SQ(:,ij) =cum_full_dTi_SQ(:,ij)+full_dTi**2.0
         enddo
      enddo
   endif

end subroutine Cumulate_Some_Variables_for_Measure

subroutine SetZeros_Cumulative_Variables
   use field_array
   implicit none

   cum_phi00       =0.0
   cum_phi00_SQ    =0.0
   cum_full_dphi    =0.0
   cum_full_dphi_SQ =0.0

   cum_zonali      =0.0
   cum_zonali_SQ   =0.0
   cum_full_dni    =0.0
   cum_full_dni_SQ =0.0

   cum_zonale      =0.0
   cum_zonale_SQ   =0.0
   cum_full_dne    =0.0
   cum_full_dne_SQ =0.0

   cum_zonal_dTe   =0.0
   cum_zonal_dTe_SQ=0.0
   cum_full_dTe    =0.0
   cum_full_dTe_SQ =0.0

   cum_zonal_dTi   =0.0
   cum_zonal_dTi_SQ=0.0
   cum_full_dTi    =0.0
   cum_full_dTi_SQ =0.0

end subroutine SetZeros_Cumulative_Variables

subroutine Do_FieldLineFollowing_for_Poincare_Plot

   use precision
   use global_parameters
   use particle_decomp
   use FieldLineFollowing

   implicit  none

   character(len=64) chdum

   real(wp) out_ga, out_grad_ga(3), out_partial_dg0dt
   real(wp) B_vec_cntv(3)

   real(wp) dgf ! toroidal step for field-line-following integration
   real(wp) a_RK(4), gq_RK(4)
   real(wp) ini_a, ini_gq, ini_gf
   real(wp) now_a, now_gq, now_gf, specific_time
   real(wp) prev_a, prev_gq, prev_gf
   integer :: maximum_tor_rotation
   integer tor_dir, aid, num_a, gqid, num_gq, PPI_id
   real(wp):: aaa, the, Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
   type(Poincare_Point_Information), dimension(:,:,:), allocatable:: PPI
   real(wp) intp_weight, dummy_real
   real(wp) t_start, t_finish
   logical ex
   real(wp), external :: Cal_FLabel_Function

   call cpu_time(t_start)



   ! if(irun/=0) then
   if(.true.) then

      ! In the case of irun!=0, if there is previously calculated FLF_Poincare data,
      ! then just return this function to skip calculating Poincare Plot
      write(chdum,'("VacFieldAnalysis/FLF_Poincare",i5.5,".dat")') myrank_toroidal
      inquire(file=chdum, exist=ex)
      if(ex) then
         return
      endif
   endif

   ! If there is no PPI data, then calculate the following routines
   maximum_tor_rotation= FLF_max_tor_num
   specific_time = 0.0 ! time is constant because we want the field Poincare Plot at specific time t=t0

   dgf = 2.0*pi/500.0 ! toroidal step for field-line-following integration
   num_a = 40
   num_gq = 7

   allocate( PPI(2*maximum_tor_rotation, num_gq, num_a))
   !PPI(:,:,:) = 0

   !write(*,*) "Field Line Integration start!", myrank_partd, myrank_toroidal

   do aid=1,num_a
      do gqid = 1, num_gq


         ini_a = a0 + (aid-1)*(a1-a0)/(num_a-1.0)
         ini_gq = (gqid-1)*(2.0*pi)/(num_gq) ! To avoid overlap between gq=0 and gq=2pi points
         ini_gf = 0.0


         aaa=ini_a
         the=ini_gq
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
         !PPI(:, gqid, aid)%Jacobian = rj*(zaj*rqj-zqj*raj); ! Jacobian at initial position.
         PPI(:, gqid, aid)%Fid = (aid-1)*num_gq + (gqid-1)

         PPI_id = 0

         do tor_dir=-1,1,2
            dgf = (tor_dir)*2.0*pi/500.0

            now_a = ini_a
            now_gq = ini_gq
            now_gf = ini_gf

            prev_a = now_a
            prev_gq = now_gq
            prev_gf = now_gf

            do while(.true.)

               ! When B-line-integrated point just passed gf=zetamax poloidal plane.
               if( floor((now_gf-zetamax)/(2.0*pi)) /= floor((prev_gf-zetamax)/(2.0*pi)) ) then
                  PPI_id = PPI_id + 1

                  !intp_weight = modulo(abs(now_gf), 2.0*pi) / dgf
                  intp_weight = abs(min(modulo((now_gf-zetamax), 2.0*pi), 2.0*pi - modulo((now_gf-zetamax), 2.0*pi))/dgf) ! interpolation weight for previous point

                  PPI(PPI_id, gqid, aid)%a  = intp_weight*prev_a + (1.0-intp_weight)*now_a
                  PPI(PPI_id, gqid, aid)%gq = intp_weight*prev_gq + (1.0-intp_weight)*now_gq
                  PPI(PPI_id, gqid, aid)%gf = intp_weight*prev_gf + (1.0-intp_weight)*now_gf
                  PPI(PPI_id, gqid, aid)%target_normalized_gf = zetamax
                  PPI(PPI_id, gqid, aid)%ini_FLabel = Cal_FLabel_Function(PPI(PPI_id,gqid,aid)%a, PPI(PPI_id,gqid,aid)%gq, PPI(PPI_id,gqid,aid)%gf)

                  call Convert_a_to_s(PPI(PPI_id, gqid, aid)%a, PPI(PPI_id, gqid, aid)%s)

                  !call esigetrzb(PPI(PPI_id, gqid, aid)%R, PPI(PPI_id, gqid, aid)%Z, dummy_real, PPI(PPI_id, gqid, aid)%a, modulo(PPI(PPI_id, gqid, aid)%gq, 2*pi), 1)

                  aaa=PPI(PPI_id, gqid, aid)%a
                  the=modulo(PPI(PPI_id, gqid, aid)%gq, 2.0*pi)
                  call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
                  PPI(PPI_id, gqid, aid)%inv_B = 1.0/Bj; ! Jacobian at initial position.
                  !PPI(PPI_id, gqid, aid)%inv_B = 1.0;
                  PPI(PPI_id, gqid, aid)%R = rj
                  PPI(PPI_id, gqid, aid)%Z = zj

                  call esiget_Intp_B_vec_contravariant(aaa, the, modulo(PPI(PPI_id, gqid, aid)%gf,2.0*pi), specific_time, B_vec_cntv)
                  PPI(PPI_id, gqid, aid)%dl =abs(dgf)*sqrt( ((B_vec_cntv(1)/B_vec_cntv(3))**2)*(raj**2+zaj**2)  &
                     +((B_vec_cntv(2)/B_vec_cntv(3))**2)*(rqj**2+zqj**2)  &
                     +(B_vec_cntv(1)*B_vec_cntv(2)/(B_vec_cntv(3)**2))*(raj*rqj+zaj*zqj) &
                     + rj**2  )
               endif

               if(abs(now_gf)>maximum_tor_rotation*(2*pi) .or. PPI_id>= 2*maximum_tor_rotation) then
                  exit
               endif

               prev_a = now_a
               prev_gq = now_gq
               prev_gf = now_gf

               ! Runge-Kutta 4th order
               call esiget_Intp_B_vec_contravariant(now_a, modulo(now_gq,2.0*pi), modulo(now_gf,2.0*pi), specific_time, B_vec_cntv)
               a_RK(1)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
               gq_RK(1)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
               call esiget_Intp_B_vec_contravariant(now_a+0.5*a_RK(1), modulo(now_gq+0.5*gq_RK(1),2.0*pi), modulo(now_gf+0.5*dgf,2.0*pi), specific_time, B_vec_cntv)
               a_RK(2)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
               gq_RK(2)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
               call esiget_Intp_B_vec_contravariant(now_a+0.5*a_RK(2), modulo(now_gq+0.5*gq_RK(2),2.0*pi), modulo(now_gf+0.5*dgf,2.0*pi), specific_time, B_vec_cntv)
               a_RK(3)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
               gq_RK(3)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
               call esiget_Intp_B_vec_contravariant(now_a+a_RK(3), modulo(now_gq+gq_RK(3),2.0*pi), modulo(now_gf+dgf,2.0*pi), specific_time, B_vec_cntv)
               a_RK(4)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
               gq_RK(4)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv

               now_a = now_a + (a_RK(1)+2.0*a_RK(2)+2.0*a_RK(3)+a_RK(4))/6.0
               now_gq = now_gq + (gq_RK(1)+2.0*gq_RK(2)+2.0*gq_RK(3)+gq_RK(4))/6.0
               now_gf = now_gf + dgf

            end do
         end do
      end do
   end do

   call cpu_time(t_finish)
   write(stdout,*) "I'm=", mype, "FLF_Poincare Finished! delta_time = ", t_finish-t_start, " seconds."
   flush(stdout)

   call system("mkdir -p " // 'VacFieldAnalysis')
   write(chdum,'("VacFieldAnalysis/FLF_Poincare",i5.5,".dat")') myrank_toroidal
   open(181125,file=chdum,status='replace',action='write')
   do aid=1,num_a
      do gqid = 1, num_gq
         do PPI_id=1,2*maximum_tor_rotation
            write(181125, "(I8,1x,20(ES16.8,1x))") PPI(PPI_id,gqid,aid)
         enddo
      enddo
   enddo
   flush(181125)
   close(181125)

   !write(*,*) "Do_FieldLineFollowing_for_Poincare_Plot OK"


   deallocate(PPI)

end subroutine Do_FieldLineFollowing_for_Poincare_Plot


subroutine Construct_gPPI_for_FieldLineAverage


   use precision
   use global_parameters
   use particle_decomp
   use FieldLineFollowing
   use esi_interface
   use field_array
   implicit  none

   character(len=64) chdum
   ! Global_PPI(num_of_maximum_intersections, 1:num_of_SampledFieldLines)
   real(wp) dgf ! toroidal step for field-line-following integration
   real(wp) a_RK(4), gq_RK(4)

   real(wp) ini_a, ini_gq, ini_gf
   real(wp),dimension(:), allocatable:: ini_a_arr, ini_gq_arr
   integer target_mode_m, target_mode_n
   !real(wp) :: target_safety_q(9) = (/0.5,0.75,1.0,1.3333,1.5,1.6666,2.0,2.5,3.0/)
   real(wp) target_safety_q(200)! Arrange enough initial length for target_safety_q. Later, effective num_of_target_saftegy_q will be calculated.
   integer num_of_target_safety_q
   integer sampling_num_at_each_rational_q, i,j
   integer MaxNum_of_gPPI_for_each_FieldLine
   real(wp) Find_a_at_specific_safety_q
   real(wp) now_a, now_gq, now_gf, specific_time
   real(wp) prev_a, prev_gq, prev_gf
   real(wp) out_ga, out_grad_ga(3), out_partial_dg0dt
   real(wp) B_vec_cntv(3)
   integer tor_dir, PPI_id, Fid, maximum_tor_rotation
   integer(kind=8) :: sobol_seed=1, sobol_dim_num=2
   real(kind=8) sobol_rand_2d(2)
   real(wp) intp_weight, dummy_real, tmp_rand
   real(wp):: aaa, the, Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
   integer, parameter:: num_of_target_zeta=1
   logical record_flag
   integer :: k, tz_id
   real(wp)::target_zeta(num_of_target_zeta)
   real(wp) t_start, t_finish
   logical ex ! existence of the reading file
   integer ios! IOSTAT of the reading file
   type(Poincare_Point_Information) :: tmp_gPPI ! For reading purpose
   real(wp) target_q, resonant_a, MP_ga, MP_wga, MP_gat ! Magnetic perturbation variables
   real(wp), external :: Cal_FLabel_Function


   ! if(irun/=0) then
   if(.true.) then
      ! First, check if there is previously calculated gPPI data, and read it
      write(chdum,'("VacFieldAnalysis/gPPI_",i5.5,".dat")') myrank_toroidal
      inquire(file=chdum, exist=ex)

      if(ex) then
         open(190704,file="VacFieldAnalysis/gPPI_Information.dat",status='old',form='formatted',action='read')
         read(190704,"(3(I8,1x))") gPPI_max_tor_num, MaxNum_of_gPPI_for_each_FieldLine, num_of_SampledFieldLines
         close(190704)

         allocate(gPPI(MaxNum_of_gPPI_for_each_FieldLine, num_of_SampledFieldLines))
         allocate(num_of_gPPI_for_each_FieldLine(num_of_SampledFieldLines))

         open(190705,file=chdum,status='old',form='formatted',action='read')

         Fid=1
         PPI_id=0
         do while(.true.)
            read(190705,"(I8,1x,11(ES16.8,1x))",IOSTAT=ios) tmp_gPPI
            if(ios>0) then ! something wrong with ios
               write(stdout,*) "ios=",ios,", something is wrong during reading ",chdum
               call flush(stdout)
               stop
            elseif(ios<0) then ! end of file reached
               num_of_gPPI_for_each_FieldLine(Fid)=PPI_id
               exit
            else
               if(Fid /= tmp_gPPI%Fid) then
                  num_of_gPPI_for_each_FieldLine(Fid)=PPI_id
                  Fid = Fid+1
                  PPI_id=0
               endif
               PPI_id=PPI_id+1
               gPPI(PPI_id,Fid) = tmp_gPPI
            endif
         enddo
         close(190705)

         ! This is required to calculate Do_FieldLineAverage_by_using_gPPI later
         call Cal_FLA_Denominator_using_gPPI
         call Get_FLA_Te0_Ti0_n0

         return
      endif
   endif

   ! If there is no gPPI data, then calculate the following routines
   call cpu_time(t_start)

   maximum_tor_rotation = gPPI_max_tor_num
   MaxNum_of_gPPI_for_each_FieldLine = 2*(num_of_target_zeta)*maximum_tor_rotation
   !num_of_SampledFieldLines = (mpsi/4)*10  ! num_of_SampledFieldLines might depend more on the radial resolution than poloidal one.

   !maximum_tor_rotation = MaxNum_of_gPPI_for_each_FieldLine/2 !  Maximum toroidal roation for each Forward integration & Backward Integration

   specific_time = 0.0 ! time is constant because we want the field Poincare Plot at specific time t=t0
   dgf = 2.0*pi/500.0 ! toroidal step for field-line-following integration

   do k=1,num_of_target_zeta
      target_zeta(k) = zetamin + k*(zetamax-zetamin)/num_of_target_zeta
   enddo
   !write(stdout,*) 'mype=', mype, '(zetamin,zetamax)=', zetamin, zetamax,', target_zeta=', target_zeta(1:num_of_target_zeta)
   !call flush(stdout)
   k=1
   do target_mode_n = 1,3
      do target_mode_m = 1,10
         if(target_mode_m <= 3*target_mode_n) then ! search for the safety_q<=3
            target_safety_q(k) = real(target_mode_m)/real(target_mode_n)
            k=k+1
         endif
         ! target_safety_q(k) = 0.95*real(target_mode_m)/real(target_mode_n)
         ! k=k+1
         ! target_safety_q(k) = 1.05*real(target_mode_m)/real(target_mode_n)
         ! k=k+1
      enddo
   enddo
   !call flush(stdout)
   call find_unique_sorted_array_and_its_num_of_elements(target_safety_q, 1, k-1, num_of_target_safety_q)
   if(myrank_toroidal==0) then
      write(stdout,*) 'num_of_target_safety_q = ', num_of_target_safety_q
      write(stdout,*) 'target q vals = ',target_safety_q(1:num_of_target_safety_q)
      call flush(stdout)
   endif

   sampling_num_at_each_rational_q = 50
   !sampling_num_at_each_rational_q = 30
   ! First, find target rational safety q positions for accurate investigation
   ! Fid = 1
   ! do target_mode_n = 1,3
   !   do target_mode_m = 1,10
   !     ini_a = Find_a_at_specific_safety_q(target_mode_m/real(target_mode_n))
   !     if(ini_a>0) then
   !       Fid = Fid + 1
   !     endif
   !   enddo
   ! enddo

   Fid=1
   do i = 1,num_of_target_safety_q
      ini_a = Find_a_at_specific_safety_q(target_safety_q(i))
      if(ini_a>0) then
         !Fid = Fid + 1
         Fid = Fid + 3 ! 3 field lines will be sampled at each target_saftey_q and their neighborhoods.
      endif
   enddo

   num_of_SampledFieldLines = 50*20 + (sampling_num_at_each_rational_q*(Fid-1))

   allocate(gPPI(MaxNum_of_gPPI_for_each_FieldLine, num_of_SampledFieldLines))
   allocate(ini_a_arr(num_of_SampledFieldLines), ini_gq_arr(num_of_SampledFieldLines))
   allocate(num_of_gPPI_for_each_FieldLine(num_of_SampledFieldLines))


   ! ! Find ini_a_arr and ini_gq_arr at rational q surface and sobol_rand position
   ! Fid = 1
   ! do target_mode_n = 1,3
   !   do target_mode_m = 1,10
   !     ini_a = Find_a_at_specific_safety_q(target_mode_m/real(target_mode_n))
   !     if(ini_a>0) then
   !       ini_a_arr(Fid:Fid+(sampling_num_at_each_rational_q-1)) = ini_a
   !       do i=0,sampling_num_at_each_rational_q-1
   !         ini_gq_arr(Fid+i) = i*(2.0*pi/real(sampling_num_at_each_rational_q))
   !       enddo
   !       Fid = Fid + sampling_num_at_each_rational_q
   !     endif
   !   enddo
   ! enddo


   !Find ini_a_arr and ini_gq_arr at rational q surface and sobol_rand position
   Fid = 1
   do i = 1, num_of_target_safety_q
      ini_a = Find_a_at_specific_safety_q(target_safety_q(i))
      if(ini_a>0) then
         ! At exact target_safety_q
         ini_a_arr(Fid:Fid+(sampling_num_at_each_rational_q-1)) = ini_a
         do j=0,sampling_num_at_each_rational_q-1
            ini_gq_arr(Fid+j) = j*(2.0*pi/real(sampling_num_at_each_rational_q))
         enddo
         Fid = Fid + sampling_num_at_each_rational_q

         ! Sligtly left from target_safety_q
         ini_a_arr(Fid:Fid+(sampling_num_at_each_rational_q-1)) = ini_a - 0.01
         do j=0,sampling_num_at_each_rational_q-1
            ini_gq_arr(Fid+j) = j*(2.0*pi/real(sampling_num_at_each_rational_q))
         enddo
         Fid = Fid + sampling_num_at_each_rational_q

         ! Slightly right from target_safety_q
         ini_a_arr(Fid:Fid+(sampling_num_at_each_rational_q-1)) = ini_a + 0.01
         do j=0,sampling_num_at_each_rational_q-1
            ini_gq_arr(Fid+j) = j*(2.0*pi/real(sampling_num_at_each_rational_q))
         enddo
         Fid = Fid + sampling_num_at_each_rational_q
      endif
   enddo


   do while(Fid<=num_of_SampledFieldLines)
      call i8_sobol(sobol_dim_num, sobol_seed, sobol_rand_2d)
      ini_a_arr(Fid) = a0 + (a1-a0)*sobol_rand_2d(1)
      ini_gq_arr(Fid) = 2.0*pi*sobol_rand_2d(2)
      Fid = Fid+1
   enddo

   !write(stdout,*) 'mype=', mype,', before quicksort_real_together'
   call flush(stdout)
   call quicksort_real_together(ini_a_arr, ini_gq_arr, 1,num_of_SampledFieldLines)
   !write(stdout,*) 'mype=', mype,', after quicksort_real_together'
   !write(stdout,*) 'mype=', mype,', ini_a_arr(1:3)', ini_a_arr(1:3), ', ini_a_arr(end-3:end)', ini_a_arr(num_of_SampledFieldLines-3:num_of_SampledFieldLines)
   !write(stdout,*) 'mype=', mype,', ini_gq_arr(1:3)', ini_gq_arr(1:3), ', ini_gq_arr(end-3:end)', ini_gq_arr(num_of_SampledFieldLines-3:num_of_SampledFieldLines)
   call flush(stdout)

   ! First Resonant q and a (only the First element of EMP)
   resonant_a = Find_a_at_specific_safety_q(real(EMP(1)%gam)/real(EMP(1)%gan))
   write(stdout,*) '1st resonant_q=',real(EMP(1)%gam)/real(EMP(1)%gan), 'resonant_a=',resonant_a
   call flush(stdout)

   do Fid = 1, num_of_SampledFieldLines

      ! ini_a = a0 + (a1-a0)*sobol_rand_2d(1)
      ! call i8_sobol(2, sobol_seed, sobol_rand_2d)
      ! ini_gq = 2.0*pi*sobol_rand_2d(2)
      ! call random_number(tmp_rand)
      ! ini_a = a0 + (a1-a0)*tmp_rand
      ! call random_number(tmp_rand)
      ! ini_gq = 2.0*pi*tmp_rand
      ini_a = ini_a_arr(Fid)
      ini_gq = ini_gq_arr(Fid)
      ini_gf = 0


      aaa=ini_a
      the=ini_gq
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      !gPPI(:, Fid)%Jacobian = rj*(zaj*rqj-zqj*raj); ! Jacobian at initial position.
      gPPI(:, Fid)%Fid = Fid


      PPI_id = 0

      do tor_dir=-1,1,2  ! Forward integration & Backward Integration

         dgf = (tor_dir)*2.0*pi/500.0

         now_a = ini_a
         now_gq = ini_gq
         now_gf = ini_gf

         prev_a = now_a
         prev_gq = now_gq
         prev_gf = now_gf

         do while(.true.)
            record_flag=.false.
            !if( modulo(now_gf,2*pi)>zetamin-abs(dgf) .and.  modulo(now_gf,2*pi)<=zetamax+abs(dgf)) then
            do k=1,num_of_target_zeta
               if(PPI_id==0 .or. floor((now_gf-target_zeta(k))/(2.0*pi)) /= floor((prev_gf-target_zeta(k))/(2.0*pi))) then
                  record_flag=.true.
                  tz_id = k
                  exit
               endif
            enddo
            !endif

            if(record_flag) then
               !if( floor((now_gf-zetamax)/(2.0*pi)) /= floor((prev_gf-zetamax)/(2.0*pi)) ) then
               PPI_id = PPI_id + 1

               gPPI(PPI_id, Fid)%Fid = Fid
               !intp_weight = modulo(abs(now_gf), 2.0*pi) / dgf
               !intp_weight = abs(min(modulo((now_gf-zetamax), 2.0*pi), 2.0*pi - modulo((now_gf-zetamax), 2.0*pi))/dgf) ! interpolation weight for previous point
               intp_weight = abs(min(modulo((now_gf-target_zeta(tz_id)), 2.0*pi), 2.0*pi - modulo((now_gf-target_zeta(tz_id)), 2.0*pi))/dgf) ! interpolation weight for previous point

               gPPI(PPI_id, Fid)%a  = intp_weight*prev_a + (1.0-intp_weight)*now_a
               gPPI(PPI_id, Fid)%gq = intp_weight*prev_gq + (1.0-intp_weight)*now_gq
               gPPI(PPI_id, Fid)%gf = intp_weight*prev_gf + (1.0-intp_weight)*now_gf
               gPPI(PPI_id, Fid)%target_normalized_gf = target_zeta(tz_id)
               gPPI(PPI_id, Fid)%ini_FLabel = Cal_FLabel_Function(gPPI(PPI_id, Fid)%a, gPPI(PPI_id, Fid)%gq, gPPI(PPI_id, Fid)%gf)

               call Convert_a_to_s(gPPI(PPI_id, Fid)%a, gPPI(PPI_id, Fid)%s)

               !call esigetrzb(gPPI(PPI_id, Fid)%R, gPPI(PPI_id, Fid)%Z, dummy_real, gPPI(PPI_id, Fid)%a, modulo(gPPI(PPI_id, Fid)%gq, 2*pi), 1)

               aaa=gPPI(PPI_id, Fid)%a
               the=modulo(gPPI(PPI_id, Fid)%gq, 2.0*pi)
               call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
               gPPI(PPI_id, Fid)%inv_B = 1.0/Bj; ! Jacobian at initial position.
               !gPPI(PPI_id, Fid)%inv_B = 1.0; ! Jacobian at initial position.
               gPPI(PPI_id, Fid)%R = rj
               gPPI(PPI_id, Fid)%Z = zj

               call esiget_Intp_B_vec_contravariant(aaa, the, modulo(gPPI(PPI_id, Fid)%gf,2.0*pi), specific_time, B_vec_cntv)
               gPPI(PPI_id, Fid)%dl = abs(dgf)*sqrt( ((B_vec_cntv(1)/B_vec_cntv(3))**2)*(raj**2+zaj**2)  &
                  +((B_vec_cntv(2)/B_vec_cntv(3))**2)*(rqj**2+zqj**2)  &
                  +(B_vec_cntv(1)*B_vec_cntv(2)/(B_vec_cntv(3)**2))*(raj*rqj+zaj*zqj) &
                  + rj**2  )
            endif

            if(abs(now_gf)>maximum_tor_rotation*(2*pi) .or.  PPI_id>=MaxNum_of_gPPI_for_each_FieldLine) then
               exit
            endif

            prev_a = now_a
            prev_gq = now_gq
            prev_gf = now_gf

            ! Runge-Kutta 4th order
            call esiget_Intp_B_vec_contravariant(now_a, modulo(now_gq,2.0*pi), modulo(now_gf,2.0*pi), specific_time, B_vec_cntv)
            a_RK(1)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
            gq_RK(1)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
            call esiget_Intp_B_vec_contravariant(now_a+0.5*a_RK(1), modulo(now_gq+0.5*gq_RK(1),2.0*pi), modulo(now_gf+0.5*dgf,2.0*pi), specific_time, B_vec_cntv)
            a_RK(2)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
            gq_RK(2)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
            call esiget_Intp_B_vec_contravariant(now_a+0.5*a_RK(2), modulo(now_gq+0.5*gq_RK(2),2.0*pi), modulo(now_gf+0.5*dgf,2.0*pi), specific_time, B_vec_cntv)
            a_RK(3)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
            gq_RK(3)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
            call esiget_Intp_B_vec_contravariant(now_a+a_RK(3), modulo(now_gq+gq_RK(3),2.0*pi), modulo(now_gf+dgf,2.0*pi), specific_time, B_vec_cntv)
            a_RK(4)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
            gq_RK(4)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv

            now_a = now_a + (a_RK(1)+2.0*a_RK(2)+2.0*a_RK(3)+a_RK(4))/6.0
            now_gq = now_gq + (gq_RK(1)+2.0*gq_RK(2)+2.0*gq_RK(3)+gq_RK(4))/6.0
            now_gf = now_gf + dgf

            if (now_a<a0 .or. now_a>a1) then
               exit
            endif
         enddo
      enddo

      num_of_gPPI_for_each_FieldLine(Fid) = PPI_id
   enddo

   ! write(stdout,*) 'mype=', mype,', before Cal_FLA_Denominator_using_gPPI'
   !call flush(stdout)
   !test to get test_FLA_a
   call Cal_FLA_Denominator_using_gPPI

   !write(stdout,*) 'mype=', mype,', before Get_FLA_Te0_Ti0_n0'
   !call flush(stdout)
   call Get_FLA_Te0_Ti0_n0


   call cpu_time(t_finish)
   write(stdout,*) "I'm=", mype, "gPPI Construction Finished! delta_time = ", t_finish-t_start, " seconds."
   flush(stdout)

   call system("mkdir -p " // 'VacFieldAnalysis')
   write(chdum,'("VacFieldAnalysis/gPPI_",i5.5,".dat")') myrank_toroidal
   open(181202,file=chdum,status='replace',action='write')
   do Fid = 1, num_of_SampledFieldLines
      do PPI_id =1,num_of_gPPI_for_each_FieldLine(Fid)
         write(181202, "(I8,1x,11(ES16.8,1x))") gPPI(PPI_id,Fid)
      enddo
   enddo
   flush(181202)
   close(181202)


   open(190704,file="VacFieldAnalysis/gPPI_Information.dat",status='replace',form='formatted',action='write')
   write(190704,"(3(I8,1x))") gPPI_max_tor_num, MaxNum_of_gPPI_for_each_FieldLine, num_of_SampledFieldLines
   close(190704)

   deallocate(ini_a_arr, ini_gq_arr)


end subroutine Construct_gPPI_for_FieldLineAverage



subroutine Construct_gPPI_for_FieldLineAverage_Parallelized

   use precision
   use global_parameters
   use particle_decomp
   use FieldLineFollowing
   use esi_interface
   use field_array
   use mpi_defs

   implicit  none

   character(len=128) chdum, chdum2
   real(wp) dgf ! toroidal step for field-line-following integration
   real(wp) a_RK(4), gq_RK(4)

   real(wp) ini_a, ini_gq, ini_gf
   real(wp),dimension(:), allocatable:: ini_a_arr, ini_gq_arr, ini_gf_arr
   integer target_mode_m, target_mode_n
   !real(wp) :: target_safety_q(9) = (/0.5,0.75,1.0,1.3333,1.5,1.6666,2.0,2.5,3.0/)
   real(wp) target_safety_q(200)! Arrange enough initial length for target_safety_q. Later, effective num_of_target_saftegy_q will be calculated.
   integer num_of_target_safety_q
   integer sampling_num_at_each_rational_q, i,j
   integer MaxNum_of_gPPI_for_each_FieldLine
   real(wp) Find_a_at_specific_safety_q
   real(wp) now_a, now_gq, now_gf, specific_time
   real(wp) prev_a, prev_gq, prev_gf
   real(wp) out_ga, out_grad_ga(3), out_partial_dg0dt
   real(wp) B_vec_cntv(3)
   integer tor_dir, Pid, Fid, maximum_tor_rotation
   integer(kind=8) :: sobol_seed=1, sobol_dim_num=2
   real(kind=8) sobol_rand(2)
   real(wp) intp_weight, dummy_real, tmp_rand
   real(wp):: aaa, the, Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
   integer :: num_of_target_zeta
   real(wp),allocatable,dimension(:)::target_zeta
   logical record_flag
   integer :: k, tzid, tor_id, ierr, length
   real(wp) t_start, t_finish
   logical ex ! existence of the reading file
   integer ios! IOSTAT of the reading file
   type(Poincare_Point_Information) :: tmp_single_gPPI ! For reading purpose
   type(Poincare_Point_Information),allocatable, dimension(:,:,:) :: tmp_gPPI ! For parallel calculation of gPPI
   integer, dimension(:,:), allocatable :: tmp_num_of_gPPI_for_each_FieldLine
   integer tmp_num_of_SampledFieldLines

   real(wp) target_q, resonant_a, MP_ga, MP_wga, MP_gat ! Magnetic perturbation variables
   real(wp), external :: Cal_FLabel_Function


   ! if(irun/=0) then
   if(.true.) then
      ! First, check if there is previously calculated gPPI data, and if so read it
      write(chdum,'("VacFieldAnalysis/gPPI_",i5.5,".dat")') myrank_toroidal
      inquire(file=chdum, exist=ex)

      if(ex) then
         open(190704,file="VacFieldAnalysis/gPPI_Information.dat",status='old',form='formatted',action='read')
         read(190704,"(3(I8,1x))") gPPI_max_tor_num, MaxNum_of_gPPI_for_each_FieldLine, num_of_SampledFieldLines
         close(190704)

         allocate(gPPI(MaxNum_of_gPPI_for_each_FieldLine, num_of_SampledFieldLines))
         allocate(num_of_gPPI_for_each_FieldLine(num_of_SampledFieldLines))

         open(190705,file=chdum,status='old',form='formatted',action='read')

         Fid=1
         Pid=0
         do while(.true.)
            read(190705,"(I8,1x,11(ES16.8,1x))",IOSTAT=ios) tmp_single_gPPI
            if(ios>0) then ! something wrong with ios
               write(stdout,*) "ios=",ios,", something is wrong during reading ",chdum
               call flush(stdout)
               stop
            elseif(ios<0) then ! end of file reached
               num_of_gPPI_for_each_FieldLine(Fid)=Pid
               exit
            else
               if(Fid /= tmp_single_gPPI%Fid) then
                  num_of_gPPI_for_each_FieldLine(Fid)=Pid
                  Fid = Fid+1
                  Pid=0
               endif
               Pid=Pid+1
               gPPI(Pid,Fid) = tmp_single_gPPI
            endif
         enddo
         close(190705)

         ! This is required to calculate Do_FieldLineAverage_by_using_gPPI later
         call Cal_FLA_Denominator_using_gPPI
         call Get_FLA_Te0_Ti0_n0

         return
      endif
   endif

   ! If there is no gPPI data, then calculate the following routines
   call cpu_time(t_start)

   num_of_target_zeta = 1*ntoroidal
   allocate(target_zeta(num_of_target_zeta))
   do k=1,num_of_target_zeta
      target_zeta(k) = 2.0*pi*k/real(num_of_target_zeta)
   enddo
   maximum_tor_rotation = gPPI_max_tor_num
   MaxNum_of_gPPI_for_each_FieldLine = 2*maximum_tor_rotation*(num_of_target_zeta/ntoroidal)

   specific_time = 0.0 ! time is constant because we want the field Poincare Plot at specific time t=t0
   dgf = 2.0*pi/500.0 ! toroidal step for field-line-following integration

   k=1
   do target_mode_n = 1,3
      do target_mode_m = 1,10
         if(target_mode_m <= 3*target_mode_n) then ! search for the safety_q<=3
            ini_a = Find_a_at_specific_safety_q(real(target_mode_m)/real(target_mode_n))
            if(ini_a>0) then
               target_safety_q(k) = real(target_mode_m)/real(target_mode_n)
               k=k+1
            endif
         endif
      enddo
   enddo
   !call flush(stdout)
   call find_unique_sorted_array_and_its_num_of_elements(target_safety_q, 1, k-1, num_of_target_safety_q)
   if(myrank_toroidal==0) then
      write(stdout,*) 'num_of_target_safety_q = ', num_of_target_safety_q
      write(stdout,*) 'target q vals = ',target_safety_q(1:num_of_target_safety_q)
      call flush(stdout)
   endif

   sampling_num_at_each_rational_q = 50/(ntoroidal)


   Fid=1
   do i = 1,num_of_target_safety_q
      ini_a = Find_a_at_specific_safety_q(target_safety_q(i))
      if(ini_a>0) then
         Fid = Fid + 3 ! 3 field lines will be sampled at each target_saftey_q and their neighborhoods.
      endif
   enddo

   ! tmp_num_of_SampledFieldLines = 50*20/(ntoroidal) + (sampling_num_at_each_rational_q*(Fid-1))
   tmp_num_of_SampledFieldLines = 50*60/(ntoroidal) + 3*num_of_target_safety_q*sampling_num_at_each_rational_q

   allocate(tmp_gPPI(MaxNum_of_gPPI_for_each_FieldLine, tmp_num_of_SampledFieldLines, num_of_target_zeta))
   allocate(tmp_num_of_gPPI_for_each_FieldLine(tmp_num_of_SampledFieldLines,num_of_target_zeta))
   allocate(ini_a_arr(tmp_num_of_SampledFieldLines), ini_gq_arr(tmp_num_of_SampledFieldLines), ini_gf_arr(tmp_num_of_SampledFieldLines))



   !Find ini_a_arr and ini_gq_arr at rational q surface and sobol_rand position
   Fid = 1
   do i = 1, num_of_target_safety_q
      ini_a = Find_a_at_specific_safety_q(target_safety_q(i))
      if(ini_a>0) then
         ! At exact target_safety_q
         ini_a_arr(Fid:Fid+(sampling_num_at_each_rational_q-1)) = ini_a
         ini_gf_arr(Fid:Fid+(sampling_num_at_each_rational_q-1)) = 1e-8
         do j=0,sampling_num_at_each_rational_q-1
            ini_gq_arr(Fid+j) = (2.0*pi/ntoroidal)*(myrank_toroidal + j/real(sampling_num_at_each_rational_q))
         enddo
         Fid = Fid + sampling_num_at_each_rational_q

         ! Sligtly left from target_safety_q
         ini_a_arr(Fid:Fid+(sampling_num_at_each_rational_q-1)) = ini_a - 0.01
         ini_gf_arr(Fid:Fid+(sampling_num_at_each_rational_q-1)) = 1e-8
         do j=0,sampling_num_at_each_rational_q-1
            ! ini_gq_arr(Fid+j) = j*(2.0*pi/real(sampling_num_at_each_rational_q))
            ini_gq_arr(Fid+j) = (2.0*pi/ntoroidal)*(myrank_toroidal + j/real(sampling_num_at_each_rational_q))
         enddo
         Fid = Fid + sampling_num_at_each_rational_q

         ! Slightly right from target_safety_q
         ini_a_arr(Fid:Fid+(sampling_num_at_each_rational_q-1)) = ini_a + 0.01
         ini_gf_arr(Fid:Fid+(sampling_num_at_each_rational_q-1)) = 1e-8
         do j=0,sampling_num_at_each_rational_q-1
            ini_gq_arr(Fid+j) = (2.0*pi/ntoroidal)*(myrank_toroidal + j/real(sampling_num_at_each_rational_q))
         enddo
         Fid = Fid + sampling_num_at_each_rational_q
      endif
   enddo


   sobol_seed = 1+ myrank_toroidal*(tmp_num_of_SampledFieldLines-3*num_of_target_safety_q*sampling_num_at_each_rational_q)
   do while(Fid<=tmp_num_of_SampledFieldLines)
      call i8_sobol(sobol_dim_num, sobol_seed, sobol_rand)
      ini_a_arr(Fid) = a0 + (a1-a0)*sobol_rand(1)
      ini_gq_arr(Fid) = 2.0*pi*sobol_rand(2)
      ini_gf_arr(Fid) = 1e-8
      Fid = Fid+1
   enddo

   call quicksort_real_together(ini_a_arr, ini_gq_arr, 1,tmp_num_of_SampledFieldLines)

   ! First Resonant q and a (only the First element of EMP)
   resonant_a = Find_a_at_specific_safety_q(real(EMP(1)%gam)/real(EMP(1)%gan))
   write(stdout,*) '1st resonant_q=',real(EMP(1)%gam)/real(EMP(1)%gan), 'resonant_a=',resonant_a
   call flush(stdout)

   tmp_num_of_gPPI_for_each_FieldLine(:,:) = 0
   do Fid=1,tmp_num_of_SampledFieldLines
      tmp_gPPI(:,Fid,:)%Fid = Fid + myrank_toroidal*tmp_num_of_SampledFieldLines
   enddo

   do Fid = 1, tmp_num_of_SampledFieldLines

      ini_a = ini_a_arr(Fid)
      ini_gq = ini_gq_arr(Fid)
      ini_gf = ini_gf_arr(Fid)

      aaa=ini_a
      the=ini_gq


      do tor_dir=-1,1,2  ! Forward integration & Backward Integration

         dgf = (tor_dir)*2.0*pi/500.0

         now_a = ini_a
         now_gq = ini_gq
         now_gf = ini_gf

         prev_a = now_a
         prev_gq = now_gq
         prev_gf = now_gf

         do while(.true.)
            record_flag=.false.
            !if( modulo(now_gf,2*pi)>zetamin-abs(dgf) .and.  modulo(now_gf,2*pi)<=zetamax+abs(dgf)) then
            do k=1,num_of_target_zeta
               Pid = tmp_num_of_gPPI_for_each_FieldLine(Fid,k)

               if(floor((now_gf-target_zeta(k))/(2.0*pi)) /= floor((prev_gf-target_zeta(k))/(2.0*pi))) then
                  record_flag=.true.
                  tzid = k
                  tor_id = 1 +(k-1)/(num_of_target_zeta/ntoroidal)
                  exit
               endif
            enddo
            !endif

            if(record_flag) then

               tmp_num_of_gPPI_for_each_FieldLine(Fid,tor_id) = tmp_num_of_gPPI_for_each_FieldLine(Fid,tor_id) +1
               Pid = tmp_num_of_gPPI_for_each_FieldLine(Fid,tor_id)

               if(Pid > MaxNum_of_gPPI_for_each_FieldLine) then
                  tmp_num_of_gPPI_for_each_FieldLine(Fid,tor_id) = tmp_num_of_gPPI_for_each_FieldLine(Fid,tor_id) -1
                  exit
               endif

               intp_weight = abs(min(modulo((now_gf-target_zeta(tzid)), 2.0*pi), 2.0*pi - modulo((now_gf-target_zeta(tzid)), 2.0*pi))/dgf) ! interpolation weight for previous point

               tmp_gPPI(Pid,Fid,tor_id)%a  = intp_weight*prev_a + (1.0-intp_weight)*now_a
               tmp_gPPI(Pid,Fid,tor_id)%gq = intp_weight*prev_gq + (1.0-intp_weight)*now_gq
               tmp_gPPI(Pid,Fid,tor_id)%gf = intp_weight*prev_gf + (1.0-intp_weight)*now_gf
               tmp_gPPI(Pid,Fid,tor_id)%target_normalized_gf = target_zeta(tzid)
               tmp_gPPI(Pid,Fid,tor_id)%ini_FLabel = Cal_FLabel_Function(tmp_gPPI(Pid,Fid,tor_id)%a, tmp_gPPI(Pid,Fid,tor_id)%gq, tmp_gPPI(Pid,Fid,tor_id)%gf)

               call Convert_a_to_s(tmp_gPPI(Pid,Fid,tor_id)%a, tmp_gPPI(Pid,Fid,tor_id)%s)


               aaa=tmp_gPPI(Pid,Fid,tor_id)%a
               the=modulo(tmp_gPPI(Pid,Fid,tor_id)%gq, 2.0*pi)
               call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
               tmp_gPPI(Pid,Fid,tor_id)%inv_B = 1.0/Bj; ! inv_B
               tmp_gPPI(Pid,Fid,tor_id)%R = rj
               tmp_gPPI(Pid,Fid,tor_id)%Z = zj

               call esiget_Intp_B_vec_contravariant(aaa, the, modulo(tmp_gPPI(Pid,Fid,tor_id)%gf,2.0*pi), specific_time, B_vec_cntv)
               tmp_gPPI(Pid,Fid,tor_id)%dl = abs(dgf)*sqrt( ((B_vec_cntv(1)/B_vec_cntv(3))**2)*(raj**2+zaj**2)  &
                                                         +((B_vec_cntv(2)/B_vec_cntv(3))**2)*(rqj**2+zqj**2)  &
                                                         +(B_vec_cntv(1)*B_vec_cntv(2)/(B_vec_cntv(3)**2))*(raj*rqj+zaj*zqj) &
                                                         + rj**2  )
            endif

            ! if(abs(now_gf-ini_gf)>maximum_tor_rotation*(2*pi) .or.  Pid>=MaxNum_of_gPPI_for_each_FieldLine) then
            if(abs(now_gf-ini_gf)>maximum_tor_rotation*(2*pi)) then
               exit
            endif

            prev_a = now_a
            prev_gq = now_gq
            prev_gf = now_gf

            ! Runge-Kutta 4th order
            call esiget_Intp_B_vec_contravariant(now_a, modulo(now_gq,2.0*pi), modulo(now_gf,2.0*pi), specific_time, B_vec_cntv)
            a_RK(1)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
            gq_RK(1)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
            call esiget_Intp_B_vec_contravariant(now_a+0.5*a_RK(1), modulo(now_gq+0.5*gq_RK(1),2.0*pi), modulo(now_gf+0.5*dgf,2.0*pi), specific_time, B_vec_cntv)
            a_RK(2)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
            gq_RK(2)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
            call esiget_Intp_B_vec_contravariant(now_a+0.5*a_RK(2), modulo(now_gq+0.5*gq_RK(2),2.0*pi), modulo(now_gf+0.5*dgf,2.0*pi), specific_time, B_vec_cntv)
            a_RK(3)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
            gq_RK(3)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
            call esiget_Intp_B_vec_contravariant(now_a+a_RK(3), modulo(now_gq+gq_RK(3),2.0*pi), modulo(now_gf+dgf,2.0*pi), specific_time, B_vec_cntv)
            a_RK(4)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
            gq_RK(4)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv

            now_a = now_a + (a_RK(1)+2.0*a_RK(2)+2.0*a_RK(3)+a_RK(4))/6.0
            now_gq = now_gq + (gq_RK(1)+2.0*gq_RK(2)+2.0*gq_RK(3)+gq_RK(4))/6.0
            now_gf = now_gf + dgf

            if (now_a<a0 .or. now_a>a1) then
               exit
            endif
         enddo
      enddo
   enddo


   num_of_SampledFieldLines = num_of_target_zeta*tmp_num_of_SampledFieldLines


   allocate(gPPI(MaxNum_of_gPPI_for_each_FieldLine, num_of_SampledFieldLines))
   allocate(num_of_gPPI_for_each_FieldLine(num_of_SampledFieldLines))

   ! write(stdout,*) "storage_size(tmp_gPPI & gPPI)=", storage_size(tmp_gPPI), storage_size(gPPI),storage_size(k), storage_size(now_a),storage_size(MPI_BYTE), storage_size(EMP)
   ! call flush(stdout)

   ! length = MaxNum_of_gPPI_for_each_FieldLine*tmp_num_of_SampledFieldLines
   length = MaxNum_of_gPPI_for_each_FieldLine*tmp_num_of_SampledFieldLines*storage_size(tmp_gPPI)/storage_size(now_a)
   call MPI_Alltoall(tmp_gPPI, length, MPI_Rsize, gPPI, length, MPI_Rsize, toroidal_comm, ierr)

   length = tmp_num_of_SampledFieldLines
   call MPI_Alltoall(tmp_num_of_gPPI_for_each_FieldLine, length, MPI_INTEGER, num_of_gPPI_for_each_FieldLine, length, MPI_INTEGER, toroidal_comm, ierr);

   call Cal_FLA_Denominator_using_gPPI
   call Get_FLA_Te0_Ti0_n0


   call cpu_time(t_finish)
   write(stdout,'(a,I5,a,ES11.3,a)') "I'm=", mype, ", gPPI Construction Finished! delta_time = ", t_finish-t_start, " seconds."
   flush(stdout)

   call system("mkdir -p " // 'VacFieldAnalysis')
   write(chdum,'("VacFieldAnalysis/gPPI_",i5.5,".dat")') myrank_toroidal
   open(181202,file=chdum,status='replace',action='write')
   do Fid = 1, num_of_SampledFieldLines
      do Pid =1,num_of_gPPI_for_each_FieldLine(Fid)
         write(181202, "(I8,1x,11(ES16.8,1x))") gPPI(Pid,Fid)
      enddo
   enddo
   flush(181202)
   close(181202)


   open(190704,file="VacFieldAnalysis/gPPI_Information.dat",status='replace',form='formatted',action='write')
   write(190704,"(3(I8,1x))") gPPI_max_tor_num, MaxNum_of_gPPI_for_each_FieldLine, num_of_SampledFieldLines
   close(190704)

   deallocate(ini_a_arr, ini_gq_arr, ini_gf_arr)
   deallocate(tmp_gPPI,tmp_num_of_gPPI_for_each_FieldLine)
   deallocate(target_zeta)

end subroutine Construct_gPPI_for_FieldLineAverage_Parallelized


function Find_a_at_specific_safety_q(target_safety_q) result(out_a)

   use global_parameters
   use field_array
   use precision

   implicit none

   integer aid
   real(wp) prev_q, this_q, prev_a, this_a
   real(wp) target_safety_q, out_a


   do aid=1,mpsi
      prev_q = q_grid(aid-1)
      this_q = q_grid(aid)

      if ( (prev_q-target_safety_q)*(this_q-target_safety_q) < 0 ) then
         ! Take linear interpolation to find solution
         prev_a = psi_grid(aid-1)
         this_a = psi_grid(aid)
         out_a = prev_a + ((this_a-prev_a)/(this_q-prev_q))*(target_safety_q-prev_q)
         return
      endif
   enddo

   !write(*,*) ("We cannot find q=",target_safety_q)
   out_a = -100.0

end function

function Find_safety_q_at_specific_a(target_a) result(out_safety_q)

   use global_parameters
   use field_array
   use precision

   implicit none

   integer aid
   real(wp) prev_q, this_q, prev_a, this_a
   real(wp) target_a, out_safety_q


   do aid=1,mpsi

      prev_a = psi_grid(aid-1)
      this_a = psi_grid(aid)

      if((target_a-prev_a)*(target_a-this_a) <0) then
         prev_q = q_grid(aid-1)
         this_q = q_grid(aid)

         out_safety_q  = prev_q + (target_a-prev_a)/(this_a-prev_a)*(this_q-prev_q)
         return
      endif
   enddo

   write(*,*) ("We cannot find safety_q at a=",target_a)
   out_safety_q = -100.0
end function


subroutine Cal_FLA_Denominator_using_gPPI
   use mpi_defs
   use precision
   use global_parameters
   use particle_decomp
   use FieldLineFollowing
   use field_array

   implicit none

   real(wp) this_a, this_gq, this_s, this_gf, sdum, gqdum,inv_ds, inv_dgq(0:mpsi), inv_gf
   real(wp) ws0, ws1, wgq00, wgq01, wgq10, wgq11, wgf0, wgf1
   real(wp) tmp_weight
   integer Fid, Pid ! (Fid=Field ID), (Pid=gPPI ID)
   integer inner_sid, outter_sid, inner_gqid, outter_gqid, inner_ij, outter_ij
   integer i,ii,jt, ierror

   real(wp) sendl(mgrid),recvr(mgrid)
   integer icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)

   if (num_of_SampledFieldLines ==0) then
      write(*,*) "Error! Please construct gPPI first before using FieldLineAverage!!"
      STOP
   endif

   ! if(myrank_toroidal==0) then
   !   write(*,*) "Do_FieldLineAverage_by_using_gPPI"
   ! endif
   inv_ds = real(mpsi)/(s1-s0)
   inv_dgq = 1.0/deltat
   inv_gf = 1.0/deltaz


   FLA_Denominator = 0.0 ! Initial zero value can result NaN (FLA_ResultData/FLA_Denominator) in the case of no touching grid nodes
   do Fid=1,num_of_SampledFieldLines
      do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

         this_s = gPPI(Pid,Fid)%s
         this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
         !this_gf = modulo(gPPI(Pid,Fid)%gf, 2.0*pi)
         this_gf = gPPI(Pid,Fid)%target_normalized_gf

         sdum = (this_s-s0)*inv_ds
         inner_sid = int(sdum)
         ws1  = sdum - inner_sid
         ws0  = 1.0 - ws1
         gqdum = inv_dgq(inner_sid)*modulo(this_gq - qtinv(inner_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
         inner_gqid = int(gqdum)
         wgq01 = gqdum - inner_gqid
         wgq00 = 1.0 - wgq01
         inner_ij = igrid(inner_sid) + inner_gqid

         outter_sid = inner_sid+1
         gqdum = inv_dgq(outter_sid)*modulo(this_gq - qtinv(outter_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
         outter_gqid = int(gqdum)
         wgq11 = gqdum - outter_gqid
         wgq10 = 1.0 - wgq11
         outter_ij = igrid(outter_sid) + outter_gqid

         !wgf1 = min(modulo(this_gf-zetamin,2.0*pi),2.0*pi-modulo(this_gf-zetamin,2.0*pi))*inv_gf
         wgf1 = (this_gf-zetamin)*inv_gf
         wgf0 = 1.0 - wgf1


         tmp_weight = gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
         !FLA_Denominator(0,inner_ij)    = FLA_Denominator(0,inner_ij)    +wgf0*ws0*wgq00*gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
         !FLA_Denominator(0,inner_ij)    = FLA_Denominator(0,inner_ij)    +wgf0*ws0*wgq00*tmp_weight
         FLA_Denominator(0,inner_ij)    = FLA_Denominator(0,inner_ij)    +wgf0*ws0*wgq00
         FLA_Denominator(0,inner_ij+1)  = FLA_Denominator(0,inner_ij+1)  +wgf0*ws0*wgq01
         FLA_Denominator(0,outter_ij)   = FLA_Denominator(0,outter_ij)   +wgf0*ws1*wgq10
         FLA_Denominator(0,outter_ij+1) = FLA_Denominator(0,outter_ij+1) +wgf0*ws1*wgq11
         FLA_Denominator(1,inner_ij)    = FLA_Denominator(1,inner_ij)    +wgf1*ws0*wgq00
         FLA_Denominator(1,inner_ij+1)  = FLA_Denominator(1,inner_ij+1)  +wgf1*ws0*wgq01
         FLA_Denominator(1,outter_ij)   = FLA_Denominator(1,outter_ij)   +wgf1*ws1*wgq10
         FLA_Denominator(1,outter_ij+1) = FLA_Denominator(1,outter_ij+1) +wgf1*ws1*wgq11

      enddo
   enddo


   ! This part is the required post-processing and/or communication to get correct FLA_Denominator
   ! poloidal end cell; need valid data for densityi(:,igrid(i))
   do i=0,mpsi
      FLA_Denominator(:,igrid(i)+mtheta(i))=FLA_Denominator(:,igrid(i)+mtheta(i))+FLA_Denominator(:,igrid(i))
      FLA_Denominator(:,igrid(i))=FLA_Denominator(:,igrid(i)+mtheta(i))
   enddo
   ! Communicate FLA_Denominator between nodes
   ! toroidal end cell
   sendl=FLA_Denominator(0,:)
   recvr=0.0_wp
   icount=mgrid
   idest=left_pe
   isource=right_pe
   isendtag=myrank_toroidal
   irecvtag=isource
   ! send densityi to left and receive from right
   call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,recvr,icount,&
      mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      ! B.C. at zeta=2*pi is shifted
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         FLA_Denominator(1,ii+1:ii+jt)=FLA_Denominator(1,ii+1:ii+jt)+ cshift(recvr(ii+1:ii+jt),itran(i))
         FLA_Denominator(1,igrid(i))=FLA_Denominator(1,igrid(i)+mtheta(i))
      enddo
   else
      ! B.C. at zeta<2*pi is continuous
      FLA_Denominator(1,:)=FLA_Denominator(1,:)+recvr
   endif

end subroutine Cal_FLA_Denominator_using_gPPI

subroutine Get_FLA_Te0_Ti0_n0

   use mpi_defs
   use precision
   use global_parameters
   use particle_decomp
   use FieldLineFollowing
   use field_array

   implicit none

   real(wp) Te0_Data(0:1,mgrid), Ti0_Data(0:1,mgrid), n0_Data(0:1,mgrid), radial_psi(0:1, mgrid)
   integer i,j, ij

   do i=0,mpsi
      Te0_Data(:,igrid(i):igrid(i)+mtheta(i)) = Te_grid(i)
      Ti0_Data(:,igrid(i):igrid(i)+mtheta(i)) = vth_grid(i)**2
      n0_Data(:,igrid(i):igrid(i)+mtheta(i)) = den_grid(i)
      radial_psi(:,igrid(i):igrid(i)+mtheta(i)) = psi_grid(i)
   enddo

   call Do_FieldLineAverage_by_using_gPPI(Te0_Data, FLA_Te0)
   call Do_FieldLineAverage_by_using_gPPI(Ti0_Data, FLA_Ti0)
   call Do_FieldLineAverage_by_using_gPPI(n0_Data, FLA_n0)
   call Do_FieldLineAverage_by_using_gPPI(radial_psi, FLA_psi0)

   call Do_FieldLineAverage_by_using_gPPI(FLabel_Function, FLA_FLabel)

end subroutine Get_FLA_Te0_Ti0_n0



subroutine Do_FieldLineAverage_by_using_gPPI(TargetData, FLA_ResultData)
   use mpi_defs
   use precision
   use global_parameters
   use particle_decomp
   use FieldLineFollowing
   use field_array

   implicit none

   real(wp) TargetData(0:1,mgrid), FLA_ResultData(0:1,mgrid)!, FLA_Denominator(0:1,mgrid)

   real(wp) each_FLA_val(num_of_SampledFieldLines), each_FLA_denominator(num_of_SampledFieldLines), intpl_val ! FLA=Field Line Average
   real(wp) this_a, this_gq, this_s, this_gf, sdum, gqdum,inv_ds, inv_dgq(0:mpsi), inv_gf
   real(wp) ws0, ws1, wgq00, wgq01, wgq10, wgq11, wgf0, wgf1
   integer Fid, Pid ! (Fid=Field ID), (Pid=gPPI ID)
   integer inner_sid, outter_sid, inner_gqid, outter_gqid, inner_ij, outter_ij
   integer i,ii,jt, ierror

   real(wp) tmp_weight
   real(wp) sendl(mgrid),recvr(mgrid)
   integer icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)


   if (num_of_SampledFieldLines ==0) then
      write(*,*) "Error! Please construct gPPI first before using FieldLineAverage!!"
      STOP
   endif

   if(myrank_toroidal==0) then
      write(*,*) "Do_FieldLineAverage_by_using_gPPI"
   endif
   inv_ds = real(mpsi)/(s1-s0)
   inv_dgq = 1.0/deltat
   inv_gf = 1.0/deltaz

   ! Get FieldLineAverages of Target_data
   do Fid=1,num_of_SampledFieldLines
      each_FLA_val(Fid)=0
      each_FLA_denominator(Fid)=0
      do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

         !this_a = gPPI(Pid,Fid)%a
         this_s = gPPI(Pid,Fid)%s
         this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
         !this_gf = modulo(gPPI(Pid,Fid)%gf, 2.0*pi)
         this_gf = gPPI(Pid,Fid)%target_normalized_gf
         ! Warning!!!
         !Interpolation 하려면 a값으로부터 s grid 값 알아야함..
         !일단 아래 같은 함수 통해서 알 수 있다고 가정하고 코딩 진행..
         !call Convert_a_to_s(this_a, this_s)

         !sid = max(0, min(mpsi, int((this_s-s0)*inv_ds)))
         !gqid = max(0, min(mtheta(sid), int(gqdum*inv_dgq(sid))))
         sdum = (this_s-s0)*inv_ds
         inner_sid = int(sdum)
         ws1  = sdum - inner_sid
         ws0  = 1.0 - ws1
         gqdum = inv_dgq(inner_sid)*modulo(this_gq - qtinv(inner_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
         inner_gqid = int(gqdum)
         wgq01 = gqdum - inner_gqid
         wgq00 = 1.0 - wgq01
         inner_ij = igrid(inner_sid) + inner_gqid

         outter_sid = inner_sid+1
         gqdum = inv_dgq(outter_sid)*modulo(this_gq - qtinv(outter_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
         outter_gqid = int(gqdum)
         wgq11 = gqdum - outter_gqid
         wgq10 = 1.0 - wgq11
         outter_ij = igrid(outter_sid) + outter_gqid


         !wgf1 = min(modulo(this_gf-zetamin,2.0*pi),2.0*pi-modulo(this_gf-zetamin,2.0*pi))*inv_gf
         wgf1 = (this_gf-zetamin)*inv_gf
         wgf0 = 1.0 - wgf1

         intpl_val = wgf0*(ws0*(wgq00*TargetData(0,inner_ij) + wgq01*TargetData(0,inner_ij+1)) + ws1*(wgq10*TargetData(0,outter_ij)+wgq11*TargetData(0,outter_ij+1))) &
            +wgf1*(ws0*(wgq00*TargetData(1,inner_ij) + wgq01*TargetData(1,inner_ij+1)) + ws1*(wgq10*TargetData(1,outter_ij)+wgq11*TargetData(1,outter_ij+1)))

         each_FLA_val(Fid) = each_FLA_val(Fid) + intpl_val*gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
         each_FLA_denominator(Fid) = each_FLA_denominator(Fid) + gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B

      enddo
   enddo

   ! Do Allreduce across toroidal_comm to get a sum of FLA
   call MPI_ALLREDUCE(each_FLA_val,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
   each_FLA_val = recvr(1:num_of_SampledFieldLines)
   call MPI_ALLREDUCE(each_FLA_denominator,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
   each_FLA_denominator = recvr(1:num_of_SampledFieldLines)

   each_FLA_val =  each_FLA_val/ each_FLA_denominator ! Take average




   ! Scatter Averaged data to grid nodes
   FLA_ResultData = 0.0
   do Fid=1,num_of_SampledFieldLines
      do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

         !this_a = gPPI(Pid,Fid)%a
         this_s = gPPI(Pid,Fid)%s
         this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
         !this_gf = modulo(gPPI(Pid,Fid)%gf, 2.0*pi)
         this_gf = gPPI(Pid,Fid)%target_normalized_gf

         ! Warning!!!
         !Interpolation 하려면 a값으로부터 s grid 값 알아야함..
         !일단 아래 같은 함수 통해서 알 수 있다고 가정하고 코딩 진행..
         !call Convert_a_to_s(this_a, this_s)

         sdum = (this_s-s0)*inv_ds
         inner_sid = int(sdum)
         ws1  = sdum - inner_sid
         ws0  = 1.0 - ws1
         gqdum = inv_dgq(inner_sid)*modulo(this_gq - qtinv(inner_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
         inner_gqid = int(gqdum)
         wgq01 = gqdum - inner_gqid
         wgq00 = 1.0 - wgq01
         inner_ij = igrid(inner_sid) + inner_gqid

         outter_sid = inner_sid+1
         gqdum = inv_dgq(outter_sid)*modulo(this_gq - qtinv(outter_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
         outter_gqid = int(gqdum)
         wgq11 = gqdum - outter_gqid
         wgq10 = 1.0 - wgq11
         outter_ij = igrid(outter_sid) + outter_gqid


         !wgf1 = min(modulo(this_gf-zetamin,2.0*pi),2.0*pi-modulo(this_gf-zetamin,2.0*pi))*inv_gf
         wgf1 = (this_gf-zetamin)*inv_gf
         wgf0 = 1.0 - wgf1

         !tmp_weight = gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
         !      FLA_ResultData(0,inner_ij)    = FLA_ResultData(0,inner_ij)    +wgf0*ws0*wgq00*each_FLA_val(Fid)*gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
         !FLA_ResultData(0,inner_ij)    = FLA_ResultData(0,inner_ij)    +wgf0*ws0*wgq00*each_FLA_val(Fid)*tmp_weight
         FLA_ResultData(0,inner_ij)    = FLA_ResultData(0,inner_ij)    +wgf0*ws0*wgq00*each_FLA_val(Fid)
         FLA_ResultData(0,inner_ij+1)  = FLA_ResultData(0,inner_ij+1)  +wgf0*ws0*wgq01*each_FLA_val(Fid)
         FLA_ResultData(0,outter_ij)   = FLA_ResultData(0,outter_ij)   +wgf0*ws1*wgq10*each_FLA_val(Fid)
         FLA_ResultData(0,outter_ij+1) = FLA_ResultData(0,outter_ij+1) +wgf0*ws1*wgq11*each_FLA_val(Fid)
         FLA_ResultData(1,inner_ij)    = FLA_ResultData(1,inner_ij)    +wgf1*ws0*wgq00*each_FLA_val(Fid)
         FLA_ResultData(1,inner_ij+1)  = FLA_ResultData(1,inner_ij+1)  +wgf1*ws0*wgq01*each_FLA_val(Fid)
         FLA_ResultData(1,outter_ij)   = FLA_ResultData(1,outter_ij)   +wgf1*ws1*wgq10*each_FLA_val(Fid)
         FLA_ResultData(1,outter_ij+1) = FLA_ResultData(1,outter_ij+1) +wgf1*ws1*wgq11*each_FLA_val(Fid)
      enddo
   enddo

   ! This part is the required post-processing and/or communication to get correct FLA_ResultData
   ! poloidal end cell; need valid data for densityi(:,igrid(i))
   do i=0,mpsi
      FLA_ResultData(:,igrid(i)+mtheta(i))=FLA_ResultData(:,igrid(i)+mtheta(i))+FLA_ResultData(:,igrid(i))
      FLA_ResultData(:,igrid(i))=FLA_ResultData(:,igrid(i)+mtheta(i))
   enddo
   ! Communicate FLA_ResultData between nodes
   ! toroidal end cell
   sendl=FLA_ResultData(0,:)
   recvr=0.0_wp
   icount=mgrid
   idest=left_pe
   isource=right_pe
   isendtag=myrank_toroidal
   irecvtag=isource
   ! send densityi to left and receive from right
   call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,recvr,icount,&
      mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      ! B.C. at zeta=2*pi is shifted
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         FLA_ResultData(1,ii+1:ii+jt)=FLA_ResultData(1,ii+1:ii+jt)+ cshift(recvr(ii+1:ii+jt),itran(i))
         FLA_ResultData(1,igrid(i))=FLA_ResultData(1,igrid(i)+mtheta(i))
      enddo
   else
      ! B.C. at zeta<2*pi is continuous
      FLA_ResultData(1,:)=FLA_ResultData(1,:)+recvr
   endif


   ! Last part is to divide FLA_ResultData by FLA_Denominator to get correct line-averaged FLA_ResultData
   FLA_ResultData = FLA_ResultData/FLA_Denominator

end subroutine Do_FieldLineAverage_by_using_gPPI

subroutine Do_FieldLineAverage_from_oriGrid_to_DG(oriG_TargetData, DG_FLA_ResultData)
   use mpi_defs
   use precision
   use global_parameters
   use particle_decomp
   use FieldLineFollowing
   use field_array
   use Diagnostic_Grid
   implicit none

   real(wp), intent(in) :: oriG_TargetData(0:1,mgrid)
   real(wp), intent(out):: DG_FLA_ResultData(DG%dim_gq, DG%dim_a,0:1)
   real(wp) tmp_DG_Data(DG%dim_gq, DG%dim_a)

   real(wp) each_FLA_val(num_of_SampledFieldLines), each_FLA_denominator(num_of_SampledFieldLines), intpl_val ! FLA=Field Line Average
   real(wp) this_gq, this_s, this_gf, sdum, gqdum,inv_ds, inv_dgq(0:mpsi), inv_gf
   real(wp) ws0, ws1, wgq00, wgq01, wgq10, wgq11, wgf0, wgf1
   integer Fid, Pid ! (Fid=Field ID), (Pid=gPPI ID)
   integer inner_sid, outter_sid, inner_gqid, outter_gqid, inner_ij, outter_ij
   integer ierror

   real(wp) recvr(mgrid)
   integer istatus(MPI_STATUS_SIZE)

   real(wp) target_gq, target_zeta
   real(wp) wgq0,wgq1,wgf,ratio
   integer ori_sid, j,k, sid, gqid, qq, ss, length
   real(wp) ff, this_safety_q



   if (num_of_SampledFieldLines ==0) then
      write(*,*) "Error! Please construct gPPI first before using FieldLineAverage!!"
      STOP
   endif

   if(myrank_toroidal==0) then
      write(*,*) "Do_FieldLineAverage_by_using_gPPI"
   endif
   inv_ds = real(mpsi)/(s1-s0)
   inv_dgq = 1.0/deltat
   inv_gf = 1.0/deltaz

   ! Get FieldLineAverages of Target_data
   do Fid=1,num_of_SampledFieldLines
      each_FLA_val(Fid)=0
      each_FLA_denominator(Fid)=0
      do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

         this_s = gPPI(Pid,Fid)%s
         this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
         this_gf = gPPI(Pid,Fid)%target_normalized_gf

         sdum = (this_s-s0)*inv_ds
         inner_sid = int(sdum)
         ws1  = sdum - inner_sid
         ws0  = 1.0 - ws1
         gqdum = inv_dgq(inner_sid)*modulo(this_gq - qtinv(inner_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
         inner_gqid = int(gqdum)
         wgq01 = gqdum - inner_gqid
         wgq00 = 1.0 - wgq01
         inner_ij = igrid(inner_sid) + inner_gqid

         outter_sid = inner_sid+1
         gqdum = inv_dgq(outter_sid)*modulo(this_gq - qtinv(outter_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
         outter_gqid = int(gqdum)
         wgq11 = gqdum - outter_gqid
         wgq10 = 1.0 - wgq11
         outter_ij = igrid(outter_sid) + outter_gqid


         wgf1 = (this_gf-zetamin)*inv_gf
         wgf0 = 1.0 - wgf1

         intpl_val = wgf0*(ws0*(wgq00*oriG_TargetData(0,inner_ij) + wgq01*oriG_TargetData(0,inner_ij+1)) + ws1*(wgq10*oriG_TargetData(0,outter_ij)+wgq11*oriG_TargetData(0,outter_ij+1))) &
            +wgf1*(ws0*(wgq00*oriG_TargetData(1,inner_ij) + wgq01*oriG_TargetData(1,inner_ij+1)) + ws1*(wgq10*oriG_TargetData(1,outter_ij)+wgq11*oriG_TargetData(1,outter_ij+1)))

         each_FLA_val(Fid) = each_FLA_val(Fid) + intpl_val*gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
         each_FLA_denominator(Fid) = each_FLA_denominator(Fid) + gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B

      enddo
   enddo

   ! Do Allreduce across toroidal_comm to get a sum of FLA
   call MPI_ALLREDUCE(each_FLA_val,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
   each_FLA_val = recvr(1:num_of_SampledFieldLines)
   call MPI_ALLREDUCE(each_FLA_denominator,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
   each_FLA_denominator = recvr(1:num_of_SampledFieldLines)

   each_FLA_val =  each_FLA_val/ each_FLA_denominator ! Take average



   ! Scatter Averaged data to Diagnostic Grid (DG) nodes
   DG_FLA_ResultData = 0.0
   do Fid=1,num_of_SampledFieldLines
      do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

         this_s = gPPI(Pid,Fid)%s
         this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
         this_gf = gPPI(Pid,Fid)%target_normalized_gf

         ori_sid = (this_s-s0)*inv_ds ! s index of original full grid
         ff  = (this_s-s_grid(ori_sid))*inv_ds
         this_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)

         sid = max(0,min(int((this_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
         ws1 = (this_s-DG%s_1D(sid))*DG%inv_ds
         ws0 = 1.0-ws1

         do k=0,1
            if(k==0) target_zeta =zetamin
            if(k==1) target_zeta =zetamax

            wgf = (1.0-abs(this_gf-target_zeta)/deltaz)

            target_gq = modulo(this_gq + (target_zeta-this_gf)/this_safety_q, 2.0*pi)
            gqid = (target_gq)*DG%inv_dgq +1
            wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
            wgq0 = 1.0-wgq1

            do j=1,4
               select case(j)
                  case(1)
                     qq=gqid; ss=sid
                     ratio=ws0*wgq0*wgf
                  case(2)
                     qq=gqid+1; ss=sid
                     ratio=ws0*wgq1*wgf
                  case(3)
                     qq=gqid; ss=sid+1
                     ratio=ws1*wgq0*wgf
                  case(4)
                     qq=gqid+1; ss=sid+1
                     ratio=ws1*wgq1*wgf
               end select
               DG_FLA_ResultData(qq,ss,k)  = DG_FLA_ResultData(qq,ss,k)  + ratio*each_FLA_val(Fid)
            enddo
         enddo
      enddo
   enddo



   ! Match poloidal end cell
   ! (gq=0 and gq=2pi is the same position physically
   ! although they have different gqid)
   do sid=1,DG%dim_a
      DG_FLA_ResultData(1,sid,:) = DG_FLA_ResultData(1,sid,:) + DG_FLA_ResultData(DG%dim_gq,sid,:)
      DG_FLA_ResultData(DG%dim_gq,sid,:) = DG_FLA_ResultData(1,sid,:)
   enddo

   ! Communicate DG_FLA_ResultData between nodes
   ! toroidal end cell
   length = DG%dim_gq*DG%dim_a
   call MPI_SENDRECV(DG_FLA_ResultData(:,:,0),length,mpi_Rsize,left_pe,1, tmp_DG_Data,length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   DG_FLA_ResultData(:,:,1) = DG_FLA_ResultData(:,:,1) + tmp_DG_Data
   call MPI_SENDRECV(DG_FLA_ResultData(:,:,1),length,mpi_Rsize,right_pe,1, DG_FLA_ResultData(:,:,0),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

   ! Last part is to divide FLA_ResultData by DG%FLA_Denominator to get correct line-averaged FLA_ResultData
   DG_FLA_ResultData = DG_FLA_ResultData/DG%FLA_Denominator

end subroutine Do_FieldLineAverage_from_oriGrid_to_DG

subroutine Do_FieldLineAverage_on_DG_by_using_gPPI(DG_TargetData, DG_FLA_ResultData)
   use mpi_defs
   use precision
   use global_parameters
   use particle_decomp
   use FieldLineFollowing
   use field_array
   use Diagnostic_Grid
   implicit none

   real(wp), intent(in) :: DG_TargetData(DG%dim_gq, DG%dim_a,0:1)
   real(wp), intent(out):: DG_FLA_ResultData(DG%dim_gq, DG%dim_a,0:1)
   real(wp) tmp_DG_Data(DG%dim_gq, DG%dim_a)

   real(wp) each_FLA_val(num_of_SampledFieldLines), each_FLA_denominator(num_of_SampledFieldLines), intpl_val ! FLA=Field Line Average
   real(wp) this_gq, this_s, this_gf, sdum, gqdum,ori_inv_ds
   real(wp) ws0, ws1, wgq00, wgq01, wgq10, wgq11, wgf0, wgf1
   integer Fid, Pid ! (Fid=Field ID), (Pid=gPPI ID)
   integer inner_sid, outter_sid, inner_gqid, outter_gqid, inner_ij, outter_ij
   integer ierror

   real(wp) recvr(mgrid)
   integer istatus(MPI_STATUS_SIZE)

   real(wp) target_gq, target_zeta
   real(wp) wgq0,wgq1,wgf,ratio
   integer ori_sid, j,k, sid, gqid, qq, ss, length
   real(wp) ff, this_safety_q



   if (num_of_SampledFieldLines ==0) then
      write(*,*) "Error! Please construct gPPI first before using FieldLineAverage!!"
      STOP
   endif

   if(myrank_toroidal==0) then
      write(*,*) "Do_FieldLineAverage_by_using_gPPI"
   endif
   ori_inv_ds = real(mpsi)/(s1-s0)



   ! Scatter Averaged data to Diagnostic Grid (DG) nodes
   DG_FLA_ResultData = 0.0
   each_FLA_val = 0.0
   each_FLA_denominator = 0.0
   do Fid=1,num_of_SampledFieldLines
      do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

         this_s = gPPI(Pid,Fid)%s
         this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
         this_gf = gPPI(Pid,Fid)%target_normalized_gf

         ori_sid = (this_s-s0)*ori_inv_ds ! s index of original full grid
         ff  = (this_s-s_grid(ori_sid))*ori_inv_ds
         this_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)

         sid = max(0,min(int((this_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
         ws1 = (this_s-DG%s_1D(sid))*DG%inv_ds
         ws0 = 1.0-ws1

         intpl_val=0
         do k=0,1
            if(k==0) target_zeta =zetamin
            if(k==1) target_zeta =zetamax

            wgf = (1.0-abs(this_gf-target_zeta)/deltaz)

            target_gq = modulo(this_gq + (target_zeta-this_gf)/this_safety_q, 2.0*pi)
            gqid = (target_gq)*DG%inv_dgq +1
            wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
            wgq0 = 1.0-wgq1

            do j=1,4
               select case(j)
                  case(1)
                     qq=gqid; ss=sid
                     ratio=ws0*wgq0*wgf
                  case(2)
                     qq=gqid+1; ss=sid
                     ratio=ws0*wgq1*wgf
                  case(3)
                     qq=gqid; ss=sid+1
                     ratio=ws1*wgq0*wgf
                  case(4)
                     qq=gqid+1; ss=sid+1
                     ratio=ws1*wgq1*wgf
               end select
               intpl_val = intpl_val + ratio*DG_TargetData(qq,ss,k)
            enddo
         enddo
         each_FLA_val(Fid) = each_FLA_val(Fid) + intpl_val*gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
         each_FLA_denominator(Fid) = each_FLA_denominator(Fid) + gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
      enddo
   enddo

   ! Do Allreduce across toroidal_comm to get a sum of FLA
   call MPI_ALLREDUCE(each_FLA_val,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
   each_FLA_val = recvr(1:num_of_SampledFieldLines)
   call MPI_ALLREDUCE(each_FLA_denominator,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
   each_FLA_denominator = recvr(1:num_of_SampledFieldLines)

   each_FLA_val =  each_FLA_val/ each_FLA_denominator ! Take average



   ! Scatter Averaged data to Diagnostic Grid (DG) nodes
   DG_FLA_ResultData = 0.0
   do Fid=1,num_of_SampledFieldLines
      do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

         this_s = gPPI(Pid,Fid)%s
         this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
         this_gf = gPPI(Pid,Fid)%target_normalized_gf

         ori_sid = (this_s-s0)*ori_inv_ds ! s index of original full grid
         ff  = (this_s-s_grid(ori_sid))*ori_inv_ds
         this_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)

         sid = max(0,min(int((this_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
         ws1 = (this_s-DG%s_1D(sid))*DG%inv_ds
         ws0 = 1.0-ws1

         do k=0,1
            if(k==0) target_zeta =zetamin
            if(k==1) target_zeta =zetamax

            wgf = (1.0-abs(this_gf-target_zeta)/deltaz)

            target_gq = modulo(this_gq + (target_zeta-this_gf)/this_safety_q, 2.0*pi)
            gqid = (target_gq)*DG%inv_dgq +1
            wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
            wgq0 = 1.0-wgq1

            do j=1,4
               select case(j)
                  case(1)
                     qq=gqid; ss=sid
                     ratio=ws0*wgq0*wgf
                  case(2)
                     qq=gqid+1; ss=sid
                     ratio=ws0*wgq1*wgf
                  case(3)
                     qq=gqid; ss=sid+1
                     ratio=ws1*wgq0*wgf
                  case(4)
                     qq=gqid+1; ss=sid+1
                     ratio=ws1*wgq1*wgf
               end select
               DG_FLA_ResultData(qq,ss,k)  = DG_FLA_ResultData(qq,ss,k)  + ratio*each_FLA_val(Fid)
            enddo
         enddo
      enddo
   enddo


   ! Match poloidal end cell
   ! (gq=0 and gq=2pi is the same position physically
   ! although they have different gqid)
   do sid=1,DG%dim_a
      DG_FLA_ResultData(1,sid,:) = DG_FLA_ResultData(1,sid,:) + DG_FLA_ResultData(DG%dim_gq,sid,:)
      DG_FLA_ResultData(DG%dim_gq,sid,:) = DG_FLA_ResultData(1,sid,:)
   enddo

   ! Communicate DG_FLA_ResultData between nodes
   ! toroidal end cell
   length = DG%dim_gq*DG%dim_a
   call MPI_SENDRECV(DG_FLA_ResultData(:,:,0),length,mpi_Rsize,left_pe,1, tmp_DG_Data,length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   DG_FLA_ResultData(:,:,1) = DG_FLA_ResultData(:,:,1) + tmp_DG_Data
   call MPI_SENDRECV(DG_FLA_ResultData(:,:,1),length,mpi_Rsize,right_pe,1, DG_FLA_ResultData(:,:,0),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

   ! Last part is to divide FLA_ResultData by DG%FLA_Denominator to get correct line-averaged FLA_ResultData
   DG_FLA_ResultData = DG_FLA_ResultData/DG%FLA_Denominator

end subroutine Do_FieldLineAverage_on_DG_by_using_gPPI

subroutine Cal_FLA_FLabel_on_DG_by_using_gPPI_iniFLabel(DG_FLA_ResultData)
   use mpi_defs
   use precision
   use global_parameters
   use particle_decomp
   use FieldLineFollowing
   use field_array
   use Diagnostic_Grid
   implicit none

   real(wp), intent(out):: DG_FLA_ResultData(DG%dim_gq, DG%dim_a,0:1)
   real(wp) tmp_DG_Data(DG%dim_gq, DG%dim_a)

   real(wp) each_FLA_val(num_of_SampledFieldLines), each_FLA_denominator(num_of_SampledFieldLines), intpl_val ! FLA=Field Line Average
   real(wp) this_gq, this_s, this_gf, sdum, gqdum,ori_inv_ds
   real(wp) ws0, ws1, wgq00, wgq01, wgq10, wgq11, wgf0, wgf1
   integer Fid, Pid ! (Fid=Field ID), (Pid=gPPI ID)
   integer inner_sid, outter_sid, inner_gqid, outter_gqid, inner_ij, outter_ij
   integer ierror

   real(wp) recvr(mgrid)
   integer istatus(MPI_STATUS_SIZE)

   real(wp) target_gq, target_zeta
   real(wp) wgq0,wgq1,wgf,ratio
   integer ori_sid, j,k, sid, gqid, qq, ss, length
   real(wp) ff, this_safety_q



   if (num_of_SampledFieldLines ==0) then
      write(*,*) "Error! Please construct gPPI first before using FieldLineAverage!!"
      STOP
   endif

   if(myrank_toroidal==0) then
      write(*,*) "Do_FieldLineAverage_by_using_gPPI"
   endif
   ori_inv_ds = real(mpsi)/(s1-s0)

   ! Scatter Averaged data to Diagnostic Grid (DG) nodes
   DG_FLA_ResultData = 0.0
   each_FLA_val = 0.0
   each_FLA_denominator = 0.0
   do Fid=1,num_of_SampledFieldLines
      do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)
         each_FLA_val(Fid) = each_FLA_val(Fid) + gPPI(Pid,Fid)%ini_FLabel*gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
         each_FLA_denominator(Fid) = each_FLA_denominator(Fid) + gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
      enddo
   enddo

   ! Do Allreduce across toroidal_comm to get a sum of FLA
   call MPI_ALLREDUCE(each_FLA_val,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
   each_FLA_val = recvr(1:num_of_SampledFieldLines)
   call MPI_ALLREDUCE(each_FLA_denominator,recvr,num_of_SampledFieldLines,mpi_Rsize, MPI_SUM, toroidal_comm,ierror)
   each_FLA_denominator = recvr(1:num_of_SampledFieldLines)

   each_FLA_val =  each_FLA_val/ each_FLA_denominator ! Take average



   ! Scatter Averaged data to Diagnostic Grid (DG) nodes
   DG_FLA_ResultData = 0.0
   do Fid=1,num_of_SampledFieldLines
      do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

         this_s = gPPI(Pid,Fid)%s
         this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
         this_gf = gPPI(Pid,Fid)%target_normalized_gf

         ori_sid = (this_s-s0)*ori_inv_ds ! s index of original full grid
         ff  = (this_s-s_grid(ori_sid))*ori_inv_ds
         this_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)

         sid = max(0,min(int((this_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
         ws1 = (this_s-DG%s_1D(sid))*DG%inv_ds
         ws0 = 1.0-ws1

         do k=0,1
            if(k==0) target_zeta =zetamin
            if(k==1) target_zeta =zetamax

            wgf = (1.0-abs(this_gf-target_zeta)/deltaz)

            target_gq = modulo(this_gq + (target_zeta-this_gf)/this_safety_q, 2.0*pi)
            gqid = (target_gq)*DG%inv_dgq +1
            wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
            wgq0 = 1.0-wgq1

            do j=1,4
               select case(j)
                  case(1)
                     qq=gqid; ss=sid
                     ratio=ws0*wgq0*wgf
                  case(2)
                     qq=gqid+1; ss=sid
                     ratio=ws0*wgq1*wgf
                  case(3)
                     qq=gqid; ss=sid+1
                     ratio=ws1*wgq0*wgf
                  case(4)
                     qq=gqid+1; ss=sid+1
                     ratio=ws1*wgq1*wgf
               end select
               DG_FLA_ResultData(qq,ss,k)  = DG_FLA_ResultData(qq,ss,k)  + ratio*each_FLA_val(Fid)
            enddo
         enddo
      enddo
   enddo


   ! Match poloidal end cell
   ! (gq=0 and gq=2pi is the same position physically
   ! although they have different gqid)
   do sid=1,DG%dim_a
      DG_FLA_ResultData(1,sid,:) = DG_FLA_ResultData(1,sid,:) + DG_FLA_ResultData(DG%dim_gq,sid,:)
      DG_FLA_ResultData(DG%dim_gq,sid,:) = DG_FLA_ResultData(1,sid,:)
   enddo

   ! Communicate DG_FLA_ResultData between nodes
   ! toroidal end cell
   length = DG%dim_gq*DG%dim_a
   call MPI_SENDRECV(DG_FLA_ResultData(:,:,0),length,mpi_Rsize,left_pe,1, tmp_DG_Data,length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   DG_FLA_ResultData(:,:,1) = DG_FLA_ResultData(:,:,1) + tmp_DG_Data
   call MPI_SENDRECV(DG_FLA_ResultData(:,:,1),length,mpi_Rsize,right_pe,1, DG_FLA_ResultData(:,:,0),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

   ! Last part is to divide FLA_ResultData by DG%FLA_Denominator to get correct line-averaged FLA_ResultData
   DG_FLA_ResultData = DG_FLA_ResultData/DG%FLA_Denominator

end subroutine Cal_FLA_FLabel_on_DG_by_using_gPPI_iniFLabel

subroutine Get_FLA_FieldLabel_by_using_gPPI(FLA_ResultData)
   use mpi_defs
   use precision
   use global_parameters
   use particle_decomp
   use FieldLineFollowing
   use field_array

   implicit none

   real(wp) FLA_ResultData(0:1,mgrid)!, FLA_Denominator(0:1,mgrid)

   real(wp) each_FLA_val(num_of_SampledFieldLines), each_FLA_denominator(num_of_SampledFieldLines), intpl_val ! FLA=Field Line Average
   real(wp) this_a, this_gq, this_s, this_gf, sdum, gqdum,inv_ds, inv_dgq(0:mpsi), inv_gf
   real(wp) ws0, ws1, wgq00, wgq01, wgq10, wgq11, wgf0, wgf1
   integer Fid, Pid ! (Fid=Field ID), (Pid=gPPI ID)
   integer inner_sid, outter_sid, inner_gqid, outter_gqid, inner_ij, outter_ij
   integer i,ii,jt, ierror

   real(wp) tmp_weight
   real(wp) sendl(mgrid),recvr(mgrid)
   integer icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)


   if (num_of_SampledFieldLines ==0) then
      write(*,*) "Error! Please construct gPPI first before using FieldLineAverage!!"
      STOP
   endif

   if(myrank_toroidal==0) then
      write(*,*) "Do_FieldLineAverage_by_using_gPPI"
   endif
   inv_ds = real(mpsi)/(s1-s0)
   inv_dgq = 1.0/deltat
   inv_gf = 1.0/deltaz

   ! Scatter Averaged data to grid nodes
   FLA_ResultData = 0.0
   do Fid=1,num_of_SampledFieldLines
      do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

         this_s = gPPI(Pid,Fid)%s
         this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
         this_gf = gPPI(Pid,Fid)%target_normalized_gf

         sdum = (this_s-s0)*inv_ds
         inner_sid = int(sdum)
         ws1  = sdum - inner_sid
         ws0  = 1.0 - ws1
         gqdum = inv_dgq(inner_sid)*modulo(this_gq - qtinv(inner_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
         inner_gqid = int(gqdum)
         wgq01 = gqdum - inner_gqid
         wgq00 = 1.0 - wgq01
         inner_ij = igrid(inner_sid) + inner_gqid

         outter_sid = inner_sid+1
         gqdum = inv_dgq(outter_sid)*modulo(this_gq - qtinv(outter_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
         outter_gqid = int(gqdum)
         wgq11 = gqdum - outter_gqid
         wgq10 = 1.0 - wgq11
         outter_ij = igrid(outter_sid) + outter_gqid


         !wgf1 = min(modulo(this_gf-zetamin,2.0*pi),2.0*pi-modulo(this_gf-zetamin,2.0*pi))*inv_gf
         wgf1 = (this_gf-zetamin)*inv_gf
         wgf0 = 1.0 - wgf1

         !tmp_weight = gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
         !      FLA_ResultData(0,inner_ij)    = FLA_ResultData(0,inner_ij)    +wgf0*ws0*wgq00*each_FLA_val(Fid)*gPPI(Pid,Fid)%dl*gPPI(Pid,Fid)%inv_B
         !FLA_ResultData(0,inner_ij)    = FLA_ResultData(0,inner_ij)    +wgf0*ws0*wgq00*each_FLA_val(Fid)*tmp_weight
         FLA_ResultData(0,inner_ij)    = FLA_ResultData(0,inner_ij)    +wgf0*ws0*wgq00*gPPI(Pid,Fid)%ini_FLabel
         FLA_ResultData(0,inner_ij+1)  = FLA_ResultData(0,inner_ij+1)  +wgf0*ws0*wgq01*gPPI(Pid,Fid)%ini_FLabel
         FLA_ResultData(0,outter_ij)   = FLA_ResultData(0,outter_ij)   +wgf0*ws1*wgq10*gPPI(Pid,Fid)%ini_FLabel
         FLA_ResultData(0,outter_ij+1) = FLA_ResultData(0,outter_ij+1) +wgf0*ws1*wgq11*gPPI(Pid,Fid)%ini_FLabel
         FLA_ResultData(1,inner_ij)    = FLA_ResultData(1,inner_ij)    +wgf1*ws0*wgq00*gPPI(Pid,Fid)%ini_FLabel
         FLA_ResultData(1,inner_ij+1)  = FLA_ResultData(1,inner_ij+1)  +wgf1*ws0*wgq01*gPPI(Pid,Fid)%ini_FLabel
         FLA_ResultData(1,outter_ij)   = FLA_ResultData(1,outter_ij)   +wgf1*ws1*wgq10*gPPI(Pid,Fid)%ini_FLabel
         FLA_ResultData(1,outter_ij+1) = FLA_ResultData(1,outter_ij+1) +wgf1*ws1*wgq11*gPPI(Pid,Fid)%ini_FLabel
      enddo
   enddo

   ! This part is the required post-processing and/or communication to get correct FLA_ResultData
   ! poloidal end cell; need valid data for densityi(:,igrid(i))
   do i=0,mpsi
      FLA_ResultData(:,igrid(i)+mtheta(i))=FLA_ResultData(:,igrid(i)+mtheta(i))+FLA_ResultData(:,igrid(i))
      FLA_ResultData(:,igrid(i))=FLA_ResultData(:,igrid(i)+mtheta(i))
   enddo
   ! Communicate FLA_ResultData between nodes
   ! toroidal end cell
   sendl=FLA_ResultData(0,:)
   recvr=0.0_wp
   icount=mgrid
   idest=left_pe
   isource=right_pe
   isendtag=myrank_toroidal
   irecvtag=isource
   ! send densityi to left and receive from right
   call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,recvr,icount,&
      mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      ! B.C. at zeta=2*pi is shifted
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         FLA_ResultData(1,ii+1:ii+jt)=FLA_ResultData(1,ii+1:ii+jt)+ cshift(recvr(ii+1:ii+jt),itran(i))
         FLA_ResultData(1,igrid(i))=FLA_ResultData(1,igrid(i)+mtheta(i))
      enddo
   else
      ! B.C. at zeta<2*pi is continuous
      FLA_ResultData(1,:)=FLA_ResultData(1,:)+recvr
   endif


   ! Last part is to divide FLA_ResultData by FLA_Denominator to get correct line-averaged FLA_ResultData
   FLA_ResultData = FLA_ResultData/FLA_Denominator

end subroutine Get_FLA_FieldLabel_by_using_gPPI

subroutine Cal_FLA_Denominator_on_DG_using_gPPI
   use mpi_defs
   use precision
   use global_parameters
   use particle_decomp
   use FieldLineFollowing
   use field_array
   use Diagnostic_Grid

   implicit none

   real(wp) tmp_Data(DG%dim_gq,DG%dim_a)
   real(wp) this_gq, this_s, this_gf, ori_inv_ds, this_safety_q
   real(wp) target_gq, target_zeta
   real(wp) ff, ws0,ws1, wgq0,wgq1, wgf,ratio ! ratios along each dimension
   integer Fid, Pid ! (Fid=Field ID), (Pid=gPPI ID)
   integer :: ori_sid, j,k, sid, gqid, qq, ss, length
   integer istatus(MPI_STATUS_SIZE),ierror


   if (num_of_SampledFieldLines ==0) then
      write(*,*) "Error! Please construct gPPI first before using FieldLineAverage!!"
      STOP
   endif

   ori_inv_ds = real(mpsi)/(s1-s0)

   ! Scatter Averaged data to grid nodes
   DG%FLA_Denominator = 0.0
   do Fid=1,num_of_SampledFieldLines
      do Pid=1,num_of_gPPI_for_each_FieldLine(Fid)

         this_s = gPPI(Pid,Fid)%s
         this_gq = modulo(gPPI(Pid,Fid)%gq, 2.0*pi)
         this_gf = gPPI(Pid,Fid)%target_normalized_gf

         ori_sid = (this_s-s0)*ori_inv_ds ! s index of original full grid
         ff  = (this_s-s_grid(ori_sid))*ori_inv_ds
         this_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)


         sid = max(0,min(int((this_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
         ws1 = (this_s-DG%s_1D(sid))*DG%inv_ds
         ws0 = 1.0-ws1

         do k=0,1
            if(k==0) target_zeta =zetamin
            if(k==1) target_zeta =zetamax

            wgf = (1.0-abs(this_gf-target_zeta)/deltaz)

            target_gq = modulo(this_gq + (target_zeta-this_gf)/this_safety_q, 2.0*pi)
            gqid = (target_gq)*DG%inv_dgq +1
            wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
            wgq0 = 1.0-wgq1

            do j=1,4
               select case(j)
                  case(1)
                     qq=gqid; ss=sid
                     ratio=ws0*wgq0*wgf
                  case(2)
                     qq=gqid+1; ss=sid
                     ratio=ws0*wgq1*wgf
                  case(3)
                     qq=gqid; ss=sid+1
                     ratio=ws1*wgq0*wgf
                  case(4)
                     qq=gqid+1; ss=sid+1
                     ratio=ws1*wgq1*wgf
               end select
               DG%FLA_Denominator(qq,ss,k) = DG%FLA_Denominator(qq,ss,k) + ratio
            enddo
         enddo
      enddo
   enddo

   ! Match poloidal end cell
   ! (gq=0 and gq=2pi is the same position physically
   ! although they have different gqid)
   do sid=1,DG%dim_a
      DG%FLA_Denominator(1,sid,:) = DG%FLA_Denominator(1,sid,:) + DG%FLA_Denominator(DG%dim_gq,sid,:)
      DG%FLA_Denominator(DG%dim_gq,sid,:) = DG%FLA_Denominator(1,sid,:)
   enddo

   ! Communicate DG%FLA_FLabel between nodes
   ! toroidal end cell
   length = DG%dim_gq*DG%dim_a
   call MPI_SENDRECV(DG%FLA_Denominator(:,:,0),length,mpi_Rsize,left_pe,1, tmp_Data,length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   DG%FLA_Denominator(:,:,1) = DG%FLA_Denominator(:,:,1) + tmp_Data
   call MPI_SENDRECV(DG%FLA_Denominator(:,:,1),length,mpi_Rsize,right_pe,1, DG%FLA_Denominator(:,:,0),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

end subroutine Cal_FLA_Denominator_on_DG_using_gPPI

recursive subroutine quicksort_integer(a, first, last)
   ! quicksort.f -*-f90-*-
   ! Author: t-nissie
   ! License: GPLv3
   ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
   !!
   use precision
   implicit none

   integer  a(*), x, t
   integer first, last
   integer i, j

   x = a( (first+last) / 2 )
   i = first
   j = last
   do
      do while (a(i) < x)
         i=i+1
      end do

      do while (x < a(j))
         j=j-1
      end do

      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      i=i+1
      j=j-1
   end do

   if (first < i-1) call quicksort_integer(a, first, i-1)
   if (j+1 < last)  call quicksort_integer(a, j+1, last)
end subroutine quicksort_integer

recursive subroutine quicksort_real(a, first, last)
   ! quicksort.f -*-f90-*-
   ! Author: t-nissie
   ! License: GPLv3
   ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
   !!
   use precision
   implicit none

   real(wp)  a(*), x, t
   integer first, last
   integer i, j

   x = a( (first+last) / 2 )
   i = first
   j = last
   do
      do while (a(i) < x)
         i=i+1
      end do

      do while (x < a(j))
         j=j-1
      end do

      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      i=i+1
      j=j-1
   end do

   if (first < i-1) call quicksort_real(a, first, i-1)
   if (j+1 < last)  call quicksort_real(a, j+1, last)
end subroutine quicksort_real


recursive subroutine quicksort_real_together(a, b, first, last)
   ! quicksort.f -*-f90-*-
   ! Author: t-nissie
   ! License: GPLv3
   ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
   !!
   use precision
   implicit none

   real(wp)  a(*), b(*), x, t
   integer first, last
   integer i, j

   x = a( (first+last) / 2 )
   i = first
   j = last
   do
      do while (a(i) < x)
         i=i+1
      end do

      do while (x < a(j))
         j=j-1
      end do

      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      t = b(i);  b(i) = b(j);  b(j) = t
      i=i+1
      j=j-1
   end do

   if (first < i-1) call quicksort_real_together(a, b, first, i-1)
   if (j+1 < last)  call quicksort_real_together(a, b, j+1, last)
end subroutine quicksort_real_together


subroutine find_unique_sorted_array_and_its_num_of_elements(a, first, last, num_of_unique_elements)

   use precision
   implicit none

   real(wp) a(*), x, t
   integer first, last, i, j, num_of_unique_elements

   call quicksort_real(a, first, last)

   i=first
   j=i+1
   do while(i<last)
      do while(j<=last .and. a(j)==a(i))
         j=j+1
      enddo
      if(j>last) exit
      if(j/=i+1) then
         a(i+1)=a(j)
      endif
      i=i+1
   enddo
   num_of_unique_elements = i
end subroutine find_unique_sorted_array_and_its_num_of_elements

subroutine Convert_a_to_s(in_a, out_s)
   use global_parameters

   implicit none

   real(wp) in_a, out_s
   integer ierr, IFAIL, minpts
   real(wp) error_8, acc
   real(wp), external :: D01AHF,dsda

   if (dsda_power == 0) then
      out_s = s0 + (in_a-a0) !!
   else
      error_8=1.0e-5_wp
      IFAIL=1
      out_s=D01AHF(0.0_wp,in_a,error_8,minpts,acc,dsda,-1,IFAIL)
      out_s=min(s1,max(s0,out_s))
   endif
end subroutine Convert_a_to_s

subroutine Cal_Weighted_Neighboring_Volume_on_each_GridNode

   use mpi_defs
   use precision
   use global_parameters
   use particle_decomp
   use field_array

   implicit none

   real(wp),external :: jacobian

   integer i,j,k
   character(len=64) chdum
   !integer, parameter:: SAMP_a_num=400, SAMP_gq_num=500, SAMP_gf_num=10 ! number of sampling points for each dimension
   !integer, parameter:: SAMP_a_num=1000, SAMP_gq_num=4000, SAMP_gf_num=20 ! number of sampling points for each dimension
   ! integer, parameter:: SAMP_a_num=mpsi*20, SAMP_gq_num=mtheta(ipsi)*20, SAMP_gf_num=20 ! number of sampling points for each dimension
   ! real(wp) SAMP_X_vec_arr(3, SAMP_a_num*SAMP_gq_num*SAMP_gf_num) ! (a, gq, gf) arr
   integer SAMP_a_num, SAMP_gq_num, SAMP_gf_num ! number of sampling points for each dimension
   !real(wp), dimension(:,:), allocatable:: SAMP_X_vec_arr! (a, gq, gf) arr

   real(wp) SAMP_da, SAMP_dgq, SAMP_dgf, SAMP_d3x
   integer SAMP_id, SAMP_NUM
   real(wp) dpsi, this_ds

   real(wp) this_a, this_s, this_gq, this_gf, this_Jacob
   real(wp) inv_ds, inv_dgq(0:mpsi), inv_dgf
   real sdum, gqdum, ws0, ws1, wgq00, wgq01, wgq10, wgq11, wgf0, wgf1
   integer inner_sid, outter_sid, inner_gqid, outter_gqid, inner_ij, outter_ij
   real(wp), external :: dsda

   integer ij,ii,jt, ierror
   real(wp) sendl(mgrid),recvr(mgrid)
   integer icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
   real(wp) t_start, t_finish

   call cpu_time(t_start)

   SAMP_a_num  = mpsi*20
   SAMP_gq_num = mtheta(ipsi)*20
   SAMP_gf_num = 20
   !allocate(SAMP_X_vec_arr(3, SAMP_a_num*SAMP_gq_num*SAMP_gf_num))

   SAMP_da  = (a1-a0)/(SAMP_a_num)
   SAMP_dgq = 2.0*pi/(SAMP_gq_num)
   SAMP_dgf = (zetamax-zetamin)/(SAMP_gf_num)

   SAMP_d3x = ((a1-a0)*2.0*pi*(zetamax-zetamin))/real(SAMP_a_num*SAMP_gq_num*SAMP_gf_num)


   ! SAMP_NUM = 0
   ! do i=0, SAMP_a_num-1
   !   do j=0, SAMP_gq_num-1
   !     do k=0, SAMP_gf_num-1
   !       SAMP_NUM=SAMP_NUM+1
   !       SAMP_X_vec_arr(1,SAMP_NUM) = a0 + i*SAMP_da
   !       SAMP_X_vec_arr(2,SAMP_NUM) = j*SAMP_dgq
   !       SAMP_X_vec_arr(3,SAMP_NUM) = zetamin + k*SAMP_dgf
   !     enddo
   !   enddo
   ! enddo

   inv_ds = real(mpsi)/(s1-s0)
   inv_dgq = 1.0/deltat
   inv_dgf = 1.0/deltaz

   Weighted_Vol = 0
   !do SAMP_id=1, SAMP_NUM
   do i=0, SAMP_a_num-1
      do j=0, SAMP_gq_num-1
         do k=0, SAMP_gf_num-1

            ! this_a  = SAMP_X_vec_arr(1,SAMP_id)
            ! this_gq = SAMP_X_vec_arr(2,SAMP_id)
            ! this_gf = SAMP_X_vec_arr(3,SAMP_id)
            this_a  = a0 + (i+0.5)*SAMP_da
            this_gq = (j+0.5)*SAMP_dgq
            this_gf = zetamin + (k+0.5)*SAMP_dgf
            this_Jacob = jacobian(this_a, this_gq)

            call Convert_a_to_s(this_a, this_s)

            sdum = (this_s-s0)*inv_ds
            inner_sid = int(sdum)
            ws1  = sdum - inner_sid
            ws0  = 1.0 - ws1
            gqdum = inv_dgq(inner_sid)*modulo(this_gq - qtinv(inner_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
            inner_gqid = int(gqdum)
            wgq01 = gqdum - inner_gqid
            wgq00 = 1.0 - wgq01
            inner_ij = igrid(inner_sid) + inner_gqid

            outter_sid = inner_sid+1
            gqdum = inv_dgq(outter_sid)*modulo(this_gq - qtinv(outter_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
            outter_gqid = int(gqdum)
            wgq11 = gqdum - outter_gqid
            wgq10 = 1.0 - wgq11
            outter_ij = igrid(outter_sid) + outter_gqid

            wgf1 = (this_gf-zetamin)*inv_dgf
            wgf0 = 1.0 - wgf1

            Weighted_Vol(0,inner_ij)    = Weighted_Vol(0,inner_ij)    +wgf0*ws0*wgq00*this_Jacob
            Weighted_Vol(0,inner_ij+1)  = Weighted_Vol(0,inner_ij+1)  +wgf0*ws0*wgq01*this_Jacob
            Weighted_Vol(0,outter_ij)   = Weighted_Vol(0,outter_ij)   +wgf0*ws1*wgq10*this_Jacob
            Weighted_Vol(0,outter_ij+1) = Weighted_Vol(0,outter_ij+1) +wgf0*ws1*wgq11*this_Jacob
            Weighted_Vol(1,inner_ij)    = Weighted_Vol(1,inner_ij)    +wgf1*ws0*wgq00*this_Jacob
            Weighted_Vol(1,inner_ij+1)  = Weighted_Vol(1,inner_ij+1)  +wgf1*ws0*wgq01*this_Jacob
            Weighted_Vol(1,outter_ij)   = Weighted_Vol(1,outter_ij)   +wgf1*ws1*wgq10*this_Jacob
            Weighted_Vol(1,outter_ij+1) = Weighted_Vol(1,outter_ij+1) +wgf1*ws1*wgq11*this_Jacob

            !enddo
         enddo
      enddo
   enddo

   Weighted_Vol = Weighted_Vol * SAMP_d3x ! Calculate the real volume by multiplying d3x, so that (Vol = sum(Jacob*da*dgq*dgf))

   do i=0,mpsi
      Weighted_Vol(:,igrid(i)+mtheta(i))=Weighted_Vol(:,igrid(i)+mtheta(i))+Weighted_Vol(:,igrid(i))
      Weighted_Vol(:,igrid(i))=Weighted_Vol(:,igrid(i)+mtheta(i))
   enddo
   ! Communicate Weighted_Vol between nodes
   ! toroidal end cell
   sendl=Weighted_Vol(0,:)
   recvr=0.0_wp
   icount=mgrid
   idest=left_pe
   isource=right_pe
   isendtag=myrank_toroidal
   irecvtag=isource
   ! send densityi to left and receive from right
   call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,recvr,icount,&
      mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      ! B.C. at zeta=2*pi is shifted
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         Weighted_Vol(1,ii+1:ii+jt)=Weighted_Vol(1,ii+1:ii+jt)+ cshift(recvr(ii+1:ii+jt),itran(i))
         Weighted_Vol(1,igrid(i))=Weighted_Vol(1,igrid(i)+mtheta(i))
      enddo
   else
      ! B.C. at zeta<2*pi is continuous
      Weighted_Vol(1,:)=Weighted_Vol(1,:)+recvr
   endif



   call cpu_time(t_finish)
   write(stdout,*) 'Im=',mype,  'Cal_Vol time =',t_finish-t_start,' seconds'
   flush(stdout)

   if (myrank_partd==0) then
      write(chdum,'("Weighted_Vol_",i5.5,".dat")') myrank_toroidal
      open(181215,file=chdum,status='replace',action='write')
      do i=0,mpsi
         do j=0, mtheta(i)
            ij=igrid(i) + j
            if(i==0) then
               dpsi=(psi_grid(1)-psi_grid(0))*0.5_wp         ! half cell
               this_ds = 0.5*(1.0/inv_ds)
            else if(i==mpsi) then
               dpsi=(psi_grid(mpsi)-psi_grid(mpsi-1))*0.5_wp ! half cell
               this_ds = 0.5*(1.0/inv_ds)
            else
               dpsi=(psi_grid(i+1)-psi_grid(i-1))*0.5_wp     ! full cell
               this_ds = (1.0/inv_ds)
            endif
            !write(181215, "(2(ES16.8,1x))") Weighted_Vol(1,ij), jacobian(psi_grid(i), modulo(j*deltat(i) + qtinv(i)*zetamax,2.0*pi))*deltat(i)*dpsi*deltaz
            write(181215, "(3(ES16.8,1x))") Weighted_Vol(1,ij),&
               jacobian(psi_grid(i), modulo(j*deltat(i) + qtinv(i)*zetamax,2.0*pi))*deltat(i)*dpsi*deltaz, &
               (1.0/dsda(psi_grid(i)))*jacobian(psi_grid(i), modulo(j*deltat(i) + qtinv(i)*zetamax,2.0*pi))*deltat(i)*this_ds*deltaz
         enddo
      enddo
      flush(181215)
      close(181215)
   endif

   !deallocate(SAMP_X_vec_arr)

end subroutine Cal_Weighted_Neighboring_Volume_on_each_GridNode


subroutine Cal_Weighted_Neighboring_Volume_on_each_GridNode_by_MCI

   use mpi_defs
   use precision
   use global_parameters
   use particle_decomp
   use field_array

   implicit none

   real(wp),external :: jacobian

   integer i,j,k
   character(len=64) chdum


   real(wp) dpsi, this_ds
   real(wp) this_a, this_s, this_gq, this_gf, this_Jacob
   real(wp) inv_ds, inv_dgq(0:mpsi), inv_dgf
   real sdum, gqdum, ws0, ws1, wgq00, wgq01, wgq10, wgq11, wgf0, wgf1
   integer inner_sid, outter_sid, inner_gqid, outter_gqid, inner_ij, outter_ij
   real(wp), external :: dsda

   integer ij,ii,jt, ierror
   real(wp) sendl(mgrid),recvr(mgrid)
   integer icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)

   integer(kind=8) :: sobol_seed=1 ! sobol sedd should not be parameter, becuase sobol_seed changes at every sampling
   integer(kind=8), parameter :: sobol_dim_num=3
   real(kind=8) sobol_rand_3d(sobol_dim_num)
   integer :: iteration=0, max_iteration = 1e9
   real(wp) prev_Weighted_Vol(0:1, mgrid), now_Weighted_Vol(0:1, mgrid), d3X_factor
   real(wp) :: max_rel_err, rel_tolerance=1e-4
   real(wp) t_start, t_finish

   call cpu_time(t_start)

   inv_ds = real(mpsi)/(s1-s0)
   inv_dgq = 1.0/deltat
   inv_dgf = 1.0/deltaz

   Weighted_Vol = 0

   ! SAMP_da  = (a1-a0)/(SAMP_a_num-1.0)
   ! SAMP_dgq = 2.0*pi/(SAMP_gq_num-1.0)
   ! SAMP_dgf = (zetamax-zetamin)/(SAMP_gf_num-1.0)

   ! SAMP_d3x = SAMP_da*SAMP_dgq*SAMP_dgf ! da*dgq*dgf


   ! SAMP_NUM = 0
   ! do i=0, SAMP_a_num-1
   !   do j=0, SAMP_gq_num-1
   !     do k=0, SAMP_gf_num-1
   !       SAMP_NUM=SAMP_NUM+1
   !       SAMP_X_vec_arr(1,SAMP_NUM) = a0 + i*SAMP_da
   !       SAMP_X_vec_arr(2,SAMP_NUM) = j*SAMP_dgq
   !       SAMP_X_vec_arr(3,SAMP_NUM) = zetamin + k*SAMP_dgf
   !     enddo
   !   enddo
   ! enddo


   Weighted_Vol = 0
   prev_Weighted_Vol = 0
   do while( iteration < max_iteration)

      call i8_sobol(sobol_dim_num, sobol_seed, sobol_rand_3d)
      this_a = a0 + (a1-a0)*sobol_rand_3d(1)
      this_gq = 2.0*pi*sobol_rand_3d(2)
      this_gf = zetamin + (zetamax-zetamin)*sobol_rand_3d(3)
      this_Jacob = jacobian(this_a, this_gq)
      call Convert_a_to_s(this_a, this_s)

      sdum = (this_s-s0)*inv_ds
      inner_sid = int(sdum)
      ws1  = sdum - inner_sid
      ws0  = 1.0 - ws1
      gqdum = inv_dgq(inner_sid)*modulo(this_gq - qtinv(inner_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
      inner_gqid = int(gqdum)
      wgq01 = gqdum - inner_gqid
      wgq00 = 1.0 - wgq01
      inner_ij = igrid(inner_sid) + inner_gqid

      outter_sid = inner_sid+1
      gqdum = inv_dgq(outter_sid)*modulo(this_gq - qtinv(outter_sid)*(this_gf), 2.0*pi) ! to consider initial gq at this poloidal plane.
      outter_gqid = int(gqdum)
      wgq11 = gqdum - outter_gqid
      wgq10 = 1.0 - wgq11
      outter_ij = igrid(outter_sid) + outter_gqid

      wgf1 = (this_gf-zetamin)*inv_dgf
      wgf0 = 1.0 - wgf1

      Weighted_Vol(0,inner_ij)    = Weighted_Vol(0,inner_ij)    +wgf0*ws0*wgq00*this_Jacob
      Weighted_Vol(0,inner_ij+1)  = Weighted_Vol(0,inner_ij+1)  +wgf0*ws0*wgq01*this_Jacob
      Weighted_Vol(0,outter_ij)   = Weighted_Vol(0,outter_ij)   +wgf0*ws1*wgq10*this_Jacob
      Weighted_Vol(0,outter_ij+1) = Weighted_Vol(0,outter_ij+1) +wgf0*ws1*wgq11*this_Jacob
      Weighted_Vol(1,inner_ij)    = Weighted_Vol(1,inner_ij)    +wgf1*ws0*wgq00*this_Jacob
      Weighted_Vol(1,inner_ij+1)  = Weighted_Vol(1,inner_ij+1)  +wgf1*ws0*wgq01*this_Jacob
      Weighted_Vol(1,outter_ij)   = Weighted_Vol(1,outter_ij)   +wgf1*ws1*wgq10*this_Jacob
      Weighted_Vol(1,outter_ij+1) = Weighted_Vol(1,outter_ij+1) +wgf1*ws1*wgq11*this_Jacob

      iteration = iteration + 1
      if( mod(iteration, (max_iteration/100))==0) then

         d3X_factor = ((a1-a0)*2.0*pi*(zetamax-zetamin))/real(iteration)
         now_Weighted_Vol = Weighted_Vol*d3X_factor
         max_rel_err =  maxval( abs((now_Weighted_Vol - prev_Weighted_Vol)/now_Weighted_Vol))

         write(stdout,*) 'iteration=', iteration, 'max_rel_err', max_rel_err
         if( max_rel_err < rel_tolerance) then
            exit
         endif
         prev_Weighted_Vol = now_Weighted_Vol
      endif
   enddo


   d3X_factor = ((a1-a0)*2.0*pi*(zetamax-zetamin))/real(iteration)
   Weighted_Vol = Weighted_Vol*d3X_factor

   if( iteration >= max_iteration) then
      write(stdout,*) 'Although upto #', max_iteration,  'iterations, max_rel_err = ', max_rel_err
      flush(stdout)
   endif

   !  Weighted_Vol = Weighted_Vol * SAMP_d3x ! Calculate the real volume by multiplying d3x, so that (Vol = sum(Jacob*da*dgq*dgf))
   do i=0,mpsi
      Weighted_Vol(:,igrid(i)+mtheta(i))=Weighted_Vol(:,igrid(i)+mtheta(i))+Weighted_Vol(:,igrid(i))
      Weighted_Vol(:,igrid(i))=Weighted_Vol(:,igrid(i)+mtheta(i))
   enddo
   ! Communicate Weighted_Vol between nodes
   ! toroidal end cell
   sendl=Weighted_Vol(0,:)
   recvr=0.0_wp
   icount=mgrid
   idest=left_pe
   isource=right_pe
   isendtag=myrank_toroidal
   irecvtag=isource
   ! send densityi to left and receive from right
   call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,recvr,icount,&
      mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
   if(myrank_toroidal == ntoroidal-1)then
      ! B.C. at zeta=2*pi is shifted
      do i=0,mpsi
         ii=igrid(i)
         jt=mtheta(i)
         Weighted_Vol(1,ii+1:ii+jt)=Weighted_Vol(1,ii+1:ii+jt)+ cshift(recvr(ii+1:ii+jt),itran(i))
         Weighted_Vol(1,igrid(i))=Weighted_Vol(1,igrid(i)+mtheta(i))
      enddo
   else
      ! B.C. at zeta<2*pi is continuous
      Weighted_Vol(1,:)=Weighted_Vol(1,:)+recvr
   endif


   call cpu_time(t_finish)
   write(stdout,*) 'Im=',mype,  'Cal_Vol_by_MCI time =',t_finish-t_start,' seconds'
   flush(stdout)
   if (myrank_partd==0) then
      write(chdum,'("Weighted_Vol_by_MCI_",i5.5,".dat")') myrank_toroidal
      open(181218,file=chdum,status='replace',action='write')
      do i=0,mpsi
         do j=0, mtheta(i)
            ij=igrid(i) + j
            if(i==0) then
               dpsi=(psi_grid(1)-psi_grid(0))*0.5_wp         ! half cell
               this_ds = 0.5*(1.0/inv_ds)
            else if(i==mpsi) then
               dpsi=(psi_grid(mpsi)-psi_grid(mpsi-1))*0.5_wp ! half cell
               this_ds = 0.5*(1.0/inv_ds)
            else
               dpsi=(psi_grid(i+1)-psi_grid(i-1))*0.5_wp     ! full cell
               this_ds = (1.0/inv_ds)
            endif
            write(181218, "(3(ES16.8,1x))") Weighted_Vol(1,ij),&
               jacobian(psi_grid(i), modulo(j*deltat(i) + qtinv(i)*zetamax,2.0*pi))*deltat(i)*dpsi*deltaz, &
               (1.0/dsda(psi_grid(i)))*jacobian(psi_grid(i), modulo(j*deltat(i) + qtinv(i)*zetamax,2.0*pi))*deltat(i)*this_ds*deltaz
         enddo
      enddo
      flush(181218)
      close(181218)
   endif


end subroutine Cal_Weighted_Neighboring_Volume_on_each_GridNode_by_MCI

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

   if(acc_newer_H_zeon) then
      !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
      acc_newer_H_zeon=.false.
      acc_newer_D_zeon=.false.
   endif
   if(plasma==1) then    ! use experimental ni(r), T(r) and omega(r) profiles
      if(linterp_from_nodes==1) then
         dela=(a1-a0)/real(jr,wp)
         dela=1.0_wp/dela

         !$omp parallel do private(i,r,ip,fip)
         !$acc parallel loop present(zeon)
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
         !$acc parallel loop present(zeon,backT)
         do i=1,me
            backT(i)=T_0*Bt*exp(-dTemi*zeon(1,i)**xt)
            backT(i)=Tc*(backT(i)/Tc)**r0
         end do
      else
         dela=(a1-a0)/real(jr,wp)
         dela=1.0_wp/dela
         !$omp parallel do private(i,r,ip,fip,tempi)
         !$acc parallel loop present(zeon,backT)
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

   !$acc update self(backT(:me))

end subroutine Ti_at_electrons

module function_gg_module

   contains

      FUNCTION SMOOTHSTEP_M(xmin,xmax,x)
         !$acc routine seq
         use precision
         !      use global_parameters,only: pi
         !
         !      SMOOTHSTEP generates a smooth step function using
         !      SMOOTHSTEP=1/2-cos(pi*x)/2
         !      SMOOTHSTEP=0 if x<xmin; =1 if x>xmax; otherwise 0<SMOOTHSTEP<1
         !      with xmax > xmin
         !

         implicit  none
         real(wp) SMOOTHSTEP_M
         real(wp) xmin,xmax,x,y
         real(wp) pi

         pi = acos(-1.0_wp)
         y=(x-xmin)/(xmax-xmin)
         y=max(0.0_wp,min(1.0_wp,y))
         SMOOTHSTEP_M=0.5_wp-0.5_wp*cos(pi*y)

         RETURN
      end function SMOOTHSTEP_M
end module function_gg_module

FUNCTION SMOOTHSTEP(xmin,xmax,x)
   use function_gg_module, only : smoothstep_m
   use precision
   !
   !      SMOOTHSTEP generates a smooth step function using
   !      SMOOTHSTEP=1/2-cos(pi*x)/2
   !      SMOOTHSTEP=0 if x<xmin; =1 if x>xmax; otherwise 0<SMOOTHSTEP<1
   !      with xmax > xmin
   !

   implicit  none
   real(wp) SMOOTHSTEP
   real(wp) xmin,xmax,x

   SMOOTHSTEP = smoothstep_m (xmin,xmax,x)

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


subroutine Init_Diagnostic_Grid

   use Diagnostic_Grid
   use DG_operations
   use global_parameters
   use field_array
   use precision
   use particle_decomp
   use mpi_defs
   implicit none

   integer i,j,k,aid,gqid,ij,Na,Nq
   integer pbid ! probe id
   real(wp), external :: temperature, temperature_e, density, Cal_FLabel_Function, dsda

   ! approximated delta_a and delta_gf for approximated calculation of the weighted cell volume
   real(wp) approx_da, approx_dgf, tmp
   real(wp), allocatable, dimension(:,:,:) :: tmp_DG_Data1, tmp_DG_Data2
   integer ierror, istatus(MPI_STATUS_SIZE)
   real(wp) uvec_B_cov(3), uvec_N_cov(3), uvec_BxN_cov(3)

   DG%dim_a  = (mpsi+1)/DG%a_jump
   DG%dim_gq = mtheta(ipsi)/DG%gq_jump +1
   DG%num_eSD = 6 +6 +6
   DG%num_eVD = 5

   DG%num_iSD = 8 +6 +6! 2 additional dni_Larmor, dTi_Larmor by considering larmor orbit of ions
   DG%num_iVD = 5+4 ! 4 additional energy, particle flows (total and only ExB)

   DG%dgq   = 2.0*pi/(DG%dim_gq-1)
   DG%inv_dgq = 1.0/DG%dgq
   DG%inv_ds  = DG%a_jump*(real(mpsi)/(s1-s0))

   DG%zetamax = zetamax
   DG%zetamin = zetamin

   Na=DG%dim_a
   Nq=DG%dim_gq

   allocate(DG%gq_1D(Nq))
   allocate(DG%a_1D(Na), DG%ori_aid(Na), DG%s_1D(Na), DG%n0_1D(Na), DG%Te0_1D(Na), DG%Ti0_1D(Na))
   allocate(DG%a_2D(Nq,Na))
   allocate(DG%gq_2D(Nq,Na), DG%Cell_dV(Nq,Na), DG%inv_Cell_dV(Nq,Na), DG%dne_adiab(Nq,Na))
   allocate(DG%Cell_dS_a(Nq,Na), DG%inv_Cell_dS_a(Nq,Na), DG%Cell_dS_gq(Nq,Na), DG%inv_Cell_dS_gq(Nq,Na), DG%Cell_dS_gf(Nq,Na), DG%inv_Cell_dS_gf(Nq,Na))
   allocate(DG%FLabel_Function(Nq,Na,0:1))
   allocate(DG%FLA_Denominator(Nq,Na,0:1), DG%FLA_FLabel(Nq,Na,0:1), DG%FLA_psi0(Nq,Na,0:1))
   allocate(DG%FLA_SM_Denominator(Nq,Na), DG%FLA_SM_FLabel(Nq,Na), DG%FLA_SM_psi0(Nq,Na))
   allocate(DG%FLA_SM_FLabel_Left(Nq,Na), DG%FLA_SM_FLabel_Right(Nq,Na))
   allocate(DG%grad_FLA_SM_FLabel(3,Nq,Na,0:1),DG%uvec_grad_FL_cov(3,Nq,Na,0:1), DG%mag_grad_FLA_SM_FLabel(Nq,Na,0:1))

   allocate(DG%eScaData(Nq,Na,DG%num_eSD,0:1), DG%eVecData(3,Nq,Na,DG%num_eVD,0:1) )
   allocate(DG%iScaData(Nq,Na,DG%num_iSD,0:1), DG%iVecData(3,Nq,Na,DG%num_iVD,0:1) )


   k=1
   do i=0,mpsi,DG%a_jump
      DG%ori_aid(k) = i
      DG%a_1D(k) = psi_grid(i)
      DG%a_2D(:,k) = psi_grid(i)
      DG%s_1D(k) = s_grid(i)

      DG%n0_1D(k)  = den_grid(i)
      DG%Te0_1D(k) = Te_grid(i)
      DG%Ti0_1D(k) = vth_grid(i)**2.0
      k=k+1
   enddo

   do i=1,Nq
      DG%gq_1D(i) = (2.0*pi)/(Nq-1)*(i-1)
      DG%gq_2D(i,:) = DG%gq_1D(i)
   enddo

   ! Calculate FLabel Function and add FLA_psi0
   do j=1,Nq
      do i=1,Na
         DG%FLabel_Function(j,i,0) = Cal_FLabel_Function(DG%a_2D(j,i), DG%gq_2D(j,i), zetamin)
         DG%FLabel_Function(j,i,1) = Cal_FLabel_Function(DG%a_2D(j,i), DG%gq_2D(j,i), zetamax)
      enddo
   enddo

   !Set and Calculate Magnetic Information (ESI variables)
   allocate(DG%Jacobian(Nq,Na), DG%F(Nq,Na),DG%dFda(Nq,Na),DG%dgFda(Nq,Na),DG%dgYda(Nq,Na),DG%R(Nq,Na),DG%dRda(Nq,Na),DG%dRdgq(Nq,Na),DG%Z(Nq,Na),DG%dZda(Nq,Na),DG%dZdgq(Nq,Na),DG%BB(Nq,Na),DG%dBBda(Nq,Na),DG%dBBdgq(Nq,Na),DG%P(Nq,Na),DG%T(Nq,Na))
   allocate(DG%esi_N(Nq,Na),DG%esi_K(Nq,Na),DG%g_aa_cov(Nq,Na),DG%g_aq_cov(Nq,Na),DG%g_qq_cov(Nq,Na),DG%g_ff_cov(Nq,Na))
   allocate(DG%D(Nq,Na),DG%g_aa_cntv(Nq,Na),DG%g_aq_cntv(Nq,Na),DG%g_qq_cntv(Nq,Na),DG%g_ff_cntv(Nq,Na))
   allocate(DG%B0_cov3(3,Nq,Na))
   allocate(DG%B0_cntv3(3,Nq,Na),DG%j0_cov3(3,Nq,Na),DG%j0_cntv3(3,Nq,Na))

   allocate(DG%g_ss_cntv(Nq,Na),DG%g_sq_cntv(Nq,Na),DG%g_ss_cov(Nq,Na),DG%g_sq_cov(Nq,Na),DG%Jacob_sqf(Nq,Na))
   allocate(DG%B0_cntv3_sqf(3,Nq,Na))

   call esiget2dthese_more(DG%F,DG%dFda,DG%dgFda,DG%dgYda,DG%R,DG%dRda,DG%dRdgq,DG%Z,DG%dZda,DG%dZdgq,DG%BB,DG%dBBda,DG%dBBdgq,DG%P,DG%T,DG%a_2D,DG%gq_2D, Nq*Na)

   DG%Jacobian = DG%R*(DG%dZda*DG%dRdgq-DG%dZdgq*DG%dRda)
   DG%esi_N = (DG%dRda*DG%dRdgq+DG%dZda*DG%dZdgq)/DG%Jacobian
   DG%esi_K = (DG%dRdgq*DG%dRdgq+DG%dZdgq*DG%dZdgq)/DG%Jacobian

   DG%B0_cov3(1,:,:) = -DG%esi_N*DG%dgYda
   DG%B0_cov3(2,:,:) = -DG%esi_K*DG%dgYda
   DG%B0_cov3(3,:,:) = DG%F

   DG%B0_cntv3(1,:,:) = 0.0
   DG%B0_cntv3(2,:,:) = -DG%dgYda/DG%Jacobian
   DG%B0_cntv3(3,:,:) = DG%F/(DG%R**2)

   DG%j0_cntv3(1,:,:) = 0.0
   DG%j0_cntv3(2,:,:) = -DG%dFda/DG%Jacobian
   DG%j0_cntv3(3,:,:) = DG%P + DG%T/(DG%R**2)

   DG%g_aa_cov = DG%dRda*DG%dRda+DG%dZda*DG%dZda
   DG%g_aq_cov = DG%dRda*DG%dRdgq+DG%dZda*DG%dZdgq
   DG%g_qq_cov = DG%dRdgq*DG%dRdgq+DG%dZdgq*DG%dZdgq
   DG%g_ff_cov = DG%R**2

   DG%D=DG%dZda*DG%dRdgq-DG%dZdgq*DG%dRda

   DG%g_aa_cntv=(DG%dRdgq**2+DG%dZdgq**2)/DG%D**2   ! DG%grada^2
   DG%g_qq_cntv=(DG%dZda**2+DG%dRda**2)/DG%D**2   ! DG%gradtheta^2
   DG%g_ff_cntv=1.0/DG%R**2                 ! DG%gradphi^2
   DG%g_aq_cntv=(-DG%dRda*DG%dRdgq-DG%dZda*DG%dZdgq)/DG%D**2  ! DG%gradagradtheta


   ! Metrics regarding s coordinate
   do aid=1,DG%dim_a
      tmp=dsda(DG%a_1D(aid))
      DG%g_ss_cov(:,aid)   = DG%g_aa_cov(:,aid)/(tmp**2.0)
      DG%g_sq_cov(:,aid)   = DG%g_aq_cov(:,aid)/tmp

      DG%g_ss_cntv(:,aid)  = DG%g_aa_cntv(:,aid)*(tmp**2.0)
      DG%g_sq_cntv(:,aid)  = DG%g_aq_cntv(:,aid)*tmp

      DG%Jacob_sqf(:,aid) = DG%Jacobian(:,aid)/tmp

      DG%B0_cntv3_sqf(1,:,aid) = DG%B0_cntv3(1,:,aid)*tmp
   enddo
   DG%B0_cntv3_sqf(2:3,:,:)=DG%B0_cntv3(2:3,:,:)

   DG%j0_cov3(1,:,:) = DG%g_aa_cov*DG%j0_cntv3(1,:,:) + DG%g_aq_cov*DG%j0_cntv3(2,:,:)
   DG%j0_cov3(2,:,:) = DG%g_aq_cov*DG%j0_cntv3(1,:,:) + DG%g_qq_cov*DG%j0_cntv3(2,:,:)
   DG%j0_cov3(3,:,:) = DG%g_ff_cov*DG%j0_cntv3(3,:,:)


   ! NOTE : Get Cell volume like dV=(Jacobian)*dpsi*dtheta*dgf
   !        Following Cell_dV is an approximation
   !        For higher accuracy, a Weighted Volume by considering hellical particle shape and varying da is required
   approx_dgf=zetamax-zetamin
   do aid=1,Na
      if(aid==1) then
         approx_da=0.5*(DG%a_1D(2)-DG%a_1D(1))    ! half cell
      else if(aid==Na) then
         approx_da=0.5*(DG%a_1D(Na)-DG%a_1D(Na-1)) ! half cell
      else
         approx_da=0.5*(DG%a_1D(aid+1)-DG%a_1D(aid-1))  ! full cell
      endif

      do gqid=1,Nq
         DG%Cell_dV(gqid,aid)=DG%Jacobian(gqid,aid)*approx_da*DG%dgq*approx_dgf
         DG%Cell_dS_a(gqid,aid) = sqrt(DG%g_aa_cntv(gqid,aid))*DG%Jacobian(gqid,aid)*DG%dgq*approx_dgf
         DG%Cell_dS_gq(gqid,aid) = sqrt(DG%g_qq_cntv(gqid,aid))*DG%Jacobian(gqid,aid)*approx_da*approx_dgf
         DG%Cell_dS_gf(gqid,aid) = sqrt(DG%g_ff_cntv(gqid,aid))*DG%Jacobian(gqid,aid)*approx_da*DG%dgq
      enddo
   enddo
   DG%inv_Cell_dV=1.0/DG%Cell_dV
   DG%inv_Cell_dS_a  =1.0/DG%Cell_dS_a
   DG%inv_Cell_dS_gq =1.0/DG%Cell_dS_gq
   DG%inv_Cell_dS_gf =1.0/DG%Cell_dS_gf

   ! Set PVS (ProbeVelocitySpace) parameters

   nPVS = nprobe
   allocate(PVS(0:nPVS)) ! 0-index PVS will be used to measure particles in a whole poloidal section.

   ! Set location and size of probes
   PVS(0)%a  = 0.5*(a0+a1)
   PVS(0)%gq = 0
   PVS(0)%gf = 0.5*(zetamax+zetamin)
   PVS(0)%L_a  = a1-a0  !To deal with a whole poloidal section
   PVS(0)%L_gq = 2.0*pi
   PVS(0)%L_gf = zetamax-zetamin

   do pbid=1,nPVS
      PVS(pbid)%a  = rprobe(pbid)
      PVS(pbid)%gq = tprobe(pbid)
      PVS(pbid)%gf = 0.5*(zetamax+zetamin)
      PVS(pbid)%L_a  = dradius
      PVS(pbid)%L_gq = dtheta
      PVS(pbid)%L_gf = zetamax-zetamin
   enddo

   ! Set other information
   do pbid=0,nPVS
      PVS(pbid)%n0  = density(PVS(pbid)%a,0.0)
      PVS(pbid)%Ti0 = temperature(PVS(pbid)%a,PVS(pbid)%gq)
      PVS(pbid)%Te0 = msm0*temperature_e(PVS(pbid)%a,PVS(pbid)%gq)

      PVS(pbid)%inv_n0  = 1.0/PVS(pbid)%n0
      PVS(pbid)%inv_Ti0 = 1.0/PVS(pbid)%Ti0
      PVS(pbid)%inv_Te0 = 1.0/PVS(pbid)%Te0

      PVS(pbid)%nv  = nv
      PVS(pbid)%dv  = delta_v
      PVS(pbid)%inv_dv = 1.0/PVS(pbid)%dv
      PVS(pbid)%max_norm_v = PVS(pbid)%nv*PVS(pbid)%dv ! or just use umax parameters (ex. 1.2*umax)

      PVS(pbid)%dim_vpara = 2*PVS(pbid)%nv+1
      PVS(pbid)%dim_vperp = PVS(pbid)%nv+1

      PVS(pbid)%num_of_fi  = 5
      PVS(pbid)%num_of_fe  = 5

      allocate(PVS(pbid)%vpara_1D(-PVS(pbid)%nv:PVS(pbid)%nv))
      allocate(PVS(pbid)%vperp_1D(0:PVS(pbid)%nv))
      allocate(PVS(pbid)%fi(-PVS(pbid)%nv:PVS(pbid)%nv, 0:PVS(pbid)%nv, PVS(pbid)%num_of_fi))
      allocate(PVS(pbid)%fe(-PVS(pbid)%nv:PVS(pbid)%nv, 0:PVS(pbid)%nv, PVS(pbid)%num_of_fe))

      do i=0,PVS(pbid)%nv
         PVS(pbid)%vperp_1D(i) = i*PVS(pbid)%dv
      enddo
      do i=-PVS(pbid)%nv,PVS(pbid)%nv
         PVS(pbid)%vpara_1D(i) = i*PVS(pbid)%dv
      enddo
   enddo

   if(myrank_partd==nproc_partd-1) then

      call Cal_FLA_Denominator_on_DG_using_gPPI

      allocate(tmp_DG_Data1(Nq, Na,0:1),tmp_DG_Data2(Nq, Na,0:1))

      DG%FLA_SM_Denominator=DG%FLA_Denominator(:,:,1)
      call Smooth_DG_Data_2D(3, DG%FLA_SM_Denominator)

      call Do_FieldLineAverage_on_DG_by_using_gPPI([DG%a_2D,DG%a_2D], DG%FLA_psi0)
      DG%FLA_SM_psi0 = DG%FLA_psi0(:,:,1)
      call Smooth_DG_Data_2D(3, DG%FLA_SM_psi0)

      ! Cal FLabel using repetitive FieldLineAverage
      call Cal_FLA_FLabel_on_DG_by_using_gPPI_iniFLabel(DG%FLA_FLabel)
      tmp_DG_Data1 = DG%FLA_FLabel
      do i=1,4
         call Smooth_DG_Data_2D(2, tmp_DG_Data1(:,:,0))
         call Smooth_DG_Data_2D(2, tmp_DG_Data1(:,:,1))
         call Do_FieldLineAverage_on_DG_by_using_gPPI(tmp_DG_Data1, tmp_DG_Data2)
         tmp_DG_Data1=tmp_DG_Data2
      enddo
      ! DG%FLA_FLabel = tmp_DG_Data1
      DG%FLA_SM_FLabel = tmp_DG_Data1(:,:,1)
      call Smooth_DG_Data_2D(1, DG%FLA_SM_FLabel)

      deallocate(tmp_DG_Data1,tmp_DG_Data2)
   endif

   ! Broadcast calculated FLA data
   call MPI_BCAST(DG%FLA_Denominator, 2*Na*Nq,mpi_Rsize, nproc_partd-1,partd_comm,ierror)
   call MPI_BCAST(DG%FLA_FLabel,      2*Na*Nq,mpi_Rsize, nproc_partd-1,partd_comm,ierror)
   call MPI_BCAST(DG%FLA_psi0,        2*Na*Nq,mpi_Rsize, nproc_partd-1,partd_comm,ierror)
   call MPI_BCAST(DG%FLA_SM_Denominator, Na*Nq,mpi_Rsize, nproc_partd-1,partd_comm,ierror)
   call MPI_BCAST(DG%FLA_SM_FLabel,      Na*Nq,mpi_Rsize, nproc_partd-1,partd_comm,ierror)
   call MPI_BCAST(DG%FLA_SM_psi0,        Na*Nq,mpi_Rsize, nproc_partd-1,partd_comm,ierror)

   ! Communicate btw toroidal_comm to get Left and Right FLA_SM_FLabel
   call MPI_SENDRECV(DG%FLA_SM_FLabel, Nq*Na, mpi_Rsize, right_pe, myrank_toroidal, &
   DG%FLA_SM_FLabel_Left, Nq*Na, mpi_Rsize, left_pe, left_pe, toroidal_comm, istatus, ierror)
   call MPI_SENDRECV(DG%FLA_SM_FLabel, Nq*Na, mpi_Rsize, left_pe, myrank_toroidal, &
   DG%FLA_SM_FLabel_Right, Nq*Na, mpi_Rsize, right_pe, right_pe, toroidal_comm, istatus, ierror)


   !Calculate normal vector (gradient) of FLA_SM_FLabel
   Na=Na
   Nq=Nq
   ! d(FL)/da
   DG%grad_FLA_SM_FLabel(1,:,1,1) = (DG%FLA_SM_FLabel(:,2)-DG%FLA_SM_FLabel(:,1))/(DG%a_1D(2)-DG%a_1D(1))
   DG%grad_FLA_SM_FLabel(1,:,2:Na-1,1) = (DG%FLA_SM_FLabel(:,3:Na)-DG%FLA_SM_FLabel(:,1:Na-2))/(DG%a_2D(:,3:Na)-DG%a_2D(:,1:Na-2))
   DG%grad_FLA_SM_FLabel(1,:,Na,1) = (DG%FLA_SM_FLabel(:,Na)-DG%FLA_SM_FLabel(:,Na-1))/(DG%a_1D(Na)-DG%a_1D(Na-1))

   ! d(FL)/dgq with periodicity
   DG%grad_FLA_SM_FLabel(2,1,:,1) = (DG%FLA_SM_FLabel(2,:)-DG%FLA_SM_FLabel(Nq,:))/(2.0*DG%dgq)
   DG%grad_FLA_SM_FLabel(2,2:Nq-1,:,1) = (DG%FLA_SM_FLabel(3:Nq,:)-DG%FLA_SM_FLabel(1:Nq-2,:))/(2.0*DG%dgq)
   DG%grad_FLA_SM_FLabel(2,Nq,:,1) = DG%grad_FLA_SM_FLabel(2,1,:,1)

   ! d(FL)/dgf with periodicity
   DG%grad_FLA_SM_FLabel(3,:,:,1) = (DG%FLA_SM_FLabel_Right - DG%FLA_SM_FLabel_Left)/(2.0*(zetamax-zetamin))

   call MPI_SENDRECV(DG%grad_FLA_SM_FLabel(:,:,:,1), 3*Nq*Na, mpi_Rsize, right_pe, myrank_toroidal, &
                     DG%grad_FLA_SM_FLabel(:,:,:,0), 3*Nq*Na, mpi_Rsize, left_pe, left_pe, toroidal_comm, istatus, ierror)



   ! Calculate magnitude of grad(FL)
   do k=0,1
      DG%mag_grad_FLA_SM_FLabel(:,:,k) = sqrt(  DG%g_aa_cntv*(DG%grad_FLA_SM_FLabel(1,:,:,k)**2.0) &
                                                +2.0*DG%g_aq_cntv*DG%grad_FLA_SM_FLabel(1,:,:,k)*DG%grad_FLA_SM_FLabel(2,:,:,k) &
                                                +DG%g_qq_cntv*(DG%grad_FLA_SM_FLabel(2,:,:,k)**2.0)  &
                                                +DG%g_ff_cntv*(DG%grad_FLA_SM_FLabel(3,:,:,k)**2.0))
   enddo

   do k=1,3
      DG%uvec_grad_FL_cov(k,:,:,:) = DG%grad_FLA_SM_FLabel(k,:,:,:)/(DG%mag_grad_FLA_SM_FLabel)
   enddo

end subroutine Init_Diagnostic_Grid


subroutine Cal_unit_vectors_tangential_normal_binormal_to_FluxSurface(in_a, in_gq, in_gf, in_s, uvec_B_cov, uvec_B_cntv, uvec_N_cov, uvec_N_cntv, uvec_BxN_cov, uvec_BxN_cntv)
   use precision
   use Diagnostic_Grid
   use global_parameters
   use field_array
   use esi_interface

   implicit none
   real(wp),intent(in) :: in_a, in_gq, in_gf, in_s

   real(wp) target_gq ! projected gq at nearby poloidal plane following the B0 field
   integer :: ori_sid !sid of original full grid
   real(wp) ori_inv_ds ! inv_ds of original full grid
   real(wp) ff, ws0,ws1, wgq0,wgq1, wgf,ratio ! ratios along each dimension
   real(wp) target_zeta, safety_q
   integer :: j,k, sid, gqid, qq, ss

   real(wp):: F_,dFda_,dgFda_,dgYda_,R_,dRda_,dRdgq_,Z_,dZda_,dZdgq_,BB_,dBBda_,dBBdgq_,P_,T_
   real(wp) Jacobian_, D_
   real(wp):: g_aa_cov_, g_aq_cov_, g_qq_cov_, g_ff_cov_, g_aa_cntv_, g_aq_cntv_, g_qq_cntv_, g_ff_cntv_

   real(wp) tot_B_vec_cov(3), tot_B_vec_cntv(3)
   real(wp) vec_BxN_cntv(3), vec_BxN_cov(3) !Note tha theses are not the unit vectors
   real(wp), intent(out) :: uvec_B_cov(3), uvec_N_cov(3), uvec_BxN_cov(3) ! unit vectors in total B direction, Normal (gradFL), and binormal
   real(wp), intent(out) :: uvec_B_cntv(3), uvec_N_cntv(3), uvec_BxN_cntv(3)

   ori_inv_ds = real(mpsi)/(s1-s0)

   !Calculate Interpolation of uvec_N_cov which uses the precalculated DG%uvec_grad_FL_cov on DG nodes
   !Get uvec_B_cov(3), uvec_N_cov(3), uvec_BxN_cov(3)
   call esiget_Intp_B_vec_cov_and_cntv(in_a, in_gq, in_gf, istep*tstep, tot_B_vec_cov, tot_B_vec_cntv)
   uvec_B_cov   = tot_B_vec_cov/sqrt(dot_product(tot_B_vec_cov,tot_B_vec_cntv))

   call esiget2dthese_more(F_,dFda_,dgFda_,dgYda_,R_,dRda_,dRdgq_,Z_,dZda_,dZdgq_,BB_,dBBda_,dBBdgq_,P_,T_,in_a,in_gq, 1)
   Jacobian_ = R_*(dZda_*dRdgq_-dZdgq_*dRda_)
   D_=dZda_*dRdgq_-dZdgq_*dRda_

   g_aa_cov_   = dRda_*dRda_+dZda_*dZda_;
   g_aq_cov_  = dRda_*dRdgq_+dZda_*dZdgq_;
   g_qq_cov_ = dRdgq_*dRdgq_+dZdgq_*dZdgq_;
   g_ff_cov_ = R_**2

   g_aa_cntv_=(dRdgq_**2+dZdgq_**2)/D_**2   ! grada^2
   g_qq_cntv_=(dZda_**2+dRda_**2)/D_**2   ! gradtheta^2
   g_ff_cntv_=1.0_wp/R_**2                 ! gradphi^2
   g_aq_cntv_=(-dRda_*dRdgq_-dZda_*dZdgq_)/D_**2  ! gradagradtheta


   ori_sid = (in_s-s0)*ori_inv_ds ! s index of original full grid
   ff  = (in_s-s_grid(ori_sid))*ori_inv_ds
   safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)
   sid = max(0,min(int((in_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
   ws1 = (in_s-DG%s_1D(sid))*DG%inv_ds
   ws0 = 1.0-ws1

   uvec_N_cov = 0
   do k=0,1
      if(k==0) target_zeta =zetamin
      if(k==1) target_zeta =zetamax

      wgf = (1.0-abs(in_gf-target_zeta)/deltaz)

      target_gq = modulo(in_gq + (target_zeta-in_gf)/safety_q, 2.0*pi)
      gqid = (target_gq)*DG%inv_dgq +1
      wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
      wgq0 = 1.0-wgq1
      do j=1,4
         select case(j)
            case(1)
               qq=gqid; ss=sid
               ratio=ws0*wgq0*wgf
            case(2)
               qq=gqid+1; ss=sid
               ratio=ws0*wgq1*wgf
            case(3)
               qq=gqid; ss=sid+1
               ratio=ws1*wgq0*wgf
            case(4)
               qq=gqid+1; ss=sid+1
               ratio=ws1*wgq1*wgf
         end select
         uvec_N_cov = uvec_N_cov+ratio*DG%uvec_grad_FL_cov(:,qq,ss,k)
      enddo
   enddo

   ! Remove a possible small error part of uvec_N_cov which is parallel to the magnetic field
   uvec_B_cntv = tot_B_vec_cntv/sqrt(dot_product(tot_B_vec_cov,tot_B_vec_cntv))
   uvec_N_cov= uvec_N_cov - dot_product(uvec_N_cov,uvec_B_cntv)*uvec_B_cov

   uvec_N_cov = uvec_N_cov/sqrt( g_aa_cntv_*uvec_N_cov(1)**2.0 + 2.0*g_aq_cntv_*uvec_N_cov(1)*uvec_N_cov(2) + g_qq_cntv_*uvec_N_cov(2)**2.0 + g_ff_cntv_*uvec_N_cov(3)**2.0)
   uvec_N_cntv(1) = g_aa_cntv_*uvec_N_cov(1) + g_aq_cntv_*uvec_N_cov(2)
   uvec_N_cntv(2) = g_aq_cntv_*uvec_N_cov(1) + g_qq_cntv_*uvec_N_cov(2)
   uvec_N_cntv(3) = g_ff_cntv_*uvec_N_cov(3)

   vec_BxN_cntv(1) = (1.0/Jacobian_)*(uvec_B_cov(2)*uvec_N_cov(3)-uvec_B_cov(3)*uvec_N_cov(2))
   vec_BxN_cntv(2) = (1.0/Jacobian_)*(uvec_B_cov(3)*uvec_N_cov(1)-uvec_B_cov(1)*uvec_N_cov(3))
   vec_BxN_cntv(3) = (1.0/Jacobian_)*(uvec_B_cov(1)*uvec_N_cov(2)-uvec_B_cov(2)*uvec_N_cov(1))

   vec_BxN_cov(1) = g_aa_cov_*vec_BxN_cntv(1) + g_aq_cov_*vec_BxN_cntv(2)
   vec_BxN_cov(2) = g_aq_cov_*vec_BxN_cntv(1) + g_qq_cov_*vec_BxN_cntv(2)
   vec_BxN_cov(3) = g_ff_cov_*vec_BxN_cntv(3)

   uvec_BxN_cov   = vec_BxN_cov/sqrt(dot_product(vec_BxN_cov,vec_BxN_cntv))
   uvec_BxN_cntv  = vec_BxN_cntv/sqrt(dot_product(vec_BxN_cov,vec_BxN_cntv))

end subroutine Cal_unit_vectors_tangential_normal_binormal_to_FluxSurface


subroutine Measure_Electrons_on_Diagnostic_Grid(in_P_data, Full_dadgqdgf)

   use precision
   use Diagnostic_Grid
   use global_parameters
   use field_array
   use mpi_defs
   use particle_array
   use particle_decomp
   use esi_interface

   implicit none

   real(wp), intent(in) :: in_P_data(3,me) ! additional particle data from input
   real(wp), intent(in) :: Full_dadgqdgf(3,me)
   real(wp) energy(me)
   !P varaiables are at particle position
   real(wp) P_s, P_safety_q !safety factor
   real(wp) P_a, P_gq, P_gf

   real(wp) delta_f(me), f_tot(me) ! delta_f and f_tot (f_tot=f0+delta_f)

   real(wp) target_gq ! projected gq at nearby poloidal plane following the B0 field
   integer :: ori_sid !sid of original full grid
   real(wp) ori_inv_ds,ori_inv_dgf, P_inv_ori_dgq ! inv_ds of original full grid
   real(wp) ff, ws0,ws1, wgq0,wgq1, wgf,ratio ! ratios along each dimension
   real(wp) target_zeta
   integer :: i,m, j,k, aid, sid, vid, gqid, qq, ss
   integer ::pbid,vxid,vyid,fid
   real(wp):: norm_vpara, norm_vperp, rvx,rvy ! norm_vpara=vpara/vth, norm_vperp=vperp/vth, where vth=sqrt(T/m)
   real(wp),allocatable, dimension(:,:,:) :: tmp_PVS_fe
   integer :: length, istatus(MPI_STATUS_SIZE),ierror

   real(wp) tmp_Scalar_Data(DG%dim_gq,DG%dim_a,DG%num_eSD,0:1), tmp_Vector_Data(3,DG%dim_gq,DG%dim_a,DG%num_eVD,0:1)
   real(wp) :: intpl_dphi(DG%dim_gq, DG%dim_a), inv_tot_ne(DG%dim_gq, DG%dim_a)

   real(wp) :: uvec_N_cov(3), uvec_B_cov(3), uvec_BxN_cov(3) ! unit vectors in total B direction, Normal (gradFL), and binormal
   real(wp) :: uvec_N_cntv(3), uvec_B_cntv(3), uvec_BxN_cntv(3) ! unit vectors in total B direction, Normal (gradFL), and binormal

   DG%eScaData=0
   DG%eVecData=0
   do pbid=0,nPVS
      PVS(pbid)%fe = 0
   enddo

   ori_inv_ds = real(mpsi)/(s1-s0)
   ori_inv_dgf= 1.0/(zetamax-zetamin)

   energy(:) = in_P_data(2,:)

   ! Calculate dne_adiabatic by interpolating dphi on original grids
   call Interpolate_OriGridData_into_DG(phi(1,:), intpl_dphi)
   do aid=1,DG%dim_a
      i = DG%ori_aid(aid)
      DG%dne_adiab(:,aid)=den_grid(i)*abs(csc0)*(intpl_dphi(:,aid)-real(1.0-ZF_SEPARATION)*phi00(i))/Te_grid(i)
   enddo


   if(Flag_Load_and_WeightDef==0)then
      delta_f(:) = zeon0(5,1:me)*eNrNm
      f_tot(:)   = eNrNm
   else
      delta_f(:) = zeon0(5,1:me)
      f_tot(:)   = zeon0(6,1:me)
   endif

   do m=1,me

      if (Flag_w_Threshold==1 .and. (zeon0(5,m)<w_Thold_min .or. zeon0(5,m)>w_Thold_max)) cycle

      P_a  = zeon0(1,m)
      P_gq = modulo(zeon0(2,m),2*pi)
      P_gf = modulo(zeon0(3,m),2*pi)
      P_s  = zeon0(7,m)

      ori_sid = (P_s-s0)*ori_inv_ds ! s index of original full grid
      ff  = (P_s-s_grid(ori_sid))*ori_inv_ds
      P_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)
      P_inv_ori_dgq = 1.0/(ff*deltat(ori_sid+1)+(1.0-ff)*deltat(ori_sid))

      sid = max(0,min(int((P_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
      ws1 = (P_s-DG%s_1D(sid))*DG%inv_ds
      ws0 = 1.0-ws1

      call Cal_unit_vectors_tangential_normal_binormal_to_FluxSurface(P_a,P_gq,P_gf,P_s, uvec_B_cov, uvec_B_cntv, uvec_N_cov, uvec_N_cntv, uvec_BxN_cov, uvec_BxN_cntv)

      ! Scatter data into neighboring DG nodes
      do k=0,1
         if(k==0) target_zeta =zetamin
         if(k==1) target_zeta =zetamax

         wgf = (1.0-abs(P_gf-target_zeta)/deltaz)

         target_gq = modulo(P_gq + (target_zeta-P_gf)/P_safety_q, 2.0*pi)
         gqid = (target_gq)*DG%inv_dgq +1
         wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
         wgq0 = 1.0-wgq1

         do j=1,4
            select case(j)
               case(1)
                  qq=gqid; ss=sid
                  ratio=ws0*wgq0*wgf
               case(2)
                  qq=gqid+1; ss=sid
                  ratio=ws0*wgq1*wgf
               case(3)
                  qq=gqid; ss=sid+1
                  ratio=ws1*wgq0*wgf
               case(4)
                  qq=gqid+1; ss=sid+1
                  ratio=ws1*wgq1*wgf
            end select

            vid=1; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m) ! Density
            vid=2; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*(2.0/3.0*energy(m)-backTe(m)*msm0)*delta_f(m) !Temperature definition 1
            vid=3; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*(2.0/3.0*energy(m)-backTe(m)*msm0)*delta_f(m) !
            vid=4; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio ! Marker Distribution
            vid=5; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*f_tot(m)*(2.0/3.0*energy(m)) ! Integrating f_tot*energy, so Te_tot

            ! Volumetric delta-f Partcle Flux (=df*v_tot) (#/(m^2*s)) along tangential, normal, binormal directions, respectively
            vid=7; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_B_cov)
            vid=8; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_N_cov)
            vid=9; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_BxN_cov)

            vid=10; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_B_cntv(1))*uvec_B_cov)
            vid=11; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_N_cntv(1))*uvec_N_cov)
            vid=12; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_BxN_cntv(1))*uvec_BxN_cov)

            ! Volumetric delta-f Energy Flux (=df*KE*v_tot) (#/(m^2*s)) along tangential, normal, binormal directions, respectively
            vid=13; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_B_cov)
            vid=14; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_N_cov)
            vid=15; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_BxN_cov)

            vid=16; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_B_cntv(1))*uvec_B_cov)
            vid=17; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_N_cntv(1))*uvec_N_cov)
            vid=18; DG%eScaData(qq,ss,vid,k) = DG%eScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_BxN_cntv(1))*uvec_BxN_cov)


            vid=1 ! Turbulent Particle Flux
            ! DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*(delta_f(m)*ghq(m)) ! Radial Flux
            ! DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*(delta_f(m)*Ta(m)) ! theta Flux
            ! DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*(delta_f(m)*dve_gf(m)) ! zeta Flux
            DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*delta_f(m)*Full_dadgqdgf(1,m) ! Radial Flux
            DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*delta_f(m)*Full_dadgqdgf(2,m) ! theta Flux
            DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*delta_f(m)*Full_dadgqdgf(3,m) ! zeta Flux


            vid=2 ! Turbulent Particle Flux due to ExB only
            DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*(delta_f(m)*(ghq(m)-vpara_dve_a(m))) ! Radial Flux
            DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*(delta_f(m)*(Ta(m)-vpara_dve_gq(m))) ! theta Flux
            DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*(delta_f(m)*(dve_gf(m)-vpara_dve_gf(m)) )! zeta Flux

            vid=3 ! Turbulent Energy Flux
            ! DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*ghq(m)) ! Radial Flux
            ! DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Ta(m)) ! theta Flux
            ! DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*dve_gf(m)) ! zeta Flux
            DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Full_dadgqdgf(1,m)) ! Radial Flux
            DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Full_dadgqdgf(2,m)) ! theta Flux
            DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Full_dadgqdgf(3,m)) ! zeta Flux

            vid=4 ! Turbulent Energy Flux due to ExB only
            DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*(ghq(m)-vpara_dve_a(m))) ! Radial Flux
            DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*(Ta(m)-vpara_dve_gq(m))) ! theta Flux
            DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*(dve_gf(m)-vpara_dve_gf(m)) )! zeta Flux

            vid=5 ! To measure CFL condition (v*dt/dx)
            DG%eVecData(1,qq,ss,vid,k) = DG%eVecData(1,qq,ss,vid,k) + ratio*abs(tstep*((Tc/backT(m))**dsda_power*Full_dadgqdgf(1,m))*ori_inv_ds)
            DG%eVecData(2,qq,ss,vid,k) = DG%eVecData(2,qq,ss,vid,k) + ratio*abs(tstep*(Full_dadgqdgf(2,m)-Full_dadgqdgf(3,m)/P_safety_q)*P_inv_ori_dgq) !time derivative of target_gq=gq-gf/safety_q
            DG%eVecData(3,qq,ss,vid,k) = DG%eVecData(3,qq,ss,vid,k) + ratio*abs(tstep*Full_dadgqdgf(3,m)*ori_inv_dgf)
         enddo
      enddo


      ! Velocity Space Measurement
      do pbid=0,nPVS
         target_gq = modulo(P_gq + (PVS(pbid)%gf-P_gf)/P_safety_q, 2.0*pi)

         if( abs(P_gf-PVS(pbid)%gf)<0.5*PVS(pbid)%L_gf &
         .and. abs(P_a-PVS(pbid)%a) <0.5*PVS(pbid)%L_a &
         .and. min(abs(target_gq-PVS(pbid)%gq), 2.0*pi-abs(target_gq-PVS(pbid)%gq)) < 0.5*PVS(pbid)%L_gq) then

            norm_vpara = zeon0(4,m)*csc0*BB(m)*sqrt(m0ms*PVS(pbid)%inv_Te0) ! norm_vpara = vpara/sqrt(T/m) = rho_para*(eB/m)/sqrt(T/m) = rho_para*(e*B)/sqrt(m*T)
            norm_vperp = zeon(6,m)*sqrt(2.0*BB(m))*sqrt(PVS(pbid)%inv_Te0) ! norm_vperp = sqrt(mu)*sqrt(2B/m)/sqrt(T/m) = sqrt(mu)*sqrt(2B)/(sqrt(T))

            vxid = floor(norm_vpara*PVS(pbid)%inv_dv)
            vyid = floor(norm_vperp*PVS(pbid)%inv_dv)
            rvx  = (norm_vpara-PVS(pbid)%vpara_1D(vxid))*PVS(pbid)%inv_dv
            rvy  = (norm_vperp-PVS(pbid)%vperp_1D(vyid))*PVS(pbid)%inv_dv

            if(vxid< -PVS(pbid)%nv) then
               vxid= -PVS(pbid)%nv
               rvx=0.0
            elseif (vxid>= PVS(pbid)%nv) then
               vxid=  PVS(pbid)%nv-1
               rvx=1.0
            endif

            if(vyid>= PVS(pbid)%nv) then
               vyid=  PVS(pbid)%nv-1
               rvy=1.0
            endif

            do j=1,4
               select case(j)
                  case(1)
                     qq=vxid; ss=vyid
                     ratio=(1.0-rvx)*(1.0-rvy)
                  case(2)
                     qq=vxid+1; ss=vyid
                     ratio=(rvx)*(1.0-rvy)
                  case(3)
                     qq=vxid; ss=vyid+1
                     ratio=(1.0-rvx)*(rvy)
                  case(4)
                     qq=vxid+1; ss=vyid+1
                     ratio=(rvx)*(rvy)
               end select

               fid=1; PVS(pbid)%fe(qq,ss,fid) = PVS(pbid)%fe(qq,ss,fid) + ratio*delta_f(m) ! Density
               fid=2; PVS(pbid)%fe(qq,ss,fid) = PVS(pbid)%fe(qq,ss,fid) + ratio*delta_f(m)*(2.0/3.0*energy(m)-backT(m)) ! Temperature definition
               fid=3; PVS(pbid)%fe(qq,ss,fid) = PVS(pbid)%fe(qq,ss,fid) + ratio ! Marker Density
               fid=4; PVS(pbid)%fe(qq,ss,fid) = PVS(pbid)%fe(qq,ss,fid) + ratio*f_tot(m) ! f_tot
               fid=5; PVS(pbid)%fe(qq,ss,fid) = PVS(pbid)%fe(qq,ss,fid) + ratio*f_tot(m)*(2.0/3.0*energy(m)) ! f_tot*energy
            enddo
         endif
      enddo

   enddo


   if(npartdom>1) then
      tmp_Scalar_Data = DG%eScaData
      tmp_Vector_Data = DG%eVecData

      length = DG%dim_gq*DG%dim_a*DG%num_eSD*2
      call MPI_ALLREDUCE(tmp_Scalar_Data, DG%eScaData, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)

      length = 3*DG%dim_gq*DG%dim_a*DG%num_eVD*2
      call MPI_ALLREDUCE(tmp_Vector_Data, DG%eVecData, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)

      allocate(tmp_PVS_fe, mold=PVS(0)%fe)
      do pbid=0,nPVS
         tmp_PVS_fe = PVS(pbid)%fe
         length = PVS(pbid)%dim_vpara*PVS(pbid)%dim_vperp*PVS(pbid)%num_of_fe
         call MPI_ALLREDUCE(tmp_PVS_fe, PVS(pbid)%fe, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      enddo
      deallocate(tmp_PVS_fe)
   endif


   ! Match poloidal end cell
   ! (gq=0 and gq=2pi is the same position physically
   ! although they have different gqid)
   do sid=1,DG%dim_a
      DG%eScaData(1,sid,:,:) = DG%eScaData(1,sid,:,:) + DG%eScaData(DG%dim_gq,sid,:,:)
      DG%eScaData(DG%dim_gq,sid,:,:) = DG%eScaData(1,sid,:,:)

      DG%eVecData(:,1,sid,:,:) = DG%eVecData(:,1,sid,:,:) + DG%eVecData(:,DG%dim_gq,sid,:,:)
      DG%eVecData(:,DG%dim_gq,sid,:,:) = DG%eVecData(:,1,sid,:,:)
   enddo


   ! toroidal end cell
   length = DG%dim_gq*DG%dim_a*DG%num_eSD
   call MPI_SENDRECV(DG%eScaData(:,:,:,0),length,mpi_Rsize,left_pe,1, tmp_Scalar_Data(:,:,:,1),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   DG%eScaData(:,:,:,1) = DG%eScaData(:,:,:,1) + tmp_Scalar_Data(:,:,:,1)



   length = 3*DG%dim_gq*DG%dim_a*DG%num_eVD
   call MPI_SENDRECV(DG%eVecData(:,:,:,:,0),length,mpi_Rsize,left_pe,1, tmp_Vector_Data(:,:,:,:,1),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   DG%eVecData(:,:,:,:,1) = DG%eVecData(:,:,:,:,1) + tmp_Vector_Data(:,:,:,:,1)


   ! Do normalization and post-processing
   ! Divide by Cell Volume to get density quantities
   do vid=1,DG%num_eSD
      DG%eScaData(:,:,vid,1) = DG%eScaData(:,:,vid,1)*DG%inv_Cell_dV
   enddo

   do aid=1,DG%dim_a
      DG%eScaData(:,aid,2,1) = DG%eScaData(:,aid,2,1)/DG%n0_1D(aid) ! Convert Energy to Temperature
   enddo


   ! Calculate inverse of tot density (including dne_adiabatic in the case of delta-h electron (ifelectron==2))
   do aid=1,DG%dim_a
      inv_tot_ne(:,aid) = 1.0/(DG%n0_1D(aid) + DG%eScaData(:,aid,1,1) + int(ifelectron/2)*DG%dne_adiab(:,aid))
   enddo

   ! Consistent Temperature by considering density fluctuation in the denominator
   DG%eScaData(:,:,3,1) = DG%eScaData(:,:,3,1)*inv_tot_ne
   vid=6; DG%eScaData(:,:,vid,1) = DG%eScaData(:,:,vid,1)*inv_tot_ne! Te_tot


   ! Get Flux Density by dividing Area in each direction
   do vid=1,4
      j=1; DG%eVecData(j,:,:,vid,1) = DG%eVecData(j,:,:,vid,1)*DG%inv_Cell_dS_a  ! Radial Flux Density
      j=2; DG%eVecData(j,:,:,vid,1) = DG%eVecData(j,:,:,vid,1)*DG%inv_Cell_dS_gq ! Theta Flux Density
      j=3; DG%eVecData(j,:,:,vid,1) = DG%eVecData(j,:,:,vid,1)*DG%inv_Cell_dS_gf ! Zeta Flux Density
   enddo

   ! To Normalize CFL condtion values, divide by number of markers (4th iSca)
   vid=5
   do j=1,3
      DG%eVecData(j,:,:,vid,1) = DG%eVecData(j,:,:,vid,1)/(DG%eScaData(:,:,4,1)*DG%Cell_dV)
   enddo

end subroutine Measure_Electrons_on_Diagnostic_Grid



subroutine Measure_Ions_on_Diagnostic_Grid(in_P_data, Full_dadgqdgf)

   use precision
   use Diagnostic_Grid
   use global_parameters
   use field_array
   use mpi_defs
   use particle_array
   use particle_decomp
   use esi_interface

   implicit none

   real(wp), intent(in) :: in_P_data(3,mi) ! additional particle data from input
   real(wp), intent(in) :: Full_dadgqdgf(3,mi)
   real(wp) energy(mi)
   !P varaiables are at particle position
   real(wp) P_s, P_safety_q !safety factor
   real(wp) P_a, P_gq, P_gf

   real(wp) delta_f(mi), f_tot(mi) ! delta_f and f_tot (f_tot=f0+delta_f)

   real(wp) target_gq ! projected gq at nearby poloidal plane following the B0 field
   integer :: ori_sid !sid of original full grid
   real(wp) ori_inv_ds,ori_inv_dgf, P_inv_ori_dgq ! inv_ds of original full grid
   real(wp) ff, ws0,ws1, wgq0,wgq1, wgf,ratio ! ratios along each dimension
   real(wp) target_zeta
   integer :: i,m, j,k, aid, sid, vid, gqid, qq, ss
   integer ::pbid,vxid,vyid,fid
   real(wp):: norm_vpara, norm_vperp, rvx,rvy ! norm_vpara=vpara/vth, norm_vperp=vperp/vth, where vth=sqrt(T/m)
   real(wp),allocatable, dimension(:,:,:) :: tmp_PVS_fi
   integer :: length, istatus(MPI_STATUS_SIZE),ierror
   integer :: larmor
   real(wp) tmp_Scalar_Data(DG%dim_gq,DG%dim_a,DG%num_iSD,0:1), tmp_Vector_Data(3,DG%dim_gq,DG%dim_a,DG%num_iVD,0:1)
   real(wp) ::inv_tot_ni(DG%dim_gq, DG%dim_a)

   real(wp) :: uvec_N_cov(3), uvec_B_cov(3), uvec_BxN_cov(3) ! unit vectors in total B direction, Normal (gradFL), and binormal
   real(wp) :: uvec_N_cntv(3), uvec_B_cntv(3), uvec_BxN_cntv(3) ! unit vectors in total B direction, Normal (gradFL), and binormal

   DG%iScaData=0
   DG%iVecData=0
   do pbid=0,nPVS
      PVS(pbid)%fi = 0
   enddo


   ori_inv_ds = real(mpsi)/(s1-s0)
   ori_inv_dgf= 1.0/(zetamax-zetamin)

   energy(:) = in_P_data(2,:)

   if(Flag_Load_and_WeightDef==0)then
      delta_f(:) = zion0(5,1:mi)*iNrNm
      f_tot(:)   = iNrNm
   else
      delta_f(:) = zion0(5,1:mi)
      f_tot(:)   = zion0(6,1:mi)
   endif

   do m=1,mi

      if (Flag_w_Threshold==1 .and. (zion0(5,m)<w_Thold_min .or. zion0(5,m)>w_Thold_max)) cycle

      P_a  = zion0(1,m)
      P_gq = modulo(zion0(2,m),2*pi)
      P_gf = modulo(zion0(3,m),2*pi)
      P_s  = zion0(7,m)

      ori_sid = (P_s-s0)*ori_inv_ds ! s index of original full grid
      ff  = (P_s-s_grid(ori_sid))*ori_inv_ds
      P_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)
      P_inv_ori_dgq = 1.0/(ff*deltat(ori_sid+1)+(1.0-ff)*deltat(ori_sid))

      sid = max(0,min(int((P_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
      ws1 = (P_s-DG%s_1D(sid))*DG%inv_ds
      ws0 = 1.0-ws1

      call Cal_unit_vectors_tangential_normal_binormal_to_FluxSurface(P_a,P_gq,P_gf,P_s, uvec_B_cov, uvec_B_cntv, uvec_N_cov, uvec_N_cntv, uvec_BxN_cov, uvec_BxN_cntv)

      ! Scatter data into neighboring DG nodes
      do k=0,1
         if(k==0) target_zeta =zetamin
         if(k==1) target_zeta =zetamax

         wgf = (1.0-abs(P_gf-target_zeta)/deltaz)

         target_gq = modulo(P_gq + (target_zeta-P_gf)/P_safety_q, 2.0*pi)
         gqid = (target_gq)*DG%inv_dgq +1
         wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
         wgq0 = 1.0-wgq1

         do j=1,4
            select case(j)
               case(1)
                  qq=gqid; ss=sid
                  ratio=ws0*wgq0*wgf
               case(2)
                  qq=gqid+1; ss=sid
                  ratio=ws0*wgq1*wgf
               case(3)
                  qq=gqid; ss=sid+1
                  ratio=ws1*wgq0*wgf
               case(4)
                  qq=gqid+1; ss=sid+1
                  ratio=ws1*wgq1*wgf
            end select

            vid=1; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m) ! Density
            vid=2; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*(2.0/3.0*energy(m)-backT(m))*delta_f(m) !Temperature definition 1
            vid=3; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*(2.0/3.0*energy(m)-backT(m))*delta_f(m) !
            vid=4; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio ! Marker Distribution
            vid=5; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*f_tot(m)*(2.0/3.0*energy(m)) ! Integrating f_tot*energy, so Te_tot

            ! Volumetric delta-f Partcle Flux (=df*v_tot) (#/(m^2*s)) along tangential, normal, binormal directions, respectively
            vid=9; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_B_cov)
            vid=10;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_N_cov)
            vid=11;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_BxN_cov)

            vid=12;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_B_cntv(1))*uvec_B_cov)
            vid=13;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_N_cntv(1))*uvec_N_cov)
            vid=14;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_BxN_cntv(1))*uvec_BxN_cov)

            ! Volumetric delta-f Energy Flux (=df*KE*v_tot) (#/(m^2*s)) along tangential, normal, binormal directions, respectively
            vid=15;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_B_cov)
            vid=16;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_N_cov)
            vid=17;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),uvec_BxN_cov)

            vid=18;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_B_cntv(1))*uvec_B_cov)
            vid=19;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_N_cntv(1))*uvec_N_cov)
            vid=20;DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*energy(m)*delta_f(m)*dot_product(Full_dadgqdgf(:,m),sign(1.0,uvec_BxN_cntv(1))*uvec_BxN_cov)

            vid=1 ! Turbulent Particle Flux
            ! DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*(delta_f(m)*ghq(m)) ! Radial Flux
            ! DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*(delta_f(m)*Ta(m)) ! theta Flux
            ! DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*(delta_f(m)*dvi_gf(m)) ! zeta Flux
            DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*delta_f(m)*Full_dadgqdgf(1,m) ! Radial Flux
            DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*delta_f(m)*Full_dadgqdgf(2,m) ! theta Flux
            DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*delta_f(m)*Full_dadgqdgf(3,m) ! zeta Flux


            vid=2 ! Turbulent Particle Flux due to ExB only
            DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*(delta_f(m)*(ghq(m)-vpara_dvi_a(m))) ! Radial Flux
            DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*(delta_f(m)*(Ta(m)-vpara_dvi_gq(m))) ! theta Flux
            DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*(delta_f(m)*(dvi_gf(m)-vpara_dvi_gf(m)) )! zeta Flux

            vid=3 ! Turbulent Energy Flux
            ! DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*ghq(m)) ! Radial Flux
            ! DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Ta(m)) ! theta Flux
            ! DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*dvi_gf(m)) ! zeta Flux
            DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Full_dadgqdgf(1,m)) ! Radial Flux
            DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Full_dadgqdgf(2,m)) ! theta Flux
            DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*Full_dadgqdgf(3,m)) ! zeta Flux



            vid=4 ! Turbulent Energy Flux due to ExB only
            DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*(ghq(m)-vpara_dvi_a(m))) ! Radial Flux
            DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*(Ta(m)-vpara_dvi_gq(m))) ! theta Flux
            DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*energy(m)*(delta_f(m)*(dvi_gf(m)-vpara_dvi_gf(m)) )! zeta Flux

            vid=5 ! To measure CFL condition (v*dt/dx)
            DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*abs(tstep*((Tc/backT(m))**dsda_power*Full_dadgqdgf(1,m))*ori_inv_ds)
            DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*abs(tstep*(Full_dadgqdgf(2,m)-Full_dadgqdgf(3,m)/P_safety_q)*P_inv_ori_dgq) !time derivative of target_gq=gq-gf/safety_q
            DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*abs(tstep*Full_dadgqdgf(3,m)*ori_inv_dgf)
         enddo
      enddo

      ! Velocity Space Measurement
      do pbid=0,nPVS
         target_gq = modulo(P_gq + (PVS(pbid)%gf-P_gf)/P_safety_q, 2.0*pi)

         if( abs(P_gf-PVS(pbid)%gf)<0.5*PVS(pbid)%L_gf &
         .and. abs(P_a-PVS(pbid)%a) <0.5*PVS(pbid)%L_a &
         .and. min(abs(target_gq-PVS(pbid)%gq), 2.0*pi-abs(target_gq-PVS(pbid)%gq)) < 0.5*PVS(pbid)%L_gq) then

            norm_vpara = zion0(4,m)*BB(m)*sqrt(PVS(pbid)%inv_Ti0) ! norm_vpara = vpara/sqrt(T/m)
            norm_vperp = zion(6,m)*sqrt(2.0*BB(m))*sqrt(PVS(pbid)%inv_Ti0) ! norm_vperp = sqrt(mu)*sqrt(2B/m)/sqrt(T/m)

            vxid = floor(norm_vpara*PVS(pbid)%inv_dv)
            vyid = floor(norm_vperp*PVS(pbid)%inv_dv)
            rvx  = (norm_vpara-PVS(pbid)%vpara_1D(vxid))*PVS(pbid)%inv_dv
            rvy  = (norm_vperp-PVS(pbid)%vperp_1D(vyid))*PVS(pbid)%inv_dv

            if(vxid< -PVS(pbid)%nv) then
               vxid= -PVS(pbid)%nv
               rvx=0.0
            elseif (vxid>= PVS(pbid)%nv) then
               vxid=  PVS(pbid)%nv-1
               rvx=1.0
            endif

            if(vyid>= PVS(pbid)%nv) then
               vyid=  PVS(pbid)%nv-1
               rvy=1.0
            endif

            do j=1,4
               select case(j)
                  case(1)
                     qq=vxid; ss=vyid
                     ratio=(1.0-rvx)*(1.0-rvy)
                  case(2)
                     qq=vxid+1; ss=vyid
                     ratio=(rvx)*(1.0-rvy)
                  case(3)
                     qq=vxid; ss=vyid+1
                     ratio=(1.0-rvx)*(rvy)
                  case(4)
                     qq=vxid+1; ss=vyid+1
                     ratio=(rvx)*(rvy)
               end select

               fid=1; PVS(pbid)%fi(qq,ss,fid) = PVS(pbid)%fi(qq,ss,fid) + ratio*delta_f(m) ! Density
               fid=2; PVS(pbid)%fi(qq,ss,fid) = PVS(pbid)%fi(qq,ss,fid) + ratio*delta_f(m)*(2.0/3.0*energy(m)-backT(m)) ! Temperature definition
               fid=3; PVS(pbid)%fi(qq,ss,fid) = PVS(pbid)%fi(qq,ss,fid) + ratio ! Marker Density
               fid=4; PVS(pbid)%fi(qq,ss,fid) = PVS(pbid)%fi(qq,ss,fid) + ratio*f_tot(m) ! f_tot
               fid=5; PVS(pbid)%fi(qq,ss,fid) = PVS(pbid)%fi(qq,ss,fid) + ratio*f_tot(m)*(2.0/3.0*energy(m)) ! f_tot*energy
            enddo
         endif
      enddo
   enddo

   !$acc update self (zion_s_on_GyroOrbit(:,mi),zion_gq_on_GyroOrbit(:,mi))
   do m=1,mi
      if (Flag_w_Threshold==1 .and. (zion0(5,m)<w_Thold_min .or. zion0(5,m)>w_Thold_max)) cycle

      do larmor=1,larmor_num
         P_a  = zion0(1,m)
         P_gq = modulo(zion_gq_on_GyroOrbit(larmor,m),2*pi)
         P_gf = modulo(zion0(3,m),2*pi)
         P_s  = zion_s_on_GyroOrbit(larmor,m)

         ! ori_sid = (P_s-s0)*ori_inv_ds ! s index of original full grid
         ori_sid = max(0,min(int((P_s-s0)*ori_inv_ds),mpsi-1)) ! To avoid segmentation fault
         ff  = (P_s-s_grid(ori_sid))*ori_inv_ds
         P_safety_q = ff*q_grid(ori_sid+1) + (1.0-ff)*q_grid(ori_sid)
         P_inv_ori_dgq = 1.0/(ff*deltat(ori_sid+1)+(1.0-ff)*deltat(ori_sid))

         sid = max(0,min(int((P_s-s0)*DG%inv_ds) +1,DG%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
         ws1 = (P_s-DG%s_1D(sid))*DG%inv_ds
         ws0 = 1.0-ws1

         do k=0,1
            if(k==0) target_zeta =zetamin
            if(k==1) target_zeta =zetamax

            wgf = (1.0-abs(P_gf-target_zeta)/deltaz)

            target_gq = modulo(P_gq + (target_zeta-P_gf)/P_safety_q, 2.0*pi)
            gqid = (target_gq)*DG%inv_dgq +1
            wgq1 = (target_gq-DG%gq_1D(gqid))*DG%inv_dgq
            wgq0 = 1.0-wgq1

            do j=1,4
               select case(j)
                  case(1)
                     qq=gqid; ss=sid
                     ratio=ws0*wgq0*wgf
                  case(2)
                     qq=gqid+1; ss=sid
                     ratio=ws0*wgq1*wgf
                  case(3)
                     qq=gqid; ss=sid+1
                     ratio=ws1*wgq0*wgf
                  case(4)
                     qq=gqid+1; ss=sid+1
                     ratio=ws1*wgq1*wgf
               end select

               ! Scalar Data
               vid=7; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*inv_larmor_num*delta_f(m) ! Density
               vid=8; DG%iScaData(qq,ss,vid,k) = DG%iScaData(qq,ss,vid,k) + ratio*inv_larmor_num*(2.0/3.0*energy(m)-backT(m))*delta_f(m) !Temperature definition 1

               ! Vector Data
               vid=6 ! Turbulent Particle Flux
               DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*inv_larmor_num*(delta_f(m)*ghq(m)) ! Radial Flux
               DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*inv_larmor_num*(delta_f(m)*Ta(m)) ! theta Flux
               DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*inv_larmor_num*(delta_f(m)*dvi_gf(m)) ! zeta Flux

               vid=7 ! Turbulent Particle Flux due to ExB only
               DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*inv_larmor_num*(delta_f(m)*(ghq(m)-vpara_dvi_a(m))) ! Radial Flux
               DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*inv_larmor_num*(delta_f(m)*(Ta(m)-vpara_dvi_gq(m))) ! theta Flux
               DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*inv_larmor_num*(delta_f(m)*(dvi_gf(m)-vpara_dvi_gf(m)) )! zeta Flux

               vid=8 ! Turbulent Energy Flux
               DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*energy(m)*inv_larmor_num*(delta_f(m)*ghq(m)) ! Radial Flux
               DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*energy(m)*inv_larmor_num*(delta_f(m)*Ta(m)) ! theta Flux
               DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*energy(m)*inv_larmor_num*(delta_f(m)*dvi_gf(m)) ! zeta Flux

               vid=9 ! Turbulent Energy Flux due to ExB only
               DG%iVecData(1,qq,ss,vid,k) = DG%iVecData(1,qq,ss,vid,k) + ratio*energy(m)*inv_larmor_num*(delta_f(m)*(ghq(m)-vpara_dvi_a(m))) ! Radial Flux
               DG%iVecData(2,qq,ss,vid,k) = DG%iVecData(2,qq,ss,vid,k) + ratio*energy(m)*inv_larmor_num*(delta_f(m)*(Ta(m)-vpara_dvi_gq(m))) ! theta Flux
               DG%iVecData(3,qq,ss,vid,k) = DG%iVecData(3,qq,ss,vid,k) + ratio*energy(m)*inv_larmor_num*(delta_f(m)*(dvi_gf(m)-vpara_dvi_gf(m)) )! zeta Flux
            enddo
         enddo
      enddo
   enddo

   if(npartdom>1) then
      tmp_Scalar_Data = DG%iScaData
      tmp_Vector_Data = DG%iVecData

      length = DG%dim_gq*DG%dim_a*DG%num_iSD*2
      call MPI_ALLREDUCE(tmp_Scalar_Data, DG%iScaData, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)

      length = 3*DG%dim_gq*DG%dim_a*DG%num_iVD*2
      call MPI_ALLREDUCE(tmp_Vector_Data, DG%iVecData, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)

      allocate(tmp_PVS_fi, mold=PVS(0)%fi)
      do pbid=0,nPVS
         tmp_PVS_fi = PVS(pbid)%fi
         length = PVS(pbid)%dim_vpara*PVS(pbid)%dim_vperp*PVS(pbid)%num_of_fi
         call MPI_ALLREDUCE(tmp_PVS_fi, PVS(pbid)%fi, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
      enddo
      deallocate(tmp_PVS_fi)
   endif


   ! Match poloidal end cell
   ! (gq=0 and gq=2pi is the same position physically
   ! although they have different gqid)
   do sid=1,DG%dim_a
      DG%iScaData(1,sid,:,:) = DG%iScaData(1,sid,:,:) + DG%iScaData(DG%dim_gq,sid,:,:)
      DG%iScaData(DG%dim_gq,sid,:,:) = DG%iScaData(1,sid,:,:)

      DG%iVecData(:,1,sid,:,:) = DG%iVecData(:,1,sid,:,:) + DG%iVecData(:,DG%dim_gq,sid,:,:)
      DG%iVecData(:,DG%dim_gq,sid,:,:) = DG%iVecData(:,1,sid,:,:)
   enddo


   ! toroidal end cell
   length = DG%dim_gq*DG%dim_a*DG%num_iSD
   call MPI_SENDRECV(DG%iScaData(:,:,:,0),length,mpi_Rsize,left_pe,1, tmp_Scalar_Data(:,:,:,1),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   DG%iScaData(:,:,:,1) = DG%iScaData(:,:,:,1) + tmp_Scalar_Data(:,:,:,1)



   length = 3*DG%dim_gq*DG%dim_a*DG%num_iVD
   call MPI_SENDRECV(DG%iVecData(:,:,:,:,0),length,mpi_Rsize,left_pe,1, tmp_Vector_Data(:,:,:,:,1),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
   DG%iVecData(:,:,:,:,1) = DG%iVecData(:,:,:,:,1) + tmp_Vector_Data(:,:,:,:,1)


   ! Do normalization and post-processing
   ! Divide by Cell Volume to get density quantities
   do vid=1,DG%num_iSD
      DG%iScaData(:,:,vid,1) = DG%iScaData(:,:,vid,1)*DG%inv_Cell_dV
   enddo

   do aid=1,DG%dim_a
      DG%iScaData(:,aid,2,1) = DG%iScaData(:,aid,2,1)/DG%n0_1D(aid) ! Convert Energy to Temperature
   enddo


   ! Calculate inverse of tot density (including dne_adiabatic in the case of delta-h electron (ifelectron==2))
   do aid=1,DG%dim_a
      inv_tot_ni(:,aid) = 1.0/(DG%n0_1D(aid) + DG%iScaData(:,aid,1,1))
   enddo

   ! Consistent Temperature by considering density fluctuation in the denominator
   DG%iScaData(:,:,3,1) = DG%iScaData(:,:,3,1)*inv_tot_ni
   vid=6; DG%iScaData(:,:,vid,1) = DG%iScaData(:,:,vid,1)*inv_tot_ni! Te_tot
   vid=8; DG%iScaData(:,:,vid,1) = DG%iScaData(:,:,vid,1)*inv_tot_ni! Te_tot

   ! Get Flux Density by dividing Area in each direction
   do vid=1,4
      j=1; DG%iVecData(j,:,:,vid,1) = DG%iVecData(j,:,:,vid,1)*DG%inv_Cell_dS_a  ! Radial Flux Density
      j=2; DG%iVecData(j,:,:,vid,1) = DG%iVecData(j,:,:,vid,1)*DG%inv_Cell_dS_gq ! Theta Flux Density
      j=3; DG%iVecData(j,:,:,vid,1) = DG%iVecData(j,:,:,vid,1)*DG%inv_Cell_dS_gf ! Zeta Flux Density
   enddo

   ! To Normalize CFL condtion values, divide by number of markers (4th iSca)
   vid=5
   do j=1,3
      DG%iVecData(j,:,:,vid,1) = DG%iVecData(j,:,:,vid,1)/(DG%iScaData(:,:,4,1)*DG%Cell_dV)
   enddo

   ! Get Flux Density (using the positions on Larmor Radius) by dividing Area in each direction
   do vid=6,9
      j=1; DG%iVecData(j,:,:,vid,1) = DG%iVecData(j,:,:,vid,1)*DG%inv_Cell_dS_a  ! Radial Flux Density
      j=2; DG%iVecData(j,:,:,vid,1) = DG%iVecData(j,:,:,vid,1)*DG%inv_Cell_dS_gq ! Theta Flux Density
      j=3; DG%iVecData(j,:,:,vid,1) = DG%iVecData(j,:,:,vid,1)*DG%inv_Cell_dS_gf ! Zeta Flux Density
   enddo

end subroutine Measure_Ions_on_Diagnostic_Grid

subroutine Interpolate_OriGridData_into_DG(oriData_1D, out_IntpData_2D)

   use precision
   use global_parameters
   use field_array
   use Diagnostic_Grid
   implicit none

   real(wp), dimension(mgrid), intent(in) :: oriData_1D
   real(wp), dimension(DG%dim_gq,DG%dim_a), intent(out) :: out_IntpData_2D

   integer i,j,ij,jt, ori_aid

   real(wp) ori_inv_dgq(0:mpsi), gqdum, wt

   ori_inv_dgq=1.0/deltat

   do i=1,DG%dim_a
      do j=1,DG%dim_gq

         ori_aid = DG%ori_aid(i)

         gqdum = modulo(DG%gq_1D(j) - qtinv(ori_aid)*zetamax, 2.0*pi) * ori_inv_dgq(ori_aid)
         jt=int(gqdum);  wt=gqdum-jt
         ij=igrid(ori_aid)+jt
         out_IntpData_2D(j,i) =(1.0-wt)*oriData_1D(ij) + wt*oriData_1D(ij+1)
      enddo
   enddo

end subroutine Interpolate_OriGridData_into_DG


! module DG_operations

!    use Diagnostic_Grid

!    contains
!       subroutine Smooth_DG_Data_2D(num_smooth, target_Data2D, Reference_Data2D)

!          use precision
!          use global_parameters
!          use Diagnostic_Grid
!          implicit none

!          integer, intent(in) :: num_smooth
!          real(wp),intent(inout) :: target_Data2D(DG%dim_gq, DG%dim_a)
!          real(wp),optional:: Reference_Data2D(DG%dim_gq, DG%dim_a)
!          real(wp) :: VV(DG%dim_gq, DG%dim_a)
!          integer :: smid, aid, gqid, Na,Ngq, vid,k
!          real(wp) :: SM_Data(0:DG%dim_gq+1, 0:DG%dim_a+1)  ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
!          real(wp) :: SM_gq_Coeff(-1:1,0:DG%dim_gq+1, 0:DG%dim_a+1), SM_a_Coeff(-1:1,0:DG%dim_gq+1, 0:DG%dim_a+1)

!          Na  = DG%dim_a
!          Ngq = DG%dim_gq

!          VV=1.0
!          if(present(Reference_Data2D)) VV=Reference_Data2D

!          ! Get Smooth Coeffs which conserves Volume Weighted Integration
!          call Get_Conserving_Smooth_Coeffs_on_DG(Reference_Data2D, SM_a_Coeff, SM_gq_Coeff)

!          SM_Data = 0.0
!          SM_Data(1:Ngq,1:Na) = target_Data2D
!          do smid=1,num_smooth
!             ! Set ghost cell periodic in theta direction
!             SM_Data(0,:)     = SM_Data(Ngq,:)
!             SM_Data(Ngq+1,:) = SM_Data(1,:)
!             ! Smooth in theta direction
!             SM_Data(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*SM_Data(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*SM_Data(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*SM_Data(2:Ngq+1,:)
!             ! Smooth in radial direction
!             SM_Data(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*SM_Data(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*SM_Data(:,1:Na) + SM_a_Coeff(1,:,1:Na)*SM_Data(:,2:Na+1)
!          enddo
!          target_Data2D = SM_Data(1:Ngq,1:Na)

!       end subroutine Smooth_DG_Data_2D
! end module DG_operations


subroutine Smooth_DG_Scalar_Vector_Data(num_smooth)

   use precision
   use global_parameters
   use Diagnostic_Grid
   implicit none

   integer, intent(in) :: num_smooth
   integer :: smid, aid, gqid, Na,Ngq, vid,k
   real(wp) :: VV(0:DG%dim_gq+1, 0:DG%dim_a+1),  target_Data2D(0:DG%dim_gq+1, 0:DG%dim_a+1) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
   real(wp) :: SM_gq_Coeff(-1:1,0:DG%dim_gq+1, 0:DG%dim_a+1), SM_a_Coeff(-1:1,0:DG%dim_gq+1, 0:DG%dim_a+1)

   Na  = DG%dim_a
   Ngq = DG%dim_gq


   target_Data2D=0

   ! Get Smooth Coeffs which conserves Volume Weighted Integration
   call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dV, SM_a_Coeff, SM_gq_Coeff)
   ! =======================================
   ! Smooth Electron Scalar Data by conserving Volume Weighted Integration
   ! =======================================
   ! In the case of electrons, dne_adiab should be smoothed as well.
   target_Data2D(1:Ngq,1:Na) = DG%dne_adiab
   do smid=1,num_smooth
      ! Set ghost cell periodic in theta direction
      target_Data2D(0,:)     = target_Data2D(Ngq,:)
      target_Data2D(Ngq+1,:) = target_Data2D(1,:)
      ! Smooth in theta direction
      target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
      ! Smooth in radial direction
      target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
   enddo
   DG%dne_adiab = target_Data2D(1:Ngq,1:Na)

   do vid=1,DG%num_eSD
      target_Data2D(1:Ngq,1:Na) = DG%eScaData(:,:,vid,1)
      do smid=1,num_smooth
         ! Set ghost cell periodic in theta direction
         target_Data2D(0,:)     = target_Data2D(Ngq,:)
         target_Data2D(Ngq+1,:) = target_Data2D(1,:)
         ! Smooth in theta direction
         target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
         ! Smooth in radial direction
         target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
      enddo
      DG%eScaData(:,:,vid,1) = target_Data2D(1:Ngq,1:Na)
   enddo
   ! =======================================
   ! Smooth Ion Scalar Data by conserving Volume Weighted Integration
   ! =======================================
   do vid=1,DG%num_iSD
      target_Data2D(1:Ngq,1:Na) = DG%iScaData(:,:,vid,1)
      do smid=1,num_smooth
         ! Set ghost cell periodic in theta direction
         target_Data2D(0,:)     = target_Data2D(Ngq,:)
         target_Data2D(Ngq+1,:) = target_Data2D(1,:)
         ! Smooth in theta direction
         target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
         ! Smooth in radial direction
         target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
      enddo
      DG%iScaData(:,:,vid,1) = target_Data2D(1:Ngq,1:Na)
   enddo

   ! =======================================
   ! Smooth Electron Vector Data by conserving Surface Area Weighted Integration in each direction
   ! =======================================
   do vid=1,DG%num_eVD
      do k=1,3
         select case(k)
            case(1)
               call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dS_a, SM_a_Coeff, SM_gq_Coeff)
            case(2)
               call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dS_gq, SM_a_Coeff, SM_gq_Coeff)
            case(3)
               call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dS_gf, SM_a_Coeff, SM_gq_Coeff)
         end select

         target_Data2D(1:Ngq,1:Na) = DG%eVecData(k,:,:,vid,1)
         do smid=1,num_smooth
            ! Set ghost cell periodic in theta direction
            target_Data2D(0,:)     = target_Data2D(Ngq,:)
            target_Data2D(Ngq+1,:) = target_Data2D(1,:)
            ! Smooth in theta direction
            target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
            ! Smooth in radial direction
            target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
         enddo
         DG%eVecData(k,:,:,vid,1) = target_Data2D(1:Ngq,1:Na)
      enddo
   enddo

   ! =======================================
   ! Smooth Ion Vector Data by conserving Surface Area Weighted Integration in each direction
   ! =======================================
   do vid=1,DG%num_iVD
      do k=1,3
         select case(k)
         case(1)
            call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dS_a, SM_a_Coeff, SM_gq_Coeff)
         case(2)
            call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dS_gq, SM_a_Coeff, SM_gq_Coeff)
         case(3)
            call Get_Conserving_Smooth_Coeffs_on_DG(DG%Cell_dS_gf, SM_a_Coeff, SM_gq_Coeff)
         end select

         target_Data2D(1:Ngq,1:Na) = DG%iVecData(k,:,:,vid,1)
         do smid=1,num_smooth
            ! Set ghost cell periodic in theta direction
            target_Data2D(0,:)     = target_Data2D(Ngq,:)
            target_Data2D(Ngq+1,:) = target_Data2D(1,:)
            ! Smooth in theta direction
            target_Data2D(1:Ngq,:) = SM_gq_Coeff(-1,1:Ngq,:)*target_Data2D(0:Ngq-1,:) + SM_gq_Coeff(0,1:Ngq,:)*target_Data2D(1:Ngq,:) + SM_gq_Coeff(1,1:Ngq,:)*target_Data2D(2:Ngq+1,:)
            ! Smooth in radial direction
            target_Data2D(:,1:Na) = SM_a_Coeff(-1,:,1:Na)*target_Data2D(:,0:Na-1) + SM_a_Coeff(0,:,1:Na)*target_Data2D(:,1:Na) + SM_a_Coeff(1,:,1:Na)*target_Data2D(:,2:Na+1)
         enddo
         DG%iVecData(k,:,:,vid,1) = target_Data2D(1:Ngq,1:Na)
      enddo
   enddo

end subroutine Smooth_DG_Scalar_Vector_Data


subroutine Get_Conserving_Smooth_Coeffs_on_DG(in_VOLorAREA, SM_a_Coeff, SM_gq_Coeff)

   use precision
   use global_parameters
   use field_array
   use Diagnostic_Grid
   implicit none

   integer :: Na,Ngq
   real(wp), intent(in) :: in_VOLorAREA(DG%dim_gq, DG%dim_a) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
   real(wp) :: VV(0:DG%dim_gq+1, 0:DG%dim_a+1) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
   real(wp), intent(out):: SM_gq_Coeff(-1:1,0:DG%dim_gq+1, 0:DG%dim_a+1), SM_a_Coeff(-1:1,0:DG%dim_gq+1, 0:DG%dim_a+1)

   Na  = DG%dim_a
   Ngq = DG%dim_gq

   VV = 0
   VV(1:Ngq,1:Na) = in_VOLorAREA

   ! Note that theta direction is periodic
   VV(0,:) = VV(Ngq,:)
   VV(Ngq+1,:) = VV(1,:)

   SM_gq_Coeff(-1,1:Ngq,:) = 0.5*VV(0:Ngq-1,:)/(VV(0:Ngq-1,:) + VV(1:Ngq,:)) ! Coeff_L = 0.5*VL/(VL+VC)
   SM_gq_Coeff( 0,1:Ngq,:) = 0.5*VV(1:Ngq,:)*(VV(0:Ngq-1,:)+2.0*VV(1:Ngq,:)+VV(2:Ngq+1,:))/((VV(0:Ngq-1,:)+VV(1:Ngq,:))*(VV(1:Ngq,:)+VV(2:Ngq+1,:))) ! 0.5*VC*(VL+2*VC+VR)/((VL+VC)*(VC+VR))
   SM_gq_Coeff( 1,1:Ngq,:) = 0.5*VV(2:Ngq+1,:)/(VV(2:Ngq+1,:) + VV(1:Ngq,:)) ! Coeff_R = 0.5*VR/(VR+VC)

   SM_a_Coeff(-1,:,1:Na) = 0.5*VV(:,0:Na-1)/(VV(:,0:Na-1) + VV(:,1:Na)) ! Coeff_L = 0.5*VL/(VL+VC)
   SM_a_Coeff( 0,:,1:Na) = 0.5*VV(:,1:Na)*(VV(:,0:Na-1)+2.0*VV(:,1:Na)+VV(:,2:Na+1))/((VV(:,0:Na-1)+VV(:,1:Na))*(VV(:,1:Na)+VV(:,2:Na+1))) ! 0.5*VC*(VL+2*VC+VR)/((VL+VC)*(VC+VR))
   SM_a_Coeff( 1,:,1:Na) = 0.5*VV(:,2:Na+1)/(VV(:,2:Na+1) + VV(:,1:Na)) ! Coeff_R = 0.5*VR/(VR+VC)

   ! To remove NaN values at the raidal ghost cells
   SM_gq_Coeff(:,:,0)   = 0
   SM_gq_Coeff(:,:,Na+1)= 0
   SM_a_Coeff(:,:,0)    = 0
   SM_a_Coeff(:,:,Na+1) = 0

end subroutine Get_Conserving_Smooth_Coeffs_on_DG