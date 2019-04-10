subroutine get_ZF_shear(LINES)
 
!      Evaluation of neoclassical Er
 
  use precision
  use mpi_defs
  use global_parameters
  use particle_array
  use field_array
  use esi_interface
  use spl_profiles

  implicit none
 
  real(wp), dimension(10) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the
  real(wp)  dlnn_da,dlnT_da,taur,vthr,da,hint1,ep,hint2,u_int,taui,qq,  &
          B_p0,rho_0,R2, &
          xe,ac,half_width,Er_amp,Yrr,dum, grad_a2,a_minor
  real(wp), external :: temperature,density,u_parallel,dT_da,dn_da
  integer i, ifail,ierror,j,LINES

!  integer Er_from_GTC_Neo,n_node
!  character(25) Er_data_name
  real(wp), dimension(:), allocatable :: Er_on_node,r_node
  real(wp), dimension(:), allocatable :: Er401_TRANSP, Er_grid_TRANSP

      
  real(wp) phia00(0:mpsi), phiaa(0:mpsi), omega_E(0:mpsi), Era_grid(0:mpsi)


      je=jr
      allocate(Yr(0:je), coty(0:je), Er401_pc(0:je), Er401_bp(0:je),Er401_bp0(0:je))
      allocate(Er401_TRANSP(0:je), Er_grid_TRANSP(0:mpsi))
      xe=2.0_wp
!      ac=0.35
      ac=0.50_wp
      half_width=0.08_wp
!      Er_amp=1.0d-4
!      Er_amp=2.0d-4
!      Er_amp=6.0d-5
      Er_amp=0.0_wp

!      Er=dPhi/da=Er_amp*exp{-[(a-ac)/half_width]^xe}

if(mype==0) then 
   print*, 'Er=dPhi/da=Er_amp*exp{-[(a-ac)/half_width]^xe}' 
   print*, 'Er_amp=', Er_amp, '  ac=', ac, '  xe=', xe, '  half_width=', half_width
end if

      aaa(1)=1.0_wp
      the(1)=0.0_wp
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      rho_0=sqrt(gFaj(1)/B_0)    ! minor radius defined in (5.35) of HH; da/dr=1/rho_0
      R2=rj(1)
      the(1)=pi
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      if(mype==0)write(stdout,*)'rho_0=',rho_0,R2-R_0,(R2-rj(1))/2.0_wp
      a_minor=(R2-rj(1))/2.0d0

      Ydr=(a1-a0)/real(je,wp)
      do 29 i=0,je
         Yr(i)=a0+Ydr*real(i,wp)
         vthr=sqrt(2.0_wp*temperature(Yr(i),0.0_wp))

         Taur=0.25_wp/pi*pmass0  &
              *(vthr*uveloci)**3.0_wp/(density(Yr(i),0.0_wp)*charg0*charg0  &
              *charg0*charg0*dlni)*pmass0   ! in unit of sec
         Taui=0.75_wp*sqrt(pi)*Taur
         Taui=Taui*sqrt(2.0_wp)/utime      ! normalized collision time 
 
         U_int=u_parallel(Yr(i),0.0_wp)

         aaa(1)=Yr(i)
         the(1)=0.0_wp
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)

         B_p0=-gYaj(1)/(R_0*rho_0)  ! poloidal field defined in (5.36) of HH

         U_int=U_int*Bj(1)   ! =B*u_i//=B_0*u(r)

         R2=rj(1)
         qq=-gFaj(1)/gYaj(1)
         the(1)=pi
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
         ep=0.5_wp*abs(R2-rj(1))/R_0  ! aspect ratio
         coty(i)=sqrt(2.0_wp)*qq/(ep**1.5_wp)*(R_0/vthr)/Taui

         hint1=(1.17_wp-0.35_wp*sqrt(coty(i)))/(1.0_wp+0.7_wp*sqrt(coty(i)))
         hint2=(hint1-2.1_wp*coty(i)*coty(i)*ep*ep*ep)/(1.0_wp+coty(i)*coty(i)*ep*ep*ep)
    
         dlnT_da=dT_da(Yr(i),0._wp)/temperature(Yr(i),0._wp)
         dlnn_da=dn_da(Yr(i),0._wp)/density(Yr(i),0._wp)
      
                 ! Er=dPhi/da
         Er401_bp(i)=((hint1-1._wp)*dlnT_da-dlnn_da)*0.5_wp*vthr**2-U_int*gFaj(1)/(qq*Fj(1))
         Er401_pc(i)=((hint2-1._wp)*dlnT_da-dlnn_da)*0.5_wp*vthr**2-U_int*gFaj(1)/(qq*Fj(1))

         if(plasma==1) call ezspline_derivative(spPHI,1,Yr(i),Er401_TRANSP(i),ierror) ! Er from TRANSP

!   given a model Er profile
!         Er401_bp(i)=Er_amp*exp(-((Yr(i)-ac)/half_width)**xe)

29    continue

      do i=0,mpsi
         Yrr=psi_grid(i)
         vthr=sqrt(2._wp*temperature(Yrr,0._wp))

         Taur=0.25_wp/pi*pmass0  &
              *(vthr*uveloci)**3._wp/(density(Yrr,0._wp)*charg0*charg0  &
              *charg0*charg0*dlni)*pmass0   ! in unit of sec
         Taui=0.75_wp*sqrt(pi)*Taur
         Taui=Taui*dsqrt(2._wp)/utime      ! normalized collision time

         U_int=u_parallel(Yrr,0._wp)

         aaa(1)=Yrr
         the(1)=0._wp
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)

         B_p0=-gYaj(1)/(R_0*rho_0)  ! poloidal field defined in (5.36) of HH

         U_int=U_int*Bj(1)   ! =B*u_i//=B_0*u(r)

         R2=rj(1)
         qq=-gFaj(1)/gYaj(1)
         the(1)=pi
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
         ep=0.5_wp*abs(R2-rj(1))/R_0  ! aspect ratio
         dum=sqrt(2._wp)*qq/(ep**1.5_wp)*(R_0/vthr)/Taui

         hint1=(1.17_wp-0.35_wp*sqrt(dum))/(1._wp+0.7_wp*sqrt(dum))
         hint2=(hint1-2.1_wp*dum*dum*ep*ep*ep)/(1._wp+dum*dum*ep*ep*ep)
   
         dlnT_da=dT_da(Yrr,0._wp)/temperature(Yrr,0._wp)
         dlnn_da=dn_da(Yrr,0._wp)/density(Yrr,0._wp)

                 ! Er=dPhi/da
         Er_grid(i)=((hint1-1._wp)*dlnT_da-dlnn_da)*0.5_wp*vthr**2-U_int*gFaj(1)/(qq*Fj(1)) ! using bp regime

         if(plasma==1) call ezspline_derivative(spPHI,1,Yr(i),Er_grid_TRANSP(i),ierror) ! =Er from TRANSP
      end do

!      do i=0,mpsi
!         Er_grid(i)=Er_amp*exp(-((psi_grid(i)-ac)/half_width)**xe)
!      end do
     


  if(Er_from_GTC_Neo==1) then  ! import Er data from GTC-Neo simulation
                               ! n_node+2: number of radial nodes of Er data
                               ! Er_on_node: Er data at nodes
                               ! r_node:     radial location of nodes
                               ! E01_K,E01_C,E01_WRK: spline coefficients and working space
     if(mype==0) then
        allocate(Er_on_node(n_node+2),r_node(n_node+2))
        open(unit=20,file=Er_data_name)
        do i=1,n_node+2
           read(20,100) r_node(i), Er_on_node(i)

! rescale Er to examine sensitivity of turbulence on Er shear
!           Er_on_node(i)=Er_on_node(i)*0.7       ! rescale Er to examine sensitivity on Er shear

        end do
100  format(10(1x,ES13.6))

        call ezspline_init(spEr0,n_node+2,(/1,0/),ierror)
        spEr0%x1=r_node
        spEr0%bcval1min=0.0d0
        call ezspline_setup(spEr0,Er_on_node,ierror)

        do i=0,mpsi
           call ezspline_interp(spEr0,psi_grid(i),Er_grid(i),ierror)
           call ezspline_derivative(spEr0,1,psi_grid(i),Era_grid(i),ierror) ! Era_grid=(d^2/da^2)phi
        end do
        deallocate(Er_on_node,r_node)
     end if

     call MPI_BCAST(Er_grid,mpsi+1,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(Era_grid,mpsi+1,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  end if

  call ezspline_free(spEr0,ierror)

!   calculate EXB shear rate of ZF

    phiaa=0.0
    if(mype==0) then

print*, 'EXB SHEAR with  vthc/a_minor= ', vthc, a_minor

       open(502,file='ZF_shearate.dat',status='replace',RECL=8192) 
       do j=1,LINES
          !  read(122,801) istep, (phia00(i),i=0,19)  ! not work after edison relocated 2016
          read(122,*) istep, (phia00(i),i=0,19)
          read(123,*) istep, (phia00(i),i=20,39)
          read(124,*) istep, (phia00(i),i=40,59)
          read(125,*) istep, (phia00(i),i=60,79)
          read(126,*) istep, (phia00(i),i=80,99)
!          read(126,801) istep, (phia00(i),i=80,mpsi)
          read(127,*) istep, (phia00(i),i=100,119)
!          read(128,801) istep, (phia00(i),i=120,mpsi)
          read(128,*) istep, (phia00(i),i=120,139)
          read(129,*) istep, (phia00(i),i=140,159)
          read(130,*) istep, (phia00(i),i=160,179)
          read(131,*) istep, (phia00(i),i=180,199)
          read(132,*) istep, (phia00(i),i=200,219)
          read(133,*) istep, (phia00(i),i=220,239)
          read(134,*) istep, (phia00(i),i=240,259)
          read(135,*) istep, (phia00(i),i=260,279)
          read(136,*) istep, (phia00(i),i=280,mpsi)
!          read(132,*) istep, (phia00(i),i=200,mpsi)
!          read(129,801) istep, (phia00(i),i=140,mpsi)
!          read(130,*) istep, (phia00(i),i=160,mpsi)
!          read(131,801) istep, (phia00(i),i=180,mpsi)

!         -------------------------------------------------------------------------
!         if calculate total shear ZF +equilibrium Er, use:
!          phia00=phia00+Er_grid     ! --> phia00= total Er=d<phi>/da +d(phi_0)/da 
!          phia00=phia00+Er_grid_TRANSP     ! --> phia00= total Er=d<phi>/da + d(phi_0)/da with d(phi_0)/da from TRANSP 
!          phia00=Er_grid_TRANSP     ! --> phia00=  d(phi_0)/da with d(phi_0)/da from TRANSP 
!         -------------------------------------------------------------------------

          call ezspline_init(spEr,mpsi+1,(/1,1/),ierror) ! 1st derivative boundary condition for both left and right
          spEr%x1=psi_grid
          spEr%bcval1min=0.0d0  ! derivative at left
          spEr%bcval1max=0.0d0  ! derivative at right
          call ezspline_setup(spEr,phia00,ierror)

!   print*, 'read(130 0000', j

          do i=0,mpsi
             call ezspline_derivative(spEr,1,psi_grid(i),phiaa(i),ierror) ! phiaa=(d^2/da^2)phi
             aaa(1)=psi_grid(i)  
             the(1)=0.0d0
             call esiget2dfunctions(aaa,the,1)
             grad_a2=(rq(1)*rq(1)+zq(1)*zq(1))/(za(1)*rq(1)-zq(1)*ra(1))**2

             omega_E(i)=-grad_a2/BB(1)*(phiaa(i)-phia00(i)*gYaa(1)/gYa(1))        ! only ZF shear or total shear depending on 
                                                                                  ! whether phia00 includes EQ Er or not  Jan., 2010

!             
!             omega_E(i)=-grad_a2/BB(1)*((phiaa(i)+Era_grid(i))-(phia00(i)+Er_grid(i))*gYaa(1)/gYa(1)) ! total shear
!                                   ! not used; total shear can be included via add Er_grid(i) to phia00(i);  Jan., 2010
!

             omega_E(i)=omega_E(i)/(vthc/a_minor)  ! normalized by c_s/a

!             omega_E(i)=abs(omega_E(i))            ! taking absolute value
          end do

          call ezspline_free(spEr,ierror)

          write(502,5) istep, (omega_E(i),i=0,mpsi)

          call flush(502)
       end do
       close(502)
    end if





!   calculate intensity gradient

    phiaa=0.0
    if(mype==0) then

       open(502,file='Intensity_gradient.dat',status='replace',RECL=8192)
       do j=1,LINES
          !  read(80,802) istep, (phia00(i),i=0,mpsi)  ! not working after edison relocated
          read(80,*) istep, (phia00(i),i=0,mpsi)

          call ezspline_init(spEr,mpsi+1,(/1,1/),ierror) ! 1st derivative boundary condition for both left and right
          spEr%x1=psi_grid
          spEr%bcval1min=0.0d0  ! derivative at left
          spEr%bcval1max=0.0d0  ! derivative at right
          call ezspline_setup(spEr,phia00,ierror)

          do i=0,mpsi
             call ezspline_derivative(spEr,1,psi_grid(i),phiaa(i),ierror) ! phiaa=(d/da)(delta_phi^2)
          end do

          call ezspline_free(spEr,ierror)

          write(502,5) istep, (phiaa(i),i=0,mpsi)

          call flush(502)
       end do
       close(502)
    end if





5   format(1x,i6, 503(1x,d13.6)) ! for maxium mpsi=500
801 format(1x,i6, 21(1x,d13.6))       ! new format for mstep>100000
! 801 format(1x,i5, 21(1x,d13.6))       ! new format for mstep>100000
! 801 format(i6, 21(1x,d13.6))

802 format(1x,i6, 503(1x,d13.6)) ! for maxium mpsi=500

  if(mype==0) then
     open(502,file='Er_estimate.dat')

     do i=0,je
        write(502,50) i, Yr(i), Er401_bp(i), Er401_pc(i)
     end do

     call flush(502)

50   format(2x, i5, 3(2x,ES13.6))
  end if

      deallocate(Er401_TRANSP, Er_grid_TRANSP)

end subroutine get_ZF_shear
