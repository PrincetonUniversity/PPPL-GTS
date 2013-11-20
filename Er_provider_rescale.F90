!subroutine Er_provider(Er_from_GTC_Neo,n_node,Er_data_name)
subroutine Er_provider
 
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
          xe,ac,half_width,Er_amp,Yrr,dum, &
          Er_scale
  real(wp), external :: temperature,density,u_parallel,dT_da,dn_da
  integer i, ifail,ierror

!  integer Er_from_GTC_Neo,n_node
!  character(25) Er_data_name
  real(wp), dimension(:), allocatable :: Er_on_node,r_node,E01_K,E01_C,E01_WRK
      


      je=jr
      allocate(Yr(0:je), coty(0:je), Er401_pc(0:je), Er401_bp(0:je),Er401_bp0(0:je))
      xe=2.0_wp
!      ac=0.35
      ac=0.50_wp
      half_width=0.08_wp
!      Er_amp=1.0d-4
!      Er_amp=2.0d-4
!      Er_amp=6.0d-5
      Er_amp=0.0_wp
      Er_scale=0.5_wp

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

!   given a model Er profile
!         Er401_bp(i)=Er_amp*exp(-((Yr(i)-ac)/half_width)**xe)

         Er401_bp(i)=Er401_bp(i)*Er_scale
         Er401_pc(i)=Er401_pc(i)*Er_scale
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
      end do

      do i=0,mpsi
!         Er_grid(i)=Er_amp*exp(-((psi_grid(i)-ac)/half_width)**xe)
      end do
     
      Er_grid=Er_grid*Er_scale

  if(Er_from_GTC_Neo==1) then  ! import Er data from GTC-Neo simulation
                               ! n_node+2: number of radial nodes of Er data
                               ! Er_on_node: Er data at nodes
                               ! r_node:     radial location of nodes
                               ! E01_K,E01_C,E01_WRK: spline coefficients and working space 
     if(mype==0) then
        allocate(Er_on_node(n_node+2),r_node(n_node+2),E01_K(n_node+2+4),E01_C(n_node+2+4),E01_WRK(6*(n_node+2)+16))
        open(unit=20,file=Er_data_name)
        do i=1,n_node+2
           read(20,100) r_node(i), Er_on_node(i)

! rescale Er to examine sensitivity of turbulence on Er shear
!           Er_on_node(i)=Er_on_node(i)*0.7       ! rescale Er to examine sensitivity on Er shear

        end do
100  format(10(1x,ES13.6))

!        ifail=0
!        call E01BAF(n_node+2,r_node,Er_on_node,E01_K,E01_C,n_node+2+4,E01_WRK,6*(n_node+2)+16,ifail)
        call ezspline_init(spEr,n_node+2,(/1,0/),ierror)
        spEr%x1=r_node
        spEr%bcval1min=0.0d0
        call ezspline_setup(spEr,Er_on_node,ierror)

        do i=0,je
!           ifail=0
!           call E02BBF(n_node+2+4,E01_K,E01_C,Yr(i),Er401_bp(i),ifail)
           call ezspline_interp(spEr,Yr(i),Er401_bp(i),ierror)
        end do
        do i=0,mpsi
!           ifail=0
!           call E02BBF(n_node+2+4,E01_K,E01_C,psi_grid(i),Er_grid(i),ifail)
           call ezspline_interp(spEr,psi_grid(i),Er_grid(i),ierror)
        end do
     end if

        call ezspline_free(spEr,ierror)

     call MPI_BCAST(Er401_bp,je+1,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(Er_grid,mpsi+1,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  end if 
 
!  Er401_bp=1.2_wp*Er401_bp
!  Er_grid=1.2_wp*Er_grid
  Er401_bp0=Er401_bp
  Er_grid0=Er_grid

  if(mype==0) then
     open(502,file='Er_estimate.dat')

     do i=0,je
        write(502,5) i, Yr(i), Er401_bp(i), Er401_pc(i) 
        !print*, ' Er_nc:', i, Yr(i), Er401_bp(i), Er401_pc(i)
     end do

     call flush(502)

5    format(2x, i5, 3(2x,ES13.6))
  end if
 
end subroutine Er_provider
