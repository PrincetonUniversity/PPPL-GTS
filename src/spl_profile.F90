module profile_array
    use precision
    real(wp), dimension(:), allocatable :: rcor,ne_50,Ti_50,Te_50,O_50,PHI_50
    integer  ndum,ier
    integer  jj1,jj2
    real(wp), dimension(:), allocatable :: Ti_node,dTi_node,ne_node,dne_node,&
                                           Te_node,dTe_node,Iomg_node,dIomg_node
    real(wp), dimension(:), allocatable :: rcor2,ne_51,Ti_51,Te_51
    real(wp), dimension(:), allocatable :: dne_51,dTi_51,dTe_51
    real(wp) :: f1, f2, f3
end module profile_array



module spl_profiles

   use ezspline
   use ezspline_obj

   type(ezspline1_r8):: spn,spTi,spO,spTe,spEr,spEr0,spPHI

end module spl_profiles



subroutine read_profiles(profile_data_name)

  use precision
  use mpi_defs
  use global_parameters, only:mype
  use spl_profiles
  use profile_array
  use control_parameters, only:ain, aout

  implicit none

  real(wp), dimension(:,:), allocatable :: worka
  integer i,j,k,ierr
  character(36) chdum
  character(len=36),dimension(6) :: profname
  character(25)  profile_data_name

  real(wp) xx, acdum, error_8, acc
  integer minpts,ifail
  real(wp), external :: D01AHF
  real(wp), external :: fln,flTe,flTi,frn,frTe,frTi

! Names of the profiles that we want to read out of the file
  profname(1)=' Zone centrs = sq_rt(norm toroidal f'   ! radii (r/a)
  profname(2)=' Electron dens                      '   ! Ne
  profname(3)=' Hydrogenic Temperature             '   ! Ti
  profname(4)=' Electron Temperature               '   ! Te
  profname(5)=' thermal plasma toroidal rotation   '   ! Toroidal rotation
  profname(6)=' Radial Potential due to rotation NC'   ! Equilibrium Radial electric potential
  !    profname(6)=' Radial Electrical Potential        '   ! Equilibrium Radial electric potential

! The root process reads the file containing the profiles.
! We first read the number of radial grid points there are for each profile.
  if(mype == 0)then
     open(unit=3,file=profile_data_name)
     read(3,'(a36)')chdum
     do while(chdum(2:7)/='Nzones')
        read(3,'(a36)')chdum
Print*, chdum
     enddo
     read(3,10) ndum
     write(*,*)'in read_profiles: Number of grid points ndum=',ndum
  endif
!         10 format(12x,I3)
10 format(11x,I3)

! The value of "ndum" is broadcast to the other processes for memory allocation
  call mpi_bcast(ndum,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  allocate(worka(ndum,6),rcor(ndum+2),ne_50(ndum+2),Ti_50(ndum+2),&
           Te_50(ndum+2),O_50(ndum+2),PHI_50(ndum+2))
  allocate(rcor2(ndum),ne_51(ndum),Ti_51(ndum),Te_51(ndum),dne_51(ndum),dTi_51(ndum),dTe_51(ndum))

! The root process continues to read the file and store the relevant profiles
! in the work array "worka".
! We read the file line by line and keep only the profiles corresponding to
! the names profname(1:6)
  if(mype == 0)then
     do i=1,6
        do while(chdum/=profname(i))
           read(3,'(a36)')chdum
        enddo
       read(3,20) worka(1:ndum,i)
     enddo
     close(3) ! close file

     open(unit=3,file='TRANSP_profs.d')   ! write original TRANSP profiles to TRANSP_profs.d
     do i=1,ndum
        write(3,30) (worka(i,j), j=1,6)
     end do
30   format(6(2x, ES13.4))
     close(3) ! close file

  endif
20 format(5(ES13.4))



!!! added to modify profiles automatically   Nov. 2014
  if(mype == 0)then

     !  ain=0.20_wp
     !  aout=0.63_wp

     f1=2.0_wp
     f2=2.0_wp
     f3=2.0_wp


     !  jj1=8
     !  jj2=27

     do j=1,ndum
        if(worka(j,1)>ain) then
           jj1=j-1
           go to 100
        endif
     enddo

100  do j=1,ndum
        if(worka(j,1)>aout) then
           jj2=j
           go to 200
        endif
     enddo

200  print*, 'jj1=',jj1, 'jj2=',jj2   ! profiles Ti, Te and ne are modified at boundary layer inside jj1 and outside jj2
                                      ! say, dTi/da=dTi/da(jj2-1)*fdown and
                                      ! fdown=exp{-f2*[(r-rc))/(f1*rc)]^f3} with rc=r_jj2-r_(jj2-1)

     !  allocate(rcor2(ndum),ne_51(ndum),Ti_51(ndum),Te_51(ndum),dne_51(ndum),dTi_51(ndum),dTe_51(ndum))
     rcor2(1:ndum)=worka(1:ndum,1)
     ne_51(1:ndum)=worka(1:ndum,2)
     Ti_51(1:ndum)=worka(1:ndum,3)
     Te_51(1:ndum)=worka(1:ndum,4)

     call spl_init2

     do j=1,ndum
        xx=rcor2(j)
        call ezspline_derivative(spn,1,xx,dne_51(j),ier)

        call ezspline_derivative(spTi,1,xx,dTi_51(j),ier)

        call ezspline_derivative(spTe,1,xx,dTe_51(j),ier)
     enddo

     error_8=1.0e-6_wp
     do j=1,jj1
        k=jj1-j+1
        xx=rcor2(k)
        acdum=rcor2(jj1+1)
        ifail=1
        ne_51(k)=ne_51(jj1+1) + D01AHF(acdum,xx,error_8,minpts,acc,fln,-1,ifail)     !1D integral over r
        ifail=1
        Te_51(k)=Te_51(jj1+1) + D01AHF(acdum,xx,error_8,minpts,acc,flTe,-1,ifail)     !1D integral over r
        ifail=1
        Ti_51(k)=Ti_51(jj1+1) + D01AHF(acdum,xx,error_8,minpts,acc,flTi,-1,ifail)     !1D integral over r
     enddo

     do j=jj2,ndum
        xx=rcor2(j)
        acdum=rcor2(jj2-1)
        ifail=1
        ne_51(j)=ne_51(jj2-1) + D01AHF(acdum,xx,error_8,minpts,acc,frn,-1,ifail)      !1D integral over r
        ifail=1
        Te_51(j)=Te_51(jj2-1) + D01AHF(acdum,xx,error_8,minpts,acc,frTe,-1,ifail)     !1D integral over r
        ifail=1
        Ti_51(j)=Ti_51(jj2-1) + D01AHF(acdum,xx,error_8,minpts,acc,frTi,-1,ifail)     !1D integral over r
     enddo

     worka(1:ndum,2)=ne_51(1:ndum)
     worka(1:ndum,3)=Ti_51(1:ndum)
     worka(1:ndum,4)=Te_51(1:ndum)

   ! Delete and deallocate the spline coefficients so that we can reuse the sp* structures
     call ezspline_free(spn,ier)
     call ezspline_free(spTi,ier)
     call ezspline_free(spTe,ier)
  endif



! The root process broadcast the work array "worka"
  call mpi_bcast(worka,6*ndum,mpi_Rsize,0,MPI_COMM_WORLD,ierr)

! Unpack the data to the relevant arrays
  rcor(2:ndum+1)=worka(1:ndum,1)
  ne_50(2:ndum+1)=worka(1:ndum,2)
  Ti_50(2:ndum+1)=worka(1:ndum,3)
  Te_50(2:ndum+1)=worka(1:ndum,4)
  O_50(2:ndum+1)=worka(1:ndum,5)
  PHI_50(2:ndum+1)=worka(1:ndum,6)

  deallocate(worka)  ! We don't need the work array anymore

  rcor(1)=0.0_wp
  rcor(ndum+2)=1.0_wp

  ne_50(1)=(rcor(3)**2*ne_50(2)-rcor(2)**2*ne_50(3))/(rcor(3)**2-rcor(2)**2)
  ne_50(ndum+2)=ne_50(ndum+1)+(rcor(ndum+2)-rcor(ndum+1))  &
                /(rcor(ndum+1)-rcor(ndum))*(ne_50(ndum+1)-ne_50(ndum))

  Ti_50(1)=(rcor(3)**2*Ti_50(2)-rcor(2)**2*Ti_50(3))/(rcor(3)**2-rcor(2)**2)
  Ti_50(ndum+2)=Ti_50(ndum+1)+(rcor(ndum+2)-rcor(ndum+1)) &
                /(rcor(ndum+1)-rcor(ndum))*(Ti_50(ndum+1)-Ti_50(ndum))

  Te_50(1)=(rcor(3)**2*Te_50(2)-rcor(2)**2*Te_50(3))/(rcor(3)**2-rcor(2)**2)
  Te_50(ndum+2)=Te_50(ndum+1)+(rcor(ndum+2)-rcor(ndum+1)) &
                /(rcor(ndum+1)-rcor(ndum))*(Te_50(ndum+1)-Te_50(ndum))

  O_50(1)=(rcor(3)**2*O_50(2)-rcor(2)**2*O_50(3))/(rcor(3)**2-rcor(2)**2)
  O_50(ndum+2)=O_50(ndum+1)+(rcor(ndum+2)-rcor(ndum+1)) &
               /(rcor(ndum+1)-rcor(ndum))*(O_50(ndum+1)-O_50(ndum))

  PHI_50(1)=(rcor(3)**2*PHI_50(2)-rcor(2)**2*PHI_50(3))/(rcor(3)**2-rcor(2)**2)
  PHI_50(ndum+2)=PHI_50(ndum+1)+(rcor(ndum+2)-rcor(ndum+1)) &
               /(rcor(ndum+1)-rcor(ndum))*(PHI_50(ndum+1)-PHI_50(ndum))

  Ti_50=Ti_50/1000._wp      !in Kev
  Te_50=Te_50/1000._wp      !in Kev

!!!!!! DEBUG !!!!!!!
!  if(mype==0)then
!     open(unit=987,file='DEBUG_PROFILES.dat',status='unknown')
!     do i=1,ndum+2
!        write(987,'(6(1pe14.5))')rcor(i),ne_50(i),Ti_50(i),Te_50(i),O_50(i),PHI_50(i)
!     enddo
!     close(987)
!  endif
!!!!!!! DEBUG !!!!!!

  deallocate(dne_51,dTi_51,dTe_51)

end subroutine read_profiles

  FUNCTION fln(pdum)
      use precision
      use profile_array
      implicit  none
      real(wp) fln,pdum, fdown

      fdown=abs((pdum-rcor2(jj1+1))/(f1*(rcor2(jj1+1)-rcor2(jj1))))
      fdown=dexp(-f2*fdown**f3)
      fln=dne_51(jj1+1)*fdown

      RETURN
   end function fln

   FUNCTION flTi(pdum)
      use precision
      use profile_array
      implicit  none
      real(wp) flTi, pdum, fdown

      fdown=abs((pdum-rcor2(jj1+1))/(f1*(rcor2(jj1+1)-rcor2(jj1))))
      fdown=dexp(-f2*fdown**f3)
      flTi=dTi_51(jj1+1)*fdown

      RETURN
   end function flTi

   FUNCTION flTe(pdum)
      use precision
      use profile_array
      implicit  none
      real(wp) flTe, pdum, fdown

      fdown=abs((pdum-rcor2(jj1+1))/(f1*(rcor2(jj1+1)-rcor2(jj1))))
      fdown=dexp(-f2*fdown**f3)
      flTe=dTe_51(jj1+1)*fdown

      RETURN
   end function flTe

   FUNCTION frn(pdum)
      use precision
      use profile_array
      implicit  none
      real(wp) frn, pdum, fdown

      fdown=abs((pdum-rcor2(jj2-1))/(f1*(rcor2(jj2-1)-rcor2(jj2))))
      fdown=dexp(-f2*fdown**f3)
      frn=dne_51(jj2-1)*fdown

      RETURN
   end function frn

   FUNCTION frTi(pdum)
      use precision
      use profile_array
      implicit  none
      real(wp) frTi, pdum, fdown

      fdown=abs((pdum-rcor2(jj2-1))/(f1*(rcor2(jj2-1)-rcor2(jj2))))
      fdown=dexp(-f2*fdown**f3)
      frTi=dTi_51(jj2-1)*fdown

      RETURN
   end function frTi

   FUNCTION frTe(pdum)
      use precision
      use profile_array
      implicit  none
      real(wp) frTe, pdum, fdown

      fdown=abs((pdum-rcor2(jj2-1))/(f1*(rcor2(jj2-1)-rcor2(jj2))))
      fdown=dexp(-f2*fdown**f3)
      frTe=dTe_51(jj2-1)*fdown

      RETURN
   end function frTe




subroutine read_profiles_NSTX(profile_data_name)

   use precision
   use mpi_defs
   use global_parameters, only:mype
   use profile_array

   implicit none

   real(wp), dimension(:,:), allocatable :: worka
   integer i,j,ierr
   character(100) chdum
   character(25)  profile_data_name

   ! Only the root process reads the file containing the profiles.
   ! We first read the number of radial grid points there are for each profile.
   if(mype == 0)then
      open(unit=3,file=profile_data_name)
      read(3,*) chdum
      read(3,*) chdum
      read(3,*) chdum
      read(3,*) chdum
      read(3,*) chdum
      read(3,*) chdum
      read(3,*) chdum

      read(3,10) ndum
      write(*,*)'in read_profiles_NSTX: Number of grid points ndum=',ndum
   endif
   10 format(41x,I3)

   ! The value of "ndum" is broadcast to the other processes for memory allocation
   call mpi_bcast(ndum,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   allocate(worka(ndum,5),rcor(ndum+2),ne_50(ndum+2),Ti_50(ndum+2),&
      Te_50(ndum+2),O_50(ndum+2))

   ! The root process continues to read the file and store the relevant profiles
   ! in the wrok array "worka".
   if(mype == 0)then
      do i=1,5
         read(3,20) worka(1:ndum,i)
         read(3,*) chdum
      end do
      close(3)
   endif
   20 format(ES16.8)

   ! The root process broadcast the work array "worka"
   call mpi_bcast(worka,5*ndum,mpi_Rsize,0,MPI_COMM_WORLD,ierr)

   rcor(1:ndum)=worka(1:ndum,1)     ! normalized radial coordinate
   Te_50(1:ndum)=worka(1:ndum,2)    ! in KeV
   ne_50(1:ndum)=worka(1:ndum,3)    ! in m^-3
   Ti_50(1:ndum)=worka(1:ndum,4)    ! in KeV
   O_50(1:ndum)=worka(1:ndum,5)     ! in rad/s

   ne_50=ne_50/1.0e6_wp      !in cm^-3
   ndum=ndum-2

   deallocate(worka)

end subroutine read_profiles_NSTX




subroutine spl_init

   use spl_profiles
   use profile_array

   implicit none

   call ezspline_init(spn,ndum+2,(/1,0/),ier)
   spn%x1=rcor
   spn%bcval1min=0.0d0
   call ezspline_setup(spn,ne_50,ier)

   call ezspline_init(spTi,ndum+2,(/1,0/),ier)
   spTi%x1=rcor
   spTi%bcval1min=0.0d0
   call ezspline_setup(spTi,Ti_50,ier)

   call ezspline_init(spTe,ndum+2,(/1,0/),ier)
   spTe%x1=rcor
   spTe%bcval1min=0.0d0
   call ezspline_setup(spTe,Te_50,ier)

   call ezspline_init(spO,ndum+2,(/1,0/),ier)
   spO%x1=rcor
   spO%bcval1min=0.0d0
   call ezspline_setup(spO,O_50,ier)

   call ezspline_init(spPHI,ndum+2,(/1,0/),ier)
   spPHI%x1=rcor
   spPHI%bcval1min=0.0d0
   call ezspline_setup(spPHI,PHI_50,ier)

end subroutine spl_init




subroutine spl_init2

   use spl_profiles
   use profile_array

   implicit none

   call ezspline_init(spn,ndum,(/1,0/),ier)
   spn%x1=rcor2
   spn%bcval1min=0.0d0
   call ezspline_setup(spn,ne_51,ier)

   call ezspline_init(spTi,ndum,(/1,0/),ier)
   spTi%x1=rcor2
   spTi%bcval1min=0.0d0
   call ezspline_setup(spTi,Ti_51,ier)

   call ezspline_init(spTe,ndum,(/1,0/),ier)
   spTe%x1=rcor2
   spTe%bcval1min=0.0d0
   call ezspline_setup(spTe,Te_51,ier)

   ! call ezspline_init(spO,ndum+2,(/1,0/),ier)
   ! spO%x1=rcor
   ! spO%bcval1min=0.0d0
   ! call ezspline_setup(spO,O_50,ier)

end subroutine spl_init2




subroutine spl_free

   use spl_profiles
   use profile_array

   call ezspline_free(spn,ier)
   call ezspline_free(spTi,ier)
   call ezspline_free(spTe,ier)
   call ezspline_free(spO,ier)
   call ezspline_free(spPHI,ier)

end subroutine spl_free




subroutine get_profiles
   ! This routine must be called by a single processor since it writes to a file
   ! that has a processor-independent name.

   use precision
   use global_parameters
   use spl_profiles
   use profile_array

   implicit  none
   integer, parameter:: np=307
   real(doubleprec), dimension(np) :: xx,yn,yT,yTe,yo
   integer i
   real(doubleprec) thedum,pdum
   real(wp), external :: density,temperature,temperature_e,domega_da


   open(105,file='show_profile.dat')

   thedum=0.0_wp
   do i=1,np
      xx(i)=1.0d0/real(np-1,doubleprec)*real(i-1,doubleprec)
      yn(i)=density(xx(i),thedum)
      yT(i)=temperature(xx(i),thedum)
      yTe(i)=temperature_e(xx(i),thedum)
      pdum=domega_da(xx(i),thedum,yo(i))
   end do

   !     call ezspline_interp(spn,np,xx,yn,ier)
   !     call ezspline_interp(spTi,np,xx,yT,ier)
   !     call ezspline_interp(spTe,np,xx,yTe,ier)
   !     call ezspline_interp(spO,np,xx,yo,ier)

   yT=yT*Tem_00             ! convert to Kev
   yTe=yTe*Tem_00/m0ms      ! convert to Kev

   do i=1,np
      write(105,1) xx(i), yn(i), yT(i), yTe(i), yo(i)
   end do
   1    format(5(2x,ES13.4))
   close(105)

end subroutine get_profiles


subroutine get_profiles_on_nodes
   ! prepare data for linear interpolation of profiles

   use precision
   use global_parameters
   use spl_profiles
   use profile_array

   implicit  none
   real(doubleprec) xx,dxx
   integer i

   allocate(Ti_node(0:jr),dTi_node(0:jr),ne_node(0:jr),dne_node(0:jr),&
      Te_node(0:jr),dTe_node(0:jr),Iomg_node(0:jr),dIomg_node(0:jr))

   dxx=(a1-a0)/real(jr,doubleprec)
   do i=0,jr
      xx=a0+dxx*real(i,doubleprec)
      call ezspline_interp(spn,xx,ne_node(i),ier)
      call ezspline_derivative(spn,1,xx,dne_node(i),ier)

      call ezspline_interp(spTi,xx,Ti_node(i),ier)
      call ezspline_derivative(spTi,1,xx,dTi_node(i),ier)

      call ezspline_interp(spTe,xx,Te_node(i),ier)
      call ezspline_derivative(spTe,1,xx,dTe_node(i),ier)

      call ezspline_interp(spO,xx,Iomg_node(i),ier)        ! Iomg=omega
      call ezspline_derivative(spO,1,xx,dIomg_node(i),ier) ! dIomg=d(omega)/da

      if(mype==0) write(699,1) xx, ne_node(i), Ti_node(i)*Tem_00, &
         Te_node(i)*Tem_00/m0ms, dne_node(i)/ne_node(i), &
         dTi_node(i)/Ti_node(i), dTe_node(i)/Te_node(i)
   end do

   1    format(7(2x,ES13.4))
   call flush(699)

end subroutine get_profiles_on_nodes
