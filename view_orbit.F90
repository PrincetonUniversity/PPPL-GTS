Subroutine view_orbit

  use precision
  use global_parameters
  use particle_array
  use field_array
  use esi_interface

  implicit none

  integer i, k, j

  real(wp), dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,en_error

  if(istep==0) then
     open(unit=10,file='orbit.dat')
     open(unit=15,file='en_error.dat')
     open(unit=17,file='weight.dat')
  end if

  do i=1,10

    aaa(i)=zion(1,mi/2+i)
    the(i)=zion(2,mi/2+i)

!    aaa(i)=zion(1,i)
!    the(i)=zion(2,i)
  end do

  call esigetrzb(rj,zj,Bj,aaa,the,10)

  write(10,100) (rj(k), zj(k), k=1,10)
100 format(1x,20(ES12.6,1x))


  en_error(1:10)=0.5_wp*(zion(4,mi/2+1:mi/2+10)*Bj(1:10))**2+zion(6,mi/2+1:mi/2+10)**2*Bj(1:10)
  en_error(1:10)=abs(enion(mi/2+1:mi/2+10)-en_error(1:10))/enion(mi/2+1:mi/2+10)


  write(15,300) istep, (en_error(k), k=1,10)

  write(17,300) istep, (zion(5,mi/2+k), k=1,10)

200 format(1x,10(ES12.6,1x))


300 format(1x,i4,1x,20(ES12.6,1x))

  call flush(10)
  call flush(15)
  call flush(17)

end subroutine view_orbit
