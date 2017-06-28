Subroutine view_distribution(zero_or_one)

  use global_parameters
  use particle_array
  use esi_interface

  implicit none
  integer zero_or_one, i
  
  if(zero_or_one==0) then
    open(unit=1,file='dis0.dat')
    open(unit=81,file='tz0.dat')
    call esiget2dfunctions(zion(1,:),zion(2,:),mi)
    do i=1,mi
      write(1,1000) rr(i), z(i)
      if(zion(3,i).lt.zetamin.or.zion(3,i).gt.zetamax) write(81,1000) zion(3,i),zion(2,i) 
    end do
  else
    open(unit=2,file='dis1.dat')
    open(unit=82,file='tz1.dat')
    call esiget2dfunctions(zion(1,:),zion(2,:),mi)
    do i=1,mi
      write(2,1000) rr(i), z(i)
      if(zion(3,i).lt.zetamin.or.zion(3,i).gt.zetamax) write(82,1000) zion(3,i),zion(2,i) 
    end do
  end if


1000 format(1x,2(ES13.6,1x))

end subroutine view_distribution 
