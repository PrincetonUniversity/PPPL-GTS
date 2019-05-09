Subroutine check_eq

  use precision
  use global_parameters
  use particle_array
  use esi_interface

  implicit none
  real(wp), dimension(100) :: Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,Jco
  integer n_tgrid, n_rgrid,n_t
  integer i, j, k, dp

!  print*,'check_eq', mype

  open(unit=11,file='eq_surface.dat')
  open(unit=12,file='eq_theta.dat')
  open(unit=16,file='jacobian.dat')

  n_tgrid=500
  n_rgrid=10
  n_t=20

  do i=0,n_tgrid
    do j=0,n_rgrid
      aaa(j+1)=1.0_wp/real(n_rgrid,wp)*real(n_rgrid-j,wp)
      the(j+1)=2.0_wp*pi/real(n_tgrid,wp)*real(i,wp)
    end do
    aaa(n_rgrid+1)=1.0e-7_wp
  
  call esigetrzb(rj,zj,Bj,aaa,the,n_rgrid+1)
  write(11,100) rj(1), zj(1),rj(2), zj(2),rj(3), zj(3),rj(4), zj(4), &
                rj(5), zj(5),rj(6), zj(6),rj(7), zj(7),rj(8), zj(8), &
                rj(9), zj(9),rj(10), zj(10),rj(11), zj(11)

! check uniformity of Jacobian
!  call esiget2dfunctions(aaa,the,n_rgrid+1)  ! replaced by:
  call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,n_rgrid+1)
     do j=1,n_rgrid+1
     Jco(j)=rj(j)*(zaj(j)*rqj(j)-zqj(j)*raj(j))
     end do
  write(16,100) (the(k),Jco(k),k=1,n_rgrid+1)

  end do

  do  i=0,n_rgrid
    do j=0,n_t
      if(i.eq.n_rgrid) then
      aaa(j+1)=1.0e-7_wp
      else
      aaa(j+1)=1.0_wp/real(n_rgrid,wp)*real(n_rgrid-i,wp)
      end if
      the(j+1)=2.0_wp*pi/real(n_t,wp)*real(j,wp)
    end do
 
  call esigetrzb(rj,zj,BB,aaa,the,n_t+1)
  write(12,200) rj(1), zj(1),rj(2), zj(2),rj(3), zj(3),rj(4), zj(4), &
                rj(5), zj(5),rj(6), zj(6),rj(7), zj(7),rj(8), zj(8), &
                rj(9), zj(9),rj(10), zj(10),&
                rj(11), zj(11),rj(12), zj(12),rj(13), zj(13),rj(14), zj(14), &
                rj(15), zj(15),rj(16), zj(16),rj(17), zj(17),rj(18), zj(18), &
                rj(19), zj(19),rj(20), zj(20),rj(21), zj(21)
  end do

100 format(1x,22(d13.6,1x))
200 format(1x,42(d13.6,1x))

  close(11)
  close(12)
  close(16)
end subroutine check_eq
  
  
