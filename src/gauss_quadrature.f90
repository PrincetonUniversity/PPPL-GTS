!       program main
!          integer :: m
!          real*8 :: a1,a2,delta,intgall
!          real*8, allocatable :: x(:), y(:)
! 
!          !integrate from a1 to a2 interval
!          a1=2
!          a2=5
! 
!          !divide [a1,a2] into m of integration subintervals
!          m=30
! 
!       !the coordinate of grid [x0,x1,x2,x3, ... ,xm-1,xm]
!       !y is the function to be interpolated
!       allocate(x(0:m), y(0:m))
!       do i=0,m
!          delta=(a2-a1)/m
!          x(i)=a1+i*delta
!          y(i)=x(i)*x(i)
!          !   write(*,'(A,i, 2(2x,f8.2))') "function at x =", i, x(i), y(i)
!       enddo
!       print *
! 
!          call gts_integration(m,a1,a2,x,y,intgall)
! 
!          print *
!          !   write(*, '(A,f6.2,2x,f6.2,2x,f8.2)') "int x^2 dx =", a1,a2,intgall
!  
!          deallocate(x, y)
!        end program main

      subroutine gts_integration(m,a1,a2,x,y,intgall)
         integer :: m
         real*8 :: a1,a2,intgall
         real*8 :: x(0:m), y(0:m)

         integer :: n,i,j
         real*8, allocatable :: wi(:), xi(:) !qudrature weight and location
         real*8, allocatable :: rf(:) !funtion to be integrated
         real*8 :: xa,xb,intg,jacob, rx,delta
         real*8, allocatable :: b(:), c(:), d(:)
         real*8 :: xx, eval_at_xx


         !9 quadrature points
         n=9
         allocate(wi(n))
         allocate(xi(n))
         allocate(rf(n))
         call gauss_quadrature_9(n,wi,xi)



         !================= =================!
         !cubic spline
      allocate(b(0:m), c(0:m), d(0:m))

      ! call spline.f to calculate b c d of cubic spline
      call spline (m+1, x, y, b, c, d)
      !   do i=1,m
      !      write(*,'(A,i2.2,3(2x,f8.2))') "b c d at x =", i, b(i), c(i), d(i)
      !   enddo
      ! print *
         !================= =================!



         intgall=0.
         do i=1,m
            delta=x(i)-x(i-1)
            xa=x(i-1)
            xb=x(i)
            jacob=(xb-xa)/2.

            ! examplef=x^2
            do j=1,n
               rx=(xb+xa)/2. + (xb-xa)/2. * xi(j)

         !================= =================!
         ! call seval.f to evaluate the interpolate function at a given x
         xx=rx
         call seval(m+1, xx, x, y, b, c, d, eval_at_xx)
               rf(j)=eval_at_xx
         !   write(*,'(A,2(2x,f8.2))') "function at xx =", xx, eval_at_xx
         !================= =================!

            enddo

            intg=0.
            do j=1,n
               intg=intg+rf(j)*wi(j)*jacob
            enddo

            !   write(*, '(A,f6.2,2x,f6.2,2x,f8.2)') "int x^2 dx =", xa,xb,intg
            intgall=intgall+intg
         enddo

         !   write(*, *)
         !   write(*, '(A,f6.2,2x,f6.2,2x,f8.2)') "int x^2 dx =", a1,a2,intgall

         deallocate(wi,xi,rf)
         deallocate(b, c, d)

      end subroutine gts_integration


      subroutine gauss_quadrature_9(n,wi,xi)
              !https://pomax.github.io/bezierinfo/legendre-gauss.html
         integer :: n
         real*8 :: wi(n), xi(n)

         !i      weight - wi     abscissa - xi
         wi(1)=0.3302393550012598;       xi(1)= 0.0000000000000000
         wi(2)=0.1806481606948574;       xi(2)=-0.8360311073266358
         wi(3)=0.1806481606948574;       xi(3)= 0.8360311073266358
         wi(4)=0.0812743883615744;       xi(4)=-0.9681602395076261
         wi(5)=0.0812743883615744;       xi(5)= 0.9681602395076261
         wi(6)=0.3123470770400029;       xi(6)=-0.3242534234038089
         wi(7)=0.3123470770400029;       xi(7)= 0.3242534234038089
         wi(8)=0.2606106964029354;       xi(8)=-0.6133714327005904
         wi(9)=0.2606106964029354;       xi(9)= 0.6133714327005904
      end subroutine gauss_quadrature_9



      subroutine gauss_quadrature_25(a,b,intg)
              !https://pomax.github.io/bezierinfo/legendre-gauss.html
         integer :: n, i
         real*8 :: a, b, intg
         real*8 :: jacob, rx
         real*8, allocatable :: wi(:), xi(:), rf(:)

         n=25
         allocate(wi(n))
         allocate(xi(n))
         allocate(rf(n))

         !i      weight - wi     abscissa - xi
         wi( 1)=0.1231760537267154;       xi(1)=0.0000000000000000
         wi( 2)=0.1222424429903100;       xi(2)=-0.1228646926107104
         wi( 3)=0.1222424429903100;       xi(3)=0.1228646926107104
         wi( 4)=0.1194557635357848;       xi(4)=-0.2438668837209884
         wi( 5)=0.1194557635357848;       xi(5)=0.2438668837209884
         wi( 6)=0.1148582591457116;       xi(6)=-0.3611723058093879
         wi( 7)=0.1148582591457116;       xi(7)=0.3611723058093879
         wi( 8)=0.1085196244742637;       xi(8)=-0.4730027314457150
         wi( 9)=0.1085196244742637;       xi(9)=0.4730027314457150
         wi(10)=0.1005359490670506;       xi(10)=-0.5776629302412229
         wi(11)=0.1005359490670506;       xi(11)=0.5776629302412229
         wi(12)=0.0910282619829637;       xi(12)=-0.6735663684734684
         wi(13)=0.0910282619829637;       xi(13)=0.6735663684734684
         wi(14)=0.0801407003350010;       xi(14)=-0.7592592630373576
         wi(15)=0.0801407003350010;       xi(15)=0.7592592630373576
         wi(16)=0.0680383338123569;       xi(16)=-0.8334426287608340
         wi(17)=0.0680383338123569;       xi(17)=0.8334426287608340
         wi(18)=0.0549046959758352;       xi(18)=-0.8949919978782753
         wi(19)=0.0549046959758352;       xi(19)=0.8949919978782753
         wi(20)=0.0409391567013063;       xi(20)=-0.9429745712289743
         wi(21)=0.0409391567013063;       xi(21)=0.9429745712289743
         wi(22)=0.0263549866150321;       xi(22)=-0.9766639214595175
         wi(23)=0.0263549866150321;       xi(23)=0.9766639214595175
         wi(24)=0.0113937985010263;       xi(24)=-0.9955569697904981
         wi(25)=0.0113937985010263;       xi(25)=0.9955569697904981

         jacob=(b-a)/2.

         ! examplef=x^2
         do i=1,n
            rx=(b+a)/2. + (b-a)/2. * xi(i)
            rf(i)=rx*rx
         enddo

         intg=0.
         do i=1,n
            intg=intg+rf(i)*wi(i)*jacob
         enddo

         deallocate(wi)
         deallocate(xi)
         deallocate(rf)
      end subroutine gauss_quadrature_25
