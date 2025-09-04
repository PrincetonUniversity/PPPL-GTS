program rewrite_snap
  
  integer,parameter :: mbin_psi=5, mbin_u=32, mpsi=64,jm=320,numberpe=32, &
                       num_mode=8,m_poloidal=9
  integer i,j,jj,jt,k,kz,m,nsnap,n1,n2,n3,n4,n5,n6,m1,m2,m3,m4,m5,m6
  real a1,a2,c1,c2
  real,dimension(0:mpsi) :: marker,fflows,dflows,ftem,dtem
  real,dimension(mbin_u,mbin_psi) :: ubin,dubin,pbin,dpbin,ebin,debin

           open(10,file='snap121.out')
           open(11,file='SAMPLE/snap121.out')
           open(60,file='snap400.out')
           open(61,file='snap500.out')

     read(10,102) a1,a2
     read(11,102) c1,c2
     write(60,102) a1,a2
     write(61,102) c1,c2

     read(10,101)n1,n2,n3,n4,n5,n6
     read(11,101)m1,m2,m3,m4,m5,m6
     write(60,101)n1,n2,n3,n4,n5,n6
     write(61,101)m1,m2,m3,m4,m5,m6
     
     read(10,101)n1
     read(11,101)m1
     write(60,101)n1
     write(61,101)m1

     do i=1,mbin_u
        read(10,102) a1
        read(11,102) c1
        write(60,102) a1
        write(61,102) c1
     enddo

     do i=1,mbin_u
        read(10,102) a1
        read(11,102) c1
        write(60,102) a1
        write(61,102) c1
     enddo

     do i=1,mbin_u
        read(10,102) a1
        read(11,102) c1
        write(60,102) a1
        write(61,102) c1
     enddo

! velocity space distribution function
     do j=1,mbin_psi
        read(10,102)(ubin(i,j),i=1,mbin_u)
        read(10,102)(dubin(i,j),i=1,mbin_u)
        read(10,102)(pbin(i,j),i=1,mbin_u)
        read(10,102)(dpbin(i,j),i=1,mbin_u)
        read(10,102)(ebin(i,j),i=1,mbin_u)
        read(10,102)(debin(i,j),i=1,mbin_u)

        write(60,102)(ubin(i,j),i=1,mbin_u)
        write(60,102)(dubin(i,j),i=1,mbin_u)
        write(60,102)(pbin(i,j),i=1,mbin_u)
        write(60,102)(dpbin(i,j),i=1,mbin_u)
        write(60,102)(ebin(i,j),i=1,mbin_u)
        write(60,102)(debin(i,j),i=1,mbin_u)

        read(11,102)(ubin(i,j),i=1,mbin_u)
        read(11,102)(dubin(i,j),i=1,mbin_u)
        read(11,102)(pbin(i,j),i=1,mbin_u)
        read(11,102)(dpbin(i,j),i=1,mbin_u)
        read(11,102)(ebin(i,j),i=1,mbin_u)
        read(11,102)(debin(i,j),i=1,mbin_u)

        write(61,102)(ubin(i,j),i=1,mbin_u)
        write(61,102)(dubin(i,j),i=1,mbin_u)
        write(61,102)(pbin(i,j),i=1,mbin_u)
        write(61,102)(dpbin(i,j),i=1,mbin_u)
        write(61,102)(ebin(i,j),i=1,mbin_u)
        write(61,102)(debin(i,j),i=1,mbin_u)
     enddo

     read(10,101)n1
     read(11,101)m1
     write(60,101)n1
     write(61,101)m1

     do i=1,mpsi,2
        read(10,102) a1
        read(11,102) c1
        write(60,102) a1
        write(61,102) c1
     enddo

     read(10,102)(marker(i),i=1,mpsi,2)
     write(60,102)(marker(i),i=1,mpsi,2)
     read(11,102)(marker(i),i=1,mpsi,2)
     write(61,102)(marker(i),i=1,mpsi,2)

     do i=1,mpsi,2
        read(10,102) a1
        read(11,102) c1
        write(60,102) a1
        write(61,102) c1
     enddo

     read(10,102)(marker(i),i=1,mpsi,2)
     write(60,102)(marker(i),i=1,mpsi,2)
     read(11,102)(marker(i),i=1,mpsi,2)
     write(61,102)(marker(i),i=1,mpsi,2)

     read(10,102)(marker(i),i=1,mpsi,2)
     write(60,102)(marker(i),i=1,mpsi,2)
     read(11,102)(marker(i),i=1,mpsi,2)
     write(61,102)(marker(i),i=1,mpsi,2)

     read(10,102)(marker(i),i=1,mpsi,2)
     write(60,102)(marker(i),i=1,mpsi,2)
     read(11,102)(marker(i),i=1,mpsi,2)
     write(61,102)(marker(i),i=1,mpsi,2)

     read(10,102)(marker(i),i=1,mpsi,2)
     write(60,102)(marker(i),i=1,mpsi,2)
     read(11,102)(marker(i),i=1,mpsi,2)
     write(61,102)(marker(i),i=1,mpsi,2)

     read(10,102)(marker(i),i=1,mpsi,2)
     write(60,102)(marker(i),i=1,mpsi,2)
     read(11,102)(marker(i),i=1,mpsi,2)
     write(61,102)(marker(i),i=1,mpsi,2)

     read(10,101)n1
     read(11,101)m1
     write(60,101)n1
     write(61,101)m1

     do j=0,jm
        jj=j
        if(j>jm/2)jj=j-jm
        do i=1,mpsi,2
        read(10,102) a1
        read(11,102) c1
        write(60,102) c1
        write(61,102) a1
        enddo
     enddo

     do j=0,jm
        jj=j
        if(j>jm/2)jj=j-jm
        do i=1,mpsi,2
        read(10,102) a1
        read(11,102) c1
        write(60,102) c1
        write(61,102) a1
        enddo
     enddo
     
     do j=0,jm
        jj=j
        if(j>jm/2)jj=j-jm
        do i=1,mpsi,2
        read(10,102) a1
        read(11,102) c1
        write(60,102) a1
        write(61,102) c1
        enddo
     enddo

     read(10,101)n1
     read(11,101)m1
     write(60,101)n1
     write(61,101)m1

     do j=1,jm
        do kz=1,numberpe
           read(10,102) a1
           read(11,102) c1
           write(60,102) a1
           write(61,102) c1
        enddo
     enddo

     read(10,101)n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
     read(11,101)m1,m2,m3,m4,m5,m6,m7,m8,m9,m10
     write(60,101)n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
     write(61,101)m1,m2,m3,m4,m5,m6,m7,m8,m9,m10

     do i=1,mpsi,2
        do kz=1,num_mode
           do j=1,m_poloidal
             read(10,102) a1
             read(11,102) c1
             write(60,102) a1
             write(61,102) c1
           enddo
        enddo
     enddo

101 format(i6)
102 format(e10.4)
  
end program rewrite_snap

