program rewrite  
  integer,parameter :: num_mode=8
  integer  i,j,k,mall,i1,i2,i3,i4,i5
  real fdum
  character(len=31) cdum

           open(10,file='history.out')            
           open(20,file='rewrite.out')            

           read(10,101)  i1,i2,i3,i4,i5
           write(20,101) i1,i2,i3,i4,i5
           read(10,102) fdum
           write(20,102) fdum

  do i1=2,2000,2
     read(10,*) cdum, i5, cdum
     do i2=1,25
     read(10,102) fdum
     write(20,102) fdum
     enddo

     do i=1,2
        do j=1,num_mode
           do k=1,2
              read(10,103)  fdum,k,j,i
              write(20,102) fdum
           enddo
        enddo
     enddo

  enddo

101 format(i6)
102 format(e10.4)
103 format(e10.4,2x,i4,2x,i4,2x,i4)

end program rewrite
