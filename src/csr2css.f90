subroutine csrcsc ( n, job, ipos, a, ja, ia, ao, jao, iao )

   !*******************************************************************************
   !
   !! CSRCSC converts Compressed Sparse Row to Compressed Sparse Column.
   !
   !  Discussion:
   !
   !    This is essentially a transposition operation.
   !
   !    It is NOT an in-place algorithm.
   !
   !    This routine transposes a matrix stored in a, ja, ia format.
   !
   !  Modified:
   !
   !    07 January 2004
   !
   !  Author:
   !
   !    Youcef Saad
   !
   !  Parameters:
   !
   !    Input, integer N, the order of the matrix.
   !
   !    Input, integer JOB, indicates whether or not to fill the values of the
   !    matrix AO or only the pattern (IA, and JA).  Enter 1 for yes.
   !
   ! ipos  = starting position in ao, jao of the transposed matrix.
   !         the iao array takes this into account (thus iao(1) is set to ipos.)
   !         Note: this may be useful if one needs to append the data structure
   !         of the transpose to that of A. In this case use
   !                call csrcsc (n,1,n+2,a,ja,ia,a,ja,ia(n+2))
   !        for any other normal usage, enter ipos=1.
   !
   !    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
   !    Compressed Sparse Row format.
   !
   !    Output, real AO(*), JAO(*), IAO(N+1), the matrix in CSC
   !    Compressed Sparse Column format.
   !
   implicit none

   integer n

   real ( kind = 8 ) a(*)
   real ( kind = 8 ) ao(*)
   integer i
   integer ia(n+1)
   integer iao(n+1)
   integer ipos
   integer j
   integer ja(*)
   integer jao(*)
   integer job
   integer k
   integer next
   !
   !  Compute lengths of rows of A'.
   !
   iao(1:n+1) = 0

   do i = 1, n
      do k = ia(i), ia(i+1)-1
         j = ja(k) + 1
         iao(j) = iao(j) + 1
      end do
   end do
   !
   !  Compute pointers from lengths.
   !
   iao(1) = ipos
   do i = 1, n
      iao(i+1) = iao(i) + iao(i+1)
   end do
   !
   !  Do the actual copying.
   !
   do i = 1, n
      do k = ia(i), ia(i+1)-1
         j = ja(k)
         next = iao(j)
         if ( job == 1 ) then
            ao(next) = a(k)
         end if
         jao(next) = i
         iao(j) = next + 1
      end do
   end do
   !
   !  Reshift IAO and leave.
   !
   do i = n, 1, -1
      iao(i+1) = iao(i)
   end do
   iao(1) = ipos

   return
   end

