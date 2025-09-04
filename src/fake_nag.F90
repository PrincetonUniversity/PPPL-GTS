!     D01FCF attempts to evaluate a multi-dimensional integral
!     with constant and finite limites, to a specified relative
!     accuracy using an adaptive subdivision strategy
!
      subroutine  d01fcf(ndim,A,B,minpts,maxpts,functn,eps,                 &
     &                   acc, lenwrk, wrkstr, finval, ifail )
      implicit none
      integer ndim, minpts, maxpts, lenwrk, ifail
      real A(ndim), B(ndim), eps, acc, wrkstr(lenwrk)
      real finval

      real functn
      external functn              

      interface
        subroutine adapt(ndim,A,B,minpts,maxpts,functn,eps,                 &
     &        relerr,lenwrk,wrkstr,finest, ifail)
        implicit none
        integer ndim,minpts,maxpts,lenwrk, ifail
        real A(ndim),B(ndim), eps, acc, wrkstr(lenwrk)
        real finest, relerr

        real functn
        external functn
        end subroutine adapt
      end interface

      call adapt(ndim,A,B,minpts,maxpts,functn,eps,                         &
     &            acc, lenwrk, wrkstr, finval, ifail )

      return
      end

!     =====================================================
!     C05ADF locates a zero of a continuous function in a given interval
!     by a combination of the methods of linear interpolation,
!     extrapolation and bisection
!     
!     abs(x-root) < eps
!     abs( f(x) ) < eta
!

      subroutine C05ADF(a,b, eps, eta, F, x, ifail )
      implicit none
      integer ifail
      real a,b,eps,eta,  x
      real F
      external F

      logical is_error
      integer iflag
      real re, ae, aa, bb

!        RE    - RELATIVE ERROR USED FOR RW IN THE STOPPING CRITERION.
!                IF THE REQUESTED RE IS LESS THAN MACHINE PRECISION,
!                THEN RW IS SET TO APPROXIMATELY MACHINE PRECISION.
!        AE    - ABSOLUTE ERROR USED IN THE STOPPING CRITERION.  IF THE
!                GIVEN INTERVAL (B,C) CONTAINS THE ORIGIN, THEN A
!                NONZERO VALUE SHOULD BE CHOSEN FOR AE.
!        IFLAG - RETURNS A STATUS OF THE RESULTS INDICATING WHICH
!                OF THE FOLLOWING CONDITIONS HOLD.
!                    A - ABS(B-C) .LE. 2.*(RW*ABS(B)+AE)
!                    B - F(B) * F(C) .LT. 0.
!                    C - ABS(F(B)) .LE. ABS(F(C))
!                    D - ABS(F(B   )) .LE. MAX(ABS(F(B  )),ABS(F(C  )))
!                               OUT                   IN          IN
!                    E - NUMBER OF EVALUATIONS OF F(X) .LE. 500
!               =1 INDICATES NORMAL CASE.  ALL CONDITIONS ABOVE HOLD.
!               =2 INDICATES F(B) = 0.  CONDITION A MAY NOT HOLD.
!               =3 INDICATES CONDITIONS A, B, C, AND E HOLD BUT D DOES
!                  NOT.  (B,C) PROBABLY CONTAINS A SINGULAR POINT OF F.
!               =4 INDICATES CONDITIONS A AND E HOLD BUT B DOES NOT.
!                  A LOCAL MINIMUM OF F(X) IN (B,C) MAY HAVE BEEN FOUND.
!               =5 INDICATES SEARCH WAS ABORTED WHEN CONDITION E FAILED.

      interface
        subroutine zeroin(F,b,c,re,ae,iflag)
        implicit none
        integer iflag
        real b,c,re,ae
        real F
        external F
        end subroutine
      end interface


      is_error = (eps < 0.0) .or. (a .eq. b).or.f(a)*f(b) > 0.0
      if (is_error) then
         ifail = 1
         return
      endif
    
      re = 0.25d0*eps/abs(b-a)
      ae = epsilon(ae)
      aa = a
      bb = b
      call zeroin(F, aa,bb, re,ae,iflag)

      if(iflag==2)then
         x = aa
      else
         x = (aa + bb)/2.0d0
      endif

!       -----------
!       normal case
!       -----------
      if (iflag.eq.1 .or. iflag.eq.2) then
        ifail = 0
        return
      endif

!     ---------------------------------
!     for all other cases
!     indicate serious error and return
!     ---------------------------------
      print*,'zeroin returns iflag = ',iflag
      print*,'aa,bb,re,ae ',aa,bb,re,ae

      ifail  = 4
      return
      end

!     =====================================================
!     E01BAF determines a cubic-spline interpolant to a 
!     given set of data
!
      subroutine E01BAF(m, x, y, lamda, c, lck, wrk, lwrk, ifail)
      implicit none
      integer m, lck, lwrk, ifail
      real x(m), y(m), lamda(lck), c(lck), wrk(lwrk)
!
!     derivative end conditions are not imposed
!     1st 4 and last 4 knots of B splines set to x(1) and x(m) 

      logical is_error, is_ok
      integer i

      integer ibcl, ibcr, n, k, kntopt
      real fbcl, fbcr

      interface
       subroutine BINT4(x,y,ndata,ibcl,ibcr,fbcl,fbcr,kntopt,               &
     &                  t,bcoef,n,k,w)
       implicit none
       integer ndata,ibcl,ibcr,kntopt,n,k
       real x(ndata),y(ndata),w(*), t(*),bcoef(*), fbcl,fbcr
       end subroutine BINT4
      end interface 

!     ------------------
!     checking arguments
!     ------------------
      is_error = (m.lt.4).or.(lck.lt.m+4).or.(lwrk.lt.6*M+16)
      if (is_error) then
          ifail = 1
          return
      endif

      is_ok = .false.
      do i=1,m-1
        is_ok = x(i).lt.x(i+1)
        if (.not.is_ok) exit
      enddo
      is_error = .not. is_ok
      if (is_error) then
         ifail = 2
         return
      endif

     





      kntopt=1
      ibcl = 0
      ibcr = 0
      fbcl = 0.0
      fbcr = 0.0
      call  bint4( x,y,m, ibcl, ibcr, fbcl,fbcr, kntopt,                    &
     &             lamda, c, n,k,wrk)

      ifail = 0

      return
      end
!     =====================================================
!     E02BBF  evaluates a cubic spline from B-spline representation
!
      subroutine E02BBF(ncap7, lamda, c, x, s, ifail)
      implicit none
      integer ncap7, ifail
      real lamda(ncap7), c(ncap7), x, s

      interface
        real function bvalu( t, bcoef, n, k, ideriv, x, inbv, work)
        implicit none
        integer n,k,ideriv,inbv
        real t(*), bcoef(*), x, work(*)
        end function bvalu
      end interface


      logical is_error, is_ok

      integer k, n, ideriv
      real work(3*7)
      integer inbv
      logical isfirst
      save inbv, isfirst

      data isfirst/.true./

      s = 0.0
      is_error = (ncap7 .lt. 8)
      if (is_error) then
        ifail = 2
        return
      endif

      is_ok = (lamda(4).le. x) .and. (x.le.lamda(ncap7-3))
      is_error = .not.is_ok
      if (is_error) then
         ifail = 1
         return
      endif





      if (isfirst) then
        isfirst = .false.
        inbv = 1
      endif



      k = 7
      n = ncap7 - 7
      ideriv = 0
      s = bvalu( lamda, c, n, k, ideriv, x, inbv, work )

      ifail = 0
      return
      end 

!     =====================================================
!     D01AHF computes a definite integral over a finite range to a specified 
!     relative accuracy using a method described by Patterson

      real function d01ahf_new(a,b, epsr, npts, relerr, F, nlimit_in, ifail )
      implicit none
      integer npts, nlimit_in, ifail
      real a, b, epsr, relerr
      real F
      external F

      integer nlimit
      logical is_error

      real eps, r, e
      integer iflag, kf
      real*8 fmin,fmax
      logical rst 
      integer nmax
      parameter(nmax=10000)
      real*8 w(nmax,6)
      integer nint

      interface
        subroutine q1dax(F,a,b,eps,r,e,nint,rst,w,nmax,                    &
     &                   fmin,fmax,kf,iflag)
        integer nint,nmax,kf,iflag
        real*8 a,b,eps,r,e,fmin,fmax
        logical rst
        real*8 F
        real*8 w(nmax,6)
        external F
        end subroutine

        subroutine q1db(F,a,b,eps,r,e,kf,iflag)
        implicit none
        integer kf, iflag
        real*8 a,b,eps,r,e
        real*8 F
        external F
        end subroutine
      end interface
      ifail = 0
      nlimit = nlimit_in
      if (nlimit_in .le. 0) then
        nlimit = 10000
      endif

      print *, "   nlimit =",nlimit
      is_error = (epsr .le. 0.0)
      if (is_error) then
         ifail = 3
         write(*,*) 'd01ahf_new: invalid epsr = ', epsr
         return
      endif
      
      eps = epsr
      nint = 1
      rst = .false.
      !!print *, 'before q1dax',a,b,eps,r
      call  q1dax(F,a,b,eps,r,e,nint,rst,w,nmax,fmin,fmax,kf,iflag)
      !!print *, 'after q1dax'
      npts = kf

      if (abs(r).gt.epsilon(r)) then
        relerr = abs(e/r)
      else
        relerr = abs(e)
      endif

      if ((0.le.iflag).and.(iflag.le.3)) then
!       ------------------
!       normal termination
!       ------------------
        ifail = 0
      else
        write(*,*) 'd01ahf_new: q1da return iflag ',iflag
        ifail = iflag
      endif

          
      d01ahf_new = r


      return
      end

!************************************************************************
!*     D01AHF simulate NAG'S D01AHF 
!*     20250902
!************************************************************************

      real FUNCTION D01AHF(A,B,EPR,NPTS,RELERR,F,NL,IFAIL)

!     INPUT ARGUMENTS
!     ----- ----------
!     A,B     -  LOWER AND UPPER INTEGRATION LIMITS.
!     EPR     -  REQUIRED RELATIVE ACCURACY.
!     NL      -  APPROXIMATE LIMIT ON NUMBER OF INTEGRAND EVALUATIONS.
!                IF SET NEGATIVE OR ZERO THE DEFAULT IS 10000.
!     F       -  THE USER NAMED AND PREPARED FUNCTION  F(X)
!                GIVES THE VALUE OF THE INTEGRAND AT X.
!     IFAIL      INTEGER VARIABLE
!             - 0  FOR HARD FAIL REPORT
!             - 1  FOR SOFT FAIL REPORT
!
!     OUTPUT ARGUMENTS
!     ------ ----------
!     NPTS    -  NUMBER OF INTEGRAND EVALUATIONS USED IN OBTAINING THE RESULT.
!     RELERR  -  ROUGH ESTIMATE OF RELATIVE ACCURACY ACHIEVED.
!     IFAIL   -  VALUE INDICATES THE OUTCOME OF THE INTEGRATION -
!                IFAIL  = 0  CONVERGED
!                IFAIL  = 1  INTEGRAND EVALUATIONS EXCEEDED  NL.
!                            THE RESULT WAS OBTAINED BY CONTINUING
!                            BUT IGNORING ANY NEED TO SUBDIVIDE.
!                            RESULT LIKELY TO BE INACCURATE.
!                IFAIL  = 2  DURING THE SUBDIVISION PROCESS, THE STACK BECAME FULL
!                            (PRESENTLY SET TO HOLD 20 LEVELS OF INFORMATION.  
!                            MAY BE INCREASED BY  ALTERING  ISMAX
!                            AND THE DIMENSIONS OF STACK AND ISTACK).
!                            RESULT IS OBTAINED BY CONTINUING BUT IGNORING CONVERGENCE
!                            FAILURES ON INTERVALS WHICH CANNOT BE ACCOMMODATED ON THE STACKS.
!                            RESULT LIKELY TO BE INACCURATE.
!                IFAIL  = 3  INVALID ACCURACY REQUEST.
!
!     .. Parameters ..
      CHARACTER*7                      SRNAME
      PARAMETER                        (SRNAME='D01AHF')

!     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B, EPR, RELERR
      INTEGER                          IFAIL, NL, NPTS
!     .. Function Arguments ..
      DOUBLE PRECISION                 F
      EXTERNAL                         F

!     Prepare to call dqag

      double precision epsabs,epsrel,qresult,abserr
      integer key,neval,ier,limit,lenw,last
      integer,dimension(:),allocatable :: iwork
      double precision,dimension(:),allocatable :: work

      interface
              subroutine dqag(f,a,b,epsabs,epsrel,key,result,abserr,neval,ier, limit,lenw,last,iwork,work)
                      double precision a,abserr,b,epsabs,epsrel,f,result,work
                      integer ier,iwork,key,last,lenw,limit,lvl,l1,l2,l3,neval
                      dimension iwork(limit),work(lenw)
                      external f
                      end subroutine dqag
      end interface

      epsrel=EPR
      epsabs=epsrel*1.d-5
      key=6
      limit=100
      lenw=4*limit
      allocate(iwork(limit))
      allocate(work(lenw))

      call dqag(F,A,B,epsabs,epsrel,key,qresult,abserr,neval,ier,limit,lenw,last,iwork,work)

      RELERR=abserr
      IFAIL=ier
      NPTS=neval
      deallocate(iwork)
      deallocate(work)

      D01AHF = qresult
      RETURN
      END


!     =====================================================
!     S15AEF returns the value of the error function erf(x)
!     
!   this replacement may cause problem at least on Franklin!!!  Sep. 08
!
      real function s15aef(x,ifail)
      integer ifail
      real x

      ifail = 0

      s15aef = derf( x )
      return
      end



!     =====================================================
!     S14AEF returns the value of the kth derivative of the 
!     psi function for real x, k=0,1,...6 

      real function s14aef(x,k,ifail)
      implicit none
      integer k, ifail
      real x

      logical is_error
      real, parameter ::  eps = epsilon(x)

      is_error = abs(x-int(x)) .lt. eps * int( abs(x) )
      is_error = is_error.or.(k .lt. 0).or.(k.gt.6)
      if (is_error) then
        ifail = 1
        return
      endif


      write(*,*) 's14aef not implemented yet'

      s14aef = 0.0
      return
      end


