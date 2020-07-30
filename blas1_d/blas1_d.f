      function dasum ( n, dx, incx )

c*********************************************************************72
c
cc DASUM takes the sum of the absolute values.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 2008
c
c  Author:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision X(*), the vector to be examined.
c
c    Input, integer INCX, the increment between successive entries of X.
c    INCX must not be negative.
c
c    Output, double precision DASUM, the sum of the absolute values of X.
c
      implicit none

      double precision dasum
      double precision dtemp
      double precision dx(*)
      integer i
      integer incx
      integer m
      integer n
      integer nincx

      dasum = 0.0D+00
      dtemp = 0.0D+00

      if ( n .le. 0 ) then
        return
      end if

      if ( incx .le. 0 ) then
        return
      end if

      if ( incx .ne. 1 ) then

        nincx = n * incx
        do i = 1, nincx, incx
          dtemp = dtemp + dabs ( dx(i) )
        end do

      else

        m = mod ( n, 6 )

        do i = 1, m
          dtemp = dtemp + dabs ( dx(i) )
        end do

        do i = m + 1, n, 6
          dtemp = dtemp 
     &      + dabs ( dx(i) ) 
     &      + dabs ( dx(i+1) ) 
     &      + dabs ( dx(i+2) )
     &      + dabs ( dx(i+3) ) 
     &      + dabs ( dx(i+4) ) 
     &      + dabs ( dx(i+5) )
        end do

      end if

      dasum = dtemp

      return
      end
      subroutine daxpy ( n, da, dx, incx, dy, incy )

c*********************************************************************72
c
cc DAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 2008
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of elements in DX and DY.
c
c    Input, double precision DA, the multiplier of DX.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of DX.
c
c    Input/output, double precision DY(*), the second vector.
c    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
c
c    Input, integer INCY, the increment between successive entries of DY.
c
      implicit none

      double precision da
      double precision dx(*)
      double precision dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n

      if ( n .le. 0 ) then
        return
      end if

      if ( da .eq. 0.0d0 ) then
        return
      end if

      if ( incx .ne. 1 .or. incy .ne. 1 ) then

        ix = 1
        iy = 1
        if ( incx .lt. 0 ) then
          ix = (-n+1) * incx + 1
        end if

        if ( incy .lt. 0 ) then
          iy = (-n+1) * incy + 1
        end if

        do i = 1, n
          dy(iy) = dy(iy) + da * dx(ix)
          ix = ix + incx
          iy = iy + incy
        end do

      else

        m = mod ( n, 4 )

        do i = 1, m
          dy(i) = dy(i) + da * dx(i)
        end do

        do i = m + 1, n, 4
          dy(i) = dy(i) + da * dx(i)
          dy(i + 1) = dy(i + 1) + da * dx(i + 1)
          dy(i + 2) = dy(i + 2) + da * dx(i + 2)
          dy(i + 3) = dy(i + 3) + da * dx(i + 3)
        end do

      end if

      return
      end
      subroutine dcopy ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DCOPY copies a vector.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    The routine uses unrolled loops for increments equal to one.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of elements in DX and DY.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of DX.
c
c    Output, double precision DY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries of DY.
c
      implicit none

      double precision dx(*)
      double precision dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n

      if ( n .le. 0 ) then
        return
      end if

      if ( incx .ne. 1 .or. incy .ne. 1 ) then

        ix = 1
        iy = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        if(incy.lt.0)iy = (-n+1)*incy + 1
        do i = 1, n
          dy(iy) = dx(ix)
          ix = ix + incx
          iy = iy + incy
        end do

      else

        m = mod(n,7)

        do i = 1,m
          dy(i) = dx(i)
        end do

        do i = m + 1, n, 7
          dy(i) = dx(i)
          dy(i + 1) = dx(i + 1)
          dy(i + 2) = dx(i + 2)
          dy(i + 3) = dx(i + 3)
          dy(i + 4) = dx(i + 4)
          dy(i + 5) = dx(i + 5)
          dy(i + 6) = dx(i + 6)
        end do

      end if

      return
      end
      function ddot ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DDOT forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c    
c  Modified:
c
c    09 February 2014
c
c  Author:
c
c    Original FORTRAN77 version by Jack Dongarra.
c    This version by John Burkardt
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries in DX.
c
c    Input, double precision DY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries in DY.
c
c    Output, double precision DDOT, the sum of the product of the 
c    corresponding entries of DX and DY.
c
      implicit none

      double precision ddot
      double precision dx(*)
      double precision dy(*)
      double precision dtemp
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n

      ddot = 0.0D+00
      dtemp = 0.0D+00

      if ( n .le. 0 ) then
        return
      end if

      if ( incx .eq. 1 .and. incy .eq. 1 ) then

        m = mod ( n, 5 )

        do i = 1, m
          dtemp = dtemp + dx(i) * dy(i)
        end do

        do i = m + 1, n, 5
          dtemp = dtemp 
     &      + dx(i)   * dy(i) 
     &      + dx(i+1) * dy(i+1) 
     &      + dx(i+2) * dy(i+2) 
     &      + dx(i+3) * dy(i+3) 
     &      + dx(i+4) * dy(i+4)
        end do
c
c  Code for unequal increments or equal increments not equal to 1
c
      else

        if ( incx .lt. 0 ) then
          ix = ( - n + 1 ) * incx + 1
        else
          ix = 1
        end if

        if ( incy .lt. 0 ) then
          iy = ( - n + 1 ) * incy + 1
        else
          iy = 1
        end if

        do i = 1, n
          dtemp = dtemp + dx(ix) * dy(iy)
          ix = ix + incx
          iy = iy + incy
        end do

      end if

      ddot = dtemp

      return
      end
      function dnrm2 ( n, x, incx )

c*********************************************************************72
c
cc DNRM2 returns the euclidean norm of a vector. 
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c      DNRM2 ( X ) = sqrt ( X' * X )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2014
c
c  Author:
c
c    Sven Hammarling
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision X(*), the vector whose norm is to be computed.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, double precision DNRM2, the Euclidean norm of X.
c
      implicit none

      integer incx, n

      double precision dnrm2 
      double precision x( * )

      double precision one         , zero
      parameter ( one = 1.0d+0, zero = 0.0d+0 )

      integer ix
      double precision absxi, norm, scale, ssq

      intrinsic abs, sqrt

      if( n.lt.1 .or. incx.lt.1 )then
         norm  = zero
      else if( n.eq.1 )then
         norm  = abs( x( 1 ) )
      else
         scale = zero
         ssq   = one
c
c  The following loop is equivalent to this call to the LAPACK
c  auxiliary routine:
c  call dlassq( n, x, incx, scale, ssq )
c
         do ix = 1, 1 + ( n - 1 )*incx, incx
            if( x( ix ).ne.zero )then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi )then
                  ssq   = one   + ssq*( scale/absxi )**2
                  scale = absxi
               else
                  ssq   = ssq   +     ( absxi/scale )**2
               end if
            end if
         end do
         norm  = scale * sqrt( ssq )
      end if

      dnrm2 = norm

      return
      end
      subroutine drot ( n, dx, incx, dy, incy, c, s )

c*********************************************************************72
c
cc DROT applies a plane rotation.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input/output, double precision X(*), one of the vectors to be rotated.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, double precision Y(*), one of the vectors to be rotated.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
c    Input, double precision C, S, parameters (presumably the cosine and
c    sine of some angle) that define a plane rotation.
c
      implicit none

      double precision c
      double precision dtemp
      double precision dx(*)
      double precision dy(*)
      double precision s
      integer i,incx,incy,ix,iy,n

      if ( n .le. 0 ) then
        return
      end if

      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0) then
        ix = (-n+1)*incx + 1
      end if

      if(incy.lt.0) then
        iy = (-n+1)*incy + 1
      end if

      do i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c  code for both increments equal to 1
c
   20 do i = 1,n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
      end do

      return
      end
      subroutine drotg ( da, db, c, s )

c*********************************************************************72
c
cc DROTG constructs a Givens plane rotation.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    Given values A and B, this routine computes
c
c    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
c          = sign ( B ) if abs ( A ) <= abs ( B );
c
c    R     = SIGMA * ( A * A + B * B );
c
c    C = A / R if R is not 0
c      = 1     if R is 0;
c
c    S = B / R if R is not 0,
c        0     if R is 0.
c
c    The computed numbers then satisfy the equation
c
c    (  C  S ) ( A ) = ( R )
c    ( -S  C ) ( B ) = ( 0 )
c
c    The routine also computes
c
c    Z = S     if abs ( A ) > abs ( B ),
c      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
c      = 1     if C is 0.
c
c    The single value Z encodes C and S, and hence the rotation:
c
c    If Z = 1, set C = 0 and S = 1;
c    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
c    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input/output, double precision SA, SB.  On input, SA and SB are the values
c    A and B.  On output, SA is overwritten with R, and SB is
c    overwritten with Z.
c
c    Output, double precision C, S, the cosine and sine of the
c    Givens rotation.
c
      implicit none

      double precision c
      double precision da
      double precision db,s,roe,scale,r,z

      roe = db
      if( dabs(da) .gt. dabs(db) ) roe = da
      scale = dabs(da) + dabs(db)

      if( scale .eq. 0.0d0 ) then
         c = 1.0d0
         s = 0.0d0
         r = 0.0d0
         z = 0.0d0
      else
        r = scale*dsqrt((da/scale)**2 + (db/scale)**2)
        r = dsign(1.0d0,roe)*r
        c = da/r
        s = db/r
        z = 1.0d0
        if( dabs(da) .gt. dabs(db) ) z = s
        if( dabs(db) .ge. dabs(da) .and. c .ne. 0.0d0 ) z = 1.0d0/c
      end if

      da = r
      db = z

      return
      end
      subroutine drotm ( n, dx, incx, dy, incy, dparam )

c*********************************************************************72
c
cc DROTM applies a modified Givens rotation matrix.
c
c  Purpose
c  =======
c
c     APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
c
c     (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
c     (DY**T)
c
c     DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
c     LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
c     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
c
c     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
c
c       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
c     H=(          )    (          )    (          )    (          )
c       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
c     SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2014
c
c  Author:
c
c    Original FORTRAN77 version by Jack Dongarra.
c    This version by John Burkardt
c
c  Arguments
c  =========
c
c  N      (input) INTEGER
c         number of elements in input vector(s)
c
c  DX     (input/output) DOUBLE PRECISION array, dimension N
c         double precision vector with N elements
c
c  INCX   (input) INTEGER
c         storage spacing between elements of DX
c
c  DY     (input/output) DOUBLE PRECISION array, dimension N
c         double precision vector with N elements
c
c  INCY   (input) INTEGER
c         storage spacing between elements of DY
c
c  DPARAM (input/output)  DOUBLE PRECISION array, dimension 5 
c     DPARAM(1)=DFLAG
c     DPARAM(2)=DH11
c     DPARAM(3)=DH21
c     DPARAM(4)=DH12
c     DPARAM(5)=DH22
c
      integer incx,incy,n

      double precision dparam(5),dx(*),dy(*)
      double precision dflag,dh11,dh12,dh21,dh22,two,w,z,zero
      integer i,kx,ky,nsteps

      data zero,two/0.d0,2.d0/

      dflag = dparam(1)
      if (n.le.0 .or. (dflag+two.eq.zero)) return
      if (incx.eq.incy.and.incx.gt.0) then

         nsteps = n*incx
         if (dflag.lt.zero) then
            dh11 = dparam(2)
            dh12 = dparam(4)
            dh21 = dparam(3)
            dh22 = dparam(5)
            do i = 1,nsteps,incx
               w = dx(i)
               z = dy(i)
               dx(i) = w*dh11 + z*dh12
               dy(i) = w*dh21 + z*dh22
            end do
         else if (dflag.eq.zero) then
            dh12 = dparam(4)
            dh21 = dparam(3)
            do i = 1,nsteps,incx
               w = dx(i)
               z = dy(i)
               dx(i) = w + z*dh12
               dy(i) = w*dh21 + z
            end do
         else
            dh11 = dparam(2)
            dh22 = dparam(5)
            do i = 1,nsteps,incx
               w = dx(i)
               z = dy(i)
               dx(i) = w*dh11 + z
               dy(i) = -w + dh22*z
            end do
         end if
      else
         kx = 1
         ky = 1
         if (incx.lt.0) kx = 1 + (1-n)*incx
         if (incy.lt.0) ky = 1 + (1-n)*incy
*
         if (dflag.lt.zero) then
            dh11 = dparam(2)
            dh12 = dparam(4)
            dh21 = dparam(3)
            dh22 = dparam(5)
            do i = 1,n
               w = dx(kx)
               z = dy(ky)
               dx(kx) = w*dh11 + z*dh12
               dy(ky) = w*dh21 + z*dh22
               kx = kx + incx
               ky = ky + incy
            end do
         else if (dflag.eq.zero) then
            dh12 = dparam(4)
            dh21 = dparam(3)
            do i = 1,n
               w = dx(kx)
               z = dy(ky)
               dx(kx) = w + z*dh12
               dy(ky) = w*dh21 + z
               kx = kx + incx
               ky = ky + incy
            end do
         else
             dh11 = dparam(2)
             dh22 = dparam(5)
             do i = 1,n
                w = dx(kx)
                z = dy(ky)
                dx(kx) = w*dh11 + z
                dy(ky) = -w + dh22*z
                kx = kx + incx
                ky = ky + incy
            end do
         end if
      end if
      return
      end
      subroutine drotmg(dd1,dd2,dx1,dy1,dparam)

c*********************************************************************72
c
cc DROTMG generates a modified Givens rotation matrix.
c
*  Purpose
*  =======
*
*     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
*     THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*
*     DY2)**T.
*     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
*
*     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
*
*       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
*     H=(          )    (          )    (          )    (          )
*       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
*     LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
*     RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
*     VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
*
*     THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
*     INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
*     OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
*
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2014
c
*  Arguments
*  =========
*
*  DD1    (input/output) DOUBLE PRECISION
*
*  DD2    (input/output) DOUBLE PRECISION
*
*  DX1    (input/output) DOUBLE PRECISION
*
*  DY1    (input) DOUBLE PRECISION
*
*  DPARAM (input/output)  DOUBLE PRECISION array, dimension 5
*     DPARAM(1)=DFLAG
*     DPARAM(2)=DH11
*     DPARAM(3)=DH21
*     DPARAM(4)=DH12
*     DPARAM(5)=DH22
*
      double precision dd1,dd2,dx1,dy1
      double precision dparam(5)
      double precision dflag,dh11,dh12,dh21,dh22,dp1,dp2,dq1,dq2,dtemp,
     $                 du,gam,gamsq,one,rgamsq,two,zero
*     ..
*     .. intrinsic functions ..
      intrinsic dabs
*     ..
*     .. data statements ..
*
      data zero,one,two/0.d0,1.d0,2.d0/
      data gam,gamsq,rgamsq/4096.d0,16777216.d0,5.9604645d-8/
*     ..

      if (dd1.lt.zero) then
*        go zero-h-d-and-dx1..
         dflag = -one
         dh11 = zero
         dh12 = zero
         dh21 = zero
         dh22 = zero
         dd1 = zero
         dd2 = zero
         dx1 = zero
      else
*        case-dd1-nonnegative
         dp2 = dd2*dy1
         if (dp2.eq.zero) then
            dflag = -two
            dparam(1) = dflag
            return
         end if 
*        regular-case..
         dp1 = dd1*dx1
         dq2 = dp2*dy1
         dq1 = dp1*dx1
*
         if (dabs(dq1).gt.dabs(dq2)) then
            dh21 = -dy1/dx1
            dh12 = dp2/dp1
*
            du = one - dh12*dh21
*
           if (du.gt.zero) then
             dflag = zero
             dd1 = dd1/du
             dd2 = dd2/du
             dx1 = dx1*du
           end if
         else

            if (dq2.lt.zero) then
*              go zero-h-d-and-dx1..
               dflag = -one
               dh11 = zero
               dh12 = zero
               dh21 = zero
               dh22 = zero
*
               dd1 = zero
               dd2 = zero
               dx1 = zero
            else
               dflag = one
               dh11 = dp1/dp2
               dh22 = dx1/dy1
               du = one + dh11*dh22
               dtemp = dd2/du
               dd2 = dd1/du
               dd1 = dtemp
               dx1 = dy1*du
            end if
         end if

*     procedure..scale-check
         if (dd1.ne.zero) then
            do while ((dd1.le.rgamsq) .or. (dd1.ge.gamsq))
               if (dflag.eq.zero) then
                  dh11 = one
                  dh22 = one
                  dflag = -one
               else
                  dh21 = -one
                  dh12 = one
                  dflag = -one
               end if
               if (dd1.le.rgamsq) then
                  dd1 = dd1*gam**2
                  dx1 = dx1/gam
                  dh11 = dh11/gam
                  dh12 = dh12/gam
               else
                  dd1 = dd1/gam**2
                  dx1 = dx1*gam
                  dh11 = dh11*gam
                  dh12 = dh12*gam
               end if
            enddo
         end if
  
         if (dd2.ne.zero) then
            do while ( (dabs(dd2).le.rgamsq) .or. (dabs(dd2).ge.gamsq) )
               if (dflag.eq.zero) then
                  dh11 = one
                  dh22 = one
                  dflag = -one
               else
                  dh21 = -one
                  dh12 = one
                  dflag = -one
               end if
               if (dabs(dd2).le.rgamsq) then
                  dd2 = dd2*gam**2
                  dh21 = dh21/gam
                  dh22 = dh22/gam
               else
                  dd2 = dd2/gam**2
                  dh21 = dh21*gam
                  dh22 = dh22*gam
               end if      
            end do
         end if
     
      end if

      if (dflag.lt.zero) then
         dparam(2) = dh11
         dparam(3) = dh21
         dparam(4) = dh12
         dparam(5) = dh22
      else if (dflag.eq.zero) then
         dparam(3) = dh21
         dparam(4) = dh12 
      else
         dparam(2) = dh11
         dparam(5) = dh22
      end if

      dparam(1) = dflag

      return
      end
      subroutine dscal ( n, da, dx, incx )

c*********************************************************************72
c
cc DSCAL scales a vector by a constant.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c    
c  Modified:
c
c    09 February 2014
c
c  Author:
c
c    Original FORTRAN77 version by Jack Dongarra.
c    This version by John Burkardt.
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision SA, the multiplier.
c
c    Input/output, double precision X(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of X.
c
      implicit none

      double precision da
      double precision dx(*)
      integer i
      integer incx
      integer m
      integer n
      integer nincx

      if ( n .le. 0 ) then
        return
      end if

      if ( incx .le. 0 ) then
        return
      end if

      if ( incx .eq. 1 ) then

        m = mod ( n, 5 )

        do i = 1, m
          dx(i) = da * dx(i)
        end do

        do i = m + 1, n, 5
          dx(i) =     da * dx(i)
          dx(i+1) = da * dx(i+1)
          dx(i+2) = da * dx(i+2)
          dx(i+3) = da * dx(i+3)
          dx(i+4) = da * dx(i+4)
        end do

      else

        nincx = n * incx
        do i = 1, nincx, incx
          dx(i) = da * dx(i)
        end do
      end if

      return
      end
      function dsdot ( n, sx, incx, sy, incy )

c*********************************************************************72
c
cc DSDOT computes the inner product of two vectors with extended precision.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
*     ..
*
*  AUTHORS
*  =======
*  Lawson, C. L., (JPL), Hanson, R. J., (SNLA), 
*  Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)
*
*  Purpose
*  =======
*  Compute the inner product of two vectors with extended
*  precision accumulation and result.
*
*  Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
*  DSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
*  where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
*  defined in a similar way using INCY.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         number of elements in input vector(s)
*
*  SX     (input) REAL array, dimension(N)
*         single precision vector with N elements
*
*  INCX   (input) INTEGER
*          storage spacing between elements of SX
*
*  SY     (input) REAL array, dimension(N)
*         single precision vector with N elements
*
*  INCY   (input) INTEGER
*         storage spacing between elements of SY
*
*  DSDOT  (output) DOUBLE PRECISION
*         DSDOT  double precision dot product (zero if N.LE.0)
*
*  Further Details
*  ===============
*
*  REFERENCES
*      
*  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
*  Krogh, Basic linear algebra subprograms for Fortran
*  usage, Algorithm No. 539, Transactions on Mathematical
*  Software 5, 3 (September 1979), pp. 308-323.
c
      double precision dsdot
      INTEGER INCX,INCY,N
      REAL SX(*),SY(*)
      integer i,kx,ky,ns

      intrinsic dble

      dsdot = 0.0d0
      if (n.le.0) return
      if (incx.eq.incy .and. incx.gt.0) then
c
c     code for equal, positive, non-unit increments.
c
         ns = n*incx
         do i = 1,ns,incx
            dsdot = dsdot + dble(sx(i))*dble(sy(i))
         end do
      else
c
c     code for unequal or nonpositive increments.
c
         kx = 1
         ky = 1
         if (incx.lt.0) kx = 1 + (1-n)*incx
         if (incy.lt.0) ky = 1 + (1-n)*incy
         do i = 1,n
            dsdot = dsdot + dble(sx(kx))*dble(sy(ky))
            kx = kx + incx
            ky = ky + incy
         end do
      end if

      return
      end
      subroutine dswap ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DSWAP interchanges two vectors.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input/output, double precision X(*), one of the vectors to swap.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, double precision Y(*), one of the vectors to swap.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
      implicit none

      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
      end do
      if( n .lt. 3 ) return
   40 continue

      do i = m+1, n, 3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
      end do

      return
      end
      function idamax ( n, dx, incx )

c*********************************************************************72
c
cc IDAMAX finds the index of element having maximum absolute value.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision X(*), the vector to be examined.
c
c    Input, integer INCX, the increment between successive entries of SX.
c
c    Output, integer IDAMAX, the index of the element of SX of maximum
c    absolute value.
c
      implicit none

      double precision dmax
      double precision dx(*)
      integer i
      integer idamax
      integer incx
      integer ix
      integer n

      idamax = 0
      if ( n .lt. 1 .or. incx .le. 0 ) then
        return
      end if

      idamax = 1
      if ( n .eq. 1 ) then
        return
      end if

      if ( incx .eq. 1 ) then

        dmax = dabs ( dx(1) )
        do i = 2, n
          if( dmax .lt. dabs ( dx(i) ) ) then
            idamax = i
            dmax = dabs ( dx(i) )
          end if
        end do

      else

        ix = 1
        dmax = dabs ( dx(1) )
        ix = ix + incx
        do  i = 2, n
          if ( dmax .lt. dabs ( dx(ix) ) ) then
            idamax = i
            dmax = dabs ( dx(ix) )
          end if
          ix = ix + incx
        end do

      end if

      return
      end
