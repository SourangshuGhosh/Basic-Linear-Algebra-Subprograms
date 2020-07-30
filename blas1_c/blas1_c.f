      subroutine caxpy ( n, ca, cx, incx, cy, incy )

c*********************************************************************72
c
cc CAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input, integer N, the number of elements in CX and CY.
c
c    Input, complex CA, the multiplier of CX.
c
c    Input, complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Input/output, complex CY(*), the second vector.
c    On output, CY(*) has been replaced by CY(*) + CA * CX(*).
c
c    Input, integer INCY, the increment between successive entries of CY.
c
      implicit none

      complex ca
      complex cx(*)
      complex cy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer n

      if ( n .le. 0 ) then
        return
      end if

      if ( abs ( real ( ca ) ) + abs ( aimag ( ca ) ) .eq. 0.0 ) then
        return
      end if
c
c  Code for both increments equal to 1.
c
      if ( incx .eq. 1 .and. incy .eq. 1 ) then

        do i = 1, n
          cy(i) = cy(i) + ca * cx(i)
        end do
c
c  Code for unequal increments or equal increments
c  not equal to 1.
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
          cy(iy) = cy(iy) + ca * cx(ix)
          ix = ix + incx
          iy = iy + incy
        end do

      end if

      return
      end
      subroutine ccopy ( n, cx, incx, cy, incy )

c*********************************************************************72
c
cc CCOPY copies a vector, x, to a vector, y.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input, integer N, the number of elements in CX and CY.
c
c    Input, complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Output, complex CY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries of CY.
c
      implicit none

      complex cx(*),cy(*)
      integer i,incx,incy,ix,iy,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments
c  not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        cy(iy) = cx(ix)
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c  code for both increments equal to 1
c
   20 do i = 1,n
        cy(i) = cx(i)
      end do

      return
      end
      function cdotc ( n, cx, incx, cy, incy )

c*********************************************************************72
c
cc CDOTC forms the conjugated  dot product of two vectors.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input, complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries in CX.
c
c    Input, complex CY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries in CY.
c
c    Output, complex CDOTC, the conjugated dot product of 
c    the corresponding entries of CX and CY.
c
      implicit none

      complex cdotc
      complex cx(*),cy(*),ctemp
      integer i,incx,incy,ix,iy,n

      ctemp = (0.0,0.0)
      cdotc = (0.0,0.0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments
c  not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        ctemp = ctemp + conjg(cx(ix))*cy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      cdotc = ctemp
      return
c
c  code for both increments equal to 1
c
   20 do i = 1,n
        ctemp = ctemp + conjg(cx(i))*cy(i)
      end do

      cdotc = ctemp
      return
      end
      function cdotu ( n, cx, incx, cy, incy )

c*********************************************************************72
c
cc CDOTU forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input, complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries in CX.
c
c    Input, complex CY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries in CY.
c
c    Output, complex CDOTU, the unconjugated dot product of 
c    the corresponding entries of CX and CY.
c
      implicit none

      complex cdotu
      complex cx(*),cy(*),ctemp
      integer i,incx,incy,ix,iy,n

      ctemp = (0.0,0.0)
      cdotu = (0.0,0.0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments
c  not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        ctemp = ctemp + cx(ix)*cy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      cdotu = ctemp
      return
c
c  code for both increments equal to 1
c
   20 do i = 1,n
        ctemp = ctemp + cx(i)*cy(i)
      end do

      cdotu = ctemp

      return
      end
      subroutine crotg ( ca, cb, c, s )

c*********************************************************************72
c
cc CROTG generates a Givens rotation.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c    Given values A and B, this routine computes:
c
c    If A = 0:
c
c      R = B
c      C = 0
c      S = (1,0).
c
c    If A /= 0:
c
c      ALPHA = A / abs ( A )
c      NORM  = sqrt ( ( abs ( A ) )^2 + ( abs ( B ) )^2 )
c      R     = ALPHA * NORM
c      C     = abs ( A ) / NORM
c      S     = ALPHA * conj ( B ) / NORM
c
c    In either case, the computed numbers satisfy the equation:
c
c    (         C    S ) * ( A ) = ( R )
c    ( -conj ( S )  C )   ( B ) = ( 0 )
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
c    Input/output, complex CA; on input, the value A.  On output,
c    the value R.
c
c    Input, complex CB, the value B.
c
c    Output, real C, the cosine of the Givens rotation.
c
c    Output, complex S, the sine of the Givens rotation.
c
      implicit none

      complex ca,cb,s
      real c
      real norm,scale
      complex alpha

      if ( cabs(ca) .eq. 0.0 ) then
         c = 0.
         s = (1.,0.)
         ca = cb
      else
         scale = cabs(ca) + cabs(cb)
         norm = scale * sqrt((cabs(ca/scale))**2 + (cabs(cb/scale))**2)
         alpha = ca /cabs(ca)
         c = cabs(ca) / norm
         s = alpha * conjg(cb) / norm
         ca = alpha * norm
      end if

      return
      end
      subroutine cscal ( n, ca, cx, incx )

c*********************************************************************72
c
cc CSCAL scales a vector by a constant.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input, complex CA, the multiplier.
c
c    Input/output, complex CX(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
      implicit none

      complex ca,cx(*)
      integer i,incx,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        cx(i) = ca*cx(i)
      end do
      return
c
c  code for increment equal to 1
c
   20 do i = 1,n
        cx(i) = ca*cx(i)
      end do

      return
      end
      subroutine csrot ( n, cx, incx, cy, incy, c, s )

c*********************************************************************72
c
cc CSROT applies a plane rotation.
c
c  Discussion:
c
c    The cos and sin (c and s) are real and the vectors cx and cy are complex.
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
c    Input/output, complex CX(*), one of the vectors to be rotated.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Input/output, complex CY(*), one of the vectors to be rotated.
c
c    Input, integer INCY, the increment between successive elements of CY.
c
c    Input, real C, S, parameters (presumably the cosine and sine of
c    some angle) that define a plane rotation.
c
      implicit none

      complex cx(1),cy(1),ctemp
      real c,s
      integer i,incx,incy,ix,iy,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        ctemp = c*cx(ix) + s*cy(iy)
        cy(iy) = c*cy(iy) - s*cx(ix)
        cx(ix) = ctemp
        ix = ix + incx
        iy = iy + incy
      end do

      return
c
c       code for both increments equal to 1
c
   20 do i = 1,n
        ctemp = c*cx(i) + s*cy(i)
        cy(i) = c*cy(i) - s*cx(i)
        cx(i) = ctemp
      end do
      return
      end
      subroutine csscal ( n, sa, cx, incx )

c*********************************************************************72
c
cc CSSCAL scales a complex vector by a real constant.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input, real SA, the multiplier.
c
c    Input/output, complex CX(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of 
c    the vector CX.
c
      implicit none

      complex cx(*)
      real sa
      integer i,incx,n,nincx

      if( n.le.0 .or. incx.le.0 )return

      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        cx(i) = cmplx ( sa*real(cx(i)), sa*aimag(cx(i)) )
      end do
      return
c
c  code for increment equal to 1
c
   20 do i = 1,n
        cx(i) = cmplx ( sa*real(cx(i)), sa*aimag(cx(i)) )
      end do

      return
      end
      subroutine cswap ( n, cx, incx, cy, incy )

c*********************************************************************72
c
cc CSWAP interchanges two vectors.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input/output, complex CX(*), one of the vectors to swap.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Input/output, complex CY(*), one of the vectors to swap.
c
c    Input, integer INCY, the increment between successive elements of CY.
c
      implicit none

      complex cx(*),cy(*),ctemp
      integer i,incx,incy,ix,iy,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal
c  to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        ctemp = cx(ix)
        cx(ix) = cy(iy)
        cy(iy) = ctemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c       code for both increments equal to 1
c
   20 do i = 1,n
        ctemp = cx(i)
        cx(i) = cy(i)
        cy(i) = ctemp
      end do

      return
      end
      function icamax ( n, cx, incx )

c*********************************************************************72
c
cc ICAMAX finds the index of element having maximum absolute value.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input, complex X(*), the vector.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, integer ICAMAX, the index of the element of maximum
c    absolute value.
c
      implicit none

      complex cx(*)
      integer icamax
      real smax
      integer i,incx,ix,n
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))

      icamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      icamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = cabs1(cx(1))
      ix = ix + incx
      do i = 2,n
         if(cabs1(cx(ix)).le.smax) go to 5
         icamax = i
         smax = cabs1(cx(ix))
    5    ix = ix + incx
      end do

      return
c
c        code for increment equal to 1
c
   20 smax = cabs1(cx(1))
      do i = 2,n
        if ( smax .lt. cabs1(cx(i)) ) then
          icamax = i
          smax = cabs1(cx(i))
        end if
      end do

      return
      end
      function scasum ( n, cx, incx )

c*********************************************************************72
c
cc SCASUM takes the sum of the absolute values of a complex vector.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input, complex X(*), the vector.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, real SCASUM, the sum of the absolute values.
c
      implicit none

      complex cx(*)
      real scasum
      real stemp
      integer i,incx,n,nincx

      scasum = 0.0e0
      stemp = 0.0e0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        stemp = stemp + abs(real(cx(i))) + abs(aimag(cx(i)))
      end do
      scasum = stemp
      return
c
c        code for increment equal to 1
c
   20 do i = 1,n
        stemp = stemp + abs(real(cx(i))) + abs(aimag(cx(i)))
      end do

      scasum = stemp

      return
      end
      function scnrm2 ( n, x, incx )

c*********************************************************************72
c
cc SCNRM2 returns the euclidean norm of a complex vector.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c    SCNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
c            = sqrt ( dot_product ( x(1:n), x(1:n) ) )
c
c  Modified:
c
c    07 July 2007
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
c    Input, complex X(*), the vector.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, real SCNRM2, the norm of the vector.
c
      implicit none

      integer incx
      integer ix
      integer n
      real one
      parameter ( one = 1.0e+0 )
      real scnrm2
      complex x( * )
      real zero
      parameter ( zero = 0.0e+0 )


      real                  norm, scale, ssq, temp
      intrinsic             abs, aimag, real, sqrt

      if( n.lt.1 .or. incx.lt.1 )then
         norm  = zero
      else
         scale = zero
         ssq   = one
c        The following loop is equivalent to this call to the LAPACK
c        auxiliary routine:
c        CALL CLASSQ( N, X, INCX, SCALE, SSQ )
c
         do ix = 1, 1 + ( n - 1 )*incx, incx
            if( real( x( ix ) ).ne.zero )then
               temp = abs( real( x( ix ) ) )
               if( scale.lt.temp )then
                  ssq   = one   + ssq*( scale/temp )**2
                  scale = temp
               else
                  ssq   = ssq   +     ( temp/scale )**2
               end if
            end if
            if( aimag( x( ix ) ).ne.zero )then
               temp = abs( aimag( x( ix ) ) )
               if( scale.lt.temp )then
                  ssq   = one   + ssq*( scale/temp )**2
                  scale = temp
               else
                  ssq   = ssq   +     ( temp/scale )**2
               end if
            end if
         end do
         norm  = scale * sqrt( ssq )
      end if

      scnrm2 = norm

      return
      end
