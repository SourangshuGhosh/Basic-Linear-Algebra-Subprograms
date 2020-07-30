      function dzasum ( n, zx, incx )

c*********************************************************************72
c
cc DZASUM takes the sum of the absolute values.
c
c  Discussion:
c
c    This routine uses double precision complex arithmetic.
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
c    Input, double complex X(*), the vector.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, double precision DZASUM, the sum of the absolute values.
c
      implicit none

      double precision dzasum
      integer i
      integer incx
      integer ix
      integer n
      double precision stemp
      double precision zabs1
      double complex zx(*)

      dzasum = 0.0d0
      stemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do i = 1,n
        stemp = stemp + zabs1(zx(ix))
        ix = ix + incx
      end do
      dzasum = stemp
      return
c
c        code for increment equal to 1
c
   20 do i = 1,n
        stemp = stemp + zabs1(zx(i))
      end do
      dzasum = stemp
      return
      end
      function dznrm2 ( n, x, incx )

c*********************************************************************72
c
cc DZNRM2 returns the euclidean norm of a vector.
c
c  Discussion:
c
c    This routine uses double precision complex arithmetic.
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
c    Input, double complex X(*), the vector.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, double precision DZNRM2, the norm of the vector.
c
      implicit none

      integer                           incx, n
      double complex x( * )


      double precision      one         , zero
      parameter           ( one = 1.0d+0, zero = 0.0d+0 )

      double precision dznrm2
      integer               ix
      double precision      norm, scale, ssq, temp

      intrinsic             abs, dimag, dble, sqrt

      if( n.lt.1 .or. incx.lt.1 )then
         norm  = zero
      else
         scale = zero
         ssq   = one

         do ix = 1, 1 + ( n - 1 )*incx, incx
            if( dble( x( ix ) ).ne.zero )then
               temp = abs( dble( x( ix ) ) )
               if( scale.lt.temp )then
                  ssq   = one   + ssq*( scale/temp )**2
                  scale = temp
               else
                  ssq   = ssq   +     ( temp/scale )**2
               end if
            end if
            if( dimag( x( ix ) ).ne.zero )then
               temp = abs( dimag( x( ix ) ) )
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

      dznrm2 = norm

      return
      end
      function izamax ( n, zx, incx )

c*********************************************************************72
c
cc IZAMAX finds the index of element having maximum absolute value.
c
c  Discussion:
c
c    This routine uses double precision complex arithmetic.
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
c    Input, double complex X(*), the vector.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, integer IZAMAX, the index of the element of maximum
c    absolute value.
c
      implicit none

      integer izamax
      double complex zx(*)
      double precision smax
      integer i,incx,ix,n
      double precision zabs1

      izamax = 0
      if( n.lt.1 .or. incx.le.0 )return
      izamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = zabs1(zx(1))
      ix = ix + incx
      do i = 2,n
        if ( smax .lt. zabs1(zx(ix)) ) then
          izamax = i
          smax = zabs1(zx(ix))
        end if
        ix = ix + incx
      end do
      return
c
c        code for increment equal to 1
c
   20 smax = zabs1(zx(1))
      do i = 2,n
        if ( smax .lt. zabs1(zx(i)) ) then
          izamax = i
          smax = zabs1(zx(i))
        end if
      end do

      return
      end
      subroutine zaxpy ( n, za, zx, incx, zy, incy )

c*********************************************************************72
c
cc ZAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses double precision complex arithmetic.
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
c    Input, double complex CA, the multiplier of CX.
c
c    Input, double complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Input/output, double complex CY(*), the second vector.
c    On output, CY(*) has been replaced by CY(*) + CA * CX(*).
c
c    Input, integer INCY, the increment between successive entries of CY.
c
      implicit none

      double complex zx(*),zy(*),za
      integer i,incx,incy,ix,iy,n
      double precision zabs1
      if(n.le.0)return
      if (zabs1(za) .eq. 0.0d0) return
      if (incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c        code for both increments equal to 1
c
   20 do i = 1,n
        zy(i) = zy(i) + za*zx(i)
      end do
      return
      end
      subroutine zcopy ( n, zx, incx, zy, incy )

c*********************************************************************72
c
cc ZCOPY copies a vector, x, to a vector, y.
c
c  Discussion:
c
c    This routine uses double precision complex arithmetic.
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
c    Input, double complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Output, double complex CY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries of CY.
c
      implicit none

      double complex zx(*),zy(*)
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        zy(iy) = zx(ix)
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c        code for both increments equal to 1
c
   20 do i = 1,n
        zy(i) = zx(i)
      end do

      return
      end
      function zdotc ( n, zx, incx, zy, incy )

c*********************************************************************72
c
cc ZDOTC forms the dot product of a vector.
c
c  Discussion:
c
c    This routine uses double precision complex arithmetic.
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
c    Input, double complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries in CX.
c
c    Input, double complex CY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries in CY.
c
c    Output, double complex ZDOTC, the conjugated dot product of 
c    the corresponding entries of CX and CY.
c
      implicit none

      double complex zx(*),zy(*),ztemp
      double complex zdotc
      integer i,incx,incy,ix,iy,n

      ztemp = (0.0d0,0.0d0)
      zdotc = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        ztemp = ztemp + dconjg(zx(ix))*zy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      zdotc = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do i = 1,n
        ztemp = ztemp + dconjg(zx(i))*zy(i)
      end do

      zdotc = ztemp

      return
      end
      function zdotu ( n, zx, incx, zy, incy )

c*********************************************************************72
c
cc ZDOTU forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses double precision complex arithmetic.
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
c    Input, double complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries in CX.
c
c    Input, double complex CY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries in CY.
c
c    Output, double complex ZDOTU, the unconjugated dot product of 
c    the corresponding entries of CX and CY.
c
      implicit none

      double complex zx(*),zy(*),ztemp
      double complex zdotu
      integer i,incx,incy,ix,iy,n

      ztemp = (0.0d0,0.0d0)
      zdotu = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        ztemp = ztemp + zx(ix)*zy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      zdotu = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do i = 1,n
        ztemp = ztemp + zx(i)*zy(i)
      end do

      zdotu = ztemp

      return
      end
      subroutine zdrot ( n, zx, incx, zy, incy, c, s )

c*********************************************************************72
c
cc ZDROT applies a plane rotation.
c
c  Discussion:
c
c    The cos and sin (c and s) are double precision and the vectors 
c    zx and zy are double complex.
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
c    Input/output, double complex CX(*), one of the vectors to be rotated.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Input/output, double complex CY(*), one of the vectors to be rotated.
c
c    Input, integer INCY, the increment between successive elements of CY.
c
c    Input, double precision C, S, parameters (presumably the cosine and 
c    sine of some angle) that define a plane rotation.
c
      implicit none

      double complex zx(1),zy(1),ztemp
      double precision c,s
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
        ztemp = c*zx(ix) + s*zy(iy)
        zy(iy) = c*zy(iy) - s*zx(ix)
        zx(ix) = ztemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c       code for both increments equal to 1
c
   20 do i = 1,n
        ztemp = c*zx(i) + s*zy(i)
        zy(i) = c*zy(i) - s*zx(i)
        zx(i) = ztemp
      end do

      return
      end
      subroutine zdscal ( n, da, zx, incx )

c*********************************************************************72
c
cc ZDSCAL scales a vector by a constant.
c
c  Discussion:
c
c    This routine uses double precision complex arithmetic.
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
c    Input, double precision SA, the multiplier.
c
c    Input/output, double complex CX(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of 
c    the vector CX.
c
      implicit none

      double precision da
      integer i
      integer incx
      integer ix
      integer n
      double complex zx(*)

      if ( n .le. 0 .or. incx .le. 0 ) then
        return
      end if

      if ( incx .eq. 1 ) then

        do i = 1, n
          zx(i) = dcmplx ( da, 0.0d0 ) * zx(i)
        end do

      else

        ix = 1
        do i = 1, n
          zx(ix) = dcmplx ( da, 0.0d0 ) * zx(ix)
          ix = ix + incx
        end do

      end if

      return
      end
      subroutine zrotg ( ca, cb, c, s )

c*********************************************************************72
c
cc ZROTG generates a Givens rotation.
c
c  Discussion:
c
c    This routine uses double precision complex arithmetic.
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
c    Input/output, double complex CA, on input, the value A.  On output,
c    the value R.
c
c    Input, double complex CB, the value B.
c
c    Output, double precision C, the cosine of the 
c    Givens rotation.
c
c    Output, double complex S, the sine of the 
c    Givens rotation.
c
      implicit none

      double complex ca,cb,s
      double precision c
      double precision norm,scale
      double complex alpha

      if (cdabs(ca) .eq. 0.0d0) then
        c = 0.0d0
        s = (1.0d0,0.0d0)
        ca = cb
      else
        scale = cdabs(ca) + cdabs(cb)
        norm = scale*dsqrt((cdabs(ca/dcmplx(scale,0.0d0)))**2 +
     &                     (cdabs(cb/dcmplx(scale,0.0d0)))**2)
        alpha = ca /cdabs(ca)
        c = cdabs(ca) / norm
        s = alpha * dconjg(cb) / norm
        ca = alpha * norm
      end if

      return
      end
      subroutine zscal ( n, za, zx, incx )

c*********************************************************************72
c
cc ZSCAL scales a vector by a constant.
c
c  Discussion:
c
c    This routine uses double precision complex arithmetic.
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
c    Input, double complex CA, the multiplier.
c
c    Input/output, double complex CX(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
      implicit none

      double complex za,zx(*)
      integer i,incx,ix,n

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
      end do

      return
c
c        code for increment equal to 1
c
   20 do i = 1,n
        zx(i) = za*zx(i)
      end do

      return
      end
      subroutine zswap ( n, zx, incx, zy, incy )

c*********************************************************************72
c
cc ZSWAP interchanges two vectors.
c
c  Discussion:
c
c    This routine uses double precision complex arithmetic.
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
c    Input/output, double complex CX(*), one of the vectors to swap.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Input/output, double complex CY(*), one of the vectors to swap.
c
c    Input, integer INCY, the increment between successive elements of CY.
c
      implicit none

      double complex zx(*),zy(*),ztemp
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
        ztemp = zx(ix)
        zx(ix) = zy(iy)
        zy(iy) = ztemp
        ix = ix + incx
        iy = iy + incy
      end do

      return
c
c       code for both increments equal to 1
c
   20 do i = 1,n
        ztemp = zx(i)
        zx(i) = zy(i)
        zy(i) = ztemp
      end do

      return
      end
