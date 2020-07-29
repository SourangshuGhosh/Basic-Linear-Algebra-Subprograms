      function c4_uniform_01 ( seed )

c*********************************************************************72
c
cc C4_UNIFORM_01 returns a unit double precision complex pseudorandom number.
c
c  Discussion:
c
c    The angle should be uniformly distributed between 0 and 2 * PI,
c    the square root of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 March 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, complex C4_UNIFORM_01, a pseudorandom complex value.
c
      implicit none

      real pi
      parameter ( pi = 3.141592653589793E+00 )

      complex c4_uniform_01
      real r
      integer k
      integer seed
      real theta

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      r = sqrt ( real ( dble ( seed ) * 4.656612875D-10 ) )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      theta = 2.0E+00 * pi
     &  * real ( dble ( seed ) * 4.656612875D-10 )

      c4_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ) )

      return
      end
      subroutine c4mat_print ( m, n, a, title )

c*********************************************************************72
c
cc C4MAT_PRINT prints a C4MAT.
c
c  Discussion:
c
c    A C4MAT is a matrix of C4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    in the matrix.
c
c    Input, complex A(M,N), the matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      complex a(m,n)
      character * ( * ) title

      call c4mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine c4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, 
     &  title )

c*********************************************************************72
c
cc C4MAT_PRINT_SOME prints some of a C4MAT.
c
c  Discussion:
c
c    A C4MAT is a matrix of C4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    in the matrix.
c
c    Input, complex A(M,N), the matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 4 )
      integer m
      integer n

      complex a(m,n)
      character * ( 20 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title
      complex zero

      zero = cmplx ( 0.0E+00, 0.0E+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      if ( m .le. 0 .or. n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (None)' 
        return
      end if
c
c  Print the columns of the matrix, in strips of INCX.
c
      do j2lo = jlo, min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i10,10x)' ) j
        end do

        write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) INCX entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( aimag ( a(i,j) ) .eq. 0.0E+00 ) then
              write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j) )
            else
              write ( ctemp(j2), '(2g10.3)' ) a(i,j)
            end if

          end do

          write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine c4mat_test ( n, a )

c*********************************************************************72
c
cc C4MAT_TEST sets up a test matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Output, complex A(N,N), the Fourier matrix.
c
      implicit none

      integer n

      complex a(n,n)
      real angle
      complex c4_i
      integer i
      integer j
      real pi
      parameter ( pi = 3.141592653589793E+00 )

      c4_i = cmplx ( 0.0E+00, 1.0E+00 )

      do i = 1, n
        do j = 1, n

          angle = 2.0E+00 * pi * real ( ( i - 1 ) * ( j - 1 ) ) 
     &      / real ( n )

          a(i,j) = exp ( c4_i * angle ) / sqrt ( real ( n ) )

        end do
      end do

      return
      end
      subroutine c4mat_test_inverse ( n, a )

c*********************************************************************72
c
cc C4MAT_TEST_INVERSE returns the inverse of the test matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Output, complex A(N,N), the matrix.
c
      implicit none

      integer n

      complex a(n,n)
      complex b(n,n)
      integer i
      integer j

      call c4mat_test ( n, b )

      do j = 1, n
        do i = 1, n
          a(i,j) = conjg ( b(j,i) )
        end do
      end do

      return
      end
      function c8_uniform_01 ( seed )

c*********************************************************************72
c
cc C8_UNIFORM_01 returns a unit double precision complex pseudorandom number.
c
c  Discussion:
c
c    The angle should be uniformly distributed between 0 and 2 * PI,
c    the square root of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 February 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double complex C8_UNIFORM_01, a pseudorandom complex value.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      double precision r
      integer k
      integer seed
      double precision theta
      double complex c8_uniform_01

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r = sqrt ( dble ( seed ) * 4.656612875D-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      theta = 2.0D+00 * pi * ( dble ( seed ) * 4.656612875D-10 )

      c8_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ) )

      return
      end
      subroutine c8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc C8MAT_PRINT prints a C8MAT.
c
c  Discussion:
c
c    A C8MAT is a matrix of C8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    in the matrix.
c
c    Input, double complex A(M,N), the matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double complex a(m,n)
      character * ( * ) title

      call c8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine c8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, 
     &  title )

c*********************************************************************72
c
cc C8MAT_PRINT_SOME prints some of a C8MAT.
c
c  Discussion:
c
c    A C8MAT is a matrix of C8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    in the matrix.
c
c    Input, double complex A(M,N), the matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 4 )
      integer m
      integer n

      double complex a(m,n)
      character * ( 20 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title
      double complex zero

      zero = dcmplx ( 0.0D+00, 0.0D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      if ( m .le. 0 .or. n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (None)' 
        return
      end if
c
c  Print the columns of the matrix, in strips of INCX.
c
      do j2lo = jlo, min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i10,10x)' ) j
        end do

        write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) INCX entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( dimag ( a(i,j) ) .eq. 0.0D+00 ) then
              write ( ctemp(j2), '(g10.3,10x)' ) dreal ( a(i,j) )
            else
              write ( ctemp(j2), '(2g10.3)' ) a(i,j)
            end if

          end do

          write ( *, '(i5,a,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine c8mat_test ( n, a )

c*********************************************************************72
c
cc C8MAT_TEST sets up a test matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Output, double complex A(N,N), the Fourier matrix.
c
      implicit none

      integer n

      double complex a(n,n)
      double precision angle
      double complex c8_i
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      c8_i = dcmplx ( 0.0D+00, 1.0D+00 )

      do i = 1, n
        do j = 1, n

          angle = 2.0D+00 * pi * dble ( ( i - 1 ) * ( j - 1 ) ) 
     &      / dble ( n )

          a(i,j) = exp ( c8_i * angle ) / sqrt ( dble ( n ) )

        end do
      end do

      return
      end
      subroutine c8mat_test_inverse ( n, a )

c*********************************************************************72
c
cc C8MAT_TEST_INVERSE returns the inverse of the test matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Output, double complex A(N,N), the matrix.
c
      implicit none

      integer n

      double complex a(n,n)
      double complex b(n,n)
      integer i
      integer j

      call c8mat_test ( n, b )

      do j = 1, n
        do i = 1, n
          a(i,j) = dconjg ( b(j,i) )
        end do
      end do

      return
      end
      function cabs1 ( z )

c*********************************************************************72
c
cc CABS1 returns the L1 norm of a single precision complex number.
c
c  Discussion:
c
c    The L1 norm of a complex number is the sum of the absolute values
c    of the real and imaginary components.
c
c    CABS1 ( Z ) = abs ( real ( Z ) ) + abs ( imaginary ( Z ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 May 2002
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
c    Input, complex Z, the number whose norm is desired.
c
c    Output, real CABS1, the L1 norm of Z.
c
      implicit none

      real cabs1
      complex z

      cabs1 = abs ( real ( z ) ) + abs ( aimag ( z ) )

      return
      end
      function cabs2 ( z )

c*********************************************************************72
c
cc CABS2 returns the L2 norm of a single precision complex number.
c
c  Discussion:
c
c    The L2 norm of a complex number is the square root of the sum of the 
c    squares of the real and imaginary components.
c
c    CABS2 ( Z ) = sqrt ( ( real ( Z ) )**2 + ( imaginary ( Z ) )**2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 April 2006
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
c    Input, complex Z, the number whose norm is desired.
c
c    Output, real CABS2, the L2 norm of Z.
c
      implicit none

      real cabs2
      complex z

      cabs2 = sqrt ( ( real ( z ) )**2 + ( aimag ( z ) )**2 )

      return
      end
      function cmach ( job )

c*********************************************************************72
c
cc CMACH computes single precision complex machine parameters.
c
c  Discussion:
c
c    Assume the computer has
c
c      B = base of arithmetic;
c      T = number of base B digits;
c      L = smallest possible exponent;
c      U = largest possible exponent;
c
c    then
c
c      EPS = B^(1-T)
c      TINY = 100.0 * B^(-L+T)
c      HUGE = 0.01 * B^(U-T)
c
c    If complex division is done by
c
c      1 / (X+i*Y) = (X-i*Y) / (X^2+Y^2)
c
c    then
c
c      TINY = sqrt ( TINY )
c      HUGE = sqrt ( HUGE )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2014
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
c    Input, integer JOB:
c    1, EPS is desired;
c    2, TINY is desired;
c    3, HUGE is desired.
c
c    Output, real CMACH, the requested value.
c
      implicit none

      real cmach
      real eps
      real huge
      integer job
      real s
      real tiny

      eps = 1.0E+00

   10 continue

      eps = eps / 2.0E+00
      s = 1.0E+00 + eps
      if ( 1.0E+00 .lt. s ) then
        go to 10
      end if
      eps = 2.0E+00 * eps
      cmach =eps
      if ( job .eq. 1 ) then
        return
      end if

      s = 1.0E+00

   20 continue

      tiny = s
      s = s/16.0E+00
      if (s*1.0E+00 .ne. 0.0E+00) go to 20
      tiny = (tiny/eps)*100.0E+00
      s = real((1.0E+00,0.0E+00)/cmplx(tiny,0.0E+00))
      if (s .ne. 1.0E+00/tiny) tiny = sqrt(tiny)
      huge = 1.0E+00/tiny

      if ( job .eq. 1 ) then
        cmach = eps
      else if ( job .eq. 2 ) then
        cmach = tiny
      else if ( job .eq. 3 ) then
        cmach = huge
      end if

      return
      end
      function csign1 ( z1, z2 )

c*********************************************************************72
c
cc CSIGN1 is a single precision complex transfer-of-sign function.
c
c  Discussion:
c
c    The L1 norm is used.
c
c  Modified:
c
c    14 May 2004
c
c  Parameters:
c
c    Input, complex Z1, Z2, the arguments.
c
c    Output, complex CSIGN1,  a complex value, with the magnitude of
c    Z1, and the argument of Z2.
c
      implicit none

      real cabs1
      complex csign1
      complex z1
      complex z2

      if ( cabs1 ( z2 ) == 0.0E+00 ) then
        csign1 = ( 0.0E+00, 0.0E+00 )
      else
        csign1 = cabs1 ( z1 ) * ( z2 / cabs1 ( z2 ) )
      end if

      return
      end
      function csign2 ( z1, z2 )

c*********************************************************************72
c
cc CSIGN2 is a single precision complex transfer-of-sign function.
c
c  Discussion:
c
c    The L2 norm is used.
c
c  Modified:
c
c    19 March 2006
c
c  Parameters:
c
c    Input, complex Z1, Z2, the arguments.
c
c    Output, complex CSIGN2,  a complex value, with the magnitude of
c    Z1, and the argument of Z2.
c
      implicit none

      real cabs2
      complex csign2
      complex z1
      complex z2

      if ( cabs2 ( z2 ) == 0.0E+00 ) then
        csign2 = ( 0.0E+00, 0.0E+00 )
      else
        csign2 = cabs2 ( z1 ) * ( z2 / cabs2 ( z2 ) )
      end if

      return
      end
      function dmach ( job )

c*********************************************************************72
c
cc DMACH computes machine parameters of floating point arithmetic.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    Assume the computer has
c
c      B = base of arithmetic;
c      T = number of base B digits;
c      L = smallest possible exponent;
c      U = largest possible exponent;
c
c    then
c
c      EPS = B^(1-T)
c      TINY = 100.0 * B^(-L+T)
c      HUGE = 0.01 * B^(U-T)
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
c    Input, integer JOB:
c    1, EPS is desired;
c    2, TINY is desired;
c    3, HUGE is desired.
c
c    Output, double precision DMACH, the requested value.
c
      implicit none

      double precision dmach
      integer job
      double precision eps,tiny,huge,s

      eps = 1.0d0

   10 continue

      eps = eps/2.0d0
      s = 1.0d0 + eps
      if (s .gt. 1.0d0) go to 10
      eps = 2.0d0*eps

      s = 1.0d0

   20 continue

      tiny = s
      s = s/16.0d0
      if (s*1.0 .ne. 0.0d0) go to 20
      tiny = (tiny/eps)*100.0
      huge = 1.0d0/tiny

      if (job .eq. 1) dmach = eps
      if (job .eq. 2) dmach = tiny
      if (job .eq. 3) dmach = huge

      return
      end
      function lsame ( ca, cb )

c*********************************************************************72
c
cc LSAME returns TRUE if CA is the same letter as CB regardless of case.
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
c    Input, character CA, CB, the character to compare.
c
c    Output, logical LSAME, is TRUE if the characters are equal,
c    disregarding case.
c
      implicit none

      character ca
      character cb
      intrinsic ichar
      integer inta
      integer intb
      logical lsame
      integer zcode
c
c  Test if the characters are equal.
c
      lsame = ca .eq. cb

      if ( lsame ) then
        return
      end if
c
c  Now test for equivalence if both characters are alphabetic.
c
      zcode = ichar ( 'Z' )
c
c  Use 'Z' rather than 'A' so that ASCII can be detected on Prime
c  machines, on which ICHAR returns a value with bit 8 set.
c  ICHAR('A') on Prime machines returns 193 which is the same as
c  ICHAR('A') on an EBCDIC machine.
c
      inta = ichar ( ca )
      intb = ichar ( cb )
c
c  ASCII is assumed.
c  ZCODE is the ASCII code of either lower or upper case 'Z'.
c
      if ( zcode .eq. 90 .or. zcode .eq. 122 ) then

        if ( inta .ge. 97 .and. inta .le. 122 ) then
          inta = inta - 32
        end if

        if ( intb .ge. 97 .and. intb .le. 122 ) then
          intb = intb - 32
        end if
c
c  EBCDIC is assumed.
c  ZCODE is the EBCDIC code of either lower or upper case 'Z'.
c
      else if ( zcode .eq. 233 .or. zcode .eq. 169 ) then

        if ( inta .ge. 129 .and. inta .le. 137 .or.
     &       inta .ge. 145 .and. inta .le. 153 .or.
     &       inta .ge. 162 .and. inta .le. 169 ) then
          inta = inta + 64
        end if

        if ( intb .ge. 129 .and. intb .le. 137 .or.
     &       intb .ge. 145 .and. intb .le. 153 .or.
     &       intb .ge. 162 .and. intb .le. 169 ) then
          intb = intb + 64
        end if
c
c  ASCII is assumed, on Prime machines.
c  ZCODE is the ASCII code plus 128 of either lower or upper case 'Z'.
c
      else if ( zcode .eq. 218 .or. zcode .eq. 250 ) then

        if ( inta .ge. 225 .and. inta .le. 250 ) then
          inta = inta - 32
        end if

        if ( intb .ge. 225 .and. intb .le. 250 ) then
          intb = intb - 32
        end if

      end if

      lsame = inta .eq. intb

      return
      end
      function r4_abs ( x )

c*********************************************************************72
c
cc R4_ABS returns the absolute value of an R4.
c
c  Discussion:
c
c    FORTRAN90 supplies the ABS function, which should be used instead
c    of this function!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the number whose absolute value is desired.
c
c    Output, real R4_ABS, the absolute value of X.
c
      implicit none

      real r4_abs
      real x

      if ( 0.0 .le. x ) then
        r4_abs = + x
      else
        r4_abs = - x
      end if

      return
      end
      function r4_sign ( x )

c*********************************************************************72
c
cc R4_SIGN returns the sign of an R4.
c
c  Discussion:
c
c    value = -1 if X < 0;
c    value =  0 if X => 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the number whose sign is desired.
c
c    Output, real R4_SIGN, the sign of X:
c
      implicit none

      real r4_sign
      real x

      if ( x .lt. 0.0E+00 ) then
        r4_sign = -1.0E+00
      else
        r4_sign = +1.0E+00
      end if

      return
      end
      function r4_uniform_01 ( seed )

c*********************************************************************72
c
cc R4_UNIFORM_01 returns a unit single precision pseudorandom number.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2^31 - 1 )
c      r4_uniform_01 = seed / ( 2^31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R4_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    P A Lewis, A S Goodman, J M Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R4_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer k
      integer seed
      real r4_uniform_01

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r4_uniform_01 = real ( dble ( seed ) * 4.656612875D-10 )

      return
      end
      function r4_uniform_ab ( a, b, seed )

c*********************************************************************72
c
cc R4_UNIFORM_AB returns a scaled pseudorandom R4.
c
c  Discussion:
c
c    The pseudorandom number should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 January 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input, real A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R4_UNIFORM_AB, a number strictly between A and B.
c
      implicit none

      real a
      real b
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      real r4_uniform_ab
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_UNIFORM_AB - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r4_uniform_ab = a + ( b - a )
     &  * real ( dble ( seed ) * 4.656612875E-10 )

      return
      end
      subroutine r4mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R4MAT_PRINT prints an R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of R4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, real A(M,N), the matrix.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer m
      integer n

      real a(m,n)
      character*(*) title

      call r4mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi,
     &  title )

c*********************************************************************72
c
cc R4MAT_PRINT_SOME prints some of an R4MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, real A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      real a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      if ( m .le. 0 .or. n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (None)'
        return
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(g14.6)' ) a(i,j)

          end do

          write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine r4mat_test ( trans, lda, m, n, a )

c*********************************************************************72
c
cc R4MAT_TEST sets up a test matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c  
c  Modified:
c
c    10 February 2014
c
c  Author:
c
c    John Burkardt.
c
c  Parameters:
c
c    Input, character * ( 1 ) TRANS, indicates whether matrix is to be 
c    transposed.
c    'N', no transpose.
c    'T', transpose the matrix.
c
c    Input, integer LDA, the leading dimension of the matrix.
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Output, real A(LDA,*), the matrix.
c    if TRANS is 'N', then the matrix is stored in LDA*N entries,
c    as an M x N matrix;
c    if TRANS is 'T', then the matrix is stored in LDA*M entries,
c    as an N x M matrix.
c
      implicit none

      integer lda

      real a(lda,*)
      integer i
      integer j
      integer m
      integer n
      character * ( 1 ) trans

      if ( trans .eq. 'N' ) then

        do j = 1, n
          do i = 1, m
            a(i,j) = real ( 10 * i + j )
          end do
        end do

      else

        do j = 1, n
          do i = 1, m
            a(j,i) = real ( 10 * i + j )
          end do
        end do

      end if

      return
      end
      subroutine r4vec_print ( n, a, title )

c*********************************************************************72
c
cc R4VEC_PRINT prints an R4VEC.
c
c  Discussion:
c
c    An R4VEC is a vector of R4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, real A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      real a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

      return
      end
      function r8_abs ( x )

c*********************************************************************72
c
cc R8_ABS returns the absolute value of an R8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 September 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the number whose absolute value is desired.
c
c    Output, double precision R8_ABS, the absolute value of X.
c
      implicit none

      double precision r8_abs
      double precision x

      if ( 0.0D+00 .le. x ) then
        r8_abs = + x
      else
        r8_abs = - x
      end if

      return
      end
      function r8_sign ( x )

c*********************************************************************72
c
cc R8_SIGN returns the sign of an R8.
c
c  Discussion:
c
c    value = -1 if X < 0;
c    value = +1 if X => 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the number whose sign is desired.
c
c    Output, double precision R8_SIGN, the sign of X.
c
      implicit none

      double precision r8_sign
      double precision x

      if ( x .lt. 0.0D+00 ) then
        r8_sign = -1.0D+00
      else
        r8_sign = +1.0D+00
      end if

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit double precision pseudorandom number.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2^31 - 1 )
c      r8_uniform_01 = seed / ( 2^31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      r8_uniform_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Philip Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
      integer k
      integer seed

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      function r8_uniform_ab ( a, b, seed )

c*********************************************************************72
c
cc R8_UNIFORM_AB returns a pseudorandom R8 scaled to [A,B].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_AB, a number strictly between A and B.
c
      implicit none

      double precision a
      double precision b
      integer k
      double precision r8_uniform_ab
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_AB - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop 1
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r8_uniform_ab = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_PRINT prints an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 May 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character ( len = * ) title

      call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi,
     &  title )

c*********************************************************************72
c
cc R8MAT_PRINT_SOME prints some of an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      double precision a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      if ( m .le. 0 .or. n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (None)'
        return
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(g14.6)' ) a(i,j)

          end do

          write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine r8mat_test ( trans, lda, m, n, a )

c*********************************************************************72
c
cc R8MAT_TEST sets up a test matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c  
c  Modified:
c
c    10 February 2014
c
c  Author:
c
c    John Burkardt.
c
c  Parameters:
c
c    Input, character * ( 1 ) TRANS, indicates whether matrix is to be 
c    transposed.
c    'N', no transpose.
c    'T', transpose the matrix.
c
c    Input, integer LDA, the leading dimension of the matrix.
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Output, double precision A(LDA,*), the matrix.
c    if TRANS is 'N', then the matrix is stored in LDA*N entries,
c    as an M x N matrix;
c    if TRANS is 'T', then the matrix is stored in LDA*M entries,
c    as an N x M matrix.
c
      implicit none

      integer lda

      double precision a(lda,*)
      integer i
      integer j
      integer m
      integer n
      character * ( 1 ) trans

      if ( trans .eq. 'N' ) then

        do j = 1, n
          do i = 1, m
            a(i,j) = dble ( 10 * i + j )
          end do
        end do

      else

        do j = 1, n
          do i = 1, m
            a(j,i) = dble ( 10 * i + j )
          end do
        end do

      end if

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

      return
      end
      function smach ( job )

c*********************************************************************72
c
cc SMACH computes machine parameters of floating point arithmetic.
c
c  Discussion:
c
c    Assume the computer has
c
c      B = base of arithmetic;
c      T = number of base B digits;
c      L = smallest possible exponent;
c      U = largest possible exponent;
c
c    then
c
c      EPS = B^(1-T)
c      TINY = 100.0 * B^(-L+T)
c      HUGE = 0.01 * B^(U-T)
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
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input, integer JOB:
c    1, EPS is desired;
c    2, TINY is desired;
c    3, HUGE is desired.
c
c    Output, real SMACH, the requested value.
c
      implicit none

      integer job

      real eps,tiny,huge,s
      real smach

      eps = 1.0

   10 continue

      eps = eps / 2.0
      s = 1.0 + eps
      if ( 1.0 .lt. s ) then
        go to 10
      end if
      eps = 2.0*eps

      s = 1.0
   20 continue
      tiny = s
      s = s/16.0
      if (s*100. .ne. 0.0) then
        go to 20
      end if
      tiny = (tiny/eps)*100.0
      huge = 1.0/tiny

      if (job .eq. 1) then
        smach = eps
      else if (job .eq. 2) then
        smach = tiny
      else if (job .eq. 3) then
        smach = huge
      end if

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
      subroutine triangle_upper_to_i4 ( i, j, k )

c*********************************************************************72
c
cc TRIANGLE_UPPER_TO_I4 converts an upper triangular coordinate to an integer.
c
c  Discussion:
c
c    Triangular coordinates are handy when storing a naturally triangular
c    array (such as the upper half of a matrix) in a linear array.
c
c    Thus, for example, we might consider storing 
c
c    (1,1) (1,2) (1,3) (1,4)
c          (2,2) (2,3) (2,4)
c                (3,3) (3,4)
c                      (4,4)
c
c    as the linear array
c
c    (1,1) (1,2) (2,2) (1,3) (2,3) (3,3) (1,4) (2,4) (3,4) (4,4)    
c
c    Here, the quantities in parenthesis represent the natural row and
c    column indices of a single number when stored in a rectangular array.
c
c    Thus, our goal is, given the row I and column J of the data,
c    to produce the value K which indicates its position in the linear
c    array.
c
c    The triangular numbers are the indices associated with the
c    diagonal elements of the original array, T(1,1), T(2,2), T(3,3)
c    and so on.
c
c    The formula is:
c
c      K = I + ( (J-1) * J ) / 2
c
c  First Values:
c
c     I  J  K
c
c     0  0  0
c     1  1  1
c     1  2  2
c     2  2  3
c     1  3  4
c     2  3  5
c     3  3  6
c     1  4  7
c     2  4  8
c     3  4  9
c     4  4 10
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2017
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, the row and column indices.  I and J must
c    be nonnegative, and I must not be greater than J.
c
c    Output, integer K, the linear index of the (I,J) element.
c
      implicit none

      integer i
      integer j
      integer k

      if ( i .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGLE_UPPER_TO_I4 - Fatal error!'
        write ( *, '(a)' ) '  I < 0.'
        write ( *, '(a,i8)' ) '  I = ', i
        stop 1
      else if ( j .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGLE_UPPER_TO_I4 - Fatal error!'
        write ( *, '(a)' ) '  J < 0.'
        write ( *, '(a,i8)' ) '  J = ', j
        stop 1
      else if ( j .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGLE_UPPER_TO_I4 - Fatal error!'
        write ( *, '(a)' ) '  J < I.'
        write ( *, '(a,i8)' ) '  I = ', i
        write ( *, '(a,i8)' ) '  J = ', j
        stop 1
      end if

      k = i + ( ( j - 1 ) * j ) / 2

      return
      end
      subroutine xerbla ( srname, info )

c*********************************************************************72
c
cc XERBLA is an error handler.
c
c  Discussion:
c
c    XERBLA is called if an input parameter has an invalid value.  
c    A message is printed and execution stops.
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
c    This FORTRAN77 version by John Burkardt.
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
c    Input, character * ( * ) SRNAME, the name of the routine 
c    which called XERBLA.
c
c    Input, integer INFO, the position of the invalid parameter in the 
c    parameter list of the calling routine.
c
      implicit none

      integer info
      character * ( * ) srname

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XERBLA - Fatal error!'
      write ( *, '(a)' ) '  On entry to routine ' // trim ( srname )
      write ( *, '(a)' ) '  an illegal value was detected for' 
      write ( *, '(a,i6)' ) '  parameter number ', info
      stop 1

      end
      function zabs1 ( z )

c*********************************************************************72
c
cc ZABS1 returns the L1 norm of a double complex number.
c
c  Discussion:
c
c    The L1 norm of a complex number is the sum of the absolute values
c    of the real and imaginary components.
c
c    ZABS1 ( Z ) = abs ( real ( Z ) ) + abs ( imaginary ( Z ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 March 2006
c
c  Author:
c
c    John Burkardt
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
c    Input, double complex Z, the number whose norm is desired.
c
c    Output, double precision ZABS1, the L1 norm of Z.
c
      implicit none

      double complex z
      double precision zabs1

      zabs1 = abs ( dble ( z ) ) + abs ( dimag ( z ) )

      return
      end
      function zabs2 ( z )

c*********************************************************************72
c
cc ZABS2 returns the L2 norm of a double complex number.
c
c  Discussion:
c
c    The L2 norm of a complex number is the square root of the sum of the 
c    squares of the real and imaginary components.
c
c    ZABS2 ( Z ) = sqrt ( ( real ( Z ) )^2 + ( imaginary ( Z ) )^2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 April 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, double complex Z, the number whose norm is desired.
c
c    Output, double precision ZABS2, the L2 norm of Z.
c
      implicit none

      double precision zabs2
      double complex z

      zabs2 = sqrt ( ( dble ( z ) )**2 + ( dimag ( z ) )**2 )

      return
      end
      function zmach ( job )

c*********************************************************************72
c
cc ZMACH computes double complex floating point arithmetic constants.
c
c  Discussion:
c
c    Assume the computer has
c
c      B = base of arithmetic;
c      T = number of base B digits;
c      L = smallest possible exponent;
c      U = largest possible exponent;
c
c    then
c
c      EPS = B^(1-T)
c      TINY = 100.0 * B^(-L+T)
c      HUGE = 0.01 * B^(U-T)
c
c    If complex division is done by
c
c      1 / (X+i*Y) = (X-i*Y) / (X^2+Y^2)
c
c    then
c
c      TINY = sqrt ( TINY )
c      HUGE = sqrt ( HUGE )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 July 2007
c
c  Author:
c
c    This FORTRAN77 version by John Burkardt
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
c    Input, integer JOB, indicates the desired quantity.
c    1: EPS
c    2: TINY
c    3: HUGE
c
c    Output, double precision ZMACH, the value of the desired quantity.
c
      implicit none

      double precision eps
      double precision huge
      integer job
      double precision s
      double precision tiny
      double precision zmach

      eps = 1.0d0

   10 continue

      eps = eps / 2.0d0
      s = 1.0d0 + eps
      if ( 1.0d0 .lt. s ) then
        go to 10
      end if
      eps = 2.0d0 * eps

      s = 1.0d0
   20 continue
      tiny = s
      s = s / 16.0d0
      if ( s * 1.0d0 .ne. 0.0d0 ) then
        go to 20
      end if
      tiny = tiny / eps
      s = ( 1.0d0, 0.0d0 ) / dcmplx ( tiny, 0.0d0 )
      if ( s .ne. 1.0d0 / tiny ) then
        tiny = dsqrt ( tiny )
      end if

      huge = 1.0d0 / tiny

      if ( job .eq. 1 ) then
        zmach = eps
      else if ( job .eq. 2 ) then
        zmach = tiny
      else if ( job .eq. 3 ) then
        zmach = huge
      end if

      return
      end
      function zsign1 ( z1, z2 )

c*********************************************************************72
c
cc ZSIGN1 is a double precision complex transfer-of-sign function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 May 2004
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
c    Input, double complex Z1, Z2, the arguments.
c
c    Output, double complex ZSIGN1,  a complex value, with the
c    magnitude of Z1, and the argument of Z2.
c
      implicit none

      double complex z1
      double complex z2
      double precision zabs1
      double complex  zsign1

      if ( zabs1 ( z2 ) == 0.0D+00 ) then
        zsign1 = ( 0.0D+00, 0.0D+00 )
      else
        zsign1 = zabs1 ( z1 ) * ( z2 / zabs1 ( z2 ) )
      end if

      return
      end
      function zsign2 ( z1, z2 )

c*********************************************************************72
c
cc ZSIGN2 is a double precision complex transfer-of-sign function.
c
c  Discussion:
c
c    The L2 norm is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 March 2006
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
c    Input, double complex Z1, Z2, the arguments.
c
c    Output, double complex ZSIGN2, a complex value, with the magnitude of
c    Z1, and the argument of Z2.
c
      implicit none

      double complex z1
      double complex z2
      double precision zabs2
      double complex zsign2

      if ( zabs2 ( z2 ) == 0.0D+00 ) then
        zsign2 = ( 0.0D+00, 0.0D+00 )
      else
        zsign2 = zabs2 ( z1 ) * ( z2 / zabs2 ( z2 ) )
      end if

      return
      end
