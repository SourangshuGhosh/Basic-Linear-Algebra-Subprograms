      program main

c*********************************************************************72
c
cc MAIN is the main program for BLAS1_D_PRB.
c
c  Discussion:
c
c    BLAS1_D_PRB tests the BLAS1 library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLAS1_D_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BLAS1_D library.'

      call dasum_test ( )
      call daxpy_test ( )
      call dcopy_test ( )
      call ddot_test ( )
      call dnrm2_test ( )
      call drot_test ( )
      call drotg_test ( )
      call dscal_test ( )
      call dswap_test ( )
      call idamax_test ( )
!
!  Terminate.
!
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLAS1_D_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine dasum_test ( )

c*********************************************************************72
c
cc DASUM_TEST tests DASUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer lda
      integer ma
      integer na
      integer nx

      parameter ( lda = 6 )
      parameter ( ma = 5 )
      parameter ( na = 4 )
      parameter ( nx = 10 )

      double precision a(lda,na)
      double precision dasum
      integer i
      integer j
      double precision x(nx)

      do i = 1, nx
        x(i) = (-1.0D+00)**i * dble ( 2 * i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DASUM_TEST'
      write ( *, '(a)' ) '  DASUM adds the absolute values of '
      write ( *, '(a)' ) '  elements of a double precision vector.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, nx
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  DASUM ( NX, X, 1 ) =   ',
     &  dasum ( nx, x, 1 )
      write ( *, '(a,g14.6)' ) '  DASUM ( NX/2, X, 2 ) = ',
     &  dasum ( nx/2, x, 2 )
      write ( *, '(a,g14.6)' ) '  DASUM ( 2, X, NX/2 ) = ',
     &  dasum ( 2, x, nx/2 )

      do i = 1, lda
        do j = 1, na
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, ma
        do j = 1, na
          a(i,j) = (-1.0D+00)**(i+j) * dble ( 10 * i + j )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Demonstrate with a matrix A:'
      write ( *, '(a)' ) ' '
      do i = 1, ma
        write ( *, '(2x,5g14.6)' ) ( a(i,j), j = 1, na )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  DASUM(MA,A(1,2),1) =   ',
     &  dasum ( ma, a(1,2), 1 )
      write ( *, '(a,g14.6)' ) '  DASUM(NA,A(2,1),LDA) = ',
     &  dasum ( na, a(2,1), lda )

      return
      end
      subroutine daxpy_test ( )

c*********************************************************************72
c
cc DAXPY_TEST tests DAXPY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      double precision da
      integer i
      double precision x(n)
      double precision y(n)

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = dble ( 100 * i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DAXPY_TEST'
      write ( *, '(a)' ) '  DAXPY adds a double precision multiple of '
      write ( *, '(a)' ) '  vector X to vector Y.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Y = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      da = 1.0D+00
      call daxpy ( n, da, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DAXPY ( N, ', da, ', X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      do i = 1, n
        y(i) = dble ( 100 * i )
      end do

      da = -2.0D+00
      call daxpy ( n, da, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DAXPY ( N, ', da, ', X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      do i = 1, n
        y(i) = dble ( 100 * i )
      end do

      da = +3.0D+00
      call daxpy ( 3, da, x, 2, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DAXPY ( 3, ', da, ', X, 2, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      do i = 1, n
        y(i) = dble ( 100 * i )
      end do

      da = -4.0D+00
      call daxpy ( 3, da, x, 1, y, 2 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DAXPY ( 3, ', da, ', X, 1, Y, 2 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      return
      end
      subroutine dcopy_test ( )

c*********************************************************************72
c
cc DCOPY_TEST tests DCOPY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a(5,5)
      integer i
      integer j
      double precision x(10)
      double precision y(10)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DCOPY_TEST'
      write ( *, '(a)' ) '  DCOPY copies one double precision vector'
      write ( *, '(a)' ) '  into another.'

      do i = 1, 10
        x(i) = dble ( i )
      end do

      do i = 1, 10
        y(i) = dble ( 10 * i )
      end do

      do i = 1, 5
        do j = 1, 5
          a(i,j) = dble ( 10 * i + j )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Y = '
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A = '
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,5f8.2)' ) ( a(i,j), j = 1, 5 )
      end do

      call dcopy ( 5, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DCOPY ( 5, X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      do i = 1, 10
        y(i) = dble ( 10 * i )
      end do

      call dcopy ( 3, x, 2, y, 3 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DCOPY ( 3, X, 2, Y, 3 )'
      write ( *, '(a)' ) ' '
      do i = 1, 10
        write ( *, '(2x,i6,g14.6)' ) i, y(i)
      end do

      call dcopy ( 5, x, 1, a, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DCOPY ( 5, X, 1, A, 1 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A = '
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,5f8.2)' ) ( a(i,j), j = 1, 5 )
      end do

      do i = 1, 5
        do j = 1, 5
          a(i,j) = dble ( 10 * i + j )
        end do
      end do

      call dcopy ( 5, x, 2, a, 5 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DCOPY ( 5, X, 2, A, 5 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A = '
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,5f8.2)' ) ( a(i,j), j = 1, 5 )
      end do

      return
      end
      subroutine ddot_test ( )

c*********************************************************************72
c
cc DDOT_TEST tests DDOT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda
      integer ldb
      integer ldc

      parameter ( n = 5 )
      parameter ( lda = 10 )
      parameter ( ldb = 7 )
      parameter ( ldc = 6 )

      double precision a(lda,lda)
      double precision b(ldb,ldb)
      double precision c(ldc,ldc)
      integer i
      integer j
      double precision ddot
      double precision sum1
      double precision x(n)
      double precision y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DDOT_TEST'
      write ( *, '(a)' ) '  DDOT computes the dot product of '
      write ( *, '(a)' ) '  double precision vectors.'

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = - dble ( i )
      end do

      do i = 1, n
        do j = 1, n
          a(i,j) = dble ( i + j )
        end do
      end do

      do i = 1, n
        do j = 1, n
          b(i,j) = dble ( i - j )
        end do
      end do
c
c  To compute a simple dot product of two vectors, use a
c  call like this:
c
      sum1 = ddot ( n, x, 1, y, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Dot product of X and Y is ', sum1
c
c  To multiply a ROW of a matrix A times a vector X, we need to
c  specify the increment between successive entries of the row of A:
c
      sum1 = ddot ( n, a(2,1), lda, x, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Product of row 2 of A and X is ', sum1
c
c  Product of a column of A and a vector is simpler:
c
      sum1 = ddot ( n, a(1,2), 1, x, 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  Product of column 2 of A and X is ', sum1
c
c  Here's how matrix multiplication, c = a*b, could be done
c  with DDOT:
c
      do i = 1, n
        do j = 1, n
          c(i,j) = ddot ( n, a(i,1), lda, b(1,j), 1 )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Matrix product computed with DDOT:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,5g14.6)' ) ( c(i,j), j = 1, n )
      end do

      return
      end
      subroutine dnrm2_test ( )

c*********************************************************************72
c
cc DNRM2_TEST tests DNRM2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lda

      parameter ( n = 5 )
      parameter ( lda = n + 5 )
c
c  These parameters illustrate the fact that matrices are typically
c  dimensioned with more space than the user requires.
c
      double precision a(lda,lda)
      integer i
      integer incx
      integer j
      double precision dnrm2
      double precision sum1
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DNRM2_TEST'
      write ( *, '(a)' ) '  DNRM2 computes the Euclidean norm of '
      write ( *, '(a)' ) '  a double precision vector.'
c
c  Compute the euclidean norm of a vector:
c
      do i = 1, n
        x(i) = dble ( i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The vector X:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,f8.4)' ) i, x(i)
      end do
      incx = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The 2-norm of X is ',
     &  dnrm2 ( n, x, incx )
c
c  Compute the euclidean norm of a row or column of a matrix:
c
      do i = 1, n
        do j = 1, n
          a(i,j) = dble ( i + j )
        end do
      end do

      incx = lda
      sum1 = dnrm2 ( n, a(2,1), incx )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The 2-norm of row 2 of A is ', sum1

      incx = 1
      sum1 = dnrm2 ( n, a(1,2), incx )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The 2-norm of column 2 of A is ', sum1

      return
      end
      subroutine drot_test ( )

c*********************************************************************72
c
cc DROT_TEST tests DROT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 6 )

      double precision c
      integer i
      double precision s
      double precision x(n)
      double precision y(n)

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = dble ( i * i - 12 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DROT_TEST'
      write ( *, '(a)' ) '  DROT carries out a double precision '
      write ( *, '(a)' ) '  Givens rotation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      c = 0.5D+00
      s = sqrt ( 1.0D+00 - c * c )
      call drot ( n, x, 1, y, 1, c, s )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a,f8.4,a)' )
     &  '  DROT ( N, X, 1, Y, 1, ', c, ',', s, ' )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = dble ( i * i - 12 )
      end do

      c = x(1) / sqrt ( x(1) * x(1) + y(1) * y(1) )
      s = y(1) / sqrt ( x(1) * x(1) + y(1) * y(1) )
      call drot ( n, x, 1, y, 1, c, s )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a,f8.4,a)' )
     &  '  DROT ( N, X, 1, Y, 1, ', c, ',', s, ' )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      return
      end
      subroutine drotg_test ( )

c*********************************************************************72
c
cc DROTG_TEST tests DROTG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      double precision a
      double precision b
      double precision c
      double precision r8_uniform_01
      double precision r
      double precision s
      double precision sa
      double precision sb
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 5 )
      double precision z

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  DROTG generates a real Givens rotation'
      write ( *, '(a)' ) '    (  C  S ) * ( A ) = ( R )'
      write ( *, '(a)' ) '    ( -S  C )   ( B )   ( 0 )'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        a = r8_uniform_01 ( seed )
        b = r8_uniform_01 ( seed )

        sa = a
        sb = b

        call drotg ( sa, sb, c, s )

        r = sa
        z = sb

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6,a,g14.6)' ) '  A =  ', a,  '  B =  ', b
        write ( *, '(a,g14.6,a,g14.6)' ) '  C =  ', c,  '  S =  ', s
        write ( *, '(a,g14.6,a,g14.6)' ) '  R =  ', r,  '  Z =  ', z
        write ( *, '(a,g14.6)' ) '   C*A+S*B = ',  c * a + s * b
        write ( *, '(a,g14.6)' ) '  -S*A+C*B = ', -s * a + c * b

      end do

      return
      end
      subroutine dscal_test ( )

c*********************************************************************72
c
cc DSCAL_TEST tests DSCAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 6 )

      double precision da
      integer i
      double precision x(n)

      do i = 1, n
        x(i) = dble ( i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DSCAL_TEST'
      write ( *, '(a)' ) '  DSCAL multiplies a double precision scalar'
      write ( *, '(a)' ) '  times a vector.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X = '
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do

      da = 5.0D+00
      call dscal ( n, da, x, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DSCAL ( N, ', da, ', X, 1 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do

      do i = 1, n
        x(i) = dble ( i )
      end do

      da = -2.0D+00
      call dscal ( 3, da, x, 2 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DSCAL ( 3, ', da, ', X, 2 )'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6)' ) i, x(i)
      end do

      return
      end
      subroutine dswap_test ( )

c*********************************************************************72
c
cc DSWAP_TEST tests DSWAP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 6 )

      integer i
      double precision x(n)
      double precision y(n)

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = dble ( 100 * i )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DSWAP_TEST'
      write ( *, '(a)' ) '  DSWAP swaps two vectors.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      call dswap ( n, x, 1, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DSWAP ( N, X, 1, Y, 1 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = dble ( 100 * i )
      end do

      call dswap ( 3, x, 2, y, 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.4,a)' ) '  DSWAP ( 3, X, 2, Y, 1 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X and Y'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,g14.6,g14.6)' ) i, x(i), y(i)
      end do

      return
      end
      subroutine idamax_test ( )

c*********************************************************************72
c
cc IDAMAX_TEST tests IDAMAX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 11 )

      integer i
      integer i1
      integer incx
      integer idamax
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'IDAMAX_TEST'
      write ( *, '(a)' ) '  IDAMAX returns the index of maximum '
      write ( *, '(a)' ) '  magnitude;'

      do i = 1, n
        x(i) = dble ( mod ( 7 * i, 11 ) ) - dble ( n / 2 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The vector X:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,f8.4)' ) i, x(i)
      end do

      incx = 1

      i1 = idamax ( n, x, incx )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The index of maximum magnitude = ', i1

      return
      end
