      program main

c*********************************************************************72
c
cc MAIN is the main program for BLAS2_C_PRB.
c
c  Discussion:
c
c    BLAS2_C_PRB tests the BLAS library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2014
c
c  Author:
c
c    Sourangshu Ghosh
c
c  License:
c  
c  Released under MIT License by Sourangshu Ghosh
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLAS2_C_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BLAS library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLAS2_C_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests CGEMV.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2014
c
c  Author:
c
c    Sourangshu Ghosh
c
c  License:
c  
c  Released under MIT License by Sourangshu Ghosh
c
      implicit none

      integer m
      parameter ( m = 5 )
      integer n
      parameter ( n = 5 )

      complex a(m,n)
      complex alpha
      complex beta
      integer i
      integer incx
      integer incy
      integer j
      integer lda
      character trans
      complex x(n)
      real x1
      real x2
      complex y(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  For a general matrix A,'
      write ( *, '(a)' ) 
     &  '  CGEMV computes y := alpha * A * x + beta * y'
    
      trans = 'N'
      alpha = ( 10.0E+00, 1.0E+00 ) 
      lda = m
      incx = 1
      beta = 3.0E+00
      incy = 1

      do i = 1, m
        do j = 1, n
          if ( i .eq. j ) then
            a(i,j) = ( 2.0E+00, 0.0E+00 )
          else if ( i .eq. j - 1 .or. i .eq. j + 1 ) then
            a(i,j) = ( -1.0E+00, 0.0E+00 )
          else
            a(i,j) = ( 0.0E+00, 0.0E+00 )
          end if
        end do
      end do

      x1 = 0.0E+00
      x2 = real ( n )
      do i = 1, n
        x(i) = cmplx ( x1, x2 )
        x1 = x1 + 1.0E+00
        x2 = x2 - 2.0E+00
      end do

      do i = 1, m
        y(i) = ( 100.0E+00, 1.0E+00 )
      end do

      call cgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Result vector Y = '
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,2g14.6)' ) y(i)
      end do

      return
      end
