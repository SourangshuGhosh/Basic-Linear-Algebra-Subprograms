      program main

c*********************************************************************72
c
cc MAIN is the main program for BLAS0_PRB.
c
c  Discussion:
c
c    BLAS0_PRB tests the BLAS0 library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 March 2017
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'BLAS0_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BLAS0 library.'

      call dmach_test ( )
      call test01 ( )
      call test015 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'BLAS0_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine dmach_tesst ( )

c*********************************************************************72
c
cc DMACH_TEST tests DMACH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision dmach

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DMACH_TEST'
      write ( *, '(a)' ) '  DMACH computes several machine-dependent'
      write ( *, '(a)' ) '  double precision arithmetic parameters.'

      write ( *, '(a)' ) ' '
      write ( *, * ) '  DMACH(1)  = machine epsilon = ', dmach ( 1 )
      write ( *, * ) '  DMACH(2)  = a tiny value    = ', dmach ( 2 )
      write ( *, * ) '  DMACH(3)  = a huge value    = ', dmach ( 3 )

      return
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests R4_ABS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 March 2014
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real r4
      real r4_abs
      real r4_absolute
      real r4_hi
      real r4_lo
      real r4_uniform_ab
      integer seed
      integer test
      integer test_num

      r4_hi = 5.0E+00
      r4_lo = -3.0E+00
      seed = 123456789
      test_num = 10

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  R4_ABS returns the absolute value of an R4.'
      write ( *, '(a)' ) ''

      do test = 1, test_num
        r4 = r4_uniform_ab ( r4_lo, r4_hi, seed )
        r4_absolute = r4_abs ( r4 )
        write ( *, '(2x,f10.6,2x,f10.6)' ) r4, r4_absolute
      end do

      return
      end
      subroutine test015 ( )

c*********************************************************************72
c
cc TEST015 tests R4_SIGN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2014
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 5 )

      real r4_sign
      integer test
      real x
      real x_test(test_num)

      save x_test

      data x_test /
     & -1.25E+00, -0.25E+00, 0.0E+00, +0.5E+00, +9.0E+00 /

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST015'
      write ( *, '(a)' ) '  R4_SIGN returns the sign of a number.'
      write ( *, '(a)' ) ''

      do test = 1, test_num
        x = x_test(test)
        write ( *, '(2x,f8.2,2x,f8.2)' ) x, r4_sign ( x )
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests R8_ABS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2014
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision r8
      double precision r8_abs
      double precision r8_absolute
      double precision r8_hi
      double precision r8_lo
      double precision r8_uniform_ab
      integer seed
      integer test
      integer test_num

      r8_hi = 5.0D+00
      r8_lo = -3.0D+00
      seed = 123456789
      test_num = 10

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  R8_ABS returns the absolute value of an R8.'
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '      X         R8_ABS(X)'
      write ( *, '(a)' ) ''

      do test = 1, test_num
        r8 = r8_uniform_ab ( r8_lo, r8_hi, seed )
        r8_absolute = r8_abs ( r8 )
        write ( *, '(2x,f10.6,2x,f1.06)' ) r8, r8_absolute
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests R8_SIGN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2014
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 5 )

      double precision r8_sign
      integer test
      double precision x
      double precision x_test(test_num)

      save x_test

      data x_test /
     & -1.25D+00, -0.25D+00, 0.0D+00, +0.5D+00, +9.0D+00 /

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  R8_SIGN returns the sign of a number.'
      write ( *, '(a)' ) ''

      do test = 1, test_num
        x = x_test(test)
        write ( *, '(2x,f8.4,2x,f8.4)' ) x, r8_sign ( x )
      end do

      return
      end
