# BLAS0
Auxilliary Functions for Basic Linear Algebra Subprograms (BLAS)
BLAS0 is a FORTRAN77 library which contains auxilliary functions for the Basic Linear Algebra Subprograms (BLAS).

## Licensing:
The computer code and data files made available on this web page are distributed under the MIT license by Sourangshu Ghosh

## Author:

Sourangshu Ghosh

## Reference:
Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
Algorithm 539: Basic Linear Algebra Subprograms for Fortran Usage,
ACM Transactions on Mathematical Software,
Volume 5, Number 3, September 1979, pages 308-323.

## Source Code:
blas0.f, the source code.
blas0.sh, BASH commands to compile the source code.

## Examples and Tests:
blas0_prb.f, a sample calling program.
blas0_prb.sh, BASH commands to compile and run the sample program.
blas0_prb_output.txt, the output file.

## List of Routines:
- C4_UNIFORM_01 returns a unit double precision complex pseudorandom number.
- C4MAT_PRINT prints a C4MAT.
- C4MAT_PRINT_SOME prints some of a C4MAT.
- C4MAT_TEST sets up a test matrix.
- C4MAT_TEST_INVERSE returns the inverse of the test matrix.
- C8_UNIFORM_01 returns a unit double precision complex pseudorandom number.
- C8MAT_PRINT prints a C8MAT.
- C8MAT_PRINT_SOME prints some of a C8MAT.
- C8MAT_TEST sets up a test matrix.
- C8MAT_TEST_INVERSE returns the inverse of the test matrix.
- CABS1 returns the L1 norm of a single precision complex number.
- CABS2 returns the L2 norm of a single precision complex number.
- CMACH computes single precision complex machine parameters.
- CSIGN1 is a single precision complex transfer-of-sign function.
- CSIGN2 is a single precision complex transfer-of-sign function.
- DMACH computes machine parameters of floating point arithmetic.
- LSAME returns TRUE if CA is the same letter as CB regardless of case.
- R4_ABS returns the absolute value of an R4.
- R4_SIGN returns the sign of an R4.
- R4_UNIFORM_01 returns a unit single precision pseudorandom number.
- R4_UNIFORM_AB returns a scaled pseudorandom R4.
- R4MAT_PRINT prints an R4MAT.
- R4MAT_PRINT_SOME prints some of an R4MAT.
- R4MAT_TEST sets up a test matrix.
- R4VEC_PRINT prints an R4VEC.
- R8_ABS returns the absolute value of an R8.
- R8_SIGN returns the sign of an R8.
- R8_UNIFORM_01 returns a unit double precision pseudorandom number.
- R8_UNIFORM_AB returns a pseudorandom R8 scaled to [A,B].
- R8MAT_PRINT prints an R8MAT.
- R8MAT_PRINT_SOME prints some of an R8MAT.
- R8MAT_TEST sets up a test matrix.
- R8VEC_PRINT prints an R8VEC.
- SMACH computes machine parameters of floating point arithmetic.
- TIMESTAMP prints out the current YMDHMS date as a timestamp.
- XERBLA is an error handler.
- ZABS1 returns the L1 norm of a double complex number.
- ZABS2 returns the L2 norm of a double complex number.
- ZMACH computes double complex floating point arithmetic constants.
- ZSIGN1 is a double precision complex transfer-of-sign function.
- ZSIGN2 is a double precision complex transfer-of-sign function.
