# BLAS1_S
## Basic Linear Algebra Subprograms
## Level 1
## Single Precision Real Arithmetic

BLAS1_S, a FORTRAN77 library which constitutes the Level 1 Basic Linear Algebra Subprograms (BLAS), for vector-vector operations using single precision real arithmetic, by Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh.

The BLAS are a small core library of linear algebra utilities, which can be highly optimized for various architectures. Software that relies on the BLAS is thus highly portable, and will typically run very efficiently.

## Licensing:
The computer code and data files described and made available on this web page are distributed under the MIT license.

## Author:
Original FORTRAN77 version by Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh.

## Reference:

1.Thomas Coleman, Charles vanLoan,
Handbook for Matrix Computations,
SIAM, 1988,
ISBN13: 978-0-898712-27-8,
LC: QA188.C65.

2.Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
LINPACK User's Guide,
SIAM, 1979,
ISBN13: 978-0-898711-72-1,
LC: QA214.L56.

3.Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
Algorithm 539: Basic Linear Algebra Subprograms for Fortran Usage,
ACM Transactions on Mathematical Software,
Volume 5, Number 3, September 1979, pages 308-323.

## Source Code:
- **blas1_s.f**, the source code;
- **blas1_s.sh**, commands to compile the source code;

## Examples and Tests:
- **blas1_s_prb.f**, a calling program;
- **blas1_s_prb.sh**, commands to compile, link, load and run the calling program;
- **blas1_s_prb_output.txt**, the output file.

## List of Routines:
- **ISAMAX** finds the index of element having maximum absolute value.
- **LSAME** returns TRUE if CA is the same letter as CB regardless of case.
- **SASUM** takes the sum of the absolute values.
- **SAXPY** computes constant times a vector plus a vector.
- **SCOPY** copies a vector, x, to a vector, y.
- **SDOT** forms the dot product of two vectors.
- **SDSDOT** computes the inner product of two vectors with extended precision.
- **SMACH** computes machine parameters of floating point arithmetic.
- **SNRM2** returns the euclidean norm of a real vector.
- **SROT** applies a plane rotation.
- **SROTG** constructs a Givens plane rotation.
- **SROTM** applies a modified Givens transformation matrix.
- **SROTMG** constructs a modified Givens transformation matrix
- **SSCAL** scales a vector by a constant.
- **SSWAP** interchanges two vectors.
- **XERBLA** is an error handler.
