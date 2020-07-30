# BLAS1_D
## Basic Linear Algebra Subprograms
## Level 1
## Double Precision Real Arithmetic

BLAS1_D a FORTRAN77 library which constitutes the Level 1 Basic Linear Algebra Subprograms (BLAS), for vector-vector operations using double precision real arithmetic, by Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh.

The BLAS are a small core library of linear algebra utilities, which can be highly optimized for various architectures. Software that relies on the BLAS is thus highly portable, and will typically run very efficiently.

## Licensing:
The computer code and data files described and made available on this web page are distributed under the MIT license.

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
- **blas1_d.f**, the source code;

## Examples and Tests:
- **blas1_d_prb.f**, a calling program;
- **blas1_d_prb_output.txt**, the output file.

## List of Routines:

- **DASUM** takes the sum of the absolute values.
- **DAXPY** computes constant times a vector plus a vector.
- **DCOPY** copies a vector.
- **DDOT** forms the dot product of two vectors.
- **DMACH** computes machine parameters of floating point arithmetic.
- **DNRM2** returns the euclidean norm of a vector.
- **DROT** applies a plane rotation.
- **DROTG** constructs a Givens plane rotation.
- **DROTM** applies a modified Givens rotation matrix.
- **DROTMG** generates a modified Givens rotation matrix.
- **DSCAL** scales a vector by a constant.
- **DSDOT** computes the inner product of two vectors with extended precision.
- **DSWAP** interchanges two vectors.
- **IDAMAX** finds the index of element having maximum absolute value.
- **LSAME** returns TRUE if CA is the same letter as CB regardless of case.
- **XERBLA** is an error handler.
