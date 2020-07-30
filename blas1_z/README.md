# BLAS1_Z
## Basic Linear Algebra Subprograms
## Level 1
## Double Precision Complex Arithmetic

BLAS1_Z a FORTRAN77 library which constitutes the Level 1 Basic Linear Algebra Subprograms (BLAS), for vector-vector operations using double precision complex arithmetic, by Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh.

The BLAS are a small core library of linear algebra utilities, which can be highly optimized for various architectures. Software that relies on the BLAS is thus highly portable, and will typically run very efficiently.

## Licensing:
The computer code and data files described and made available on this web page are distributed under the MIT license.

## Author:

Sourangshu Ghosh

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
- **blas1_z.f**, the source code;
- **blas1_z.sh**, commands to compile the source code;

## Examples and Tests:
- **blas1_z_prb.f**, a calling program;
- **blas1_z_prb.sh**, commands to compile, link, load and run the calling program;
- **blas1_z_prb_output.txt**, the output file.

## List of Routines:
- **DZASUM** takes the sum of the absolute values.
- **DZNRM2** returns the euclidean norm of a vector.
- **IZAMAX** finds the index of element having maximum absolute value.
- **LSAME** returns TRUE if CA is the same letter as CB regardless of case.
- **XERBLA** is an error handler.
- **ZAXPY** computes constant times a vector plus a vector.
- **ZCOPY** copies a vector, x, to a vector, y.
- **ZDOTC** forms the dot product of a vector.
- **ZDOTU** forms the dot product of two vectors.
- **ZDSCAL** scales a vector by a constant.
- **ZROTG** generates a Givens rotation.
- **ZSCAL** scales a vector by a constant.
- **ZSWAP** interchanges two vectors.
