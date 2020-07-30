# BLAS1_C
## Basic Linear Algebra Subprograms
## Level 1
## Single Precision Complex Arithmetic

BLAS1_C, a FORTRAN77 library which constitutes the Level 1 Basic Linear Algebra Subprograms (BLAS), for vector-vector operations using single precision complex arithmetic, by Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh.

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
- **blas1.f**, the source code;
- **blas1_c.sh**, commands to compile the source code;

## Examples and Tests:
- **blas1_c_prb.f**, a calling program;
- **blas1_c_prb.sh**, commands to compile, link, load and run the calling program;
- **blas1_c_prb_output.txt**, the output file.

## List of Routines:

- **CAXPY** computes constant times a vector plus a vector.
- **CCOPY** copies a vector, x, to a vector, y.
- **CDOTC** forms the conjugated dot product of two vectors.
- **CDOTU** forms the dot product of two vectors.
- **CROTG** generates a Givens rotation.
- **CSCAL** scales a vector by a constant.
- **CSSCAL** scales a complex vector by a real constant.
- **ICAMAX** finds the index of element having maximum absolute value.
- **LSAME** returns TRUE if CA is the same letter as CB regardless of case.
- **SCASUM** takes the sum of the absolute values of a complex vector.
- **SCNRM2** returns the euclidean norm of a complex vector.
- **XERBLA** is an error handler.
