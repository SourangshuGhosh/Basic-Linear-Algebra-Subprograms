# BLAS2_C

## Single Precision Complex Matrix-Vector Basic Linear Algebra Subprograms

BLAS2_C, a FORTRAN77 library which constitutes the Level 2 Basic Linear Algebra Subprograms (BLAS), for matrix-vector operations using single precision complex arithmetic.

The BLAS are a small core library of linear algebra utilities, which can be highly optimized for various architectures.

## Licensing:
The computer code and data files described and made available on this web page are distributed under the MIT license.

## Reference:
1.Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, Sven Hammarling, Alan McKenney, Danny Sorensen,
LAPACK User's Guide,
Third Edition,
SIAM, 1999,
ISBN: 0898714478,
LC: QA76.73.F25L36.

2.Thomas Coleman, Charles vanLoan,
Handbook for Matrix Computations,
SIAM, 1988,
ISBN13: 978-0-898712-27-8,
LC: QA188.C65.

3.Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
LINPACK User's Guide,
SIAM, 1979,
ISBN13: 978-0-898711-72-1,
LC: QA214.L56.

4.Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
Algorithm 539: Basic Linear Algebra Subprograms for Fortran Usage,
ACM Transactions on Mathematical Software,
Volume 5, Number 3, September 1979, pages 308-323.

## Source Code:
- **blas2_c.f**, the source code;
- **blas2_c.sh**, commands to compile the source code;

## Examples and Tests:
- **blas2_c_prb.f**, a calling program;
- **blas2_c_prb.sh**, commands to compile, link, load and run the calling program;
- **blas2_c_prb_output.txt**, the output file.

## List of Routines:
- **CGBMV** computes y := alpha * A * x + beta * y, A a complex band matrix.
- **CGEMV** computes y := alpha * A * x + beta * y, A a general complex matrix.
- **CGERC** performs the rank 1 operation A := A + alpha * x * conjg ( y' ).
- **CGERU** performs the rank 1 operation A := A + alpha * x * y'.
- **CHBMV** performs the matrix-vector operation y := alpha * A * x + beta * y.
- **CHEMV** performs the matrix-vector operation y := alpha * A * x + beta * y.
- **CHER2** performs the hermitian rank 2 operation
- **CHER** performs the hermitian rank 1 operation A := A + alpha*x*conjg( x' ).
- **CHPMV** performs the matrix-vector operation y := alpha*A*x + beta*y.
- **CHPR2** performs the hermitian rank 2 operation
- **CHPR** performs the hermitian rank 1 operation A := A + alpha*x*conjg( x' ).
- **CTBMV** performs one of the matrix-vector operations
- **CTBSV** solves one of the systems of equations
- **CTPMV** performs one of the matrix-vector operations
- **CTPSV** solves one of the systems of equations
- **CTRMV** performs one of the matrix-vector operations
- **CTRSV** solves one of the systems of equations
