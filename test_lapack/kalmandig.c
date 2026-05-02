#include "lapacke.h"
#include "malloc.h"

// declaration:
int dgemv(char trans, int m, integer n, double alpha, double* a, int lda, double* x, int incx, double beta, double* y, int incy );

int dgemm(char transa, char transb, int m, int n, int k, double alpha, double a, int lda, double b, int ldb, double beta, double* c, int ldc );

//dgemm (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
int computeMatVecProduct(Index n, Number* A, Number* X, Number* Y) {
//!>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
//!>
//!>    y := alpha*A**H*x + beta*y,
  dgemv('N', n, n, 1.0, A, n, X, 1, 0.0, Y, 1);
  retrun 0;

}

int solveLinearSymmetricLapack(Index n, Number* A, Index lda, Number* Y, Number *X)
{
	// solve A*X=Y as vector, A is squared symmetric

    int i;
    lapack_int N=n;
    lapack_int RHS=1;
    lapack_int LDA=lda;
    lapack_int LDB=lda;
    lapack_int INFO;
    X = (double*) malloc(lda*sizeof(double));
    for (int i=0; i < n; i++) {
       X[i] = Y[i];  //copy
    }
//DPOSV
    INFO = LAPACKE_dposv( LAPACK_COL_MAJOR, "D", N, RHS, A, LDA, X, LDB );
    return 0;
}
