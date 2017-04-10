//#include "MathFunctions.h"
#include <stdio.h>
#include <malloc.h>
#include "/home/frederic/lapack-3.7.0/LAPACKE/include/lapacke.h"


//#define SIZE 9
int mySolveLinearLapack(int n, int rhs, double * A, int lda, int * ipiv, double * B, int ldb, int info) {

 // char JOBU;
 // char JOBVT;

  int i;
 // integer JOB=n;
  lapack_int N=n;
  lapack_int RHS=rhs;
  lapack_int LDA=lda;
  lapack_int * IPIV;
  lapack_int LDB=ldb;
  lapack_int INFO;
 
  IPIV=malloc(n*sizeof(lapack_int));

  for (i=0;i<n;i++) {
    IPIV[i]=ipiv[i];
  }
//  dgesv_(&N, &RHS, A, &LDA, IPIV, B, &LDB, &INFO);
  INFO=LAPACKE_dgesv(LAPACK_COL_MAJOR,N,RHS,A,LDA, IPIV, B,LDB);
//  for (i=0; i<n; i++) {
//     printf("\n b[ %d ] = %f", i, B[i] );
//  }
  if (INFO==1) {
    printf("LAPACK ERROR: dgesv has not found a solution\n");
    return -1; 
 }
 else {
    return 0;
 }
}
