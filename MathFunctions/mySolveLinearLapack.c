#include "MathFunctions.h"
#include "/home/frederic/CLAPACK/INCLUDE/f2c.h"
#include <stdio.h>
#include <malloc.h>
#include "/home/frederic/CLAPACK/INCLUDE/clapack.h"


//#define SIZE 9
int mySolveLinearLapack(int n, int rhs, double * A, int lda, int * ipiv, double * B, int ldb, int info) {

 // char JOBU;
 // char JOBVT;

  int i;
 // integer JOB=n;
  integer N=n;
  integer RHS=rhs;
  integer LDA=lda;
  integer * IPIV;
  integer LDB=ldb;
  integer INFO;
 
  IPIV=malloc(n*sizeof(int));

  for (i=0;i<n;i++) {
    IPIV[i]=ipiv[i];
  }
  dgesv_(&N, &RHS, A, &LDA, IPIV, B, &LDB, &INFO);
  
//  for (i=0; i<n; i++) {
//     printf("\n b[ %d ] = %f", i, B[i] );
//  }
  return 0;
}
