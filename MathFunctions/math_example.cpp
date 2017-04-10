//#include "MathFunctions.h"
extern "C" int mySolveLinearLapack(int n, int rhs, double * A, int lda, int * ipiv, double * B, int ldb, int info);

//#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <istream>

int main() {
 double a[4*4]={1, 2, 3, 4, 1, 4, 9, 16, 1, 8, 27, 64, 1, 16, 81, 256};
 double b[4]={1, 0, 0, 0};
 int n=4;
 int lda=4;
 int ldb=4;
 int rhs=1;
 int ipiv[4]={1,1,1,1};
 int info=0;
 int i;

 if (!mySolveLinearLapack(n, rhs, a, lda, ipiv, b, ldb, info)) {
  
    for (i=0;i<4;i++) {
      std::cout << "b[" << i << "]= " <<  b[i] << "\n";
    }
 }  
 else {
  printf("error of solve linear with lapack\n");
 } 
 return 0;
}
