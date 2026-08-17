#include <iostream>
#include "decl_matrixop.h"

void print_matrix(const matrixop Mat) {
   for (int i=0; i<Mat.n; i++) {
      for (int j=0; j<Mat.m; j++) {
	      std::cout << Mat(i,j);
      }
      std::cout << "\n";
   }
   std::cout << "\n";
}

int main(int argc, char **argv)  {

   matrixop A(2,3);
   A(1,1)=1.0; A(1,2)=-0.5; A(1,3)=0.;
   A(2,1)=0.5; A(2,2)=1.0 ; A(2,3)=-0.25;
   print_matrix(A);

   matrixop B(3,2);
   B(1,1)=0.0; B(1,2)=0;
   B(2,1)=1.0; B(2,2)=0.;
   B(3,1)=0.0; B(3,2)=1.0;
   print_matrix(B);

 //  matrixop C(A);
 //  matrixop C=A;
 //  C *= B;
 //  print_matrix(C);
   print_matrix(A);
   return 0;
}
