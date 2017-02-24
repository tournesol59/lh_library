#include <iostream>

#include "../include/lhTypes.hpp"
#include "../galerkin/ident05_gakm.hpp"
#include <malloc.h>

using namespace std;

int main() {
  ifstream a_file ( "../optim/ident06_data" );
  double * t1,x1;
  int i;

  t1=malloc(10*sizeof(Number));
  x1=malloc(10*sizeof(Number));
 
  cout << "Hello World!\n" ;
  
  IDENT05_GAKM classgakm = IDENT05_GAKM(6);
  for (i=0;i < 10; i++) {
       t1[i]=double(i)/10.0;
  }
  if (classgakm::eval_leg_pol(6,10,t1, x1) == 0 ) {
 
      for (i=0; i < 10; i++) {
            printf("%f  %f\n", t1[i], x1[i] );
      }
  }
  else {
      printf("error evaluation");
      exit(1);
  }

  a_file.close();

}
