#include <iostream>
#include "./galerkin/ident05_gakm.hpp"

using namespace std;


int main() {
  ifstream a_file ( "ident06_data" );
  
  cout << "Hello World!\n" ;
  
  IDENT05_GAKM classgakm = IDENT05_GAKM(6);
  float x1= classgakm.eval_leg_pol(6, 1.0);  // order goes up to 6

  a_file.close();

}
