/********************************************************************************************
* The aim is to test the collocation method, this time not through Ipopt on the equations
*    y'' + 1/2*y' + 1/4*y = 0  (1) with y(0)=0 y'(0)=1
* then
*    y'' + 4*y = u(t)          (2) with y(0)=1 y'(0)=0
*
*********************************************************************************************/

#include "ident05_coll.hpp"

using namespace lhlib;

int main() {
  char * strFileName[15];
  int lenFileName = 12;
  // open a class as Smart Pointer
  strFileName = "ThisFileInp.in";  
  SmartPtr<IDENT05_COLL> cnClassIns(new IDENT05_COLL(6, strFileName));

  cnClassInst->read_parse_file();
  
  
}

