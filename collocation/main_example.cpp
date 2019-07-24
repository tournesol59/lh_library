/********************************************************************************************
* The aim is to test the collocation method, this time not through Ipopt on the equations
*    y'' + 1/2*y' + 16*y = 0  (1) with y(0)=0 y'(0)=1
* then
*    ty''+ y' + t*y = u(t)          (2) with y(0)=1 y'(0)=0
*
*********************************************************************************************/

#include "ident05_coll.hpp"
#include <memory>  // for Smart Ptr: unique_ptr.h

using namespace lhlib;

int main(int argc, char ** argv) {
  char strFileName[14];
  char strCodeName[14];
  int lenFileName=14;

  strncpy(strFileName, argv[1], 8);  // 8 < 14 but works do not touch anything!
//  strncpy(strCodeName, argv[2], 8);  // same thing
//  strFileName = "TheFile";  //prohibited, call $ ./main_example "TheFile" "TheCode"  instead
  std::cout << "Have string " << strFileName << " of length " << lenFileName << " as arg[1] and " << strCodeName << " as arg[2]\n";

  IDENT05_COLL cnClassInst=IDENT05_COLL(6,strFileName,strCodeName);
  cnClassInst.read_parse_file();
  cnClassInst.read_parse_code();

// If Testing: Decompose:
  cnClassInst.ExpandSeriesLinearSys_ref1();
  cnClassInst.SolveSeriesLinearSys_ref1();

  cnClassInst.SolveNumRangesSys_ref1();
}

