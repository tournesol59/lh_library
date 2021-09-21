/********************************************************************************************
* The aim is to test the collocation method, this time not through Ipopt on the equations
*    y'' + 1/2*y' + 16*y = 0  (1) with y(0)=0 y'(0)=1
* then
*    ty''+ y' + t*y = u(t)          (2) with y(0)=1 y'(0)=0
*TO DO: after simplification of NumRangesLinearSys
*       link as a library
*********************************************************************************************/

#include "../include/ident05_coll.hpp"
#include "../include/ident05_fftdata.hpp"
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

  // list of import export vectors
  li_doubles datafft;  // li_doubles is a list<double> container
  li_doubles datatimesol;  // idem

  // solver class
  IDENT05_COLL cnClassInst=IDENT05_COLL(6,strFileName,strCodeName);

  // data solution container and class list<T>
  char fftFileName[14];
  strncpy(fftFileName,"fft.out", 8); //..
  IDENT05_IODATA fftClassInst=IDENT05_IODATA(256, fftFileName);

  char expFileName[14];
  strncpy(expFileName,"out.bat", 8); //..
  IDENT05_IODATA expClassInst=IDENT05_IODATA(2000, expFileName); // two signals, two thousands points
	  
  // read problem formulation:
  cnClassInst.read_parse_file();
  cnClassInst.read_parse_code();


// If Testing: Decompose:
  cnClassInst.ExpandSeriesLinearSys_ref1();
  cnClassInst.SolveSeriesLinearSys_ref1();
// or complete auto solving:
  cnClassInst.SolveNumRangesSys_ref1();

    // data exchange between classes
  cnClassInst.pass_dataarray_col( 2000, datatimesol);
  //Index select[2]; //={1,2}
  expClassInst.exportToDisk(datatimesol);
// import fft after being verified
  fftClassInst.read_extern_output(256, datafft);

  return 0;
}

