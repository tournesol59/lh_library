#include <iostream>
#include "../include/ident05_fftdata.hpp" 
#include "decl_relsq.hpp"
#include <vector>

//using namespace std;
using namespace lhlib;

int main() {
  const int Npty = 30;
  int k;
  dpair singleton;
  std::vector<double> ctldata, outdata;
  li_doubles datalist, estoutlist;  // exchange list of doubles

  // zero ctl
  for (k=0; k<Npty; k++) {
     ctldata.push_back(0.0);
  }

  // generate Npty points over line 0.8+0.0533*t, t=1...n, added with randoms numbers
  // between -0.5 and 0.5
  generate_rdom_example(outdata, Npty, 0.8, 0.0833, 0.5);

  // export class
  char expFileName[14];
  strncpy(expFileName, "a.rand", 7);
  IDENT05_IODATA expClassInst=IDENT05_IODATA(Npty, expFileName);

  // copy data into list and save to Disk
  k=0; 
  for (auto it=outdata.begin(); it != outdata.end(); it++) {
     singleton.x = (double (k))/10.0;
     singleton.y = (*it);
     datalist.push_back( singleton );
     std::cout << (*it) << " ";
     k++;
  }
  std::cout << "\n";
  expClassInst.exportToDisk(datalist);

  //create recursive intrument estimation class (FVE)
  reestClassInst=IDENT05_REEST(ydata, udata, dsize, dn1, dn2, dr, fTs, fvar_phi, fvar_y);
     // create child class of transfer functions

     // ...
   
  // begin while loop:
  
  // first: initiate generation of output of model from LSE

  // second: continuous filer to generate component of instr. (psi)

  // third: find an ARMA model for noise (stubb this phase with identity matrix) 
  // fouth: from the last filter generate component of psi

  // fifth: perform LSE-like inversion-based algorithm
  //
  // end while loop

  //create another export class
  char expvarFileName[14];
  strncpy(expvarFileName, "a.lsqo", 7);
  IDENT05_IODATA expvarClassInst=IDENT05_IODATA(Npty, expvarFileName);

  // copy output data into list for export
  std::string str_data="yst0";
  reestClassInst.pass_iodata(estoutlist, str_data);
  explsqClassInst.exportToDisk(estoutlist);
	
  return 0;
}



