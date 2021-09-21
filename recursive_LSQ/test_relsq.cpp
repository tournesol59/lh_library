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
  li_doubles datalist, lsqoutlist;  // exchange list of doubles

  // zero ctl
  for (k=0; k<Npty; k++) {
     ctldata.push_back(0.0);
  }

  // generate Npty points over line 0.8+0.0*t, t=1...n, added with randoms numbers
  // between -0.5 and 0.5, (slope 0.0 to test the variance of the Cpp random generator (see line 57)
  generate_rdom_example(outdata, Npty, 0.8, 0.0, 0.5);

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

  // create a instance of class which manages to approx "data"
  // with size 30 points that are given,...
  IDENT05_RELSQ lsqClassInst=IDENT05_RELSQ(outdata, ctldata, Npty, 1, 0.1, 0.5);

  // perform least square 1st order, 30 iteration
  for (k=1; k<Npty; k++) { // NEVER BEGIN at k=0 (div by zero)
      lsqClassInst.recursive_algorithm(k,1);
  }
  
  //create another export class
  char explsqFileName[14];
  strncpy(explsqFileName, "a.lsqo", 7);
  IDENT05_IODATA explsqClassInst=IDENT05_IODATA(Npty, explsqFileName);

  // copy output data into list for export
  std::string str_data="myvar";  // varianz of y to test the Cpp random generator
  lsqClassInst.pass_iodata(lsqoutlist, str_data);
  // "approx value of varianz"
  explsqClassInst.exportToDisk(lsqoutlist);
	
  return 0;
}


