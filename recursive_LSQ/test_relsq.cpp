#include <iostream>
#include "../include/ident05_fftdata.hpp" 
#include "decl_relsq.hpp"
#include <vector>

//using namespace std;
using namespace lhlib;

int main() {
  ////////////////////////////////////
  // Test an simple constant slope noise-corrupted signal: determine by least square the dynamic estimation of fit-values (slope or acceleration) 	
  int Npty = 30;
  int k;
  std::pair<double,double> singleton; // TBRD
  std::vector<double> ctldata, outdata, ar_params;
  li_doubles datalist, lsqoutlist;  // exchange list of doubles // check TBD

  // select type of orig. data: slope(0) or AR1(1)
  int dselect;
  std::cout << "Please select random source source (0: slope, 1: AR(1), 2: from file)";
  std::cin >> dselect;
  // zero ctl
  for (k=0; k<Npty; k++) {
     ctldata.push_back(0.0);
  }

  // generate Npty points over line 0.8+0.0*t, t=1...n, added with randoms numbers
  // between -0.5 and 0.5, (slope 0.0 to test the variance of the Cpp random generator (see line 57)
  if (dselect==0) {
     generate_rdom_example(outdata, Npty, 0.8, 0.0, 0.5);
  } else if (dselect==1) {
     ar_params.push_back(0.875);
     generate_arN_example(outdata, Npty, 1, ar_params, 0.5);
  } else if (dselect==2) {
     generate_from_file(outdata, Npty, "secondorder.dat", 3, 2);
  }

  // export class
  char expFileName[14];
  strncpy(expFileName, "a.rand", 7);
  IDENT05_IODATA expClassInst=IDENT05_IODATA(Npty, 0.1, expFileName); //do not forget 0.1=sampling

  // copy data into list and save to Disk
  k=0; 
  for (auto it=outdata.begin(); it != outdata.end(); it++) {
     singleton.first = (double (k))/10.0;
     singleton.second = (*it);
     datalist.push_back( singleton );
     std::cout << (*it) << " ";
     k++;
  }
  std::cout << "\n";
  expClassInst.exportToDisk(datalist);

  // create a instance of class which manages to approx "data"
  // with size 30 points that are given,...
  IDENT05_RELSQ lsqClassInst=IDENT05_RELSQ(outdata, ctldata, Npty, 1, 0.1, 0.5);
  // perform differentiation with small filter, not ideal but
  for (k=1; k<30; k++) { // 
      lsqClassInst.apply_derivative_filter(k, 0.9355, 1);
  }
  // perform least square 1st order, 30 iteration
  for (k=1; k<Npty; k++) { // NEVER BEGIN at k=0 (div by zero)
      lsqClassInst.recursive_algorithm(k,1);
  }
  
  //create another export class
  char explsqFileName[14];
  strncpy(explsqFileName, "a.lsqo", 7);
  IDENT05_IODATA explsqClassInst=IDENT05_IODATA(Npty, 0.1, explsqFileName);

  // copy output data into list for export
  std::string str_data="myvar";  // varianz of y to test the Cpp random generator
  lsqClassInst.pass_iodata(lsqoutlist, str_data);  //! here a modification of dpair has TBC
  // "approx value of varianz"
  explsqClassInst.exportToDisk(lsqoutlist);

  ////////////////////////////////////
  // Test derivated classes: 
  return 0;
}


