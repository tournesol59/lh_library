/* Scope: parse the first of the two input file
 * LIEBHERR TOULOUSE
 *******************************************************/
#include "../include/ident05_coll.hpp"
#include <iomanip>
#include <sstream>

#ifndef __TEST_COLL_ONLY__
#define __TEST_COLL_ONLY__
#endif
//using namespace std;
using namespace lhlib;

bool IDENT05_COLL::subSaveRangesRes() {
   Index k=0;
   Number numval1, numval2;
   Number kcoeffs[10];
   Number x;
   std::ofstream fspec ("TheSpec", std::ofstream::out);

   k=0;
   while (k < num_ranges) {
      fspec << "INT:RNGORD\n";
      fspec << k << "\t" << order << "\n"; // order is already known but useful for reimport	   
      // recopy spectral coeffs:
      for (int i=0; i < order+1; i++) {
         kcoeffs[i] = coeffarray[(order+1)*k+i];
      }
      // export these coeffs to the file
      for (int i=0; i < (order/2); i++) {
        // prepare to output two coefficients per lines
        numval1 = kcoeffs[2*i];
        numval2 = kcoeffs[2*i+1];
        fspec << "DOU:LDCAL" << (2*i) << "\n";
        fspec << numval1;
        fspec << "\t" << numval2 << "\n"; 
      }
      // evaluate solution at end point
      x=-1.;
      numval1 = evalCollocation(k, x, kcoeffs);
      x=1.;
      numval2 = evalCollocation(k, x, kcoeffs);
      fspec << "DOU:BVALUY" << "\n";
      fspec <<  numval1;
      fspec << "\t" << numval2 << "\n";
      // evaluate derivative at end point
      x=-1.;
      numval1 = evalDerivCollocation(k, x, kcoeffs);
      boundary_derivall[2*k] = numval1;  // it must be here
      x=1.;
      numval2 = evalDerivCollocation(k, x, kcoeffs);
      boundary_derivall[2*k+1] = numval2; // fred: same thing
      fspec << "DOU:BVALUD" << "\n";
      fspec << numval1;
      fspec << "\t" << numval2 << "\n";         //
      numval1 = 0.0;
      numval2 = 0.0;
      // zero to complete when iterative approach is taken in the future
      fspec << "DOU:BVJACA" << "\n";
      fspec <<  numval1 << "\t" << numval2 << "\n";
      fspec << "DOU:BVJACB" << "\n";
      fspec <<  numval1 << "\t" << numval2 << "\n";

      k++; // next interval
   } // loop over intervals

   fspec.close(); 
   return 1;
}

bool IDENT05_COLL::subSaveGraphRes() {
  std::ofstream lh_test;
  lh_test.open("TheTest", std::ofstream::out);  // test w/o variable Name
  int row=0;
  while (row<num_ranges*num_points) {
//  while (lh_test.good()) {
      lh_test << solarray[row][0] << "\t";
      lh_test << solarray[row][1] << "\n";
      row++;
  }
  lh_test.close();
  return 1;
}
