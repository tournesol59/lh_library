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
   char strSpecFileNAME[14];

   if (Miter==0) {
      strncpy(strSpecFileNAME, "TheSpec0", 8);
   }
   else if (Miter==1) {
      strncpy(strSpecFileNAME, "TheSpec1", 8);
   }
   std::ofstream fspec (strSpecFileNAME, std::ofstream::out);
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
  while (row < num_ranges*num_points) {
//  while (lh_test.good()) {
      lh_test << solarray[row][0] << "\t";
      lh_test << solarray[row][1] << "\n";
      row++;
  }
  lh_test.close();
  return 1;
}

bool IDENT05_COLL::subEraseCodeInp(Index k) {
   Index validx1, validx2;
   Number numval1, numval2;
   Index row=0;
   Index maxrow=8*2;
   char strCodeFileNAME[14];
   
   // manually set the name of the file to erase
   if (k==0) {
      strncpy(strCodeFileNAME, "TheCode1", 8);
   }
   else if (k==1) {
      strncpy(strCodeFileNAME, "TheCode2", 8);
   }
   else if (k==2) {
      strncpy(strCodeFileNAME, "TheCode3", 8);
   }
   // open the file in writing mode
   std::ofstream fspec (strCodeFileNAME, std::ofstream::out);

   //first line: type of the problem
   validx1 = type_ovp;
   validx2 = type_eqn;
   fspec << "INT:IBVEQN:END" << "\n";
   fspec << validx1 << "\t" << validx2 << "\n";

   //second line: logic how to use the predictparams
   validx1 = type_predict;
   validx2 = repeat_predict;
   fspec << "INT:PRELOG:END" << "\n";
   fspec << validx1 << "\t" << validx2 << "\n";
   
   //third line: times of the interval
   numval1 = tinit;
   numval2 = tend;
   fspec << "INT:INIEND:END" << "\n";
   fspec << numval1 << "\t" << numval2 << "\n";

   //fourth line: boundary values
   numval1 = boundary_all[2*k];
   numval2 = boundary_all[2*k+1];
   fspec << "DOU:BVALUE:END" << "\n";
   fspec << numval1 << "\t" << numval2 << "\n";

   //fifth line: num of intervals + num of points per interval
   validx1 = num_ranges;
   validx2 = num_points;
   fspec << "INT:NUMPTS:END" << "\n";
   fspec << validx1 << "\t" << validx2 << "\n";

   //sixth line: predict parameters
   numval1 = predictparams[0];
   numval2 = predictparams[1];
   fspec << "DOU:PREVAL:END" << "\n";
   fspec << numval1 << "\t" << numval2 << "\n";

   //seventh line
   numval1 = predictparams[2];
   numval2 = predictparams[3];
   fspec <<  "DOU:PREEXP:END" << "\n";
   fspec << numval1 << "\t" << numval2 << "\n";

   //eigth line: logic to import output spectral coefficients
   validx1 = enableload;
   validx2 = numload;
   fspec << "INT:ENALDC:END" << "\n";
   fspec << validx1 << "\t" << validx2 << "\n";

   //ninth line: values of the coeff 0 and 1
   numval1 = coeffarray[k*7+0];
   numval2 = coeffarray[k*7+1];
   fspec << "DOU:LDCAL0:END" << "\n";
   fspec << numval1 << "\t" << numval2 << "\n";

   //tenth line: values of the coeff 2 and 3
   numval1 = coeffarray[k*7+2];
   numval2 = coeffarray[k*7+3];
   fspec << "DOU:LDCAL2:END" << "\n";
   fspec << numval1 << "\t" << numval2 << "\n";

   //eleventh line: values of the coeff 4 and 5
   numval1 = coeffarray[k*7+4];
   numval2 = coeffarray[k*7+5];
   fspec << "DOU:LDCAL4:END" << "\n";
   fspec << numval1 << "\t" << numval2 << "\n";

   //twelvth line: values of the coeff 6 and 7
   numval1 = coeffarray[k*7+6];
   numval2 = coeffarray[k*7+7];
   fspec << "DOU:LDCAL6:END" << "\n";
   fspec << numval1 << "\t" << numval2 << "\n";

   fspec.close(); 
   return 1;

}

