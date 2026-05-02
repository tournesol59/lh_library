/* Scope: parse the first of the two input file
 * LIEBHERR TOULOUSE
 *******************************************************/
#include "../include/ident05_coll.hpp"
//#include <cstring>

#ifndef __TEST_COLL_ONLY__
#define __TEST_COLL_ONLY__
#endif
//using namespace std;
using namespace lhlib;

bool IDENT05_COLL::read_parse_file() 
{
   Index row,col,row_guess;

// Open the input file, which is a text file with 6 columns
  std::ifstream lh_file;
  lh_file.open(strFileNameInp, std::ifstream::in);  // test w/o variable Name
  
  row=0;
  row_guess=0; // prescaler output of row % num_points
  while ((!lh_file.eof()) && (row<41)) {
  /**
  * complete here
  **/
         std::string line;
	 std::getline(lh_file, line);
         std::stringstream ss(line);
         col = 0;	 
	 while (ss >> dataarray[row][col]) col++; // fait la detection automatique chaine vide pour s'arreter
#ifdef __TEST_COLL_ONLY__
	 std::cout << "tf= " << dataarray[row][0] << " ";
	 std::cout << "utf= " << dataarray[row][1] << " ";
	 std::cout << "ytf= " << dataarray[row][2] << " ";
	 std::cout << "c2f= " << dataarray[row][3] << " ";
	 std::cout << "c1f= " << dataarray[row][4] << " ";
	 std::cout << "c0f= " << dataarray[row][5] << "\n";
	 
#endif
         if (row==row_guess*(num_points)-1) {
             yguessed[row_guess][0]=dataarray[row][0];
             yguessed[row_guess][1]=dataarray[row][2];
	     row_guess++;
         }
	 row++;
   }
   num_rows=row;
   lh_file.close();

   return 0;
}

inline double lh_atof(const char *str)
{
      while ((*str) && isspace(*str))
       ++str;
      return strcmp(str, "nan") ? atof(str) : NAN;
}

// Exercice: do a little bit differently for the next function
bool IDENT05_COLL::read_parse_specfile() {
   //std::ifstream lh_fspec;
   char head[256];
   std::string headstr;
   char flag[256];
   std::string flagstr;
   char fl[3];
   char value[256];
   std::string valuestr;
   
   Index validx1;
   Index validx2;
   Number valued1;
   Number valued2;
   Index row=0;
   Index kint=-1; // needed since first operation is to incremment it
   Index ksave=0;
   Index ccal=0;
   Index maxrow=num_ranges*(8*2);  // TBC if more params coming to be exported

   std::cout << "****** test reread export specfile ****** \n";
   //lh_fspec.open(strSpecNameInp, std::ios::in);
   std::ifstream lh_fspec(strSpecNameInp, std::ios::in);
 
   //lh_fspec.getline(headstr, 256, '\n');
   while ((getline(lh_fspec, headstr)) && (row < maxrow)) {
      
      if ((row % 2)==0) { // even line
	      //copy the string into flag for matching wildcard
         strncpy(flag, headstr.c_str(), 10);
	 if (!strncmp(flag, "INT:RNGORD", 10)) {
            kint++; // needed below for wildcard "DOU:LDCAL0"
	 }
      }
      else {  // odd line
           //copy the string into value for parsing int or double
	 if ((!strncmp(flag, "INT", 3)) && (getline(lh_fspec, valuestr))) {
            std::stringstream iss(valuestr);
            iss >> validx1 >> validx2; // test
	 }
         else if (!strncmp(flag, "DOU", 3)) {
            std::stringstream iss(valuestr);
	    iss >> valued1 >> valued2; // test
	 }
      }
      // NOW store the imported value in the right place depending on wildcard
      if (!strncmp(flag, "INT:RNGORD", 10)) {
	 // here simple verfication of order
	 if (validx2 != order) {
	    std::cout<< "mismatch of parameter: order=" << order << " in output file\n";
	 }
      }
      else if (!strncmp(flag, "DOU:LDCAL", 9)) {
	 coeffarray_prev[(order+1)*kint+2*ccal] = valued1;
	 coeffarray_prev[(order+1)*kint+2*ccal+1] = valued2;
         ccal++; // there are order%2 lines of two coefficients (per lines)
	 ccal=(2*(ccal+1)==order)? 0 : ccal;  // preview for next interval
      }
      else if (!strncmp(flag, "DOU:BVALUY", 10)) {
         boundry_prev[0] = valued1;
	 boundry_prev[1] = valued2;
      }
      else if (!strncmp(flag, "DOU:BVALUD", 10)) {
         boundary_derivprev[0] = valued1;
	 boundary_derivprev[1] = valued2;
      }
      else if (!strncmp(flag, "DOU:BVJACA", 10)) {
         boundary_jacobian[0] = valued1;
	 boundary_jacobian[1] = valued2;
      }
      else if (!strncmp(flag, "DOU:BVJACB", 10)) {
         boundary_jacobian[2] = valued1;
	 boundary_jacobian[3] = valued2;
      }
 }
#ifdef __TEST_COLL_ONLY__
      if (kint == (ksave+1)) {
	      // we are just at the beginning of a new interval => display the prev:
         std::cout << "Reread output spec (coeff) file for interv k=" << kint <<" :\n";
         std::cout << "y(a)=" << boundry_prev[0] << ", y(b)=" << boundry_prev[1] << "\n";
	 std::cout << "y'(a)=" << boundary_derivprev[0] << ", y'(b)=" << boundary_derivprev[1] << "\n";
	 //std::cout << "Ddy/da(a)=" << boundary_jacobian[0] << ",Ddy/db(a)=" << boundary_jacobian[1] << ",Ddy/da(b)=" << boundary_jacobian[3] << ",Ddy/db(b)=" << boundary_jacobian[0] << "\n";
      }
#endif
   lh_fspec.close();
   return 0;
}
