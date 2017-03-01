/* Global Scope: Parameter Estimation, ML
 * Scope:
 * Solve a second-order linear (for the moment) equation
 * coefficients could be a polynom up to degree two.
 * The collocation method (with Chebyshev interpolation)
 * and Lanczos method is used
 * Ref: Wright (60's), Lanczos (1957), Clenshaw (50's)
 *
 * LIEBHERR TOULOUSE
 *******************************************************/
#include "../include/lhTypes.hpp"
#include <cassert>
#include <iostream>
#include <fstream>
#include <malloc.h>
#include <cstdio>
#include <cstring>

#ifndef __IDENT05_COLL_HPP__
#define __IDENT05_COLL_HPP__

using namespace lhlib;
//using namespace std;

class IDENT05_COLL {

   public:
	//constructor:
      IDENT05_COLL(Index typode, Index ord, const char * FileNameInp);
      // destructor:
      ~IDENT05_COLL();

      //Parse the data text file:
      bool read_parse_file();

      // sets the linear system of relation:
      bool ExpandSeriesLinearSys_ref1();
      
   private:
      	   // Ordinary Differential equation
      Index type_eqn;  // 1. Linear2 single coeff 2. Linear2 with polynom up to deg 2; 
      char * strFileNameInp;
      Index lenFileNameInp;
	   //Chebyshev polynoms
      Number t0[7];
      Number t1[7];
      Number t2[7];
      Number t3[7];
      Number t4[7];
      Number t5[7];
      Number t6[7];
	  //dimension of the problem
      Index order;      //order of interpolation, also called N in comments
      Index num_ranges;  //also called ma
      Index num_total_coeffs; //shall be initialized to N*ma

      // Content from InputData read
      Number * d_t;
      Number * d_ut;
      Number * d_yt;

};
#endif
