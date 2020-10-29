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
#include "../MathFunctions/MathFunctions.hpp"
#include "ident05_fftdata.hpp"
//extern "C" int mySolveLinearLapack(int n, int rhs, double * A, int lda, int * ipiv, double * B, int ldb, int info);

//#include "../../../CLAPACK/INCLUDE/clapack.h"
//#include <stdio.h>
#include <cassert>
#include <iostream>
#include <istream>
#include <sstream>
#include <fstream>
#include <malloc.h>
#include <cstdio>
#include <vector>  // for read_parse_file
#include <array>
#include <cstring>
#include <math.h>

#ifndef __IDENT05_COLL_HPP__
#define __IDENT05_COLL_HPP__

//constants macros
#define COL_TYP_INITIAL 0
//#define COL_TYP_BVP 1
#define COL_TYP_BVP_PRED_TRICKED 0
//#define COL_TYP_BVP_PRED_NORM 1
#define COL_TYP_BVP_NREPEAT 0
//#define COL_TYP_BVP_NSAMEONE 1

using namespace lhlib;
//using namespace std;

typedef struct ident05_options {
// duplets lines (same type per line) will be read from options code file
// in method read_parse_code 
      Index type_eqn;  // 1. Linear2 single coeff 2. Linear2 with polynom up to deg 2; 
      Index type_ovp; // 1 if boundary value problem 0, if initial
      Index type_predict;  //as boolean, use the B*sin(wt) func to generate boundry[1] value
      Index repeat_predict; //(1) or use the same value[1] for all ranges
      Index num_ranges;  //also called ma, should be equal to 4
      Index num_points;  // points per range      
      Number boundry[2]; // same struct for ivp or bvp
      Number tinit;
      Number tend;
      Number predictparams[2];      

} IDENT05_OPTIONS;

class IDENT05_COLL {

   public:
	//constructor:
      IDENT05_COLL(Index ord, const char * FileNameInp, const char * CodeNameInp);
      // destructor:
      ~IDENT05_COLL();

      //Parse the data text file ("TheFile"):
      bool read_parse_file();
      //Parse the other data text file ("TheCode"):
      bool read_parse_code();

      // sets the linear system of relation:
      bool ExpandSeriesLinearSys_ref1();
      // solve it with LAPACKE (for the moment) after: IPOPT
      bool SolveSeriesLinearSys_ref1();
      // repeat the operation ma times and adjusting for each interval:
	// - coeff of equ1[] because of changing t variable
	// - boudary conditions boundry[](tk+1)=solution(boundry[](tk), t=tk+1)
      bool NumRangesCalcBoundary(Number* boundary_all, Number* times_end, int k);
      bool SolveNumRangesSys_ref1();
	// only for test purposes: evaluate the  solution and display the solution
     bool pass_dataarray_col(Index dim, li_doubles &exports);
      Number evalCollocation(Number t, Number * coeff);
      Number evalChebyshevPolynom(Number t, Index i);

   private:

      Number equ1[3];  // equ[0]*y + equ[1]*y' + equ[2]*y''=0
      Number equ2[9];  // (equ[0]+equ[1]*t+equ[2]*t^2)*y + (equ[3]+equ[4]*t+equ[5]*t^2)*y' + (equ[6]+equ[7]*t+equ[8]*t^2)*y'' = 0  not implemented up to the present

      char strFileNameInp[14];
      char strCodeNameInp[14];
      Index lenFileNameInp;
	   //Chebyshev polynoms
      Number t0[7];
      Number t1[7];
      Number t2[7];
      Number t3[7];
      Number t4[7];
      Number t5[7];
      Number t6[7];

      Number t7[8]; //for completion
	  //dimension of the problem
      Index order;      //order of interpolation, also called N in comments, order=6 is implemented

      // Ordinary Differential equation parameters
      IDENT05_OPTIONS code_opts;

      Index type_eqn;  // 1. Linear2 single coeff 2. Linear2 with polynom up to deg 2; 
      Index type_predict;  //as boolean, use the B*sin(wt) func to generate boundry[1] value
      Index repeat_predict; //(1) or use the same value[1] for all ranges      
      Index type_ovp; // 1 if boundary value problem 0, if initial
      Number boundry[2]; // if bvp: xa and xb
      Number initial[2]; // or if ivp: x0 and dx0
      Number tinit;
      Number tend;
      Index num_ranges;  //also called ma, should be equal to 4
      Index num_points;  // points per range   
      Index num_rows;
      Index num_total_coeffs; //shall be initialized to N*ma
      Number predictparams[2];

      Number A_l1[81];  // matrix (order+3)*(order+3)
      Number B_l1[9];   // matrix RHS (order+3)
	// this shall be solved ma times (the number of intervals)

      Number coeffarray[4*(6+1)];  //Chebyshev coefficients for the principal solution dim=[num_total_coeffs]
//      Number coeffdistarray[4*(6+1)];  //Chebyshev coefficients for the perturbation solution (in the future Van der Pol)

      // Content from InputData read
      Number dataarray[2000][6]; //  0 = t, 1 = ut, 2 = yt, 3=c2, 4=c1, 5=c0
      Number solarray[2000][2]; //yt_approx
};

#endif