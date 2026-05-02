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

#include "../../lapack_inst/include/lapacke.h" // fred: path shall now be indicated in makefile
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

#define PEXIT_ERR_IMPORTCODE 1
#define PEXIT_ERR_IMPORTFILE 2
#define PEXIT_ERR_EXPORTTEST 3
#define PEXIT_ERR_EXPORTSPEC 4
#define PEXIT_ERR_EXACT_NULLDIV 5

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
      Number tinit;
      Number tend;      
      Number boundry[2]; // same struct for ivp or bvp;
      Number predictparams[4];  // FRED, WAS: [2], IS:[4] one extra params for exponential decay 
      Index enableload;
      Index numload;     
      Number coeffload[(7+1)];

} IDENT05_OPTIONS;

class IDENT05_COLL {

   public:
	//constructor:
      IDENT05_COLL(Index ord, const char * FileNameInp, const char * CodeNameInp , const char * SpecNameInp);
      // destructor:
      ~IDENT05_COLL();

      //Parse the data text file ("TheFile"):
      bool read_parse_file();
      //Parse the other data text file ("TheCode"):
      bool read_parse_code();
      //Parse the result data text outputed from subSaveRangesRes
      bool read_parse_specfile();
      // set the rhs by inversion of Taylor Expansion

      bool InverseRHSLinearSys_ref1();
      // sets the linear system of relation:
      bool ExpandSeriesLinearSys_ref1();
      // solve it with LAPACKE (for the moment) after: IPOPT
      bool SolveSeriesLinearSys_ref1();
      // repeat the operation ma times and adjusting for each interval:
	// - coeff of equ1[] because of changing t variable
	// - boudary conditions boundry[](tk+1)=solution(boundry[](tk), t=tk+1)
      bool NumRangesCalcBoundary(Index k);
      bool SolveNumRangesSys_ref1();
	// only for test purposes: evaluate the  solution and display the solution

      bool UpdateBoundaryIterative();  // deprecated
      bool UpdateBoundaryIterativeSlope(); // new
      bool UpdateBoundaryIterativePhase(); // new
	bool exploreForExtremumFirstOdd(int k, int &found1, int &found2, double &yabs, double &dyabs, double &dphi); 
      bool exploreForExtremumSecEven(int k, int &found1, int &found2, double &yabs, double &dyabs, double &dphi); 

      bool subSaveRangesRes();  // subroutine that is called by SolveNumRangesSys_ref1()
      bool subSaveGraphRes();  // subroutine that is called by SolveNumRangesSys_ref1()
      bool subEraseCodeInp(Index k);

      bool pass_dataarray_col(Index dim, li_doubles &exports);
      Number evalCollocation(Index k, Number t, Number * coeff);
      Number evalChebyshevPolynom(Number t, Index i);
      Number evalDerivCollocation(Index k, Number t, Number * coeff);
      Number evalDerivChebyshevPolynom(Number t, Index i);      
      bool evalExactSolution(Index kitv);

   private:

      Number equ1[3];  // equ[0]*y + equ[1]*y' + equ[2]*y''=0
      Number equ2[9];  // (equ[0]+equ[1]*t+equ[2]*t^2)*y + (equ[3]+equ[4]*t+equ[5]*t^2)*y' + (equ[6]+equ[7]*t+equ[8]*t^2)*y'' = 0  not implemented up to the present

      char strFileNameInp[14];
      char strCodeNameInp[14];
      Index lenFileNameInp;
      // char strTestNameOut[14];
      char strSpecNameInp[14];

      // Chebyshev polynoms
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
      Index Miter;
      // Ordinary Differential equation parameters
      IDENT05_OPTIONS code_opts;

      Index type_eqn;  // 1. Linear2 single coeff 2. Linear2 with polynom up to deg 2;
      Index type_rhs; // 0=cosinus, 1=sinus
      Number omega_rhs; 
      Index type_predict;  //as boolean, use the B*sin(wt) func to generate boundry[1] value
      Index repeat_predict; //(1) or use the same value[1] for all ranges      
      Index type_ovp; // 1 if boundary value problem 0, if initial
      Number times_end[10];
      Number boundry[2]; // if bvp: xa and xb
      Number boundry_prev[2]; // same as boundry but imported value
      Number boundary_deriv[2]; // if bvp: x'(a) and x'(b)
      Number boundary_derivprev[2]; // same as boundary_deriv but imported value
      Number boundary_jacobian[4]; // for future dev: jacobian of boundary_deriv % boundry 
      Number boundary_all[10*2]; // load or predict boundaries, convention is:
             // _all[2*k] -> y(a) as begin [k*T,(k+1)*T]
	     // _all[2*k+1] -> y(b) as end [k*T,(k+1)*T]
	     // _all[2*(k+1)] -> y(b) as begin [(k+1)*T,(k+2)*T] redondant
	     // _all[2*(k+1)+1] -> y(c) as end [(k+1)*T,(k+2)*T] iterate ..
      Number boundary_derivall[10*2]; // store result derivatives computed convention:
             // _derivall[2*k] -> y'(a) as begin [k*T, (k+1)*T]
	     // _derivall[2*k+1] -> y'(b) as end of [k*T, (k+1)*T]
             // derivall[2*(k+1)] -> y'(b) as begin of [(k+1)*T, (k+2)*T] => MATCH condition
	     // derivall[2*(k+1)+1] -> y'(c) as end of [(k+1)*T, (k+2)*T]

      Index extremum_alllocus[10*2];  // location (index in of two first extremums (one for y, one for y' in solarray) per interval
      Number extremum_allval[10*2];  // found extremum values (y and y') per interval
      Number extremum_alldphi[10*2];  // analysis of phase increase per interval
      Number initial[2]; // or if ivp: x0 and dx0
      Number tinit;
      Number tend;
      Index num_ranges;  //also called ma, should be equal to 4
      Index num_rows_file; // total rows in "TheFile" problem definition data
      Index num_points;  // points per range to display   
      Index num_rows;
      Index num_total_coeffs; //shall be initialized to N*ma
      Number predictparams[4];

      Number A_l1[81];  // matrix (order+3)*(order+3)
      Number B_l1[9];   // matrix RHS (order+3)
      Number B_l1_rhs[9]; // matrix RHS for an individual source
	// this shall be solved ma times (the number of intervals)

      Number coeffarray[8*(6+1)];  //Chebyshev coefficients for the principal solution dim=[num_total_coeffs]
      Number coeffarray_prev[8*(6+1)];
//      Number coeffdistarray[4*(6+1)];  //Chebyshev coefficients for the perturbation solution (in the future Van der Pol)

      Index enableload;
      Index numload;     
      Number coeffload[7+1]; // Chebyshev coefficients for the right hand side for the first interval only, loaded at runtime

      // Content from InputData read
      Number dataarray[2000][6]; //  0 = t, 1 = ut, 2 = yt, 3=c2, 4=c1, 5=c0
      Number yguessed[100][2];      
      Number solarray[2000][3]; //yt_approx
};

#endif
