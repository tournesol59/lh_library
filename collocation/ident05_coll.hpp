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
#include <cstring>
#include <math.h>

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
      // solve it with LAPACKE (for the moment) after: IPOPT
      bool SolveSeriesLinearSys_ref1();
      // repeat the operation ma times and adjusting for each interval:
	// - coeff of equ1[] because of changing t variable
	// - boudary conditions boundry[](tk+1)=solution(boundry[](tk), t=tk+1)
      bool SolveNumRangesSys_ref1();
	// only for test purposes: evaluate the  solution and display the solution
      Number evalCollocation(Number t, Number * coeff);
      Number evalChebyshevPolynom(Number t, Index i);

   private:
      	   // Ordinary Differential equation
      Index type_eqn;  // 1. Linear2 single coeff 2. Linear2 with polynom up to deg 2; 
      Number equ1[3];  // equ[0]*y + equ[1]*y' + equ[2]*y''=0
      Number equ2[9];  // (equ[0]+equ[1]*t+equ[2]*t^2)*y + (equ[3]+equ[4]*t+equ[5]*t^2)*y' + (equ[6]+equ[7]*t+equ[8]*t^2)*y'' = 0  not implemented up to the present
      Number boundry[2]; // x0 and dx0
      
      char strFileNameInp[14];
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
      Number tinit,tend;
      Index num_ranges;  //also called ma, should be equal to 4
      Index num_points;  // points per range
      Index num_total_coeffs; //shall be initialized to N*ma

      Number A_l1[81];  // matrix (order+3)*(order+3)
      Number B_l1[9];   // matrix RHS (order+3)
	// this shall be solved ma times (the number of intervals)

      Number coeffarray[4*(6+1)];  //Chebyshev coefficients for the principal solution dim=[num_total_coeffs]
      Number coeffdistarray[4*(6+1)];  //Chebyshev coefficients for the perturbation solution (in the future Van der Pol)

      // Content from InputData read
      Number dataarray[2000][6]; //  0 = t, 1 = ut, 2 = yt, 3=c2, 4=c1, 5=c0
      Number solarray[2000][2]; //yt_approx
};
#endif
