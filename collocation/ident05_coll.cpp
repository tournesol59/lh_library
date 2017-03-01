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
#include "ident05_coll.hpp"

#ifndef __TEST_COLL_ONLY__
#define __TEST_COLL_ONLY__
#endif

//using namespace std;
using namespace lhlib;

/** attemps to access optimization parameters via macros */

//#define P1_EPS (ident05_nlp::get_p1_eps)
//#define P2_OM (ident05_nlp::get_p2_om)

// constructor
IDENT05_COLL::IDENT05_COLL(Index typode, Index ord, const char * FileNameInp )
{
  Index i;

/*eqn Problem type: */
  type_eqn=typode;

/* Input File Name */
  lenFileNameInp = 15;
  std::strcpy(FileNameInp, strFileNameInp, lenFileNameInp);

/*  t0=(Number)malloc(7*sizeof(Number)); */
  t0[0]=1.0;t0[1]=0.0;t0[2]=0.0;t0[3]=0.0;t0[4]=0.0;t0[5]=0.0;t0[6]=0.0;
  t1[0]=0.0;t1[1]=1.0;t1[2]=0.0;t1[3]=0.0;t1[4]=0.0;t1[5]=0.0;t1[6]=0.0;
  t2[0]=-1.0;t2[1]=0.0;t2[2]=2.0;t2[3]=0.0;t2[4]=0.0;t2[5]=0.0;t2[6]=0.0;
  t3[0]=0.0;t3[1]=-3.0;t3[2]=0.0;t3[3]=4.0;t3[4]=0.0;t3[5]=0.0;t3[6]=0.0;
  t4[0]=1.0;t4[1]=0.0;t4[2]=-8.0;t4[3]=0.0;t4[4]=8.0;t4[5]=0.0;t4[6]=0.0;
  t5[0]=0.0;t5[1]=5.0;t5[2]=0.0;t5[3]=-20.0;t5[4]=0.0;t5[5]=16.0;t5[6]=0.0;
  t6[0]=-1.0;t6[1]=0.0;t6[2]=18.0;t6[3]=0.0;t6[4]=-48.0;t6[5]=0.0;t6[6]=32.0;

  ord=order;
  num_ranges=1;
  num_total_coeffs=(order+1)*num_ranges;

}

IDENT05_COLL::~IDENT05_COLL() 
{
}
bool IDDENT05_COLL::read_parse_file() 
{
   Index i;
   char * ReadLine[39];
/* Open the input file, which is a text file with 3 columns */
  ifstream lh_file("ThisFileInp.in");  // test w/o variable Name
  lh_file.open();
  for (i=0; i<43; i++) {
     lh_file >> ReadLine;
     cout << "i= " << i << "  " << ReadLine;
  }
   lh_file.close();
}

bool IDENT05_COLL::ExpandSeriesLinearSys_ref1()
{
   Index i,j,k;

   return 0;
}

