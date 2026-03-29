#include "../include/lhTypes.hpp"
#include "../include/ident05_fftdata.hpp"

#if 0
#include "lapackpp.h"
#endif

#if 1
#include "lafnames.h"       /* macros for LAPACK++ filenames */
#include LA_GEN_MAT_DOUBLE_H
#include LA_VECTOR_DOUBLE_H
#include LA_GEN_QR_FACT_DOUBLE_H
#include "blaspp.h"
#include LA_SOLVE_DOUBLE_H
#include LA_GENERATE_MAT_DOUBLE_H
#include LA_EXCEPTION_H
#include LA_UTIL_H
#include "lasvd.h"
#endif


#include <iostream>
#include <cassert>
#include <vector>

// #include "../../boost_1_72_0/libs/??"
#ifndef __SUBDETERM_HPP__
#define __SUBDETERM_HPP__

#define _

class IDENT05_PROJECTIVE
// contains no matrix data nor result, a matrix simple C-array is passed
{
  public:
    Index dimy;      // the size of a vector-like measurement
    Index dimcol;            // common to Yp,Yf,Up,Uf number of column=size of delayed in/outputs
    Index dimrow_y_p; // past measurement inputs number of rows
    Index dimrow_y_f; // future measurement outputs number of rows

    // Index dimrow_all_p; // past masurements inputs+outputs for oblique projection

// constructor    
    IDENT05_PROJECTIVE(Index dy, Index nc, Index ny_1, Index ny);
// destructor
    ~IDENT05_PROJECTIVE();
// estimate extended observability matrix for deterministic case:
    bool LUprojectY_U(Number *Ypast, Number *Upast, Number *Yfuture, Number *Ufuture, Number *Res, Number* L, Number *Q); 
     /* In details: performs the normal projection of future outputs onto
       * past outputs and measurements: it is an error it must be pro
       * jected oblique from future inputs. It uses QR decomp in lieu 
       * of pseudo inverse
       */ 

    // bool LUobl_project_Y_YU_1();

  private:


};

class IDENT05_DETERMSUB 
{
  public:
    Index istart, irange;
    Index nmeasurement;

    IDENT05_DETERMSUB(Index start, Index end, Index nmeas, Number *Ym, Number * Um);

    ~IDENT05_DETERMSUB();

    bool WriteData2Hankel(Number *Y, Number *U);

  private:
    Number *Ypast;
    Number *Upast;
    Number *Yfuture;
    Number *Ufuture;
    Number *Y_on_W;
    Number *Y_on_W_c;

};

#endif

