#include "mex.h"
#include "matrix.h"
//#include "collocation/ident05_coll.hpp"
#include "collocation/calc.hpp"

//double *mxGetPr(const mxArray *pm);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
/*
** ex:[outvector1, scalardiff]=mexFunction('filename', invector1, size1)
*/


{
   int i;
// get input parameters of the mex function
   std::string datfilename=readMATLABString(prhs[0]);

   double *invec = mxgetPr(prhs[1]);

   const int *dim = mxGetDimensions(prhs[1]);

   const int *order = mxGetPr(prhs[2]);

// prepare outputs
   mwSize outdim[1] = {dim[1]};

   plhs[0] = mxCreateNumericArray(1, outdim, mxSINGLE_CLASS, mxREAL);

   double *data = (double *)mxGetPr(plhs[0]);

   double *ptr=&data[0];
//void demo_mexfun(double * out, const double * in, const int size) 
   demo_mexfunc(ptr, invec, dim[1]);

   for (i=1;i<dim[1];i++) {
       
       data[i]=ptr[i];
 
   }
}

