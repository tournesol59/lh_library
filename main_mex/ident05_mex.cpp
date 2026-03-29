#include "mex.h"
#include "matrix.h"
#include "ident05_wrap.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
// Example:
// [y_solve,t_s]=mexFunction(str_file, str_code,n)
//
   int i;
   const int cdim=40;

// get input parameters of the mex function
   std::string datfilename=readMATLABString(prhs[0]);
   std::string codfilename=readMATLABString(prhs[1]);
   const int *order = mxGetPr(prhs[2]);

   (INARG_S*) invec = (INARG_S *) malloc(1*sizeof(INARG_S));
   invec->order = order;
   
   if (strcpy(invec->filenameinp, datfilename,8)) {
      std::cout << "error copy of string datfilename\n";
   }
   
   if (strcpy(invec->codenameinp, codfilename,8)) {
      std::cout << "error copy of string codfilename\n";
   }

   //future dev:
//   double *invec = mxgetPr(prhs[1]);
//   const int *dim = mxGetDimensions(prhs[1]);
// prepare outputs
   // mwSize outdim[1] = {dim[1]};

   mwSize outdim[1]={cdim};
   plhs[0] = mxCreateNumericArray(1, outdim, mxSINGLE_CLASS, mxREAL);
   plhs[1] = mxCreateNumericArray(1, outdim, mxSINGLE_CLASS, mxREAL);   
// forget   plhs[2]
   double *data = (double *)mxGetPr(plhs[0]);
   double *times = (double *)mxGetPr(plhs[1]);

   double *ptr_y=&data[0];
   double *ptr_t=&times[0];
   //void demo_mexfun(double * out, const double * in, const int size) 

   wrapper_mexfunc_coll(ptr_y, ptr_t, inarg);

   for (i=1;i<cdim;i++) {
       
       data[i]=ptr_y[i];
       times[i]=ptr_t[i];
 
   }
}


