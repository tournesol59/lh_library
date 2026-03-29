#include "decl_relsq.hpp"

#ifndef  _RELSQ_TEST_FIXED_MATDIM_
#define  _RELSQ_TEST_FIXED_MATDIM_
#endif

/* method of class IDENT05_EXTLSQ */

    IDENT05_EXTLSQ::IDENT05_EXTLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dna, int dnb, int dnd, double fTs, double fvarian) : IDENT05_ABSLSQ(ydata, udata, dsize, dna, fTs, fvarian) 
{
  nb=dnb;
  nd=dnd;
  estvarian=0.0;
  for (int j=0; j<nd; j++) {
     werr.push_back(0.0); // initialization of vector copied
  }
/* to be completed by comparing with the latest version of IDENT05_BASICLSQ, yes */

}
    IDENT05_EXTLSQ::~IDENT05_EXTLSQ()
{

}

bool IDENT05_EXTLSQ::innovation(double &epsilon) 
{
   int i;
#ifdef _RELSQ_TEST_FIXED_MATDIM_
// calculation of the iteration with the procedures above for 3x3 matrices
   for (i=0; i<(na+nb+nd); i++) {
       theta.at(i) = theta.at(i) + matrixK.at(i)*epsilon;
   }
#else
// calculation of the iteration with the boost lib (general)

#endif
	return 0;
}

 
bool IDENT05_EXTLSQ::predict(double yiter, double &epsilon) 
{
   int j,flag;
   double ypred=0.0;
#ifdef _RELSQ_TEST_FIXED_MATDIM_
// calculation of the iteration with the procedures above for 3x3 matrices
   flag = multiply_scalprod(phi, theta, ypred, na+nb);
#else
// calculation of the iteration with the boost lib (general)
//   epsilon = inner_prod (phi, theta);
#endif
   for (j=0; j<nd-1; j++) { // shift data
	   werr.at(j)=werr.at(j+1);
   }
   werr.at(nd-1)=(yiter-ypred);  
   for (j=0; j<nd; j++) {
         epsilon=epsilon + theta.at(na+nb+j)*werr.at(nd-j); //calculate e(k)=D(q^^-1)*werr(k)
      }
   yest.push_back(ypred);// needed to monitor the procedure

   return flag;
}

bool  IDENT05_EXTLSQ::pass_iodata(std::vector<double> &dataoutvec)
{
    int row=0; 
    // TBC
// skeleton
  //  auto it=dataoutvec.begin();
     while (row<size)  {
        dataoutvec.push_back(yhat[row]); // warning about what is really calculte in yhat
        row++;
//	it++;
   }//end while
  return 0;
}

