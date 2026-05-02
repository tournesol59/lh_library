#include "decl_basiclsq.hpp"

/******* IDENT05_ACFLSQ class ********/
 IDENT05_ACFLSQ :: IDENT05_ACFLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dnp, double fTs, double fvarian, int dnlag) : 
   IDENT05_ABSLSQ(ydata, udata, dsize, dnp, fTs, fvarian),
   rxx_acf(std::vector<double> ((dsize/4),0.0))
   nlag(dnlag)
{

}


 IDENT05_ACFLSQ :: IDENT05_ACFLSQ(const IDENT05_ACFLSQ &source) : IDENT05_ABSLSQ(source), rxx_acf(source.rxx_acf), nlag(source.dnlag)
{
}

 IDENT05_ACFLSQ & :: IDENT05_ACFLSQ operator=(const IDENT05_ACFLSQ &source)
{
}      

 IDENT05_ACFLSQ :: ~IDENT05_ACFLSQ(void)
{
}

 bool IDENT05_ACFLSQ :: sumacf()
{
  double acf_k;
  for (int k=0; k < nlag; k++) {
     acf_k = 0.0;
     for (int i=0; i < (size-k); i++) {
        acf_k = acf_k + y[i+k]*y[i];
     }
     rxx_acf[k] = acf_k;
  }
  return 1;
}


/******* IDENT05_YWLSQ class ********/
 IDENT05_YWLSQ :: IDENT05_YWLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dnp, double fTs, double fvarian, int dnlag) : 
   IDENT05_ACFLSQ(ydata, udata, dsize, dnp, fTs, fvarian, dnlag),
   matrixYw(std::vector<double> ((dnp)*(dnp), 0.0)),
   rhxYw(std::vector<double> (dnp,0.0))
{

}

 IDENT05_YWLSQ :: IDENT05_YWLSQ(const IDENT05_YWLSQ &source) : IDENT05_ACFLSQ(source), matrixYw(source.matrixYw), rhsYw(source.rhsYw)
{
}

 IDENT05_YWLSQ & :: IDENT05_YWLSQ operator=(const IDENT05_YWLSQ &source)
{
}      

 IDENT05_YWLSQ :: ~IDENT05_YWLSQ(void)
{
}

 IDENT05_YWLSQ & :: expandYuleW() 
{
   double rho_k;
   for (int k=0; k < np; k++) {
      rho_k = rxx_acf[k]/rxx_acf[0];
      for (int i=0; i < (np-k); i++) {
         matrixYw[i,i+k] = rho_k;
      }
      for (int i=(np-k); i < np; i++) {
         matrixYw[i,i-np+k] = rho_k; 
      }
      if (k==0) {
        rhsYw[k]=1.0 - varian/rxx_acf[0];
      }
      else {
        rhsYw[k]=rho_k;
      }
   }
}

