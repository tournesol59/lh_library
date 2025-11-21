// complement of impl_basicrelsq.cpp, implements methods of the class
// IDENT05_MANLSQ

#include "impl_basiclsq.hpp"

// this constructor inherits from two different base classes
 IDENT05_MANLSQ::IDENT05_MANLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dna, double fTs, double fvarian, std::vector<double> dxinit, std::vector<double> drhs, ddim, ftol) : IDENT05_ABSLSQ(ydata, udata, dsize, dna, fTs, fvarian) , IDENT05_OPTLSQ(dxinit, drhs, ddim, ftol)
 
{
   for (int i=0; i<dim; i++) {
      theta[i]=dxinit[i];
      xd[i]=dxinit[i];
   }
}

 IDENT05_MANLSQ::IDENT05_MANLSQ(const IDENT05_MANLSQ &source) : IDENT05_ABSLSQ(source),  IDENT05_OPTLSQ(source)
{
	
}

IDENT05_MANLSQ &IDENT05_MANLSQ :: operator =(const IDENT05_MANLSQ &source)
{

}

IDENT05_MANLSQ:: ~IDENT05_MANLSQ(void)
{

}

 bool IDENT05_MANLSQ:: algorithm() {

   std::vector<double> xd;
   std::vector<double> res;
   bool flag=true;
   flag=fres(xd, res);
   
   while ((iter <= Itermax) && (fabs(res[0]) > tol) && (fabs(res[1]) > tol) && (flag == true)) {

     flag=fres(xd, res);

     flag=Jac(xd, jacobian);
   
     flag=update(xd, jacobian, res, theta);
// prepare next step
     for (int i=0; i<dim; i++) {
        xd[i]=theta[i];
     }
     iter++;
   }
   return 1;
}

  
