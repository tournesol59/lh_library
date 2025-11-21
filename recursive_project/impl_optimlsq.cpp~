#include "decl_optimlsq.hpp"

IDENT05_OPTLSQ::IDENT05_OPTLSQ(std::vector<double> xinit, std::vector<double> rhs, int ddim, double ftol) :
	theta(std::vector<double>(ddim, 0.0)),
	xd(dxinit),
	xinit(dxinit),
	jacobian(std::vector<double>(ddim^2,0.0)),
	rhs(drhs),
	dim(ddim),
	tol(ftol)
{
   iter=0;
   Itermax=100;
}	

IDENT05_OPTLSQ:: ~IDENT05_OPTLSQ(void)
{
}

bool IDENT05_OPTLSQ:: fres(std::vector<double> xd, std::vector<double> &res) {
// sum elevated to square only if needed: two counter convex func, here it is not the case
   if (dim==2) {
      res[0] = ( rhs[0]*(xd[0]^2 + xd[1]^2 + 1) - xd[0]-xd[0]*xd[1] );
      res[1] =  (rhs[1]*(xd[0]^2 + xd[1]^2 + 1) - xd[1]);
   }
   return 1;
};
   
bool IDENT05_OPTLSQ:: Jac(std::vector<double> xd, std::vector<double> &jac) {
   if (dim==2) {
      jac[0] = (rhs[0])*(2*xd[0])-xd[1]-1.0;
      jac[1] = (rhs[0])*(2*xd[1])-xd[0];
      jac[2] = (rhs[1])*(2*xd[0]);
      jac[0] = (rhs[1])*(2*xd[1])-1.0;
   }
   return 1;
}

bool IDENT05_OPTLSQ::update(std::vector<double> xd, std::vector<double> jac, std::vector<double> res, std::vector<double> &newxd) {
   // newxd = xd - res./(Df)  voir si ca marche 
   std::vector<double> invjac;
   double det = (jac[0]*jac[3]-jac[1]*jac[2]);
   if (fabs(det) > 1e-5) {
	   // compute inverse
      invjac.push_back(jac[3]/det);
      invjac.push_back(-jac[1]/det);
      invjac.push_back(-jac[2]/det);
      invjac.push_back(jac[0]/det);
           // matrix wise operation
      newxd[0] = xd[0] - invjac[0]*res[0] - invjac[1]*res[1];
      newxd[1] = xd[1] - invjac[2]*res[0] - invjac[3]*res[1];
      return 1;
   }
   else {
      return -1;
   }
   
}
