/*
 *
*/
#include <iostream>
#include <cstring>
#include <vector>
#include <list>
#include <iterator>
#include "../include/lhTypes.hpp"

using namespace lhlib;
////////////////
// Class Definition
 class IDENT05_OPTLSQ {
      public:
   IDENT05_OPTLSQ(std::vector<double> dxinit, std::vector<double> drhs, int ddim, double ftol);
   IDENT05_OPTLSQ(const IDENT05_OPTLSQ &source);
   IDENT05_OPTLSQ &operator=(const IDENT05_OPTLSQ &source);
   ~IDENT05_OPTLSQ(void);

   bool fres(std::vector<double> xd, double &res);
   bool Jac(std::vector<double> xd, std::vector<double> &jac);
   bool update(std::vector<double> xd, std::vector<double> jac, std::vector<double> res, std::vector<double> &newxd);

      protected:
   std::vector<double> theta;
   std::vector<double> xd;
   std::vector<double> jacobian;
   std::vector<double> xinit;
   std::vector<double> rhs;
   int dim;
   double tol;
   int Itermax;
   int iter;
};
