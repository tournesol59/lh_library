#include <iostream>
#include <math.h>
#include <vector>

class matrixop {

   int n,m;
   std::vector<double> coeff;

 public:
   // Constructeurs et operateur de copie:
   matrixop(int dn, int dm, std::vector<double> dcoeff);
   matrixop(const matrixop &);
   matrixop &operator=(const matrixop &);
   ~matrixop();
   
   //fonctions allowing to read the double values coeff
   std::vector<double> getcoeff(void) const;
   int getrows(void) const;  // return n
   int getcols(void) const; // return m

   // basic operators. This particular overriding of op={+,-,*} when the
   //  code appear as:
   //    A = A op B
   //  is to be understood as
   //  A.opOperand(B)
   //  that is the result of the operation will be placed in the left operand
   //   for matrix mult, an internal copy operation will be performed before

   matrixop &operator+=(const matrixop & );
   matrixop &operator-=(const matrixop & );
   matrixop &operator*=(const matrixop & );
   
};

