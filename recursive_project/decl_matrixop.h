#include <iostream>
#include <math.h>
#include <vector>

class matrixop {

   typedef double *ligne;
   ligne *lignes;

 public:
   unsigned short int n;
   unsigned short int m;
      // Constructeurs et operateur de copie:
   matrixop(unsigned short int dl, unsigned short int nc);
   matrixop(const matrixop &source);
   matrixop &operator=(const matrixop &m1);
   ~matrixop(void);
   
   //functions to write access the values
   double &operator() (unsigned short int i, unsigned short int j);
   // other function to read access function
   double operator() (unsigned short int i, unsigned short int j) const; 
   std::vector<double> getValues();  // instead of a casting
   void setValues(std::vector<double> values);

   // basic operators. This particular overriding of op={+,-,*} when the
   //  code appear as:
   //    A = A op B
   //  is to be understood as
   //  A.opOperand(B)
   //  that is the result of the operation will be placed in the left operand
   //   for matrix mult, an internal copy operation will be performed before

   matrixop &operator+=(const matrixop & );
   matrixop &operator-=(const matrixop & );
   void multiply(const matrixop & matB);
   matrixop &operator*=(const matrixop & );
   void multscal(double scalar);
   double vec_scal_product(std::vector<double> coeffs);

};

