#include <iostream>
#include "decl_matrixop.h"
#include <vector>

//----- constructors -----

matrixop::matrixop(unsigned short int nl, unsigned short int nc)
{
   n=nl;
   m=nc;
   lignes = new ligne[n];
   for (unsigned short int i=0; i < n; i++) {
      lignes[i] = new double[m];
   } 
   return;
}

matrixop::matrixop(const matrixop &source)
{
   n=source.n;
   m=source.m;
   lignes = new ligne[n];
   for (int i=0; i < n; i++) {
      lignes[i] = new double[m];
      for (int j=0; j < m; j++) {
         lignes[i][j]=source.lignes[i][j];
      }
   } 
   return;
}


matrixop &matrixop::operator=(const matrixop &source)
{
   if (&source != this)
   {
      if ((source.n!=n) || (source.m!=m)) {
         for (unsigned short int i=0; i<n; i++) {
	    delete[] lignes[i];
         }
         delete[] lignes;
         n=source.n;
         m=source.m;
	 lignes = new ligne[n];
         for (int i=0; i<n; i++) {
            lignes = new ligne[n];
         }
      }
      for (int i=0; i < n; i++) {
         for (int j=0; j < m; j++) {
            lignes[i][j]=source.lignes[i][j];
         }
      }
  } 
   return *this;
}

matrixop::~matrixop()
{
}

   //functions to write access the values
double &matrixop::operator() (unsigned short int i, unsigned short int j) {
    return lignes[i][j];
}
   // other function to read access function
double matrixop::operator() (unsigned short int i, unsigned short int j) const {
    return lignes[i][j];
}

// tip to transform a vector to a matrixop, not a constructor
void matrixop::setValues(std::vector<double> values) {
   if (n*m==values.size()) {
      for (int i=0; i<n; i++) {
          for (int j=0; j<m; j++) {
              lignes[i][j]=values.at(i*m+j);
	  }
      }
   }
}

// tips: permit the export to a vector<double> container
std::vector<double> matrixop::getValues() {
	std::vector<double> temp = std::vector<double>(n*m, 0.0);
    for (int i=0; i<n; i++) {
       for (int j=0; j<m; j++) {
           temp[i*m+j]=lignes[i][j];
       }
    }
    return temp;
}


// ------ override internal operators ------//

matrixop &matrixop::operator+=(const matrixop & matB)
{
  if ((n != matB.n) || (m != matB.m)) {
     std::cout << "Error not matching dimensions for op addition";
  }
  else {
   for (int i=0; i<n; i++) {
     for (int j=0; j<m; j++) {
        lignes[i][j]  += matB(i,j);
     }
   }
  }
  return *this;
}

matrixop &matrixop::operator-=(const matrixop & matB)
{
  if ((n != matB.n) || (m != matB.m)) {
     std::cout << "Error not matching dimensions for op addition";
  }
  else {
   for (int i=0; i<n; i++) {
     for (int j=0; j<m; j++) {
        lignes[i][j]  -= matB(i,j);
     }
   }
  }
  return *this;
}

// following makes a cpy firstly
void matrixop::multiply(const matrixop & matB) {
  matrixop matC=matrixop(*this); // FRED: a tester

  int i,j,k;
  double sum;

  if (m != matB.n) {
     std::cout << "Error dimension mismatch in multiplication of matrices";
  }
  else {
   for (i=0; i<n; i++) {
     for (j=0; j<n; j++) {
	sum=0.0;
	for (k=0; k<n; k++) {
            sum=sum+lignes[i][k] + matB(k,j);
	}
        matC(i,j)=sum;
     }
    }
   // free matC
   *this = matC;
   }
}

matrixop &matrixop::operator*=(const matrixop & matB)
 {
  this->multiply(matB);
  return *this;
 }

// the basic: allow a multiplication by a scalar
void matrixop::multscal(double scalar) {
   for (int i=0; i<n; i++) {
      for (int j=0; j<m; j++) {
	 lignes[i][j] *= scalar;
      }
   }
}

// tip: allow a scalar product with a vector container

double matrixop::vec_scal_product(std::vector<double> coeffs) {
   double sum=0.0;
   for (int i=0; i<n; i++) {
      sum += lignes[i][0] * coeffs.at(i);
   }
   return sum;
}

