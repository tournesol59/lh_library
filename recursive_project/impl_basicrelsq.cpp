#include "decl_basiclsq.hpp"

#ifndef _RELSQ_TEST_FIXED_MATDIM_
#define _RELSQ_TEST_FIXED_MATDIM_
#endif

/**********************************************************
 * procedures
 * *******************************************************/

int copy_matrix(std::vector<double> matA, std::vector<double> &matcpy, int n, int m) {
  for (int i=0; i<n; i++) {
     for (int j=0; j<m; j++) {
        matcpy.at(i*m+j)  = matA.at(i*m+j);
     }
  }
  return 0;
}

int inverse_nthree_matrix(std::vector<double> &mat, std::vector<double> &invmat){
  int i,j;
  double factor, determ;

  for (i=0; i<3; i++) {
     for (j=0; j<3; j++) {
         factor=0.0;
    	 if ((i==0) && (j==0)) {
	    factor=mat.at(3*i+3+j+1)*mat.at(3*i+6+j+2) -
		    mat.at(3*i+3+j+2)*mat.at(3*i+6+j+1);
	 }
	 else if ((i==0) && (j<2)) {
	    factor=mat.at(3*i+3+j-1)*mat.at(3*i+6+j+1) -
		    mat.at(3*i+3+j+1)*mat.at(3*i+6+j-1);
	 }
	 else if ((i==0) && (j==2)) {
	    factor=mat.at(3*i+3+j-2)*mat.at(3*i+6+j-1) -
		    mat.at(3*i+3+j-1)*mat.at(3*i+6+j-2);
	 }
	 else if ((i<2) && (j==0)) {
	    factor=mat.at(3*i-3+j+1)*mat.at(3*i+3+j+2) -
		    mat.at(3*i+3+j+1)*mat.at(3*i-3+j+2);
	 }
	 else if ((i<2) && (j<2)) {
	    factor=mat.at(3*i-3+j-1)*mat.at(3*i+3+j+1) -
		    mat.at(3*i+3+j-1)*mat.at(3*i-3+j+1);
	 }
	 else if ((i<2) && (j==2)) {
	    factor=mat.at(3*i-3+j-2)*mat.at(3*i+3+j-1) -
		    mat.at(3*i+3+j-2)*mat.at(3*i-3+j-1);
	 }
	 else if ((i==2) && (j==0)) {
	    factor=mat.at(3*i-6+j+1)*mat.at(3*i-3+j+2) -
		    mat.at(3*i-3+j+2)*mat.at(3*i-6+j+1);
	 }
	 else if ((i==2) && (j<2)) {
	    factor=mat.at(3*i-6+j-1)*mat.at(3*i-3+j+1) -
		    mat.at(3*i-3+j-1)*mat.at(3*i-6+j+1);
	 }
	 else if ((i==2) && (j==2)) {
	    factor=mat.at(3*i-6+j-2)*mat.at(3*i-3+j-1) -
		    mat.at(3*i-3+j-2)*mat.at(3*i-6+j-1);
	 }
         invmat.at(3*j+i)=factor;

     }
  }
  determ = mat.at(0)*invmat.at(0) + mat.at(1)*invmat.at(3) + mat.at(2)*invmat.at(6);
  for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	 if (i!=j) {
		 invmat.at(3*i+j)=invmat.at(3*i+j)*(-1.0)/determ;
	 }
	 else {
		 invmat.at(3*i+j)=invmat.at(3*i+j)/determ;
	 }
      }
  }

  return 0;
}

int sum_matrices(std::vector<double> matA, std::vector<double> matB, std::vector<double> res, int n, int m) {
  for (int i=0; i<n; i++) {
     for (int j=0; j<m; j++) {
        res.at(i*m+j) = matA.at(i*m+j)+matB.at(i*m+j);
     }
  }
  return 0;
}

int multiply_tvec_vec(std::vector<double> &tvec, std::vector<double> vec, std::vector<double> res, int n) {  // matrice colonne fois une matrice ligne
  int i,j;
  for (i=0; i<n; i++) {
     for (j=0; j<n; j++) {
         res.at(n*i+j) = tvec.at(i)*vec.at(j);
     }
  }
      
  return 0;
}

int multiply_matrix_vec(std::vector<double> &mat, std::vector<double> vec, std::vector<double> res, int n) {
  int i,j;
  double sum;
  for (i=0; i<n; i++) {
     sum=0.0;
     for (j=0; j<n; j++) {
         sum=sum+mat.at(n*i+j)*vec.at(j);
     }
     res.at(i)=sum;
  }
      
  return 0;
}


int multiply_tvec_matrix(std::vector<double> &mat, std::vector<double> vec, std::vector<double> res, int n) {
  int i,j;
  double sum;
  for (i=0; i<n; i++) {
     sum=0.0;
     for (j=0; j<n; j++) {
         sum=sum+mat.at(n*j+i)*vec.at(j);
     }
     res.at(i)=sum;
  }
      
  return 0;
}

int multiply_matrix_matrix(std::vector<double> &matA, std::vector<double> &matB, std::vector<double> &matC, int n) {
  int i,j,k;
  double sum;
  for (i=0; i<n; i++) {
     for (j=0; j<n; j++) {
	sum=0.0;
	for (k=0; k<n; k++) {
            sum=sum+matA.at(n*i+k)*matB.at(n*k+j);
	}
     }
     matC.at(n*i+j)=sum;
  }
      
  return 0;
}

int multiply_scal_matrix(std::vector<double> &mat, double scal, int n) {
    int i,j;
    for (i=0; i<n; i++) {
       for (j=0; j<n; j++) {
          mat.at(n*i+j)=mat.at(n*i+j)*scal;     
       }
    }
    return 0;
}

int multiply_scalprod(std::vector<double> &row, std::vector<double> &vec, double &res, int n) {
  int i;
  res=0.0;
  for (i=0; i<n; i++) {
	  res=res+row.at(i)*vec.at(i);
  }
  return 0;
}

int multiply_crossprod(std::vector<double> vec, std::vector<double> row, std::vector<double> mat, int n) {
  int i,j;
  for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
          mat.at(n*i+j)=vec.at(i)*row.at(j);
      }
  }
  return 0;
}


/***********************************************************
 * mother class IDENT05_ABSLSQ
 * ********************************************************/
 IDENT05_ABSLSQ :: IDENT05_ABSLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dnp, double fTs, double fvarian ) :
	  y(std::vector<double> (dsize, 0.0)),
	  u(std::vector<double> (dsize, 0.0)),
	  size(dsize),
	  np(dnp),
	  Ts(fTs),
	  varian(fvarian)
{
 // nothing here
  int k;
  for (k=0; k<dsize; k++) {
     y[k]=ydata[k];
     u[k]=udata[k];
  }
} //end constructor

//copy constructor
 IDENT05_ABSLSQ:: IDENT05_ABSLSQ(const IDENT05_ABSLSQ &source) 
{
}  
// afffectation constructor
 IDENT05_ABSLSQ &IDENT05_ABSLSQ :: operator=(const IDENT05_ABSLSQ &source) 
{
}

//destructor
 IDENT05_ABSLSQ:: ~IDENT05_ABSLSQ(void) 
{
}     

/***********************************************************
 * derived class IDENT05_BASICLSQ
 *********************************************************/
 IDENT05_BASICLSQ :: IDENT05_BASICLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dna, double fTs, double fvarian) : IDENT05_ABSLSQ(ydata, udata, dsize, dna, fTs, fvarian)
// IDENT05_BASICLSQ::IDENT05_BASICLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dna, double fTs, double fvarian ) :
//       y(std::vector<double> (dsize, 0.0)), // constructor this is important
//       u(std::vector<double> (dsize, 0.0)),
//       size(dsize),
//       na(dna),
//       Ts(fTs),
//       varian(fvarian)
{
 // nothing here
  na=dna;
  int k;
  for (k=0; k<dsize; k++) {
     y[k]=ydata[k];
     u[k]=udata[k];
  }
} //end constructor

//copy constructor
//
 IDENT05_BASICLSQ :: IDENT05_BASICLSQ(const IDENT05_BASICLSQ &source) : IDENT05_ABSLSQ(source)
{
}  
// afffectation constructor
 IDENT05_BASICLSQ &IDENT05_BASICLSQ :: operator=(const IDENT05_BASICLSQ &source) //: IDENT05_ABSLSQ(source) 
{
}

//destructor
 IDENT05_BASICLSQ:: ~IDENT05_BASICLSQ(void) 
{
}   

bool IDENT05_BASICLSQ::innovation(double &epsilon) 
{
#ifdef _RELSQ_TEST_FIXED_MATDIM_
// calculation of the iteration with the procedures above for 3x3 matrices

#else
// calculation of the iteration with the boost lib (general)


#endif
	return 0;
}
 
bool IDENT05_BASICLSQ::predict(int k, double yiter, double &epsilon) 
{
#ifdef _RELSQ_TEST_FIXED_MATDIM_
// calculation of the iteration with the procedures above for 3x3 matrices

#else
// calculation of the iteration with the boost lib (general)


#endif
	return 0;
}

bool IDENT05_BASICLSQ::update()
{
#ifdef _RELSQ_TEST_FIXED_MATDIM_
// calculation of the iteration with the procedures above for 3x3 matrices

#else
// calculation of the iteration with the boost lib (general)


#endif
	return 0;
}

bool IDENT05_BASICLSQ::getParams() 
{
   std::cout<< " Result from statlsq (AR type) params:\n";
   std::cout<< "\t autoregressive part:\n";
   for (int j=0; j<(na); j++) {
       std::cout<< "a[" << j << "]="<< theta[j] <<", ";
   }
   return 0;
}

/******************************************************
 * derived-derived class IDENT05_STATLSQ
 *****************************************************/
/*
 IDENT05_STATLSQ::IDENT05_STATLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dna, int dnb, double fTs, double fvarian ) : IDENT05_ABSLSQ( ydata, udata, dsize, dna, fTs, fvarian ) // first constructor as derived class: must be present even if the args (np vs na,nb) does not correspond to the derived class
{
  nb=dnb;
  iter=0;
 // nothing here
  int k;
  for (k=0; k<dsize; k++) {
     y[k]=ydata[k];
     u[k]=udata[k];
  }
} //end constructor

//destructor
 IDENT05_STATLSQ:: ~IDENT05_STATLSQ() 
{
}  

bool IDENT05_STATLSQ::innovation(double &epsilon) 
{
   int i;
#ifdef _RELSQ_TEST_FIXED_MATDIM_
// calculation of the iteration with the procedures above for 3x3 matrices
   for (i=0; i<(na+nb); i++) {
       theta.at(i) = theta.at(i) + matrixK.at(i)*epsilon;
   }
#else
// calculation of the iteration with the boost lib (general)

#endif
	return 0;
}
 
bool IDENT05_STATLSQ::predict(double yiter, double &epsilon) 
{
   int flag;
   double ypred;
#ifdef _RELSQ_TEST_FIXED_MATDIM_
// calculation of the iteration with the procedures above for 3x3 matrices
   flag = multiply_scalprod(phi, theta, ypred, na+nb);
#else
// calculation of the iteration with the boost lib (general)
//   epsilon = inner_prod (phi, theta);
#endif
   epsilon = yiter - ypred;
   yhat.push_back(ypred);

   return flag;
}

bool IDENT05_STATLSQ::update()
{
   int flag;
   std::vector<double> Fphi, Kphi, F_add;
   double res, invS;
      
#ifdef _RELSQ_TEST_FIXED_MATDIM_
// calculation of the iteration with the procedures above for 3x3 matrices
   flag = multiply_matrix_vec(matrixF, phi, Fphi, na+nb);
   flag = multiply_scalprod(phi, Fphi, res, na+nb);
   invS = 1.0/(lamb1/lamb2 + res);
   flag = copy_matrix(Fphi, matrixK, na+nb, na+nb);
   flag = multiply_scal_matrix(matrixK, invS, na+nb);  // should work with copy container as a scalar
   flag = multiply_tvec_vec(matrixK, phi, Kphi, na+nb);
   flag = multiply_matrix_matrix(Kphi, matrixF, F_add, na+nb);
   flag = sum_matrices(matrixF, F_add, matrixF, na+nb, na+nb);
   flag = multiply_scal_matrix(matrixF, 1.0/lamb1, na+nb);
#else
// calculation of the iteration with the boost lib (general)


#endif
     return flag;
}

bool IDENT05_STATLSQ::pass_iodata(std::list<dpair> &list_yh, std::string str_data) 
{
   dpair singleton;
   for (int j=0; j<size; j++) {
       singleton.x = j*Ts;
       singleton.y = yhat.at(j);
       list_yh.push_back(singleton);
   }
   return 0;
}

bool IDENT05_STATLSQ::getParams() 
{
   std::cout<< " Result from statlsq (AR type control) params:\n";
   std::cout<< "\t autoregressive part:\n";
   for (int j=0; j<(na); j++) {
       std::cout<< "a[" << j << "]="<< theta[j] <<", ";
   }
   std::cout<< "\n \t and control part:\n";
   for (int j=na; j<(na+nb); j++) {
       std::cout<< "b[" << j << "]="<< theta[j] <<", ";
   }
   return 0;
}

*/
