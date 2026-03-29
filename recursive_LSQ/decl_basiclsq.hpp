/*
 * Class IDENT05_BASICLSQ implements actually time-vector of data 
 * for being identified parametrically. For the moment only declaration are
 * exposed for inheritance mechanism
*/
#include <iostream>
#include <cstring>
#include <vector>
#include <list>
#include <iterator>
#include "../include/lhTypes.hpp"


//using namespace std;
using namespace lhlib;

///////////////////
// Prototypes
int copy_matrix(std::vector<double> matA, std::vector<double> &matcpy, int n, int m) ;
int inverse_nthree_matrix(std::vector<double> &mat, std::vector<double> &invmat) ;
int sum_matrices(std::vector<double> matA, std::vector<double> matB, std::vector<double> res, int n, int m) ;
int multiply_tvec_vec(std::vector<double> &tvec, std::vector<double> vec, std::vector<double> res, int n) ;
int multiply_matrix_vec(std::vector<double> &mat, std::vector<double> vec, std::vector<double> res, int n) ;
int multiply_matrix_matrix(std::vector<double> &matA, std::vector<double> &matB, std::vector<double> &matC, int n) ;
}
int multiply_scal_matrix(std::vector<double> &mat, double scal, int n) ;
int multiply_scalprod(std::vector<double> &row, std::vector<double> &vec, double &res, int n) ;
int multiply_crossprod(std::vector<double> vec, std::vector<double> row, std::vector<double> mat, int n) ;

////////////////
// Class Definition

   class IDENT05_ABSLSQ {
         public:
    IDENT05_ABSLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dnp, double fTs, double fvarian );
    ~IDENT05_ABSLSQ();

	 private:
       std::vector<double> y;
       std::vector<double> u;
       Index size;
       Index np;
       Number Ts;
       Number varian;

} // end class IDENT05_ABSLSQ


/////////////////
// Inheritance for BASICLSQ from ABSLSQ
 class IDENT05_BASICLSQ : IDENT05_ABSLSQ {
         public:
    IDENT05_BASICLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dna, double fTs, double fvarian );
    ~IDENT05_BASICLSQ();

	 private:
       std::vector<double> y;
       std::vector<double> u;
       Index size;
       Index na;
       Number Ts;
       Number varian;

} // end class IDENT05_BASICLSQ

