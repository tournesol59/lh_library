/*
 * Program that imports (t,u,y) data and perform a model estimation of params
 *  y- a1*y-1 + ... + an*y-n = b1*u + ... + bm*u-m-1 + e
 * this time with the stage method from J. Sabiti (1994) An efficient estimating method for ARMA Process
 * F Peugny LTS
*/
#include <math.h>
#include "decl_iogenerate.hpp"
#include "decl_basiclsq.hpp"
#include "decl_optimlsq.hpp"
#include "../include/ident05_fftdata.hpp"
//#include "../recursive_LSQ/decl_relsq.hpp";

//immport func generate_from_file:
int generate_more_from_file_second(std::vector<double> &sig, int dsize, const char* fileName, int m) {

int main(int argc, char **argv) {

   char inFile[20];
   int Npty;
   std::vector<double> u_indata;
   std::vector<double> y_indata;

   strcpy(inFile, (const char*) argv[1]);
   Npty=std::atoi(argv[2]);
   generate_more_from_file_second(u_indata, int dsize, inFile, int 1);
   return 0;
}
