/*
 * Program that imports (t,u,y) data and perform a model estimation of params
 * y- a1*y-1 + ... + an*y-n = b1*u + ... + bm*u-m-1 + e
 * n=na, m=nb
 * F Peugny LTS
 *
 * still an error at compilation unit impl_basicrelsq.hpp
 * decl_basiclsq.hpp:56:8: error: candidates are: IDENT05_BASICLSQ::IDENT05_BASICLSQ(const IDENT05_BASICLSQ&)
 */
#include "decl_basiclsq.hpp"
#include "../include/ident05_fftdata.hpp"
//#include "../recursive_LSQ/decl_relsq.hpp"; // generate_from_file:
int generate_from_file(std::vector<double> &sig, int &Npty, const char * filename, int n, int m);

int main(int argc, char **argc) {

   char inFileName[14];
   int Npty;
   std::vector<double> u_indata;
   std::vector<double> y_indata;
   
   //import data
   strcpy(inFileName, (const char*) argv[2]);
   Npty=500;  //argv[3]; prefere mettre en dur pour secondorder.dat
	   // assume fsample=0.1 (not imported)
   generate_from_file(u_indata, Npty, inFileName, 3, 2);
   generate_from_file(y_indata, Npty, inFileName, 3, 3); //y is-third sigin file
   
   // instanciates ident_..lsq class, assumes varianz=0.5
   IDENT05_BASICLSQ inst_basiclsq(y_indata, u_indata, Npty, 2, 0.1, 0.5 );

   double epsilon=0.;
   // loop over samples (recursive) 
   // for (int i=3; i <Npty; i++) { // 3 because 0,1,2 are used for init
     i=3;
        inst_basiclsq.predict(i, y_indata[i], epsilon);
        inst_basiclsq.innovation(epsilon);
        inst_basiclsq.update();
   // }

   // program terminates correctly, vector ressources free-ed automaatically
   return 0;
}
