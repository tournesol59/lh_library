#include "decl_relsq.hpp"

#ifndef _RELSQ_TEST_FIXED_MATDIM_
#define _RELSQ_TEST_FIXED_MATDIM_
#endif


/*
 * method of class IDENT05_ABSGAUSSLSQ: abstract class !
 */

//  IDENT05_ABSGAUSSLSQ::IDENT05_ABSGAUSSLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dna, int dnb, int dnc, double fTs, double fvarian)
//{
//}

//  IDENT05_ABSGAUSSLSQ::~IDENT05_ABSGAUSSLSQ() 
//{ 
//}

/*
 * method of daughter class IDENT05_MAGAUSSLSQ
 */
 //(Comments to separate the future code from the actual)

   // constructor
    IDENT05_MAGAUSSLSQ :: IDENT05_MAGAUSSLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dnc, double fTs, double fvarian) : IDENT05_ABSGAUSSLSQ(ydata, udata, size, dnc, fTs, fvarian)
{
     // nothing here
  int k;
  for (k=0; k<dsize; k++) {
     y[k]=ydata[k];
     u[k]=udata[k];
     wk[k]=0.0;
     zk[k]=0.0;
  }
} //end constructor

    IDENT05_MAGAUSSLSQ :: ~IDENT05_MAGAUSSLSQ()
{
}

// exact copy of the method of BASICLSQ (see "impl_basiclsq.cpp")
bool  IDENT05_MAGAUSSLSQ::pass_iodata(std::vector<double> &dataoutvec)
{
    int row=0; 
    // TBC
// skeleton
  //  auto it=dataoutvec.begin();
    while (row<size)  {
        dataoutvec.push_back(yhat[row]); // warning about what is really calculte in yhat
        row++;
//	it++;
   }//end while
  return 0;
}

bool IDENT05_MAGAUSSLSQ::getParams() 
{
   std::cout<< " Result from gausslsq (MA type control) params:\n";
   std::cout<< "\t moving average part:\n";
   for (int j=0; j<(nc); j++) {
       std::cout<< "c[" << j << "]="<< theta[j] <<", ";
   }
   return 0;
}


   // method for calling a method from another class? we will make a copy if no mean found


/*
 * method of daughter class IDENT05_ARMAGENLSQ
 */
   // constructor
    IDENT05_ARMAGENLSQ :: IDENT05_ARMAGENLSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dna, int dnb, int dnc, double fTs, double fvarian) : IDENT05_ABSGAUSSLSQ(ydata, udata, size, dnc, fTs, fvarian)
{
  na = dna;
  nb = dnb;
     // nothing here
  int k;
  for (k=0; k<dsize; k++) {
     y[k]=ydata[k];
     u[k]=udata[k];
     wk[k]=0.0;
     zk[k]=0.0;
  }
} //end constructor

    IDENT05_ARMAGENLSQ :: ~IDENT05_ARMAGENLSQ()
{
}




