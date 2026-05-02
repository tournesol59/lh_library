#include "tseries.h"
#include <vector>
#include <cmath>

int main() {

   int Nt=20;
   int Na=2; int Nb=1; int Nc=2;
   std::vector<double> data; // size Nt input data: outputs
   std::vector<double> cldata; // size Nt input data: controls
   std::vector<double> inCa; // size Na, regression coefficients
   std::vector<double> inCb;
   std::vector<double> inCc;

   for (int i=0; i<Nt; i++) {
      int c = i*i % 7;
      int d = i*i % 12;
      int e = i*i % 3;
      double x = (double) d*pow(-1,e);
      double y = (double) c*pow(-1,e);
      data.push_back(x); // random
      cldata.push_back(y); // random      
   }
   inCa.push_back(1.0);
   for (int i=1; i<Na; i++) {
      inCa.push_back(0.0);
   }
   inCb.push_back(1.0);
   for (int i=1; i<Nb ; i++) {
      inCb.push_back(0.0);
   }
   inCc.push_back(1.0);   
   for (int i=1; i<Nc ; i++) {
      inCc.push_back(0.0);
   }
   BASIC_MODEL iBasicmodel = BASIC_MODEL( data,  inCa, 0.1, 1.0, Nt, Na) ;
   AR_MODEL iARmodel = AR_MODEL( data, inCa, 0.1, 1.0, Nt, Na) ;
   ARX_MODEL iARXmodel = ARX_MODEL( data, cldata, inCa, inCb, 0.1, 1.0, Nt, Na-1, Nb);
   // test base class methods of class BASIC:
  // iBasicmodel.predict(3);
   iBasicmodel.update(3);
   iBasicmodel.export_ts();
   std::cout << std::endl;

   // test derived class AR methods are called (cout << "this is B method..."
  // iARmodel.predict(3);
   iARmodel.update(3);
   iARmodel.export_ts();
   std::cout << std::endl;
   
   // test derived class ARX methods are called (cout << "this is C method..."
   iARXmodel.update(0); // indep. from k
   iARXmodel.export_ts();
   std::cout << std::endl;
//   delete iBasicmodel;
//  delete iARmodel;
//  delete iARXmodel;

   return 0;
}
