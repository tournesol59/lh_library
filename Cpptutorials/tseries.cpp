#include "tseries.h"

#include <iostream>

// A: base class w/o algorithm, only purpose is data initialization
BASIC_MODEL:: BASIC_MODEL(std::vector<double> data,std::vector<double> incoeffa,double fTs,double fvarian,int dnt,int dna) :
  Ts(fTs),   // order must follows declaration order
  varian(fvarian),
  nt(dnt),
  na(dna),
  ydata(std::vector<double> (dnt,0.0)),
  coeffa(std::vector<double> (dna,0.0)),	
  tspred(std::vector<double> (dnt,0.0))
{
   for (int i=0; i<nt; i++) { 
      ydata[i]=data[i]; // fill imported values 
   }

   for (int i=0; i<na; i++) { 
      coeffa[i]=incoeffa[i]; // fills initial values of coeffa 
   }
}

BASIC_MODEL:: ~BASIC_MODEL()
{
}

int BASIC_MODEL:: predict(int k)
{
// empty function, compute mean for test

   long double sum=0.;

   for (int i=0; i<nt; i++) {
      sum=sum+ydata[i] ;
   }
   tspred[k]= sum/nt ;
   std::cout << "this is our base class predict method: " << tspred[k] << std::endl;
   return 0;
}

int BASIC_MODEL:: update(   int k  )
{
// again empty function, only to test the inheritance
   int stat = predict(k);  // call the previous function, here only to experiment the inheritence recipies
 if (stat==0) {  
   double meas = ydata[k];
   double pred= tspred[k];

   for (int i=0; i<na ;i++) { // looks a bit simple but no outer loop k shall be in this function
      coeffa[i]=coeffa[i]+1.0/(k+1)*(meas-pred) ;
   }
   std::cout << "this is our base class update method: " << coeffa[0] << std::endl;
 }
   return stat;
}

int BASIC_MODEL:: export_ts(  )
{
// again empty function, only
    std::cout << "the data values " << std::endl;

    for (int k=0; k<nt; k++) {
      std::cout << ydata[k] << "  ";
    }
    std::cout << std::endl;
    std::cout << "the predicted values " << std::endl;

    for (int k=0; k<nt; k++) {
      std::cout << tspred[k] << "  ";
    }
    std::cout << "This is our base class export method! " << std::endl; 
    return 0;
}

// --------- B --------------
// 
AR_MODEL:: AR_MODEL(std::vector<double> data,std::vector<double> incoeffa,double fTs,double fvarian,int dnt,int dna) : BASIC_MODEL(data,incoeffa,fTs,fvarian,dnt,dna) 
{
 /*
 * for (int i=0; i<nt; i++) { 
      ydata[i]=(data[i]); // fills initial values with import
   }

   for (int i=0; i<na ; i++) { 
      coeffa[i]=(incoeffa[i]); // fills initial values of coeffa 
   }
 */
}

AR_MODEL:: ~AR_MODEL()
{
}

int AR_MODEL:: predict(int k) 
{
// here modification for derived class: compute ponderated mean for test

 long double sum=0.;
 //  // call to test the inheritance of a virtual method
 int stat = -1;
 if (stat==0) {
   std::cout << "this is our derived (AR) class predict method with call to base class, value is: " << tspred[k] << std::endl;
 }
 else {
  for (int i=0; i<na; i++) {
      double x;
      if ((k-na+i) < 0) { x=0.0; }
      else { x=ydata[k-na+i]; }

      sum=sum+coeffa[i]*x ;
   }
   tspred[k]= sum ;
   std::cout << "this is our derived (AR) class predict method: " << tspred[k] << std::endl;
 }
   return 0;
}

int AR_MODEL::update( int k )
{
  int flag=BASIC_MODEL::update(k);
  
  std::cout << "this is our derived (AR) class predict method: " << tspred[k] << std::endl;

  return flag;
}

/*
int AR_MODEL:: update( int k ) 
{
// test the inheritance with linear reg (valable only if na=1)
   int stat = predict(k);
 if (stat==0) {
   double meas = ydata[k];
   double pred= tspred[k];

   for (int i=0; i<na ;i++) { // looks a bit simple but no outer loop k shall be in this function
      coeffa[i]=coeffa[i]+1.0/(k+1)*(meas-pred) ;
   }
   std::cout << "this is our derived (AR) class update method: " << coeffa[0] << std::endl;
 }
   return stat;
}
*/

// -------- C ---------
  ARX_MODEL:: ARX_MODEL(std::vector<double> data,std::vector<double> cldata,std::vector<double> incoeffa,std::vector<double> incoeffb,double fTs,double fvarian,int dnt,int dna,int dnb) : AR_MODEL(data,incoeffa,fTs,fvarian,dnt,dna),
  coeffb(std::vector<double> (dnb,0.0)),
  udata(std::vector<double> (dnt,0.0))
{
   nb=dnb;
   /*
   udata.alloc();
   coeffb.alloc();
   */
   for (int i=0; i<nt ; i++) { 
      udata[i] = (cldata[i]); // fills control values 
   }
   
   for (int i=0; i<nb ; i++) { 
      coeffb[i] = (incoeffb[i]); // fills initial values of coeffa 
   }

}

ARX_MODEL:: ~ARX_MODEL()
{
}

int ARX_MODEL:: update(   int k  ) 
{
// test the inheritance with linear reg (valable only if na=1)

   double corr_xx=0.;
   double corr_yx=0.;
   
   if (na==1) {
      for (int i=0; i<(nt-1) ;i++) { // works only for na==1
         corr_xx += ydata[i]*ydata[i];
	 corr_yx += ydata[i+1]*ydata[i];
      }
      coeffa[0]= corr_yx / corr_xx; // corr_xx >0 guaranted if nt>=1
   }
   else { std::cout << "ARX_MODEL::update not currently work for na>1" << std::endl;
   }
   std::cout << "this is our derived (ARX) class update method: " << coeffa[0] << std::endl;
   return 0;
}

int ARX_MODEL:: export_ts( )
{
// again empty function, only
    std::cout << "the udata values " << std::endl;

    for (int k=0; k<nt; k++) {
      std::cout << udata[k] << "  ";
    }
    std::cout << std::endl;
    std::cout << "the predicted values " << std::endl;

    for (int k=0; k<nt; k++) {
      std::cout << tspred[k] << "  ";
    }
    std::cout << "This is our derived (ARX) class export method! " << std::endl; 
    return 0;
}




