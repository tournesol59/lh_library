#ifndef __TSERIES_H__
#define __TSERIES_H__
#endif
/*
 * BASIC_MODEL |=|  IDENT05_ABSLSQ // IDENT05_BASICLSQ
 * AR_MODEL    |=|  IDENT05_STATLSQ
 * ARX_MODEL   |=|  NO EQUIVALENT
*/

#include <iostream>
#include <vector>
// A
class BASIC_MODEL {
   protected:

      double Ts;
      double varian;
      int nt, na;  //nt is the number of time samples, na is the order A(q1)x(t)=eps(t)
      std::vector<double> ydata;
      std::vector<double> coeffa;
      std::vector<double> tspred;

   public:

      BASIC_MODEL(std::vector<double> data,
		  std::vector<double> incoeffa,
		  double fTs, double fvarian,
		  int dnt, int dna) ;
      ~BASIC_MODEL();
      // virtual
     virtual int predict(int k); // for a simple sample
     virtual int update(int k); // think of a linear regression for autoregressive model. This is not the case for arma models
     int export_ts();

};

// B
class AR_MODEL : public BASIC_MODEL
{
    protected:

    public:
      AR_MODEL(std::vector<double> data,
	       std::vector<double> incoeffa,
	       double fTs, double fvarian,
	       int dnt,int dna) ;
      ~AR_MODEL();
     int predict(int k); 
     int update(int k);      
};

// C
class ARX_MODEL : public AR_MODEL
{
   protected:
      int nb;  // na, nb are the order A(q1)x(t)=B(q1)u(t)+eps(t)
      std::vector<double> coeffb ;
      std::vector<double> udata;
   public:
      ARX_MODEL(std::vector<double> data,
	        std::vector<double> cldata,
		std::vector<double> incoeffa,
		std::vector<double> incoeffb,
		double fTs,double fvarian,
		int dnt,int dna,int dnb);
      ~ARX_MODEL();
     int update(int k);      
     int export_ts();
};

// D 
/*
class ARMAX_MODEL : public AR_MODEL
{
   protected:
     int nb;  // na, nb are the order A(q1)x(t)=B(q1)u(t)+eps(t)
     vector<double> coeffb ;

     int nc;  // na, nb are the order A(q1)x(t)=B(q1)u(t)+eps(t)
     vector<double> coeffc ;

   public:
      ARMAX_MODEL(Ts,data,incoeffa,incoeffb,incoeffc,varian,nt,na,nb,nc);

};

*/
