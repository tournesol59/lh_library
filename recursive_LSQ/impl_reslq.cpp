#include "decl_reslq.h"

void IDENT05_RELSQ::IDENT05_RELSQ(int dsize, int dn, double fTs ) 
{
	size=dsize;
	n=dn;
	Ts=fTs;
}

void IDENT05_RELSQ::~IDENT05_RELSQ() 
{
}

bool IDENT05_RESLQ::recursive_algorithm(int k, double &y[], double &u[],
		             double varian,
                             double y_h_1, double y_hd_1,
                             double &y_h[], double &y_hd[], int order,
			     double &pe_y[], double &pe_yd[]) {

  double K1k=2*(2*k-1)/k/(k+1);
  double K2k=6/k/(k+1)/Ts;

  double Xres = y[k] - y_h_1-y_hd_1*Ts;
  double y_h_p=y_h_1+y_hd_1*Ts+K1k*Xres;
  double y_hd_p=y_hd_1+K2k*Xres;

  y_h[k]=y_h_p;
  y_hd[k]=y_hd_p;

  Pe_y[k]=2*(2*k-1)/k/(k+1)*varian*varian;
  Pe_yd[k]=12/k/(k*k-1)/Ts*varian*varian;

  return 1;
}

