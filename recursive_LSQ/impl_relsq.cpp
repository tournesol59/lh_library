#include <vector>
#include "decl_relsq.hpp"
#include <math.h>

using namespace lhlib;

IDENT05_RELSQ::IDENT05_RELSQ(std::vector<double> ydata, std::vector<double> udata, int dsize, int dn, double fTs, double fvarian ) :
       y(std::vector<double> (dsize, 0.0)), // constructor this is important
	u(std::vector<double> (dsize, 0.0)),
	y_h(std::vector<double> (dsize, 0.0)),
	y_hd(std::vector<double> (dsize, 0.0)),
	y_hdd(std::vector<double> (dsize, 0.0)),
	y_f(std::vector<double> (dsize, 0.0)),
	y_fd_con(std::vector<double> (dsize, 0.0)),	
	pe_y(std::vector<double> (dsize, 0.0)),
	pe_yd(std::vector<double> (dsize, 0.0)),
	pe_ydd(std::vector<double> (dsize, 0.0)),
        myvar(std::vector<double> (dsize, 0.0)),
       size(dsize),
       n(dn),
       Ts(fTs),
       varian(fvarian)
{
 // nothing here
  int k;
  for (k=0; k<dsize; k++) {
     y[k]=ydata[k];
     u[k]=udata[k];
  }
}//end constructor

IDENT05_RELSQ::~IDENT05_RELSQ() 
{
}

/*
bool IDENT05_RELSQ::test_assign() 
{
   auto it=y_h.begin();
   std::cout << *it;
   y_h[2]=3.1415;

   return 1;
}
*/
bool IDENT05_RELSQ::pass_iodata(std::vector<std::pair<double,double>> &list_yh, std::string str_data)
{
  int k, sel=1;
  std::pair<double,double> singleton;

  if (!(str_data.compare("yori"))) {
	  std::cout << "match yori\n";
	  sel=1;
  }
  else if (!(str_data.compare("uori"))) {
	  sel=2;
  }
  else if (!(str_data.compare("yst0"))) {
	  sel=3;
  }
  else if (!(str_data.compare("yst1"))) {
	  std::cout << "yst1 match\n";	  
	  sel=4;
  }
  else if (!(str_data.compare("yst2"))) {
          sel=5;
  }
  else if (!(str_data.compare("pys0"))) {
	  sel=6;
  }
  else if (!(str_data.compare("pys1"))) {
	  sel=7;
  }
  else if (!(str_data.compare("pys2"))) {
	  sel=8;
  }
  else if (!(str_data.compare("myvar"))) {
	  sel=9;
  }
  for (k=0; k<size; k++) {
      singleton.first = (double (k))*Ts;
      switch (sel) {
         case 1:  singleton.second = y[k];
		  break;
	 case 2:  singleton.second = u[k];
		  break;
         case 3:  singleton.second = y_h[k];
		  break;
	 case 4:  singleton.second = y_hd[k];
		  break;
	 case 5:  singleton.second = y_hdd[k];
		  break;
	 case 6:  singleton.second = pe_y[k];
		  break;
	 case 7:  singleton.second = pe_yd[k];
		  break;
	 case 8:  singleton.second = pe_ydd[k];
 		  break;
         case 9:  singleton.second = myvar[k];
		  break;
      }
      std::cout << singleton.second << " ";
      list_yh.push_back( singleton );
  }

  return 1;
}

bool IDENT05_RELSQ::apply_derivative_filter(int k, double poleT1, int manage) {
  
  if (k==0) {
     y_f[k]=y[k];
  }
  else {
     y_f[k]=y_f[k-1]*poleT1+y[k]*(1-poleT1);
     y_fd_con[k] = (y_f[k]-y_f[k-1])/Ts;
  }
  // manage signal to operate with follwing recursive_algorithm
  if (manage==1) {  // if =1, erase the original data with the small filtered one
     y[k]=y_f[k];
  }
  else if (manage==2) {
     y[k]=y_fd_con[k];
  }
  return 0;
}


bool IDENT05_RELSQ::recursive_algorithm(int k, int order) {
  double y_h_1 =0.0;
  double y_hd_1 =0.0;
  double y_hdd_1 =0.0;
  double y_h_p =0.0;
  double y_hd_p =0.0;
  double y_hdd_p =0.0;
  double K1k =0.0;
  double K2k =0.0;
  double K3k =0.0;
  double Xres =0.0;

        if (k<1) {
         std::cout << "Warning divide by zero!\n";
         return 1;
	}
        if (k==1) {
         y_h_1 = 1.0;
         y_hd_1 = 0.0;
         y_hdd_1 = 0.0;
	}
	else {
         y_h_1 = y_h[k-1];
         y_hd_1 = y_hd[k-1];
         y_hdd_1 = y_hdd[k-1];
	}

	if (order == 1)    // case first order truncation, lsq tracks a line
	{

         K1k=((double) 2*(2*k-1))/((double) k*(k+1));
         K2k=((double) (6))/((double) k*(k+1)*Ts);
         Xres = y[k] - y_h_1 - y_hd_1*Ts;
         y_h_p = y_h_1 + y_hd_1*Ts + K1k*Xres;
         y_hd_p = y_hd_1 + K2k*Xres;
     
	 y_h[k]=y_h_p;
 	 y_hd[k]=y_hd_p;
         y_hdd[k]=0.0;
	 
         pe_y[k]=((double) 2*(2*k-1))/((double) (k*(k+1)))*varian*varian;
	 
         if (k>1) {
		 pe_yd[k]=12.0/((double) k*(k*k-1)*Ts)*varian*varian;
	 }
	 else {
		 pe_yd[k]=12.0/Ts*varian*varian;
	 }
         pe_ydd[k]=0.0;
	 // user variable to troubleshoot quality of rand generator:
	 myvar[k]=sqrt(k*pe_y[k]);

	} // endif
	if (order == 2)    // case second order truncation, lsq tracks a 2nd O. Curve
	{
        K1k=((double) (3*(3*k*k-3*k+2)))/((double) (k*(k+1)*(k+2)));
        K2k= ((double) (18*(2*k-1)))/((double) (k*(k+1)*(k+2)*Ts));
        K3k= ((double) (60/k/(k+1)))/((double) ((k+2)*Ts*Ts));

        Xres = y[k] - y_h_1 - y_hd_1*Ts - 0.5*y_hdd_1*Ts*Ts;
        y_h_p = y_h_1 + y_hd_1*Ts + 0.5*y_hdd_1*Ts*Ts + K1k*Xres;
        y_hd_p = y_hd_1 + y_hdd_1*Ts*Ts + K2k*Xres;
        y_hdd_p = y_hdd_1 + K3k*Xres;
	
        y_h[k]=y_h_p;
        y_hd[k]=y_hd_p;
        y_hdd[k]=y_hdd_p;

        pe_y[k]=((double) 3*(3*k*k-3*k+2))/((double) (k*(k+1)*(k+2)*varian*varian));
	if (k>1) {
          pe_yd[k]=((double) 12*(16*k*k-30*k+11))/((double) (k*(k*k-1)*(k*k-4)*Ts*Ts))*varian*varian;
          pe_ydd[k]=((double) 720)/((double) (k*(k*k-1)*(k*k-4)*Ts*Ts*Ts*Ts))*varian*varian;
	}
	else {
	  pe_yd[k]=12.0/Ts/Ts*varian*varian;
          pe_ydd[k]=720.0/Ts/Ts/Ts/Ts*varian*varian;
	}

	}//endif order
//       std::cout << y_h[k] << " " << y_hd[k] << " " << y_hdd[k] << " "
//	       << pe_y[k] << " " << pe_yd[k] << " " << pe_ydd[k] << "\n";
  return 1;
}

/******************************************
 * class IDENT05_TFcont_OP
 * ***************************************/
IDENT05_TFcont_OP::IDENT05_TFcont_OP(int dsize, int dna, int dnb, double dTf, double dTs, std::vector<double> dcoeffs) :
	yinit(std::vector<double> ((dna+dnb), 0.0)), // constructor
	yout(std::vector<double> (dsize*(dna+dnb), 0.0))
	
{
	na=dna;
	nb=dnb;
	maxsize=(dna+dnb)*dsize;
	Tf=dTf;
	Ts=dTs;
	for (int i=0; i< (dna+dnb); i++) {
           Gab[i]=dcoeffs[i];
	}
	for (int i=0; i< dsize; i++) {
           for (int j=0; j<(dna+dnb); j++) {
	      yout[i*(dna+dnb)+j]=0.0;  // if it there is no other common mean..
	   }
	}
  // Nothing here
}//end constructor

IDENT05_TFcont_OP::~IDENT05_TFcont_OP()
{
}

bool IDENT05_TFcont_OP:: advanceintime(int k) 
{      
	return 1;
}	// knowing actual time and time step, steps k time integration forwards
	  
//double IDENT05_TFcont_OP::printoutput(int k, int order j) {}

IDENT05_TFdisc_OP::IDENT05_TFdisc_OP(int dsize, int dnc, int dnd, double dTs, std::vector<double> dcoeffs) :
	yinit(std::vector<double> ((dnc+dnd), 0.0)), // constructor
	yout(std::vector<double> (dsize*(dnc+dnd), 0.0))
	
{
	nc=dnc;
	nd=dnd;
	maxsize=(dnc+dnd)*dsize;
	Ts=dTs;
	for (int i=0; i< (dnc+dnd); i++) {
           Hcd[i]=dcoeffs[i];
	}
	for (int i=0; i< dsize; i++) {
           for (int j=0; j<(dnc+dnd); j++) {
	      yout[i*(dnc+dnd)+j]=0.0;  // if it there is no other common mean..
	   }
	}
  // Nothing here
}//end constructor

IDENT05_TFdisc_OP::~IDENT05_TFdisc_OP()
{
}

/******************************************
 * inherited class IDENT05_REEST
 * ***************************************/

IDENT05_REEST::IDENT05_REEST(std::vector<double> ydata, std::vector<double> udata, int dsize, int dn1, int dn2, int dr, double fTs, double fvarian_var_phi, double fvarian_var_y) :
       y(std::vector<double> (dsize, 0.0)), // constructor this is important
	u(std::vector<double> (dsize, 0.0)),
	y_h(std::vector<double> (dsize, 0.0)),
	y_hd(std::vector<double> (dsize, 0.0)),
	y_hdd(std::vector<double> (dsize, 0.0)),
	pe_y(std::vector<double> (dsize, 0.0)),
	pe_yd(std::vector<double> (dsize, 0.0)),
	pe_ydd(std::vector<double> (dsize, 0.0)), 
       theta(std::vector<double> ((dn1+dn2)*dsize, 0.0)),
       phi(std::vector<double> (dsize*2, 0.0)),
       y_fc(std::vector<double> (dsize, 0.0)),
       y_fd(std::vector<double> (dsize, 0.0)),
       u_fc(std::vector<double> (dsize, 0.0)),
       u_fd(std::vector<double> (dsize, 0.0)),
       psi(std::vector<double> (dsize*(n1+n2), 0.0)),
       Rcov_var_phi(std::vector<double> (dsize*1, 0.0)),
       Rcov_var_y(std::vector<double> (dsize*1, 0.0)),
       L_EIV(std::vector<double> (dsize*4*4, 0.0)),
       g(std::vector<double> (dsize*(4+1), 0.0)),
       LA(std::vector<double> (dsize*4, 0.0)),
       size(dsize),
       n1(dn1),
       n2(dn2),
       r(dr),
       Ts(fTs),
       varian_var_phi(fvarian_var_phi),
       varian_var_y(fvarian_var_y)
    {
 // nothing here
  int k;
  for (k=0; k<dsize; k++) {
     y[k]=ydata[k];
     u[k]=udata[k];
  }
}//end constructor                          double fvarian_var_y);

IDENT05_REEST::~IDENT05_REEST() 
{
}
bool  IDENT05_REEST::block_algorithm(int k)
{
	return 1;
}

bool  IDENT05_REEST::firststep_getvals(int k)
{
	return 1;
}

bool  IDENT05_REEST::twostep_filtercvals(int k)
{
	return 1;
}

bool  IDENT05_REEST::thirdstep_getdisc()
{
	return 1;
}

bool  IDENT05_REEST::fourthstep_filterdvals(int k)
{
	return 1;
}

bool  IDENT05_REEST::fifthstep_lseinstr()
{
	return 1;
}

bool IDENT05_REEST::pass_iodata(std::vector<std::pair<double,double>> &list_yh, std::string str_data)
{
  int k, sel=1;
  std::pair<double, double> singleton;

  if (!(str_data.compare("yori"))) {
	  std::cout << "match yori\n";
	  sel=1;
  }
  for (k=0; k<size; k++) {
      singleton.first = (double (k))*Ts;
      switch (sel) {
         case 1:  singleton.second = y[k];  // more case when complete
		  break;
      }
      std::cout << singleton.second << " ";
      list_yh.push_back( singleton );
  }

  return 1;
}      
