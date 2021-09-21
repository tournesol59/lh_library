/*
 * Class IDENT05_RELSQ implements actually time-vector of data 
 * for being identified parametrically. For the moment only non regressive
 * linear (a*t+b) are implemented
 */
#include <iostream>
#include <cstring>
#include <vector>
#include <list>
#include <iterator>
#include "../include/lhTypes.hpp"


//using namespace std;
using namespace lhlib;

int generate_rdom_example(std::vector<double> &sig, int n, double ca, double cb, double fvar);

class IDENT05_RELSQ {

  public:
    IDENT05_RELSQ(std::vector<double> ydata, std::vector<double> udata,
		    int dsize, int dn, double fTs, double fvarian);
   ~IDENT05_RELSQ();

    bool recursive_algorithm(int k, int order);
    bool test_assign();
    bool pass_iodata(std::list<dpair> &list_yh, std::string str_data);

    std::vector<double> y;  // known size at init
    std::vector<double> u;
    std::vector<double> y_h; // size equal to length(y) or u
    std::vector<double> y_hd;
    std::vector<double> y_hdd;
    std::vector<double> pe_y;
    std::vector<double> pe_yd;
    std::vector<double> pe_ydd;
    std::vector<double> myvar;
  private:
    int size;  // size of the original data
    int n;  // order for future recursive filter algorithms
    double Ts; // sample time
    double varian; // estimation of the variance (y,original - y,mean (1st or 2nd order))

};

  class IDENT05_TFcont_OP { // operative transfer function is continuous time form

  public:
	  IDENT05_TFcont_OP(int dsize, int dna, int dnb, double dTf, double dTs, std::vector<double> dcoeffs);
	  ~IDENT05_TFcont_OP();
	  int stepactual;
	  std::vector<double> yinit;  // init value for y and derivative must be filled by constructor
	  std::vector<double> yout;  // output value ! caution: y and na-1 derivatives follows each other and then other values for time step
          // void setintegcolloc(IDENT05_COLLOC cnfirst);
	  bool advanceintime(int k);  // knowing actual time and time step, steps k time integration forwards
	  //double printoutput(int k, int order j);

  private:
          double Tf; // computing step time size
	  double Ts; // sampling time
	  int maxsize;
	  int na,nb;
	  std::vector<double> Gab;             // model in transfer func (na numerator values first, nb denom follows the na values)

  };

  class IDENT05_TFdisc_OP {
  
  public:
	  IDENT05_TFdisc_OP(int dsize, int dnc, int dnd, double dTs, std::vector<double> dcoeffs);
	  ~IDENT05_TFdisc_OP();
	  int stepactual;
	  std::vector<double> yinit;  // init value for y and derivative must be filled by constructor
	  std::vector<double> yout;  // output value ! caution: y and na-1 derivatives follows each other and then other values for time step

	  bool advanceintime(int k);  // knowing actual time and time step, steps k time integration forwards

  private:
	  double Ts; // sampling time
	  int maxsize;
	  int nc,nd;

	  std::vector<double> Hcd;             // ARMA filter for white noise <- (y-yh) (nc discrete numerator first, nd denom follows the nc values)
  };


  class IDENT05_REEST { // : public IDENT05_RELSQ{
//*************** CURRENT FORM ************************
 

//******************* TO USE LATER: *******************
// l'algorithme iteratif decrit au tout début de la thèse soit avec / l'amélioration
// proposée utilisant un vecteur indpt de l'erreur construit avec des variables
// instumentales
//
// KE = L_t_1 * psi_t*(LA+ psi_t' L_t_1 psi_t)^(-1)
// psi_t = [[Rcov_var,phi_t_1*var, phi_t]
//          [ y_t ]]
// g(t) = [[ var'Rcov_var,y_t_1 ]
//         [ y(t) ]]
// LA = [[-var*var, id],
//       [id,       0]]
// estim_t=estim_t_1+KE(g(t) - psi_t'*estim_t_1)
//
// Rcov_var,phi_t=Rcov_var,phi_t_1 + var*phi_t'
// Rcov,var,y_t=Rcov,var,y_t_1+var*y'
// L_t=L_t_1 - KE*psi*L_t_1
// var(t)=[-y(t-r) ... y(t-n1-r) u(t) ... u(t-n2)]
// phi_t=[y(t) y(t-1) y(t-2)]  // discrete delayed outputs (subject to noise)

//*************************************************

  public:
    IDENT05_REEST(std::vector<double> ydata, std::vector<double> udata, int dsize, int dn1, int dn2, int dr, double fTs, double fvarian_var_phi, double fvarian_var_y);
    ~IDENT05_REEST();

    std::vector<double> y;      // known size at init
    std::vector<double> u;      // input
    std::vector<double> y_h;    // estim, size equal to length(y) or u
    std::vector<double> y_hd;   // estim of deriv
    std::vector<double> y_hdd;  // estim of second deriv
    std::vector<double> pe_y;   // covariance of estim(y)-y
    std::vector<double> pe_yd;  //
    std::vector<double> pe_ydd; //

    std::vector<double> theta;           // the params vector,size n1+n2 prev. named "estim"
    std::vector<double> phi;             // the components for ident (delayed outputs) vector
    std::vector<double> y_fc;            // by fc filtered values
    std::vector<double> y_fd;            // by fd filtered values
    std::vector<double> u_fc;            // by fc filtered input
    std::vector<double> u_fd;            // by fd filtered input
    std::vector<double> psi;             // the instrumental variable vector (step 2 and 4: filtered output and inputs)
    std::vector<double> Rcov_var_phi;
    std::vector<double> Rcov_var_y;

     // for future recursive version:
    std::vector<double> L_EIV;
    std::vector<double> g;               
    std::vector<double> LA;

    bool block_algorithm(int k);
    bool firststep_getvals(int k);
    bool twostep_filtercvals(int k);
    bool thirdstep_getdisc();
    bool fourthstep_filterdvals(int k);
    bool fifthstep_lseinstr();
    bool pass_iodata(std::list<dpair> &list_yh, std::string str_data);

  private:
    int size;
    int n1,n2,r; // order in numerator of TF, in denominator of TF, r=discret delay
    double Ts; // sampling time
    double varian_var_phi, varian_var_y; // known variances Cov(X,Y)=Sum (Xi-Exi)*(Yi-Eyi)

//    IDENT05_TFcont_OP link
    std::vector<double> y_gab_first;  // output value auxiliary init model (from simple LSE)
//    IDENT05_TFcont_OP link
    std::vector<double> y_ac_two;     // output value continuous filter used in filtering second step of the algorithm
//    IDENT05_TFdisc_OP link
    std::vector<double> y_hcd_third; // output value discrete filter of noise used in filtering third step of the algorithm

  };

