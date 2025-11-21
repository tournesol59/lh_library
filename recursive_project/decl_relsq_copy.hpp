/*
 * Class IDENT05_RELSQ implements actually time-vector of data 
 * for being identified parametrically. For the moment only non regressive
 * linear (a*t+b) are implemented
 *
 * However, second part of this file are tentative of project
 * declaration about the experimental variable method (not continued)
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
int generate_arN_example(std::vector<double> &sig, int n, int order, double x0, std::vector<double> &ar, double fvar) ;
int generate_from_file(std::vector<double> &sig, int &n, const char * filename) ;

class IDENT05_RELSQ {

  public:
    IDENT05_RELSQ(std::vector<double> ydata, std::vector<double> udata,
		    int dsize, int dn, double fTs, double fvarian);
   ~IDENT05_RELSQ();

    bool apply_derivative_filter(int k, double poleT1, int manage); // this must always be called firstly if you want to operate on a filtered signal
    bool recursive_algorithm(int k, int order);
    
    bool pass_iodata(std::vector<std::pair<double,double>> &list_yh, std::string str_data);

    std::vector<double> y;  // known size at init
    std::vector<double> u;
    std::vector<double> y_h; // size equal to length(y) or u
    std::vector<double> y_hd; // identified derivative of y_h signal
    std::vector<double> y_hdd; // identified second derivative of y_h signal
    std::vector<double> pe_y; // evolution of the covariance of (y-y_h)
    std::vector<double> pe_yd; // ... covariance of derivative
    std::vector<double> pe_ydd; // .. covariance of second derivative
    // note that the last covariance should decrease to zero
    std::vector<double> myvar;
    std::vector<double> y_f; // output from T1 filter
    std::vector<double> y_fd_con; // output from derivate
  protected:
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




