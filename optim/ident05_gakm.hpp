/************************************************
* LIEBHERR AEROSPACE TLSE SAS
* Function declaration to implement the Galerkin method
* for ODE 2nd order resolution
* Context is the Ipopt ident05 problem:
*
* Problem ident05 looks like this:
 * defined a set parameter, states and slack variables:
 * p1,p2,p3,p4,p5,p6,p7,p8.., x0..xN, v1..vN
 *
 * min sum v1/q*v1 + v2/q*v2 .. + vN/q*vN + log(Sk)
 * s.t. discrete equality from ode:
 *      dx/dt=g(x,u,t)		discretized in: xk+1=g_(x1..xk,u1..uk)
 *      y=h(x,t) + v(t)		discretized in: yk=h(xk)+vk
 * AND: at each ipopt iteration
 *      x~k+1=Ak*x~k+Bk*udata,k +EKF*(ydata,k-Ck*x~k)
 *      Ak,Bk,Ck are from a linearization routine
 *      Extended Kalman Filter furnish EKF matrix and S
*************************************************/


#ifndef __IDENT05_GAKM_HPP__
#define __IDENT05_GAKM_HPP__

//#include "IpTNLP.hpp"

class IDENT05_GAKM
{
public:
  /** default constructor */
  IDENT05_GAKM(fstream& filedata, int order);

  /** default destructor */
  virtual ~IDENT05_GAKM();
   
  /** Methods to calculate the Legendre Polynoms */
  float eval_leg_pol(int order, float t);  // order goes up to 6
  /** Method to prepare the calculus of integrals and interpolate the input u */
  bool set_gakm_data();

  /** Methods to solve the VanderPol problem and only in that scope
   * shall be called in eval_g() Ipopt method: */
  bool get_gakm_invec_trans method(int order, int ma, int N,
		                  numbers* u_t,
				  numbers* u_phi); //uphi has dims 1x(ma*order)
  bool set_gakm_matvec_VdP(int order, int ma, int,
		              numbers* mat_phi, numbers* vec_phi); //specific to VanderPol Problem
  bool eval_gakm_coeff_trans_xi(int order, int ma
		              number* evalm_leg,
			      numbers* time_vec);

//  bool prepare_gakm_pb();
//  bool inverse_gakm_mat();
//
  /** Methods to display output solution in state vector form + coeff at each time quantiz unit */
  bool print_gakm_out(std::ftream fileout);


private:
  /** Data that should be returned by methods only */
  /** first a copy of the parameters from the NLP program */
  float p1_eps_copy;
  float p2_om_copy;

  /** second the 6order Legendre polynoms should be init */
  float * l0, l1, l2, l3, l4, l5, l6;

  /** an array for each state of the VanderPol model */
  float final_time=2.0; // from 0 to 2 sec
  int vec_st_len=200; // for 2sec, also named N in comments
  float tsample=0.01;
  float * states1;
  float * states2;
  float * input1;

  /** an array for the 6th order approximation of x(t) and dx(t)/dt */
  /**  first a quantized time that is a multiple of tsample */
  float tquantiz=0.1;
  int num_ranges = 20; // also named ma in comments
  int order = 6;
  int vec_phi_len = 6*20; // order*ma gives the length of the Legendre coeff vector
  /** the integrated solution coefficients, corresponding to a Legendre approximation */
  float * integPol1;
  float * integPol2;

  /** The upper triangular matrix*/
  int mat_g_len; //format column index first
  float * integ_mat;

    /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
  //@{
  //  HS071_NLP();
  IDENT05_GAKM(const IDENT05_GAKM&);
  IDENT05_GAKM& operator=(const IDENT05_GAKM&);
  //@}
#endif


