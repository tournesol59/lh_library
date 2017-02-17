/************************************************
* LIEBHERR AEROSPACE TLSE SAS
* Function implementation to implement the Galerkin method
* for ODE 2nd order resolution
* Context is the Ipopt ident05 problem:
*
* Problem ident05 looks like this:
 * defined a set parameter, states and slack variables:
 * p1,p2,p3,p4,p5,p6,p7,p8.., x0..xN, v1..vN
 *
 * min sum v1/q*v1 + v2/q*v2 .. + vN/q*vN + log(Sk)
 * s.t. discrete equality from ode OR galerkin/collocation method:
 *      dx/dt=g(x,u,t)		discretized in: g(x_phi_j)=d(u_phi_j)
 *      			in galerkin method: xk=x(tk)=sum,j cj*phi_j(tk)
 *      y=h(x,t) + v(t)		discretized in: yk=h(xk)+vk
 * AND: at each ipopt iteration
 *      x~k+1=Ak*x~k+Bk*udata,k +EKF*(ydata,k-Ck*x~k)
 *      Ak,Bk,Ck are from a linearization routine
 *      Extended Kalman Filter furnish EKF matrix and S
*************************************************//


#ifndef __IDENT05_GAKM_HPP__
#define __IDENT05_GAKM_HPP__

#include "../../Ipopt-3.12.6/include/coin/IpTypes.hpp"  // define Number as double
#include "./MathFunctions.h"
#include <fstream>
#include <iostream>
#include <math>
#include <malloc>  //<map>?

class IDENT05_GAKM
{
public:
  /** default constructor */
  IDENT05_GAKM();

  /** default destructor */
  virtual ~IDENT05_GAKM();
   
  /** Methods to calculate the Legendre Polynoms */
  bool eval_leg_pol(int order, int len, Number* tvec, Number* out);  // order goes up to 6
  /** Method to prepare the calculus of integrals and interpolate the input u */
  bool set_gakm_data();

  /** Methods to solve the VanderPol problem and only in that scope
   * shall be called in eval_g() Ipopt method: */
  bool get_gakm_invec_trans method(int order, int ma, int N, Number* u_t,
				  Number* u_phi); //uphi has dims 1x(ma*order)
  bool get_gakm_iterative_avg_and_corr(int order, int ma, int N,
		                  Number* xi, Number* zi, 
				  number* avgX2, number* corrXz);

  bool set_gakm_matvec_VdP(int order, int ma, int N,
		              Number* mat_phi, Number* vec_phi); //specific to VanderPol Problem
  bool eval_gakm_coeff_trans_xi(int order, int ma
		              number* evalm_leg,
			      Number* time_vec);

//  bool prepare_gakm_pb();
//  bool inverse_gakm_mat();
//
  /** Methods to display output solution in state vector form + coeff at each time quantiz unit */
  bool print_gakm_out(std::ifstream fileout);


private:
  /** Data that should be returned by methods only */
  /** first a copy of the parameters from the NLP program */
  float p1_eps_copy;
  float p2_om_copy;

  /** second the 6order Legendre polynoms should be init */
  float * l0, l1, l2, l3, l4, l5, l6;

  /** an array for each state of the VanderPol model */
  float final_time; // from 0 to 2 sec
  int vec_st_len; // for 2sec, also named N in comments
  float tsample;

  float * states1;
  float * states2;
  float * input1;
  float * output1;
  float * parameters;  //  eps and wO


  /** an array for the 6th order approximation of x(t) and dx(t)/dt */
  /**  first a quantized time that is a multiple of tsample */
  float tquantiz;
  int num_ranges; // also named ma in comments
  int order = 6;
  int vec_phi_len; // order*ma gives the length of the Legendre coeff vector

  /** the fermeture equation coefficients to approx VdP eqn with linear eqn */
  float vdp_alpha;
  float vdp_gamma;

  /**  coefficients a Legendre approximation */
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


