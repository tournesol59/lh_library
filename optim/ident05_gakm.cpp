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
 * s.t. discrete equality from ode:
 *      dx/dt=g(x,u,t)		discretized in: xk+1=g_(x1..xk,u1..uk)
 *      y=h(x,t) + v(t)		discretized in: yk=h(xk)+vk
 * AND: at each ipopt iteration
 *      x~k+1=Ak*x~k+Bk*udata,k +EKF*(ydata,k-Ck*x~k)
 *      Ak,Bk,Ck are from a linearization routine
 *      Extended Kalman Filter furnish EKF matrix and S
*************************************************/
#include "ident05_gakm.hpp"

#include "ident05_nlp.hpp"
#include <iostream>

using namespace Ipopt;

/** attemps to access optimization parameters via macros */
#define P1_EPS (ident05_nlp::get_p1_eps)
#define P2_OM (ident05_nlp::get_p2_om)

// constructor
IDENT05_GAKM::IDENT05_GAKM(fstream& filedata, int order);
{
   p1_eps_copy=P1_EPS;
   p2_om_copy= P2_OM;

   l0=malloc(7*sizeof(float));
   l1=malloc(7*sizeof(float));
   l2=malloc(7*sizeof(float));
   l3=malloc(7*sizeof(float));
   l4=malloc(7*sizeof(float));
   l5=malloc(7*sizeof(float));
   l6=malloc(7*sizeof(float));

   l0[0]=1.0; l0[1]=0.0; l0[2]=0.0; l0[3]=0.0; l0[4]=0.0; l0[5]=0.0; l0[6]=0.0;
   l1[0]=0.0; l1[1]=1;   l1[2]=0.0; l1[3]=0.0; l1[4]=0.0; l1[5]=0.0; l1[6]=0.0;
   l2[0]=-0.5;l2[1]=0;   l2[2]=1.5; l2[3]=0.0; l2[4]=0.0; l2[5]=0.0; l2[6]=0.0;
   l3[0]=0.0; l3[1]=-1.5;l3[2]=0.0;l3[3]=2.25;l3[4]=0.0; l3[5]=0.0; l3[6]=0.0;
}

//destructor
IDENT05_GAKM::~IDENT05_GAKM()
{}

