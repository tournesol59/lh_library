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
IDENT05_GAKM::IDENT05_GAKM(fstream& filedata, int order)
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
   l3[0]=0.0; l3[1]=-1.5;l3[2]=0.0;l3[3]=2.5;l3[4]=0.0; l3[5]=0.0; l3[6]=0.0;
   l4[0]=0.375;l4[1]=0.0;l4[2]=-3.75;l4[3]=0.0;l4[4]=4.375;l4[5]=0.0;l4[6]=0.0;
   l5[0]=0.0; l5[1]=1.875; l5[2]=0.0; l5[3]=-8.75; l5[4]=0.0; l5[5]=7.875; l5[6]=0.0;   
   l6[0]=-0.3125; l6[1]=0.0; l6[2]=6.5625; l6[3]=0.0; l6[4]=-19.6875; l6[5]=0.0; l6[6]=14.4375;
}

IDENT05_GAKM::eval_leg_pol(int order, float t)   // order goes up to 6
{
   switch case order {
	   0: {
		   return l0[0];
	      };
   	   1: {
		   return l1[0]+t*l1[1];
	      };
	   2: {
		   return l2[0]+t*l2[1]+t*t*l2[2];
	      };
		   return l3[0]+t*l3[1]+t*t*l3[2]+t*t*t*l3[3];
	      };
	   4: {
		   return l4[0]+t*l4[1]+t*t*l4[2]+t*t*t*l4[3]+t*t*t*t*l4[4];
	      };
	   5: {
		   return l5[0]+t*l5[1]+t*t*l5[2]+t*t*t*l5[3]+t*t*t*t*l5[4]+t*t*t*t*t*l5[5];
	      };
	   6: {
		   return l6[0]+t*l6[1]+t*t*l6[2]+t*t*t*l6[3]+t*t*t*t*l6[4]+t*t*t*t*t*l6[5]+t*t*t*t*t*t*l6[6];
	      };

   }
}

//destructor
IDENT05_GAKM::~IDENT05_GAKM()
{}
  /** Methods to solve the VanderPol problem and only in that scope
   * shall be called in eval_g() Ipopt method: */
IDENT05_GAKM::get_gakm_invec_trans method(int order, int ma, int N,
		                  numbers* u_t,
				  numbers* u_phi)
{ //uphi has dims 1x(ma*order)
// algorithm of trapezoidal rule to compute <u(t),leg[k](t)>=c_phi(k)
  float* ti, tii, dt, dt_1, a, b;
  numbers*  sum;
  sum=malloc(order*sizeof(numbers));
  dt=tsample;
  n=floor(tquantiz/tsample);

  for (k=0;k<=order;k++) {
	for (j=0;j<ma;j++) {
		a=j*tquantiz;
		b=(j+1)*tquantiz;
		dt_1=dt*(b-a)/2;
		sum[k]=0.0;

		for (i=0;i<n;i++) {
			ti=-1+float(i)/float(n)*2;
			tii=-1+float(i+1)/float(n)*2;
			sum[k]+=dt_1*0.5*((u_t[j*n+i]*eval_leg_pol(k,ti))+
					u_t[j*n+i+1]*eval_leg_pol(k,tii));
		}
		switch case k {
			0: u_phi[j*order+k]=sum[0]*2/5;

			1: u_phi[j*order+k]=sum[1]*2/5;			   

			2: u_phi[j*order+k]=sum[2]*2/5;			   

			3: u_phi[j*order+k]=sum[3]*2/5;			   

			4: u_phi[j*order+k]=sum[4]*2/5;			   

			5: u_phi[j*order+k]=sum[5]*2/5;			   

			6: u_phi[j*order+k]=sum[6]*2/5;			   
		}
	}//end range
  }//end order

   return 0;   
}
IDENT05_GAKM::set_gakm_matvec_VdP(int order, int ma, int,
		              numbers* mat_phi, numbers* vec_phi) 
{
       	//specific to VanderPol Problem
	//
	//
};

IDENT05_GAKM::eval_gakm_coeff_trans_xi(int order, int ma
		              number* evalm_leg,
			      numbers* time_vec) 
{
   // VanderMonde Matrix of evaluation coeff of legendre polynomials
};

