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
*************************************************/
#include "ident05_gakm.hpp"
#include <cassert>
#include <iostream>

#ifndef __TEST_GAKM_ONLY__
#define __TEST_GAKM_ONLY__
#endif

using namespace std;
using namespace lhlib;

/** attemps to access optimization parameters via macros */

//#define P1_EPS (ident05_nlp::get_p1_eps)
//#define P2_OM (ident05_nlp::get_p2_om)

// constructor
IDENT05_GAKM::IDENT05_GAKM()
{
#ifndef __TEST_GAKM_ONLY__
   p1_eps_copy=P1_EPS;
   p2_om_copy= P2_OM;
#else
   p1_eps_copy=0.01;
   p2_om_copy=4;
#endif

  l0=malloc(7*sizeof(Number));
  l1=malloc(7*sizeof(Number));
  l2=malloc(7*sizeof(Number));
  l3=malloc(7*sizeof(Number));
  l4=malloc(7*sizeof(Number));
  l5=malloc(7*sizeof(Number));
  l6=malloc(7*sizeof(Number));

   l0[0]=1.0; l0[1]=0.0; l0[2]=0.0; l0[3]=0.0; l0[4]=0.0; l0[5]=0.0; l0[6]=0.0;
   l1[0]=0.0; l1[1]=1.0; l1[2]=0.0; l1[3]=0.0; l1[4]=0.0; l1[5]=0.0; l1[6]=0.0;
   l2[0]=-0.5;l2[1]=0.0; l2[2]=1.5; l2[3]=0.0; l2[4]=0.0; l2[5]=0.0; l2[6]=0.0;
   l3[0]=0.0; l3[1]=-1.5;l3[2]=0.0; l3[3]=2.5; l3[4]=0.0; l3[5]=0.0; l3[6]=0.0;
   l4[0]=0.375;l4[1]=0.0;l4[2]=-3.75;l4[3]=0.0;l4[4]=4.375;l4[5]=0.0;l4[6]=0.0;
   l5[0]=0.0; l5[1]=1.875;l5[2]=0.0;l5[3]=-8.75;l5[4]=0.0;l5[5]=7.875;l5[6]=0.0;
   l6[0]=-0.3125;l6[1]=0.0;l6[2]=6.5625;l6[3]=0.0;l6[4]=-19.6875;l6[5]=0.0;l6[6]=14.4375;

   final_time=2.0;
   vec_st_len=200;
   tsample=0.01;

   tquantiz=0.1;
   num_ranges = 20; // also named ma in comments   
   order = 6;
   vec_phi_len =120; // order*ma gives the length of the Legendre coeff vector

   vdp_alpha = -0.3;
   vdp_gamma = 0.01;
}

IDENT05_GAKM::eval_leg_pol(int order, int len, Number* tvec, Number* out)   // order goes up to 6
{
   int i;
   double result;
   for (i=0; i<=len; i++) {
       switch case (order) {
	  case 0:	   result= l0[0];
	     		   break;
   	  case 1: 	   result= l1[0]+pow(tvec[i],1)*l1[1];
	      		   break;
	  case 2: 	   result= l2[0]+pow(tvec[i],1)*l2[1]+pow(tvec[i],2)*l2[2];
	      		   break;
	  case 3:      result= l3[0]+pow(tvec[i],1)*l3[1]+pow(tvec[i],2)*l3[2]+pow(tvec[i],3)*l3[3];
	     	           break;
	  case 4:	   result= l4[0]+pow(tvec[i],1)*l4[1]+pow(tvec[i],2)*l4[2]+pow(tvec[i],3)*l4[3]+pow(tvec[i],4)*l4[4];
	     	           break;
	  case 5: 	   result= l5[0]+pow(tvec[i],1)*l5[1]+pow(tvec[i],2)*l5[2]+pow(tvec[i],3)*l5[3]+pow(tvec[i],4)*l5[4]+pow(tvec[i],5)*l5[5];
	      		   break;
	  case 6:          result= l6[0]+pow(tvec[i],1)*l6[1]+pow(tvec[i],2)*l6[2]+pow(tvec[i],3)*l6[3]+pow(tvec[i],4)*l6[4]+pow(tvec[i],5)*l6[5]+pow(tvec[i],6)*l6[6];
			   break;
       }
       out[i]=result;
   }
}

//destructor
IDENT05_GAKM::~IDENT05_GAKM()
{}
  /** Methods to solve the VanderPol problem and only in that scope
   * shall be called in eval_g() Ipopt method: */
IDENT05_GAKM::get_gakm_invec_trans_method(int order, int ma, int N,
		                  Number* u_t,
				  Number* u_phi)
{ //uphi has dims 1x(ma*order)
// algorithm to compute <u(t),leg[k](t)>=c_phi(k)
  int i,j,k,n;
  Number* t, dt, dt_1, a, b;
  Number* t_vec,u_vec, phi_vec; // time-vectors time and legendre evaluation
  Number*  sum;

  dt=tsample;
  n=N/ma;  // recalculated

  sum=malloc(order*sizeof(Number));
  t_vec=malloc(n*sizeof(Number));
  u_vec=malloc(n*sizeof(Number));  
  phi_vec=malloc(n*sizeof(Number));


  for (j=0;j<ma;j++) {
      a=j*tquantiz;
      b=(j+1)*tquantiz;

      sum=0.0;
     
      for (i=0; i<n; i++) {
	  t=-a+i*dt;  // goes from a (i=0) to b (i=n)
	  t_vec[i]= (t-(a+b)/2)*2/(b-a);  // spanned to [-1,1]
	  u_vec[i]=u_t[j*n+i];  // copy the range (a,b) from u_t
      }

      for (k=0;k<=order;k++) {
	      // compute phi_vec, evaluation of Legendre,k at t_vec times in [-1,1]
	  eval_leg_pol(k, n, &tvec, &phi_vec);
	  sum = myComputeScalarProduct(n, tvec, &u_vec, &phi_vec);
          switch(k) {
			case 0: u_phi[j*order+k]=sum*0.5;
				break;   
			// normed function int(l0(x)^2.dx, x=-1..1)
			case 1: u_phi[j*order+k]=sum*1.5;   
				break;
			case 2: u_phi[j*order+k]=sum*2.5;
				break;
			case 3: u_phi[j*order+k]=sum*3.5;		   
				break;
			case 4: u_phi[j*order+k]=sum*4.5;
				break;
			case 5: u_phi[j*order+k]=sum[5]*5.5;
				break;
			case 6: u_phi[j*order+k]=sum[6]*6.5;
				break;
	}
      }//end order
  }//end range

   return 0;   
}

IDENT05_GAKM::get_gakm_iterative_avg_and_corr(int order, int ma, int N,
		                  Number* xi, Number* zi, 
				  Number* avgX2, Number* corrXz)
{

  Number eps, ohm0;
  int i,j,k;
  double sumX2, sumXz;
  //------------------- explanation:-------------------
  // set the average and correlation Number needed to linearize the VdP equation system
  // IPOPT take an array of Number simply called 'x' (no matter xi, zi in our functions) 
  //    x(1:p)           :	set of Parameters: p=2 for VdP problem
  //    x(p+1:p+2*N)     :	set of states (xi), 2*N because dim of VdP states =2
  //    x(p+2*N+1:p+4*N) :	set of complement states (zi)
  //    x(p+4*N+1:p+5*N) :	set of slack variable used for calc the cost function (the noise)
  //  and this is all. Note that zero indexing ius used in the implementation
  // hence the values xi, zi that are inputs of this function are part of the 'x' array

  sumX2=0.0;
  sumXz=0.0;
  j=0;
  k=0;
  for (i=0;i<N;i++) {
	sumX2+= (xi[i])^2;
	sumXz+= (xi[i])*(zi[k]);
	j++; // then dxi/dt[]
	j++; // then xi[+1]
	k++; // then dzi/dt[]
	k++; // then zi[+1]
  }
  avgX2=double(sumX2)/double(N);  // conv to Number can also be used
  corrXz=double(sumXz)/double(N);

  if (j<2*N) || (k<2*N) {
	  return 0; }
  else {
	  return -1; } // detect segment. error

}

IDENT05_GAKM::set_gakm_matvec_VdP(int order, int ma, int N,
		              Number* mat_phi, Number* vec_phi) 
{
   //specific to VanderPol Problem
   // equation with linear solution: d2X/dt2 - wo*X(t) = 0 (+u)
   // equation with approx. deviation: 
   //    d2z/dt2 - eps.wo* (1 - 2*avg{X*z}-avg{X^2})*dz/dt - wo^2*z(t) = eps.wo*(1-avg{X^2})dX/dt
   //    this should be augmented by following iteration depdt. equation (for the moment): 
   //    d(avg{X*z})/dk=Alpha*avg{X*z} + Gamma**avg{X^2}, i.e. Alpha=-0.3 <0, Gamma=0.01 >0



};

IDENT05_GAKM::eval_gakm_coeff_trans_xi(int order, int ma,
		              Number* evalm_leg,
			      Number* time_vec) 
{
   // VanderMonde Matrix of evaluation coeff of legendre polynomials

};

