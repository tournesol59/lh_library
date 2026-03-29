// NLP problem

#include "nlp_ctr01.hpp"
//extern "C" int mySolveLinearLapack(int n, int rhs, double * A, int lda, int * ipiv, double * B, int ldb, int info);
#include <cassert>

using namespace Ipopt;

/* Constructor. */
NLP_CTR01::NLP_CTR01()
{}

NLP_CTR01::~NLP_CTR01()
{}

bool NLP_CTR01::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style) {

  n=15;
   
  m=12;

  nnz_jac_g=2;

  nnz_h_lag=0;

  index_style=C_STYLE;
 
  return true;

}

bool NLP_CTR01:: get_bounds_info(Index n, Number* x_l, Number* x_u, 
                     Index m, Number* g_l, Number* g_u) {

  assert(n==15);

  assert(m==12);

 for(Index i=0; i<15; i++) {

    x_l[i]=-5.0;
  }
 for(Index i=0; i<11; i++) {
    g_l[i]=0.0;
    g_u[i]=1.0E19; // ineauqlity superior constraint 
  }
   g_l[11]=-1.0E19;
   g_u[11]=10.0;
  
  return true;
}

bool NLP_CTR01:: get_starting_point(Index n, bool init_x, Number* x, 
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
                                    Number* lambda)  {
  assert(init_x==true);
  assert(init_z==false);
  assert(init_lambda==false);
// a1, a2, a3, a4, k1, k2
 x[0]=1.0;
 x[1]=1.0;
 x[2]=1.0;
 x[3]=1.0;
 x[4]=1.0;
 x[5]=1.0;
//gamma hat 1,2,3,4; gamma bar 1,2,3,4
 x[6]=1.0;
 x[7]=0.0;
 x[8]=0.0;
 x[9]=1.0;
 x[10]=0.10;
 x[11]=0.0;
 x[12]=0.0;
 x[13]=0.10;
//ro:
 x[14]=0.001;
return true;
}

bool NLP_CTR01:: eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {

 assert(n==15);

 obj_value=x[7]+x[11];  //NO yes only " elements counts for the trace;

 return true;

}

bool NLP_CTR01:: eval_grad_f(Index n, const Number* x, bool new_x, Number * grad_f)  {

 assert(n==15);

 grad_f[0]=0.0;
 grad_f[1]=0.0;
 grad_f[2]=0.0;
 grad_f[3]=0.0;
 grad_f[4]=0.0;
 grad_f[5]=0.0;

 grad_f[6]=0.0;

 grad_f[7]=1.0;  //count in trace 
 grad_f[8]=0.0;
 grad_f[9]=0.0;
 grad_f[10]=0.0;
 grad_f[11]=1.0;  //count in trace 
 grad_f[12]=0.0;
 grad_f[13]=0.0;
 
 grad_f[14]=0.0;

 return true;

}


bool NLP_CTR01:: eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {

 assert(n==15);
 assert(m==11);
 
 //construct the matrix %Cons% (constraint) and determine its smallest proper value
 // to see if the %Cons% matrix is definite positive or not
 // the proper value are solved by a C to Fortran Call so it 
 int i,j,k;
 int n12,lda, info, work; 

 double Fi[4]={0.0, 1.0, -0.5, 1.0}; // ! 2x2
 double Gi[2]={-6.0, 1.0};                // 2x1
 double Hi[2]={-100.0, 10.0};	     // 1x2
 double Ci1[2]= {0.0, 1.0};	     // 2x1
 double Ci2[2]= {0.0, 0.0};	     // 2x1
 double Zi[1]={1.0};		     // 1x1
 double Ei[2]={0.0, 0.03};	     // 1x2
 double Qi[1]={1.0};		     // 1x1
 double Ri[1]={0.1};		     // 1x1

 double Ai[4];  // 2x2
 double Ki[2];  // 2x1

 double Fb[16]; // 4x4
 double Cb[4];  // 4x1
 double Eb[4];
 double Wb[4];  // 2x2
 double Gb[8];  // 4x2
 double Gam[16];// 4x4
 double Sig[16];// 4x4
 
 double Ba[100];//10x10
 double N[20];  //10x2
 double NtN[100];//10x10
 double Ma[10]; //1x10
 double Id[1];  //1x1

 double Cons[121]; //11x11
 double tCons[121]; //11x11
 double Diag[11];  // Eigenvalues of tCons
 //filling the matrix Fb:
 for (i=0; i<4; i++) {
    Ai[i]=x[i]; //a1,2,3,4
 }
 for (i=0; i<1; i++) {//yes only one iteration
   Ki[i]=x[i+4]; //k1
 }
 Ki[1]=0.0;
 for (i=2; i<3; i++) {
   Ki[i]=x[i+4]; //k2
 }
 Ki[3]=0.0;

 for (i=0; i<16; i++) {
    Fb[i]=0.0; //malloc->to zero
 }
 for (i=0; i<2; i++) {
    Fb[i]=Fi[i]; //f11,f12
 }
 for (i=4; i<6; i++) {
    Fb[i]=Fi[i-2]; //f21,f22
 }
 for (i=8; i<9; i++) {
    Fb[i]=Ki[0]; //k1
 } 
 for (i=10; i<12; i++) {   
   Fb[i]=Ai[i-10]-Ki[i-10]; //a11-k1, a12
 }
 for (i=12; i<13; i++) {   
   Fb[i]=Ki[1]; //0
 }
 for (i=14; i<16; i++) {   
   Fb[i]=Ai[i-12]-Ki[i-12]; //a21-k2, a22
 }
// filling the matrix Wb:
 Wb[0]=Qi[1];
 Wb[1]=0.0;
 Wb[2]=0.0;
 Wb[3]=Ri[1];

//filling the matrix Gb:
 for (i=0; i<8; i++) {
    Gb[i]=0.0;
 }
//Gi,Ki
 for (i=0; i<8; i++) {
    Gb[i]=0.0; //malloc->to zero
 }
 for (i=0; i<1; i++) {
    Gb[i]=Gi[i]; //g1
 }
 for (i=2; i<3; i++) {
    Gb[i]=Gi[i-1]; //g2
 }
 for (i=5; i<6; i++) {
    Gb[i]=Ki[i-5]; //k1
 }
 for (i=7; i<8; i++) {
    Gb[i]=Ki[i-5]; //k2
 }
// filling matrix Cb:
 for (i=0; i<4; i++) {
    Cb[i]=0.0;
 }
 Cb[0]=Ci1[0];
 Cb[1]=Ci1[1];
 Cb[2]=Ki[0]*Ci2[0]+Ki[1]*Ci2[1];  // Ki*Ci2
 Cb[3]=Ki[2]*Ci2[0]+Ki[3]*Ci2[1];
 // filling matrix Eb:
 for (i=0; i<4; i++) {
    Eb[i]=0.0;
 }
 Eb[0]=Ei[0];
 Eb[1]=Ei[1];
//filling the matrix Sig:
  for (i=0; i<8; i++) {
    Sig[i]=0.0;
 }
 Sig[0]=2.0;
 Sig[5]=2.0;
 Sig[10]=2.0;
 Sig[16]=2.0;

//filling the matrix Gam with the actual upper bound for covqriqnce matrix:
  for (i=0; i<16; i++) {
    Gam[i]=0.0; //malloc->to zero
 }
 Gam[0]=x[7]+x[11]; //gam,bar+gam,hat ,11
 Gam[5]=x[8]+x[12]; //" ", 12
 Gam[10]=x[9]+x[13]; //" ",21
 Gam[16]=x[10]+x[14]; // " ",22


//filling the matrix Cons then transpose it in tCons:
//filling the matrix Ba; NtN, Ma:
  for (i=0; i<100; i++) {
     Ba[i]=0.0; //malloc to zero
     NtN[i]=0.0;
  }
  for (i=0; i<20; i++) {
     N[i]=0.0; //malloc to zero
  }
  for (i=0; i<10; i++) {
     Ma[i]=0.0; //malloc to zero
  }
// Ba: this will be long:
// [Sig,   0,   Fb+Cb*Z*Eb]
// [0,    diag(Q,R)^-1,tGb]
// [t(Fb+Cb*Z*Eb),Gb,  Gam]
  for (i=0; i<4; i++) {
     Ba[i]=Sig[i];
  }
  for (i=10; i<14; i++) {
     Ba[i]=Sig[i-6];  // Sig[4] begins
  }
  for (i=20; i<24; i++) {
     Ba[i]=Sig[i-12];    // Sig[8] begins
  }
  for (i=30; i<34; i++) {
     Ba[i]=Sig[i-18];      // Sig[12] begins
  }
 // beginning of line 5 and 6 are zeros
 Ba[44]=1.0/Qi[0];

 Ba[55]=1.0/Ri[0];
 // first tiers of block line three: Fb+Cb*Zi*Ei
  for (i=60; i<64; i++) {
     Ba[i]=Fb[i-60] + Cb[0]*Zi[0]*Eb[i-60];
  }
  for (i=70; i<74; i++) {
     Ba[i]=Fb[i-66] + Cb[1]*Zi[0]*Eb[i-70];
  } 
  for (i=80; i<84; i++) {
     Ba[i]=Fb[i-72] + Cb[2]*Zi[0]*Eb[i-80];
  } 
  for (i=90; i<94; i++) {
     Ba[i]=Fb[i-78] + Cb[3]*Zi[0]*Eb[i-90];
  } 
  // second tiers of block line three: Gb
  for (i=64; i<66; i++) {
     Ba[i]=Gb[i-64];
  } 
  for (i=74; i<76; i++) {
     Ba[i]=Gb[i-72];
  }
  for (i=84; i<86; i++) {
     Ba[i]=Gb[i-80];
  } 
  for (i=94; i<96; i++) {
     Ba[i]=Gb[i-88];
  }
  // third tiers of block line three: Gam
   for (i=66; i<70; i++) {
     Ba[i]=Gam[i-66];
  } 
   for (i=76; i<80; i++) {
     Ba[i]=Gam[i-72];
  }
   for (i=86; i<90; i++) {
     Ba[i]=Gam[i-78];
  }
   for (i=96; i<100; i++) {
     Ba[i]=Gam[i-84];
  }

  // then the last two blocks (1,3) and (2,3) of matrix Ba are deduced from
  // the transpose of blocks (3,1) and (3,2) respectively 
  // because Ba is symmetric
  // so transpose block (3,1)
  for (i=6; i<10; i++) {
     Ba[i]=Ba[60+i*10];
  }
  for (i=16; i<20; i++) {
     Ba[i]=Ba[61+i*10];
  }
  for (i=26; i<30; i++) {
     Ba[i]=Ba[62+i*10];
  }
  for (i=36; i<40; i++) {
     Ba[i]=Ba[63+i*10];
  }
// transpose block (3,2)
  for (i=46; i<50; i++) {
     Ba[i]=Ba[64+i*10];
  }
  for (i=56; i<60; i++) {
     Ba[i]=Ba[65+i*10];
  }
// THE MATRIX BA IS COMPLETED (100x100)
// fill the matrix N and NtN:
  for (i=0; i<20; i++) {
     N[i]=0.0; //malloc to zero
  }
  N[0]=Ei[0];
  N[1]=Ei[1];
  for (i=0; i<100; i++) {
     NtN[i]=0.0; //malloc to zero
  }
  NtN[0]=N[0]*N[0];
  NtN[1]=N[0]*N[1];
  NtN[10]=N[1]*N[0];
  NtN[11]=N[1]*N[1];
// fill the matrix M:
  for (i=0; i<10; i++) {
     Ma[i]=0.0; //malloc to zero
  }
  Ma[6]=Cb[0];
  Ma[7]=Cb[1];
  Ma[8]=Cb[2];
  Ma[9]=Cb[3];
// NOW FILL THE MATRIX Cons
  for (i=0; i<121; i++) {
     Cons[i]=0.0; //malloc to zero
     tCons[i]=0.0;
  }
  for (i=0; i<10; i++) {
     Cons[i]=Ba[i]-x[14]*NtN[i];
  }
  Cons[10]=Ma[0];
  for (i=10; i<20; i++) {
     Cons[i+1]=Ba[i]-x[14]*NtN[i];
  }
  Cons[20+1]=Ma[1];
  for (i=20; i<30; i++) { 
     Cons[i+2]=Ba[i]-x[14]*NtN[i];
  }
  Cons[30+2]=Ma[2];
  for (i=30; i<40; i++) {
     Cons[i+3]=Ba[i]-x[14]*NtN[i];
  }
  Cons[40+3]=Ma[3];
  for (i=40; i<50; i++) {
     Cons[i+4]=Ba[i]-x[14]*NtN[i];
  }
  Cons[50+4]=Ma[4];  
  for (i=50; i<60; i++) {
     Cons[i+5]=Ba[i]-x[14]*NtN[i];
  }
  Cons[60+5]=Ma[5];
  for (i=60; i<70; i++) {
     Cons[i+6]=Ba[i]-x[14]*NtN[i];
  }
  Cons[70+6]=Ma[6];
  for (i=70; i<80; i++) {
     Cons[i+7]=Ba[i]-x[14]*NtN[i];
  }
  Cons[80+7]=Ma[7];
  for (i=80; i<90; i++) {
     Cons[i+8]=Ba[i]-x[14]*NtN[i];
  }
  Cons[90+8]=Ma[8];
  for (i=90; i<100; i++) {
     Cons[i+9]=Ba[i]-x[14]*NtN[i];
  }
  Cons[100+9]=Ma[9];

  for (i=110; i<120; i++) {
     Cons[i]=Ma[i-110];
  }
  Cons[120]=1.0*x[14];

  for (i=0; i<11; i++) {//transpose Cons into tCons for fortran convention in Lapack
     for (j=0; j<11; j++) {
         tCons[j*11+i]=Cons[i*11+j];
     }
  }
// THE PROBLEM STATEMENT (MATRIX FILLING) IN NOW TERMINATED
// NOW THE PROBLEM INSIDE THIS FUNCTION IS THE FOLLOWING:
// ASSURE THAT THE MATRIX Ba IS POSITIVE DEFINITE
// or: ALL THEIR EIGENVALUES ARE POSITIVE
// FOR THIS A LAPACK SUBROUTINE MUST BE CALLED
// 1/ DSYTRD: to factorize tCons into QTQ'
// 2/ DSTEQR: to compute EVal of T which are also EVal of tCons
 
   n12=12; // <>n!
   lda=12;
   info=1;
   work=0;
   myCompEigValLapack(n12,tCons,lda,Diag,info,work);

   for (i=0; i<11; i++) {
     g[i]=Diag[i]; 
   }
 // to complete: 12th constraint: g[11]=trace(Gamma_i) <= b
  //g[11]=x[7]+x[11];
}//end eval_g

bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values)
{//implementation
  bool status;
  Index i,j,k;
  Number eval_g_param0[12], eval_g_param1[12] ;
  Number var_params[15]={1E-3, 1E-3, 1E-3, 1E-3, 1E-3, 1E-3, 1E-3, 1E-3, 1E-3, 1E-3, 1E-3, 1E-3, 1E-3, 1E-3, 1E-3};
  Number param0[15], param1[15];

/* 1) if 'values' is NULL: return structure of the jacobian
   it is assumed in all cases to be full matrix 12x15 */ 
 
  if (values==NULL) {
     iRow[0]=0;
     jCol[0]=0;
     iRow[1]=0;
     jCol[1]=1;
     iRow[2]=0;
     jCol[2]=2;
     iRow[3]=0;
     jCol[3]=3;   
     iRow[4]=0;
     jCol[4]=4;
     iRow[5]=0;
     jCol[5]=5;
     iRow[6]=0;
     jCol[6]=6;
     iRow[7]=0;
     jCol[7]=7;
     iRow[8]=0;
     jCol[8]=8;
     iRow[9]=0;
     jCol[9]=9;
     iRow[10]=0;
     jCol[10]=10;
     iRow[11]=0;
     jCol[11]=11;
  // Now let's make a loop, at least!
    for (i=1; i<15; i++) {
        for (j=0; j<12; j++) {
            iRow[i*15+j]=i;
            jCol[i*15+j]=j;
        }
    } 
  } //if values==NULL
  else {
 /*2) if 'values' is not NULL: return the jacobian values
 hint: which can here only be approximated by finite difference (similar as the hessian with quasi newton) 
 **/ 
    for (i=0; i<15; i++) { 
       for (k=0; k<15; k++) {
         param0[k]=x[k];
         param1[k]=x[k];
       }
       param1[i] = param1[i] + var_params[i];  //small variation
  // INTERNAL CALL TO A FUNCTION eval_g !!! should normally work as it is defined in this same file nlp_ctr01.cpp::

 //      status=eval_g(n, (const Number*) param0, new_x, m, eval_g_param0); // eval constraints at param0 values (EVL computation, slow);
        status=eval_g(n, param0, new_x, m, eval_g_param0); // eval constraints at param0 values (EVL computation, slow);
       status=eval_g(n, param1, new_x, m, eval_g_param1); // eval constraints at param1 values 
      // finite difference for the all 12th components of the constraints: 
       for (j=0; j<12; j++) {
          values[i*12+j] = (eval_g_param1[j] - eval_g_param0[j]) /var_params[i]; 
       }
    }
  } // end if values==NULL

}

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   *
 bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values)  {

     // should not be implemented so that the Ipopt Solver automatically does a Quasi Newton Iteration
  }
******/
  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
				 const IpoptData* ip_data,
				 IpoptCalculatedQuantities* ip_cq) {
 // what to do?
}
