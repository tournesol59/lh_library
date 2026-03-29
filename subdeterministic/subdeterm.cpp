#include "subdeterm.hpp"
#include <malloc.h>

//useful procedure to test the linear QR-based inversion of a square matrix
double residual(const LaGenMatDouble &A, const LaGenMatDouble &X)
//                const LaGenMatDouble& b) // is equal to identity(N)
{
    int M = A.size(0);
    int N = A.size(1);
    LaVectorDouble xi(N);

    for (int i=0; i<N; i++) {
       for (int j=0; j<M; j++) { // extract a first col vector from X->xi
          xi(j+1)=X(j,i);
       }    
       LaVectorDouble Ei(M); Ei(i)=1.0;  // i-th canonical vector
       std::cout << "\tNorm_Inf(A*x-b)" << Norm_Inf(A * xi - Ei) << std::endl;
       std::cout << "\tNorm_Inf(A) " << Norm_Inf(A) << std::endl;
       std::cout << "\tNorm_Inf(x) " << Norm_Inf(xi) << std::endl;
       std::cout << "\tMacheps :" << Mach_eps_double() << std::endl;

      if (M > N)
      {
          LaVectorDouble Axb = A * xi - Ei;
          LaVectorDouble R(N);

          Blas_Mat_Trans_Vec_Mult(A, Axb, R);
          return Norm_Inf(R) /
               (Norm_Inf(A) * Norm_Inf(xi) * N * Mach_eps_double());

      }
      else
      {
          return Norm_Inf(A * xi - Ei ) /
               ( Norm_Inf(A) * Norm_Inf(xi) * N * Mach_eps_double());
      }
    }
}

// constructor
IDENT05_DETERMSUB::IDENT05_DETERMSUB(Index start, Index end, Index nmeas, Number *Ym, Number *Um)
{
   Index ir = (end-start)/2;  // --a verification here would be fine
   if ((ir*2) != (end-start)) {
	   std::cout << "difference between end and start should be div by two\n";
   }
   else {
     istart=start;
     irange=ir;
     nmeasurement=nmeas;

     Ypast = new Number[ir*nmeasurement];
     Upast = new Number[ir*nmeasurement];
     Yfuture = new Number[ir*nmeasurement];
     Ufuture = new Number[ir*nmeasurement];
     Y_on_W = new Number[ir*nmeasurement];
     Y_on_W_c = new Number[ir*nmeasurement];

     for (Index i=0; i<irange; i++) {
         for (Index j=0; j<nmeasurement; j++) {
            Y_on_W[i*nmeasurement+j]= 0.0; // past-part the Hankel matrix is not square
            Y_on_W_c[i*nmeasurement+j]=0.0;
	 }
     }

   }
}

//destructor
 IDENT05_DETERMSUB::~IDENT05_DETERMSUB()
{
   
  free(Ypast);
  free(Upast);
  free(Yfuture);
  free(Ufuture);
  free(Y_on_W);
  free(Y_on_W_c);

}

// construct matrices
bool IDENT05_DETERMSUB::WriteData2Hankel(Number *Y, Number *U) {
// ROW major matrix format

     for (Index i=0; i<2*irange; i++) {  // 2*irange hence 'past' and 'future' rows
         for (Index j=0; j<nmeasurement; j++) {
            Ypast[i*nmeasurement+j] = Y[istart+i*nmeasurement+j]; // past-part the Hankel matrix is not square
            Upast[i*nmeasurement+j] = U[istart+i*nmeasurement+j];
	 }
         for (Index j=0; j<nmeasurement; j++) {
            Yfuture[i*nmeasurement+j] = Y[istart+(irange+i)*nmeasurement+j]; // past-part the Hankel matrix is not square
            Ufuture[i*nmeasurement+j] = U[istart+(irange+i)*nmeasurement+j];
         }
    }
 return 0;
}


//----------------------------------
//class for the projective method
 IDENT05_PROJECTIVE::IDENT05_PROJECTIVE(Index dy, Index nc,Index ny_1, Index ny)
{
    dimy=dy;
    dimcol=nc;
    dimrow_y_p=ny_1;
    dimrow_y_f=ny;
}

 IDENT05_PROJECTIVE::~IDENT05_PROJECTIVE()
{

}


 bool IDENT05_PROJECTIVE::LUprojectY_U(Number *Ypast, Number *Upast, Number *Yfuture, Number *Ufuture, Number *Res, Number *L, Number *Q)
{
// Number *Wp = new Number[2*dimcol*dimrow_y_p*dimy];
  int N=2*dimrow_y_p*dimy; // N1=N2/2; N2=0; N3=N/2;  // TBD
  int M=dimcol; // ok

  // lapack++ library objects creation: create the transpose of [[Up][Yp]]
  // since upper QR factorisation is available but lower LQ is needed

  LaGenMatDouble AR(N,M);

  for (int i=1; i<N/2+1; i++) {      //first part of columns: Upast
	 // Fortran index based (start by one)
     for (int j=1; j<M; j++) {
        AR(j,i)=Upast[(i-1)*M+j-1];    // but C arrays start by zero
     }
  }
  for (int i=N/2+1; i<(N+1); i++) {  // second part of colums: Ypast
     for (int j=1; j<M; j++) {
        AR(j,i)=Ypast[(i-N/2-1)*M+j-1];
     }
  }
  /* for (int i=N+1; i<(3*N/2+1); i++) { // TBD third part of columns
   *   for (int j=1; j<M; j++) {
   *      AR(j,i)=Upast[(i-N-1)*M+j-1];
   *   }
   * }
   */
  for (int i=N+1; i<(3*N/2+1); i++) {  // TBD, 3*N/2+1, fourth part of colums: Yfuture
     for (int j=1; j<M; j++) {
        AR(j,i)=Yfuture[(i-N-1)*M+j-1];
     }
  }

 // create the QRFactorization class instance (null object):
   LaGenQRFactDouble QR;

     /** Ref liblapack++: Calculate the QR decomposition of A.
     * in "gfqrd.h"
     * This is in-place, i.e. it destroys the input matrix A and
     * keeps a reference to its memory around. In other words, you
     * cannot do anything with your input matrix A anymore. You can
     * safely delete any references to A because this object will
     * keep its own references still around.
     */
  QR.decomposeQR_IP(AR);     

  // copy results with transpose, careful: AR is equal to R now M*(3N/2)
  // but only the (3N/2) first rows of AR make sense
  LaGenMatDouble L11(N/2,N/2);
  LaGenMatDouble L21(N/2,N/2);
  LaGenMatDouble L22(N/2,N/2);
  LaGenMatDouble L31(N/2,N/2);
  LaGenMatDouble L32(N/2,N/2);
  LaGenMatDouble L33(N/2,N/2);

  for (int i=1; i<N/2+1; i++) {
     for (int j=1; j<N/2+1; j++) {
         L11(i,j)=AR(j,i);       // L11=tr(R11), TBD (1..N,1..N)
     }
  }
  for (int i=N/2+1; i<N+1; i++) {
     for (int j=1; j<N/2+1; j++) {
         L21(i-N/2,j)=AR(j,i);       // L21=tr(R12), TBD(N+1..3N/2,1..N)
     }
  }
  for (int i=N+1; i<3*N/2+1; i++) {
     for (int j=1; j<N/2+1; j++) {
         L31(i-N,j)=AR(j,i);       // L31=tr(R13), TBD (3N/2+1..2N, 1..N)
     }
  }
  for (int i=N/2+1; i<N+1; i++) {
     for (int j=N/2+1; j<N+1; j++) {
         L22(i-N/2,j-N/2)=AR(j,i);    // L22=tr(R22), TBD (N+1..3N/2,N+1..3N/2)
     }
  }
  for (int i=N+1; i<3*N/2+1; i++) {
     for (int j=N/2+1; j<N+1; j++) {
         L32(i-N,j-N/2)=AR(j,i);     // L32=tr(R23), TBD (3N/2+1..2N,N+1..3N/2)
     }
  }
  for (int i=N+1; i<3*N/2+1; i++) {
     for (int j=N+1; j<3*N/2+1; j++) {
         L33(i-N,j-N)=AR(j,i);       // L33=tr(R33), TBD (3N/2+1..2N,3N/2+1..N)s
     }
  }

  // generate Q:
  QR.generateQ(AR);
  LaGenMatDouble Q1(N/2,N/2);
  LaGenMatDouble Q2(N/2,N/2);
  LaGenMatDouble Q3(N/2,N/2);

  for (int i=1; i<N/2+1; i++) {
     for (int j=1; j<N/2+1; j++) {
         Q1(i,j)=AR(j,i);       // Q1j=tr(Q1i)
     }
 }

  for (int i=N/2+1; i<N+1; i++) {
     for (int j=1; j<N/2+1; j++) {
         Q2(i-N/2,j)=AR(j,i);       // Q2j=tr(Q2i)
     }
 }

  for (int i=N+1; i<3*N/2+1; i++) {
     for (int j=1; j<N/2+1; j++) {
         Q3(i-N,j)=AR(j,i);       // Q3j=tr(Q3i)
     }
 }
/*
 * preview a fourth Q matrix)
 */

  // complete the post processing to output projection part
  LaGenMatDouble Yf_Yp(N/2,N/2);  // previewed for Result
  LaGenMatDouble Yf_Up(N/2,N/2);  // previewed for Result
  LaGenMatDouble IN2(N/2,N/2); // previewed for identity matrix
  IN2.eye(N/2);  
  LaGenMatDouble L11_1(N/2,N/2) ; // previewed to inverse (L31-L32*L21)*(L11)^(-1) column by column
  LaGenMatDouble L22_1(N/2,N/2) ; // previewed to inverse L32*(L22)^(-1)

  LaVectorDouble LXi(N/2);  // previewed for each (X*L22=Ei) column vector system
  LaVectorDouble LDXi(N/2);  // previewed for each (X*L11=Ei) column vector system

  
  LaGenMatDouble transL21(N/2,N/2);   // used to transpose
  Blas_Mat_Trans_Mat_Mult(L21, IN2, transL21);  // trans21=L21'*IN2=L21'
  LaGenMatDouble L31_L32cL21(N/2,N/2);
  L31_L32cL21 = Blas_Mat_Mat_Mult(L32,L21);
  L31_L32cL21 = L31_L32cL21 = L31 - L31_L32cL21; // test if it is correct

  
  for (int i=1; i<N/2; i++) {  // loop over the columns of tr(L22) and tr(L11) as an inversion must be done from right and only VectorProcedure linear solving is available
     // firstly, each transpose of column (each row) of L31_L32cL21 and Q1:
    LaIndex LDXInd(1,N/2);
    LaVectorDouble Ei(N/2);  // previewed for right hand side
    Ei= L31_L32cL21( 1, LDXInd );
    Blas_Mat_Trans_Mat_Mult(Q1, IN2, L11_1) ; // trans error
    // secondly, perform sys solve for this column (tr):
    std::cout << "LUprojectY_U" << ": LaQRLinearSolve: Matrix A=" << L31
                  << "  Right hand side B=" << Ei;
    LaQRLinearSolveIP(L11_1, LDXi, Ei);
    std::cout << "  Solution X=" << LDXi;

    //   thirdly copy the result in Yf_Up:
    Yf_Up(i,LDXInd) = LDXi;

    // lastly but not least Yf_Up has to be transposed 

 //   //we could also have used  LaLULinearSolve(L22,Yf_inv,eye(N/2)) 
 //   ?=Blas_Mat_Mat_Mult(L32,Q2,false,false); //simple
 //
 //   Then repeat same procedure for L21*inv(L22)
    LaIndex LXInd(1,N/2);
    LaVectorDouble Ei2(N/2);
    Ei2=L21(1,LXInd);
    Blas_Mat_Trans_Mat_Mult(L22, IN2, L22_1) ; // error trans L22_1=L22
    // secondly, perform sys solve for this column (tr):
    std::cout << "LUprojectY_U" << ": LaQRLinearSolve: Matrix A=" << L22_1
                  << "  Right hand side B=" << Ei2;
    LaQRLinearSolveIP(L22_1, LXi, Ei2);
    std::cout << "  Solution X=" << LXi;
    // thirdly copy the result in Yf_Yp:
    Yf_Yp(i,LXInd) = LXi;

    // lastly but not least Yf_Yp has to be transposed 
 
 //   remind that Fortran lapack procedure call:
 // void  Blas_Mat_Mat_Mult (const LaGenMatDouble &A, const LaGenMatDouble &B, LaGenMatDouble &C, bool transpose_A, bool transpose_B=false, double alpha=1.0, double beta=0.0) 
 
  } // end for each column

 //      //    std::cout << "Residual " << residual(L31_L32cL21, Yf_Up, IN2) << std::endl;
 //
 // deallocate
 // free(Wp);

 return 0;
}


