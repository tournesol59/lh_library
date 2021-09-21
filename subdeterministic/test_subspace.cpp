#include "./subdeterm.hpp"

#include "lafnames.h"       /* macros for LAPACK++ filenames */
#include LA_GEN_MAT_DOUBLE_H
#include LA_VECTOR_DOUBLE_H
#include LA_SYMM_MAT_DOUBLE_H
#include "blaspp.h"
#include LA_SOLVE_DOUBLE_H
#include LA_GENERATE_MAT_DOUBLE_H
#include LA_EXCEPTION_H
#include LA_UTIL_H

const double val_Y[9] = {0.1, 0.4, 1.15, 1.45, 1.9, 2.65, 3.0, 3.35, 3.95};
const double val_U[9] = {0.9, 1.15, 1.0, 0.95, 0.8, 1.025, 1.1, 0.95, 1.0};

double eig_residual(const LaSymmMatDouble &A, double lambda,
                    const LaVectorDouble &x)
{
    int N = A.size(0);

    return Norm_Inf(A * x - lambda * x) /
           (Norm_Inf(A) * Norm_Inf(x) * N * Mach_eps_double());


}

void TestGenEigSolve(int N)
{
//    LaSymmMatDouble A(N, N);
      LaGenMatDouble  A(N,N), Bp;
      LaVectorDouble  v(N);      
    if (N==3) {
      LaIndex  A_Ind(0,2,1);
      Bp.ref(A(A_Ind,A_Ind)); 
      Bp(0,0)=val_Y[0];   Bp(0,1)=val_Y[1];   Bp(0,2)=val_Y[2];
      Bp(1,0)=val_Y[1];   Bp(1,1)=val_Y[2];   Bp(1,2)=val_Y[3];
      Bp(2,0)=val_Y[2];   Bp(2,1)=val_Y[3];   Bp(2,2)=val_Y[4];
      v(0)=1.0;
      v(1)=0.0;
      v(2)=0.0;
    }
    else {
	    std::cout << "TestEigSolve: nothing done if N <> 3\n";
    }

#ifndef HPPA
    const char fname[] = "TestGenEigSolve() ";
#else
    char *fname = NULL;
#endif

    //char e = 'e';

    LaGenerateMatDouble(A);

    LaGenMatDouble Eigenvectors(N, N);

    std::cerr << fname <<
              ": testing LaEigSolve(LaSymmMat, eig_value, eig_vectors) \n";

    LaEigSolve(A, v, Eigenvectors);
// see laslv.h


    for (int i = 0; i < A.size(0); i++)
    {
        LaIndex I(0, N - 1);
        double res = eig_residual(A, v(i), Eigenvectors(I, i));
        if (res > 1)
        {
            std::cerr << fname << " residual " <<  res << " is too high.\n";
            exit(1);
        }
    }

    // if we've made this far, all eigenvalue/vector pairs have
    // been tested.

    std::cerr << fname << ": LaEigSolve() success.\n\n";

}


int main(int argc, char **argv)
{

    std::cout.precision(4);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);

    if (argc < 2)
    {
        std::cerr << "Usage " << argv[0] << " N " << std::endl;
        exit(1);
    }
    int N = atoi(argv[1]);

    std::cout << "Testing " << N << " x " << N << " system." << std::endl;

    TestGenEigSolve(N);

}

