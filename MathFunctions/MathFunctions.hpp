extern "C" double mysqrt(double x);
extern "C" double myComputeScalarProduct(int len, double * t, double * f, double * g);
extern "C" int mySolveLinearLapack(int n, int rhs, double * A, int lda, int * ipiv, double * B, int ldb, int info);
