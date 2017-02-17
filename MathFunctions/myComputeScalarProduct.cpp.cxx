#include "MathFunctions.h"
#include <math.h>

// a hack square root calculation using simple operations
double myComputeScalarProduct(int len, double* t, double* f, double* g)
{

  double result;
  result = 0.0;

  // do len iterations
  int i;
  for (i = 0; i < len-1; ++i) {
	result+= ((f[i]*g[i])+(f[i+1]*g[i+1]))/2*(t[i+1]-t[i]);
  }

  return result;
}
