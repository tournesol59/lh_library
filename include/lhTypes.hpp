// allow a different namespace as Ipopt
#ifndef __LHTYPES_HPP__
#define __LHTYPES_HPP__

namespace lhlib
{
  /** Type of all numbers */
  typedef double Number;
  /** Type of all indices of vectors, matrices etc */
  typedef int Index;
  /** Type of default integer */
  typedef int Int;

  typedef struct dpair {
     Number x; 
     Number y;
  }* p_dpair;

} // namespace lhlib

#endif
