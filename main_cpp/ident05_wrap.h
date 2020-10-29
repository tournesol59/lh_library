#include <iostream>
#include <cstring>
#include <vector>
#include <list>

#include "../include/lhTypes.hpp"
#include "../include/ident05_coll.hpp"
#include "../include/ident05_fftdata.hpp"
#include "../recursive_LSQ/decl_relsq.hpp"

typedef struct inarg_s {
  Index order;

  Index dim;

  std::string filenameinp;
  
  std::string codenameinp;

  std::string lsqnameoutp;

} INARG_S;

//void wrapper_mexfunc_coll(Number *ptr_y, Number *ptr_t, INARG_S *inarg);

void wrapper_mexfunc_lsq(Number *ptr_u, Number *ptr_t, INARG_S *inarg);

