#include <iostream>
#include "collocation/ident05_coll.hpp"

typedef struct inarg_s {
  Index order;

  std::string filenameinp;
  
  std::string codenameinp;

} INARG_S;

void  wrapper_mexfunc_coll(Number *ptr_y, Number *ptr_t, INARG_S *inarg);

