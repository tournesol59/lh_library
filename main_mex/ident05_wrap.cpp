#include "ident05_wrap.h"

void wrapper_mexfunc_coll(Number *ptr_y, Number *ptr_t, INARG_S *inarg) {
  char strFileName[14];
  char strCodeName[14];
  int lenFileName=14;
  int i;

  strncpy(strFileName, inarg->filenameinp, 8); 
  strncpy(strCodeName, inarg->codenameinp, 8); 
  std::cout << "Have string " << strFileName << " of length " << lenFileName << " as arg[1] and " << strCodeName << " as arg[2]\n";

  // call to constructor
  IDENT05_COLL cnClassInst=IDENT05_COLL(6,strFileName, strCodeName);

  cnClassInst.read_parse_file();
  cnClassInst.read_parse_code();

// If Testing: Decompose:
  cnClassInst.ExpandSeriesLinearSys_ref1();
  cnClassInst.SolveSeriesLinearSys_ref1();

   for (i=1; i<cnClassInst.num_points; i++) {
       ptr_y[1] = solarray[i][0];
       ptr_y[i] = solarray[i][1];
   }
// Else:
//  cnClassInst.SolveNumRangesSys_ref1();
} 
