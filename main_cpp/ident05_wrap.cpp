#include "ident05_wrap.h"

/************** COLLOCATION ***************
 *
void wrapper_mexfunc_coll(Number *ptr_y, Number *ptr_t, INARG_S *inarg) {
  char strFileName[14];
  char strCodeName[14];
  int lenFileName=14;
  int i;

  strncpy(strFileName, inarg->filenameinp, 8); 
  strncpy(strCodeName, inarg->codenameinp, 8); 
//  std::cout << "Have string " << strFileName << " of length " << lenFileName << " as arg[1] and " << strCodeName << " as arg[2]\n";

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
************ END COLLOCATION ****************/

void wrapper_mexfunc_lsq(Number *ptr_u, Number *ptr_t, INARG_S *inarg) {

  int i;
  int matlab_emb=1;
  std::vector<double> data={};
  std::vector<double> dtime={};
  lhlib::dpair singleton;
  li_doubles datalist;  // exchange list of doubles
  

  // generate 30 points over line 0.8+0.533*t added with randoms numbers
  // between -0.5 and 0.5
  if (matlab_emb) {
     for (i=1; i<inarg->dim; i++) {
        data.push_back(ptr_u[i]);
        dtime.push_back(ptr_t[i]);
     }
  }
  else {
     for (i=0; i<30; i++) {
        dtime.push_back( ((double) i)/10.0 );
     }
       generate_rdom_example(data, 30, 0.8, 0.533);
  }
  // export class
  char expFileName[14];
  strncpy(expFileName, ".randin", 8);
  IDENT05_IODATA expClassInst=IDENT05_IODATA(inarg->dim, expFileName);

  // copy data into list but not found an optimal iterator
  for (i=0; i< inarg->dim; i++) {
     singleton.x = dtime.at(i);
     singleton.y = data.at(i); 
     datalist.push_back( singleton );
  }
  i=1;
  // so not optimal..
  //
  expClassInst.exportToDisk(datalist);

  // create a instance of class which manages to approx "data"
  // with size 30 points that are given, order is 1, sample Time is 1.0
  IDENT05_RELSQ lsqClassInst=IDENT05_RELSQ(30,1,1.0);


// Else
}
