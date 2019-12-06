#include "decl_reslq.h"

int generate_rdom_example( std::vector<double> &sig, int n, double ca, double cb);

int main() {
  int k;
  std::vector<double> data={};

  generate_rdom_example(&data, 10, 3.0, 1.0);

//  IDENT05_RELSQ slqClassInst=IDENT05_RELSQ(10,1,1.0);

  return 0;
}
