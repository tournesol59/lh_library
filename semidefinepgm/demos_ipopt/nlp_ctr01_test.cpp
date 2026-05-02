//#Include "/home/frederic/Ipopt-3.12.6/include/coin/IpIpoptApplication.hpp"
#include "IpIpoptApplication.hpp"
#include "nlp_ctr01.hpp"
//#include "IpTNLP.hpp"  //needed after some debug TBD
#include <iostream>

// defined in ..hpp
 
using namespace Ipopt;
  IndexStyleEnum the_index_style=TNLP::C_STYLE;
 
int main() {
  //vars to be displayed for tests
  Index the_n, the_m;
  Index the_nnz_jac_g, the_nnz_h_lag;
//  Index the_nele_jac;
//  Index the_iRow[15];
//  Index the_jCol[12];
//  Number the_x_l[15];
//  Number the_x_u[15];
//  Number the_x[15];
//  Number the_g_l[12];
//  Number the_g_u[12];
//  Number the_g[12];
//  Number the_z_l[15];
//  Number the_z_u[15];

  Number the_lambda[15];
//  Number obj_value;
//  Number grad_f[15];
//  Number values[12*15];
  
  bool status;
//  bool init_x, init_z, init_lambda;
//  bool new_x;

//  NLP_CTR01 mynlp=new NLP_CTR01();
  NLP_CTR01 mynlp;
  status=mynlp.get_nlp_info(&the_n, &the_m, &the_nnz_jac_g,
                             &the_nnz_h_lag, &the_index_style);
  std::cout << "the m: " << the_m << " the n: " << the_n << "\n";
//  printf("the m: %d, the n: %d \n", the_m, the_n);
//  status=mynlp.get_bounds_info(the_n, the_x_l, the_x_u,
//                              the_m, the_g_l, the_g_u);

//  status=mynlp.get_starting_point(the_n, init_x, the_x,
//                                   init_z, the_z_l, the_z_u,
//                                   the_m, init_lambda, 
//                                   the_lambda);

//  status=mynlp.eval_f(the_n, the_x, new_x, &obj_value);.

//  status=mynlp.eval_grad_f(the_n, the_x, new_x, grad_f);

//  status=mynlp.eval_g(the_n, the_x, new_x, the_m, the_g); 

//  status=mynlp.eval_jac_g(the_n, (const Number*) the_x, new_x,
//                         the_m, the_nele_jac, the_iRow, the_jCol, 
//                         values);


} 
