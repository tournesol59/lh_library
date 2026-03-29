#include <iostream>

class IDENT05_RELSQ {

  public:
    bool IDENT05_RELSQ(int dsize, int dn double fTs);
    bool ~IDENT05_RELSQ();
    bool recursive_algorithm(int k, double &y[], double &u[],
		             double varian,
                             double y_h_1, double y_hd_1,
                             double &y_h[], double &y_hd[], int order,
			     double &pe_y[], double &pe_yd[]);
  private:
    int size;
    int n;
    double Ts;

};
