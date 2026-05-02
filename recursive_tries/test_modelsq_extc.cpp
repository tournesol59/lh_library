#include "decl_modelsq_extc.h"
#include "../include/ident05_fftdata.hpp" 
#include "../recursive_project/decl_relsq_copy.hpp"

// puis le main:
int main(int argc, char **argv) {
// variables au debut: y_data, inputfilename

    char inFileName[14];
    int Npty=30;
    std::vector<double> y_data;
    vectorop y_dataop(Npty);
    std::vector<double> u_data=std::vector<double>(Npty, 0.0);
    vectorop u_dataop(Npty);

    std::vector<std::pair<double,double>> lsqoutlist;  // exchange list of doubles // check TBD
    std::string str_data="yst0";
    generate_from_file(y_data, Npty, inFileName);
 // recopy
    for (int i=0; i<Npty; i++) {
       y_dataop.values[i]=y_data.at(i);
       u_dataop.values[i]=u_data.at(i);
    }
// puis instances de classe
    double fvarian=0.5;
    double fTs=0.1;
    int order=2;
    
    IDENT05_MODELSQ_EXTC instlsqExtC(y_dataop, Npty, order);
// puis appels de methodes
/* afterwards
    double epsilon;
    instlsqExtC.predict(2, epsilon);
    instlsqExtC.innovation(2, epsilon);
    instlsqExtC.update(2);
    */
// export

      //create another export class
    char explsqFileName[14];
    strncpy(explsqFileName, "a.lsqo", 7);
    IDENT05_IODATA explsqClassInst=IDENT05_IODATA(Npty, 0.1, explsqFileName);
    instlsqExtC.pass_iodata(lsqoutlist, str_data);
  
  // recalculated ydata
    explsqClassInst.exportToDisk(lsqoutlist);
}

