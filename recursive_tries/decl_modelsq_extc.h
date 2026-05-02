#include "decl_modelsq.h"
//une autre classe: celle en cpp aussi mais n utilisant que des objets compatibles C. Attention maintenant c est la classe vectorop qui est utilisee comme vecteur

/* ------- class IDENT05_MODELSQ_EXTC ------*/

  class IDENT05_MODELSQ_EXTC : public IDENT05_MODELSQ {
     public:
   IDENT05_MODELSQ_EXTC(vectorop yvalues, int dsize, int dn);
   IDENT05_MODELSQ_EXTC(const IDENT05_MODELSQ_EXTC &source);
   IDENT05_MODELSQ_EXTC &operator=(const IDENT05_MODELSQ &source);
   ~IDENT05_MODELSQ_EXTC();
   // ajouter une methode read 
   void setphivalues(int ind);
   void predict(int ind, double &epsilon);
   void innovation(int ind, double epsilon);  
   void update(int ind);
   bool pass_iodata(std::vector<std::pair<double,double>> &list_yh, std::string str_data); // one must redeclare it

     private:
   int size; // size of vector y
   int n; // size of vector coeffs and unique column-matrixK
   vectorop y;
   vectorop phi;
   vectorop coeffs;
   vectorop matrixK;
   vectorop matrixF;
};

