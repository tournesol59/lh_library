#include <iostream>
#include <vector>
#include "../include/lhTypes.hpp"

using namespace lhlib;

// demonstration d'une fonction virtuelle

class MaClasseBase {

   public:
     MaClasseBase(std::vector<double> entrees, int dsize, int dnp);
     ~MaClasseBase(void);
     bool entre(void);
     bool recalcule(void);
     // virtual bool recalcule(void);
     std::vector<double> valeurs;
     std::vector<double> newvaleurs;

   protected:
     std::vector<double> coeffs;
     int size;
     int np;  
};

class MaClasseDerivee : public MaClasseBase {

   public:
     MaClasseDerivee(std::vector<double> entrees, int dsize, int dnp);
     ~MaClasseDerivee(void);
     
     // bool entre(void); // do not declare as it is identical to base class
     bool recalcule(void);
/*
     std::vector<double> valeurs;
     std::vector<double> newvaleurs;
*/
   protected:
/*     
     std::vector<double> coeffs;
     int size;
     int np;
*/
};

class MaDeuxClasseDerivee : public MaClasseBase {
  // prend un signal de commande supplementaire
   public:
     MaDeuxClassDerivee(std::vector<double> entrees, std::vector<double> com, int dsize, int dna, int dnb);
     ~MaDeuxClasseDerivee(void);

     // bool entre(voi);  // from base class
     bool recalcule(void);
     std::vector<double> commande;
/*
     std::vector<double> valeurs;
     std::vector<double> newvaleurs;
*/
   protected:
/*     
     std::vector<double> coeffs;
     int size;
     int np;
*/
     
     std::vector<double> states;  //same length as coeff, but it is only a demo
     int na;
     int nb;
};

