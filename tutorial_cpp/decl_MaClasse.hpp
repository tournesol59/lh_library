#include <iostream>
#include <vector>
#include "../include/lhTypes.hpp"

using namespace lhlib;

// demonstration d'une fonction virtuelle: une classe de base plus deux classes
// derivee public.
// Nota: il serait utile de rajoutee une classe derivee supplementaire (meme si pas reprise dans a suite) derivee de MaClasseDerivee qui elle serait virtuelle
// pour tester les dynamic_cast<>

class MaClasseBase {

   public:
     MaClasseBase(std::vector<double> entrees, int dsize, int dnp);
     ~MaClasseBase(void);
  //   MaClasseBase(const MaClasseBase &source);
  //   MaClasseBase &operator=(const MaClasseBase &source);

     virtual bool entre(void);
     bool recalcule(void);
     bool pass_iodata(std::vector<std::pair<double,double>> listy, int select);
     
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
 //    MaClasseDerivee(const MaClasseDerivee &source);
 //    MaClasseDerivee &operator=(const MaClasseDerivee &source);
     
     bool entre(void); // do not declare as it is identical to base class
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

/*
class MaDeuxClasseDerivee : public MaClasseDerivee {
  // prend un signal de commande supplementaire
   public:
     MaDeuxClasseDerivee(std::vector<double> entrees, std::vector<double> com, int dsize, int dna, int dnb);	   
    // MaDeuxClasseDerivee(std::vector<double> entrees, std::vector<double> com, int dsize, int dna, int dnb) ;
     MaDeuxClasseDerivee(const MaDeuxClasseDerivee &source);
     MaDeuxClasseDerivee &operator=(const MaDeuxClasseDerivee &source);

     ~MaDeuxClasseDerivee(void);

     bool entre(void);  // from base class
     bool recalcule(void);
     std::vector<double> commande;

//     std::vector<double> valeurs;
//     std::vector<double> newvaleurs;

   protected:
     
//     std::vector<double> coeffs;
//     int size;
//     int np;
     
     std::vector<double> states;  //same length as coeff, but it is only a demo
     int na;
     int nb;
};


class MaTroisClasseDerivee : public MaDeuxClasseDerivee {
  // prend un delay sur la commande supplementaire
   public:
     MaTroisClasseDerivee(std::vector<double> entrees, std::vector<double> com, int dsize, int dna, int dnb, int idelay);	   
     MaTroisClasseDerivee(const MaTroisClasseDerivee &);
     MaTroisClasseDerivee &operator=(const MaTroisClasseDerivee &);

     ~MaTroisClasseDerivee(void);

     bool entre(void);  // from base class
     bool recalcule(void);
     std::vector<double> commande;

//     std::vector<double> valeurs;
//     std::vector<double> newvaleurs;

   protected:
     
//     std::vector<double> coeffs;
//     int size;
//     int np;

     std::vector<double> states;  //same length as coeff, but it is only a demo
     int na;
     int nb;
     int delay;
};
*/


