#include <iostream>
#include <cstring>
#include <vector>

/*---- la classe mere ----*/
class MaClasseBase {

   public:
     MaClasseBase(std::vector<double> entrees, int dsize, int dnp);
     ~MaClasseBase(void);
  //   MaClasseBase(const MaClasseBase &source);
  //   MaClasseBase &operator=(const MaClasseBase &source);

     virtual void entre(void)=0;
     void recalcule(void); // va appeler entre()
    // pass_iodata(std::vector<std::pair<double,double>> &list_yh, int select);
     std::vector<double> valeurs;
     std::vector<double> newvaleurs;
     
   protected:
     std::vector<double> coeffs;
     int size;
     int np;  
};

MaClasseBase::MaClasseBase(std::vector<double> entrees, int dsize, int dnp) :
	valeurs(std::vector<double>(dsize,0.)),
	newvaleurs(std::vector<double>(dsize,0.)),
	coeffs(std::vector<double>(dnp,1.)),
	np(dnp)
{
	coeffs[0]=1.0;
	for (int i=0; i<dsize; i++) {
		valeurs[i]=entrees[i];
	}
}

MaClasseBase::~MaClasseBase()
{
}

void MaClasseBase::entre(void) 
{
  double val;
  std::cout << "Entrez une valeur\n";
  std::cin >> val; 
  valeurs.push_back(val);
}

void MaClasseBase::recalcule(void) 
{
  entre();
  newvaleurs[0] = newvaleurs[0] + (valeurs.at(valeurs.size())-newvaleurs[0]) * coeffs[0];
  std::cout << "La valeur recalculee: " << newvaleurs[0] << "\n";
}

/*
MaClasseBase::pass_iodata(std::vector<std::pair<double,double>> &list_yh, int select) {

  int k, sel=1;
  double Ts=0.05;
  std::pair<double, double> singleton;
  for (k=0; k<size; k++) {
      singleton.first = (double (k))*Ts;
      switch (sel) {
         case 1:  singleton.second = newvaleurs[k];  // more case when complete
		  break;
      }
      std::cout << singleton.second << " ";
      list_yh.push_back( singleton );
  }
}
*/

/*---- la classe fille ----*/
class MaClasseDerivee : public MaClasseBase {

   public:
     MaClasseDerivee(std::vector<double> entrees, int dsize, int dnp);
     ~MaClasseDerivee(void);
  //   MaClasseDerivee(const MaClasseDerivee &source);
  //   MaClasseDerivee &operator=(const MaClasseDerivee &source);

     virtual void entre(void);
     void recalcule(void); // va appeler entre()
     //bool pass_iodata(std::vector<std::pair<double,double>> listy, int select);
     std::vector<double> valeurs;
     std::vector<double> newvaleurs;
     
   protected:
     std::vector<double> coeffs;
     int size;
     int np;  
     
     double biases;

};

MaClasseDerivee::MaClasseDerivee(std::vector<double> entrees, int dsize, int dnp) : MaClasseBase(entrees, dsize, dnp)
{
   biases = -0.5;
}

MaClasseDerivee::~MaClasseDerivee()
{
}

void MaClasseDerivee::entre(void) 
{
  double val;
  for (int i=1; i<2; i++) {
  std::cout << "Entrez une deuxieme valeur\n";
  std::cin >> val; 
  valeurs.push_back(val);
  }
  MaClasseBase::entre();

}

void MaClasseDerivee::recalcule(void) 
{
  entre();
  for (int i=1; i<size; i++) {
    newvaleurs[i] = newvaleurs[i-1] + (valeurs[i]-newvaleurs[i-1]) * coeffs[0];
  }
  std::cout << "La valeur recalculee: " << newvaleurs[1] << "\n";
}

int main() {
   std::vector<double> entrees;
   entrees.push_back(1.5);
   entrees.push_back(0.45);
   entrees.push_back(0.2);
//   MaClasseBase instA(entrees, 3, 2) ;
//   instA.entre();
   //   instA.recalcule();	

   MaClasseDerivee instB(entrees, 3, 2) ;
   instB.entre();
//   instB.recalcule();	
   return 0;
}
