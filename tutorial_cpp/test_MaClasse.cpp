#include "./decl_MaClasse.hpp"
#include <cstring>
// idee: mettre en commentaire tous les objets de MaClasseDerivee dans decl_MaClasse.hpp sauf les methodes.

int main(int argc, char **argv) {
   std::vector<double> y(5,0.0);
   for (int i=0; i<5; i++) {
      std::cout << "Entrez valeur i= " << i << std::endl;
      std::cin >> y[i];
   }
   
   // demo of derived wih virtual member function
   MaClasseBase inst_maClBase(y,5,2);  // 5=size(y), 2=no params
   std::cout << "entre() and recalcule de la classe de base: " << std::endl; 
   inst_maClBase.recalcule();
   MaClasseDerivee inst_maClDer(y,5,2);
   std::cout << "normally entre() and recalcule de la classe derivee: " << std::endl;
   inst_maClDer.recalcule();

   // demo of derived contructor with call to base constructor
   std::vector<double> u(5,0.0);
   for (int i=0; i<5; i++) {
      std::cout << "Entrez valeur i= " << i << std::endl;
      std::cin >> u[i];
   }
   MaDeuxClasseDerivee inst_ma2ClD(y,u,5,2,1);
   inst_ma2ClD.recalcule();

   // demo of derived copy constructor
   MaDeuxClasseDerivee insttwo_ma2ClD= MaDeuxClasseDerivee(inst_ma2ClD);
   inst_ma2ClD.recalcule();

   // another derived class
   //MaTroisClasseDerivee instthree_ma3ClD(y,u,5,2,1,1);
   //instthree_ma3ClD.recalcule();

   // rules of pointers:
   MaClasseBase *pB, *pB0;
   MaClasseDerivee *pD1;
   MaDeuxClasseDerivee *pD2All, *pD2two;

   pB = &inst_maClBase;  // no inheritance, ok
   pD1 = &inst_maClDer; // same type, ok
   pD2two = &insttwo_ma2ClD;
   MaDeuxClasseDerivee *pD2= new MaDeuxClasseDerivee(y,u,5,2,1);; // same time as object instanciation ok
//   MaTroisClasseDerivee *pD3 = &instthree_ma3ClD;

   // conversion
   pB0 = (MaClasseBase *) pD2;  // target is basis class, ok
   //pD2All = (MaDeuxClasseDerivee *) pB;  // incorrect?or inc if D2 : virtual B
   //pD2All = dynamic_cast<MaDeuxClasseDerivee *>(pB); // correct vertic if D2 : B
   //std::cout << pB->newvaleurs[3] <<" : "<< pD2All->newvaleurs[3] << std::endl;
   //pD2All = dynamic_cast<MaDeuxClasseDerivee *>(pD2two); // correct, honrizont

   std::cout<< "Program terminated normally" << std::endl;
   return 0;
}