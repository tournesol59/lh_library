#include "decl_MonContainer.hpp"
#include <string.h>
#include <cstring>
#include <iostream>
#include <sstream>

int main(int argc, char **argv) {
   const int taille=10;
   char chainecible[10];
   const char *chaine;
   std::string maligne;

   std::cout<< "Entrez une chaine de "<< taille <<" 10 caracteres" << std::endl;
   std::cin >> maligne;
   strncpy(chainecible, maligne.c_str(), taille);
   std::cout << "La chaine copiee dans cible : " << chainecible;
   MonObjetDeriv instMonObj1; // = MonObjetDeriv(taille, chainecible);

   std::cout<< "Entrez une chaine de "<< taille <<" 10 caracteres" << std::endl;
   std::cin >> maligne;
   strncpy(chainecible, maligne.c_str(), taille);
   std::cout << "La chaine copiee dans cible : " << chainecible;
   MonObjetDeriv instMonObj2; // = MonObjetDeriv(taille, chainecible);


   MonContainer instMonCont = MonContainer();
   instMonCont.add( instMonObj1 );
   instMonCont.add( instMonObj2 );
   instMonCont.print(); // print the list of Objects in the "Bag" (container) 
   instMonCont.remove( instMonObj2 );
   instMonCont.print();

   return 0;
}
