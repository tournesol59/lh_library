#include "decl_iogenerate.hpp"

int generate_more_from_file(std::vector<double> &sig, int &Npty, const char * filename, int n, int m);
int generate_more_from_file_second(std::vector<double> &sig, int dsize, const char* fileName, int m);

IDENT05_IOGENERATE::IDENT05_IOGENERATE(const char *sdirName) {
   std::strcpy(dirName, sdirName);
}

IDENT05_IOGENERATE::~IDENT05_IOGENERATE() 
{
}

int IDENT05_IOGENERATE::choiceInterface() {
  // basic for the moment, select two files that must be present in dirName
   std::cout << "Importez deux fichiers \n";

   std::cout << "Veuillez entrer un nom de fichier a importer \n";
   std::cin >> rawDataName1;
   std::cout << "Veuillez entre le nombre de ligne du fichier \n";
   std::cin >> dsize1;
   
   std::cout << "Veuillez entrer un nom de fichier a importer \n";
   std::cin >> rawDataName2;
   std::cout << "Veuillez entre le nombre de ligne du fichier \n";
   std::cin >> dsize2;

   return 0;
}


int IDENT05_IOGENERATE::loadvectors(std::vector<double> &sig, int m) {
    // TBC
    // shall use the method generate_more_from_file_second(vector<double &sig,...) in the file ../recursive_LSQ/generate.cpp
   if (m==1) {
      generate_more_from_file_second(sig, dsize1, (const char*) rawDataName1, 1);
   }
   else if (m==2) {
      generate_more_from_file_second(sig, dsize2, (const char*) rawDataName2, 1);
   }
   return 1;
}
