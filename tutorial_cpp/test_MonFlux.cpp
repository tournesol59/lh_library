#include "decl_MonFlux.hpp"
#include <exception>

void myEx_handler(void) throw () {
  std::cout << "Error occurred .. doing clean up before termination" << std::endl;
  throw 1;
}

int main(int argc, char **argv) {

   // read argument for the name of the file that will be our input stream
   char filename[40]; //="InFile";
   strcpy(filename, argv[1]);  // prefer this, file can be easily changed
   std::cout << "first init ok: " << filename << std::endl;

   // string data to add at the end of the stream
   char DataToAdd[40];
   const char text[40]="Monday,27      February      XXXY";
   int taille=38;
   strncpy(DataToAdd, text, taille);
   DataToAdd[38]='\0';  // copy filename,since a char is input for null character manually added
   std::cout << "second step ok: " << DataToAdd << std::endl;
   
   // set exception handler, which is called if all defined exception are not detected
   std::set_unexpected( myEx_handler);
   std::cout << "third step ok: handler" << std::endl;
   
   // input stream treatment, under exception braces
   try {
      MyStreamFile strmfile=MyStreamFile(filename, taille, DataToAdd);
      std::cout << "fourth step: object created" << std::endl;
      strmfile.ReadFromFile();
      std::cout << "fifth step: file read" << std::endl;
   }
   catch (int i) {
      std::cout << "Exception of type int caught :" << i << std::endl;
   }
   catch (...) {
      std::cout << "Exception of some oyher type caught :" << std::endl;
   }
   
   return 0;
}
