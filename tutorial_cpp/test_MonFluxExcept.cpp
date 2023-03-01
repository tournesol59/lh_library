#include "decl_MonFlux.hpp"
#include <exception>

void myEx_handler(void) throw () {
   std::cout << "Error occurred .. doing clean up before termination" << std::endl;
  throw 1;
}

// demo of a file openin (as least) but with throwing exception
int fdemo_ex(const char* fname) {
   char filename[10];
   strcpy(filename, fname);
   std::fstream infile;
   int res;
   try {
      infile.open(filename, std::fstream::in | std::fstream::out);
      if ( (infile.rdstate() & std::ifstream::failbit ) != 0 ) {
        // std::cerr << "Error opening file\n";
	MyFileNameExcept except=MyFileNameExcept(1);
	throw except;  // some amout of code, but this idea is to demonstrate this
      }
      else {
        std::cout << "File normally opened" << std::endl;
      }
      infile.close();
   }
   catch (int e) {
      res=e;
      std::cout << "An error exception of type int: " << e << " was thrown " << std::endl;
   }
   catch (MyFileNameExcept except) {
      res=except.getcode();
      except.printex();
   }
   catch (...) {
      res=-1;
      std::cout << "An other type of exception was thrown" << std::endl;
   }
   return res;
}

/**
 * main
 */

int main(int argc, char **argv) {
   // set exception handler, which is called if all defined exception are not detected
   std::set_unexpected( myEx_handler);
   int res;
   try {	   
     res= fdemo_ex((const char*) argv[1]);
   }
   catch (MyException except) {
      std::cout << "In main from fdemo(), code=" << res << "   ";
      except.printex();
   }
   return res;
}
