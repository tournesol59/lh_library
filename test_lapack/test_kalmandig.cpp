/*
 * import Lapack
 */

#include "lapacke.h"
#include "../include/lhTypes.hpp"
#include <cstring>
#include <stdlib.h>
#include <iostream>
#include <istream>
#include <sstream>
#include <fstream>

using namespace lhlib;

int importDynamicProblem(int n, Number *Fi) //, Number *Gi, Number *Hi, Number *Qi, Number *Ri)
{
   //Fi = (Number*) malloc(n*n*sizeof(Nummber));
   Fi = new Number[n*n];
   
   char strFileName[14];
   const char *FileName ="inputtwo.txt";
   std::string headstr;
//   std::string valuestr;
   strcpy(strFileName, FileName);
   // open file
   std::ifstream kinp_file;
   kinp_file.open(strFileName,std::ios::in);
   
   int nmat=0;
   int irow, icol, i, j;
   bool finished=true;
   int line=0;
   int linemax=1*(n+1);
   Number value;
   
   while ((getline(kinp_file, headstr)) && (line < linemax)) {
      // while iterating each line of the file and counter 'nmat' and bool finished helps to situate oneself where we are in the file:
      if (finished == true) {
         finished == false;
         // head line
         std::stringstream iss(headstr);
         iss >> irow >> icol;
	 std::cout << "no rows: " << irow << ", no cols: " << icol << "\n";
	 i=0;
      }
      else if (nmat==0) {
     // import one row of Fi matrix
         if (i < irow) {
            std::stringstream iss(headstr);
	//    for (int j=0; j < icol; j++) {
        //       iss >> value;
	//       Fi[i*icol+j]=value;
	//    }
            iss >> Fi[i*icol+0] >> Fi[i*icol+1];
	    i++;
	 }
	 if (i == irow) {
           finished ==true;
	   nmat++;
	   i=0;
	 }
      }
      line++;
   }
   kinp_file.close();
   // control:
   std::cout << "imported Fi matrix: \n";
   for (int i=0; i < irow; i++) {
      for (int j=0; j < icol; j++) {
        std::cout << "  " << Fi[i*icol+j];
      }
      std::cout << " end \n";
   }
   std::cout << "\n";
   return 0;
}
   
int main(int argc, char** argv) {

   Number* Fi, Gi, Hi, Qi, Ri;
   int flag;

   flag = importDynamicProblem(2, Fi);

   return 0;
}
