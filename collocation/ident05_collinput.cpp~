/* Scope: parse the first of the two input file
 * LIEBHERR TOULOUSE
 *******************************************************/
#include "../include/ident05_coll.hpp"

#ifndef __TEST_COLL_ONLY__
#define __TEST_COLL_ONLY__
#endif
//using namespace std;
using namespace lhlib;

bool IDENT05_COLL::read_parse_file() 
{
   Index row,col,row_guess;

// Open the input file, which is a text file with 6 columns
  std::ifstream lh_file;
  lh_file.open("TheFile", std::ifstream::in);  // test w/o variable Name
  
  row=0;
  row_guess=0; // prescaler output of row % num_points
  while ((!lh_file.eof()) && (row<41)) {
  /**
  * complete here
  **/
         std::string line;
	 std::getline(lh_file, line);
         std::stringstream ss(line);
         col = 0;	 
	 while (ss >> dataarray[row][col]) col++; // fait la detection automatique chaine vide pour s'arreter
#ifdef __TEST_COLL_ONLY__
	 std::cout << "tf= " << dataarray[row][0] << " ";
	 std::cout << "utf= " << dataarray[row][1] << " ";
	 std::cout << "ytf= " << dataarray[row][2] << " ";
	 std::cout << "c2f= " << dataarray[row][3] << " ";
	 std::cout << "c1f= " << dataarray[row][4] << " ";
	 std::cout << "c0f= " << dataarray[row][5] << " ";
	 
#endif
         if (row==row_guess*(num_points)-1) {
             yguessed[row_guess][0]=dataarray[row][0];
             yguessed[row_guess][1]=dataarray[row][2];
	     row_guess++;
         }
	 row++;
   }
   num_rows=row;
   lh_file.close();

   return 0;
}
