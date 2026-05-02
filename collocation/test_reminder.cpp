#include <iostream>
#include <stdlib.h>

using namespace std;

int main(int argc, const char **argv) {
   int row=3;
   struct { int quot; int rem; } dv;

   dv.quot = row;

   dv.rem = row % 2;
   std::cout << "result of division of row by two: " << dv.rem << "\n";
    
   return 0;
   }
