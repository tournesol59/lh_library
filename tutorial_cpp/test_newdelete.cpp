#include <iostream>

int main() {

   int* montable;
   montable=new int[3];
   montable[0]=1;
   std::cout << "premier element: " << montable[0] << "\n";
   int*p=montable;
   std::cout << "premier element: " << p[0] << "\n";
   delete montable; 
   std::cout << "premier element: " << p[0] << "\n";

// p n’est plus valide
  return 0;
}
