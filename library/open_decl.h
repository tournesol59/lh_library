#include <iostream>

class MyClass {

   public:
      MyClass(const char* input, int ia, int ib);
      ~MyClass();
      int add();
   protected:
      int a,b,c;
      char inputFileName[14];
};

