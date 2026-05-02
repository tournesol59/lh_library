#include "open_decl.h"

MyClass::MyClass(const char* input, int ia, int ib) : 
   a(ia),
   b(ib)
{
   strncpy(inputFileName, input, 14);
}

MyClass::~MyClass() {}

int MyClass::add() {
   return (a+b);
}

